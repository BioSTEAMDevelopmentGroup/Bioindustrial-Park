#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Enzyme-burden (proteome-allocation) constraint for the kinetic Bayesian
optimization (kinetic_optimization.py).

Every capacity k_* of the nskinetics fermentation model is a specific
activity per unit biomass (enzyme amount x turnover), so raising it means
expressing more enzyme out of a finite proteome. This module converts the
sampled capacities to enzyme mass fractions (g enzyme / g dry cell weight),
sums them into the modeled metabolic pool Phi_M, closes the cell's
FLEXIBLE protein sector F_flex against Phi_M plus the translation sector
phi_T (which scales with the sampled growth capacities k_7 / k_8), and
derates growth linearly to zero as the flexible sector fills:

    native steps (r1, r2, r3, r5, r6):
        pool_i = pool_wt,i * max_j(k_ij / k_ij,ref)        (ratio route)
    r4: pool_r4 = pool_wt,r4                               (k_4 burden-free)
    Ehrlich steps (r13-r16):
        pool_i = k_i * MW_enz / (kcat * SIGMA_EFF * 3600 * MW_sub)
    Phi_M = sum of pools
    g = max(k_7/k_7,ref, k_8/k_8,ref);  phi_T = PHI_T_WT * g
    d = clip((F_flex - Phi_M) / phi_T, 0, 1)
    k_7_eff = d * k_7;  k_8_eff = d * k_8
    feasible  <=>  d > 0  <=>  Phi_M < F_flex

Exactly inert at the reference the model is built from (d = 1; apply is
the identity on k_7 / k_8), so a burden-off run and the smoke baselines
are unchanged. Design record: docs/superpowers/specs/
2026-09-05-enzyme-burden-design.md; wild-type pools from
docs/reports/enzyme-abundance-g-per-gDCW.md (both gitignored, so the
tables are hardcoded here with their sources).

Stdlib-only: imports neither numpy, biosteam nor nskinetics, so
analyses/test_enzyme_burden_offline.py and the stdlib-only BO supervisor
can load it. SCENARIO_B_EHRLICH is read from nskinetics' scenarios.py BY
FILE PATH (scenario_b_ehrlich), never by importing the package (which
pulls tellurium/roadrunner/biosteam).
"""
from dataclasses import dataclass
import importlib.util
import math
import os

__all__ = ('PROTEIN_CONTENT', 'POOL_TABLE_PROTEIN_CONTENT', 'HOUSEKEEPING_FRACTION',
           'TRANSLATION_FRACTION_WT', 'SIGMA_EFF', 'F_FLEX', 'PHI_T_WT',
           'NATIVE_STEPS', 'EHRLICH_STEPS', 'ANCHOR_STEPS',
           'GROWTH_CAPACITIES', 'STEP_ORDER', 'BURDEN_COLUMNS',
           'ehrlich_unit_cost', 'anchor_sigma',
           'BurdenResult', 'BurdenModel',
           'scenarios_path', 'scenario_b_ehrlich')

#%% Sector constants (spec 4.1; never change without asking)

#: g protein / gDCW. 0.49 = the crude-protein content of S. cerevisiae
#: biomass from the beer manufacturing process (Jach, Serefko, Ziaja &
#: Kieliszek 2022, Metabolites 12:63, doi:10.3390/metabo12010063, Table 1:
#: 49 % of dry weight) -- the upper end of the 0.40-0.50 range of the
#: enzyme-abundance report 3, whose central 0.45 was used until
#: 2026-09-06. The wild-type pool table NATIVE_STEPS was tabulated at
#: POOL_TABLE_PROTEIN_CONTENT and is rescaled to this value (report 3:
#: every pool is a proteome mass fraction x the protein content, so it
#: scales linearly with it).
PROTEIN_CONTENT = 0.49
#: g protein / gDCW at which the NATIVE_STEPS pools were tabulated
#: (enzyme-abundance report 3: proteome mass fraction x 0.45).
POOL_TABLE_PROTEIN_CONTENT = 0.45
#: Fraction of protein unavailable to the modeled + translation sectors
#: (Scott et al. 2010 zero-growth intercept ~ half the proteome;
#: Metzl-Raz et al. 2017; Xia et al. 2022).
HOUSEKEEPING_FRACTION = 0.50
#: phi_T,wt / P: the ACTIVE ribosomal/translation sector of fast-growing
#: yeast at the wild-type growth capacity. Metzl-Raz, Kafri, Yaari,
#: Shreberk-Shaked, Eldar & Barkai 2017, eLife 6:e28034: the ribosomal
#: proteins are ~30 % of the proteome of rapidly growing S. cerevisiae
#: (~8 % in the slowest-growing cells; Warner 1999, Trends Biochem. Sci.
#: 24:437 gives the same ~30 %), of which ~25 % does not contribute to
#: translation (an ~8 %-of-proteome reserve, constant across growth
#: rates) -- 0.30 x (1 - 0.25) = 0.225 is the translating part, the
#: sector the proportional ribosome law scales with k_7 (0.30, the whole
#: ribosomal fraction, until 2026-09-06). The housekeeping fraction was
#: deliberately left at 0.50 (the reserve is NOT re-booked there), so the
#: change frees 0.037 g/gDCW of reallocation slack at wild-type growth
#: (0.0286 -> 0.0654) and leaves the hard cap F_flex unchanged.
TRANSLATION_FRACTION_WT = 0.225
#: Average in-vivo enzyme saturation of the kcat route -- the fraction of
#: its in-vitro kcat an enzyme delivers in vivo. GECKO's fitted value for
#: S. cerevisiae (Sanchez et al. 2017, Mol. Syst. Biol. 13:935; kept as
#: the GECKO 3 default, Chen et al. 2024). Scales every Ehrlich pool
#: inversely, so it is the single largest lever on the Ehrlich cost. The
#: two native single-enzyme steps give sigma_r3 ~ 2.0 and sigma_r6 ~ 0.08
#: (geometric mean 0.41; anchor_sigma) -- a diagnostic, not a calibration.
SIGMA_EFF = 0.50

#: Flexible protein sector (g/gDCW): the zero-growth point of the modeled
#: pool. ~ 0.245 (0.225 at the pre-2026-09-06 PROTEIN_CONTENT of 0.45).
F_FLEX = PROTEIN_CONTENT*(1.0 - HOUSEKEEPING_FRACTION)
#: Translation sector at the wild-type growth capacity (g/gDCW). ~ 0.110
#: (0.147 at TRANSLATION_FRACTION_WT = 0.30, 0.135 at 0.30 and P = 0.45).
PHI_T_WT = PROTEIN_CONTENT*TRANSLATION_FRACTION_WT

#%% Wild-type pool table (spec 4.2; enzyme-abundance report 1)
# step -> (pool_wt in g enzyme / gDCW, model capacities charged by ratio).
# r4's pool is fixed: the model already carries X_AcDH as a literal pool
# synthesized by r9 at the cost of active biomass, so k_4 is a turnover
# number here (Q10). r7/r8 are growth (translation sector); r9-r11 are
# the AcDH synthesis/decay law itself. Sources: Ho, Baryshnikova & Brown
# 2018 (per-cell medians), PaxDB, Kulak 2014; tabulated at the report's
# protein content POOL_TABLE_PROTEIN_CONTENT (0.45) and rescaled to
# PROTEIN_CONTENT below.

_NATIVE_STEPS_TABLE = {
    'r1': (0.044, ('k_1h', 'k_1l', 'k_1e')),   # glycolysis lump (17 genes)
    'r2': (0.0032, ('k_2',)),                   # PDH complex + TCA
    'r3': (0.0085, ('k_3',)),                   # Pdc1 (+Pdc5/6)
    'r4': (0.0032, ()),                         # Ald6 (AcDH pool; fixed)
    'r5': (0.0008, ('k_5', 'k_5e')),            # Acs2 (+Acs1)
    'r6': (0.0040, ('k_6',)),                   # Adh1 (+Adh2-5)
}

#: The pool table at PROTEIN_CONTENT: pool_wt = tabulated pool x
#: PROTEIN_CONTENT / POOL_TABLE_PROTEIN_CONTENT (the identity at 0.45;
#: x 1.089 at 0.49 -- r1 0.0479, r2 0.00348, r3 0.00926, r4 0.00348,
#: r5 0.000871, r6 0.00436 g/gDCW; Phi_M,wt 0.0694).
NATIVE_STEPS = {
    step: (pool*PROTEIN_CONTENT/POOL_TABLE_PROTEIN_CONTENT, capacities)
    for step, (pool, capacities) in _NATIVE_STEPS_TABLE.items()
}

#%% Ehrlich enzyme table (spec 4.3; UniProt reviewed S288C masses)
# step -> (capacity, MW of the substrate the model writes k_i on [g/mol],
#          ((enzyme, MW_enzyme [Da], kcat [s^-1 per mole of that
#            substrate]), ...)).
# Specific activities are converted at the listed subunit mass:
# kcat [s^-1] = SA [umol min^-1 mg^-1] * MW [g/mol] * 1e-3 / 60.

def _kcat_from_specific_activity(U_per_mg, MW):
    """kcat (s^-1) from a specific activity (U/mg = umol min^-1 mg^-1)
    and the enzyme (subunit) mass MW (g/mol)."""
    return U_per_mg*MW*1e-3/60.0

EHRLICH_STEPS = {
    # Ilv2 (P07342, 74 937 Da) + one Ilv6 regulatory subunit (P25605,
    # 33 987 Da); Pang & Duggleby 1999 (Biochemistry 38:5222): reconstituted
    # Ilv2+Ilv6 49.0 U/mg Ilv2 -> 61 acetolactate s^-1 per Ilv2 subunit;
    # two pyruvate per acetolactate; BRENDA EC 2.2.1.6.
    'r13': ('k_13', 88.06,
            (('Ilv2+Ilv6', 108924.0,
              2.0*_kcat_from_specific_activity(49.0, 74937.0)),)),
    # Ilv5 KARI (P06168). Fungal surrogate: N. crassa reductoisomerase
    # 18.5 U/mg (Kiritani, Narise & Wagner 1966, JBC 241:2047; BRENDA EC
    # 1.1.1.86 ref 639171) at the Ilv5 mass; no S. cerevisiae entry exists
    # (E. coli IlvC 1.7 s^-1 would be ~8x more expensive).
    'r14': ('k_14', 132.11,
            (('Ilv5', 44368.0, _kcat_from_specific_activity(18.5, 44368.0)),)),
    # Ilv3 DHAD (P39522). Fungal surrogate: A. fumigatus Ilv3A 18 U/mg at
    # pH 8, 22 C (Oliver et al. 2012, PLoS ONE 7:e43559) at the Ilv3 mass;
    # no S. cerevisiae entry (E. coli IlvD 69 s^-1, Flint et al. 1993).
    'r15': ('k_15', 134.13,
            (('Ilv3', 62861.0, _kcat_from_specific_activity(18.0, 62861.0)),)),
    # KDC Aro10 (Q06408): kcat 19 s^-1, Km 8.5 mM on 2-ketoisovalerate,
    # pH 7.0, 30 C (Kneen et al. 2011, FEBS J 278:1842; BRENDA EC
    # 4.1.1.43). ADH Adh6 (Q04894): 296 s^-1 on 2-methylpropanal + NADPH,
    # pH 7.0, 25 C (Larroy et al. 2002, Biochem J 361:163; BRENDA EC
    # 1.1.1.2). Both sized on the KIV molar flux, k_16 / 116.12, because
    # the model writes k_16 in g KIV gDCW^-1 h^-1.
    'r16': ('k_16', 116.12,
            (('Aro10', 71384.0, 19.0),
             ('Adh6', 39618.0, 296.0))),
}

#: Native single-enzyme steps used as the sigma DIAGNOSTIC (spec 4.4):
#: step -> (capacity, MW_substrate [g/mol], MW_enzyme [Da], kcat [s^-1]).
#: Pdc1 (P06169) 60 s^-1 per subunit on pyruvate, pH 6.0, 30 C
#: (Balakrishnan et al. 2012, JACS 134:3873; BRENDA EC 4.1.1.1); Adh1
#: (P00330) 1800 s^-1 acetaldehyde reduction, pH 7.3, 30 C (Ganzhorn et
#: al. 1987, JBC 262:3754).
ANCHOR_STEPS = {
    'r3': ('k_3', 88.06, 61495.0, 60.0),
    'r6': ('k_6', 44.05, 36849.0, 1800.0),
}

#: Growth capacities derated by the burden (both draw on the translation
#: sector; one machinery sized by the larger demand, Q5a/Q8).
GROWTH_CAPACITIES = ('k_7', 'k_8')

#: Steps whose pools sum to Phi_M, in trajectory-column order.
STEP_ORDER = (*NATIVE_STEPS, *EHRLICH_STEPS)

#: Trajectory-CSV columns a burden study records (BurdenResult.as_record).
#: violation = Phi_M - F_flex is derivable, so it is not a column.
BURDEN_COLUMNS = (*(f'pool_{step}' for step in STEP_ORDER),
                  'Phi_M', 'phi_T', 'F_flex', 'burden_factor',
                  'k_7_eff', 'k_8_eff')

def ehrlich_unit_cost(step, sigma_eff=SIGMA_EFF):
    """g enzyme / gDCW needed per unit capacity (1 g gDCW^-1 h^-1) of an
    Ehrlich step: sum over its enzymes of
    MW_enz / (kcat * sigma_eff * 3600 * MW_sub)."""
    _, mw_sub, enzymes = EHRLICH_STEPS[step]
    return sum(mw_enz/(kcat*sigma_eff*3600.0*mw_sub)
               for _, mw_enz, kcat in enzymes)

def anchor_sigma(k_ref):
    """Implied in-vivo saturation of the two native single-enzyme steps,
    sigma_i = k_ref,i * MW_enz / (kcat * 3600 * MW_sub * pool_wt,i), and
    their geometric mean -- {'r3': ..., 'r6': ..., 'geomean': ...}. With
    the model's reference capacities: r3 ~ 2.0 (the fitted k_3 exceeds
    what the tag-biased-low Pdc1 pool delivers at its in-vitro kcat), r6
    ~ 0.08 (Adh1 is expressed in ten-fold excess); geometric mean ~ 0.41
    (0.45 at the pre-2026-09-06 protein content of 0.45), consistent with
    SIGMA_EFF. A diagnostic, not a calibration."""
    sigmas = {}
    for step, (capacity, mw_sub, mw_enz, kcat) in ANCHOR_STEPS.items():
        pool_wt = NATIVE_STEPS[step][0]
        sigmas[step] = float(k_ref[capacity])*mw_enz/(kcat*3600.0*mw_sub*pool_wt)
    sigmas['geomean'] = math.exp(sum(math.log(v) for v in sigmas.values())
                                 /len(ANCHOR_STEPS))
    return sigmas

#%% Burden model

@dataclass(frozen=True)
class BurdenResult:
    """One evaluated decision point. `pools` maps every step of
    STEP_ORDER to its enzyme mass fraction (g/gDCW); `burden_factor` is
    the growth derating d in [0, 1]; `k_7_eff` / `k_8_eff` are the
    derated growth capacities the model receives. `violation` =
    Phi_M - F_flex (<= 0 feasible) is what the sampler constraint
    reads; `feasible` <=> d > 0 <=> Phi_M < F_flex."""
    pools: dict
    Phi_M: float
    phi_T: float
    F_flex: float
    burden_factor: float
    k_7_eff: float
    k_8_eff: float

    @property
    def violation(self):
        return self.Phi_M - self.F_flex

    @property
    def feasible(self):
        return self.burden_factor > 0.0

    def as_record(self):
        """{column: value} over BURDEN_COLUMNS (the trajectory-CSV
        columns of a burden study)."""
        record = {f'pool_{step}': self.pools[step] for step in STEP_ORDER}
        record.update(Phi_M=self.Phi_M, phi_T=self.phi_T, F_flex=self.F_flex,
                      burden_factor=self.burden_factor,
                      k_7_eff=self.k_7_eff, k_8_eff=self.k_8_eff)
        return record


class BurdenModel:
    """The burden evaluated against a snapshot of the model's reference
    capacities `k_ref` ({name: value}; every name of
    required_capacities() must be present -- the engine passes its
    kinetic_baselines, i.e. the live values after the scenario workbook
    was loaded). Native steps are charged by ratio to the reference (so
    the burden is exactly inert there), Ehrlich steps by kcat + MW at
    `sigma_eff` (module default SIGMA_EFF).

    Construction raises KeyError on a missing capacity, ValueError on a
    non-positive native or growth reference (the ratio route is
    undefined), on a constant set that makes the wild type itself
    infeasible (Phi_M,wt + phi_T,wt > F_flex), and on a reference point
    that is infeasible under the burden (e.g. the scenario-B Ehrlich
    constants, spec 4.5 / Q11): a burden study must start from a
    burden-feasible reference -- use burden=False / --no-burden for a
    burden-free study from such a point."""

    def __init__(self, k_ref, sigma_eff=SIGMA_EFF):
        self.sigma_eff = float(sigma_eff)
        if self.sigma_eff <= 0.0:
            raise ValueError('sigma_eff must be positive.')
        required = self.required_capacities()
        missing = [name for name in required if name not in k_ref]
        if missing:
            raise KeyError(f'k_ref lacks the capacities {missing} needed by '
                           'the enzyme burden (pass the engine\'s full '
                           'kinetic_baselines).')
        self.reference = {name: float(k_ref[name]) for name in required}
        nonpositive = [name for name in (*self._native_capacities(),
                                         *GROWTH_CAPACITIES)
                       if self.reference[name] <= 0.0]
        if nonpositive:
            raise ValueError('reference capacities charged by ratio must be '
                             f'positive; got {nonpositive}.')
        self.F_flex = F_FLEX
        self.phi_T_wt = PHI_T_WT
        self.Phi_M_wt = sum(pool_wt for pool_wt, _ in NATIVE_STEPS.values())
        if self.Phi_M_wt + self.phi_T_wt > self.F_flex:
            raise ValueError(
                'the sector constants make the wild type infeasible: '
                f'Phi_M,wt {self.Phi_M_wt:.4f} + phi_T,wt {self.phi_T_wt:.4f} '
                f'> F_flex {self.F_flex:.4f} g/gDCW.')
        self.sigma_diagnostic = anchor_sigma(self.reference)
        self.reference_result = self.evaluate(self.reference)
        if not self.reference_result.feasible:
            ehrlich = sum(self.reference_result.pools[s] for s in EHRLICH_STEPS)
            raise ValueError(
                'the reference point is infeasible under the enzyme burden '
                f'(Phi_M {self.reference_result.Phi_M:.4f} >= F_flex '
                f'{self.F_flex:.4f} g/gDCW; its Ehrlich capacities alone need '
                f'{ehrlich:.4f}). A burden study must start from a '
                'burden-feasible reference (the scenario-A baseline); pass '
                'burden=False / --no-burden for a burden-free study from '
                'this point.')

    @classmethod
    def from_reference(cls, k_ref, sigma_eff=SIGMA_EFF):
        """Snapshot `k_ref` (see the class docstring)."""
        return cls(k_ref, sigma_eff=sigma_eff)

    @staticmethod
    def _native_capacities():
        return tuple(c for _, caps in NATIVE_STEPS.values() for c in caps)

    @classmethod
    def required_capacities(cls):
        """Every capacity the burden reads, in table order: the native
        ratio-route capacities, the Ehrlich capacities, then k_7, k_8."""
        return (*cls._native_capacities(),
                *(cap for cap, _, _ in EHRLICH_STEPS.values()),
                *GROWTH_CAPACITIES)

    def _value(self, values, name):
        """`values[name]`, or the reference when the key is absent (a
        study that excludes some capacities still evaluates)."""
        return float(values[name]) if name in values else self.reference[name]

    def pools(self, values):
        """{step: g enzyme/gDCW} over STEP_ORDER at the decision point
        `values` ({name: value}; missing names -> reference)."""
        pools = {}
        for step, (pool_wt, capacities) in NATIVE_STEPS.items():
            multiplier = max((self._value(values, c)/self.reference[c]
                              for c in capacities), default=1.0)
            pools[step] = pool_wt*multiplier
        for step, (capacity, _, _) in EHRLICH_STEPS.items():
            pools[step] = (self._value(values, capacity)
                           *ehrlich_unit_cost(step, self.sigma_eff))
        return pools

    def evaluate(self, values):
        """BurdenResult at the decision point `values` (see pools)."""
        pools = self.pools(values)
        Phi_M = sum(pools.values())
        k_7 = self._value(values, 'k_7')
        k_8 = self._value(values, 'k_8')
        g = max(k_7/self.reference['k_7'], k_8/self.reference['k_8'])
        phi_T = self.phi_T_wt*g
        if phi_T > 0.0:
            d = min(1.0, max(0.0, (self.F_flex - Phi_M)/phi_T))
        else:  # both growth capacities sampled at zero: nothing to derate
            d = 1.0 if Phi_M < self.F_flex else 0.0
        return BurdenResult(pools=pools, Phi_M=Phi_M, phi_T=phi_T,
                            F_flex=self.F_flex, burden_factor=d,
                            k_7_eff=d*k_7, k_8_eff=d*k_8)

    def apply(self, values):
        """Copy of `values` with k_7 and k_8 set to their derated
        effective values (both keys are always present in the result,
        from the reference when not in `values`, so the engine's setattr
        loop writes the derating even when growth is not a decision
        variable); nothing else is touched."""
        result = self.evaluate(values)
        applied = dict(values)
        applied['k_7'] = result.k_7_eff
        applied['k_8'] = result.k_8_eff
        return applied

    def describe_point(self, values, label=''):
        """Printable burden report of the decision point `values`: every
        step's pool and its multiple of the wild type, Phi_M, F_flex,
        phi_T, the burden factor d with the effective growth capacities,
        a FEASIBLE / INFEASIBLE status line and the sigma diagnostic.
        The driver prints it for the scenario-A reference and for the
        scenario-B Ehrlich constants (the Q11 sanity report)."""
        result = self.evaluate(values)
        title = 'Enzyme burden report' + (f' -- {label}' if label else '')
        lines = [f'{title} (g enzyme / gDCW)',
                 f'  {"step":<5}{"pool":>10}{"x wild type":>13}   capacities']
        for step, (pool_wt, capacities) in NATIVE_STEPS.items():
            caps = ', '.join(f'{c}={self._value(values, c):.4g}'
                             for c in capacities) or 'fixed (k_4 burden-free)'
            lines.append(f'  {step:<5}{result.pools[step]:>10.4f}'
                         f'{result.pools[step]/pool_wt:>13.2f}   {caps}')
        for step, (capacity, _, enzymes) in EHRLICH_STEPS.items():
            names = '+'.join(name for name, _, _ in enzymes)
            lines.append(f'  {step:<5}{result.pools[step]:>10.4f}{"-":>13}   '
                         f'{capacity}={self._value(values, capacity):.4g} '
                         f'({names}; {ehrlich_unit_cost(step, self.sigma_eff):.4f} '
                         f'per unit at sigma_eff {self.sigma_eff:g})')
        lines.append(f'  Phi_M = {result.Phi_M:.4f} (wild type {self.Phi_M_wt:.4f}); '
                     f'F_flex = {result.F_flex:.4f}; phi_T = {result.phi_T:.4f} '
                     f'(wild type {self.phi_T_wt:.4f}, '
                     f'g = {result.phi_T/self.phi_T_wt:.2f})')
        k_7 = self._value(values, 'k_7')
        k_8 = self._value(values, 'k_8')
        lines.append(f'  burden factor d = {result.burden_factor:.4f} -> '
                     f'k_7_eff = {result.k_7_eff:.4g} (sampled {k_7:.4g}), '
                     f'k_8_eff = {result.k_8_eff:.4g} (sampled {k_8:.4g})')
        if result.feasible:
            lines.append(f'  status: FEASIBLE (Phi_M - F_flex = {result.violation:+.4f}; '
                         'reallocation slack at undiminished growth '
                         f'{result.F_flex - result.Phi_M - result.phi_T:+.4f})')
        else:
            lines.append(f'  status: INFEASIBLE (Phi_M exceeds F_flex by '
                         f'{result.violation:.4f}; growth capacity zero -- such '
                         'a trial is pruned before simulating)')
        sd = self.sigma_diagnostic
        lines.append('  sigma diagnostic (native anchors): '
                     f'sigma_r3 = {sd["r3"]:.2f}, sigma_r6 = {sd["r6"]:.3f}, '
                     f'geometric mean {sd["geomean"]:.2f} vs SIGMA_EFF = {SIGMA_EFF:g}')
        return '\n'.join(lines)

#%% nskinetics scenario presets (read by file path)

#: Path of scenarios.py relative to the nskinetics package directory.
_SCENARIOS_RELPATH = ('models', 's_cerevisiae_ferm_fb_inhib_mod_ibo',
                      'scenarios.py')

def scenarios_path():
    """Absolute path of nskinetics' scenarios.py for the shipped
    S. cerevisiae ethanol/isobutanol model, located via the package's
    spec WITHOUT importing it (importlib.util.find_spec only resolves
    the file) -- the kinetic_optimization.kinetic_parameter_roles_path
    pattern."""
    spec = importlib.util.find_spec('nskinetics')
    if spec is None or not spec.origin:
        raise ImportError('nskinetics is not installed (needed for the '
                          'scenario-B Ehrlich constants).')
    return os.path.join(os.path.dirname(spec.origin), *_SCENARIOS_RELPATH)

def scenario_b_ehrlich(path=None):
    """A copy of nskinetics' SCENARIO_B_EHRLICH ({'k_13': 5.81, 'k_14':
    4.8, 'k_15': 4.8, 'k_16': 2.82, 'k_16r': 0.0125}) -- the Ehrlich
    constants of scenario B, whose implied burden the driver reports (Q11).
    The file is executed BY PATH (spec_from_file_location + exec_module):
    it imports nothing, whereas importing it through the package pulls
    tellurium/roadrunner/biosteam, which the offline tests and the
    stdlib-only supervisor must never do. `path` overrides the located
    file (tests)."""
    spec = importlib.util.spec_from_file_location(
        '_nskinetics_scenarios', scenarios_path() if path is None else path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return dict(module.SCENARIO_B_EHRLICH)
