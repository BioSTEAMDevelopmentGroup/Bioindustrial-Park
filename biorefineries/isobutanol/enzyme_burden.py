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

__all__ = ('PROTEIN_CONTENT', 'HOUSEKEEPING_FRACTION',
           'TRANSLATION_FRACTION_WT', 'SIGMA_EFF', 'F_FLEX', 'PHI_T_WT',
           'NATIVE_STEPS', 'EHRLICH_STEPS', 'ANCHOR_STEPS',
           'GROWTH_CAPACITIES', 'STEP_ORDER', 'BURDEN_COLUMNS',
           'ehrlich_unit_cost', 'anchor_sigma',
           'BurdenResult', 'BurdenModel',
           'scenarios_path', 'scenario_b_ehrlich')

#%% Sector constants (spec 4.1; never change without asking)

#: g protein / gDCW (enzyme-abundance report 3; central of 0.40-0.50).
PROTEIN_CONTENT = 0.45
#: Fraction of protein unavailable to the modeled + translation sectors
#: (Scott et al. 2010 zero-growth intercept ~ half the proteome;
#: Metzl-Raz et al. 2017; Xia et al. 2021).
HOUSEKEEPING_FRACTION = 0.50
#: phi_T,wt / P: the ribosomal/translation sector of fast-growing yeast
#: at the wild-type growth capacity (Metzl-Raz 2017; Xia 2021;
#: Bjorkeroth 2020).
TRANSLATION_FRACTION_WT = 0.30
#: Average in-vivo enzyme saturation of the kcat route -- the fraction of
#: its in-vitro kcat an enzyme delivers in vivo. GECKO's fitted value for
#: S. cerevisiae (Sanchez et al. 2017, Mol. Syst. Biol. 13:935; kept as
#: the GECKO 3 default, Chen et al. 2024). Scales every Ehrlich pool
#: inversely, so it is the single largest lever on the Ehrlich cost. The
#: two native single-enzyme steps give sigma_r3 ~ 2.2 and sigma_r6 ~ 0.09
#: (geometric mean 0.45; anchor_sigma) -- a diagnostic, not a calibration.
SIGMA_EFF = 0.50

#: Flexible protein sector (g/gDCW): the zero-growth point of the modeled
#: pool. ~ 0.225.
F_FLEX = PROTEIN_CONTENT*(1.0 - HOUSEKEEPING_FRACTION)
#: Translation sector at the wild-type growth capacity (g/gDCW). ~ 0.135.
PHI_T_WT = PROTEIN_CONTENT*TRANSLATION_FRACTION_WT

#%% Wild-type pool table (spec 4.2; enzyme-abundance report 1)
# step -> (pool_wt in g enzyme / gDCW, model capacities charged by ratio).
# r4's pool is fixed: the model already carries X_AcDH as a literal pool
# synthesized by r9 at the cost of active biomass, so k_4 is a turnover
# number here (Q10). r7/r8 are growth (translation sector); r9-r11 are
# the AcDH synthesis/decay law itself. Sources: Ho, Baryshnikova & Brown
# 2018 (per-cell medians), PaxDB, Kulak 2014; protein content 0.45.

NATIVE_STEPS = {
    'r1': (0.044, ('k_1h', 'k_1l', 'k_1e')),   # glycolysis lump (17 genes)
    'r2': (0.0032, ('k_2',)),                   # PDH complex + TCA
    'r3': (0.0085, ('k_3',)),                   # Pdc1 (+Pdc5/6)
    'r4': (0.0032, ()),                         # Ald6 (AcDH pool; fixed)
    'r5': (0.0008, ('k_5', 'k_5e')),            # Acs2 (+Acs1)
    'r6': (0.0040, ('k_6',)),                   # Adh1 (+Adh2-5)
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
    the model's reference capacities: r3 ~ 2.2 (the fitted k_3 exceeds
    what the tag-biased-low Pdc1 pool delivers at its in-vitro kcat), r6
    ~ 0.09 (Adh1 is expressed in ten-fold excess); geometric mean ~ 0.45,
    consistent with SIGMA_EFF. A diagnostic, not a calibration."""
    sigmas = {}
    for step, (capacity, mw_sub, mw_enz, kcat) in ANCHOR_STEPS.items():
        pool_wt = NATIVE_STEPS[step][0]
        sigmas[step] = float(k_ref[capacity])*mw_enz/(kcat*3600.0*mw_sub*pool_wt)
    sigmas['geomean'] = math.exp(sum(math.log(v) for v in sigmas.values())
                                 /len(ANCHOR_STEPS))
    return sigmas
