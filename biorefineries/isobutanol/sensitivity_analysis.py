#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Variance-based global sensitivity analysis over a kinetic-optimization
campaign space (spec docs/superpowers/specs/2026-09-20-sobol-sensitivity-
split12d-design.md). SIM-FREE and package-free: stdlib + numpy / scipy /
sklearn only, loadable by file path; callers pass the kinetic_optimization
module (`ko`) and the enzyme_burden module (`eb`) in where needed.

The campaign space is a box with a closed-form feasibility constraint
(enzyme burden Phi_M < F_flex; fed-batch volume-ratio bound <= cap), so the
inputs are DEPENDENT on the feasible domain. The quantities estimated here
are well defined under dependence: the closed subset index
S_u = Var(E[Y | X_u]) / Var(Y) (pick-freeze with conditional rejection
sampling) and the exact Shapley effects built from all 2^d of them."""
import dataclasses
import math
import warnings

import numpy as np

__all__ = ('campaign_engine_kwargs', 'screening_search_space',
           'custom_search_space', 'search_space_tag', 'design_record', 'unit_to_values',
           'feasible_sobol_stream', 'VectorizedFeasibility',
           'sample_feasible', 'conditional_partners', 'closed_index',
           'all_closed_indices', 'first_order', 'total_order',
           'shapley_effects', 'best_subsets', 'smallest_subset_reaching',
           'mask_names', 'tail_variance_share',
           'RELIABLE_Q2', 'GP_HYPER_N', 'SURROGATE_CANDIDATES', 'Surrogate',
           'fit_surrogates')

#%% Campaign space

_PRESET_ENGINE_KEYS = ('include_params', 'exclude_params', 'multiplier_bounds',
                       'rate_multiplier_bounds', 'rate_params',
                       'parameter_multiplier_bounds', 'stage_1_max_x_bounds',
                       'parameter_groups', 'group_multiplier_bounds',
                       'spike_delta_bounds', 'group_references',
                       'param_bounds_override')

def campaign_engine_kwargs(ko, study_target_products, study_type):
    """(preset, engine_kwargs) exactly as analyses/optimize_kinetics_BO.run
    derives them for a named preset: the preset's keys, with
    param_bounds_override = the bounds workbook's bands
    (ko.workbook_kinetic_bounds of preset['kinetic_bounds_scenario'])
    UPDATED by the preset's own override (the scenario-A-anchored IBO
    pathway bands win). engine_kwargs feeds ko.build_search_space and
    ko._prepare_optimization unchanged."""
    preset = ko.resolve_study_preset(study_target_products, study_type)
    engine_kwargs = {key: preset[key] for key in _PRESET_ENGINE_KEYS}
    derived = ko.workbook_kinetic_bounds(
        preset['kinetic_bounds_scenario'],
        multiplier_bounds=engine_kwargs['multiplier_bounds'],
        rate_multiplier_bounds=engine_kwargs['rate_multiplier_bounds'],
        rate_params=engine_kwargs['rate_params'],
        parameter_multiplier_bounds=engine_kwargs['parameter_multiplier_bounds'])
    derived.update(engine_kwargs['param_bounds_override'] or {})
    engine_kwargs['param_bounds_override'] = derived
    return preset, engine_kwargs

#: The pre-optimization SCREENING measure (spec
#: docs/superpowers/specs/2026-09-21-sobol-screening-measure-design.md): the
#: three from-zero pathway axes are sampled LINEAR-uniform on [0, ceiling]
#: (a log band on an axis whose baseline is 0 hangs off an arbitrary
#: numerical floor and starves the co-production corner: 0.03 % of the
#: campaign measure), and the three native ethanol-branch axes on a
#: 0.1x-4x band (the campaign's 1e-3x knock-outs are the failure tail that
#: carries 100 % of Var(PI); a knock-out is a discrete decision, not a dial).
SCREENING_LINEAR_AXES = ('k_13', 'k_17', 'ehrlich_downstream')
SCREENING_NATIVE_AXES = ('k_3', 'k_6', 'glycolysis')
SCREENING_NATIVE_BAND = (0.1, 4.0)

def screening_search_space(search_space, kinetic_baselines):
    """A NEW search-space dict (same insertion order) carrying the screening
    measure: SCREENING_LINEAR_AXES -> {'low': 0.0, 'high': <campaign
    ceiling>, 'log': False}; 'k_3' / 'k_6' -> SCREENING_NATIVE_BAND x their
    live baseline, log; 'glycolysis' (a multiplier, baseline 1) ->
    SCREENING_NATIVE_BAND itself, log; every other entry copied. The
    measure is then carried by the existing linear branch of unit_to_values
    / ko.unit_to_external, so VectorizedFeasibility, ko.external_to_unit and
    design_record need no change. KeyError for a missing axis; ValueError
    for a nonpositive native baseline or an int-typed axis. `search_space`
    is not mutated."""
    m_lo, m_hi = SCREENING_NATIVE_BAND
    out = {}
    for name, sp in search_space.items():
        if name in SCREENING_LINEAR_AXES or name in SCREENING_NATIVE_AXES:
            if sp.get('int'):
                raise ValueError(f'{name!r} is an integer axis; the screening '
                                 'measure re-specifies float axes only')
        if name in SCREENING_LINEAR_AXES:
            out[name] = {'low': 0.0, 'high': float(sp['high']), 'log': False}
        elif name == 'glycolysis':
            out[name] = {'low': float(m_lo), 'high': float(m_hi), 'log': True}
        elif name in SCREENING_NATIVE_AXES:
            base = float(kinetic_baselines[name])
            if not base > 0.0:
                raise ValueError(f'{name!r}: the screening band is a multiplier '
                                 f'of a positive live baseline; got {base!r}')
            out[name] = {'low': m_lo*base, 'high': m_hi*base, 'log': True}
        else:
            out[name] = dict(sp)
    missing = [n for n in SCREENING_LINEAR_AXES + SCREENING_NATIVE_AXES if n not in out]
    if missing:
        raise KeyError(f'screening measure: axes {missing} are not in the search space')
    return out

def custom_search_space(search_space, bands):
    """A NEW search-space dict (same insertion order) carrying a CUSTOM
    measure: every axis of `search_space` must appear in `bands` (and no
    other), each as {'low', 'high', 'log'} in ABSOLUTE decision-variable
    units (a group axis = its multiplier / its referenced absolute value, as
    in the search space). An integer axis keeps 'int': True and needs
    integer bounds; a log axis needs low > 0; low < high always. The
    measure is carried by the existing branches of unit_to_values /
    ko.unit_to_external, so nothing downstream changes. `search_space` is
    not mutated."""
    extra = sorted(set(bands) - set(search_space))
    missing = [n for n in search_space if n not in bands]
    if extra or missing:
        raise KeyError(f'custom measure: bands must cover exactly the search-space axes; '
                       f'missing {missing}, unknown {extra}')
    out = {}
    for name, sp in search_space.items():
        b = bands[name]
        lo, hi, log = float(b['low']), float(b['high']), bool(b['log'])
        if not lo < hi:
            raise ValueError(f'{name!r}: low {lo!r} must be < high {hi!r}')
        if log and not lo > 0.0:
            raise ValueError(f'{name!r}: a log band needs low > 0; got {lo!r}')
        if sp.get('int'):
            if log or lo != int(lo) or hi != int(hi):
                raise ValueError(f'{name!r} is an integer axis: linear, integer bounds only')
            out[name] = {'low': int(lo), 'high': int(hi), 'log': False, 'int': True}
        else:
            out[name] = {'low': lo, 'high': hi, 'log': log}
    return out

def search_space_tag(search_space, n_hex=8):
    """Short content hash of a search space (sha1 of its canonical JSON) for
    study names: a custom band set is too long to spell out under the
    260-character path limit; the design record carries the full bands."""
    import hashlib
    import json
    blob = json.dumps(search_space, sort_keys=True, separators=(',', ':'))
    return hashlib.sha1(blob.encode()).hexdigest()[:n_hex]

def design_record(*, search_space, parameter_groups, group_references,
                  kinetic_baselines, burden_model, eb, baseline_model_kwargs,
                  baseline_max_n_spikes, volume_cap, target_conc_max,
                  meta=None):
    """Everything stage 2 needs to rebuild the feasible domain WITHOUT the
    model, as a JSON-able dict (written by stage 1 as <study>_design.json).
    The burden block is the closed form BurdenModel.pools evaluates: native
    steps pool_wt x max_c(value_c / reference_c), Ehrlich steps value x
    unit cost; feasible <=> Phi_M < F_flex."""
    for legacy in ('spike_delta', 'target_conc', 'threshold_delta', 'spike_conc'):
        if legacy in search_space:
            raise NotImplementedError(
                f'{legacy!r} is sampled: VectorizedFeasibility covers only a '
                'pinned spike feed on the threshold-anchored scheme.')
    return dict(
        search_space={name: {k: (bool(v) if k in ('log', 'int') else float(v))
                             for k, v in sp.items()}
                      for name, sp in search_space.items()},
        parameter_groups={g: list(m) for g, m in (parameter_groups or {}).items()},
        group_references={g: {m: float(r) for m, r in refs.items()}
                          for g, refs in (group_references or {}).items()},
        kinetic_baselines={n: float(v) for n, v in kinetic_baselines.items()},
        burden=dict(
            F_flex=float(burden_model.F_flex),
            reference={n: float(v) for n, v in burden_model.reference.items()},
            native_steps={step: [float(pool_wt), list(caps)]
                          for step, (pool_wt, caps) in eb.NATIVE_STEPS.items()},
            ehrlich_costs={cap: float(eb.ehrlich_unit_cost(step, burden_model.sigma_eff))
                           for step, (cap, _, _) in eb.EHRLICH_STEPS.items()}),
        volume=dict(spike_conc=float(baseline_model_kwargs['spike_conc']),
                    baseline_max_n_spikes=int(baseline_max_n_spikes),
                    volume_cap=float(volume_cap),
                    target_conc_max=float(target_conc_max)),
        meta=dict(meta or {}))

#%% Unit cube <-> decision values (vectorized twin of ko.unit_to_external)

def unit_to_values(U, search_space):
    """{name: ndarray} for unit-cube rows U (n, d), columns in search-space
    order: log floats exp(ln lo + u (ln hi - ln lo)) clipped into [lo, hi],
    ints lo + floor(u (hi - lo + 1)) clipped to hi (int64), linear floats
    lo + u (hi - lo) -- ko.unit_to_internal's measure, elementwise."""
    U = np.clip(np.asarray(U, dtype=float), 0.0, 1.0)
    values = {}
    for j, (name, sp) in enumerate(search_space.items()):
        u, lo, hi = U[:, j], sp['low'], sp['high']
        if sp.get('int'):
            lo, hi = int(lo), int(hi)
            values[name] = np.minimum(hi, lo + np.floor(u*(hi - lo + 1))).astype(np.int64)
        elif sp['log']:
            x = np.exp(math.log(lo) + u*(math.log(hi) - math.log(lo)))
            values[name] = np.clip(x, lo, hi)
        else:
            values[name] = lo + u*(hi - lo)
    return values

def feasible_sobol_stream(search_space, is_feasible, seed, unit_to_external,
                          block=1024):
    """Endless iterator of (candidate_index, values): scrambled-Sobol'
    candidates of the unit cube in blocks of `block` (a power of 2),
    mapped by `unit_to_external` (pass ko.unit_to_external -- the ONE
    definition of the campaign measure) and kept when `is_feasible(values)`.
    Rejection of a uniform stream = uniform on the feasible set. Fully
    deterministic in `seed`, so the accepted sequence is NESTED (its first n
    points are a valid design for every n) and resuming = regenerating and
    skipping."""
    from scipy.stats import qmc
    engine = qmc.Sobol(d=len(search_space), scramble=True, seed=seed)
    candidate_index = 0
    while True:
        for u in engine.random(block):
            values = unit_to_external(u, search_space)
            if is_feasible(values):
                yield candidate_index, values
            candidate_index += 1

#%% Vectorized feasibility (stage 2 needs ~1e7-1e8 checks)

class VectorizedFeasibility:
    """Batch twin of ko.feasibility_predicate(burden_on=True, volume_on=True)
    over unit-cube arrays, rebuilt from a design_record dict. The scalar
    predicate stays the authority: the offline test and stage 2 both
    cross-check this class against it / against the recorded Phi_M."""

    def __init__(self, design):
        self.search_space = design['search_space']
        self.names = list(self.search_space)
        self.d = len(self.names)
        self.parameter_groups = design['parameter_groups']
        self.group_references = design['group_references']
        self.kinetic_baselines = design['kinetic_baselines']
        self.burden = design['burden']
        self.volume = design['volume']

    @classmethod
    def from_design(cls, design):
        return cls(design)

    def _applied(self, values):
        """ko.expand_grouped_values, vectorized: member = basis x multiplier
        (basis = the group's reference when it has one, else the live
        baseline); individually sampled kinetics pass through."""
        applied = {}
        for name, x in values.items():
            if name in self.parameter_groups:
                refs = self.group_references.get(name)
                for member in self.parameter_groups[name]:
                    basis = (self.kinetic_baselines[member] if refs is None
                             else refs[member])
                    applied[member] = basis*x
            elif name in self.kinetic_baselines:
                applied[name] = x
        return applied

    def _phi_M(self, values):
        applied, reference = self._applied(values), self.burden['reference']
        n = len(next(iter(values.values())))
        def val(c):
            return applied[c] if c in applied else np.full(n, reference[c])
        phi = np.zeros(n)
        for pool_wt, caps in self.burden['native_steps'].values():
            # A native step with no capacities (r4, the fixed Ald6 pool) has a
            # constant unit multiplier -- BurdenModel.pools uses default=1.0 --
            # so guard the empty reduction (np.max([]) has no identity).
            multiplier = (np.max([val(c)/reference[c] for c in caps], axis=0)
                          if caps else 1.0)
            phi += pool_wt*multiplier
        for cap, cost in self.burden['ehrlich_costs'].items():
            phi += val(cap)*cost
        return phi

    def phi_M(self, U):
        return self._phi_M(unit_to_values(U, self.search_space))

    def __call__(self, U):
        values = unit_to_values(U, self.search_space)
        ok = self._phi_M(values) < self.burden['F_flex']
        vol = self.volume
        thr = values['threshold_conc']
        tgt = np.minimum(vol['target_conc_max'], thr + values['target_delta'])
        spk = vol['spike_conc']
        n_spikes = (values['max_n_spikes'] if 'max_n_spikes' in values
                    else np.full(len(thr), vol['baseline_max_n_spikes']))
        with np.errstate(divide='ignore', invalid='ignore'):
            log_ratio = n_spikes*np.log((spk - thr)/(spk - tgt))
        log_ratio = np.where(spk <= tgt, np.inf, log_ratio)
        log_ratio = np.where(n_spikes <= 0, 0.0, log_ratio)
        return ok & (log_ratio <= math.log(vol['volume_cap']))

#%% Index estimators on the feasible domain

def sample_feasible(n, d, is_feasible_U, rng):
    """n uniform rows of the unit cube satisfying `is_feasible_U` (a
    vectorized predicate U -> bool array), by batch rejection."""
    out, have = [], 0
    while have < n:
        cand = rng.random((max(1024, 2*(n - have)), d))
        cand = cand[is_feasible_U(cand)]
        out.append(cand)
        have += len(cand)
    return np.concatenate(out)[:n]

def _mask_columns(mask, d):
    return np.array([(mask >> j) & 1 for j in range(d)], dtype=bool)

def conditional_partners(X, mask, is_feasible_U, rng, max_tries=200):
    """Pick-freeze partners under DEPENDENCE: Xp shares the columns in
    `mask` with X and redraws the others uniformly until the whole row is
    feasible, i.e. Xp_~u ~ P(X_~u | X_u). A row still infeasible after
    `max_tries` rounds keeps its own X_~u (always feasible; a slight upward
    bias) -- returns (Xp, that count)."""
    free = ~_mask_columns(mask, X.shape[1])
    Xp, todo = X.copy(), np.arange(len(X))
    for _ in range(max_tries):
        if todo.size == 0:
            break
        cand = X[todo].copy()
        cand[:, free] = rng.random((todo.size, int(free.sum())))
        ok = is_feasible_U(cand)
        Xp[todo[ok]] = cand[ok]
        todo = todo[~ok]
    return Xp, int(todo.size)

def closed_index(y, yp):
    """Var(E[Y | X_u]) / Var(Y) from a pick-freeze pair (Janon-Monod
    estimator: mean and variance pooled over both samples)."""
    mu = 0.5*(y.mean() + yp.mean())
    var = 0.5*((y*y).mean() + (yp*yp).mean()) - mu*mu
    return float(((y*yp).mean() - mu*mu)/var) if var > 0.0 else float('nan')

def all_closed_indices(predictors, is_feasible_U, d, *, n_base=4096,
                       n_replicates=3, seed=0, max_tries=200, progress=None):
    """Closed indices of EVERY subset (bitmask 0 .. 2^d - 1; bit j <-> input
    j) for every metric in `predictors` ({metric: callable(U) -> y}). One
    base sample and one partner sample per (replicate, subset) are shared by
    all metrics. S[..., 0] = 0 and S[..., 2^d - 1] = 1 by definition.
    Returns (S {metric: (n_replicates, 2^d)}, fallback_fraction (2^d,))."""
    n_masks, full = 2**d, 2**d - 1
    S = {m: np.zeros((n_replicates, n_masks)) for m in predictors}
    fallback = np.zeros(n_masks)
    for r in range(n_replicates):
        rng = np.random.default_rng([seed, r])
        X = sample_feasible(n_base, d, is_feasible_U, rng)
        y = {m: np.asarray(f(X), dtype=float) for m, f in predictors.items()}
        for mask in range(1, full):
            Xp, n_fb = conditional_partners(X, mask, is_feasible_U, rng, max_tries)
            fallback[mask] += n_fb/(n_base*n_replicates)
            for m, f in predictors.items():
                S[m][r, mask] = closed_index(y[m], np.asarray(f(Xp), dtype=float))
            if progress is not None:
                progress(r, mask, full - 1)
        for m in predictors:
            S[m][r, full] = 1.0
    return S, fallback

def first_order(S, d):
    return np.stack([S[..., 1 << i] for i in range(d)], axis=-1)

def total_order(S, d):
    full = 2**d - 1
    return np.stack([1.0 - S[..., full ^ (1 << i)] for i in range(d)], axis=-1)

def _popcounts(d):
    return np.array([bin(m).count('1') for m in range(2**d)])

def shapley_effects(S, d):
    """Exact Shapley effects from all closed indices: Sh_i = sum over
    u not containing i of |u|! (d - |u| - 1)! / d! x (S[u + i] - S[u]).
    Sums to S[full] - S[empty] = 1 by construction."""
    sizes, masks = _popcounts(d), np.arange(2**d)
    weight = np.array([math.factorial(k)*math.factorial(d - k - 1)/math.factorial(d)
                       for k in range(d)])
    out = []
    for i in range(d):
        bit = 1 << i
        without = masks[(masks & bit) == 0]
        out.append(np.sum(weight[sizes[without]]
                          *(S[..., without | bit] - S[..., without]), axis=-1))
    return np.stack(out, axis=-1)

def best_subsets(S_mean, d, sizes=(1, 2, 3, 4)):
    """{k: (mask, S)} -- the size-k subset with the largest closed index."""
    pop, best = _popcounts(d), {}
    for k in sizes:
        masks = np.flatnonzero(pop == k)
        if masks.size:
            j = masks[np.argmax(S_mean[masks])]
            best[k] = (int(j), float(S_mean[j]))
    return best

def smallest_subset_reaching(S_mean, d, share=0.8):
    """(mask, S) of the best subset of the SMALLEST size whose closed index
    reaches `share` (the full set always does)."""
    for k, (mask, value) in best_subsets(S_mean, d, sizes=range(1, d + 1)).items():
        if value >= share:
            return mask, value
    return 2**d - 1, 1.0

def mask_names(mask, names):
    return tuple(name for j, name in enumerate(names) if (mask >> j) & 1)

def tail_variance_share(y, threshold=0.0):
    """Share of Var(y) that would vanish if the tail (y < threshold) were
    absent from the variance budget: 1 - p_hi Var(y | y >= threshold) /
    Var(y) (law of total variance: everything except the within-variance of
    the non-tail group). 0 when there is no tail."""
    y = np.asarray(y, dtype=float)
    hi = y[y >= threshold]
    if hi.size == y.size:
        return 0.0
    within_hi = (hi.size/y.size)*np.var(hi) if hi.size else 0.0
    return float(1.0 - within_hi/np.var(y))

#%% Surrogates

#: Below this cross-validated Q2 a metric's indices are still computed but
#: flagged unreliable in every output table.
RELIABLE_Q2 = 0.8
#: Rows used to optimize the GP's ARD hyper-parameters (then held fixed).
GP_HYPER_N = 1500
#: Surrogate families, in the order they are built and scored (also the
#: tie-break order of the best-Q2 pick).
SURROGATE_CANDIDATES = ('gp', 'hgb')

@dataclasses.dataclass
class Surrogate:
    name: str          # 'gp' | 'hgb'
    model: object
    q2: dict           # {candidate: Q2} (n_folds-fold CV)
    reliable: bool

    def predict(self, U, chunk=4096):
        U = np.asarray(U, dtype=float)
        if self.name != 'gp':
            return self.model.predict(U)
        return np.concatenate([self.model.predict(U[i:i + chunk])
                               for i in range(0, len(U), chunk)])

def _q2(y, y_hat):
    return float(1.0 - np.sum((y - y_hat)**2)/np.sum((y - y.mean())**2))

def fit_surrogates(U, y, *, seed=0, n_folds=5, candidates=SURROGATE_CANDIDATES):
    """Fit a GP (constant x Matern-5/2 ARD + white noise, normalized target)
    and gradient-boosted trees on unit-cube inputs U -> y; return the one
    with the better `n_folds`-fold CV Q2, refit on all rows. The GP kernel
    is optimized ONCE on a <= GP_HYPER_N-row subsample and held fixed
    (optimizer=None) in the folds and the final fit: optimizing ARD
    hyper-parameters on thousands of rows inside CV for every metric would
    take hours. Its CV Q2 therefore carries a mild hyper-parameter leak; the
    trees' does not.

    `candidates` is a non-empty subset of SURROGATE_CANDIDATES: only those
    families are built, scored and eligible, and `Surrogate.q2` has exactly
    those keys. Dropping 'gp' skips the kernel tuning entirely -- a GP
    PREDICTS thousands of times more slowly than the trees, so the index
    estimation of a metric that does not need one is the cheap path."""
    candidates = tuple(candidates)
    unknown = [c for c in candidates if c not in SURROGATE_CANDIDATES]
    if not candidates or unknown:
        raise ValueError(f'candidates must be a non-empty subset of '
                         f'{SURROGATE_CANDIDATES}; got {candidates!r}')
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.model_selection import KFold
    U, y = np.asarray(U, dtype=float), np.asarray(y, dtype=float)
    makers = {}
    if 'gp' in candidates:
        from sklearn.exceptions import ConvergenceWarning
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
        d = U.shape[1]
        rng = np.random.default_rng(seed)
        sub = rng.choice(len(U), size=min(GP_HYPER_N, len(U)), replace=False)
        kernel = (ConstantKernel(1.0, (1e-3, 1e3))
                  *Matern(length_scale=np.ones(d), length_scale_bounds=(1e-2, 1e2), nu=2.5)
                  + WhiteKernel(1e-3, (1e-8, 1e1)))
        tuned = GaussianProcessRegressor(kernel, normalize_y=True,
                                         n_restarts_optimizer=1, random_state=seed)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', ConvergenceWarning)
            tuned.fit(U[sub], y[sub])
        makers['gp'] = lambda: GaussianProcessRegressor(tuned.kernel_, optimizer=None,
                                                        normalize_y=True)
    if 'hgb' in candidates:
        makers['hgb'] = lambda: HistGradientBoostingRegressor(
            max_iter=500, learning_rate=0.05, early_stopping=True,
            random_state=seed)
    folds = list(KFold(n_folds, shuffle=True, random_state=seed).split(U))
    q2 = {}
    for name, make in makers.items():
        y_hat = np.empty_like(y)
        for train, test in folds:
            y_hat[test] = make().fit(U[train], y[train]).predict(U[test])
        q2[name] = _q2(y, y_hat)
    name = max(q2, key=q2.get)
    return Surrogate(name=name, model=makers[name]().fit(U, y), q2=q2,
                     reliable=q2[name] >= RELIABLE_Q2)
