#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Publication figure comparing the 15 decision variables of the
`metabolic_minimal_subset` kinetic-optimization preset across the
scenario-A baseline and a handful of hand-picked campaign trials (each
trial optionally from a different campaign -- one trajectory CSV per
objective). Three stacked panels:

  a  outcome incumbent trajectories -- IRR, ethanol titer, isobutanol
     titer, batch time. One step line per set: the metric evaluated at
     that set's running optimization incumbent vs trial_number, where the
     incumbent is the arg-best of the set's own selection metric (its
     objective for a "best" trial, COL for "best:COL"). So the IRR panel
     shows best-so-far IRR for an IRR study and the at-best-titer IRR for
     a titer study. The scenario-A baseline (no study) is a dashed
     reference line.
  b  the 15 parameter values, grouped by nskinetics role and arranged in
     pathway order (glycolysis/fermentation r1->r3->r6, Ehrlich branch
     r13->r16, product-inhibition effector multipliers, feeding), one
     bar per set with the searched band shaded.
  c  enzyme burden -- one horizontal stacked bar per set with the SAME
     five fixed pathway categories on every bar (a partition of the
     Phi_M pools: glycolysis r1; TCA cycle r2; acetate / acetyl-CoA
     production r4->r5; ethanol production r3+r6; isobutanol production
     r13->r16), plus the translation sector phi_T as a dotted final-
     category tail; the F_flex cap and the growth-derating factor d are
     marked.

"campaign" is this figure's word for a kinetic-optimization study; code,
CSV, and study names keep "study".

Sim-safe: kinetic_optimization.py (ko) and enzyme_burden.py (eb) are
loaded by file path (no biosteam import, no load()); campaign CSVs and
the parameter-distribution workbooks are plain pandas reads. Runnable
while a campaign is in flight. Run:

    python plots/plot_kin_opt_parameter_sets.py \
        --set "Best IRR" <campaign> best \
        --set "Best ethanol titer" <campaign> "best:EtOH titer" \
        --set "Best isobutanol titer" <campaign> "best:IBO titer"

With no --set arguments it plots the "best" trial of each of the five most
recent minimal-subset campaigns -- one per objective (IRR / ethanol titer /
isobutanol titer / ethanol yield / isobutanol yield) -- against the
baseline. Writes <stem>_<stamp>.png and .pdf to --out-dir.
"""
import os
import argparse
import importlib.util
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.ticker import (AutoMinorLocator, FixedLocator, FuncFormatter,
                               LogLocator, MultipleLocator, NullFormatter,
                               NullLocator, PercentFormatter)

# --- sim-safe module loads (by file path; never import the package) ----------
PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')
# The three most recent minimal-subset studies share one search space
# (rb0.001-10_ib0.2-2_burden) and differ only in the optimized objective,
# so each set is that study's own "best" trial (max of its objective).
_MINIMAL_SUBSET_STUDY = ('kin_opt_ethanol_isobutanol_metabolic_minimal_subset'
                         '_%s_rb0.001-10_ib0.2-2_burden')
DEFAULT_STUDY = _MINIMAL_SUBSET_STUDY % 'irr'          # financial (IRR) optimum
ETOH_TITER_STUDY = _MINIMAL_SUBSET_STUDY % 'etoh_titer'  # ethanol-titer optimum
IBO_TITER_STUDY = _MINIMAL_SUBSET_STUDY % 'ibo_titer'    # isobutanol-titer optimum
IBO_YIELD_STUDY = _MINIMAL_SUBSET_STUDY % 'ibo_yield'    # isobutanol-yield optimum
ETOH_YIELD_STUDY = _MINIMAL_SUBSET_STUDY % 'etoh_yield'  # ethanol-yield optimum


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(PKG_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ko = _load('ko', 'kinetic_optimization.py')
eb = _load('eb', 'enzyme_burden.py')

# --- the 15 decision variables, in figure order -----------------------------
RATE_VARS = list(ko.METABOLIC_MINIMAL_SUBSET_RATES)          # 9
# Ehrlich-branch capacities: genuinely zero (branch off) at the scenario-A
# baseline, so their baseline "0" bar labels are suppressed (see bar_cell)
EHRLICH_RATE_VARS = ('k_13', 'k_14', 'k_15', 'k_16')
GROUP_VARS = list(ko.METABOLIC_MINIMAL_SUBSET_GROUPS)        # 3
FEED_VARS = ['threshold_conc', 'target_delta', 'max_n_spikes']  # 3 CSV columns
DECISION_VARS = RATE_VARS + GROUP_VARS + FEED_VARS
# feeding cells DRAW derived / realized quantities in place of two raw decision
# columns: the applied target sugar concentration (target_conc = min(
# TARGET_CONC_MAX, threshold_conc + target_delta), via ko._applied_feeding) for
# target_delta, and the realized spike count (n_glu_spikes, a tracked metric)
# for the max_n_spikes cap.
FEED_DRAW_VARS = ['target_conc', 'threshold_conc', 'n_glu_spikes']

BANDS = [
    ('Glycolysis (r1) + ethanol production (r3 → r6)',
     ['k_1l', 'k_1h', 'k_1e', 'k_3', 'k_6']),
    ('Isobutanol production (Ehrlich pathway, r13 → r14 → r15 → r16)',
     ['k_13', 'k_14', 'k_15', 'k_16']),
    ('Product inhibition relative to baseline '
     '(applied to r1, r4, r6, r7, r10, r16)',
     GROUP_VARS),
    ('Feeding strategy', FEED_DRAW_VARS),
]

# cell titles: parameter on the first line, reaction index + enzyme below
REACTION_LABELS = {
    'k_1l': 'k_1l\nr1 glycolysis\n(low affinity)',
    'k_1h': 'k_1h\nr1 glycolysis\n(high affinity)',
    'k_1e': 'k_1e\nr1 glycolysis\n(ethanol phase)',
    'k_3': 'k_3\nr3 PDC\n(Pdc1)',
    'k_6': 'k_6\nr6 ADH\n(Adh1)',
    'k_13': 'k_13\nr13 ALS\n(Ilv2+Ilv6)',
    'k_14': 'k_14\nr14 KARI\n(Ilv5)',
    'k_15': 'k_15\nr15 DHAD\n(Ilv3)',
    'k_16': 'k_16\nr16 KDC+ADH\n(Aro10+Adh6)',
}
GROUP_LABELS = {'inhib_ethanol': 'Ethanol\ninhibition',
                'inhib_isobutanol': 'Isobutanol\ninhibition',
                'inhib_acetate': 'Acetate\ninhibition'}
# value-axis units. Rate constants: every sampled capacity carries the
# Antimony unit g_per_l_per_h (gram/(litre*hour)) in the shipped model, so
# g·L^-1·h^-1 for all nine. Effector multipliers are dimensionless fold-changes
# relative to the per-family baseline (1 = baseline).
RATE_UNIT = 'g·L$^{-1}$·h$^{-1}$'
# feeding cell: (title, (shaded-range low, high)) -- engine default bounds
# target_conc shares the threshold cell's 0-300 g/L axis (it is clipped at
# TARGET_CONC_MAX = 300, and always sits at or above the threshold), so the two
# feeding concentrations read on the same scale.
FEED_LABELS = {'threshold_conc': ('Thresh. sugar\nconc. [g·L$^{-1}$]', (0, 300)),
               'target_conc': ('Target sugar\nconc. [g·L$^{-1}$]', (0, 300)),
               'n_glu_spikes': ('No. of spikes', (0, 50))}

# the seven study steps -> enzyme name and charging parameter(s); read
# against the eb tables so a table drift here raises at import
STEP_ENZYME = {'r1': 'Glycolysis lump', 'r3': 'Pdc1', 'r6': 'Adh1',
               'r13': 'Ilv2+Ilv6', 'r14': 'Ilv5', 'r15': 'Ilv3',
               'r16': 'Aro10+Adh6'}
STEP_PARAMS = {'r1': 'k_1l, k_1h, k_1e', 'r3': 'k_3', 'r6': 'k_6',
               'r13': 'k_13', 'r14': 'k_14', 'r15': 'k_15', 'r16': 'k_16'}
for _s in STEP_ENZYME:
    if _s not in eb.NATIVE_STEPS and _s not in eb.EHRLICH_STEPS:
        raise KeyError(f'STEP_ENZYME step {_s!r} not in eb tables')

# outcomes panel: (CSV column, cell title, (y-low, y-high)), in draw order.
# Titles are rotated y-axis labels in narrow cells -- keep each to <=2 lines
# so the label stays out of the neighbouring cell's plot box.
OUTCOMES = (('IRR', 'Financial attractiveness\nas IRR [%]', (0, 0.3)),
            ('IBO titer', 'Isobutanol titer\n[g·L$^{-1}$]', (0, 100)),
            ('IBO yield', 'Isobutanol yield\n[g·g$^{-1}$]', (0, 0.4)),
            ('EtOH titer', 'Ethanol titer\n[g·L$^{-1}$]', (0, 300)),
            ('EtOH yield', 'Ethanol yield\n[g·g$^{-1}$]', (0, 0.5)))

# per-outcome override for the major-tick step (else _linear_cap's step).
# IRR is stored as a fraction shown in percent, so 0.05 -> 0/5/.../30 %;
# isobutanol titer reads cleaner on 0/40/80/120 than the default 0/20/.../120.
OUTCOME_TICK_STEP = {'IRR': 0.05, 'IBO titer': 40.0}

# outcomes whose negative values are drawn at zero: a loss-making IRR (finite
# negative) and an unsolvable IRR (-inf) both read as 0 rather than diving
# below the axis floor / being omitted. Only nan (never solved) stays
# non-finite and is omitted.
CLAMP_NEG_TO_ZERO = {'IRR'}

# one color per set: baseline dark grey, campaigns from the hue palette
BASELINE_COLOR = '0.25'
HUE_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#8c564b']
MAX_SETS = 1 + len(HUE_COLORS)   # 6

FONTS = {'band': 12, 'cell': 10, 'tick': 9, 'callout': 9,
         'legend': 10, 'axis': 11, 'panel': 14}

# --- scenario-A baseline outcomes: HARD-CODED (spec decision b).
# Simulated 2026-09-07 (smoke_test_1 protocol, IBO_2026). A cached one-off
# baseline simulation (option a) is a later change; when done, replace this
# dict with the cached read.
BASELINE_A = {
    'IRR': 0.1230, 'EtOH titer': 118.4, 'IBO titer': 0.0, 'tau': 45.4,
    'n_glu_spikes': 10, 'EtOH yield': 0.455, 'IBO yield': 0.0,
    'Cell density': 15.7, 'TCI': 139.6, 'threshold_conc': 217.125,
    'target_delta': 4.125, 'max_n_spikes': 16,
}


# --- data assembly -----------------------------------------------------------
def baseline_set():
    """Scenario-A baseline as a flat record shaped like a campaign row.

    Kinetics: the scenario-A workbook rows (the four Ehrlich rates are
    absent from the A workbook == zero in the live A model); effector
    multipliers 1.0; feeding from BASELINE_A. Burden pools from
    eb.BurdenModel at that point (r13-r16 == 0). Outcomes are the
    hard-coded BASELINE_A constants (spec decision b).
    """
    A = ko.workbook_kinetic_baselines('A')
    # k_ref for the burden model: every required capacity, Ehrlich at 0
    k_ref = dict(A)
    for k in ('k_13', 'k_14', 'k_15', 'k_16'):
        k_ref.setdefault(k, 0.0)
    missing = [c for c in eb.BurdenModel.required_capacities()
               if c not in k_ref]
    if missing:
        raise KeyError('scenario-A workbook is missing burden capacities '
                       f'{missing}')
    res = eb.BurdenModel(k_ref).evaluate(k_ref)

    rec = {'label': 'Baseline', 'campaign': None,
           'trial_number': None, 'is_baseline': True, 'extra_sampled': []}
    for k in RATE_VARS:
        rec[k] = float(k_ref[k])
    for g in GROUP_VARS:
        rec[g] = 1.0
    rec['threshold_conc'] = BASELINE_A['threshold_conc']
    rec['target_delta'] = BASELINE_A['target_delta']
    rec['max_n_spikes'] = BASELINE_A['max_n_spikes']
    # applied target sugar concentration (clip at TARGET_CONC_MAX), drawn
    # in place of target_delta
    rec['target_conc'] = float(ko._applied_feeding(rec)[0])
    for col in ('IRR', 'TCI', 'EtOH titer', 'IBO titer', 'EtOH yield',
                'IBO yield', 'tau', 'n_glu_spikes'):
        rec[col] = BASELINE_A[col]
    for st, pool in res.pools.items():
        rec[f'pool_{st}'] = float(pool)
    rec['Phi_M'] = float(res.Phi_M)
    rec['phi_T'] = float(res.phi_T)
    rec['F_flex'] = float(res.F_flex)
    rec['burden_factor'] = float(res.burden_factor)
    # no study -> no incumbent trajectory; drawn as a dashed reference line
    rec['sel_col'] = None
    rec['sel_dir'] = None
    rec['traj_x'] = None
    rec['traj'] = None
    return rec


# objective slugs, longest first so "IBO yield x titer" beats "IBO yield"
_OBJ_SLUGS = sorted(ko.OBJECTIVE_REGISTRY, key=len, reverse=True)


def resolve_campaign_csv(campaign):
    """A study name -> RESULTS_DIR/<name>_trajectory.csv, or a literal
    path to a trajectory CSV used as-is."""
    if campaign.lower().endswith('.csv'):
        path = campaign
    else:
        path = os.path.join(RESULTS_DIR, f'{campaign}_trajectory.csv')
    if not os.path.isfile(path):
        raise FileNotFoundError(f'campaign trajectory CSV not found: {path}')
    return path


def campaign_objective(campaign):
    """The objective slug embedded in a study name (the token the driver
    slugged in), or None. Matches an OBJECTIVE_REGISTRY key with spaces
    turned to underscores, e.g. ..._minimal_subset_irr_... -> 'IRR'."""
    stem = os.path.basename(campaign)
    if stem.endswith('.csv'):
        return None
    low = stem.lower()
    for obj in _OBJ_SLUGS:
        if f'_{obj.lower().replace(" ", "_")}_' in low or \
           low.endswith('_' + obj.lower().replace(' ', '_')):
            return obj
    return None


def selection_spec(df, trial, campaign, warn=False):
    """(selection column, direction) defining the running incumbent for a
    set -- the same rule resolve_trial uses to pick its representative row:
    "best:COL" -> (COL, 'maximize'); "best" and an explicit integer trial
    (a placeholder while its study runs) -> the 'objective' column with the
    campaign objective's registry direction ('maximize' if the slug is
    unknown). Factored out so the representative row (resolve_trial) and the
    panel-a incumbent trajectory (incumbent_trajectory) cannot drift."""
    if isinstance(trial, str) and trial.startswith('best') and ':' in trial:
        col = trial.split(':', 1)[1]
        if col not in df.columns:
            raise ValueError(
                f'campaign {campaign}: no metric column {col!r} for '
                f'"best:{col}"')
        return col, 'maximize'
    obj = campaign_objective(campaign)
    if obj is None or obj not in ko.OBJECTIVE_REGISTRY:
        if warn:
            print(f'  WARNING campaign {campaign}: objective slug not in the '
                  'registry; "best" maximizes the objective column')
        return 'objective', 'maximize'
    return 'objective', ko.OBJECTIVE_REGISTRY[obj]['direction']


def resolve_trial(df, trial, campaign):
    """Return the requested COMPLETE row as a Series."""
    ok = df[df['state'] == 'COMPLETE']
    if isinstance(trial, str) and trial.startswith('best'):
        sel_col, direction = selection_spec(df, trial, campaign, warn=True)
        idx = ok[sel_col].idxmin() if direction == 'minimize' \
            else ok[sel_col].idxmax()
        return ok.loc[idx]
    # explicit integer trial_number
    n = int(trial)
    hit = df[df['trial_number'] == n]
    if hit.empty:
        raise ValueError(f'campaign {campaign}: trial {n} not found')
    row = hit.iloc[0]
    if row['state'] != 'COMPLETE':
        raise ValueError(f'campaign {campaign}: trial {n} is '
                         f'{row["state"]}, not COMPLETE (no outcomes/pools)')
    return row


def incumbent_trajectory(df, sel_col, direction):
    """Running-incumbent trajectory over a study's COMPLETE trials.

    Walking COMPLETE trials in trial_number order, keep the arg-best of
    `sel_col` under `direction` (non-finite selection values never win;
    before the first finite one the incumbent is undefined). Return
    (x, {outcome_col: y}): x is the trial_number of every COMPLETE trial
    (plus the study's final trial_number so the step reaches the right
    edge), and each y is that outcome column evaluated at the running
    incumbent -- non-finite -> NaN, so the line breaks rather than diving
    to a floor (e.g. an unsolvable-IRR trial that is a titer incumbent)."""
    ok = df[df['state'] == 'COMPLETE'].sort_values('trial_number')
    n = len(ok)
    x = ok['trial_number'].to_numpy(dtype=float)
    sel = ok[sel_col].to_numpy(dtype=float)
    inc = np.full(n, -1, dtype=int)
    best_pos, best_val = -1, None
    for i in range(n):
        v = sel[i]
        if np.isfinite(v) and (best_pos < 0 or
                               (v > best_val if direction == 'maximize'
                                else v < best_val)):
            best_pos, best_val = i, v
        inc[i] = best_pos
    mask = inc >= 0
    traj = {}
    for col, _, _ in OUTCOMES:
        vals = ok[col].to_numpy(dtype=float)
        y = np.full(n, np.nan)
        y[mask] = vals[inc[mask]]
        y[~np.isfinite(y)] = np.nan
        traj[col] = y
    if n:
        x_end = float(df['trial_number'].max())
        if x_end > x[-1]:
            x = np.append(x, x_end)
            for col in traj:
                traj[col] = np.append(traj[col], traj[col][-1])
    return x, traj


def load_set(label, campaign, trial):
    """A campaign trial as a flat record (CSV row + metadata + the panel-a
    incumbent trajectory)."""
    df = ko.load_trajectory(resolve_campaign_csv(campaign))
    missing = [c for c in DECISION_VARS if c not in df.columns]
    if missing:
        raise ValueError(
            f'campaign {campaign}: not a metabolic_minimal_subset campaign '
            f'(missing decision columns {missing}). A metabolic / '
            'metabolic_protein campaign samples the inhibition coefficients '
            'individually and has no inhib_* group columns.')
    missing_b = [c for c in eb.BURDEN_COLUMNS if c not in df.columns]
    if missing_b:
        raise ValueError(f'campaign {campaign}: missing burden columns '
                         f'{missing_b}')
    row = resolve_trial(df, trial, campaign)
    rec = row.to_dict()
    rec['label'] = label
    rec['campaign'] = campaign
    rec['trial_number'] = int(row['trial_number'])
    rec['is_baseline'] = False
    # applied target sugar concentration (clip at TARGET_CONC_MAX), drawn
    # in place of the raw target_delta decision column
    rec['target_conc'] = float(ko._applied_feeding(rec)[0])
    # decision columns beyond the 15 (e.g. a metabolic_minimal campaign's
    # extra rates + stage_1_max_x); reported, not drawn
    known = set(DECISION_VARS) | {'trial_number', 'state', 'objective',
                                  'error'}
    known |= set(eb.BURDEN_COLUMNS)
    known |= set(ko.TRACKED_METRICS)
    known |= {c for c in df.columns if c.startswith('applied_')}
    rec['extra_sampled'] = [c for c in df.columns if c not in known]
    sel_col, sel_dir = selection_spec(df, trial, campaign)
    tx, traj = incumbent_trajectory(df, sel_col, sel_dir)
    rec['sel_col'] = sel_col
    rec['sel_dir'] = sel_dir
    rec['traj_x'] = tx
    rec['traj'] = traj
    # full per-trial cloud for panel a, drawn only on the cell whose metric
    # this campaign optimized: every trial's number, its completion flag and
    # each outcome value (non-COMPLETE / unsolved -> NaN via the CSV).
    rec['objective'] = campaign_objective(campaign)
    rec['scatter_x'] = df['trial_number'].to_numpy(dtype=float)
    rec['scatter_complete'] = (df['state'] == 'COMPLETE').to_numpy()
    rec['scatter'] = {c: df[c].to_numpy(dtype=float) for c, _, _ in OUTCOMES}
    return rec


# --- searched-band shading, from the driver's own preset -> search-space
# construction -----------------------------------------------------------
_DEFAULT_AXES = ('ethanol_isobutanol', 'metabolic_minimal_subset')
_TARGET_PRODUCTS = ('ethanol_isobutanol', 'ethanol_only')
_STUDY_TYPES = ('metabolic_minimal_subset', 'metabolic_minimal',
                'metabolic_protein', 'metabolic')   # longest-first


def campaign_axes_from_name(campaign):
    """(study_target_products, study_type) parsed from a study name;
    default to the minimal-subset ethanol+isobutanol preset (with a
    warning) for a bare CSV path or an unrecognised name."""
    stem = os.path.basename(campaign)
    if stem.endswith('.csv') or 'kin_opt_' not in stem:
        print(f'  WARNING campaign {campaign}: band source not inferrable '
              'from the name; assuming the metabolic_minimal_subset preset')
        return _DEFAULT_AXES
    tp = next((t for t in _TARGET_PRODUCTS if t in stem), _DEFAULT_AXES[0])
    ty = next((t for t in _STUDY_TYPES if t in stem), _DEFAULT_AXES[1])
    return tp, ty


def campaign_band(campaign):
    """{var: (lo, hi)} searched band for the 15 decision vars, exactly as
    the driver built it: ko.resolve_study_preset -> ko.build_search_space
    on the scenario-A baselines. Feeding vars from the engine defaults.

    ko.build_search_space returns (space, excluded_parameter_names), and
    it needs a `kinetic_baselines` dict that already carries a value for
    EVERY name it should place in the space -- it iterates
    kinetic_baselines.items(), never preset['include_params'] on its own
    -- and, for a `parameter_groups` member, a value that is POSITIVE
    (the group-forming loop checks kinetic_baselines[member] directly,
    before param_bounds_override is even consulted). The scenario-A
    workbook (ko.workbook_kinetic_baselines('A')) alone is not enough:
    - the Ehrlich-branch rate constants (k_13..k_16, role `capacity`)
      are absent from it because they are genuinely ZERO in the live
      scenario-A model (the branch is off); a zero/absent baseline with
      no override is EXCLUDED by build_search_space, so the driver's own
      real run (which discovers every parameter off the LIVE model,
      still 0 there) supplies a `param_bounds_override` derived from the
      scenario-B workbook instead (`kinetic_bounds_from_scenario` in
      analyses/optimize_kinetics_BO.py, ko.workbook_kinetic_bounds here)
      so those rates still get a positive, B-baseline-derived band.
    - the inhibition-coefficient group members (k_1ie, k_4ie, ...) are
      likewise absent from the A workbook, but -- unlike the Ehrlich
      rates -- they are NOT zero in the live A model: A's workbook simply
      doesn't curate them, so they carry the shipped Antimony default,
      which is numerically the scenario-B workbook value (see
      plot_kin_opt_best_vs_baseline_heatmap.py's `build_rows`, same
      fallback). A group member needs that positive value in
      `kinetic_baselines` itself, not just an override.
    Every name shared by both workbooks (k_1l, k_1h, k_1e, k_3, k_6) has
    the identical baseline in each, so folding in the B-derived override
    changes nothing for those -- the shared-baseline rate band is exact.
    """
    tp, ty = campaign_axes_from_name(campaign)
    preset = ko.resolve_study_preset(tp, ty)
    kb_scenario = preset['kinetic_bounds_scenario']
    A = ko.workbook_kinetic_baselines('A')
    A_bounds_wb = ko.workbook_kinetic_baselines(kb_scenario)
    kinetic_baselines = dict(A)
    group_members = {m for members in (preset['parameter_groups'] or {}).values()
                     for m in members}
    for name in group_members:
        if name not in kinetic_baselines:
            # inherited (nonzero) Antimony-default value; see docstring
            kinetic_baselines[name] = A_bounds_wb[name]
    for name in preset['include_params']:
        if name not in kinetic_baselines:
            # genuinely zero in the live scenario-A model (e.g. k_13..k_16)
            kinetic_baselines[name] = 0.0
    # the scenario-B-derived absolute bounds that let a zero-baseline
    # Ehrlich rate stay in the space instead of being excluded; identical
    # in value to the A-baseline band for every name the two workbooks share
    param_bounds_override = ko.workbook_kinetic_bounds(
        kb_scenario,
        multiplier_bounds=preset['multiplier_bounds'],
        rate_multiplier_bounds=preset['rate_multiplier_bounds'],
        rate_params=preset['rate_params'],
        parameter_multiplier_bounds=preset['parameter_multiplier_bounds'],
    )
    space, _excluded = ko.build_search_space(
        kinetic_baselines,
        multiplier_bounds=preset['multiplier_bounds'],
        param_bounds_override=param_bounds_override,
        include_params=preset['include_params'],
        exclude_params=preset['exclude_params'],
        rate_multiplier_bounds=preset['rate_multiplier_bounds'],
        rate_params=preset['rate_params'],
        parameter_multiplier_bounds=preset['parameter_multiplier_bounds'],
        parameter_groups=preset['parameter_groups'],
        group_multiplier_bounds=preset['group_multiplier_bounds'],
        spike_delta_bounds=preset['spike_delta_bounds'],
        stage_1_max_x_bounds=preset['stage_1_max_x_bounds'],
    )
    band = {}
    for v in DECISION_VARS:
        if v in space:
            band[v] = (space[v]['low'], space[v]['high'])
        elif v in GROUP_VARS:
            # fallback, in case a future preset groups differently and the
            # group itself doesn't end up as its own entry in `space`
            band[v] = tuple(preset['group_multiplier_bounds'])
    # feeding fallbacks (build_search_space always emits these, but be safe)
    band.setdefault('threshold_conc', (0.0, 300.0))
    band.setdefault('target_delta', (5.0, 500.0))
    band.setdefault('max_n_spikes', (0, 50))
    return band


# --- drawing -----------------------------------------------------------------
def apply_fonts():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['xtick.labelsize'] = FONTS['tick']
    plt.rcParams['ytick.labelsize'] = FONTS['tick']
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['hatch.linewidth'] = 0.6
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['mathtext.bf'] = 'Arial:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


def set_colors(sets):
    colors, hue = {}, 0
    for s in sets:
        if s.get('is_baseline'):
            colors[id(s)] = BASELINE_COLOR
        else:
            if hue >= len(HUE_COLORS):
                n_campaign = len([x for x in sets if not x.get('is_baseline')])
                raise ValueError(
                    f'too many campaign sets ({n_campaign}); at most '
                    f'{len(HUE_COLORS)} plus the baseline (six total)')
            colors[id(s)] = HUE_COLORS[hue]
            hue += 1
    return colors


def style_cell_axes(ax):
    ax.tick_params(axis='y', which='major', direction='inout', right=False,
                   length=4)
    ax.tick_params(axis='y', which='minor', direction='inout', right=False,
                   length=2.2)
    ax.tick_params(axis='x', which='both', top=False, bottom=False,
                   labelbottom=False)
    for sp in ('right', 'top'):
        ax.spines[sp].set_visible(False)


def _baseline_value(sets, var):
    for s in sets:
        if s.get('is_baseline'):
            return s.get(var)
    return None


# --- nice, tick-aligned axis caps: every quantitative axis is bounded on
# both ends by a labeled tick (like the IRR cell's 0 -> 30) ------------------
_NICE_MANTISSAS = (1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)


def _nice_ceiling(x):
    """Smallest 'nice' number (mantissa in _NICE_MANTISSAS) >= x."""
    if x <= 0:
        return 1.0
    e = np.floor(np.log10(x))
    scale = 10.0 ** e
    m = x / scale
    for cand in _NICE_MANTISSAS:
        if m <= cand * (1 + 1e-9):
            return cand * scale
    return 10.0 * scale


def _linear_cap(vmax, floor_hi=0.0):
    """(hi, step) for a linear axis based at 0 with both 0 and hi as major
    ticks. hi is the nice ceiling of max(vmax, floor_hi); step = hi / k for
    the fewest-tick k in (3, 4, 5, 6, 8) that lands on a clean 1/2/2.5/5
    mantissa (falling back to hi / 5). hi is an integer multiple of step, so
    a MultipleLocator(step) labels both ends."""
    hi = _nice_ceiling(max(vmax, floor_hi))
    for k in (3, 4, 5, 6, 8):
        step = hi / k
        e = np.floor(np.log10(step))
        m = step / 10.0 ** e
        if min(abs(m - t) for t in (1.0, 2.0, 2.5, 5.0)) < 0.03:
            return hi, step
    return hi, hi / 5.0


def _log_stop_below(x):
    """Largest 1/2/5-decade stop strictly below x (e.g. 0.2 -> 0.1)."""
    e = np.floor(np.log10(x))
    for m in (5.0, 2.0, 1.0):
        v = m * 10.0 ** e
        if v < x * (1 - 1e-9):
            return v
    return 10.0 ** (e - 1) * 5.0


def _log_stop_above(x):
    """Smallest 1/2/5-decade stop strictly above x (e.g. 2 -> 5)."""
    e = np.floor(np.log10(x))
    for m in (1.0, 2.0, 5.0):
        v = m * 10.0 ** e
        if v > x * (1 + 1e-9):
            return v
    return 10.0 ** (e + 1)


def _decade_floor(lo):
    """Enclosing decade at or below lo, dropped one more when the band would
    otherwise touch the bottom edge (keeps ~>=0.3 decade headroom for the
    '0'/absent bar labels); the result is a decade, so it carries a tick."""
    a = np.floor(np.log10(lo))
    if 10.0 ** a > lo * 0.5:
        a -= 1
    return 10.0 ** a


def _decade_ceil(hi):
    """Enclosing decade at or above hi, raised one more when the band would
    otherwise touch the top edge; the result is a decade (a tick)."""
    b = np.ceil(np.log10(hi))
    if 10.0 ** b < hi * 2:
        b += 1
    return 10.0 ** b


# hand-tuned major y-ticks for a few rate cells (readability): an explicit
# tick list overrides the auto _linear_cap ticks and sets the axis top to its
# largest tick. A uniform list keeps half-step minor ticks (so e.g. k_14's
# dropped 0.5/1.5 survive as unlabeled minors); a non-uniform list (k_1h)
# drops the removed values outright.
RATE_YTICKS = {
    'k_1h': [0.0, 0.5, 1.0],        # even 0.5 steps (drop 0.25, 0.75)
    'k_14': [0.0, 1.0, 2.0],        # drop 0.5, 1.5
    'k_15': [0.0, 1.0, 2.0],        # drop 0.5, 1.5
    'k_16': [0.0, 1.0, 2.0, 3.0],   # re-cap 2.5 -> 3.0, whole-number steps
}


def bar_cell(ax, sets, colors, var, kind, ylabel, subtitle=None, ylim=None,
             band=None):
    """One parameter cell: colored bars per set with the searched band
    shaded and the baseline dashed. `ylabel` names the value axis (the
    parameter symbol / group / feeding quantity); `subtitle`, if given, is
    a smaller descriptor above the cell (reaction + enzyme). Value axes get
    conventional linear ticks (major + one minor between) on the rate,
    feeding and group cells alike."""
    n = len(sets)
    base = _baseline_value(sets, var)
    if kind == 'rate':
        # linear value axis from 0 to a nice ceiling just above the largest
        # bar (both edges on labeled ticks); no search-band shading
        data_max = 0.0
        for s in sets:
            v = s.get(var)
            if v is not None and np.isfinite(v):
                data_max = max(data_max, float(v))
        if base and base > 0:
            data_max = max(data_max, base)
        bottom = 0
        if var in RATE_YTICKS:
            ticks = RATE_YTICKS[var]
            ax.set_ylim(bottom, ticks[-1])
            ax.yaxis.set_major_locator(FixedLocator(ticks))
            uniform = len(set(np.round(np.diff(ticks), 6))) == 1
            ax.yaxis.set_minor_locator(
                AutoMinorLocator(2) if uniform else NullLocator())
        else:
            cap, step = _linear_cap(data_max)
            ax.set_ylim(bottom, cap)
            ax.yaxis.set_major_locator(MultipleLocator(step))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        if base and base > 0:
            ax.axhline(base, color='k', lw=0.8, ls='--', zorder=1)
    elif kind == 'group':
        # linear value axis from 0 to a nice ceiling just above the largest
        # bar (both edges on labeled ticks), baseline 1.0 dashed; no
        # search-band shading -- one minor tick between majors like the rates
        data_max = 0.0
        for s in sets:
            v = s.get(var)
            if v is not None and np.isfinite(v):
                data_max = max(data_max, float(v))
        data_max = max(data_max, 1.0)   # keep the baseline reference framed
        bottom = 0
        cap, _ = _linear_cap(data_max)
        # exactly three major ticks (0, cap/2, cap) with one minor between,
        # e.g. 0 / 0.5 / 1.0 on the ethanol cell, 0 / 1.0 / 2.0 elsewhere
        ax.set_ylim(bottom, cap)
        ax.yaxis.set_major_locator(FixedLocator([0.0, cap / 2.0, cap]))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.axhline(1.0, color='k', lw=0.8, ls='--', zorder=1)
    else:  # feed
        lo, _ = ylim
        # linear axis from the declared floor to a nice ceiling just above the
        # largest bar (both edges on labeled ticks); no search-band shading
        data_max = 0.0
        for s in sets:
            v = s.get(var)
            if v is not None and np.isfinite(v):
                data_max = max(data_max, float(v))
        if base is not None and np.isfinite(base):
            data_max = max(data_max, base)
        cap, step = _linear_cap(data_max)
        ax.set_ylim(lo, cap)
        if base is not None:
            ax.axhline(base, color='k', lw=0.8, ls='--', zorder=1)
        ax.yaxis.set_major_locator(MultipleLocator(step))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        bottom = 0
    ax.set_ylabel(ylabel, fontsize=FONTS['cell'], labelpad=3)
    if subtitle:
        ax.set_title(subtitle, fontsize=FONTS['tick'], pad=4)
    for j, s in enumerate(sets):
        v = s.get(var)
        c = colors[id(s)]
        if v is not None and np.isfinite(v) and v > 0:
            ax.bar(j, v - bottom, bottom=bottom, width=0.72, color=c, zorder=2)
        else:
            top = ax.get_ylim()[1]
            y = bottom * 1.25 if bottom else 0.01 * top
            label = 'n/a' if (v is None or (isinstance(v, float)
                                            and not np.isfinite(v))) else '0'
            # the baseline's Ehrlich branch is genuinely off (zero flux) and a
            # batch set runs zero glucose spikes; those "0" labels only clutter
            # the four Ehrlich cells and the spikes cell, so drop them
            drop = label == '0' and (
                (s.get('is_baseline') and var in EHRLICH_RATE_VARS)
                or var == 'n_glu_spikes')
            if not drop:
                ax.text(j, y, label, ha='center', va='bottom',
                        fontsize=FONTS['tick'], color=c)
    ax.set_xlim(-0.6, n - 0.4); style_cell_axes(ax)


def outcome_cell(ax, sets, colors, col, title, ylim, xmax):
    """One outcome metric as incumbent trajectories over trial_number: one
    step line per campaign set (the metric at that set's running incumbent,
    in its panel-b/c color), the scenario-A baseline as a dashed reference.
    """
    lo, hi = ylim

    def _yv(y):
        # every value is drawn at the axis floor lo when it has no finite
        # position: failures (no finite value: nan / +-inf) drop to lo instead
        # of being omitted (cloud) or gapping (line), and for the flagged
        # outcomes (IRR) finite losses (< lo) are drawn at lo too. So both the
        # cloud and the incumbent line dip to zero on a failure, no hole.
        y = np.asarray(y, dtype=float)
        out = np.where(np.isfinite(y), y, lo)
        if col in CLAMP_NEG_TO_ZERO:
            out = np.where(out < lo, lo, out)
        return out

    base = _baseline_value(sets, col)
    finite = [v for s in sets if s.get('traj') is not None
              for v in s['traj'][col] if np.isfinite(v)]
    if base is not None and np.isfinite(base):
        finite.append(base)
    vmax = max(finite) if finite else hi
    # bound the value axis top and bottom by labeled ticks (like IRR's 0->30):
    # a nice ceiling >= the data, no ~8% auto-overshoot past the last tick
    cap, step = _linear_cap(vmax, floor_hi=hi)
    step = OUTCOME_TICK_STEP.get(col, step)
    ax.set_ylim(lo, cap)
    ax.set_xlim(0, xmax * 1.02)
    if base is not None and np.isfinite(base):
        ax.axhline(base, color=BASELINE_COLOR, lw=0.9, ls='--', zorder=1)
    # individual trial cloud from the ONE campaign that optimized THIS metric:
    # every trial as a translucent dot in the campaign's own color, behind the
    # incumbent lines (zorder 1). Failed / unsolved trials (no finite value)
    # are drawn at the axis floor rather than omitted.
    for s in sets:
        if s.get('scatter') is None or s.get('objective') != col:
            continue
        sx, sy, cc = s['scatter_x'], _yv(s['scatter'][col]), colors[id(s)]
        ax.scatter(sx, sy, s=9, color=cc, alpha=0.25, linewidths=0, zorder=1)
    for s in sets:
        if s.get('traj') is None:
            continue
        ax.step(s['traj_x'], _yv(s['traj'][col]), where='post',
                color=colors[id(s)], lw=1.4, zorder=2)
    # metric name next to the value axis itself (not a title above the cell)
    ax.set_ylabel(title, fontsize=FONTS['cell'], labelpad=3)
    ax.set_xlabel('Trial', fontsize=FONTS['tick'], labelpad=2)
    # the wide IRR cell fits five majors (0..2000 by 500); the narrow cells
    # take three (0, 1000, 2000). Four minor ticks sit between each major pair.
    ax.xaxis.set_major_locator(MultipleLocator(500 if col == 'IRR' else 1000))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(step))
    if col == 'IRR':   # fraction stored; show the value axis in percent
        ax.yaxis.set_major_formatter(
            PercentFormatter(xmax=1.0, decimals=0, symbol=''))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='y', which='major', direction='inout', right=False,
                   length=4)
    ax.tick_params(axis='y', which='minor', direction='inout', right=False,
                   length=2.2)
    ax.tick_params(axis='x', which='major', top=False, bottom=True,
                   labelbottom=True, direction='inout', length=4)
    ax.tick_params(axis='x', which='minor', top=False, bottom=True,
                   direction='inout', length=2.2)
    for sp in ('right', 'top'):
        ax.spines[sp].set_visible(False)


def draw_outcomes(fig, gs_cell, sets, colors):
    xmax = 1.0
    for s in sets:
        tx = s.get('traj_x')
        if tx is not None and len(tx):
            xmax = max(xmax, float(tx[-1]))
    # 2x4 grid: IRR (the headline outcome) fills the left 2x2 block and is the
    # largest cell; the remaining four outcomes fill the right 2x2 block, one
    # per cell (top row IBO titer / IBO yield, bottom row EtOH titer / yield).
    sub_gs = gs_cell.subgridspec(2, 4, wspace=0.62, hspace=0.62)
    axes = []
    big, rest = OUTCOMES[0], OUTCOMES[1:]
    ax_big = fig.add_subplot(sub_gs[0:2, 0:2])
    outcome_cell(ax_big, sets, colors, big[0], big[1], big[2], xmax)
    axes.append(ax_big)
    for (col, label, yl), (r, c) in zip(rest,
                                        ((0, 2), (0, 3), (1, 2), (1, 3))):
        ax = fig.add_subplot(sub_gs[r, c])
        outcome_cell(ax, sets, colors, col, label, yl, xmax)
        axes.append(ax)
    return axes


def draw_parameters(fig, gs_rows, sets, colors, band):
    axes = []
    for gs_row, (title, params) in zip(gs_rows, BANDS):
        sub_gs = gs_row.subgridspec(1, 5, wspace=0.75)
        last_ax = None
        for i, p in enumerate(params):
            ax = fig.add_subplot(sub_gs[0, i]); last_ax = ax
            axes.append(ax)
            if p in RATE_VARS:
                # reaction/enzyme descriptors go in the figure caption, not
                # above each cell -- keep only the rate symbol on the ylabel
                sym = REACTION_LABELS[p].partition('\n')[0]
                bar_cell(ax, sets, colors, p, 'rate', f'{sym}\n[{RATE_UNIT}]',
                         band=band)
            elif p in GROUP_VARS:
                # the "relative to baseline" unit lives in the band title now
                bar_cell(ax, sets, colors, p, 'group',
                         GROUP_LABELS[p], band=band)
            else:
                t, rng = FEED_LABELS[p]
                bar_cell(ax, sets, colors, p, 'feed', t, ylim=rng)
        # A band title is a header for the band BELOW it, so it should hug its
        # own cells and leave the larger gap to the band above. The feeding
        # cells carry long two-line vertical ylabels ("Threshold sugar
        # conc. ...") that overflow the short cell and poke above its top box,
        # so that band takes a middling offset to clear them; the rate and
        # group cells have short ylabels and take the tight offset that pins
        # the title to its band.
        if any(p in FEED_DRAW_VARS for p in params):
            offset = 0.022
        else:
            offset = 0.010
        # left-aligned with the plot boxes (the gridspec left margin)
        fig.text(0.083, last_ax.get_position().y1 + offset, title,
                 fontsize=FONTS['band'], fontweight='bold', va='bottom')
    return axes


def tint(color, t):
    r, g, b = to_rgb(color)
    return (r + (1 - r) * t, g + (1 - g) * t, b + (1 - b) * t)


_STUDY_STEPS = ['r1', 'r3', 'r6', 'r13', 'r14', 'r15', 'r16']
_UNSAMPLED_STEPS = ['r2', 'r4', 'r5']   # never sampled; folded into "other"

if set(_STUDY_STEPS) | set(_UNSAMPLED_STEPS) != set(eb.STEP_ORDER):
    raise AssertionError(
        'burden step partition drift: _STUDY_STEPS | _UNSAMPLED_STEPS != '
        'eb.STEP_ORDER (%r vs %r)'
        % (sorted(set(_STUDY_STEPS) | set(_UNSAMPLED_STEPS)),
           sorted(eb.STEP_ORDER)))

# Panel c: a FIXED set of pathway categories drawn as the same stacked
# segments on every bar (a partition of the Phi_M pools), with the
# translation sector phi_T as the dotted final-category tail. r2 is the
# PDH complex (TCA cycle); r4 (Ald6) -> r5 (Acs2) are the acetate bypass
# that regenerates cytosolic acetyl-CoA (its own category). The partition
# is asserted against eb.STEP_ORDER so a table drift raises at import.
BURDEN_CATEGORIES = (
    ('Glycolysis', ['r1']),
    ('TCA cycle', ['r2']),
    ('Acetate / acetyl-CoA\nproduction', ['r4', 'r5']),
    ('Ethanol production', ['r3', 'r6']),
    ('Isobutanol production', ['r13', 'r14', 'r15', 'r16']),
)
_CAT_STEPS = [st for _, steps in BURDEN_CATEGORIES for st in steps]
if sorted(_CAT_STEPS) != sorted(eb.STEP_ORDER) \
        or len(_CAT_STEPS) != len(set(_CAT_STEPS)):
    raise AssertionError(
        'burden category partition drift: BURDEN_CATEGORIES steps != '
        'eb.STEP_ORDER (%r vs %r)' % (sorted(_CAT_STEPS),
                                      sorted(eb.STEP_ORDER)))


def draw_burden(ax, sets, colors):
    # the four modeled categories share one tint ramp per bar colour
    # (darkest = glycolysis, lightest = isobutanol production); translation
    # is the dotted 5th-category tail
    cats = BURDEN_CATEGORIES
    tints = np.linspace(0.0, 0.66, len(cats))
    n = len(sets)
    F = sets[0]['F_flex']
    h = 0.6
    ypos = {id(s): n - i for i, s in enumerate(sets)}
    seg_centers = {}
    for s in sets:
        y = ypos[id(s)]; x = 0.0; c = colors[id(s)]; centers = []
        for i, (_, steps) in enumerate(cats):
            w = sum(s[f'pool_{st}'] for st in steps)
            ax.barh(y, w, left=x, height=h, color=tint(c, tints[i]),
                    edgecolor=c, lw=0.5, zorder=2)
            centers.append(x + w / 2); x += w
        # translation sector phi_T -- the 5th fixed category (dotted)
        ax.barh(y, s['phi_T'], left=x, height=h, color='none',
                edgecolor=c, lw=0.6, hatch='....', zorder=2)
        centers.append(x + s['phi_T'] / 2)
        end = x + s['phi_T']
        ax.text(max(end, F) + 0.005, y, f'd = {s["burden_factor"]:.2f}',
                va='center', fontsize=FONTS['tick'])
        seg_centers[id(s)] = centers
    # category callouts once, leaders to the FIRST campaign set (all study
    # pools non-zero there; the baseline's Ehrlich pools are zero)
    campaign_sets = [s for s in sets if not s.get('is_baseline')]
    anchor = campaign_sets[0] if campaign_sets else sets[0]
    y = ypos[id(anchor)]
    labels = [name for name, _ in cats] + ['Translation (ribosomes);\n'
                                           'growth derated by d']
    centers = seg_centers[id(anchor)]
    xs = np.linspace(0.015, 0.30, len(labels))
    for i, (lab, cx) in enumerate(zip(labels, centers)):
        row = n + 0.45 if i % 2 == 0 else n + 1.25
        ax.annotate(lab, xy=(cx, y + h / 2), xytext=(xs[i], row),
                    fontsize=FONTS['callout'], ha='center', va='bottom',
                    arrowprops=dict(arrowstyle='-', lw=0.5, color='0.35',
                                    shrinkA=0, shrinkB=0))
    ax.axvline(F, color='k', ls='--', lw=0.9, zorder=3)
    # F_flex rides its own dashed cap line, rotated 90 deg to run along it,
    # just below the lowest bar -- no separate empty row for it any more
    ax.text(F, 0.62, 'F$_{flex}$', rotation=90, ha='center', va='top',
            fontsize=FONTS['callout'])
    ax.set_ylim(0.05, n + 1.9)
    # rows are keyed by color through the legend, so the categorical y axis
    # carries no labels of its own
    ax.set_yticks([])
    # bounded on both ends by a labeled tick: 0 -> 0.35 in 0.05 steps
    # (0.35 clears the tallest bar + phi_T tail and the right-hand callouts)
    ax.set_xlim(0, 0.35)
    ax.set_xlabel('Enzyme burden Φ$_M$ [g enzyme·(g DCW)$^{-1}$]',
                  fontsize=FONTS['axis'])
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='y', right=False, length=0)
    ax.tick_params(axis='x', which='major', direction='inout', top=False,
                   length=4)
    ax.tick_params(axis='x', which='minor', direction='inout', top=False,
                   length=2.2)
    for sp in ('left', 'right', 'top'):
        ax.spines[sp].set_visible(False)


def plot(sets, band, out_stem, dpi=300):
    apply_fonts()
    fig = plt.figure(figsize=(9.5, 12.4))
    # Top-anchored vertical layout, figure fractions. Panel a and the wide
    # a -> b gap (which carries panel b's letter/title and the first band
    # title) are unchanged. Panel b's four band cells are ~40% shorter than
    # before (cell height ~0.0335 vs ~0.0558), so the band block ends higher;
    # the band -> band gaps stay ~0.05 (hspace 1.49 x the shorter cell) to keep
    # room for the rate-band titles. Panel c keeps its height and rides up just
    # below the bands, freeing space at the bottom of the canvas. Each region
    # is its own gridspec so the three vertical positions are set directly.
    LEFT, RIGHT = 0.083, 0.97
    a_gs = fig.add_gridspec(1, 1, left=LEFT, right=RIGHT, top=0.945, bottom=0.775)
    band_gs = fig.add_gridspec(4, 1, left=LEFT, right=RIGHT, top=0.689,
                               bottom=0.405, hspace=1.49)
    c_gs = fig.add_gridspec(1, 1, left=LEFT, right=RIGHT, top=0.359, bottom=0.134)
    colors = set_colors(sets)
    a_axes = draw_outcomes(fig, a_gs[0], sets, colors)
    b_axes = draw_parameters(fig, [band_gs[0], band_gs[1], band_gs[2],
                                   band_gs[3]], sets, colors, band)
    axc = fig.add_subplot(c_gs[0]); draw_burden(axc, sets, colors)
    # each panel gets a bold letter and a descriptive title on the same
    # baseline; the panel title (13 pt) outranks the band sub-titles (12 pt).
    # Panel b's letter is lifted into the panel-a -> b gap so its title clears
    # the first band title below it.
    panels = ((a_axes[0].get_position().y1 + 0.012, 'A',
               'Optimization incumbent trajectories'),
              (b_axes[0].get_position().y1 + 0.035, 'B',
               'Optimized kinetic and process parameters'),
              (axc.get_position().y1 + 0.005, 'C', 'Enzyme-burden allocation'))
    for y, letter, title in panels:
        fig.text(0.03, y, letter, fontsize=FONTS['panel'], fontweight='bold',
                 va='baseline')
        fig.text(0.055, y, title, fontsize=FONTS['panel'] - 1,
                 fontweight='bold', va='baseline')
    handles = [plt.Rectangle((0, 0), 1, 1, fc=colors[id(s)], label=s['label'])
               for s in sets]
    # the legend sits in the empty right columns of panel b's lower bands and
    # doubles as the row key for panel c (whose categorical y axis is
    # unlabelled). Raised above panel c -- which now rides higher after the
    # bands were shortened -- so the two no longer overlap.
    leg = fig.legend(handles=handles, loc='center left',
                     bbox_to_anchor=(0.635, 0.450), ncol=1, frameon=True,
                     fontsize=FONTS['legend'] + 1, title='Optimization campaign',
                     labelspacing=0.5, handlelength=1.7, handleheight=1.1,
                     borderpad=0.6, edgecolor='0.6', fancybox=False)
    leg.get_title().set_fontweight('bold')
    leg.get_title().set_fontsize(FONTS['legend'] + 2)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out_stem}.{ext}', dpi=dpi)
    plt.close(fig)
    return out_stem


def console_report(sets, band_campaign):
    print(f'band source (campaign): {band_campaign}')
    rank = sorted(_STUDY_STEPS,
                  key=lambda st: -max(s[f'pool_{st}'] for s in sets))[:5]
    for s in sets:
        tag = 'baseline' if s.get('is_baseline') \
            else (f'{s["campaign"]} trial {int(s["trial_number"])}; '
                  f'panel-a incumbent by {s["sel_col"]} ({s["sel_dir"]})')
        print(f'[{s["label"]}] {tag}')
        print(f'    IRR {s.get("IRR")!s:>8}  EtOH {s.get("EtOH titer")!s:>7}'
              f'  IBO {s.get("IBO titer")!s:>7}  tau {s.get("tau")!s:>6}')
        print(f'    Phi_M {s["Phi_M"]:.4f}  d {s["burden_factor"]:.2f}  pools: '
              + ', '.join(f'{st} {s[f"pool_{st}"]:.4f}' for st in rank))
        if s.get('extra_sampled'):
            print(f'    campaign {s["campaign"]}: also sampled '
                  f'{s["extra_sampled"]} (not shown)')


def build_sets(set_specs, include_baseline):
    sets = [baseline_set()] if include_baseline else []
    band = band_campaign = None
    for label, campaign, trial in set_specs:
        sets.append(load_set(label, campaign, trial))
        if band is None:
            band = campaign_band(campaign); band_campaign = campaign
        else:
            other = campaign_band(campaign)
            if other != band:
                print(f'  WARNING campaign {campaign}: searched band differs '
                      f'from the first campaign ({band_campaign}); keeping '
                      'the first campaign band for shading')
    if not sets:
        raise ValueError('no sets to plot (use --set or drop --no-baseline)')
    return sets, band, band_campaign


def main(argv=None):
    ap = argparse.ArgumentParser(description='Kinetic-optimization '
                                 'parameter-set comparison figure.')
    ap.add_argument('--set', dest='sets', action='append', nargs=3,
                    metavar=('LABEL', 'CAMPAIGN', 'TRIAL'), default=None,
                    help='a set to plot; repeatable, in bar order after the '
                         'baseline. TRIAL is an int, "best", or "best:COL".')
    ap.add_argument('--no-baseline', action='store_true',
                    help='drop the scenario-A baseline row')
    ap.add_argument('--out-dir', default=RESULTS_DIR)
    ap.add_argument('--stem', default=None)
    ap.add_argument('--dpi', type=int, default=300)
    args = ap.parse_args(argv)

    def norm_trial(t):
        if isinstance(t, str) and t.startswith('best'):
            return t
        return int(t)

    if args.sets:
        specs = [(lab, camp, norm_trial(tr)) for lab, camp, tr in args.sets]
    else:  # default: each optimum from the study that optimized it -- the
        # five most recent minimal-subset campaigns, one per objective
        specs = [('Financial attractiveness', DEFAULT_STUDY, 'best'),
                 ('Isobutanol titer', IBO_TITER_STUDY, 'best'),
                 ('Ethanol titer', ETOH_TITER_STUDY, 'best'),
                 ('Isobutanol yield', IBO_YIELD_STUDY, 'best'),
                 ('Ethanol yield', ETOH_YIELD_STUDY, 'best')]

    if len(specs) + (0 if args.no_baseline else 1) > MAX_SETS:
        raise ValueError(f'at most {MAX_SETS} sets (baseline + '
                         f'{len(HUE_COLORS)} campaign sets)')

    sets, band, band_campaign = build_sets(specs, not args.no_baseline)
    stem = args.stem or f'{os.path.basename(specs[0][1]).replace(".csv", "")}' \
                        '_parameter_sets'
    stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
    out_stem = os.path.join(args.out_dir, f'{stem}_{stamp}')
    plot(sets, band, out_stem, dpi=args.dpi)
    console_report(sets, band_campaign)
    print(f'wrote {out_stem}.png / .pdf')
    return out_stem


if __name__ == '__main__':
    main()
