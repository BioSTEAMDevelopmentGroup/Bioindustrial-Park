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
  c  enzyme burden -- one horizontal stacked bar per set: the five
     reaction steps with the largest pool (named by callouts), one
     hatched "other enzymes" lump (every remaining step, sampled or
     not: r2, r4, r5 and any study step outside the top five), and the
     translation sector phi_T as a dotted tail; the F_flex cap and the
     growth-derating factor d are marked.

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

With no --set arguments it plots the three best/best:EtOH titer/
best:IBO titer trials of the default minimal-subset IRR campaign against
the baseline. Writes <stem>_<stamp>.png and .pdf to --out-dir.
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
from matplotlib.ticker import (AutoMinorLocator, FuncFormatter, LogLocator,
                               MultipleLocator, NullFormatter, PercentFormatter)

# --- sim-safe module loads (by file path; never import the package) ----------
PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')
DEFAULT_STUDY = ('kin_opt_ethanol_isobutanol_metabolic_minimal_subset_irr'
                 '_rb0.001-10_ib0.2-2_burden')


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
GROUP_VARS = list(ko.METABOLIC_MINIMAL_SUBSET_GROUPS)        # 3
FEED_VARS = ['threshold_conc', 'target_delta', 'max_n_spikes']  # 3
DECISION_VARS = RATE_VARS + GROUP_VARS + FEED_VARS

BANDS = [
    ('Glycolysis / fermentation capacities   r1 → r3 → r6',
     ['k_1l', 'k_1h', 'k_1e', 'k_3', 'k_6']),
    ('Ehrlich-branch capacities   r13 → r14 → r15 → r16',
     ['k_13', 'k_14', 'k_15', 'k_16']),
    ('Product-inhibition effector multipliers (× baseline family)',
     GROUP_VARS),
    ('Feeding strategy', FEED_VARS),
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
# g/L/h for all nine. Effector multipliers are dimensionless fold-changes
# relative to the per-family baseline (1 = baseline).
RATE_UNIT = 'g/L/h'
GROUP_UNIT = '× baseline'
# feeding cell: (title, (shaded-range low, high)) -- engine default bounds
FEED_LABELS = {'threshold_conc': ('Feed threshold\n(g/L)', (0, 300)),
               'target_delta': ('Target − threshold\n(g/L)', (5, 500)),
               'max_n_spikes': ('Max. glucose\nspikes\n(count)', (0, 50))}

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

# outcomes panel: (CSV column, cell title, (y-low, y-high)), in draw order
OUTCOMES = (('IRR', 'Financial attractiveness,\nas IRR [%]', (0, 0.3)),
            ('IBO titer', 'Isobutanol titer\n(g/L)', (0, 100)),
            ('IBO yield', 'Isobutanol yield\n(g/g)', (0, 0.4)),
            ('EtOH titer', 'Ethanol titer\n(g/L)', (0, 300)),
            ('EtOH yield', 'Ethanol yield\n(g/g)', (0, 0.5)),
            ('tau', 'Batch time\n(h)', (0, 80)))

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

    rec = {'label': 'Scenario A baseline', 'campaign': None,
           'trial_number': None, 'is_baseline': True, 'extra_sampled': []}
    for k in RATE_VARS:
        rec[k] = float(k_ref[k])
    for g in GROUP_VARS:
        rec[g] = 1.0
    rec['threshold_conc'] = BASELINE_A['threshold_conc']
    rec['target_delta'] = BASELINE_A['target_delta']
    rec['max_n_spikes'] = BASELINE_A['max_n_spikes']
    for col in ('IRR', 'EtOH titer', 'IBO titer', 'EtOH yield', 'IBO yield',
                'tau', 'n_glu_spikes'):
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


def bar_cell(ax, sets, colors, var, kind, ylabel, subtitle=None, ylim=None,
             band=None):
    """One parameter cell: colored bars per set with the searched band
    shaded and the baseline dashed. `ylabel` names the value axis (the
    parameter symbol / group / feeding quantity); `subtitle`, if given, is
    a smaller descriptor above the cell (reaction + enzyme). Value axes get
    conventional ticks -- log decade majors + minor subdivisions on the
    rate/group cells, linear major + minor on the feeding cells."""
    n = len(sets)
    base = _baseline_value(sets, var)
    if kind == 'rate':
        lo, hi = band[var]
        floor = lo / 6
        ax.set_yscale('log')
        ax.set_ylim(floor, hi * 2.5)
        ax.axhspan(lo, hi, color='0.92', zorder=0)
        if base and base > 0:
            ax.axhline(base, color='k', lw=0.8, ls='--', zorder=1)
        bottom = floor
    elif kind == 'group':
        lo, hi = band[var]
        ax.set_yscale('log'); ax.set_ylim(lo * 0.7, hi * 1.4)
        ax.axhspan(lo, hi, color='0.92', zorder=0)
        ax.axhline(1.0, color='k', lw=0.8, ls='--', zorder=1)
        # narrow (<2 decade) log axis: plain-number majors at nice log stops,
        # the remaining log subdivisions as unlabeled minor ticks (default
        # log labeling would sci-notate every minor here)
        ax.yaxis.set_major_locator(
            LogLocator(base=10.0, subs=(0.2, 0.5, 1.0, 2.0), numticks=12))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
        ax.yaxis.set_minor_locator(LogLocator(
            base=10.0, subs=(0.3, 0.4, 0.6, 0.7, 0.8, 0.9), numticks=12))
        ax.yaxis.set_minor_formatter(NullFormatter())
        bottom = lo * 0.7
    else:  # feed
        lo, hi = ylim
        ax.set_ylim(lo, hi * 1.06); ax.axhspan(lo, hi, color='0.92', zorder=0)
        if base is not None:
            ax.axhline(base, color='k', lw=0.8, ls='--', zorder=1)
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_minor_locator(AutoMinorLocator())
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
            ax.text(j, y, label, ha='center', va='bottom',
                    fontsize=FONTS['tick'], color=c)
    ax.set_xlim(-0.6, n - 0.4); style_cell_axes(ax)


def outcome_cell(ax, sets, colors, col, title, ylim, xmax):
    """One outcome metric as incumbent trajectories over trial_number: one
    step line per campaign set (the metric at that set's running incumbent,
    in its panel-b/c color), the scenario-A baseline as a dashed reference.
    """
    lo, hi = ylim
    base = _baseline_value(sets, col)
    finite = [v for s in sets if s.get('traj') is not None
              for v in s['traj'][col] if np.isfinite(v)]
    if base is not None and np.isfinite(base):
        finite.append(base)
    vmax = max(finite + [hi]) if finite else hi
    ax.set_ylim(lo, vmax * 1.08 if vmax > hi else hi)
    ax.set_xlim(0, xmax * 1.02)
    if base is not None and np.isfinite(base):
        ax.axhline(base, color=BASELINE_COLOR, lw=0.9, ls='--', zorder=1)
    for s in sets:
        if s.get('traj') is None:
            continue
        ax.step(s['traj_x'], s['traj'][col], where='post',
                color=colors[id(s)], lw=1.4, zorder=2)
    # metric name next to the value axis itself (not a title above the cell)
    ax.set_ylabel(title, fontsize=FONTS['cell'], labelpad=3)
    ax.set_xlabel('Trial', fontsize=FONTS['tick'], labelpad=2)
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if col == 'IRR':   # fraction stored; show the value axis in percent
        ax.yaxis.set_major_locator(MultipleLocator(0.1))   # 10% steps
        ax.yaxis.set_major_formatter(
            PercentFormatter(xmax=1.0, decimals=0, symbol=''))
    else:
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
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
    sub_gs = gs_cell.subgridspec(1, len(OUTCOMES), wspace=0.85)
    axes = []
    for i, (col, label, yl) in enumerate(OUTCOMES):
        ax = fig.add_subplot(sub_gs[0, i])
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
                sym, _, sub = REACTION_LABELS[p].partition('\n')
                bar_cell(ax, sets, colors, p, 'rate', f'{sym}\n({RATE_UNIT})',
                         subtitle=sub, band=band)
            elif p in GROUP_VARS:
                bar_cell(ax, sets, colors, p, 'group',
                         f'{GROUP_LABELS[p]}\n({GROUP_UNIT})', band=band)
            else:
                t, rng = FEED_LABELS[p]
                bar_cell(ax, sets, colors, p, 'feed', t, ylim=rng)
        fig.text(0.19, last_ax.get_position().y1 + 0.043, title,
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


def draw_burden(ax, sets, colors):
    rank = sorted(_STUDY_STEPS,
                  key=lambda st: -max(s[f'pool_{st}'] for s in sets))
    top5, rest = rank[:5], rank[5:]
    other_steps = rest + _UNSAMPLED_STEPS
    tints = [0.0, 0.25, 0.45, 0.62, 0.78]   # darkest = largest pool
    n = len(sets)
    F = sets[0]['F_flex']
    h = 0.6
    ypos = {id(s): n - i for i, s in enumerate(sets)}
    seg_centers = {}
    for s in sets:
        y = ypos[id(s)]; x = 0.0; c = colors[id(s)]; centers = []
        segs = [(s[f'pool_{st}'], tint(c, tints[i]), None)
                for i, st in enumerate(top5)]
        segs.append((sum(s[f'pool_{st}'] for st in other_steps), c, '////'))
        for w, fc, hatch in segs:
            ax.barh(y, w, left=x, height=h,
                    color=('none' if hatch else fc),
                    edgecolor=c, lw=0.5, hatch=hatch, zorder=2)
            centers.append(x + w / 2); x += w
        ax.barh(y, s['phi_T'], left=x, height=h, color='none',
                edgecolor=c, lw=0.6, hatch='....', zorder=2)
        end = x + s['phi_T']
        ax.text(max(end, F) + 0.005, y, f'd = {s["burden_factor"]:.2f}',
                va='center', fontsize=FONTS['tick'])
        seg_centers[id(s)] = centers
    # callouts once, leaders to the FIRST campaign set (all seven study
    # pools non-zero there; the baseline's Ehrlich pools are zero)
    campaign_sets = [s for s in sets if not s.get('is_baseline')]
    anchor = campaign_sets[0] if campaign_sets else sets[0]
    y = ypos[id(anchor)]
    labels = [f'{STEP_ENZYME[st]}\n({STEP_PARAMS[st]})' for st in top5] + \
             ['other enzymes']
    xs = np.linspace(0.01, 0.235, len(labels))
    for i, (lab, cx) in enumerate(zip(labels, seg_centers[id(anchor)])):
        row = n + 1.0 if i % 2 == 0 else n + 2.1
        ax.annotate(lab, xy=(cx, y + h / 2), xytext=(xs[i], row),
                    fontsize=FONTS['callout'], ha='center', va='bottom',
                    arrowprops=dict(arrowstyle='-', lw=0.5, color='0.35',
                                    shrinkA=0, shrinkB=0))
    xt = anchor['Phi_M'] + anchor['phi_T'] * 0.8
    ax.annotate('translation sector φ$_T$ (ribosomes):\ngrowth derated '
                'by d where the\nbar crosses the F$_{flex}$ cap',
                xy=(xt, y + h / 2), xytext=(0.328, n + 2.1),
                fontsize=FONTS['callout'], ha='right', va='bottom',
                arrowprops=dict(arrowstyle='-', lw=0.5, color='0.35',
                                shrinkA=0, shrinkB=0))
    ax.axvline(F, color='k', ls='--', lw=0.9, zorder=3)
    ax.text(F + 0.004, -0.15, 'F$_{flex}$ = %.3f' % F, ha='left',
            va='center', fontsize=FONTS['callout'])
    ax.set_ylim(-0.7, n + 3.3)
    ax.set_yticks([ypos[id(s)] for s in sets])
    ax.set_yticklabels([s['label'] for s in sets], fontsize=FONTS['tick'])
    ax.set_xlim(0, 0.33)
    ax.set_xlabel('Enzyme burden Φ$_M$ (g enzyme / g DCW)',
                  fontsize=FONTS['axis'])
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
    fig = plt.figure(figsize=(9.5, 14.0))
    gs = fig.add_gridspec(6, 1,
                          height_ratios=[1.45, 1.15, 1.15, 1.15, 1.15, 2.2],
                          hspace=1.0, left=0.19, right=0.97, top=0.94,
                          bottom=0.05)
    colors = set_colors(sets)
    a_axes = draw_outcomes(fig, gs[0], sets, colors)
    b_axes = draw_parameters(fig, [gs[1], gs[2], gs[3], gs[4]], sets, colors,
                             band)
    axc = fig.add_subplot(gs[5]); draw_burden(axc, sets, colors)
    fig.text(0.03, a_axes[0].get_position().y1 + 0.012, 'a',
             fontsize=FONTS['panel'], fontweight='bold')
    fig.text(0.03, b_axes[0].get_position().y1 + 0.03, 'b',
             fontsize=FONTS['panel'], fontweight='bold')
    fig.text(0.03, axc.get_position().y1 + 0.005, 'c',
             fontsize=FONTS['panel'], fontweight='bold')
    handles = []
    for s in sets:
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=colors[id(s)],
                                     label=s['label']))
    fig.legend(handles=handles, loc='upper center',
               bbox_to_anchor=(0.56, 0.99), ncol=len(sets), frameon=False,
               fontsize=FONTS['legend'], columnspacing=1.2, handlelength=1.4)
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
    else:  # default: three trials of the default minimal-subset campaign
        specs = [('Financial attractiveness optimum', DEFAULT_STUDY, 'best'),
                 ('Ethanol titer optimum', DEFAULT_STUDY, 'best:EtOH titer'),
                 ('Isobutanol titer optimum', DEFAULT_STUDY, 'best:IBO titer')]

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
