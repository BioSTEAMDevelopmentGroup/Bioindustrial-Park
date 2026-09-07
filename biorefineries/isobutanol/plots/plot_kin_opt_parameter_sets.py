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

  a  outcomes -- IRR, ethanol titer, isobutanol titer, batch time.
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
import sys
import argparse
import importlib.util
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgb

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
# feeding cell: (title, (shaded-range low, high)) -- engine default bounds
FEED_LABELS = {'threshold_conc': ('Feed threshold\n(g/L)', (0, 300)),
               'target_delta': ('Target − threshold\n(g/L)', (5, 500)),
               'max_n_spikes': ('Max. glucose\nspikes', (0, 50))}

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

# outcomes panel: (CSV column, axis label, (y-low, y-high))
OUTCOMES = (('IRR', 'IRR', (0, 0.3)),
            ('EtOH titer', 'Ethanol titer (g/L)', (0, 300)),
            ('IBO titer', 'Isobutanol titer (g/L)', (0, 100)),
            ('tau', 'Batch time (h)', (0, 80)))

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
    'n_glu_spikes': 10, 'EtOH yield': 0.455, 'Cell density': 15.7,
    'TCI': 139.6, 'threshold_conc': 217.125, 'target_delta': 4.125,
    'max_n_spikes': 16,
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
    for col in ('IRR', 'EtOH titer', 'IBO titer', 'tau', 'n_glu_spikes'):
        rec[col] = BASELINE_A[col]
    for st, pool in res.pools.items():
        rec[f'pool_{st}'] = float(pool)
    rec['Phi_M'] = float(res.Phi_M)
    rec['phi_T'] = float(res.phi_T)
    rec['F_flex'] = float(res.F_flex)
    rec['burden_factor'] = float(res.burden_factor)
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


def resolve_trial(df, trial, campaign):
    """Return the requested COMPLETE row as a Series."""
    ok = df[df['state'] == 'COMPLETE']
    if isinstance(trial, str) and trial.startswith('best'):
        if ':' in trial:
            col = trial.split(':', 1)[1]
            if col not in df.columns:
                raise ValueError(
                    f'campaign {campaign}: no metric column {col!r} for '
                    f'"best:{col}"')
            return ok.loc[ok[col].idxmax()]
        obj = campaign_objective(campaign)
        if obj is None or obj not in ko.OBJECTIVE_REGISTRY:
            print(f'  WARNING campaign {campaign}: objective slug not in the '
                  'registry; "best" maximizes the objective column')
            direction = 'maximize'
        else:
            direction = ko.OBJECTIVE_REGISTRY[obj]['direction']
        idx = ok['objective'].idxmin() if direction == 'minimize' \
            else ok['objective'].idxmax()
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


def load_set(label, campaign, trial):
    """A campaign trial as a flat record (CSV row + metadata)."""
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


def sub(name):
    if '_' in name and name[0] in 'kK':
        head, tail = name.split('_', 1)
        return rf'$\mathit{{{head}}}_{{\mathrm{{{tail}}}}}$'
    return name


def fmt(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return 'n/a'
    return f'{v:.2g}' if v < 1000 else f'{v:.0f}'


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
    ax.tick_params(axis='y', which='both', direction='inout', right=False,
                   length=4)
    ax.tick_params(axis='x', which='both', top=False, bottom=False,
                   labelbottom=False)
    for sp in ('right', 'top'):
        ax.spines[sp].set_visible(False)


def _baseline_value(sets, var):
    for s in sets:
        if s.get('is_baseline'):
            return s.get(var)
    return None


def bar_cell(ax, sets, colors, var, kind, title, ylim=None, band=None):
    n = len(sets)
    base = _baseline_value(sets, var)
    if kind == 'rate':
        lo, hi = band[var]
        floor = lo / 6
        ax.set_yscale('log')
        ax.set_ylim(floor, hi * 2.5)
        ax.axhspan(lo, hi, color='0.92', zorder=0)
        ticks, labels = [lo, hi], [fmt(lo), fmt(hi)]
        if base and base > 0:
            ax.axhline(base, color='k', lw=0.8, ls='--', zorder=1)
            ticks.insert(1, base); labels.insert(1, fmt(base))
        ax.set_yticks(ticks); ax.set_yticklabels(labels)
        ax.set_yticks([], minor=True)
        bottom = floor
    elif kind == 'group':
        lo, hi = band[var]
        ax.set_yscale('log'); ax.set_ylim(lo * 0.7, hi * 1.4)
        ax.axhspan(lo, hi, color='0.92', zorder=0)
        ax.axhline(1.0, color='k', lw=0.8, ls='--', zorder=1)
        ax.set_yticks([lo, 1, hi]); ax.set_yticklabels([fmt(lo), '1', fmt(hi)])
        ax.set_yticks([], minor=True)
        bottom = lo * 0.7
    elif kind == 'feed':
        lo, hi = ylim
        ax.set_ylim(lo, hi * 1.06); ax.axhspan(lo, hi, color='0.92', zorder=0)
        if base is not None:
            ax.axhline(base, color='k', lw=0.8, ls='--', zorder=1)
        bottom = 0
    else:  # outcome
        lo, hi = ylim
        vmax = max([s.get(var, 0) or 0 for s in sets] + [hi])
        ax.set_ylim(lo, vmax * 1.1 if vmax > hi else hi)
        if base:
            ax.axhline(base, color='k', lw=0.8, ls='--', zorder=1)
        bottom = 0
    ax.set_title(title, fontsize=FONTS['cell'], pad=4)
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


def draw_outcomes(fig, gs_cell, sets, colors):
    sub_gs = gs_cell.subgridspec(1, len(OUTCOMES), wspace=0.55)
    for i, (col, label, yl) in enumerate(OUTCOMES):
        ax = fig.add_subplot(sub_gs[0, i])
        bar_cell(ax, sets, colors, col, 'outcome', label, ylim=yl)


def draw_parameters(fig, gs_rows, sets, colors, band):
    for gs_row, (title, params) in zip(gs_rows, BANDS):
        sub_gs = gs_row.subgridspec(1, 5, wspace=0.55)
        last_ax = None
        for i, p in enumerate(params):
            ax = fig.add_subplot(sub_gs[0, i]); last_ax = ax
            if p in RATE_VARS:
                bar_cell(ax, sets, colors, p, 'rate',
                         REACTION_LABELS[p], band=band)
            elif p in GROUP_VARS:
                bar_cell(ax, sets, colors, p, 'group',
                         GROUP_LABELS[p], band=band)
            else:
                t, rng = FEED_LABELS[p]
                bar_cell(ax, sets, colors, p, 'feed', t, ylim=rng)
        fig.text(0.19, last_ax.get_position().y1 + 0.043, title,
                 fontsize=FONTS['band'], fontweight='bold', va='bottom')
