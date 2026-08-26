#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Per-metric scatter grids against every uncertain parameter, from an
uncertainty-analysis results workbook (*_1_full_evaluation.xlsx as written by
analyses/full/uncertainties_IBO_EtOH.py).

One figure per metric: a grid of scatter panels (one per uncertain
parameter), each with density-colored points, a binned-median trend line, a
Spearman rho annotation (bold when p < 0.05), and a grey-diamond baseline
marker when the companion *_0_baseline.xlsx and the scenario's
parameter-distributions workbook are available.

Pure pandas/matplotlib/scipy -- does NOT import biorefineries, so it runs in
seconds with no simulation.

Usage:
    python plot_uncertainty_scatters.py <results.xlsx>
        [--metrics NAME [NAME ...]] [--baseline-file PATH]
        [--distributions-file PATH] [--output-dir DIR]
        [--dpi N] [--format EXT]
"""

import os
import re
import argparse
import textwrap
from math import ceil

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from scipy import stats

__all__ = ('plot_uncertainty_scatters', 'DEFAULT_METRICS', 'MetricSpec')


#%% Metric specifications

class MetricSpec:
    """Display label plus case-insensitive substring matching rules applied
    to the 'name [units]' level of the workbook's metric columns. `include`
    (then each of `fallbacks`, in order) must be a substring; none of
    `exclude` may be. `units` overrides the display units parsed from the
    matched column (used where the workbook's label is wrong, e.g. the
    ethanol MPSP column historically labeled '$/kg IBO')."""
    def __init__(self, label, include, exclude=(), fallbacks=(), units=None):
        self.label = label
        self.include = include
        self.exclude = exclude
        self.fallbacks = fallbacks
        self.units = units

DEFAULT_METRICS = [
    MetricSpec('IRR', 'irr'),
    MetricSpec('Ethanol MPSP', 'ethanol mpsp',
               fallbacks=('adjusted minimum selling price',), units='$/kg'),
    MetricSpec('Isobutanol MPSP', 'isobutanol mpsp', units='$/kg'),
    MetricSpec('Ethanol yield', 'et oh yield', exclude=('combined',)),
    MetricSpec('Ethanol titer', 'et oh titer'),
    MetricSpec('Ethanol productivity', 'et oh productivity'),
    MetricSpec('Isobutanol yield', 'ibo yield', exclude=('combined',)),
    MetricSpec('Isobutanol titer', 'ibo titer'),
    MetricSpec('Isobutanol productivity', 'ibo productivity'),
    MetricSpec('Cell density', 'cell loading', exclude=('active',)),
    MetricSpec('Active cell density', 'active cell loading'),
]


#%% Name/units helpers

def _strip_units(name):
    """'k 1e [g_per_l_per_h]' -> 'k 1e'."""
    return re.sub(r'\s*\[[^\]]*\]\s*$', '', str(name)).strip()

def _units_of(name):
    """'k 1e [g_per_l_per_h]' -> 'g_per_l_per_h'; '' if absent/'nan'."""
    m = re.search(r'\[([^\]]*)\]\s*$', str(name))
    units = m.group(1).strip() if m else ''
    return '' if units.lower() == 'nan' else units

def _pretty_units(units):
    """'g_per_l_per_h' -> 'g/L/h'; leaves e.g. '$/wet-kg' untouched."""
    u = units.replace('_per_', '/')
    return re.sub(r'/l\b', '/L', u)

def _normalize_name(name):
    """Match a results-workbook parameter name to a distributions-workbook
    one: biosteam turns '_' into ' ' and prefixes '.' oddly (distributions
    '._k_13' appears as '. k 13 [g_per_l_per_h]'), so strip the units
    bracket and remove '.', '_', and spaces. Case is PRESERVED: it is what
    distinguishes rate constant 'k 13' from saturation constant 'K 13'."""
    return re.sub(r'[.\s_]', '', _strip_units(name))

def _param_display(col):
    """Wrapped x-axis label for a parameter column (element, 'name [units]')."""
    name = _strip_units(col[1])
    name = re.sub(r'^\.\s*', '', name)  # '. k 1e' -> 'k 1e'
    label = '\n'.join(textwrap.wrap(name, 18))
    units = _pretty_units(_units_of(col[1]))
    if units:
        label += f'\n[{units}]'
    return label


#%% Colormap (constants copied from the kinetic-sweep scripts /
#   biosteam.utils.colors so this module needs no biosteam import)

def JBEI_UCB_colormap(N_levels=90, reverse=False):
    JBEI_orange = (233/255, 83/255, 39/255)
    UCB_blue = (0/255, 38/255, 118/255)
    UCB_yellow = (253/255, 181/255, 21/255)
    grey_dark = (0.282, 0.282, 0.278)  # biosteam colors.grey_dark.RGBn
    cmap_colors = [UCB_yellow, JBEI_orange, UCB_blue, grey_dark]
    if reverse:
        cmap_colors.reverse()
    return LinearSegmentedColormap.from_list('JBEI_UCB', cmap_colors, N_levels)


#%% Workbook loaders

def load_uncertainty_results(results_file):
    """Read the results workbook.

    Returns (raw, param_cols, metric_cols): `raw` is the 'Raw data' sheet
    (parameters + metrics on aligned rows, failed-simulation rows already
    dropped by the uncertainty script); `param_cols`/`metric_cols` are lists
    of its (element, 'name [units]') column tuples. Parameters are identified
    as the columns that also appear in the 'Parameters' sheet (ignoring its
    inserted 'Probability*' columns); the 'Blank parameter' placeholder is
    dropped from `param_cols`. If 'Raw data' is empty (a scenario whose
    all-NaN isobutanol MPSP column emptied the dropna'd table), fall back to
    joining the 'Parameters' values with 'TEA results' (a warning is printed).
    """
    pars_header = pd.read_excel(results_file, sheet_name='Parameters',
                                header=[0, 1], index_col=0, nrows=0)
    param_cols_sheet = [c for c in pars_header.columns
                        if not str(c[1]).startswith('Probability')]
    param_col_set = set(param_cols_sheet)

    raw = pd.read_excel(results_file, sheet_name='Raw data',
                        header=[0, 1], index_col=0)
    if raw.empty:
        print("WARNING: 'Raw data' sheet is empty; falling back to "
              "'Parameters' + 'TEA results'.")
        pars = pd.read_excel(results_file, sheet_name='Parameters',
                             header=[0, 1], index_col=0)
        pars = pars[param_cols_sheet]
        tea = pd.read_excel(results_file, sheet_name='TEA results',
                            header=[0, 1], index_col=0)
        raw = pars.join(tea)

    param_cols = [c for c in raw.columns if c in param_col_set
                  and 'blank parameter' not in str(c[1]).lower()]
    metric_cols = [c for c in raw.columns if c not in param_col_set]
    return raw, param_cols, metric_cols


def find_metric_column(spec, metric_cols):
    """First metric column whose lowercase 'name [units]' contains
    spec.include (else each fallback, in order) and none of spec.exclude.
    Returns the column tuple or None; warns when a needle is ambiguous."""
    for needle in (spec.include, *spec.fallbacks):
        matches = [
            c for c in metric_cols
            if needle in str(c[1]).lower()
            and not any(x in str(c[1]).lower() for x in spec.exclude)
        ]
        if matches:
            if len(matches) > 1:
                print(f"WARNING: '{spec.label}' matcher '{needle}' matched "
                      f"{len(matches)} columns; using {matches[0]}.")
            return matches[0]
    return None


def load_param_baselines(distributions_file, param_cols):
    """{results param column tuple: baseline value} from the distributions
    workbook's 'Baseline' column, matched on (element lowercased,
    case-preserving normalized name)."""
    df = pd.read_excel(distributions_file)
    lookup = {}
    for _, row in df.iterrows():
        key = (str(row['Element']).lower(), _normalize_name(row['Parameter name']))
        lookup[key] = row['Baseline']
    out = {}
    for c in param_cols:
        key = (str(c[0]).lower(), _normalize_name(c[1]))
        if key in lookup and pd.notna(lookup[key]):
            out[c] = float(lookup[key])
    return out


def load_metric_baselines(baseline_file):
    """The 'initial' row of the *_0_baseline.xlsx workbook, as a Series
    indexed by (element, 'name [units]') column tuples."""
    df = pd.read_excel(baseline_file, header=[0, 1], index_col=0)
    return df.loc['initial']


def detect_scenario(results_file):
    """'A' or 'B' from the results filename; None if undetectable."""
    name = os.path.basename(results_file)
    m = re.search(r"\['([AB])'\]", name)
    if m:
        return m.group(1)
    m = re.search(r'_([AB])_[01]_', name)
    if m:
        return m.group(1)
    return None


#%% Panel ingredients

def _point_densities(x, y, max_fit=1000, seed=0):
    """Gaussian-KDE local density per point, fit on z-scored coordinates
    (subsampled to <= max_fit points for speed) and evaluated on all points.
    Zeros when either axis is constant or the KDE is singular."""
    sx, sy = x.std(), y.std()
    if not (np.isfinite(sx) and sx > 0 and np.isfinite(sy) and sy > 0):
        return np.zeros(len(x))
    z = np.column_stack([(x - x.mean())/sx, (y - y.mean())/sy])
    fit = z
    if len(z) > max_fit:
        rng = np.random.default_rng(seed)
        fit = z[rng.choice(len(z), max_fit, replace=False)]
    try:
        return stats.gaussian_kde(fit.T)(z.T)
    except Exception:
        return np.zeros(len(x))


def _binned_median(x, y, n_bins=15, min_per_bin=5):
    """(bin-median x, bin-median y) arrays over ~equal-count x-bins, or None
    when there is not enough spread (< 4 distinct bin edges or < 3 usable
    bins)."""
    edges = np.unique(np.quantile(x, np.linspace(0, 1, n_bins + 1)))
    if len(edges) < 4:
        return None
    idx = np.clip(np.searchsorted(edges, x, side='right') - 1, 0, len(edges) - 2)
    xs, ys = [], []
    for b in range(len(edges) - 1):
        mask = idx == b
        if mask.sum() >= min_per_bin:
            xs.append(np.median(x[mask]))
            ys.append(np.median(y[mask]))
    if len(xs) < 3:
        return None
    return np.array(xs), np.array(ys)


#%% Figure

def plot_one_metric(label, units, y, params_df, param_baselines,
                    metric_baseline, out_path, dpi=600, ncols=9):
    """One scatter-grid figure: `y` (Series, aligned with params_df rows)
    against every column of params_df. `param_baselines` maps param column ->
    baseline x; `metric_baseline` is the baseline y (NaN to skip markers)."""
    n = params_df.shape[1]
    nrows = ceil(n/ncols)
    fig, axs = plt.subplots(nrows, ncols, squeeze=False, sharey=True,
                            figsize=(2.0*ncols, 1.7*nrows),
                            constrained_layout=True)
    y_all = y.to_numpy(float)
    y_finite = y_all[np.isfinite(y_all)]
    if y_finite.size:
        lo, hi = np.percentile(y_finite, [1, 99])
        if hi <= lo:
            lo, hi = lo - 1.0, hi + 1.0
        pad = 0.05*(hi - lo)
        axs[0][0].set_ylim(lo - pad, hi + pad)
    cmap = JBEI_UCB_colormap(reverse=True)
    n_baseline_markers = 0
    for i, col in enumerate(params_df.columns):
        ax = axs[i//ncols][i % ncols]
        x = params_df[col].to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(y_all)
        xx, yy = x[mask], y_all[mask]
        if xx.size:
            dens = _point_densities(xx, yy)
            order = np.argsort(dens)  # high-density points drawn last
            ax.scatter(xx[order], yy[order], c=dens[order], cmap=cmap,
                       s=3, linewidths=0, rasterized=True)
            trend = _binned_median(xx, yy)
            if trend is not None:
                ax.plot(*trend, color='0.15', lw=1.2, zorder=4)
            rho, p = stats.spearmanr(xx, yy)
            if np.isfinite(rho):
                significant = np.isfinite(p) and p < 0.05
                ax.text(0.04, 0.95, f'ρ = {rho:.2f}',
                        transform=ax.transAxes, fontsize=7, va='top',
                        fontweight='bold' if significant else 'normal',
                        color='k' if significant else '0.45', zorder=6)
        if col in param_baselines and np.isfinite(metric_baseline):
            ax.plot(param_baselines[col], metric_baseline, ls='none',
                    marker='D', ms=5, mfc='0.6', mec='k', mew=0.7, zorder=5)
            n_baseline_markers += 1
        ax.set_xlabel(_param_display(col), fontsize=6.5)
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(which='both', direction='in', labelsize=6)
    for j in range(n, nrows*ncols):
        axs[j//ncols][j % ncols].set_axis_off()
    title = f'{label} [{units}]' if units else label
    fig.suptitle(title, fontweight='bold', fontsize=13)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return n_baseline_markers


#%% Orchestration

def plot_uncertainty_scatters(results_file,
                              metrics=None,
                              baseline_file=None,
                              distributions_file=None,
                              output_dir=None,
                              dpi=600,
                              file_format='png'):
    """Generate one scatter-grid figure per metric from an uncertainty
    results workbook. `metrics` is None (-> DEFAULT_METRICS), or a list of
    MetricSpec and/or name strings (a string is matched as a
    case-insensitive substring). Returns the list of saved figure paths;
    requested metrics absent from the workbook are warned and skipped."""
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

    results_file = os.path.abspath(results_file)
    if output_dir is None:
        output_dir = os.path.dirname(results_file)

    raw, param_cols, metric_cols = load_uncertainty_results(results_file)

    if metrics is None:
        specs = DEFAULT_METRICS
    else:
        specs = [s if isinstance(s, MetricSpec) else MetricSpec(s, s.lower())
                 for s in metrics]
    resolved, missing = [], []
    for spec in specs:
        col = find_metric_column(spec, metric_cols)
        (resolved.append((spec, col)) if col is not None
         else missing.append(spec.label))
    if missing:
        print(f'WARNING: metrics not found in this workbook (skipped): '
              f'{missing}')

    # Baseline metric values (optional)
    if baseline_file is None:
        candidate = results_file.replace('_1_full_evaluation', '_0_baseline')
        baseline_file = (candidate if candidate != results_file
                         and os.path.isfile(candidate) else None)
    metric_baselines = None
    if baseline_file is not None:
        try:
            metric_baselines = load_metric_baselines(baseline_file)
        except Exception as e:
            print(f'NOTE: could not read baseline workbook '
                  f'({baseline_file}): {e}; baseline markers skipped.')
    else:
        print('NOTE: no companion *_0_baseline.xlsx found; '
              'baseline markers skipped.')

    # Baseline parameter values (optional)
    if distributions_file is None:
        scenario = detect_scenario(results_file)
        if scenario is not None:
            candidate = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'full', 'parameter_distributions',
                f'parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx')
            distributions_file = (candidate if os.path.isfile(candidate)
                                  else None)
    param_baselines = {}
    if distributions_file is not None:
        try:
            param_baselines = load_param_baselines(distributions_file,
                                                   param_cols)
        except Exception as e:
            print(f'NOTE: could not read distributions workbook '
                  f'({distributions_file}): {e}; baseline markers skipped.')
    else:
        print('NOTE: no parameter-distributions workbook found; '
              'baseline markers skipped.')

    os.makedirs(output_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(results_file))[0]
    saved = []
    for spec, col in resolved:
        units = spec.units or _pretty_units(_units_of(col[1]))
        baseline_val = np.nan
        if metric_baselines is not None:
            try:
                baseline_val = float(metric_baselines[col])
            except (KeyError, TypeError, ValueError):
                pass
        out_path = os.path.join(
            output_dir,
            f"{stem}_scatter_{spec.label.replace(' ', '-')}.{file_format}")
        n_markers = plot_one_metric(spec.label, units, raw[col],
                                    raw[param_cols], param_baselines,
                                    baseline_val, out_path, dpi=dpi)
        print(f'Saved {out_path} ({n_markers}/{len(param_cols)} '
              f'baseline markers).')
        saved.append(out_path)
    return saved


#%% CLI

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Per-metric scatter grids vs. every uncertain parameter '
                    'from an uncertainty-analysis results workbook.')
    parser.add_argument('results_file',
                        help='path to a *_1_full_evaluation.xlsx workbook')
    parser.add_argument('--metrics', nargs='+', default=None,
                        help='metric names (case-insensitive substrings); '
                             'default: the standard set')
    parser.add_argument('--baseline-file', default=None)
    parser.add_argument('--distributions-file', default=None)
    parser.add_argument('--output-dir', default=None)
    parser.add_argument('--dpi', type=int, default=600)
    parser.add_argument('--format', dest='file_format', default='png')
    args = parser.parse_args(argv)
    plot_uncertainty_scatters(args.results_file,
                              metrics=args.metrics,
                              baseline_file=args.baseline_file,
                              distributions_file=args.distributions_file,
                              output_dir=args.output_dir,
                              dpi=args.dpi,
                              file_format=args.file_format)

if __name__ == '__main__':
    main()
