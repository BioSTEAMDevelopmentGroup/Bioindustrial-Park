#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Plots for the scenario-B one-at-a-time kinetic-parameter sweeps produced by
evaluate_OAT_kinetics_B.py. Reads only the CSVs in
analyses/results/OAT_kinetics_B/ -- no simulation, safe to rerun anytime.

Outputs (same directory):
* OAT_B_<param>.png -- per-parameter figure, 3x2 panels: TEA metrics (IRR,
  ethanol MPSP, isobutanol MPSP) in the left column; fermentation metrics
  (product/cell titers, yields, productivities) in the right column.
* OAT_B_overview_<metric>.png -- one grid of mini-panels per metric, all
  parameters on a shared normalized x axis (value / baseline, 0-10).
"""

import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'results', 'OAT_kinetics_B')

baseline_df = pd.read_csv(os.path.join(results_dir,
                                       'kinetic_params_baseline_B.csv'))
baseline = pd.read_csv(os.path.join(results_dir,
                                    'baseline_TEA_B.csv')).iloc[0]

UNIT_LABELS = {'g_per_l_per_h': r'g L$^{-1}$ h$^{-1}$',
               'g_per_l': r'g L$^{-1}$',
               'l_per_g': r'L g$^{-1}$',
               'dimensionless': 'dimensionless'}

# column -> (label, units, color); colors: UCB blue / JBEI orange / green
# for TEA, and product-consistent hues for fermentation metrics
METRICS = {
    'IRR':               ('IRR', '', '#002676'),
    'MPSP_ethanol':      ('Ethanol MPSP', r'\$ kg$^{-1}$', '#E95327'),
    'MPSP_isobutanol':   ('Isobutanol MPSP', r'\$ kg$^{-1}$', '#3B7D23'),
    'EtOH_titer':        ('EtOH titer', r'g L$^{-1}$', '#E95327'),
    'IBO_titer':         ('IBO titer', r'g L$^{-1}$', '#3B7D23'),
    'cell_titer':        ('Cell mass titer', r'g L$^{-1}$', '#6B4C9A'),
    'EtOH_yield':        ('EtOH yield', r'g g-glucose$^{-1}$', '#E95327'),
    'IBO_yield':         ('IBO yield', r'g g-glucose$^{-1}$', '#3B7D23'),
    'EtOH_productivity': ('EtOH productivity', r'g L$^{-1}$ h$^{-1}$',
                          '#E95327'),
    'IBO_productivity':  ('IBO productivity', r'g L$^{-1}$ h$^{-1}$',
                          '#3B7D23'),
}

TEA_COLS = ['IRR', 'MPSP_ethanol', 'MPSP_isobutanol']
# right-hand column of the per-parameter figure: grouped fermentation panels
FERM_PANELS = [
    ('Titers (g L$^{-1}$)', ['EtOH_titer', 'IBO_titer', 'cell_titer']),
    ('Yields (g g-glucose$^{-1}$)', ['EtOH_yield', 'IBO_yield']),
    ('Productivities (g L$^{-1}$ h$^{-1}$)',
     ['EtOH_productivity', 'IBO_productivity']),
]


def case_safe_stem(p):
    """Windows filenames are case-insensitive (k_1l vs K_1l would collide),
    so uppercase (affinity K_*) parameters get a '_cap' suffix -- must match
    evaluate_OAT_kinetics_B.sweep_csv_name."""
    return f'OAT_B_{p}_cap' if p[0] == 'K' else f'OAT_B_{p}'


sweeps = {}  # param -> DataFrame
for _, prow in baseline_df.iterrows():
    p = prow['parameter']
    csv_path = os.path.join(results_dir, case_safe_stem(p) + '.csv')
    if not os.path.exists(csv_path):
        print(f'No results yet for {p}; skipping.')
        continue
    sweeps[p] = pd.read_csv(csv_path)

baseline_by_param = baseline_df.set_index('parameter')

#%% Per-parameter figures
for p, df in sweeps.items():
    b = baseline_by_param.loc[p, 'baseline_B']
    unit = UNIT_LABELS.get(baseline_by_param.loc[p, 'units'],
                           str(baseline_by_param.loc[p, 'units']))
    x = df['value'].values
    fig, axs = plt.subplots(3, 2, figsize=(10.5, 7.5), sharex=True)

    # Left column: TEA metrics
    for ax, col in zip(axs[:, 0], TEA_COLS):
        label, yunits, color = METRICS[col]
        y = df[col].values
        ax.axvline(b, color='0.55', ls='--', lw=1,
                   label=f'baseline = {b:.4g}')
        base_y = baseline.get(col, np.nan)
        if np.isfinite(base_y):
            ax.axhline(base_y, color='0.75', ls=':', lw=1)
            ax.plot([b], [base_y], 'o', mfc='none', mec='0.3', ms=7, zorder=5)
        ax.plot(x, y, '-o', color=color, ms=4, lw=1.3)
        failed = ~np.isfinite(y)
        if failed.any():
            ax.plot(x[failed], np.full(failed.sum(), ax.get_ylim()[0]),
                    'x', color='crimson', ms=5, clip_on=False,
                    label=f'{failed.sum()} NaN/failed')
        ax.set_ylabel(label + (f' ({yunits})' if yunits else ''), fontsize=10)
        ax.tick_params(labelsize=9)
        ax.legend(fontsize=8, loc='best', frameon=False)

    # Right column: fermentation metrics, grouped
    for ax, (panel_label, cols) in zip(axs[:, 1], FERM_PANELS):
        ax.axvline(b, color='0.55', ls='--', lw=1)
        for col in cols:
            label, _, color = METRICS[col]
            if col not in df.columns: continue
            y = df[col].values
            ax.plot(x, y, '-o', color=color, ms=3.5, lw=1.2, label=label)
            base_y = baseline.get(col, np.nan)
            if np.isfinite(base_y):
                ax.plot([b], [base_y], 'o', mfc='none', mec=color, ms=7,
                        zorder=5)
        ax.set_ylabel(panel_label, fontsize=10)
        ax.tick_params(labelsize=9)
        ax.legend(fontsize=8, loc='best', frameon=False)

    for ax in axs[-1, :]:
        ax.set_xlabel(f'{p} ({unit})', fontsize=11)
    fig.suptitle(f'OAT sweep of {p} (scenario B)', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(results_dir, case_safe_stem(p) + '.png'), dpi=300)
    plt.close(fig)
print(f'Saved {len(sweeps)} per-parameter figures.')

#%% Overview grids (one per metric); shared normalized x = value/baseline
params_sorted = [p for p in baseline_df['parameter'] if p in sweeps]
n = len(params_sorted)
ncols = 7
nrows = int(np.ceil(n/ncols))

for col, (label, yunits, color) in METRICS.items():
    if not any(col in df.columns for df in sweeps.values()):
        print(f'No data recorded for {col}; skipping overview grid.')
        continue
    fig, axs = plt.subplots(nrows, ncols,
                            figsize=(2.05*ncols, 1.75*nrows),
                            sharex=True)
    axs = np.atleast_2d(axs)
    base_y = baseline.get(col, np.nan)
    for i, p in enumerate(params_sorted):
        ax = axs[i//ncols, i % ncols]
        df = sweeps[p]
        b = baseline_by_param.loc[p, 'baseline_B']
        xn = df['value'].values / b if b else df['value'].values
        y = df[col].values
        if np.isfinite(base_y):
            ax.axhline(base_y, color='0.8', ls=':', lw=0.8)
        ax.axvline(1.0, color='0.8', ls='--', lw=0.8)
        ax.plot(xn, y, '-', color=color, lw=1.1)
        ax.plot(xn, y, '.', color=color, ms=2.5)
        ax.set_title(p, fontsize=8, pad=2)
        ax.tick_params(labelsize=6.5)
    for j in range(n, nrows*ncols):
        axs[j//ncols, j % ncols].set_visible(False)
    sup_units = f' ({yunits})' if yunits else ''
    fig.suptitle(f'{label}{sup_units} vs kinetic parameter '
                 '(x = value / scenario-B baseline, 0–10)',
                 fontsize=12)
    fig.supxlabel('parameter value / baseline', fontsize=10)
    fig.tight_layout(rect=(0, 0.01, 1, 0.97))
    fig.savefig(os.path.join(results_dir, f'OAT_B_overview_{col}.png'),
                dpi=300)
    plt.close(fig)
    print(f'Saved overview grid for {col}.')
