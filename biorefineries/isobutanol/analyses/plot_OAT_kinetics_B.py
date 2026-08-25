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
* OAT_B_<param>.png -- per-parameter figure, 3 stacked panels (IRR,
  ethanol MPSP, isobutanol MPSP) vs the parameter value, baseline marked.
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
baseline_TEA = pd.read_csv(os.path.join(results_dir,
                                        'baseline_TEA_B.csv')).iloc[0]

UNIT_LABELS = {'g_per_l_per_h': r'g L$^{-1}$ h$^{-1}$',
               'g_per_l': r'g L$^{-1}$',
               'l_per_g': r'L g$^{-1}$',
               'dimensionless': 'dimensionless'}

METRICS = [('IRR', 'IRR', ''),
           ('MPSP_ethanol', 'Ethanol MPSP', r'\$ kg$^{-1}$'),
           ('MPSP_isobutanol', 'Isobutanol MPSP', r'\$ kg$^{-1}$')]

COLORS = {'IRR': '#002676',            # UCB blue
          'MPSP_ethanol': '#E95327',   # JBEI orange
          'MPSP_isobutanol': '#3B7D23'}

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
    fig, axs = plt.subplots(3, 1, figsize=(5.5, 7.5), sharex=True)
    for ax, (col, label, yunits) in zip(axs, METRICS):
        y = df[col].values
        x = df['value'].values
        ax.axvline(b, color='0.55', ls='--', lw=1,
                   label=f'baseline = {b:.4g}')
        base_y = baseline_TEA[col]
        if np.isfinite(base_y):
            ax.axhline(base_y, color='0.75', ls=':', lw=1)
            ax.plot([b], [base_y], 'o', mfc='none', mec='0.3', ms=7, zorder=5)
        ax.plot(x, y, '-o', color=COLORS[col], ms=4, lw=1.3)
        failed = ~np.isfinite(y)
        if failed.any():
            ax.plot(x[failed], np.full(failed.sum(), ax.get_ylim()[0]),
                    'x', color='crimson', ms=5, clip_on=False,
                    label=f'{failed.sum()} NaN/failed')
        ylabel = f'{label}' + (f' ({yunits})' if yunits else '')
        ax.set_ylabel(ylabel, fontsize=10)
        ax.tick_params(labelsize=9)
        ax.legend(fontsize=8, loc='best', frameon=False)
    axs[-1].set_xlabel(f'{p} ({unit})', fontsize=11)
    axs[0].set_title(f'OAT sweep of {p} (scenario B)', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, case_safe_stem(p) + '.png'), dpi=300)
    plt.close(fig)
print(f'Saved {len(sweeps)} per-parameter figures.')

#%% Overview grids (one per metric); shared normalized x = value/baseline
params_sorted = [p for p in baseline_df['parameter'] if p in sweeps]
n = len(params_sorted)
ncols = 7
nrows = int(np.ceil(n/ncols))

for col, label, yunits in METRICS:
    fig, axs = plt.subplots(nrows, ncols,
                            figsize=(2.05*ncols, 1.75*nrows),
                            sharex=True)
    axs = np.atleast_2d(axs)
    base_y = baseline_TEA[col]
    for i, p in enumerate(params_sorted):
        ax = axs[i//ncols, i % ncols]
        df = sweeps[p]
        b = baseline_by_param.loc[p, 'baseline_B']
        xn = df['value'].values / b if b else df['value'].values
        y = df[col].values
        if np.isfinite(base_y):
            ax.axhline(base_y, color='0.8', ls=':', lw=0.8)
        ax.axvline(1.0, color='0.8', ls='--', lw=0.8)
        ax.plot(xn, y, '-', color=COLORS[col], lw=1.1)
        ax.plot(xn, y, '.', color=COLORS[col], ms=2.5)
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
