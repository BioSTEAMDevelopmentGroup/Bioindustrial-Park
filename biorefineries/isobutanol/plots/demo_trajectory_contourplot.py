#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Headless demo for plot_metric_with_trajectories.

Loads the plotting module standalone (no biorefinery import / no simulation),
finds a real k_13 x k_7ii sweep CSV set under analyses/results/ if present
(falling back to a synthetic grid otherwise), draws the MPSP contour with
greedy trajectories for MPSP / IRR / IBO MPSP, and saves a PNG next to this
file. Run:  & "$py" biorefineries/isobutanol/plots/demo_trajectory_contourplot.py
"""
import os
import glob
import importlib.util

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, '..', 'analyses', 'results')


def _load_module():
    path = os.path.join(HERE, 'trajectory_contourplot.py')
    spec = importlib.util.spec_from_file_location('trajectory_contourplot', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_metric_csv(prefix, metric):
    # sweep CSVs are saved as '<prefix>_<metric>.csv' with a leading index col;
    # the 2-D block is the remaining columns.
    matches = glob.glob(os.path.join(RESULTS_DIR, f'*{metric}.csv'))
    matches = [m for m in matches if metric in os.path.basename(m)]
    if not matches:
        return None
    df = pd.read_csv(matches[0])
    arr = df[df.columns[1:]].to_numpy(dtype='float', na_value=np.nan)
    return arr[None, :, :]


def _real_results():
    metrics = ['MPSP', 'IRR', 'IBO MPSP']
    results = {}
    for m in metrics:
        arr = _load_metric_csv('*k_13*k_7ii*', m)
        if arr is None:
            return None
        results[m] = arr
    ny, nx = results['MPSP'].shape[1:]
    spec_1 = np.linspace(0.0, 40.0, nx)     # k_13
    spec_2 = np.linspace(0.0001, 0.5, ny)   # k_7ii
    return results, spec_1, spec_2


def _synthetic_results():
    spec_1 = np.linspace(0.0, 40.0, 20)
    spec_2 = np.linspace(0.0001, 0.5, 20)
    iy, ix = np.mgrid[0:20, 0:20]
    results = {
        'MPSP': (0.4 + 0.02 * ((iy - 0) ** 2 + (ix - 0) ** 2) ** 0.5)[None, :, :],
        'IRR': (0.05 + 0.01 * (iy + ix))[None, :, :],
        'IBO MPSP': (0.9 + 0.03 * ((iy - 19) ** 2 + (ix - 0) ** 2) ** 0.5)[None, :, :],
    }
    return results, spec_1, spec_2


def main():
    tc = _load_module()
    real = _real_results()
    if real is None:
        print("No real k_13 x k_7ii sweep CSVs found; using synthetic data.")
        results, s1, s2 = _synthetic_results()
    else:
        print("Loaded real sweep CSVs.")
        results, s1, s2 = real
    # k_13 > 0 so IBO MPSP is defined (the k_13 = 0 column is NaN: no isobutanol at zero formation rate)
    fig, ax, traj = tc.plot_metric_with_trajectories(
        results, s1, s2,
        color_metric='MPSP',
        color_metric_units=r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$",
        baseline_point=(10.0, 0.04),
        trajectory_metrics=['MPSP', 'IRR', 'IBO MPSP'],
        senses={'MPSP': 'min', 'IRR': 'max', 'IBO MPSP': 'min'},
        x_label=r"$\mathbf{k}_{13}$", y_label=r"$\mathbf{k}_{7ii}$",
    )
    out = os.path.join(HERE, 'demo_trajectory_contourplot.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print("Saved:", out)
    for m, d in traj.items():
        print(f"{m}: n_steps={d['n_steps']}, optimum_xy={d['optimum_xy']}, "
              f"optimum_value={d['optimum_value']:.4g}")


if __name__ == '__main__':
    main()
