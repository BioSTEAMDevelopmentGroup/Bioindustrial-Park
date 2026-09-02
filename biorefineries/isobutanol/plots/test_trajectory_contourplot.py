#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Standalone tests for plots/trajectory_contourplot.py.

Run:  & "$py" biorefineries/isobutanol/plots/test_trajectory_contourplot.py
Loads the module directly via importlib so it never imports the biorefinery
package (no import-time simulation, no numba cache). Uses the Agg backend.
"""
import os
import importlib.util

import matplotlib
matplotlib.use('Agg')
import numpy as np


def _load_module():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'trajectory_contourplot.py')
    spec = importlib.util.spec_from_file_location('trajectory_contourplot', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


tc = _load_module()


def test_squeeze_grid_accepts_2d_and_3d():
    a2 = np.arange(6.0).reshape(2, 3)
    assert tc._squeeze_grid(a2).shape == (2, 3)
    a3 = np.arange(6.0).reshape(1, 2, 3)
    assert tc._squeeze_grid(a3).shape == (2, 3)
    assert np.array_equal(tc._squeeze_grid(a3), a2)


def test_snap_index_picks_nearest():
    spec = np.array([0.0, 10.0, 20.0, 30.0, 40.0])
    assert tc._snap_index(spec, 0.0) == 0
    assert tc._snap_index(spec, 12.0) == 1
    assert tc._snap_index(spec, 26.0) == 3


def test_greedy_climb_walks_diagonally_to_min():
    # bowl with min at (4,4); 8-neighbor climb from (0,0) goes diagonally.
    iy, ix = np.mgrid[0:5, 0:5]
    grid = (iy - 4.0) ** 2 + (ix - 4.0) ** 2
    path = tc._greedy_climb(grid, (0, 0), 'min')
    assert path[0] == (0, 0)
    assert path[-1] == (4, 4)
    assert path == [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]
    vals = [grid[c] for c in path]
    assert all(b < a for a, b in zip(vals, vals[1:]))  # strictly improving


def test_greedy_climb_max_sense():
    iy, ix = np.mgrid[0:5, 0:5]
    grid = -((iy - 4.0) ** 2 + (ix - 4.0) ** 2)  # max at (4,4)
    path = tc._greedy_climb(grid, (0, 0), 'max')
    assert path[-1] == (4, 4)


def test_greedy_climb_skips_nan_neighbors():
    iy, ix = np.mgrid[0:5, 0:5]
    grid = (iy - 4.0) ** 2 + (ix - 4.0) ** 2
    grid[1, 1] = np.nan  # block the diagonal; climb must route around it
    path = tc._greedy_climb(grid, (0, 0), 'min')
    assert (1, 1) not in path
    assert path[-1] == (4, 4)


def test_greedy_climb_start_at_optimum_is_single_point():
    iy, ix = np.mgrid[0:5, 0:5]
    grid = (iy - 4.0) ** 2 + (ix - 4.0) ** 2
    path = tc._greedy_climb(grid, (4, 4), 'min')
    assert path == [(4, 4)]


def test_greedy_climb_nan_start_raises():
    grid = np.zeros((3, 3))
    grid[0, 0] = np.nan
    try:
        tc._greedy_climb(grid, (0, 0), 'min')
    except ValueError as e:
        assert 'NaN' in str(e)
    else:
        raise AssertionError("expected ValueError for NaN start cell")


def test_greedy_climb_invalid_sense_raises():
    try:
        tc._greedy_climb(np.zeros((2, 2)), (0, 0), 'minimize')
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for invalid sense")


def test_greedy_climb_min_rel_improvement_ignores_round_off_plateau():
    # a flat column (differences ~1e-8 relative, like a converged simulation's
    # noise) next to a genuinely better cell two columns away that the climb
    # cannot reach without first taking a noise-level step.
    grid = np.full((3, 3), 100.0)
    grid[:, 0] = [100.0, 100.0 + 1e-6, 100.0 + 2e-6]   # plateau, rising "up"
    grid[2, 1] = 100.0 + 3e-6
    grid[2, 2] = 150.0
    # default (any strict improvement): walks the plateau and reaches 150
    path = tc._greedy_climb(grid, (0, 0), 'max')
    assert path[-1] == (2, 2)
    # with a threshold above the noise the start is already a local optimum
    # ((1, 1) is exactly equal, the plateau cells are within tolerance)
    assert tc._greedy_climb(grid, (0, 0), 'max', min_rel_improvement=1e-4) == [(0, 0)]
    # a real improvement still passes the threshold
    assert tc._greedy_climb(grid, (2, 1), 'max', min_rel_improvement=1e-4) == [(2, 1), (2, 2)]
    # zero current value: any strict improvement counts
    grid0 = np.zeros((1, 2)); grid0[0, 1] = 1e-9
    assert tc._greedy_climb(grid0, (0, 0), 'max', min_rel_improvement=1e-4) == [(0, 0), (0, 1)]


HELPER_TESTS = [
    test_squeeze_grid_accepts_2d_and_3d,
    test_greedy_climb_min_rel_improvement_ignores_round_off_plateau,
    test_snap_index_picks_nearest,
    test_greedy_climb_walks_diagonally_to_min,
    test_greedy_climb_max_sense,
    test_greedy_climb_skips_nan_neighbors,
    test_greedy_climb_start_at_optimum_is_single_point,
    test_greedy_climb_nan_start_raises,
    test_greedy_climb_invalid_sense_raises,
]

def _synthetic_results():
    # 6x7 grid (ny=6 over spec_2, nx=7 over spec_1). MPSP-like bowl (min at a
    # corner), IRR-like plane (max at the opposite corner), plus a constant.
    spec_1 = np.linspace(0.0, 30.0, 7)   # x
    spec_2 = np.linspace(0.0, 0.5, 6)    # y
    iy, ix = np.mgrid[0:6, 0:7]
    mpsp = (iy - 0.0) ** 2 + (ix - 0.0) ** 2 + 1.0   # min at (0,0)
    irr = iy.astype(float) + ix.astype(float)         # max at (5,6)
    const = np.full((6, 7), 3.0)
    results = {
        'MPSP': mpsp[None, :, :],
        'IRR': irr[None, :, :],
        'Const': const[None, :, :],
    }
    return results, spec_1, spec_2


def test_derive_levels_ranges_and_errors():
    lv, ct = tc._derive_levels(np.array([1.0, np.nan, 3.0]))
    assert lv[0] == 1.0 and lv[-1] == 3.0
    assert ct[0] == 1.0 and ct[-1] == 3.0
    for bad in (np.array([np.nan, np.nan]), np.array([2.0, 2.0])):
        try:
            tc._derive_levels(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("expected ValueError from _derive_levels")


def test_plot_returns_fig_ax_and_trajectory_data():
    results, s1, s2 = _synthetic_results()
    fig, ax, traj = tc.plot_metric_with_trajectories(
        results, s1, s2, color_metric='MPSP', baseline_point=(30.0, 0.5),
        trajectory_metrics=['MPSP', 'IRR'], senses={'MPSP': 'min', 'IRR': 'max'})
    assert set(traj) == {'MPSP', 'IRR'}
    # baseline snapped to the far corner (ix=6, iy=5)
    assert traj['MPSP']['path_ij'][0] == (5, 6)
    assert traj['MPSP']['optimum_xy'] == (0.0, 0.0)      # min bowl corner
    assert traj['IRR']['n_steps'] == 0                    # baseline already max
    # each path strictly improves per its sense
    mpsp = tc._squeeze_grid(results['MPSP'])
    vals = [mpsp[c] for c in traj['MPSP']['path_ij']]
    assert all(b < a for a, b in zip(vals, vals[1:]))
    # at least one trajectory polyline was drawn, plus a legend
    assert len(ax.lines) >= 1
    assert ax.get_legend() is not None
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_missing_sense_raises():
    results, s1, s2 = _synthetic_results()
    try:
        tc.plot_metric_with_trajectories(
            results, s1, s2, color_metric='MPSP', baseline_point=(0.0, 0.0),
            trajectory_metrics=['IRR'], senses={})
    except KeyError:
        pass
    else:
        raise AssertionError("expected KeyError for missing sense")


def test_plot_constant_color_metric_raises():
    results, s1, s2 = _synthetic_results()
    try:
        tc.plot_metric_with_trajectories(
            results, s1, s2, color_metric='Const', baseline_point=(0.0, 0.0),
            trajectory_metrics=['MPSP'], senses={'MPSP': 'min'})
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for constant color metric")


def test_plot_trajectory_styles_and_legend_kwargs():
    results, s1, s2 = _synthetic_results()
    fig, ax, traj = tc.plot_metric_with_trajectories(
        results, s1, s2, color_metric='MPSP', baseline_point=(30.0, 0.5),
        trajectory_metrics=['MPSP', 'IRR'], senses={'MPSP': 'min', 'IRR': 'max'},
        trajectory_linestyles={'MPSP': '--'}, trajectory_linewidths={'MPSP': 3.0},
        legend_kwargs={'loc': 'lower left'})
    # the MPSP polyline (multi-point line) carries the requested style/width
    polylines = [l for l in ax.lines if len(l.get_xdata()) > 1]
    assert len(polylines) == 1                       # IRR has n_steps == 0
    assert polylines[0].get_linestyle() == '--'
    assert polylines[0].get_linewidth() == 3.0
    # legend entry mirrors the style; unspecified metric keeps the defaults
    by_label = {h.get_label(): h for h in ax.get_legend().legend_handles}
    assert by_label['MPSP'].get_linestyle() == '--'
    assert by_label['MPSP'].get_linewidth() == 3.0
    assert by_label['IRR'].get_linestyle() == '-'
    assert by_label['IRR'].get_linewidth() == 1.6
    assert ax.get_legend()._loc == 3                 # matplotlib code for 'lower left'
    import matplotlib.pyplot as plt
    plt.close(fig)


PLOT_TESTS = [
    test_plot_trajectory_styles_and_legend_kwargs,
    test_derive_levels_ranges_and_errors,
    test_plot_returns_fig_ax_and_trajectory_data,
    test_plot_missing_sense_raises,
    test_plot_constant_color_metric_raises,
]


def _run(tests):
    failures = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except Exception as e:
            failures += 1
            print(f"FAIL {t.__name__}: {type(e).__name__}: {e}")
    return failures


if __name__ == '__main__':
    import sys
    n = _run(HELPER_TESTS + PLOT_TESTS)
    print(f"\n{'ALL PASS' if n == 0 else str(n) + ' FAILED'}")
    sys.exit(1 if n else 0)
