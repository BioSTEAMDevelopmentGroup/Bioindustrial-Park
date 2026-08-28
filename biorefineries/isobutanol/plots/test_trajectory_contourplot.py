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


HELPER_TESTS = [
    test_squeeze_grid_accepts_2d_and_3d,
    test_snap_index_picks_nearest,
    test_greedy_climb_walks_diagonally_to_min,
    test_greedy_climb_max_sense,
    test_greedy_climb_skips_nan_neighbors,
    test_greedy_climb_start_at_optimum_is_single_point,
    test_greedy_climb_nan_start_raises,
    test_greedy_climb_invalid_sense_raises,
]

PLOT_TESTS = []  # populated in Task 2


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
