#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Contour plot of a 2-D kinetic-parameter sweep for one metric, overlaid with
greedy (8-neighbor hill-climb) optimization trajectories from a baseline point
to each of several metrics' local optima.

Consumes existing sweep results only (a {metric_name: array} dict plus the two
1-D spec grids); it runs no biorefinery simulation.
"""
import numpy as np

__all__ = ('plot_metric_with_trajectories',)


def _squeeze_grid(arr):
    """Return a 2-D (ny, nx) float array from a metric result that is either
    (ny, nx) or (z, ny, nx). A z axis of any length uses its first slice."""
    a = np.asarray(arr, dtype=float)
    if a.ndim == 3:
        a = a[0]
    elif a.ndim != 2:
        raise ValueError(f"metric array must be 2-D or 3-D, got shape {a.shape}")
    return a


def _snap_index(spec, value):
    """Index of the entry in 1-D `spec` nearest to `value`."""
    spec = np.asarray(spec, dtype=float)
    return int(np.argmin(np.abs(spec - value)))


def _greedy_climb(grid, start, sense, max_steps=None):
    """8-neighbor (Moore) hill-climb over a 2-D `grid`, returning the list of
    visited `(iy, ix)` cells from `start` to the greedy local optimum.

    At each step the up-to-8 in-bounds neighbors are considered; NaN neighbors
    are skipped; the strictly-best-improving neighbor (per `sense`) is taken;
    the climb stops when no neighbor improves. `sense` is 'min' or 'max'. Ties
    are broken toward the smallest `(iy, ix)` for determinism. Raises
    ValueError if `sense` is invalid or the start cell is NaN.
    """
    grid = np.asarray(grid, dtype=float)
    ny, nx = grid.shape
    if sense not in ('min', 'max'):
        raise ValueError(f"sense must be 'min' or 'max', got {sense!r}")
    iy, ix = int(start[0]), int(start[1])
    if np.isnan(grid[iy, ix]):
        raise ValueError(
            f"baseline cell {(iy, ix)} is NaN; cannot start a trajectory there")
    if max_steps is None:
        max_steps = ny * nx
    improves = (lambda cand, ref: cand < ref) if sense == 'min' \
        else (lambda cand, ref: cand > ref)
    path = [(iy, ix)]
    for _ in range(max_steps):
        best = None
        best_val = grid[iy, ix]
        for diy in (-1, 0, 1):
            for dix in (-1, 0, 1):
                if diy == 0 and dix == 0:
                    continue
                jy, jx = iy + diy, ix + dix
                if not (0 <= jy < ny and 0 <= jx < nx):
                    continue
                v = grid[jy, jx]
                if np.isnan(v):
                    continue
                if improves(v, best_val):
                    best_val = v
                    best = (jy, jx)
        if best is None:
            break
        iy, ix = best
        path.append((iy, ix))
    return path
