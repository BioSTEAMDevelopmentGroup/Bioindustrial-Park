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
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import contourplots
from biosteam.utils import colors

__all__ = ('plot_metric_with_trajectories', 'JBEI_UCB_colormap')

#: default trajectory line/marker colors, cycled across trajectory metrics
_DEFAULT_TRAJECTORY_COLORS = (
    '#33ccff', '#ff5555', '#ffcc00', '#66cc66',
    '#cc66ff', '#ff9933', '#00cccc', '#ff66b3',
)


def JBEI_UCB_colormap(N_levels=90, reverse=False):
    """The JBEI/UCB contour colormap used across the isobutanol kinetic plots."""
    JBEI_orange = (233 / 255, 83 / 255, 39 / 255)
    UCB_blue = (0 / 255, 38 / 255, 118 / 255)
    UCB_yellow = (253 / 255, 181 / 255, 21 / 255)
    cmap_colors = [UCB_yellow, JBEI_orange, UCB_blue, colors.grey_dark.RGBn]
    if reverse:
        cmap_colors = cmap_colors[::-1]
    return LinearSegmentedColormap.from_list('CABBI', cmap_colors, N_levels)


def _derive_levels(data, n_levels=80, n_cbar_ticks=5):
    """Filled-contour levels and colorbar ticks spanning the non-NaN range of
    `data`. Raises ValueError if there is nothing to contour (all NaN) or no
    range (constant)."""
    d = np.asarray(data, dtype=float)
    d = d[~np.isnan(d)]
    if d.size == 0:
        raise ValueError("color metric has no non-NaN values to contour")
    lo, hi = float(d.min()), float(d.max())
    if lo == hi:
        raise ValueError("color metric is constant; no range to contour")
    w_levels = np.linspace(lo, hi, n_levels + 1)
    cbar_ticks = np.linspace(lo, hi, n_cbar_ticks + 1)
    return w_levels, cbar_ticks


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


def plot_metric_with_trajectories(
        results, spec_1, spec_2, color_metric, baseline_point,
        trajectory_metrics, senses,
        x_label='', y_label='', x_units='', y_units='',
        x_ticks=None, y_ticks=None,
        color_metric_units='', color_metric_label=None,
        cmap=None, w_levels=None, cbar_ticks=None, w_ticks=None,
        extend_cmap='max', cmap_over_color=None,
        trajectory_colors=None,
        trajectory_linestyles=None, trajectory_linewidths=None,
        baseline_marker=('D', 'gray', 6),
        optimum_marker='*', optimum_marker_size=12,
        show_legend=True, legend_kwargs=None, fig_ax=None,
        **contourplot_kwargs):
    """Filled contour of `color_metric` over (`spec_1`, `spec_2`), with a
    baseline marker and one greedy 8-neighbor trajectory (+ optimum marker)
    per metric in `trajectory_metrics`.

    Parameters
    ----------
    results : dict
        {metric_name: array of shape (z, ny, nx) or (ny, nx)}; z is squeezed.
    spec_1, spec_2 : 1-D array_like
        x-axis (len nx) and y-axis (len ny) grid values.
    color_metric : str
        Metric drawn as the filled contour + colorbar.
    baseline_point : (x_value, y_value)
        In data units; snapped to the nearest grid cell (shared start for all
        trajectories).
    trajectory_metrics : list of str
        Metrics to hill-climb and draw.
    senses : dict
        {metric_name: 'min' | 'max'} for every metric in `trajectory_metrics`.
    trajectory_colors, trajectory_linestyles, trajectory_linewidths : dict, optional
        {metric_name: value} overrides for each trajectory's polyline (and its
        legend entry); defaults cycle a fixed color list, solid, 1.6 pt.
        Distinct styles/widths keep coincident path segments distinguishable.
    legend_kwargs : dict, optional
        Overrides for `ax.legend` (default loc='upper right', fontsize=8,
        framealpha=0.9), e.g. {'loc': 'center right'}.

    Returns
    -------
    (fig, ax, trajectory_data) where trajectory_data[metric] =
        {'path_xy': [(x,y),...], 'path_ij': [(iy,ix),...],
         'optimum_xy': (x,y), 'optimum_value': float, 'n_steps': int}.
    """
    spec_1 = np.asarray(spec_1, dtype=float)
    spec_2 = np.asarray(spec_2, dtype=float)

    if color_metric not in results:
        raise KeyError(f"color_metric {color_metric!r} not in results")
    color_grid = _squeeze_grid(results[color_metric])

    for m in trajectory_metrics:
        if m not in results:
            raise KeyError(f"trajectory metric {m!r} not in results")
        if m not in senses:
            raise KeyError(f"no sense ('min'/'max') given for metric {m!r}")

    # Filled-contour levels / colorbar ticks (validates all-NaN / constant).
    derived_levels, derived_cbar = _derive_levels(color_grid)
    if w_levels is None:
        w_levels = derived_levels
    if cbar_ticks is None:
        cbar_ticks = derived_cbar
    if w_ticks is None:
        interior = list(np.asarray(cbar_ticks, dtype=float)[1:-1])
        w_ticks = interior if interior else [float(np.asarray(cbar_ticks)[0])]
    if cmap is None:
        cmap = JBEI_UCB_colormap()
    if cmap_over_color is None:
        cmap_over_color = colors.grey_dark.shade(8).RGBn
    if x_ticks is None:
        x_ticks = [float(spec_1[0]), float(spec_1[-1])]
    if y_ticks is None:
        y_ticks = [float(spec_2[0]), float(spec_2[-1])]
    if color_metric_label is None:
        color_metric_label = color_metric

    if fig_ax is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = fig_ax

    # Base contour + colorbar via contourplots (approach A). include_top_bar
    # MUST be False for fig_ax_to_use to be honored (else it makes its own axes).
    contourplots.animated_contourplot(
        w_data_vs_x_y_at_multiple_z=color_grid[None, :, :],
        x_data=spec_1, y_data=spec_2, z_data=[0],
        x_label=x_label, y_label=y_label, z_label='',
        w_label=color_metric_label,
        x_ticks=x_ticks, y_ticks=y_ticks, z_ticks=[0],
        w_levels=w_levels, w_ticks=w_ticks,
        x_units=x_units, y_units=y_units, z_units='',
        w_units=color_metric_units,
        cmap=cmap, extend_cmap=extend_cmap, cmap_over_color=cmap_over_color,
        cbar_ticks=cbar_ticks,
        include_top_bar=False, include_cbar=True,
        fig_ax_to_use=(fig, ax),
        **contourplot_kwargs,
    )

    ix0 = _snap_index(spec_1, baseline_point[0])
    iy0 = _snap_index(spec_2, baseline_point[1])

    if trajectory_colors is None:
        trajectory_colors = {
            m: _DEFAULT_TRAJECTORY_COLORS[i % len(_DEFAULT_TRAJECTORY_COLORS)]
            for i, m in enumerate(trajectory_metrics)}

    trajectory_data = {}
    legend_handles = []
    for m in trajectory_metrics:
        grid = _squeeze_grid(results[m])
        path_ij = _greedy_climb(grid, (iy0, ix0), senses[m])
        path_xy = [(float(spec_1[ix]), float(spec_2[iy])) for (iy, ix) in path_ij]
        color = trajectory_colors[m]
        linestyle = (trajectory_linestyles or {}).get(m, '-')
        linewidth = (trajectory_linewidths or {}).get(m, 1.6)
        xs = [p[0] for p in path_xy]
        ys = [p[1] for p in path_xy]
        if len(path_xy) > 1:
            ax.plot(xs, ys, linestyle=linestyle, color=color,
                    linewidth=linewidth, zorder=400, clip_on=False)
        ax.plot(xs[-1], ys[-1], linestyle='None', marker=optimum_marker,
                markerfacecolor=color, markeredgecolor='k', markeredgewidth=0.6,
                markersize=optimum_marker_size, zorder=600, clip_on=False)
        opt_iy, opt_ix = path_ij[-1]
        trajectory_data[m] = {
            'path_xy': path_xy,
            'path_ij': path_ij,
            'optimum_xy': path_xy[-1],
            'optimum_value': float(grid[opt_iy, opt_ix]),
            'n_steps': len(path_xy) - 1,
        }
        legend_handles.append(Line2D(
            [0], [0], color=color, marker=optimum_marker, markerfacecolor=color,
            markeredgecolor='k', linestyle=linestyle, linewidth=linewidth,
            label=m))

    b_shape, b_color, b_size = baseline_marker
    ax.plot(float(spec_1[ix0]), float(spec_2[iy0]), linestyle='None',
            marker=b_shape, markerfacecolor=b_color, markeredgecolor='k',
            markeredgewidth=0.8, markersize=b_size, zorder=550, clip_on=False)
    baseline_handle = Line2D(
        [0], [0], color='none', marker=b_shape, markerfacecolor=b_color,
        markeredgecolor='k', linestyle='None', label='baseline')

    if show_legend:
        legend_kw = dict(loc='upper right', fontsize=8, framealpha=0.9)
        if legend_kwargs:
            legend_kw.update(legend_kwargs)
        ax.legend(handles=[baseline_handle, *legend_handles], **legend_kw)

    return fig, ax, trajectory_data
