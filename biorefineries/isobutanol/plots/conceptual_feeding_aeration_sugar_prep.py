#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Publication-ready conceptual diagram (Nature Communications format) of how
sugar-solution preparation, fed-batch feeding, and aeration are handled
together in the corn ethanol-isobutanol biorefinery model.

Three panels, concepts rather than code:

* **a** — process + control schematic: the splitter feeding parallel
  initial-feed and spike-feed conditioning trains (multi-effect evaporator,
  dilution-water mixer, cooler), the fed-batch fermentor driven by a
  kinetic model, the demand-pull compressed-air aeration loop, and the
  fed-batch strategy specification that tunes all of it.
* **b** — in-silico fed-batch control logic within one batch: threshold-
  triggered glucose spikes toward a target concentration (capped at
  N_max), the stage-1 aerobic growth window that sizes aeration, and the
  glucose-depletion policy that sets the batch time tau (capped at
  tau_max). Trajectories are schematic, not simulation output.
* **c** — simulation orchestration: the load-specifications + simulate
  fixed-point sweep (actuator root-solves, kinetic integration with the
  spike-retry loop, feed-split and air-supply closure, full-flowsheet
  convergence) and the failure-recovery cascade.

Standalone on purpose: imports only numpy/matplotlib (NOT
``biorefineries.isobutanol``), so running it triggers no baseline
simulation and touches no numba cache. Fonts embed as TrueType (fonttype
42) as required for journal submission; the figure is sized to the
double-column width (180 mm). Colors are the Okabe-Ito colorblind-safe
palette (Wong 2011, Nature Methods), assigned by role: blue = sugar
preparation, green = fermentation, orange = aeration, reddish purple =
control/specification, vermillion = failure recovery.

Run (any Python with matplotlib; headless-safe)::

    python plots/conceptual_feeding_aeration_sugar_prep.py

Outputs ``conceptual_feeding_aeration_sugar_prep.pdf`` (vector, for
submission) and ``.png`` (600 dpi, for review) under ``plots/results/``.
"""

import os

import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

__all__ = ('make_figure',)

# --------------------------------------------------------------------------
# Style constants
# --------------------------------------------------------------------------

# Okabe-Ito, assigned by role (fixed order, never cycled)
C_PREP = '#0072B2'      # blue — sugar-solution preparation
C_FERM = '#009E73'      # bluish green — fermentation
C_AIR = '#E69F00'       # orange — aeration
C_CTRL = '#CC79A7'      # reddish purple — control / specification
C_FAIL = '#D55E00'      # vermillion — failure recovery
C_INK = '#1A1A1A'       # primary text
C_INK2 = '#4D4D4D'      # secondary text / stream arrows
C_GRAY = '#999999'      # biomass line, dashed guides
C_AIR_TXT = '#8A6100'   # darker orange for small text (contrast)

TINT = {C_PREP: '#E9F1F8', C_FERM: '#E7F5F0', C_AIR: '#FBF1DD',
        C_CTRL: '#F8ECF3', C_FAIL: '#FBE9DF', C_INK2: '#F0F0F0'}

FS_BASE = 6.2       # box titles (pt)
FS_SUB = 5.4        # box sublines / annotations
FS_PANEL = 8.0      # bold panel letters
FS_ZONE = 6.4       # zone headers

MM = 1 / 25.4
FIG_W = 180 * MM    # Nature Communications double-column width
FIG_H = 168 * MM


def _set_style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.default': 'regular',   # math in the text font, not Computer Modern
        'pdf.fonttype': 42,              # embed TrueType (journal requirement)
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
        'axes.linewidth': 0.6,
        'text.color': C_INK,
        'axes.edgecolor': C_INK2,
        'axes.labelcolor': C_INK,
        'xtick.color': C_INK2,
        'ytick.color': C_INK2,
    })


# --------------------------------------------------------------------------
# Drawing helpers
# --------------------------------------------------------------------------

def box(ax, x, y, w, h, title, sub=None, color=C_PREP, lw=0.9,
        title_fs=FS_BASE, sub_fs=FS_SUB, dashed=False, rounding=1.1,
        stack=False):
    """Rounded process box with a bold title and optional sublines.

    ``stack=True`` anchors the title to the box top and the sublines to the
    box bottom (for panels whose y-scale is finer than panel a's).
    """
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f'round,pad=0,rounding_size={rounding}',
        facecolor=TINT.get(color, '#FFFFFF'), edgecolor=color,
        linewidth=lw, linestyle=(0, (3.2, 1.8)) if dashed else 'solid',
        zorder=3)
    ax.add_patch(patch)
    cx, cy = x + w / 2, y + h / 2
    if sub and stack:
        ax.text(cx, y + h - 1.1, title, ha='center', va='top',
                fontsize=title_fs, fontweight='bold', color=C_INK, zorder=4)
        ax.text(cx, y + 0.9, sub, ha='center', va='bottom',
                fontsize=sub_fs, color=C_INK2, zorder=4, linespacing=1.25)
    elif sub:
        ax.text(cx, cy + 0.16 * h, title, ha='center', va='center',
                fontsize=title_fs, fontweight='bold', color=C_INK, zorder=4)
        ax.text(cx, cy - 0.21 * h, sub, ha='center', va='center',
                fontsize=sub_fs, color=C_INK2, zorder=4, linespacing=1.25)
    else:
        ax.text(cx, cy, title, ha='center', va='center', fontsize=title_fs,
                fontweight='bold', color=C_INK, zorder=4)
    return patch


def arrow(ax, p0, p1, color=C_INK2, lw=1.0, dashed=False, zorder=2,
          head=3.2):
    """Stream (solid) or control/information (dashed) arrow."""
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle=f'-|>,head_width={head*0.55},head_length={head}',
        color=color, linewidth=lw,
        linestyle=(0, (3.0, 1.8)) if dashed else 'solid',
        shrinkA=0.0, shrinkB=0.0, zorder=zorder, mutation_scale=1.0))


def elbow(ax, pts, color=C_INK2, lw=1.0, dashed=False, zorder=2, head=3.2):
    """Polyline through ``pts`` with an arrowhead on the final segment."""
    xs, ys = zip(*pts)
    ax.plot(xs[:-1], ys[:-1], color=color, lw=lw, zorder=zorder,
            linestyle=(0, (3.0, 1.8)) if dashed else 'solid',
            solid_capstyle='round')
    arrow(ax, pts[-2], pts[-1], color=color, lw=lw, dashed=dashed,
          zorder=zorder, head=head)


def slabel(ax, x, y, text, color=C_INK2, fs=FS_SUB, ha='center', va='center',
           style='normal', weight='normal'):
    ax.text(x, y, text, ha=ha, va=va, fontsize=fs, color=color,
            style=style, fontweight=weight, zorder=5, linespacing=1.25)


# --------------------------------------------------------------------------
# Panel a — process + control schematic
# --------------------------------------------------------------------------

def draw_panel_a(ax):
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 57)
    ax.axis('off')

    # ---- zone headers -----------------------------------------------------
    slabel(ax, 39, 55.0, 'SUGAR-SOLUTION PREPARATION', color=C_PREP,
           fs=FS_ZONE, weight='bold')
    slabel(ax, 84.5, 55.0, 'FED-BATCH FERMENTATION', color=C_FERM,
           fs=FS_ZONE, weight='bold')
    slabel(ax, 43, 10.4, 'AERATION (stage 1 only)', color=C_AIR,
           fs=FS_ZONE, weight='bold')

    # ---- feed + splitter --------------------------------------------------
    box(ax, 1.5, 23.5, 11.5, 9, 'Saccharified\ncorn mash',
        sub='glucose-rich\nslurry', color=C_PREP)
    box(ax, 17.5, 23.5, 8.5, 9, 'Splitter', sub='split $\\phi$', color=C_PREP)
    arrow(ax, (13, 28), (17.5, 28))

    # splitter -> the two trains
    elbow(ax, [(26, 30.0), (28, 30.0), (28, 42.0), (31, 42.0)], lw=1.0)
    elbow(ax, [(26, 26.0), (28, 26.0), (28, 15.5), (31, 15.5)], lw=1.0)

    # ---- initial-feed train (top) ----------------------------------------
    yi = 38.0   # box bottom
    box(ax, 31, yi, 11.5, 8, 'Evaporator',
        sub='multi-effect\nknob: vapor frac.', color=C_PREP)
    box(ax, 46.5, yi, 11.5, 8, 'Dilution mixer',
        sub='knob: water :\nsugar ratio', color=C_PREP)
    box(ax, 62, yi, 8.5, 8, 'Cooler', sub='to 32 °C', color=C_PREP)
    arrow(ax, (42.5, yi + 4), (46.5, yi + 4))
    arrow(ax, (58, yi + 4), (62, yi + 4))
    # condensate out (up), dilution water in (down)
    arrow(ax, (36.75, yi + 8), (36.75, yi + 12.2))
    slabel(ax, 36.75, yi + 13.6, 'condensate  →  WWT')
    arrow(ax, (52.25, yi + 12.2), (52.25, yi + 8))
    slabel(ax, 52.25, yi + 13.6, 'dilution water')

    # initial charge into fermentor
    arrow(ax, (70.5, 42), (74.5, 42), lw=1.1)
    slabel(ax, 71, 47.6, 'initial charge\nat $C_{target}$', fs=FS_SUB)

    # ---- spike train (bottom) --------------------------------------------
    ys = 11.5
    box(ax, 31, ys, 11.5, 8, 'Evaporator',
        sub='multi-effect\nknob: vapor frac.', color=C_PREP)
    box(ax, 46.5, ys, 11.5, 8, 'Dilution mixer',
        sub='knob: water :\nsugar ratio', color=C_PREP)
    box(ax, 62, ys, 8.5, 8, 'Cooler', sub='to 32 °C', color=C_PREP)
    arrow(ax, (42.5, ys + 4), (46.5, ys + 4))
    arrow(ax, (58, ys + 4), (62, ys + 4))

    # spike feed into fermentor
    arrow(ax, (70.5, 15.5), (74.5, 15.5), lw=1.1)
    slabel(ax, 68, 10.2, 'spike feed at $C_{spike}$', fs=FS_SUB)

    # ---- fermentor --------------------------------------------------------
    fx, fy, fw, fh = 74.5, 13.5, 15.5, 28.5
    box(ax, fx, fy, fw, fh,
        'Fed-batch\nfermentor',
        sub='S. cerevisiae\nkinetic model (ODEs)\nethanol + isobutanol\nco-production',
        color=C_FERM, lw=1.2, rounding=1.6)
    # vent + broth
    arrow(ax, (83, fy + fh), (83, fy + fh + 4.3), lw=1.0)
    slabel(ax, 83, 48.9, 'CO$_2$-rich vent\n→ scrubber')
    arrow(ax, (90, 21.5), (98.5, 21.5), lw=1.1)
    slabel(ax, 94.3, 19.3, 'broth →\nseparation')
    # nutrient feed-forward triggers (bottom right)
    slabel(ax, 93.5, 6.9,
           'seed yeast $\\propto$ mash flow\n'
           'glucoamylase $\\propto$ dry mash\n'
           'NH$_3$ $\\propto$ yeast grown\n(flows re-set every run)')
    arrow(ax, (93.5, 10.9), (88.5, 13.4), dashed=True, color=C_CTRL, lw=0.9)

    # ---- aeration loop ----------------------------------------------------
    box(ax, 59.5, 1.0, 11.5, 6.8, 'Compressor', sub='air, isothermal',
        color=C_AIR)
    box(ax, 74.5, 1.0, 8.5, 6.8, 'Valve', sub='let-down', color=C_AIR)
    slabel(ax, 56, 4.4, 'air')
    arrow(ax, (57.3, 4.4), (59.5, 4.4))
    arrow(ax, (71, 4.4), (74.5, 4.4))
    elbow(ax, [(83, 4.4), (86.5, 4.4), (86.5, 13.5)], lw=1.0)
    # demand-pull re-simulation (control)
    arrow(ax, (77.5, 13.5), (77.5, 7.9), color=C_AIR, lw=0.9, dashed=True)
    slabel(ax, 43, 4.6,
           'sized to 2× the stage-1 O$_2$ demand\n'
           '(aerobic growth, ends at cell-density\n'
           'target or time cap); valve → compressor\n'
           're-simulated after every fermentation run',
           fs=FS_SUB, color=C_AIR_TXT)

    # ---- fed-batch strategy specification (control hub) -------------------
    hx, hy, hw, hh = 31, 24.5, 39.5, 8.5
    box(ax, hx, hy, hw, hh, 'Fed-batch strategy specification',
        sub='$C_{threshold}$ < $C_{target}$ < $C_{spike}$ ;  '
            '$\\tau_{max}$ ;  spike cap $N_{max}$',
        color=C_CTRL, lw=1.1, dashed=True)
    # control arrows: split, both trains (actuator pairs), kinetic model
    arrow(ax, (hx, 28.75), (26.2, 28.75), color=C_CTRL, lw=0.9, dashed=True)
    for xa in (36.75, 52.25):
        arrow(ax, (xa, hy + hh), (xa, yi), color=C_CTRL, lw=0.9, dashed=True)
        arrow(ax, (xa, hy), (xa, ys + 8), color=C_CTRL, lw=0.9, dashed=True)
    slabel(ax, 44.5, 35.7, 'root-solved to hit $C_{target}$ / $C_{spike}$',
           color=C_CTRL, fs=FS_SUB)
    arrow(ax, (hx + hw, 28.75), (fx, 28.75), color=C_CTRL, lw=0.9,
          dashed=True)
    slabel(ax, 72.5, 31.0, 'set-\npoints', color=C_CTRL, fs=FS_SUB)

    # legend for arrow types
    ax.plot([2.5, 6.5], [7.2, 7.2], color=C_INK2, lw=1.0)
    slabel(ax, 7.5, 7.2, 'material stream', ha='left', fs=FS_SUB)
    ax.plot([2.5, 6.5], [4.2, 4.2], color=C_CTRL, lw=1.0,
            linestyle=(0, (3.0, 1.8)))
    slabel(ax, 7.5, 4.2, 'control / specification', ha='left', fs=FS_SUB)
    ax.plot([2.5, 6.5], [1.2, 1.2], color=C_AIR, lw=1.0,
            linestyle=(0, (3.0, 1.8)))
    slabel(ax, 7.5, 1.2, 'demand-pull sizing', ha='left', fs=FS_SUB)


# --------------------------------------------------------------------------
# Panel b — in-batch fed-batch control logic (schematic trajectories)
# --------------------------------------------------------------------------

def _glucose_sawtooth(c_target, c_threshold, spike_times, tau):
    """Piecewise-exponential decline segments between spikes (schematic)."""
    knots_t, knots_c = [], []
    seg_starts = [0.0] + list(spike_times)
    seg_ends = list(spike_times) + [tau]
    end_concs = [c_threshold] * len(spike_times) + [0.0]
    for t0, t1, c1 in zip(seg_starts, seg_ends, end_concs):
        t = np.linspace(t0, t1, 60)
        # convex decline from c_target to c1 (consumption accelerates
        # early in the batch as biomass grows; near-linear later)
        frac = (t - t0) / (t1 - t0)
        shape = 1 - (1 - frac) ** 1.6 if t0 < 25 else frac
        c = c_target - (c_target - c1) * shape
        knots_t.append(t)
        knots_c.append(c)
    return knots_t, knots_c


def draw_panel_b(ax):
    c_target, c_threshold = 1.0, 0.55
    spike_times = [18.0, 31.0, 42.0]
    tau, tau_max, t_air = 62.0, 74.0, 12.0

    ax.set_xlim(0, 80)
    ax.set_ylim(0, 1.45)

    # stage-1 aerobic window
    ax.axvspan(0, t_air, color=TINT[C_AIR], zorder=0)
    ax.plot([t_air, t_air], [0, 1.45], color=C_AIR, lw=0.7,
            linestyle=(0, (3.0, 1.8)), zorder=1)
    slabel(ax, t_air / 2, 1.36, 'stage 1\naerobic\ngrowth', color=C_AIR_TXT,
           fs=FS_SUB)
    slabel(ax, t_air + 2, 1.40, 'stage 2 · anaerobic production',
           color=C_INK2, fs=FS_SUB, ha='left')

    # target / threshold guides (labeled at left, inside the plot)
    for c_val in (c_target, c_threshold):
        ax.plot([0, tau], [c_val, c_val], color=C_GRAY, lw=0.6,
                linestyle=(0, (1.2, 1.6)), zorder=1)
    slabel(ax, 1.2, 1.03, '$C_{target}$', color=C_INK2, ha='left',
           va='bottom', fs=FS_SUB)
    slabel(ax, 1.2, 0.52, '$C_{threshold}$', color=C_INK2, ha='left',
           va='top', fs=FS_SUB)

    # glucose sawtooth
    ts, cs = _glucose_sawtooth(c_target, c_threshold, spike_times, tau)
    for i, (t, c) in enumerate(zip(ts, cs)):
        ax.plot(t, c, color=C_PREP, lw=1.3, zorder=3)
        if i < len(spike_times):   # the vertical spike jump
            ax.plot([t[-1], t[-1]], [c_threshold, c_target], color=C_PREP,
                    lw=1.3, zorder=3)
    ax.plot(spike_times, [c_threshold] * len(spike_times), 'o', ms=3.4,
            markerfacecolor='white', markeredgecolor=C_PREP,
            markeredgewidth=1.0, zorder=4)

    # annotation band (top): spike rule (left), batch-end rule (right)
    slabel(ax, 1.5, 1.175,
           'glucose $\\leq C_{threshold}$ → spike restores $C_{target}$\n'
           '(at most $N_{max}$ spikes; then glucose runs out)',
           color=C_PREP, fs=FS_SUB, ha='left')
    ax.annotate('', xy=(17.6, 0.60), xytext=(15.5, 1.10),
                arrowprops=dict(arrowstyle='-|>', color=C_PREP, lw=0.7,
                                shrinkA=0, shrinkB=0))
    slabel(ax, 47.5, 1.175,
           'batch ends: $\\tau$ = first\nglucose minimum (≤ $\\tau_{max}$)',
           color=C_INK, fs=FS_SUB, ha='left')
    ax.annotate('', xy=(61.5, 1.06), xytext=(57.5, 1.115),
                arrowprops=dict(arrowstyle='-|>', color=C_INK2, lw=0.7,
                                shrinkA=0, shrinkB=0))

    # biomass, ethanol, isobutanol (schematic)
    t = np.linspace(0, tau, 300)
    biomass = 0.04 + 0.36 / (1 + np.exp(-(t - 6.0) / 2.6))
    ethanol = 0.9 * (t / tau) ** 1.25 * (1 - 0.25 * np.exp(-t / 8))
    ibo = 0.30 * (t / tau) ** 1.45
    ax.plot(t, biomass, color=C_GRAY, lw=1.1, zorder=2)
    ax.plot(t, ethanol, color=C_FERM, lw=1.3, zorder=2)
    ax.plot(t, ibo, color=C_FAIL, lw=1.3, zorder=2)

    # direct labels (right of curve ends / beside the curve)
    slabel(ax, tau + 1.2, float(biomass[-1]), 'cells', color=C_INK2,
           ha='left', fs=FS_SUB)
    slabel(ax, tau + 1.2, float(ethanol[-1]), 'ethanol', color=C_FERM,
           ha='left', fs=FS_SUB, weight='bold')
    slabel(ax, tau + 1.2, float(ibo[-1]), 'isobutanol', color=C_FAIL,
           ha='left', fs=FS_SUB, weight='bold')
    slabel(ax, 13, 0.44, 'glucose', color=C_PREP, fs=FS_SUB, ha='left',
           weight='bold')

    # tau and tau_max
    ax.plot([tau, tau], [0, 1.05], color=C_INK2, lw=0.8,
            linestyle=(0, (4.0, 2.0)), zorder=2)
    ax.plot([tau_max, tau_max], [0, 1.05], color=C_GRAY, lw=0.7,
            linestyle=(0, (1.2, 1.6)), zorder=2)
    slabel(ax, tau_max + 1.2, 0.99, '$\\tau_{max}$\ncap', color=C_INK2,
           fs=FS_SUB, ha='left')

    # axes cosmetics: schematic — no numeric ticks
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.set_xlabel('Time in batch (schematic)', fontsize=FS_BASE)
    ax.set_ylabel('Concentration (schematic)', fontsize=FS_BASE)


# --------------------------------------------------------------------------
# Panel c — simulation orchestration flowchart
# --------------------------------------------------------------------------

def draw_panel_c(ax):
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')

    bx, bw = 8, 56       # main-chain box geometry
    steps = [
        # (y_bottom, height, title, subline, color)
        (89, 11, 'Specify feeding strategy',
         '$C_{target}$, $C_{threshold}$, $C_{spike}$, $\\tau_{max}$, '
         '$N_{max}$\n(call arguments override stored values)', C_CTRL),
        (77, 8.5, 'Impose set-points on kinetic model',
         'spike trigger / target / concentration + cap', C_CTRL),
        (62.5, 11, 'Tune feed-train actuators',
         'root-solve evaporation, then dilution,\n'
         'to hit $C_{target}$ (initial) and $C_{spike}$ (spike)', C_PREP),
        (48, 11, 'Simulate fermentation kinetics',
         'integrate ODEs · retry with fewer spikes\n'
         'until glucose exhausted · select $\\tau$', C_FERM),
        (33.5, 11, 'Close feed and air balances',
         'split $\\phi$ from realized spike volume;\n'
         'air chain re-sized at fresh O$_2$ demand', C_AIR),
        (21.5, 8.5, 'Simulate full flowsheet',
         'recycles · wastewater · boiler · HX network', C_INK2),
    ]
    for y, h, title, sub, color in steps:
        box(ax, bx, y, bw, h, title, sub=sub, color=color, rounding=1.6,
            stack=True)
    # connecting arrows, edge-to-edge
    for i in range(len(steps) - 1):
        y_top = steps[i][0]
        y_next = steps[i + 1][0] + steps[i + 1][1]
        arrow(ax, (bx + bw / 2, y_top), (bx + bw / 2, y_next), lw=0.9)

    # decision + converged
    dy, dh = 12.5, 6
    box(ax, bx + 6, dy, bw - 12, dh, 'flowsheet state drift ≤ 0.01%?',
        color=C_INK2, rounding=2.4)
    arrow(ax, (bx + bw / 2, 21.5), (bx + bw / 2, dy + dh), lw=0.9)
    box(ax, bx + 6, 1, bw - 12, 8, 'Converged',
        sub='≤ 5 sweeps · then read the TEA', color=C_FERM,
        rounding=2.4, stack=True)
    arrow(ax, (bx + bw / 2, dy), (bx + bw / 2, 9), lw=0.9)
    slabel(ax, bx + bw / 2 + 2, 10.7, 'yes', ha='left', fs=FS_SUB)
    # "no" loop back to step 2
    elbow(ax, [(bx + 6, dy + dh / 2), (2.5, dy + dh / 2), (2.5, 81.25),
               (bx, 81.25)], lw=0.9)
    slabel(ax, 4.8, 23, 'no:\nsweep\nagain', ha='center', fs=FS_SUB)

    # failure-recovery cascade (right column)
    rx, rw = 70, 28
    slabel(ax, rx + rw / 2, 63.2, 'on failure', color=C_FAIL, fs=FS_ZONE,
           weight='bold')
    slabel(ax, rx + rw / 2, 60.3, '(any simulation error)', color=C_FAIL,
           fs=FS_SUB)
    recovery = [
        (49, 'retry at reduced $\\tau_{max}$,\nthen restore'),
        (38, 'tighten integrator rtol;\nswitch integrator if needed'),
        (27, 'reset flowsheet caches,\nempty recycles'),
        (16, 'switch flowsheet solver\n(fixed-point → Aitken)'),
    ]
    for y, text in recovery:
        box(ax, rx, y, rw, 8.5, '', color=C_FAIL, rounding=1.6)
        slabel(ax, rx + rw / 2, y + 4.25, text, color=C_INK, fs=FS_SUB)
    for (y1, _), (y0, _) in zip(recovery[:-1], recovery[1:]):
        arrow(ax, (rx + rw / 2, y1), (rx + rw / 2, y0 + 8.5), color=C_FAIL,
              lw=0.8)
    # trigger from the fermentation-simulation step; return to the sweep
    arrow(ax, (bx + bw, 53.5), (rx, 53.5), color=C_FAIL, lw=0.8,
          dashed=True)
    elbow(ax, [(rx + rw / 2, 16), (rx + rw / 2, dy + dh / 2),
               (bx + bw - 6, dy + dh / 2)],
          color=C_FAIL, lw=0.8, dashed=True)
    slabel(ax, 68, 13.0, 're-enter sweep', color=C_FAIL, fs=FS_SUB)


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------

def make_figure():
    _set_style()
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    ax_a = fig.add_axes([0.005, 0.555, 0.99, 0.425])
    ax_b = fig.add_axes([0.065, 0.055, 0.44, 0.42])
    ax_c = fig.add_axes([0.565, 0.02, 0.43, 0.485])

    draw_panel_a(ax_a)
    draw_panel_b(ax_b)
    draw_panel_c(ax_c)

    fig.text(0.005, 0.99, 'a', fontsize=FS_PANEL, fontweight='bold',
             va='top')
    fig.text(0.005, 0.515, 'b', fontsize=FS_PANEL, fontweight='bold',
             va='top')
    fig.text(0.53, 0.515, 'c', fontsize=FS_PANEL, fontweight='bold',
             va='top')
    return fig


def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'results')
    os.makedirs(out_dir, exist_ok=True)
    fig = make_figure()
    stem = os.path.join(out_dir, 'conceptual_feeding_aeration_sugar_prep')
    fig.savefig(stem + '.pdf')
    fig.savefig(stem + '.png', dpi=600, facecolor='white')
    plt.close(fig)
    print(f'Saved {stem}.pdf and {stem}.png')


if __name__ == '__main__':
    main()
