#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""k_13 x isobutanol-inhibition sweep figure for the opt_PI_TRY_informed scenario.

IRR at default prices over the k_13 (ALS capacity) x grouped inhib_isobutanol
multiplier grid of analyses/evaluate_EtOH_k13_inhib_isobutanol.py
(opt_PI_TRY_informed baseline = trial #1912 of the 2026-09-24 TRY-informed PI
(log-tail) relay campaign _rl15c111dc; enzyme burden ON, the scenario's
feeding strategy (1-spike cap / 189.19 / 194.19), no feeding-strategy
optimization), with the grid optimum of IRR, TCI, cell density and the ethanol
and isobutanol titer / yield / productivity marked and annotated. Format and
style follow plots/plot_glycolysis_inhib_ethanol_A.py (its IRR twin); the
drawing lives in plots/_kinetic_sweep_figure.py (loaded by file path).

Renamed 2026-10-08 from plot_k13_inhib_isobutanol_opt_IRR.py when the opt_IRR
scenario (split_12d trial 1602) was removed in favour of opt_PI_TRY_informed.
STALE INPUTS: the sweep CSVs of the 2026-10-01 run are of the OLD opt_IRR point
(max_n=50 prefix); re-run the sweep on opt_PI_TRY_informed (its CSVs carry the
max_n=1 prefix read below) before rendering. The annotation offsets, marker
nudges and IRR colour range below were tuned on the old opt_IRR grid and will
likely need retuning on the regenerated one.

Sim-safe: reads the sweep CSVs only and never imports the package.

Output: analyses/results/publication/Kinetic-sweeps/k13_inhib_isobutanol_opt_PI_TRY_informed_IRR.{png,pdf}
"""

import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    '_kinetic_sweep_figure', os.path.join(HERE, '_kinetic_sweep_figure.py'))
ksf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ksf)

# Grid of the sweep to read: a 40 x 40 run of
# analyses/evaluate_EtOH_k13_inhib_isobutanol.py on opt_PI_TRY_informed (NOT
# yet generated; the 2026-10-01 run, hensmith 6d4776f, was on the old opt_IRR
# point), k_13 over 1e-3x-4x the baseline (4.0 g/L/h for both opt_IRR and
# opt_PI_TRY_informed) and the inhib_isobutanol multiplier over 0.75x-1.5x
# (the metabolic_14d default group band). SPEC_1 / SPEC_2 must match that
# script's linspaces. (The 2026-09-20 opt_IRR run over a 1e-3x-1.5x multiplier
# band is kept under
# results/_superseded/k13_inhib_isobutanol_ib0.001-1.5_2026-09-20/.)
SWEEP_STEPS = (40, 40, 1)
BASELINE_K13 = 4.0
SPEC_1 = np.linspace(1e-3 * BASELINE_K13, 4.0 * BASELINE_K13, SWEEP_STEPS[0])  # k_13 [g/L/h]
SPEC_2 = np.linspace(0.75, 1.5, SWEEP_STEPS[1])  # inhib_isobutanol multiplier [-]

# Optimum markers: shape = metric (as in the glycolysis figure), face colour =
# product (ethanol white, isobutanol grey; economic optima cyan). Offsets were
# set by eye for the old opt_IRR sweep (2026-10-01); retune them on the
# regenerated opt_PI_TRY_informed grid if the optima move.
IBO_FACE = '#9a9a9a'
ETOH_TITER = 'ethanol\ntiter'
IBO_YIELD = 'isobutanol\nyield'
OPTIMA = [
    ('EtOH Titer',        'max', ETOH_TITER,                 '^', 'white',   10, (4, 52),     0.),
    ('EtOH Yield',        'max', 'ethanol yield',            'p', 'white',   10, (57, 27),    0.),
    ('EtOH Productivity', 'max', 'ethanol productivity',     's', 'white',    9, (22, 10),   -0.1),
    ('IBO Titer',         'max', 'isobutanol titer',         '^', IBO_FACE,  10, (9, 4),      0.),
    ('IBO Yield',         'max', IBO_YIELD,                  'p', IBO_FACE,  10, (-24, 44),   0.2),
    ('IBO Productivity',  'max', 'isobutanol productivity',  's', IBO_FACE,   9, (28, 17),    0.),
    ('Cell loading',      'max', 'cell density',             'o', 'white',   10, (10, 130),   0.),
    ('TCI',               'min', 'TCI',                      'p', '#33ccff', 10, (-30, -16), -0.2),
    ('IRR',               'max', 'profitability',            '*', '#33ccff', 14, (40, 30),   -0.2),
]
# the TCI label sits on the dark (IRR < 0) cells: it and its arrow in white
LABEL_COLOR = {'TCI': 'white'}

FIGURE = ksf.SweepFigure(
    output_stem='k13_inhib_isobutanol_opt_PI_TRY_informed_IRR',
    # max_n=1: opt_PI_TRY_informed's 1-spike cap (the old opt_IRR CSVs carry
    # max_n=50), so stale opt_IRR CSVs are never read by mistake
    csv_prefix=f'ibo_{SWEEP_STEPS}_k_13_inhib_Spike_opt=False_max_n=1_',
    spec_1=SPEC_1, spec_2=SPEC_2,
    xlabel=r'$\bf{ALS\ capacity}$, $k_{\mathrm{13}}$' + f' [{ksf.G_PER_L_PER_H}]',
    ylabel=r'$\bf{Isobutanol\ inhibition}$ [× baseline]',
    optima=OPTIMA,
    xlim=(0., 16.), ylim=(0.75, 1.5),
    x_major=4., x_minor_div=4, y_major=0.25, y_minor_div=5,
    # co-located optima fanned apart [pt]: ethanol titer and yield share the
    # cell next to the cell-density corner; isobutanol titer and productivity
    # sit in neighbouring bottom-row cells
    marker_nudge={ETOH_TITER: (0., 8.), 'ethanol yield': (9., -2.),
                  'isobutanol productivity': (-6., 0.), 'isobutanol titer': (6., 0.)},
    label_color=LABEL_COLOR,
    # each arrow leaves its label's near edge, so arrows do not cross the
    # stacked labels: the column right of the bottom-left cluster from the
    # left edge, the upper-left labels from the left / bottom-left
    arrow_relpos={'isobutanol titer': (0., 0.5), 'isobutanol productivity': (0., 0.5),
                  'ethanol yield': (0., 0.5), ETOH_TITER: (0., 0.),
                  'cell density': (0., 0.5), IBO_YIELD: (1., 0.5),
                  'TCI': (1., 0.5), 'baseline': (1., 0.5)},
    color_axis='IRR',
    # IRR axis -15 to 30 %: set for the OLD opt_IRR grid (it topped out at
    # 28.1 %, that baseline 28.0 %); may need revisiting after the re-run on
    # opt_PI_TRY_informed (baseline IRR 0.2792)
    levels=np.arange(-15., 30.00001, 0.5),
    cbar_ticks=np.arange(-15., 30.00001, 5.),
    cbar_minor_step=1.,
    # opt_PI_TRY_informed baseline (k_13 = 4.0 g/L/h, multiplier 1x)
    baseline=(BASELINE_K13, 1.),
    baseline_callout=('baseline', (-10, 34), 0.2),
)

if __name__ == '__main__':
    ksf.main(FIGURE)
