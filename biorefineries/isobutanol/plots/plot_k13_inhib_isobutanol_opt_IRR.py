#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""k_13 x isobutanol-inhibition sweep figure for the opt_IRR scenario.

IRR at default prices over the k_13 (ALS capacity) x grouped inhib_isobutanol
multiplier grid of analyses/evaluate_EtOH_k13_inhib_isobutanol.py (opt_IRR
baseline = split_12d trial 1602, enzyme burden ON, opt_IRR's feeding strategy,
no feeding-strategy optimization), with the grid optimum of IRR, TCI, cell
density and the ethanol and isobutanol titer / yield / productivity marked and
annotated. Format and style follow plots/plot_glycolysis_inhib_ethanol_A.py
(its IRR twin); the drawing lives in plots/_kinetic_sweep_figure.py (loaded by
file path).

Sim-safe: reads the sweep CSVs only and never imports the package.

Output: analyses/results/publication/Kinetic-sweeps/k13_inhib_isobutanol_opt_IRR_IRR.{png,pdf}
"""

import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    '_kinetic_sweep_figure', os.path.join(HERE, '_kinetic_sweep_figure.py'))
ksf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ksf)

# Grid of the sweep to read: the 2026-09-20 40 x 40 run of
# analyses/evaluate_EtOH_k13_inhib_isobutanol.py (commit 248f75f4), k_13 over
# 1e-3x-4x the opt_IRR baseline (4.0 g/L/h) and the inhib_isobutanol
# multiplier over 1e-3x-1.5x. The script has since narrowed the multiplier
# band to 0.75x-1.5x (2fb99458) but was not re-run; SPEC_2 must be changed
# with the CSVs if it is.
SWEEP_STEPS = (40, 40, 1)
BASELINE_K13 = 4.0
SPEC_1 = np.linspace(1e-3 * BASELINE_K13, 4.0 * BASELINE_K13, SWEEP_STEPS[0])  # k_13 [g/L/h]
SPEC_2 = np.linspace(1e-3, 1.5, SWEEP_STEPS[1])  # inhib_isobutanol multiplier [-]

# Optimum markers: shape = metric (as in the glycolysis figure), face colour =
# product (ethanol white, isobutanol grey; economic optima cyan). Offsets are
# set by eye for the current sweep; retune them if the optima move.
IBO_FACE = '#9a9a9a'
OPTIMA = [
    ('EtOH Titer',        'max', 'ethanol titer',            '^', 'white',   10, (22, 24),   -0.2),
    ('EtOH Yield',        'max', 'ethanol yield',            'p', 'white',   10, (16, 32),   -0.3),
    ('EtOH Productivity', 'max', 'ethanol productivity',     's', 'white',    9, (22, 10),   -0.1),
    ('IBO Titer',         'max', 'isobutanol\ntiter',        '^', IBO_FACE,  10, (-40, 40),   0.2),
    ('IBO Yield',         'max', 'isobutanol\nyield',        'p', IBO_FACE,  10, (-30, 44),   0.2),
    ('IBO Productivity',  'max', 'isobutanol productivity',  's', IBO_FACE,   9, (56, -4),   -0.1),
    ('Cell loading',      'max', 'cell density',             'o', 'white',   10, (-46, 26),   0.3),
    ('TCI',               'min', 'TCI',                      'p', '#33ccff', 10, (-40, 32),   0.2),
    ('IRR',               'max', 'profitability',            '*', '#33ccff', 14, (36, 8),    -0.2),
]
# labels over the dark (IRR < 0) cells, and their arrows, in white
LABEL_COLOR = {label: 'white' for label in
               ('isobutanol\ntiter', 'isobutanol\nyield', 'isobutanol productivity',
                'cell density', 'TCI')}

FIGURE = ksf.SweepFigure(
    output_stem='k13_inhib_isobutanol_opt_IRR_IRR',
    csv_prefix=f'ibo_{SWEEP_STEPS}_k_13_inhib_Spike_opt=False_max_n=50_',
    spec_1=SPEC_1, spec_2=SPEC_2,
    xlabel=r'$\bf{ALS\ capacity}$, $k_{\mathrm{13}}$' + f' [{ksf.G_PER_L_PER_H}]',
    ylabel=r'$\bf{Isobutanol\ inhibition}$ [× baseline]',
    optima=OPTIMA,
    xlim=(0., 16.), ylim=(0., 1.5),
    x_major=4., x_minor_div=4, y_major=0.5, y_minor_div=5,
    # co-located optima fanned apart [pt]: isobutanol titer and TCI share a cell
    marker_nudge={'isobutanol\ntiter': (0., 6.), 'TCI': (0., -6.)},
    label_color=LABEL_COLOR,
    # stacked labels: each arrow leaves its label's near edge, so it does not
    # cross the label below (left column from the left edge, right column
    # from the right edge)
    arrow_relpos={**{label: (0., 0.5) for label in
                     ('ethanol yield', 'ethanol titer', 'ethanol productivity')},
                  **{label: (1., 0.5) for label in
                     ('isobutanol\ntiter', 'isobutanol\nyield', 'TCI', 'cell density')}},
    color_axis='IRR',
    # IRR axis -15 to 35 % (the grid tops out at 31.1 %, the baseline 27.3 %)
    levels=np.arange(-15., 35.00001, 0.5),
    cbar_ticks=np.arange(-15., 35.00001, 5.),
    cbar_minor_step=1.,
    # opt_IRR baseline (k_13 = 4.0 g/L/h, multiplier 1x)
    baseline=(BASELINE_K13, 1.),
    baseline_callout=('baseline', (14, 22), -0.3),
)

if __name__ == '__main__':
    ksf.main(FIGURE)
