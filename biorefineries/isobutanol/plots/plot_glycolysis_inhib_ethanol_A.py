#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Glycolysis x ethanol-inhibition sweep figure for scenario A.

Minimum ethanol selling price (MESP = the sweep's purity-adjusted ethanol MPSP
at a 15 % IRR, converted from $/kg to $/GGE exactly as in
plots/plot_uncertainty_MPSP_vs_TCI.py) over the grouped glycolysis-capacity
multiplier (k_1l / k_1h / k_1e) x inhib_ethanol-multiplier grid of
analyses/evaluate_EtOH_glycolysis_inhib_ethanol.py (both axes on their wide
bands; scenario A, enzyme burden ON, A's fed-batch feeding strategy, no
feeding-strategy optimization), with the grid optimum of cell density,
ethanol titer, productivity, yield, TCI and MESP marked and annotated. Format
and style follow plots/plot_feeding_strategy_baseline_A.py (panel A); the
drawing lives in plots/_kinetic_sweep_figure.py (loaded by file path).

The lowest-glycolysis column (0.001x) makes essentially no ethanol (titer
~1e-5 g/L), so its MESP is undefined and it is left blank; the cell-density
and TCI optima nevertheless sit on it (no ethanol plant, all glucose to
growth).

Sim-safe: reads the sweep CSVs only and never imports the package.

Output: analyses/results/publication/Kinetic-sweeps/glycolysis_inhib_ethanol_A.{png,pdf}
"""

import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    '_kinetic_sweep_figure', os.path.join(HERE, '_kinetic_sweep_figure.py'))
ksf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ksf)

# Grid of the sweep to read (analyses/evaluate_EtOH_glycolysis_inhib_ethanol.py,
# the both-axes-wide run tagged _rb0.001-4_ib0.001-2); SPEC_1 / SPEC_2 must
# match that script's linspaces.
SWEEP_STEPS = (40, 40, 1)
SPEC_1 = np.linspace(1e-3, 4.0, SWEEP_STEPS[0])  # glycolysis multiplier [-]
SPEC_2 = np.linspace(1e-3, 2.0, SWEEP_STEPS[1])  # inhib_ethanol multiplier [-]

# Optimum markers, as in the feeding-strategy figure. 'Combined Yield' is the
# ethanol yield here (scenario A makes no isobutanol). Offsets are set by eye
# for the current sweep; retune them if the optima move.
OPTIMA = [
    ('EtOH Titer',        'max', 'titer',        '^', 'white',   10, (10, 20),   -0.3),
    ('Cell loading',      'max', 'cell density', 'o', 'white',   10, (28, 40),   -0.3),
    ('EtOH Productivity', 'max', 'productivity', 's', 'white',    9, (18, -4),    0.3),
    ('Combined Yield',    'max', 'yield',        'p', 'white',   10, (-14, 16),   0.3),
    ('TCI',               'min', 'TCI',          'p', '#33ccff', 10, (58, 14),   -0.2),
    ('MPSP',              'min', 'MESP',         '*', '#33ccff', 14, (12, 24),   -0.3),
]
# cell density and TCI sit on the blank no-ethanol column, so their labels
# are carried onto the orange plateau; productivity's label sits in the dark
# high-MESP wall, hence white
LABEL_COLOR = {'productivity': 'white'}

FIGURE = ksf.SweepFigure(
    output_stem='glycolysis_inhib_ethanol_A',
    csv_prefix=(f'ibo_{SWEEP_STEPS}_glyco_inhib_Spike_rb0.001-4_ib0.001-2'
                '_opt=False_max_n=16_'),
    spec_1=SPEC_1, spec_2=SPEC_2,
    xlabel=r'$\bf{Glycolysis\ capacity}$ [× baseline]',
    ylabel=r'$\bf{Ethanol\ inhibition}$ [× baseline]',
    optima=OPTIMA,
    xlim=(0., 4.), ylim=(0., 2.),
    x_major=1., x_minor_div=4, y_major=0.5, y_minor_div=5,
    label_color=LABEL_COLOR,
)

if __name__ == '__main__':
    ksf.main(FIGURE)
