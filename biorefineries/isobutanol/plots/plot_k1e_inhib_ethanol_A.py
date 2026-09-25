#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""k_1e x ethanol-inhibition sweep figure for scenario A.

Minimum ethanol selling price (MESP = the sweep's purity-adjusted ethanol MPSP
at a 15 % IRR, converted from $/kg to $/GGE exactly as in
plots/plot_uncertainty_MPSP_vs_TCI.py) over the k_1e x inhib_ethanol-multiplier
grid of analyses/evaluate_EtOH_k1e_inhib_ethanol.py (scenario A, enzyme burden
ON, A's fed-batch feeding strategy, no feeding-strategy optimization), with the
grid optimum of cell density, ethanol titer, productivity, yield, TCI and MESP
marked and annotated. Format and style follow
plots/plot_feeding_strategy_baseline_A.py (panel A); the drawing lives in
plots/_kinetic_sweep_figure.py (loaded by file path).

Sim-safe: reads the sweep CSVs only and never imports the package.

Output: analyses/results/publication/Kinetic-sweeps/k1e_inhib_ethanol_A.{png,pdf}
"""

import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    '_kinetic_sweep_figure', os.path.join(HERE, '_kinetic_sweep_figure.py'))
ksf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ksf)

# Grid of the sweep to read (analyses/evaluate_EtOH_k1e_inhib_ethanol.py);
# SPEC_1 / SPEC_2 must match that script's linspaces: k_1e over the
# metabolic_14d glycolysis band 0.2x-4x of scenario A's k_1e (47.1 g/L/h, the
# nskinetics antimony value), the inhib_ethanol multiplier over the default
# group band 0.75x-1.5x.
SWEEP_STEPS = (40, 40, 1)
BASELINE_K1E = 47.1
SPEC_1 = np.linspace(0.2 * BASELINE_K1E, 4.0 * BASELINE_K1E, SWEEP_STEPS[0])  # k_1e [g/L/h]
SPEC_2 = np.linspace(0.75, 1.5, SWEEP_STEPS[1])  # inhib_ethanol multiplier [-]

# Optimum markers, as in the feeding-strategy figure. The sweep has no
# separate ethanol-yield column: 'Combined Yield' (ethanol + isobutanol per g
# sugars added) is the ethanol yield here, since scenario A makes no
# isobutanol. Offsets are set by eye for the current sweep; retune them if the
# optima move.
OPTIMA = [
    ('EtOH Titer',        'max', 'titer',        '^', 'white',   10, (-10, 22),  0.3),
    ('Cell loading',      'max', 'cell density', 'o', 'white',   10, (6, 32),   -0.3),
    ('EtOH Productivity', 'max', 'productivity', 's', 'white',    9, (0, 22),   -0.2),
    ('Combined Yield',    'max', 'yield',        'p', 'white',   10, (-14, 16),  0.3),
    ('TCI',               'min', 'TCI',          'p', '#33ccff', 10, (24, 6),    0.3),
    ('MPSP',              'min', 'MESP',         '*', '#33ccff', 14, (10, 34),  -0.3),
]

FIGURE = ksf.SweepFigure(
    output_stem='k1e_inhib_ethanol_A',
    csv_prefix=f'ibo_{SWEEP_STEPS}_k_1e_inhib_Spike_opt=False_max_n=16_',
    spec_1=SPEC_1, spec_2=SPEC_2,
    xlabel=r'$\bf{Glycolytic\ capacity}$, $k_{\mathrm{1e}}$' + f' [{ksf.G_PER_L_PER_H}]',
    ylabel=r'$\bf{Ethanol\ inhibition}$ [× baseline]',
    optima=OPTIMA,
    x_major=50., x_minor_div=2, y_major=0.25, y_minor_div=5,
    label_ha={'productivity': 'center'},
    # co-located optima fanned apart [pt]: cell density and TCI both sit at
    # the lowest-k_1e, lowest-multiplier corner
    marker_nudge={'cell density': (0., 5.), 'TCI': (5., 0.)},
)

if __name__ == '__main__':
    ksf.main(FIGURE)
