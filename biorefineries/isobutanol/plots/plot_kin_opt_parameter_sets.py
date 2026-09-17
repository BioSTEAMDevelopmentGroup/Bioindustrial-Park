#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Publication figure comparing the 15 decision variables of the
`metabolic_minimal_subset` kinetic-optimization preset across the
scenario-A baseline and a handful of hand-picked campaign trials (each
trial optionally from a different campaign -- one trajectory CSV per
objective). Three stacked panels:

  a  outcome incumbent trajectories -- IRR, ethanol titer, isobutanol
     titer, batch time. One step line per set: the metric evaluated at
     that set's running optimization incumbent vs trial_number, where the
     incumbent is the arg-best of the set's own selection metric (its
     objective for a "best" trial, COL for "best:COL"). So the IRR panel
     shows best-so-far IRR for an IRR study and the at-best-titer IRR for
     a titer study. The scenario-A baseline (no study) is a dashed
     reference line.
  b  the 15 parameter values, grouped by nskinetics role and arranged in
     pathway order (glycolysis/fermentation r1->r3->r6, Ehrlich branch
     r13->r16, product-inhibition effector multipliers, feeding), one
     bar per set with the searched band shaded.
  c  full proteome allocation -- one horizontal stacked bar per set, each
     summing to the proteome cap eb.PROTEIN_CONTENT (0.49 g protein/gDCW).
     Sectors left to right: the fixed housekeeping block; the modeled
     metabolic pool Phi_M split into four fixed pathway categories
     (glycolysis r1; TCA cycle + acetate / acetyl-CoA production r2+r4+r5;
     ethanol production r3+r6; isobutanol production r13->r16); the
     unallocated flexible slack; and the growth-derated translation sector
     phi_T, flush against the cap on the right. Categories share the set's
     colour and are told apart by hatch. A single vertical dashed line marks
     the un-derated translation demand (its extent measured right-to-left
     from the cap).

"campaign" is this figure's word for a kinetic-optimization study; code,
CSV, and study names keep "study".

Sim-safe: kinetic_optimization.py (ko) and enzyme_burden.py (eb) are
loaded by file path (no biosteam import, no load()); campaign CSVs and
the parameter-distribution workbooks are plain pandas reads. Runnable
while a campaign is in flight. Run:

    python plots/plot_kin_opt_parameter_sets.py \
        --set "Best IRR" <campaign> best \
        --set "Best ethanol titer" <campaign> "best:EtOH titer" \
        --set "Best isobutanol titer" <campaign> "best:IBO titer"

With no --set arguments it plots the "best" trial of each of the five most
recent minimal-subset campaigns -- one per objective (IRR / ethanol titer /
isobutanol titer / ethanol yield / isobutanol yield) -- against the
baseline. Writes <stem>_<stamp>.png and .pdf to --out-dir.
"""
import os
import argparse
import importlib.util
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D, TICKDOWN, TICKLEFT
from matplotlib.ticker import (AutoMinorLocator, FixedLocator, FuncFormatter,
                               LogLocator, MultipleLocator, NullFormatter,
                               NullLocator, PercentFormatter)

# --- sim-safe module loads (by file path; never import the package) ----------
PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')
# The three most recent minimal-subset studies share one search space
# (rb0.001-10_ib0.2-2_burden) and differ only in the optimized objective,
# so each set is that study's own "best" trial (max of its objective).
_MINIMAL_SUBSET_STUDY = ('kin_opt_ethanol_isobutanol_metabolic_minimal_subset'
                         '_%s_rb0.001-10_ib0.2-2_burden')
DEFAULT_STUDY = _MINIMAL_SUBSET_STUDY % 'irr'          # financial (IRR) optimum
ETOH_TITER_STUDY = _MINIMAL_SUBSET_STUDY % 'etoh_titer'  # ethanol-titer optimum
IBO_TITER_STUDY = _MINIMAL_SUBSET_STUDY % 'ibo_titer'    # isobutanol-titer optimum
IBO_YIELD_STUDY = _MINIMAL_SUBSET_STUDY % 'ibo_yield'    # isobutanol-yield optimum
ETOH_YIELD_STUDY = _MINIMAL_SUBSET_STUDY % 'etoh_yield'  # ethanol-yield optimum


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(PKG_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ko = _load('ko', 'kinetic_optimization.py')
eb = _load('eb', 'enzyme_burden.py')

# --- the 15 decision variables, in figure order -----------------------------
RATE_VARS = list(ko.METABOLIC_MINIMAL_SUBSET_RATES)          # 9
# Ehrlich-branch capacities: genuinely zero (branch off) at the scenario-A
# baseline, so their baseline "0" bar labels are suppressed (see bar_cell)
EHRLICH_RATE_VARS = ('k_13', 'k_14', 'k_15', 'k_16')
GROUP_VARS = list(ko.METABOLIC_MINIMAL_SUBSET_GROUPS)        # 3
FEED_VARS = ['threshold_conc', 'target_delta', 'max_n_spikes']  # 3 CSV columns
DECISION_VARS = RATE_VARS + GROUP_VARS + FEED_VARS
# feeding cells DRAW derived / realized quantities in place of two raw decision
# columns: the applied target sugar concentration (target_conc = min(
# TARGET_CONC_MAX, threshold_conc + target_delta), via ko._applied_feeding) for
# target_delta, and the realized spike count (n_glu_spikes, a tracked metric)
# for the max_n_spikes cap.
FEED_DRAW_VARS = ['target_conc', 'threshold_conc', 'n_glu_spikes']

BANDS = [
    ('Glycolysis ($r_1$) + ethanol production ($r_3$ → $r_6$)',
     ['k_1l', 'k_1h', 'k_1e', 'k_3', 'k_6']),
    ('Isobutanol production (Ehrlich pathway, '
     '$r_{13}$ → $r_{14}$ → $r_{15}$ → $r_{16}$)',
     ['k_13', 'k_14', 'k_15', 'k_16']),
    ('Product inhibition relative to baseline '
     '(applied to $r_1$, $r_4$, $r_6$, $r_7$, $r_{10}$, $r_{16}$)',
     GROUP_VARS),
    ('Feeding strategy', FEED_DRAW_VARS),
]

# cell titles: parameter on the first line, reaction index + enzyme below
REACTION_LABELS = {
    'k_1l': 'k_1l\nr1 glycolysis\n(low affinity)',
    'k_1h': 'k_1h\nr1 glycolysis\n(high affinity)',
    'k_1e': 'k_1e\nr1 glycolysis\n(ethanol phase)',
    'k_3': 'k_3\nr3 PDC\n(Pdc1)',
    'k_6': 'k_6\nr6 ADH\n(Adh1)',
    'k_13': 'k_13\nr13 ALS\n(Ilv2+Ilv6)',
    'k_14': 'k_14\nr14 KARI\n(Ilv5)',
    'k_15': 'k_15\nr15 DHAD\n(Ilv3)',
    'k_16': 'k_16\nr16 KDC\n(Aro10)',
    'k_17': 'k_17\nr17 ADH\n(Adh6)',
}
GROUP_LABELS = {'inhib_ethanol': 'Ethanol\ninhibition',
                'inhib_isobutanol': 'Isobutanol\ninhibition',
                'inhib_acetate': 'Acetate\ninhibition'}
# value-axis units. Rate constants: every sampled capacity carries the
# Antimony unit g_per_l_per_h (gram/(litre*hour)) in the shipped model, so
# g·L^-1·h^-1 for all nine. Effector multipliers are dimensionless fold-changes
# relative to the per-family baseline (1 = baseline).
RATE_UNIT = 'g·L$^{-1}$·h$^{-1}$'
# per-rate-cell ylabel override (else the cell uses _mathify(var)). Empty for the
# minimal_subset layout; the split_12d layout gives ehrlich_downstream -- an
# absolute-capacity bar standing for three rates -- a reaction-range symbol
# instead of the mangled "$ehrlich_{downstream}$" _mathify would produce.
RATE_CELL_LABEL = {}
# derived "× baseline" fold-change cells: {display_var: raw_var}. A raw rate
# with a nonzero scenario-A baseline can be drawn as a fold-change group cell
# (value / baseline, baseline 1.0) instead of an absolute g/L/h bar. Empty for
# the minimal_subset layout; the split_12d layout registers k_3 / k_6 (the
# glycolysis-adjacent fermentation rates) and k_16 / k_17 (Aro10 / Adh6, whose
# fitted scenario-A rates are small but nonzero). load_set / baseline_set fill
# the display_var from the raw value divided by its scenario-A baseline.
REL_RATE_VARS = {}
# scenario-A "× baseline" divisors for REL_RATE_VARS raw rates that are NOT in
# the A workbook. k_3 / k_6 come from ko.workbook_kinetic_baselines('A'); k_16
# (Aro10) and k_17 (Adh6) are absent from the A workbook (it omits the Ehrlich
# block and Adh6), so their fitted-antimony scenario-A values are pinned by the
# split_12d layout (g/L/h). Empty for the minimal_subset layout.
REL_RATE_BASELINES_A = {}
# raw rates whose per-rate value lives in a CSV `applied_*` column rather than a
# decision column: metabolic_split_12d samples k_14/k_15/k_16 as the single
# ehrlich_downstream multiplier, so load_set copies applied_k_14/_k_15/_k_16 to
# plain k_14/k_15/k_16 keys. Empty for the minimal_subset layout.
APPLIED_RATE_COLS = {}
# feeding cell: (title, (shaded-range low, high)) -- engine default bounds
# target_conc shares the threshold cell's 0-300 g/L axis (it is clipped at
# TARGET_CONC_MAX = 300, and always sits at or above the threshold), so the two
# feeding concentrations read on the same scale.
FEED_LABELS = {'threshold_conc': ('Thresh. sugar\nconc. [g·L$^{-1}$]', (0, 300)),
               'target_conc': ('Target sugar\nconc. [g·L$^{-1}$]', (0, 300)),
               'n_glu_spikes': ('No. of spikes', (0, 50))}

# --- metabolic_split_12d layout (the 12-variable Ehrlich-split preset) --------
# The three newest campaigns (pi_log-tail / ibo_titer / ibo_yield, 2026-09-16)
# sample metabolic_split_12d: individual rates k_3, k_6, k_13, k_17; two
# capacity GROUPS -- glycolysis (0.2-4x of the live k_1l/k_1h/k_1e) and
# ehrlich_downstream (an ABSOLUTE g/L/h band on the anchor k_14, stoichiometric
# weights -> k_14/k_15/k_16); three inhibition-effector multipliers; and
# feeding. Panel B draws the underlying rates, NOT the raw ehrlich_downstream
# multiplier: k_13/k_14/k_15 (zero at the scenario-A baseline) as absolute g/L/h
# bars, and k_3/k_6/k_16/k_17 (nonzero baseline) plus glycolysis and the three
# inhibition families as × baseline fold-change bars (baseline 1.0). k_14/k_15/
# k_16 come from the CSV applied_* columns. Panels A (outcomes) and C (proteome
# allocation) are preset-agnostic. use_split_12d_layout() swaps the panel-B
# layout globals in place; the default minimal_subset layout is untouched.
SPLIT_12D_DECISION_VARS = ['k_3', 'k_6', 'k_13', 'k_17', 'glycolysis',
                           'ehrlich_downstream', 'inhib_ethanol',
                           'inhib_isobutanol', 'inhib_acetate',
                           'threshold_conc', 'target_delta', 'max_n_spikes']
# rate-kind cells (absolute g/L/h): the Ehrlich rates that are ZERO at the
# scenario-A baseline -- k_13 and the individually-plotted k_14 / k_15. The
# preset samples k_14/k_15/k_16 as the single ehrlich_downstream multiplier;
# their per-rate applied values live in the CSV applied_k_14 / _k_15 / _k_16
# columns, exposed under plain keys by load_set (SPLIT_12D_APPLIED_RATES).
SPLIT_12D_RATE_VARS = ['k_13', 'k_14', 'k_15']
# × baseline fold-change cells {display: raw}: the fermentation rates k_3 / k_6
# and the small-but-nonzero-baseline Ehrlich rates k_16 (Aro10) / k_17 (Adh6).
SPLIT_12D_REL_RATE_VARS = {'k_3_rel': 'k_3', 'k_6_rel': 'k_6',
                           'k_16_rel': 'k_16', 'k_17_rel': 'k_17'}
# scenario-A divisors for the raw rates absent from the A workbook: the fitted
# antimony values (g/L/h), confirmed in the nskinetics antimony 2026-09-16 --
# k_16 the Aro10 decarboxylase leak, k_17 the Adh6 rate (corrected 44 ->
# 0.1077). Re-pin if the antimony moves. (k_3 / k_6 come from the A workbook.)
SPLIT_12D_REL_RATE_BASELINES_A = {'k_16': 0.02115, 'k_17': 0.1077}
# k_14 / k_15 / k_16 are collapsed into the ehrlich_downstream decision column;
# their per-rate applied values are the applied_* CSV columns. load_set copies
# them to plain keys so the separated cells (and the k_16 fold-change) read them.
SPLIT_12D_APPLIED_RATES = {'k_14': 'applied_k_14', 'k_15': 'applied_k_15',
                           'k_16': 'applied_k_16'}
# glycolysis and the k_3/k_6/k_16/k_17 fold-change cells are group-kind too
# (multiplier bars, baseline 1.0), so they join the inhibition families in
# GROUP_VARS for the draw-kind dispatch; the inhibition BAND lists only the
# three multipliers.
SPLIT_12D_INHIB_GROUPS = ['inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate']
SPLIT_12D_GROUP_VARS = (['glycolysis', 'k_3_rel', 'k_6_rel',
                         'k_16_rel', 'k_17_rel'] + SPLIT_12D_INHIB_GROUPS)
# genuinely zero (branch off) at the scenario-A baseline -> suppress the "0"
# baseline bar label. k_3/k_6/k_16/k_17 are nonzero (drawn as × baseline group
# cells, baseline 1.0), so only k_13/k_14/k_15 need the suppression.
SPLIT_12D_EHRLICH_RATE_VARS = ('k_13', 'k_14', 'k_15')
SPLIT_12D_BANDS = [
    ('Glycolysis capacity ($r_1$) + ethanol production ($r_3$, $r_6$)',
     ['glycolysis', 'k_3_rel', 'k_6_rel']),
    ('Isobutanol production (Ehrlich pathway $r_{13}$–$r_{16}$; Adh6 $r_{17}$)',
     ['k_13', 'k_14', 'k_15', 'k_16_rel', 'k_17_rel']),
    ('Product inhibition relative to baseline '
     '(applied to $r_1$, $r_4$, $r_6$, $r_7$, $r_{10}$, $r_{17}$)',
     SPLIT_12D_INHIB_GROUPS),
    ('Feeding strategy', FEED_DRAW_VARS),
]


def use_split_12d_layout():
    """Reassign the panel-B layout globals to the metabolic_split_12d preset.
    Idempotent; called once from main() when the plotted campaigns are 12d.
    Panels A and C are preset-agnostic, so only the panel-B constants change."""
    global DECISION_VARS, RATE_VARS, GROUP_VARS, EHRLICH_RATE_VARS, BANDS
    global REL_RATE_VARS, REL_RATE_BASELINES_A, APPLIED_RATE_COLS
    DECISION_VARS = SPLIT_12D_DECISION_VARS
    RATE_VARS = SPLIT_12D_RATE_VARS
    GROUP_VARS = SPLIT_12D_GROUP_VARS
    EHRLICH_RATE_VARS = SPLIT_12D_EHRLICH_RATE_VARS
    BANDS = SPLIT_12D_BANDS
    REL_RATE_VARS = dict(SPLIT_12D_REL_RATE_VARS)
    REL_RATE_BASELINES_A = dict(SPLIT_12D_REL_RATE_BASELINES_A)
    APPLIED_RATE_COLS = dict(SPLIT_12D_APPLIED_RATES)
    GROUP_LABELS['glycolysis'] = '$k_{1e}$, $k_{1h}$, $k_{1l}$\n(× baseline)'
    GROUP_LABELS['k_3_rel'] = '$k_3$\n(× baseline)'
    GROUP_LABELS['k_6_rel'] = '$k_6$\n(× baseline)'
    GROUP_LABELS['k_16_rel'] = '$k_{16}$\n(× baseline)'
    GROUP_LABELS['k_17_rel'] = '$k_{17}$\n(× baseline)'
    # the minimal_subset figure hard-capped IRR at 25 % (its incumbents topped
    # out at 23.5 %); the 12d PI campaign reaches ~27 %, so let the IRR cell
    # auto-scale (to 30 % at the 0.05 step) instead of clipping the line.
    OUTCOME_TICK_CAP.pop('IRR', None)

# the seven study steps -> enzyme name and charging parameter(s); read
# against the eb tables so a table drift here raises at import
STEP_ENZYME = {'r1': 'Glycolysis lump', 'r3': 'Pdc1', 'r6': 'Adh1',
               'r13': 'Ilv2+Ilv6', 'r14': 'Ilv5', 'r15': 'Ilv3',
               'r16': 'Aro10', 'r17': 'Adh6'}   # r17 (Adh6) is the native step of the 2026-09-15 split; sampled by metabolic_split_12d/14d campaigns
STEP_PARAMS = {'r1': 'k_1l, k_1h, k_1e', 'r3': 'k_3', 'r6': 'k_6',
               'r13': 'k_13', 'r14': 'k_14', 'r15': 'k_15', 'r16': 'k_16',
               'r17': 'k_17'}
for _s in STEP_ENZYME:
    if _s not in eb.NATIVE_STEPS and _s not in eb.EHRLICH_STEPS:
        raise KeyError(f'STEP_ENZYME step {_s!r} not in eb tables')

# outcomes panel: (CSV column, cell title, (y-low, y-high)), in draw order.
# Titles are rotated y-axis labels in narrow cells -- keep each to <=2 lines
# so the label stays out of the neighbouring cell's plot box. Order:
# IRR (the big left cell), then the isobutanol metrics (yield, titer,
# productivity) down one column, then the ethanol metrics down the next.
OUTCOMES = (('IRR', 'Financial attractiveness\nas IRR [%]', (0, 0.3)),
            ('IBO yield', 'Isobutanol yield\n[g·g$^{-1}$]', (0, 0.4)),
            ('IBO titer', 'Isobutanol titer\n[g·L$^{-1}$]', (0, 60)),
            ('IBO productivity',
             'Isobutanol productivity\n[g·L$^{-1}$·h$^{-1}$]', (0, 2)),
            ('EtOH yield', 'Ethanol yield\n[g·g$^{-1}$]', (0, 0.5)),
            ('EtOH titer', 'Ethanol titer\n[g·L$^{-1}$]', (0, 200)),
            ('EtOH productivity',
             'Ethanol productivity\n[g·L$^{-1}$·h$^{-1}$]', (0, 8)))

# per-outcome override for the major-tick step (else _linear_cap's step).
# IRR is stored as a fraction shown in percent, so 0.05 -> 0/5/.../30 %;
# isobutanol titer reads cleaner on 0/20/40/60 than the default 0/15/30/45/60.
OUTCOME_TICK_STEP = {'IRR': 0.05, 'IBO titer': 20.0}

# per-outcome value-axis ceiling override (fraction/native units), used in place
# of the auto nice-ceiling. IRR tops out at 25 % rather than the auto 30 % so the
# top major tick is 25; the best incumbent (23.5 %) still fits below it.
OUTCOME_TICK_CAP = {'IRR': 0.25}

# outcomes whose negative values are drawn at zero: a loss-making IRR (finite
# negative) and an unsolvable IRR (-inf) both read as 0 rather than diving
# below the axis floor / being omitted. Only nan (never solved) stays
# non-finite and is omitted.
CLAMP_NEG_TO_ZERO = {'IRR'}

# which campaign "owns" each outcome panel -- its incumbent line is drawn thick
# and its per-trial cloud appears there. By default a panel is owned by the
# campaign whose objective slug equals the panel's own column. The financial-
# attractiveness panel is plotted on the IRR axis but its owning campaign may
# optimize a profitability-index objective instead (the split_12d campaigns
# maximize PI / PI (log-tail), whose incumbent is shown against the same IRR
# axis), so that panel is owned by any financial-attractiveness objective.
OUTCOME_OWN_OBJECTIVES = {'IRR': frozenset({'IRR', 'PI', 'PI (log-tail)'})}


def _owns_outcome(objective, col):
    """True if a campaign with this objective slug owns the `col` outcome
    panel (thick incumbent line + trial cloud)."""
    return objective in OUTCOME_OWN_OBJECTIVES.get(col, frozenset({col}))

# one color per set: baseline grey, campaigns from the hue palette. Fixed
# 7-color palette (blue, orange, green, purple, red, grey, yellow); the baseline
# takes the grey and the five campaigns take blue/orange/green/purple/yellow in
# order (red is the spare sixth hue).
BASELINE_COLOR = '#90918e'
HUE_COLORS = ['#18C4DC', '#f98f60', '#79bf82', '#a280b9', '#f3c354', '#ED586F']
MAX_SETS = 1 + len(HUE_COLORS)   # 7

FONTS = {'band': 12, 'cell': 10, 'tick': 9, 'callout': 9,
         'legend': 10, 'axis': 11, 'panel': 14}

# light-grey backing for panel b, so it reads as a distinct block from the
# white panels a and c
PANEL_B_BG = '0.93'

# --- scenario-A baseline outcomes: HARD-CODED (spec decision b).
# Re-simulated 2026-09-13 (smoke_test_1 protocol, IBO_2026: load(),
# scenarios.load_scenario('A'), model_specification, solve_TEA; ethanol MPSP
# 0.86867) on the current nskinetics dev (anaerobic_growth_mult 0.75, r13 /
# r14-r16 rate laws -- the latter inert for A). A cached one-off baseline
# simulation (option a) is a later change; when done, replace this dict with
# the cached read.
# History: fermentation entries first recorded 2026-09-07 (EtOH titer 118.4,
# tau 45.4, 10 spikes, EtOH yield 0.455, cell density 15.7) under
# anaerobic_growth_mult 1.0; TEA entries (IRR, TCI) re-recorded 2026-09-12 for
# the DDGS-dryer acid split (the 2026-09-07 values were 0.1230 / 139.6) and
# unchanged on 2026-09-13.
BASELINE_A = {
    'IRR': 0.1260, 'EtOH titer': 114.5, 'IBO titer': 0.0, 'tau': 61.6,
    'n_glu_spikes': 7, 'EtOH yield': 0.460, 'IBO yield': 0.0,
    'Cell density': 13.0, 'TCI': 139.3, 'threshold_conc': 217.125,
    'target_delta': 4.125, 'max_n_spikes': 16,
}


# --- data assembly -----------------------------------------------------------
def baseline_set():
    """Scenario-A baseline as a flat record shaped like a campaign row.

    Kinetics: the scenario-A workbook rows (the four Ehrlich rates are
    absent from the A workbook == zero in the live A model); effector
    multipliers 1.0; feeding from BASELINE_A. Burden pools from
    eb.BurdenModel at that point (r13-r16 == 0). Outcomes are the
    hard-coded BASELINE_A constants (spec decision b).
    """
    A = ko.workbook_kinetic_baselines('A')
    # k_ref for the burden model: every required capacity, Ehrlich at 0
    k_ref = dict(A)
    for k in ('k_13', 'k_14', 'k_15', 'k_16'):
        k_ref.setdefault(k, 0.0)
    # r17 (Adh6) is a NATIVE constitutive step, not an Ehrlich one: it is
    # present in scenario A at its live rate, so it does NOT default to 0. The
    # A workbook has no k_17 row; the live scenario-A value is the antimony
    # default, corrected 44 -> 0.1077 g/L/h on 2026-09-16 (CLAUDE.md k_17,ref
    # note; confirmed live off r_te 2026-09-16, matching the enzyme-burden
    # reference). Its native pool is self-referential (multiplier
    # k_17/reference == 1), so pool_r17 equals the wild-type value for any
    # nonzero k_17 -- panel C is unaffected; only the panel-B baseline k_17 bar
    # reads this value. NB the split_12d/14d campaigns' k_17 SEARCH band anchors
    # on the B workbook's still-baked 44.0, so a plotted campaign's k_17 may sit
    # far above this corrected baseline. Re-pin if the antimony k_17 moves.
    k_ref.setdefault('k_17', 0.1077)
    missing = [c for c in eb.BurdenModel.required_capacities()
               if c not in k_ref]
    if missing:
        raise KeyError('scenario-A workbook is missing burden capacities '
                       f'{missing}')
    res = eb.BurdenModel(k_ref).evaluate(k_ref)

    rec = {'label': 'Baseline (no optimization)', 'campaign': None,
           'trial_number': None, 'is_baseline': True, 'extra_sampled': []}
    # rate-kind cells read their live scenario-A value; a rate absent from the
    # workbook is 0 in the live A model (the Ehrlich rates k_13..k_16, and the
    # split_12d ehrlich_downstream capacity, which is off at baseline). Group
    # multipliers (inhibition families, and split_12d glycolysis) are 1.0.
    for k in RATE_VARS:
        rec[k] = float(k_ref.get(k, 0.0))
    for g in GROUP_VARS:
        rec[g] = 1.0
    rec['threshold_conc'] = BASELINE_A['threshold_conc']
    rec['target_delta'] = BASELINE_A['target_delta']
    rec['max_n_spikes'] = BASELINE_A['max_n_spikes']
    # applied target sugar concentration (clip at TARGET_CONC_MAX), drawn
    # in place of target_delta
    rec['target_conc'] = float(ko._applied_feeding(rec)[0])
    for col in ('IRR', 'TCI', 'EtOH titer', 'IBO titer', 'EtOH yield',
                'IBO yield', 'tau', 'n_glu_spikes'):
        rec[col] = BASELINE_A[col]
    # productivity = water-basis titer / batch time (the registry definition);
    # baseline IBO titer is 0 so IBO productivity is 0
    rec['IBO productivity'] = BASELINE_A['IBO titer'] / BASELINE_A['tau']
    rec['EtOH productivity'] = BASELINE_A['EtOH titer'] / BASELINE_A['tau']
    for st, pool in res.pools.items():
        rec[f'pool_{st}'] = float(pool)
    rec['Phi_M'] = float(res.Phi_M)
    rec['phi_T'] = float(res.phi_T)
    rec['F_flex'] = float(res.F_flex)
    rec['burden_factor'] = float(res.burden_factor)
    # no study -> no incumbent trajectory; drawn as a dashed reference line
    rec['sel_col'] = None
    rec['sel_dir'] = None
    rec['traj_x'] = None
    rec['traj'] = None
    return rec


# objective slugs, longest first so "IBO yield x titer" beats "IBO yield"
_OBJ_SLUGS = sorted(ko.OBJECTIVE_REGISTRY, key=len, reverse=True)


def resolve_campaign_csv(campaign):
    """A study name -> RESULTS_DIR/<name>_trajectory.csv, or a literal
    path to a trajectory CSV used as-is."""
    if campaign.lower().endswith('.csv'):
        path = campaign
    else:
        path = os.path.join(RESULTS_DIR, f'{campaign}_trajectory.csv')
    if not os.path.isfile(path):
        raise FileNotFoundError(f'campaign trajectory CSV not found: {path}')
    return path


def campaign_objective(campaign):
    """The objective slug embedded in a study name (the token the driver
    slugged in), or None. Matches an OBJECTIVE_REGISTRY key with spaces
    turned to underscores, e.g. ..._minimal_subset_irr_... -> 'IRR'."""
    stem = os.path.basename(campaign)
    if stem.endswith('.csv'):
        return None
    low = stem.lower()
    for obj in _OBJ_SLUGS:
        if f'_{obj.lower().replace(" ", "_")}_' in low or \
           low.endswith('_' + obj.lower().replace(' ', '_')):
            return obj
    return None


def selection_spec(df, trial, campaign, warn=False):
    """(selection column, direction) defining the running incumbent for a
    set -- the same rule resolve_trial uses to pick its representative row:
    "best:COL" -> (COL, 'maximize'); "best" and an explicit integer trial
    (a placeholder while its study runs) -> the 'objective' column with the
    campaign objective's registry direction ('maximize' if the slug is
    unknown). Factored out so the representative row (resolve_trial) and the
    panel-a incumbent trajectory (incumbent_trajectory) cannot drift."""
    if isinstance(trial, str) and trial.startswith('best') and ':' in trial:
        col = trial.split(':', 1)[1]
        if col not in df.columns:
            raise ValueError(
                f'campaign {campaign}: no metric column {col!r} for '
                f'"best:{col}"')
        return col, 'maximize'
    obj = campaign_objective(campaign)
    if obj is None or obj not in ko.OBJECTIVE_REGISTRY:
        if warn:
            print(f'  WARNING campaign {campaign}: objective slug not in the '
                  'registry; "best" maximizes the objective column')
        return 'objective', 'maximize'
    return 'objective', ko.OBJECTIVE_REGISTRY[obj]['direction']


def resolve_trial(df, trial, campaign):
    """Return the requested COMPLETE row as a Series."""
    ok = df[df['state'] == 'COMPLETE']
    if isinstance(trial, str) and trial.startswith('best'):
        sel_col, direction = selection_spec(df, trial, campaign, warn=True)
        idx = ok[sel_col].idxmin() if direction == 'minimize' \
            else ok[sel_col].idxmax()
        return ok.loc[idx]
    # explicit integer trial_number
    n = int(trial)
    hit = df[df['trial_number'] == n]
    if hit.empty:
        raise ValueError(f'campaign {campaign}: trial {n} not found')
    row = hit.iloc[0]
    if row['state'] != 'COMPLETE':
        raise ValueError(f'campaign {campaign}: trial {n} is '
                         f'{row["state"]}, not COMPLETE (no outcomes/pools)')
    return row


def incumbent_trajectory(df, sel_col, direction):
    """Running-incumbent trajectory over a study's COMPLETE trials.

    Walking COMPLETE trials in trial_number order, keep the arg-best of
    `sel_col` under `direction` (non-finite selection values never win;
    before the first finite one the incumbent is undefined). Return
    (x, {outcome_col: y}): x is the trial_number of every COMPLETE trial
    (plus the study's final trial_number so the step reaches the right
    edge), and each y is that outcome column evaluated at the running
    incumbent -- non-finite -> NaN, so the line breaks rather than diving
    to a floor (e.g. an unsolvable-IRR trial that is a titer incumbent)."""
    ok = df[df['state'] == 'COMPLETE'].sort_values('trial_number')
    n = len(ok)
    x = ok['trial_number'].to_numpy(dtype=float)
    sel = ok[sel_col].to_numpy(dtype=float)
    inc = np.full(n, -1, dtype=int)
    best_pos, best_val = -1, None
    for i in range(n):
        v = sel[i]
        if np.isfinite(v) and (best_pos < 0 or
                               (v > best_val if direction == 'maximize'
                                else v < best_val)):
            best_pos, best_val = i, v
        inc[i] = best_pos
    mask = inc >= 0
    traj = {}
    for col, _, _ in OUTCOMES:
        vals = ok[col].to_numpy(dtype=float)
        y = np.full(n, np.nan)
        y[mask] = vals[inc[mask]]
        y[~np.isfinite(y)] = np.nan
        traj[col] = y
    if n:
        x_end = float(df['trial_number'].max())
        if x_end > x[-1]:
            x = np.append(x, x_end)
            for col in traj:
                traj[col] = np.append(traj[col], traj[col][-1])
    return x, traj


def load_set(label, campaign, trial):
    """A campaign trial as a flat record (CSV row + metadata + the panel-a
    incumbent trajectory)."""
    df = ko.load_trajectory(resolve_campaign_csv(campaign))
    missing = [c for c in DECISION_VARS if c not in df.columns]
    if missing:
        raise ValueError(
            f'campaign {campaign}: not a metabolic_minimal_subset campaign '
            f'(missing decision columns {missing}). A metabolic / '
            'metabolic_protein campaign samples the inhibition coefficients '
            'individually and has no inhib_* group columns.')
    missing_b = [c for c in eb.BURDEN_COLUMNS if c not in df.columns]
    if missing_b:
        raise ValueError(f'campaign {campaign}: missing burden columns '
                         f'{missing_b}')
    row = resolve_trial(df, trial, campaign)
    rec = row.to_dict()
    rec['label'] = label
    rec['campaign'] = campaign
    rec['trial_number'] = int(row['trial_number'])
    rec['is_baseline'] = False
    # applied target sugar concentration (clip at TARGET_CONC_MAX), drawn
    # in place of the raw target_delta decision column
    rec['target_conc'] = float(ko._applied_feeding(rec)[0])
    # metabolic_split_12d collapses k_14/k_15/k_16 into the ehrlich_downstream
    # multiplier; expose each rate's applied value under a plain key so the
    # separated panel-B cells (and the k_16 fold-change) read them uniformly.
    for plain, applied in APPLIED_RATE_COLS.items():
        if applied in rec:
            rec[plain] = float(rec[applied])
    # × baseline fold-change cells: raw value / scenario-A baseline. k_3/k_6 use
    # the A workbook; k_16/k_17 use the pinned antimony divisors.
    if REL_RATE_VARS:
        A_ref = dict(ko.workbook_kinetic_baselines('A'))
        A_ref.update(REL_RATE_BASELINES_A)
        for disp, raw in REL_RATE_VARS.items():
            b = A_ref.get(raw)
            rec[disp] = float(rec[raw]) / b if b else np.nan
    # decision columns beyond the 15 (e.g. a metabolic_minimal campaign's
    # extra rates + stage_1_max_x); reported, not drawn
    known = set(DECISION_VARS) | {'trial_number', 'state', 'objective',
                                  'error'}
    known |= set(eb.BURDEN_COLUMNS)
    known |= set(ko.TRACKED_METRICS)
    known |= {c for c in df.columns if c.startswith('applied_')}
    rec['extra_sampled'] = [c for c in df.columns if c not in known]
    sel_col, sel_dir = selection_spec(df, trial, campaign)
    tx, traj = incumbent_trajectory(df, sel_col, sel_dir)
    rec['sel_col'] = sel_col
    rec['sel_dir'] = sel_dir
    rec['traj_x'] = tx
    rec['traj'] = traj
    # full per-trial cloud for panel a, drawn only on the cell whose metric
    # this campaign optimized: every trial's number, its completion flag and
    # each outcome value (non-COMPLETE / unsolved -> NaN via the CSV).
    rec['objective'] = campaign_objective(campaign)
    rec['scatter_x'] = df['trial_number'].to_numpy(dtype=float)
    rec['scatter_complete'] = (df['state'] == 'COMPLETE').to_numpy()
    rec['scatter'] = {c: df[c].to_numpy(dtype=float) for c, _, _ in OUTCOMES}
    return rec


# --- searched-band shading, from the driver's own preset -> search-space
# construction -----------------------------------------------------------
_DEFAULT_AXES = ('ethanol_isobutanol', 'metabolic_minimal_subset')
_TARGET_PRODUCTS = ('ethanol_isobutanol', 'ethanol_only')
_STUDY_TYPES = ('metabolic_minimal_subset', 'metabolic_minimal',
                'metabolic_split_12d', 'metabolic_split_14d',
                'metabolic_protein', 'metabolic')   # longest-first


def campaign_axes_from_name(campaign):
    """(study_target_products, study_type) parsed from a study name;
    default to the minimal-subset ethanol+isobutanol preset (with a
    warning) for a bare CSV path or an unrecognised name."""
    stem = os.path.basename(campaign)
    if stem.endswith('.csv') or 'kin_opt_' not in stem:
        print(f'  WARNING campaign {campaign}: band source not inferrable '
              'from the name; assuming the metabolic_minimal_subset preset')
        return _DEFAULT_AXES
    tp = next((t for t in _TARGET_PRODUCTS if t in stem), _DEFAULT_AXES[0])
    ty = next((t for t in _STUDY_TYPES if t in stem), _DEFAULT_AXES[1])
    return tp, ty


def _split_12d_band(campaign):
    """{var: (lo, hi)} searched band for the 12 metabolic_split_12d decision
    variables. Cosmetic: bar_cell does not shade the band, so this is used only
    for the cross-campaign equality check in build_sets. Built deterministically
    from the preset so the three 12d campaigns share one band (no spurious
    warning). The individual rates and the glycolysis multiplier take their
    multiplier bands; ehrlich_downstream is the absolute anchor band."""
    tp, ty = campaign_axes_from_name(campaign)
    preset = ko.resolve_study_preset(tp, ty)
    rb = tuple(preset['rate_multiplier_bounds'])
    band = {v: rb for v in ('k_3', 'k_6', 'k_13', 'k_17')}
    band['ehrlich_downstream'] = tuple(ko.IBO_PATHWAY_ZERO_A_RATE_BOUNDS)
    gmb = preset['group_multiplier_bounds']
    for g in ('glycolysis', 'inhib_ethanol', 'inhib_isobutanol',
              'inhib_acetate'):
        try:
            band[g] = ko.group_bounds_for(g, gmb)
        except Exception:
            band[g] = (0.2, 4.0) if g == 'glycolysis' else (0.75, 1.5)
    band['threshold_conc'] = (0.0, 300.0)
    band['target_delta'] = (5.0, 500.0)
    band['max_n_spikes'] = (0, 50)
    return band


def campaign_band(campaign):
    """{var: (lo, hi)} searched band for the decision vars, exactly as
    the driver built it: ko.resolve_study_preset -> ko.build_search_space
    on the scenario-A baselines. Feeding vars from the engine defaults.
    A metabolic_split_12d campaign takes the dedicated _split_12d_band (its
    ehrlich_downstream referenced capacity group is built differently).

    ko.build_search_space returns (space, excluded_parameter_names), and
    it needs a `kinetic_baselines` dict that already carries a value for
    EVERY name it should place in the space -- it iterates
    kinetic_baselines.items(), never preset['include_params'] on its own
    -- and, for a `parameter_groups` member, a value that is POSITIVE
    (the group-forming loop checks kinetic_baselines[member] directly,
    before param_bounds_override is even consulted). The scenario-A
    workbook (ko.workbook_kinetic_baselines('A')) alone is not enough:
    - the Ehrlich-branch rate constants (k_13..k_16, role `capacity`)
      are absent from it because they are genuinely ZERO in the live
      scenario-A model (the branch is off); a zero/absent baseline with
      no override is EXCLUDED by build_search_space, so the driver's own
      real run (which discovers every parameter off the LIVE model,
      still 0 there) supplies a `param_bounds_override` derived from the
      scenario-B workbook instead (`kinetic_bounds_from_scenario` in
      analyses/optimize_kinetics_BO.py, ko.workbook_kinetic_bounds here)
      so those rates still get a positive, B-baseline-derived band.
    - the inhibition-coefficient group members (k_1ie, k_4ie, ...) are
      likewise absent from the A workbook, but -- unlike the Ehrlich
      rates -- they are NOT zero in the live A model: A's workbook simply
      doesn't curate them, so they carry the shipped Antimony default,
      which is numerically the scenario-B workbook value (see
      plot_kin_opt_best_vs_baseline_heatmap.py's `build_rows`, same
      fallback). A group member needs that positive value in
      `kinetic_baselines` itself, not just an override.
    Every name shared by both workbooks (k_1l, k_1h, k_1e, k_3, k_6) has
    the identical baseline in each, so folding in the B-derived override
    changes nothing for those -- the shared-baseline rate band is exact.
    """
    tp, ty = campaign_axes_from_name(campaign)
    if ty == 'metabolic_split_12d':
        return _split_12d_band(campaign)
    preset = ko.resolve_study_preset(tp, ty)
    kb_scenario = preset['kinetic_bounds_scenario']
    A = ko.workbook_kinetic_baselines('A')
    A_bounds_wb = ko.workbook_kinetic_baselines(kb_scenario)
    kinetic_baselines = dict(A)
    group_members = {m for members in (preset['parameter_groups'] or {}).values()
                     for m in members}
    for name in group_members:
        if name not in kinetic_baselines:
            # inherited (nonzero) Antimony-default value; see docstring
            kinetic_baselines[name] = A_bounds_wb[name]
    for name in preset['include_params']:
        if name not in kinetic_baselines:
            # genuinely zero in the live scenario-A model (e.g. k_13..k_16)
            kinetic_baselines[name] = 0.0
    # the scenario-B-derived absolute bounds that let a zero-baseline
    # Ehrlich rate stay in the space instead of being excluded; identical
    # in value to the A-baseline band for every name the two workbooks share
    param_bounds_override = ko.workbook_kinetic_bounds(
        kb_scenario,
        multiplier_bounds=preset['multiplier_bounds'],
        rate_multiplier_bounds=preset['rate_multiplier_bounds'],
        rate_params=preset['rate_params'],
        parameter_multiplier_bounds=preset['parameter_multiplier_bounds'],
    )
    space, _excluded = ko.build_search_space(
        kinetic_baselines,
        multiplier_bounds=preset['multiplier_bounds'],
        param_bounds_override=param_bounds_override,
        include_params=preset['include_params'],
        exclude_params=preset['exclude_params'],
        rate_multiplier_bounds=preset['rate_multiplier_bounds'],
        rate_params=preset['rate_params'],
        parameter_multiplier_bounds=preset['parameter_multiplier_bounds'],
        parameter_groups=preset['parameter_groups'],
        group_multiplier_bounds=preset['group_multiplier_bounds'],
        spike_delta_bounds=preset['spike_delta_bounds'],
        stage_1_max_x_bounds=preset['stage_1_max_x_bounds'],
    )
    band = {}
    for v in DECISION_VARS:
        if v in space:
            band[v] = (space[v]['low'], space[v]['high'])
        elif v in GROUP_VARS:
            # fallback, in case a future preset groups differently and the
            # group itself doesn't end up as its own entry in `space`
            # (group_bounds_for resolves the per-group dict / shared tuple and
            # falls back to the default band for a group absent from the dict)
            band[v] = ko.group_bounds_for(v, preset['group_multiplier_bounds'])
    # feeding fallbacks (build_search_space always emits these, but be safe)
    band.setdefault('threshold_conc', (0.0, 300.0))
    band.setdefault('target_delta', (5.0, 500.0))
    band.setdefault('max_n_spikes', (0, 50))
    return band


# --- drawing -----------------------------------------------------------------
def apply_fonts():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['xtick.labelsize'] = FONTS['tick']
    plt.rcParams['ytick.labelsize'] = FONTS['tick']
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['hatch.linewidth'] = 0.6
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['mathtext.bf'] = 'Arial:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


def set_colors(sets):
    colors, hue = {}, 0
    for s in sets:
        if s.get('is_baseline'):
            colors[id(s)] = BASELINE_COLOR
        else:
            if hue >= len(HUE_COLORS):
                n_campaign = len([x for x in sets if not x.get('is_baseline')])
                raise ValueError(
                    f'too many campaign sets ({n_campaign}); at most '
                    f'{len(HUE_COLORS)} plus the baseline (six total)')
            colors[id(s)] = HUE_COLORS[hue]
            hue += 1
    return colors


def _inward_top_right_ticks(ax, do_x=True, do_y=True):
    """Retarget the top/right tick lines to inward-only markers (bottom/left
    keep their in+out markers). matplotlib's tick `direction` is per-axis, so
    the secondary side (tick2: top for x, right for y) is given an inward
    marker whose reach (half the tick length) matches the inner half of the
    in+out ticks."""
    if do_x:
        for ticks, half in ((ax.xaxis.get_major_ticks(), 2.0),
                            (ax.xaxis.get_minor_ticks(), 1.1)):
            for t in ticks:
                t.tick2line.set_marker(TICKDOWN)
                t.tick2line.set_markersize(half)
    if do_y:
        for ticks, half in ((ax.yaxis.get_major_ticks(), 2.0),
                            (ax.yaxis.get_minor_ticks(), 1.1)):
            for t in ticks:
                t.tick2line.set_marker(TICKLEFT)
                t.tick2line.set_markersize(half)


def style_cell_axes(ax):
    # value (y) axis: ticks in+out on the left, mirrored inward-only on the
    # right; the box is closed on all four sides, matching panel A.
    ax.tick_params(axis='y', which='major', direction='inout',
                   left=True, right=True, labelright=False, length=4)
    ax.tick_params(axis='y', which='minor', direction='inout',
                   left=True, right=True, length=2.2)
    # categorical x axis (one bar per campaign, keyed by the legend): no x
    # ticks or labels on either the bottom or the mirrored top axis.
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis='x', which='both', top=False, bottom=False,
                   labelbottom=False)
    for sp in ('right', 'top'):
        ax.spines[sp].set_visible(True)
    _inward_top_right_ticks(ax, do_x=False, do_y=True)


def _baseline_value(sets, var):
    for s in sets:
        if s.get('is_baseline'):
            return s.get(var)
    return None


# --- nice, tick-aligned axis caps: every quantitative axis is bounded on
# both ends by a labeled tick (like the IRR cell's 0 -> 30) ------------------
_NICE_MANTISSAS = (1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)


def _nice_ceiling(x):
    """Smallest 'nice' number (mantissa in _NICE_MANTISSAS) >= x."""
    if x <= 0:
        return 1.0
    e = np.floor(np.log10(x))
    scale = 10.0 ** e
    m = x / scale
    for cand in _NICE_MANTISSAS:
        if m <= cand * (1 + 1e-9):
            return cand * scale
    return 10.0 * scale


def _linear_cap(vmax, floor_hi=0.0):
    """(hi, step) for a linear axis based at 0 with both 0 and hi as major
    ticks. hi is the nice ceiling of max(vmax, floor_hi); step = hi / k for
    the fewest-tick k in (3, 4, 5, 6, 8) that lands on a clean 1/2/2.5/5
    mantissa (falling back to hi / 5). hi is an integer multiple of step, so
    a MultipleLocator(step) labels both ends."""
    hi = _nice_ceiling(max(vmax, floor_hi))
    for k in (3, 4, 5, 6, 8):
        step = hi / k
        e = np.floor(np.log10(step))
        m = step / 10.0 ** e
        if min(abs(m - t) for t in (1.0, 2.0, 2.5, 5.0)) < 0.03:
            return hi, step
    return hi, hi / 5.0


def _log_stop_below(x):
    """Largest 1/2/5-decade stop strictly below x (e.g. 0.2 -> 0.1)."""
    e = np.floor(np.log10(x))
    for m in (5.0, 2.0, 1.0):
        v = m * 10.0 ** e
        if v < x * (1 - 1e-9):
            return v
    return 10.0 ** (e - 1) * 5.0


def _log_stop_above(x):
    """Smallest 1/2/5-decade stop strictly above x (e.g. 2 -> 5)."""
    e = np.floor(np.log10(x))
    for m in (1.0, 2.0, 5.0):
        v = m * 10.0 ** e
        if v > x * (1 + 1e-9):
            return v
    return 10.0 ** (e + 1)


def _decade_floor(lo):
    """Enclosing decade at or below lo, dropped one more when the band would
    otherwise touch the bottom edge (keeps ~>=0.3 decade headroom for the
    '0'/absent bar labels); the result is a decade, so it carries a tick."""
    a = np.floor(np.log10(lo))
    if 10.0 ** a > lo * 0.5:
        a -= 1
    return 10.0 ** a


def _decade_ceil(hi):
    """Enclosing decade at or above hi, raised one more when the band would
    otherwise touch the top edge; the result is a decade (a tick)."""
    b = np.ceil(np.log10(hi))
    if 10.0 ** b < hi * 2:
        b += 1
    return 10.0 ** b


# hand-tuned major y-ticks for a few rate cells (readability): an explicit
# tick list overrides the auto _linear_cap ticks and sets the axis top to its
# largest tick. A uniform list keeps half-step minor ticks (so e.g. k_14's
# dropped 0.5/1.5 survive as unlabeled minors); a non-uniform list (k_1h)
# drops the removed values outright.
RATE_YTICKS = {
    'k_1h': [0.0, 0.5, 1.0],        # even 0.5 steps (drop 0.25, 0.75)
    'k_13': [0.0, 2.0, 4.0],        # whole-even steps (minor at 1, 3)
    'k_14': [0.0, 1.0, 2.0],        # drop 0.5, 1.5
    'k_15': [0.0, 1.0, 2.0],        # drop 0.5, 1.5
    'k_16': [0.0, 1.0, 2.0, 3.0],   # re-cap 2.5 -> 3.0, whole-number steps
}

# same idea for feeding cells: an explicit major-tick list overrides the auto
# _linear_cap ticks and sets the axis top to its largest tick (uniform list ->
# half-step minors kept).
FEED_YTICKS = {
    'n_glu_spikes': [0.0, 4.0, 8.0],   # 0/4/8 (minor at 2, 6)
}


def _mathify(token):
    """Render a parameter / reaction token ('k_1l', 'k_13', 'r1', 'r13') as a
    matplotlib mathtext symbol: an italic base letter with the trailing
    identifier as a subscript. mathtext styles each subscript character by its
    class on its own -- digits upright, letters italic -- so '$k_{1l}$' yields
    an italic k with an upright '1' and an italic 'l'."""
    base, _, sub = token.partition('_')
    if not sub:                       # no underscore ('r13'): split off the digits
        i = 0
        while i < len(token) and not token[i].isdigit():
            i += 1
        base, sub = token[:i], token[i:]
    return f'${base}_{{{sub}}}$'


def bar_cell(ax, sets, colors, var, kind, ylabel, subtitle=None, ylim=None,
             band=None):
    """One parameter cell: colored bars per set with the searched band
    shaded and the baseline dashed. `ylabel` names the value axis (the
    parameter symbol / group / feeding quantity); `subtitle`, if given, is
    a smaller descriptor above the cell (reaction + enzyme). Value axes get
    conventional linear ticks (major + one minor between) on the rate,
    feeding and group cells alike."""
    n = len(sets)
    base = _baseline_value(sets, var)
    if kind == 'rate':
        # linear value axis from 0 to a nice ceiling just above the largest
        # bar (both edges on labeled ticks); no search-band shading
        data_max = 0.0
        for s in sets:
            v = s.get(var)
            if v is not None and np.isfinite(v):
                data_max = max(data_max, float(v))
        if base and base > 0:
            data_max = max(data_max, base)
        bottom = 0
        if var in RATE_YTICKS:
            ticks = RATE_YTICKS[var]
            ax.set_ylim(bottom, ticks[-1])
            ax.yaxis.set_major_locator(FixedLocator(ticks))
            uniform = len(set(np.round(np.diff(ticks), 6))) == 1
            ax.yaxis.set_minor_locator(
                AutoMinorLocator(2) if uniform else NullLocator())
        else:
            cap, step = _linear_cap(data_max)
            ax.set_ylim(bottom, cap)
            ax.yaxis.set_major_locator(MultipleLocator(step))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        if base and base > 0:
            ax.axhline(base, color='k', lw=0.8, ls='--', zorder=1)
    elif kind == 'group':
        # linear value axis from 0 to a nice ceiling just above the largest
        # bar (both edges on labeled ticks), baseline 1.0 dashed; no
        # search-band shading -- one minor tick between majors like the rates
        data_max = 0.0
        for s in sets:
            v = s.get(var)
            if v is not None and np.isfinite(v):
                data_max = max(data_max, float(v))
        data_max = max(data_max, 1.0)   # keep the baseline reference framed
        bottom = 0
        cap, _ = _linear_cap(data_max)
        # exactly three major ticks (0, cap/2, cap) with one minor between,
        # e.g. 0 / 0.5 / 1.0 on the ethanol cell, 0 / 1.0 / 2.0 elsewhere
        ax.set_ylim(bottom, cap)
        ax.yaxis.set_major_locator(FixedLocator([0.0, cap / 2.0, cap]))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.axhline(1.0, color='k', lw=0.8, ls='--', zorder=1)
    else:  # feed
        lo, _ = ylim
        # linear axis from the declared floor to a nice ceiling just above the
        # largest bar (both edges on labeled ticks); no search-band shading
        data_max = 0.0
        for s in sets:
            v = s.get(var)
            if v is not None and np.isfinite(v):
                data_max = max(data_max, float(v))
        if base is not None and np.isfinite(base):
            data_max = max(data_max, base)
        if var in FEED_YTICKS:
            ticks = FEED_YTICKS[var]
            ax.set_ylim(lo, ticks[-1])
            ax.yaxis.set_major_locator(FixedLocator(ticks))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        else:
            cap, step = _linear_cap(data_max)
            ax.set_ylim(lo, cap)
            ax.yaxis.set_major_locator(MultipleLocator(step))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        if base is not None:
            ax.axhline(base, color='k', lw=0.8, ls='--', zorder=1)
        bottom = 0
    ax.set_ylabel(ylabel, fontsize=FONTS['cell'], labelpad=3)
    if subtitle:
        ax.set_title(subtitle, fontsize=FONTS['tick'], pad=4)
    for j, s in enumerate(sets):
        v = s.get(var)
        c = colors[id(s)]
        if v is not None and np.isfinite(v) and v > 0:
            ax.bar(j, v - bottom, bottom=bottom, width=0.72, color=c, zorder=2)
        else:
            top = ax.get_ylim()[1]
            y = bottom * 1.25 if bottom else 0.01 * top
            label = 'n/a' if (v is None or (isinstance(v, float)
                                            and not np.isfinite(v))) else '0'
            # the baseline's Ehrlich branch is genuinely off (zero flux) and a
            # batch set runs zero glucose spikes; those "0" labels only clutter
            # the four Ehrlich cells and the spikes cell, so drop them
            drop = label == '0' and (
                (s.get('is_baseline') and var in EHRLICH_RATE_VARS)
                or var == 'n_glu_spikes')
            if not drop:
                ax.text(j, y, label, ha='center', va='bottom',
                        fontsize=FONTS['tick'], color=c)
    ax.set_xlim(-0.6, n - 0.4); style_cell_axes(ax)


def outcome_cell(ax, sets, colors, col, title, ylim, xmax, point_size=9):
    """One outcome metric as incumbent trajectories over trial_number: one
    step line per campaign set (the metric at that set's running incumbent,
    in its panel-b/c color), the scenario-A baseline as a dashed reference.
    """
    lo, hi = ylim

    def _yv(y):
        # every value is drawn at the axis floor lo when it has no finite
        # position: failures (no finite value: nan / +-inf) drop to lo instead
        # of being omitted (cloud) or gapping (line), and for the flagged
        # outcomes (IRR) finite losses (< lo) are drawn at lo too. So both the
        # cloud and the incumbent line dip to zero on a failure, no hole.
        y = np.asarray(y, dtype=float)
        out = np.where(np.isfinite(y), y, lo)
        if col in CLAMP_NEG_TO_ZERO:
            out = np.where(out < lo, lo, out)
        return out

    base = _baseline_value(sets, col)
    finite = [v for s in sets if s.get('traj') is not None
              for v in s['traj'][col] if np.isfinite(v)]
    if base is not None and np.isfinite(base):
        finite.append(base)
    vmax = max(finite) if finite else hi
    # bound the value axis top and bottom by labeled ticks (like IRR's 0->30):
    # a nice ceiling >= the data, no ~8% auto-overshoot past the last tick
    cap, step = _linear_cap(vmax, floor_hi=hi)
    step = OUTCOME_TICK_STEP.get(col, step)
    cap = OUTCOME_TICK_CAP.get(col, cap)
    ax.set_ylim(lo, cap)
    ax.set_xlim(0, xmax * 1.02)
    if base is not None and np.isfinite(base):
        ax.axhline(base, color=BASELINE_COLOR, lw=0.9, ls='--', zorder=1)
    # individual trial cloud from the ONE campaign that optimized THIS metric:
    # every trial as a translucent dot in the campaign's own color, behind the
    # incumbent lines (zorder 1). Failed / unsolved trials (no finite value)
    # are drawn at the axis floor rather than omitted.
    for s in sets:
        if s.get('scatter') is None or not _owns_outcome(s.get('objective'), col):
            continue
        sx, sy, cc = s['scatter_x'], _yv(s['scatter'][col]), colors[id(s)]
        ax.scatter(sx, sy, s=point_size, color=cc, alpha=0.35, linewidths=0,
                   zorder=1)
    for s in sets:
        if s.get('traj') is None:
            continue
        # the campaign that optimized THIS metric gets a ~2.25x-thick SOLID line
        # so its own trajectory stands out among the cross-plotted campaigns; the
        # thinner lines are drawn on top of it so none is hidden underneath, and
        # DASHED to mark that those campaigns were not optimizing this metric
        own = _owns_outcome(s.get('objective'), col)
        ax.step(s['traj_x'], _yv(s['traj'][col]), where='post',
                color=colors[id(s)], lw=3.15 if own else 1.4,
                ls='-' if own else (0, (2, 1.5)),
                zorder=2 if own else 3)
    # metric name next to the value axis itself (not a title above the cell)
    ax.set_ylabel(title, fontsize=FONTS['cell'], labelpad=3)
    ax.set_xlabel('Trial', fontsize=FONTS['tick'], labelpad=2)
    # the wide IRR cell fits five majors (0..2000 by 500); the narrow cells
    # take three (0, 1000, 2000). Four minor ticks sit between each major pair.
    ax.xaxis.set_major_locator(MultipleLocator(500 if col == 'IRR' else 1000))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(step))
    if col == 'IRR':   # fraction stored; show the value axis in percent
        ax.yaxis.set_major_formatter(
            PercentFormatter(xmax=1.0, decimals=0, symbol=''))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # ticks mirrored on all four sides (top/right mirror bottom/left); the
    # value/trial labels stay on the left/bottom only.
    ax.tick_params(axis='y', which='major', direction='inout',
                   left=True, right=True, labelright=False, length=4)
    ax.tick_params(axis='y', which='minor', direction='inout',
                   left=True, right=True, length=2.2)
    ax.tick_params(axis='x', which='major', top=True, bottom=True,
                   labeltop=False, labelbottom=True, direction='inout', length=4)
    ax.tick_params(axis='x', which='minor', top=True, bottom=True,
                   direction='inout', length=2.2)
    for sp in ('right', 'top'):
        ax.spines[sp].set_visible(True)
    # the mirrored top/right ticks point INWARD only (bottom/left stay in+out)
    _inward_top_right_ticks(ax, do_x=True, do_y=True)


def draw_outcomes(fig, gs_cell, sets, colors):
    xmax = 1.0
    for s in sets:
        tx = s.get('traj_x')
        if tx is not None and len(tx):
            xmax = max(xmax, float(tx[-1]))
    # 3x4 grid: IRR (the headline outcome) fills the left 3x2 block and is the
    # largest cell; the remaining six outcomes fill the right two columns
    # product-by-product -- the isobutanol metrics (yield, titer, productivity)
    # down column 3, the ethanol metrics down column 4, each top-to-bottom.
    sub_gs = gs_cell.subgridspec(3, 4, wspace=0.62, hspace=0.62)
    axes = []
    big, rest = OUTCOMES[0], OUTCOMES[1:]
    ax_big = fig.add_subplot(sub_gs[0:3, 0:2])
    outcome_cell(ax_big, sets, colors, big[0], big[1], big[2], xmax)
    axes.append(ax_big)
    big_w = ax_big.get_position().width
    for (col, label, yl), (r, c) in zip(rest,
                                        ((0, 2), (1, 2), (2, 2),
                                         (0, 3), (1, 3), (2, 3))):
        ax = fig.add_subplot(sub_gs[r, c])
        # these cells are ~0.38x the linear size of the IRR cell (1 of 4 columns
        # vs 2 columns + a wspace), so shrink the trial-cloud marker to match:
        # scale its AREA by (cell width / IRR width)^2, i.e. its diameter by the
        # linear ratio, from the actual rendered widths (robust to layout tweaks)
        ratio = ax.get_position().width / big_w
        outcome_cell(ax, sets, colors, col, label, yl, xmax,
                     point_size=9 * ratio ** 2)
        axes.append(ax)
    return axes


def draw_parameters(fig, gs_rows, sets, colors, band):
    axes = []
    for gs_row, (title, params) in zip(gs_rows, BANDS):
        sub_gs = gs_row.subgridspec(1, 5, wspace=0.75)
        last_ax = None
        for i, p in enumerate(params):
            ax = fig.add_subplot(sub_gs[0, i]); last_ax = ax
            ax.set_facecolor(PANEL_B_BG)     # panel-b block colour (see backing)
            axes.append(ax)
            if p in RATE_VARS:
                # reaction/enzyme descriptors go in the figure caption, not
                # above each cell -- keep only the rate symbol on the ylabel
                # (mathtext: italic base, subscript numbers upright / letters
                # italic)
                sym = RATE_CELL_LABEL.get(p) or _mathify(p)
                bar_cell(ax, sets, colors, p, 'rate', f'{sym}\n[{RATE_UNIT}]',
                         band=band)
            elif p in GROUP_VARS:
                # the "relative to baseline" unit lives in the band title now
                bar_cell(ax, sets, colors, p, 'group',
                         GROUP_LABELS[p], band=band)
            else:
                t, rng = FEED_LABELS[p]
                bar_cell(ax, sets, colors, p, 'feed', t, ylim=rng)
        # A band title is a header for the band BELOW it, so it should hug its
        # own cells and leave the larger gap to the band above. The feeding
        # cells carry long two-line vertical ylabels ("Threshold sugar
        # conc. ...") that overflow the short cell and poke above its top box,
        # so that band takes a middling offset to clear them; the rate and
        # group cells have short ylabels and take the tight offset that pins
        # the title to its band.
        if any(p in FEED_DRAW_VARS for p in params):
            offset = 0.022
        else:
            offset = 0.010
        # left-aligned with the plot boxes (the gridspec left margin)
        fig.text(0.083, last_ax.get_position().y1 + offset, title,
                 fontsize=FONTS['band'], fontweight='bold', va='bottom')
    return axes


# r17 (Adh6, k_17; nskinetics 2026-09-15 split) is now a sampled step -- the
# metabolic_split_12d/14d campaigns optimize k_17 -- so it carries an enzyme
# label (STEP_ENZYME / STEP_PARAMS / REACTION_LABELS) and enters the console
# pool ranking. A minimal_subset campaign simply leaves it at its baseline.
_STUDY_STEPS = ['r1', 'r3', 'r6', 'r13', 'r14', 'r15', 'r16', 'r17']
# never sampled by any campaign; folded into "other".
_UNSAMPLED_STEPS = ['r2', 'r4', 'r5']

if set(_STUDY_STEPS) | set(_UNSAMPLED_STEPS) != set(eb.STEP_ORDER):
    raise AssertionError(
        'burden step partition drift: _STUDY_STEPS | _UNSAMPLED_STEPS != '
        'eb.STEP_ORDER (%r vs %r)'
        % (sorted(set(_STUDY_STEPS) | set(_UNSAMPLED_STEPS)),
           sorted(eb.STEP_ORDER)))

# Panel c: a FIXED set of pathway categories drawn as the same stacked
# segments on every bar (a partition of the Phi_M pools). r2 is the PDH
# complex (TCA cycle); r4 (Ald6) -> r5 (Acs2) is the acetate bypass that
# regenerates cytosolic acetyl-CoA -- merged with the TCA category here.
# r17 (Adh6, the isobutyraldehyde -> isobutanol reduction of the 2026-09-15
# split) belongs to isobutanol production.
# The partition is asserted against eb.STEP_ORDER so a table drift raises
# at import.
BURDEN_CATEGORIES = (
    ('Glycolysis', ['r1']),
    ('TCA cycle + acetate /\nacetyl-CoA production', ['r2', 'r4', 'r5']),
    ('Ethanol production', ['r3', 'r6']),
    ('Isobutanol production', ['r13', 'r14', 'r15', 'r16', 'r17']),
)
_CAT_STEPS = [st for _, steps in BURDEN_CATEGORIES for st in steps]
if sorted(_CAT_STEPS) != sorted(eb.STEP_ORDER) \
        or len(_CAT_STEPS) != len(set(_CAT_STEPS)):
    raise AssertionError(
        'burden category partition drift: BURDEN_CATEGORIES steps != '
        'eb.STEP_ORDER (%r vs %r)' % (sorted(_CAT_STEPS),
                                      sorted(eb.STEP_ORDER)))

# sector styling for panel c. The four metabolic categories are filled in the
# campaign colour; a hatch tells them apart, EXCEPT the TCA category, which is
# drawn fill-only (hatch None). Housekeeping and translation are instead drawn
# HOLLOW -- no fill, with the bar's outline AND hatch both in the campaign
# colour (matplotlib ties a patch's hatch colour to its edge colour), so they
# read as the fixed non-metabolic anchors while still keying to the campaign.
# The unallocated flexible slack stays empty room (no fill, no hatch).
_METABOLIC_HATCHES = ('///', None, '...', 'xxx')
_TRANSLATION_HATCH = 'oo'
_HOUSEKEEPING_HATCH = '++'
if len(_METABOLIC_HATCHES) != len(BURDEN_CATEGORIES):
    raise AssertionError('one hatch per BURDEN_CATEGORIES entry required')


def draw_burden(fig, gs_cell, sets, colors):
    # full proteome allocation: each bar sums to the proteome cap PC (0.49).
    # Sectors left to right -- four metabolic categories | the unallocated
    # flexible slack | the growth-derated translation sector | housekeeping,
    # flush against the cap on the right. The four metabolic sectors are filled
    # in the set's colour and told apart by hatch (TCA fill-only); housekeeping
    # and translation are drawn hollow -- no fill, campaign-coloured outline and
    # hatch -- as the fixed non-metabolic anchors. Translation is
    # derated: metabolism fills the flexible sector F_flex first and the cell
    # builds only d * phi_T,demand = min(phi_T,demand, F_flex - Phi_M) of it,
    # so a single dashed line at F_flex - phi_T,demand marks how far translation
    # would reach un-derated, measured leftward from the housekeeping edge.
    #
    # Housekeeping is the SAME fixed 0.245-wide block on every bar (half the
    # proteome), now on the RIGHT, flush against the cap: it spans [0.245, 0.49].
    # A broken x-axis cuts a chunk out of the MIDDLE of it: the left (wide)
    # window keeps [0, BREAK_L], the right (stub) window resumes at BREAK_R and
    # runs to the cap, with BREAK_L and BREAK_R both strictly inside
    # housekeeping. So the reader still sees where housekeeping begins -- its
    # boundary with translation at 0.245 sits in the left window, a little
    # before BREAK_L -- and that it fills to the cap, in the right stub. Both
    # windows share ONE scale: their column width ratios equal their data
    # ranges, so a 0.05 tick step is the same physical distance on each side.
    # Housekeeping straddles the cut, so the full stack is drawn on BOTH windows
    # and each one clips it to its own range.
    cats = BURDEN_CATEGORIES
    n = len(sets)
    PC = float(eb.PROTEIN_CONTENT)
    housekeeping = PC * float(eb.HOUSEKEEPING_FRACTION)
    h = 0.62
    ypos = {id(s): n - i for i, s in enumerate(sets)}

    # the cut, strictly inside housekeeping (0.245 < BREAK_L < BREAK_R < 0.49):
    # the left window runs from the origin through the start of housekeeping,
    # the right stub shows housekeeping filling to the cap.
    BREAK_L, BREAK_R = 0.28, 0.46
    xmax = 0.5                       # end the right window on the 0.5 major tick
    # equal scale on both windows <=> width ratios == their data ranges
    sub = gs_cell.subgridspec(1, 2, width_ratios=[BREAK_L, xmax - BREAK_R],
                              wspace=0.025)
    axL = fig.add_subplot(sub[0, 0])                     # housekeeping stub
    axR = fig.add_subplot(sub[0, 1], sharey=axL)         # rest of the proteome

    def seg(ax, y, x, w, c, hatch=None, fill=True, outline=None, zorder=2):
        # fill=False draws empty room (white, faint grey outline). `outline`,
        # when given, sets the edge colour explicitly -- with fill=False this
        # draws a sector HOLLOW but keyed to the campaign: no fill, campaign-
        # coloured outline and (the hatch colour follows the edge) campaign-
        # coloured hatch. The solid metabolic sectors are drawn at a higher
        # zorder than the hollow ones so their black outlines render ON TOP of
        # a neighbouring hollow sector's campaign/grey outline at the shared
        # edge (not the other way round).
        ax.barh(y, w, left=x, height=h,
                facecolor=(c if fill else 'white'),
                edgecolor=(outline if outline is not None
                           else ('0.15' if fill else '0.6')),
                lw=0.5, hatch=hatch, zorder=zorder)

    # k_7/k_8 are pinned in this study, so the translation demand phi_T is the
    # same for every set and the un-derated marker is one vertical line.
    demands = [float(s['phi_T']) for s in sets]
    if max(demands) - min(demands) > 1e-4:
        raise ValueError('panel c assumes a shared translation demand phi_T '
                         '(k_7/k_8 pinned); sets differ: %r' % demands)
    # translation now sits against the housekeeping edge (F_flex = PC -
    # housekeeping), so the un-derated demand line falls F_flex - phi_T,demand
    # from the origin -- in the left window.
    demand_x = PC - housekeeping - demands[0]

    def draw_stack(ax, s, y, c):
        # the whole proteome bar; the axis window clips it to its own range.
        # order left->right: metabolic | slack | translation | housekeeping
        x = 0.0
        for (_, steps), hatch in zip(cats, _METABOLIC_HATCHES):  # metabolic
            w = sum(s[f'pool_{st}'] for st in steps)
            if w > 0:
                seg(ax, y, x, w, c, hatch, zorder=3)   # black edge on top
            x += w
        Phi_M = float(s['Phi_M'])
        phi_T_built = float(s['burden_factor']) * float(s['phi_T'])
        slack = max(0.0, PC - housekeeping - Phi_M - phi_T_built)
        if slack > 0:                                            # empty slack
            seg(ax, y, x, slack, c, fill=False)
        x += slack
        if phi_T_built > 0:                                      # translation
            seg(ax, y, x, phi_T_built, c, _TRANSLATION_HATCH,    # (hollow)
                fill=False, outline=c)
        x += phi_T_built
        seg(ax, y, x, housekeeping, c, _HOUSEKEEPING_HATCH,      # housekeeping
            fill=False, outline=c)                               # (hollow, at cap)

    for s in sets:
        y = ypos[id(s)]; c = colors[id(s)]
        draw_stack(axL, s, y, c)
        draw_stack(axR, s, y, c)

    # un-derated translation demand (in the left window), spanning just the
    # bar rows -- not the legend headroom above them
    axL.plot([demand_x, demand_x], [0.5, n + 0.6], color='0.15', ls='--',
             lw=1.0, zorder=4)   # above the raised metabolic sectors (zorder 3)
    # a second dashed line where the housekeeping sector begins (the flexible
    # sector edge, F_flex = PC - housekeeping); also in the left window
    hk_start = PC - housekeeping
    axL.plot([hk_start, hk_start], [0.5, n + 0.6], color='0.15', ls='--',
             lw=1.0, zorder=4)
    # the proteome cap at PC (the right end of every bar) -- a dashed line in
    # the right stub, labelled parallel to it (rotated 90 deg) in the empty gap
    # between the cap and the 0.5 edge
    axR.plot([PC, PC], [0.5, n + 0.6], color='0.15', ls='--', lw=1.0, zorder=4)
    axR.text(0.5 * (PC + xmax), 0.5 * (0.5 + n + 0.6), 'proteome cap',
             ha='center', va='center', rotation=90, rotation_mode='anchor',
             fontsize=11, color='0.15', clip_on=False, zorder=5)

    # sector-demand brackets, hovering just above the top bar: double-headed
    # arrows with square end caps and a centred label (as |<--- ... --->|).
    #  * penalty-free metabolic budget -- metabolism grows from the start of
    #    glycolysis (the origin) and only derates translation once it passes the
    #    un-derated demand line, so that span (F_flex - phi_T,demand) is free;
    #  * translation demand -- from the demand line to the housekeeping edge.
    by = n + 0.68                          # above the top border spine (n+0.50)
    cap = 0.10                             # end-cap half-height

    def bracket(ax, x0, x1, label):
        ax.annotate('', xy=(x1, by), xytext=(x0, by), annotation_clip=False,
                    arrowprops=dict(arrowstyle='<->', color='0.15', lw=1.0,
                                    shrinkA=0.0, shrinkB=0.0), zorder=5)
        for bx in (x0, x1):                 # vertical end caps ('|')
            ax.plot([bx, bx], [by - cap, by + cap], color='0.15', lw=1.0,
                    solid_capstyle='butt', clip_on=False, zorder=5)
        ax.text(0.5 * (x0 + x1), by + cap + 0.06, label, ha='center',
                va='bottom', fontsize=11, color='0.15', clip_on=False, zorder=5)

    bracket(axL, 0.0, demand_x, 'penalty-free metabolic budget')
    bracket(axL, demand_x, hk_start, 'translation demand')

    top_edge = n + 0.50                              # top border, below the brackets
    axL.set_ylim(0.4, n + 1.05)                      # shared: sets both windows
    for ax in (axL, axR):
        # rows are keyed by colour through the campaign legend, so the
        # categorical y axis carries no labels of its own
        ax.set_yticks([])
        ax.tick_params(axis='y', right=False, length=0)
        ax.tick_params(axis='x', which='major', direction='inout', top=False,
                       length=4)
        ax.tick_params(axis='x', which='minor', direction='inout', top=False,
                       length=2.2)
        # equal scale, so one tick cadence matches physically across the break
        ax.xaxis.set_major_locator(MultipleLocator(0.05))
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        # box the panel: the OUTER edges of the broken axis carry the left and
        # right spines (axL's left, axR's right), while the inner edges stay open
        # to read as the axis break. The top border + its mirrored ticks are a
        # secondary x-axis added below (at top_edge, just above the bars and
        # BELOW the sector-demand brackets that hover in the headroom); the
        # native top spine stays hidden and the left/right spines are bounded to
        # end at top_edge rather than run up through the bracket region. The
        # categorical y axis carries no ticks.
        ax.spines['top'].set_visible(False)
    axL.spines['left'].set_visible(True)
    axL.spines['left'].set_bounds(0.4, top_edge)
    axL.spines['right'].set_visible(False)
    axR.spines['right'].set_visible(True)
    axR.spines['right'].set_bounds(0.4, top_edge)
    axR.spines['left'].set_visible(False)
    axL.set_xlim(0, BREAK_L)
    axR.set_xlim(BREAK_R, xmax)
    # mirror the bottom x-axis onto the top border: a secondary x-axis at
    # top_edge gives a real spine + ticks there, and because each window keeps
    # its own xlim the break gap is reproduced on top. Ticks point inward only
    # (down), half-length like panels A/B; labels off.
    y0, y1 = axL.get_ylim()
    top_frac = (top_edge - y0) / (y1 - y0)
    for ax in (axL, axR):
        sax = ax.secondary_xaxis(top_frac)
        sax.xaxis.set_major_locator(MultipleLocator(0.05))
        sax.xaxis.set_minor_locator(AutoMinorLocator())
        sax.tick_params(axis='x', which='major', direction='in', length=2.0,
                        labeltop=False, labelbottom=False)
        sax.tick_params(axis='x', which='minor', direction='in', length=1.1,
                        labeltop=False, labelbottom=False)
    # if a break edge happens to land on a 0.05 major tick, its label would sit
    # right at the cut and crowd the break marks (or collide across the narrow
    # gap), so blank any label falling exactly on a break edge while the tick
    # marks stay. At the current 0.28/0.46 edges neither is on a tick, so this is
    # a no-op and the axis reads ... 0.25 // 0.50.
    def _blank_at(edge):
        return FuncFormatter(lambda x, pos:
                             '' if abs(x - edge) < 1e-9 else f'{x:.2f}')
    axL.xaxis.set_major_formatter(_blank_at(BREAK_L))
    axR.xaxis.set_major_formatter(_blank_at(BREAK_R))
    # diagonal break marks at the cut, fixed physical size (point markers) so
    # the unequal panel widths do not skew them
    mk = dict(marker=[(-1, -3.2), (1, 3.2)], markersize=7, linestyle='none',
              color='0.15', mec='0.15', mew=1.1, clip_on=False)
    axL.plot([1], [0], transform=axL.transAxes, **mk)          # bottom cut
    axR.plot([0], [0], transform=axR.transAxes, **mk)
    axL.plot([1], [top_frac], transform=axL.transAxes, **mk)   # top cut (at top_edge)
    axR.plot([0], [top_frac], transform=axR.transAxes, **mk)

    # one x-axis label and one sector legend, both centred over the whole
    # (broken) panel rather than either sub-panel
    box_l = axL.get_position(); box_r = axR.get_position()
    mid = 0.5 * (box_l.x0 + box_r.x1)
    fig.text(mid, box_l.y0 - 0.030,
             'Proteome allocation [g protein·(g DCW)$^{-1}$]',
             ha='center', va='top', fontsize=FONTS['axis'])

    # sector key, in the margin to the left of the narrowed panel: neutral
    # swatches so the hatches read independent of the campaign colours (metabolic
    # sectors grey-filled; the hollow housekeeping/translation sectors white with
    # a grey outline + hatch, matching how they are drawn on the bars). Grouped
    # under a bold 'Metabolic' header, then Translation and Housekeeping, with
    # the two flexible-sector annotations (slack, un-derated demand) last.
    handles = [Line2D([], [], linestyle='none', marker='none',
                      label='Metabolic')]                     # section header
    for (name, steps), hatch in zip(cats, _METABOLIC_HATCHES):  # indented members
        plain_ids = '(%s)' % ', '.join(steps)                  # width test (unrendered)
        ids = '(%s)' % ', '.join(_mathify(st) for st in steps)  # rendered (mathtext)
        lines = name.split('\n')
        if len(lines[-1]) + 1 + len(plain_ids) <= 28:          # ids fit inline
            lines[-1] += ' ' + ids
        else:                                                  # else wrap below
            lines.append(ids)
        label = '\n'.join('  ' + ln for ln in lines)           # 2-space indent
        handles.append(plt.Rectangle((0, 0), 1, 1, facecolor='0.72',
                                     edgecolor='0.15', lw=0.5, hatch=hatch,
                                     label=label))
    # the metabolic members, unallocated, translation and housekeeping run at
    # the plain labelspacing (no blank spacer rows): the bold, un-indented
    # headers already mark the group breaks, and dropping the two spacers pulls
    # the box bottom up to sit roughly level with the x-axis title. Unallocated
    # sits before Translation to follow the new bar order (metabolic |
    # unallocated | translation | housekeeping).
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor='white',
                                 edgecolor='0.6', lw=0.5,
                                 label='Unallocated'))
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor='white',
                                 edgecolor='0.15', lw=0.5,
                                 hatch=_TRANSLATION_HATCH,
                                 label='Translation'))
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor='white',
                                 edgecolor='0.15', lw=0.5,
                                 hatch=_HOUSEKEEPING_HATCH,
                                 label='Housekeeping'))
    # at panel-b's larger font the wrapped key is tall, so hang it from the top
    # of the panel (just under the title) rather than centring it -- centring
    # pushed its head into the panel title. It runs down the free left column.
    leg_c = fig.legend(handles=handles, loc='upper left',
                       bbox_to_anchor=(0.02, box_l.y1), ncol=1, frameon=True,
                       fontsize=FONTS['legend'] + 1, title='Sectors',
                       handlelength=1.6, handleheight=1.1, labelspacing=0.5,
                       borderpad=0.6, edgecolor='0.6', fancybox=False)
    leg_c.get_title().set_fontweight('bold')
    leg_c.get_title().set_fontsize(FONTS['legend'] + 2)       # match panel-b title
    for t in leg_c.get_texts():                               # bold the headers
        if t.get_text() in ('Metabolic', 'Translation', 'Housekeeping'):
            t.set_fontweight('bold')
    return axR


def plot(sets, band, out_stem, dpi=300, include_parameters=True):
    apply_fonts()
    LEFT, RIGHT = 0.083, 0.97
    if include_parameters:
        fig = plt.figure(figsize=(9.5, 13.454))
        # Top-anchored vertical layout, figure fractions. Panel a now holds a
        # 3x4 grid (was 2x4), so its box is 1.5x taller; the canvas grew by that
        # one extra row-height (12.4 -> 13.454 in) and panels b and c keep their
        # absolute inch heights, positions and gaps (their fractions rescaled by
        # 12.4/13.454). Panel a's bottom edge is unchanged in inches, so the wide
        # a -> b gap (which carries panel b's letter/title and the first band
        # title) is preserved. Panel b's four band cells are ~40% shorter than
        # its single-cell height (~0.0335 vs ~0.0558 of the OLD canvas); the
        # band -> band gaps stay ~0.05 (hspace 1.49 x the shorter cell) to keep
        # room for the rate-band titles. Panel c rides up just below the bands.
        # Each region is its own gridspec so the three vertical positions are
        # set directly.
        a_gs = fig.add_gridspec(1, 1, left=LEFT, right=RIGHT,
                                top=0.9493, bottom=0.7143)
        band_gs = fig.add_gridspec(4, 1, left=LEFT, right=RIGHT, top=0.6351,
                                   bottom=0.3733, hspace=1.49)
        # panel c is narrowed on the left (left=0.32 vs LEFT) to clear a column
        # for its framed sector legend, which sits in that margin rather than
        # above the bars (wide enough for the box + the longest label, clear of
        # the 0.00 tick)
        c_gs = fig.add_gridspec(1, 1, left=0.32, right=RIGHT,
                                top=0.3309, bottom=0.1235)
        colors = set_colors(sets)
        a_axes = draw_outcomes(fig, a_gs[0], sets, colors)
        b_axes = draw_parameters(fig, [band_gs[0], band_gs[1], band_gs[2],
                                       band_gs[3]], sets, colors, band)
        # light-grey backing behind the whole of panel b (band cells, their
        # titles, the panel-b letter/title and the campaign legend), so the
        # panel reads as one block set off from the white panels a and c. Drawn
        # at zorder 0 so the cells, bars, text and legend all sit on top; the
        # band cells share its colour (set in draw_parameters) so the fill is
        # seamless across the gaps.
        b_top = b_axes[0].get_position().y1
        b_bot = b_axes[-1].get_position().y0
        fig.add_artist(plt.Rectangle(
            (0.02, b_bot - 0.030), 0.985 - 0.02,
            (b_top + 0.058) - (b_bot - 0.030),
            transform=fig.transFigure, facecolor=PANEL_B_BG, edgecolor='none',
            zorder=0))
        axc = draw_burden(fig, c_gs[0], sets, colors)
        # each panel gets a bold letter and a descriptive title on the same
        # baseline; the panel title (13 pt) outranks the band sub-titles
        # (12 pt). Panel b's letter is lifted into the panel-a -> b gap so its
        # title clears the first band title below it.
        panels = ((a_axes[0].get_position().y1 + 0.012, 'A',
                   'Optimization incumbent trajectories'),
                  (b_axes[0].get_position().y1 + 0.035, 'B',
                   'Final kinetic and process parameters'),
                  (axc.get_position().y1 + 0.005, 'C',
                   'Final proteome allocation'))
        # the campaign legend sits in the empty right columns of panel b's lower
        # bands and doubles as the row key for panel c (whose categorical y axis
        # is unlabelled).
        legend_loc, legend_anchor, legend_ncol = 'center left', (0.635, 0.4147), 1
    else:
        # two-panel variant (panel B omitted): panel A on top, the proteome-
        # allocation panel below it -- relettered B. A shorter canvas keeps both
        # panels at ~their three-panel absolute heights; the freed middle band
        # becomes the a -> b gap that carries panel b's letter/title and the
        # campaign legend.
        # panel a holds a 3x4 grid (was 2x4), 1.5x taller; the canvas grew by
        # that one extra row-height (8.8 -> 9.856 in) and panel c keeps its
        # absolute inch height, position and the a -> b gap (fractions rescaled
        # by 8.8/9.856; panel a's bottom edge is unchanged in inches).
        fig = plt.figure(figsize=(9.5, 9.856))
        a_gs = fig.add_gridspec(1, 1, left=LEFT, right=RIGHT,
                                top=0.9330, bottom=0.6116)
        c_gs = fig.add_gridspec(1, 1, left=0.32, right=RIGHT,
                                top=0.4598, bottom=0.1741)
        colors = set_colors(sets)
        a_axes = draw_outcomes(fig, a_gs[0], sets, colors)
        axc = draw_burden(fig, c_gs[0], sets, colors)
        panels = ((a_axes[0].get_position().y1 + 0.012, 'A',
                   'Optimization incumbent trajectories'),
                  (axc.get_position().y1 + 0.005, 'B',
                   'Final proteome allocation'))
        # the campaign legend sits in the a -> b gap, doubling as the row key for
        # the proteome panel below. A single vertical column is too tall and
        # clips panel a's lower cells; a single horizontal row of all campaign
        # items runs off the figure once there are 6 campaigns (7 sets, the
        # widest being "Baseline (no optimization)"), so it wraps to at most four
        # columns -- two rows under the title.
        legend_loc = 'center'
        legend_anchor, legend_ncol = (0.527, 0.5357), min(len(sets), 4)
    for y, letter, title in panels:
        fig.text(0.03, y, letter, fontsize=FONTS['panel'], fontweight='bold',
                 va='baseline')
        fig.text(0.055, y, title, fontsize=FONTS['panel'] - 1,
                 fontweight='bold', va='baseline')
    handles = [plt.Rectangle((0, 0), 1, 1, fc=colors[id(s)], label=s['label'])
               for s in sets]
    leg = fig.legend(handles=handles, loc=legend_loc,
                     bbox_to_anchor=legend_anchor, ncol=legend_ncol, frameon=True,
                     fontsize=FONTS['legend'] + 1, title='Optimization campaign',
                     labelspacing=0.5, handlelength=1.7, handleheight=1.1,
                     borderpad=0.6, edgecolor='0.6', fancybox=False)
    leg.get_title().set_fontweight('bold')
    leg.get_title().set_fontsize(FONTS['legend'] + 2)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out_stem}.{ext}', dpi=dpi)
    plt.close(fig)
    return out_stem


def console_report(sets, band_campaign):
    print(f'band source (campaign): {band_campaign}')
    rank = sorted(_STUDY_STEPS,
                  key=lambda st: -max(s[f'pool_{st}'] for s in sets))[:5]
    for s in sets:
        tag = 'baseline' if s.get('is_baseline') \
            else (f'{s["campaign"]} trial {int(s["trial_number"])}; '
                  f'panel-a incumbent by {s["sel_col"]} ({s["sel_dir"]})')
        print(f'[{s["label"]}] {tag}')
        print(f'    IRR {s.get("IRR")!s:>8}  EtOH {s.get("EtOH titer")!s:>7}'
              f'  IBO {s.get("IBO titer")!s:>7}  tau {s.get("tau")!s:>6}')
        print(f'    Phi_M {s["Phi_M"]:.4f}  d {s["burden_factor"]:.2f}  pools: '
              + ', '.join(f'{st} {s[f"pool_{st}"]:.4f}' for st in rank))
        if s.get('extra_sampled'):
            print(f'    campaign {s["campaign"]}: also sampled '
                  f'{s["extra_sampled"]} (not shown)')


def build_sets(set_specs, include_baseline):
    sets = [baseline_set()] if include_baseline else []
    band = band_campaign = None
    for label, campaign, trial in set_specs:
        sets.append(load_set(label, campaign, trial))
        if band is None:
            band = campaign_band(campaign); band_campaign = campaign
        else:
            other = campaign_band(campaign)
            if other != band:
                print(f'  WARNING campaign {campaign}: searched band differs '
                      f'from the first campaign ({band_campaign}); keeping '
                      'the first campaign band for shading')
    if not sets:
        raise ValueError('no sets to plot (use --set or drop --no-baseline)')
    return sets, band, band_campaign


def main(argv=None):
    ap = argparse.ArgumentParser(description='Kinetic-optimization '
                                 'parameter-set comparison figure.')
    ap.add_argument('--set', dest='sets', action='append', nargs=3,
                    metavar=('LABEL', 'CAMPAIGN', 'TRIAL'), default=None,
                    help='a set to plot; repeatable, in bar order after the '
                         'baseline. TRIAL is an int, "best", or "best:COL".')
    ap.add_argument('--no-baseline', action='store_true',
                    help='drop the scenario-A baseline row')
    ap.add_argument('--no-panel-b', action='store_true',
                    help='omit panel B (final kinetic and process parameters); '
                         'the proteome-allocation panel then becomes panel B')
    ap.add_argument('--out-dir', default=RESULTS_DIR)
    ap.add_argument('--stem', default=None)
    ap.add_argument('--dpi', type=int, default=300)
    args = ap.parse_args(argv)

    def norm_trial(t):
        if isinstance(t, str) and t.startswith('best'):
            return t
        return int(t)

    if args.sets:
        specs = [(lab, camp, norm_trial(tr)) for lab, camp, tr in args.sets]
    else:  # default: each optimum from the study that optimized it -- the
        # five most recent minimal-subset campaigns, one per objective
        specs = [('Financial attractiveness', DEFAULT_STUDY, 'best'),
                 ('Isobutanol titer', IBO_TITER_STUDY, 'best'),
                 ('Ethanol titer', ETOH_TITER_STUDY, 'best'),
                 ('Isobutanol yield', IBO_YIELD_STUDY, 'best'),
                 ('Ethanol yield', ETOH_YIELD_STUDY, 'best')]

    if len(specs) + (0 if args.no_baseline else 1) > MAX_SETS:
        raise ValueError(f'at most {MAX_SETS} sets (baseline + '
                         f'{len(HUE_COLORS)} campaign sets)')

    # panel-B layout follows the plotted preset. The figure is homogeneous
    # (one decision-variable set per figure), so the first --set campaign
    # selects it; the default no-arg path stays on the minimal_subset layout.
    stems = [os.path.basename(c) for _, c, _ in specs]
    if any('metabolic_split_12d' in st for st in stems):
        if not all('metabolic_split_12d' in st for st in stems):
            raise ValueError('mixing metabolic_split_12d and other presets in '
                             'one figure is not supported (decision variables '
                             'differ); plot them separately')
        use_split_12d_layout()

    sets, band, band_campaign = build_sets(specs, not args.no_baseline)
    stem = args.stem or f'{os.path.basename(specs[0][1]).replace(".csv", "")}' \
                        '_parameter_sets'
    stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
    out_stem = os.path.join(args.out_dir, f'{stem}_{stamp}')
    plot(sets, band, out_stem, dpi=args.dpi,
         include_parameters=not args.no_panel_b)
    console_report(sets, band_campaign)
    print(f'wrote {out_stem}.png / .pdf')
    return out_stem


if __name__ == '__main__':
    main()
