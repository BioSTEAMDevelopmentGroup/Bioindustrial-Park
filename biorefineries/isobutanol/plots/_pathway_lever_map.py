#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Enzyme lever-map figure -- the fourth variant of the kinetic-optimization
parameter-set comparison figure (companion to plot_kin_opt_parameter_sets.py).

Re-presents the metabolic_split_12d decision variables ENZYME-FIRST, laid on
the glucose -> ethanol / isobutanol carbon backbone, each lever shown as a
fold-change from wild type (or an absolute newly-expressed rate), and closes
to profitability with an outcome ribbon (max-product != max-profit).

Sim-safe: this module imports only matplotlib and numpy. The parent passes in
its already-file-path-loaded ko / eb modules and a bundle of its own
callables/tables (`helpers`); nothing here loads a module by path or imports
the biorefineries package. See plot_pathway_lever_map."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

# --- pathway node -> lever mapping (metabolic_split_12d layout) --------------
# step id -> (record key, kind). 'fold' = a x-baseline fold-change lever
# (baseline 1.0: glycolysis group + the k_*_rel fermentation/Aro10/Adh6 cells);
# 'absolute' = a newly-expressed Ehrlich enzyme, off (0) in wild type, drawn in
# g/L/h. Keys match what load_set / baseline_set expose under use_split_12d_layout
# (SPLIT_12D_REL_RATE_VARS -> k_*_rel; SPLIT_12D_RATE_VARS -> k_13/k_14/k_15).
NODE_LEVERS = {
    'r1':  ('glycolysis', 'fold'),
    'r3':  ('k_3_rel', 'fold'),
    'r6':  ('k_6_rel', 'fold'),
    'r13': ('k_13', 'absolute'),
    'r14': ('k_14', 'absolute'),
    'r15': ('k_15', 'absolute'),
    'r16': ('k_16_rel', 'fold'),
    'r17': ('k_17_rel', 'fold'),
}
# product-inhibition multipliers (baseline 1.0, fold-change). Larger multiplier
# => larger exp(-k*[product]) exponent in the nskinetics rate laws (r1/r4/r6/
# r7/r10/r17) => MORE inhibition => LESS tolerant, so a bar above 1x reads as
# less tolerant and engineered tolerance rises as the multiplier drops below 1x.
TOLERANCE_LEVERS = ('inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate')
# outcome-ribbon columns (the punchline): financial + the two product titers.
OUTCOME_COLS = ('IRR', 'IBO titer', 'EtOH titer')


# --- pure-logic value helpers (unit-tested offline) -------------------------
def lever_value(rec, key):
    """The lever's plotted value for set record `rec`: rec[key] (a fold-change
    for a 'fold' lever, g/L/h for an 'absolute' lever), or nan when the key is
    absent or the value is non-finite."""
    v = rec.get(key)
    try:
        v = float(v)
    except (TypeError, ValueError):
        return np.nan
    return v if np.isfinite(v) else np.nan


def is_off_in_wildtype(baseline_rec, key):
    """True if an absolute-lever enzyme is genuinely off (0 or absent) in the
    wild type, so its 'off in wild type' marker is drawn and its baseline bar
    suppressed."""
    v = lever_value(baseline_rec, key)
    return not (np.isfinite(v) and v > 0.0)


def clamp_outcome(col, value, floor, clamp_neg):
    """Outcome-ribbon value with the parent's CLAMP_NEG_TO_ZERO treatment: for
    a clamp column (IRR) a finite loss or -inf is drawn at `floor`; nan (never
    solved) stays nan so the caller omits that bar."""
    v = float(value)
    if np.isnan(v):
        return np.nan
    if col in clamp_neg:
        if not np.isfinite(v) or v < floor:   # -inf or a finite loss -> floor
            return floor
    if not np.isfinite(v):
        return np.nan
    return v


def proteome_segments(rec, eb, categories):
    """Reduce a record's burden pools to the slim stacked-strip segments:
    the isobutanol-pathway pool, the rest of the modeled metabolic pool, the
    flexible slack, the derated translation sector, and housekeeping. Sums to
    eb.PROTEIN_CONTENT. `categories` is helpers['BURDEN_CATEGORIES'] (the parent's
    ((name, [steps]), ...) partition of Phi_M)."""
    PC = float(eb.PROTEIN_CONTENT)
    housekeeping = PC * float(eb.HOUSEKEEPING_FRACTION)
    ibo_steps = next(steps for name, steps in categories
                     if name == 'Isobutanol production')
    ibo = sum(float(rec[f'pool_{st}']) for st in ibo_steps)
    Phi_M = float(rec['Phi_M'])
    rest = max(0.0, Phi_M - ibo)
    phi_T_built = float(rec['burden_factor']) * float(rec['phi_T'])
    slack = max(0.0, PC - housekeeping - Phi_M - phi_T_built)
    return {'isobutanol': ibo, 'rest': rest, 'slack': slack,
            'translation': phi_T_built, 'housekeeping': housekeeping}
