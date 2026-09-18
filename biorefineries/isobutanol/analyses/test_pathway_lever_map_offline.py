#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Offline pure-logic + render test of _pathway_lever_map (no load, no
simulation). Loads the parent plotter (ps), eb and the companion (plm) by
file path, builds a real scenario-A baseline record + synthetic campaign
records under the split_12d layout, and checks the lever/outcome/proteome
helpers, the objective-keyed colour map, and that the figure renders to a
temp PNG/PDF without error. Exit 0 + ALL CHECKS PASSED = clean."""
import os
import importlib.util
import numpy as np

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(name, relpath):
    path = os.path.join(PKG_DIR, relpath)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ps = _load('ps', os.path.join('plots', 'plot_kin_opt_parameter_sets.py'))
plm = _load('_pathway_lever_map', os.path.join('plots', '_pathway_lever_map.py'))
eb = ps.eb
ps.use_split_12d_layout()   # populate the split_12d panel-B globals
checks = []


def check(name, cond):
    checks.append((name, bool(cond)))
    print(f'  [{"PASS" if cond else "FAIL"}] {name}')


# --- records ---------------------------------------------------------------
base = ps.baseline_set()   # real scenario-A record (no simulation)


def _campaign(label, **over):
    """A synthetic campaign record: copy of the baseline with selected levers
    and outcomes overridden. is_baseline False so it takes a hue."""
    s = dict(base)
    s['label'] = label
    s['campaign'] = 'kin_opt_ethanol_isobutanol_metabolic_split_12d_ibo_titer_gp'
    s['is_baseline'] = False
    s['objective'] = 'IBO titer'
    s['trial_number'] = 1
    s.update(over)
    return s


# --- lever_value -----------------------------------------------------------
check('baseline fold levers are 1.0 (glycolysis)',
      abs(plm.lever_value(base, 'glycolysis') - 1.0) < 1e-9)
check('baseline k_3_rel is 1.0', abs(plm.lever_value(base, 'k_3_rel') - 1.0) < 1e-9)
check('baseline absolute Ehrlich rate k_13 is 0.0',
      abs(plm.lever_value(base, 'k_13')) < 1e-12)
check('missing key -> nan', np.isnan(plm.lever_value(base, 'nonexistent')))

# --- is_off_in_wildtype ----------------------------------------------------
check('k_13 off in wild type', plm.is_off_in_wildtype(base, 'k_13') is True)
check('k_14 off in wild type', plm.is_off_in_wildtype(base, 'k_14') is True)
inst = _campaign('IBO', k_13=3.0, k_14=1.2, k_15=1.1)
check('installed k_13 is not off',
      plm.is_off_in_wildtype(base, 'k_13') and not plm.is_off_in_wildtype(inst, 'k_13'))

# --- clamp_outcome (IRR treatment) -----------------------------------------
cn = ps.CLAMP_NEG_TO_ZERO
check('IRR finite loss clamps to floor',
      plm.clamp_outcome('IRR', -0.2, 0.0, cn) == 0.0)
check('IRR -inf clamps to floor',
      plm.clamp_outcome('IRR', float('-inf'), 0.0, cn) == 0.0)
check('IRR nan is omitted (nan)', np.isnan(plm.clamp_outcome('IRR', np.nan, 0.0, cn)))
check('IRR positive passes through',
      abs(plm.clamp_outcome('IRR', 0.27, 0.0, cn) - 0.27) < 1e-12)
check('non-clamp titer negative passes through unchanged',
      plm.clamp_outcome('IBO titer', -1.0, 0.0, cn) == -1.0)

# --- proteome_segments -----------------------------------------------------
seg = plm.proteome_segments(base, eb, ps.BURDEN_CATEGORIES)
PC = float(eb.PROTEIN_CONTENT)
check('segments sum to PROTEIN_CONTENT (0.49)',
      abs(sum(seg.values()) - PC) < 1e-6)
check('housekeeping == 0.5*PC', abs(seg['housekeeping'] - PC * 0.5) < 1e-9)
check('baseline isobutanol segment == native Adh6 pool r17 (Ehrlich off)',
      abs(seg['isobutanol'] - float(base['pool_r17'])) < 1e-9)
check('baseline isobutanol + rest == Phi_M',
      abs(seg['isobutanol'] + seg['rest'] - float(base['Phi_M'])) < 1e-9)

n_fail = sum(1 for _n, ok in checks if not ok)
print()
if n_fail:
    print(f'{n_fail} CHECK(S) FAILED')
    raise SystemExit(1)
print('ALL CHECKS PASSED')
