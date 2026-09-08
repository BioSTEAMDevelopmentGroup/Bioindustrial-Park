#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Offline pure-logic test of plot_proteome_voronoi (no load, no
simulation). Loads the Stage-1 plotter by file path and checks the
baseline tile's allocation against the eb constants and the spec's
verification numbers. Exit 0 + ALL CHECKS PASSED = clean."""
import os
import importlib.util

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pv = _load('plot_proteome_voronoi',
           os.path.join(PKG_DIR, 'plots', 'plot_proteome_voronoi.py'))
eb = pv.eb
checks = []


def check(name, cond):
    checks.append((name, bool(cond)))
    print(f'  [{"PASS" if cond else "FAIL"}] {name}')


# --- baseline tile ----------------------------------------------------------
rec = pv.ps.baseline_set()
tile = pv.tile_from_record(rec)
kids = {c['piece']: c for c in tile['children']}
PC = eb.PROTEIN_CONTENT

check('four top-level pieces', len(tile['children']) == 4)
check('housekeeping == 0.5*PROTEIN_CONTENT',
      abs(kids['housekeeping']['value'] - PC * eb.HOUSEKEEPING_FRACTION) < 1e-9)
check('housekeeping == 0.245', abs(kids['housekeeping']['value'] - 0.245) < 1e-6)
check('translation == record phi_T',
      abs(kids['translation']['value'] - float(rec['phi_T'])) < 1e-9)
check('baseline translation ~= 0.110 (PHI_T_WT, g=1)',
      abs(kids['translation']['value'] - eb.PHI_T_WT) < 1e-3)

metabolic = kids['metabolic']
mkids = {c['piece']: c for c in metabolic['children']}
check('metabolic has five categories', len(metabolic['children']) == 5)
check('metabolic value == sum of five categories',
      abs(metabolic['value'] - sum(c['value'] for c in metabolic['children']))
      < 1e-12)
check('metabolic value == record Phi_M',
      abs(metabolic['value'] - float(rec['Phi_M'])) < 1e-4)
check('baseline isobutanol category == 0 (Ehrlich branch off at A)',
      abs(mkids['cat_isobutanol']['value']) < 1e-9)

slack = kids['slack']['value']
check('slack == F_flex - Phi_M - phi_T',
      abs(slack - (eb.F_FLEX - metabolic['value']
                   - kids['translation']['value'])) < 1e-9)
check('slack >= 0 at baseline', slack > -1e-9)

total = (kids['housekeeping']['value'] + metabolic['value']
         + kids['translation']['value'] + slack)
check('all cells sum to PROTEIN_CONTENT (0.49)', abs(total - PC) < 1e-6)
check('no warning on the feasible baseline', 'warning' not in tile)


def _synthetic(phi_m, phi_t_demand, label):
    """A record shaped like a set row with all metabolic pool in Glycolysis
    (r1) so category-sum Phi_M == phi_m, and burden_factor = the model's
    derating d = clip((F_flex - Phi_M)/phi_T, 0, 1)."""
    s = dict(rec)
    for st in eb.STEP_ORDER:
        s[f'pool_{st}'] = 0.0
    s['pool_r1'] = phi_m
    s['Phi_M'] = phi_m
    s['phi_T'] = phi_t_demand
    s['F_flex'] = eb.F_FLEX
    d = (eb.F_FLEX - phi_m) / phi_t_demand
    s['burden_factor'] = max(0.0, min(1.0, d))
    s['label'] = label
    s['is_baseline'] = False
    return s


# --- derated-growth trial (feasible: Phi_M < F_flex, but Phi_M + demand >
#     F_flex so d < 1): translation cell is the derated allocation, slack 0,
#     no warning, cells still sum to 0.49 ---------------------------------------
dr = pv.tile_from_record(_synthetic(0.20, 0.10, 'Derated'))
dk = {c['piece']: c for c in dr['children']}
check('derated translation == F_flex - Phi_M (d*phi_T,demand)',
      abs(dk['translation']['value'] - (eb.F_FLEX - 0.20)) < 1e-9)
check('derated translation < demand', dk['translation']['value'] < 0.10 - 1e-9)
check('derated slack == 0', abs(dk['slack']['value']) < 1e-9)
check('derated tile has no warning (feasible)', 'warning' not in dr)
check('derated cells sum to PROTEIN_CONTENT (0.49)',
      abs(sum(c['value'] for c in dr['children']) - PC) < 1e-6)

# --- truly infeasible trial (Phi_M >= F_flex, d = 0): no room for any
#     translation; slack clamped to 0 + warning ---------------------------------
inf = pv.tile_from_record(_synthetic(0.30, 0.10, 'Infeasible'))
ik = {c['piece']: c for c in inf['children']}
check('infeasible translation == 0', abs(ik['translation']['value']) < 1e-12)
check('infeasible slack clamped to 0', abs(ik['slack']['value']) < 1e-12)
check('infeasible tile carries a warning', 'warning' in inf)

# --- document shape ---------------------------------------------------------
doc = pv.build_document([rec], band_campaign=None)
check('meta.protein_content == 0.49', abs(doc['meta']['protein_content'] - 0.49) < 1e-9)
check('meta.F_flex == eb.F_FLEX', abs(doc['meta']['F_flex'] - eb.F_FLEX) < 1e-9)
check('one tile', len(doc['tiles']) == 1)
check('piece_order has 8 entries', len(doc['meta']['piece_order']) == 8)

# --- per-tile base color (matches the parameter-sets figure) ----------------
check('baseline tile base_color == ps.BASELINE_COLOR',
      doc['tiles'][0]['base_color'] == pv.ps.BASELINE_COLOR)
two = pv.build_document([pv.ps.baseline_set(),
                         _synthetic(0.20, 0.10, 'Campaign')], band_campaign=None)
check('baseline tile keeps the grey base_color',
      two['tiles'][0]['base_color'] == pv.ps.BASELINE_COLOR)
check('first campaign tile base_color == ps.HUE_COLORS[0]',
      two['tiles'][1]['base_color'] == pv.ps.HUE_COLORS[0])

# --- CLI --no-render writes a valid JSON document ---------------------------
import json as _json
import tempfile

with tempfile.TemporaryDirectory() as td:
    # baseline-only, no campaigns -> no CSV reads required
    argv = ['--no-baseline', '--out-dir', td, '--stem', 'unit', '--no-render']
    # a baseline-only run needs --set-free defaults dropped; force one tile via
    # a direct document write instead of the study-CSV default path:
    doc = pv.build_document([pv.ps.baseline_set()], band_campaign=None)
    path = pv.write_document(doc, td, 'unit', 'stamp')
    with open(path, encoding='utf-8') as fh:
        loaded = _json.load(fh)
    check('write_document round-trips JSON',
          loaded['tiles'][0]['label'] == pv.ps.baseline_set()['label'])
    check('JSON tile cells sum to 0.49',
          abs(sum((c['value'] for c in loaded['tiles'][0]['children']))
              - 0.49) < 1e-6)

n_fail = sum(1 for _n, ok in checks if not ok)
print()
if n_fail:
    print(f'{n_fail} CHECK(S) FAILED')
    raise SystemExit(1)
print('ALL CHECKS PASSED')
