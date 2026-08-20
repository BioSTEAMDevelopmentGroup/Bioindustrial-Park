#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Verification of the standalone IBO/EtOH/water separation system
(separations.py) against its success criteria:

1. >90% recovery of ethanol to the ethanol product stream.
2. >99 wt% purity of the (pre-denaturant) ethanol product stream.
3. >90% recovery of isobutanol.
4. >95 wt% purity of the isobutanol product stream.
5. Material & energy flows stable across 3 consecutive re-simulations.
6. Criteria 1-5 hold for the scenario B ``P301-0`` baseline feed and for
   variations with higher/lower isobutanol and ethanol flows.

Loads separations.py by file path (fresh kernel; does not trigger the
package's import-time baseline simulation). Exits non-zero on any failure.
"""

import os
import sys
import importlib.util
import traceback
import numpy as np
import biosteam as bst
import thermosteam as tmo

#%% Load the separations module standalone (no package __init__)
_here = os.path.dirname(os.path.abspath(__file__))
_sep_path = os.path.join(os.path.dirname(_here), 'separations.py')
_spec = importlib.util.spec_from_file_location('ibo_separations', _sep_path)
separations = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(separations)

#%% Build the system once
tmo.settings.set_thermo(separations.create_separation_chemicals())
bst.main_flowsheet.set_flowsheet('IBO_EtOH_sep_test')

feed = separations.create_scenario_B_feed()
sep_sys = separations.create_IBO_EtOH_separation_system(ins=feed)
sep_sys.set_tolerance(rmol=1e-5, mol=1e-3, subsystems=True)

f = bst.main_flowsheet
ref_flows = dict(separations.scenario_B_P301_flows)

def set_feed(EtOH_mult, IBO_mult):
    for ID, m in ref_flows.items():
        feed.imass[ID] = m
    feed.imass['Ethanol'] *= EtOH_mult
    feed.imass['Isobutanol'] *= IBO_mult
    feed.T, feed.P = 305.15, 607950.0

def snapshot():
    mat = np.array([s.mol for s in sep_sys.streams], float).flatten()
    duties = []
    for u in sorted(sep_sys.units, key=lambda u: u.ID):
        duties.extend(hu.duty for hu in u.heat_utilities)
        duties.append(u.power_utility.rate)
    return mat, np.array(duties, float)

#%% Criteria checker
def check(name, EtOH_mult, IBO_mult):
    print(f'\n===== {name} (EtOH x{EtOH_mult}, IBO x{IBO_mult}) =====',
          flush=True)
    set_feed(EtOH_mult, IBO_mult)
    feed_EtOH = feed.imass['Ethanol']
    feed_IBO = feed.imass['Isobutanol']
    snaps = []
    for _ in range(3):
        sep_sys.simulate()
        snaps.append(snapshot())

    drift_mat = drift_duty = 0.0
    for a, b in ((0, 1), (1, 2), (0, 2)):
        ma, da = snaps[a]
        mb, db = snaps[b]
        drift_mat = max(drift_mat, np.abs(ma - mb).max() / np.abs(ma).max())
        drift_duty = max(drift_duty,
                         np.abs(da - db).max() / max(np.abs(da).max(), 1.0))

    etoh = f.stream.ethanol_product
    ibo = f.stream.isobutanol_product
    rec_EtOH = etoh.imass['Ethanol'] / feed_EtOH
    pur_EtOH = etoh.imass['Ethanol'] / etoh.F_mass
    rec_IBO = ibo.imass['Isobutanol'] / feed_IBO
    pur_IBO = ibo.imass['Isobutanol'] / ibo.F_mass

    checks = [
        ('EtOH recovery > 0.90', rec_EtOH, rec_EtOH > 0.90),
        ('EtOH purity   > 0.99', pur_EtOH, pur_EtOH > 0.99),
        ('IBO  recovery > 0.90', rec_IBO, rec_IBO > 0.90),
        ('IBO  purity   > 0.95', pur_IBO, pur_IBO > 0.95),
        ('mol drift    < 1e-3', drift_mat, drift_mat < 1e-3),
        ('energy drift < 1e-3', drift_duty, drift_duty < 1e-3),
    ]
    ok = all(c[2] for c in checks)
    for label, val, passed in checks:
        print(f'  {label}: {val:.6g}  [{"ok" if passed else "FAIL"}]')
    print(f'  => {"PASS" if ok else "FAIL"}')
    return ok

#%% Run all feed variants
def main():
    variants = [
        ('baseline scenario B', 1.0, 1.0),
        ('higher IBO (x2)', 1.0, 2.0),
        ('lower IBO (x0.5)', 1.0, 0.5),
        ('lower EtOH (x0.5)', 0.5, 1.0),
        ('higher EtOH (x1.5)', 1.5, 1.0),
        ('higher IBO, lower EtOH', 0.5, 2.0),
    ]
    results = {}
    for name, em, im in variants:
        try:
            results[name] = check(name, em, im)
        except Exception:
            traceback.print_exc()
            results[name] = False
    print('\n===== SUMMARY =====')
    for k, v in results.items():
        print(f'  {k:<26s} {"PASS" if v else "FAIL"}')
    return int(not all(results.values()))

if __name__ == '__main__':
    sys.exit(main())
