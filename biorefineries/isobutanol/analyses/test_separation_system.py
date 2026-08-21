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
   variations of the isobutanol and ethanol flows spanning ~0-200 g/L
   titer equivalents. At very low titers (multiplier <= 0.05) reduced
   recovery of the diminished alcohol is acceptable (floor of 60%), but
   product purities must hold whenever the product stream exists, and at
   exactly zero titer the corresponding product stream must be empty.

Variants are simulated SEQUENTIALLY on one system (no cache resets), like
a titer sweep would, so state carried between operating points is part of
what is being tested.

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

#: Product stream counts as "empty" below this total mass flow [kg/hr];
#: the reference feed carries ~1e5 kg/hr of water.
EMPTY_F_MASS = 1e-2

def set_feed(EtOH_mult, IBO_mult):
    for ID, m in ref_flows.items():
        feed.imass[ID] = m
    feed.imass['Ethanol'] *= EtOH_mult
    feed.imass['Isobutanol'] *= IBO_mult
    feed.T, feed.P = 305.15, 607950.0

def snapshot():
    mat = np.array([s.mol for s in sep_sys.streams], float).flatten()
    # Aggregate duty per unit (heat-utility lists can change length when a
    # unit's design is skipped at zero throughput; per-unit sums keep the
    # vector alignment fixed).
    duties = []
    for u in sorted(sep_sys.units, key=lambda u: u.ID):
        duties.append(sum([hu.duty for hu in u.heat_utilities]))
        duties.append(u.power_utility.rate)
    return mat, np.array(duties, float)

#%% Criteria checker
def check(name, EtOH_mult, IBO_mult,
          rec_EtOH_min=0.90, pur_EtOH_min=0.99,
          rec_IBO_min=0.90, pur_IBO_min=0.95):
    """Simulate the variant thrice and check its criteria. A recovery
    minimum of None means the corresponding product stream must be empty
    (zero-titer feed); purity is only checked on nonempty products."""
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
        # 1e3 kJ/hr (~0.3 kW) normalization floor: sub-kW absolute duty
        # differences on an (almost) all-zero duty vector are noise, not
        # instability (normal-operation duties here are >=1e7 kJ/hr).
        drift_duty = max(drift_duty,
                         np.abs(da - db).max() / max(np.abs(da).max(), 1e3))

    etoh = f.stream.ethanol_product
    ibo = f.stream.isobutanol_product
    checks = []
    for label, product, product_chem, feed_mass, rec_min, pur_min in (
            ('EtOH', etoh, 'Ethanol', feed_EtOH, rec_EtOH_min, pur_EtOH_min),
            ('IBO ', ibo, 'Isobutanol', feed_IBO, rec_IBO_min, pur_IBO_min),
            ):
        if rec_min is None:
            checks.append((f'{label} product empty ', product.F_mass,
                           product.F_mass < EMPTY_F_MASS))
        else:
            rec = product.imass[product_chem] / feed_mass
            checks.append((f'{label} recovery > {rec_min:.2f}', rec,
                           rec > rec_min))
            if product.F_mass > EMPTY_F_MASS:
                pur = product.imass[product_chem] / product.F_mass
                checks.append((f'{label} purity   > {pur_min:.2f}', pur,
                               pur > pur_min))
    checks.extend((
        ('mol drift    < 1e-3', drift_mat, drift_mat < 1e-3),
        ('energy drift < 1e-3', drift_duty, drift_duty < 1e-3),
    ))
    ok = all(c[2] for c in checks)
    for label, val, passed in checks:
        print(f'  {label}: {val:.6g}  [{"ok" if passed else "FAIL"}]')
    print(f'  => {"PASS" if ok else "FAIL"}')
    return ok

#%% Run all feed variants
def main():
    low = dict(rec_min=0.60)  # relaxed recovery floor for very low titers
    variants = [
        # name, EtOH_mult, IBO_mult, criteria overrides
        ('baseline scenario B', 1.0, 1.0, {}),
        ('higher IBO (x2)', 1.0, 2.0, {}),
        ('lower IBO (x0.5)', 1.0, 0.5, {}),
        ('lower EtOH (x0.5)', 0.5, 1.0, {}),
        ('higher EtOH (x1.5)', 1.5, 1.0, {}),
        ('higher IBO, lower EtOH', 0.5, 2.0, {}),
        # ~200 g/L titer equivalents
        ('EtOH 200 g/L (x3.33)', 3.33, 1.0, {}),
        ('IBO 200 g/L (x4.2)', 1.0, 4.2, {}),
        ('both 200 g/L', 3.33, 4.2, {}),
        # low-titer continuity range (reduced IBO recovery acceptable)
        ('low IBO (x0.05)', 1.0, 0.05, dict(rec_IBO_min=low['rec_min'])),
        ('low IBO (x0.01)', 1.0, 0.01, dict(rec_IBO_min=low['rec_min'])),
        ('low IBO (x0.001)', 1.0, 0.001, dict(rec_IBO_min=low['rec_min'])),
        ('low EtOH (x0.01)', 0.01, 1.0, dict(rec_EtOH_min=low['rec_min'])),
        ('low EtOH (x0.001)', 0.001, 1.0, dict(rec_EtOH_min=low['rec_min'])),
        # zero titers: corresponding product stream must be empty
        ('zero IBO', 1.0, 0.0, dict(rec_IBO_min=None)),
        ('zero EtOH', 0.0, 1.0, dict(rec_EtOH_min=None)),
        ('zero both', 0.0, 0.0, dict(rec_EtOH_min=None, rec_IBO_min=None)),
    ]
    results = {}
    for name, em, im, criteria in variants:
        try:
            results[name] = check(name, em, im, **criteria)
        except Exception:
            traceback.print_exc()
            results[name] = False
    print('\n===== SUMMARY =====')
    for k, v in results.items():
        print(f'  {k:<26s} {"PASS" if v else "FAIL"}')
    return int(not all(results.values()))

if __name__ == '__main__':
    sys.exit(main())
