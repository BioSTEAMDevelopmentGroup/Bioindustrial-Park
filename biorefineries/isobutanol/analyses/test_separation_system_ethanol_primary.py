#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Verification of the ethanol-primary separation branch
(``create_EtOH_primary_separation_system``) and the process-gated wrapper
(``create_separation_system``) in separations.py.

Branch (primarily-ethanol feeds: scenario B composition, reduced IBO):

1. >90% recovery of ethanol to the (denatured) ethanol product.
2. >99 wt% ethanol purity of the product on a denaturant-free basis.
3. ALL isobutanol leaves via the bottoms-water outlet (the rectifier
   bottoms; the integrated flowsheet sends it to the WWT mixer): >99% of
   feed IBO there, <0.1% in the ethanol product, <0.1% in the stillage.
   Physical routing basis: dilute IBO (gamma ~ 50, alpha ~ 26 vs. water)
   travels overhead in the beer column; near the ethanol azeotrope in the
   rectifier its activity coefficient collapses, so it is retained in the
   rectifier bottoms, far below the ~10.7 wt% miscibility gap (recovery
   by decantation infeasible -> divert to WWT).
4. Material & energy flows stable across 3 consecutive re-simulations.
5. Zero-alcohol and fully-empty feeds simulate without error with empty
   products (the integrated baseline parks this branch at zero flow).

Wrapper:

6. processes=('IBO_EtOH', 'ethanol') builds the gating splitter S201;
   ONE built system is re-gated in place through split = 1.0 -> 0.0 ->
   0.5 -> 1.0: recoveries hold on the active side, the train that just
   went inactive drains ALL of its outlets to exactly empty (specs/
   guards + recycle-loop drain below min_key_flow), and alcohol mass
   closes across all outlets at every operating point.
7. Single-process calls build without a splitter and trim the outs list.
8. Invalid ``processes`` raise ValueError.

Runs standalone by file path (fresh kernel). Exits non-zero on any
failure.
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

tmo.settings.set_thermo(separations.create_separation_chemicals())

#: Product stream counts as "empty" below this total mass flow [kg/hr];
#: the reference feed carries ~1e5 kg/hr of water.
EMPTY_F_MASS = 1e-2

results = {}

#%% Shared checking machinery

def snapshot(sep_sys):
    mat = np.array([s.mol for s in sep_sys.streams], float).flatten()
    duties = []
    for u in sorted(sep_sys.units, key=lambda u: u.ID):
        duties.append(sum([hu.duty for hu in u.heat_utilities]))
        duties.append(u.power_utility.rate)
    return mat, np.array(duties, float)

def simulate_thrice_and_drift(sep_sys):
    snaps = []
    for _ in range(3):
        sep_sys.simulate()
        snaps.append(snapshot(sep_sys))
    drift_mat = drift_duty = 0.0
    for a, b in ((0, 1), (1, 2), (0, 2)):
        ma, da = snaps[a]
        mb, db = snaps[b]
        drift_mat = max(drift_mat,
                        np.abs(ma - mb).max() / max(np.abs(ma).max(), 1e-9))
        drift_duty = max(drift_duty,
                         np.abs(da - db).max() / max(np.abs(da).max(), 1e3))
    return drift_mat, drift_duty

def report(name, checks):
    ok = all(c[2] for c in checks)
    print(f'\n===== {name} =====', flush=True)
    for label, val, passed in checks:
        print(f'  {label}: {val:.6g}  [{"ok" if passed else "FAIL"}]')
    print(f'  => {"PASS" if ok else "FAIL"}')
    results[name] = ok
    return ok

def record_crash(name):
    traceback.print_exc()
    results[name] = False

#%% Part 1: ethanol-primary branch, standalone

def run_branch_tests():
    bst.main_flowsheet.set_flowsheet('EtOH_primary_branch_test')
    feed = separations.create_scenario_B_feed()
    branch_sys = separations.create_EtOH_primary_separation_system(ins=feed)
    branch_sys.set_tolerance(rmol=1e-5, mol=1e-3, subsystems=True)
    fs = bst.main_flowsheet
    ref_flows = dict(separations.scenario_B_P301_flows)

    def set_feed(EtOH_mult, IBO_mult, empty=False):
        for ID, m in ref_flows.items():
            feed.imass[ID] = 0.0 if empty else m
        feed.imass['Ethanol'] *= EtOH_mult
        feed.imass['Isobutanol'] *= IBO_mult
        feed.T, feed.P = 305.15, 607950.0

    def check_branch(name, EtOH_mult, IBO_mult, rec_EtOH_min=0.90):
        set_feed(EtOH_mult, IBO_mult)
        feed_EtOH = feed.imass['Ethanol']
        feed_IBO = feed.imass['Isobutanol']
        drift_mat, drift_duty = simulate_thrice_and_drift(branch_sys)
        etoh = fs.stream.ethanol_product
        bw = fs.stream.bottoms_water
        stl = fs.stream.stillage
        checks = []
        if feed_EtOH:
            rec = etoh.imass['Ethanol'] / feed_EtOH
            checks.append((f'EtOH recovery > {rec_EtOH_min:.2f}', rec,
                           rec > rec_EtOH_min))
            denat_free = etoh.F_mass - etoh.imass['Octane']
            pur = etoh.imass['Ethanol'] / denat_free
            checks.append(('EtOH purity (denat.-free) > 0.99', pur,
                           pur > 0.99))
        else:
            checks.append(('EtOH product empty', etoh.F_mass,
                           etoh.F_mass < EMPTY_F_MASS))
        if feed_IBO:
            frac_bw = bw.imass['Isobutanol'] / feed_IBO
            frac_pr = etoh.imass['Isobutanol'] / feed_IBO
            frac_st = stl.imass['Isobutanol'] / feed_IBO
            checks.extend((
                ('IBO to bottoms_water > 0.99', frac_bw, frac_bw > 0.99),
                ('IBO in EtOH product < 1e-3', frac_pr, frac_pr < 1e-3),
                ('IBO in stillage < 1e-3', frac_st, frac_st < 1e-3),
            ))
        checks.extend((
            ('mol drift    < 1e-3', drift_mat, drift_mat < 1e-3),
            ('energy drift < 1e-3', drift_duty, drift_duty < 1e-3),
        ))
        return report(f'branch: {name}', checks)

    variants = [
        # name, EtOH_mult, IBO_mult
        ('primarily EtOH, IBO x0.1', 1.0, 0.1),
        ('primarily EtOH, IBO x0.3', 1.0, 0.3),
        ('primarily EtOH, trace IBO x0.01', 1.0, 0.01),
        ('no IBO at all', 1.0, 0.0),
        ('lower EtOH x0.5, IBO x0.1', 0.5, 0.1),
    ]
    for name, em, im in variants:
        try:
            check_branch(name, em, im)
        except Exception:
            record_crash(f'branch: {name}')

    # Zero-alcohol and fully-empty feeds: no crash, empty products. The
    # integrated biorefinery parks this branch at exactly this state.
    for name, empty in (('zero alcohols', False), ('fully empty feed', True)):
        try:
            set_feed(0.0, 0.0, empty=empty)
            branch_sys.simulate()
            etoh = fs.stream.ethanol_product
            bw = fs.stream.bottoms_water
            report(f'branch: {name}', [
                ('EtOH product empty', etoh.F_mass,
                 etoh.F_mass < EMPTY_F_MASS),
                ('bottoms_water IBO-free', bw.imass['Isobutanol'],
                 bw.imass['Isobutanol'] < EMPTY_F_MASS),
            ])
        except Exception:
            record_crash(f'branch: {name}')

#%% Part 2: wrapper with gating splitter

def build_wrapper(flowsheet_ID, processes):
    bst.main_flowsheet.set_flowsheet(flowsheet_ID)
    wfeed = separations.create_scenario_B_feed()
    wsys, udct = separations.create_separation_system(
        ins=wfeed, processes=processes, udct=True)
    wsys.set_tolerance(rmol=1e-5, mol=1e-3, subsystems=True)
    return wsys, udct, wfeed

def check_wrapper_split(wsys, udct, wfeed, split, tag):
    name = f'wrapper both, split={split} ({tag})'
    udct['S201'].split = split
    feed_EtOH = wfeed.imass['Ethanol']
    feed_IBO = wfeed.imass['Isobutanol']
    drift_mat, drift_duty = simulate_thrice_and_drift(wsys)
    fs = bst.main_flowsheet
    e1 = fs.stream.ethanol_product
    i1 = fs.stream.isobutanol_product
    st1 = fs.stream.stillage
    db1 = fs.stream.D103_bottoms
    e2 = fs.stream.ethanol_product_2
    st2 = fs.stream.stillage_2
    bw2 = fs.stream.bottoms_water_2
    outs_all = (e1, i1, st1, db1, e2, st2, bw2)
    EtOH_rec = (e1.imass['Ethanol'] + e2.imass['Ethanol'])/feed_EtOH
    EtOH_closure = sum([s.imass['Ethanol'] for s in outs_all])/feed_EtOH
    IBO_closure = sum([s.imass['Isobutanol'] for s in outs_all])/feed_IBO
    checks = [
        ('total EtOH recovery > 0.90', EtOH_rec, EtOH_rec > 0.90),
        ('EtOH closure within 2%', EtOH_closure,
         abs(1.0 - EtOH_closure) < 0.02),
        ('IBO closure within 2%', IBO_closure,
         abs(1.0 - IBO_closure) < 0.02),
        ('mol drift    < 1e-3', drift_mat, drift_mat < 1e-3),
        ('energy drift < 1e-3', drift_duty, drift_duty < 1e-3),
    ]
    if split == 1.0:
        rec_ibo = i1.imass['Isobutanol']/feed_IBO
        # Every branch-2 outlet must be EXACTLY drained -- after the 0.5
        # visit this verifies the inactive train's specs/guards empty all
        # of its streams (recycle loops included) rather than parking
        # stale flows.
        idle_2 = e2.F_mass + st2.F_mass + bw2.F_mass
        checks.extend((
            ('IBO recovery (branch 1) > 0.90', rec_ibo, rec_ibo > 0.90),
            ('branch-2 outs fully drained', idle_2,
             idle_2 < EMPTY_F_MASS),
        ))
    elif split == 0.0:
        rec_e2 = e2.imass['Ethanol']/feed_EtOH
        frac_bw2 = bw2.imass['Isobutanol']/feed_IBO
        # Branch 1 was active on the previous visit (split=1.0): all four
        # of its outlets must drain to exactly empty.
        idle_1 = e1.F_mass + i1.F_mass + st1.F_mass + db1.F_mass
        checks.extend((
            ('branch-2 EtOH recovery > 0.90', rec_e2, rec_e2 > 0.90),
            ('branch-1 outs fully drained', idle_1,
             idle_1 < EMPTY_F_MASS),
            ('all IBO to bottoms_water_2 > 0.99', frac_bw2,
             frac_bw2 > 0.99),
        ))
    else:
        share = i1.imass['Isobutanol']/(split*feed_IBO)
        checks.append(('branch-1 IBO recovery of its share > 0.90',
                       share, share > 0.90))
    return report(name, checks)

def check_single_process(processes, expected_n_outs, product_ID):
    name = f'wrapper processes={processes}'
    wsys, udct, wfeed = build_wrapper(
        'wrapper_single_' + '_'.join(processes), processes)
    checks = [
        ('no splitter built', 0.0, 'S201' not in udct),
        (f'outs trimmed to {expected_n_outs}', float(len(wsys.outs)),
         len(wsys.outs) == expected_n_outs),
    ]
    wsys.simulate()
    product = getattr(bst.main_flowsheet.stream, product_ID)
    rec = product.imass['Ethanol']/wfeed.imass['Ethanol']
    checks.append(('EtOH recovery > 0.90', rec, rec > 0.90))
    return report(name, checks)

def check_invalid_processes():
    bst.main_flowsheet.set_flowsheet('wrapper_errors')
    ok_empty = ok_bogus = False
    try:
        separations.create_separation_system(
            ins=separations.create_scenario_B_feed('broth_err1'),
            processes=())
    except ValueError:
        ok_empty = True
    except Exception:
        traceback.print_exc()
    try:
        separations.create_separation_system(
            ins=separations.create_scenario_B_feed('broth_err2'),
            processes=('bogus',))
    except ValueError:
        ok_bogus = True
    except Exception:
        traceback.print_exc()
    return report('wrapper invalid processes', [
        ('processes=() raises ValueError', float(ok_empty), ok_empty),
        ("('bogus',) raises ValueError", float(ok_bogus), ok_bogus),
    ])

def run_wrapper_tests():
    # ONE system, re-gated in place (1.0 -> 0.0 -> 0.5 -> 1.0): exercises
    # active -> inactive draining and reactivation of BOTH trains, exactly
    # how the integrated flowsheet is re-gated (sep_udct['S201'].split = x).
    # Stale state carried between operating points is part of what is
    # being tested: every inactive-train unit must leave its outlets empty
    # (its spec/guard empties them; recycle loops drain below
    # min_key_flow, which sits above the flow tolerance, so they reach
    # exactly zero).
    try:
        wsys, udct, wfeed = build_wrapper('wrapper_both',
                                          ('IBO_EtOH', 'ethanol'))
    except Exception:
        record_crash('wrapper both (build)')
    else:
        for split, tag in ((1.0, 'initial'), (0.0, 'drain branch 1'),
                           (0.5, 'shared'), (1.0, 'drain branch 2')):
            try:
                check_wrapper_split(wsys, udct, wfeed, split, tag)
            except Exception:
                record_crash(f'wrapper both, split={split} ({tag})')
    try:
        check_single_process(('IBO_EtOH',), 4, 'ethanol_product')
    except Exception:
        record_crash("wrapper processes=('IBO_EtOH',)")
    try:
        check_single_process(('ethanol',), 3, 'ethanol_product_2')
    except Exception:
        record_crash("wrapper processes=('ethanol',)")
    try:
        check_invalid_processes()
    except Exception:
        record_crash('wrapper invalid processes')

#%% Run everything

def main():
    try:
        run_branch_tests()
    except Exception:
        record_crash('branch tests (setup)')
    try:
        run_wrapper_tests()
    except Exception:
        record_crash('wrapper tests (setup)')
    print('\n===== SUMMARY =====')
    for k, v in results.items():
        print(f'  {k:<44s} {"PASS" if v else "FAIL"}')
    return int(not (results and all(results.values())))

if __name__ == '__main__':
    sys.exit(main())
