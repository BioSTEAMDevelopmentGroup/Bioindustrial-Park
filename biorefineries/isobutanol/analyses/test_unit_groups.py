#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Sim-based regression of the unit-group restructure
(docs/superpowers/specs/2026-09-19-facility-unit-groups-design.md): the
boiler / turbogenerator / cooling-utility groups, the four cost
pseudo-groups, the 'Material cost' -> 'Operating cost' rename and the shared
'alcohol recovery' + two purification groups. One load() (default
both-trains build), scenario B (isobutanol flows; dryer and turbine live).
Checks 1-9 are the spec's section 7; 10-11 pin the section-5 reconciliation
rule (each priced feed in exactly one group; the boiler row + the steam-fuel
pseudo-group close on BT801). Every check runs even if an earlier one fails,
so a single run names everything that is wrong. Prints the section-5 probe
numbers and the Operating-cost column. Exit 0 + ALL 11 CHECKS PASSED =
clean. Fresh kernel."""
import os
import sys
import threading
import numpy as np
import biosteam as bst
from biorefineries import isobutanol
isobutanol.load()
from biorefineries.isobutanol import scenarios, system

# A hang (a non-converging flowsheet) must FAIL the test, not block it.
_watchdog = threading.Timer(
    600.0, lambda: (print('\nWATCHDOG: 600 s exceeded -> exit 2', flush=True),
                    os._exit(2)))
_watchdog.daemon = True
_watchdog.start()

EXPECTED_GROUPS = [
    'feedstock acquisition', 'feedstock saccharification',
    'sugar solution preparation', 'fermentation',
    'alcohol recovery', 'ethanol purification', 'isobutanol purification',
    'storage and handling', 'DDGS recovery', 'wastewater treatment',
    'heat exchanger network',
    'boiler', 'turbogenerator', 'cooling utility facilities',
    'other facilities',
    'natural gas (for steam generation)', 'natural gas (for product drying)',
    'fixed operating cost', 'excess electricity',
    ]
EXPECTED_METRICS = ['Installed equipment cost', 'Cooling duty', 'Heating duty',
                    'Electricity consumption', 'Operating cost']
INST, COOL, HEAT, ELEC, OPER = EXPECTED_METRICS

n_pass = 0
failures = []
def CHECK(label, fn):
    """Run one check; record (never raise) a failure so every check runs."""
    global n_pass
    try:
        detail = fn()
    except Exception as e:
        failures.append(label)
        print(f'FAIL: {label}\n      {type(e).__name__}: {e}', flush=True)
    else:
        n_pass += 1
        print(f'PASS {n_pass}: {label}' + (f' ({detail})' if detail else ''),
              flush=True)

def close(a, b, rel):
    return abs(a - b) <= rel*max(abs(a), abs(b), 1e-300)

b = scenarios.load_scenario('B')   # runs the baseline model_specification
sys_ = system.corn_EtOH_IBO_sys
tea = system.corn_EtOH_IBO_sys_tea
groups = system.unit_groups
ugd = system.unit_groups_dict
sep_udct = system.sep_udct
u = sys_.flowsheet.unit
BT, M510 = u.BT801, system.M510

def M(group_name, metric_name):
    """Value of the named metric of the named group (by NAME, never index)."""
    for m in ugd[group_name].metrics:
        if m.name == metric_name: return m()
    raise KeyError(f'group {group_name!r} has no metric {metric_name!r}')

def stream_utility_cost(unit):
    """USD/hr of a unit's define_utility stream cash flows (fuel, ash
    disposal, RO water): utility cost less heat- and power-utility costs."""
    return ((unit.utility_cost or 0.)
            - sum([i.cost for i in unit.heat_utilities])
            - unit.power_utility.cost)

def traced_feeds(group):
    """The feed set UnitGroup.get_material_cost sums over (same tracing)."""
    inlets = set(bst.utils.feeds_from_units(group.units))
    bst.utils.filter_out_missing_streams(inlets)
    if group.extend_feed_ends:
        inlets = [bst.utils.get_inlet_origin(i) for i in inlets]
    return set(bst.utils.feeds(inlets))

# (1) installed cost closes on the TEA
def check_1():
    total = sum([M(g.name, INST) for g in groups])
    ref = tea.installed_equipment_cost/1e6
    assert close(total, ref, 1e-9), f'groups {total!r} != TEA {ref!r} MM$'
    return f'{total:.4f} MM$'
CHECK('sum of the groups\' installed cost == tea.installed_equipment_cost', check_1)

# (2) boiler + turbogenerator == BT801 + M510
def check_2():
    bo, tg = M('boiler', INST), M('turbogenerator', INST)
    ref = (BT.installed_cost + M510.installed_cost)/1e6
    assert bo > 0. and tg > 0., f'boiler {bo!r}, turbogenerator {tg!r} must both be > 0'
    assert close(bo + tg, ref, 1e-9), f'{bo!r} + {tg!r} != BT801 + M510 {ref!r} MM$'
    return f'boiler {bo:.4f} + turbogenerator {tg:.4f} MM$'
CHECK('boiler + turbogenerator installed cost == BT801 + M510, both > 0', check_2)

# (3) electricity consumption
def check_3():
    bo, tg = M('boiler', ELEC), M('turbogenerator', ELEC)
    ref = (BT.power_utility.consumption + M510.power_utility.consumption)/1e3
    assert tg == 0., f'turbogenerator electricity consumption {tg!r} != 0'
    assert close(bo + tg, ref, 1e-9), f'{bo!r} + {tg!r} != BT801 consumption {ref!r} MW'
    return f'boiler {bo:.4f} MW'
CHECK('boiler + turbogenerator electricity == BT801 consumption; turbogenerator 0', check_3)

# (4) cooling duty
def check_4():
    bo, tg = M('boiler', COOL), M('turbogenerator', COOL)
    ref = abs(BT.cooling_duty)/1e6
    assert bo == 0., f'boiler cooling duty {bo!r} != 0'
    assert close(tg, ref, 1e-12), f'turbogenerator cooling duty {tg!r} != {ref!r} GJ/hr'
    return f'turbogenerator {tg:.4f} GJ/hr'
CHECK('boiler cooling duty 0; turbogenerator == abs(BT801.cooling_duty)/1e6', check_4)

# (5) the section-5 rule: the Operating-cost column closes on tea.AOC
def check_5():
    column = sum([M(g.name, OPER) for g in groups])
    ref = tea.AOC/tea.operating_hours
    if not close(column, ref, 1e-6):
        cu = list(sys_.cost_units)
        print('      --- AOC decomposition [USD/hr] ---')
        print(f'      FOC/h              {tea.FOC/tea.operating_hours:.6f}')
        print(f'      priced feeds       {sum([s.cost for s in sys_.feeds if s.price]):.6f}')
        print(f'      inlet fees         {sum([i._inlet_cost for i in cu]):.6f}')
        print(f'      heat utilities     {sum([sum([h.cost for h in i.heat_utilities]) for i in cu]):.6f}')
        print(f'      power (net)        {sum([i.power_utility.cost for i in cu]):.6f}')
        for i in cu:
            v = stream_utility_cost(i)
            if abs(v) > 1e-9: print(f'      stream utility {i.ID:<8} {v:.6f}')
        print('      --- Operating-cost column [USD/hr] ---')
        for g in groups: print(f'      {g.name:<38} {M(g.name, OPER):.6f}')
    assert close(column, ref, 1e-6), \
        f'column {column!r} != tea.AOC/operating_hours {ref!r} USD/hr (diff {column - ref!r})'
    return f'{column*tea.operating_hours/1e6:.4f} MM$/yr'
CHECK('sum of Operating cost x operating_hours == tea.AOC', check_5)

# (6) exclusive, exhaustive membership
def check_6():
    bad = []
    for unit in sys_.units:
        homes = [g.name for g in groups if any(unit is j for j in g.units)]
        if len(homes) != 1: bad.append((unit.ID, homes))
    assert not bad, f'units not in exactly one group: {bad}'
    return f'{len(sys_.units)} units'
CHECK('every in-system unit is in exactly one group', check_6)

# (7) names and metric shape
def check_7():
    names = [g.name for g in groups]
    assert names == EXPECTED_GROUPS, f'group names {names}'
    assert list(ugd) == EXPECTED_GROUPS, f'unit_groups_dict keys {list(ugd)}'
    for g in groups:
        mnames = [m.name for m in g.metrics]
        assert mnames == EXPECTED_METRICS, f'{g.name!r} metrics {mnames}'
    for stale in ('ethanol separation', 'isobutanol separation'):
        assert stale not in ugd, f'stale group {stale!r} still registered'
    return f'{len(names)} groups x {len(EXPECTED_METRICS)} metrics'
CHECK('the 19 expected groups, five metrics each ending in Operating cost, nothing stale', check_7)

# (8) separation membership, by object identity
def ids(units): return sorted([i.ID for i in units])
def check_8():
    got = set(ugd['alcohol recovery'].units)
    # system.P301 is system.py's broth pump feeding MX8 -- NOT
    # sep_udct['P301'], the ethanol-primary train's own beer pump.
    expected = {system.P301, system.MX8, u.V409, u.P410, system.MX7, u.PX, u.MX5,
                u.P406, u.P407, u.MX,
                sep_udct['S201'], sep_udct['D101'], sep_udct['M201'], sep_udct['D102']}
    assert got == expected, (f"'alcohol recovery' strays {ids(got - expected)}, "
                             f"missing {ids(expected - got)}")
    got = set(ugd['isobutanol purification'].units)
    expected = {sep_udct[i] for i in ('D103', 'M301', 'H301', 'S301', 'D104', 'H302')}
    assert got == expected, (f"'isobutanol purification' strays {ids(got - expected)}, "
                             f"missing {ids(expected - got)}")
    got = set(ugd['ethanol purification'].units)
    required = ({sep_udct[i] for i in ('H202', 'MS201', 'H201',
                                        'D302', 'D303', 'U301', 'H304')}
                | {system.MX6, u.P512})
    assert required <= got, f"'ethanol purification' missing {ids(required - got)}"
    assert sep_udct['P301'] in got, "ethanol-primary beer pump sep_udct['P301'] not in 'ethanol purification'"
    return f"ethanol purification = {ids(got)}"
CHECK('separation-group membership matches spec section 4a', check_8)

# (9) scenario B: all three separation groups live; the models metric evaluates
def check_9():
    for name in ('alcohol recovery', 'ethanol purification', 'isobutanol purification'):
        v = M(name, INST)
        assert v > 0., f'{name!r} installed cost {v!r} not > 0'
    h = M('isobutanol purification', HEAT)
    assert h > 0., f"'isobutanol purification' heating duty {h!r} not > 0"
    model = isobutanol.models.models_EtOH_IBO_corn.model
    mnames = [m.name for m in model.metrics]
    assert 'Ethanol separation operating cost' not in mnames, 'stale models metric still present'
    hits = [m for m in model.metrics if m.name == 'Separation operating cost']
    assert len(hits) == 1, f"{len(hits)} 'Separation operating cost' metrics"
    assert hits[0].element == 'Separation', f'element {hits[0].element!r}'
    v = hits[0]()
    assert np.isfinite(v) and v > 0., f"'Separation operating cost' = {v!r}"
    return f'Separation operating cost {v:.5f} $/kg total alcohol'
CHECK('scenario B: three live separation groups + the Separation operating cost metric', check_9)

# (10) each priced feed is counted by exactly one group
def check_10():
    feeds_by_group = {g.name: traced_feeds(g) for g in groups}
    bad = []
    for s in sys_.feeds:
        if not s.price: continue
        homes = [name for name, fs in feeds_by_group.items() if s in fs]
        if len(homes) != 1: bad.append((s.ID, homes))
    assert not bad, f'priced feeds not counted exactly once: {bad}'
    return f'{len([s for s in sys_.feeds if s.price])} priced feeds'
CHECK('every priced feed is traced by exactly one group', check_10)

# (11) the boiler row + the steam-fuel pseudo-group close on BT801
def check_11():
    row = M('boiler', OPER) + M('natural gas (for steam generation)', OPER)
    ref = ugd['boiler'].get_material_cost() + stream_utility_cost(BT)
    assert close(row, ref, 1e-9), f'boiler + steam fuel {row!r} != native + BT801 stream utilities {ref!r}'
    assert M('natural gas (for steam generation)', OPER) > 0., 'no boiler fuel cost in scenario B'
    assert M('natural gas (for product drying)', OPER) > 0., 'no dryer fuel cost in scenario B'
    return f'{row:.4f} USD/hr'
CHECK('boiler + natural gas (for steam generation) == BT801 materials + fuel + ash', check_11)

# ---- section-5 probe numbers (reported to the user; informational) ----
print('\n--- section-5 probes ---')
ng, ash = BT.natural_gas, BT.ash_disposal
print(f'BT801 natural gas: stream price {ng.price!r}, stream cost {ng.cost!r}, '
      f'utility price {BT.natural_gas_price!r}, F_mass {ng.F_mass:.3f} kg/hr')
print(f'BT801 ash disposal: utility price {BT.ash_disposal_price!r}, F_mass {ash.F_mass:.3f} kg/hr')
for i in sys_.units:
    if isinstance(i, bst.DrumDryer):
        print(f'{i.ID} natural gas: stream price {i.natural_gas.price!r}, '
              f'stream-utility cost {stream_utility_cost(i):.4f} USD/hr')
chems = u.CT901.ins[1]
print(f'CT901 chemicals {chems.ID}: price {chems.price!r}, cost {chems.cost:.4f} USD/hr')
print(f'PWC901 stream-utility cost (RO-water charge): {stream_utility_cost(u.PWC901):.4f} USD/hr')
if not failures or OPER in [m.name for m in groups[0].metrics]:
    print('\n--- Operating cost by group ---')
    for g in groups:
        try: v = M(g.name, OPER)
        except KeyError: continue
        print(f'{g.name:<38} {v:12.4f} USD/hr  {v*tea.operating_hours/1e6:9.4f} MM$/yr')

_watchdog.cancel()
if failures:
    print(f'\n{len(failures)} CHECK(S) FAILED: {failures}')
    sys.exit(1)
print(f'ALL {n_pass} CHECKS PASSED')
