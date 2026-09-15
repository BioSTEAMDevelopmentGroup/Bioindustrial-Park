#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Sim-based regression of the kinetic optimization's V406-effluent
(g/L-BROTH) and acetate metrics -- kinetic_optimization._broth_conc and the
five TRACKED_METRICS columns added 2026-09-14 -- read exactly as the
kinetic BO reads them: ko.get_handles() on the live flowsheet, every
TRACKED_METRICS getter evaluated as the engine's per-trial recording loop
does. One load(), scenario B (nonzero ethanol, isobutanol, cells AND
acetate in the broth; a batch, so no fed-batch volume rescaling). Checks:
every getter evaluates without raising and the five new metrics are
finite and positive; the stream the getters read is the broth effluent
(P406's feed, not the vent) and each broth value IS imass/F_vol of it;
and the physics of the two bases: the fermentation CO2 strips a
species-specific fraction of each volatile into the vent (outs[0]) BEFORE
the broth leaves, so broth/water = (water-volume fraction) x (1 - vent
fraction) -- dividing the vent fraction out (mass conservation over the two
outlets) must collapse the three product ratios onto ONE water-volume
fraction in (0.5, 1] (no inlet carries ethanol, isobutanol or acetate),
with the cell ratio within 1 % of it. Exit 0 + ALL CHECKS PASSED = clean."""
import math
import os
import threading
import numpy as np
from biorefineries import isobutanol
isobutanol.load(separation_processes=('IBO_EtOH',))
from biorefineries.isobutanol import scenarios, system
from biorefineries.isobutanol import kinetic_optimization as ko

# A hang (a non-converging flowsheet) must FAIL the test, not block it.
_watchdog = threading.Timer(
    600.0, lambda: (print('\nWATCHDOG: 600 s exceeded -> exit 2', flush=True),
                    os._exit(2)))
_watchdog.daemon = True
_watchdog.start()

n_pass = 0
def PASS(m):
    global n_pass; n_pass += 1; print(f'PASS {n_pass}: {m}')

NEW = ('EtOH titer (broth)', 'IBO titer (broth)', 'Cell density (broth)',
       'Acetate titer', 'Acetate titer (broth)')
#: chemical -> (kinetic g/L-water metric, effluent g/L-broth metric)
PAIRS = {'Ethanol':    ('EtOH titer',    'EtOH titer (broth)'),
         'Isobutanol': ('IBO titer',     'IBO titer (broth)'),
         'AceticAcid': ('Acetate titer', 'Acetate titer (broth)'),
         'Yeast':      ('Cell density',  'Cell density (broth)')}

b = scenarios.load_scenario('B')
b['model_specification'](**b['feeding_kwargs'])          # the B baseline run
handles = ko.get_handles()                                # the BO's own handles
handles['latest_TEA_solution'].update(
    handles['solve_TEA'](stream_IDs=('ethanol', 'isobutanol')))

# 1. The engine's per-trial recording loop, on the live flowsheet.
record = {name: getter(handles) for name, getter in ko.TRACKED_METRICS.items()}
for name in NEW:
    v = float(record[name])
    assert math.isfinite(v) and v > 0.0, (name, v)
print('  ' + '\n  '.join(f'{k}: {float(record[k]):.6g}' for k in
                          [*sum(PAIRS.values(), ()), 'tau', 'n_glu_spikes']))
PASS('every TRACKED_METRICS getter evaluates on the live handles; the five '
     'new metrics are finite and positive on the scenario-B broth')

# 2. The stream the getters read is the broth effluent (P406's feed), and
#    each broth value IS imass/F_vol of it (kg/hr over m3/hr = g/L-broth).
V406 = handles['V406']
effluent = ko._effluent(handles)
assert effluent is V406.outs[1] and effluent is system.f.P406.ins[0]
assert effluent is not V406.outs[0]                       # not the vent
for chem, (_, broth_name) in PAIRS.items():
    direct = effluent.imass[chem]/effluent.F_vol
    assert np.isclose(record[broth_name], direct, rtol=1e-12, atol=0.0), (
        chem, record[broth_name], direct)
PASS('the getters read V406.outs[1] (the broth, P406 feed; not the vent) and '
     'each broth metric equals imass/F_vol of that stream')

# 3. Physical self-consistency of the two bases. NSKBatchReactor builds the
#    effluent as imass = [s]_water x V_water for every mapped species, then
#    the fermentation CO2 strips a species-specific fraction of each
#    volatile into the vent (outs[0]; scenario B ~1.2 % of the ethanol and
#    ~1.9 % of the isobutanol -- what V409 recovers) BEFORE the broth
#    leaves as outs[1]. So broth/water = (V_water/F_vol_broth) x (1 -
#    f_vent), f_vent = vent/(vent + broth) of that species by mass
#    conservation over the two outlets: dividing the vent factor out must
#    collapse the three product ratios onto ONE water-volume fraction in
#    (0.5, 1] (no inlet carries ethanol, isobutanol or acetate; the broth
#    is the water plus the solutes). Cells also receive the inoculum feed,
#    so their corrected ratio is pinned only within 1 % of it.
vent = V406.outs[0]
f_vent = {chem: vent.imass[chem]/(vent.imass[chem] + effluent.imass[chem])
          for chem in PAIRS}
ratios = {chem: float(record[broth])/float(record[water])
          for chem, (water, broth) in PAIRS.items()}
corrected = {chem: ratios[chem]/(1.0 - f_vent[chem]) for chem in PAIRS}
print('  broth/water:    ' + ', '.join(f'{c} {r:.6f}' for c, r in ratios.items()))
print('  vent fraction:  ' + ', '.join(f'{c} {f:.4%}' for c, f in f_vent.items()))
print('  vent-corrected: ' + ', '.join(f'{c} {r:.6f}' for c, r in corrected.items()))
assert 0.0 < f_vent['Ethanol'] < 0.05 and 0.0 < f_vent['Isobutanol'] < 0.05, f_vent
products = [corrected[c] for c in ('Ethanol', 'Isobutanol', 'AceticAcid')]
assert all(0.5 < r <= 1.0 + 1e-9 for r in products), products
assert np.ptp(products) < 1e-6*max(products), products
assert abs(corrected['Yeast'] - products[0]) < 1e-2*products[0], (
    corrected['Yeast'], products[0])
PASS('broth/water = one water-volume fraction in (0.5, 1] x (1 - vent '
     'fraction) for ethanol, isobutanol and acetate; cells within 1 % of it')

print(f'\nALL {n_pass} CHECKS PASSED')
