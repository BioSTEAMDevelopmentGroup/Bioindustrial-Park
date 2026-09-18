#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Sim-based regression of the DDGS-dryer acid routing
(docs/superpowers/specs/2026-09-12-ddgs-dryer-acid-split-design.md): the
broth acids that reach the DDGS dryer D610 (system.DDGS_DRYER_OVERHEAD_ACIDS,
acetic acid today) are split 1.0 to the dryer exhaust and burned in the
thermal oxidizer X611, like ethanol, instead of being sold as DDGS mass.
One load() (default both-trains build), scenario B (the largest acetate
flow). Checks: (1) D610 splits every listed acid present in the chemical
set, and Ethanol, 1.0 overhead; (2) the DDGS product carries exactly zero
acetic acid; (3) the X611 inlet carries all of the dryer-feed acetic acid;
(4) the WWT feed's acetic acid (the Ev607 vapor share, MX5 -> M501) is
unchanged from the pre-change capture. In the current model that Ev607
vapor share is exactly 0.0 kg/hr -- no acetic acid evaporates in the DDGS
evaporator, so 100 % of the dryer-feed acetic acid was sold as DDGS mass
today and (after this change) 100 % goes to the oxidizer; the WWT feed is
unchanged because it was already zero. Prints the WWT : oxidizer acid split.
Exit 0 + ALL CHECKS PASSED = clean. Fresh kernel."""
from biorefineries import isobutanol
isobutanol.load()
from biorefineries.isobutanol import scenarios, system

#: Acetic acid in the WWT mixer M501 feed [kg/hr] on the scenario-B
#: baseline BEFORE the dryer routing change (captured 2026-09-12 by the
#: implementation plan's Task 0). The Ev607 vapor share of acetic acid,
#: which the dryer split must leave unchanged. It is exactly 0.0 in the
#: current model: no acetic acid evaporates in the DDGS evaporator, so all
#: of it went to the DDGS solids today and none reached WWT.
M501_ACETIC_ACID_PRE_KG_H = 0.0

n_pass = 0
def PASS(m):
    global n_pass; n_pass += 1; print(f'PASS {n_pass}: {m}')

b = scenarios.load_scenario('B')   # runs the baseline model_specification
f = isobutanol.models.models_EtOH_IBO_corn.model.system.flowsheet
D610, X611, M501, DDGS = f.D610, f.X611, f.M501, f.DDGS

# (1) the dryer split table
assert 'AceticAcid' in system.DDGS_DRYER_OVERHEAD_ACIDS, \
    f'AceticAcid missing from DDGS_DRYER_OVERHEAD_ACIDS {system.DDGS_DRYER_OVERHEAD_ACIDS}'
for ID in system.DDGS_DRYER_OVERHEAD_ACIDS:
    if ID in D610.chemicals:
        s = float(D610.isplit[ID])
        assert s == 1.0, f'D610 split of {ID} is {s}, expected 1.0 (to the exhaust)'
assert float(D610.isplit['Ethanol']) == 1.0, "corn's Ethanol=1.0 dryer split was lost"
PASS('D610 splits every listed acid and ethanol 1.0 to the exhaust')

# (2) no acid sold as DDGS
ddgs_acid = float(DDGS.imass['AceticAcid'])
assert ddgs_acid == 0.0, f'DDGS product carries {ddgs_acid} kg/hr acetic acid'
PASS('DDGS product carries exactly zero acetic acid')

# (3) all of the dryer-feed acid reaches the oxidizer
feed_acid = float(D610.ins[0].imass['AceticAcid'])
ox_acid = float(X611.ins[0].imass['AceticAcid'])
assert feed_acid > 0.0, 'dryer feed carries no acetic acid (scenario B should)'
assert abs(ox_acid - feed_acid) <= 1e-9*feed_acid, \
    f'X611 inlet acetic acid {ox_acid} != dryer feed {feed_acid} kg/hr'
PASS(f'X611 inlet carries all {feed_acid:.3f} kg/hr of the dryer-feed acetic acid')

# (4) WWT feed unchanged (nothing upstream of the dryer moved). The pre value
# is 0.0 (no Ev607 vapor share), so this asserts the WWT feed stays ~0 with an
# absolute tolerance robust to the zero baseline.
wwt_acid = float(M501.outs[0].imass['AceticAcid'])
assert abs(wwt_acid - M501_ACETIC_ACID_PRE_KG_H) <= 1e-6 + 1e-4*M501_ACETIC_ACID_PRE_KG_H, \
    f'M501 feed acetic acid {wwt_acid} vs pre-change {M501_ACETIC_ACID_PRE_KG_H} kg/hr'
PASS(f'WWT feed acetic acid unchanged from the pre-change value ({wwt_acid:.6g} kg/hr)')

print(f'WWT : oxidizer acetic acid split = {wwt_acid:.3f} : {ox_acid:.3f} kg/hr '
      f'({100*wwt_acid/(wwt_acid + ox_acid):.1f} % to WWT)')
print(f'ALL {n_pass} CHECKS PASSED')
