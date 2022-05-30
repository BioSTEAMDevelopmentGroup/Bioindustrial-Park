#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

# =============================================================================
# Setup
# =============================================================================

from biosteam import HeatUtility, Facility
from biosteam.units.decorators import cost
from biorefineries.lactic import _facilities as la_facilities
from biorefineries.lactic._facilities import *
from biorefineries.lactic._utils import CEPCI

# The lactic acid biorefinery does not have CWP
__all__ = (*la_facilities.__all__, 'CWP')


# %%

# =============================================================================
# Chilled water package
# =============================================================================

@cost('Duty', 'Chilled water package', units= 'kJ/hr',
      # Original basis is 14 in Gcal/hr
      kW=2535.38, cost=1275750, S=14*4184000, CE=CEPCI[2010], n=0.6, BM=1.6)
class CWP(Facility):
    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 1
    _units= {'Duty': 'kJ/hr'}

    network_priority = 1
    line = 'Chilled water package'

    def __init__(self, ID='', ins=None, outs=()):
        Facility.__init__(self, ID, ins, outs)
        self.agent = HeatUtility.get_cooling_agent('chilled_water')

    def _run(self):
        system_chilled_water_utilities = self.system_chilled_water_utilities = {}

        total_duty = 0
        agent = self.agent
        number = 1
        for u in self.system.units:
            if u is self: continue
            if hasattr(u, 'heat_utilities'):
                for hu in u.heat_utilities:
                    if hu.agent is agent:
                        system_chilled_water_utilities[f'#{number}: {u.ID} - {hu.ID}'] = hu
                        number += 1
                        total_duty -= hu.duty

        hu_chilled = self.heat_utilities[0]
        hu_chilled.mix_from([i for i in system_chilled_water_utilities.values()])
        hu_chilled.reverse()
        self.system_chilled_water_duty = -hu_chilled.duty

        # Total amount of chilled water needed in the whole system
        total_chilled_water = self.total_chilled_water = \
            - hu_chilled.flow * self.chemicals.H2O.MW
        self.ins[0].imass['H2O'] = self.outs[0].imass['H2O'] = total_chilled_water
        self.ins[0].T = hu_chilled.agent.T_limit
        self.outs[0].T = hu_chilled.agent.T

        self.design_results['Duty'] = hu_chilled.duty