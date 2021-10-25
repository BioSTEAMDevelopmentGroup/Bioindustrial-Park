# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

__all__ = ('SLLECentrifuge', 'SolidLiquidsSplitCentrifuge',)

import biosteam as bst
from thermosteam import separations as sep

@bst.units.decorators.cost('Flow rate', units='m^3/hr',
    CE=525.4, cost=28100, n=0.574, kW=1.4, ub=400, BM=2.03,
    N='Number of centrifuges')
class SLLECentrifuge(bst.Unit):
    _N_ins = 1
    _N_outs = 3
    _N_heat_utilities = 0
    
    @property
    def solids_split(self):
        return self._solids_isplit.data
    @property
    def solids_isplit(self):
        return self._solids_isplit
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 solids_split, top_chemical=None, efficiency=1.0,
                 moisture_content=0.5):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        
        # [ChemicalIndexer] Splits to 0th outlet stream.
        self._solids_isplit = self.thermo.chemicals.isplit(solids_split)
        
        #: [str] Identifier of chemical that will be favored in the "liquid" phase.
        #: If none given, the "liquid" phase will the lightest and the "LIQUID"
        #: phase will be the heaviest.
        self.top_chemical = top_chemical
        
        #: Fraction of feed in liquid-liquid equilibrium.
        #: The rest of the feed is divided equally between phases.
        self.efficiency = efficiency 
        
        #: Moisture content of retentate
        self.moisture_content = moisture_content
        assert self._solids_isplit['7732-18-5'] == 0, 'cannot define water split, only moisture content'

    def _run(self):
        top, bottom, solids = self.outs
        sep.lle(self.ins[0], top, bottom, self.top_chemical, self.efficiency)
        bottom.split_to(solids, bottom, self.solids_split)
        sep.adjust_moisture_content(solids, bottom, self.moisture_content)


@bst.units.decorators.cost('Flow rate', units='m^3/hr',
    CE=525.4, cost=28100, n=0.574, kW=1.4, ub=200, BM=2.03,
    N='Number of centrifuges')
class SolidLiquidsSplitCentrifuge(bst.Unit):
    
    _N_ins = 1
    _N_outs = 3
    _N_heat_utilities = 0
    
    @property
    def solids_split(self):
        return self._solids_isplit.data
    @property
    def solids_isplit(self):
        return self._solids_isplit
    
    @property
    def liquids_split(self):
        return self._liquids_isplit.data
    @property
    def liquids_isplit(self):
        return self._liquids_isplit
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 liquids_split, solids_split, 
                 top_chemical=None,
                 moisture_content=0.5):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        
        # [ChemicalIndexer] Splits to 2th outlet stream (solids) from heavier phase.
        self._solids_isplit = self.thermo.chemicals.isplit(solids_split)
        
        # [ChemicalIndexer] Splits to 0th outlet stream (lighter phase).
        self._liquids_isplit = self.thermo.chemicals.isplit(liquids_split)
        
        #: [str] Identifier of chemical that will be favored in the "liquid" phase.
        #: If none given, the "liquid" phase will the lightest and the "LIQUID"
        #: phase will be the heaviest.
        self.top_chemical = top_chemical
        
        #: Moisture content of retentate
        self.moisture_content = moisture_content
        assert self._solids_isplit['7732-18-5'] == 0, 'cannot define water split to solids, only moisture content'

    def _run(self):
        top, bottom, solids = self.outs
        self.ins[0].split_to(top, bottom, self.liquids_split)
        bottom.split_to(solids, bottom, self.solids_split)
        sep.adjust_moisture_content(solids, bottom, self.moisture_content)