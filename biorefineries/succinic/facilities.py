#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat
"""


# %% 

# =============================================================================
# Setup
# =============================================================================

# import biosteam as bst
from biosteam import HeatUtility, Facility
from biosteam.units.decorators import cost
from thermosteam import Stream
from biorefineries.make_a_biorefinery.utils import CEPCI

__all__ = ('CIP', 'ADP', 'CT', 'PWC',)

# !!! Add facilities as needed

# %% 

# =============================================================================
# Clean-in-place system
# =============================================================================
@cost(basis='Flow rate', ID='System', units='kg/hr',
      cost=421000, S=63, CE=CEPCI[2009], n=0.6, BM=1.8)
class CIP(Facility):
    network_priority = 3
    line = 'Clean-in-place system'


# %% 

# =============================================================================
# Air distribution package, size based on stream 950 in Humbird et al.
# =============================================================================
@cost(basis='Flow rate', ID='Plant air compressor', units='kg/hr',
      kW=111.855, cost=28000, S=83333, CE=CEPCI[2010], n=0.6, BM=1.6)
@cost(basis='Flow rate', ID='Plant air reciever', units='kg/hr',
      cost=16000, S=83333, CE=CEPCI[2009], n=0.6, BM=3.1)
@cost(basis='Flow rate', ID='Instrument air dryer', units='kg/hr',
      cost=15000, S=83333, CE=CEPCI[2009], n=0.6, BM=1.8)
class ADP(Facility): 
    network_priority = 3
    line = 'Air distribution package'
    
    def __init__(self, ID='', ins=None, outs=(), ratio=None):
        Facility.__init__(self, ID, ins, outs)
        self.ratio = ratio
    
    def _design(self):
        self.design_results['Flow rate'] = 83333 * self.ratio


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
        chilled_water_utilities = self.chilled_water_utilities = {}
        
        total_duty = 0
        agent = self.agent
        number = 1
        for u in self.system.units:
            if u is self: continue
            if hasattr(u, 'heat_utilities'):
                for hu in u.heat_utilities:
                    if hu.agent is agent:
                        chilled_water_utilities[f'#{number}: {u.ID} - {hu.ID}'] = hu
                        number += 1
                        total_duty -= hu.duty
        hu_chilled = self.heat_utilities[0]
        hu_chilled.mix_from([i for i in chilled_water_utilities.values()])
        if total_duty != 0:
            hu_chilled.reverse() # will trigger an error if hu_chilled is None
            self.ins[0].T = hu_chilled.agent.T_limit
            self.outs[0].T = hu_chilled.agent.T
        self.system_chilled_water_duty = -hu_chilled.duty
        
        # Total amount of chilled water needed in the whole system
        total_chilled_water = self.total_chilled_water = \
            - hu_chilled.flow * self.chemicals.H2O.MW
        self.ins[0].imass['H2O'] = self.outs[0].imass['H2O'] = total_chilled_water

        self.design_results['Duty'] = hu_chilled.duty


# %% 

# =============================================================================
# Cooling tower
# =============================================================================

@cost('Flow rate', 'Cooling tower', units= 'kg/hr',
      kW=559.275, cost=1375000, S=10037820, CE=CEPCI[2010], n=0.6, BM=1.5)
@cost('Flow rate', 'Cooling water pump', units='kg/hr',
      kW=1118.55, cost=283671, S=10982556,  CE=CEPCI[2010], n=0.8, BM=3.1)
class CT(Facility):
    """
    Create a cooling tower with capital cost and power based on the flow rate 
    of cooling water as in [1]_.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    network_priority = 1
    _N_ins = 3
    _N_outs = 2    
    _N_heat_utilities = 1
    _units= {'Flow rate': 'kg/hr'}

    # Page 55 of Humbird et al., including windage
    blowdown = 0.00005+0.0015
    
    def __init__(self, ID='', ins=None, outs=()):
        Facility.__init__(self, ID, ins, outs)
        self.agent = HeatUtility.get_cooling_agent('cooling_water')
        
    def _run(self):
        return_cw, ct_chems, makeup_water = self.ins
        process_cw, blowdown = self.outs
        system_cooling_water_utilities = self.system_cooling_water_utilities = {}

        # Based on stream 945 in Humbird et al.
        return_cw.T = 37 + 273.15
        # Based on streams 940/944 in Humbird et al.
        process_cw.T = blowdown.T = 28 + 273.15
        
        total_duty = 0
        number = 1
        agent = self.agent
        for u in self.system.units:
            if u is self: continue
            if hasattr(u, 'heat_utilities'):
                for hu in u.heat_utilities:
                    if hu.agent is agent:
                        system_cooling_water_utilities[f'#{number}: {u.ID} - {hu.ID}'] = hu
                        number += 1
                        total_duty -= hu.duty
        
        hu_cooling = self.heat_utilities[0]
        hu_cooling.mix_from(system_cooling_water_utilities.values())
        hu_cooling.reverse()        
        self.system_cooling_water_duty = -hu_cooling.duty
        
        # Total amount of cooling water needed in the whole system
        total_cooling_water = self.total_cooling_water = \
            - hu_cooling.flow * self.chemicals.H2O.MW
        return_cw.imass['H2O'] = process_cw.imass['H2O'] = total_cooling_water
        makeup_water.imass['H2O'] = total_cooling_water * self.blowdown
        blowdown.imass['H2O'] = makeup_water.imass['H2O']
        
        # 2 kg/hr from Table 30 on Page 63 of Humbird et al., 4.184 is kcal to kJ,
        # 97.401 MMkcal/hr is the cooling duty on Page 134 of Humbird et al.
        ct_chems.imass['CoolingTowerChems'] = 2 * (hu_cooling.duty/(97.401*4184000))
        self.design_results['Flow rate'] = total_cooling_water


# %% 

# =============================================================================
# Process water center
# =============================================================================

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=250000, S=451555, CE=CEPCI[2009], n=0.7, BM=1.7)
@cost(basis='Flow rate', ID='Circulating pump', units='kg/hr',
      kW=55.9275, cost=15292, S=518924, CE=CEPCI[2010], n=0.8, BM=3.1)
@cost(basis='Flow rate', ID='Makeup water pump', units='kg/hr',
      kW=14.914, cost=6864, S=155564, CE=CEPCI[2010], n=0.8, BM=3.1)
class PWC(Facility):
    _N_ins = 2
    _N_outs = 2
    _units= {'Flow rate': 'kg/hr'}
    
    network_priority = 2
    line = 'Process water center'
    
    def __init__(self, ID='', ins=None, outs=(), process_water_streams=None,
                 recycled_blowdown_streams=None):
        Facility.__init__(self, ID, ins, outs)
        self.process_water_streams = process_water_streams
        self.recycled_blowdown_streams = recycled_blowdown_streams

    def _run(self):
        makeup, RO_water = self.ins
        process_water, discharged = self.outs
        
        water_demand = sum(i.imol['Water'] for i in self.process_water_streams)
        water_needs = water_demand - RO_water.imol['Water']
        self.recycled_water = RO_water.imass['Water']
        
        if self.recycled_blowdown_streams:
            water_needs -= sum(i.imol['Water'] for i in self.recycled_blowdown_streams)
            self.recycled_water += sum(i.imass['Water'] for i in self.recycled_blowdown_streams)
        
        if water_needs > 0:
            makeup.imol['Water'] = water_needs
            discharged.empty()
        else:
            discharged.imol['Water'] = - water_needs
            makeup.empty()

        process_water.mol = makeup.mol + RO_water.mol - discharged.mol

        self.design_results['Flow rate'] = self.F_mass_in

