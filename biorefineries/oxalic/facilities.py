#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Oxalic acid biorefineries.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for 2,3-Butanediol instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

All units are explicitly defined here for transparency and easy reference

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
from biorefineries.HP.utils import CEPCI

__all__ = ('CIP', 'ADP', 'CT', 'PWC', 'BT')


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


# %% 

# =============================================================================
# Boiler and turbogenerator
# =============================================================================

@cost(basis='Flow rate', ID='Boiler', units='kg/hr',
      kW=2752, cost=28550000, S=238686, CE=CEPCI[2010], n=0.6, BM=1.8)
@cost(basis='Work', ID='Turbogenerator', units='kW',
      kW=23.6, cost=9500000, S=42200, CE=CEPCI[2010], n=0.6, BM=1.8)
@cost(basis='Flow rate', ID='Hot process water softener system', units='kg/hr',
      cost=78000, S=235803, CE=CEPCI[2010], n=0.6, BM=1.8)
@cost(basis='Flow rate', ID='Deaerator', units='kg/hr',
      cost=305000, S=235803, CE=CEPCI[2010], n=0.6, BM=3)
@cost(basis='Flow rate', ID='Amine addition pkg', units='kg/hr',
      cost=40000, S=235803, CE=CEPCI[2010], n=0, BM=1.8)
class BT(Facility):
    """
    Create a BoilerTurbogenerator object that will calculate electricity
    generation from burning the feed. It also takes into account how much
    steam is being produced, and the required cooling utility of the turbo
    generator. Combustion reactions are populated based on molecular formula
    of the ins. Purchase natural gas if BT cannot meet system heating/steam demand.
    
    Parameters
    ----------
    ins : 
        [0] Liquid/solid wastes to be burned.
        [1] Gas waste to be burned.        
        [2] FGD lime.        
        [3] Boiler chemicals.        
        [4] Baghouse bag.                
        [5] Supplementary natural gas in the case that waste energy is not enough
            to meet system steam demand
        [6] Makeup water.

    outs :
        [0] Gas emission.        
        [1] Ash residues.        
        [2] Blowdown water.
   
    boiler_efficiency : float
        Fraction of heat transfered to steam.        
    turbo_generator_efficiency : float
        Fraction of steam heat converted to electricity.        
    agent : UtilityAgent
        Steam produced.
        
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
        
    .. [2] Cortes-Peña, Y.; Kumar, D.; Singh, V.; Guest, J. S. 
        BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, 
        and Techno-Economic Analysis of Biorefineries under Uncertainty. 
        ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
        https://doi.org/10.1021/acssuschemeng.9b07040.
        
    """
    network_priority = 0
    _N_ins = 7
    _N_outs = 3
    _units= {'Flow rate': 'kg/hr',
             'Work': 'kW'} 
    
    blowdown = 0.03
    
    def __init__(self, ID='', ins=None, outs=(), *, B_eff=0.8,
                 TG_eff=0.85, combustibles, side_streams_to_heat=()):
        Facility.__init__(self, ID, ins, outs)
        self.B_eff = B_eff
        self.TG_eff = TG_eff
        self.combustibles = combustibles
        self.side_streams_to_heat = side_streams_to_heat
        self.side_streams_lps = None

    def _run(self): pass

    def _design(self):
        feed_solids, feed_gases, lime, boiler_chems, bag, natural_gas, makeup_water \
            = self.ins
        emission, ash, blowdown_water = self.outs
        side_streams_to_heat = self.side_streams_to_heat
        side_streams_lps = self.side_streams_lps
        system_heating_utilities = self.system_heating_utilities = {}
        lps = HeatUtility.get_heating_agent('low_pressure_steam')
        hps = HeatUtility.get_heating_agent('high_pressure_steam')
        
        # Use combustion reactions to create outs
        combustion_rxns = self.chemicals.get_combustion_reactions()
        combustable_feeds = Stream(None)
        emission.mol = feed_solids.mol + feed_gases.mol
        combustable_feeds.copy_flow(emission, tuple(self.combustibles), remove=True)
        combustion_rxns.force_reaction(combustable_feeds.mol)
        emission.mol += combustable_feeds.mol

        # FGD lime scaled based on SO2 generated,
        # 20% stoichiometetric excess based on P52 of Humbird et al.
        # 895 kg/hr based on Table 30 in Humbird et al.
        lime.imol['Lime'] = emission.imol['SO2'] * 1.2
        
        # 92% of SO2 removed by lime
        CaSO4_mol = emission.imol['SO2'] * 0.92
        emission.imol['SO2'] *= (1 - 0.92)

        # Air/O2 usage not rigorously modeled
        emission.imol['O2'] = 0
        
        ash.empty()
        for chemical in emission.chemicals:
            if chemical.ID not in ('Water', 'H2O') and chemical.locked_state != 'g':
                ash.imol[chemical.ID] = emission.imol[chemical.ID]
                emission.imol[chemical.ID] = 0
        emission.imol['Water'] = feed_solids.imol['Water'] + feed_gases.imol['Water']
                
        ash.mol += boiler_chems.mol
        ash.imol['CaSO4'] = CaSO4_mol
        ash.imol['Lime'] += lime.F_mol - CaSO4_mol
        
        emission.phase = 'g'
        ash.phase = 's'
        # Assume T of emission and ash are the same as hps, which
        # has highest T amont all heating agents
        emission.T = ash.T = hps.T

        # Total heat generated by the boiler (kJ/hr)
        H_in = feed_solids.H + feed_gases.H
        H_out = emission.H + ash.H
        # Water evaporation energy is already accounted for in emission enthalpy
        heat_from_combustion = -(feed_solids.HHV+feed_gases.HHV)
        heat_generated = self.heat_generated = \
            (H_in+heat_from_combustion)*self.B_eff - H_out
        
        number = 1
        for u in self.system.units:
            if u is self: continue
            if hasattr(u, 'heat_utilities'):
                for hu in u.heat_utilities:
                    # Including low/medium/high_pressure_steam
                    if hu.flow*hu.duty > 0:
                        system_heating_utilities[f'#{number}: {u.ID} - {hu.ID}'] = hu
                        number += 1
        
        # Use lps to account for the energy needed for the side steam
        if side_streams_to_heat:
            if not side_streams_lps:
                side_streams_lps = self.side_streams_lps = HeatUtility()
                side_streams_lps.load_agent(lps)
            side_streams_lps(duty=sum([i.H for i in side_streams_to_heat]), 
                             T_in=298.15)
            system_heating_utilities['BT - side_streams_lps'] = side_streams_lps

        system_heating_demand = self.system_heating_demand = \
            sum([i.duty for i in system_heating_utilities.values()])
            
        BT_heat_surplus = self.BT_heat_surplus = heat_generated - system_heating_demand
        self.system_steam_demand = sum([i.flow for i in system_heating_utilities.values()])
        
        hu_cooling = HeatUtility()
        
        # BT can meet system heating/steam demand
        if BT_heat_surplus >0:
            # 3600 is conversion of kJ/hr to kW (kJ/s)
            electricity_generated = self.electricity_generated = \
                BT_heat_surplus * self.TG_eff / 3600
            
            # Take the opposite for cooling duty (i.e., cooling duty should be negative)
            # this is to condense the unused steam
            cooling_need = self.cooling_need = -(BT_heat_surplus-electricity_generated*3600)
            # hu_cooling = HeatUtility()
            hu_cooling(duty=cooling_need, T_in=lps.T)
            natural_gas.empty()
            # hu_BT.add(hu_cooling)
            
        # BT cannot meet system heating/steam demand, supplement with natural gas
        else:
            CH4_LHV = natural_gas.chemicals.CH4.LHV
            natural_gas.imol['CH4'] = BT_heat_surplus / (CH4_LHV*self.B_eff)
            emission.imol['CO2'] += natural_gas.imol['CH4']
            emission.imol['H2O'] += 2 * natural_gas.imol['CH4']
            electricity_generated = self.electricity_generated = 0
            
        heating_utilities = HeatUtility.sum_by_agent(system_heating_utilities.values())
        for i in heating_utilities:
            i.reverse()
        
        if hu_cooling.duty != 0:
            self.heat_utilities = tuple([hu_cooling, *heating_utilities])
        else:
            self.heat_utilities = tuple(heating_utilities)
        
        total_steam_mol = sum([i.flow for i in system_heating_utilities.values()])
        total_steam = total_steam_mol * self.chemicals.H2O.MW
        blowdown_water.imass['H2O'] = total_steam * self.blowdown
        blowdown_water.T = 373.15

        # Additional need for making lime slurry
        makeup_water.imol['H2O'] = blowdown_water.imol['H2O'] + lime.F_mol/0.2*0.8
        
        # 1.23 is $2007/hour and 2.2661 is $/lb from Table 30 in Humbird et al.,
        # 144.629 is the total duty from Page 131 in Humbird et al.,
        # 2.20462 is kg to lb, 4.184 is kcal to kJ, 
        ratio = system_heating_demand/(144.629*1e6*4.184)
        boiler_chems.imass['BoilerChems'] = 1.23/2.2661/2.20462 * ratio
        bag.imass['BaghouseBag'] = ratio

        Design = self.design_results        
        Design['Flow rate'] = total_steam
        Design['Work'] = electricity_generated

    def _cost(self):
        # If no kW setting in decorated cost, then power_utility won't be reset
        self.power_utility(0)
        self._decorated_cost()
        self.electricity_CHP_used = self.power_utility.rate
        self.power_utility.production = self.electricity_generated

