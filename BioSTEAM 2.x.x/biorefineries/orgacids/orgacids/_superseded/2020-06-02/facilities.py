#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 09:53:03 2020

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

@author: yalinli_cabbi
"""


# %% Setup

from biosteam import HeatUtility, Facility
from biosteam.units.decorators import cost
from thermosteam import Stream

__all__ = ('OrganicAcidsCIP', 'OrganicAcidsADP', 'OrganicAcidsCT', 'OrganicAcidsPWC', 
           'OrganicAcidsBT')

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
# Year  1997    1998    2009    2010    2016
# CE    386.5   389.5   521.9   550.8   541.7


# %% Clean-in-place system

@cost(basis='Flow rate', ID='System', units='kg/hr',
      cost=421000, S=63, CE=521.9, BM=1.8, n=0.6)
class OrganicAcidsCIP(Facility):
    network_priority = 2
    line = 'Clean-in-place system'


# %% Air distribution package

# Size based on stream 950 in Humbird et al.
@cost(basis='Flow rate', ID='Plant air compressor', units='kg/hr',
      kW=111.855, cost=28000, S=1372608, CE=550.8, n=0.6, BM=1.6)
@cost(basis='Flow rate', ID='Plant air reciever', units='kg/hr',
      cost=16000, S=1372608, CE=521.9, n=0.6, BM=3.1)
@cost(basis='Flow rate', ID='Instrument air dryer', units='kg/hr',
      cost=15000, S=1372608, CE=521.9, n=0.6, BM=1.8)
class OrganicAcidsADP(Facility): 
    line = 'Air distribution package'
    network_priority = 2
    

# %% Cooling tower

@cost('Flow rate', 'Cooling tower', units= 'kg/hr',
      kW=559.275, cost=1375000, S=10037820, CE=550.8, n=0.6, BM=1.5)
@cost('Flow rate', 'Cooling water pump', units='kg/hr',
      kW=1118.55, cost=283671, S=10982556,  CE=550.8, n=0.8, BM=3.1)
class OrganicAcidsCT(Facility):
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

    # Page 55 of Humbird et al.
    blowdown = 0.00005+0.0015
    
    def __init__(self, ID='', ins=None, outs=(), system_cooling_utilities={}, ratio=0):
        Facility.__init__(self, ID, ins, outs)
        self.system_cooling_utilities = system_cooling_utilities
        self.ratio = ratio
        
    def _run(self):
        return_cw, makeup_water, ct_chems = self.ins
        process_cw, blowdown = self.outs
        system_cooling_utilities = self.system_cooling_utilities

        # Based on stream 945 in Humbird et al.
        return_cw.T = 37 + 273.15
        # Based on streams 940/944 in Humbird et al.
        process_cw.T = blowdown.T = 28 + 273.15
        
        for u in self.system.units:
            if u is self: continue
            if hasattr(u, 'heat_utilities'):
                for cw in u.heat_utilities:
                    # hu.flow>0 is to avoid having BT's steam counted here
                    if cw.duty <0 and cw.flow >0:
                        # Note that technically, CT can only produce cooling water,
                        # but here includes all cool_utilities to check if
                        # cooling water is the only cooling utility in the system
                        system_cooling_utilities[f'{u.ID} - {cw.ID}'] = cw
        
        hu_cooling = HeatUtility().sum_by_agent(system_cooling_utilities.values())[0]
        hu_cooling.reverse()
        self.hu_cooling = hu_cooling
        self.heat_utilities = (hu_cooling,)
        self.system_cooling_demand = -hu_cooling.duty
        
        # Total amount of cooling water needed in the whole system (kmol/hr)
        total_cooling_water = self.total_cooling_water = -hu_cooling.flow
        return_cw.imol['H2O'] = process_cw.imol['H2O'] = total_cooling_water
        makeup_water.imol['H2O'] = total_cooling_water * self.blowdown
        blowdown.imol['H2O'] = makeup_water.imol['H2O']
        
        # 2 kg/hr from Table 30 on Page 63 of Humbird et al.
        ct_chems.imass['CoolingTowerChems'] = 2 * self.ratio

        self.design_results['Flow rate'] = total_cooling_water


# %% Process water center

# Assume all process water streams use fresh water
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      # Size basis changed from 451555 as the original design is not sufficient
      CE=521.9, cost=250000, S=106453, n=0.7, BM=1.7)
@cost(basis='Flow rate', ID='Circulating pump', units='kg/hr',
      CE=550.8, kW=55.9275, cost=15292, S=518924, n=0.8, BM=3.1)
@cost(basis='Flow rate', ID='Makeup water pump', units='kg/hr',
      CE=550.8, kW=14.914, cost=6864, S=155564, n=0.8, BM=3.1)
class OrganicAcidsPWC(Facility):
    network_priority = 1
    _N_ins = 1
    _N_outs = 1
    _units= {'Flow rate': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), process_water_streams=None):
        Facility.__init__(self, ID, ins, outs)
        self.process_water_streams = process_water_streams

    def _run(self):
        makeup_in = self.ins[0]
        makeup_out = self.outs[0]
        makeup_out.link_with(makeup_in)
        
        makeup_in.imol['Water'] = sum(i.imol['Water'] for i in self.process_water_streams)           
            
        self.design_results['Flow rate'] = makeup_in.F_mass
        

# %% Boiler and turbogenerator

@cost(basis='Flow rate', ID='Boiler', units='kg/hr',
      kW=2752, cost=28550000, S=238686, CE=550.8, n=0.6, BM=1.8)
@cost(basis='Work', ID='Turbogenerator', units='kW',
      cost=9500000, S=42200, CE=550.8, n=0.6, BM=1.8)
@cost(basis='Flow rate', ID='Hot process water softener system', units='kg/hr',
      cost=78000, S=235803, CE=550.8, n=0.6, BM=1.8)
@cost(basis='Flow rate', ID='Deaerator', units='kg/hr',
      cost=305000, S=235803, CE=550.8, n=0.6, BM=3)
@cost(basis='Flow rate', ID='Amine addition pkg', units='kg/hr',
      cost=40000, S=235803, CE=550.8, n=0, BM=1.8)
class OrganicAcidsBT(Facility):
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
                 TG_eff=0.85, combustibles, ratio, side_streams_to_heat=(),
                 system_heating_utilities={}):
        Facility.__init__(self, ID, ins, outs)
        self.B_eff = B_eff
        self.TG_eff = TG_eff
        self.combustibles = combustibles
        self.ratio = ratio
        self.side_streams_to_heat = side_streams_to_heat
        self.side_streams_lps = None
        self.system_heating_utilities = system_heating_utilities

    def _run(self): pass

    # Below cannot be put in the _run function because ins will be empty when simulated,
    # _design is run after _run so it'll work this way
    def _design(self):
        feed_solids, feed_gases, lime, boiler_chems, bag, natural_gas, makeup_water \
            = self.ins
        emission, ash, blowdown_water = self.outs
        ratio = self.ratio
        side_streams_to_heat = self.side_streams_to_heat
        side_streams_lps = self.side_streams_lps
        system_heating_utilities = self.system_heating_utilities
        lps = HeatUtility.get_heating_agent('low_pressure_steam')
        
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
        lime.imol['CalciumDihydroxide'] = emission.imol['SO2'] * 1.2
        # boiler_chems and bag are scaled based on plant size,
        # 1.23 is $2007/hour and 2.2661 is $/lb from Table 30 in Humbird et al.
        # 2.20462 is kg to lb
        boiler_chems.imass['BoilerChems'] = 1.23 / 2.2661 / 2.20462 * ratio
        bag.imass['BaghouseBag'] = ratio
        
        # 92% of SO2 removed by lime
        CaSO4_mol = emission.imol['SO2'] * 0.92
        emission.imol['SO2'] *= (1 - 0.92)

        # Air/O2 usage not rigorously modeled
        emission.imol['O2'] = 0
        
        for chemical in emission.chemicals:
            if chemical.locked_state == 'l' or chemical.locked_state == 's':
                ash.imol[chemical.ID] = emission.imol[chemical.ID]
                emission.imol[chemical.ID] = 0
                
        ash.mol += boiler_chems.mol
        ash.imol['CaSO4'] = CaSO4_mol
        ash.imol['Lime'] += lime.F_mol - CaSO4_mol
        
        emission.phase = 'g'
        ash.phase = 's'
        # Assume T of emission and ash are the same as lps, which
        # has lowest T amont all heating agents
        emission.T = ash.T = lps.T

        # Total heat generated by the boiler (kJ/hr)
        H_in = feed_solids.H + feed_gases.H
        H_out = emission.H + ash.H
        heat_from_combustion = -(feed_solids.LHV+feed_gases.LHV)
        heat_generated = self.heat_generated = \
            (H_in+heat_from_combustion)*self.B_eff - H_out
        
        for u in self.system.units:
            if u is self: continue
            if hasattr(u, 'heat_utilities'):
                for hu in u.heat_utilities:
                    # hu.flow>0 is to avoid having CT's cooling water counted here
                    if hu.duty >0 and hu.flow >0:
                        system_heating_utilities[f'{u.ID} - {hu.ID}'] = hu
        
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
        
        hu_BT = set()
        # BT can meet system heating/steam demand
        if BT_heat_surplus >0:
            # 3600 is conversion of kJ/hr to kW (kJ/s)
            electricity_generated = self.electricity_generated = \
                BT_heat_surplus * self.TG_eff / 3600
            
            # Take the opposite for cooling duty (i.e., cooling duty should be negative)
            # this is to condense the unused steam
            cooling_need = self.cooling_need = -(BT_heat_surplus-electricity_generated)
            hu_cooling = HeatUtility()
            hu_cooling(duty=cooling_need, T_in=lps.T)
            hu_BT.add(hu_cooling)
            
        # BT cannot meet system heating/steam demand, supplement with natural gas
        else:
            CH4_LHV = natural_gas.chemicals.CH4.LHV
            natural_gas.imol['CH4'] = BT_heat_surplus / (CH4_LHV*self.B_eff)
            emission.imol['CO2'] += natural_gas.imol['CH4']
            emission.imol['H2O'] += 2 * natural_gas.imol['CH4']
            electricity_generated = self.electricity_generated = 0
            

        for i in system_heating_utilities.values():
            j = HeatUtility()
            j.copy_like(i)
            j.reverse()
            hu_BT.add(j)
        
        total_steam = sum([i.flow for i in system_heating_utilities.values()])
        blowdown_water.imol['H2O'] = total_steam * self.blowdown
        blowdown_water.T = 373.15

        # Additional need from making lime slurry
        makeup_water.imol['H2O'] = blowdown_water.imol['H2O'] + lime.F_mol/0.2*0.8

        BT_utilities = self.BT_utilities = HeatUtility().sum_by_agent(hu_BT)
        self.heat_utilities = tuple(BT_utilities)
        Design = self.design_results        
        Design['Flow rate'] = total_steam
        Design['Work'] = electricity_generated

    def _end_decorated_cost_(self):
        self.power_utility(self.power_utility.rate - self.electricity_generated)
        
                

