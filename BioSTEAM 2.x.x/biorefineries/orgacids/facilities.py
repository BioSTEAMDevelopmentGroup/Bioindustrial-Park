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
from biosteam.units._hx import HXutility
from biosteam.units.decorators import cost
from thermosteam import Stream

__all__ = ('OrganicAcidsCIP', 'OrganicAcidsADP', 'OrganicAcidsCT', 'OrganicAcidsPWC', 
           'OrganicAcidsBT')

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
# Year  1997    1998    2009    2010    2016
# CE    386.5   389.5   521.9   550.8   541.7


# %% Clean-in-place package

@cost(basis='Flow rate', ID='Clean-in-place package', units='kg/hr',
      cost=421000, S=63, CE=521.9, BM=1.8, n=0.6)
class OrganicAcidsCIP(Facility):
    network_priority = 2
    line = 'Clean-in-place package'


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
    _N_ins = 2
    _N_outs = 2    
    _N_heat_utilities = 1
    _units = {'Flow rate': 'kg/hr'}

    # Page 55 of Humbird et al.
    blowdown = 0.00005+0.0015
    
    def __init__(self, ID='', ins=None, outs=()):
        Facility.__init__(self, ID, ins, outs)
        self.cooling_water_utilities = set()
        
    def _run(self):
        return_cw, makeup_water = self.ins
        process_cw, blowdown = self.outs

        # Based on stream 945 in Humbird et al.
        return_cw.T = 37 + 273.15
        # Based on stream 941 in Humbird et al.
        makeup_water.T = 33 + 273.15
        # Based on streams 940/944 in Humbird et al.
        process_cw.T = blowdown.T = 28 + 273.15
        
        cwu = self.cooling_water_utilities
        if not cwu:
            for u in self.system.units:
                if u is self: continue
                for hu in u.heat_utilities:
                    if hu.ID == 'cooling_water':
                        cwu.add(hu)
        
        # Total amount of cooling water needed in the whole system (kmol/hr)
        total_cooling_water = self.total_cooling_water = sum(i.flow for i in cwu)
        return_cw.imol['H2O'] = process_cw.imol['H2O'] = total_cooling_water
        makeup_water.imol['H2O'] = total_cooling_water * self.blowdown
        
        hu = self.heat_utilities[0]
        hu.get_cooling_agent('cooling_water')
        hu.ID = 'cooling_water'
        # Regeneration cost is modeled by this unit
        hu.cost = 0

        self.design_results['Flow rate'] = total_cooling_water


# %% Process water center

@cost(basis='Total water flow rate', ID='Tank', units='kg/hr',
      # Size basis changed from 451555 as the original design is not sufficient
      CE=521.9, cost=250000, S=106453, n=0.7, BM=1.7)
@cost(basis='Total water flow rate', ID='PWC circulating pump', units='kg/hr',
      CE=550.8, kW=55.9275, cost=15292, S=518924, n=0.8, BM=3.1)
@cost(basis='Balance/discharged water flow rate', ID='Balance/discharged water pump', units='kg/hr',
      CE=550.8, kW=14.914, cost=6864, S=155564, n=0.8, BM=3.1)
class OrganicAcidsPWC(Facility):
    """
    Create a ProcessWaterCenter object that takes care of balancing the amount
    of water (and energy associated with the water) required for the process.
    The capital cost and power are based on the flow rate of process and makeup water as in [1]_.
    
    Parameters
    ----------
    ins :
        [0] Recycled water.
        
        [1] Balance water (>0 when recycled water < process water).
    outs :
        [0] Process water.
        
        [1] Discharged water (>0 when recycled water > process water).
    process_water_streams : streams, optional
        All process water streams.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    network_priority = 1
    _N_ins = 2
    _N_outs = 2
    _N_heat_utilities = 1
    _units = {'Total water flow rate': 'kg/hr',
              'Balance/discharged water flow rate': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), process_water_streams=None):
        Facility.__init__(self, ID, ins, outs)
        self.process_water_streams = process_water_streams
        self.HX = HXutility(None)

    def _run(self):
        recycled, balance = self.ins
        process, discharged = self.outs
        process_water_streams = self.process_water_streams        
        
        process.mix_from(process_water_streams)
        discharged.imol['Water'] = recycled.imol['Water'] - process.imol['Water']
        if discharged.imol['Water'] < 0:
            balance.imol['Water'] -= discharged.imol['Water']
            discharged.imol['Water'] = 0
            
        Design = self.design_results
        total_water = recycled.imass['Water'] + balance.imass['Water']
        Design['Total water flow rate'] = total_water
        Design['Balance/discharged water flow rate'] = max(balance.imass['Water'],
                                                           discharged.imass['Water'])
        
        HX = self.HX
        #!!! Why need a stream and should self.outs[0] be used?
        H_net = process.H - recycled.H
        HX.simulate_as_auxiliary_exchanger(H_net, process)
        self.purchase_costs['Heat exchanger'] = HX.purchase_costs['Heat exchanger']


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
    of the ins.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Liquid/solid feed that will be burned.
        
        [1] Gas that will be burned.
        
        [2] FGD lime.
        
        [3] Boiler chemicals.
        
        [4] Baghouse bag.
                
        [5] Makeup water.
    outs : stream sequence
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
    _N_ins = 6
    _N_outs = 3
    _N_heat_utilities = 2
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW'}
    
    blowdown = 0.03
    
    def __init__(self, ID='', ins=None, outs=(), *, B_eff=0.8,
                 TG_eff=0.85, combustables, ratio):
        Facility.__init__(self, ID, ins, outs)
        
        self.B_eff = B_eff
        self.TG_eff = TG_eff
        self.combustables = combustables
        self.ratio = ratio
        # All needed steam streams in the system, using set to avoid adding duplicates
        self.steam_utilities = set()
    
    def _run(self): pass

    # Below cannot be put in the _run function because ins will be empty when simulated,
    # _design is run after _run so it'll work this way
    def _design(self):
        feed_solids, feed_gases, lime, boiler_chemicals, bag, makeup_water = self.ins
        emission, ash, blowdown_water = self.outs
        # hu_cooling used to condense blowdown
        hu_steam, hu_cooling = self.heat_utilities
        ratio = self.ratio
        steam_utilities = self.steam_utilities
        
        # Use combustion reactions to create outs
        combustion_rxns = self.chemicals.get_combustion_reactions()
        combustable_feeds = Stream(None)
        emission.mol = feed_solids.mol + feed_gases.mol
        combustable_feeds.copy_flow(emission, tuple(self.combustables), remove=True)
        combustion_rxns(combustable_feeds.mol)
        emission.mol += combustable_feeds.mol

        # FGD lime scaled based on SO2 generated,
        # 20% stoichiometetric excess based on P52 of Humbird et al.
        # 895 kg/hr based on Table 30 in Humbird et al.
        lime.imol['CalciumDihydroxide'] = emission.imol['SO2'] * 1.2
        # boiler_chemicals and bag are scaled based on plant size,
        # 1.23 is $2007/hour and 2.2661 is $/lb from Table 30 in Humbird et al.
        # 2.20462 is kg to lb
        boiler_chemicals.imass['BoilerChemicals'] = 1.23 / 2.2661 / 2.20462 * ratio
        bag.imass['BaghouseBag'] = ratio
        
        # 92% of SO2 removed by lime
        CaSO4_mol = emission.imol['SO2'] * 0.92
        emission.imol['SO2'] *= (1 - 0.92)

        # Air usage not rigorously modeled
        emission.imol['O2'] = 0
        
        ash_mass = 0
        for chemical in emission.chemicals:
            if chemical.locked_state == 'l' or chemical.locked_state == 's':
                ash_mass += emission.imass[chemical.ID]
                emission.imass[chemical.ID] = 0
        ash.imass['Ash'] = ash_mass + CaSO4_mol*136.14 + (lime.F_mol-CaSO4_mol)*56.0774 \
                           + boiler_chemicals.F_mass

        # Total heat generated by the boiler (kJ/hr), 
        # LHV is initially negative so take the opposite here
        heat_generated = self.heat_generated = \
            -(feed_solids.LHV+feed_gases.LHV-emission.H-ash.H)*self.B_eff
        
        #!!! For trouble-shotting
        # print('Generated heat is '+str(heat_generated))

        # To get steam demand of the whole system
        if not steam_utilities:
            for u in self.system.units:
                if u is self: continue
                for hu in u.heat_utilities:
                    if hu.ID in (i.ID for i in HeatUtility.heating_agents):
                        steam_utilities.add(hu)
                        #!!! For trouble-shotting
                        # print(u.ID)
        
        # Total demand of steam by other units in the whole system (kmol/hr)           
        steam_demand = self.steam_demand = sum(i.flow for i in steam_utilities)
        # Heat needed to generate the steam
        heat_demand = self.heat_demand = sum(i.agent.H*i.flow for i in steam_utilities)
        
        #!!! For trouble-shotting
        # print('Total heat demand is '+str(heat_demand)) 
        
        # Use low_pressure_steam as the agent
        lps = HeatUtility.get_heating_agent('low_pressure_steam')
        hu_steam.agent = lps
        # Regeneration cost is modeled by this unit
        hu_steam.cost = 0
        
        # heat_generated - heat_demand
        heat_surplus = self.heat_surplus = max(0, heat_generated - heat_demand)
        # Steam generated by the boiler with the surplus energy (kmol/hr)
        steam_surplus = self.steam_surplus = heat_surplus / lps.H
        blowdown_water.imol['H2O'] = steam_surplus * self.blowdown
        # Additional need from making lime slurry
        makeup_water.imol['H2O'] = blowdown_water.imol['H2O'] + lime.F_mol/0.2*0.8

        # 3600 is conversion of kJ/hr to kW (kJ/s)
        electricity = heat_surplus * self.TG_eff / 3600
        # Take the opposite to for cooling
        cooling = -(heat_surplus - electricity)
        
        hu_cooling.get_cooling_agent('cooling_water')
        hu_cooling(duty=cooling, T_in=lps.T)

        ash.phase = 's'
        emission.phase = 'g'
        makeup_water.phase = 'l'
        
        for i in range(3):
            self.outs[i].T = 373.15
            self.outs[i].P = 101325

        Design = self.design_results
        Design['Flow rate'] = (steam_demand+steam_surplus) * 18.01528
        Design['Work'] = electricity

    def _end_decorated_cost_(self):
        self.power_utility(self.power_utility.rate - self.design_results['Work'])
                

