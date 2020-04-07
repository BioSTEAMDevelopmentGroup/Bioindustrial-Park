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

__all__ = ('CIPpackage', 'OrganicAcidsADP', 'OrganicAcidsPWC', 'OrganicAcidsBT')

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
# Year  1997    1998    2009    2010    2016
# CE    386.5   389.5   521.9   550.8   541.7


# %% Clean-in-place package

@cost(basis='Flow rate', ID='Clean-in-place package', units='kg/hr',
      cost=421000, S=63, CE=521.9, BM=1.8, n=0.6)
class CIPpackage(Facility):
    network_priority = 2
    line = 'Clean in place package'


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
    The capital cost and power are based on the flow rate of process and make-up water as in [1]_.
    
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
    
    def __init__(self, ID='', ins=None, outs=(),
                 process_water_streams=None):
        Facility.__init__(self, ID, ins, outs)
        self.process_water_streams = process_water_streams
        self.HX = HXutility(None)
    
    # def _assert_compatible_property_package(self): pass

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
                
        [5] Make-up water.
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
    
    duty_over_mol = 0
    boiler_blowdown = 0.03
    RO_rejection = 0    
    generated_heat = 0
    total_steam = 0
    
    def __init__(self, ID='', ins=None, outs=(), *,
                 boiler_efficiency=0.80,
                 turbogenerator_efficiency=0.85, 
                 combustables, ratio,
                 agent=HeatUtility.get_heating_agent('low_pressure_steam')):
        Facility.__init__(self, ID, ins, outs)
        
        self.agent = agent
        self.makeup_water = makeup_water = agent.to_stream('boiler_makeup_water')
        loss = makeup_water.flow_proxy()
        loss.ID = 'rejected_water_and_blowdown'
        self.ins[-1] = makeup_water
        self.outs[-1] = loss
        self.boiler_efficiency = boiler_efficiency
        self.turbogenerator_efficiency = turbogenerator_efficiency
        self.steam_utilities = set()
        self.steam_demand = agent.to_stream()
        self.steam_demand.chemicals.set_synonym('Water', 'H2O')
        self.combustables = combustables
        self.ratio = ratio
        self.duty_over_mol = agent.H * agent.heat_transfer_efficiency
    
    def _run(self):
        B_eff = self.boiler_efficiency
        ratio = self.ratio
        feed_solids, feed_gases, lime, boiler_chemicals, bag, _ = self.ins
        emission, ash, _ = self.outs
        combustables = self.combustables
        emission.mol = feed_solids.mol + feed_gases.mol

        # Use combustion reactions to create outs
        combustion_rxns = self.chemicals.get_combustion_reactions()
        combustable_feeds = Stream(None)
        combustable_feeds.copy_flow(emission, tuple(combustables), remove=True)
        combustion_rxns(combustable_feeds.mol)
        emission.mol += combustable_feeds.mol
        # Air usage not rigorously modeled, but can easily be added
        emission.imol['O2'] = 0

        # FGD lime scaled based on SO2 generated,
        # 20% stoichiometetric excess based on P52 of Humbird et al.
        # 895 kg/hr based on Table 30 in Humbird et al.
        lime.imol['CalciumDihydroxide'] = emission.imol['SO2'] * 1.2
        # boiler_chemicals and bag are scaled based on plant size,
        # 1.23 is $2007/hour and 2.2661 is $/lb from Table 30 in Humbird et al.
        # 2.20462 is kg to lb
        boiler_chemicals.imass['BoilerChemicals'] = 1.23 / 2.2661 / 2.20462 * ratio
        bag.imass['BaghouseBag'] = ratio
        
        ash_mass = 0
        for chemical in emission.chemicals:
            if chemical.locked_state == 'l' or chemical.locked_state == 's':
                ash_mass += emission.imass[chemical.ID]
                emission.imass[chemical.ID] = 0
        ash.imass['Ash'] = ash_mass
        ash.phase = 's'
        emission.phase = 'g'
        self.outs[2].phase = 'l'

        for i in range(0, 3):
            self.outs[i].T = 373.15
            self.outs[i].P = 101325

        # Generated heat
        #!!! See if need to consider situations where not all water can be evaporated
        self.generated_heat = emission.H - (feed_solids.LHV+feed_gases.LHV)*B_eff

    def _design(self):
        TG_eff = self.turbogenerator_efficiency
        steam = self.steam_demand
        steam_utilities = self.steam_utilities
        generated_heat = self.generated_heat
        duty_over_mol = self.duty_over_mol
        #: [float] Total steam produced by the boiler (kmol/hr)
        total_steam = self.total_steam = self.generated_heat / duty_over_mol
        
        if not steam_utilities:
            for u in self.system.units:
                if u is self: continue
                for hu in u.heat_utilities:
                    if hu.ID == self.agent.ID:
                        steam_utilities.add(hu)
        steam.imol['H2O'] = steam_mol = sum([i.flow for i in steam_utilities])
        hu_cooling, hu_steam = self.heat_utilities
        H_steam = steam_mol * duty_over_mol
        #: [float] Total steam produced by the boiler (kmol/hr)
        # Note: A portion of the steam produced is at milder conditions,
        #       so it does not consume as much energy.
        #       This is a really vague approximation, a more rigorous 
        #       model is needed (i.e. simulate whole system).

        self.makeup_water.imol['H2O'] = total_steam * self.boiler_blowdown * \
                                        1/(1-self.RO_rejection)

        # Heat available for the turbogenerator, in kJ/hr
        H_electricity = generated_heat - H_steam
        
        Design = self.design_results
        Design['Flow rate'] = total_steam * 18.01528
        
        if H_electricity < 0:
            H_steam = generated_heat
            cooling = electricity = 0
        else:
            electricity = H_electricity * TG_eff
            cooling = electricity - H_electricity
        hu_cooling(cooling, steam.T)
        hu_steam.agent = self.agent
        hu_steam.cost = -sum([i.cost for i in steam_utilities])
        # 3600 is conversion of kJ/hr to kW (kJ/s)
        Design['Work'] = electricity/3600

    def _end_decorated_cost_(self):
        self.power_utility(self.power_utility.rate - self.design_results['Work'])


