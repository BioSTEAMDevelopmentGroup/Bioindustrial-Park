#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 09:53:03 2020

@author: yalinli_cabbi
"""

# %% Setup

import biosteam as bst
from biosteam import HeatUtility, Unit
from biosteam.units.decorators import cost
from thermosteam import Stream

__all__ = ('CIPpackage', 'OrganicAcidsBT')


# %% Clean-in-place package

@cost('Flow rate', units='kg/hr', ub=False,
      S=63, cost=421e3, CE=521.9, BM=1.8, n=0.6)
class CIPpackage(bst.Facility):
    line = 'Clean in place package'
    _N_ins = 1
    _N_outs = 1
    network_priority = 0


# %% Boiler and turbogenerator

@cost('Work', 'Turbogenerator',
      CE=551, S=42200, kW=0, cost=9500e3, n=0.60, BM=1.8)
@cost('Flow rate', 'Hot process water softener system', 
      CE=551, cost=78e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Amine addition pkg', 
      CE=551, cost=40e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Deaerator',
      CE=551, cost=305e3, S=235803, n=0.6, BM=3.0)
@cost('Flow rate', 'Boiler',
      CE=551, cost=28550e3, kW=1000, S=238686, n=0.6, BM=1.8)
class OrganicAcidsBT(bst.units.Facility):
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
        [0] Solid and gas emissions.
        
        [1] Blowdown water.
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

    _N_ins = 6
    _N_outs = 2
    _N_heat_utilities = 2
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW'}
    
    network_priority = 0
    duty_over_mol = 40000 # Superheat steam with 40000 kJ/kmol
    boiler_blowdown = 0.03
    RO_rejection = 0    
    generated_heat = 0
    total_steam = 0
    
    plant_size_ratio = 1
    combustables = []

    
    def __init__(self, ID='', ins=None, outs=(), *,
                 boiler_efficiency=0.80,
                 turbogenerator_efficiency=0.85, 
                 combustables=[], plant_size_ratio=1,
                 agent=HeatUtility.get_heating_agent('low_pressure_steam')):
        Unit.__init__(self, ID, ins, outs)
        
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
    
    def _run(self):
        B_eff = self.boiler_efficiency
        plant_size_ratio = self.plant_size_ratio
        feed_solids, feed_gases, lime, boiler_chemicals, bag, _ = self.ins
        emissions, _ = self.outs
        combustables = self.combustables
        emissions.mol = feed_solids.mol + feed_gases.mol

        # Use combustion reactions to create outs
        combustion_rxns = self.chemicals.get_combustion_reactions()
        combustable_feeds = Stream(None)
        combustable_feeds.copy_flow(emissions, tuple(combustables), remove=True)
        combustion_rxns(combustable_feeds.mol)
        emissions.mol += combustable_feeds.mol
        # Air usage not rigorously modeled, but can easily be added
        emissions.imol['O2'] = 0
        emissions.T = 373.15
        emissions.P = 101325

        # FGD lime scaled based on SO2 generated,
        # 20% stoichiometetric excess based on P52 of Humbird et al.
        # 895 kg/hr based on Table 30 in Humbird et al.
        lime.imol['CalciumDihydroxide'] = emissions.imol['SO2'] * 1.2
        # boiler_chemicals and bag are scaled based on plant size,
        # 1.23 is $2007/hour and 2.2661 is $/lb from Table 30 in Humbird et al.
        # 2.20462 is kg to lb
        boiler_chemicals.imass['BoilerChemicals'] = 1.23 / 2.2661 / 2.20462 * plant_size_ratio
        bag.imass['BaghouseBag'] = plant_size_ratio

        # Generated heat
        #!!! See if need to consider situations where not all water can be evaporated
        self.generated_heat = emissions.H - (feed_solids.LHV+feed_gases.LHV)*B_eff

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

        # Heat available for the turbogenerator
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
        Design['Work'] = electricity/3600 #!!! what is this 3600?

    def _end_decorated_cost_(self):
        self.power_utility(self.power_utility.rate - self.design_results['Work'])

