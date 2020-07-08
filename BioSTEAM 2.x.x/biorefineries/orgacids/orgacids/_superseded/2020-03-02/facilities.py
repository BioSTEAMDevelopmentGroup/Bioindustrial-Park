#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 09:53:03 2020

@author: yalinli_cabbi
"""

# %% Setup

import biosteam as bst
import thermosteam as tmo
from biosteam import Unit
from biosteam.units.decorators import cost
from thermosteam import Stream


# %% Clean-in-place package

@cost('Flow rate', units='kg/hr', ub=False,
      S=63, cost=421e3, CE=521.9, BM=1.8, n=0.6)
class CIPpackage(bst.Facility):
    line = 'Clean in place package'
    _N_ins = 1
    _N_outs = 1


# %% Boiler and turbogenerator


#!!! STILL NEEDS TO FIX THE DUTY THING



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
    generator. No emissions or mass balances are taken into account.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Liquid/solid feed that will be burned.
        
        [1] Gas that will be burned.
        
        [2] Make-up water. 
    outs : stream sequence
        [0] Emission from burned liquid/solid feed.
        
        [1] Emission from burned gas feed.
        
        [2] Blowdown water.
    boiler_efficiency : float
        Fraction of heat transfered to steam.
    turbo_generator_efficiency : float
        Fraction of steam heat converted to electricity.
    """

    duty_over_mol = 45000 # Superheat steam with 45000 kJ/kmol
    boiler_blowdown = 0.03
    RO_rejection = 0
    _N_outs = _N_ins = 3
    _N_heat_utilities = 2
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW'}
    
    def __init__(self, ID='', ins=None, outs=(), *,
                 boiler_efficiency=0.80,
                 turbogenerator_efficiency=0.85,
                 side_steam=None):
        Unit.__init__(self, ID, ins, outs)
        
        lps = bst.HeatUtility.heating_agents['Low pressure steam']
        self.makeup_water = makeup_water = Stream('boiler_makeup_water', 
                                                  thermo=lps.thermo)
        loss = makeup_water.flow_proxy()
        loss.ID = 'rejected_water_and_blowdown'
        self.ins[-1] = makeup_water
        self.outs[-1] = loss
        self.boiler_efficiency = boiler_efficiency
        self.turbogenerator_efficiency = turbogenerator_efficiency
        self.steam_utilities = set()
        self.steam_demand = Stream(None, thermo=lps.thermo, 
                                   P=lps.P, T=lps.T, phase='g')
        self.side_steam = side_steam
    
    def _run(self): pass
    
    def _design(self):
        B_eff = self.boiler_efficiency
        TG_eff = self.turbogenerator_efficiency
        steam = self.steam_demand
        steam_utilities = self.steam_utilities
        if not steam_utilities:
            for u in self.system.units:
                if u is self: continue
                for hu in u.heat_utilities:
                    if hu.ID == 'Low pressure steam':
                        steam_utilities.add(hu)
        steam.imol['Water'] = steam_mol = sum([i.flow for i in steam_utilities])
        duty_over_mol = self.duty_over_mol
        feed_solids, feed_gas, _ = self.ins
        emission_1, emission_2, _ = self.outs
        emission_1.phase = 's'
        emission_2.phase = 'g'
        hu_cooling, hu_steam = self.heat_utilities
        H_steam = steam_mol * duty_over_mol
        if self.side_steam: 
            H_steam += self.side_steam.H
        Hc_feed = 0
        for feed in (feed_solids, feed_gas):
            if feed: Hc_feed += feed.Hc
        
        # This is simply for the mass balance (no special purpose)
        # TODO: In reality, this should be modeled based on elemental composition
        # Turn all chemicals into ash
        if feed_solids: 
            ash = feed_solids.F_mass
            emission_1.empty()
            emission_1.imass['Ash'] = ash
        if feed_gas: emission_2.mol[:] = feed_gas.mol
        
        # A percent of solids is water, so remove latent heat of vaporization
        if feed_solids:
            F_mass_feed = feed_solids.F_mass
            moisture_content = feed_solids.imass['Water']/F_mass_feed
            H_moisture_evap = F_mass_feed*2300*moisture_content
        else:
            H_moisture_evap = 0
        H_content = max(0, (Hc_feed*B_eff-H_moisture_evap))
        
        if H_content == 0:
            evaporated_water = Hc_feed*B_eff/2300
            emission_2.imass['Water'] = feed_solids.imass['Water']-evaporated_water
        else: emission_2.imass['Water'] = 0
        #: [float] Total steam produced by the boiler (kmol/hr)
        self.total_steam = H_content / 40000 
        # Note: A portion of the steam produced is at milder conditions,
        #       so it does not consume as much energy.
        #       This is a really vague approximation, a more rigorous 
        #       model is needed (i.e. simulate whole system).
        
        self.makeup_water.imol['Water'] = makeup_mol = (
            self.total_steam * self.boiler_blowdown * 1/(1-self.RO_rejection)
        )
        # Heat available for the turbogenerator
        H_electricity = H_content - H_steam - makeup_mol * 18000
        
        Design = self.design_results
        Design['Flow rate'] = sum([i.get_total_flow('kg/hr') for i in self._outs])
        
        if H_electricity < 0:
            H_steam = H_content
            cooling = electricity = 0
        else:
            electricity = H_electricity * TG_eff
            cooling = electricity - H_electricity
        hu_cooling(cooling, steam.T)
        hu_steam.ID = 'Low pressure steam'
        hu_steam.cost = -sum([i.cost for i in steam_utilities])
        Design['Work'] = electricity/3600

    def _end_decorated_cost_(self):
        self.power_utility(self.power_utility.rate - self.design_results['Work'])

# Simulation of ethanol production from sugarcane
# in Brazil: economic study of an autonomous
# distillery
# Marina O.S. Diasa,b, Marcelo P. Cunhaa, Charles D.F. Jesusa, Mirna I.G.
# Scandiffioa, Carlos E.V. Rossella,b, Rubens Maciel Filhob, Antonio Bonomia
# a CTBE – Bioethanol Science and Technology National Laboratory, PO Box 6170 –
# CEP 13083-970, Campinas – SP, Brazil, marina.dias@bioetanol.org.br
# b School of Chemical Engineering, University of Campinas, PO Box 6066 – CEP
# 13083-970, Campinas – SP, Brazil
# LHV = 7565 # Bagasse lower heating value (kJ/kg)