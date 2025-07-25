# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 14:17:48 2024

@author: IGB
"""


import numpy as np
from math import exp, pi, log, ceil
import biosteam as bst
import thermosteam as tmo
from biosteam import Stream, Unit, BinaryDistillation
from biosteam.units import HXutility, Mixer, SolidsSeparator, Compressor
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch
from biosteam.exceptions import DesignError
from biosteam.units.design_tools.cost_index import CEPCI_by_year
from thermosteam import MultiStream
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

#%% 
# =============================================================================
# H2 production through Liquid Alkaline (LA) electrolysis 
# =============================================================================
# Model H2 production by H2O => H2 + 0.5O2
@cost(basis='Hydrogen flow rate', ID='Stack', S=500000/24, cost=346378000, 
      CE=CEPCI_by_year[2020], n=0.6, BM=1.01, lifetime=10)
@cost('Hydrogen flow rate', 'Mechanical BOP', S=500000/24, cost=181172000, 
      CE=CEPCI_by_year[2020], n=0.6, BM=1.38)
@cost('Hydrogen flow rate', 'Electrical BOP', S=500000/24, cost=155982000, 
      CE=CEPCI_by_year[2020], n=0.6, BM=1.05)
class Electrolyzer(Unit): 
    _N_ins = 1
    _N_outs = 2
    
    _units= {'Hydrogen flow rate': 'kg/hr'}
    
    def __init__(self, ID="", ins=None, outs=(),
                 base_kWh_per_kg_hydrogen=1208000/500000*24, # 1208 MW for 500 MT H2 per day; based on https://www.osti.gov/biblio/2203367 P17
                 **kwargs):
        bst.Unit.__init__(self, ID, ins, outs)
        self.base_kWh_per_kg_hydrogen = base_kWh_per_kg_hydrogen
        
    def _run(self):
        water, = self.ins
        hydrogen, oxygen = self.outs
        hydrogen.phase = 'g'
        hydrogen.P = 30 * 101325 # Model H2 purification so that output is 99.99% purity and 30 bar
        oxygen.phase = 'g'
        hydrogen.imol['H2'] = water.imol['Water']
        oxygen.imol['O2'] = water.imol['Water'] / 2
    
    def _design(self):
        Design = self.design_results
        Design['Hydrogen flow rate'] = self.outs[0].imass['H2']
        self.power_utility.power = self.base_kWh_per_kg_hydrogen * self.outs[0].imass['H2']

# =============================================================================
# MeOH synthesis
# =============================================================================

# Reactor is modelled as an adiabatic ideal plug flow reactor (PFR)
class MeOH_SynthesisReactor(bst.units.design_tools.PressureVessel, bst.Unit):
    _N_ins = _N_outs = 2
    
    _units = {**bst.design_tools.PressureVessel._units,
              'Volume': 'ft^3',
              'Catalyst weight': 'kg'}
    
    _F_BM_default = {**bst.design_tools.PressureVessel._F_BM_default,
                     'Catalyst': 1.0} # Cu/ZnO/Al2O3
    
    def __init__(self, ID="", ins=None, outs=(),
                 T=210+273.15,
                 P=7600000/6894.76,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 length_to_diameter=0.15/0.016, # 10.1016/j.jclepro.2013.06.008
                 WHSV = 0.028 * 3600 / 34.8, # Residence time in kg/hr feed / kg catalyst; 10.1016/j.jclepro.2013.06.008
                 catalyst_density = 1775, # In kg/m^3; 10.1016/j.jclepro.2013.06.008
                 catalyst_price = 30, # reference in process_settings
                 catalyst_longevity = 7884, # 1 year; from Perez report
                 porosity = 0.5,
                 **kwargs):
            bst.Unit.__init__(self, ID, ins, outs)
            self.T = T
            self.P = P
            self.vessel_material = vessel_material # Vessel material
            self.vessel_type = vessel_type # 'Horizontal' or 'Vertical'
            self.length_to_diameter = length_to_diameter #: Length to diameter ratio
            self.WHSV = WHSV
            self.catalyst_density = catalyst_density
            self.catalyst_price = 13 # In USD/kg
            self.catalyst_longevity = catalyst_longevity
            self.porosity = 0.5 #  fraction of the bed volume occupied by void space
            self.reactions = ParallelRxn([
                # Reaction definition                 Reactant    Conversion
                Rxn('CO2 + 3H2 -> CH3OH + H2O',       'CO2',      0.210),
                Rxn('CO2 + H2 -> CO + H2O',           'CO2',      0.004),])
                
    
    def _run(self):
        feed, fresh_catalyst = self.ins
        feed.phase = 'g'
        feed.T = self.T
        effluent, spent_catalyst = self.outs
        effluent.copy_like(feed)
        effluent.phase = 'g'
        fresh_catalyst.imass['CaO'] = self.catalyst_weight = feed.F_mass / self.WHSV / self.catalyst_longevity
        self.reactions.adiabatic_reaction(effluent)
        spent_catalyst.copy_like(fresh_catalyst)
        
    def _design(self):
        feed, fresh_catalyst = self.ins
        effluent, spent_catalyst = self.outs
        length_to_diameter = self.length_to_diameter
        P = effluent.get_property('P', 'psi')
        results = self.design_results
        
        results['Catalyst weight'] = self.catalyst_weight
        results['Volume'] = volume = self.catalyst_weight / self.catalyst_density / (1-self.porosity) * 35.3147 # to ft3
        D = cylinder_diameter_from_volume(volume, length_to_diameter)
        L = length_to_diameter * D
        results.update(self._vessel_design(P, D, L))
        
    def _cost(self):
        D = self.design_results
        self.purchase_costs.update(
            self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
        )
        self.purchase_costs['Catalyst'] = self.catalyst_price * D['Catalyst weight']

# =============================================================================
# Formic acid synthesis
# =============================================================================
class Reactor(bst.units.design_tools.PressureVessel, bst.Unit, isabstract=True):
    _N_ins = 5
    _N_outs = 1
    
    _units = {**bst.units.design_tools.PressureVessel._units,
              'Total volume': 'm3',
              'Reactor volume': 'm3',
              'Single reactor volume': 'm3'}
    
    # For a single reactor, based on diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _V_max = pi/4*(20**2)*40/35.3147
    
    def _setup(self):
        super()._setup()
        
    def _design(self):
        Design = self.design_results
        Design['Catalyst_weight'] = catalyst_weight = self.catalyst_weight
        
        V_total = (self.ins[-2].F_vol + self.ins[-1].F_vol) * self.residence_time / self.liquid_volume_frac
        P = self.P * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        wall_thickness_factor = self.wall_thickness_factor
        
        N = ceil(V_total/self._V_max)
        if N == 0:
            V_reactor = 0
            D = 0
            L = 0
        else:
            V_reactor = V_total / N
            D = (4*V_reactor/pi/length_to_diameter)**(1/3)
            D *= 3.28084 # convert from m to ft
            L = D * length_to_diameter

        Design['Total volume'] = V_total
        Design['Single reactor volume'] = V_reactor
        Design['Number of reactors'] = N
        P, D, L = float(P), float(D), float(L)
        Design.update(self._vessel_design(P, D, L))
        if wall_thickness_factor == 1: pass
        elif wall_thickness_factor < 1:
            raise DesignError('wall_thickness_factor must be larger than 1')
        else:
              Design['Wall thickness'] *= wall_thickness_factor
              # Weight is proportional to wall thickness in PressureVessel design
              Design['Weight'] = round(Design['Weight']*wall_thickness_factor,2)

    def _cost(self):    
        Design = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        
        if Design['Total volume'] == 0:
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] = 0
        
        else:
            baseline_purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] *= Design['Number of reactors']
    


# class HCOOH_SynthesisReactor(Reactor):
#     _N_ins = 5
#     _N_outs = 1
    
#     _F_BM_default = {**Reactor._F_BM_default}
    
#     auxiliary_unit_names = ('heat_exchanger',)
    
#     reaction = Rxn(
#         'CO2 + H2 + C6H15N -> TREAHCOOH',   'CO2',   0.84) # No CO2 recycled; from Accelerating the net-zero economy with CO2-hydrogenated formic acid production
    
#     def __init__(self, ID="", ins=None, outs=(),
#                  P=12000000, T=120+273.15, # Park et.al
#                  residence_time=0.5,
#                  liquid_volume_frac=0.2, # 0.05-0.25 from rules of thumb
#                  length_to_diameter=3,
#                  wall_thickness_factor=1,
#                  vessel_material='Stainless steel 316',
#                  vessel_type='Vertical'):
#         bst.Unit.__init__(self, ID, ins, outs)
#         self.P = P
#         self.T = T
#         self.residence_time = residence_time
#         self.liquid_volume_frac = liquid_volume_frac
#         self.length_to_diameter = length_to_diameter
#         self.wall_thickness_factor = wall_thickness_factor
#         self.vessel_material = vessel_material
#         self.vessel_type = vessel_type
#         self.heat_exchanger = self.auxiliary(
#             'heat_exchanger',
#             bst.HXutility,
#             ins=(bst.Stream(),),
#             T=self.T,)

#     def _run(self):
#         feed_CO2, feed_H2, recycled_gas, makeup_TREA, recycled_TREA = self.ins
#         effluent, = self.outs
#         required_makeup_TREA = 0.48 * (feed_CO2.imol['CO2'] + recycled_gas.imol['CO2']) - recycled_TREA.imol['triethylamine']
#         if required_makeup_TREA > 0:
#             makeup_TREA.imol['triethylamine'] = required_makeup_TREA
#         else:
#             makeup_TREA.imol['triethylamine'] = 0
        
#         mixed = bst.Stream()
#         mixed.mix_from([feed_CO2, feed_H2, recycled_gas, makeup_TREA, recycled_TREA], conserve_phases=True)
        
#         self.heat_exchanger.ins[0].copy_like(mixed)
#         self.heat_exchanger._run()
        
#         heated = self.heat_exchanger.outs[0]
#         self.reaction(heated)
        
#         effluent.copy_like(heated)
#         effluent.T = self.T
#         effluent.P = self.P
#         effluent.vle(T=effluent.T, P=effluent.P)
    
#     def _design(self):
#         self.catalyst_weight = (self.ins[0].imol['CO2'] +
#                                 self.ins[1].imol['CO2'] +
#                                 self.ins[2].imol['CO2'] ) * 0.48 * 46 / 669 * 24 # 669.0 gform. gcat−1 d−1
#         super()._design()
#         self.heat_exchanger._design()
        
#     def _cost(self):
#         super()._cost()
#         self.heat_exchanger._cost()


        
# class Amine_Exchange_Reactor(bst.Unit):
#     _N_ins = 3 
#     _N_outs = 1 
    
#     def _setup(self):
#         super()._setup()
#         self.reaction = Rxn(
#             'TREAHCOOH + nBIM -> nBIMHCOOH + C6H15N',  'TREAHCOOH', 1.)
    
#     def _run(self):
#         adduct, makeup_nBIM, recycled_nBIM = self.ins
#         effluent, = self.outs
#         required_makeup_nBIM = adduct.imol['TREAHCOOH'] - recycled_nBIM.imol['nBIM']
#         if required_makeup_nBIM > 0:
#             makeup_nBIM.imol['nBIM'] = required_makeup_nBIM
#         else:
#             makeup_nBIM.imol['nBIM'] = 0
#         effluent.mix_from([adduct, makeup_nBIM, recycled_nBIM], energy_balance=False)
#         self.reaction(effluent)
    
#     def _design(self):
#         pass
    
#     def _cost(self):
#         pass



            

# class nBIM_Exchange_Reactor(bst.Unit):
#     _N_ins = 1 
#     _N_outs = 1 
    
#     def _setup(self):
#         super()._setup()
#         self.reaction = Rxn(
#             'nBIMHCOOH -> nBIM + HCOOH',  'nBIMHCOOH', 1.)
    
#     def _run(self):
#         influent, = self.ins
#         effluent, = self.outs
#         effluent.copy_like(influent)
#         self.reaction(effluent)
    
#     def _design(self):
#         pass
    
#     def _cost(self):
#         pass