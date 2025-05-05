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
from thermosteam import MultiStream
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume
from qsdsan._sanunit import SanUnit

from biorefineries.ccu._chemicals import chems

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

tmo.settings.set_thermo(chems, cache=True)

#%% 
# =============================================================================
# H2 production through Liquid Alkaline (LA) electrolysis 
# =============================================================================
# Model H2 production by H2O => H2 + 0.5O2
class Electrolyzer(SanUnit): 
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID="", ins=None, outs=(),
                 stack_lifetime=10,
                 renewable_electricity_price=0.03, # currently possible in U.S. markets with plentiful  wind
                 **kwargs):
        bst.Unit.__init__(self, ID, ins, outs)
        self.stack_lifetime = stack_lifetime
        self.renewable_electricity_price = renewable_electricity_price
        
    def _run(self):
        water, = self.ins
        hydrogen, oxygen = self.outs
        hydrogen.phase = 'g'
        hydrogen.P = 30 * 101325 # Model H2 purification so that output is 99.99% purity and 30 bar
        oxygen.phase = 'g'
        hydrogen.imol['H2'] = water.F_mol
        oxygen.imol['O2'] = water.F_mol / 2
    
    def _cost(self): # based on https://www.osti.gov/biblio/2203367 P17
        self.power_utility.rate = self.outs[0].imass['H2'] * 122000 * 24 / 50000 # base case 50 MTD / 122 MW
        # self.power_utility.cost = self.renewable_electricity_price * self.power_utility.rate
        self.baseline_purchase_costs['Stack'] = self.power_utility.rate * 292 # process equipment, piping, valves, and  instrumentation including temperature, pressure, flow, and level indicators
        self.baseline_purchase_costs['Mechanical BOP'] = self.power_utility.rate * 280
        self.baseline_purchase_costs['Electrical BOP'] = self.power_utility.rate * 155 # for AC-> DC, transformer, power substation
        self.add_OPEX = self._calculate_replacement_cost()
    
    def _calculate_replacement_cost(self):
        stack_replacement_cost  = self.baseline_purchase_costs['Stack'] / self.stack_lifetime / (365*24) # convert from USD/year to USD/hour
        return stack_replacement_cost
    


    

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
                 catalyst_price = 13,
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
        fresh_catalyst.imass['CaO'] = self.catalyst_weight = feed.F_mass / self.WHSV
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
    _N_ins = 2
    _N_outs = 1
    
    auxiliary_unit_names = ('heat_exchanger',)
    
    _units = {**bst.units.design_tools.PressureVessel._units,
              'Total volume': 'm3',
              'Reactor volume': 'm3',
              'Single reactor volume': 'm3'}
    
    # For a single reactor, based on diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _V_max = pi/4*(20**2)*40/35.3147
    
    def __init__(self, ID="", ins=None, outs=(),
                 P=12000000, T=120+273.15, bulk_density=0.5, # 0.5 g/ml = 0.5 kg/L
                 length_to_diameter=3,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        bst.Unit.__init__(self, ID, ins, outs)
        self.P = P
        self.T = T
        self.bulk_density = bulk_density
        self.length_to_diameter = length_to_diameter
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
        heat_exchanger = self.auxiliary('heat exchanger',
                                        bst.HXutility,
                                        ins=bst.Stream(),
                                        T=self.T)
    def _setup(self):
        super()._setup()
        self.heat_exchanger = self.auxiliary(
            'heat_exchanger',
            bst.HXutility,
            ins=(bst.Stream(),),
            T=self.T,
        )
    
    def _run(self):
        feed, recycle = self.ins
        product, = self.outs
        
        # Reactor reaction
        product.mix_from([feed, recycle], energy_balance=False)
        self.reaction(product)
        product.T = self.T
        product.P = self.P
        
        # Send product to heat exchanger
        self.heat_exchanger.ins[0].copy_like(product)
        self.heat_exchanger._run()
    
    def _design(self):
        self.heat_exchanger._design()
        
        Design = self.design_results
        catalyst_weight = self.catalyst_weight
        V_total = catalyst_weight / self.bulk_density / 1000 # from L to m3
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
        self.heat_exchanger._cost()
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
    


class HCOOH_SynthesisReactor(Reactor):
    _N_ins = 3
    _N_outs = 1
    _F_BM_default = {**Reactor._F_BM_default}
    
    reaction = bst.Reaction(
        'CO2 + H2 + C6H15N -> TREAHCOOH',   'CO2',   0.48)
    
    def __init__(self, ID="", ins=None, outs=(), **kwargs):
        super().__init__(ID, ins, outs, **kwargs)
        
    def _run(self):
        feed, makeup_TREA, recycled_TREA = self.ins
        effluent, = self.outs
        effluent.mix_from([feed, makeup_TREA, recycled_TREA], energy_balance=False)
        self.reaction(effluent)
        effluent.T = self.T
        effluent.P = self.P
        self.heat_exchanger.ins[0].copy_like(effluent)
        self.heat_exchanger._run()
    
    def _design(self):
        self.catalyst_weight = (self.ins[0].imass['CO2'] + self.ins[0].imass['H2'] +
                                self.ins[1].imass['CO2'] + self.ins[1].imass['H2']) / 669 # 669.0 gform. gcat−1 d−1
        super()._design()
        
    def _cost(self):
        super()._cost()


        
class Amine_Exchange_Reactor(bst.Unit):
    _N_ins = 2 
    _N_outs = 1 
    
    def _setup(self):
        super()._setup()
        self.reaction = bst.Reaction(
            'TREAHCOOH + nBIM -> nBIMHCOOH + C6H15N',  'TREAHCOOH', 1.)
    
    def _run(self):
        adduct, nBIM = self.ins
        effluent, = self.outs
        nBIM.imol['nBIM'] = adduct.imol['TREAHCOOH']
        effluent.mix_from([adduct, nBIM], energy_balance=False)
        self.reaction(effluent)
    
    def _design(self):
        pass
    
    def _cost(self):
        pass



            

class nBIM_Exchange_Reactor(bst.Unit):
    _N_ins = 1 
    _N_outs = 1 
    
    def _setup(self):
        super()._setup()
        self.reaction = bst.Reaction(
            'nBIMHCOOH -> nBIM + HCOOH',  'nBIMHCOOH', 1.)
    
    def _run(self):
        influent, = self.ins
        effluent, = self.outs
        effluent.copy_like(influent)
        self.reaction(effluent)
    
    def _design(self):
        pass
    
    def _cost(self):
        pass