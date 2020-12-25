# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 15:31:18 2020

@author: yrc2
"""
import flexsolve as flx
from biosteam.exceptions import InfeasibleRegion

class TiterAndInhibitorSpecification:
    __slots__ = (
        'evaporator',
        'mixer',
        'reactor',
        'product',
        'products',
        'sugars',
        'inhibitors',
        'target_titer',
        'maximum_inhibitor_concentration',
    )
    
    max_sugar_concentration = 600 # g / L
    
    def __init__(self, evaporator, mixer, reactor,
                 target_titer, product, products,
                 maximum_inhibitor_concentration=1.,
                 sugars= ('Glucose', 'Xylose'),
                 inhibitors = ('AceticAcid', 'HMF', 'Furfural')):
        self.evaporator = evaporator
        self.mixer = mixer
        self.reactor = reactor
        self.product = product
        self.products = products
        self.sugars = sugars
        self.inhibitors = inhibitors
        self.target_titer = target_titer
        self.maximum_inhibitor_concentration = maximum_inhibitor_concentration
        
    @property
    def feed(self):
        return self.evaporator.ins[0]
    
    @property
    def evaporated_sugar_solution(self):
        return self.evaporator.outs[0]
    
    @property
    def dilution_water(self):
        return self.mixer.ins[1]
    
    def run_units(self):
        self.evaporator._run()
        self.mixer._run()
        self.reactor._run()
    
    def calculate_titer(self):
        product = self.product
        return product.imass[self.products].sum() / product.F_vol
    
    def calculate_inhibitors(self): # g / L
        product = self.product
        return product.imass[self.inhibitors].sum() / product.F_vol
    
    def calculate_sugar_concentration(self): # g / L
        s = self.evaporated_sugar_solution
        return s.imass[sugars].sum() / s.F_vol 
    
    def check_sugar_concentration(self):
        if self.calculate_sugar_concentration() > self.max_sugar_concentration:
            raise InfeasibleRegion('sugar concentration')
    
    def titer_objective_function(self, V):
        self.evaporator.V = V
        self.run_units()
        return self.calculate_titer() - self.target_titer
    
    def inhibitor_objective_function(self, V):
        self.evaporator.V = V
        self.run_units()
        self.update_dilution_water()
        self.run_units()
        return self.calculate_inhibitors() - self.maximum_inhibitor_concentration
    
    def run(self):
        self.dilution_water.empty()
        self.evaporator.V = 0.
        self.run_units()
        x_titer = self.calculate_titer()
        V_min = 1e-6 
        if x_titer < self.target_titer: # Evaporate
            self.evaporator.V = V_min = flx.IQ_interpolation(self.titer_objective_function,
                                                             V_min, 1., ytol=1e-4, maxiter=100) 
        elif x_titer > self.target_titer: # Dilute
            self.update_dilution_water(x_titer)
            self.mixer._run()
            self.reactor._run()
        
        self.check_sugar_concentration()
        x_inhibitor = self.calculate_inhibitors()
        if x_inhibitor > self.maximum_inhibitor_concentration:
            self.evaporator.V = flx.IQ_interpolation(self.inhibitor_objective_function,
                                                     V_min, 1., ytol=1e-4, maxiter=100) 
        else:
            self.check_sugar_concentration()
    
    def update_dilution_water(self, x_titer=None):
        if x_titer is None: x_titer = self.calculate_titer()
        water = self.water_required_to_dilute_to_set_titer(x_titer)
        molar_volume = product.chemicals.Water.V(product.T, 101325) # m3 / mol
        self.dilution_water.imol['Water'] += water / molar_volume / 1000
        
    def water_required_to_dilute_to_set_titer(self, x_titer):
        return (1./x_titer - 1./self.target_titer) * self.product.imass[self.products].sum()
      
        