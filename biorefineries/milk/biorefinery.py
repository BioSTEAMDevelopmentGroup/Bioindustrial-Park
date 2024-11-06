# -*- coding: utf-8 -*-
"""
"""

__all__ = (
        
)
import os
import biosteam as bst
import numpy as np
import chaospy as cp
from warnings import catch_warnings
from .chemicals import create_chemicals
from .process_settings import (
    load_process_settings, GWP as GWPkey,
    characterization_factors as CFs
)
from .systems import create_system
from biorefineries.cane.data import price_distributions_2023 as dist
from biorefineries.tea import (
    create_cellulosic_ethanol_tea as create_tea
)
from inspect import signature

__all__ = (
    'Biorefinery',
    'ConfigurationKey',
    'ConfigurationComparison',
)

# DO NOT DELETE:
# natural_gas = bst.Chemical('CH4')
# natural_gas.phase = 'g'
# natural_gas.set_property('T', 60, 'degF')
# natural_gas.set_property('P', 14.73, 'psi')
# original_value = natural_gas.imol['CH4']
# natural_gas.imass['CH4'] = 1 
# V_ng = natural_gas.get_total_flow('m3/hr')
# natural_gas.imol['CH4'] = original_value
V_ng = 1.473318463076884 # Natural gas volume at 60 F and 14.73 psi [m3 / kg]

results_folder = os.path.join(os.path.dirname(__file__), 'results')

class ConfigurationKey:
    __slots__ = ('feed', 'product')
    
    def __init__(self, feed, product):
        self.feed = feed
        self.product = product
        
    def __sub__(self, other):
        return ConfigurationComparison(self, other)
        
    def __repr__(self):
        return f"ConfigurationKey({self.feed}, {self.product})"


class ConfigurationComparison:
    __slots__ = ('left', 'right')
    
    def __init__(self, left, right):
        self.left = left
        self.right = right
        
    def __repr__(self):
        return f"ConfigurationComparison(left={self.left}, right={self.right})"


class Biorefinery(bst.ProcessModel):
    """
    Examples
    --------
    >>> from biorefineries import milk
    >>> br = milk.Biorefinery(simulate=False)
    >>> br.evaluate_scenario()
    MSP [USD/kg]        2.27
    GWP [kg*CO2e/kg]   0.218
    TCI [MMUSD]         73.1
    dtype: float64
    
    """
    
    def __new__(
            cls,
            feed=None,
            product=None,
            simulate=True,
            cache={},
        ):
        if feed is None: feed = 'AcidWhey'
        if product is None: product = 'Dodecanol'
        key = (feed, product)
        if key in cache: return cache[key]
        self = super().__new__(cls)
        self.config = ConfigurationKey(*key)
        bst.settings.set_thermo(create_chemicals())
        load_process_settings()
        system = create_system(product=product, feed=feed)
        self.tea = create_tea(
            system,
        )
        self.tea.income_tax = 0.21
        self.load_system(system)
        model = bst.Model(system)
        parameter = model.parameter
        metric = model.metric
        
        def uniform(baseline, *args, **kwargs):
            bounds = [0.8 * baseline, 1.2 * baseline]
            return parameter(*args, **kwargs, baseline=baseline, bounds=bounds)
        
        @uniform(units='USD/kg', element='hexane', baseline=0.73) 
        def set_hexane_price(price):
            self.hexane.price = price
        
        @parameter(units='USD/kg', element='dodecylacetate', bounds=[3, 6],
                   baseline=3) # https://www.alibaba.com/product-detail/Factory-direct-sale-DODECYL-ACETATE-CAS_1601041319372.html
        def set_product_price(price):
            self.product.price = price
        
        @parameter(units='g/L', element=self.fermentation, 
                   bounds=[10, 60], baseline=30)
        def set_titer(titer):
            self.fermentation.titer = titer
            
        @parameter(units='g/L/h', element=self.fermentation,
                   bounds=[0.1, 1.3], baseline=1)
        def set_productivity(productivity):
            self.fermentation.productivity = productivity
        
        @parameter(units='% theoretical', name='yield', element=self.fermentation, 
                   bounds=[30, 90], baseline=50)
        def set_yield(X):
            X = X / 100.
            for reaction_series in self.fermentation.reactions: reaction_series.X[0] = X
        
        # https://www.eia.gov/energyexplained/natural-gas/prices.php
        # @parameter(distribution=dist.natural_gas_price_distribution, element='Natural gas', units='USD/m3',
        #            baseline=4.73 * 35.3146667/1e3)
        # def set_natural_gas_price(price): 
        #     self.BT.natural_gas_price = price * V_ng
    
        # @parameter(distribution=dist.electricity_price_distribution, units='USD/kWh',
        #            element='electricity', baseline=dist.mean_electricity_price)
        # def set_electricity_price(price): 
        #     bst.settings.electricity_price = price
        # Capacity based on U.S. production in 2015 - https://cen.acs.org/articles/95/i6/Acid-whey-waste-product-untapped.html
        @parameter(units='MT/yr', element=self.feedstock, bounds=[1e5, 1e6], baseline=5e5)
        def set_processing_capacity(processing_capacity):
            self.feedstock.F_mass = processing_capacity / system.operating_hours * 1000 # kg / hr
            
        self.BT.fuel.set_CF(GWPkey, CFs['Miscanthus'])
        self.BT.fuel.price = 0.08
        self.hexane.set_CF(GWPkey, CFs['Hexane'])
        
        @metric(units='USD/kg')
        def MSP():
            return self.tea.solve_price(self.product)
        
        @metric(units='kg*CO2e/kg')
        def GWP():
            return (
                self.system.get_net_impact(GWPkey)
                / self.system.get_mass_flow(self.product)
            )
        
        @metric(units='MMUSD')
        def TCI():
            return self.tea.TCI / 1e6
        
        for i in model.parameters: i.setter(i.baseline)
        for i in model.parameters:
            if i.distribution is None: i.distribution = cp.Uniform(*i.bounds)
        
        system.__class__.strict_convergence = False
        system.set_tolerance(rmol=1e-4, mol=1e-3, maxiter=50, 
                             method='wegstein', subsystems=True)
        if simulate:
            system.simulate()
            self.load_system(system)
        self.load_model(model)
        cache[key] = self
        return self

    def evaluate_scenario(self,
            hexane_price=None,
            dodecylacetate_price=None,
            titer=None,
            productivity=None,
            conversion=None,
            processing_capacity=None,
        ):
        sample = self.model.get_baseline_sample()
        values = [
            hexane_price,
            dodecylacetate_price,
            titer,
            productivity,
            conversion,
            processing_capacity,
        ]
        for n, value in enumerate(values):
            if value is None: continue
            sample.iloc[n] = value
        return self.model(sample)['-']
        