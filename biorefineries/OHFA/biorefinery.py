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
from .units import AeratedCoFermentation
from .systems import create_system
from biorefineries.cane.data import price_distributions_2023 as dist
from biorefineries.tea import (
    create_cellulosic_ethanol_tea as create_tea
)

__all__ = (
    'Biorefinery',
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

class Biorefinery(bst.ProcessModel):
    """
    Examples
    --------
    >>> import OHFA
    >>> model = OHFA.Biorefinery(simulate=False)
    >>> model.system.simulate()
    >>> model.system.diagram()
    >>> (model.MSP(), model.tea.TCI / 1e6)
    (1.66, 631.57)
    
    >>>
    
    """
    name = 'OHFA'
    
    def optimize(self):
        results = self.model.optimize(
            self.MSP,  
            convergence_model='linear regressor', 
        )
        for p, x in zip(self.model.optimized_parameters, results.x): p.setter(x)
        return results
    
    def __init__(
            self,
            simulate=True
        ):
        bst.settings.set_thermo(create_chemicals())
        load_process_settings()
        system = create_system()
        self.tea = create_tea(system)
        self.load_system(system)
        model = bst.Model(system)
        parameter = model.parameter
        metric = model.metric
        for OHFA_production in system.units:
            if OHFA_production.__class__ is AeratedCoFermentation: break
        OHFA_production.register_alias('OHFA_production')
        self.OHFA_production = self.R202
        def uniform(baseline, *args, **kwargs):
            bounds = [0.8 * baseline, 1.2 * baseline]
            return parameter(*args, **kwargs, baseline=baseline, bounds=bounds)
        
        @uniform(units='USD/kg', element='hexane', baseline=0.73) 
        def set_hexane_price(price):
            self.hexane.price = price
        
        @parameter(units='g/L', element=self.OHFA_production, 
                   bounds=[40, 100], baseline=70)
        def set_OHFA_titer(titer):
            self.OHFA_production.titer = titer
            
        @parameter(units='g/L/h', element=self.OHFA_production,
                   bounds=[0.2, 2], baseline=1)
        def set_OHFA_productivity(productivity):
            self.OHFA_production.productivity = productivity
        
        @parameter(units='% theoretical', name='yield', element=self.OHFA_production, 
                   bounds=[40, 90], baseline=50)
        def set_OHFA_yield(X):
            self.OHFA_production.cofermentation[0].X[:] = X / 100.
        
        # https://www.eia.gov/energyexplained/natural-gas/prices.php
        @parameter(distribution=dist.natural_gas_price_distribution, element='Natural gas', units='USD/m3',
                   baseline=4.73 * 35.3146667/1e3)
        def set_natural_gas_price(price): 
            self.BT.natural_gas_price = price * V_ng
    
        @parameter(distribution=dist.electricity_price_distribution, units='USD/kWh',
                   element='electricity', baseline=dist.mean_electricity_price)
        def set_electricity_price(price): 
            bst.settings.electricity_price = price
        
        chemicals = bst.settings.chemicals
        self.direct_nonbiogenic_emissions = lambda: (
            sum([i.get_atomic_flow('C') for i in self.system.products if i.phase == 'g'])
            * chemicals.CO2.MW * system.operating_hours
        )
        
        self.system.define_process_impact(
            key=GWPkey,
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=self.direct_nonbiogenic_emissions,
            CF=1.,
        )
        self.BT.natural_gas.set_CF(GWPkey, CFs['Natural gas'])
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
        
        for i in model.parameters: i.setter(i.baseline)
        for i in model.parameters:
            if i.distribution is None: i.distribution = cp.Uniform(*i.bounds)
        
        # system.__class__.strict_convergence = False
        system.set_tolerance(rmol=1e-4, mol=1e-3, maxiter=50, 
                             method='wegstein', subsystems=True)
        if simulate:
            system.simulate()
            self.load_system(system)
        self.load_model(model)
