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
from .chemicals import create_acetate_ester_chemicals
from .process_settings import (
    load_process_settings, GWP as GWPkey,
    characterization_factors as CFs
)
from .systems import create_acetyl_ester_system
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
    __slots__ = ('carbon_capture', 'dewatering')
    name = 'AcEster'
    
    def __init__(self, carbon_capture, dewatering):
        self.carbon_capture = carbon_capture
        self.dewatering = dewatering
        
    def __sub__(self, other):
        return ConfigurationComparison(self, other)
        
    def __repr__(self):
        return f"ConfigurationKey(carbon_capture={self.carbon_capture}, dewatering={self.dewatering})"


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
    >>> import biorefineries.acester as ace
    >>> br = ace.Biorefinery(simulate=False)
    >>> br.system.simulate()
    >>> br.system.diagram()
    >>> (br.MSP(), br.GWP(), br.tea.TCI / 1e6)
    (7.30, 2.80, 438.53)
    
    """
    name = ConfigurationKey.name
    
    def optimize(self):
        results = self.model.optimize(
            self.MSP,  
            convergence_model='linear regressor', 
        )
        for p, x in zip(self.model.optimized_parameters, results.x): p.setter(x)
        return results
    
    def __new__(
            cls,
            simulate=True,
            dewatering=True,
            carbon_capture=True,
            cache={},
        ):
        key = (dewatering, carbon_capture)
        if key in cache: return cache[key]
        self = super().__new__(cls)
        self.config = ConfigurationKey(carbon_capture, dewatering)
        bst.settings.set_thermo(create_acetate_ester_chemicals())
        bst.settings.chemicals.define_group(
            'Biomass', 
            ['Water', 'Extract', 'Acetate', 'Ash', 'Lignin', 
             'Protein', 'Glucan', 'Xylan', 'Arabinan', 'Galactan'],
            [0.2, 0.031964705882352916, 0.004917647058823527, 
             0.015647058823529406, 0.1616, 0.008270588235294111, 
             0.3472, 0.2048, 0.0136, 0.012],
            wt=True,
        )
        load_process_settings()
        system = create_acetyl_ester_system(
            carbon_capture=carbon_capture,
            dewatering=dewatering
        )
        self.tea = create_tea(
            system,
        )
        self.tea.steam_power_depreciation = 'MACRS7'
        self.tea.income_tax = 0.21
        self.load_system(system)
        model = bst.Model(system)
        parameter = model.parameter
        metric = model.metric
        
        # repeated_names = [
        #     ('U_g', ('U_G', 'V_G'))
        # ]
        
        # fermentation_variables = {
        #     'epsilon_g': (0.5, 1.2), 
        # }
        # for name, others in repeated_names: 
        #     value = fermentation_variables[name]
        #     for other in others:
        #         fermentation_variables[other] = value
        
        # KLa_correlations = [
        #     f1, f2, f3
        # ]
        
        # for f in KLa_correlations:
        #     names = list(signature(f).parameters)
        #     for name in names:
        #         @parameter(element='KLa', bounds=fermentation_variables[name]) 
        #         def param(value):
        #             setattr(self.AcEster_production, name, value)
        
        def uniform(baseline, *args, **kwargs):
            bounds = [0.8 * baseline, 1.2 * baseline]
            return parameter(*args, **kwargs, baseline=baseline, bounds=bounds)
        
        # @uniform(units='USD/MT', element='Carbon capture', baseline=100) 
        # def set_carbon_capture_cost(price):
        #     self.CC.b = 4.230769230769226 + price
        
        @parameter(units='wt %', element='Boiler flue gas', baseline=5.5, bounds=(5, 10)) 
        def set_boiler_flue_gas_CO2_content(CO2_content):
            self.BT.CO2_emissions_concentration = CO2_content / 100
        
        @uniform(units='USD/kg', element='EtAc', baseline=1.57) 
        def set_ethyl_acetate_price(price):
            self.ethyl_acetate.price = price
        
        @uniform(units='USD/kg', element='hexane', baseline=0.73) 
        def set_hexane_price(price):
            self.hexane.price = price
        
        @parameter(units='USD/kg', element='dodecylacetate', bounds=[3, 6],
                   baseline=3) # https://www.alibaba.com/product-detail/Factory-direct-sale-DODECYL-ACETATE-CAS_1601041319372.html
        def set_dodecylacetate_price(price):
            self.dodecylacetate.price = price
        
        @parameter(units='USD/kg', element='H2', bounds=[3, 7],
                   baseline=4.5) # Natural gas 1.5 - 5; Electrolysis 3 - 7
        def set_H2_price(price):
            self.H2.price = price
        
        @parameter(units='g/L', element=self.AcOH_production,
                   bounds=[50, 90], baseline=60)
        def set_AcOH_titer(titer):
            self.AcOH_production.titer['AceticAcid'] = titer
        
        @parameter(units='g/L/h', element=self.AcOH_production,
                   bounds=[1, 3], baseline=2.25)
        def set_AcOH_productivity(productivity):
            self.AcOH_production.productivity = productivity
        
        @parameter(units='g/L', element=self.AcEster_production, 
                   bounds=[40, 90], baseline=70)
        def set_AcEster_titer(titer):
            self.AcEster_production.titer = titer
            
        @parameter(units='g/L/h', element=self.AcEster_production,
                   bounds=[0.2, 2], baseline=1)
        def set_AcEster_productivity(productivity):
            self.AcEster_production.productivity = productivity
        
        @parameter(units='% theoretical', name='yield', element=self.AcEster_production, 
                   bounds=[40, 90], baseline=50)
        def set_AcEster_yield(X):
            self.AcEster_production.reactions.X[0] = X / 100.
        
        # https://www.eia.gov/energyexplained/natural-gas/prices.php
        # @parameter(distribution=dist.natural_gas_price_distribution, element='Natural gas', units='USD/m3',
        #            baseline=4.73 * 35.3146667/1e3)
        # def set_natural_gas_price(price): 
        #     self.BT.natural_gas_price = price * V_ng
    
        # @parameter(distribution=dist.electricity_price_distribution, units='USD/kWh',
        #            element='electricity', baseline=dist.mean_electricity_price)
        # def set_electricity_price(price): 
        #     bst.settings.electricity_price = price
        
        @parameter(units='MT/yr', element='AcEster', bounds=[5000, 50000], baseline=30000)
        def set_production_capacity(production_capacity):
            capacity = production_capacity / system.operating_hours * 1000 # kg / hr
            original = self.AcOH_media.F_mass * self.AcOH_production.titer['AceticAcid'] * self.AcEster_production.reactions[0].product_yield('DodecylAcetate', basis='wt') / 1000
            self.system.rescale(self.AcOH_media, capacity / original)
        
        if carbon_capture:
            self.AcOH_production.length_to_diameter = 8
        else:
            self.AcOH_production.length_to_diameter = 12
            
        if dewatering:
            optimized_parameter = model.optimized_parameter
            
            @optimized_parameter(bounds=[5, 20], baseline=12, name='extractor stages')
            def set_stages(N):
                N = int(N)
                self.extractor.N_stages = N
            
            @optimized_parameter(bounds=[1, 2], baseline=1.5, name='solvent to feed ratio')
            def set_solvent_to_feed(solvent_feed_ratio):
                self.extractor.solvent_feed_ratio = solvent_feed_ratio
            
            @optimized_parameter(bounds=[1.01, 2.5], baseline=1.2, element=self.extract_distiller)
            def extract_distiller_reflux(k):
                self.extract_distiller.k = k
            
            @optimized_parameter(bounds=[90, 99.9], baseline=95, units='%', element=self.extract_distiller)
            def extract_distiller_heavy_key_recovery(Hr):
                self.extract_distiller.Hr = Hr / 100
            
            @optimized_parameter(bounds=[1.01, 2.5], baseline=1.2, element=self.raffinate_distiller)
            def raffinate_distiller_reflux(k):
                self.raffinate_distiller.k = k
            
            @optimized_parameter(bounds=[80, 99.9], baseline=99, units='%', element=self.raffinate_distiller)
            def raffinate_distiller_light_key_recovery(Lr):
                self.raffinate_distiller.Lr = Lr / 100
            
            @optimized_parameter(bounds=[80, 99.9], baseline=99, units='%', element=self.raffinate_distiller)
            def raffinate_distiller_heavy_key_recovery(Hr):
                self.raffinate_distiller.Hr = Hr / 100
            
            @optimized_parameter(bounds=[0, 1], baseline=0.5, element=self.extract_heater)
            def extract_heater_vapor_fraction(V):
                self.extract_heater.V = V
            
            @optimized_parameter(bounds=[2, 12], baseline=8, element='AcOH bioreactor', name='length to diameter')
            def set_AcOH_bioreactor_length_to_diameter(length_to_diameter):
                self.AcOH_production.length_to_diameter = length_to_diameter
            
            @optimized_parameter(bounds=[2, 12], baseline=4, element='AcEster bioreactor', name='length to diameter')
            def set_AcEster_bioreactor_length_to_diameter(length_to_diameter):
                self.AcEster_production.length_to_diameter = length_to_diameter
            
            # optimal_values = [15, 1.5e+00, 1.02e+00, 9.99e01, 1.02e+00, 9.99e01, 9.99e01, 0]
            # for i, j in zip(model.optimized_parameters, optimal_values): 
            #     i.baseline = j
            #     i.setter(j)
        # chemicals = bst.settings.chemicals
        # self.direct_nonbiogenic_emissions = lambda: (
        #     sum([i.get_atomic_flow('C') for i in self.system.products if i.phase == 'g'])
        #     * chemicals.CO2.MW * system.operating_hours
        # )
        
        # self.system.define_process_impact(
        #     key=GWPkey,
        #     name='Direct non-biogenic emissions',
        #     basis='kg',
        #     inventory=self.direct_nonbiogenic_emissions,
        #     CF=1.,
        # )
        if not dewatering: self.ethyl_acetate = bst.MockStream('ethyl_acetate')
        self.BT.fuel.set_CF(GWPkey, CFs['Miscanthus'])
        self.BT.fuel.price = 0.08
        self.hydrogen.set_CF(GWPkey, CFs['H2'])
        self.hexane.set_CF(GWPkey, CFs['Hexane'])
        self.ethyl_acetate.set_CF(GWPkey, CFs['Ethyl acetate'])
        
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

    def key_flows(self):
        def total(s):
            return s.get_total_flow('kg/s')
        
        def comps(s, components=None):
            if components is None: components = [i.ID for i in s.available_chemicals]
            total = s.get_total_flow('kg/s')
            dct = {'total': total}
            for i in components: 
                flow = s.get_flow('kg/s', i)
                if flow: dct[i] = round(flow / total, 3) * 100
            return dct
        
        def biomass(s):
            return {'total': s.get_total_flow('kg/s'),
                    'Moisture content': s.imass['Water'] / s.F_mass}
        
        return {
            'Flue gas': comps(self.emissions, ['CO2', 'H2O', 'O2', 'N2']),
            'Exhaust': comps(self.CC.outs[1], ['H2O', 'O2', 'N2']),
            'CO2': total(self.CO2),
            'Recycled CO2': comps(self.AcEster_production.outs[0], ['CO2', 'H2O', 'O2', 'N2']),
            'H2': total(self.system.flowsheet.stream.H2),
            'Unreacted H2': comps(self.AcOH_production.outs[0], ['H2', 'CO2', 'H2O', 'O2', 'N2']),
            'Dilute AcOH': comps(self.AcOH_production.outs[1], ['AceticAcid', 'H2O']),
            'Ethyl acetate': total(self.ethyl_acetate),
            'AcOH-Wastewater': comps(self.raffinate_distiller.outs[1]),
            'Glacial AcOH': comps(self.extract_distiller.outs[1]),
            'Air': total(self.air),
            'Wastewater': comps(self.M5.outs[1] + self.C2.outs[1]),
            'Biogas': comps(self.biogas),
            'Biogenic emissions': comps(self.biogenic_emissions),
            'Biomass': biomass(self.BT.fuel),
        }
