# -*- coding: utf-8 -*-
"""
"""
import os
import biosteam as bst
import numpy as np
import chaospy as cp
from warnings import catch_warnings
from .property_package import create_acetate_ester_chemicals
from .process_settings import (
    load_process_settings, GWP as GWPkey,
    characterization_factors as CFs
)
from .systems import create_hydrogen_oleochemical_system, create_substrate_oleochemical_system
from biorefineries.cane.data import price_distributions_2023 as dist
from biorefineries.tea import (
    create_cellulosic_ethanol_tea as create_tea
)
from warnings import catch_warnings
from inspect import signature
from chaospy import distributions as shape

__all__ = (
    'Biorefinery',
    'TRYBiorefinery',
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
    >>> from biorefineries.gas_fermentation import Biorefinery
    >>> br = Biorefinery(simulate=False, scenario='glucose growth')
    >>> br.system.simulate()
    >>> assumptions, results = br.baseline()
    >>> br.system.diagram() # View diagram
    >>> print(assumptions)
    EtAc                     Price [USD/kg]                                1.57
    hexane                   Price [USD/kg]                                0.73
    glucose                  Price [USD/kg]                               0.413
    oleochemical             Price [USD/kg]                                   3
    H2                       Price [USD/kg]                                   3
    AcOH production          Titer [g/L]                                     60
                             Productivity [g/L/h]                           1.5
    Oleochemical production  Titer [g/L]                                      5
                             Productivity [g/L/h]                             1
                             Bioreactor yield [% theoretical]                36
                             Specific yield [g_{Dodecanol}/g_{cell}]       1.57
    Flue gas                 Processing capacity [MT/yr]               2.32e+05
    biomass                  Price [USD/MT]                                54.7
    dtype: float64
    
    >>> print(results)
    -             MSP [USD/kg]                                    7.43
                  Carbon intensity [kg*CO2e/kg]                   2.26
                  TCI [10^6 USD]                                   646
                  Product yield to biomass [wt %]                0.675
                  Product yield to hydrogen [% theoretical]       0.14
                  Biomass burned [10^3 MT/yr]                     80.7
                  Hydrogen consumption [10^3 MT/yr]               54.6
                  Electricity demand [kWh/kg-H2]                  13.7
    oleochemical  Production capacity [10^3 MT/yr]            5.45e+04
    dtype: float64
    
    """
    class Scenario:
        name: str = 'acetate/glucose-seed'
        product: str = 'Dodecanol'
        carbon_capture: bool = False
        dewatering: bool = False
        glucose_growth: bool = True
        carbon_source: str = 'BFG' # Alternatively 'biomass' or 'glucose'
        biomass: str = 'cornstover' # Either 'miscanthus' or 'cornstover'Al
        hydrogen_price: str = 'all' # Alternatively 'min', 'max', or 'baseline'
        fermentation_performance: str = 'all' # Alternatively 'baseline yield' or 'potential yield'
    
    @property
    def name(self):
        return self.scenario.name
    
    @classmethod
    def as_scenario(cls, scenario):
        if scenario == 'acetate':
            carbon_source = 'BFG'
            glucose_growth = False
        elif scenario == 'glucose':
            carbon_source = 'glucose'
            glucose_growth = True
        elif scenario == 'acetate/glucose-seed':
            carbon_source = 'BFG'
            glucose_growth = True
        else:
            raise ValueError("invalid scenario")
        return cls.Scenario(
            product='Dodecanol',
            carbon_capture=False,
            dewatering=False,
            glucose_growth=glucose_growth,
            carbon_source=carbon_source,
            biomass='cornstover',
            hydrogen_price='all',
            fermentation_performance='all',
            name=scenario,
        )
    
    def optimize(self):
        with catch_warnings(action="ignore"):
            results, convergence_model = self.model.optimize(
                self.MSP,  
                convergence_model='linear regressor', 
                # method='differential evolution',
            )
        for p, x in zip(self.model.optimized_parameters, results.x): p.setter(x)
        return results
    
    def create_thermo(self):
        return create_acetate_ester_chemicals()
    
    def create_system(self):
        scenario = self.scenario
        if scenario.biomass == 'miscanthus':
            self.chemicals.define_group(
                'Biomass', 
                ['Water', 'Extract', 'Acetate', 'Ash', 'Lignin', 
                 'Protein', 'Glucan', 'Xylan', 'Arabinan', 'Galactan'],
                [0.2, 0.031964705882352916, 0.004917647058823527, 
                 0.015647058823529406, 0.1616, 0.008270588235294111, 
                 0.3472, 0.2048, 0.0136, 0.012],
                wt=True,
            )
        elif scenario.biomass == 'cornstover':
            self.chemicals.define_group(
                'Biomass', 
                ['Water',
                 'Sucrose',
                 'Extract',
                 'Acetate',
                 'Ash',
                 'Lignin',
                 'Protein',
                 'Glucan',
                 'Xylan',
                 'Arabinan',
                 'Mannan',
                 'Galactan'],
                [0.2022,
                 0.0062,
                 0.1185,
                 0.0146,
                 0.0399,
                 0.1274,
                 0.0251,
                 0.283,
                 0.1579,
                 0.0192,
                 0.0049,
                 0.0012],
                wt=True,
            )
        if scenario.carbon_source == 'BFG':
            # Composition
            # https://www.researchgate.net/publication/342635556_Advanced_Packed-Bed_Ca-Cu_Looping_Process_for_the_CO2_Capture_From_Steel_Mill_Off-Gases
            bfg = bst.Stream(CO=23, CO2=23, H2=4.5, N2=49.5, units='m3/hr', phase='g')
            gases = 'CO', 'CO2', 'H2', 'N2'
            gas_flows = bfg.imass[gases] 
            self.chemicals.define_group(
                'BFG', gases, gas_flows, wt=True
            )
        load_process_settings()
        if scenario.carbon_source == 'BFG':
            system = create_hydrogen_oleochemical_system(
                carbon_capture=scenario.carbon_capture,
                dewatering=scenario.dewatering,
                glucose_growth=scenario.glucose_growth,
                product=scenario.product,
                carbon_source=scenario.carbon_source,
            )
            areas = ('Acetate prod.',
                     'Oleochem prod',
                     'Facilities')
        elif scenario.carbon_source == 'glucose':
            system = create_substrate_oleochemical_system(
                specific_yield=2.7, 
                substrate='glucose',
                product=scenario.product,
                glucose_growth=scenario.glucose_growth,
                facilities=True,
            )
            system.flowsheet.emissions.ID = 'flue_gas'
            self.feedstock = feedstock = system.ins[1]
            feedstock.empty()
            glucose = 2.84e+04 # Typical amount produced by corn-dry grind
            feedstock.imass['Water', 'Glucose'] = [9 * glucose, glucose]
            self.feedstock = feedstock
            feedstock.set_CF('GWP', 0.1 * 0.9375) # Sugar GREET 2023
            feedstock.price = 0.413 * 0.10
            feedstock.ID = 'feedstock'
            areas = ('Oleochem prod',
                     'Facilities')
        self.unit_groups = unit_groups = bst.UnitGroup.group_by_area(system.units)
        for i, j in zip(unit_groups, areas): i.name = j
        self.tea = create_tea(
            system,
        )
        self.tea.steam_power_depreciation = 'MACRS7'
        self.tea.income_tax = 0.21
        system.__class__.strict_convergence = False
        system.set_tolerance(rmol=1e-4, mol=1e-3, maxiter=50, 
                             method='wegstein', subsystems=True)
        return system
    
    def create_model(self):
        if not hasattr(self, 'BT'): self.BT = self.B
        system = self.system
        scenario = self.scenario
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
        #             setattr(self.oleochemical_production, name, value)
        
        def uniform(baseline, *args, **kwargs):
            bounds = [0.8 * baseline, 1.2 * baseline]
            return parameter(*args, **kwargs, baseline=baseline, bounds=bounds)
        
        @uniform(units='USD/kg', element='EtAc', baseline=1.57) 
        def set_ethyl_acetate_price(price):
            self.ethyl_acetate.price = price
        
        @uniform(units='USD/kg', element='hexane', baseline=0.73) 
        def set_hexane_price(price):
            self.hexane.price = price
        
        # https://tradingeconomics.com/commodity/sugar
        @parameter(units='USD/kg', element='glucose', bounds=[0.215, 0.611]) 
        def set_glucose_price(price): # https://scijournals.onlinelibrary.wiley.com/doi/epdf/10.1002/bbb.1976?saml_referrer
            if scenario.glucose_growth:
                self.seedtrain_feed.price = price * 0.12
            if scenario.carbon_source == 'glucose':
                self.feedstock.price = price * 0.10
        
        @parameter(units='USD/kg', element='oleochemical', bounds=[3, 8],
                   baseline=5, distribution='uniform') # https://www.alibaba.com/product-detail/Factory-direct-sale-DODECYL-ACETATE-CAS_1601041319372.html
        def set_product_price(price):
            self.product.price = price
        
        # https://www.hydrogen.energy.gov/docs/hydrogenprogramlibraries/pdfs/20004-cost-electrolytic-hydrogen-production.pdf?Status=Master
        # Actual cost of hydrogen between 4-6 USD/kg according to DOE
        @parameter(units='USD/kg', element='H2', bounds=[2, 6],
                   baseline=3, distribution='uniform') # Natural gas 1.5 - 5; Electrolysis 3 - 7
        def set_H2_price(price):
            if scenario.carbon_source == 'glucose': return
            self.hydrogen.price = price
        
        @parameter(units='g/L', element='AcOH production',
                   bounds=[40, 90], baseline=60, distribution='uniform')
        def set_AcOH_titer(titer):
            if scenario.carbon_source == 'glucose': return
            self.AcOH_production.titer['AceticAcid'] = titer
        
        @parameter(units='g/L/h', element='AcOH production',
                   bounds=[1, 2], baseline=1.5, distribution='uniform')
        def set_AcOH_productivity(productivity):
            if scenario.carbon_source == 'glucose': return
            self.AcOH_production.productivity = productivity
        
        @parameter(units='g/L', element='Oleochemical production', 
                   bounds=[1, 10], baseline=5, distribution='uniform')
        def set_oleochemical_titer(titer):
            self.oleochemical_production.titer = titer
            
        @parameter(units='g/L/h', element='Oleochemical production',
                   bounds=[0.1, 2.0], baseline=1.0, distribution='uniform')
        def set_oleochemical_productivity(productivity):
            self.oleochemical_production.productivity = productivity
        
        @parameter(units='% theoretical', name='bioreactor yield', element='Oleochemical production', 
                   bounds=[35, 85], baseline=36, distribution='uniform')
        def set_oleochemical_bioreactor_yield(X):
            self.oleochemical_production.reactions.X[0] = X / 100.
        
        @parameter(units='g_{' + scenario.product + '}/g_{cell}', element='Oleochemical production', 
                   bounds=[0.45, 2.7], baseline=1.57, distribution='uniform')
        def set_oleochemical_specific_yield(specific_yield):
            self.oleochemical_production.specific_yield = specific_yield
        
        # https://www.eia.gov/energyexplained/natural-gas/prices.php
        # @parameter(distribution=dist.natural_gas_price_distribution, element='Natural gas', units='USD/m3',
        #            baseline=4.73 * 35.3146667/1e3)
        # def set_natural_gas_price(price): 
        #     self.BT.natural_gas_price = price * V_ng
    
        # @parameter(distribution=dist.electricity_price_distribution, units='USD/kWh',
        #            element='electricity', baseline=dist.mean_electricity_price)
        # def set_electricity_price(price): 
        #     bst.settings.electricity_price = price
        
        # https://pmc.ncbi.nlm.nih.gov/articles/PMC3947793/
        # https://www.capturemap.no/the-biogenic-co2-breakdown/
        self.BT.CO2_emissions_concentration = 15.0 / 100 # 
        
        if scenario.carbon_source == 'biomass':
            @parameter(units='MT/yr', element='oleochemical', bounds=[20000, 50000], baseline=35000)
            def set_production_capacity(production_capacity):
                self.production_capacity = production_capacity
                
            # @uniform(units='USD/MT', element='Carbon capture', baseline=100) 
            # def set_carbon_capture_cost(price):
            #     self.CC.b = 4.230769230769226 + price
            
            # @parameter(
            #     units='wt %', element='Boiler flue gas', 
            #     distribution=shape.Trunc(shape.Normal(6.3 * 1.55, 0.5 * 1.55), 5.3 * 1.55, 7.3 * 1.55), 
            #     baseline=6.3 * 1.55, 
            #     bounds=(5.3 * 1.55, 7.3 * 1.55)
            # ) 
            # def set_boiler_flue_gas_CO2_content(CO2_content):
            #     self.BT.CO2_emissions_concentration = CO2_content / 100
            
            @system.add_specification(simulate=True)
            def adjust_production_capacity():
                capacity = self.production_capacity / system.operating_hours * 1000 # kg / hr
                self.system.simulate()
                self.system.rescale(self.AcOH_media, capacity / self.product.F_mass) 
        elif scenario.carbon_source == 'BFG':
            total_pig_iron_produced_US = 22.3e6 # MT / yr
            N_facilities_US = 12
            BFG_per_ton_pig_iron = 2.5 # Blast furnace gas (2.5 to 3.5 BFG / steel by wt)
            pig_iron_per_facility = total_pig_iron_produced_US / N_facilities_US # MT / yr
            BFG_per_facility = pig_iron_per_facility * BFG_per_ton_pig_iron  
            
            @parameter(units='MT/yr', element='Flue gas',
                       bounds=[BFG_per_facility / 30, BFG_per_facility / 15]) # Only a fraction is used to prevent TCI > 1 billion
            def set_flue_gas_processing_capacity(processing_capacity):
                self.flue_gas.imass[scenario.carbon_source] = processing_capacity / system.operating_hours * 1000 # kg / hr
        else:
            # Must be glucose, but no uncertainty is accounted for at this level
            pass
            
        if scenario.carbon_capture:
            baseline_length_to_diameter = 8
        else:
            baseline_length_to_diameter = 12
        
        optimized_parameter = model.optimized_parameter
        # if dewatering:
        #     @optimized_parameter(bounds=[5, 20], baseline=12, name='extractor stages')
        #     def set_stages(N):
        #         N = int(N)
        #         self.extractor.N_stages = N
            
        #     @optimized_parameter(bounds=[1, 2], baseline=1.5, name='solvent to feed ratio')
        #     def set_solvent_to_feed(solvent_feed_ratio):
        #         self.extractor.solvent_feed_ratio = solvent_feed_ratio
            
        #     @optimized_parameter(bounds=[1.01, 2.5], baseline=1.2, element=self.extract_distiller)
        #     def extract_distiller_reflux(k):
        #         self.extract_distiller.k = k
            
        #     @optimized_parameter(bounds=[90, 99.9], baseline=95, units='%', element=self.extract_distiller)
        #     def extract_distiller_heavy_key_recovery(Hr):
        #         self.extract_distiller.Hr = Hr / 100
            
        #     @optimized_parameter(bounds=[1.01, 2.5], baseline=1.2, element=self.raffinate_distiller)
        #     def raffinate_distiller_reflux(k):
        #         self.raffinate_distiller.k = k
            
        #     @optimized_parameter(bounds=[80, 99.9], baseline=99, units='%', element=self.raffinate_distiller)
        #     def raffinate_distiller_light_key_recovery(Lr):
        #         self.raffinate_distiller.Lr = Lr / 100
            
        #     @optimized_parameter(bounds=[80, 99.9], baseline=99, units='%', element=self.raffinate_distiller)
        #     def raffinate_distiller_heavy_key_recovery(Hr):
        #         self.raffinate_distiller.Hr = Hr / 100
            
        #     @optimized_parameter(bounds=[0, 1], baseline=0.5, element=self.extract_heater)
        #     def extract_heater_vapor_fraction(V):
        #         self.extract_heater.V = V
        
        @optimized_parameter(bounds=[2, 12], baseline=baseline_length_to_diameter, element='AcOH bioreactor', name='length to diameter')
        def set_AcOH_bioreactor_length_to_diameter(length_to_diameter):
            if scenario.carbon_source == 'glucose': return
            self.AcOH_production.length_to_diameter = length_to_diameter
        
        @optimized_parameter(bounds=[2, 12], baseline=baseline_length_to_diameter, element='oleochemical bioreactor', name='length to diameter')
        def set_oleochemical_bioreactor_length_to_diameter(length_to_diameter):
            self.oleochemical_production.length_to_diameter = length_to_diameter
        
        @optimized_parameter(bounds=[0.2, 0.6], baseline=0.5, element='oleochemical bioreactor', name='agitation power')
        def set_oleochemical_bioreactor_agitation_power(kW_per_m3):
            self.oleochemical_production.kW_per_m3 = kW_per_m3
        
        if not scenario.dewatering: self.ethyl_acetate = bst.MockStream('ethyl_acetate')
        
        chemicals = bst.settings.chemicals
        self.credited_carbon_intake = lambda: (
            self.product.get_atomic_flow('C') * chemicals.CO2.MW * system.operating_hours
        )
        self.system.define_process_impact(
            key=GWPkey,
            name='Credited carbon intake',
            basis='kg',
            inventory=self.credited_carbon_intake,
            CF=-1.,
        )
        key = scenario.biomass
        try:
            self.BT.fuel.set_CF(GWPkey, CFs[key.capitalize()])
        except:
            breakpoint()
        if key == 'miscanthus':
            price_ub = 59 / 907.185 * 0.8 # https://www.sciencedirect.com/science/article/pii/S096195340700205X
            price_lb = 61.98 / 907.185 * 0.8 # https://farmdoc.illinois.edu/fast-tools/biomass-crop-budget-tool-miscanthus-and-switchgrass
        elif key == 'cornstover':
            price_lb = 59 / 907.185 * 0.8 # Humbird NREL 2011 cellulosic ethanol
            price_ub = 64.96 / 907.185 * 0.8 # https://www.extension.purdue.edu/extmedia/ec/re-3-w.pdf
        else:
            raise ValueError('invalid carbon source')
        self.BT.fuel.price = price_baseline = 0.5 * (price_lb + price_ub)
        
        @parameter(units='USD/MT', element='biomass', 
                   bounds=[price_lb * 1000, price_ub * 1000], 
                   baseline=price_baseline * 1000)
        def set_biomass_price(price):
            self.BT.fuel.price = price / 1000
        
        # bst.settings.electricity_price = 0.060 # Maryland solar REC (renewable energy certificates) # https://escholarship.org/uc/item/80n4q8xc
        if scenario.carbon_source != 'glucose':
            self.hydrogen.set_CF(GWPkey, CFs['H2'])
        self.hexane.set_CF(GWPkey, CFs['Hexane'])
        self.ethyl_acetate.set_CF(GWPkey, CFs['Ethyl acetate'])
        
        @metric(units='USD/kg')
        def MSP():
            return self.tea.solve_price(self.product)
        
        @metric(units='kg*CO2e/kg')
        def carbon_intensity():
            return (
                self.system.get_net_impact(GWPkey)
                / self.system.get_mass_flow(self.product)
            )
        
        @metric(units='10^6 USD')
        def TCI():
            return self.tea.TCI / 1e6
        
        @metric(units='% theoretical')
        def product_yield_to_hydrogen():
            if scenario.carbon_source == 'glucose': return 0
            return self.product.get_atomic_flow('H') / self.hydrogen.get_atomic_flow('H')
        
        @metric(units='10^3 MT/yr')
        def biomass_burned(): # 825 MT / y for NREL's cornstover model
            return self.BT.fuel.F_mass * self.system.operating_hours / 1e6
        
        @metric(units='10^3 MT/yr')
        def hydrogen_consumption(): 
            if scenario.carbon_source == 'glucose': return 0
            return self.hydrogen.F_mass * self.system.operating_hours / 1e6
        
        @metric(units='10^3 MT/yr')
        def glucose_consumption(): 
            if scenario.carbon_source == 'glucose':
                glucose = self.feedstock.imass['Glucose'] + self.seedtrain_feed.imass['Glucose']
            elif scenario.glucose_growth:
                glucose = self.seedtrain_feed.imass['Glucose']
            else:
                return 0
            return glucose * self.system.operating_hours / 1e6
        
        @metric(units='kWh/kg-H2')
        def electricity_demand(): 
            return self.system.get_electricity_consumption() / (self.product.F_mass * self.system.operating_hours)
        
        if scenario.carbon_source != 'biomass':
            @metric(units='10^3 MT/yr', element='oleochemical')
            def production_capacity():
                return self.system.get_mass_flow(self.product) / 1e6
        
        for i in model.parameters:
            if i.distribution is None: i.distribution = cp.Uniform(*i.bounds)
        key = scenario.carbon_capture
        if key in self.pseudo_optimal_design_decisions:
            optimal_values = self.pseudo_optimal_design_decisions[key]
            for i, j in zip(model.optimized_parameters, optimal_values): 
                i.baseline = j
                i.setter(j)
        else:
            self.pseudo_optimal_design_decisions[key] = self.optimize().x
        
        match scenario.hydrogen_price:
            case 'min':
                set_H2_price.active = False
                set_H2_price.baseline = set_H2_price.bounds[0]
            case 'baseline':
                set_H2_price.active = False
            case 'max':
                set_H2_price.active = False
                set_H2_price.baseline = set_H2_price.bounds[1]
            case 'all':
                pass
            case _:
                raise ValueError('invalid hydrogen price')
                
        self.fermentation_parameters = (
            set_oleochemical_bioreactor_yield,
            set_oleochemical_specific_yield,
            set_oleochemical_productivity,
            set_oleochemical_titer,
            set_AcOH_productivity,
            set_AcOH_titer,
        )
        match scenario.fermentation_performance:
            case 'potential yield':
                for i in (set_oleochemical_bioreactor_yield,
                          set_oleochemical_specific_yield):
                    i.active = False
                    i.baseline = i.bounds[1]
            case 'baseline yield':
                set_oleochemical_bioreactor_yield.active = False
                set_oleochemical_specific_yield.active = False
            case 'all':
                pass
            case _:
                raise ValueError('invalid hydrogen price')
        
        return model

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
            'Recycled CO2': comps(self.oleochemical_production.outs[0], ['CO2', 'H2O', 'O2', 'N2']),
            'H2': total(self.system.flowsheet.stream.hydrogen),
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
    
    def MSP_contributions(self):
        total_cost = self.MSP() * self.product.F_mass * self.tea.operating_hours
        cost = self.tea.total_production_cost(self.product)
        f = cost / total_cost
        return {
            'OPEX': f,
            'CAPEX': 1 - f,
        }
    
    def to_experimental_conditions(self):
        if 'acetate' in self.scenario.name:
            # -Titer = 600 mg/L alcohols (dodecanol)
            # -Yield = 0.11 g alcohols (dodecanol) / g acetate (0.265 max)
            # -Specific yield = 0.49 g alcohols / g biomass
            # -Productivity = 0.01 g alcohols / L / h, or 0.008 g alcohols / gDCW/h
            self.set_oleochemical_titer.setter(0.6)
            self.set_oleochemical_bioreactor_yield.setter(42)
            self.set_oleochemical_specific_yield.setter(0.49)
            self.set_oleochemical_productivity.setter(0.01)
            self.system.simulate()
        elif self.scenario.name == 'glucose':
            # -Titer = 1,551 mg/L
            # -Yield = 0.126 g alcohols / g glucose (0.32 max)
            # -Specific yield = 1.5 g alcohol / gDCW
            # -Productivity = 0.016 g alcohol / gDCW / h 
            self.set_oleochemical_titer.setter(1.55)
            self.set_oleochemical_bioreactor_yield.setter(39)
            self.set_oleochemical_specific_yield.setter(1.5)
            self.set_oleochemical_productivity.setter(0.016)
            self.system.simulate()
        else:
            raise NotImplementedError(f'experimental conditions for scenario {self.scenario.name!r}')
    
    def H2_price_breakeven_configurations(self):
        if self.scenario.glucose_growth:
            br_glucose = self
            if hasattr(self, 'acetate_growth_biorefinery'):
                br_acetate = self.acetate_growth_biorefinery
            else:
                self.acetate_growth_biorefinery = br_acetate = type(self)(scenario=self.scenario, glucose_growth=False)
        else:
            br_acetate = self
            if hasattr(self, 'glucose_growth_biorefinery'):
                br_glucose = self.glucose_growth_biorefinery
            else:
                self.glucose_growth_biorefinery = br_glucose = type(self)(scenario=self.scenario, glucose_growth=True)
        
        H2 = br_acetate.hydrogen.F_mass
        product = br_acetate.product.F_mass
        contribution = H2 / product
        original_price = br_acetate.hydrogen.price
        br_acetate.hydrogen.price = 0
        MSP_acetate = br_acetate.MSP()
        br_acetate.hydrogen.price = original_price
        MSP_glucose = br_glucose.MSP()
        return (MSP_glucose - MSP_acetate) / contribution
    
    def MSP_CI_vs_specific_yield(self):
        specific_yields = [0.7, 1.0, 1.5, 2.0, 2.5, 3, 3.5]
        CIs = []
        MSPs = []
        for specific_yield in specific_yields:
            self.oleochemical_production.specific_yield = specific_yield
            self.system.simulate()
            CIs.append(self.carbon_intensity())
            MSPs.append(self.MSP())
        return CIs, MSPs
    
    def MSP_CI_vs_H2_over_C(self):
        H2_over_Cs = [1.0, 1.5, 2.0]
        CIs = []
        MSPs = []
        CC = []
        for H2_over_C in H2_over_Cs:
            self.AcOH_production.H2_over_C = H2_over_C
            self.system.simulate()
            CIs.append(self.carbon_intensity())
            MSPs.append(self.MSP())
            CC.append(self.credited_carbon_intake())
        return CIs, MSPs, CC
        
Biorefinery.pseudo_optimal_design_decisions = {
    False: np.array([12. , 12. ,  0.6]),
}


class TRYBiorefinery(bst.ProcessModel):
    
    class Scenario:
        name: str = 'H2|Dodecanol'
        substrate: str = 'H2' # Either H2, Glucose, Corn
        product: str = 'Dodecanol'
        glucose_growth: bool = False
        
    @property
    def name(self):
        return self.scenario.name
    
    @classmethod
    def as_scenario(cls, scenario):
        try:
            substrate, product = scenario.split('|')
        except:
            substrate, product = scenario.split('-')
        return cls.Scenario(
            substrate=substrate,
            product=product,
            glucose_growth=False,
            name=scenario,
        )
    
    def create_thermo(self):
        return create_acetate_ester_chemicals()
    
    def create_system(self):
        scenario = self.scenario
        self.chemicals.define_group(
            'Biomass', 
            ['Water',
             'Sucrose',
             'Extract',
             'Acetate',
             'Ash',
             'Lignin',
             'Protein',
             'Glucan',
             'Xylan',
             'Arabinan',
             'Mannan',
             'Galactan'],
            [0.2022,
             0.0062,
             0.1185,
             0.0146,
             0.0399,
             0.1274,
             0.0251,
             0.283,
             0.1579,
             0.0192,
             0.0049,
             0.0012],
            wt=True,
        )
        match scenario.substrate:
            case 'H2':
                br = Biorefinery(glucose_growth=scenario.glucose_growth) # H2 to oleochemical; baseline biorefinery
                br.set_AcOH_titer(90)
                br.set_oleochemical_specific_yield(2.7)
                br.set_H2_price(2)
                system = br.system
            case 'Glucose' | 'Corn':
                system = create_substrate_oleochemical_system(
                    specific_yield=2.7, 
                    substrate='Glucose',
                    product=scenario.product,
                    glucose_growth=scenario.glucose_growth,
                    facilities=True,
                )
                if scenario.substrate == 'Corn':
                    from biorefineries import corn
                    corn.load()
                    corn.corn.set_CF('GWP', 0.2610) # GREET 2023
                    corn.E316.disconnect(join_ends=True)
                    corn.E401.disconnect(join_ends=True)
                    corn.MH103.outs[1].price = 1e-9
                    corn.S1.disconnect()
                    @corn.E402.add_specification
                    def remove_components():
                        corn.E402.run()
                        IDs = ('Water', 'Glucose')
                        slurry = corn.E402.outs[0]
                        mol = slurry.imol[IDs]
                        slurry.empty()
                        slurry.imol[IDs] = mol
                    new_units = [i for i in corn.system.units[:31] if i not in (corn.V403, corn.P404)]
                    self.slurry = system.ins[1] = new_units[-1].outs[0]
                    bst.settings.set_thermo(system.ins[0].thermo)
                    system = bst.System.from_units(units=[*new_units, *system.units])
                elif scenario.substrate == 'Glucose':
                    self.feedstock = feedstock = system.ins[1]
                    feedstock.empty()
                    glucose = 2.84e+04
                    feedstock.imass['Water', 'Glucose'] = [9 * glucose, glucose]
                    self.feedstock = feedstock
                    feedstock.set_CF('GWP', 0.1 * 0.9375) # Sugar GREET 2023
                    feedstock.price = 0.413 * 0.10
                    feedstock.ID = 'feedstock'
                # system.diagram(kind='cluster', number=True)
        load_process_settings()
        self.tea = create_tea(
            system,
        )
        self.tea.steam_power_depreciation = 'MACRS7'
        self.tea.income_tax = 0.21
        system.__class__.strict_convergence = False
        system.set_tolerance(rmol=1e-4, mol=1e-3, maxiter=50, 
                             method='wegstein', subsystems=True)
        return system
    
    def create_model(self):
        system = self.system
        scenario = self.scenario
        model = bst.Model(system)
        parameter = model.parameter
        metric = model.metric
        
        def uniform(baseline, *args, **kwargs):
            bounds = [0.8 * baseline, 1.2 * baseline]
            return parameter(*args, **kwargs, baseline=baseline, bounds=bounds)
        
        # https://tradingeconomics.com/commodity/sugar
        @parameter(units='USD/kg', element='glucose', bounds=[0.215, 0.611]) 
        def set_glucose_price(price): # https://scijournals.onlinelibrary.wiley.com/doi/epdf/10.1002/bbb.1976?saml_referrer
            if scenario.glucose_growth:
                self.seedtrain_feed.price = price * 0.12
            if scenario.substrate == 'Glucose':
                self.feedstock.price = 0.1 * price
        
        @parameter(units='USD/kg', element='oleochemical', bounds=[3, 6],
                   baseline=3, distribution='uniform') # https://www.alibaba.com/product-detail/Factory-direct-sale-DODECYL-ACETATE-CAS_1601041319372.html
        def set_product_price(price):
            self.product.price = price
        
        @parameter(units='g/L', element='Oleochemical production', 
                   bounds=[1, 10], baseline=5, distribution='uniform')
        def set_oleochemical_titer(titer):
            self.oleochemical_production.titer = titer
            
        @parameter(units='g/L/h', element='Oleochemical production',
                   bounds=[0.1, 1.0], baseline=1.0, distribution='uniform')
        def set_oleochemical_productivity(productivity):
            self.oleochemical_production.productivity = productivity
        
        # Start from optimistic 40% for this conceptual analysis
        @parameter(units='% theoretical', name='bioreactor yield', element='Oleochemical production', 
                   bounds=[40, 85], baseline=40, distribution='uniform')
        def set_oleochemical_bioreactor_yield(X):
            self.oleochemical_production.reactions.X[0] = X / 100.
        
        @metric(units='USD/kg')
        def MSP():
            return self.tea.solve_price(self.product)
        
        chemicals = bst.settings.chemicals
        self.credited_carbon_intake = lambda: (
            self.product.get_atomic_flow('C') * chemicals.CO2.MW * system.operating_hours
        )
        if scenario.substrate != 'H2':
            self.system.define_process_impact(
                key=GWPkey,
                name='Credited carbon intake',
                basis='kg',
                inventory=self.credited_carbon_intake,
                CF=-1.,
            )
        self.BT.fuel.set_CF(GWPkey, CFs['Cornstover'])
        if hasattr(self, 'hydrogen'): self.hydrogen.set_CF(GWPkey, CFs['H2'])
        if hasattr(self, 'hexane'): self.hexane.set_CF(GWPkey, CFs['Hexane'])
        if hasattr(self, 'ethyl_acetate'): self.ethyl_acetate.set_CF(GWPkey, CFs['Ethyl acetate'])

        @metric(units='kg*CO2e/kg')
        def carbon_intensity():
            return (
                self.system.get_net_impact(GWPkey)
                / self.system.get_mass_flow(self.product)
            )
        
        return model