# -*- coding: utf-8 -*-
"""
"""
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
from .systems import create_oleochemical_system
from biorefineries.cane.data import price_distributions_2023 as dist
from biorefineries.tea import (
    create_cellulosic_ethanol_tea as create_tea
)
from warnings import catch_warnings
from inspect import signature
from chaospy import distributions as shape

__all__ = (
    'Biorefinery',
    'Scenario',
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
    >>> br = Biorefinery(simulate=False, scenario='all fermentation-glucose growth')
    >>> br.system.simulate()
    >>> br.system.diagram()
    >>> (br.MSP(), br.carbon_intensity(), br.tea.TCI / 1e6)
    (7.30, 2.80, 438.53)
    
    """
    class Scenario:
        name: str = 'all fermentation|glucose growth'
        product: str = 'Dodecanol'
        carbon_capture: bool = False
        dewatering: bool = False
        glucose_growth: bool = True
        carbon_source: str = 'BFG' # Alternatively 'biomass'
        biomass: str = 'cornstover' # Either 'miscanthus' or 'cornstover'Al
        hydrogen_price: str = 'min' # Alternatively 'min' or 'max'
        fermentation_performance: str = 'all' # Alternatively 'min' or 'max'
    
    
    @property
    def name(self):
        return self.scenario.name
    
    @classmethod
    def as_scenario(cls, scenario):
        try:
            fermentation, glucose_growth = scenario.split('|')
        except:
            fermentation, glucose_growth = scenario.split('-')
        if glucose_growth == 'glucose growth':
            glucose_growth = True
        elif glucose_growth == 'acetate growth':
            glucose_growth = False
        else:
            raise ValueError("invalid scenario")
        hydrogen_price = 'min'
        match fermentation:
            case 'all fermentation':
                fermentation_performance = 'all'
            case 'conservative fermentation':
                fermentation_performance = 'min'
            case 'optimistic fermentation':
                fermentation_performance = 'max'
            case _:
                raise ValueError('invalid scenario')
        return Scenario(
            product='Dodecanol',
            carbon_capture=False,
            dewatering=False,
            glucose_growth=glucose_growth,
            carbon_source='BFG',
            biomass='cornstover',
            hydrogen_price=hydrogen_price,
            fermentation_performance=fermentation_performance,
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
        system = create_oleochemical_system(
            carbon_capture=scenario.carbon_capture,
            dewatering=scenario.dewatering,
            glucose_growth=scenario.glucose_growth,
            product=scenario.product,
            carbon_source=scenario.carbon_source,
        )
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
                self.seedtrain_feed.price = price * 0.1
        
        @parameter(units='USD/kg', element='oleochemical', bounds=[3, 6],
                   baseline=3, distribution='uniform') # https://www.alibaba.com/product-detail/Factory-direct-sale-DODECYL-ACETATE-CAS_1601041319372.html
        def set_product_price(price):
            self.product.price = price
        
        @parameter(units='USD/kg', element='H2', bounds=[3, 7],
                   baseline=3, distribution='uniform') # Natural gas 1.5 - 5; Electrolysis 3 - 7
        def set_H2_price(price):
            self.hydrogen.price = price
        
        @parameter(units='g/L', element='AcOH production',
                   bounds=[40, 90], baseline=60, distribution='uniform')
        def set_AcOH_titer(titer):
            self.AcOH_production.titer['AceticAcid'] = titer
        
        @parameter(units='g/L/h', element='AcOH production',
                   bounds=[1, 2], baseline=1.5, distribution='uniform')
        def set_AcOH_productivity(productivity):
            self.AcOH_production.productivity = productivity
        
        @parameter(units='g/L', element='Oleochemical production', 
                   bounds=[10, 50], baseline=30, distribution='uniform')
        def set_oleochemical_titer(titer):
            self.oleochemical_production.titer = titer
            
        @parameter(units='g/L/h', element='Oleochemical production',
                   bounds=[0.1, 2.0], baseline=1.0, distribution='uniform')
        def set_oleochemical_productivity(productivity):
            self.oleochemical_production.productivity = productivity
        
        @parameter(units='% theoretical', name='bioreactor yield', element='Oleochemical production', 
                   bounds=[40, 90], baseline=65, distribution='uniform')
        def set_oleochemical_bioreactor_yield(X):
            self.oleochemical_production.reactions.X[0] = X / 100.
        
        @parameter(units='g_{' + scenario.product + '}/g_{cell}', element='Oleochemical production', 
                   bounds=[0.7, 3.5], baseline=0.7, distribution='uniform')
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
            raise ValueError('invalid carbon source')
            
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
        self.BT.fuel.set_CF(GWPkey, CFs[key.capitalize()])
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
        
        @metric(units='wt %')
        def product_yield_to_biomass():
            try: 
                return self.product.F_mass / self.BT.fuel.F_mass
            except:
                return 0
        
        @metric(units='% theoretical')
        def product_yield_to_hydrogen():
            return self.product.get_atomic_flow('H') / self.hydrogen.get_atomic_flow('H')
        
        @metric(units='10^3 MT/yr')
        def biomass_burned(): # 825 MT / y for NREL's cornstover model
            return self.BT.fuel.F_mass * self.system.operating_hours / 1e6
        
        @metric(units='10^3 MT/yr')
        def hydrogen_consumption(): 
            return self.hydrogen.F_mass * self.system.operating_hours / 1e6
        
        @metric(units='kWh/kg-H2')
        def electricity_demand(): 
            return self.system.get_electricity_production() / (self.product.F_mass * self.system.operating_hours)
        
        if scenario.carbon_source != 'biomass':
            @metric(units='MT/yr', element='oleochemical')
            def production_capacity():
                return self.system.get_mass_flow(self.product) / 1e3
        
        for i in model.parameters:
            if i.distribution is None: i.distribution = cp.Uniform(*i.bounds)
        key = (scenario.carbon_capture, scenario.dewatering)
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
            case 'max':
                set_H2_price.active = False
                set_H2_price.baseline = set_H2_price.bounds[1]
            case 'all':
                pass
            case _:
                raise ValueError('invalid hydrogen price')
                
        fermentation_parameters = (
            set_oleochemical_bioreactor_yield,
            set_oleochemical_specific_yield,
            set_oleochemical_productivity,
            set_oleochemical_titer,
            set_AcOH_productivity,
            set_AcOH_titer,
        )
        match scenario.fermentation_performance:
            case 'min':
                for i in fermentation_parameters:
                    i.active = False
                    i.baseline = i.bounds[0]
            case 'max':
                for i in fermentation_parameters:
                    i.active = False
                    i.baseline = i.bounds[1]
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
        
Scenario = Biorefinery.Scenario
Biorefinery.pseudo_optimal_design_decisions = {
    (False, False): np.array([12. , 12. ,  0.6]),
}