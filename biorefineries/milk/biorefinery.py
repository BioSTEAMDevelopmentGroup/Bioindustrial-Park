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
from .chemicals import create_chemicals, create_galacto_oligosaccharide_chemicals
from .process_settings import (
    load_process_settings, GWP as GWPkey,
    characterization_factors as CFs
)
from .systems import create_system
from biorefineries.cane.data import price_distributions_2023 as dist
from biorefineries import cellulosic
from biorefineries.tea import (
    create_cellulosic_ethanol_tea as create_tea,
    create_conventional_ethanol_tea as create_galacto_oligosaccharide_tea,
)
from inspect import signature

__all__ = (
    'GalactoOligoSaccharideProcess',
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

class GalactoOligoSaccharideProcess(bst.ProcessModel):
    """
    Examples
    --------
    >>> from biorefineries import milk
    >>> br = milk.GalactoOligoSaccharideProcess(simulate=False)
    >>> scenario, result = br.baseline()
    >>> result
    
    """
    class Scenario:
        processing_capacity: float = 30e9 * 0.06 * 0.453592e-3 / 100, '# Processing capacity [MT/y]' # 30 billion lb/yr of whey permeate * 0.06 % solids * 0.453592e-3 lb/MT * 1/100 supply available Initial estimate based on https://www.dairyfoods.com/articles/95527-its-time-to-rethink-whey-permeate
        dry_whey_permeate_price: float = 0.23 * 0.453592e-3, '# Price of dry whey permeate [USD/MT]' # Initial estimate of 0.23 USD/lb based on the price of sugar https://tradingeconomics.com/commodity/sugar
        H2SO4_price: float = 0.11, '# Price of H2SO4 [USD/kg]' # Defaults to 2021 93 w % H2SO4 price from https://catcost.chemcatbio.org/materials-library.
        CaCO3_price: float = 0.43, '# Price of CaCO3 [USD/kg]' # Defaults to 2021 price from https://catcost.chemcatbio.org/materials-library.
        CaSO4_price: float = 0, '# Price of CaSO4 [USD/kg]' # Defaults to 0 as a waste co-product.
    
    def create_thermo(self):
        return create_galacto_oligosaccharide_chemicals()
    
    def create_system(self):
        scenario = self.scenario
        operating_hours = 330 * 24
        dry_whey_permeate = bst.Stream(
            'dry_whey_permeate', 
            DryWheyPermeate=(
                1000 * scenario.processing_capacity / operating_hours
            ),
            price=scenario.dry_whey_permeate_price
        )
        dry_whey_storage = bst.StorageTank(ins=dry_whey_permeate, tau=7*24)
        sulfuric_acid = bst.Stream(
            'sulfuric_acid', 
            H2SO4=93,
            Water=7,
            price=scenario.H2SO4_price
        )
        tank = cellulosic.SulfuricAcidTank(ins=sulfuric_acid)
        pump = bst.Pump(ins=tank-0, P=101325 * 10)
        hx = bst.HXutility(ins=pump-0, T=400, V=0)
        mixer = bst.Mixer(ins=[dry_whey_storage-0, hx-0])
        @mixer.add_specification(run=True, impacted_units=[tank])
        def adjust_flow_and_temperature(hx=hx):
            # Adjust flow
            z_target = 0.76 / 2
            z_H2SO4 = sulfuric_acid.imass['H2SO4'] / sulfuric_acid.F_mass
            flow_feed = dry_whey_permeate.F_mass
            # z_H2SO4 * flow_acid + flow_feed = (flow_acid + flow_feed) * z_target
            # flow_acid = ((flow_acid + flow_feed) * z_target - flow_feed) / z_H2SO4
            # flow_acid - flow_acid * z_target / z_H2SO4 = (flow_feed * z_target - flow_feed) / z_H2SO4
            # flow_acid * (1 - z_target / z_H2SO4) = (flow_feed * z_target - flow_feed) / z_H2SO4
            flow_acid = flow_feed * (1 - z_target) / (z_H2SO4 - z_target)
            sulfuric_acid.F_mass = flow_acid
        
            # Adjust temperature
            outlet = sulfuric_acid + dry_whey_permeate
            outlet.T = 70 + 273.15
            H_final = outlet.H
            H_feed = dry_whey_permeate.H
            sulfuric_acid_feed = hx.outs[0]
            sulfuric_acid_feed.H = H_final - H_feed
            hx.T = sulfuric_acid_feed.T
            
        reactor = bst.SinglePhaseReactor(
            ins=mixer-0, 
            reaction=bst.Reaction('Lactose -> GalactoOligosaccharide + Water ', X=0.99, reactant='Lactose'),
            T=70 + 273.15,
            batch=False,
            tau=20 / 60,
        )
        dilution_water = bst.Stream('dilution_water', Water=1)
        pump = bst.Pump(ins=dilution_water, P=101325 * 10)
        hx = bst.HXutility(ins=pump-0, T=300, V=0)
        mixer = bst.Mixer(ins=[reactor-0, hx-0])
        @mixer.add_specification(run=True, impacted_units=[pump])
        def adjust_flow_and_temperature(hx=hx, feed=reactor-0):
            flow_feed = feed.F_mass
            z_target = 0.04
            z_H2SO4 = feed.imass['H2SO4'] / flow_feed
            # z_H2SO4 * flow_feed = (flow_water + flow_feed) * z_target
            # - flow_water * z_target = flow_feed * z_target - z_H2SO4 * flow_feed
            flow_water = flow_feed * (z_target - z_H2SO4) / (-z_target) 
            if flow_water <= 0 :
                dilution_water.empty()
                return
            dilution_water.F_mass = flow_water
            
            # Adjust temperature
            outlet = feed + dilution_water
            outlet.T = 25 + 273.15
            H_final = outlet.H
            H_feed = feed.H
            dilution_water_feed = hx.outs[0]
            dilution_water_feed.T = 290
            if dilution_water_feed.H < H_final - H_feed:
                dilution_water_feed.H = H_final - H_feed
            hx.T = dilution_water_feed.T
            
        hx = bst.HXutility(ins=mixer-0, T=25 + 273, V=0)
        calcium_carbonate = bst.Stream(
            'calcium_carbonate', 
            CaCO3=1,
            price=scenario.CaCO3_price
        )
        calcium_carbonate_storage = bst.StorageTank(ins=calcium_carbonate, tau=7*24)
        mixer = bst.Mixer(ins=[hx-0, calcium_carbonate_storage-0])
        @mixer.add_specification(run=True, impacted_units=[calcium_carbonate_storage])
        def adjust_neutralization(feed=reactor-0):
            calcium_carbonate.imol['CaCO3'] = feed.imol['H2SO4']
        
        reactor = bst.SinglePhaseReactor(
            ins=mixer-0, 
            reaction=bst.Reaction('CaCO3 + H2SO4 -> CO2 + CaSO4 + H2O', X=1, reactant='H2SO4'),
            T=25 + 273.15,
            batch=False,
            tau=30 / 60,
        )
        centrifuge = bst.SolidsCentrifuge(ins=reactor-0, moisture_content=0.5, split=dict(CaSO4=1))
        calcium_sulfate = bst.Stream('calcium_sulfate', price=scenario.CaSO4_price)
        calcium_sulfate_storage = bst.StorageTank(ins=centrifuge-0, outs=calcium_sulfate, tau=7*24)
        dryer = bst.SprayDryer(ins=centrifuge-1, moisture_content=0.90)
        product = bst.Stream('product')
        product_storage = bst.StorageTank(ins=dryer-1, outs=product, tau=7*24)
        
    def create_model(self):
        scenario = self.scenario
        system = self.system
        self.tea = create_galacto_oligosaccharide_tea(
            system, lang_factor=5.03, # Solids/fluid processing plant
        )
        self.tea.income_tax = 0.21
        model = bst.Model(system)
        parameter = model.parameter
        metric = model.metric
        
        @metric(units='USD/kg')
        def MSP():
            return self.tea.solve_price(self.product)
        
        return model

class Biorefinery(bst.ProcessModel):
    """
    Examples
    --------
    >>> from biorefineries import milk
    >>> br = milk.Biorefinery(simulate=False)
    >>> scenario, result = br.baseline()
    >>> result
    MSP [USD/kg]        2.27
    GWP [kg*CO2e/kg]   0.218
    TCI [MMUSD]         73.1
    dtype: float64
    
    >>> from biorefineries import milk
    >>> br = milk.Biorefinery(simulate=False, feed='UFPermeate')
    >>> scenario, result = br.baseline()
    >>> result
    MSP [USD/kg]         2.08
    GWP [kg*CO2e/kg]   0.0435
    TCI [MMUSD]          37.7
    dtype: float64
    
    """
    
    class Scenario:
        feed: str = 'AcidWhey'
        product: str = 'Dodecanol'
    
    def create_thermo(self):
        return create_chemicals()
    
    def create_system(self):
        scenario = self.scenario
        load_process_settings()
        system = create_system(product=scenario.product, feed=scenario.feed)
        self.tea = create_tea(
            system,
        )
        self.tea.income_tax = 0.21
        return system
        
    def create_model(self):
        scenario = self.scenario
        model = bst.Model(self.system)
        parameter = model.parameter
        metric = model.metric
        
        def uniform(baseline, *args, **kwargs):
            bounds = [0.8 * baseline, 1.2 * baseline]
            return parameter(*args, **kwargs, baseline=baseline, bounds=bounds)
        
        match scenario.feed:
            case 'UFPermeate':
                baseline = 0.3 * 0.05 / 0.45359237 # aboud 0.11 USD per kg sugar
                lb = 0.8 * baseline
                ub = 1.2 * baseline
            case 'AcidWhey':
                # TODO: Update based on industry expert
                baseline = 0.055 * 0.05 / 0.45359237
                lb = 0.8 * baseline
                ub = 1.2 * baseline
            case 'SaccharifiedCornSlurry':
                # https://scijournals.onlinelibrary.wiley.com/doi/epdf/10.1002/bbb.1976?saml_referrer
                baseline = 0.33 # USD per kg sugar
                baseline *= 0.21 # USD per kg slurry
                lb = 0.8 * baseline
                ub = 1.2 * baseline
            case 'GlucoseMedia':
                # https://tradingeconomics.com/commodity/sugar
                baseline = 0.50 # USD per kg sugar
                baseline *= 0.207 # USD per kg sugar
                lb = 0.8 * baseline
                ub = 1.2 * baseline
            case _:
                raise ValueError('invalid feed')
        
        @parameter(units='USD/kg', element='feedstock', 
                   baseline=baseline, bounds=[lb, ub]) 
        def set_feedstock_price(price):
            self.feedstock.price = price
        
        @uniform(units='USD/kg', element='hexane', baseline=0.73) 
        def set_hexane_price(price):
            self.hexane.price = price
        
        @parameter(units='USD/kg', element='product', bounds=[3, 6],
                   baseline=3) # https://www.alibaba.com/product-detail/Factory-direct-sale-DODECYL-ACETATE-CAS_1601041319372.html
        def set_product_price(price):
            self.product.price = price
        
        @parameter(units='g/L', element=self.fermentation, 
                   bounds=[1, 10], baseline=5) # baseline is really 0.155, but this is infeasible
        def set_titer(titer):
            self.fermentation.titer = titer
            
        @parameter(units='g/L/h', element=self.fermentation,
                   bounds=[0.1, 1], baseline=0.5) # baseline is really 0.002, but this is infeasible
        def set_productivity(productivity):
            self.fermentation.productivity = productivity
        
        @parameter(units='% theoretical', name='yield', element=self.fermentation, 
                   bounds=[30, 60], baseline=35.2)
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
        
        if scenario.feed == 'AcidWhey':
            # Capacity based on U.S. production in 2015 - https://cen.acs.org/articles/95/i6/Acid-whey-waste-product-untapped.html
            bounds=[1e5, 1e6]
            baseline=5e5
        else: # Everything the same as UFPermeate
            # https://ehlenbachscheese.com/cheese-facts.php
            # https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/whey-permeate
            ub = 0.90 * 500e6 / 1000
            lb = ub / 10
            bounds = [lb, ub]
            baseline = 0.5 * (lb + ub)
        # elif feed == 'UFPermeate':
        #     # https://ehlenbachscheese.com/cheese-facts.php
        #     # https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/whey-permeate
        #     ub = 0.90 * 500e6 / 1000
        #     lb = ub / 10
        #     bounds = [lb, ub]
        #     baseline = 0.5 * (lb + ub)
        # elif feed == 'SaccharifiedCornSlurry':
        #     # TODO: Think about the validity of this assumption
        #     # baseline = 1.20384e6
        #     # ub = 1.20 * baseline
        #     # lb = 0.80 * baseline
        #     # bounds = [lb, ub]
        #     ub = 0.90 * 500e6 / 1000
        #     lb = ub / 10
        #     bounds = [lb, ub]
        #     baseline = 0.5 * (lb + ub)
        # else:
        #     raise ValueError('invalid feed')
        
        @parameter(units='MT/yr', element=self.feedstock, bounds=bounds, baseline=baseline)
        def set_processing_capacity(processing_capacity):
            self.feedstock.F_mass = processing_capacity / self.system.operating_hours * 1000 # kg / hr
            
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
        
        self.system.__class__.strict_convergence = False
        self.system.set_tolerance(
            rmol=1e-4, mol=1e-3, maxiter=50, 
            method='wegstein', subsystems=True
        )
        return model

    def evaluate_sample(self,
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
        return self.model(sample)
        
Scenario = Biorefinery.Scenario