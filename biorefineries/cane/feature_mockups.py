# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:44:17 2021

@author: yrc2
"""
from biosteam import MockFeature
from thermosteam.units_of_measure import format_units
from thermosteam.utils import roundsigfigs
import pandas as pd
import numpy as np
import os

(set_juicing_oil_recovery, set_microbial_oil_recovery, set_bagasse_oil_recovery, 
 set_cane_operating_days, set_sorghum_operating_days, 
 set_available_land, set_dry_biomass_yield, set_crude_oil_price, 
 set_baseline_feedstock_price,
 set_cellulosic_ethanol_price, set_advanced_ethanol_price,
 set_biomass_based_diesel_price, set_cellulosic_based_diesel_price,
 set_natural_gas_price, set_electricity_price, 
 set_IRR, set_crude_glycerol_price, set_pure_glycerol_price,
 set_saccharification_reaction_time, set_cellulase_price, 
 set_cellulase_loading, set_reactor_base_cost,
 set_cane_glucose_yield, set_cane_xylose_yield, 
 set_sorghum_glucose_yield, set_sorghum_xylose_yield, 
 set_glucose_to_ethanol_yield, set_xylose_to_ethanol_yield,
 set_cofermentation_ethanol_titer, set_cofermentation_ethanol_productivity,
 set_glucose_to_microbial_oil_yield, set_xylose_to_microbial_oil_yield,
 set_fermentation_microbial_oil_titer, set_fermentation_microbial_oil_productivity,
 set_cane_PL_content, set_sorghum_PL_content, set_cane_FFA_content,
 set_sorghum_FFA_content,  set_cane_oil_content, set_relative_sorghum_oil_content,
 set_TAG_to_FFA_conversion, set_feedstock_GWP, set_methanol_GWP, 
 set_pure_glycerine_GWP, set_cellulase_GWP, set_natural_gas_GWP,
 set_income_tax,
) = all_parameter_mockups = (
    MockFeature('Juicing oil recovery', '%', '-'),
    MockFeature('Microbial oil recovery', '%', '-'),
    MockFeature('Bagasse oil recovery', '%', '-'),
    MockFeature('Cane operating days', 'day/y', '-'),
    MockFeature('Sorghum operating days', 'day/y', '-'),
    MockFeature('Available land', 'ha', 'Feedstock'),
    MockFeature('Dry biomass yield', 'DMT/ha/y', 'Feedstock'),
    MockFeature('Price', 'USD/L', 'Crude oil'),
    MockFeature('Price', 'USD/kg', 'Feedstock'),
    MockFeature('Price', 'USD/L', 'Cellulosic ethanol'),
    MockFeature('Price', 'USD/L', 'Advanced ethanol'),
    MockFeature('Price', 'USD/L', 'Biomass based diesel'),
    MockFeature('Price', 'USD/L', 'Cellulosic based diesel'),
    MockFeature('Price', 'USD/m3', 'Natural gas'),
    MockFeature('Price', 'USD/kWh', 'Electricity'),
    MockFeature('IRR', '%', '-'),
    MockFeature('Price', 'USD/kg', 'Crude glycerol'),
    MockFeature('Price', 'USD/kg', 'Pure glycerine'),
    MockFeature('Reaction time', 'h', 'Saccharification'),
    MockFeature('Price', 'USD/kg', 'Cellulase'),
    MockFeature('Cellulase loading', 'wt. % cellulose', 'Cellulase'),
    MockFeature('Base cost', 'million USD', 'Pretreatment reactor system'),
    MockFeature('Cane glucose yield', '%', 'Pretreatment and saccharification'),
    MockFeature('Sorghum glucose yield', '%', 'Pretreatment and saccharification'),
    MockFeature('Cane xylose yield', '%', 'Pretreatment and saccharification'),
    MockFeature('Sorghum xylose yield', '%', 'Pretreatment and saccharification'),
    MockFeature('Glucose to ethanol yield', '%', 'Cofermenation'),
    MockFeature('Xylose to ethanol yield', '%', 'Cofermenation'),
    MockFeature('Ethanol titer', 'g/L', 'Cofermentation'),
    MockFeature('Ethanol productivity', 'g/L/h', 'Cofermentation'),
    MockFeature('Glucose to microbial oil yield', '%', 'Cofermenation'),
    MockFeature('Xylose to microbial oil yield', '%', 'Cofermenation'),
    MockFeature('Microbial oil titer', 'g/L', 'Fermentation'),
    MockFeature('Microbial oil productivity', 'g/L/h', 'Fermentation'),
    MockFeature('Cane PL content', '% oil', 'Oilcane'),
    MockFeature('Sorghum PL content', '% oil', 'Oilsorghum'),
    MockFeature('Cane FFA content', '% oil', 'Oilcane'),
    MockFeature('Sorghum FFA content', '% oil', 'Oilsorghum'),
    MockFeature('Cane oil content', 'dry wt. %', 'Oilcane'),
    MockFeature('Relative sorghum oil content', 'dry wt. %', 'Oilsorghum'),
    MockFeature('TAG to FFA conversion', '% oil', '-'),
    MockFeature('GWP', 'kg*CO2e/kg', 'Sugarcane'),
    MockFeature('GWP', 'kg*CO2e/kg', 'Methanol'),
    MockFeature('GWP', 'kg*CO2e/kg', 'Pure glycerine'),
    MockFeature('GWP', 'kg*CO2e/kg', 'Cellulase'),
    MockFeature('GWP', 'kg*CO2e/kg', 'Natural gas'),
    MockFeature('Income tax', '%', '-'),
)
     
(MFPP, MESP, MBSP, feedstock_consumption, biodiesel_production, biodiesel_yield, ethanol_production, 
 electricity_production, net_energy_production, natural_gas_consumption, TCI, 
 heat_exchanger_network_error, GWP_economic, GWP_ethanol, GWP_biodiesel, 
 GWP_crude_glycerol, GWP_electricity, GWP_ethanol_displacement, GWP_biodiesel_displacement,
 GWP_biofuel_allocation, GWP_ethanol_allocation,
 GWP_biodiesel_allocation, GWP_crude_glycerol_allocation,
 MFPP_derivative, biodiesel_production_derivative, ethanol_production_derivative, 
 electricity_production_derivative, natural_gas_consumption_derivative, 
 TCI_derivative, GWP_economic_derivative, 
 GWP_ethanol_derivative, GWP_biodiesel_derivative,
 GWP_crude_glycerol_derivative, GWP_electricity_derivative,
 ROI, competitive_biomass_yield, energy_competitive_biomass_yield, IRR
 # competitive_microbial_oil_yield,
 # energy_competitive_microbial_oil_yield,
 ) = all_metric_mockups = (
    MockFeature('MFPP', 'USD/MT', '-'),
    MockFeature('MESP', 'USD/L', '-'),
    MockFeature('MBSP', 'USD/L', '-'),
    MockFeature('Feedstock consumption', 'MT/y', '-'),
    MockFeature('Biodiesel production', 'L/MT', '-'),
    MockFeature('Biodiesel yield', 'L/ha', '-'),
    MockFeature('Ethanol production', 'L/MT', '-'),
    MockFeature('Electricity production', 'kWh/MT', '-'),
    MockFeature('Net energy production', 'GGE/MT', '-'),
    MockFeature('Natural gas consumption', 'm3/MT', '-'),
    MockFeature('TCI', '10^6*USD', '-'),
    MockFeature('Heat exchanger network error', '%', '-'),
    MockFeature('GWP', 'kg*CO2e / USD', 'Economic allocation'),
    MockFeature('Ethanol GWP', 'kg*CO2e / L', 'Economic allocation'),
    MockFeature('Biodiesel GWP', 'kg*CO2e / L', 'Economic allocation'),
    MockFeature('Crude glycerol GWP', 'kg*CO2e / kg', 'Economic allocation'),
    MockFeature('Electricity GWP', 'kg*CO2e / MWh', 'Economic allocation'),
    MockFeature('Ethanol GWP', 'kg*CO2e / L', 'Displacement allocation'),
    MockFeature('Biodiesel GWP', 'kg*CO2e / L', 'Displacement allocation'),
    MockFeature('Biofuel GWP', 'kg*CO2e / GGE', 'Energy allocation'),
    MockFeature('Ethanol GWP', 'kg*CO2e / L', 'Energy allocation'),
    MockFeature('Biodiesel GWP', 'kg*CO2e / L', 'Energy allocation'),
    MockFeature('Crude-glycerol GWP', 'kg*CO2e / kg', 'Energy allocation'),
    MockFeature('MFPP derivative', 'USD/MT', '-'),
    MockFeature('Biodiesel production derivative', 'L/MT', '-'),
    MockFeature('Ethanol production derivative', 'L/MT', '-'),
    MockFeature('Electricity production derivative', 'kWh/MT', '-'),
    MockFeature('Natural gas consumption derivative', 'cf/MT', '-'),
    MockFeature('TCI derivative', '10^6*USD', '-'),
    MockFeature('GWP derivative', 'kg*CO2e / USD', 'Economic allocation'),
    MockFeature('GWP derivative', 'kg*CO2e / L', 'Ethanol'),
    MockFeature('GWP derivative', 'kg*CO2e / L', 'Biodiesel'),
    MockFeature('GWP derivative', 'kg*CO2e / kg', 'Crude glycerol'),
    MockFeature('GWP derivative', 'kg*CO2e / MWh', 'Electricity'),
    MockFeature('ROI', '%', '-'),
    MockFeature('Competitive biomass yield', 'dry MT/ha', 'Feedstock'),
    MockFeature('Energy competitive biomass yield', 'dry MT/ha', 'Feedstock'),
    MockFeature('Breakeven IRR', '%', '-'),
)

tea_monte_carlo_metric_mockups = (
    MFPP, 
    TCI,
    ethanol_production,
    biodiesel_production,
    electricity_production,
    natural_gas_consumption
)

tea_monte_carlo_derivative_metric_mockups = (
    MFPP_derivative, 
    TCI_derivative,
    ethanol_production_derivative,
    biodiesel_production_derivative,
    electricity_production_derivative,
    natural_gas_consumption_derivative,
)

lca_monte_carlo_metric_mockups = (
    GWP_economic,
    GWP_ethanol,
    GWP_biodiesel,
    GWP_electricity,
    GWP_crude_glycerol,
)

lca_monte_carlo_derivative_metric_mockups = (
    GWP_economic_derivative,
    GWP_ethanol_derivative,
    GWP_biodiesel_derivative, 
    GWP_electricity_derivative,
    GWP_crude_glycerol_derivative,
)

def get_YRCP2023_spearman_names(configuration, kind=None):
    from biorefineries.cane import Biorefinery, YRCP2023
    YRCP2023()
    br = Biorefinery('O1', simulate=False)
    name = 'name'
    full_name = 'full_name'
    tea_spearman_labels = {
        br.set_juicing_oil_recovery: name,
        br.set_microbial_oil_recovery: name,
        br.set_bagasse_oil_recovery: name,
        br.set_cane_operating_days: name,
        br.set_available_land: name,
        br.set_dry_biomass_yield: name,
        br.set_crude_oil_price: full_name, 
        br.set_baseline_feedstock_price: 'Baseline feedstock price',
        # br.set_cellulosic_ethanol_price: full_name,
        # br.set_advanced_ethanol_price: full_name,
        br.set_biomass_based_diesel_price: full_name,
        br.set_cellulosic_based_diesel_price: full_name,
        br.set_natural_gas_price: full_name,
        br.set_electricity_price: full_name, 
        br.set_IRR: name,
        br.set_crude_glycerol_price: full_name,
        br.set_pure_glycerol_price: full_name,
        br.set_saccharification_reaction_time: full_name,
        br.set_cellulase_price: full_name, 
        br.set_cellulase_loading: name,
        br.set_reactor_base_cost: ('PTRS base cost', 'MMUSD'),
        br.set_cane_glucose_yield: 'Glucan to glucose yield',
        br.set_cane_xylose_yield: 'Xylan to xylose yield', 
        # br.set_glucose_to_ethanol_yield: ('Glucose to ethanol yield', '% theoretical'),
        # br.set_xylose_to_ethanol_yield: ('Xylose to ethanol yield', '% theoretical'),
        # br.set_cofermentation_ethanol_titer: full_name,
        # br.set_cofermentation_ethanol_productivity: full_name,
        br.set_glucose_to_microbial_oil_yield: name, 
        br.set_xylose_to_microbial_oil_yield: name, 
        br.set_fermentation_microbial_oil_titer: name, 
        br.set_fermentation_microbial_oil_productivity: name, 
        br.set_cane_PL_content: name, 
        br.set_cane_FFA_content: name,
        br.set_cane_oil_content: name, 
        br.set_TAG_to_FFA_conversion: name,
    }
    GWP_spearman_labels = {
        br.set_feedstock_GWP: 'Baseline feedstock GWP',
        br.set_methanol_GWP: full_name, 
        br.set_pure_glycerine_GWP: full_name,
        br.set_cellulase_GWP: full_name,
        br.set_natural_gas_GWP: full_name,
    }
    if configuration == 'O7':
        del GWP_spearman_labels[br.set_cellulase_GWP]
        del GWP_spearman_labels[br.set_natural_gas_GWP]
        del tea_spearman_labels[br.set_cellulosic_based_diesel_price]
        del tea_spearman_labels[br.set_saccharification_reaction_time]
        del tea_spearman_labels[br.set_cellulase_price]
        del tea_spearman_labels[br.set_cellulase_loading]
        del tea_spearman_labels[br.set_reactor_base_cost]
        del tea_spearman_labels[br.set_cane_glucose_yield]
        del tea_spearman_labels[br.set_cane_xylose_yield]
        del tea_spearman_labels[br.set_xylose_to_microbial_oil_yield]
        del tea_spearman_labels[br.set_bagasse_oil_recovery]
    elif configuration == 'O9':
        pass
    else:
        raise ValueError(configuration)
    def with_units(f, name, units=None):
        d = f.distribution
        dname = type(d).__name__
        if units is None: units = f.units
        if dname == 'Triangle':
            distribution = ', '.join([format(j, '.3g')
                                      for j in d._repr.values()])
        elif dname == 'Uniform':
            distribution = ' $-$ '.join([format(j, '.3g')
                                         for j in d._repr.values()])
        return f"{name}\n[{distribution} {format_units(units)}]"
        
    def get_full_name(f):
        a = f.element_name
        if a == 'Cofermentation':
            a = 'Co-Fermentation'
        b = f.name
        if b == 'GWP': 
            return f"{a} {b}"
        else:
            return f"{a} {b.lower()}"
        
    for dct in (GWP_spearman_labels, tea_spearman_labels):
        for i, j in tuple(dct.items()):
            if j == name:
                dct[i.index] = with_units(i, i.name)
            elif j == full_name:
                dct[i.index] = with_units(i, get_full_name(i))
            elif isinstance(j, tuple):
                dct[i.index] = with_units(i, *j)
            elif isinstance(j, str):
                dct[i.index] = with_units(i, j)
            else:
                raise TypeError(str(j))
            del dct[i]
    
    lca_spearman_labels = {
        i: j for i, j in
        tea_spearman_labels.items()
        if not any([k in ''.join(i) for k in ('price', 'cost', 'days', 'land', 'IRR')])
    }
    lca_spearman_labels.update(GWP_spearman_labels)
    if kind == 'TEA':
        return tea_spearman_labels
    elif kind == 'LCA':
        return lca_spearman_labels
    else:
        return tea_spearman_labels, lca_spearman_labels

def get_YRCP2023_correlated_distribution_table():
    from biorefineries.cane import Biorefinery, YRCP2023, results_folder
    YRCP2023()
    br = Biorefinery('O1', simulate=False)
    models = br.price_distribution_module.models
    residuals = br.price_distribution_module.residuals
    scores = br.price_distribution_module.scores
    parameter_residuals = [
        br.set_cellulosic_ethanol_price,
        br.set_advanced_ethanol_price,
        br.set_biomass_based_diesel_price,
        br.set_cellulosic_based_diesel_price,
        br.set_natural_gas_price,
        br.set_electricity_price, 
    ]
    rows = []
    for p in parameter_residuals:
        name = p.element_name
        model = models[name]
        rows.append({
            'Price': f"{name}",
            'Units': p.units,
            'Intercept': model.intercept_,
            'Slope': model.coef_[0],
            'Residual standard error': np.std(residuals[name], ddof=2),
            'R2': round(scores[name], 3),
        })
    table = pd.DataFrame(
        rows, 
        index=list(range(35, len(rows) + 35)),
    )
    table.index.name = '#'
    file = os.path.join(results_folder, 'correlated_distributions.xlsx')
    table.to_excel(file)
    return table

def get_YRCP2023_distribution_table(kind=None):
    from biorefineries.cane import Biorefinery, YRCP2023, results_folder
    YRCP2023()
    br = Biorefinery('O1', simulate=False)
    name = 'name'
    full_name = 'full_name'
    parameters = {
        br.set_juicing_oil_recovery: name,
        br.set_microbial_oil_recovery: name,
        br.set_bagasse_oil_recovery: name,
        br.set_cane_operating_days: name,
        br.set_available_land: name,
        br.set_dry_biomass_yield: name,
        br.set_crude_oil_price: full_name, 
        br.set_baseline_feedstock_price: 'Baseline feedstock price',
        br.set_IRR: name,
        br.set_crude_glycerol_price: full_name,
        br.set_pure_glycerol_price: full_name,
        br.set_saccharification_reaction_time: full_name,
        br.set_cellulase_price: full_name, 
        br.set_cellulase_loading: name,
        br.set_reactor_base_cost: ('PTRS base cost', 'MMUSD'),
        br.set_cane_glucose_yield: 'Glucan to glucose yield',
        br.set_cane_xylose_yield: 'Xylan to xylose yield', 
        br.set_glucose_to_ethanol_yield: ('Glucose to ethanol yield', '% theoretical'),
        br.set_xylose_to_ethanol_yield: ('Xylose to ethanol yield', '% theoretical'),
        br.set_cofermentation_ethanol_titer: full_name,
        br.set_cofermentation_ethanol_productivity: full_name,
        br.set_glucose_to_microbial_oil_yield: name, 
        br.set_xylose_to_microbial_oil_yield: name, 
        br.set_fermentation_microbial_oil_titer: name, 
        br.set_fermentation_microbial_oil_productivity: name, 
        br.set_cane_PL_content: name, 
        br.set_cane_FFA_content: name,
        br.set_cane_oil_content: name, 
        br.set_TAG_to_FFA_conversion: name,
        br.set_feedstock_GWP: 'Baseline feedstock GWP',
        br.set_methanol_GWP: full_name, 
        br.set_pure_glycerine_GWP: full_name,
        br.set_cellulase_GWP: full_name,
        br.set_natural_gas_GWP: full_name,
    }
    
    def with_units(f, name, units=None):
        if units is None: units = f.units
        return name, units
        
    def get_distribution_dict(f):
        d = f.distribution
        dname = type(d).__name__
        values = [roundsigfigs(j, 3) for j in d._repr.values()]
        if dname == 'Triangle':
            return {
                'Shape': 'Triangular',
                'Lower': values[0],
                'Upper': values[2],
                'Mode': values[1],
            }
            
        elif dname == 'Uniform':
            return {
                'Shape': 'Uniform',
                'Lower': values[0],
                'Upper': values[1],
                'Mode': '-',
            }
    
    def get_full_name(f):
        a = f.element_name
        if a == 'Cofermentation':
            a = 'Co-Fermentation'
        b = f.name
        if b == 'GWP': 
            return f"{a} {b}"
        else:
            return f"{a} {b.lower()}"
    
    rows = []
    for i, j in parameters.items():
        if j == name:
            parameter_name, units = with_units(i, i.name)
        elif j == full_name:
            parameter_name, units = with_units(i, get_full_name(i))
        elif isinstance(j, tuple):
            parameter_name, units = with_units(i, *j)
        elif isinstance(j, str):
            parameter_name, units = with_units(i, j)
        else:
            raise TypeError(str(j))
        rows.append({
            'Parameter': parameter_name,
            'Units': units,
            'Baseline': roundsigfigs(i.baseline, 3),
            **get_distribution_dict(i),
        })
    table = pd.DataFrame(
        rows, 
        index=list(
            range(1, len(parameters) + 1)
        )
    )
    file = os.path.join(results_folder, 'distributions.xlsx')
    table.to_excel(file)
    table.index.name = '#'
    return table