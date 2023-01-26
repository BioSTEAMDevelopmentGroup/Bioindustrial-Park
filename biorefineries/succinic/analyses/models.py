# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:52:07 2023

Modified from the biorefineries constructed in [1], [2], and [3] for the production of
[1] 3-hydroxypropionic acid, [2] lactic acid, and [3] ethanol from lignocellulosic feedstocks

[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040


@author: sarangbhagwat
"""


# %% 

# =============================================================================
# Setup
# =============================================================================

# import numpy as np
import biosteam as bst
from chaospy import distributions as shape
# from biosteam import main_flowsheet as find
from biosteam.evaluation import Model, Metric
# from biosteam.evaluation.evaluation_tools import Setter
from biorefineries.succinic.system_sc import succinic_sys, succinic_tea, succinic_LCA, u, s, unit_groups, unit_groups_dict, spec, price, TEA_breakdown

# get_annual_factor = lambda: succinic_tea._annual_factor


get_annual_factor = lambda: succinic_tea.operating_days*24 # hours per year

_kg_per_ton = 907.18474


system_feeds = [i for i in succinic_sys.feeds if i.price]
system_products = [i for i in succinic_sys.products if i.price]
    
# gypsum = find.stream.gypsum
# system_products.append(gypsum)

baseline_yield, baseline_titer, baseline_productivity =\
    spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity

# %% 

# =============================================================================
# Overall biorefinery metrics
# =============================================================================

feedstock = s.sugarcane
product_stream = s.SuccinicAcid
CSL = s.CSL


R302 = u.R302
R303 = u.R303

BT = u.BT701

_feedstock_factor = feedstock.F_mass / (feedstock.F_mass-feedstock.imass['Water'])
# Minimum selling price of succinic stream
def get_MSP():
    for i in range(3):
        product_stream.price = succinic_tea.solve_price(product_stream)
    # return product_stream.price*product_stream.F_mass/sum(product_stream.imass['Octyl_5_hydroxyhexanoate','Octyl_3_5_dihydroxyhexanoate', 'DHL'])
    return product_stream.price

# Mass flow rate of succinic stream
get_yield = lambda: product_stream.imass['SuccinicAcid']*get_annual_factor()/1e6
# Purity (%) of succinic in the final product
get_purity = lambda: product_stream.imass['SuccinicAcid']/product_stream.F_mass
# Adjust for purity
get_adjusted_MSP = lambda: get_MSP() / get_purity()
get_adjusted_yield = lambda: get_yield() * get_purity()
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: product_stream.imol['SuccinicAcid']\
    /(R302.outs[0].imol['SuccinicAcid', 'CalciumSuccinate'].sum())
get_overall_TCI = lambda: succinic_tea.TCI/1e6

get_overall_installed_cost = lambda: succinic_tea.installed_equipment_cost/1e6

# Annual operating cost, note that AOC excludes electricity credit
get_overall_AOC = lambda: succinic_tea.AOC/1e6
get_material_cost = lambda: succinic_tea.material_cost/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: succinic_tea.sales/1e6
# System power usage, individual unit power usage should be positive
excess_power = lambda: (succinic_sys.power_utility.production-succinic_sys.power_utility.consumption)
electricity_price = bst.PowerUtility.price
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*electricity_price*get_annual_factor())/1e6

metrics = [Metric('Minimum selling price', get_MSP, '$/kg', 'Biorefinery'),
           Metric('Product yield', get_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product purity', get_purity, '%', 'Biorefinery'),
           Metric('Adjusted minimum selling price', get_adjusted_MSP, '$/kg', 'Biorefinery'),
           Metric('Adjusted product yield', get_adjusted_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product recovery', get_recovery, '%', 'Biorefinery'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $', 'Biorefinery'),
           Metric('Total installed equipment cost', get_overall_installed_cost, '10^6 $', 'Biorefinery'),
           Metric('Annual operating cost', get_overall_AOC, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual material cost', get_material_cost, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual product sale', get_annual_sale, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr', 'Biorefinery')
           ]

# To see if TEA converges well for each simulation
get_NPV = lambda: succinic_tea.NPV
metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))


metrics_labels_dict = {
    'Installed cost':(0, '10^6 $'), 
    'Material cost':(4,'USD/h'), 
    'Cooling duty':(1,'GJ/h'), 
    'Heating duty':(2,'GJ/h'), 
    'Electricity usage':(3, 'MW'), 
    }

for m, u_i in metrics_labels_dict.items():
    for ug in unit_groups:
        metrics.append(Metric(ug.name, ug.metrics[u_i[0]], u_i[1], m))
#%%
# =============================================================================
# Construct base model and add parameters
# =============================================================================

model = succinic_model = Model(succinic_sys, metrics)

param = model.parameter

def baseline_uniform(baseline, ratio):
    return shape.Uniform(baseline*(1-ratio), baseline*(1+ratio))

def baseline_triangle(baseline, ratio):
    return shape.Triangle(baseline*(1-ratio), baseline, baseline*(1+ratio))

# A fake parameter serving as a "blank" in sensitivity analysis to capture
# fluctuations due to converging errors
D = baseline_uniform(1, 0.1)
@param(name='Blank parameter', element='TEA', kind='coupled', units='',
       baseline=1, distribution=D)
def set_blank_parameter(anything):
    # This does nothing
    feedstock.T = feedstock.T

#%% ######################## TEA parameters ########################
# U101 = SSCF.U101

D = shape.Triangle(0.5479*0.9, 0.5479, 0.5479*1.1)
@param(name='Plant uptime', element='TEA', kind='isolated', units='%',
       baseline=0.5479, distribution=D)
def set_plant_uptime(uptime):
    succinic_tea.operating_days = 365. * uptime



# Only include materials that account for >5% of annual material cost,
D = shape.Triangle(0.8*price['Feedstock'], price['Feedstock'], 1.2*price['Feedstock'])
@param(name='Feedstock unit price', element='TEA', kind='isolated', units='$/wet-kg',
       baseline=price['Feedstock'], distribution=D)
def set_feedstock_price(f_price):
    # feedstock.price = f_price / _feedstock_factor
    feedstock.price = f_price

F404 = u.F404
D = shape.Triangle(0.198, 0.2527, 0.304)
@param(name='Natural gas unit price', element='TEA', kind='isolated', units='$/kg',
        baseline=0.2527, distribution=D)
def set_natural_gas_price(gas_price):
    BT.natural_gas_price = F404.ins[2].price = gas_price


 
D = shape.Triangle(0.067, 0.070, 0.074)
@param(name='Electricity unit price', element='TEA', kind='isolated', units='$/kWh',
       baseline=0.070, distribution=D)
def set_electricity_price(e_price):
    bst.PowerUtility.price = e_price


makeup_MEA_A401 = s.makeup_MEA_A401
D = shape.Triangle(1.032, 1.426, 1.819)
@param(name='Monoethanolamine unit price', element='TEA', kind='isolated', units='$/kg',
       baseline=1.426, distribution=D)
def set_MEA_price(MEA_price):
    makeup_MEA_A401.price = MEA_price
    
CO2_fermentation = s.CO2_fermentation
CO2_seedtrain = s.CO2_seedtrain
D = shape.Triangle(0.401, 0.409, 0.418)
@param(name='Liquid CO2 unit price', element='TEA', kind='isolated', units='$/kg',
       baseline=0.409, distribution=D)
def set_CO2_price(CO2_price):
    CO2_fermentation.price = CO2_seedtrain.price = CO2_price
    
#%% ######################## Feedstock parameters ########################

D = shape.Triangle(96000*0.8, 96000, 96000*1.2)
@param(name='Feedstock capacity', element='Conversion', kind='coupled', units='kg/h',
       baseline=96000, distribution=D)
def set_feedstock_capacity(feedstock_capacity):
    feedstock.F_mass = feedstock_capacity
    
#%% ######################## Conversion parameters ########################

# Fermentation
D = shape.Triangle(5, 10, 15)
@param(name='CSL loading', element='Conversion', kind='coupled', units='g/L',
       baseline=10, distribution=D)
def set_CSL_loading(loading):
    R302.CSL_loading = loading

R302 = u.R302
# 1e-6 is to avoid generating tiny negative flow (e.g., 1e-14)
D = shape.Triangle(0.81, 0.9, 0.99)
@param(name='Seed train fermentation ratio', element='Conversion', kind='coupled', units='%',
       baseline=0.95, distribution=D)
def set_ferm_ratio(ratio):
    R303.ferm_ratio = ratio

### Fermentation

D = shape.Triangle(baseline_yield*0.8, baseline_yield, baseline_yield*1.2)
@param(name='Succinic acid yield', element='Conversion', kind='coupled', units='g/g',
       baseline=baseline_yield, distribution=D)
def set_succinic_yield(succinic_yield):
    # spec.load_specifications(succinic_yield,
    #                           spec.spec_2,
    #                           spec.spec_3)
    spec.spec_1 = succinic_yield

D = shape.Triangle(baseline_titer*0.8, baseline_titer, baseline_titer*1.2)
@param(name='Succinic acid  titer', element='Conversion', kind='coupled', units='g/L',
       baseline=baseline_titer, distribution=D)
def set_succinic_titer(succinic_titer):
    # spec.load_specifications(spec.spec_1,
    #                           succinic_titer,
    #                           spec.spec_3)
    spec.spec_2 = succinic_titer
    
D = shape.Triangle(baseline_productivity*0.8, baseline_productivity, baseline_productivity*1.2)
@param(name='Succinic acid  productivity', element='Conversion', kind='coupled', units='g/L/hr',
       baseline=baseline_productivity, distribution=D)
def set_succinic_productivity(succinic_prod):
    # spec.load_specifications(spec.spec_1,
    #                           spec.spec_2,
    #                           succinic_prod)
    spec.spec_3 = succinic_prod
###

D = shape.Triangle(0.8*0.066, 0.066, 1.2*0.066)
@param(name='A. succinogenes yield', element='Conversion', kind='coupled', units='% theoretical',
       baseline=0.066, distribution=D)
def set_microbe_yield(yield_):
    R302.glucose_to_microbe_rxn.X = yield_
    R303.glucose_to_microbe_rxn.X = R303.ferm_ratio*yield_
   

baseline_neutralization_safety_factor = R302.neutralization_safety_factor
D = shape.Triangle(0.8*baseline_neutralization_safety_factor, 
                   1.*baseline_neutralization_safety_factor, 
                   1.2*baseline_neutralization_safety_factor)
@param(name='Neutralization or pH control safety factor', element='Conversion', kind='coupled', units='N/A',
       baseline=baseline_neutralization_safety_factor, distribution=D)
def set_safety_factor(safety_factor):
    R302.neutralization_safety_factor = safety_factor

#%%
######################## Separation parameters ########################

gypsum = s.gypsum
D = shape.Uniform(-0.0288, 0.00776)
@param(name='Gypsum unit price', element='TEA', kind='isolated', units='$/kg',
       baseline=0, distribution=D)
def set_gypsum_price(price):
    gypsum.price = price



C401, C402, C403 = u.C401, u.C402, u.C403
baseline_crystallization_output_conc_multiplier = 1.
D = shape.Triangle(0.8*baseline_crystallization_output_conc_multiplier, 
                   baseline_crystallization_output_conc_multiplier,
                   1.2*baseline_crystallization_output_conc_multiplier) # assumed
@param(name='Crystallization output concentration multiplier', element='Separation', kind='coupled', units='N/A',
       baseline=baseline_crystallization_output_conc_multiplier, distribution=D)
def set_multiplier(multiplier):
    C401.output_conc_multiplier = C402.output_conc_multiplier = C403.output_conc_multiplier = multiplier
    

S402, S403, S404 = u.S402, u.S403, u.S404
baseline_crystallization_recovery = S402.recovery
D = shape.Triangle(0.8*baseline_crystallization_recovery, 
                   baseline_crystallization_output_conc_multiplier,
                   1.2*baseline_crystallization_recovery) # assumed
@param(name='Crystallization pressure filter recovery', element='Separation', kind='coupled', units='N/A',
       baseline=baseline_crystallization_recovery, distribution=D)
def set_recovery(recovery):
    S402.recovery = S403.recovery = S404.recovery = recovery
    
    

#%%
######################## Facility parameters ########################
D = baseline_uniform(0.8, 0.1)
@param(name='boiler efficiency', element='Co-heat and power', units='%',
       baseline=0.8, distribution=D)
def set_boiler_efficiency(efficiency):
    BT.boiler_efficiency = efficiency

T601 = u.T601
D = shape.Triangle(0.8*7*24, 7*24, 1.2*7*24)
@param(name='Product succinic storage time', element='Storage', kind='coupled', units='h',
       baseline=7*24, distribution=D)
def set_product_storage_time(storage_time):
    T601.tau = storage_time
    

parameters = model.get_parameters()


