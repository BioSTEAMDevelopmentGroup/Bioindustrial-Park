#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

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
from biorefineries.TAL.system_ethyl_esters import TAL_sys, TAL_tea, u, s, unit_groups, spec, price

# get_annual_factor = lambda: TAL_tea._annual_factor


get_annual_factor = lambda: TAL_tea.operating_days*24 # hours per year

_kg_per_ton = 907.18474


system_feeds = [i for i in TAL_sys.feeds if i.price]
system_products = [i for i in TAL_sys.products if i.price]
    
# gypsum = find.stream.gypsum
# system_products.append(gypsum)

baseline_yield, baseline_titer, baseline_productivity =\
    spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity

# %% 

# =============================================================================
# Overall biorefinery metrics
# =============================================================================

feedstock = s.feedstock
Mixed_esters = s.Mixed_esters
ethanol_fresh = s.ethanol_fresh
CSL_fresh = s.CSL_fresh


R302 = u.R302
R303 = u.R303

AC401 = u.AC401

# U401 = u.U401

F401 = u.F401
M403 = u.M403

# F402 = u.F402
# H403 = u.H403

# H402 = u.H402

T620 = u.T620

BT = u.BT701

_feedstock_factor = feedstock.F_mass / (feedstock.F_mass-feedstock.imass['Water'])
# Minimum selling price of TAL stream
def get_MSP():
    for i in range(3):
        Mixed_esters.price = TAL_tea.solve_price(Mixed_esters)
    # return Mixed_esters.price*Mixed_esters.F_mass/sum(Mixed_esters.imass['Octyl_5_hydroxyhexanoate','Octyl_3_5_dihydroxyhexanoate', 'DHL'])
    return Mixed_esters.price

# Mass flow rate of TAL stream
get_yield = lambda: sum(Mixed_esters.imass['Ethyl_5_hydroxyhexanoate','Ethyl_3_5_dihydroxyhexanoate', 'DHL'])*get_annual_factor()/1e6
# Purity (%) of TAL in the final product
get_purity = lambda: sum(Mixed_esters.imass['Ethyl_5_hydroxyhexanoate','Ethyl_3_5_dihydroxyhexanoate', 'DHL'])/Mixed_esters.F_mass
# Adjust for purity
get_adjusted_MSP = lambda: get_MSP() / get_purity()
get_adjusted_yield = lambda: get_yield() * get_purity()
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: sum(Mixed_esters.imol['Ethyl_5_hydroxyhexanoate','Ethyl_3_5_dihydroxyhexanoate', 'DHL']) \
    /(R302.outs[0].imol['TAL'])
get_overall_TCI = lambda: TAL_tea.TCI/1e6
# Annual operating cost, note that AOC excludes electricity credit
get_overall_AOC = lambda: TAL_tea.AOC/1e6
get_material_cost = lambda: TAL_tea.material_cost/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: TAL_tea.sales/1e6
# System power usage, individual unit power usage should be positive
excess_power = lambda: (TAL_sys.power_utility.production-TAL_sys.power_utility.consumption)
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
           Metric('Annual operating cost', get_overall_AOC, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual material cost', get_material_cost, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual product sale', get_annual_sale, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr', 'Biorefinery')
           ]


# To see if TEA converges well for each simulation
get_NPV = lambda: TAL_tea.NPV
metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))

# for ug in unit_groups:
#     ug_metrics = ug.metrics
#     metrics.append(Metric(ug.name, ug.metrics[0], '10^6 $', 'Installed cost'))
#     metrics.append(Metric(ug.name, ug.metrics[5], 'USD/h', 'Material cost'))
#     metrics.append(Metric(ug.name, ug.metrics[1], 'GJ/h', 'Cooling duty'))
#     metrics.append(Metric(ug.name, ug.metrics[2], 'GJ/h', 'Heating duty'))
#     metrics.append(Metric(ug.name, ug.metrics[3], 'MW', 'Electricity usage'))

metrics_labels_dict = {
    'Installed cost':(0, '10^6 $'), 
    'Material cost':(5,'USD/h'), 
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

model = TAL_model = Model(TAL_sys, metrics)

param = model.parameter

def baseline_uniform(baseline, ratio):
    return shape.Uniform(baseline*(1-ratio), baseline*(1+ratio))

def baseline_triangle(baseline, ratio):
    return shape.Triangle(baseline*(1-ratio), baseline, baseline*(1+ratio))

# A fake parameter serving as a "blank" in sensitivity analysis to capture
# fluctuations due to converging errors
D = baseline_uniform(1, 0.1)
@param(name='Blank parameter', element=feedstock, kind='coupled', units='',
       baseline=1, distribution=D)
def set_blank_parameter(anything):
    # This does nothing
    feedstock.T = feedstock.T

#%% ######################## TEA parameters ########################
# U101 = SSCF.U101
# D = baseline_uniform(2205, 0.1)
# @param(name='Feedstock flow rate', element=feedstock, kind='coupled', units='dry-ton/day',
#        baseline=2205, distribution=D)
# def set_feedstock_flow_rate(rate):
#     feedstock.mass *= rate / U101._cached_flow_rate
#     U101._cached_flow_rate = rate

D = shape.Triangle(0.84, 0.9, 0.96)
@param(name='Plant uptime', element='TEA', kind='isolated', units='%',
       baseline=0.9, distribution=D)
def set_plant_uptime(uptime):
    TAL_tea.operating_days = 365. * uptime


# D = baseline_triangle(1, 0.25)
# @param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
#         baseline=1, distribution=D)
# def set_TCI_ratio(new_ratio):
#     old_ratio = TAL_tea._TCI_ratio_cached
#     # old_ratio = TAL_no_CHP_tea._TCI_ratio_cached
#     for unit in TAL_sys.units:
#         if hasattr(unit, 'cost_items'):
#             for item in unit.cost_items:
#                 unit.cost_items[item].cost /= old_ratio
#                 unit.cost_items[item].cost *= new_ratio
#     TAL_tea._TCI_ratio_cached = new_ratio
#     # TAL_no_CHP_tea._TCI_ratio_cached = new_ratio


# Only include materials that account for >5% of total annual material cost,
# enzyme not included as it's cost is more affected by the loading (considered later)
D = shape.Triangle(0.8*price['Glucose'], price['Glucose'], 1.2*price['Glucose'])
@param(name='Feedstock unit price', element='TEA', kind='isolated', units='$/dry-kg',
       baseline=price['Glucose'], distribution=D)
def set_feedstock_price(f_price):
    feedstock.price = f_price / _feedstock_factor


D = shape.Triangle(0.198, 0.2527, 0.304)
@param(name='Natural gas unit price', element='TEA', kind='isolated', units='$/kg',
        baseline=0.2527, distribution=D)
def set_natural_gas_price(gas_price):
    BT.natural_gas_price = gas_price

H2_fresh = s.H2_fresh
baseline_H2_fresh_price = H2_fresh.price
D = shape.Triangle(0.8*baseline_H2_fresh_price, baseline_H2_fresh_price, 1.2*baseline_H2_fresh_price)
@param(name='H2 price', element='TEA', kind='isolated', units='$/kg',
        baseline=baseline_H2_fresh_price, distribution=D)
def set_H2_price(H2_fresh_price):
    H2_fresh.price = H2_fresh_price
    

PdC_fresh = s.PdC_fresh
baseline_PdC_fresh_price = PdC_fresh.price
D = shape.Triangle(0.8*baseline_PdC_fresh_price, baseline_PdC_fresh_price, 1.2*baseline_PdC_fresh_price)
@param(name='PdC price', element='TEA', kind='isolated', units='$/kg',
        baseline=baseline_PdC_fresh_price, distribution=D)
def set_H2_price(PdC_fresh_price):
    PdC_fresh.price = PdC_fresh_price
    
    
D = shape.Triangle(0.067, 0.070, 0.074)
@param(name='Electricity unit price', element='TEA', kind='isolated', units='$/kWh',
       baseline=0.070, distribution=D)
def set_electricity_price(e_price):
    bst.PowerUtility.price = e_price


D = shape.Triangle(0.8*price['Activated carbon'], price['Activated carbon'], 1.2*price['Activated carbon'])
@param(name='Activated carbon unit price', element='TEA', kind='isolated', units='$/ft^3',
       baseline=41, distribution=D)
def set_adsorbent_price(ac_price):
    AC401.adsorbent_cost['Activated carbon'] = ac_price

# 2.2 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal
# based on Annual Energy Outlook (AEO) from Energy Information Adiministration (EIA)
# (https://www.eia.gov/outlooks/aeo/), which is $0.7328/gal and similar to the
# 2.2/(2988/1e3) = $0.736/gal based on a density of 2988 g/gal from H2 Tools
# Lower and upper bounds are $1.37/gal and $2.79/gal, or $0.460/kg and $0.978/kg
D = shape.Triangle(0.460, 0.7328, 0.978)
@param(name='Ethanol unit price', element='TEA', kind='isolated', units='$/kg',
       baseline=0.7328, distribution=D)
def set_ethanol_price(etoh_price):
    ethanol_fresh.price = etoh_price
    

#%% ######################## Feedstock parameters ########################
U101 = u.U101
D = shape.Triangle(29728.5*0.8, 29728.5, 29728.5*1.2)
@param(name='Feedstock capacity', element=U101, kind='coupled', units='kg/h',
       baseline=19301.5, distribution=D)
def set_feedstock_capacity(feedstock_capacity):
    feedstock.F_mass = feedstock_capacity
    
#%% ######################## Conversion parameters ########################

# Fermentation
D = shape.Triangle(5, 10, 15)
@param(name='CSL loading', element=R302, kind='coupled', units='g/L',
       baseline=10, distribution=D)
def set_CSL_loading(loading):
    R302.CSL_loading = loading

R302 = u.R302
# 1e-6 is to avoid generating tiny negative flow (e.g., 1e-14)
D = shape.Triangle(0.9, 0.95, 1-1e-6)
@param(name='Seed train fermentation ratio', element=R302, kind='coupled', units='%',
       baseline=0.95, distribution=D)
def set_ferm_ratio(ratio):
    R303.ferm_ratio = ratio

### Fermentation

D = shape.Triangle(baseline_yield*0.8, baseline_yield, baseline_yield*1.2)
@param(name='TAL yield', element=R302, kind='coupled', units='g/g',
       baseline=baseline_yield, distribution=D)
def set_TAL_yield(TAL_yield):
    # spec.load_specifications(TAL_yield,
    #                           spec.spec_2,
    #                           spec.spec_3)
    spec.spec_1 = TAL_yield

D = shape.Triangle(baseline_titer*0.8, baseline_titer, baseline_titer*1.2)
@param(name='TAL titer', element=R302, kind='coupled', units='g/L',
       baseline=baseline_titer, distribution=D)
def set_TAL_titer(TAL_titer):
    # spec.load_specifications(spec.spec_1,
    #                           TAL_titer,
    #                           spec.spec_3)
    spec.spec_2 = TAL_titer
    
D = shape.Triangle(baseline_productivity*0.8, baseline_productivity, baseline_productivity*1.2)
@param(name='TAL productivity', element=R302, kind='coupled', units='g/L/hr',
       baseline=baseline_productivity, distribution=D)
def set_TAL_productivity(TAL_prod):
    # spec.load_specifications(spec.spec_1,
    #                           spec.spec_2,
    #                           TAL_prod)
    spec.spec_3 = TAL_prod
###

# D = shape.Uniform(0.8*0.005, 0.005, 1.2*0.005)
# @param(name='VitaminA yield', element=R302, kind='coupled', units='% theoretical',
#        baseline=0.005, distribution=D)
# def set_vitaminA_yield(yield_):
#     R302.glucose_to_VitaminA_rxn.X = yield_
#     R303.glucose_to_VitaminA_rxn.X = R303.ferm_ratio*yield_

# D = shape.Uniform(0.8*0.005, 0.005, 1.2*0.005)
# @param(name='VitaminD2 yield', element=R302, kind='coupled', units='% theoretical',
#        baseline=0.005, distribution=D)
# def set_vitaminD2_yield(yield_):
#     R302.glucose_to_VitaminD2_rxn.X = yield_
#     R303.glucose_to_VitaminD2_rxn.X = R303.ferm_ratio*yield_

D = shape.Triangle(0.8*0.05, 0.05, 1.2*0.05)
@param(name='Microbe yield', element=R302, kind='coupled', units='% theoretical',
       baseline=0.05, distribution=D)
def set_microbe_yield(yield_):
    R302.glucose_to_microbe_rxn.X = yield_
    R303.glucose_to_microbe_rxn.X = R303.ferm_ratio*yield_
   

#%%
######################## Separation parameters ########################

D = shape.Triangle(0.0739, 0.0910, 0.2474) # experimental data from Singh group
@param(name='Adsorbent unsaturated capacity', element=AC401, kind='coupled', units='g/g',
       baseline=0.0910, distribution=D)
def set_adsorbent_cap(cap):
    AC401.adsorbent_capacity = cap

D = shape.Uniform(0.4, 0.6) # Seader et al., Table 15.2
@param(name='Adsorbent void volume fraction', element=AC401, kind='coupled', units='L/L',
       baseline=0.5, distribution=D)
def set_adsorbent_vvf(frac):
    AC401.void_fraction = frac

D = shape.Uniform(500, 900) # Seader et al., Table 15.2
@param(name='Adsorbent solid particle density', element=AC401, kind='coupled', units='kg/m^3',
       baseline=700, distribution=D)
def set_adsorbent_solid_rho(rho):
    AC401.rho_adsorbent_solid = rho

D = shape.Triangle(0.8*0.07795, 0.07795, 1.2*0.07795) # experimental data from Singh group
@param(name='Desorption single-wash partition coefficient', element=AC401, kind='coupled', units='(g/L)/(g/L)',
       baseline=0.07795, distribution=D)
def set_desorption_K(K):
    AC401.K = K
    
D = shape.Uniform(0.1, 1.9) # assumed
@param(name='Adsorbent replacement period', element=AC401, kind='coupled', units='y',
       baseline=1., distribution=D)
def set_adsorbent_lifetime(lt):
    AC401._default_equipment_lifetime['Activated carbon'] = lt

D = shape.Uniform(0.05, 0.95) # assumed
@param(name='Regeneration fluid retention in column', element=AC401, kind='coupled', units='L-ethanol/L-void',
       baseline=0.5, distribution=D)
def set_adsorption_ethanol_retention(wr):
    AC401.wet_retention = wr

# D = shape.Uniform(0.01, 0.09) # assumed
# @param(name='Ethanol retention in product after drying', element=F402, kind='coupled', units='g-ethanol/g-TAL',
#         baseline=0.05, distribution=D)
# def set_drying_ethanol_retention_in_product_stream(ethanol_retention_drying):
#     F402.product_ethanol_content = ethanol_retention_drying
    
D = shape.Triangle(0.144945, 0.166880, 0.187718) # experimental data from Singh group
@param(name='TAL solubility in ethanol', element=F401, kind='coupled', units='g-TAL/g-solution',
       baseline=0.166880, distribution=D)
def set_TAL_solubility_ethanol(solubility):
    F401.TAL_solubility_in_ethanol_ww = solubility
    M403.TAL_solubility_in_ethanol_ww = solubility

R401 = u.R401
baseline_mono_di_hydroxy_esters_conversion_ratio = R401.mono_di_hydroxy_esters_conversion_ratio
D = shape.Triangle(0.8*baseline_mono_di_hydroxy_esters_conversion_ratio, 
                   baseline_mono_di_hydroxy_esters_conversion_ratio, 
                   1.2*baseline_mono_di_hydroxy_esters_conversion_ratio) # experimental data from Huber group
@param(name='mono_di_hydroxy_esters_conversion_ratio', element=R401, kind='coupled', units='%',
       baseline=baseline_mono_di_hydroxy_esters_conversion_ratio, distribution=D)
def set_mono_di_hydroxy_esters_conversion_ratio(mono_di_hydroxy_esters_conversion_ratio):
    R401.mono_di_hydroxy_esters_conversion_ratio = mono_di_hydroxy_esters_conversion_ratio
    R401.hydrogenation_rxns[0].X = R401.TAL_to_esters_conversion * mono_di_hydroxy_esters_conversion_ratio
    R401.hydrogenation_rxns[1].X = R401.TAL_to_esters_conversion * (1.-mono_di_hydroxy_esters_conversion_ratio)
    # R401.TAL_to_DHL_rxn.X = 1. - sum([i.X for i in list(R401.hydrogenation_rxns[:5]) if not i==R401.TAL_to_DHL_rxn and i.reactant=='TAL'])



baseline_esters_DHL_conversion_ratio = R401.esters_DHL_conversion_ratio
D = shape.Triangle(0.9*baseline_esters_DHL_conversion_ratio, 
                   baseline_esters_DHL_conversion_ratio, 
                   1.1*baseline_esters_DHL_conversion_ratio) # experimental data from Huber group
@param(name='esters_DHL_conversion_ratio', element=R401, kind='coupled', units='%',
       baseline=baseline_esters_DHL_conversion_ratio, distribution=D)
def set_esters_DHL_conversion_ratio(esters_DHL_conversion_ratio):
    TAL_to_esters_and_DHL_conversion = R401.TAL_to_esters_and_DHL_conversion
    R401.esters_DHL_conversion_ratio = esters_DHL_conversion_ratio
    R401.TAL_to_esters_conversion = TAL_to_esters_conversion =\
        TAL_to_esters_and_DHL_conversion * esters_DHL_conversion_ratio - 5e-5
    R401.TAL_to_DHL_rxn.X =\
        TAL_to_esters_and_DHL_conversion * (1-esters_DHL_conversion_ratio) - 5e-5
    R401.hydrogenation_rxns[0].X = TAL_to_esters_conversion * R401.mono_di_hydroxy_esters_conversion_ratio
    R401.hydrogenation_rxns[1].X = TAL_to_esters_conversion * (1.-R401.mono_di_hydroxy_esters_conversion_ratio)
    
    R401.hydrogenation_rxns[4].X = 5e-5
    

baseline_TAL_to_esters_and_DHL_conversion = R401.TAL_to_esters_and_DHL_conversion
D = shape.Triangle(0.9*baseline_TAL_to_esters_and_DHL_conversion, 
                    baseline_TAL_to_esters_and_DHL_conversion, 
                   1.1*baseline_TAL_to_esters_and_DHL_conversion) # experimental data from Huber group
@param(name='TAL_to_esters_and_DHL_conversion', element=R401, kind='coupled', units='%',
       baseline=baseline_TAL_to_esters_and_DHL_conversion, distribution=D)
def set_TAL_to_esters_and_DHL_conversion(TAL_to_esters_and_DHL_conversion):
    R401.TAL_to_esters_conversion = TAL_to_esters_conversion =\
        TAL_to_esters_and_DHL_conversion * R401.esters_DHL_conversion_ratio - 5e-5
    R401.TAL_to_DHL_rxn.X =\
        TAL_to_esters_and_DHL_conversion * (1-R401.esters_DHL_conversion_ratio) - 5e-5
    R401.hydrogenation_rxns[0].X = TAL_to_esters_conversion * R401.mono_di_hydroxy_esters_conversion_ratio
    R401.hydrogenation_rxns[1].X = TAL_to_esters_conversion * (1.-R401.mono_di_hydroxy_esters_conversion_ratio)
    
    R401.hydrogenation_rxns[4].X = 5e-5
    # R401.TAL_to_DHL_rxn.X = min(0.499, max(1e-6, 1. - sum([i.X for i in R401.TAL_conversion_rxns if i.reactant=='TAL']) +R401.TAL_to_DHL_rxn.X))


S405=u.S405
baseline_PdC_recovery_over_project_period = S405.PdC_recovery_over_project_period
D = shape.Uniform(0.8*baseline_PdC_recovery_over_project_period, 
                   1.2*baseline_PdC_recovery_over_project_period) # assumed
@param(name='Pd|C catalyst recovery over project period', element=S405, kind='coupled', units='%',
       baseline=baseline_PdC_recovery_over_project_period, distribution=D)
def set_catalyst_recovery(PdC_rec):
    S405.PdC_recovery_over_project_period = PdC_rec
    
    
# baseline_catalyst_replacements_per_year = R401.catalyst_replacements_per_year
# D = shape.Uniform(0.5*baseline_catalyst_replacements_per_year, 
#                    # baseline_catalyst_replacements_per_year, 
#                    1.5*baseline_catalyst_replacements_per_year) # assumed
# @param(name='Pd|C catalyst replacement frequency', element=R401, kind='coupled', units='y^-1',
#        baseline=baseline_catalyst_replacements_per_year, distribution=D)
# def set_catalyst_deactivation(catalyst_replacements_per_year):
#     R401.catalyst_replacements_per_year = catalyst_replacements_per_year
    

# D = shape.Uniform(0.01, 0.09) # assumed
# @param(name='Ethanol retention in product after drying', element=U401, kind='coupled', units='g-ethanol/g-product',
#        baseline=0.05, distribution=D)
# def set_drying_ethanol_retention_in_product_stream(ethanol_retention_drying):
#     U401.moisture_content = ethanol_retention_drying
    
  
    
# D = shape.Uniform(0.991, 0.999) # assumed
# @param(name='Crystallization and centrifugation recovery', element=S404, kind='coupled', units='g-recovered/g-influent',
#        baseline=0.995, distribution=D)
# def set_crystallization_and_centrifugation_combined_recovery(TAL_recovery):
#     S404.split = np.array([0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
#            0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
#            0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
#            0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
#            0.  , 0.  , 0.  , 0.  , 0.  , TAL_recovery, 0.  , 0.  , 0.  , 0.  , 0.  ,
#            0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
#            0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ])

#%%
######################## Facility parameters ########################
D = baseline_uniform(0.8, 0.1)
@param(name='boiler efficiency', element=BT, kind='coupled', units='%',
       baseline=0.8, distribution=D)
def set_boiler_efficiency(efficiency):
    BT.boiler_efficiency = efficiency

D = shape.Triangle(0.8*7*24, 7*24, 1.2*7*24)
@param(name='Product TAL storage time', element=BT, kind='coupled', units='h',
       baseline=7*24, distribution=D)
def set_product_storage_time(storage_time):
    T620.tau = storage_time
    

parameters = model.get_parameters()


