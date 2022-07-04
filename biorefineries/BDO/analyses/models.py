#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for 2,3-Butanediol instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

All units are explicitly defined here for transparency and easy reference

@author: yalinli_cabbi
"""


# %% 

# =============================================================================
# Setup
# =============================================================================

import numpy as np
import biosteam as bst
from chaospy import distributions as shape
from biosteam import main_flowsheet as flowsheet
from biosteam.evaluation import Model, Metric
from biorefineries.BDO.system_MS3 import BDO_sys, BDO_tea, spec, feedstock, MEK, BT, BDO_lca, unit_groups_dict
from biorefineries.BDO.process_settings import CFs

_kg_per_ton = 907.18474
_feedstock_factor = _kg_per_ton / 0.8

# =============================================================================
# Construction base model
# =============================================================================

BDO_model = Model(BDO_sys)
param = BDO_model.parameter

def baseline_uniform(baseline, ratio):
    lb, ub = baseline*(1-ratio), baseline*(1+ratio)
    if lb > ub: ub, lb = lb, ub
    return shape.Uniform(lb, ub)

def baseline_triangle(baseline, ratio):
    lb, mid, ub = baseline*(1-ratio), baseline, baseline*(1+ratio)
    if lb > ub: ub, lb = lb, ub
    return shape.Triangle(lb, mid, ub)


# =============================================================================
# TEA parameters
# =============================================================================

D = shape.Triangle(0.84, 0.9, 0.96)
@param(name='Plant uptime', element='TEA', kind='isolated', units='%',
       baseline=0.9, distribution=D)
def set_operating_days(uptime):
    BDO_tea.operating_days = 365. * uptime

D = shape.Triangle(60, 71.3, 83.7)
@param(name='Feedstock unit price', element='TEA', kind='isolated', units='$/dry-ton',
       baseline=71.3, distribution=D)
def set_feedstock_price(price):
    feedstock.price = price / _feedstock_factor

# Impactful parameters are set to triangular distribution based on literature,
# less important ones are set to ±10% of baseline value
special_price = {
#     stream           distribution     min  mid   max
    # 'feedstock':        ('Triangle',   (60, 71.3, 83.7)),
    # 'CSL_fresh':        ('Uniform',   (0.0673,    0.112)),
    'lime_fresh':       ('Triangle',   (0.160, 0.262, 0.288)),
    # 'ethanol_fresh':    ('Triangle',  (0.460,     0.978)),
    # 'hexanol_fresh':    ('Uniform',   (0.551*0.9,     0.551*1.1)),
    # 'natural_gas':      ('Triangle',  (0.198,     0.304)),
    'gypsum':           ('Uniform',   (-0.0288,   0.00776))
    }



# Prices for boiler_chems, baghouse_bag, and cooling_tower_chems are not included
# as they are tied to BT/CT duties
default_price_streams = ('sulfuric_acid_T201', 'enzyme', 
                         'makeup_water', 'caustic', 'ash_disposal', 'oleyl_alcohol',
                         'makeup_TCP_catalyst', 'makeup_KieCNi_catalyst', 'H2_fresh')
# add isobutanol
def add_stream_price_param(stream, D):
    param(setter=lambda price: setattr(stream, 'price', price),
          name=f'{stream.ID} price', element='TEA', kind='isolated', units='$/kg',
          baseline=stream.price, distribution=D)

for ID in special_price:
    stream = flowsheet.stream[ID]
    distribution, values = special_price[ID]
    if distribution == 'Triangle':
        D = shape.Triangle(*values)
    elif distribution == 'Uniform':
        D = shape.Uniform(*values)
    add_stream_price_param(stream, D)

for ID in default_price_streams:
    stream = flowsheet.stream[ID]
    baseline = stream.price
    D = baseline_triangle(baseline, 0.1)
    add_stream_price_param(stream, D)

D = shape.Triangle(0.067, 0.070, 0.073)
@param(name='Electricity price', element='TEA', kind='isolated', units='$/kWh',
       baseline=0.070, distribution=D)
def set_electricity_price(price): 
    bst.PowerUtility.price = price

natural_gas_price = BT.natural_gas_price
D = shape.Triangle(natural_gas_price*0.9, natural_gas_price, natural_gas_price*1.1)
@param(name='Natural gas price', element='TEA', kind='isolated', units='$/kWh',
       baseline=natural_gas_price, distribution=D)
def set_natural_gas_price(price): 
    BT.natural_gas_price = price
    


# =============================================================================
# LCA parameters
# =============================================================================

D = shape.Uniform(0.09646, 0.12894) # see Feedstock_impacts_YL
@param(name='Feedstock GHTP GWP100 CF', element='LCA', kind='isolated', units='kg-CO2eq. / dry-kg',
       baseline=0.10945, distribution=D)
def set_feedstock_GWP_CF(CF):
    CFs['GWP_CFs']['FGHTP Corn stover'] = CF
    
D = shape.Uniform(1.32576, 1.75741) # see Feedstock_impacts_YL
@param(name='Feedstock GHTP FEC CF', element='LCA', kind='isolated', units='MJeq. / dry-kg',
       baseline=1.68, distribution=D)
def set_feedstock_FEC_CF(CF):
    CFs['FEC_CFs']['FGHTP Corn stover'] = CF
    
# =============================================================================
# Pretreatment parameters
# =============================================================================

M202 = flowsheet.unit.M202
D = shape.Triangle(0.25, 0.3, 0.4)
@param(name='Pretreatment solid loading', element=M202, kind='coupled', units='%', 
       baseline=0.3, distribution=D)
def set_pretreatment_solid_loading(loading): 
    M202.solid_loading = loading

R201 = flowsheet.unit.R201
D = shape.Triangle(0.06, 0.099, 0.12)
@param(name='Pretreatment glucan-to-glucose', element=R201, kind='coupled', units='%',
       baseline=0.099, distribution=D)
def set_R201_glucan_conversion(X): R201.pretreatment_rxns[0].X = X    

D = shape.Triangle(0.8, 0.9, 0.92)
@param(name='Pretreatment xylan-to-xylose', element=R201, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_R201_xylan_conversion(X): R201.pretreatment_rxns[4].X = X        


# =============================================================================
# Conversion parameters
# =============================================================================

M301 = flowsheet.unit.M301
R302 = flowsheet.unit.R302
R303 = flowsheet.unit.R303

D = shape.Triangle(0.175, 0.2, 0.25)
@param(name='Enzymatic hydrolysis solid loading', element=M301, kind='coupled', units='%',
       baseline=0.2, distribution=D)
def set_R301_hydrolysis_solid_loading(loading): M301.solid_loading = loading

D = shape.Triangle(10, 20, 30)
@param(name='Enzyme loading', element=M301, kind='coupled', units='mg/g glucan',
        baseline=20, distribution=D)
def set_M301_enzyme_loading(loading): M301.enzyme_loading = loading

# Enzymatic hydrolysis
D = shape.Triangle(0, 24, 56)
R301 = flowsheet('R301')
@param(name='Enzymatic hydrolysis time', element=R301, kind='coupled', units='hr',
       baseline=24, distribution=D)
def set_R301_hydrolysis_time(tau): R301.tau_saccharification = tau

D = shape.Triangle(0.75, 0.9, 0.948-1e-6)
@param(name='Enzymatic hydrolysis glucan-to-glucose', element=R301, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_R301_glucan_conversion(X): R301.saccharification_rxns[2].X = X

# Fermentation
# D = shape.Triangle(76, 120, 145)
# @param(name='Fermentation time', element=R302, kind='coupled', units='hr',
#        baseline=120, distribution=D)
# def set_R302_fermentation_time(tau): R302.tau_cofermentation = tau


D = shape.Triangle(0.9, 1, 1.1) # Arbitrary +- 10%
@param(name='Productivity', element=R302, kind='coupled', units='g/L/hr',
       baseline=1., distribution=D)
def set_BDO_productivity(productivity):
    spec.spec_3 = productivity
    
D = shape.Triangle(5, 10, 15)
@param(name='CSL loading', element=R301, kind='coupled', units='g/L',
       baseline=10, distribution=D)
def set_CSL_loading(loading): R302.CSL_loading = loading



########################### COMMENT Y,T OUT FOR TRY P_MC ANALYSIS ############################
baseline_R302_BDO_yield = spec.spec_1
D = shape.Triangle(baseline_R302_BDO_yield*0.8, baseline_R302_BDO_yield, baseline_R302_BDO_yield*1.2) # +/- 20% of baseline
@param(name='2,3-BDO yield', element=R302, kind='coupled', units='% theoretical',
        baseline=baseline_R302_BDO_yield, distribution=D)
def set_R302_BDO_yield(X):
    spec.spec_1 = X

D = shape.Triangle(109.9*0.8, 109.9, 109.9*1.2) # +/- 20% of baseline
@param(name='2,3-BDO titer', element=R302, kind='coupled', units='g/L',
        baseline=109.9, distribution=D)
def set_R302_BDO_titer(X):
    spec.spec_2 = X



##############################################################################################

### Only produces acetoine, but its recycled and used up anyway
### Maybe doesn't matter, but TBD
# D = shape.Triangle(0.032, 0.040, 0.048)
# @param(name='Acetic acid and glycerol yield', element=R303, kind='coupled', units='% theoretical',
#         baseline=0.040, distribution=D)
# def set_R301_acetic_acid_yield(X): 
#     # 1e6 is to avoid generating tiny negative flow (e.g., 1e-14) in R301
#     # R302_X = R302.cofermentation_rxns.X
#     ferm_ratio = R303.ferm_ratio
    
#     X1 = min(X, 1-1e-6-R302.glucose_to_BDO_rxn.X-R302.glucose_to_biomass_rxn.X)
#     X2 = min(X, 1-1e-6-R302.xylose_to_BDO_rxn.X-R302.xylose_to_biomass_rxn.X)
    
#     R302.glucose_to_acetic_acid_rxn.X = X1
#     R303.glucose_to_acetic_acid_rxn.X = X1 * ferm_ratio
    
#     R302.xylose_to_acetic_acid_rxn.X = X2
#     R303.xylose_to_acetic_acid_rxn.X = X2 * ferm_ratio
    
#     X1_glycerol = X1 if X1==X else 0.
#     X2_glycerol = X2 if X2==X else 0.  
    
#     R302.glucose_to_glycerol_rxn.X = X1_glycerol
#     R303.glucose_to_glycerol_rxn.X = X1_glycerol * ferm_ratio
    
#     R302.xylose_to_glycerol_rxn.X = X2_glycerol
#     R303.xylose_to_glycerol_rxn.X = X2_glycerol * ferm_ratio

# D = shape.Uniform(0.4, 0.6)
# @param(name='Unfermented sugars routed to CO2 generation (balance is routed to cell mass generation)', element=R303, kind='coupled', units='% theoretical',
#        baseline=0.5, distribution=D)
# def set_sugars_to_CO2_gen(X): 
#    R302.CO2_generation_rxns.X[0] = R302.CO2_generation_rxns.X[1] = X
#    R303.CO2_generation_rxns.X[0] = R303.CO2_generation_rxns.X[1] = X * R303.ferm_ratio
   
# D = shape.Triangle(0.05, 0.07, 0.1)
# @param(name='Innoculum ratio', element=R301, kind='coupled', units='%',
#        baseline=0.07, distribution=D)
# def set_innoculum_ratio(ratio): R301.inoculum_ratio = ratio

# Seed train fermentation yield as a ratio of the main fermenter
# D = baseline_uniform(36, 0.1)
# @param(name='Seed train time', element=R303, kind='coupled', units='hr',
#        baseline=36, distribution=D)
# def set_R303_fermentation_time(tau): R303.tau_batch = tau

# D = shape.Triangle(0.8, 0.9, 1)
# @param(name='Seed train yield', element=R302, kind='coupled', units='% of R301',
#        baseline=0.9, distribution=D)
# def set_R303_ratio(ratio):
#     R302_X = R302.cofermentation_rxns.X
#     R303_X = R303.cofermentation_rxns.X
#     ratio = min(ratio, (1-1e-6-R303_X[2])/(R302_X[0]+R302_X[1]))
#     R303.ferm_ratio = ratio

# =============================================================================
# Separation parameters
# =============================================================================

# S402 = flowsheet.unit.S402
# D = shape.Triangle(0.95, 0.995, 1.)
# @param(name='Gypsum split', element=S402, kind='coupled', units='',
#        baseline=0.995, distribution=D)
# def set_S402_gypsum_split(split):
#     gypsum_index = S402.chemicals.index('Gypsum')
#     S402.split[gypsum_index] = split

# R401 = flowsheet.unit.R401
# D = baseline_triangle(1., 0.1)
# @param(name='Acidulation time', element=R401, kind='coupled', units='hr',
#        baseline=1., distribution=D)
# def set_R401_tau(tau):
#     R401.tau = tau

# R402 = flowsheet.unit.R402
# # D = baseline_triangle(0.95, 0.05)
# D = shape.Triangle(0.64, 0.8, 0.96)
# @param(name='Dehydration conversion', element=R402, kind='coupled', units='',
#        baseline=0.8, distribution=D)
# def set_R402_conversion(X):
#     R402.dehydration_reactions[0].X = X
    
# R403 = flowsheet.unit.R403
# D = shape.Triangle(0.72, 0.8, 0.88)
# @param(name='Hydrolysis conversion', element=R403, kind='coupled', units='%',
#        baseline=0.8, distribution=D)
# def set_R403_conversion_factor(X):
#     R403.hydrolysis_rxns.X[:] = X
    
# =============================================================================
# Separation parameters
# =============================================================================

# Dehydration (BDO to MEK and IBA)
R401 = flowsheet.unit.R401
D = baseline_triangle(0.9, 0.1)
@param(name='dehydration overall 2,3-BDO conversion', element=R401, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_dehydration_overall_conversion(X):
    prev_overall_conversion = float(R401.overall_BDO_conversion)
    R401.overall_BDO_conversion = X
    R401.BDO_to_MEK_rxn.X = X * R401.BDO_to_MEK_selectivity
    R401.BDO_to_MEK_rxn.X = X * R401.BDO_to_IBA_selectivity

D = shape.Triangle(0.813, 0.9, 0.96) # lb: https://doi.org/10.1039/C5RA23251A ; m: https://doi.org/10.1021/acs.iecr.6b03678; ub: https://doi.org/10.1021/i300007a025
@param(name='dehydration 2,3-BDO -> MEK selectivity', element=R401, kind='coupled', units='%',
   baseline=0.9, distribution=D)
def set_dehydration_MEK_selectivity(S):
    R401.BDO_to_MEK_selectivity = S
    R401.BDO_to_MEK_rxn.X = R401.overall_BDO_conversion * S

D = baseline_triangle(0.1, 0.1) # m: https://doi.org/10.1021/acs.iecr.6b03678;
@param(name='dehydration 2,3-BDO -> IBA selectivity', element=R401, kind='coupled', units='%',
       baseline=0.1, distribution=D)
def set_dehydration_IBA_selectivity(S):
    R401.BDO_to_IBA_selectivity = S
    R401.BDO_to_IBA_rxn.X = min(R401.overall_BDO_conversion - 1e-6 - R401.BDO_to_MEK_rxn.X, 
                                R401.overall_BDO_conversion * S)
    
D = baseline_triangle(0.5, 0.1) # m: https://doi.org/10.1021/i300007a025
@param(name='dehydration residence time', element=R401, kind='coupled', units='h',
       baseline=0.5, distribution=D)
def set_dehydration_residence_time(tau):
    R401.tau = tau

D = shape.Triangle(180., 200., 200.*1.1) # m: https://doi.org/10.1021/acs.iecr.6b03678 & https://doi.org/10.1039/C5RA23251A; lb: https://doi.org/10.1021/i300007a025
@param(name='dehydration temperature', element=R401, kind='coupled', units='Celsius',
       baseline=200., distribution=D)
def set_dehydration_temperature(T):
    R401.T = T + 273.15
    
# Hydrogenation (IBA to IBO)
R402 = flowsheet.unit.R402
D = baseline_triangle(0.86, 0.1)
@param(name='hydrogenation IBA -> IBO conversion', element=R402, kind='coupled', units='%',
       baseline=0.86, distribution=D)
def set_hydrogenation_conversion(X):
    R402.hydrogenation_rxns[0].X = X

D = baseline_triangle(16., 0.1) 
@param(name='dehydration residence time', element=R402, kind='coupled', units='h',
       baseline=16., distribution=D)
def set_hydrogenation_residence_time(tau):
    R402.tau = tau

D = baseline_triangle(25., 0.1) 
@param(name='hydrogenation temperature', element=R402, kind='coupled', units='Celsius',
       baseline=25., distribution=D)
def set_hydrogenation_temperature(T):
    R402.T = T + 273.15

# =============================================================================
# Facilities parameters
# =============================================================================

# Facilities are currently not in unit path. thus not set the element here
D = baseline_uniform(0.8, 0.1)
@param(name='BT combustion efficiency', element=BT, kind='coupled', units='%',
       baseline=0.8, distribution=D)
def set_BT_combustion_efficiency(efficiency):
    BT.boiler_efficiency = efficiency

# All parameters
parameters = BDO_model.get_parameters()

# %%

# =============================================================================
# Model to evalute system across feedstock price and carbohydate content
# =============================================================================

get_NPV = lambda: BDO_tea.NPV

def create_feedstock_price_metrics(price):
    def get_price_based_MSP():
        price_per_kg = price / _kg_per_ton * 0.8
        feedstock.price = price_per_kg
        return BDO_tea.solve_price(MEK)
    return [Metric('Minimum selling price', get_price_based_MSP, '$/kg', f'Price={price:.0f} [$/dry-ton]'),
            Metric('Net present value', get_NPV, '$', f'Price={price:.0f} [$/dry-ton]')]

prices = np.linspace(50, 300, 26)
prices = np.concatenate((prices, np.array([71.26])))
feedstock_price_metrics = sum([create_feedstock_price_metrics(price) 
                                for price in prices],[])

def set_feedstock_carbs(carbs_content):
    carbs = ('Glucan', 'Xylan', 'Arabinan', 'Galactan', 'Mannan')
    dry_mass = feedstock.F_mass.copy() - feedstock.imass[('H2O',)].copy()
    old_carbs_mass_total = feedstock.imass[carbs].sum().copy()
    ratio = feedstock.get_normalized_mass(carbs)
    new_carbs_mass = dry_mass * carbs_content * ratio
    feedstock.set_flow(new_carbs_mass, 'kg/hr', carbs)
    mass_diff = new_carbs_mass.sum() - old_carbs_mass_total
    feedstock.imass['Extract'] -= mass_diff
    if any(feedstock.mass < 0):
        raise ValueError(f'Carbohydrate content of {carbs_content*100:.0f}% dry weight is infeasible')

BDO_model_feedstock = Model(BDO_sys, feedstock_price_metrics)

param = BDO_model_feedstock.parameter

# %% Create metrics

metrics = [Metric('Minimum selling price', lambda: BDO_tea.solve_price(MEK), '$/kg'),
           Metric('Total capital investment', lambda: BDO_tea.TCI / 1e6, '10^6 $'),
           Metric('Annual operating cost', lambda: BDO_tea.AOC / 1e6, '10^6 $/yr'),
           Metric('Annual material cost', lambda: BDO_tea.material_cost / 1e6, '10^6 $/yr'),
           Metric('Annual product sale', lambda: BDO_tea.sales / 1e6, '10^6 $/yr'),
           Metric('Annual electricity credit', lambda: min(-BDO_tea.utility_cost, 0) * BDO_tea.operating_hours / 1e6, '10^6 $/yr')
]

# # %% Breakdowns by process groups
# def get_group_heating_demand(group):
#     return sum([sum([hu.duty for hu in unit.heat_utilities if hu.duty*hu.flow>0.]) for unit in group.units])

# def get_group_cooling_demand(group):
#     return sum([sum([hu.duty for hu in unit.heat_utilities if hu.duty*hu.flow<0.]) for unit in group.units])

# # Heating duty

# metrics.extend((Metric('feedstock_group - heating demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow>0.]) for unit in \
#                                           unit_groups_dict['feedstock_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('pretreatment_group - heating demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow>0.]) for unit in \
#                                           unit_groups_dict['pretreatment_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('conversion_group - heating demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow>0.]) for unit in \
#                                           unit_groups_dict['conversion_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('separation_group - heating demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow>0.]) for unit in \
#                                           unit_groups_dict['separation_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('WWT_group - heating demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow>0.]) for unit in \
#                                           unit_groups_dict['WWT_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('HXN_group - heating demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow>0.]) for unit in \
#                                           unit_groups_dict['HXN_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('BT_group - heating demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow>0.]) for unit in \
#                                           unit_groups_dict['BT_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('CT_group - heating demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow>0.]) for unit in \
#                                           unit_groups_dict['CT_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('facilities_no_hu_group - heating demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow>0.]) for unit in \
#                                           unit_groups_dict['facilities_no_hu_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# # Cooling duty

# metrics.extend((Metric('feedstock_group - cooling demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow<0.]) for unit in \
#                                           unit_groups_dict['feedstock_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('pretreatment_group - cooling demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow<0.]) for unit in \
#                                           unit_groups_dict['pretreatment_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('conversion_group - cooling demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow<0.]) for unit in \
#                                           unit_groups_dict['conversion_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('separation_group - cooling demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow<0.]) for unit in \
#                                           unit_groups_dict['separation_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('WWT_group - cooling demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow<0.]) for unit in \
#                                           unit_groups_dict['WWT_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('HXN_group - cooling demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow<0.]) for unit in \
#                                           unit_groups_dict['HXN_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('BT_group - cooling demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow<0.]) for unit in \
#                                           unit_groups_dict['BT_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('CT_group - cooling demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow<0.]) for unit in \
#                                           unit_groups_dict['CT_group'].units])/AA.F_mass,
#                        'MJ/kg'),))

# metrics.extend((Metric('facilities_no_hu_group - cooling demand',
#                        lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
#                                                if hu.duty*hu.flow<0.]) for unit in \
#                                           unit_groups_dict['facilities_no_hu_group'].units])/AA.F_mass,
#                        'MJ/kg'),))  


# # Installed equipment cost

# metrics.extend((Metric('feedstock_group - installed equipment cost',
#                        lambda:unit_groups_dict['feedstock_group'].get_installed_cost(),
#                        '10^6 $'),))

# metrics.extend((Metric('pretreatment_group - installed equipment cost',
#                        lambda:unit_groups_dict['pretreatment_group'].get_installed_cost(),
#                        '10^6 $'),))

# metrics.extend((Metric('conversion_group - installed equipment cost',
#                        lambda:unit_groups_dict['conversion_group'].get_installed_cost(),
#                        '10^6 $'),))

# metrics.extend((Metric('separation_group - installed equipment cost',
#                        lambda:unit_groups_dict['separation_group'].get_installed_cost(),
#                        '10^6 $'),))

# metrics.extend((Metric('WWT_group - installed equipment cost',
#                        lambda:unit_groups_dict['WWT_group'].get_installed_cost(),
#                        '10^6 $'),))

# metrics.extend((Metric('HXN_group - installed equipment cost',
#                        lambda:unit_groups_dict['HXN_group'].get_installed_cost(),
#                        '10^6 $'),))

# metrics.extend((Metric('BT_group - installed equipment cost',
#                        lambda:unit_groups_dict['BT_group'].get_installed_cost(),
#                        '10^6 $'),))

# metrics.extend((Metric('CT_group - installed equipment cost',
#                        lambda:unit_groups_dict['CT_group'].get_installed_cost(),
#                        '10^6 $'),))

# metrics.extend((Metric('facilities_no_hu_group - installed equipment cost',
#                        lambda:unit_groups_dict['facilities_no_hu_group'].get_installed_cost(),
#                        '10^6 $'),))  

# # Power utility demand

# metrics.extend((Metric('feedstock_group - power utility demand',
#                        lambda:unit_groups_dict['feedstock_group'].get_electricity_consumption()/AA.F_mass,
#                        'MW/kg'),))

# metrics.extend((Metric('pretreatment_group - power utility demand',
#                        lambda:unit_groups_dict['pretreatment_group'].get_electricity_consumption()/AA.F_mass,
#                        'MW/kg'),))

# metrics.extend((Metric('conversion_group - power utility demand',
#                        lambda:unit_groups_dict['conversion_group'].get_electricity_consumption()/AA.F_mass,
#                        'MW/kg'),))

# metrics.extend((Metric('separation_group - power utility demand',
#                        lambda:unit_groups_dict['separation_group'].get_electricity_consumption()/AA.F_mass,
#                        'MW/kg'),))

# metrics.extend((Metric('WWT_group - power utility demand',
#                        lambda:unit_groups_dict['WWT_group'].get_electricity_consumption()/AA.F_mass,
#                        'MW/kg'),))

# metrics.extend((Metric('HXN_group - power utility demand',
#                        lambda:unit_groups_dict['HXN_group'].get_electricity_consumption()/AA.F_mass,
#                        'MW/kg'),))

# metrics.extend((Metric('BT_group - power utility demand',
#                        lambda:unit_groups_dict['BT_group'].get_electricity_consumption()/AA.F_mass,
#                        'MW/kg'),))

# metrics.extend((Metric('CT_group - power utility demand',
#                        lambda:unit_groups_dict['CT_group'].get_electricity_consumption()/AA.F_mass,
#                        'MW/kg'),))

# metrics.extend((Metric('facilities_no_hu_group - power utility demand',
#                        lambda:unit_groups_dict['facilities_no_hu_group'].get_electricity_consumption()/AA.F_mass,
#                        'MW/kg'),))  

# # Material cost

# metrics.extend((Metric('feedstock_group - material cost',
#                        lambda:get_material_cost_breakdown()['feedstock_group'],
#                        '$/kg'),))

# metrics.extend((Metric('pretreatment_group - material cost',
#                        lambda:get_material_cost_breakdown()['pretreatment_group'],
#                        '$/kg'),))

# metrics.extend((Metric('conversion_group - material cost',
#                        lambda:get_material_cost_breakdown()['conversion_group'],
#                        '$/kg'),))

# metrics.extend((Metric('separation_group - material cost',
#                        lambda:get_material_cost_breakdown()['separation_group'],
#                        '$/kg'),))

# metrics.extend((Metric('WWT_group - material cost',
#                        lambda:get_material_cost_breakdown()['WWT_group'],
#                        '$/kg'),))

# metrics.extend((Metric('HXN_group - material cost',
#                        lambda:get_material_cost_breakdown()['HXN_group'],
#                        '$/kg'),))

# metrics.extend((Metric('BT_group - material cost',
#                        lambda:get_material_cost_breakdown()['BT_group'],
#                        '$/kg'),))

# metrics.extend((Metric('CT_group - material cost',
#                        lambda:get_material_cost_breakdown()['CT_group'],
#                        '$/kg'),))

# metrics.extend((Metric('facilities_no_hu_group - material cost',
#                        lambda:get_material_cost_breakdown()['facilities_no_hu_group'],
#                        '$/kg'),))
# # # =============================================================================
# # # Capital cost breakdown
# # # =============================================================================

# # def get_installed_equipment_cost(group):
# #     # return lambda: sum(i.installed_equipment_cost for i in unit_groups_dict[system])/1e6
# #     return unit_groups_dict[group].get_installed_cost()/AA.F_mass
# # for system in unit_groups_dict.keys():
# #     # if system == 'feedstock_sys': continue
# #     metrics.extend(
# #         (Metric(system, get_installed_equipment_cost(system), '10^6 $', 'Installed cost'),))

# # # All checks should be ~0
# # check_installed_equipment_cost = \
# #     lambda: sum(get_installed_equipment_cost(system)() 
# #                 for system in unit_groups_dict.keys()) - HP_tea.installed_equipment_cost/1e6
# # metrics.extend((Metric('Check', check_installed_equipment_cost, '10^6 $', 'Installed cost'),))


# # =============================================================================
# # Material cost breakdown
# # =============================================================================

# def get_material_cost(feed):
#     return lambda: feed.price*feed.F_mass*get_annual_factor()/1e6
# for feed in system_feeds:
#     metrics.extend((Metric(feed.ID, get_material_cost(feed), '10^6 $/yr', 'Material cost'),))
# fermentation_lime = find.stream.fermentation_lime
# FGD_lime = find.stream.FGD_lime
# # get_fermentation_lime_ratio = lambda: fermentation_lime.imol['Lime'] \
# #     / (fermentation_lime.imol['Lime']+FGD_lime.imol['Lime']) 
# # S601 = find.unit.S601
# # get_separation_sulfuric_acid_ratio = lambda: S601.outs[1].imol['H2SO4']/S601.ins[0].imol['H2SO4']
# # get_separation_sulfuric_acid_ratio = 
# check_material_cost = lambda: sum(get_material_cost(feed)()
#                                   for feed in system_feeds) - HP_tea.material_cost/1e6

# # metrics.extend((
# #     Metric('Fermentation lime ratio', get_fermentation_lime_ratio, 
# #            '%', 'Material cost'),
# #     # Metric('Separation sulfuric acid ratio', get_separation_sulfuric_acid_ratio, 
# #     #        '%', 'Material cost'),
# #     Metric('Check', check_material_cost, '10^6 $/yr', 'Material cost')))

# def get_product_sale(stream):
#     return lambda: stream.price*stream.F_mass*get_annual_factor()/1e6
# for product in system_products:
#     metrics.extend((Metric(product.ID, get_product_sale(product), '10^6 $/yr', 'Product sale'),))
# check_product_sale= \
#     lambda: sum(get_product_sale(product)() for product in system_products) \
#         - HP_tea.sales/1e6
# metrics.extend((Metric('Check', check_product_sale, '10^6 $/yr', 'Product sale'),))


# get_NPV = lambda: HP_tea.NPV
# metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))

# # To check HXN energy balance error
# # metrics.append(Metric('HXN energy balance error', lambda: HXN.energy_balance_percent_error))

# metrics.extend((Metric('HXN energy balance error', lambda: HXN.energy_balance_percent_error, '%', 'TEA'), ))



##### LCA #####
metrics.extend((
    Metric('Total GWP', lambda: BDO_lca.GWP, 'kg CO2-eq/kg', 'LCA'),
    Metric('Total FEC', lambda: BDO_lca.FEC, 'MJ/kg', 'LCA')
    ))


# Material GWP
# metrics.extend((Metric('GWP - H2SO4',
#                        lambda:get_GWP_by_ID('H2SO4'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('GWP - NaOH',
#                        lambda:get_GWP_by_ID('NaOH'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('GWP - AmmoniumHydroxide',
#                        lambda:get_GWP_by_ID('AmmoniumHydroxide'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('GWP - CalciumDihydroxide',
#                        lambda:get_GWP_by_ID('CalciumDihydroxide'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('GWP - Hexanol',
#                        lambda:get_GWP_by_ID('Hexanol'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('GWP - Enzyme',
#                        lambda:get_GWP_by_ID('Enzyme'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('GWP - TiO2',
#                        lambda:get_GWP_by_ID('TiO2'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('GWP - CSL',
#                        lambda:get_GWP_by_ID('CSL'),
#                        'kg CO2/kg', 'LCA'),))

# # Material FEC
# metrics.extend((Metric('FEC - H2SO4',
#                        lambda:get_FEC_by_ID('H2SO4'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - NaOH',
#                        lambda:get_FEC_by_ID('NaOH'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - AmmoniumHydroxide',
#                        lambda:get_FEC_by_ID('AmmoniumHydroxide'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - CalciumDihydroxide',
#                        lambda:get_FEC_by_ID('CalciumDihydroxide'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - Hexanol',
#                        lambda:get_FEC_by_ID('Hexanol'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - Enzyme',
#                        lambda:get_FEC_by_ID('Enzyme'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - TiO2',
#                        lambda:get_FEC_by_ID('TiO2'),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - CSL',
#                        lambda:get_FEC_by_ID('CSL'),
#                        'kg CO2/kg', 'LCA'),))
# # Natural gas
# metrics.extend((Metric('GWP - natural gas',
#                        lambda:get_ng_GWP(),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - natural gas',
#                        lambda:get_ng_FEC(),
#                        'MJ/kg', 'LCA'),))

# # Natural gas
# metrics.extend((Metric('GWP - natural gas',
#                        lambda:get_ng_GWP(),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - natural gas',
#                        lambda:get_ng_FEC(),
#                        'MJ/kg', 'LCA'),))

# # Electricity
# metrics.extend((Metric('GWP - electricity, net',
#                        lambda:get_net_electricity_GWP(),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - electricity, net',
#                        lambda:get_net_electricity_FEC(),
#                        'MJ/kg', 'LCA'),))

# # Feedstock growth, harvesting, transportation, and preprocessing
# metrics.extend((Metric('GWP - Feedstock GHTP',
#                        lambda:get_FGHTP_GWP(),
#                        'kg CO2/kg', 'LCA'),))
# metrics.extend((Metric('FEC - Feedstock GHTP',
#                        lambda:get_feedstock_FEC(),
#                        'MJ/kg', 'LCA'),))

# # Direct non-biogenic emissions GWP
# metrics.extend((Metric('GWP - Other direct non-bio emmissions',
#                        lambda:get_direct_emissions_GWP(),
#                        'kg CO2/kg', 'LCA'),))

# metrics.extend((Metric('GWP - Other direct non-bio emmissions',
#                        lambda:get_direct_emissions_GWP(),
#                        'kg CO2/kg', 'LCA'),))

# # Demand LCA contributions
# metrics.extend((Metric('cGWP - System heating demand',
#                        lambda:get_heating_demand_GWP()/get_GWP(),
#                        'frac', 'LCA'),))

# metrics.extend((Metric('cGWP - System cooling demand',
#                        lambda:get_cooling_demand_GWP()/get_GWP(),
#                        'frac', 'LCA'),))

# metrics.extend((Metric('cGWP - System non-cooling electricity demand',
#                        lambda:get_electricity_demand_non_cooling_GWP()/get_GWP(),
#                        'frac', 'LCA'),))

# metrics.extend((Metric('cFEC - System heating demand',
#                        lambda:get_heating_demand_FEC()/get_FEC(),
#                        'frac', 'LCA'),))

# metrics.extend((Metric('cFEC - System cooling demand',
#                        lambda:get_cooling_demand_FEC()/get_FEC(),
#                        'frac', 'LCA'),))

# metrics.extend((Metric('cFEC - System non-cooling electricity demand',
#                        lambda:get_electricity_demand_non_cooling_FEC()/get_FEC(),
#                        'frac', 'LCA'),))

# # Demand TEA contributions
# metrics.extend((Metric('cVOC - System heating demand',
#                        lambda:get_heating_demand_VOC()/get_VOC(),
#                        'frac', 'LCA'),))

# metrics.extend((Metric('cVOC - System cooling demand',
#                        lambda:get_cooling_demand_VOC()/get_VOC(),
#                        'frac', 'LCA'),))

# metrics.extend((Metric('cVOC - System non-cooling electricity demand',
#                        lambda:get_electricity_demand_non_cooling_VOC()/get_VOC(),
#                        'frac', 'LCA'),))

BDO_model.metrics = metrics