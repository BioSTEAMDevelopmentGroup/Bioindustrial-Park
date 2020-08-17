#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Thu Jul  9 11:54:08 2020

References:
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of 
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic 
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764; 
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf

[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234

[3] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040
    
@author: yalinli_cabbi
"""


# %%

# =============================================================================
# Simulate pretreatment efficacy
# =============================================================================

import numpy as np
import pandas as pd
from biosteam.utils import TicToc

timer_efficacy = TicToc('timer_efficacy')
timer_efficacy.tic()
 
lignin = np.arange(0, 0.41, 0.01)
conversion_max = np.ones(1000)
conversion_min = np.zeros(1000)

# Liquid hot water (LHW)
np.random.seed(3221)
intercept_LHW_1 = np.random.normal(0.84, 0.04, 1000)
intercept_LHW_2 = np.random.normal(1.32, 0.07, 1000)
slope_LHW_2 = np.random.normal(-2.33, 0.33, 1000)
df_LHW = pd.DataFrame()

# Acid
intercept_acid = np.random.normal(1.04, 0.04, 1000)
slope_acid = np.random.normal(-1.37, 0.18, 1000)
df_acid = pd.DataFrame()

# Explosion (EXP)
intercept_EXP = np.random.normal(0.83, 0.07, 1000)
df_EXP = pd.DataFrame()

# Base
intercept_base = np.random.normal(0.82, 0.09, 1000)
df_base = pd.DataFrame()

# Inoic liquid (IL)
intercept_IL = np.random.normal(1.52, 0.16, 1000)
slope_IL = np.random.normal(-2.87, 0.61, 1000)
df_IL = pd.DataFrame()

# Organic acid (ORG)
intercept_ORG = np.random.normal(0.90, 0.09, 1000)
df_ORG = pd.DataFrame()

# Oxidative (OXD)
intercept_OXD = np.random.normal(0.93, 0.04, 1000)
df_OXD = pd.DataFrame()

for i in range(1, 42):
    LHW_individual_1 = intercept_LHW_1
    LHW_individual_2 = intercept_LHW_2 + lignin[i-1]*slope_LHW_2
    # Select the smaller of two LHW estimations
    LHW_individual = np.minimum(LHW_individual_1, LHW_individual_2)
    acid_individual = intercept_acid + lignin[i-1]*slope_acid
    EXP_individual = intercept_EXP
    base_individual = intercept_base
    IL_individual = intercept_IL + lignin[i-1]*slope_IL
    ORG_individual = intercept_ORG
    OXD_individual = intercept_OXD
    
    # Constrict conversion to [0, 100%]
    df_LHW[round(lignin[i-1],2)] = np.minimum(conversion_max, 
                                              np.maximum(conversion_min, 
                                                          LHW_individual))
    df_acid[round(lignin[i-1],2)] = np.minimum(conversion_max, 
                                                np.maximum(conversion_min, 
                                                          acid_individual))
    df_EXP[round(lignin[i-1],2)] = np.minimum(conversion_max, 
                                              np.maximum(conversion_min, 
                                                          EXP_individual))
    df_base[round(lignin[i-1],2)] = np.minimum(conversion_max,
                                                np.maximum(conversion_min,
                                                          base_individual))
    df_IL[round(lignin[i-1],2)] = np.minimum(conversion_max, 
                                              np.maximum(conversion_min,
                                                        IL_individual))
    df_ORG[round(lignin[i-1],2)] = np.minimum(conversion_max, 
                                              np.maximum(conversion_min,
                                                          ORG_individual))
    df_OXD[round(lignin[i-1],2)] = np.minimum(conversion_max, 
                                              np.maximum(conversion_min,
                                                          OXD_individual))

# Obtain conversion quantiles
dfs =(df_LHW, df_acid, df_EXP, df_base, df_IL, df_ORG, df_OXD)
indices = ['LHW', 'Acid', 'EXP', 'Base', 'IL', 'ORG', 'OXD']
quantiles = [[df.quantile(q=i) for i in (0.05, 0.5, 0.95)] for df in dfs]

df_stats = pd.concat([pd.concat(quantiles[indices.index(index)], axis=1) 
                      for index in indices], 
                      axis=1, keys=indices)
df_stats.rename_axis('Lignin content')

# =============================================================================
# # Save simulated pretreatment efficacies in Excel
# =============================================================================
with pd.ExcelWriter('Pretreatment efficacies.xlsx') as writer:
    df_stats.to_excel(writer, sheet_name='Plotting stats')
    df_LHW.to_excel(writer, sheet_name='LHW')
    df_acid.to_excel(writer, sheet_name='Acid')
    df_EXP.to_excel(writer, sheet_name='EXP')
    df_base.to_excel(writer, sheet_name='Base')
    df_IL.to_excel(writer, sheet_name='IL')
    df_ORG.to_excel(writer, sheet_name='ORG')
    df_OXD.to_excel(writer, sheet_name='OXD')

print(f'\nSimulation time: {timer_efficacy.elapsed_time/60:.1f} min')
print('\n-------- Pretreatment Efficacy Simulation Completed --------\n')


# %%

# =============================================================================
# Biorefinery setup
# =============================================================================

from biosteam.process_tools import UnitGroup
from ethanol_adipic import system_acid as acid
from ethanol_adipic import system_base as base
from ethanol_adipic.chemicals import chems
from ethanol_adipic import utils
# from ethanol_adipic.utils import baseline_feedflow

# Set feedstock flow rate
baseline_feedflow = utils.baseline_feedflow.copy()
# baseline_feedflow = utils.baseline_feedflow.copy() / 2
# baseline_feedflow = utils.baseline_feedflow.copy() * 2
acid.feedstock.mass = baseline_feedflow
base.feedstock.mass = baseline_feedflow

_ethanol_kg_2_gal = acid._ethanol_kg_2_gal
_feedstock_factor = acid._feedstock_factor

simulated_composition = pd.read_excel('Feedstock compositions.xlsx', \
                                      sheet_name='Compositions')
# simulated_composition = simulated_composition[0:5] # for debugging
    
market_ethanol_price = 2.2 / _ethanol_kg_2_gal
# From 80% moisture $/kg to $/dry-U.S. ton to $/kg with 20% moiture
default_feedstock_price = 71.3 / _feedstock_factor

# Function to adjust feedstock flows
index_dict = {'Cellulose': [chems.index('Glucan')],
              'Hemicellulose': chems.indices(('Xylan', 'Arabinan', 'Galactan', 'Mannan')),
              'Lignin': [chems.index('Lignin')]}

indices = sum([index_dict[i] for i in index_dict.keys()], [])
total_CHL_flow = baseline_feedflow[indices].sum()
baseline_hemicellulose_flow = baseline_feedflow[index_dict['Hemicellulose']].sum()
hemicellulose_ratio = {i: baseline_feedflow[i]/baseline_hemicellulose_flow
                       for i in index_dict['Hemicellulose']}

def update_feedstock_flows(feedstock, composition):
    feedstock.imass['Glucan'] = total_CHL_flow * composition['Cellulose']
    feedstock.imass['Lignin'] = total_CHL_flow * composition['Lignin']
    total_hemicellulose_flow = total_CHL_flow * composition['Hemicellulose']
    for i in index_dict['Hemicellulose']:
        feedstock.mass[i] = total_hemicellulose_flow * hemicellulose_ratio[i]

# Function to compute the amount of electricity generated in 10^6 kWh/y,
# positive indicates net production and negative indicates net consumption,
# annual factor is operating hours per year
def compute_electricity(group, annual_factor):
    power_consumption = group.get_electricity_consumption()
    power_production = group.get_electricity_production()
    net_electricity = (power_production-power_consumption)*annual_factor / 1e3
    return net_electricity
        
# Function to compute ethanol yield in gal/dry-ton feedstock
feedstock_dry_ton = baseline_feedflow.sum() / _feedstock_factor
def compute_ethanol_yield(ethanol, feedstock):
    ethanol_gal = ethanol.F_mass / _ethanol_kg_2_gal
    return ethanol_gal/feedstock_dry_ton

# Function to compute minimum ethanol selling price (MESP) in $/gal
def compute_MESP(ethanol, tea):
    return tea.solve_price(ethanol)*_ethanol_kg_2_gal

# Function to compute maximum feedstock payment price (MFPP) in $/dry-ton
def compute_MFPP(feedstock, tea):
    return tea.solve_price(feedstock)*_feedstock_factor


# %%

# =============================================================================
# Run the acid-pretreatment biorefinery for different feedstock compositions
# =============================================================================

timer_acid = TicToc('timer_acid')
timer_acid.tic()

# This is for double-checking, simulated total flow rates should be the same 
# as the default value (default_total_flow, 104180 kg/hr)
acid_total_flow = []

# These lists store output results
acid_conversions_C6 = []
acid_conversions_C5 = []
acid_ethanol_yields = []
acid_produced_electricity = []
acid_MESPs = []
acid_MFPPs = []

# Run assumed compositions (varying cellulose, hemicellulose, and lignin compositions
# while keeping compositions of other components unchanged).
# The first composition in the file is the default one as in refs [1-3]
feedstock_dry_mass = acid.feedstock.F_mass - acid.feedstock.imass['Water']
acid_group = UnitGroup('Acid pretreatment', acid.ethanol_sys.units)
acid_factor = acid.ethanol_no_CHP_tea._annual_factor
for i in range(0, simulated_composition.shape[0]):
    # Update feedstock flow
    update_feedstock_flows(acid.feedstock, simulated_composition.iloc[i])
    
    # Adjust cellulose conversion 
    lignin_percent = acid.feedstock.imass['Lignin'] / feedstock_dry_mass
    # 0.04 and 0.012 are cellulose (glucan) conversion to other products,
    # subtract 1e-6 to avoid getting tiny negatives,
    # 1.04-1.37*lignin_percent is the developed correlation
    C6_conversion = min(1-0.04-0.012-1e-6, max(0, (1.04-1.37*lignin_percent)))
    acid.R301.saccharification_rxns_C6[2].X = C6_conversion

    # Adjust hemicellulose conversion
    # 0.05 and 0.024 are cellulose (glucan) conversion to other products,
    # subtract 1e-6 to avoid getting tiny negatives
    C5_conversion = min(1-0.05-0.024-1e-6, C6_conversion)
    acid.R201.pretreatment_rxns[4].X = C5_conversion # xylan
    acid.R201.pretreatment_rxns[9].X = C5_conversion # mannan
    acid.R201.pretreatment_rxns[12].X = C5_conversion # galactan
    acid.R201.pretreatment_rxns[15].X = C5_conversion # arabinanan
    
    # Simulate system and log results
    acid.ethanol_sys.simulate()
    acid_total_flow.append(acid.feedstock.F_mass)
    acid_conversions_C6.append(C6_conversion)    
    acid_conversions_C5.append(C5_conversion)
    acid_produced_electricity.append(compute_electricity(acid_group, acid_factor))
    acid_ethanol_yields.append(compute_ethanol_yield(acid.ethanol, acid.feedstock))
    
    acid.feedstock.price = default_feedstock_price
    acid_MESPs.append(compute_MESP(acid.ethanol, acid.ethanol_tea))
   
    acid.ethanol.price = market_ethanol_price
    acid_MFPPs.append(compute_MFPP(acid.feedstock, acid.ethanol_tea))
    print(f'Run #{len(acid_MESPs)}: {timer_acid.elapsed_time/60:.1f} min')

print(f'\nSimulation time: {timer_acid.elapsed_time/60:.1f} min')
print('\n-------- Acid-Pretreatment Biorefinery Simulation Completed --------\n\n')


# %%

# =============================================================================
# Run the base-pretreatment biorefinery for different feedstock compositions
# =============================================================================

timer_base = TicToc('timer_base')
timer_base.tic()

base.R502.set_titer_limit = False

# This is for double-checking, simulated total flow rates should be the same 
# as the default value (104180 kg/hr)
base_total_flow = []

# These lists store output results, conversions for C6 and C5 are the same
# (no constrains from other products)
base_conversions = []
base_muconic_titers = []
base_produced_electricity = []
base_ethanol_yields = []
base_MESPs = []
base_MFPPs = []

# Run assumed compositions (varying cellulose, hemicellulose, and lignin compositions
# while keeping compositions of other components unchanged).
# The first composition in the file is the default one as in refs [1-3]
base_group = UnitGroup('Base pretreatment', base.ethanol_adipic_sys.units)
base_factor = base.ethanol_adipic_no_CHP_tea._annual_factor
for i in range(0, simulated_composition.shape[0]):
    update_feedstock_flows(base.feedstock, simulated_composition.iloc[i])
    
    # Adjust cellulose and hemicellulose conversions, 0.82 based on developed correlation
    conversion = 0.82
    # Cellulose
    base.R301.saccharification_rxns_C6.X[2] = conversion
    # Xylan and arabinan, no mention of other carbohydrates conversion in the
    # baseline based on ref [2]
    base.R301.saccharification_rxns_C5.X[:] = conversion
    
    # Simulate system and log results
    base.ethanol_adipic_sys.simulate()
    base_total_flow.append(base.feedstock.F_mass)
    base_conversions.append(conversion)
    base_muconic_titers.append(base.R502.effluent_titer)
    base_produced_electricity.append(compute_electricity(base_group, base_factor))
    base_ethanol_yields.append(compute_ethanol_yield(base.ethanol, base.feedstock))

    base.feedstock.price = default_feedstock_price
    base_MESPs.append(compute_MESP(base.ethanol, base.ethanol_adipic_tea))
    
    base.ethanol.price = market_ethanol_price
    base_MFPPs.append(compute_MFPP(base.feedstock, base.ethanol_adipic_tea))
    print(f'Run #{len(base_MESPs)}: {timer_base.elapsed_time/60:.1f} min')
    
print(f'\nSimulation time: {timer_base.elapsed_time/60:.1f} min')
print('\n-------- Base-Pretreatment Biorefinery Simulation Completed --------\n\n')


# %%

# =============================================================================
# Run the base-pretreatment biorefinery for correlation between muconic acid titer
# (by changing lignin conversion during muconic acid fermentation) and MESP,
# for baseline feedstock composition only.
# =============================================================================

timer_varied = TicToc('timer_varied')
timer_varied.tic()

base.R502.set_titer_limit = False

# These lists store output results
varied_conversion_muconic_titers = []
varied_conversion_produced_electricity = []
# Ethanol yields should not change for different titers
varied_conversion_ethanol_yields = []
varied_conversion_MESPs = []
varied_conversion_MFPPs = []

lignin_conversions = np.arange(0.01, 101, 0.01)
for i in lignin_conversions:
    base.feedstock.mass = baseline_feedflow
    # Adjust muconic acid titer by changing lignin conversion during fermentation
    base.R502.main_fermentation_rxns.X[-1] = i

    # Simulate system and log results
    base.ethanol_adipic_sys.simulate()
    varied_conversion_muconic_titers.append(base.R502.effluent_titer)
    varied_conversion_produced_electricity.append(compute_electricity(base_group, base_factor))
    varied_conversion_ethanol_yields.append(compute_ethanol_yield(base.ethanol, base.feedstock))
    
    base.feedstock.price = default_feedstock_price
    varied_conversion_MESPs.append(compute_MESP(base.ethanol, base.ethanol_adipic_tea))
    
    base.ethanol.price = market_ethanol_price
    varied_conversion_MFPPs.append(compute_MFPP(base.feedstock, base.ethanol_adipic_tea)) 
    print(f'Run #{len(varied_conversion_MESPs)}: {timer_varied.elapsed_time/60:.1f} min')

print(f'\nSimulation time: {timer_varied.elapsed_time/60:.1f} min')
print('\n-------- Varied-conversion Simulation Completed --------\n\n')


# %%

# =============================================================================
# Save biorefinery simulation results in Excel
# =============================================================================

df_varied_composition_results = pd.DataFrame(
    {'Cellulose': simulated_composition['Cellulose'],
     'Hemicellulose': simulated_composition['Hemicellulose'],
     'Lignin': simulated_composition['Lignin'],
     'Acid total flow': acid_total_flow,
     'Acid cellulose conversion': acid_conversions_C6,
     'Acid hemicellulose conversion': acid_conversions_C5,
     'Acid produced electricity [10^6 kWh/y]': acid_produced_electricity,
     'Acid ethanol yield [gal/dry-ton]': acid_ethanol_yields,
     'Acid minimum ethanol selling price [2016$/gal]': acid_MESPs,
     'Acid maximum feedstock payment price [2016$/dry-ton]': acid_MFPPs,
     'Base total flow': base_total_flow,
     'Base carbohydrate conversion': base_conversions,
     'Base muconic acid titer [g/L]': base_muconic_titers,
     'Base produced electricity [10^6 kWh/y]': base_produced_electricity,
     'Base ethanol yield [gal/dry-ton]': base_ethanol_yields,
     'Base minimum ethanol selling price [2016$/gal]': base_MESPs,
     'Base maximum feedstock payment price [2016$/dry-ton]': base_MFPPs,
     'Electricity difference (acid-base)':
         np.asarray(acid_produced_electricity)-np.asarray(base_produced_electricity),
     'Ethanol yield difference (acid-base)':
         np.asarray(acid_ethanol_yields)-np.asarray(base_ethanol_yields),
     'MESP difference (acid-base)': np.asarray(acid_MESPs)-np.asarray(base_MESPs),
     'MFPP difference (acid-base)': np.asarray(acid_MFPPs)-np.asarray(base_MFPPs)
    })

df_varied_conversion_results = pd.DataFrame(
    {'Lignin conversion': lignin_conversions,
     'Muconic acid titer [g/L]': varied_conversion_muconic_titers,
     'Produced electricity [10^6 kWh/y]': varied_conversion_produced_electricity,
     'Ethanol yield [gal/dry-ton]': varied_conversion_ethanol_yields,
     'Minimum ethanol selling price [2016$/gal]': varied_conversion_MESPs,
     'Maximum feedstock payment price [2016$/dry-ton]': varied_conversion_MFPPs
     })

with pd.ExcelWriter('Biorefinery results.xlsx') as writer:
    df_varied_composition_results.to_excel(writer, sheet_name='Varied composition')
    df_varied_conversion_results.to_excel(writer, sheet_name='Varied lignin conversion (base)')




