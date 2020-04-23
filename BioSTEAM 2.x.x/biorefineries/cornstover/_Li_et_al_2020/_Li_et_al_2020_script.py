#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 08:11:18 2019

@author: yalinli_cabbi

This script was used to generate results for a manuscript in preparation (Li et al., 2020),
running this script will save an Excel file with all results in the same directory path 
as this script

Correlations between feedstock lignin content and pretreatment efficacy 
    (as sugar released from cellulose and hemicellulose) were developed 
    for seven different pretreatment technologies as described in the manuscript
    
Monte Carlo simulation was first conducted to get pretreatment efficacy 
    for feedstocks with lignin content ranging from 0-40% with a step size of 1%,
    1000 simulation was conducted for each lignin content

Then the simulated pretreatment efficacy for acid pretreatment was used 
    in a biorefinery model to calculate minimum ethanol selling price (MESP) 
    and maximum feedstock payment price (MFPP)
    
    The baseline biorefinery was described in Humbird et al.;
    the biorefinery model (python package biorefineries) and process modeling tool 
    (python package biosteam) were described in Cortes-Peña et al.
    
    MESP in Humbird et al. is $2.15/gal using a feedstock price of $58.5/dry-ton,
    both in 2007 U.S. dollars ($2.49/gal and $67.9/dry-ton in 2016$, respectively,
    ton is U.S. ton), results were also converted to 2016$ using 
    gross domestic product (GDP) chain-type index in Annual Energy Outlook 
    (https://www.eia.gov/outlooks/aeo/)
    
    Uncertainties of MESP were also calculated using the 1000 different 
    pretreatment efficacies simulated in the previous step
    
Finally, the biorefinery model was run for different feedstock compositions 
    (feedstock composition imported from an Excel file) to simulate MESP 
    (with given feedstock price) and MFPP (with given ethanol price)
    for these different compositions

Requirements:
    (1) Python biorefineries (https://pypi.org/project/biorefineries/)
    (2) Excel file named "_feedstock_composition_for_simulation.xlsx" in the same
        directory path as this script

Note: This script is compatible with biosteam v2.12.10, biorefineries v2.9.1, 
      and thermosteam v0.12.17. Results used in the manuscript were generated using
      biosteam v2.4.1, biorefineries v2.4.1, and thermosteam v0.4.1.
      
      Use pip install package_name==version to install the specific version
      
      If a "FloatingPointError: divide by zero encountered in double_scalars"
      was prompted relating to cool_water, then in the cornstover.unit module,
      replace the following line:
          hot_water.link_with(cool_water)
      with:
          hot_water.link_with(cool_water, TP=False)

References:
    (1) Humbirdet al., Technical Report NREL/TP-5100-47764; DOE: NREL, 2011
    (2) Cortes-Peña et al., ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310
"""


# %% System setup

import numpy as np
import pandas as pd
# R201: pretreatment unit
# R301: saccharification and co-fermentation unit
from biorefineries.cornstover.system import cornstover_sys, cornstover, ethanol, \
    R201, R301

simulated_composition = pd.read_excel('_feedstock_composition_for_simulation.xlsx', \
                                      sheet_name='composition')

cornstover_tea = cornstover_sys.TEA

_ethanol_density_kggal = 2.98668849
_feedstock_conversion = 907.185 / 0.8 # from 80% moisture $/kg to $/dry-ton (U.S. ton)
# For conversion of 2007$ to 2016$, based on GDP index in Annual Energy Outlook 
# (https://www.eia.gov/outlooks/aeo/)
_2007_to_2016 = 1.16

default_ethanol_price = 2.15 / _ethanol_density_kggal
default_feedstock_price = 58.5 / _feedstock_conversion
default_total_flow = 2205 * 365 * 0.96 * 907.185 / 8410 / (1-0.2) # kg/hr including water

# Function to calculate ethanol yield in gal/dry-ton feedstock
calculate_yield = lambda: ethanol.F_mass/_ethanol_density_kggal/ \
    (cornstover.F_mass/_feedstock_conversion)
# Function to calculate minimum ethanol selling price (MESP) in $/gal
calculate_MESP = lambda: cornstover_tea.solve_price(stream=ethanol) * \
    _ethanol_density_kggal
# Function to calculate maximum feedstock payment price (MFPP) in $/dry-ton
calculate_MFPP = lambda: cornstover_tea.solve_price(stream=cornstover)* \
    _feedstock_conversion


# %% Calculate cellulose and hemicellulose conversions 
# based on correlations developed in the manuscript

lignin = np.arange(0, 0.41, 0.01)
conversion_max = np.ones(1000)
conversion_min = np.zeros(1000)

# Liquid hot water (LHW)
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


# %% Calculate minimum ethanol selling price (MESP) and maximum feedstock payment price (MFPP)
# using composition of corn stover in the baseline biorefinery described in Humbird et al.

# Set up baseline feedstock composition as in Table 4 on Page 14 of Humbird et al.
water_flow = 0.2 * default_total_flow # 20% is the default water content
dw_flow = default_total_flow - water_flow
dw_frac = [0.0493, 0.3505, 0.1953, 0.1576, 0.0181, 0.0310, 
           0.1465, 0.0238, 0.0143, 0.0060, 0.0077]
chemical_IDs = ('Ash', 'Glucan', 'Xylan', 'Lignin', 'Acetate', 'Protein',
                'Extract', 'Arabinan', 'Galactan', 'Mannan', 'Sucrose')

cornstover.empty()
chemical_dw_flow = [i*dw_flow for i in dw_frac]
cornstover.set_flow(water_flow, 'kg/hr', 'Water')
cornstover.set_flow(chemical_dw_flow, 'kg/hr', chemical_IDs)

# Cellulose conversion
adjusted_conversion = min(1, max(0, (1.04 - 1.37*0.1576)))
R301.saccharification[2].X = adjusted_conversion
# Adjust hemicellulose conversion
R201.reactions[3].X = adjusted_conversion # galactan
R201.reactions[5].X = adjusted_conversion # mannan
R201.reactions[8].X = adjusted_conversion # xylan
R201.reactions[11].X = adjusted_conversion # arabinanan
# Simulate the system
cornstover_sys.simulate()

# Calculate ethanol yield
baseline_yield = calculate_yield()
# Calculate MESP using default feedstock price
cornstover.price = default_feedstock_price
baseline_MESP_2007 = calculate_MESP()
# Calculate MFPP using default ethanol price
ethanol.price = default_ethanol_price
baseline_MFPP_2007 = calculate_MFPP()

df_Humbird_baseline = pd.DataFrame(
    data=(baseline_yield, baseline_MESP_2007, baseline_MESP_2007*_2007_to_2016,
          baseline_MFPP_2007, baseline_MFPP_2007*_2007_to_2016),
    index=('Yield [gal/dry-ton]', 'MESP [2007$/gal]', 'MESP [2016$/gal]',
           'MFPP [2007$/dry-ton]', 'MFPP [2016$/dry-ton]'),
    columns=('Results',))


# %% Propagate uncertainties of pretreatment efficacy to ethanol yield and MESP

# Pretreatment efficacies for the lignin content (0.1576) of corn stover as in Humbird et al.
Humbird_conversion =  intercept_acid + 0.1576*slope_acid

# Calculate MESP for each of these 1000 pretreatment efficacies using default feedstock price
cornstover.price = default_feedstock_price
Humbird_yield = []
Humbird_MESP_2007 = []
for conversion in Humbird_conversion:
    # Adjust cellulose conversion
    R301.saccharification[2].X = conversion
    
    # Adjust hemicellulose conversion
    R201.reactions[3].X = conversion # galactan
    R201.reactions[5].X = conversion # mannan
    R201.reactions[8].X = conversion # xylan
    R201.reactions[11].X = conversion # arabinanan
    
    # Simulate system and get MESP
    cornstover_sys.simulate()
    Humbird_yield.append(calculate_yield())
    Humbird_MESP_2007.append(calculate_MESP())

Humbird_MESP_2016 = [i * _2007_to_2016 for i in Humbird_MESP_2007]

df_Humbird_uncertainties = pd.DataFrame(
    {'Conversion': Humbird_conversion,
     'Yield [gal/dry-ton]': Humbird_yield,
     'MESP [2007$/gal]': Humbird_MESP_2007,
     'MESP [2016$/gal]': Humbird_MESP_2016}
    )


# %% Run the baseline biorefinery for different feedstock compositions

# This is for double-checking, simulated total flow rates should be the same 
# as the default value (default_total_flow, 104179.5720749108 kg/hr)
simulated_total_flow = []

# These store output results
simulated_conversions = []
simulated_yields = []
simulated_MESP_2007 = []
simulated_MFPP_2007 = []

#               glucan, xylan, arabinan, galactan, mannan, lignin
default_fracs = [0.3505, 0.1953, 0.0238, 0.0143, 0.0060, 0.1576]
default_chl_fracs = [i/sum(default_fracs) for i in default_fracs]
default_xagm_fracs = [i/sum(default_fracs[1:5]) for i in default_fracs[1:5]]

default_cellulose_flow = cornstover.imass['Glucan']
default_hemicellulose_flow = cornstover.imass['Xylan'] + cornstover.imass['Arabinan'] \
                             + cornstover.imass['Galactan'] + cornstover.imass['Mannan']
default_lignin_flow = cornstover.imass['Lignin']
default_chl_flow = default_cellulose_flow + default_hemicellulose_flow + default_lignin_flow

# Run assumed compositions, the assumed compositions are corrected to
# 100% sum of cellulose, hemicellulose, and lignin,
# compositions of other components were assumed to be the same as in Humbird et al.
# The first composition in the file is the default one as in Humbird et al.
for i in range(0, len(simulated_composition['Cellulose'])):
    # Adjust feedstock flows
    adjusted_cellulose_flow =  default_chl_flow * simulated_composition['Cellulose'][i]
    adjusted_hemicellulose_flow =  default_chl_flow * simulated_composition['Hemicellulose'][i]
    adjusted_lignin_flow =  default_chl_flow * simulated_composition['Lignin'][i]
    cornstover.imass['Glucan'] = adjusted_cellulose_flow
    cornstover.imass['Xylan'] = adjusted_hemicellulose_flow * default_xagm_fracs[0]
    cornstover.imass['Arabinan'] = adjusted_hemicellulose_flow * default_xagm_fracs[1]
    cornstover.imass['Galactan'] = adjusted_hemicellulose_flow * default_xagm_fracs[2]
    cornstover.imass['Mannan'] = adjusted_hemicellulose_flow * default_xagm_fracs[3]    
    cornstover.imass['Lignin'] = adjusted_lignin_flow

    simulated_total_flow.append(cornstover.F_mass)
    
    # Calculate updated conversion
    lignin = simulated_composition['Lignin'][i] * default_fracs[-1] / default_chl_fracs[-1]
    conversion = min(1, max(0, (1.04-1.37*lignin)))
    simulated_conversions.append(conversion)
    
    # Adjust cellulose conversion
    R301.saccharification[2].X = conversion
    
    # Adjust hemicellulose conversion
    R201.reactions[3].X = conversion # galactan
    R201.reactions[5].X = conversion # mannan
    R201.reactions[8].X = conversion # xylan
    R201.reactions[11].X = conversion # arabinanan
    
    # Simulate system
    cornstover_sys.simulate()
    
    # Calculate ethanol yield
    simulated_yields.append(calculate_yield())
    
    # Calculate MESP based on default feedstock price
    cornstover.price = default_feedstock_price
    simulated_MESP_2007.append(calculate_MESP())  
   
    # Calculate MFPP based on market price of ethanol
    # (1. is the average of ethanol price of 2010-2019 in 2007$ from Annual Energy Outlook)
    ethanol.price = 1.88 / _ethanol_density_kggal
    simulated_MFPP_2007.append(calculate_MFPP())

simulated_MESP_2016 = [i * _2007_to_2016 for i in simulated_MESP_2007]
simulated_MFPP_2016 = [i * _2007_to_2016 for i in simulated_MFPP_2007]


df_simulated_results = pd.DataFrame(
    {'Cellulose': simulated_composition['Cellulose'],
      'Hemicellulose': simulated_composition['Hemicellulose'],
      'Lignin': simulated_composition['Lignin'],
      'Conversion': simulated_conversions,
      'Ethanol yield [gal/dry-ton]': simulated_yields,
      'MESP [2007$/gal]': simulated_MESP_2007,
      'MESP [2016$/gal]': simulated_MESP_2016,
      'MFPP [2007$/dry-ton]': simulated_MFPP_2007,
      'MFPP [2016$/dry-ton]': simulated_MFPP_2016}
    )


# %% Output all results

with pd.ExcelWriter('_Li_et_al_2020_results.xlsx') as writer:
    df_stats.to_excel(writer, sheet_name='Stats')
    df_LHW.to_excel(writer, sheet_name='LHW')
    df_acid.to_excel(writer, sheet_name='Acid')
    df_EXP.to_excel(writer, sheet_name='EXP')
    df_base.to_excel(writer, sheet_name='Base')
    df_IL.to_excel(writer, sheet_name='IL')
    df_ORG.to_excel(writer, sheet_name='ORG')
    df_OXD.to_excel(writer, sheet_name='OXD')
    df_Humbird_baseline.to_excel(writer, sheet_name='Humbird baseline')
    df_Humbird_uncertainties.to_excel(writer, sheet_name='Humbird uncertainties')
    df_simulated_results.to_excel(writer, sheet_name='Simulated feedstocks')

