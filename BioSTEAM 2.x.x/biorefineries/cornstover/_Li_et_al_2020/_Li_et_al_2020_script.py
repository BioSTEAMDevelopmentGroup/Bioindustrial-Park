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
    for feedstocks with lignin content ranging from 0-40% with a step of 1%,
    1000 simulation was conducted for each lignin content

Then the simulated pretreatment efficacy for acid pretreatment was used 
    in a biorefinery model to calculate minimum ethanol selling price (MESP) 
    and maximum feedstock payment price (MFPP)
    
    The baseline biorefinery was described in Humbird et al.;
    the biorefinery model (python package biorefineries) and process modeling tool 
    (python package biosteam) were described in Cortes-Peña et al.
    
    MESP in Humbird et al. is $2.15/gal using a feedstock price of $58.5/dry-ton,
    both in 2007 U.S. dollars ($2.49/gal and $67.9/dry-ton in 2016$, respectively),
    results were also converted to 2016$ using GDP chain index in Annual Energy Outlook 
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

Note: This script was develoepd for biosteam v2.4.1, biorefineries v2.4.1, 
      and thermosteam v0.4.1.

References:
    (1) Humbirdet al., Technical Report NREL/TP-5100-47764; DOE: NREL, 2011
    (2) Cortes-Peña et al., ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310
"""


# %% System setup

import numpy as np
import pandas as pd
from biosteam import find
from biorefineries.cornstover import system

simulated_composition = pd.read_excel('_feedstock_composition_for_simulation.xlsx', \
                                      sheet_name='composition')

cornstover_sys = system.cornstover_sys
cornstover_tea = cornstover_sys.TEA
cornstover = find.stream.cornstover
ethanol = find.stream.ethanol
# R201 is the preatreatment reactor
R201 = find.unit.R201
# R301 is the saccharification and co-fermentation reactor
R301 = find.unit.R301

_ethanol_density_kggal = 2.98668849
_feedstock_conversion = 907.185 / 0.8 # from 80% moisture $/kg to $/dry-ton (U.S. ton)
# For conversion of 2007$ to 2016$, based on GDP index in Annual Energy Outlook 
# (https://www.eia.gov/outlooks/aeo/)
_2007_to_2016 = 1.16

default_ethanol_price = 2.15 / _ethanol_density_kggal
default_feedstock_price = 58.5 / _feedstock_conversion
default_total_flow = 2205 * 365 * 0.96 * 907.185 / 8410 / (1-0.2) # kg/hr including water


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

# Obtain conversion percentiles  
LHW_5 = np.percentile(df_LHW, 5, axis=0)
LHW_50 = np.percentile(df_LHW, 50, axis=0)
LHW_95 = np.percentile(df_LHW, 95, axis=0)

acid_5 = np.percentile(df_acid, 5, axis=0)
acid_50 = np.percentile(df_acid, 50, axis=0)
acid_95 = np.percentile(df_acid, 95, axis=0)

EXP_5 = np.percentile(df_EXP, 5, axis=0)
EXP_50 = np.percentile(df_EXP, 50, axis=0)
EXP_95 = np.percentile(df_EXP, 95, axis=0)

base_5 = np.percentile(df_base, 5, axis=0)
base_50 = np.percentile(df_base, 50, axis=0)
base_95 = np.percentile(df_base, 95, axis=0)

IL_5 = np.percentile(df_IL, 5, axis=0)
IL_50 = np.percentile(df_IL, 50, axis=0)
IL_95 = np.percentile(df_IL, 95, axis=0)

ORG_5 = np.percentile(df_ORG, 5, axis=0)
ORG_50 = np.percentile(df_ORG, 50, axis=0)
ORG_95 = np.percentile(df_ORG, 95, axis=0)

OXD_5 = np.percentile(df_OXD, 5, axis=0)
OXD_50 = np.percentile(df_OXD, 50, axis=0)
OXD_95 = np.percentile(df_OXD, 95, axis=0)

df_stats = pd.DataFrame([LHW_5, LHW_50, LHW_95,
                         acid_5, acid_50, acid_95,
                         EXP_5, EXP_50, EXP_95,
                         base_5, base_50, base_95,
                         IL_5, IL_50, IL_95,
                         ORG_5, ORG_50, ORG_95,
                         OXD_5, OXD_50, OXD_95], 
                        index=['LHW 5%', 'LHW 50%', 'LHW 95%',
                               'Acid 5%', 'Acid 50%', 'Acid 95%',
                               'EXP 5%', 'EXP 50%', 'EXP 95%',
                               'Base 5%', 'Base 50%', 'Base 95%',
                               'IL 5%', 'IL 50%', 'IL 95%',
                               'ORG 5%', 'ORG 50%', 'ORG 95%',
                               'OXD 5%', 'OXD 50%', 'OXD 95%'
                               ],
                        columns=lignin)


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
cornstover.set_flow(data=water_flow, units='kg/hr', IDs='Water')
cornstover.set_flow(data=chemical_dw_flow, units='kg/hr', IDs=chemical_IDs)

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

# Calculate MESP using default feedstock price
cornstover.price = default_feedstock_price
baseline_MESP_2007 = cornstover_tea.solve_price(stream=ethanol)*_ethanol_density_kggal
# Calculate MFPP using default ethanol price
ethanol.price = default_ethanol_price
baseline_MFPP_2007 = cornstover_tea.solve_price(stream=cornstover)*_feedstock_conversion

df_Humbird_baseline = pd.DataFrame(
    {'2007$': (baseline_MESP_2007, baseline_MFPP_2007),
     '2016$': (baseline_MESP_2007*_2007_to_2016, baseline_MFPP_2007*_2007_to_2016)},
    index=('MESP [$/gal]', 'MFPP [$/dry-ton]')
    )


# %% Propagate uncertainties of pretreatment efficacy to MESP

# Pretreatment efficacies for the lignin content (0.1576) of corn stover as in Humbird et al.
Humbird_conversion =  intercept_acid + 0.1576*slope_acid

# Calculate MESP for each of these 1000 pretreatment efficacies using default feedstock price
cornstover.price = default_feedstock_price
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
    MESP = cornstover_tea.solve_price(stream=ethanol)*_ethanol_density_kggal
    Humbird_MESP_2007.append(MESP.copy())

Humbird_MESP_2016 = [i * _2007_to_2016 for i in Humbird_MESP_2007]

df_Humbird_uncertainties = pd.DataFrame(
    {'Conversion': Humbird_conversion,
     'MESP [2007$/gal]': Humbird_MESP_2007,
     'MESP [2016$/gal]': Humbird_MESP_2016}
    )


# %% Run the baseline biorefinery for different feedstock compositions

simulated_MESP_2007 = []
simulated_MFPP_2007 = []
simulated_conversions = []
# This is for double-checking, simulated total flow rates should be the same 
# as the default value (default_total_flow, 104179.5720749108 kg/hr)
simulated_total_flow = []

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
    
    # Get MESP based on default feedstock price
    cornstover.price = default_feedstock_price
    MESP = cornstover_tea.solve_price(stream=ethanol) * _ethanol_density_kggal   
    simulated_MESP_2007.append(MESP.copy())
   
    # Gest MFPP based on market price of ethanol 
    # (average of ethanol price in 2010-2019 based on Annual Energy Outlook)
    ethanol.price = 1.90 / _ethanol_density_kggal
    MFPP = cornstover_tea.solve_price(stream=cornstover) * _feedstock_conversion
    simulated_MFPP_2007.append(MFPP.copy())

simulated_MESP_2016 = [i * _2007_to_2016 for i in simulated_MESP_2007]
simulated_MFPP_2016 = [i * _2007_to_2016 for i in simulated_MFPP_2007]


df_simulated_results = pd.DataFrame(
    {'Cellulose': simulated_composition['Cellulose'],
      'Hemicellulose': simulated_composition['Hemicellulose'],
      'Lignin': simulated_composition['Lignin'],
      'Conversion': simulated_conversions,
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
    df_Humbird_baseline.to_excel(writer, sheet_name='Humbird TEA baseline')
    df_Humbird_uncertainties.to_excel(writer, sheet_name='Humbird TEA uncertainties')
    df_simulated_results.to_excel(writer, sheet_name='Simulated TEA')

