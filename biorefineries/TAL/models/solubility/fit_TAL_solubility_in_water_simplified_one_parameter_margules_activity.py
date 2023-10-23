#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2023-2024, Sarang Bhagwat <sarangb2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Created on Fri Aug 18 17:32:13 2023

@author: sarangbhagwat
"""
import numpy as np
from math import exp as math_exp
from scipy.optimize import minimize
from biorefineries.TAL.models.solubility.plot_utils import plot_solubility_model
import pandas as pd
from biorefineries import TAL

#%% Load experimental 

TAL_filepath = TAL.__file__.replace('\\__init__.py', '')


experimental_data_df = pd.read_excel(TAL_filepath+'\\models\\solubility\\'+'experimental_data_TAL_solubility_in_water.xlsx')

experimental_Ts = 273.15 + np.array(list(experimental_data_df['Temperature (degrees C)']))
experimental_solubilities = np.array(list(experimental_data_df['TAL solubility in water (g-TAL/L-water)']))

# %% Utils
R = 8.3145
TAL_molar_mass = 126.11004 # g/mol
H2O_molar_volume = 18.01528 # 18 mL/mol
mL_per_L = 1000. # mL/L

TAL_Hm = 30883.66976 # by Dannenfelser-Yalkowsky method as given in Poling, Prausnitz, O'Connel - The properties of gases and liquids
TAL_Tm = 458.15 # K
TAL_Hm_by_R = TAL_Hm/R

def get_TAL_solubility_in_water_for_fit(T, TAL_c_fit): # mol TAL : mol (TAL+water)
    return math_exp(-(TAL_Hm_by_R) * (1/T - 1/TAL_Tm) -\
                    (TAL_c_fit/(R*T)))

def get_mol_TAL_dissolved_for_fit(T, mol_water, TAL_c_fit):
    TAL_x = get_TAL_solubility_in_water_for_fit(T, TAL_c_fit)
    return mol_water*TAL_x/(1-TAL_x)

def get_TAL_solubility_in_water_gpL_for_fit(T, TAL_c_fit):
    return get_mol_TAL_dissolved_for_fit(T, mL_per_L/H2O_molar_volume, TAL_c_fit)*TAL_molar_mass

def get_negative_Rsq(c):
    TAL_c_fit = c[0]
    TSS, RSS = 0, 0
    y_mean = experimental_solubilities.mean()
    for i in range(len(experimental_Ts)):
        T = experimental_Ts[i]
        yi = experimental_solubilities[i]
        y_pred = get_TAL_solubility_in_water_gpL_for_fit(T, TAL_c_fit)
        RSS += (yi - y_pred)**2
        TSS += (yi - y_mean)**2
    return - (1 - RSS/TSS)

#%% Fit
res = minimize(fun=get_negative_Rsq, 
               x0=np.array([1000]), # initial guess
               )
TAL_c = res.x

print(f'\nSolubility model fit to experimental data with R^2 = {round(-get_negative_Rsq(TAL_c), 3)}.\n')


def get_TAL_solubility_in_water(T, TAL_c_fit=TAL_c): # mol TAL / mol (TAL+water)
    return math_exp(-(TAL_Hm_by_R) * (1/T - 1/TAL_Tm) -\
                    (TAL_c_fit/(R*T)))

def get_mol_TAL_dissolved(T, mol_water, TAL_c_fit=TAL_c): # mol TAL dissolved in given mol water
    TAL_x = get_TAL_solubility_in_water(T, TAL_c_fit)
    return mol_water*TAL_x/(1-TAL_x)

def get_TAL_solubility_in_water_gpL(T, TAL_c_fit=TAL_c): # g TAL / L water
    return get_mol_TAL_dissolved(T, mL_per_L/H2O_molar_volume, TAL_c_fit)*TAL_molar_mass

#%% Plot


plot_solubility_model(experimental_Ts=experimental_Ts,
                      experimental_solubilities=experimental_solubilities,
                      get_TAL_solubility_in_water_gpL=get_TAL_solubility_in_water_gpL,
                      R_squared=-get_negative_Rsq(TAL_c),
                      filename='plot_TAL_solubility_in_water_simplified_one_parameter_margules_activity.png',
                      )
