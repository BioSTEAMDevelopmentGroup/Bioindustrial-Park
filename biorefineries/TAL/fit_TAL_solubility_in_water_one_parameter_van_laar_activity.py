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
from flexsolve import IQ_interpolation
from matplotlib import pyplot as plt

experimental_Ts = 273.15 + np.array([
                                    0,
                                    22,
                                    42,
                                    
                                    53,
                                    54,
                                    56,
                                    
                                    75,
                                    76,
                                    77,
                                    
                                    88,
                                    90,
                                    93,

                                    ])

experimental_Xs = np.array([
                            3.52,
                            8.92,
                            16.46,
                            
                            25.03280709,
                            25.80465597,
                            26.17148061,
                            
                            57.15360402,
                            59.00479116,
                            64.24225554,
                            
                            113.5747233,
                            112.28405,
                            130.6517676,

                            ])

# %% Utils
R = 8.314
TAL_molar_mass = 126.11004 # g/mol
H2O_molar_volume = 18.01528 # mL/mol
mL_per_L = 1000. # mL/L

TAL_Hm = 30883.66976 # by Dannenfelser-Yalkowsky method
TAL_Tm = 458.15 # K
TAL_Hm_by_R = TAL_Hm/R

# V_TAL = 93.55e-6 # m3/mol # from CAS Scifinder
V_TAL = (1*33.5 + 2*13.5 + 1*18.0 + 1*10.0) * 1e-6 # m3/mol # using the method proposed in A Method for Estimating Both the Solubility Parameters and Molar Volumes of liquids, ROBERT F. FEDORS
# 1 CH3 + 2 CH=C + 1 CO2 (i.e., COO) + 1 OH

V_H2O = H2O_molar_volume * 1e-6 # m3/mol

def TAL_solubility_in_water_obj_fn_for_fit(x, T, TAL_c_fit):
    return -x + math_exp(-(TAL_Hm_by_R) * (1/T - 1/TAL_Tm) -\
                         (TAL_c_fit/(R*T))*(1 + (V_TAL*x)/(V_H2O*(1-x)))**-2)

def get_TAL_solubility_in_water_for_fit(T, TAL_c_fit): # mol TAL : mol (TAL+water)
    obj_fn = lambda x: TAL_solubility_in_water_obj_fn_for_fit(x, T=T, TAL_c_fit=TAL_c_fit)
    return IQ_interpolation(obj_fn, 1e-6, 1-1e-6, ytol=1e-6)

def get_mol_TAL_dissolved_for_fit(T, mol_water, TAL_c_fit):
    TAL_x = get_TAL_solubility_in_water_for_fit(T, TAL_c_fit)
    return mol_water*TAL_x/(1-TAL_x)

def get_TAL_solubility_in_water_gpL_for_fit(T, TAL_c_fit):
    return get_mol_TAL_dissolved_for_fit(T, mL_per_L/H2O_molar_volume, TAL_c_fit)*TAL_molar_mass

def get_negative_Rsq(c):
    TAL_c_fit = c[0]
    TSS, RSS = 0, 0
    y_mean = experimental_Xs.mean()
    for i in range(len(experimental_Ts)):
        T = experimental_Ts[i]
        yi = experimental_Xs[i]
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


def TAL_solubility_in_water_obj_fn(x, T, TAL_c_fit=TAL_c):
    return -x + math_exp(-(TAL_Hm_by_R) * (1/T - 1/TAL_Tm) -\
                         (TAL_c_fit/(R*T))*(1 + (V_TAL*x)/(V_H2O*(1-x)))**-2)

def get_TAL_solubility_in_water(T, TAL_c_fit=TAL_c): # mol TAL : mol (TAL+water)
    obj_fn = lambda x: TAL_solubility_in_water_obj_fn(x, T=T, TAL_c_fit=TAL_c_fit)
    return IQ_interpolation(obj_fn, 1e-6, 1-1e-6, ytol=1e-6)

def get_mol_TAL_dissolved(T, mol_water, TAL_c_fit=TAL_c): # mol TAL dissolved in given mol water
    TAL_x = get_TAL_solubility_in_water(T, TAL_c_fit)
    return mol_water*TAL_x/(1-TAL_x)

def get_TAL_solubility_in_water_gpL(T, TAL_c_fit=TAL_c): # g TAL / L water
    return get_mol_TAL_dissolved(T, mL_per_L/H2O_molar_volume, TAL_c_fit)*TAL_molar_mass

#%% Plot

plot_Ts = 273.15 + np.linspace(0., 99., 50)

plot_solubilities = [get_TAL_solubility_in_water_gpL(T) for T in plot_Ts]

plt.plot(plot_Ts-273.15, plot_solubilities, color='blue')

plt.scatter(experimental_Ts-273.15, experimental_Xs, marker='x', color='orange')

plt.text(5, 180, f'R^2 = {round(-get_negative_Rsq(TAL_c), 3)}')

plt.xlabel('Temperature (degrees Celsius)')
plt.ylabel('TAL solubility in water (g-TAL/L-water)')

plt.xlim(0, 100)
plt.ylim(0, 200)

plt.xticks([0 + 10*i for i in range(0,11)])
plt.show()