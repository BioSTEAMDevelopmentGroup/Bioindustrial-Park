#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 08:26:45 2019

Fermentation of organic acids using environmental biotechnology (EB) model.

Structure of EB model based on:
    (1) Tchobanoglous, G.; Burton, F. L.; Stensel, H. D.; Inc, M. & E.; Burton, F. Wastewater Engineering: Treatment and Reuse; McGraw-Hill Education, 2003.
    (2) Rittmann, B. E.; McCarty, P. L. Environmental Biotechnology: Principles and Applications; McGraw-Hill, 2001.

    as:
        dX/dt = mu_max * X * (S/(K_S+S)) * (1-(P/P_max)^n)
        dS/dt = -dX/dt * (1/Y_X/S)
        dP/dt = -dX/dt * (Y_P/S/Y_X/S)
    
    where:
        X is biomass
        S is substrate, g for glucose and x for xylose
        P is product, P_max is maximum product titer
        mu is specific cell growth rate, mu = dX/dt/X, mu_max is the maximum specific cell growth rate
        Y_X/S is apparent yield of biomass per substrate consumed, Y_X/S = -dX/dS
        Y_P/S is apparent yield of product per substrate consumed, Y_P/S = -dP/dS
        n is toxic power
        K_S is half-velocity constnat, the value of S when μ/μ_max = 0.5
        
@author: yalinli_cabbi
"""

#%% Things to do in the future

'''
Create a separate pandas dataframe to store kinetic parameters, so do not to make a module for each acid
'''



#%% Import modules

import numpy as np
from scipy.integrate import odeint


#%% Kinetic modeling

# Fitted kinetic parameters, g for glucose and x for xylose
n = 3 # [-]
P_max = 100 # [kg/m3]
Y_XS_g = 0.08 # [kg/kg]
Y_XS_x = 0.11 # [kg/kg]
Y_PS_g = 0.33 # [kg/kg]
Y_PS_x = 0.66 # [kg/kg]
mu_max_g = 0.21 # [kg/kg]
mu_max_x = 0.087 # [kg/kg]
K_S_g = 30 # [kg/m3]
K_S_x = 7.39 # [kg/m3]

def kinetic_model(C, t) -> '(dXdt, dPdt, dS_gdt, dS_xdt)':
    X, P, S_g, S_x = C

    # Compute coefficients
    mu_g = max(0, mu_max_g * S_g/(K_S_g + K_S_x + S_g + S_x))
    mu_x = mu_max_x * S_x/(K_S_g + K_S_x + S_g + S_x)

    # Compute derivatives
    dXdt = X * (1 - (P/P_max)**n) * (mu_g + mu_x)
    dS_gdt = -dXdt * (1/Y_XS_g) if S_g>0 else 0
    dS_xdt = -dXdt * (1/Y_XS_x) if S_x>0 else 0
    dPdt = dXdt * ((Y_PS_g/Y_XS_g) + (Y_PS_x/Y_XS_x))
    return (dXdt, dPdt, dS_gdt, dS_xdt)

def cal_eff(C0=(0.5,0,50,50), tau=100, step_number=1000):
    X0, P0, S_g0, S_x0 = C0

    # Integrate to get final concentration
    t = np.linspace(0, tau, step_number)
    Ct = odeint(kinetic_model, y0=C0, t=t)
    
    # Remove <0 and cache data
    Xt = Ct[:, 0]
    Pt = Ct[:, 1]
    S_gt = Ct[:, 2]
    S_xt = Ct[:, 3]
 
    # Calculate final fermentation efficiency @ t = tau
    Xf = Xt[-1]
    Pf = Pt[-1]
    S_gf = S_gt[-1]
    S_xf = S_xt[-1]
    S_gf = S_gf if S_gf >0 else 0
    S_xf = S_xf if S_xf >0 else 0
    #eff_g = (S_g0 - S_gf)/S_g0 * Y_PS_g/1. # 100% g-to-la under ideal conditions
    #eff_x = (S_x0 - S_xf)/S_x0 * Y_PS_x/1. # 100% x-to-la under ideal conditions
    #return eff_g, eff_x, Xt, Pt, S_gt, S_xt
    Cf = (Xf, Pf, S_gf, S_xf)
    return Cf


'''
To get some trial data:

import numpy as np
from orgacids import EB
import matplotlib.pyplot as plt

a = EB.cal_eff(C0=(0.5,0,50,50), tau=100)

t = np.linspace(0, 100, 1000)

plt.plot(t,sol[2])
plt.plot(t,sol[3])
plt.plot(t,sol[4])
'''











