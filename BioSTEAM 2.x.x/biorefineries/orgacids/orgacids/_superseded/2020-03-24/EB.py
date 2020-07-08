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
        X is fermentation microbes, g for glucose and x for xylose
        S is substrate, g for glucose and x for xylose
        P is product, P_max is maximum product titer
        mu is specific cell growth rate, mu = dX/dt/X, mu_max is the maximum specific cell growth rate
        Y_X/S is apparent yield of biomass per substrate consumed, Y_X/S = -dX/dS
        Y_P/S is apparent yield of product per substrate consumed, Y_P/S = -dP/dS
        n is toxic power
        K_S is half-velocity constnat, the value of S when μ/μ_max = 0.5
        
@author: yalinli_cabbi
"""


#%% Import modules

import numpy as np
from scipy.integrate import odeint


#%% Kinetic modeling

def EB_model(C, t, kinetic_parameters_g, kinetic_parameters_x):
    X_g, X_x, P, S_g, S_x = C
    n_g, P_max_g, Y_XS_g, Y_PS_g, mu_max_g, K_S_g = kinetic_parameters_g
    n_x, P_max_x, Y_XS_x, Y_PS_x, mu_max_x, K_S_x = kinetic_parameters_x

    # Compute coefficients
    toxic_g = 1 - (P/P_max_g)**n_g
    toxic_x = 1 - (P/P_max_x)**n_x
    mu_g = max(0, mu_max_g * S_g/(K_S_g + S_g + S_x))
    mu_x = max(0, mu_max_x * S_x/(K_S_x + S_g + S_x))

    # Compute derivatives
    dX_gdt = X_g * toxic_g * mu_g
    dX_xdt = X_x * toxic_x * mu_x
    dP_gdt = dX_gdt * (Y_PS_g/Y_XS_g)
    dP_xdt = dX_xdt * (Y_PS_x/Y_XS_x)
    dPdt = dP_gdt + dP_xdt
    dS_gdt = -dX_gdt * (1/Y_XS_g) if S_g>0 else 0
    dS_xdt = -dX_xdt * (1/Y_XS_x) if S_x>0 else 0
    return (dX_gdt, dX_xdt, dPdt, dS_gdt, dS_xdt)

def EB_simulation(tau, feed,
                  kinetic_parameters_g, kinetic_parameters_x):

    # Feed concentration should be in kg/m3 (equivalent of g/L)
    t = np.linspace(0, tau, 1000)
    Ct = odeint(EB_model, y0=feed, t=t, 
                args=(kinetic_parameters_g, kinetic_parameters_x))
    
    # Cache data
    X_gt = Ct[:, 0]
    X_xt = Ct[:, 1]    
    Pt = Ct[:, 2]
    S_gt = Ct[:, 3]
    S_xt = Ct[:, 4]
 
    # Calculate final fermentation efficiency @ t = tau
    X_gf = X_gt[-1]
    X_xf = X_xt[-1]
    Pf = Pt[-1]    
    S_gf = S_gt[-1]
    S_gf = S_gf if S_gf >0 else 0
    S_xf = S_xt[-1]
    S_xf = S_xf if S_xf >0 else 0
    
    # feed is (X_g0, X_x0, P0, S_g0, S_x0)
    P_gf = (feed[3]-S_gf)*kinetic_parameters_g[3]
    P_xf = (feed[4]-S_xf)*kinetic_parameters_x[3]
    
    results = (X_gf, X_xf, Pf, S_gf, S_xf, P_gf, P_xf)
    
    # if round(Pf, 3) == round(P_gf, 3) + round(P_xf, 3):
    #     results = (X_gf, X_xf, Pf, S_gf, S_xf, P_gf, P_xf)
    # else: results = 'Error in kinetic model!'
    return results


# %% To get some trial data:

# results = EB_simulation(tau=36, feed=(0.5535745419467963, 0.2627455016212306, 0.048855664517955824, 62.27719190497697, 39.41185258356007),
#                         kinetic_parameters_g=(3, 61.5, 0.08, 0.27, 0.228, 0),
#                         kinetic_parameters_x=(3, 61.5, 0.06, 0.66, 0.115, 45))

# print(results)
