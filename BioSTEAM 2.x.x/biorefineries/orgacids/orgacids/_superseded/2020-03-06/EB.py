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
Create a separate pandas dataframe to store kinetic parameters, \
    so do not need to make a module for each scenario
'''



#%% Import modules

import numpy as np
from scipy.integrate import odeint


#%% Kinetic modeling

n = 3 # [-]
P_max = 61.5 # [kg/m3]
Y_XS_g = 0.08 # [kg/kg]
Y_XS_x = 0.06 # [kg/kg]
Y_PS_g = 0.62 # [kg/kg]
Y_PS_x = 0.66 # [kg/kg]
mu_max_g = 0.101 # [kg/kg]
mu_max_x = 0.125 # [kg/kg]
K_S_g = 24.35 # [kg/m3]
K_S_x = 45 # [kg/m3]
EtOH_over_LA = 0.24 #[g/g]
AceA_over_LA = 0 #[g/g]

def kinetic_model(C, t) -> '(dXdt, dPdt, dS_gdt, dS_xdt)':
    X, P, S_g, S_x = C

    # Compute coefficients
    toxic = 1 - (P/P_max)**n
    mu_g = max(0, mu_max_g * S_g/(K_S_g + S_g + S_x))
    mu_x = max(0, mu_max_x * S_x/(K_S_x + S_g + S_x))

    # Compute derivatives
    dX_gdt = X * toxic * mu_g
    dX_xdt = X * toxic * mu_x
    dXdt = dX_gdt +dX_xdt
    dP_gdt = dX_gdt * (Y_PS_g/Y_XS_g)
    dP_xdt = dX_xdt * (Y_PS_x/Y_XS_x)
    dPdt = dP_gdt + dP_xdt
    dS_gdt = -dX_gdt * (1/Y_XS_g) if S_g>0 else 0
    dS_xdt = -dX_xdt * (1/Y_XS_x) if S_x>0 else 0
    return (dXdt, dPdt, dS_gdt, dS_xdt)

def EB_simulate(feed, tau):
    # Get initial concentrations
    mass_in = feed.imass['FermentationMicrobe', 'LacticAcid', 'Glucose', 'Xylose']
    vol = feed.F_vol
    #concentration_in = mass_in/vol
    X0, P0, S_g0, S_x0 = mass_in/vol
    
    C0 = [X0, P0, S_g0, S_x0]

    # Integrate to get final concentration
    t = np.linspace(0, tau, 1000)
    Ct = odeint(kinetic_model, y0=C0, t=t)
    
    # Cache data
    Xt = Ct[:, 0]    
    Pt = Ct[:, 1]
    S_gt = Ct[:, 2]
    S_xt = Ct[:, 3]
 
    # Calculate final fermentation efficiency @ t = tau
    Xf = Xt[-1]
    Pf = Pt[-1]
    S_gf = S_gt[-1]
    S_gf = S_gf if S_gf >0 else 0
    S_xf = S_xt[-1]
    S_xf = S_xf if S_xf >0 else 0
    Cf = [Xf, Pf, S_gf, S_xf]
    EtOH = (S_g0 - S_gf)* Y_PS_g * EtOH_over_LA
    AceA = (S_x0 - S_xf)* Y_PS_x * AceA_over_LA
    loss = sum(C0) - sum(Cf) - EtOH - AceA
    Cf.append(EtOH)
    Cf.append(AceA)
    Cf.append(loss)
    mass_out = tuple([i * vol for i in Cf])
    return mass_out


# %% To get some trial data:

'''
import biosteam as bst
from orgacids.EB import EB_simulate
from orgacids.species import species

bst.Stream.species = species
sugar = bst.Stream('sugar',
                   H2O=1000,
                   Z_mobilis=11,
                   LacticAcid = 53,
                   Glucose=182,
                   Xylose=205
                   )
sugar.mass[sugar.index('Z_mobilis')]

loss = 0


for i in range(0, 1):
    simulated_flow = EB_simulate(sugar, 24)
    print(simulated_flow, loss)
    sugar.mass[sugar.index('Z_mobilis')] = simulated_flow[0]
    sugar.mass[sugar.index('LacticAcid')] = simulated_flow[1]
    sugar.mass[sugar.index('Glucose')] = simulated_flow[2]
    sugar.mass[sugar.index('Xylose')] = simulated_flow[3]
    loss += (simulated_flow[2] + simulated_flow[3]) - \
        (simulated_flow[0] + simulated_flow[1])
'''


'''
import numpy as np
t = np.linspace(0, 100, 1000)

import matplotlib.pyplot as plt
plt.plot(t,sol[2])
plt.plot(t,sol[3])
plt.plot(t,sol[4])
'''











