# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 23:56:39 2022

@author: sarangbhagwat
"""
from biorefineries.TAL.system_TAL_adsorption_glucose import u, get_SA_MPSP
from matplotlib import pyplot as plt 
import numpy as np

column = u.AC401

MPSPs, aics = [], []
ads_caps = np.linspace(0.0739, 0.2474, 30)

for ac in ads_caps:
    column.adsorbent_capacity = ac
    MPSPs.append(get_SA_MPSP())
    aics.append(column.installed_cost)

plt.plot(ads_caps, MPSPs)
plt.plot(ads_caps, aics)