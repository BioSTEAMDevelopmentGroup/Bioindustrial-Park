# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 21:50:47 2021

@author: yrc2
"""

import numpy as np
from warnings import filterwarnings
import pandas as pd 
filterwarnings('ignore')
from biorefineries.HP.system_light_lle_vacuum_distillation import (
    HP_sys, process_groups, spec, get_AA_MPSP, get_GWP, get_FEC, AA, 
    HXN, get_material_cost_breakdown, feedstock
)

carbs = np.linspace(0.1, 0.7, 20)

def mass(x):
    spec.load_feedstock_carbohydrate_content(x)
    return feedstock.imass['Xylan', 'Glucan'].value
    
def comp(x):
    spec.load_feedstock_carbohydrate_content(x)
    return feedstock.imass['Xylan', 'Glucan'].sum() / feedstock.F_mass

# x = np.array([mass(i) for i in carbs])
# y = np.array([mass(i) for i in carbs])
# assert np.abs(x-y).sum() <= 1e-6