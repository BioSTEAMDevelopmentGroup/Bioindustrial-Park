# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 13:48:18 2025

@author: IGB
"""


import numpy as np
import pandas as pd
from chaospy import distributions as shape
import biosteam as bst
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools.parameter import Setter
from biorefineries.nitric._system import sys_plasma
from biorefineries.nitric._tea import create_plasma_tea
from biorefineries.nitric._process_settings import load_preferences_and_process_settings
from warnings import warn
from warnings import filterwarnings; filterwarnings('ignore')

#%%

load_preferences_and_process_settings(T='K',
                                      flow_units='kg/hr',
                                      N=100,
                                      P_units='Pa',
                                      CE=798, # Average 2023 https://toweringskills.com/financial-analysis/cost-indices/
                                      indicator='GWP100',
                                      electricity_EI=0.48,
                                      electricity_price=0.003,)
                                      
sys_plasma.set_tolerance(rmol=1e-6, mol=1e-5, maxiter=400)

tea_plasma = create_plasma_tea(sys=sys_plasma)

sys_plasma.operating_hours = tea_plasma.operating_days * 24

s = sys_plasma.flowsheet.stream

def set_price_of_streams():
    s.water_in.price = 0.002 # DI water or not
                

def set_GWP_of_streams(indicator):
    s.water_in.set_CF(key='GWP100', value=0.001)
    

set_prices = set_price_of_streams()
set_GWP = set_GWP_of_streams(indicator='GWP100')