#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Oxalic acid biorefineries.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>, Wenjun Guo <wenjung2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file is to generate results for varying productivities under the certain titer, yield, and feedstock price.
"""

#%%
# Run this cell first
from warnings import filterwarnings
filterwarnings('ignore')

import contourplots
get_rounded_str = contourplots.utils.get_rounded_str

from biosteam.utils import  colors
import numpy as np

from biorefineries import oxalic
# from biorefineries.oxalic.systems.system_sc_light_lle_vacuum_distillation import oxalic_tea, oxalic_lca, R302, spec, AA, simulate_and_print, get_AA_MPSP
# from biorefineries.oxalic.systems.corn.system_corn_improved_separations import oxalic_tea, oxalic_lca, R302, spec, AA, simulate_and_print, get_AA_MPSP

from biorefineries.oxalic.systems.sugarcane.system_sc_broth import oxalic_tea, oxalic_lca, R302, spec, AA, simulate_and_print, get_AA_MPSP, theoretical_max_g_oxalic_per_g_glucose, oxalic_chemicals

from biorefineries.oxalic.models.sugarcane import models_sc_broth as models

from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd

from math import floor, ceil
from datetime import datetime

from math import log

import os


import biosteam as bst

f = bst.main_flowsheet

chdir = os.chdir

dateTimeObj = datetime.now()

ig = np.seterr(invalid='ignore')

product = AA

AA_market_range=np.array([0.75, 3]) 

fossilbased_GWPs = [
                9.8411 + 2*oxalic_chemicals.CO2.MW/oxalic_chemicals.OxalicAcid.MW, # ecoinvent 3.8 (oxalic acid production, RoW) cradle-to-gate + EOL
                ]
fossilbased_FECs = [
                172.91, # ecoinvent 3.8 (oxalic acid production, RoW)
                ]

#%% Filepaths
oxalic_filepath = oxalic.__file__.replace('\\__init__.py', '')

# ## Change working directory to biorefineries\\oxalic\\analyses\\results
# chdir(oxalic.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
oxalic_results_filepath = oxalic_filepath + '\\analyses\\results\\'


#%% Load baseline

spec.reactor.neutralization = True # !!! set neutralization here

model = models.oxalic_model
system = oxalic_sys = models.oxalic_sys

simulate_and_print()

feedstock_tag = 'sugarcane'
product_tag = 'oxalic-broth'

mode = '100mL'

dist_filename = f'parameter-distributions_{feedstock_tag}_{product_tag}_' + mode + '.xlsx'

product_folder = 'oxalic_broth_product' if product_tag=='oxalic-broth' else None

parameter_distributions_filename = oxalic_filepath+\
    f'\\analyses\\full\\parameter_distributions\\{product_folder}\\'+dist_filename


print(f'\n\nLoading parameter distributions ({mode}) ...')
model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename)

# load_additional_params()
print(f'\nLoaded parameter distributions ({mode}).')

parameters = model.get_parameters()

print('\n\nLoading samples ...')
samples = model.sample(N=2000, rule='L')
model.load_samples(samples)
print('\nLoaded samples.')

# ## Change working directory to biorefineries\\oxalic\\analyses\\results
# chdir(oxalic.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##

model.exception_hook = 'warn'
print('\n\nSimulating baseline ...')
baseline_initial = model.metrics_at_baseline()

print(get_AA_MPSP())
simulate_and_print()

#%% Bugfix barrage

def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)
    spec.set_production_capacity(spec.desired_annual_production)
    # system.simulate()
    print('Loading and simulating with required specifications ...')
    spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
    # spec.set_production_capacity(spec.desired_annual_production)
    system.simulate()
    
def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    # spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
    # spec.set_production_capacity(spec.desired_annual_production)
    system.simulate()

def run_bugfix_barrage():
    try:
        reset_and_reload()
    except Exception as e:
        print(str(e))
        if 'length' in str(e).lower():
            # raise e
            system.reset_cache()
            system.empty_recycles()
            system.simulate()
        else:
            try:
                reset_and_switch_solver('fixedpoint')
            except Exception as e:
                print(str(e))
                try:
                    reset_and_switch_solver('aitken')
                except Exception as e:
                    print(str(e))
                    # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                    print("Bugfix barrage failed.\n")
                    # breakpoint()
                    raise e

#%%  Metrics
broth = R302.outs[1]
# SA_price_range = [6500, 7500]

product_chemical_IDs = ['OxalicAcid',]
get_product_MPSP = lambda: oxalic_tea.solve_price(product) / get_product_purity() # USD / pure-kg
get_product_purity = lambda: sum([product.imass[i] for i in product_chemical_IDs])/product.F_mass
# get_production = lambda: sum([product.imass[i] for i in product_chemical_IDs])

# get_product_recovery = lambda: sum([product.imol[i] for i in product_chemical_IDs])/sum([broth.imol[i] for i in ['oxalic', 'CalciumLactate']])
# get_oxalic_AOC = lambda: oxalic_tea.AOC / 1e6 # million USD / y
# get_oxalic_TCI = lambda: oxalic_tea.TCI / 1e6 # million USD

# HXN = f.HXN1001
# oxalic_metrics = [get_product_MPSP, 
              
#                 lambda: oxalic_lca.GWP,
#                 lambda: oxalic_lca.GWP - oxalic_lca.net_electricity_GWP, 
                
#                 lambda: oxalic_lca.FEC, 
#                 lambda: oxalic_lca.FEC - oxalic_lca.net_electricity_FEC,
                
#                 # lambda: len(HXN.original_heat_utils), 
                
#                # get_oxalic_AOC, get_oxalic_TCI, 
#                get_product_purity]

# %% Generate a range of productivity under TY

yield_titer_pairs = np.array([
    (0.2, 28.23), # our work
    (0.818, 122.72), # highest titer reported; possible yield
    # (0.1, 122.72) # highest titer reported; possible yield
], dtype=float)

# ---- productivity range ----
# for example: 50 points between 0.2*baseline and 10*baseline
n_prod_steps = 50
productivities = np.linspace(
    0.2 * spec.baseline_productivity,
    10.0 * spec.baseline_productivity,
    n_prod_steps
)

# ---------- example usage ----------
results = []
for (yld, tit) in yield_titer_pairs:
    for prod in productivities:
        spec.load_specifications(yld, tit, prod)
        for i in range(3): oxalic_sys.simulate()
        msp = get_product_MPSP()
        print(f"{yld:6.3f}  {tit:8.2f}  {prod:15.3f}  {msp:10.3f}")

# %% Generate a range of productivity under TY

yield_titer_pairs = np.array([
    (0.2, 28.23), # our work
    (0.818, 122.72), # highest titer reported; possible yield
    # (0.1, 122.72) # highest titer reported; possible yield
], dtype=float)

# ---- productivity range ----
# for example: 50 points between 0.2*baseline and 10*baseline
n_prod_steps = 50
productivities = np.linspace(
    0.2 * spec.baseline_productivity,
    10.0 * spec.baseline_productivity,
    n_prod_steps
)

# ---------- example usage ----------
results = []
for (yld, tit) in yield_titer_pairs:
    for prod in productivities:
        spec.load_specifications(yld, tit, prod)
        for i in range(3): oxalic_sys.simulate()
        msp = get_product_MPSP()
        print(f"{yld:6.3f}  {tit:8.2f}  {prod:15.3f}  {msp:10.3f}")
#%% Cost & CI vs transport distance

spec.load_specifications(0.2, 28.23, 0.24)
distance = np.linspace(0, 350)
results = []
for dis in distance:
    cost_transport = 0.109 * dis / 1000
    oxalic_sys.flowsheet.sugarcane.price = 0.03592 + cost_transport
    CI_transport = 0.049 * dis / 1000
    oxalic_sys.LCA.CFs['GWP_100']['Sugarcane'] = 0.044535 + CI_transport
    get_CI = lambda: oxalic_lca.GWP
    get_GWP_before_offset = lambda: oxalic_lca.GWP - oxalic_lca.net_electricity_GWP
    for i in range(3): oxalic_sys.simulate()
    msp = get_product_MPSP()
    ci = get_CI()
    ci_wo_e = get_GWP_before_offset()
    print(f"{dis:6.3f}  {msp:10.3f}  {ci:15.3f}  {ci_wo_e:15.3f}")
        

