# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 12:50:22 2021

@author: sarangbhagwat
"""

import numpy as np
from warnings import filterwarnings
import pandas as pd 
filterwarnings('ignore')
from biorefineries.HP.system_light_lle_vacuum_distillation import HP_sys, process_groups, spec, get_AA_MPSP, AA, HXN
from biorefineries.HP.analyses.models import HP_model as model

from biosteam.plots import plot_montecarlo_across_coordinate

import pandas as pd
from datetime import datetime
_kg_per_ton = 907.18474
system = HP_sys
_red_highlight_white_text = '\033[1;47;41m'
_yellow_text = '\033[1;33m'
_reset_text = '\033[1;0m'

# %% Bugfix barrage

def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    spec.load_yield(0.49)
    spec.load_titer(54.8)
    system.simulate()
    print('Loading and simulating with required specifications ...')
    spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2)
    system.simulate()

def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    system.simulate()

def run_bugfix_barrage():
    try:
        reset_and_reload()
    except Exception as e:
        print(str(e))
        try:
            reset_and_switch_solver('fixedpoint')
        except Exception as e:
            print(str(e))
            try:
                reset_and_switch_solver('aitken')
            except Exception as e:
                print(str(e))
                print("Bugfix barrage failed.")
    finally:
        system.converge_method = 'wegstein'
        print('\n')

# %% Setup

spec.load_spec_1 = spec.load_yield
spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity

spec.load_productivity(0.79)
spec.load_yield(0.49)
spec.load_titer(54.8)

####################
steps = 35
N_simulation = 150
####################

samples = model.sample(N=N_simulation, rule='L')
model.load_samples(samples)

titers = np.linspace(5., 150., steps)
yields = 100.*np.linspace(0.1, 0.99, steps)

ys = []


parameters = titers

baseline_sample = model.get_baseline_sample()

# %% Utils

def model_specification():
    try:
        model._system._converge()
        # for unit in pre_evaporator_units_path:
        #     unit._run()
        spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2)
        model._system.simulate()
    except Exception as e:
        str_e = str(e)
        print('Error in model spec: %s'%str_e)
        if 'sugar concentration' in str_e:
            raise e
        else:
            # try:
            run_bugfix_barrage()
            # except:
            #     print('An iteration in Monte Carlo failed.')
                
model.specification = model_specification

# ###############################
# # Bugfix barrage
# ###############################

# system = HP_sys

# def reset_and_reload():
#     print('Resetting cache and emptying recycles ...')
#     system.reset_cache()
#     system.empty_recycles()
#     print('Loading and simulating with baseline specifications ...')
#     spec_1, spec_2 = spec.spec_1, spec.spec_2
#     spec.load_yield(0.49)
#     spec.load_titer(54.8)
#     system.simulate()
#     print('Loading and simulating with required specifications ...')
#     spec.load_specifications(spec_1=spec_1, spec_2=spec_2)
#     system.simulate()
    
# def reset_and_switch_solver(solver_ID):
#     system.reset_cache()
#     system.empty_recycles()
#     system.converge_method = solver_ID
#     print(f"Trying {solver_ID} ...")
#     system.simulate()
    
# def run_bugfix_barrage():
#     try:
#         reset_and_reload()
#     except Exception as e:
#         print(str(e))
#         try:
#             reset_and_switch_solver('fixedpoint')
#         except Exception as e:
#             print(str(e))
#             try:
#                 reset_and_switch_solver('aitken')
#             except Exception as e:
#                 print(str(e))
#                 print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
#                 raise e
###############################


# %% Run Monte Carlo across one parameter
ys_dict = {}
ys_dict['MPSP'] = []
ys_dict['GWP'] = []
ys_dict['FEC'] = []
# for titer in titers:
for parameter in parameters:
    print(f"{parameter} ...")
    y_vectors = list(np.ones(3))
    try:
        spec.load_titer(parameter)
        # spec.load_yield(parameter/100.)
        # spec.load_titer(54.8)
    except:
        y_vectors = np.nan  *np.ones(3)
    try:
        model.evaluate()
        y_vectors[0] = model.table.Biorefinery['Minimum selling price [$/kg]'].values
        y_vectors[1] = model.table.LCA['Total GWP [kg CO2-eq/kg]'].values
        y_vectors[2] = model.table.LCA['Total FEC [MJ/kg]'].values
    except Exception as e:
        print(str(e))
        print(_red_highlight_white_text+'Point failed.\n'+_reset_text)
        y_vectors = np.nan  *np.ones(3)
    finally:
        setters = model._setters
        for f,s in zip(setters, baseline_sample): f(s)
        HP_sys.converge_method = 'wegstein'
    # ys_dict['MPSP'].append(_kg_per_ton* y_vectors[0])
    # ys_dict['GWP'].append(y_vectors[1])
    # ys_dict['FEC'].append(y_vectors[2])
    
    ys_dict['MPSP'].append(_kg_per_ton* model.table.Biorefinery['Minimum selling price [$/kg]'].values)
    ys_dict['GWP'].append(model.table.LCA['Total GWP [kg CO2-eq/kg]'].values)
    ys_dict['FEC'].append(model.table.LCA['Total FEC [MJ/kg]'].values)

# MPSPs_dict = {'MPSP':MPSPs}
# MPSPs_df = pd.DataFrame(MPSPs_dict, index = parameters)
# xs = []

# for i in range(len(parameters)):
#     xs.append(parameters)
MPSPs_df = pd.DataFrame(ys_dict['MPSP'], index = parameters)
GWPs_df = pd.DataFrame(ys_dict['GWP'], index = parameters)
FECs_df = pd.DataFrame(ys_dict['FEC'], index = parameters)

# %% Plot montecarlo across parameter
# ys = np.array(ys_dict['MPSP'])
# ys = np.array(ys_dict['GWP'])
ys_dict['MPSP'] = np.array(ys_dict['MPSP']).transpose()
ys_dict['GWP'] = np.array(ys_dict['GWP']).transpose()
ys_dict['FEC'] = np.array(ys_dict['FEC']).transpose()
# ys = np.transpose(ys)

# %% MPSP
R, G, B = 254, 221, 80
percentiles_MPSP = plot_montecarlo_across_coordinate(parameters, ys_dict['MPSP'],
                                  light_color = [R/255., G/255., B/255.])
percentiles_MPSP_df = pd.DataFrame(percentiles_MPSP.transpose())
# %% GWP
R, G, B = 0, 169, 150
percentiles_GWP = plot_montecarlo_across_coordinate(parameters, ys_dict['GWP'],
                                  light_color = [R/255., G/255., B/255.])
percentiles_GWP_df = pd.DataFrame(percentiles_GWP.transpose())
# %% FEC
R, G, B = 152, 135, 110
percentiles_FEC = plot_montecarlo_across_coordinate(parameters, ys_dict['FEC'],
                                  light_color = [R/255., G/255., B/255.])
percentiles_FEC_df = pd.DataFrame(percentiles_FEC.transpose())
# %% Save as excel file
dateTimeObj = datetime.now()
file_to_save = 'HP_montecarlo_across_single_parameter_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)

with pd.ExcelWriter(file_to_save+'.xlsx') as writer:
    MPSPs_df.to_excel(writer, sheet_name='MPSP')
    percentiles_MPSP_df.to_excel(writer, sheet_name='MPSP percentiles')
    GWPs_df.to_excel(writer, sheet_name='GWP')
    percentiles_GWP_df.to_excel(writer, sheet_name='GWP percentiles')
    FECs_df.to_excel(writer, sheet_name='FEC')
    percentiles_FEC_df.to_excel(writer, sheet_name='FEC percentiles')
    