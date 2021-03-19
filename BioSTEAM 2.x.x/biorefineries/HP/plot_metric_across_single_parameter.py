# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 14:33:03 2021

@author: sarangbhagwat
"""

import numpy as np
from warnings import filterwarnings
import pandas as pd 
filterwarnings('ignore')
from biorefineries.HP.system_light_lle_vacuum_distillation import (
    HP_sys, 
    process_groups, 
    spec, 
    get_AA_MPSP, 
    get_GWP, 
    get_FEC, 
    AA, 
    HXN, 
    F301,
    get_material_cost_breakdown,
    mock_hxn
)
import pandas as pd
from datetime import datetime
_kg_per_ton = 907.18474
system = HP_sys

# %% Bugfix barrage

def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    spec_1 = spec.spec_1
    spec_2 = spec.spec_2
    spec.load_yield(0.49)
    spec.load_titer(54.8)
    system.simulate()
    print('Loading and simulating with required specifications ...')
    spec.load_specifications(spec_1=spec_1, spec_2=spec_2)
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
spec.load_productivity(0.79)
spec.load_yield(0.49)
spec.load_titer(54.8)

steps = 3

# titers = np.linspace(15., 150., steps)
# yields = 100.*np.linspace(0.1, 0.99, steps)

titers = np.array([54.8, 54.8])

carbs = np.linspace(0.1, 0.7, 20)
hu_group_contributions = {}
cu_group_contributions = {}
pu_group_contributions = {}
ic_group_contributions = {}
mc_group_contributions = {}
gwp_group_contributions = {}
fec_group_contributions = {}

# process_groups[0].name = 'Preprocessing'
# process_groups[1].name = 'Pretreatment'
# process_groups[2].name = 'Conversion'
# process_groups[3].name = 'Separation'
# process_groups[4].name = 'Wastewater treatment'
# process_groups[5].name = 'Heat exchanger network'
# process_groups[6].name = 'Non-HX facilities'


parameters = carbs

# %% Utils
def get_group_heating_demand(group):
    return sum([sum([hu.duty for hu in unit.heat_utilities if hu.duty*hu.flow>0.]) for unit in group.units])

def get_group_cooling_demand(group):
    return sum([sum([hu.duty for hu in unit.heat_utilities if hu.duty*hu.flow<0.]) for unit in group.units])

# %% Run contributions analysis across one parameter

upstream, downstream = HP_sys.split(F301.ins[0])

def get_metrics(carbs):
    MPSPs = []
    GWPs = []
    FECs = []
    for x in carbs:
        spec.load_feedstock_carbohydrate_content(x)
        upstream._converge()
        # spec.load_yield(parameter/100.)
        
        spec.load_titer(54.8)
        HP_sys.simulate()
                       
        HP_sys.converge_method = 'wegstein'
        try:
            MPSPs.append(get_AA_MPSP())
        except:
            MPSPs.append(0)
        try:
            FECs.append(get_FEC())
        except:
            FECs.append(0)
        try:
            GWPs.append(get_GWP())
        except:
            GWPs.append(0)
    return (MPSPs, GWPs, FECs)
    

import matplotlib.pyplot as plt

carbs = sorted([*np.linspace(0.1, 0.7, 10), *[0.3, 0.4, 0.5, 0.7]])
carbs_a = [0.1, 0.4, 0.7]
carbs_b = [0.1, 0.3, 0.4, 0.5, 0.7]
MPSPs, GWPs, FECs = get_metrics(carbs)
MPSPs, GWPs, FECs_a = get_metrics(carbs_a)
MPSPs, GWPs, FECs_b = get_metrics(carbs_b)

HXN._cost = mock_hxn
MPSPs, GWPs, FECs_nhxn = get_metrics(carbs)
MPSPs, GWPs, FECs_a_nhxn = get_metrics(carbs_a)
MPSPs, GWPs, FECs_b_nhxn = get_metrics(carbs_b)

plt.figure()
plt.title('No HXN')
plt.scatter(carbs_a, FECs_a_nhxn, c='b', s=60)
plt.scatter(carbs_b, FECs_b_nhxn, c='r', s=30)
plt.scatter(carbs, FECs_nhxn, c='g', s=15)

plt.figure()
plt.title('With HXN')
plt.scatter(carbs_a, FECs_a, c='b', s=60)
plt.scatter(carbs_b, FECs_b, c='r', s=30)
plt.scatter(carbs, FECs, c='g', s=15)

# MPSPs_dict = {'MPSP':MPSPs}
# MPSPs_df = pd.DataFrame(MPSPs_dict, index = parameters)
# GWPs_dict = {'GWP':GWPs}
# GWPs_df = pd.DataFrame(GWPs_dict, index = parameters)
# FECs_dict = {'FEC':FECs}
# FECs_df = pd.DataFrame(FECs_dict, index = parameters)
# # hu_group_contributions['MPSP'] = MPSPs
# # pu_group_contributions_df['MPSP'] = MPSPs
# # ic_group_contributions_df['MPSP'] = MPSPs


# # %% Plot contributions

# MPSPs_df.plot.line(color = 'red')
# GWPs_df.plot.line(color = 'green')
# FECs_df.plot.line(color = 'blue')

# # %% Save as excel file
# dateTimeObj = datetime.now()
# file_to_save = 'HP_metrics_across_single_parameter_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)

# with pd.ExcelWriter(file_to_save+'.xlsx') as writer:
    
#     MPSPs_df.transpose().to_excel(writer, sheet_name='MPSP')
    