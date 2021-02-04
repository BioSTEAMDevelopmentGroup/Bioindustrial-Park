# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 14:33:03 2021

@author: sarangbhagwat
"""

import numpy as np
from warnings import filterwarnings
import pandas as pd 
filterwarnings('ignore')
from biorefineries.HP.system_light_lle_vacuum_distillation import HP_sys, process_groups, spec, get_AA_MPSP, AA, HXN

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

steps = 20

titers = np.linspace(15., 150., steps)
yields = 100.*np.linspace(0.1, 0.99, steps)

hu_group_contributions = {}
pu_group_contributions = {}
ic_group_contributions = {}
MPSPs = []

process_groups[0].name = 'Preprocessing'
process_groups[1].name = 'Pretreatment'
process_groups[2].name = 'Conversion'
process_groups[3].name = 'Separation'
process_groups[4].name = 'Wastewater treatment'
process_groups[5].name = 'Heat exchanger network'
process_groups[6].name = 'Non-HX facilities'


parameters = titers

# %% Utils
def get_group_heating_demand(group):
    return sum([sum([hu.duty for hu in unit.heat_utilities if hu.duty*hu.flow>0.]) for unit in group.units])

# %% Run contributions analysis across one parameter
for group in process_groups:
    hu_group_contributions[group.name] = []
    pu_group_contributions[group.name] = []
    ic_group_contributions[group.name] = []

# for titer in titers:
for parameter in parameters:
    spec.load_titer(parameter)
    # spec.load_yield(parameter/100.)
    
    try:
        HP_sys.simulate()
    except:
        run_bugfix_barrage()
                   
    HP_sys.converge_method = 'wegstein'
    try:
        MPSPs.append(_kg_per_ton*get_AA_MPSP())
    except:
        MPSPs.append(0)
    for group in process_groups:
        # hu_group_contributions[group.name].append(group.get_heating_duty()/AA.F_mass)
        hu_group_contributions[group.name].append(0.001*get_group_heating_demand(group)/AA.F_mass)
        pu_group_contributions[group.name].append(group.get_electricity_consumption()/AA.F_mass)
        ic_group_contributions[group.name].append(group.get_installed_cost()/AA.F_mass)

MPSPs_dict = {'MPSP':MPSPs}
MPSPs_df = pd.DataFrame(MPSPs_dict, index = parameters)

# hu_group_contributions['MPSP'] = MPSPs
# pu_group_contributions_df['MPSP'] = MPSPs
# ic_group_contributions_df['MPSP'] = MPSPs

hu_group_contributions_df = pd.DataFrame(hu_group_contributions, index = parameters)
pu_group_contributions_df = pd.DataFrame(pu_group_contributions, index = parameters)
ic_group_contributions_df = pd.DataFrame(ic_group_contributions, index = parameters)


# %% Plot contributions
contributions_to_plot = hu_group_contributions_df

ax = contributions_to_plot.plot.area()
MPSPs_df.plot.line(ax=ax, secondary_y  = ['MPSP'], color = 'blue')

# %% Save as excel file
dateTimeObj = datetime.now()
file_to_save = 'HP_stacked_area_contributions_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)

with pd.ExcelWriter(file_to_save+'.xlsx') as writer:
    hu_group_contributions_df.to_excel(writer, sheet_name='Heat utility')
    pu_group_contributions_df.to_excel(writer, sheet_name='Power utility')
    ic_group_contributions_df.to_excel(writer, sheet_name='Installed utility')
    MPSPs_df.to_excel(writer, sheet_name='MPSP')
    