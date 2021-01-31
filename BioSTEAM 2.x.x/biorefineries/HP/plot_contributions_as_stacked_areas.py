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

# %% Setup
spec.load_productivity(0.79)
spec.load_yield(0.49)
spec.load_titer(54.8)

titers = np.linspace(10, 180, 20)
yields = 100.*np.linspace(0.1, 0.99, 20)

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
        try:
            HP_sys.reset_cache()
            HP_sys.empty_recycles()
            HP_sys.simulate()
        except:
            try:
                HP_sys.reset_cache()
                HP_sys.empty_recycles()
                HP_sys.converge_method = 'fixedpoint'
                HP_sys.simulate()
            except:
                print('Failed point.')
                   
    HP_sys.converge_method = 'wegstein'
    try:
        MPSPs.append(get_AA_MPSP())
    except:
        MPSPs.append(0)
    for group in process_groups:
        # hu_group_contributions[group.name].append(group.get_heating_duty()/AA.F_mass)
        hu_group_contributions[group.name].append(0.001*get_group_heating_demand(group)/AA.F_mass)
        pu_group_contributions[group.name].append(group.get_electricity_consumption()/AA.F_mass)
        ic_group_contributions[group.name].append(group.get_installed_cost()/AA.F_mass)


hu_group_contributions_df = pd.DataFrame(hu_group_contributions, index = parameters)
pu_group_contributions_df = pd.DataFrame(pu_group_contributions, index = parameters)
ic_group_contributions_df = pd.DataFrame(ic_group_contributions, index = parameters)


# %% Plot contributions
contributions_to_plot = hu_group_contributions_df

ax = contributions_to_plot.plot.area()