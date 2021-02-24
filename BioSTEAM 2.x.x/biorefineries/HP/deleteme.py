# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 17:23:11 2021

@author: saran
"""
import biosteam as bst
import numpy as np
from biorefineries.HP.system_light_lle_vacuum_distillation import spec, HP_sys, get_AA_MPSP, simulate_and_print
from biorefineries.HP.system_light_lle_vacuum_distillation import get_GWP, get_feedstock_FEC, get_feedstock_GWP, get_material_GWP, get_material_FEC, get_non_bio_GWP, get_electricity_FEC, get_electricity_GWP, get_FEC, get_GWP 

from biorefineries.HP.system_light_lle_vacuum_distillation import process_groups
f = bst.main_flowsheet
f('SYS2').converge_method = 'wegstein'
# f('SYS2').molar_tolerance
f('SYS2').maxiter = 100

# %% baseline
# spec.load_yield(0.1)
# spec.load_titer(40.136)

# # %% infeasible material flow SYS2
# spec.load_yield(0.1306896551724138)
# spec.load_titer(76.55172413793103)

# # %% infeasible negative temperature D401
# spec.load_yield(0.1920689655172414)
# spec.load_titer(85.86305365071001)

# # %% maximum number of iterations exceeded; root could not be solved
# spec.load_yield(0.1306896551724138)
# spec.load_titer(178.96799091089662)

# # %% HXN maximum number of iterations exceeded; root could not be solved
# spec.load_yield(0.52)
# spec.load_titer(36.74157402306099)

# # %% HXN maximum number of iterations exceeded; root could not be solved
# spec.load_yield(0.37)
# spec.load_titer(30.00000090515838)

# # %% HXN maximum number of iterations exceeded; root could not be solved
# spec.load_yield(0.98)
# spec.load_titer(33.37078696254725)

# %% simulate
simulate_and_print()

# %% trial specifications barrage
steps = 5
yields = np.linspace(0.1, 0.99, steps) # yield
titers = np.linspace(10, 330, steps) # titer
MPSPs = []

FEC = {}
GWP = {}
shape = (len(yields),len(titers))

FEC['Total'] = np.zeros(shape)
GWP['Total'] = np.zeros(shape)

FEC['Feedstock'] = np.zeros(shape)
GWP['Feedstock'] = np.zeros(shape)

FEC['Electricity'] = np.zeros(shape)
GWP['Electricity'] = np.zeros(shape)

FEC['Material'] = np.zeros(shape)
GWP['Material'] = np.zeros(shape)

GWP['Non_bio'] = np.zeros(shape)

Heat_contribs = {}
Cool_contribs = {}
Elec_contribs = {}
for group in process_groups:
    Heat_contribs[group.name] = np.zeros(shape)
    Cool_contribs[group.name] = np.zeros(shape)
    Elec_contribs[group.name] = np.zeros(shape)

# F301_H = 
for i in range(len(yields)):
    MPSPs.append([])
    yield_ = yields[i]
    for j in range(len(titers)):
        titer = titers[j]
        try:
            spec.load_yield(yield_)
            spec.load_titer(titer)
            MPSPs[i].append(get_AA_MPSP())
            
            FEC['Total'][i,j] = get_FEC()
            GWP['Total'][i,j] = get_GWP()
            
            FEC['Feedstock'][i,j] = get_feedstock_FEC()
            GWP['Feedstock'][i,j] = get_feedstock_GWP()
            
            FEC['Electricity'][i,j] = get_electricity_FEC()
            GWP['Electricity'][i,j] = get_electricity_GWP()
            
            FEC['Material'][i,j] = get_material_FEC()
            GWP['Material'][i,j] = get_material_GWP()
            
            GWP['Non_bio'][i,j] = get_non_bio_GWP()
            
            for group in process_groups:
                Heat_contribs[group.name][i,j] = group.get_heating_duty()
                Cool_contribs[group.name][i,j] = group.get_cooling_duty()
                Elec_contribs[group.name][i,j] = group.get_electricity_consumption()
        except Exception as e:
            print('\n\n\n' + str(e) + '\n\n\n')
            MPSPs[i].append(np.nan)