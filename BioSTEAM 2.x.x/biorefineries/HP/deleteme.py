# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 17:23:11 2021

@author: saran
"""
import biosteam as bst
import numpy as np
from biorefineries.HP.system_light_lle_vacuum_distillation import spec, HP_sys, get_AA_MPSP, simulate_and_print



f = bst.main_flowsheet
f('SYS2').converge_method = 'wegstein'
# f('SYS2').molar_tolerance
f('SYS2').maxiter = 100

# %% baseline
spec.load_yield(0.1)
spec.load_titer(40.136)

# %% infeasible material flow SYS2
spec.load_yield(0.1306896551724138)
spec.load_titer(76.55172413793103)

# %% infeasible negative temperature D401
spec.load_yield(0.1920689655172414)
spec.load_titer(85.86305365071001)

# %% maximum number of iterations exceeded; root could not be solved
spec.load_yield(0.1306896551724138)
spec.load_titer(178.96799091089662)

# %% HXN maximum number of iterations exceeded; root could not be solved
spec.load_yield(0.52)
spec.load_titer(36.74157402306099)

# %% HXN maximum number of iterations exceeded; root could not be solved
spec.load_yield(0.37)
spec.load_titer(30.00000090515838)

# %% HXN maximum number of iterations exceeded; root could not be solved
spec.load_yield(0.98)
spec.load_titer(33.37078696254725)

# %% simulate
simulate_and_print()

# %% trial specifications barrage
steps = 5
yields = np.linspace(0.1, 0.99, steps) # yield
titers = np.linspace(30, 320, steps) # titer
MPSPs = []

for i in range(len(yields)):
    MPSPs.append([])
    yield_ = yields[i]
    for titer in titers:
        try:
            spec.load_yield(yield_)
            spec.load_titer(titer)
            MPSPs[i].append(get_AA_MPSP())
        except Exception as e:
            print(e)
            MPSPs[i].append(np.nan)