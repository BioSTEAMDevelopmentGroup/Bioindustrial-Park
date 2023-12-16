# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 14:03:16 2023

@author: sarangbhagwat
"""

import scipy
import numpy as np

from biorefineries.TAL.systems.system_TAL_solubility_exploit_ethanol_sugarcane import\
    get_TAL_MPSP, u


#%% Setup
np_nan = np.nan

S403, U402, M401, F401 = u.S403, u.U402, u.M401,  u.F401

# @np.vectorize

U402.decarboxylation_conversion_basis = 'equilibrium-based'

def simulate_function(x):
    try:
        # S403.split, M401.mol_acetylacetone_per_mol_TAL, F401.V = x
        S403.split, M401.mol_acetylacetone_per_mol_TAL = x
        return get_TAL_MPSP()
    except:
        return 20 # arbitrary number higher than highest possible MPSP in explored space
    
bounds = [
            (0., 0.99), # S403.split
            (0., 0.3), # M401.mol_acetylacetone_per_mol_TAL
            # (0., 0.5), # F401.V
          ]

#%% Differential evolution

# constraints_nonlinear = (
#         lambda x: abs(U402.outs[0].imol['Acetylacetone']/U402.outs[0].imol['TAL'] - 0.2903),
#         )

# NonlinearConstraint = scipy.optimize.NonlinearConstraint
# LinearConstraint = scipy.optimize.LinearConstraint


print('\n\nStarting optimization ...\n\n')
result = scipy.optimize.differential_evolution(
            simulate_function, 
            bounds,
            # constraints=NonlinearConstraint(constraints_nonlinear[0], -0.01, 0.01),
            maxiter=30,
            popsize=10,
            init='latinhypercube',
            disp=True,
            polish=False,
            seed=3221,
            )

# S403.split, M401.mol_acetylacetone_per_mol_TAL, F401.V = result.x
S403.split, M401.mol_acetylacetone_per_mol_TAL = result.x
#%% 
print(result)

#%% Enumeration

# steps = 5
# splits, addns, Vs = np.linspace(bounds[0][0], bounds[0][1], steps),\
#                     np.linspace(bounds[1][0], bounds[1][1], steps),\
#                     np.linspace(bounds[2][0], bounds[2][1], steps)

# MPSPs = []
# tot_sims = steps**3
# sim_ct = 0

# min_MPSP, opt_design = 10, []
# for i in splits:
#     MPSPs.append([])
#     for j in addns:
#         MPSPs[-1].append([])
#         for k in Vs:
#             S403.split, M401.mol_acetylacetone_per_mol_TAL, F401.V = i, j, k
#             try:
#                 MPSPs[-1][-1].append(get_TAL_MPSP())
#             except:
#                 MPSPs[-1][-1].append(get_TAL_MPSP())
#             sim_ct+=1
#             print(f'{sim_ct}/{tot_sims}\n')
#             print(f'S403.split = {i}, M401 addn = {j}, F401.V = {k}\n')
#             print(f'MPSP = {MPSPs[-1][-1][-1]}\n')
#             print('--------------------------------------------------')
#             if min_MPSP>MPSPs[-1][-1][-1]:
#                 min_MPSP = MPSPs[-1][-1][-1]
#                 opt_design = [i, j, k]

