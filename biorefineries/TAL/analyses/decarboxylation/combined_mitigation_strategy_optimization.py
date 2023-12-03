# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 14:03:16 2023

@author: sarangbhagwat
"""

import scipy
import numpy as np

from biorefineries.TAL.systems.system_TAL_solubility_exploit_ethanol_sugarcane import\
    get_TAL_MPSP, u


#%% 
np_nan = np.nan

S403, U402, M401, F401 = u.S403, u.U402, u.M401,  u.F401

# @np.vectorize

U402.decarboxylation_conversion_basis = 'equilibrium-based'
def demo_func(x):
    try:
        S403.split, M401.mol_acetylacetone_per_mol_TAL, F401.V = x
        return get_TAL_MPSP()
    except:
        return 20 # arbitrary number higher than highest possible MPSP in explored space


# constraints_nonlinear = (
#         lambda x: abs(U402.outs[0].imol['Acetylacetone']/U402.outs[0].imol['TAL'] - 0.2903),
#         )

# NonlinearConstraint = scipy.optimize.NonlinearConstraint
# LinearConstraint = scipy.optimize.LinearConstraint

bounds = [
            (0., 0.99), # S403.split
            (0., 0.3), # M401.mol_acetylacetone_per_mol_TAL
            (0., 0.5), # F401.V
          ]

print('\n\nStarting optimization ...\n\n')
result = scipy.optimize.differential_evolution(
            demo_func, 
            bounds,
            # constraints=NonlinearConstraint(constraints_nonlinear[0], -0.01, 0.01),
            maxiter=20,
            popsize=20,
            init='latinhypercube',
            disp=True,
            polish=False,
            seed=3221,
            )

#%% 
print(result)