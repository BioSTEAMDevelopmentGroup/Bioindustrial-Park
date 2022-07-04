# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 20:21:32 2021

@author: sarangbhagwat
"""

from warnings import filterwarnings
import biosteam as bst
import numpy as np

filterwarnings('ignore')
ig = np.seterr(invalid='ignore')
bst.speed_up()

print ('\n\n1. Baseline:\n')
print ('Loading ...\n')
from biorefineries.HP.system_light_lle_vacuum_distillation import *

# yield
print ('\n\n2. Yield improvement:\n')
print ('Loading ...\n')
spec.load_specifications(0.75, 54.8, 0.76)
simulate_and_print()

# yield + titer
print ('\n\n3. Yield and titer improvement:\n')
print ('Loading ...\n')
spec.load_specifications(0.75, 2*54.8, 0.76)
simulate_and_print()


# yield + titer + ssl
print ('\n\n4. Yield, titer, and saccharification solids loading improvement:\n')
print ('Loading ...\n')
M301.solid_loading = 0.4


full_path = HP_sys.path
evaporator_index = full_path.index(spec.titer_inhibitor_specification.evaporator)
pre_evaporator_units_path = full_path[0:evaporator_index]


spec.pre_conversion_units.simulate()
for unit in pre_evaporator_units_path:
    unit._run()
spec.titer_inhibitor_specification.run_units()
        
spec.load_specifications(0.75, 2*54.8, 0.76)
simulate_and_print()
