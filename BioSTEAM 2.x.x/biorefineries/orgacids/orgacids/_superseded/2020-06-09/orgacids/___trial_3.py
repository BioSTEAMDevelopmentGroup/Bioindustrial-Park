#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 11:38:47 2020

@author: yalinli_cabbi
"""


import numpy as np
import pandas as pd
from orgacids.models import orgacids_model

# Set seed to make sure each time the same set of random numbers will be used
# Note that change N_simulation to 1 won't get the same first sample as N_simulation = 6
np.random.seed(3221)
N_simulation = 6
samples = orgacids_model.sample(N=N_simulation, rule='L')

# Samples must be 2 dimensional
first_sample = samples[0]
first_sample_2d = first_sample[np.newaxis, :]

# This will save the parameters as an Excel file
parameter_index = [i.index for i in orgacids_model.get_parameters()]
samples_df = pd.DataFrame(data=samples, columns=parameter_index)
first_sample_2d_df = pd.DataFrame(data=first_sample_2d, columns=parameter_index)

with pd.ExcelWriter('Samples.xlsx') as writer:
    samples_df.to_excel(writer, sheet_name='All samples')
    first_sample_2d_df.to_excel(writer, sheet_name='First sample')

# This will trigger the HXN error, note that EtOH's critical point is around 514 K,
# so it does not have Psat beyond that point
orgacids_model.load_samples(first_sample_2d)
orgacids_model.evaluate()

# You can still use 
# from orgacids.system import *
# You can also use the old unit names (e.g., S4ex)
# Also, if you run lactic_acid.show(N=100), you can see there's significant SuccinicAcid
# in the products (I added 0-5% SA in feedstock for uncertainty analysis)
# so we might need to modify the system to add in EthylSuccinate and related rxns


# Note: Yalin found that simulating the system for three times is sufficient
# to reduce the MPSP error to below $0.001/kg with a molar_tolerance of 0.1
# Yalin has modified biosteam.evaluation._model.py (in Dropbox) to allow 
# the system to simulate for a given time (set by simulation_number)
# In the future, can run orgacids_model.evaluate(simulation_number=3)