#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from warnings import filterwarnings
filterwarnings('ignore')


from biorefineries import TAL

from biorefineries.TAL.systems.system_TAL_solubility_exploit_ethanol_sugarcane import TAL_tea, simulate_and_print

from biorefineries.TAL.models import models_TAL_solubility_exploit as models

import pandas as pd

from datetime import datetime

import os

from biorefineries.tea.cellulosic_ethanol_tea import capex_table, foc_table

chdir = os.chdir

dateTimeObj = datetime.now()


TAL_filepath = TAL.__file__.replace('\\__init__.py', '')

# ## Change working directory to biorefineries\\TAL\\analyses\\results
# chdir(TAL.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
TAL_results_filepath = TAL_filepath + '\\analyses\\results\\'


#%% Load baseline

#%% Filepaths
TAL_filepath = TAL.__file__.replace('\\__init__.py', '')

# ## Change working directory to biorefineries\\TAL\\analyses\\results
# chdir(TAL.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
TAL_results_filepath = TAL_filepath + '\\analyses\\results\\'


#%% Load baseline
model = models.TAL_model
system = TAL_sys = models.TAL_sys

modes = [
            # 'A_FGI_sugarcane',
            'A',
         ]


parameter_distributions_filenames = [
                                    'parameter-distributions_TAL_' + mode +'.xlsx' 
                                    for mode in modes
                                    ]
mode = modes[0]


parameter_distributions_filename = TAL_filepath+\
    '\\analyses\\full\\parameter_distributions\\'+parameter_distributions_filenames[0]
print(f'\n\nLoading parameter distributions ({mode}) ...')
model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename, models.namespace_dict)

# load_additional_params()
print(f'\nLoaded parameter distributions ({mode}).')

parameters = model.get_parameters()

print('\n\nLoading samples ...')
samples = model.sample(N=2000, rule='L')
model.load_samples(samples)
print('\nLoaded samples.')

# ## Change working directory to biorefineries\\TAL\\analyses\\results
# chdir(TAL.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##

model.exception_hook = 'warn'
print('\n\nSimulating baseline ...')
baseline_initial = model.metrics_at_baseline()

#%%
simulate_and_print()

#%% Save tables
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
file_to_save = f'_baseline'

foc_table(TAL_tea).to_excel('FOC_breakdown'+file_to_save+'.xlsx')
capex_table(TAL_tea).to_excel('TCI_breakdown'+file_to_save+'.xlsx')
