#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

# Run this cell first
from warnings import filterwarnings
filterwarnings('ignore')

import contourplots
get_rounded_str = contourplots.utils.get_rounded_str

from biosteam.utils import  colors
import numpy as np

from biorefineries import HP
# from biorefineries.HP.systems.system_sc_light_lle_vacuum_distillation import HP_tea, HP_lca, R302, spec, AA, simulate_and_print, get_AA_MPSP
# from biorefineries.HP.systems.corn.system_corn_improved_separations import HP_tea, HP_lca, R302, spec, AA, simulate_and_print, get_AA_MPSP

from biorefineries.HP.systems.glucose.system_glucose_hexanol import HP_tea, HP_lca, R302, spec, AA, simulate_and_print, get_AA_MPSP

from biorefineries.HP.models.glucose import models_glucose_hexanol as models

from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd

from math import floor, ceil
from datetime import datetime

from math import log

import os


import biosteam as bst


from biorefineries.tea.cellulosic_ethanol_tea import capex_table, foc_table

f = bst.main_flowsheet

chdir = os.chdir

dateTimeObj = datetime.now()

ig = np.seterr(invalid='ignore')

product = AA

AA_market_range=np.array([1.40, 1.65]) 

fossilbased_GWPs = [
                2.2414 + HP_lca.EOL_GWP, # ecoinvent 3.8 (acrylic acid production, RoW) cradle-to-gate + EOL
                5.9589 # GREET 2023 (acrylic acid from fossil energy) cradle-to-grave
                ]
fossilbased_FECs = [
                49.013, # ecoinvent 3.8 (acrylic acid production, RoW)
                116. # GREET 2023 (acrylic acid from fossil energy)
                ]

#%% Filepaths
HP_filepath = HP.__file__.replace('\\__init__.py', '')

# ## Change working directory to biorefineries\\HP\\analyses\\results
# chdir(HP.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
HP_results_filepath = HP_filepath + '\\analyses\\results\\'


#%% Load baseline

spec.reactor.neutralization = False # !!! set neutralization here

model = models.HP_model
system = HP_sys = models.HP_sys

get_AA_MPSP()

feedstock_tag = 'glucose'
product_tag = 'Acrylic'

# mode = '300L_FGI' # !!! remember to switch if running FGI analyses

mode = 'DASbox'

if mode == 'DASbox':
    spec.baseline_productivity = 0.548

dist_filename = f'parameter-distributions_{feedstock_tag}_{product_tag}_' + mode + '.xlsx'

product_folder = 'acrylic_acid_product' if product_tag=='Acrylic' else 'HP_salt_product'

parameter_distributions_filename = HP_filepath+\
    f'\\analyses\\full\\parameter_distributions\\{product_folder}\\'+dist_filename


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

# ## Change working directory to biorefineries\\HP\\analyses\\results
# chdir(HP.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##

model.exception_hook = 'warn'
print('\n\nSimulating baseline ...')
baseline_initial = model.metrics_at_baseline()

spec.load_specifications(spec_1=0.694, spec_2=92, spec_3=0.548)

print(get_AA_MPSP())
get_AA_MPSP()
simulate_and_print()

#%% Save tables
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
file_to_save = f'_baseline_glucose_neutral={R302.neutralization}_'+'HP_TRY_%s.%s.%s-%s.%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute, dateTimeObj.second)

foc_table(HP_tea).to_excel('FOC_breakdown'+file_to_save+'.xlsx')
capex_table(HP_tea).to_excel('TCI_breakdown'+file_to_save+'.xlsx')
