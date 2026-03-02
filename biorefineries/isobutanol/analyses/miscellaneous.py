# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 18:10:51 2026

@author: saran
"""

import numpy as np
import nskinetics as nsk
import biosteam as bst
from biorefineries import isobutanol
from matplotlib import pyplot as plt

model = isobutanol.models.models_EtOH_IBO_corn.model
namespace_dict = isobutanol.models.namespace_dict
fbs_spec = isobutanol.models.fbs_spec
optimize_1D_feeding_strategy_for_MPSP = isobutanol.models.optimize_1D_feeding_strategy_for_MPSP
optimize_max_n_glu_spikes_for_MPSP =  isobutanol.models.optimize_max_n_glu_spikes_for_MPSP
plot_kinetic_results = isobutanol.models.plot_kinetic_results
unit_groups_dict = isobutanol.models.unit_groups_dict
model_specification = model.specification
f = model.system.flowsheet
V406 = f.V406

IBO_filepath = isobutanol.__file__.replace('\\__init__.py', '')
parameter_distributions_filename = IBO_filepath+\
    '\\analyses\\full\\parameter_distributions\\'+\
    'parameter-distributions_corn_IBO_EtOH_A.xlsx'
        
model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
baseline_initial = model.metrics_at_baseline()
