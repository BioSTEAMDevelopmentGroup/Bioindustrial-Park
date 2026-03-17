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
optimize_stage_1_time_and_max_n_glu_spikes_for_MPSP =  isobutanol.models.optimize_stage_1_time_and_max_n_glu_spikes_for_MPSP
optimize_max_n_glu_spikes_for_MPSP = isobutanol.models.optimize_max_n_glu_spikes_for_MPSP
plot_kinetic_results = isobutanol.models.plot_kinetic_results
unit_groups_dict = isobutanol.models.unit_groups_dict
model_specification = model.specification
f = model.system.flowsheet
V406 = f.V406
r = V406.kinetic_reaction_system

scenario = 'A'
IBO_filepath = isobutanol.__file__.replace('\\__init__.py', '')
parameter_distributions_filename = IBO_filepath+\
    '\\analyses\\full\\parameter_distributions\\'+\
    f'parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx'
        
model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
baseline_initial = model.metrics_at_baseline()

model_specification()

# for forced batch mode:
# V406.kinetic_reaction_system._te.max_n_glu_spikes = 0
# V406.kinetic_reaction_system.default_max_n_glu_spikes = 0  

if scenario=='A':
    V406.kinetic_reaction_system._te.max_n_glu_spikes = 1
    V406.kinetic_reaction_system.default_max_n_glu_spikes = 1  
    model_specification(threshold_conc_sugars=167.25, target_conc_sugars=172.5)
elif scenario=='B':
    V406.kinetic_reaction_system._te.max_n_glu_spikes = 16
    V406.kinetic_reaction_system.default_max_n_glu_spikes = 16  
    model_specification(threshold_conc_sugars=316.875, target_conc_sugars=318.75)
else:
    raise ValueError(f'Scenario {scenario} not found.')
    
print(f.ethanol.price * f.ethanol.F_mass/f.ethanol.imass['Ethanol'])

plot_kinetic_results()

#%% Plot conc v time
isobutanol_filepath = isobutanol.__file__.replace('\\__init__.py', '')
isobutanol_results_pub_filepath = isobutanol_filepath + '\\analyses\\results\\publication\\'
subfolder_name='Feed-strat\\'

fig, ax = plot_kinetic_results(xlim=(0,80), ylim=(0,250))

# ax.set_xlim(0, 100)
# ax.set_ylim(0, 250)

plt.savefig(isobutanol_results_pub_filepath+subfolder_name+f'conc_v_time_{fbs_spec.threshold_conc_sugars}_{fbs_spec.target_conc_sugars}_{fbs_spec.conc_sugars_feed_spike}_{scenario}.png', 
            transparent = False,  
            facecolor = 'white',
            bbox_inches='tight',
            dpi=600,
            )
