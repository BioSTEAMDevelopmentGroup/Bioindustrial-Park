"""
Created on Sat Jan 16 19:09:17 2021

@author: sarangbhagwat
"""

#%% Import stream to analyze
from biorefineries.TAL.system_TAL_adsorption_glucose import R302
from biorefineries.make_a_biorefinery.analyses.solvents_barrage import run_solvents_barrage

#%% Run
run_solvents_barrage(stream=R302.outs[0],
                     solute_ID='TAL',
                     impurity_IDs=['VitaminA', 'VitaminD2'],
                     T=30+273.15,
                     stream_modifiers='baseline_stream',
                     criterion='Partition of solute into extract',
                     print_result_with_optimal_criterion=False,
                     show_all_mixer_settlers=False)