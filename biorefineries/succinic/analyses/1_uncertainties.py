# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:52:38 2023

Modified from the biorefineries constructed in [1], [2], and [3] for the production of
[1] 3-hydroxypropionic acid, [2] lactic acid, and [3] ethanol from lignocellulosic feedstocks

[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat

"""


# %% 

# =============================================================================
# Setup
# =============================================================================
from warnings import filterwarnings
filterwarnings('ignore')
import numpy as np
import pandas as pd
from biosteam import plots as bst_plots
# import biosteam as bst
from biosteam.utils import TicToc
# from biorefineries.succinic.system_succinic_adsorption_glucose import (
#     spec, succinic_sys,
# )
from biorefineries.succinic.analyses import models
from datetime import datetime
import os

# R301 = flowsheet('R301')
percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]


model = models.succinic_model
system = succinic_sys = models.succinic_sys
spec = models.spec
unit_groups_dict = models.unit_groups_dict
TEA_breakdown = models.TEA_breakdown

print('\n\n')

# %%

# =============================================================================
# Evaluate and organize results for Monte Carlo analysis
# =============================================================================
# Initiate a timer
timer = TicToc('timer')
timer.tic()

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221) # 3221
N_simulation = 10 # 2000

samples = model.sample(N=N_simulation, rule='L')
model.load_samples(samples)


###############################
# Bugfix barrage
###############################


baseline_spec = {'spec_1': spec.baseline_yield,
                 'spec_2': spec.baseline_titer,
                 'spec_3': spec.baseline_productivity,}

def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    spec.load_specifications(**baseline_spec)
    system.simulate()
    print('Loading and simulating with required specifications ...')
    spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
    system.simulate()
    
def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
    system.simulate()
    
def run_bugfix_barrage():
    try:
        reset_and_reload()
    except Exception as e:
        print(str(e))
        try:
            reset_and_switch_solver('fixedpoint')
        except Exception as e:
            print(str(e))
            try:
                reset_and_switch_solver('aitken')
            except Exception as e:
                print(str(e))
                # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                print("Bugfix barrage failed.\n")
                raise e
###############################

spec.load_spec_1 = spec.load_yield
# spec.load_spec_2 = spec.load_titer # defined in system
spec.load_spec_3 = spec.load_productivity

full_path = succinic_sys.path
# fermenter_index = full_path.index(spec.titer_inhibitor_specification.reactor)
pre_fermenter_units_path = list(spec.reactor.get_upstream_units())
pre_fermenter_units_path.reverse()
def model_specification():
    # !!!: bugfix barrage was removed for speed up, no failed evaluations found for now
    try:
        # for i in pre_fermenter_units_path: i._run()
        spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
        model._system.simulate()
    
    # print(spec.spec_1, spec.spec_2, spec.spec_3)
    # print(spec.reactor.effluent_titer)
    
    except Exception as e:
        str_e = str(e).lower()
        print('Error in model spec: %s'%str_e)
        # raise e
        if 'sugar concentration' in str_e:
            # flowsheet('AcrylicAcid').F_mass /= 1000.
            raise e
        else:
            run_bugfix_barrage()
            
model.specification = model_specification

model.exception_hook = 'warn'

baseline_initial = model.metrics_at_baseline()
baseline = pd.DataFrame(data=np.array([[i for i in baseline_initial.values],]), 
                        columns=baseline_initial.keys())

model.evaluate(notify=5, autoload=True, autosave=20, file='unfinished_evaluation')
model.table.to_excel('all_results.xlsx')

# Baseline results
baseline_end = model.metrics_at_baseline()
dateTimeObj = datetime.now()
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
file_to_save = 'succinic_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
    + '_' + str(N_simulation) + 'sims'

baseline = baseline.append(baseline_end, ignore_index=True)
baseline.index = ('initial', 'end')
baseline.to_excel(file_to_save+'_0_baseline.xlsx')

# Parameters
parameters = model.get_parameters()
index_parameters = len(model.get_baseline_sample())
parameter_values = model.table.iloc[:, :index_parameters].copy()

#%%
# TEA results
for index_TEA, i in enumerate(models.metrics):
    if i.element == 'LCA': break
index_TEA = index_parameters + index_TEA
TEA_results = model.table.iloc[:, index_parameters:index_TEA].copy()
TEA_percentiles = TEA_results.quantile(q=percentiles)

# LCA_results
LCA_results = \
    model.table.iloc[:, index_TEA::].copy()
LCA_percentiles = LCA_results.quantile(q=percentiles)

# # Spearman's rank correlation

table = model.table

model.table = model.table.dropna()

spearman_results = model.spearman()
spearman_results.columns = pd.Index([i.name_with_units for i in model.metrics])

model.table = table

# Calculate the cumulative probabilitie of each parameter
probabilities = {}
for i in range(index_parameters):
    p = parameters[i]
    p_values = parameter_values.iloc[:, 2*i]
    probabilities[p.name] = p.distribution.cdf(p_values)
    parameter_values.insert(loc=2*i+1, 
                      column=(parameter_values.iloc[:, 2*i].name[0], 'Probability'), 
                      value=probabilities[p.name],
                      allow_duplicates=True)

run_number = samples.shape[0]


#%%
'''Output to Excel'''
with pd.ExcelWriter(file_to_save+'_1_full_evaluation.xlsx') as writer:
    parameter_values.to_excel(writer, sheet_name='Parameters')
    TEA_results.to_excel(writer, sheet_name='TEA results')
    TEA_percentiles.to_excel(writer, sheet_name='TEA percentiles')
    LCA_results.to_excel(writer, sheet_name='LCA results')
    LCA_percentiles.to_excel(writer, sheet_name='LCA percentiles')
    spearman_results.to_excel(writer, sheet_name='Spearman')
    # one_p_df.to_excel(writer, sheet_name='One-parameter')
    model.table.to_excel(writer, sheet_name='Raw data')


os.remove("unfinished_evaluation") 


# %% Plot
import contourplots
# Results under uncertainty

# MPSP
MPSP_baseline = 1.40
MPSP_uncertainty = model.table.Biorefinery['Adjusted minimum selling price [$/kg]']
market_range = (2.57, 2.94)


contourplots.box_and_whiskers_plot(uncertainty_data=MPSP_uncertainty, 
                          baseline_value=MPSP_baseline, 
                          range_for_comparison=market_range,
                          values_for_comparison=[],
                          n_minor_ticks=1,
                          y_label=r"$\bfMPSP$",
                          y_units=r"$\mathrm{\$} \cdot \mathrm{kg}^{-1}$",
                          y_ticks=np.arange(0., 5., 0.5),
                          save_file=True,
                          filename=file_to_save+'_uncertainty_MPSP',
                          dpi=600,)

# GWP100a

#%%
# Spearman's rank order correlation coefficients
df_rho, df_p = model.spearman_r()
# print(df_rho['Biorefinery', 'Adjusted minimum selling price [$/kg]'])
bst_plots.plot_spearman_1d(df_rho['Biorefinery', 'Adjusted minimum selling price [$/kg]'],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='MPSP [$/kg]')


