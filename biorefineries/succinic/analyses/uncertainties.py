# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 00:40:41 2023

Modified from the biorefineries constructed in [1], [2], and [3] for the production of
[1] 3-hydroxypropionic acid, [2] lactic acid, and [3] ethanol from lignocellulosic feedstocks

[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat
"""

from warnings import filterwarnings
filterwarnings('ignore')
import numpy as np
import pandas as pd
import biosteam as bst
import thermosteam as tmo
import contourplots
print('\n\nLoading system ...')
from biorefineries.succinic.analyses import models
print('\nLoaded system.')
from datetime import datetime
from biosteam.utils import TicToc
import os

model = models.succinic_model

system = succinic_sys = models.succinic_sys
spec = models.spec
unit_groups = models.unit_groups

tea = models.succinic_tea
lca = models.succinic_LCA
get_adjusted_MSP = models.get_adjusted_MSP

#%% Bugfix barrage
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

#%% Model specification
pre_fermenter_units_path = list(spec.reactor.get_upstream_units())
pre_fermenter_units_path.reverse()
def model_specification():
    try:
        for i in pre_fermenter_units_path: i.simulate()
        spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
        model._system.simulate()
    

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

# %% 

N_simulations_per_mode = 1000 # 2000

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]

notification_interval = 5

results_dict = {'Baseline':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}, 
                            'GWP Breakdown':{}, 'FEC Breakdown':{},},
                'Uncertainty':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}},
                'Sensitivity':{'Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}}},}

modes = ['lab_batch', 'lab_fed-batch', 'pilot_batch']

parameter_distributions_filenames = ['parameter-distributions_lab-scale_batch.xlsx',
                                    'parameter-distributions_lab-scale_fed-batch.xlsx',
                                    'parameter-distributions_pilot-scale_batch.xlsx',
                                    ]

#%%

timer = TicToc('timer')
timer.tic()

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221) # 3221


for i in range(len(modes)):
    mode = modes[i]
    parameter_distributions_filename = parameter_distributions_filenames[i]
    print(f'\n\nLoading parameter distributions ({mode}) ...')
    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename)
    print(f'\nLoaded parameter distributions ({mode}).')
    
    parameters = model.get_parameters()
    
    print('\n\nLoading samples ...')
    samples = model.sample(N=N_simulations_per_mode, rule='L')
    model.load_samples(samples)
    print('\nLoaded samples.')
    
    
    model.exception_hook = 'warn'
    print('\n\nSimulating baseline ...')
    baseline_initial = model.metrics_at_baseline()
    baseline = pd.DataFrame(data=np.array([[i for i in baseline_initial.values],]), 
                            columns=baseline_initial.keys())
    
    results_dict['Baseline']['MPSP'][mode] = get_adjusted_MSP()
    results_dict['Baseline']['GWP100a'][mode] = tot_GWP = lca.GWP
    results_dict['Baseline']['FEC'][mode] = tot_FEC = lca.FEC
    
    material_GWP_breakdown = lca.material_GWP_breakdown
    
    results_dict['Baseline']['GWP Breakdown'][mode] = {
        'feedstock*': lca.FGHTP_GWP,
        'lime': material_GWP_breakdown['CalciumDihydroxide'],
        'sulfuric acid': material_GWP_breakdown['H2SO4'],
        'ammonium sulfate': material_GWP_breakdown['DiammoniumSulfate'],
        'magnesium sulfate': material_GWP_breakdown['MagnesiumSulfate'],
        'corn steep liquor': material_GWP_breakdown['CSL'],
        'other materials': material_GWP_breakdown['MEA'] + material_GWP_breakdown['NaOH'] + material_GWP_breakdown['H3PO4'],
        'natural gas (for steam generation)': lca.ng_GWP,
        'natural gas (for product drying)': material_GWP_breakdown['CH4'],
        'direct non-biogenic emissions': lca.direct_emissions_GWP,
        'net electricity production': lca.net_electricity_GWP,
        }
    
    for k, v in results_dict['Baseline']['GWP Breakdown'][mode].items():
        results_dict['Baseline']['GWP Breakdown'][mode][k] = v/tot_GWP
      
    
    material_FEC_breakdown = lca.material_FEC_breakdown
    
    results_dict['Baseline']['FEC Breakdown'][mode] = {
        'feedstock': lca.feedstock_FEC,
        'lime': material_FEC_breakdown['CalciumDihydroxide'],
        'sulfuric acid': material_FEC_breakdown['H2SO4'],
        'ammonium sulfate': material_FEC_breakdown['DiammoniumSulfate'],
        'magnesium sulfate': material_FEC_breakdown['MagnesiumSulfate'],
        'corn steep liquor': material_FEC_breakdown['CSL'],
        'other materials': material_FEC_breakdown['MEA'] + material_FEC_breakdown['NaOH'] + material_FEC_breakdown['H3PO4'],
        'natural gas (for steam generation)': lca.ng_GWP,
        'natural gas (for product drying)': material_FEC_breakdown['CH4'],
        'net electricity production': lca.net_electricity_FEC,
        }
    
    for k, v in results_dict['Baseline']['FEC Breakdown'][mode].items():
        results_dict['Baseline']['FEC Breakdown'][mode][k] = v/tot_FEC
    
    print(f"\nSimulated baseline. MPSP = ${round(results_dict['Baseline']['MPSP'][mode],2)}/kg.")
    print('\n\nEvaluating ...')
    model.evaluate(notify=notification_interval, autoload=None, autosave=None, file=None)
    print('\nFinished evaluation.')
    
    # Baseline results
    print('\n\nRe-simulating baseline ...')
    baseline_end = model.metrics_at_baseline()
    print(f"\nRe-simulated baseline. MPSP = ${round(results_dict['Baseline']['MPSP'][mode],2)}/kg.")
    dateTimeObj = datetime.now()
    minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
    file_to_save = '_succinic_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
        + '_' + str(N_simulations_per_mode) + 'sims'
    
    baseline = baseline.append(baseline_end, ignore_index=True)
    baseline.index = ('initial', 'end')
    baseline.to_excel(mode+'_'+file_to_save+'_0_baseline.xlsx')
    
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
    with pd.ExcelWriter(mode+'_'+file_to_save+'_1_full_evaluation.xlsx') as writer:
        parameter_values.to_excel(writer, sheet_name='Parameters')
        TEA_results.to_excel(writer, sheet_name='TEA results')
        TEA_percentiles.to_excel(writer, sheet_name='TEA percentiles')
        LCA_results.to_excel(writer, sheet_name='LCA results')
        LCA_percentiles.to_excel(writer, sheet_name='LCA percentiles')
        spearman_results.to_excel(writer, sheet_name='Spearman')
        # one_p_df.to_excel(writer, sheet_name='One-parameter')
        model.table.to_excel(writer, sheet_name='Raw data')
    
    
    results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price [$/kg]']
    results_dict['Uncertainty']['GWP100a'][mode] = model.table.Biorefinery['Total GWP100a [kg-CO2-eq/kg]']
    results_dict['Uncertainty']['FEC'][mode] = model.table.Biorefinery['Total FEC [kg-CO2-eq/kg]']
    
    df_rho, df_p = model.spearman_r()
    
    results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price [$/kg]']
    results_dict['Sensitivity']['Spearman']['GWP100a'][mode] = df_rho['Biorefinery', 'Total GWP100a [kg-CO2-eq/kg]']
    results_dict['Sensitivity']['Spearman']['FEC'][mode] = df_rho['Biorefinery', 'Total FEC [kg-CO2-eq/kg]']

# %% Plots
import contourplots

#%% Uncertainty

## MPSP
MPSP_uncertainty = [results_dict['Uncertainty']['MPSP'][modes[0]],
                    results_dict['Uncertainty']['MPSP'][modes[1]],
                    results_dict['Uncertainty']['MPSP'][modes[2]]]
market_range = (2.53, 2.89)
biobased_lit_MPSP_range = (1.08, 3.63)

contourplots.box_and_whiskers_plot(uncertainty_data=MPSP_uncertainty, 
                          baseline_values=[results_dict['Baseline']['MPSP'][mode] for mode in modes], 
                          baseline_locations=[1,2,3],
                          baseline_marker_colors=['w','w','w'],
                          ranges_for_comparison=[biobased_lit_MPSP_range, market_range],
                          ranges_for_comparison_colors=['#c0c1c2', '#646464'],
                          values_for_comparison=[],
                          n_minor_ticks=1,
                          show_x_ticks=True,
                          x_tick_labels=['Lab scale [batch]', 'Lab scale [fed-batch]', 'Pilot scale [batch]'],
                          y_label=r"$\bfMPSP$",
                          y_units=r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$",
                          y_ticks=np.arange(0., 5., 0.5),
                          save_file=True,
                          fig_width = 5.,
                          box_width=0.65,
                          filename=file_to_save+'_uncertainty_MPSP',
                          dpi=600,)

## GWP100a
GWP_uncertainty = [results_dict['Uncertainty']['GWP100a'][modes[0]],
                    results_dict['Uncertainty']['GWP100a'][modes[1]],
                    results_dict['Uncertainty']['GWP100a'][modes[2]]]


biobased_lit_GWP_values = [1, 2, 3] #!!!
contourplots.box_and_whiskers_plot(uncertainty_data=GWP_uncertainty, 
                          baseline_values=[results_dict['Baseline']['GWP100a'][mode] for mode in modes], 
                          baseline_locations=[1,2,3],
                          baseline_marker_colors=['w','w','w'],
                          # ranges_for_comparison=[biobased_lit_MPSP_range, market_range],
                          # ranges_for_comparison_colors=['#c0c1c2', '#646464'],
                          values_for_comparison=biobased_lit_GWP_values,
                          n_minor_ticks=1,
                          show_x_ticks=True,
                          x_tick_labels=['Lab scale [batch]', 'Lab scale [fed-batch]', 'Pilot scale [batch]'],
                          y_label=r"$\bfGWP-100a$",
                          y_units=r"$\mathrm{kg CO}^{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$",
                          y_ticks=np.arange(0., 5., 0.5),
                          save_file=True,
                          fig_width = 5.,
                          box_width=0.65,
                          filename=file_to_save+'_uncertainty_GWP100a',
                          dpi=600,)

## FEC
FEC_uncertainty = [results_dict['Uncertainty']['FEC'][modes[0]],
                    results_dict['Uncertainty']['FEC'][modes[1]],
                    results_dict['Uncertainty']['FEC'][modes[2]]]


biobased_lit_FEC_values = [1, 2, 3] #!!!
contourplots.box_and_whiskers_plot(uncertainty_data=FEC_uncertainty, 
                          baseline_values=[results_dict['Baseline']['FEC'][mode] for mode in modes], 
                          baseline_locations=[1,2,3],
                          baseline_marker_colors=['w','w','w'],
                          # ranges_for_comparison=[biobased_lit_MPSP_range, market_range],
                          # ranges_for_comparison_colors=['#c0c1c2', '#646464'],
                          values_for_comparison=biobased_lit_FEC_values,
                          n_minor_ticks=1,
                          show_x_ticks=True,
                          x_tick_labels=['Lab scale [batch]', 'Lab scale [fed-batch]', 'Pilot scale [batch]'],
                          y_label=r"$\bfFEC$",
                          y_units=r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$",
                          y_ticks=np.arange(-10, 10., 1.),
                          save_file=True,
                          fig_width = 5.,
                          box_width=0.65,
                          filename=file_to_save+'_uncertainty_FEC',
                          dpi=600,)


#%% TEA breakdown figure
df_TEA_breakdown = bst.UnitGroup.df_from_groups(
    unit_groups, fraction=True,
    scale_fractions_to_positive_values=False,
)


# df_TEA_breakdown['Net electricity production']*=-1
# df_TEA_breakdown = df_TEA_breakdown.rename(columns={'Net electricity production': 'Net electricity demand'})

contourplots.stacked_bar_plot(dataframe=df_TEA_breakdown, 
                 # y_ticks=[-200, -175, -150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175], 
                 y_ticks=[-150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250], 
                 y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
                 y_units = "%", 
                 colors=['#7BBD84', '#F7C652', '#63C6CE', '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#F8858A', 'red', 'magenta'],
                 filename='TEA_breakdown_stacked_bar_plot')

#%%

#%% LCA breakdown figures
# GWP
df_GWP_breakdown = pd.DataFrame([list(results_dict['Baseline']['FEC Breakdown'][modes[2]].keys()), list(results_dict['Baseline']['FEC Breakdown'][modes[2]].values())])


# df_GWP_breakdown['Net electricity production']*=-1
# df_GWP_breakdown = df_GWP_breakdown.rename(columns={'Net electricity production': 'Net electricity demand'})

contourplots.stacked_bar_plot(dataframe=df_GWP_breakdown, 
                 # y_ticks=[-200, -175, -150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175], 
                 y_ticks=[-150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250], 
                 y_label=r"$\bfGWP100a$" + r"$\bfBreakdown$",  
                 y_units = "%", 
                 colors=['#7BBD84', '#F7C652', '#63C6CE', '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#F8858A', 'magenta'],
                 filename='GWP_breakdown_stacked_bar_plot')

# FEC
df_FEC_breakdown = pd.DataFrame([list(results_dict['Baseline']['FEC Breakdown'][modes[2]].keys()), list(results_dict['Baseline']['FEC Breakdown'][modes[2]].values())])


# df_FEC_breakdown['Net electricity production']*=-1
# df_FEC_breakdown = df_FEC_breakdown.rename(columns={'Net electricity production': 'Net electricity demand'})

contourplots.stacked_bar_plot(dataframe=df_FEC_breakdown, 
                 # y_ticks=[-200, -175, -150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175], 
                 y_ticks=[-150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250], 
                 y_label=r"$\bfFEC$" + r"$\bfBreakdown$", 
                 y_units = "%", 
                 colors=['#7BBD84', '#F7C652', '#63C6CE', '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#F8858A', 'magenta'],
                 filename='FEC_breakdown_stacked_bar_plot')

#%%
# Spearman's rank order correlation coefficients
bst_plots = bst.plots

bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['MPSP'][modes[0]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='Lab scale [batch] - MPSP [$/kg]')
bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['GWP100a'][modes[1]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='Lab scale [batch] - GWP100a [kg-CO2-eq./kg]')
bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['FEC'][modes[2]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='Lab scale [batch] - FEC [MJ/kg]')

bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['MPSP'][modes[0]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='Lab scale [fed-batch] - MPSP [$/kg]')
bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['GWP100a'][modes[1]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='Lab scale [fed-batch] - GWP100a [kg-CO2-eq./kg]')
bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['FEC'][modes[2]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='Lab scale [fed-batch] - FEC [MJ/kg]')

bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['MPSP'][modes[0]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='Pilot scale [batch] - MPSP [$/kg]')
bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['GWP100a'][modes[1]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='Pilot scale [batch] - GWP100a [kg-CO2-eq./kg]')
bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['FEC'][modes[2]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='Pilot scale [batch] - FEC [MJ/kg]')

