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
# from biorefineries
# from biorefineries import succinic
from biorefineries import succinic

models = succinic.get_models()
# from . import models

print('\nLoaded system.')
from datetime import datetime
from biosteam.utils import TicToc
import os

chdir = os.chdir
succinic_filepath = succinic.__file__.replace('\\__init__.py', '')
succinic_results_filepath = succinic_filepath + '\\analyses\\results\\'
model = models.succinic_model

system = succinic_sys = models.succinic_sys
spec = models.spec
unit_groups = models.unit_groups

tea = models.succinic_tea
lca = models.succinic_LCA
get_adjusted_MSP = models.get_adjusted_MSP

# %% 

N_simulations_per_mode = 2000 # 2000

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]

notification_interval = 50

results_dict = {'Baseline':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}, 
                            'GWP Breakdown':{}, 'FEC Breakdown':{},},
                'Uncertainty':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}},
                'Sensitivity':{'Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}}},}

modes = [
            'lab_batch',
           'lab_fed-batch', 
           'pilot_batch',
         ]

parameter_distributions_filenames = [
                                    'parameter-distributions_lab-scale_batch.xlsx',
                                    'parameter-distributions_lab-scale_fed-batch.xlsx',
                                    'parameter-distributions_pilot-scale_batch.xlsx',
                                    ]

#%%

timer = TicToc('timer')
timer.tic()

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221) # 3221


for i in range(len(modes)):
    # ## Change working directory to biorefineries\\succinic
    # chdir(succinic.__file__.replace('\\__init__.py', ''))
    # ##
    mode = modes[i]
    parameter_distributions_filename = succinic_filepath+\
        '\\analyses\\parameter_distributions\\'+parameter_distributions_filenames[i]
    print(f'\n\nLoading parameter distributions ({mode}) ...')
    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename)
    print(f'\nLoaded parameter distributions ({mode}).')
    
    parameters = model.get_parameters()
    
    print('\n\nLoading samples ...')
    samples = model.sample(N=N_simulations_per_mode, rule='L')
    model.load_samples(samples)
    print('\nLoaded samples.')
    
    # ## Change working directory to biorefineries\\succinic\\analyses\\results
    # chdir(succinic.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
    # ##
    
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
        # 'sulfuric acid': material_GWP_breakdown['H2SO4'],
        'ammonium sulfate': material_GWP_breakdown['DiammoniumSulfate'],
        'magnesium sulfate': material_GWP_breakdown['MagnesiumSulfate'],
        'corn steep liquor': material_GWP_breakdown['CSL'],
        'other materials': material_GWP_breakdown['MEA'] + material_GWP_breakdown['NaOH'] + material_GWP_breakdown['H3PO4'],
        'natural gas\n(for steam generation)': lca.ng_GWP,
        'natural gas\n(product drying)': material_GWP_breakdown['CH4'],
        'net electricity': lca.net_electricity_GWP,
        'direct non-biogenic\nemissions': lca.direct_emissions_GWP,
        }
    
    tot_positive_GWP = sum([v for v in results_dict['Baseline']['GWP Breakdown'][mode].values() if v>0])
    for k, v in results_dict['Baseline']['GWP Breakdown'][mode].items():
        results_dict['Baseline']['GWP Breakdown'][mode][k] = v/tot_positive_GWP
      
    
    material_FEC_breakdown = lca.material_FEC_breakdown
    
    results_dict['Baseline']['FEC Breakdown'][mode] = {
        'feedstock': lca.feedstock_FEC,
        'lime': material_FEC_breakdown['CalciumDihydroxide'],
        # 'sulfuric acid': material_FEC_breakdown['H2SO4'],
        'ammonium sulfate': material_FEC_breakdown['DiammoniumSulfate'],
        'magnesium sulfate': material_FEC_breakdown['MagnesiumSulfate'],
        'corn steep liquor': material_FEC_breakdown['CSL'],
        'other materials': material_FEC_breakdown['MEA'] + material_FEC_breakdown['NaOH'] + material_FEC_breakdown['H3PO4'],
        'natural gas\n(for steam generation)': lca.ng_GWP,
        'natural gas\n(for product drying)': material_FEC_breakdown['CH4'],
        'net electricity': lca.net_electricity_FEC,
        }
    tot_positive_FEC = sum([v for v in results_dict['Baseline']['FEC Breakdown'][mode].values() if v>0])
    for k, v in results_dict['Baseline']['FEC Breakdown'][mode].items():
        results_dict['Baseline']['FEC Breakdown'][mode][k] = v/tot_positive_FEC
    
    if spec.reactor.base_neutralizes_product: # sulfuric acid for acidulation
        results_dict['Baseline']['GWP Breakdown'][mode]['sulfuric acid'] = material_GWP_breakdown['H2SO4']
        results_dict['Baseline']['FEC Breakdown'][mode]['sulfuric acid'] = material_FEC_breakdown['H2SO4']
        
        
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
    file_to_save = succinic_results_filepath+\
        '_succinic_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
        + '_' + str(N_simulations_per_mode) + 'sims'
    
    baseline = baseline.append(baseline_end, ignore_index=True)
    baseline.index = ('initial', 'end')
    baseline.to_excel(file_to_save+'_'+mode+'_0_baseline.xlsx')
    
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
    with pd.ExcelWriter(file_to_save+'_'+mode+'_1_full_evaluation.xlsx') as writer:
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

#%% Clean up NaN values for plotting
metrics = ['MPSP', 'GWP100a', 'FEC']
tot_NaN_vals_dict = results_dict['Errors'] = {metric: {mode: 0 for mode in modes} for metric in metrics}
for mode in modes:
    for metric in metrics:
        # median_val = np.median(results_dict['Uncertainty'][metric][mode])
        median_val = 1.5
        for i in range(len(results_dict['Uncertainty'][metric][mode])):
            if np.isnan(results_dict['Uncertainty'][metric][mode][i]):
                results_dict['Uncertainty'][metric][mode][i] = median_val
                tot_NaN_vals_dict[metric][mode] += 1
# %% Plots
import contourplots

MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"
GWP_units = r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$"
FEC_units = r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$"
#%% Uncertainty

scenario_name_labels = ['Lab. batch', 
                        'Lab. fed-batch', 
                        'Pilot batch']

def get_small_range(num, offset):
    return(num-offset, num+offset)
#%% MPSP
MPSP_uncertainty = [results_dict['Uncertainty']['MPSP'][modes[0]],
                    results_dict['Uncertainty']['MPSP'][modes[1]],
                    results_dict['Uncertainty']['MPSP'][modes[2]]]
market_range = (2.53, 2.89)
biobased_lit_MPSP_range = (1.08, 3.63)

contourplots.box_and_whiskers_plot(uncertainty_data=MPSP_uncertainty, 
                          baseline_values=[results_dict['Baseline']['MPSP'][mode] for mode in modes],
                          baseline_marker_shapes=["p", "s", "D"],
                          baseline_marker_sizes=[10, 6, 6,],
                          baseline_locations=[1,2,3],
                          baseline_marker_colors=['w','w','w'],
                          boxcolor="#A97802",
                          ranges_for_comparison=[biobased_lit_MPSP_range, market_range],
                          ranges_for_comparison_colors=['#c0c1c2', '#646464'],
                          values_for_comparison=[],
                          n_minor_ticks=1,
                          show_x_ticks=True,
                          x_tick_labels=scenario_name_labels,
                          x_tick_wrap_width=6,
                          y_label=r"$\bfMPSP$",
                          y_units=MPSP_units,
                          y_ticks=np.arange(0., 5.5, 0.5),
                          save_file=True,
                          fig_height=5.5,
                          fig_width = 3.,
                          box_width=0.65,
                          filename=file_to_save+'_uncertainty_MPSP',
                          dpi=600,)

#%% GWP100a

biobased_GWPs = [1.83, 2.20, 2.37]
fossilbased_GWPs = [3.27, 3.43, 10.3, 12.1]

GWP_uncertainty = [results_dict['Uncertainty']['GWP100a'][modes[0]],
                    results_dict['Uncertainty']['GWP100a'][modes[1]],
                    results_dict['Uncertainty']['GWP100a'][modes[2]]]


biobased_lit_GWP_values = [1, 2, 3] #!!!
contourplots.box_and_whiskers_plot(uncertainty_data=GWP_uncertainty, 
                          baseline_values=[results_dict['Baseline']['GWP100a'][mode] for mode in modes], 
                          baseline_marker_shapes=["p", "s", "D"],
                          baseline_marker_sizes=[10, 6, 6,],
                          baseline_locations=[1,2,3],
                          baseline_marker_colors=['w','w','w'],
                          boxcolor='#607429',
                          ranges_for_comparison=[get_small_range(i, 0.005) for i in biobased_GWPs+fossilbased_GWPs],
                          ranges_for_comparison_colors=['#c0c1c2' for i in range(len(biobased_GWPs))] +\
                                                       ['#646464' for i in range(len(fossilbased_GWPs))],
                          values_for_comparison=biobased_lit_GWP_values,
                          n_minor_ticks=1,
                          show_x_ticks=True,
                          x_tick_labels=scenario_name_labels,
                          x_tick_wrap_width=6,
                          # y_label=r"$\bfGWP-100a$",
                          y_label=r"$\mathrm{\bfGWP}_{\bf100}$",
                          y_units=GWP_units,
                          y_ticks=np.arange(0., 14., 1.),
                          save_file=True,
                          fig_height=5.5,
                          fig_width = 3.,
                          box_width=0.65,
                          filename=file_to_save+'_uncertainty_GWP100a',
                          dpi=600,)

#%% FEC

biobased_FECs = [26, 27.7, 32.7]
fossilbased_FECs = [59.2, 60.8, 112, 124]

FEC_uncertainty = [results_dict['Uncertainty']['FEC'][modes[0]],
                    results_dict['Uncertainty']['FEC'][modes[1]],
                    results_dict['Uncertainty']['FEC'][modes[2]]]


biobased_lit_FEC_values = [1, 2, 3] #!!!
contourplots.box_and_whiskers_plot(uncertainty_data=FEC_uncertainty, 
                          baseline_values=[results_dict['Baseline']['FEC'][mode] for mode in modes], 
                          baseline_marker_shapes=["p", "s", "D"],
                          baseline_marker_sizes=[10, 6, 6,],
                          baseline_locations=[1,2,3],
                          baseline_marker_colors=['w','w','w'],
                          boxcolor='#A100A1',
                          ranges_for_comparison=[get_small_range(i, 0.061) for i in biobased_FECs+fossilbased_FECs],
                          ranges_for_comparison_colors=['#c0c1c2' for i in range(len(biobased_FECs))] +\
                                                       ['#646464' for i in range(len(fossilbased_FECs))],
                          values_for_comparison=biobased_lit_FEC_values,
                          n_minor_ticks=1,
                          show_x_ticks=True,
                          x_tick_labels=scenario_name_labels,
                          x_tick_wrap_width=6,
                          y_label=r"$\bfFEC$",
                          y_units=FEC_units,
                          y_ticks=np.arange(-30, 140., 10.),
                          save_file=True,
                          fig_height=5.5,
                          fig_width = 3.,
                          box_width=0.65,
                          filename=file_to_save+'_uncertainty_FEC',
                          dpi=600,)


#%% TEA breakdown figure
df_TEA_breakdown = bst.UnitGroup.df_from_groups(
    unit_groups, fraction=True,
    scale_fractions_to_positive_values=True,
)


# df_TEA_breakdown['Net electricity production']*=-1
# df_TEA_breakdown = df_TEA_breakdown.rename(columns={'Net electricity production': 'Net electricity demand'})

df_TEA_breakdown = df_TEA_breakdown.rename(columns={'Material cost': 'Operating cost'})
contourplots.stacked_bar_plot(dataframe=df_TEA_breakdown, 
                 # y_ticks=[-200, -175, -150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175], 
                 y_ticks=[-40, -20, 0, 20, 40, 60, 80, 100], 
                 y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
                 y_units = "%", 
                 colors=['#7BBD84', '#F7C652', '#63C6CE', '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#D1C0E1', '#F8858A', '#F8858A', ],
                 # 'red', 'magenta'],
                 filename=file_to_save+'TEA_breakdown_stacked_bar_plot',
                 fig_height=5.5*1.1777*0.94*1.0975)

#%%

#%% LCA breakdown figures
# GWP
temp_GWP_breakdown_dict = results_dict['Baseline']['GWP Breakdown'][modes[2]]
GWP_breakdown_dict = {
                        # 'areas': list(temp_GWP_breakdown_dict.keys()), 
                      'contributions': [100*i for i in list(temp_GWP_breakdown_dict.values())]}
GWP_breakdown_list = [100*v for k, v in temp_GWP_breakdown_dict.items()]
df_GWP_breakdown = pd.DataFrame(GWP_breakdown_list,
                                            index=list(temp_GWP_breakdown_dict.keys()),
                                          # orient='index',
                                           # columns=['contributions'],
                                          )
# df_GWP_breakdown['contributions']=df_GWP_breakdown['contributions'].astype(float)


# df_GWP_breakdown['Net electricity production']*=-1
# df_GWP_breakdown = df_GWP_breakdown.rename(columns={'Net electricity production': 'Net electricity demand'})

contourplots.stacked_bar_plot(dataframe=df_GWP_breakdown, 
                  y_ticks=[-40, -20, 0, 20, 40, 60, 80, 100], 
                  # y_ticks=[-400, -300, -200, -100, 0, 100, 200, 300, 400], 
                 # y_ticks = []
                 # y_label=r"$\bfGWP-100a $" +" "+ r"$\bfBreakdown$",  
                 y_label =r"$\mathrm{\bfGWP}_{\bf100}$" +" "+ r"$\bfBreakdown$",
                 y_units = "%", 
                  colors=['#E1F8C0', '#8FAE3E', '#607429', 
                          ],
                  hatch_patterns=('\\', '//', 'x',  '|',),
                  # '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#F8858A', 'red', 'magenta'],
                 filename=file_to_save+'GWP_breakdown_stacked_bar_plot',
                 fig_width=2,
                 fig_height=5.5*1.1777*0.94)

# FEC
temp_FEC_breakdown_dict = results_dict['Baseline']['FEC Breakdown'][modes[2]]
FEC_breakdown_dict = {
                        # 'areas': list(temp_FEC_breakdown_dict.keys()), 
                      'contributions': [100*i for i in list(temp_FEC_breakdown_dict.values())]}
FEC_breakdown_list = [100*v for k, v in temp_FEC_breakdown_dict.items()]
df_FEC_breakdown = pd.DataFrame(FEC_breakdown_list,
                                            index=list(temp_FEC_breakdown_dict.keys()),
                                          # orient='index',
                                           # columns=['contributions'],
                                          )

# df_FEC_breakdown['Net electricity production']*=-1
# df_FEC_breakdown = df_FEC_breakdown.rename(columns={'Net electricity production': 'Net electricity demand'})

contourplots.stacked_bar_plot(dataframe=df_FEC_breakdown, 
                 # y_ticks=[-200, -175, -150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175], 
                 y_ticks=[-100, -75, -50, -25, 0, 25, 50, 75, 100], 
                 y_label=r"$\bfFEC$" +" "+ r"$\bfBreakdown$", 
                 y_units = "%", 
                 # colors=['#7BBD84', '#F7C652', '#63C6CE', '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#F8858A', 'magenta'],
                 colors=['#FEC1FE', '#FF80FF', '#A100A1', 
                         ],
                 hatch_patterns=('\\', '//', 'x',  '|',),
                 filename=file_to_save+'FEC_breakdown_stacked_bar_plot',
                 fig_width=2,
                 fig_height=5.5*1.1777*0.94)

#%% Spearman's rank order correlation coefficients
from matplotlib import pyplot as plt
chdir(succinic_results_filepath)
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = "7.5"

bst_plots = bst.plots

rho = r"$\mathrm{\rho}}$"

mode = modes[0]
fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['MPSP'][modes[0]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='MPSP - Lab scale [batch] '+"["+MPSP_units+"]", color="#A97802",
                           xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i)
fig[0].savefig(mode+'MPSP-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)

fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['GWP100a'][modes[1]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='GWP100 - Lab scale [batch] '+"["+GWP_units+"]", color='#607429',
                           xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i)
fig[0].savefig(mode+'GWP-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)
fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['FEC'][modes[2]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='FEC - Lab scale [batch] '+"["+FEC_units+"]", color='#A100A1',
                           xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i)
fig[0].savefig(mode+'FEC-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)

mode = modes[1]
fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['MPSP'][modes[0]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='MPSP - Lab scale [fed-batch] '+"["+MPSP_units+"]", color="#A97802",
                           xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i)
fig[0].savefig(mode+'MPSP-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)

fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['GWP100a'][modes[1]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='GWP100 - Lab scale [fed-batch] '+"["+GWP_units+"]", color='#607429',
                           xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i)
fig[0].savefig(mode+'GWP-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)

fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['FEC'][modes[2]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='FEC - Lab scale [fed-batch] '+"["+FEC_units+"]", color='#A100A1',
                           xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i)
fig[0].savefig(mode+'FEC-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)

mode = modes[2]
fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['MPSP'][modes[0]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='MPSP - Pilot scale [batch] '+"["+MPSP_units+"]", color="#A97802",
                           xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i)
fig[0].savefig(mode+'MPSP-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)

fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['GWP100a'][modes[1]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='GWP100 - Pilot scale [batch]a '+"["+GWP_units+"]", color='#607429',
                           xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i)
fig[0].savefig(mode+'GWP-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)

fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['FEC'][modes[2]],
                           index=[i.element_name + ': ' + i.name for i in model.parameters],
                           name='FEC - Pilot scale [batch] '+"["+FEC_units+"]", color='#A100A1',
                           xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i)
fig[0].savefig(mode+'FEC-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)

#%% Dunn et al. 2015 data - extracted using WebPlotDigitizer

# GHG-100
# Bar0, 12.126849894291754
# Bar1, 1.8266384778012683
# Bar2, 3.247357293868922

# FEC
# Bar0, 112.33870967741936
# Bar1, 42.05645161290324
# Bar2, 27.66129032258065


#%% GHG-100

## GREET 2022
# succinic acid, succinic acid bioproduct from sugars (corn stover): 0.7036 kg-CO2-eq/kg-SA
# incl. EOL GWP (1.4937 kg-CO2-eq/kg-SA): 2.20 kg-CO2-eq/kg-SA

## Dunn et al. 2015 (all cradle-to-grave)
# fossil-derived: 12.13 kg-CO2-eq/kg-SA
# corn stover-based, EDI: 3.25 kg-CO2-eq/kg-SA
# corn stover-based, LLE: 1.83 kg-CO2-eq/kg-SA

## Cok et al. 2014
# corn-based w/ dir. crystallization: 0.88 kg-CO2-eq/kg-SA
# incl. EOL GWP (1.4937 kg-CO2-eq/kg-SA): 2.37 kg-CO2-eq/kg-SA

# fossil-based (3 alternative processes): 1.78, 1.94, 8.82 kg-CO2-eq/kg-SA
# incl. EOL GWP (1.4937 kg-CO2-eq/kg-SA): 3.27,  3.43, 10.31 kg-CO2-eq/kg-SA


#%% FEC
## GREET 2022
# succinic acid, succinic acid bioproduct from sugars (corn stover): 26 MJ-eq/kg-SA (of which all 26 MJ-eq/kg-SA from natural gas)

## Dunn et al. 2015 (all cradle-to-gate)
# fossil-derived: 112 MJ-eq/kg-SA
# corn stover-based, EDI: 42.1 MJ-eq/kg-SA
# corn stover-based, LLE: 27.7 MJ-eq/kg-SA

## Cok et al. 2014
# corn-based w/ dir. crystallization: 32.7 MJ-eq/kg-SA

# fossil-based (3 alternative processes): 60.8, 59.2, 124.3 MJ-eq/kg-SA



## ecoinvent 3.8
# succinic acid production, GLO: 67.463 MJ-eq/kg-SA
# market for succinic acid, GLO: 68.19 MJ-eq/kdg-SA

#%% IPCC 2013 GWP100a

## ecoinvent 3.8
# succinic acid production, GLO: 3.182 kgCO2-eq/kg-SA
# incl. EOL GWP (1.4937 kg-CO2-eq/kg-SA): 4.6757 kg-CO2-eq/kg-SA

# market for succinic acid, GLO: 3.2322 kgCO2-eq/kg-SA
# incl. EOL GWP (1.4937 kg-CO2-eq/kg-SA): 4.7259 kg-CO2-eq/kg-SA

