#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import pandas as pd

#%%
# def run_IBO_uncertainty_analysis(
                        # modes=['A', 
                        #         # 'B', 'C', 'D',
                        #        ],
                        # N_simulations_per_mode=20,
                        # notification_interval=1,
                        # plot_TOC_fig=True,
                         # ):

from warnings import filterwarnings
filterwarnings('ignore')
import numpy as np
import pandas as pd
import contourplots
import biosteam as bst
print('\n\nLoading system ...')
# from biorefineries
# from biorefineries import isobutanol
from biorefineries import isobutanol
from biorefineries.isobutanol.models import models_EtOH_IBO_corn as models
# models = isobutanol.models
# from . import models

print('\nLoaded system.')
from datetime import datetime
from biosteam.utils import TicToc
import os

dateTimeObj = datetime.now()

chdir = os.chdir
IBO_filepath = isobutanol.__file__.replace('\\__init__.py', '')
IBO_results_filepath = IBO_filepath + '\\analyses\\results\\'
model = models.IBO_model

system = IBO_sys = models.IBO_sys
unit_groups = models.unit_groups

tea = models.IBO_tea
get_adjusted_MSP = models.get_adjusted_MSP
# per_kg_KSA_to_per_kg_SA = models.per_kg_KSA_to_per_kg_SA

f = bst.main_flowsheet

#%%
modes=[
       'A', 
        # 'B', 
        # 'C', 'D',
       ]
N_simulations_per_mode=1000
notification_interval=20
plot_TOC_fig=True

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]

results_dict = {'Baseline':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}, 
                            'GWP Breakdown':{}, 'FEC Breakdown':{},},
                'Uncertainty':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}},
                'Sensitivity':{'Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}},
                               'p-val Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}}},}


scenario_names =\
                {
                   'A': 'current state-of-technology',
                   'B': 'potential isobutanol production',
                }
                
parameter_distributions_filenames = {i: 'parameter-distributions_corn_IBO_EtOH_' + i + '.xlsx' for i in modes}


#%%

timer = TicToc('timer')
timer.tic()

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221) # 3221

# def load_additional_params():
   
for i in range(len(modes)):
    # ## Change working directory to biorefineries\\isobutanol
    # chdir(isobutanol.__file__.replace('\\__init__.py', ''))
    # ##
    mode = modes[i]
    parameter_distributions_filename = IBO_filepath+\
        '\\analyses\\full\\parameter_distributions\\'+parameter_distributions_filenames[mode]
    print(f'\n\nLoading parameter distributions ({mode}) ...')
    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename, models.namespace_dict)
    
    # load_additional_params()
    print(f'\nLoaded parameter distributions ({mode}).')
    
    parameters = model.get_parameters()
    
    print('\n\nLoading samples ...')
    samples = model.sample(N=N_simulations_per_mode, rule='L')
    model.load_samples(samples)
    print('\nLoaded samples.')
    
    # ## Change working directory to biorefineries\\isobutanol\\analyses\\results
    # chdir(isobutanol.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
    # ##
    
    model.exception_hook = 'warn'
    print('\n\nSimulating baseline ...')
    
    baseline_initial = model.metrics_at_baseline()
    model.specification()
    baseline_initial = model.metrics_at_baseline()
    
    baseline = pd.DataFrame(data=np.array([[i for i in baseline_initial.values],]), 
                            columns=baseline_initial.keys())
    
    results_dict['Baseline']['MPSP'][mode] = get_adjusted_MSP()
        
    print(f"\nSimulated baseline. MPSP = ${round(results_dict['Baseline']['MPSP'][mode],2)}/kg.")
    
    #%%
    print('\n\nEvaluating ...')
    model.evaluate(notify=notification_interval, autoload=None, autosave=None, file=None)
    print('\nFinished evaluation.')
    
    # Baseline results
    print('\n\nRe-simulating baseline ...')
    
    baseline_end = model.metrics_at_baseline()
    model.specification()
    baseline_end = model.metrics_at_baseline()
    
    print(f"\nRe-simulated baseline. MPSP = ${round(get_adjusted_MSP(),2)}/kg.")
    
    minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
    file_to_save = IBO_results_filepath+\
        '_IBO_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
        + '_' + str(modes) + '_' + str(N_simulations_per_mode) + 'sims'
    
    # baseline = baseline.append(baseline_end, ignore_index=True)
    baseline = pd.concat([baseline, pd.DataFrame([baseline_end])], ignore_index=True)
    # print(baseline)
    
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
    index_TEA = index_parameters + index_TEA\
                                          + 1
    TEA_results = model.table.iloc[:, index_parameters:index_TEA].copy()
    TEA_percentiles = TEA_results.quantile(q=percentiles)
 
    table = model.table
    
    model.table = model.table.dropna()
    
    spearman_results, spearman_p_values = model.spearman_r()
    
    spearman_results.columns = pd.Index([i.name_with_units for i in model.metrics])
    spearman_p_values.columns = pd.Index([i.name_with_units for i in model.metrics])
    
    # model.table = table
    
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
    
    print('\n\nSaving raw results ...')
    with pd.ExcelWriter(file_to_save+'_'+mode+'_1_full_evaluation.xlsx') as writer:
        parameter_values.to_excel(writer, sheet_name='Parameters')
        TEA_results.to_excel(writer, sheet_name='TEA results')
        TEA_percentiles.to_excel(writer, sheet_name='TEA percentiles')
        spearman_results.to_excel(writer, sheet_name='Spearman')
        spearman_p_values.to_excel(writer, sheet_name='Spearman p-values')
        # one_p_df.to_excel(writer, sheet_name='One-parameter')
        model.table.to_excel(writer, sheet_name='Raw data')
    
    
    results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price [$/kg IBO]']
    
    df_rho, df_p = model.spearman_r()
    
    results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price [$/kg IBO]']
    
    results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price [$/kg IBO]']
    
    print('\n\nSaved raw results.')
    print('---------------------------------\n\n')
    
#%% Clean up NaN values for plotting
metrics = ['MPSP', 
           # 'GWP100a', 
           # 'FEC',
           ]
tot_NaN_vals_dict = results_dict['Errors'] = {metric: {mode: 0 for mode in modes} for metric in metrics}
for mode in modes:
    for metric in metrics:
        # median_val = np.median(results_dict['Uncertainty'][metric][mode])
        median_val = np.median(results_dict['Uncertainty']['MPSP'][mode][np.where(~np.isnan(results_dict['Uncertainty']['MPSP'][mode]))[0]])
        for i in range(len(results_dict['Uncertainty'][metric][mode])):
            if np.isnan(results_dict['Uncertainty'][metric][mode][i]):
                results_dict['Uncertainty'][metric][mode][i] = median_val
                tot_NaN_vals_dict[metric][mode] += 1

# %% Plots
print('\n\nCreating and saving plots ...')

MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"
GWP_units = r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$"
FEC_units = r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$"
#%% More plot utils
def get_small_range(num, offset):
    return(num-offset, num+offset)

baseline_marker_shapes=["D", "^", "s", "h"]
baseline_marker_sizes=[6, 8, 6, 10]
# baseline_marker_shapes=["D",]
# baseline_marker_sizes=[6,]
n_cols_subplots=2 if len(modes)==1 else 1
#%% MPSP uncertainty and cost & utility breakdown

# modes = ['A',]
file_to_save = IBO_results_filepath+\
    '_IBO_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
    + '_' + str(modes) + '_' + str(N_simulations_per_mode) + 'sims'


MPSP_uncertainty = [results_dict['Uncertainty']['MPSP'][mode]
                    for mode in modes
                    ]

# # search page for high end: https://www.alibaba.com/trade/search?spm=a2700.galleryofferlist.0.0.2a995827YzqZVg&fsb=y&IndexArea=product_en&assessmentCompany=true&keywords=590-00-1+sorbate&productTag=1200000228&ta=y&tab=all&
# SA_market_range=np.array([
#                           6.74, # 2019 global high end from Sorbic Acid Market, Transparency Market Research
#                           6.50 * 1.3397087, # $6.50/kg-potassium-sorbate from https://www.alibaba.com/product-detail/Lifecare-Supply-Potassium-Sorbate-High-Quality_1600897125355.html?spm=a2700.galleryofferlist.p_offer.d_title.1bc15827eAs1TL&s=p
#                           ]) 

# IBO_maximum_viable_market_range = SA_market_range / theoretical_max_g_IBO_per_g_SA

market_range = np.array([0.64, 1.00])


# biobased_lit_MPSP_range = (1.08, 3.63)
MPSP_box_width = 1.05 if len(modes)==1 else 0.45
fig, axs = contourplots.box_and_whiskers_plot(uncertainty_data=MPSP_uncertainty, 
                          baseline_values=[results_dict['Baseline']['MPSP'][mode] for mode in modes],
                          baseline_marker_shapes=baseline_marker_shapes,
                          baseline_marker_sizes=baseline_marker_sizes,
                          baseline_locations=[i+1 for i in range(len(modes))],
                          baseline_marker_colors=['w' for mode in modes],
                          boxcolor="#A97802",
                          ranges_for_comparison=[market_range,],
                          ranges_for_comparison_colors=['#c0c1c2', 
                                                        '#646464',
                                                        ],
                          values_for_comparison=[],
                          n_minor_ticks=4,
                          show_x_ticks=True,
                          x_tick_labels=[scenario_names[i] for i in modes],
                          x_tick_wrap_width=14,
                          y_label=r"$\bfMPSP$",
                          y_units=MPSP_units,
                          y_ticks=np.arange(0., 1.26, 0.25),
                          save_file=False,
                          fig_height=5.5,
                           # fig_width = 2.,
                            fig_width = 5.5,
                            box_width=MPSP_box_width,
                            # box_width=1.05,
                          filename=file_to_save+'_uncertainty_MPSP',
                          dpi=600,
                          xticks_fontsize = 17,
                          ylabel_fontsize = 18,
                          yticks_fontsize = 17,
                          default_fontsize = 17,
                            n_cols_subplots = n_cols_subplots,
                            # n_cols_subplots = 1,
                          width_ratios = [1*6.72/8.15,6],
                           xticks_fontcolor = 'w',
                          # xlabelpad = 5,
                          # ylabelpad=5,
                          )
# 3.59, 7.18
# 3.75, 7.5

if len(modes)==1:
    ###### change operating cost unit labels $/h to MM$/y
    for i in unit_groups:
        for j in i.metrics:
            if j.name == 'Operating cost':
                j.units = r"$\mathrm{MM\$}$" + '\u00b7y\u207b\u00b9'
        if i.name=='natural gas (for product drying)':
            i.name='natural gas\n(for product drying)'
    ######
    
    df_TEA_breakdown = bst.UnitGroup.df_from_groups(
        unit_groups, fraction=True,
        scale_fractions_to_positive_values=True,
    )
    
    # totals=[sum([ui.metrics[i]() for ui in unit_groups])
    #         for i in range(len(unit_groups[0].metrics))]
    
    totals=[]
    metrics = unit_groups[0].metrics
    for i in range(len(metrics)):
        curr_total = 0.
        for ui in unit_groups:
            curr_total += ui.metrics[i]()
        if metrics[i].name=='Operating cost':
            # change total operating cost from $/h to MM$/y
            curr_total *= system.TEA.operating_hours/1e6
        totals.append(curr_total)
    
    
    
    
    
    contourplots.stacked_bar_plot(dataframe=df_TEA_breakdown, 
                     y_ticks = [-50, -25, 0, 25, 50, 75, 100],
                     y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
                     y_units = "%", 
                     colors=['#7BBD84', 
                             '#E58835', 
                             '#F7C652', 
                             '#63C6CE', 
                             # '#b00000', 
                             '#94948C', 
                             '#734A8C', 
                             '#D1C0E1', 
                             '#648496', 
                             # '#B97A57', 
                             '#D1C0E1', 
                             # '#F8858A', 
                             '#F8858A', 
                             # '#63C6CE', 
                             '#94948C', 
                             # '#7BBD84', 
                             '#b6fcd5', 
                             '#E58835', 
                             # '#648496',
                             '#b6fcd5',
                             ],
                     hatch_patterns=('\\', '//', '|', 'x',),
                     filename=file_to_save+'IBO_TEA_breakdown_stacked_bar_plot',
                     n_minor_ticks=4,
                      # fig_height=5.5*1.1777*0.94*1.0975,
                     # fig_width=10,
                     fig_height = 6.5,
                     fig_width=16.5,
                     show_totals=True,
                     totals=totals,
                     sig_figs_for_totals=3,
                     units_list=[i.units for i in unit_groups[0].metrics],
                     totals_label_text=r"$\bfsum:$",
                     xticks_fontsize = 17,
                     ylabel_fontsize = 18,
                     yticks_fontsize = 17,
                     default_fontsize = 17,
                      ax=axs[1],
                     subplot_padding = 5,
                     bar_width=0.5,
                      xticks_fontcolor = 'black',
                     )


    #%% Spearman's rank order correlation coefficients - 1D
    from matplotlib import pyplot as plt
    chdir(IBO_results_filepath)
    plt.rcParams['font.sans-serif'] = "Arial Unicode"
    plt.rcParams['font.size'] = "7.5"
    
    file_to_save = IBO_results_filepath+\
        '_IBO_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
        + '_' + str(modes) + '_' + str(N_simulations_per_mode) + 'sims'
    
    bst_plots = bst.plots
    
    rho = r"$\mathrm{\rho}}$"
    
    mode = modes[0]
    
    
    
    fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['MPSP'][modes[0]],
                               index=[i.element_name + ': ' + i.name for i in model.parameters],
                               name='MPSP '+"["+MPSP_units+"]", color="#A97802",
                               # xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i,
                               )
    
    fig[0].set_figwidth(6)
    fig[0].set_figheight(10)
    
    fig[0].savefig(file_to_save+'_MPSP-Spearman.png', dpi=600, bbox_inches='tight',
                facecolor=fig[0].get_facecolor(),
                transparent=False)
    
    fig[0].show()

    #%% #%% Spearman's rank order correlation coefficients - 2D
    from matplotlib import pyplot as plt
    from matplotlib import colors
    
    chdir(IBO_results_filepath)
    plt.rcParams['font.sans-serif'] = "Arial Unicode"
    plt.rcParams['font.size'] = "7.5"
    
    file_to_save = IBO_results_filepath+\
        '_IBO_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
        + '_' + str(modes) + '_' + str(N_simulations_per_mode) + 'sims'
    
    bst_plots = bst.plots
    
    rho = r"$\mathrm{\rho}}$"
    
    mode = modes[0]
    
    # def sort_dfs_by_index(dfs, 
    #                       key,
    #                       ):
    #     for df in dfs:
    #         df.sort_index(key=key)
    
    # sort_dfs_by_index(dfs = [results_dict['Sensitivity']['Spearman']['MPSP'][mode],
    #                                   results_dict['Sensitivity']['Spearman']['GWP100a'][mode],
    #                                   results_dict['Sensitivity']['Spearman']['FEC'][mode],],
    #                   key = lambda x: print(x)
    #                   )
    
    fig = bst_plots.plot_spearman_2d([results_dict['Sensitivity']['Spearman']['MPSP'][mode],
                                      # results_dict['Sensitivity']['Spearman']['GWP100a'][mode],
                                      # results_dict['Sensitivity']['Spearman']['FEC'][mode],
                                      ],
                               index=[i.element_name + ': ' + i.name for i in model.parameters],
                               # name='MPSP '+"["+MPSP_units+"]", 
                               color_wheel=(colors.to_rgb("#A97802"),
                                            colors.to_rgb('#607429'),
                                            colors.to_rgb('#A100A1'),),
                               # xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i,
                               sort=False,
                               )
    
    fig[0].set_figwidth(2.5)
    fig[0].set_figheight(10)
    
    fig[0].savefig(file_to_save+'_Spearman-rho-2D.png', dpi=600, bbox_inches='tight',
                facecolor=fig[0].get_facecolor(),
                transparent=False)
    
    fig[0].show()
    
    #%% #%% Spearman's p-values - 2D
    from matplotlib import pyplot as plt
    from matplotlib import colors
    
    chdir(IBO_results_filepath)
    plt.rcParams['font.sans-serif'] = "Arial Unicode"
    plt.rcParams['font.size'] = "7.5"
    
    file_to_save = IBO_results_filepath+\
        '_IBO_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
        + '_' + str(modes) + '_' + str(N_simulations_per_mode) + 'sims'
    
    bst_plots = bst.plots
    
    rho = r"$\mathrm{\rho}}$"
    
    mode = modes[0]
    
    
    
    fig = bst_plots.plot_spearman_2d([results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode],
                                      results_dict['Sensitivity']['p-val Spearman']['GWP100a'][mode],
                                      results_dict['Sensitivity']['p-val Spearman']['FEC'][mode],],
                               index=[i.element_name + ': ' + i.name for i in model.parameters],
                               # name='MPSP '+"["+MPSP_units+"]", 
                               color_wheel=(colors.to_rgb("#A97802"),
                                            colors.to_rgb('#607429'),
                                            colors.to_rgb('#A100A1'),),
                               # xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i,
                               sort=False,
                               )
    fig[1].vlines(0.05, fig[1].get_yticks()[0], fig[1].get_yticks()[-1], 
                  color='gray', 
                  linewidth=1, 
                  linestyles='dashed')
    
    fig[0].set_figwidth(2.5)
    fig[0].set_figheight(10)
    
    fig[0].savefig(file_to_save+'_Spearman-p-vals-2D.png', dpi=600, bbox_inches='tight',
                facecolor=fig[0].get_facecolor(),
                transparent=False)
    
    fig[0].show()

#%% Sensitivity analysis take-aways

if len(modes)==1:
    file_to_save = IBO_results_filepath+\
        'significant_parameters_from_sensitivity_analysis_IBO_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
        + '_' + str(modes) + '_' + str(N_simulations_per_mode) + 'sims'
        
    spearman_results, spearman_p_values = model.spearman_r()
    
    cutoff_p_value = 0.05 # p-value < cutoff_p_value => significantly sensitive
    
    cutoff_rho_value = 0. # rho >= cutoff_rho_value => significantly sensitive
    round_rho_to = 2
    round_p_value_to = 2
    
    metrics = list(spearman_results.columns)
    
    sig_sens_parameters = {}
    sig_sens_parameters_dict_of_dicts = {}
    # sig_sens_parameters_dict_of_dicts = {"Spearman's rho":{}, "p-value":{}}
    
    for i in metrics:
        str_i_lower = str(i).lower()
        if 'total' in str_i_lower or 'adjusted minimum selling price' in str_i_lower:
            print(f"\n\nThe parameters to which the metric {i[1]} is most significantly sensitive (i.e., p-value < {cutoff_p_value} and Spearman's rho >= {cutoff_rho_value}) are as follows:")
            print("Parameter\t\t\t\t\t\t\t\tSpearman's rho\t\t\t\t\t\t\t\tp-value")
            sig_sens_parameters[i] = []
            sig_sens_parameters_dict_of_dicts[i[1]]={}
            sig_sens_parameters_dict_of_dicts[i[1]]["Spearman's rho"] = {}
            sig_sens_parameters_dict_of_dicts[i[1]]["p-value"] = {}
            for j in spearman_results[i].index:
                spearman_rho_val, spearman_p_val = spearman_results[i][j], spearman_p_values[i][j]
                if spearman_p_val < cutoff_p_value and abs(spearman_rho_val) > cutoff_rho_value:
                    sig_sens_parameters[i].append((j, spearman_rho_val, spearman_p_val))
                    
            sig_sens_parameters[i].sort(key=lambda m: abs(m[1]), reverse=True)
            for n in range(len(sig_sens_parameters[i])):
                (k, spearman_rho_val, spearman_p_val) = sig_sens_parameters[i][n]
                sig_sens_parameters_dict_of_dicts[i[1]]["Spearman's rho"][k[1]] = spearman_rho_val
                sig_sens_parameters_dict_of_dicts[i[1]]["p-value"][k[1]] = spearman_p_val
                # print(f"{n+1}. {k[1]}, with a Spearman's rho of {np.round(spearman_rho_val, round_rho_to)} (p-value = {np.round(spearman_p_val, round_p_value_to)})")
                print(f"\n{n+1}. {k[1]}\t\t\t\t\t\t\t\t{np.round(spearman_rho_val, round_rho_to)}\t\t\t\t\t\t\t\t{np.round(spearman_p_val, round_p_value_to)}")
    
    sheet_dfs = []
    writer = pd.ExcelWriter(file_to_save+'.xlsx')
    
    def remove_units(metric_str):
        return metric_str[:metric_str.index('[')]
    
    for i in list(sig_sens_parameters_dict_of_dicts.keys()):
        sheet_dfs.append(pd.DataFrame(sig_sens_parameters_dict_of_dicts[i]))
        # sheet_dfs[-1].sort_values(by="Spearman's rho", 
        #                           axis='index',
        #                            key=abs,
        #                            ascending=False,
        #                           )
        try:
            sheet_dfs[-1].to_excel(writer, remove_units(i))
        except:
            pass
    
    # writer.save()
    writer.close()
    
    

#%%

def get_uncertainty_data_from_saved_results(filenames, columns, 
                                            sheet_names='TEA results',
                                            skiprows=[0,2]):
    unc_data = []
    for n, s, c in zip(filenames, sheet_names, columns):
        unc_data.append(pd.read_excel(n, sheet_name=s, usecols=c, skiprows=skiprows))
    return unc_data