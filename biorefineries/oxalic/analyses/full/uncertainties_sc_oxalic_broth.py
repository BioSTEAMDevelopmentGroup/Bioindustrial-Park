#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Oxalic acid biorefineries.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040
"""

from warnings import filterwarnings
filterwarnings('ignore')
import numpy as np
import pandas as pd
import contourplots
import biosteam as bst
print('\n\nLoading system ...')
# from biorefineries
# from biorefineries import oxalic
from biorefineries import oxalic
from biorefineries.oxalic.models.sugarcane import models_sc_broth as models
from biorefineries.oxalic.process_settings import chem_index
# models = oxalic.models
# from . import models

print('\nLoaded system.')
from datetime import datetime
from biosteam.utils import TicToc
import os

from biorefineries.oxalic.analyses.full.plot_utils import plot_kde_formatted
from matplotlib.colors import hex2color

chdir = os.chdir
oxalic_filepath = oxalic.__file__.replace('\\__init__.py', '')
oxalic_results_filepath = oxalic_filepath + '\\analyses\\results\\'
model = models.oxalic_model
oxalic_chemicals = models.oxalic_chemicals

system = oxalic_sys = models.oxalic_sys
spec = models.spec
unit_groups = models.unit_groups

tea = models.oxalic_tea
lca = models.oxalic_lca
get_adjusted_MSP = models.get_adjusted_MSP
# per_kg_AA_to_per_kg_SA = models.per_kg_AA_to_per_kg_SA

# %% 

N_simulations_per_mode = 6000 # 6000

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]

notification_interval = 10

results_dict = {'Baseline':{'MPSP':{}, 'GWP100a':{}, 'FEC':{},
                               'GWP100a no electricity credit':{},
                               'FEC no electricity credit':{}, 
                            'GWP Breakdown':{}, 'FEC Breakdown':{},},
                'Uncertainty':{'MPSP':{}, 'GWP100a':{}, 'FEC':{},
                               'GWP100a no electricity credit':{},
                               'FEC no electricity credit':{},},
                'Sensitivity':{'Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}},
                               'p-val Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}}},}


feedstock_tag = 'sugarcane'
product_tag = 'oxalic-broth'

modes = ['100mL']

scenario_names = modes

#%%

timer = TicToc('timer')
timer.tic()

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221) # 3221

# def load_additional_params():
   
for i in range(len(modes)):
    # ## Change working directory to biorefineries\\oxalic
    # chdir(oxalic.__file__.replace('\\__init__.py', ''))
    # ##
    mode = modes[i]
    
    
    dist_filename = f'parameter-distributions_{feedstock_tag}_{product_tag}_' + mode + '.xlsx'
    
    product_folder = 'oxalic_broth_product' if product_tag=='oxalic-broth' else None
    
    parameter_distributions_filename = oxalic_filepath+\
        f'\\analyses\\full\\parameter_distributions\\{product_folder}\\'+dist_filename


    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename)
    
    # load_additional_params()
    print(f'\nLoaded parameter distributions ({mode}).')
    
    parameters = model.get_parameters()
    
    print('\n\nLoading samples ...')
    samples = model.sample(N=N_simulations_per_mode, rule='L')
    model.load_samples(samples)
    print('\nLoaded samples.')
    
    # ## Change working directory to biorefineries\\oxalic\\analyses\\results
    # chdir(oxalic.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
    # ##
    
    model.exception_hook = 'warn'
    print('\n\nSimulating baseline ...')
    
    baseline_initial = model.metrics_at_baseline()
    spec.set_production_capacity()
    baseline_initial = model.metrics_at_baseline()
    
    baseline = pd.DataFrame(data=np.array([[i for i in baseline_initial.values],]), 
                            columns=baseline_initial.keys())
    
    results_dict['Baseline']['MPSP'][mode] = get_adjusted_MSP()
    results_dict['Baseline']['GWP100a'][mode] = tot_GWP = lca.GWP
    results_dict['Baseline']['FEC'][mode] = tot_FEC = lca.FEC
    results_dict['Baseline']['GWP100a no electricity credit'][mode] = lca.GWP - lca.net_electricity_GWP
    results_dict['Baseline']['FEC no electricity credit'][mode] = lca.FEC - lca.net_electricity_FEC
    
    materials_to_include_in_impact_breakdowns = {
                            'CalciumDihydroxide': 'lime', 
                            'CSL': 'corn steep liquor',
                            # 'CH4': 'natural gas\n(for steam generation)',
                            # 'AmmoniumSulfate': 'ammonium sulfate',
                            } # materials shown as distinct contributions in breakdown plot;
                              # others, except feedstock, are lumped in 'other materials' for figure clarity
    
    material_GWP_breakdown = lca.material_GWP_breakdown

    results_dict['Baseline']['GWP Breakdown'][mode] = {
        'feedstock*': lca.FGHTP_GWP,
        }
    sum_GWP_included_materials = 0.
    for k,v in materials_to_include_in_impact_breakdowns.items():
        results_dict['Baseline']['GWP Breakdown'][mode][v] = mat_GWP = material_GWP_breakdown[k]
        sum_GWP_included_materials += mat_GWP
    
    results_dict['Baseline']['GWP Breakdown'][mode]['other materials'] = lca.material_GWP - sum_GWP_included_materials
    
    results_dict['Baseline']['GWP Breakdown'][mode].update(
        {
        'net electricity': lca.net_electricity_GWP,
        'direct non-biogenic\nemissions': lca.direct_non_biogenic_emissions_GWP,
        }
    )
    
    try:
        assert(round(tot_GWP,3)==round(sum([v for v in results_dict['Baseline']['GWP Breakdown'][mode].values()]), 3))
    except:
        breakpoint()
        
    tot_positive_GWP = sum([v for v in results_dict['Baseline']['GWP Breakdown'][mode].values() if v>0])
    for k, v in results_dict['Baseline']['GWP Breakdown'][mode].items():
        results_dict['Baseline']['GWP Breakdown'][mode][k] = v/tot_positive_GWP
      
    
    material_FEC_breakdown = lca.material_FEC_breakdown

    results_dict['Baseline']['FEC Breakdown'][mode] = {
        'feedstock': lca.feedstock_FEC,
        }
    sum_FEC_included_materials = 0.
    for k,v in materials_to_include_in_impact_breakdowns.items():
        results_dict['Baseline']['FEC Breakdown'][mode][v] = mat_FEC = material_FEC_breakdown[k]
        sum_FEC_included_materials += mat_FEC
    
    results_dict['Baseline']['FEC Breakdown'][mode]['other materials'] = lca.material_FEC - sum_FEC_included_materials
    
    results_dict['Baseline']['FEC Breakdown'][mode].update(
        {
        'net electricity': lca.net_electricity_FEC,
        }
    )
    
    try:
        assert(round(tot_FEC,3)==round(sum([v for v in results_dict['Baseline']['FEC Breakdown'][mode].values()]), 3))
    except:
        breakpoint()
        
    tot_positive_FEC = sum([v for v in results_dict['Baseline']['FEC Breakdown'][mode].values() if v>0])
    for k, v in results_dict['Baseline']['FEC Breakdown'][mode].items():
        results_dict['Baseline']['FEC Breakdown'][mode][k] = v/tot_positive_FEC
    
    # if spec.reactor.base_neutralizes_product: # sulfuric acid for acidulation
    #     results_dict['Baseline']['GWP Breakdown'][mode]['sulfuric acid'] = material_GWP_breakdown['H2SO4']
    #     results_dict['Baseline']['FEC Breakdown'][mode]['sulfuric acid'] = material_FEC_breakdown['H2SO4']
        
        
    print(f"\nSimulated baseline. MPSP = ${round(results_dict['Baseline']['MPSP'][mode],2)}/kg.")
    #%%
    print('\n\nEvaluating ...')
    model.evaluate(notify=notification_interval, autoload=None, autosave=None, file=None)
    print('\nFinished evaluation.')
    
    # Baseline results
    print('\n\nRe-simulating baseline ...')
    
    baseline_end = model.metrics_at_baseline()
    spec.set_production_capacity()
    baseline_end = model.metrics_at_baseline()
    
    print(f"\nRe-simulated baseline. MPSP = ${round(get_adjusted_MSP(),2)}/kg.")
    #%%
    dateTimeObj = datetime.now()
    minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
    file_to_save = oxalic_results_filepath\
        +'_' + product_tag + '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
        + '_' + str(modes) + '_' + feedstock_tag + '_' + str(N_simulations_per_mode) + 'sims'
    
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
    
    # LCA_results
    LCA_results = \
        model.table.iloc[:, index_TEA::].copy()
    LCA_percentiles = LCA_results.quantile(q=percentiles)
    
    # # Spearman's rank correlation
    
    table = model.table
    
    model.table = model.table.dropna()
    
    spearman_results, spearman_p_values = model.spearman_r()
    
    spearman_results.columns = pd.Index([i.name_with_units for i in model.metrics])
    spearman_p_values.columns = pd.Index([i.name_with_units for i in model.metrics])
    
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
    
    print('\n\nSaving raw results ...')
    with pd.ExcelWriter(file_to_save+'_'+mode+'_1_full_evaluation.xlsx') as writer:
        parameter_values.to_excel(writer, sheet_name='Parameters')
        TEA_results.to_excel(writer, sheet_name='TEA results')
        TEA_percentiles.to_excel(writer, sheet_name='TEA percentiles')
        LCA_results.to_excel(writer, sheet_name='LCA results')
        LCA_percentiles.to_excel(writer, sheet_name='LCA percentiles')
        spearman_results.to_excel(writer, sheet_name='Spearman')
        spearman_p_values.to_excel(writer, sheet_name='Spearman p-values')
        # one_p_df.to_excel(writer, sheet_name='One-parameter')
        model.table.to_excel(writer, sheet_name='Raw data')
    
    
    # results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
    results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price [$/kg AA]']
    results_dict['Uncertainty']['GWP100a'][mode] = model.table.Biorefinery['Total gwp100a [kg-CO2-eq/kg]'] # GWP or gwp
    results_dict['Uncertainty']['FEC'][mode] = model.table.Biorefinery['Total FEC [MJ/kg]']
    results_dict['Uncertainty']['GWP100a no electricity credit'][mode] = model.table.Biorefinery['Total gwp100a excl. net electricity [kg-CO2-eq/kg]'] # GWP or gwp
    results_dict['Uncertainty']['FEC no electricity credit'][mode] = model.table.Biorefinery['Total FEC excl. net electricity [MJ/kg]']
    
    df_rho, df_p = model.spearman_r()
    
    # results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
    results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price [$/kg AA]']
    results_dict['Sensitivity']['Spearman']['GWP100a'][mode] = df_rho['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
    results_dict['Sensitivity']['Spearman']['FEC'][mode] = df_rho['Biorefinery', 'Total FEC [MJ/kg]']
    
    # results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
    results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price [$/kg AA]']
    results_dict['Sensitivity']['p-val Spearman']['GWP100a'][mode] = df_p['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
    results_dict['Sensitivity']['p-val Spearman']['FEC'][mode] = df_p['Biorefinery', 'Total FEC [MJ/kg]']
    
#%% Clean up NaN values for plotting
metrics = ['MPSP', 
            'GWP100a', 
            'FEC',
            'GWP100a no electricity credit', 
            'FEC no electricity credit',
           ]
tot_NaN_vals_dict = results_dict['Errors'] = {metric: {mode: 0 for mode in modes} for metric in metrics}
for mode in modes:
    for metric in metrics:
        median_val = np.median(results_dict['Uncertainty'][metric][mode][~np.isnan(results_dict['Uncertainty'][metric][mode])])
        # median_val = 1.5
        for i in range(len(results_dict['Uncertainty'][metric][mode])):
            if np.isnan(results_dict['Uncertainty'][metric][mode][i]):
                results_dict['Uncertainty'][metric][mode][i] = median_val
                tot_NaN_vals_dict[metric][mode] += 1

# %% Plots
print('\n\nCreating and saving plots ...')

MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"
GWP_units = r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$"
FEC_units = r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$"
#%% Uncertainty
def get_small_range(num, offset):
    return(num-offset, num+offset)
baseline_marker_shapes=["D", "D", "s", "^", "s", "h", "h"]
baseline_marker_sizes=[6, 8, 6, 10]*2
baseline_marker_colors = ['w', '#F8858A']*4
#%% MPSP
# modes = ['A',]
file_to_save = oxalic_results_filepath+\
    '_' + feedstock_tag + '_' + product_tag + \
    '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
    + '_' + str(modes) + '_' + str(N_simulations_per_mode) + 'sims'


MPSP_uncertainty = [results_dict['Uncertainty']['MPSP'][mode]
                    for mode in modes
                    ]

# # search page for high end: https://www.alibaba.com/trade/search?spm=a2700.galleryofferlist.0.0.2a995827YzqZVg&fsb=y&IndexArea=product_en&assessmentCompany=true&keywords=590-00-1+sorbate&productTag=1200000228&ta=y&tab=all&
# SA_market_range=np.array([
#                           6.74, # 2019 global high end from Sorbic Acid Market, Transparency Market Research
#                           6.50 * 1.3397087, # $6.50/kg-potassium-sorbate from https://www.alibaba.com/product-detail/Lifecare-Supply-Potassium-Sorbate-High-Quality_1600897125355.html?spm=a2700.galleryofferlist.p_offer.d_title.1bc15827eAs1TL&s=p
#                           ]) 

# oxalic_maximum_viable_market_range = SA_market_range / theoretical_max_g_AA_per_g_SA
# lower end from https://businessanalytiq.com/procurementanalytics/index/oxalic-acid-price-index/#:~:text=North%20America:US$0.67/KG,:US$0.44/KG%2C%20unchanged;
# higher end from https://www.alibaba.com/product-detail/High-quality-Oxalic-acid-CAS-144_1600467879003.html?spm=a2700.7724857.0.0.737f3dddErELOM
market_range = np.array([
                          0.75, 3
                          ]) 



# biobased_lit_MPSP_range = (1.08, 3.63)

biobased_price = 2.688 * chem_index[2019] / chem_index[2015] # Taylor et al. 2015 report

contourplots.box_and_whiskers_plot(uncertainty_data=MPSP_uncertainty, 
                          baseline_values=[results_dict['Baseline']['MPSP'][mode] for mode in modes],
                          baseline_marker_shapes=baseline_marker_shapes,
                          baseline_marker_sizes=baseline_marker_sizes,
                          baseline_marker_colors=baseline_marker_colors,
                          baseline_locations=[i+1 for i in range(len(modes))],
                          boxcolor="#A97802",
                          ranges_for_comparison=[
                                                 market_range,
                                                 # [biobased_price*0.995, biobased_price*1.005],
                                                 ],
                          ranges_for_comparison_colors=[
                                                        '#c0c1c2', 
                                                        # '#646464',
                                                        # '#c0c1c2', 
                                                        
                                                        ],
                          # values_for_comparison=[biobased_price],
                          n_minor_ticks=4,
                          show_x_ticks=True,
                          x_tick_labels=scenario_names,
                          x_tick_wrap_width=14,
                          y_label=r"$\bfMSP$",
                          y_units=MPSP_units,
                          y_ticks=np.arange(0., 4.01, 0.5),
                          save_file=True,
                          fig_height=5.5,
                          fig_width = 2.5,
                          box_width=0.9,
                          filename=file_to_save+'_uncertainty_MPSP',
                          dpi=600,)


#%% LCA

#%% GWP100a

fossilbased_GWPs = [
                9.8411 + 2*oxalic_chemicals.CO2.MW/oxalic_chemicals.OxalicAcid.MW, # ecoinvent 3.8 (oxalic acid production, RoW) cradle-to-gate + EOL
                ]

GWP_uncertainty = [results_dict['Uncertainty']['GWP100a'][mode] 
                    for mode in modes
                    ]
GWP_uncertainty += [results_dict['Uncertainty']['GWP100a no electricity credit'][mode] 
                    for mode in modes
                    ]

contourplots.box_and_whiskers_plot(uncertainty_data=GWP_uncertainty, 
                           baseline_values=[results_dict['Baseline']['GWP100a'][mode] for mode in modes] + [results_dict['Baseline']['GWP100a no electricity credit'][mode] for mode in modes],
                          baseline_marker_sizes=baseline_marker_sizes,
                          baseline_marker_colors=baseline_marker_colors,
                          baseline_locations=[i+1 for i in range(len(modes))],
                          boxcolor='#607429',
                           # ranges_for_comparison=[get_small_range(i, 0.005) for i in fossilbased_GWPs],
                           # ranges_for_comparison_colors=['#c0c1c2' for i in range(len(fossilbased_GWPs))],
                            ranges_for_comparison=[[fossilbased_GWPs[0]*0.98, 
                                                   fossilbased_GWPs[0]*1.02]],
                            ranges_for_comparison_colors=['#c0c1c2'],
                          # values_for_comparison=fossilbased_GWPs,
                          n_minor_ticks=1,
                          show_x_ticks=True,
                          x_tick_labels=[
                                           'with electricity credit', 
                                           'without electricity credit',
                                         ],
                          x_tick_wrap_width=6,
                          # y_label=r"$\bfGWP-100a$",
                          y_label=r"$\mathrm{\bfCI}$",
                          y_units=GWP_units,
                          y_ticks=np.arange(-6., 12.01, 2),
                          save_file=True,
                          # fig_height=5.5,
                          # fig_width = 3.,
                          # box_width=0.65,
                          fig_height=5.5,
                          fig_width = 5,
                          box_width=0.4,
                          filename=file_to_save+'_uncertainty_GWP100a',
                          dpi=600,)

#%% FEC

biobased_FECs = [26, 27.7, 32.7]
# fossilbased_FECs = [59.2, 60.8, 112, 124]

fossilbased_FECs = [
                    49.013, # ecoinvent 3.8 (acrylic acid production, RoW)
                    116. # GREET 2023 (acrylic acid from fossil energy)
                    ]

FEC_uncertainty = [results_dict['Uncertainty']['FEC'][mode]
                    for mode in modes
                    # results_dict['Uncertainty']['FEC'][modes[1]],
                    # results_dict['Uncertainty']['FEC'][modes[2]],
                    ]


biobased_lit_FEC_values = [1, 2, 3] #!!!
contourplots.box_and_whiskers_plot(uncertainty_data=FEC_uncertainty, 
                          baseline_values=[results_dict['Baseline']['FEC'][mode] for mode in modes], 
                          baseline_marker_shapes=baseline_marker_shapes,
                          baseline_marker_sizes=baseline_marker_sizes,
                          baseline_marker_colors=baseline_marker_colors,
                          baseline_locations=[i+1 for i in range(len(modes))],
                          boxcolor='#A100A1',
                          # ranges_for_comparison=[fossilbased_FECs],
                          ranges_for_comparison_colors=['#c0c1c2'],
                          # ranges_for_comparison=[get_small_range(i, 0.061) for i in biobased_FECs+fossilbased_FECs],
                          # ranges_for_comparison_colors=['#c0c1c2' for i in range(len(biobased_FECs))] +\
                          #                               ['#646464' for i in range(len(fossilbased_FECs))],
                          # values_for_comparison=biobased_lit_FEC_values,
                          n_minor_ticks=4,
                          show_x_ticks=True,
                          x_tick_labels=scenario_names,
                          x_tick_wrap_width=6,
                          y_label=r"$\bfFEC$",
                          y_units=FEC_units,
                          y_ticks=np.arange(-0, 150.1, 25),
                          save_file=True,
                          # fig_height=5.5,
                          # fig_width = 3.,
                          # box_width=0.65,
                          fig_height=5.5,
                          fig_width = 6,
                          box_width=0.3,
                          filename=file_to_save+'_uncertainty_FEC',
                          dpi=600,)

#%% TEA breakdown figure



###### change operating cost unit labels $/h to MM$/y
for i in unit_groups:
    ### change upgrading name
    if i.name=='upgrading': i.name = 'catalytic upgrading'

    for j in i.metrics:
        if j.name == 'Operating cost':
            j.units = r"$\mathrm{MM\$}$" + '\u00b7y\u207b\u00b9'
    # change natural gas names for spacing
    if i.name=='natural gas (for product drying)':
        i.name='natural gas\n(for product drying)'
    if i.name=='natural gas (for steam generation)':
        i.name='natural gas\n(for steam generation)'
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
                 y_ticks = [-50, -25,  0, 25, 50, 75, 100],
                 y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
                 y_units = "%", 
                 colors=['#7BBD84', 
                         '#E58835', 
                         '#F7C652', 
                         '#63C6CE', 
                         '#F8858A', 
                         '#94948C', 
                         '#734A8C', 
                         '#D1C0E1', 
                         '#648496', 
                         # '#B97A57', 
                         '#D1C0E1', 
                         # '#F8858A', 
                           '#b00000', 
                         # '#63C6CE', 
                         '#94948C', 
                         # '#7BBD84', 
                         '#b6fcd5', 
                         '#E58835', 
                         # '#648496',
                         '#b6fcd5',
                         ],
                 hatch_patterns=('\\', '//', '|', 'x',),
                 filename=file_to_save+'_'+mode+'_'+'AA_TEA_breakdown_stacked_bar_plot',
                 n_minor_ticks=4,
                 fig_height=5.5*1.1777*0.94*1.0975,
                 fig_width=10,
                 show_totals=True,
                 totals=totals,
                 sig_figs_for_totals=3,
                 units_list=[i.units for i in unit_groups[0].metrics],
                 totals_label_text=r"$\bfsum:$",
                 rotate_xticks=45.,
                 )


#%% LCA breakdown figures

mode = modes[0]

# GWP
temp_GWP_breakdown_dict = results_dict['Baseline']['GWP Breakdown'][mode]
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

df_GWP_breakdown = df_GWP_breakdown.rename(columns={0: ''})
contourplots.stacked_bar_plot(dataframe=df_GWP_breakdown, 
                  # y_ticks=[-60, -40, -20, 0, 20, 40, 60, 80, 100], 
                  y_ticks = [-200, -150, -100, -50, 0, 50, 100],
                  # y_ticks=[-400, -300, -200, -100, 0, 100, 200, 300, 400], 
                  # y_ticks = []
                  # y_label=r"$\bfGWP-100a $" +" "+ r"$\bfBreakdown$",  
                  y_label =r"$\mathrm{\bfCarbon}$" + " " + r"$\mathrm{\bfIntensity}$" +" "+ r"$\bfBreakdown$",
                  y_units = "%", 
                  colors=['#E1F8C0', 
                   '#8FAE3E', 
                          '#607429', 
                          ],
                  hatch_patterns=('\\', '//', 'x', '|',  ),
                  # '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#F8858A', 'red', 'magenta'],
                  filename=file_to_save+'_'+mode+'_''AA_GWP_breakdown_stacked_bar_plot',
                  fig_width=2,
                  fig_height=5.5*1.1777*0.94)

# FEC
temp_FEC_breakdown_dict = results_dict['Baseline']['FEC Breakdown'][mode]
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
df_FEC_breakdown = df_FEC_breakdown.rename(columns={0: ''})
contourplots.stacked_bar_plot(dataframe=df_FEC_breakdown, 
                  # y_ticks=[-200, -175, -150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175], 
                   y_ticks=[0, 20, 40, 60, 80, 100], 
                  # y_ticks=[-150, -100, -50, 0, 50, 100], 
                  y_label=r"$\bfFEC$" +" "+ r"$\bfBreakdown$", 
                  y_units = "%", 
                  # colors=['#7BBD84', '#F7C652', '#63C6CE', '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#F8858A', 'magenta'],
                  colors=['#FEC1FE', '#FF80FF', '#A100A1', 
                          ],
                  hatch_patterns=('\\', '//', '|', 'x',),
                  filename=file_to_save+'_'+mode+'_'+'AA_FEC_breakdown_stacked_bar_plot',
                  fig_width=2,
                  fig_height=5.5*1.1777*0.94)

#%% Spearman's rank order correlation coefficients - 1D
from matplotlib import pyplot as plt
chdir(oxalic_results_filepath)
plt.rcParams['font.sans-serif'] = "Arial Unicode"
plt.rcParams['font.size'] = "7.5"

file_to_save = oxalic_results_filepath+\
    '_AA_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
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

fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['GWP100a'][modes[0]],
                            index=[i.element_name + ': ' + i.name for i in model.parameters],
                            name='GWP100'+"["+GWP_units+"]", color='#607429',
                            # xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i,
                            )
fig[0].set_figwidth(6)
fig[0].set_figheight(10)

fig[0].savefig(file_to_save+'_GWP-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)
fig = bst_plots.plot_spearman_1d(results_dict['Sensitivity']['Spearman']['FEC'][modes[0]],
                            index=[i.element_name + ': ' + i.name for i in model.parameters],
                            name='FEC'+"["+FEC_units+"]", color='#A100A1',
                            # xlabel_fn=lambda i: "Spearman's "+rho+ " with "+i,
                            )
fig[0].set_figwidth(6)
fig[0].set_figheight(10)

fig[0].savefig(file_to_save+'_FEC-Spearman.png', dpi=600, bbox_inches='tight',
            facecolor=fig[0].get_facecolor(),
            transparent=False)

fig[0].show()



#%% #%% Spearman's rank order correlation coefficients - 2D
from matplotlib import pyplot as plt
from matplotlib import colors

chdir(oxalic_results_filepath)
plt.rcParams['font.sans-serif'] = "Arial Unicode"
plt.rcParams['font.size'] = "7.5"

file_to_save = oxalic_results_filepath+\
    '_AA_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
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

fig = bst_plots.plot_spearman_2d([results_dict['Sensitivity']['Spearman']['MPSP'][mode][1:],
                                  results_dict['Sensitivity']['Spearman']['GWP100a'][mode][1:],
                                  results_dict['Sensitivity']['Spearman']['FEC'][mode][1:],],
                           index=[i.element_name + ': ' + i.name for i in model.parameters[1:]],
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

chdir(oxalic_results_filepath)
plt.rcParams['font.sans-serif'] = "Arial Unicode"
plt.rcParams['font.size'] = "7.5"

file_to_save = oxalic_results_filepath+\
    '_AA_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
    + '_' + str(modes) + '_' + str(N_simulations_per_mode) + 'sims'

bst_plots = bst.plots

rho = r"$\mathrm{\rho}}$"

mode = modes[0]



fig = bst_plots.plot_spearman_2d([results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode][1:],
                                  results_dict['Sensitivity']['p-val Spearman']['GWP100a'][mode][1:],
                                  results_dict['Sensitivity']['p-val Spearman']['FEC'][mode][1:],],
                           index=[i.element_name + ': ' + i.name for i in model.parameters[1:]],
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
#!!! TODO: make compatible with multiple modes

file_to_save = oxalic_results_filepath+\
    'significant_parameters_from_sensitivity_analysis_AA_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
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
    sheet_dfs[-1].to_excel(writer, remove_units(i))

# writer.save()
writer.close()

#%% Bivariate kernel density plots
file_to_save = oxalic_results_filepath+\
    '_AA_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
    + '_' + str(modes) + '_' + str(N_simulations_per_mode) + 'sims'
    
#%% MPSP-GWP
plot_kde_formatted(
                    xdata = np.array(list(MPSP_uncertainty[0])),
                    ydata=np.array(list(GWP_uncertainty[0])),
                    
                    xlabel = r"$\bfMPSP$" + " [" + r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$" + "]",
                    
                    ylabel = r"$\mathrm{\bfGWP}_{\bf100}$" + " [" + r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$" + "]",
                    
                    
                    xticks = [i for i in range(0,9) if not i%1],
                    yticks = [i for i in range(-8,6) if not i%2],
                    
                    n_minor_ticks = 3,
                    
                    fig_width = 5,
                    fig_height = 5,
                    
                    show_x_ticks = True,
                    xbox_kwargs=dict(light=hex2color("#A97802"), dark=(0,0,0)),
                    ybox_kwargs=dict(light=hex2color('#607429'), dark=(0,0,0)),
                    
                    save_fig = True,
                    filename = file_to_save+'_Bivariate-KDE_MPSP-GWP'+'.png',
                    
                    )

#%% MPSP-FEC
plot_kde_formatted(
                    xdata = np.array(list(MPSP_uncertainty[0])),
                    ydata=np.array(list(FEC_uncertainty[0])),
                    
                    xlabel = r"$\bfMPSP$" + " [" + r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$" + "]",
                    
                    ylabel = r"$\bfFEC$" + " [" + r"$\mathrm{kg}$"+" "+ r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$" + "]",
                    
                    
                    xticks = [i for i in range(0,9) if not i%1],
                    yticks = [i for i in range(-50,51) if not i%25],
                    
                    n_minor_ticks = 3,
                    
                    fig_width = 5,
                    fig_height = 5,
                    
                    show_x_ticks = True,
                    xbox_kwargs=dict(light=hex2color("#A97802"), dark=(0,0,0)),
                    ybox_kwargs=dict(light=hex2color('#A100A1'), dark=(0,0,0)),
                    
                    save_fig = True,
                    filename = file_to_save+'_Bivariate-KDE_MPSP-FEC'+'.png',
                    
                    )

#%% GWP-FEC
plot_kde_formatted(
                    xdata = np.array(list(GWP_uncertainty[0])),
                    ydata=np.array(list(FEC_uncertainty[0])),
                    
                    xlabel = r"$\mathrm{\bfGWP}_{\bf100}$" + " [" + r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$" + "]",
                    
                    ylabel = r"$\bfFEC$" + " [" + r"$\mathrm{kg}$"+" "+ r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$" + "]",
                    
                    
                    xticks = [i for i in range(0,15) if not i%2],
                    yticks = [i for i in range(-50,51) if not i%25],
                    
                    n_minor_ticks = 3,
                    
                    fig_width = 5,
                    fig_height = 5,
                    
                    show_x_ticks = True,
                    xbox_kwargs=dict(light=hex2color('#607429'), dark=(0,0,0)),
                    ybox_kwargs=dict(light=hex2color('#A100A1'), dark=(0,0,0)),
                    
                    save_fig = True,
                    filename = file_to_save+'_Bivariate-KDE_GWP-FEC'+'.png',
                    
                    )

#%% Bivariate KDE with seaborn

# df_MPSP_GWP_results = pd.DataFrame({'MPSP': MPSP_uncertainty[0], 'GWP100a': GWP_uncertainty[0]})
# # ax = seaborn.kdeplot(
# #     data=df_MPSP_GWP_results, 
# #     x='MPSP',
# #     y='GWP100a',
# #     fill=True,
# #     cbar=True, 
# #     levels = 10,
# #     # levels=[0.05*i for i in range(1,21)],
# #     thresh=0.01,
# #     )


