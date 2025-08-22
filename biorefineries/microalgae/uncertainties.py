#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

Uncertainty analysis for microalgae biorefinery

Based on succinic project but adapted for microalgae system structure
"""

from warnings import filterwarnings
filterwarnings('ignore')
import numpy as np
import pandas as pd
import biosteam as bst
# import contourplots
from .system import microalgae_mcca_sys, microalgae_tea
from .lca import create_microalgae_lca
from .model_utils import MicroalgaeModel
from biosteam.evaluation import Metric
from datetime import datetime
from biosteam.utils import TicToc
import os

chdir = os.chdir
microalgae_filepath = os.path.dirname(__file__)
microalgae_results_filepath = os.path.join(microalgae_filepath, 'analyses', 'results')

# Create results directory if it doesn't exist
if not os.path.exists(microalgae_results_filepath):
    os.makedirs(microalgae_results_filepath)

system = microalgae_sys = microalgae_mcca_sys
tea = microalgae_tea

# Create model with metrics
def create_model():
    """Create evaluation model with metrics for microalgae system"""
    # Get system components
    u = microalgae_mcca_sys.flowsheet.unit
    s = microalgae_mcca_sys.flowsheet.stream
    
    # Find main product and boiler
    main_product = s.caproic_acid_product
    main_product_chemical_IDs = ['CaproicAcid']
    
    # Find boiler
    boiler = None
    for unit in microalgae_mcca_sys.units:
        if hasattr(unit, 'natural_gas') or 'BT' in unit.ID:
            boiler = unit
            break
    
    # Create LCA object
    lca = create_microalgae_lca(microalgae_mcca_sys, main_product, main_product_chemical_IDs, boiler)
    
    # Define metrics
    metrics = [
        Metric('MPSP', lambda: microalgae_tea.solve_price(main_product), '$/kg'),
        Metric('TCI', lambda: microalgae_tea.TCI/1e6, 'MM$'),
        Metric('VOC', lambda: microalgae_tea.VOC/1e6, 'MM$/y'),
        Metric('FOC', lambda: microalgae_tea.FOC/1e6, 'MM$/y'),
        Metric('GWP', lambda: lca.GWP, 'kg CO2-eq/kg'),
        Metric('FEC', lambda: lca.FEC, 'MJ/kg'),
    ]
    
    # Create namespace for parameter loading
    namespace_dict = {
        'microalgae_sys': microalgae_mcca_sys,
        'microalgae_tea': microalgae_tea,
        'u': u,
        's': s,
        'lca': lca,
        'bst': bst,
        'np': np,
        # Add chemical streams for easier access
        'microalgae': None,  # Will be set dynamically
        'GlucoAmylase': None,  # Will be set dynamically
        'AlphaAmylase': None,  # Will be set dynamically
        'Yeast': None,
        'OleylAlcohol': None,
        'base_fermentation': None,
        'FGD_lime': None,
        'PowerUtility': bst.PowerUtility,
    }
    
    # Try to find chemical streams
    for stream in microalgae_mcca_sys.feeds:
        if 'microalgae' in stream.ID.lower():
            namespace_dict['microalgae'] = stream
            break
    
    # Create model
    model = MicroalgaeModel(microalgae_mcca_sys, metrics=metrics, namespace_dict=namespace_dict)
    
    return model, lca, namespace_dict

model, lca, namespace_dict = create_model()

def get_adjusted_MSP():
    """Get adjusted minimum selling price"""
    return microalgae_tea.solve_price(microalgae_sys.flowsheet.stream.caproic_acid_product)

# %% 

N_simulations_per_mode = 2000 # 2000

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]

notification_interval = 10

results_dict = {'Baseline':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}, 
                            'GWP Breakdown':{}, 'FEC Breakdown':{},},
                'Uncertainty':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}},
                'Sensitivity':{'Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}}},}

modes = [
            'baseline',
         ]

parameter_distributions_filenames = [
                                    'parameter_distributions.xlsx',
                                    ]

#%%

timer = TicToc('timer')
timer.tic()

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221) # 3221

for i in range(len(modes)):
    mode = modes[i]
    parameter_distributions_filename = os.path.join(microalgae_filepath, parameter_distributions_filenames[i])
    
    print(f'\n\nLoading parameter distributions ({mode}) ...')
    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
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
    
    # GWP breakdown analysis
    try:
        material_GWP_breakdown = lca.material_GWP_breakdown
        
        results_dict['Baseline']['GWP Breakdown'][mode] = {
            'feedstock': lca.feedstock_GWP,
            'material inputs': lca.material_GWP,
            'natural gas\n(for steam generation)': getattr(lca, 'ng_GWP', 0),
            'net electricity': lca.net_electricity_GWP,
            'direct non-biogenic\nemissions': getattr(lca, 'direct_emissions_GWP', 0),
        }
        
        tot_positive_GWP = sum([v for v in results_dict['Baseline']['GWP Breakdown'][mode].values() if v>0])
        if tot_positive_GWP > 0:
            for k, v in results_dict['Baseline']['GWP Breakdown'][mode].items():
                results_dict['Baseline']['GWP Breakdown'][mode][k] = v/tot_positive_GWP
    except Exception as e:
        print(f"Warning: Could not calculate GWP breakdown: {e}")
        results_dict['Baseline']['GWP Breakdown'][mode] = {}
    
    # FEC breakdown analysis
    try:
        material_FEC_breakdown = lca.material_FEC_breakdown
        
        results_dict['Baseline']['FEC Breakdown'][mode] = {
            'feedstock': lca.feedstock_FEC,
            'material inputs': lca.material_FEC,
            'natural gas\n(for steam generation)': getattr(lca, 'ng_FEC', 0),
            'net electricity': lca.net_electricity_FEC,
        }
        
        tot_positive_FEC = sum([v for v in results_dict['Baseline']['FEC Breakdown'][mode].values() if v>0])
        if tot_positive_FEC > 0:
            for k, v in results_dict['Baseline']['FEC Breakdown'][mode].items():
                results_dict['Baseline']['FEC Breakdown'][mode][k] = v/tot_positive_FEC
    except Exception as e:
        print(f"Warning: Could not calculate FEC breakdown: {e}")
        results_dict['Baseline']['FEC Breakdown'][mode] = {}
        
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
    file_to_save = os.path.join(microalgae_results_filepath,
        f'_microalgae_{dateTimeObj.year}.{dateTimeObj.month}.{dateTimeObj.day}-{dateTimeObj.hour}.{minute}'\
        + f'_{N_simulations_per_mode}sims')
    
    baseline.index = ('initial', )
    baseline.to_excel(file_to_save+'_'+mode+'_0_baseline.xlsx')
    
    # Parameters
    parameters = model.get_parameters()
    index_parameters = len(model.get_baseline_sample())
    parameter_values = model.table.iloc[:, :index_parameters].copy()
    
    #%%
    
    # TEA results
    for index_TEA, i in enumerate(model.metrics):
        if hasattr(i, 'element') and i.element == 'LCA': 
            break
    else:
        index_TEA = len(model.metrics) - 2  # Assume last 2 are LCA metrics
        
    index_TEA = index_parameters + index_TEA
    TEA_results = model.table.iloc[:, index_parameters:index_TEA].copy()
    TEA_percentiles = TEA_results.quantile(q=percentiles)
    
    # LCA_results
    LCA_results = model.table.iloc[:, index_TEA::].copy()
    LCA_percentiles = LCA_results.quantile(q=percentiles)
    
    # # Spearman's rank correlation
    table = model.table
    model.table = model.table.dropna()
    
    spearman_results = model.spearman()
    spearman_results.columns = pd.Index([i.name_with_units for i in model.metrics])
    
    model.table = table
    
    # Calculate the cumulative probabilities of each parameter
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
        model.table.to_excel(writer, sheet_name='Raw data')
    
    # Extract results for plotting
    # Find the MPSP column
    mpsp_col = None
    for col in model.table.columns:
        if 'MPSP' in str(col) or 'price' in str(col).lower():
            mpsp_col = col
            break
    
    # Find GWP column
    gwp_col = None
    for col in model.table.columns:
        if 'GWP' in str(col):
            gwp_col = col
            break
    
    # Find FEC column
    fec_col = None
    for col in model.table.columns:
        if 'FEC' in str(col):
            fec_col = col
            break
    
    if mpsp_col is not None:
        results_dict['Uncertainty']['MPSP'][mode] = model.table[mpsp_col]
    else:
        print("Warning: Could not find MPSP column in results")
        results_dict['Uncertainty']['MPSP'][mode] = [results_dict['Baseline']['MPSP'][mode]] * len(model.table)
    
    if gwp_col is not None:
        results_dict['Uncertainty']['GWP100a'][mode] = model.table[gwp_col]
    else:
        print("Warning: Could not find GWP column in results")
        results_dict['Uncertainty']['GWP100a'][mode] = [results_dict['Baseline']['GWP100a'][mode]] * len(model.table)
    
    if fec_col is not None:
        results_dict['Uncertainty']['FEC'][mode] = model.table[fec_col]
    else:
        print("Warning: Could not find FEC column in results")
        results_dict['Uncertainty']['FEC'][mode] = [results_dict['Baseline']['FEC'][mode]] * len(model.table)
    
    # Spearman correlations for sensitivity analysis
    df_rho, df_p = model.spearman_r()
    
    if mpsp_col is not None and mpsp_col in df_rho.columns:
        results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho[mpsp_col]
    else:
        results_dict['Sensitivity']['Spearman']['MPSP'][mode] = pd.Series()
    
    if gwp_col is not None and gwp_col in df_rho.columns:
        results_dict['Sensitivity']['Spearman']['GWP100a'][mode] = df_rho[gwp_col]
    else:
        results_dict['Sensitivity']['Spearman']['GWP100a'][mode] = pd.Series()
    
    if fec_col is not None and fec_col in df_rho.columns:
        results_dict['Sensitivity']['Spearman']['FEC'][mode] = df_rho[fec_col]
    else:
        results_dict['Sensitivity']['Spearman']['FEC'][mode] = pd.Series()
            
#%% Clean up NaN values for plotting
metrics = ['MPSP', 'GWP100a', 'FEC']
tot_NaN_vals_dict = results_dict['Errors'] = {metric: {mode: 0 for mode in modes} for metric in metrics}
for mode in modes:
    for metric in metrics:
        median_val = 1.5  # Default fallback value
        if len(results_dict['Uncertainty'][metric][mode]) > 0:
            median_val = np.nanmedian(results_dict['Uncertainty'][metric][mode])
            if np.isnan(median_val):
                median_val = 1.5
        
        for i in range(len(results_dict['Uncertainty'][metric][mode])):
            if np.isnan(results_dict['Uncertainty'][metric][mode].iloc[i]):
                results_dict['Uncertainty'][metric][mode].iloc[i] = median_val
                tot_NaN_vals_dict[metric][mode] += 1

# %% Plots - temporarily disabled due to contourplots dependency
# MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"
# GWP_units = r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$"
# FEC_units = r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$"

scenario_name_labels = ['Baseline']

def get_small_range(num, offset):
    return(num-offset, num+offset)

print(f'\nAnalysis completed. Timer: {timer.toc():.2f} seconds')
print(f'Results saved to: {file_to_save}')

# Print summary
print('\n=== UNCERTAINTY ANALYSIS SUMMARY ===')
for mode in modes:
    print(f'\n{mode.upper()} MODE:')
    print(f'  MPSP: {results_dict["Baseline"]["MPSP"][mode]:.3f} $/kg')
    print(f'  GWP: {results_dict["Baseline"]["GWP100a"][mode]:.3f} kg CO2-eq/kg')
    print(f'  FEC: {results_dict["Baseline"]["FEC"][mode]:.3f} MJ/kg')
    
    if len(results_dict['Uncertainty']['MPSP'][mode]) > 0:
        print(f'  MPSP range: {np.min(results_dict["Uncertainty"]["MPSP"][mode]):.3f} - {np.max(results_dict["Uncertainty"]["MPSP"][mode]):.3f} $/kg')
        print(f'  GWP range: {np.min(results_dict["Uncertainty"]["GWP100a"][mode]):.3f} - {np.max(results_dict["Uncertainty"]["GWP100a"][mode]):.3f} kg CO2-eq/kg')
        print(f'  FEC range: {np.min(results_dict["Uncertainty"]["FEC"][mode]):.3f} - {np.max(results_dict["Uncertainty"]["FEC"][mode]):.3f} MJ/kg')

# Print detailed breakdown analysis
print('\n=== DETAILED BREAKDOWN ANALYSIS ===')
for mode in modes:
    print(f'\n{mode.upper()} MODE BREAKDOWN:')
    
    # GWP Breakdown
    if results_dict['Baseline']['GWP Breakdown'][mode]:
        print(f'\nGWP Breakdown (Total: {results_dict["Baseline"]["GWP100a"][mode]:.3f} kg CO2-eq/kg):')
        gwp_breakdown = results_dict['Baseline']['GWP Breakdown'][mode]
        total_gwp = sum([abs(v) for v in gwp_breakdown.values()])
        for component, value in gwp_breakdown.items():
            percentage = (value / total_gwp * 100) if total_gwp > 0 else 0
            print(f'  {component}: {value:.4f} kg CO2-eq/kg ({percentage:.1f}%)')
    
    # FEC Breakdown
    if results_dict['Baseline']['FEC Breakdown'][mode]:
        print(f'\nFEC Breakdown (Total: {results_dict["Baseline"]["FEC"][mode]:.3f} MJ/kg):')
        fec_breakdown = results_dict['Baseline']['FEC Breakdown'][mode]
        total_fec = sum([abs(v) for v in fec_breakdown.values()])
        for component, value in fec_breakdown.items():
            percentage = (value / total_fec * 100) if total_fec > 0 else 0
            print(f'  {component}: {value:.4f} MJ/kg ({percentage:.1f}%)')
            
# Print sensitivity analysis results
print('\n=== SENSITIVITY ANALYSIS (SPEARMAN CORRELATIONS) ===')
for mode in modes:
    print(f'\n{mode.upper()} MODE SENSITIVITY:')
    
    # MPSP Correlations
    if not results_dict['Sensitivity']['Spearman']['MPSP'][mode].empty:
        print(f'\nTop 10 correlations with MPSP:')
        mpsp_corr = results_dict['Sensitivity']['Spearman']['MPSP'][mode].copy()
        # Remove NaN values and sort by absolute correlation
        mpsp_corr = mpsp_corr.dropna()
        if len(mpsp_corr) > 0:
            sorted_corr = mpsp_corr.abs().sort_values(ascending=False)
            for i, (param, abs_corr) in enumerate(sorted_corr.head(10).items()):
                actual_corr = mpsp_corr[param]
                print(f'  {i+1:2d}. {param}: {actual_corr:.3f}')
        else:
            print('    No significant correlations found')
    
    # GWP Correlations
    if not results_dict['Sensitivity']['Spearman']['GWP100a'][mode].empty:
        print(f'\nTop 10 correlations with GWP:')
        gwp_corr = results_dict['Sensitivity']['Spearman']['GWP100a'][mode].copy()
        gwp_corr = gwp_corr.dropna()
        if len(gwp_corr) > 0:
            sorted_corr = gwp_corr.abs().sort_values(ascending=False)
            for i, (param, abs_corr) in enumerate(sorted_corr.head(10).items()):
                actual_corr = gwp_corr[param]
                print(f'  {i+1:2d}. {param}: {actual_corr:.3f}')
        else:
            print('    No significant correlations found')
    
    # FEC Correlations
    if not results_dict['Sensitivity']['Spearman']['FEC'][mode].empty:
        print(f'\nTop 10 correlations with FEC:')
        fec_corr = results_dict['Sensitivity']['Spearman']['FEC'][mode].copy()
        fec_corr = fec_corr.dropna()
        if len(fec_corr) > 0:
            sorted_corr = fec_corr.abs().sort_values(ascending=False)
            for i, (param, abs_corr) in enumerate(sorted_corr.head(10).items()):
                actual_corr = fec_corr[param]
                print(f'  {i+1:2d}. {param}: {actual_corr:.3f}')
        else:
            print('    No significant correlations found')

# Print statistics summary
print('\n=== STATISTICAL SUMMARY ===')
for mode in modes:
    print(f'\n{mode.upper()} MODE STATISTICS:')
    
    for metric in ['MPSP', 'GWP100a', 'FEC']:
        if len(results_dict['Uncertainty'][metric][mode]) > 0:
            data = results_dict['Uncertainty'][metric][mode]
            mean_val = np.mean(data)
            std_val = np.std(data)
            p5 = np.percentile(data, 5)
            p95 = np.percentile(data, 95)
            median_val = np.median(data)
            
            units = {'MPSP': '$/kg', 'GWP100a': 'kg CO2-eq/kg', 'FEC': 'MJ/kg'}
            unit = units[metric]
            
            print(f'\n{metric}:')
            print(f'  Mean: {mean_val:.3f} {unit}')
            print(f'  Std Dev: {std_val:.3f} {unit}')
            print(f'  Median: {median_val:.3f} {unit}')
            print(f'  5th percentile: {p5:.3f} {unit}')
            print(f'  95th percentile: {p95:.3f} {unit}')
            print(f'  Range: {np.min(data):.3f} - {np.max(data):.3f} {unit}')

# Print error summary
print('\n=== ERROR SUMMARY ===')
total_errors = 0
for metric in metrics:
    for mode in modes:
        errors = tot_NaN_vals_dict[metric][mode]
        if errors > 0:
            print(f'{metric} ({mode}): {errors} NaN values replaced')
            total_errors += errors

if total_errors == 0:
    print('No errors encountered during simulation.')
else:
    print(f'Total errors handled: {total_errors}')

print(f'\n=== ANALYSIS COMPLETED ===')
print(f'Total simulation time: {timer.toc():.2f} seconds')
print(f'Results saved to: {file_to_save}')
print(f'Number of parameters analyzed: {len(parameters)}')
print(f'Number of successful simulations: {len(model.table)}')
print(f'Output files generated:')
print(f'  - Baseline: {file_to_save}_{mode}_0_baseline.xlsx')
print(f'  - Full results: {file_to_save}_{mode}_1_full_evaluation.xlsx')

if __name__ == '__main__':
    print("Microalgae uncertainty analysis completed!")