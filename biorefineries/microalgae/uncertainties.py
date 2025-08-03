#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microalgae biorefinery uncertainty and sensitivity analysis

This module performs uncertainty and sensitivity analysis for the microalgae biorefinery
system, including Monte Carlo simulations and Spearman correlation analysis.

@author: Xingdong Shi
@version: 0.0.1
"""

from warnings import filterwarnings
filterwarnings('ignore')
import numpy as np
import pandas as pd
import biosteam as bst
import thermosteam as tmo
from datetime import datetime
from biosteam.utils import TicToc
import os

# Import microalgae system components
from . import system as microalgae_system
from . import lca
from . import tea
from ._chemicals import chems

# Set up the system
bst.settings.set_thermo(chems)
microalgae_sys = microalgae_system.microalgae_mcca_sys
microalgae_tea = microalgae_system.microalgae_tea

# Create LCA object
from . import analysis_utils
microalgae_lca = analysis_utils.create_microalgae_lca_simple(microalgae_sys, microalgae_tea)

# Create model
# Get main product for TEA calculations
main_product = analysis_utils.get_main_product_stream(microalgae_sys)

model = bst.Model(microalgae_sys, metrics=[
    bst.Metric('MPSP', lambda: microalgae_tea.solve_price(main_product), 'USD/kg'),
    bst.Metric('GWP100a', lambda: microalgae_lca.GWP, 'kg CO2-eq/kg'),
    bst.Metric('FEC', lambda: microalgae_lca.FEC, 'MJ/kg'),
])

print('\n\nLoaded microalgae system for uncertainty analysis.')

# Analysis parameters
N_simulations_per_mode = 2000
percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]
notification_interval = 10

results_dict = {
    'Baseline': {
        'MPSP': {}, 
        'GWP100a': {}, 
        'FEC': {}, 
        'GWP Breakdown': {}, 
        'FEC Breakdown': {}
    },
    'Uncertainty': {
        'MPSP': {}, 
        'GWP100a': {}, 
        'FEC': {}
    },
    'Sensitivity': {
        'Spearman': {
            'MPSP': {}, 
            'GWP100a': {}, 
            'FEC': {}
        }
    }
}

# Parameter distribution file
parameter_distributions_filename = 'parameter_distributions.xlsx'

def run_uncertainty_analysis():
    """Run uncertainty and sensitivity analysis for microalgae biorefinery."""
    
    timer = TicToc('timer')
    timer.tic()
    
    # Set seed for reproducibility
    np.random.seed(3221)
    
    print(f'\n\nLoading parameter distributions...')
    
    # Get current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    param_file_path = os.path.join(current_dir, parameter_distributions_filename)
    
    # Load parameter distributions
    model.parameters = ()
    
    # Create namespace dictionary for microalgae system
    namespace_dict = {
        'microalgae_sys': microalgae_sys,
        'microalgae_tea': microalgae_tea,
        'microalgae_lca': microalgae_lca,
    }
    
    # Add all units to namespace
    for unit in microalgae_sys.units:
        namespace_dict[unit.ID] = unit
    
    # Add all streams to namespace  
    for stream in microalgae_sys.streams:
        namespace_dict[stream.ID] = stream
        
    # Add chemicals to namespace
    for chemical in chems:
        namespace_dict[chemical.ID] = chemical
    
    model.load_parameter_distributions(param_file_path, namespace_dict)
    print(f'\nLoaded parameter distributions.')
    
    parameters = model.get_parameters()
    
    # Baseline analysis
    print(f'\n\nRunning baseline analysis...')
    model.metrics_at_baseline()
    baseline_metrics = model.metrics_at_baseline()
    
    results_dict['Baseline']['MPSP'] = baseline_metrics[0]
    results_dict['Baseline']['GWP100a'] = baseline_metrics[1] 
    results_dict['Baseline']['FEC'] = baseline_metrics[2]
    
    # Get baseline breakdowns
    microalgae_sys.simulate()
    
    # GWP breakdown
    material_GWP_breakdown = microalgae_lca.material_GWP_breakdown
    gwp_breakdown = {
        'feedstock*': microalgae_lca.feedstock_GWP,
        'H2SO4': material_GWP_breakdown.get('H2SO4', 0),
        'NaOH': material_GWP_breakdown.get('NaOH', 0),
        'NH4OH': material_GWP_breakdown.get('NH4OH', 0),
        'CalciumDihydroxide': material_GWP_breakdown.get('CalciumDihydroxide', 0),
        'Ethanol': material_GWP_breakdown.get('Ethanol', 0),
        'Octanol': material_GWP_breakdown.get('Octanol', 0),
        'GlucoAmylase': material_GWP_breakdown.get('GlucoAmylase', 0),
        'AlphaAmylase': material_GWP_breakdown.get('AlphaAmylase', 0),
        'CH4': material_GWP_breakdown.get('CH4', 0),
        'net electricity': microalgae_lca.net_electricity_GWP,
        'direct non-biogenic emissions': microalgae_lca.direct_non_biogenic_emissions_GWP,
    }
    
    # Normalize GWP breakdown
    tot_positive_GWP = sum([v for v in gwp_breakdown.values() if v > 0])
    if tot_positive_GWP > 0:
        for k, v in gwp_breakdown.items():
            gwp_breakdown[k] = v / tot_positive_GWP
    
    results_dict['Baseline']['GWP Breakdown'] = gwp_breakdown
    
    # FEC breakdown (simplified)
    fec_breakdown = {
        'feedstock*': microalgae_lca.feedstock_FEC if hasattr(microalgae_lca, 'feedstock_FEC') else 0,
        'materials': sum(material_GWP_breakdown.values()),  # Approximation
        'net electricity': microalgae_lca.net_electricity_FEC if hasattr(microalgae_lca, 'net_electricity_FEC') else 0,
    }
    
    results_dict['Baseline']['FEC Breakdown'] = fec_breakdown
    
    print(f'\nBaseline MPSP: {baseline_metrics[0]:.3f} USD/kg')
    print(f'Baseline GWP100a: {baseline_metrics[1]:.3f} kg CO2-eq/kg')
    print(f'Baseline FEC: {baseline_metrics[2]:.3f} MJ/kg')
    
    # Monte Carlo simulation
    print(f'\n\nRunning Monte Carlo simulation with {N_simulations_per_mode} samples...')
    
    def notification(i):
        if i % notification_interval == 0:
            print(f'Completed {i}/{N_simulations_per_mode} simulations')
    
    samples = model.sample(N_simulations_per_mode, rule='L')
    model.load_samples(samples)
    model.evaluate(notification=notification, autoload=False, autosave=False)
    
    results = model.table
    
    # Process uncertainty results
    metrics = ['MPSP', 'GWP100a', 'FEC']
    for metric in metrics:
        metric_values = results[metric].values
        percentile_values = np.percentile(metric_values, [p*100 for p in percentiles])
        results_dict['Uncertainty'][metric] = dict(zip(percentiles, percentile_values))
    
    # Sensitivity analysis (Spearman correlation)
    print(f'\n\nCalculating Spearman correlations...')
    
    for metric in metrics:
        correlations = {}
        metric_values = results[metric].values
        
        for param in parameters:
            param_values = results[param.name].values
            # Calculate Spearman correlation
            correlation = np.corrcoef(
                np.argsort(np.argsort(param_values)),
                np.argsort(np.argsort(metric_values))
            )[0, 1]
            correlations[param.name] = correlation
        
        # Sort by absolute correlation
        sorted_correlations = dict(sorted(correlations.items(), 
                                        key=lambda x: abs(x[1]), 
                                        reverse=True))
        results_dict['Sensitivity']['Spearman'][metric] = sorted_correlations
    
    timer.toc()
    
    return results_dict, results

def save_results(results_dict, results, filename_prefix='microalgae_uncertainty'):
    """Save analysis results to files."""
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Save results dictionary
    results_df = pd.DataFrame()
    
    # Baseline results
    baseline_data = {
        'Metric': ['MPSP', 'GWP100a', 'FEC'],
        'Baseline': [
            results_dict['Baseline']['MPSP'],
            results_dict['Baseline']['GWP100a'], 
            results_dict['Baseline']['FEC']
        ]
    }
    
    # Add uncertainty percentiles
    for percentile in percentiles:
        baseline_data[f'P{int(percentile*100)}'] = [
            results_dict['Uncertainty']['MPSP'][percentile],
            results_dict['Uncertainty']['GWP100a'][percentile],
            results_dict['Uncertainty']['FEC'][percentile]
        ]
    
    results_df = pd.DataFrame(baseline_data)
    
    # Save to Excel
    excel_filename = f'{filename_prefix}_results_{timestamp}.xlsx'
    with pd.ExcelWriter(excel_filename) as writer:
        results_df.to_excel(writer, sheet_name='Summary', index=False)
        
        # Save full Monte Carlo results
        results.to_excel(writer, sheet_name='Monte_Carlo', index=False)
        
        # Save sensitivity analysis
        for metric in ['MPSP', 'GWP100a', 'FEC']:
            sensitivity_data = pd.DataFrame(
                list(results_dict['Sensitivity']['Spearman'][metric].items()),
                columns=['Parameter', 'Spearman_Correlation']
            )
            sensitivity_data.to_excel(writer, sheet_name=f'Sensitivity_{metric}', index=False)
        
        # Save GWP breakdown
        gwp_breakdown_df = pd.DataFrame(
            list(results_dict['Baseline']['GWP Breakdown'].items()),
            columns=['Component', 'Fraction']
        )
        gwp_breakdown_df.to_excel(writer, sheet_name='GWP_Breakdown', index=False)
    
    print(f'\nResults saved to {excel_filename}')
    
    # Save summary text report
    txt_filename = f'{filename_prefix}_summary_{timestamp}.txt'
    with open(txt_filename, 'w') as f:
        f.write("Microalgae Biorefinery Uncertainty Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("Baseline Results:\n")
        f.write(f"MPSP: {results_dict['Baseline']['MPSP']:.3f} USD/kg\n")
        f.write(f"GWP100a: {results_dict['Baseline']['GWP100a']:.3f} kg CO2-eq/kg\n")
        f.write(f"FEC: {results_dict['Baseline']['FEC']:.3f} MJ/kg\n\n")
        
        f.write("Uncertainty Analysis (percentiles):\n")
        for metric in ['MPSP', 'GWP100a', 'FEC']:
            f.write(f"\n{metric}:\n")
            for percentile in [0.05, 0.25, 0.5, 0.75, 0.95]:
                value = results_dict['Uncertainty'][metric][percentile]
                f.write(f"  P{int(percentile*100)}: {value:.3f}\n")
        
        f.write("\nTop 10 Most Influential Parameters (by absolute Spearman correlation):\n")
        for metric in ['MPSP', 'GWP100a', 'FEC']:
            f.write(f"\n{metric}:\n")
            correlations = results_dict['Sensitivity']['Spearman'][metric]
            for i, (param, corr) in enumerate(list(correlations.items())[:10]):
                f.write(f"  {i+1}. {param}: {corr:.3f}\n")
    
    print(f'Summary saved to {txt_filename}')

if __name__ == '__main__':
    print("Starting microalgae biorefinery uncertainty analysis...")
    
    try:
        results_dict, results = run_uncertainty_analysis()
        save_results(results_dict, results)
        
        print("\nAnalysis completed successfully!")
        
    except Exception as e:
        print(f"\nError during analysis: {str(e)}")
        import traceback
        traceback.print_exc() 