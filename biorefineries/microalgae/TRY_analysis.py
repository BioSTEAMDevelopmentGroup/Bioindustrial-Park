#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microalgae biorefinery TRY (Titer, Rate, Yield) analysis

This module performs TRY analysis for the microalgae biorefinery system,
analyzing the effects of fermentation titer, productivity rate, and yield
on techno-economic and environmental metrics.

@author: Xingdong Shi
@version: 0.0.1
"""

from warnings import filterwarnings
filterwarnings('ignore')
import copy
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

print('\n\nLoaded microalgae system for TRY analysis.')

# TRY analysis parameters
titer_range = np.linspace(10, 100, 10)  # g/L, fermentation titer
rate_range = np.linspace(0.5, 5.0, 10)  # g/L/h, productivity rate
yield_range = np.linspace(0.1, 0.8, 10)  # g product/g substrate, yield

# Results storage
try_results = {
    'Titer': [],
    'Rate': [], 
    'Yield': [],
    'MPSP': [],
    'GWP100a': [],
    'FEC': [],
    'Production_Rate': [],  # kg/h
    'Annual_Production': [],  # kg/year
    'Conversion': [],  # fraction
}

def get_fermentation_unit():
    """Get the fermentation unit from the system."""
    for unit in microalgae_sys.units:
        if 'ferment' in unit.ID.lower() or 'MCCA' in unit.ID:
            return unit
    return None

def run_try_analysis():
    """Run TRY analysis for microalgae biorefinery."""
    
    timer = TicToc('timer')
    timer.tic()
    
    print(f'\n\nStarting TRY analysis...')
    print(f'Titer range: {titer_range[0]:.1f} - {titer_range[-1]:.1f} g/L')
    print(f'Rate range: {rate_range[0]:.1f} - {rate_range[-1]:.1f} g/L/h')
    print(f'Yield range: {yield_range[0]:.2f} - {yield_range[-1]:.2f} g/g')
    
    # Get fermentation unit
    fermentation_unit = get_fermentation_unit()
    if fermentation_unit is None:
        print("Warning: Could not find fermentation unit")
        return None
    
    print(f'Found fermentation unit: {fermentation_unit.ID}')
    
    total_simulations = len(titer_range) * len(rate_range) * len(yield_range)
    simulation_count = 0
    
    # Store baseline values
    baseline_titer = getattr(fermentation_unit, 'titer', 50.0)  # Default 50 g/L
    baseline_rate = getattr(fermentation_unit, 'productivity', 2.0)  # Default 2 g/L/h
    baseline_yield = getattr(fermentation_unit, 'fermentation_yield', 0.4)  # Default 0.4 g/g
    
    for titer in titer_range:
        for rate in rate_range:
            for yield_val in yield_range:
                simulation_count += 1
                
                if simulation_count % 100 == 0:
                    print(f'Completed {simulation_count}/{total_simulations} simulations')
                
                try:
                    # Set TRY parameters
                    if hasattr(fermentation_unit, 'titer'):
                        fermentation_unit.titer = titer
                    if hasattr(fermentation_unit, 'productivity'):
                        fermentation_unit.productivity = rate
                    if hasattr(fermentation_unit, 'fermentation_yield'):
                        fermentation_unit.fermentation_yield = yield_val
                    
                    # Alternative parameter names
                    if hasattr(fermentation_unit, 'product_titer'):
                        fermentation_unit.product_titer = titer
                    if hasattr(fermentation_unit, 'volumetric_productivity'):
                        fermentation_unit.volumetric_productivity = rate
                    if hasattr(fermentation_unit, 'conversion_efficiency'):
                        fermentation_unit.conversion_efficiency = yield_val
                    
                    # Simulate system
                    microalgae_sys.simulate()
                    
                    # Calculate metrics
                    main_product = analysis_utils.get_main_product_stream(microalgae_sys)
                    mpsp = microalgae_tea.solve_price(main_product) if main_product else np.nan
                    gwp = microalgae_lca.GWP
                    fec = microalgae_lca.FEC
                    
                    # Calculate production metrics
                    product_stream = None
                    for stream in microalgae_sys.products:
                        if stream.get_total_flow('kg/hr') > 0:
                            product_stream = stream
                            break
                    
                    if product_stream:
                        production_rate = product_stream.get_total_flow('kg/hr')
                        annual_production = production_rate * 24 * 365  # kg/year
                    else:
                        production_rate = 0
                        annual_production = 0
                    
                    # Calculate conversion (simplified)
                    conversion = yield_val  # Approximation
                    
                    # Store results
                    try_results['Titer'].append(titer)
                    try_results['Rate'].append(rate)
                    try_results['Yield'].append(yield_val)
                    try_results['MPSP'].append(mpsp)
                    try_results['GWP100a'].append(gwp)
                    try_results['FEC'].append(fec)
                    try_results['Production_Rate'].append(production_rate)
                    try_results['Annual_Production'].append(annual_production)
                    try_results['Conversion'].append(conversion)
                    
                except Exception as e:
                    print(f'Error in simulation {simulation_count}: {str(e)}')
                    # Store NaN values for failed simulations
                    try_results['Titer'].append(titer)
                    try_results['Rate'].append(rate)
                    try_results['Yield'].append(yield_val)
                    try_results['MPSP'].append(np.nan)
                    try_results['GWP100a'].append(np.nan)
                    try_results['FEC'].append(np.nan)
                    try_results['Production_Rate'].append(np.nan)
                    try_results['Annual_Production'].append(np.nan)
                    try_results['Conversion'].append(np.nan)
    
    # Restore baseline values
    if hasattr(fermentation_unit, 'titer'):
        fermentation_unit.titer = baseline_titer
    if hasattr(fermentation_unit, 'productivity'):
        fermentation_unit.productivity = baseline_rate
    if hasattr(fermentation_unit, 'fermentation_yield'):
        fermentation_unit.fermentation_yield = baseline_yield
    
    timer.toc()
    
    return pd.DataFrame(try_results)

def analyze_try_results(results_df):
    """Analyze TRY results and generate summary statistics."""
    
    analysis = {}
    
    # Basic statistics
    metrics = ['MPSP', 'GWP100a', 'FEC', 'Production_Rate', 'Annual_Production']
    
    for metric in metrics:
        valid_data = results_df[metric].dropna()
        if len(valid_data) > 0:
            analysis[metric] = {
                'min': valid_data.min(),
                'max': valid_data.max(),
                'mean': valid_data.mean(),
                'std': valid_data.std(),
                'median': valid_data.median(),
                'q25': valid_data.quantile(0.25),
                'q75': valid_data.quantile(0.75)
            }
    
    # Find optimal conditions
    valid_results = results_df.dropna(subset=['MPSP'])
    if len(valid_results) > 0:
        # Minimum MPSP
        min_mpsp_idx = valid_results['MPSP'].idxmin()
        analysis['optimal_MPSP'] = {
            'titer': valid_results.loc[min_mpsp_idx, 'Titer'],
            'rate': valid_results.loc[min_mpsp_idx, 'Rate'],
            'yield': valid_results.loc[min_mpsp_idx, 'Yield'],
            'MPSP': valid_results.loc[min_mpsp_idx, 'MPSP'],
            'GWP100a': valid_results.loc[min_mpsp_idx, 'GWP100a'],
            'FEC': valid_results.loc[min_mpsp_idx, 'FEC']
        }
        
        # Minimum GWP
        if 'GWP100a' in valid_results.columns:
            min_gwp_idx = valid_results['GWP100a'].idxmin()
            analysis['optimal_GWP'] = {
                'titer': valid_results.loc[min_gwp_idx, 'Titer'],
                'rate': valid_results.loc[min_gwp_idx, 'Rate'],
                'yield': valid_results.loc[min_gwp_idx, 'Yield'],
                'MPSP': valid_results.loc[min_gwp_idx, 'MPSP'],
                'GWP100a': valid_results.loc[min_gwp_idx, 'GWP100a'],
                'FEC': valid_results.loc[min_gwp_idx, 'FEC']
            }
    
    # Correlation analysis
    correlation_metrics = ['Titer', 'Rate', 'Yield', 'MPSP', 'GWP100a', 'FEC']
    correlation_data = results_df[correlation_metrics].dropna()
    if len(correlation_data) > 10:
        analysis['correlations'] = correlation_data.corr()
    
    return analysis

def save_try_results(results_df, analysis, filename_prefix='microalgae_try'):
    """Save TRY analysis results to files."""
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Save to Excel
    excel_filename = f'{filename_prefix}_results_{timestamp}.xlsx'
    with pd.ExcelWriter(excel_filename) as writer:
        # Save full results
        results_df.to_excel(writer, sheet_name='TRY_Results', index=False)
        
        # Save summary statistics
        if 'MPSP' in analysis:
            summary_data = []
            for metric in ['MPSP', 'GWP100a', 'FEC', 'Production_Rate', 'Annual_Production']:
                if metric in analysis:
                    row = {'Metric': metric}
                    row.update(analysis[metric])
                    summary_data.append(row)
            
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(writer, sheet_name='Summary_Statistics', index=False)
        
        # Save optimal conditions
        if 'optimal_MPSP' in analysis:
            optimal_data = []
            for opt_type in ['optimal_MPSP', 'optimal_GWP']:
                if opt_type in analysis:
                    row = {'Optimization': opt_type.replace('optimal_', '')}
                    row.update(analysis[opt_type])
                    optimal_data.append(row)
            
            optimal_df = pd.DataFrame(optimal_data)
            optimal_df.to_excel(writer, sheet_name='Optimal_Conditions', index=False)
        
        # Save correlations
        if 'correlations' in analysis:
            analysis['correlations'].to_excel(writer, sheet_name='Correlations')
    
    print(f'\nTRY results saved to {excel_filename}')
    
    # Save summary text report
    txt_filename = f'{filename_prefix}_summary_{timestamp}.txt'
    with open(txt_filename, 'w') as f:
        f.write("Microalgae Biorefinery TRY Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Analysis Parameters:\n")
        f.write(f"Titer range: {titer_range[0]:.1f} - {titer_range[-1]:.1f} g/L\n")
        f.write(f"Rate range: {rate_range[0]:.1f} - {rate_range[-1]:.1f} g/L/h\n")
        f.write(f"Yield range: {yield_range[0]:.2f} - {yield_range[-1]:.2f} g/g\n")
        f.write(f"Total simulations: {len(results_df)}\n\n")
        
        # Summary statistics
        if 'MPSP' in analysis:
            f.write("Summary Statistics:\n")
            for metric in ['MPSP', 'GWP100a', 'FEC']:
                if metric in analysis:
                    stats = analysis[metric]
                    f.write(f"\n{metric}:\n")
                    f.write(f"  Min: {stats['min']:.3f}\n")
                    f.write(f"  Max: {stats['max']:.3f}\n")
                    f.write(f"  Mean: {stats['mean']:.3f}\n")
                    f.write(f"  Std: {stats['std']:.3f}\n")
                    f.write(f"  Median: {stats['median']:.3f}\n")
        
        # Optimal conditions
        if 'optimal_MPSP' in analysis:
            f.write("\nOptimal Conditions:\n")
            
            f.write("\nMinimum MPSP:\n")
            opt_mpsp = analysis['optimal_MPSP']
            f.write(f"  Titer: {opt_mpsp['titer']:.1f} g/L\n")
            f.write(f"  Rate: {opt_mpsp['rate']:.1f} g/L/h\n")
            f.write(f"  Yield: {opt_mpsp['yield']:.3f} g/g\n")
            f.write(f"  MPSP: {opt_mpsp['MPSP']:.3f} USD/kg\n")
            f.write(f"  GWP: {opt_mpsp['GWP100a']:.3f} kg CO2-eq/kg\n")
            
            if 'optimal_GWP' in analysis:
                f.write("\nMinimum GWP:\n")
                opt_gwp = analysis['optimal_GWP']
                f.write(f"  Titer: {opt_gwp['titer']:.1f} g/L\n")
                f.write(f"  Rate: {opt_gwp['rate']:.1f} g/L/h\n")
                f.write(f"  Yield: {opt_gwp['yield']:.3f} g/g\n")
                f.write(f"  MPSP: {opt_gwp['MPSP']:.3f} USD/kg\n")
                f.write(f"  GWP: {opt_gwp['GWP100a']:.3f} kg CO2-eq/kg\n")
    
    print(f'Summary saved to {txt_filename}')

def run_single_parameter_analysis():
    """Run single parameter analysis (varying one parameter at a time)."""
    
    print('\n\nRunning single parameter analysis...')
    
    # Get fermentation unit
    fermentation_unit = get_fermentation_unit()
    if fermentation_unit is None:
        print("Warning: Could not find fermentation unit")
        return None
    
    # Store baseline values
    baseline_titer = getattr(fermentation_unit, 'titer', 50.0)
    baseline_rate = getattr(fermentation_unit, 'productivity', 2.0)
    baseline_yield = getattr(fermentation_unit, 'fermentation_yield', 0.4)
    
    single_param_results = {
        'Parameter': [],
        'Value': [],
        'MPSP': [],
        'GWP100a': [],
        'FEC': []
    }
    
    # Vary titer
    for titer in titer_range:
        try:
            # Reset to baseline
            if hasattr(fermentation_unit, 'titer'):
                fermentation_unit.titer = titer
            if hasattr(fermentation_unit, 'productivity'):
                fermentation_unit.productivity = baseline_rate
            if hasattr(fermentation_unit, 'fermentation_yield'):
                fermentation_unit.fermentation_yield = baseline_yield
            
            microalgae_sys.simulate()
            
            main_product = analysis_utils.get_main_product_stream(microalgae_sys)
            single_param_results['Parameter'].append('Titer')
            single_param_results['Value'].append(titer)
            single_param_results['MPSP'].append(microalgae_tea.solve_price(main_product) if main_product else np.nan)
            single_param_results['GWP100a'].append(microalgae_lca.GWP)
            single_param_results['FEC'].append(microalgae_lca.FEC)
            
        except Exception as e:
            print(f'Error varying titer to {titer}: {str(e)}')
    
    # Vary rate
    for rate in rate_range:
        try:
            # Reset to baseline
            if hasattr(fermentation_unit, 'titer'):
                fermentation_unit.titer = baseline_titer
            if hasattr(fermentation_unit, 'productivity'):
                fermentation_unit.productivity = rate
            if hasattr(fermentation_unit, 'fermentation_yield'):
                fermentation_unit.fermentation_yield = baseline_yield
            
            microalgae_sys.simulate()
            
            main_product = analysis_utils.get_main_product_stream(microalgae_sys)
            single_param_results['Parameter'].append('Rate')
            single_param_results['Value'].append(rate)
            single_param_results['MPSP'].append(microalgae_tea.solve_price(main_product) if main_product else np.nan)
            single_param_results['GWP100a'].append(microalgae_lca.GWP)
            single_param_results['FEC'].append(microalgae_lca.FEC)
            
        except Exception as e:
            print(f'Error varying rate to {rate}: {str(e)}')
    
    # Vary yield
    for yield_val in yield_range:
        try:
            # Reset to baseline
            if hasattr(fermentation_unit, 'titer'):
                fermentation_unit.titer = baseline_titer
            if hasattr(fermentation_unit, 'productivity'):
                fermentation_unit.productivity = baseline_rate
            if hasattr(fermentation_unit, 'fermentation_yield'):
                fermentation_unit.fermentation_yield = yield_val
            
            microalgae_sys.simulate()
            
            main_product = analysis_utils.get_main_product_stream(microalgae_sys)
            single_param_results['Parameter'].append('Yield')
            single_param_results['Value'].append(yield_val)
            single_param_results['MPSP'].append(microalgae_tea.solve_price(main_product) if main_product else np.nan)
            single_param_results['GWP100a'].append(microalgae_lca.GWP)
            single_param_results['FEC'].append(microalgae_lca.FEC)
            
        except Exception as e:
            print(f'Error varying yield to {yield_val}: {str(e)}')
    
    # Restore baseline values
    if hasattr(fermentation_unit, 'titer'):
        fermentation_unit.titer = baseline_titer
    if hasattr(fermentation_unit, 'productivity'):
        fermentation_unit.productivity = baseline_rate
    if hasattr(fermentation_unit, 'fermentation_yield'):
        fermentation_unit.fermentation_yield = baseline_yield
    
    return pd.DataFrame(single_param_results)

if __name__ == '__main__':
    print("Starting microalgae biorefinery TRY analysis...")
    
    try:
        # Run full TRY analysis
        results_df = run_try_analysis()
        
        if results_df is not None:
            # Analyze results
            analysis = analyze_try_results(results_df)
            
            # Save results
            save_try_results(results_df, analysis)
            
            # Run single parameter analysis
            single_param_df = run_single_parameter_analysis()
            
            if single_param_df is not None:
                # Save single parameter results
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                single_param_filename = f'microalgae_single_param_{timestamp}.xlsx'
                single_param_df.to_excel(single_param_filename, index=False)
                print(f'Single parameter results saved to {single_param_filename}')
            
            print("\nTRY analysis completed successfully!")
        
    except Exception as e:
        print(f"\nError during TRY analysis: {str(e)}")
        import traceback
        traceback.print_exc() 