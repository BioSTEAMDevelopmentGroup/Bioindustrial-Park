#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Titer, Rate, Yield (TRY) analysis for microalgae biorefinery

Based on succinic project but adapted for microalgae system structure
"""

import numpy as np
import pandas as pd
import os
from .uncertainties import create_model

def run_TRY_analysis():
    """Run Titer, Rate, Yield analysis for microalgae system"""
    print("Creating model for TRY analysis...")
    model, namespace_dict = create_model()
    
    # Load parameter distributions  
    current_dir = os.path.dirname(__file__)
    param_file = os.path.join(current_dir, 'parameter_distributions.xlsx')
    
    if not os.path.exists(param_file):
        raise FileNotFoundError(f"Parameter distributions file not found: {param_file}")
    
    print(f"Loading parameter distributions from {param_file}...")
    model.load_parameter_distributions(param_file, namespace_dict)
    
    # Define TRY parameter ranges based on microalgae system
    # These ranges should be adjusted based on your specific system
    titer_range = np.linspace(1.0, 4.0, 11)  # g/L - caproic acid titer
    yield_range = np.linspace(0.15, 0.35, 11)  # g/g - C6 yield
    
    results = []
    
    # Find TRY parameters from loaded parameters
    titer_param = None
    yield_param = None
    
    for param in model.parameters:
        param_name_lower = param.name.lower()
        if 'titer' in param_name_lower and 'c6' in param_name_lower:
            titer_param = param
            print(f"Found titer parameter: {param.name}")
        elif 'yield' in param_name_lower and 'c6' in param_name_lower:
            yield_param = param
            print(f"Found yield parameter: {param.name}")
    
    if titer_param is None:
        print("Warning: Could not find titer parameter. Available parameters:")
        for param in model.parameters:
            print(f"  - {param.name}")
        # Use first parameter with 'titer' in name as fallback
        for param in model.parameters:
            if 'titer' in param.name.lower():
                titer_param = param
                print(f"Using fallback titer parameter: {param.name}")
                break
    
    if yield_param is None:
        print("Warning: Could not find yield parameter. Available parameters:")
        for param in model.parameters:
            print(f"  - {param.name}")
        # Use first parameter with 'yield' in name as fallback
        for param in model.parameters:
            if 'yield' in param.name.lower():
                yield_param = param
                print(f"Using fallback yield parameter: {param.name}")
                break
    
    if titer_param is None or yield_param is None:
        print("Error: Could not find both titer and yield parameters")
        return None
    
    # Run TRY analysis
    print(f"\nRunning TRY analysis...")
    print(f"Titer range: {min(titer_range):.2f} - {max(titer_range):.2f}")
    print(f"Yield range: {min(yield_range):.3f} - {max(yield_range):.3f}")
    
    baseline_sample = model.get_baseline_sample()
    total_combinations = len(titer_range) * len(yield_range)
    completed = 0
    
    for i, titer in enumerate(titer_range):
        for j, yield_val in enumerate(yield_range):
            sample = baseline_sample.copy()
            sample[titer_param.index] = titer
            sample[yield_param.index] = yield_val
            
            try:
                result = model(sample)
                result['Titer'] = titer
                result['Yield'] = yield_val
                results.append(result)
                completed += 1
                
                if completed % 10 == 0:
                    print(f"Completed {completed}/{total_combinations} combinations")
                    
            except Exception as e:
                print(f"Error at titer={titer:.2f}, yield={yield_val:.3f}: {e}")
                continue
    
    if results:
        results_df = pd.DataFrame(results)
        print(f"\nTRY analysis completed! {len(results)} successful evaluations out of {total_combinations} total combinations")
        
        # Print summary statistics
        print("\nSummary statistics:")
        for col in ['MFSP', 'GWP', 'FEC', 'TCI']:
            if col in results_df.columns:
                print(f"{col}: min={results_df[col].min():.4f}, max={results_df[col].max():.4f}, mean={results_df[col].mean():.4f}")
        
        return results_df
    else:
        print("No successful evaluations in TRY analysis")
        return None

def analyze_TRY_results(results_df, save_path=None):
    """Analyze TRY results and identify optimal regions"""
    if results_df is None:
        print("No results to analyze")
        return
    
    print("\nAnalyzing TRY results...")
    
    # Find optimal conditions (minimum MFSP)
    if 'MFSP' in results_df.columns:
        min_mfsp_idx = results_df['MFSP'].idxmin()
        optimal_conditions = results_df.loc[min_mfsp_idx]
        print(f"\nOptimal conditions (minimum MFSP):")
        print(f"  Titer: {optimal_conditions['Titer']:.2f}")
        print(f"  Yield: {optimal_conditions['Yield']:.3f}")
        print(f"  MFSP: {optimal_conditions['MFSP']:.4f} $/kg")
        if 'GWP' in results_df.columns:
            print(f"  GWP: {optimal_conditions['GWP']:.4f} kg CO2-eq/kg")
    
    # Save results if path provided
    if save_path:
        results_df.to_csv(save_path, index=False)
        print(f"\nResults saved to: {save_path}")

if __name__ == '__main__':
    try:
        results = run_TRY_analysis()
        
        if results is not None:
            # Analyze results
            analyze_TRY_results(results)
            
            # Save results
            results_dir = 'results'
            if not os.path.exists(results_dir):
                os.makedirs(results_dir)
            
            save_path = os.path.join(results_dir, 'TRY_analysis_results.csv')
            analyze_TRY_results(results, save_path)
        
    except Exception as e:
        print(f"Error in TRY analysis: {e}")
        import traceback
        traceback.print_exc()