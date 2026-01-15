# -*- coding: utf-8 -*-
"""

Uncertainty and Sensitivity Analysis for LegHemoglobin Production - COMPREHENSIVE ANALYSIS

Strategy:
1. Monte Carlo WITHOUT Production Scale → Spearman correlations, KDE plots, contour plots
2. Monte Carlo WITH Production Scale → Scale effects analysis
3. Single-Point Sensitivity WITH ALL parameters → Tornado diagrams


Created on 2025-11-17 10:16:04

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""


from warnings import filterwarnings
filterwarnings('ignore')

import numpy as np
import pandas as pd
import biosteam as bst
from ..models import create_model
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import thermosteam as tmo
from datetime import datetime
import multiprocessing as mp
from functools import partial

# =============================================================================
# Worker Function for Parallel Execution
# =============================================================================

def evaluate_single_sample(sample_index_and_data, baseline_production_kg_hr, exclude_production_scale=False):
    """
    Worker function to evaluate a single Monte Carlo sample.
    
    Parameters
    ----------
    sample_index_and_data : tuple
        (index, sample_array) where sample_array contains parameter values
    baseline_production_kg_hr : float
        Baseline production rate
    exclude_production_scale : bool
        If True, fix production scale to baseline
    """
    index, sample = sample_index_and_data
    
    try:
        model = create_model(baseline_production_kg_hr=baseline_production_kg_hr)
        
        # Dictionary to store ACTUAL parameter values used in simulation
        param_values = {}
        
        # Apply the sample to the model parameters
        for i, param in enumerate(model.parameters):
            # Check if this is the production scale parameter
            is_production_scale = 'production scale' in param.name.lower()
            
            if exclude_production_scale and is_production_scale:
                # Use fixed baseline value
                param.setter(baseline_production_kg_hr)
                param_values[param.index] = baseline_production_kg_hr  # ← FIXED: Store actual value
            else:
                # Use sampled value
                param.setter(sample[i])
                param_values[param.index] = sample[i]  # ← Store sampled value
        
        model.system.simulate()
        
        # param_values now contains ACTUAL values used (not sample array)
        metric_values = {metric.index: metric() for metric in model.metrics}
        
        is_valid = all(np.isfinite(v) for v in metric_values.values())
        
        return (index, param_values, metric_values, is_valid)
        
    except Exception as e:
        # On error, still need to report what values were SUPPOSED to be used
        param_values = {}
        for i, param in enumerate(model.parameters):
            is_production_scale = 'production scale' in param.name.lower()
            if exclude_production_scale and is_production_scale:
                param_values[param.index] = baseline_production_kg_hr
            else:
                param_values[param.index] = sample[i]
        
        metric_values = {metric.index: np.nan for metric in model.metrics}
        return (index, param_values, metric_values, False)


# =============================================================================
# Robust Parallel Monte Carlo Function
# =============================================================================

def run_monte_carlo(model, N_target, baseline_production_kg_hr, exclude_production_scale=False, 
                   batch_size=100, scenario_name=""):
    """
    Run robust parallel Monte Carlo simulation.
    
    Parameters
    ----------
    model : biosteam.Model
        The model to evaluate
    N_target : int
        Target number of valid samples
    baseline_production_kg_hr : float
        Baseline production rate
    exclude_production_scale : bool
        If True, fix production scale to baseline
    batch_size : int
        Number of samples per batch
    scenario_name : str
        Description of the scenario for logging
        
    Returns
    -------
    results_table : pd.DataFrame
        Table of parameter and metric values
    simulation_stats : dict
        Statistics about the simulation
    """
    print(f"\n{'='*80}")
    print(f"RUNNING MONTE CARLO: {scenario_name}")
    print(f"{'='*80}")
    
    if exclude_production_scale:
        print(f"Production scale FIXED at {baseline_production_kg_hr} kg/hr")
    else:
        print(f"Production scale VARYING (baseline: {baseline_production_kg_hr} kg/hr)")
    
    # Determine number of CPU cores
    n_cores = mp.cpu_count()
    n_workers = max(1, n_cores - 2)
    print(f"\nUsing {n_workers} parallel workers (of {n_cores} cores)")
    
    # Storage for valid results
    valid_results = []
    total_attempts = 0
    
    worker_func = partial(evaluate_single_sample, 
                         baseline_production_kg_hr=baseline_production_kg_hr,
                         exclude_production_scale=exclude_production_scale)
    
    # Robust sampling loop
    while len(valid_results) < N_target:
        remaining = N_target - len(valid_results)
        batch_to_generate = min(batch_size, int(remaining * 1.2) + 10)
        
        print(f"\n--- Batch {len(valid_results)//batch_size + 1} ---")
        print(f"Generating {batch_to_generate} samples...")
        
        # Sample parameters
        batch_samples = model.sample(batch_to_generate, rule='L')
        
        sample_data = [(total_attempts + i, batch_samples[i]) 
                      for i in range(batch_to_generate)]
        
        with mp.Pool(processes=n_workers) as pool:
            batch_results = []
            for result in pool.imap_unordered(worker_func, sample_data):
                batch_results.append(result)
            pool.close()
            pool.join()
        
        batch_valid_count = 0
        batch_invalid_count = 0
        
        for idx, param_vals, metric_vals, is_valid in batch_results:
            total_attempts += 1
            
            if is_valid:
                valid_results.append((len(valid_results), param_vals, metric_vals))
                batch_valid_count += 1
            else:
                batch_invalid_count += 1
        
        success_rate = (len(valid_results) / total_attempts * 100) if total_attempts > 0 else 0
        print(f"  Batch results: {batch_valid_count} valid, {batch_invalid_count} invalid")
        print(f"  Progress: {len(valid_results)}/{N_target} valid samples " +
              f"({success_rate:.1f}% success rate)")
        
        if total_attempts > 100 and success_rate < 50:
            print(f"\n⚠️  WARNING: Low success rate ({success_rate:.1f}%)")
    
    print(f"\n✓ Monte Carlo simulation complete!")
    print(f"  Valid samples: {len(valid_results)}")
    print(f"  Total attempts: {total_attempts}")
    print(f"  Success rate: {len(valid_results)/total_attempts*100:.1f}%")
    
    # Reconstruct results table
    param_indices = [param.index for param in model.parameters]
    metric_indices = [metric.index for metric in model.metrics]
    all_indices = param_indices + metric_indices
    
    data = []
    for idx, param_vals, metric_vals in valid_results:
        row = {}
        row.update(param_vals)  # ← Now contains ACTUAL values used
        row.update(metric_vals)
        data.append(row)
    
    results_table = pd.DataFrame(data, columns=all_indices)
    
    # VERIFICATION STEP: Check if production scale is actually fixed
    if exclude_production_scale:
        production_scale_idx = ('Design', 'Production scale [kg/hr]')
        if production_scale_idx in results_table.columns:
            scale_values = results_table[production_scale_idx]
            unique_scales = scale_values.nunique()
            print(f"\n  Verification: Production scale has {unique_scales} unique value(s)")
            if unique_scales == 1:
                print(f"  ✓ Production scale correctly fixed at {scale_values.iloc[0]:.2f} kg/hr")
            else:
                print(f"  ⚠️  WARNING: Production scale varies despite exclude_production_scale=True!")
                print(f"     Range: {scale_values.min():.2f} - {scale_values.max():.2f} kg/hr")
    
    simulation_stats = {
        'n_valid': len(valid_results),
        'n_attempts': total_attempts,
        'success_rate': len(valid_results)/total_attempts*100,
        'exclude_production_scale': exclude_production_scale
    }
    
    return results_table, simulation_stats


# =============================================================================
# Contour Plot Generation Function
# =============================================================================

def generate_2d_contour_plots(results_table, param_indices, metric_indices, timestamp):
    """Generate 2D scatter plots with color-coded z-values showing joint impact of fermentation parameters."""
    print("\n" + "="*80)
    print("GENERATING 2D SCATTER PLOTS WITH COLOR-CODED METRICS")
    print("="*80)
    
    # Find parameter indices
    titer_idx = next((idx for idx in param_indices if 'titer' in str(idx).lower()), None)
    yield_idx = next((idx for idx in param_indices if 'yield' in str(idx).lower()), None)
    productivity_idx = next((idx for idx in param_indices if 'productivity' in str(idx).lower()), None)
    
    # Find metric indices
    msp_idx = next((idx for idx in metric_indices if 'MSP' in str(idx)), None)
    gwp_idx = next((idx for idx in metric_indices if 'GWP' in str(idx)), None)
    legh_content_idx = next((idx for idx in metric_indices if 'Leghemoglobin content' in str(idx)), None)
    protein_purity_idx = next((idx for idx in metric_indices if 'Protein purity' in str(idx)), None)
    
    if not all([titer_idx, yield_idx, productivity_idx]):
        print("  ⚠️  Could not find all fermentation parameters")
        return
    
    # Define plot configurations
    plot_configs = [
        # MSP scatter plots
        {'x': titer_idx, 'y': yield_idx, 'z': msp_idx, 
         'xlabel': 'Titer [g/L]', 'ylabel': 'Yield [%]', 'zlabel': 'MSP [$/kg]',
         'title': 'MSP_vs_Titer_and_Yield', 'cmap': 'YlOrRd', 'reverse_cmap': False},
        
        {'x': titer_idx, 'y': productivity_idx, 'z': msp_idx,
         'xlabel': 'Titer [g/L]', 'ylabel': 'Productivity [g/L/hr]', 'zlabel': 'MSP [$/kg]',
         'title': 'MSP_vs_Titer_and_Productivity', 'cmap': 'YlOrRd', 'reverse_cmap': False},
        
        {'x': yield_idx, 'y': productivity_idx, 'z': msp_idx,
         'xlabel': 'Yield [%]', 'ylabel': 'Productivity [g/L/hr]', 'zlabel': 'MSP [$/kg]',
         'title': 'MSP_vs_Yield_and_Productivity', 'cmap': 'YlOrRd', 'reverse_cmap': False},
        
        # GWP scatter plots
        {'x': titer_idx, 'y': yield_idx, 'z': gwp_idx,
         'xlabel': 'Titer [g/L]', 'ylabel': 'Yield [%]', 'zlabel': 'GWP [kg CO2-eq/kg]',
         'title': 'GWP_vs_Titer_and_Yield', 'cmap': 'RdYlGn', 'reverse_cmap': True},
        
        {'x': titer_idx, 'y': productivity_idx, 'z': gwp_idx,
         'xlabel': 'Titer [g/L]', 'ylabel': 'Productivity [g/L/hr]', 'zlabel': 'GWP [kg CO2-eq/kg]',
         'title': 'GWP_vs_Titer_and_Productivity', 'cmap': 'RdYlGn', 'reverse_cmap': True},
        
        {'x': yield_idx, 'y': productivity_idx, 'z': gwp_idx,
         'xlabel': 'Yield [%]', 'ylabel': 'Productivity [g/L/hr]', 'zlabel': 'GWP [kg CO2-eq/kg]',
         'title': 'GWP_vs_Yield_and_Productivity', 'cmap': 'RdYlGn', 'reverse_cmap': True},
        
        # Product quality scatter plots
        {'x': titer_idx, 'y': yield_idx, 'z': legh_content_idx,
         'xlabel': 'Titer [g/L]', 'ylabel': 'Yield [%]', 'zlabel': 'Leghemoglobin Content [%]',
         'title': 'LegH_Content_vs_Titer_and_Yield', 'cmap': 'viridis', 'reverse_cmap': False},
        
        {'x': titer_idx, 'y': yield_idx, 'z': protein_purity_idx,
         'xlabel': 'Titer [g/L]', 'ylabel': 'Yield [%]', 'zlabel': 'Protein Purity [%]',
         'title': 'Protein_Purity_vs_Titer_and_Yield', 'cmap': 'plasma', 'reverse_cmap': False},
    ]
    
    # Generate each scatter plot
    for i, config in enumerate(plot_configs, 1):
        if config['z'] is None:
            print(f"  Skipping plot {i}/{len(plot_configs)}: metric not found")
            continue
        
        try:
            print(f"  Generating plot {i}/{len(plot_configs)}: {config['title']}")
            
            # Extract data
            x_data = results_table[config['x']].values
            y_data = results_table[config['y']].values
            z_data = results_table[config['z']].values
            
            # Remove NaN/inf values
            mask = np.isfinite(x_data) & np.isfinite(y_data) & np.isfinite(z_data)
            x_data = x_data[mask]
            y_data = y_data[mask]
            z_data = z_data[mask]
            
            if len(x_data) < 10:
                print(f"    ⚠️  Insufficient valid data ({len(x_data)} points)")
                continue
            
            # Create figure with appropriate size
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Reverse colormap if needed (for GWP where lower is better)
            cmap = config['cmap'] + '_r' if config['reverse_cmap'] else config['cmap']
            
            # Create scatter plot with color-coded z-values
            scatter = ax.scatter(
                x_data, 
                y_data, 
                c=z_data,
                s=30,  # Slightly larger points for visibility
                cmap=cmap,
                alpha=0.7,
                edgecolors='black',
                linewidths=0.3,
                marker='o'
            )
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
            cbar.set_label(config['zlabel'], rotation=270, labelpad=25, fontsize=11)
            cbar.ax.tick_params(labelsize=10)
            
            # Add statistics text box
            stats_text = (
                f'Statistics (n={len(x_data):,}):\n'
                f'{config["zlabel"].split("[")[0].strip()}:\n'
                f'  Mean: {z_data.mean():.4f}\n'
                f'  Median: {np.median(z_data):.4f}\n'
                f'  Std: {z_data.std():.4f}\n'
                f'  Range: [{z_data.min():.4f}, {z_data.max():.4f}]'
            )
            ax.text(
                0.02, 0.98, stats_text,
                transform=ax.transAxes,
                fontsize=9,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'),
                family='monospace'
            )
            
            # Labels and title
            ax.set_xlabel(config['xlabel'], fontsize=12, fontweight='bold')
            ax.set_ylabel(config['ylabel'], fontsize=12, fontweight='bold')
            ax.set_title(config['title'].replace('_', ' '), fontsize=14, fontweight='bold', pad=15)
            
            # Grid
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            ax.set_axisbelow(True)
            
            # Add minor ticks for better readability
            from matplotlib.ticker import AutoMinorLocator
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            ax.tick_params(which='both', direction='in', top=True, right=True)
            ax.tick_params(which='major', length=6, width=1.2)
            ax.tick_params(which='minor', length=3, width=0.8)
            
            # Adjust layout to prevent label cutoff
            plt.tight_layout()
            
            # Save figure
            filename = f"LegH_scatter_{config['title']}_{timestamp}.png"
            plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"    ✓ Saved: {filename}")
            print(f"       Data range: {config['xlabel']} [{x_data.min():.2f}, {x_data.max():.2f}], "
                  f"{config['ylabel']} [{y_data.min():.2f}, {y_data.max():.2f}]")
            print(f"       {config['zlabel']} range: [{z_data.min():.4f}, {z_data.max():.4f}]")
            plt.close()
            
        except Exception as e:
            print(f"    ⚠️  Failed: {e}")
            import traceback
            traceback.print_exc()
            plt.close()
    
    print(f"\n✓ Scatter plot generation complete")
    print(f"  Generated {len([c for c in plot_configs if c['z'] is not None])} color-coded scatter plots")


# =============================================================================
# Main Execution Block
# =============================================================================

if __name__ == '__main__':
    mp.freeze_support()
    np.random.seed(1234)
    
    now = datetime.now()
    timestamp = now.strftime('%Y.%m.%d-%H.%M')
    
    print("="*80)
    print("LEGHEMOGLOBIN PRODUCTION - COMPREHENSIVE UNCERTAINTY & SENSITIVITY ANALYSIS")
    print("PreFerS (Precision Fermentation System)")
    print("="*80)
    print("\nAnalysis Strategy:")
    print("  1. Monte Carlo WITHOUT Production Scale → Spearman, KDE, Contours")
    print("  2. Monte Carlo WITH Production Scale → Scale effects")
    print("  3. Single-Point Sensitivity → Tornado diagrams")
    print("="*80)
    
    # =============================================================================
    # CONFIGURATION PARAMETERS
    # =============================================================================
    baseline_production_kg_hr = 275      # Baseline production rate [kg/hr]
    N_target = 120000                      # Number of valid samples per scenario
    batch_size = 30000                     # Number of samples per batch max 30000
    
    # %%    
    print(f"\nConfiguration:")
    print(f"  Baseline production: {baseline_production_kg_hr} kg/hr")
    print(f"  Target samples: {N_target} per scenario")
    print(f"  Batch size: {batch_size} samples per batch")
    print("="*80)
    
    # Create model
    print("\nCreating model...")
    model = create_model(baseline_production_kg_hr=baseline_production_kg_hr)
    
    print("\nModel Parameters:")
    model.show()
    
    # Evaluate baseline
    print("\nEvaluating baseline scenario...")
    baseline = model.metrics_at_baseline()
    print("\nBaseline Metrics:")
    for name, value in baseline.items():
        if 'GWP' in str(name):
            print(f"  {name}: {value:.6f}")
        else:
            print(f"  {name}: {value:.4f}")
    
    # %%
    # =============================================================================
    # SCENARIO 1: MONTE CARLO WITHOUT PRODUCTION SCALE
    # =============================================================================
    
    results_no_scale, stats_no_scale = run_monte_carlo(
        model=model,
        N_target=N_target,
        baseline_production_kg_hr=baseline_production_kg_hr,
        exclude_production_scale=True,
        batch_size=batch_size,  # ← USE CONFIGURED BATCH SIZE
        scenario_name="Fermentation Parameters Only (Fixed Scale)"
    )
    
    # Save results
    file_no_scale = f'LegH_MC_no_scale_{timestamp}_{N_target}sims.xlsx'
    print(f"\nSaving results to {file_no_scale}...")
    results_no_scale.to_excel(file_no_scale)

    # %%
    # =============================================================================
    # SCENARIO 2: MONTE CARLO WITH PRODUCTION SCALE
    # =============================================================================
    
    results_with_scale, stats_with_scale = run_monte_carlo(
        model=model,
        N_target=N_target,
        baseline_production_kg_hr=baseline_production_kg_hr,
        exclude_production_scale=False,
        batch_size=batch_size,  # ← USE CONFIGURED BATCH SIZE
        scenario_name="All Parameters Including Production Scale"
    )
    
    # Save results
    file_with_scale = f'LegH_MC_with_scale_{timestamp}_{N_target}sims.xlsx'
    print(f"\nSaving results to {file_with_scale}...")
    results_with_scale.to_excel(file_with_scale)
    
    # Get indices (use model from first scenario)
    param_indices = [param.index for param in model.parameters]
    metric_indices = [metric.index for metric in model.metrics]
    
    # Get key metric indices
    msp_idx = ('PreFerS', 'MSP [$/kg]')
    tci_idx = ('PreFerS', 'TCI [10^6 $]')
    aoc_idx = ('PreFerS', 'AOC [10^6 $/yr]')
    gwp_idx = ('PreFerS', 'GWP [kg CO2-eq/kg]')

    # %%
    # =============================================================================
    # SPEARMAN'S RANK CORRELATION (Using NO SCALE data)
    # =============================================================================
    
    print("\n" + "="*80)
    print("SPEARMAN'S RANK CORRELATION (Fermentation Parameters Only)")
    print("="*80)
    
    from scipy.stats import spearmanr
    
    param_data = results_no_scale[param_indices]
    metric_data = results_no_scale[metric_indices]
    
    # =============================================================================
    # FILTER OUT FIXED PARAMETERS (e.g., production scale when fixed)
    # =============================================================================
    
    print("\nChecking for fixed parameters...")
    variable_param_indices = []
    fixed_param_indices = []
    
    for param_idx in param_indices:
        param_values = param_data[param_idx]
        # Check if parameter varies (more than 1 unique value, considering floating point tolerance)
        unique_count = param_values.nunique()
        std_dev = param_values.std()
        
        # Consider parameter as fixed if:
        # 1. Only 1 unique value, OR
        # 2. Standard deviation is extremely small (< 1e-10)
        if unique_count == 1 or std_dev < 1e-10:
            fixed_param_indices.append(param_idx)
            print(f"  ⚠️  Fixed parameter detected: {param_idx}")
            print(f"      Value: {param_values.iloc[0]:.6f}, Std: {std_dev:.2e}")
        else:
            variable_param_indices.append(param_idx)
            print(f"  ✓ Variable parameter: {param_idx}")
            print(f"      Range: [{param_values.min():.4f}, {param_values.max():.4f}], Std: {std_dev:.4f}")
    
    print(f"\nSummary:")
    print(f"  Variable parameters: {len(variable_param_indices)}")
    print(f"  Fixed parameters: {len(fixed_param_indices)}")
    
    if len(variable_param_indices) == 0:
        print("\n⚠️  WARNING: No variable parameters found! Skipping Spearman analysis.")
    else:
        # =============================================================================
        # CALCULATE SPEARMAN CORRELATION (ONLY FOR VARIABLE PARAMETERS)
        # =============================================================================
        
        # Use only variable parameters
        param_data_variable = results_no_scale[variable_param_indices]
        
        df_rho = pd.DataFrame(index=pd.MultiIndex.from_tuples(variable_param_indices), 
                              columns=pd.MultiIndex.from_tuples(metric_indices))
        df_p = pd.DataFrame(index=pd.MultiIndex.from_tuples(variable_param_indices), 
                            columns=pd.MultiIndex.from_tuples(metric_indices))
        
        for param_idx in variable_param_indices:
            for metric_idx in metric_indices:
                try:
                    rho, p_val = spearmanr(param_data_variable[param_idx], metric_data[metric_idx])
                    df_rho.at[param_idx, metric_idx] = rho
                    df_p.at[param_idx, metric_idx] = p_val
                except Exception as e:
                    print(f"  Warning: Failed to calculate Spearman for {param_idx} vs {metric_idx}: {e}")
                    df_rho.at[param_idx, metric_idx] = np.nan
                    df_p.at[param_idx, metric_idx] = np.nan
        
        print("\nSpearman's Correlation Coefficients (ρ) for Variable Parameters:")
        print(df_rho)
        
        # Extract correlations for key metrics
        msp_rho = df_rho[msp_idx]
        tci_rho = df_rho[tci_idx]
        aoc_rho = df_rho[aoc_idx]
        gwp_rho = df_rho[gwp_idx]
        
        # Create parameter descriptions (ONLY for variable parameters)
        parameter_descriptions = []
        for param_idx in variable_param_indices:
            # Find the parameter object
            param_obj = next((p for p in model.parameters if p.index == param_idx), None)
            if param_obj:
                parameter_descriptions.append(param_obj.element_name + ': ' + param_obj.name)
            else:
                parameter_descriptions.append(str(param_idx))
        
        # =============================================================================
        # GENERATE SPEARMAN PLOTS
        # =============================================================================
        
        print("\nGenerating Spearman correlation plots...")
        
        # Individual 1D Spearman plots
        for metric_name, rho_series, color in [
            ('MSP', msp_rho, '#A97802'),
            ('TCI', tci_rho, '#607429'),
            ('AOC', aoc_rho, '#A100A1'),
            ('GWP', gwp_rho, '#E74C3C'),
        ]:
            try:
                # Remove NaN values
                rho_series_clean = rho_series.dropna()
                
                if len(rho_series_clean) == 0:
                    print(f"  ⚠️  No valid correlations for {metric_name}, skipping plot")
                    continue
                
                # Get corresponding descriptions for non-NaN entries
                valid_indices = [i for i, idx in enumerate(variable_param_indices) 
                               if idx in rho_series_clean.index]
                descriptions_clean = [parameter_descriptions[i] for i in valid_indices]
                
                fig, ax = bst.plots.plot_spearman_1d(
                    rho_series_clean,
                    index=descriptions_clean,
                    name=metric_name,
                    color=color,
                    sort=True,
                )
                plt.tight_layout()
                plt.savefig(f'LegH_spearman_{metric_name}_{timestamp}.png', dpi=300, bbox_inches='tight')
                print(f"  ✓ Saved {metric_name} Spearman plot ({len(rho_series_clean)} parameters)")
                plt.close()
            except Exception as e:
                print(f"  ⚠️  Failed {metric_name} Spearman plot: {e}")
        
        # Combined 2D Spearman plot
        try:
            # Check if we have enough variable parameters
            if len(variable_param_indices) < 2:
                print(f"  ⚠️  Too few variable parameters ({len(variable_param_indices)}) for 2D Spearman plot")
            else:
                # Remove NaN values from all series
                msp_rho_clean = msp_rho.dropna()
                tci_rho_clean = tci_rho.dropna()
                aoc_rho_clean = aoc_rho.dropna()
                gwp_rho_clean = gwp_rho.dropna()
                
                # Find common indices (parameters with valid correlations for all metrics)
                common_indices = (set(msp_rho_clean.index) & 
                                set(tci_rho_clean.index) & 
                                set(aoc_rho_clean.index) & 
                                set(gwp_rho_clean.index))
                
                if len(common_indices) < 2:
                    print(f"  ⚠️  Too few parameters with valid correlations ({len(common_indices)}) for 2D Spearman plot")
                else:
                    # Reorder to match common indices
                    common_indices = sorted(list(common_indices), 
                                          key=lambda x: variable_param_indices.index(x))
                    
                    msp_rho_final = msp_rho_clean[common_indices]
                    tci_rho_final = tci_rho_clean[common_indices]
                    aoc_rho_final = aoc_rho_clean[common_indices]
                    gwp_rho_final = gwp_rho_clean[common_indices]
                    
                    # Get corresponding descriptions
                    valid_indices = [i for i, idx in enumerate(variable_param_indices) 
                                   if idx in common_indices]
                    descriptions_final = [parameter_descriptions[i] for i in valid_indices]
                    
                    fig_2d, ax_2d = bst.plots.plot_spearman_2d(
                        rhos=[msp_rho_final, tci_rho_final, aoc_rho_final, gwp_rho_final],
                        index=descriptions_final,
                        color_wheel=(
                            mcolors.to_rgb('#A97802'),
                            mcolors.to_rgb('#607429'),
                            mcolors.to_rgb('#A100A1'),
                            mcolors.to_rgb('#E74C3C'),
                        ),
                        sort=True,
                    )
                    fig_2d.set_figwidth(10)
                    fig_2d.set_figheight(max(6, len(common_indices) * 0.5))  # Adjust height based on parameters
                    plt.tight_layout()
                    plt.savefig(f'LegH_spearman_2D_{timestamp}.png', dpi=300, bbox_inches='tight')
                    print(f"  ✓ Saved combined 2D Spearman plot ({len(common_indices)} parameters)")
                    plt.close()
        except Exception as e:
            print(f"  ⚠️  Failed 2D Spearman plot: {e}")
            import traceback
            traceback.print_exc()
    
    # =============================================================================
    # KDE PLOTS (Using BOTH scenarios)
    # =============================================================================
    
    print("\n" + "="*80)
    print("KERNEL DENSITY ESTIMATE (KDE) PLOTS")
    print("="*80)
    
    msp_values_no_scale = results_no_scale[msp_idx]
    gwp_values_no_scale = results_no_scale[gwp_idx]
    
    msp_values_with_scale = results_with_scale[msp_idx]
    gwp_values_with_scale = results_with_scale[gwp_idx]
    
    sns.set_style("whitegrid")
    
    # =============================================================================
    # 1D KDE PLOTS (Individual distributions)
    # =============================================================================
    
    # KDE for MSP (No Scale)
    print("\nGenerating 1D KDE plots...")
    print("  MSP distribution (fixed scale)...")
    fig_msp_1d, ax_msp_1d = plt.subplots(figsize=(10, 6))
    sns.kdeplot(data=msp_values_no_scale, fill=True, color='#A97802', alpha=0.6, ax=ax_msp_1d, label='Fixed Scale')
    ax_msp_1d.axvline(baseline[msp_idx], color='red', linestyle='--', linewidth=2, 
                   label=f'Baseline: ${baseline[msp_idx]:.4f}/kg')
    ax_msp_1d.axvline(msp_values_no_scale.median(), color='blue', linestyle='--', linewidth=2,
                   label=f'Median: ${msp_values_no_scale.median():.4f}/kg')
    ax_msp_1d.set_xlabel('Minimum Selling Price [$/kg]')
    ax_msp_1d.set_ylabel('Probability Density')
    ax_msp_1d.set_title('Distribution of MSP (Fixed Production Scale)')
    ax_msp_1d.legend()
    ax_msp_1d.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'LegH_kde_1D_MSP_no_scale_{timestamp}.png', dpi=300, bbox_inches='tight')
    print(f"    ✓ Saved")
    plt.close()
    
    # KDE for GWP (No Scale)
    print("  GWP distribution (fixed scale)...")
    fig_gwp_1d, ax_gwp_1d = plt.subplots(figsize=(10, 6))
    sns.kdeplot(data=gwp_values_no_scale, fill=True, color='#E74C3C', alpha=0.6, ax=ax_gwp_1d, label='Fixed Scale')
    ax_gwp_1d.axvline(baseline[gwp_idx], color='red', linestyle='--', linewidth=2,
                   label=f'Baseline: {baseline[gwp_idx]:.4f} kg CO2-eq/kg')
    ax_gwp_1d.axvline(gwp_values_no_scale.median(), color='blue', linestyle='--', linewidth=2,
                   label=f'Median: {gwp_values_no_scale.median():.4f} kg CO2-eq/kg')
    ax_gwp_1d.set_xlabel('Global Warming Potential [kg CO2-eq/kg]')
    ax_gwp_1d.set_ylabel('Probability Density')
    ax_gwp_1d.set_title('Distribution of GWP (Fixed Production Scale)')
    ax_gwp_1d.legend()
    ax_gwp_1d.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'LegH_kde_1D_GWP_no_scale_{timestamp}.png', dpi=300, bbox_inches='tight')
    print(f"    ✓ Saved")
    plt.close()
    
    # =============================================================================
    # 2D KDE PLOTS (Joint distributions - BioSTEAM tutorial style)
    # =============================================================================
    
    print("\n  Generating comprehensive 2D KDE plots (joint GWP-MSP distributions)...")
    
    # Plot 1: Fixed Production Scale (No Scale variation)
    print("    Plot 1: Fixed production scale scenario...")
    try:
        fig_kde_2d_no_scale, ax_kde_no_scale = plt.subplots(figsize=(10, 8))
        
        # Extract data
        x_no_scale = gwp_values_no_scale.values
        y_no_scale = msp_values_no_scale.values
        
        # Remove NaN/inf
        mask_no_scale = np.isfinite(x_no_scale) & np.isfinite(y_no_scale)
        x_no_scale = x_no_scale[mask_no_scale]
        y_no_scale = y_no_scale[mask_no_scale]
        
        # Create 2D KDE contour plot
        sns.kdeplot(
            x=x_no_scale, 
            y=y_no_scale, 
            ax=ax_kde_no_scale,
            fill=True,
            levels=10,
            cmap='YlOrRd',
            alpha=0.6,
            thresh=0.05,
        )
        
        # Add scatter points (subsample for clarity)
        subsample_idx = np.random.choice(len(x_no_scale), size=min(500, len(x_no_scale)), replace=False)
        ax_kde_no_scale.scatter(
            x_no_scale[subsample_idx], 
            y_no_scale[subsample_idx],
            c='#A97802',
            s=10,
            alpha=0.3,
            edgecolors='black',
            linewidths=0.3,
            label='Monte Carlo samples'
        )
        
        # Add baseline point
        ax_kde_no_scale.scatter(
            baseline[gwp_idx], 
            baseline[msp_idx],
            c='red',
            s=200,
            marker='*',
            edgecolors='black',
            linewidths=1.5,
            label='Baseline',
            zorder=10
        )
        
        # Add median point
        ax_kde_no_scale.scatter(
            gwp_values_no_scale.median(), 
            msp_values_no_scale.median(),
            c='blue',
            s=150,
            marker='D',
            edgecolors='black',
            linewidths=1.5,
            label='Median',
            zorder=10
        )
        
        # Labels and formatting
        ax_kde_no_scale.set_xlabel('Global Warming Potential [kg CO2-eq/kg]', fontsize=12)
        ax_kde_no_scale.set_ylabel('Minimum Selling Price [$/kg]', fontsize=12)
        ax_kde_no_scale.set_title('Joint Distribution: GWP vs MSP\n(Fixed Production Scale)', 
                                  fontsize=14, fontweight='bold')
        ax_kde_no_scale.legend(loc='upper right', fontsize=10)
        ax_kde_no_scale.grid(True, alpha=0.3)
        
        # Add statistics box
        stats_text = (
            f'Statistics (n={len(x_no_scale)}):\n'
            f'GWP: {gwp_values_no_scale.mean():.3f} ± {gwp_values_no_scale.std():.3f}\n'
            f'MSP: ${msp_values_no_scale.mean():.3f} ± ${msp_values_no_scale.std():.3f}'
        )
        ax_kde_no_scale.text(
            0.02, 0.98, stats_text,
            transform=ax_kde_no_scale.transAxes,
            fontsize=9,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        )
        
        plt.tight_layout()
        plt.savefig(f'LegH_kde_2D_GWP_MSP_no_scale_{timestamp}.png', dpi=300, bbox_inches='tight')
        print(f"      ✓ Saved: LegH_kde_2D_GWP_MSP_no_scale_{timestamp}.png")
        plt.close()
        
    except Exception as e:
        print(f"      ⚠️  Failed: {e}")
        plt.close()
    
    # Plot 2: Variable Production Scale
    print("    Plot 2: Variable production scale scenario...")
    try:
        fig_kde_2d_with_scale, ax_kde_with_scale = plt.subplots(figsize=(10, 8))
        
        # Extract data
        x_with_scale = gwp_values_with_scale.values
        y_with_scale = msp_values_with_scale.values
        
        # Remove NaN/inf
        mask_with_scale = np.isfinite(x_with_scale) & np.isfinite(y_with_scale)
        x_with_scale = x_with_scale[mask_with_scale]
        y_with_scale = y_with_scale[mask_with_scale]
        
        # Create 2D KDE contour plot
        sns.kdeplot(
            x=x_with_scale, 
            y=y_with_scale, 
            ax=ax_kde_with_scale,
            fill=True,
            levels=10,
            cmap='RdPu',
            alpha=0.6,
            thresh=0.05,
        )
        
        # Add scatter points (subsample for clarity)
        subsample_idx = np.random.choice(len(x_with_scale), size=min(500, len(x_with_scale)), replace=False)
        ax_kde_with_scale.scatter(
            x_with_scale[subsample_idx], 
            y_with_scale[subsample_idx],
            c='#9B59B6',
            s=10,
            alpha=0.3,
            edgecolors='black',
            linewidths=0.3,
            label='Monte Carlo samples'
        )
        
        # Add baseline point
        ax_kde_with_scale.scatter(
            baseline[gwp_idx], 
            baseline[msp_idx],
            c='red',
            s=200,
            marker='*',
            edgecolors='black',
            linewidths=1.5,
            label='Baseline',
            zorder=10
        )
        
        # Add median point
        ax_kde_with_scale.scatter(
            gwp_values_with_scale.median(), 
            msp_values_with_scale.median(),
            c='blue',
            s=150,
            marker='D',
            edgecolors='black',
            linewidths=1.5,
            label='Median',
            zorder=10
        )
        
        # Labels and formatting
        ax_kde_with_scale.set_xlabel('Global Warming Potential [kg CO2-eq/kg]', fontsize=12)
        ax_kde_with_scale.set_ylabel('Minimum Selling Price [$/kg]', fontsize=12)
        ax_kde_with_scale.set_title('Joint Distribution: GWP vs MSP\n(Variable Production Scale)', 
                                    fontsize=14, fontweight='bold')
        ax_kde_with_scale.legend(loc='upper right', fontsize=10)
        ax_kde_with_scale.grid(True, alpha=0.3)
        
        # Add statistics box
        stats_text = (
            f'Statistics (n={len(x_with_scale)}):\n'
            f'GWP: {gwp_values_with_scale.mean():.3f} ± {gwp_values_with_scale.std():.3f}\n'
            f'MSP: ${msp_values_with_scale.mean():.3f} ± ${msp_values_with_scale.std():.3f}'
        )
        ax_kde_with_scale.text(
            0.02, 0.98, stats_text,
            transform=ax_kde_with_scale.transAxes,
            fontsize=9,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        )
        
        plt.tight_layout()
        plt.savefig(f'LegH_kde_2D_GWP_MSP_with_scale_{timestamp}.png', dpi=300, bbox_inches='tight')
        print(f"      ✓ Saved: LegH_kde_2D_GWP_MSP_with_scale_{timestamp}.png")
        plt.close()
        
    except Exception as e:
        print(f"      ⚠️  Failed: {e}")
        plt.close()
    
    # Plot 3: Comparison overlay (both scenarios on same plot)
    print("    Plot 3: Comparison of both scenarios...")
    try:
        fig_kde_2d_compare, ax_kde_compare = plt.subplots(figsize=(12, 9))
        
        # Plot No Scale scenario (red/orange)
        sns.kdeplot(
            x=x_no_scale, 
            y=y_no_scale, 
            ax=ax_kde_compare,
            fill=False,
            levels=5,
            color='#E67E22',
            linewidths=2,
            alpha=0.8,
            label='Fixed Scale'
        )
        
        # Plot With Scale scenario (purple)
        sns.kdeplot(
            x=x_with_scale, 
            y=y_with_scale, 
            ax=ax_kde_compare,
            fill=False,
            levels=5,
            color='#9B59B6',
            linewidths=2,
            alpha=0.8,
            linestyles='--',
            label='Variable Scale'
        )
        
        # Add scatter points for both (smaller samples)
        subsample_no = np.random.choice(len(x_no_scale), size=min(200, len(x_no_scale)), replace=False)
        subsample_with = np.random.choice(len(x_with_scale), size=min(200, len(x_with_scale)), replace=False)
        
        ax_kde_compare.scatter(
            x_no_scale[subsample_no], 
            y_no_scale[subsample_no],
            c='#E67E22',
            s=15,
            alpha=0.3,
            edgecolors='black',
            linewidths=0.2,
        )
        
        ax_kde_compare.scatter(
            x_with_scale[subsample_with], 
            y_with_scale[subsample_with],
            c='#9B59B6',
            s=15,
            alpha=0.3,
            edgecolors='black',
            linewidths=0.2,
            marker='s'
        )
        
        # Add baseline point
        ax_kde_compare.scatter(
            baseline[gwp_idx], 
            baseline[msp_idx],
            c='red',
            s=250,
            marker='*',
            edgecolors='black',
            linewidths=2,
            label='Baseline',
            zorder=10
        )
        
        # Labels and formatting
        ax_kde_compare.set_xlabel('Global Warming Potential [kg CO2-eq/kg]', fontsize=12)
        ax_kde_compare.set_ylabel('Minimum Selling Price [$/kg]', fontsize=12)
        ax_kde_compare.set_title('Comparison: Impact of Production Scale Variation on GWP-MSP Space', 
                                 fontsize=14, fontweight='bold')
        ax_kde_compare.legend(loc='upper right', fontsize=10)
        ax_kde_compare.grid(True, alpha=0.3)
        
        # Add comparison statistics box
        stats_text = (
            f'Fixed Scale (n={len(x_no_scale)}):\n'
            f'  GWP: {gwp_values_no_scale.mean():.3f} ± {gwp_values_no_scale.std():.3f}\n'
            f'  MSP: ${msp_values_no_scale.mean():.3f} ± ${msp_values_no_scale.std():.3f}\n\n'
            f'Variable Scale (n={len(x_with_scale)}):\n'
            f'  GWP: {gwp_values_with_scale.mean():.3f} ± {gwp_values_with_scale.std():.3f}\n'
            f'  MSP: ${msp_values_with_scale.mean():.3f} ± ${msp_values_with_scale.std():.3f}'
        )
        ax_kde_compare.text(
            0.02, 0.98, stats_text,
            transform=ax_kde_compare.transAxes,
            fontsize=8,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7)
        )
        
        plt.tight_layout()
        plt.savefig(f'LegH_kde_2D_GWP_MSP_comparison_{timestamp}.png', dpi=300, bbox_inches='tight')
        print(f"      ✓ Saved: LegH_kde_2D_GWP_MSP_comparison_{timestamp}.png")
        plt.close()
        
    except Exception as e:
        print(f"      ⚠️  Failed: {e}")
        plt.close()
    
    print(f"\n✓ KDE plot generation complete")
    print(f"  Generated 5 plots:")
    print(f"    - 2 × 1D KDE (MSP, GWP)")
    print(f"    - 3 × 2D KDE (No Scale, With Scale, Comparison)")
    
    # =============================================================================
    # 2D CONTOUR PLOTS (Using NO SCALE data)
    # =============================================================================
    
    generate_2d_contour_plots(results_no_scale, param_indices, metric_indices, timestamp)
    
    # =============================================================================
    # PRODUCTION SCALE EFFECTS ANALYSIS
    # =============================================================================
    
    print("\n" + "="*80)
    print("PRODUCTION SCALE EFFECTS ANALYSIS")
    print("="*80)
    
    production_scale_idx = ('Design', 'Production scale [kg/hr]')
    
    if production_scale_idx in results_with_scale.columns:
        scale_values = results_with_scale[production_scale_idx]
        msp_vs_scale = results_with_scale[msp_idx]
        tci_vs_scale = results_with_scale[tci_idx]
        gwp_vs_scale = results_with_scale[gwp_idx]
        aoc_vs_scale = results_with_scale[aoc_idx]
        
        # =============================================================================
        # Plot 1: Scatter plots (existing)
        # =============================================================================
        fig_scale, axes_scale = plt.subplots(1, 2, figsize=(14, 6))
        
        # MSP vs Scale
        axes_scale[0].scatter(scale_values, msp_vs_scale, alpha=0.5, c='#A97802', s=20)
        axes_scale[0].set_xlabel('Production Scale [kg/hr]')
        axes_scale[0].set_ylabel('MSP [$/kg]')
        axes_scale[0].set_title('MSP vs Production Scale')
        axes_scale[0].grid(True, alpha=0.3)
        
        # TCI vs Scale
        axes_scale[1].scatter(scale_values, tci_vs_scale, alpha=0.5, c='#607429', s=20)
        axes_scale[1].set_xlabel('Production Scale [kg/hr]')
        axes_scale[1].set_ylabel('TCI [10^6 $]')
        axes_scale[1].set_title('TCI vs Production Scale')
        axes_scale[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'LegH_scale_effects_scatter_{timestamp}.png', dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved scale effects scatter plot")
        plt.close()
        
        # =============================================================================
        # Plot 2: Percentile band plots (NEW)
        # =============================================================================
        print("\n  Generating percentile band plots...")
        
        # Bin the data by production scale
        n_bins = 20
        scale_min = scale_values.min()
        scale_max = scale_values.max()
        scale_bins = np.linspace(scale_min, scale_max, n_bins + 1)
        scale_centers = (scale_bins[:-1] + scale_bins[1:]) / 2
        
        # Function to calculate percentiles per bin
        def calculate_percentiles_by_bin(x_values, y_values, bins):
            """Calculate percentiles of y for each x bin."""
            bin_indices = np.digitize(x_values, bins) - 1
            
            p05_list = []
            p25_list = []
            p50_list = []
            p75_list = []
            p95_list = []
            valid_centers = []
            
            for i in range(len(bins) - 1):
                mask = bin_indices == i
                if mask.sum() > 0:
                    y_in_bin = y_values[mask]
                    if len(y_in_bin) >= 3:  # Need at least 3 points
                        p05_list.append(np.percentile(y_in_bin, 5))
                        p25_list.append(np.percentile(y_in_bin, 25))
                        p50_list.append(np.percentile(y_in_bin, 50))
                        p75_list.append(np.percentile(y_in_bin, 75))
                        p95_list.append(np.percentile(y_in_bin, 95))
                        valid_centers.append((bins[i] + bins[i+1]) / 2)
            
            return (np.array(valid_centers), 
                    np.array(p05_list), np.array(p25_list), 
                    np.array(p50_list), 
                    np.array(p75_list), np.array(p95_list))
        
        # Calculate percentiles for each metric
        scale_clean = scale_values.values
        msp_clean = msp_vs_scale.values
        tci_clean = tci_vs_scale.values
        gwp_clean = gwp_vs_scale.values
        aoc_clean = aoc_vs_scale.values
        
        # Remove NaN values
        mask = (np.isfinite(scale_clean) & np.isfinite(msp_clean) & 
                np.isfinite(tci_clean) & np.isfinite(gwp_clean) & np.isfinite(aoc_clean))
        scale_clean = scale_clean[mask]
        msp_clean = msp_clean[mask]
        tci_clean = tci_clean[mask]
        gwp_clean = gwp_clean[mask]
        aoc_clean = aoc_clean[mask]
        
        # Calculate percentiles
        scale_c_msp, p05_msp, p25_msp, p50_msp, p75_msp, p95_msp = calculate_percentiles_by_bin(
            scale_clean, msp_clean, scale_bins)
        
        scale_c_tci, p05_tci, p25_tci, p50_tci, p75_tci, p95_tci = calculate_percentiles_by_bin(
            scale_clean, tci_clean, scale_bins)
        
        scale_c_gwp, p05_gwp, p25_gwp, p50_gwp, p75_gwp, p95_gwp = calculate_percentiles_by_bin(
            scale_clean, gwp_clean, scale_bins)
        
        scale_c_aoc, p05_aoc, p25_aoc, p50_aoc, p75_aoc, p95_aoc = calculate_percentiles_by_bin(
            scale_clean, aoc_clean, scale_bins)
        
        # Create 2x2 subplot for percentile bands
        fig_bands, axes_bands = plt.subplots(2, 2, figsize=(16, 12))
        
        # Plot configurations
        plot_configs = [
            {
                'ax': axes_bands[0, 0],
                'scale': scale_c_msp,
                'p05': p05_msp, 'p25': p25_msp, 'p50': p50_msp, 'p75': p75_msp, 'p95': p95_msp,
                'baseline_val': baseline[msp_idx],
                'ylabel': 'MSP [$/kg]',
                'title': 'Minimum Selling Price vs Production Scale',
                'color': '#A97802'
            },
            {
                'ax': axes_bands[0, 1],
                'scale': scale_c_tci,
                'p05': p05_tci, 'p25': p25_tci, 'p50': p50_tci, 'p75': p75_tci, 'p95': p95_tci,
                'baseline_val': baseline[tci_idx],
                'ylabel': 'TCI [10^6 $]',
                'title': 'Total Capital Investment vs Production Scale',
                'color': '#607429'
            },
            {
                'ax': axes_bands[1, 0],
                'scale': scale_c_gwp,
                'p05': p05_gwp, 'p25': p25_gwp, 'p50': p50_gwp, 'p75': p75_gwp, 'p95': p95_gwp,
                'baseline_val': baseline[gwp_idx],
                'ylabel': 'GWP [kg CO2-eq/kg]',
                'title': 'Global Warming Potential vs Production Scale',
                'color': '#E74C3C'
            },
            {
                'ax': axes_bands[1, 1],
                'scale': scale_c_aoc,
                'p05': p05_aoc, 'p25': p25_aoc, 'p50': p50_aoc, 'p75': p75_aoc, 'p95': p95_aoc,
                'baseline_val': baseline[aoc_idx],
                'ylabel': 'AOC [10^6 $/yr]',
                'title': 'Annual Operating Cost vs Production Scale',
                'color': '#A100A1'
            }
        ]
        
        for config in plot_configs:
            ax = config['ax']
            scale_x = config['scale']
            
            # Plot 25-75 percentile band (semi-transparent)
            ax.fill_between(scale_x, config['p25'], config['p75'],
                           color=config['color'], alpha=0.3, 
                           label='25th-75th percentile')
            
            # Plot 5-95 percentile lines (boundary lines)
            ax.plot(scale_x, config['p05'], '--', 
                   color=config['color'], linewidth=1.5, alpha=0.7,
                   label='5th percentile')
            ax.plot(scale_x, config['p95'], '--', 
                   color=config['color'], linewidth=1.5, alpha=0.7,
                   label='95th percentile')
            
            # Plot median line (50th percentile)
            ax.plot(scale_x, config['p50'], '-', 
                   color=config['color'], linewidth=2.5, alpha=0.9,
                   label='Median (50th percentile)')
            
            # Plot baseline as horizontal line
            ax.axhline(y=config['baseline_val'], color='red', 
                      linestyle='-.', linewidth=2,
                      label=f'Baseline ({config["baseline_val"]:.4g})')
            
            # Formatting
            ax.set_xlabel('Production Scale [kg/hr]', fontsize=11)
            ax.set_ylabel(config['ylabel'], fontsize=11)
            ax.set_title(config['title'], fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle=':')
            ax.legend(loc='best', fontsize=9, framealpha=0.9)
            
            # Add minor ticks
            from matplotlib.ticker import AutoMinorLocator
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        
        plt.suptitle('Scale Effects with Uncertainty Bands\n(Monte Carlo Results)', 
                    fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()
        plt.savefig(f'LegH_scale_effects_percentile_bands_{timestamp}.png', 
                   dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved scale effects percentile band plot")
        plt.close()
        
        # =============================================================================
        # Print scale statistics
        # =============================================================================
        print(f"\nProduction Scale Statistics:")
        print(f"  Range: {scale_values.min():.1f} - {scale_values.max():.1f} kg/hr")
        print(f"  Mean: {scale_values.mean():.1f} kg/hr")
        print(f"  Median: {scale_values.median():.1f} kg/hr")
        
        print(f"\nMSP vs Scale:")
        print(f"  At min scale ({scale_values.min():.1f} kg/hr): ${msp_vs_scale[scale_values.idxmin()]:.4f}/kg")
        print(f"  At max scale ({scale_values.max():.1f} kg/hr): ${msp_vs_scale[scale_values.idxmax()]:.4f}/kg")
        print(f"  Median MSP: ${msp_vs_scale.median():.4f}/kg")
        
        print(f"\nGWP vs Scale:")
        print(f"  At min scale: {gwp_vs_scale[scale_values.idxmin()]:.4f} kg CO2-eq/kg")
        print(f"  At max scale: {gwp_vs_scale[scale_values.idxmax()]:.4f} kg CO2-eq/kg")
        print(f"  Median GWP: {gwp_vs_scale.median():.4f} kg CO2-eq/kg")
    
    # =============================================================================
    # SINGLE-POINT SENSITIVITY (ALL PARAMETERS)
    # =============================================================================
    
    print("\n" + "="*80)
    print("SINGLE-POINT SENSITIVITY ANALYSIS (ALL PARAMETERS)")
    print("="*80)
    
    # Reset model
    model = create_model(baseline_production_kg_hr=baseline_production_kg_hr)
    
    print("\nRunning single-point sensitivity...")
    print("This evaluates each parameter at its lower and upper bounds individually.")
    
    try:
        baseline_sp, lower_sp, upper_sp = model.single_point_sensitivity()
        
        print("  ✓ Single-point sensitivity complete")
        
        # Check for NaN values
        has_nan = (baseline_sp.isnull().any() or 
                  lower_sp.isnull().any().any() or 
                  upper_sp.isnull().any().any())
        
        if not has_nan:
            index_sp = [p.describe(distribution=False) for p in model.parameters]
            
            # Plot tornado diagram for MSP
            fig_sp, ax_sp = bst.plots.plot_single_point_sensitivity(
                baseline_sp[msp_idx],
                lower_sp[msp_idx],
                upper_sp[msp_idx],
                name='MSP [$/kg]',
                index=index_sp,
                sort=True,
            )
            plt.tight_layout()
            plt.savefig(f'LegH_tornado_MSP_{timestamp}.png', dpi=300, bbox_inches='tight')
            print(f"  ✓ Saved MSP tornado diagram")
            plt.close()
            
            # Plot tornado diagram for GWP
            fig_sp_gwp, ax_sp_gwp = bst.plots.plot_single_point_sensitivity(
                baseline_sp[gwp_idx],
                lower_sp[gwp_idx],
                upper_sp[gwp_idx],
                name='GWP [kg CO2-eq/kg]',
                index=index_sp,
                sort=True,
            )
            plt.tight_layout()
            plt.savefig(f'LegH_tornado_GWP_{timestamp}.png', dpi=300, bbox_inches='tight')
            print(f"  ✓ Saved GWP tornado diagram")
            plt.close()
        else:
            print("  ⚠️  Some single-point simulations failed (NaN detected)")
    
    except Exception as e:
        print(f"  ⚠️  Single-point sensitivity failed: {e}")
    
    # =============================================================================
    # SAVE COMPREHENSIVE EXCEL REPORT
    # =============================================================================
    
    print("\n" + "="*80)
    print("SAVING COMPREHENSIVE EXCEL REPORT")
    print("="*80)
    
    comprehensive_file = f'LegH_comprehensive_{timestamp}.xlsx'
    print(f"\nSaving to {comprehensive_file}...")
    
    with pd.ExcelWriter(comprehensive_file) as writer:
        # Monte Carlo data
        results_no_scale.to_excel(writer, sheet_name='MC_No_Scale')
        results_with_scale.to_excel(writer, sheet_name='MC_With_Scale')
        
        # Spearman correlation
        df_rho.to_excel(writer, sheet_name='Spearman_rho')
        df_p.to_excel(writer, sheet_name='Spearman_p')
        
        # Single-point sensitivity (if available)
        if 'baseline_sp' in locals() and not has_nan:
            try:
                baseline_sp.to_excel(writer, sheet_name='SP_Baseline')
                lower_sp.to_excel(writer, sheet_name='SP_Lower')
                upper_sp.to_excel(writer, sheet_name='SP_Upper')
            except:
                pass
        
        # Summary statistics
        summary_data = {
            'Metric': ['MSP', 'TCI', 'AOC', 'GWP'],
            'Baseline': [baseline[msp_idx], baseline[tci_idx], baseline[aoc_idx], baseline[gwp_idx]],
            'Mean (No Scale)': [msp_values_no_scale.mean(), results_no_scale[tci_idx].mean(), 
                               results_no_scale[aoc_idx].mean(), gwp_values_no_scale.mean()],
            'Std (No Scale)': [msp_values_no_scale.std(), results_no_scale[tci_idx].std(), 
                              results_no_scale[aoc_idx].std(), gwp_values_no_scale.std()],
            '5th %ile': [np.percentile(msp_values_no_scale, 5), np.percentile(results_no_scale[tci_idx], 5),
                        np.percentile(results_no_scale[aoc_idx], 5), np.percentile(gwp_values_no_scale, 5)],
            '95th %ile': [np.percentile(msp_values_no_scale, 95), np.percentile(results_no_scale[tci_idx], 95),
                         np.percentile(results_no_scale[aoc_idx], 95), np.percentile(gwp_values_no_scale, 95)],
        }
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
        
        # Simulation statistics
        stats_data = {
            'Scenario': ['No Production Scale', 'With Production Scale'],
            'Valid Samples': [stats_no_scale['n_valid'], stats_with_scale['n_valid']],
            'Total Attempts': [stats_no_scale['n_attempts'], stats_with_scale['n_attempts']],
            'Success Rate (%)': [stats_no_scale['success_rate'], stats_with_scale['success_rate']],
        }
        stats_df = pd.DataFrame(stats_data)
        stats_df.to_excel(writer, sheet_name='Simulation_Stats', index=False)
    
    print(f"  ✓ Saved comprehensive report")
    
    # =============================================================================
    # FINAL SUMMARY
    # =============================================================================
    
    print("\n" + "="*80)
    print("COMPREHENSIVE ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nFiles Generated:")
    print(f"  1. MC without scale: {file_no_scale}")
    print(f"  2. MC with scale: {file_with_scale}")
    print(f"  3. Comprehensive report: {comprehensive_file}")
    print(f"  4. Spearman plots: LegH_spearman_*_{timestamp}.png (5 plots)")
    print(f"  5. KDE plots: LegH_kde_*_{timestamp}.png (5 plots)")
    print(f"  6. Contour plots: LegH_contour_*_{timestamp}.png (8 plots)")
    print(f"  7. Scale effects: LegH_scale_effects_{timestamp}.png")
    print(f"  8. Tornado diagrams: LegH_tornado_*_{timestamp}.png (2 plots)")
    print("\n" + "="*80)