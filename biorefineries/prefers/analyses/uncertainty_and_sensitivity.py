# -*- coding: utf-8 -*-
"""
Created on 2025-01-XX

Uncertainty and Sensitivity Analysis for LegH Biorefinery - PARALLELIZED WITH ROBUST SAMPLING

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
from biorefineries.prefers.models import create_model
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import thermosteam as tmo
from datetime import datetime
import multiprocessing as mp
from functools import partial

# =============================================================================
# Worker Function for Parallel Execution
# =============================================================================

def evaluate_single_sample(sample_index_and_data, baseline_production_kg_hr):
    """
    Worker function to evaluate a single Monte Carlo sample.
    This function runs in a separate process.
    
    Parameters
    ----------
    sample_index_and_data : tuple
        (index, sample_array) where sample_array contains parameter values
    baseline_production_kg_hr : float
        Baseline production rate for system initialization
    
    Returns
    -------
    tuple
        (index, parameter_values_dict, metric_values_dict, is_valid)
    """
    index, sample = sample_index_and_data
    
    try:
        # Create a fresh model for this worker process
        model = create_model(baseline_production_kg_hr=baseline_production_kg_hr)
        
        # Apply the sample to the model parameters
        for i, param in enumerate(model.parameters):
            param.setter(sample[i])
        
        # Simulate the system
        model.system.simulate()
        
        # Collect parameter values
        param_values = {param.index: sample[i] for i, param in enumerate(model.parameters)}
        
        # Collect metric values
        metric_values = {metric.index: metric() for metric in model.metrics}
        
        # Validate that all metrics are finite
        is_valid = all(np.isfinite(v) for v in metric_values.values())
        
        return (index, param_values, metric_values, is_valid)
        
    except Exception as e:
        # Return invalid result with NaN values
        param_values = {param.index: sample[i] for i, param in enumerate(model.parameters)}
        metric_values = {metric.index: np.nan for metric in model.metrics}
        return (index, param_values, metric_values, False)


# =============================================================================
# Main Execution Block (Must be protected by if __name__ == '__main__')
# =============================================================================

if __name__ == '__main__':
    # Protect multiprocessing code to ensure cross-platform compatibility
    mp.freeze_support()  # Required for Windows
    
    # Set random seed for reproducibility
    np.random.seed(1234)
    
    # Get timestamp for file naming
    now = datetime.now()
    timestamp = now.strftime('%Y.%m.%d-%H.%M')
    
    # Create model with baseline production
    print("="*80)
    print("LEGHEMOGLOBIN BIOREFINERY - UNCERTAINTY & SENSITIVITY ANALYSIS")
    print("="*80)
    print("\nCreating model...")
    
    baseline_production_kg_hr = 275
    model = create_model(baseline_production_kg_hr=baseline_production_kg_hr)
    
    # Display model structure
    print("\nModel Parameters:")
    model.show()
    
    # Set sample size
    N_target = 1000  # Number of VALID samples we want
    batch_size = 100  # Process samples in batches
    
    print(f"\nTarget: {N_target} valid samples")
    print("Using robust sampling strategy with automatic retry for failed simulations...")
    
    # Evaluate baseline (serial, only once)
    print("\nEvaluating baseline scenario...")
    baseline = model.metrics_at_baseline()
    print("\nBaseline Metrics:")
    for name, value in baseline.items():
        print(f"  {name}: {value:.4f}")
    
    # =============================================================================
    # ROBUST PARALLEL MONTE CARLO SIMULATION
    # =============================================================================
    
    print(f"\n{'='*80}")
    print(f"RUNNING ROBUST PARALLEL MONTE CARLO SIMULATION")
    print(f"{'='*80}")
    
    # Determine number of CPU cores to use
    n_cores = mp.cpu_count()
    n_workers = max(1, n_cores - 1)  # Leave one core free
    print(f"\nSystem has {n_cores} CPU cores")
    print(f"Using {n_workers} parallel workers")
    
    # Storage for valid results
    valid_results = []
    total_attempts = 0
    
    # Create partial function with baseline_production_kg_hr fixed
    worker_func = partial(evaluate_single_sample, 
                         baseline_production_kg_hr=baseline_production_kg_hr)
    
    # Robust sampling loop
    while len(valid_results) < N_target:
        # Calculate how many more valid samples we need
        remaining = N_target - len(valid_results)
        
        # Generate batch of samples (with 20% buffer for failed simulations)
        batch_to_generate = min(batch_size, int(remaining * 1.2) + 10)
        
        print(f"\n--- Batch {len(valid_results)//batch_size + 1} ---")
        print(f"Generating {batch_to_generate} samples...")
        
        # Sample parameters using Latin Hypercube sampling
        # Use different seed for each batch to ensure diversity
        batch_samples = model.sample(batch_to_generate, rule='L')
        
        # Prepare sample data for parallel processing
        sample_data = [(total_attempts + i, batch_samples[i]) 
                      for i in range(batch_to_generate)]
        
        # Create multiprocessing pool and execute in parallel
        with mp.Pool(processes=n_workers) as pool:
            batch_results = []
            
            # Process samples in parallel
            for result in pool.imap_unordered(worker_func, sample_data):
                batch_results.append(result)
            
            pool.close()
            pool.join()
        
        # Filter valid results from this batch
        batch_valid_count = 0
        batch_invalid_count = 0
        
        for idx, param_vals, metric_vals, is_valid in batch_results:
            total_attempts += 1
            
            if is_valid:
                # Store valid result with sequential index
                valid_results.append((len(valid_results), param_vals, metric_vals))
                batch_valid_count += 1
            else:
                batch_invalid_count += 1
        
        # Progress update
        success_rate = (len(valid_results) / total_attempts * 100) if total_attempts > 0 else 0
        print(f"  Batch results: {batch_valid_count} valid, {batch_invalid_count} invalid")
        print(f"  Collected {len(valid_results)} / {N_target} valid samples " +
              f"after {total_attempts} attempts ({success_rate:.1f}% success rate)")
        
        # Safety check: if success rate is too low, warn user
        if total_attempts > 100 and success_rate < 50:
            print(f"\n⚠️  WARNING: Success rate is low ({success_rate:.1f}%). " +
                  "This may indicate systematic issues with parameter ranges.")
    
    print(f"\n✓ Robust Monte Carlo simulation complete!")
    print(f"  Valid samples collected: {len(valid_results)}")
    print(f"  Total simulation attempts: {total_attempts}")
    print(f"  Overall success rate: {len(valid_results)/total_attempts*100:.1f}%")
    
    # =============================================================================
    # RECONSTRUCT RESULTS TABLE FROM VALID SAMPLES ONLY
    # =============================================================================
    
    print(f"\nReconstructing results table from valid samples...")
    
    # Build results table structure
    param_indices = [param.index for param in model.parameters]
    metric_indices = [metric.index for metric in model.metrics]
    all_indices = param_indices + metric_indices
    
    # Create DataFrame from valid results only
    data = []
    for idx, param_vals, metric_vals in valid_results:
        row = {}
        row.update(param_vals)
        row.update(metric_vals)
        data.append(row)
    
    results_table = pd.DataFrame(data, columns=all_indices)
    
    print(f"  Results table shape: {results_table.shape}")
    print(f"  All {len(results_table)} samples are guaranteed to be valid (no NaN/inf)")
    
    # Verify no NaN or inf values remain
    has_nan = results_table.isnull().any().any()
    has_inf = np.isinf(results_table.select_dtypes(include=[np.number])).any().any()
    
    if has_nan or has_inf:
        print("\n⚠️  WARNING: Unexpected NaN or inf values detected in final table!")
    else:
        print("  ✓ Verification passed: No NaN or inf values in results table")
    
    # Save results to Excel
    output_file = f'legH_uncertainty_analysis_{timestamp}_{N_target}sims_parallel_robust.xlsx'
    print(f"\nSaving results to {output_file}...")
    results_table.to_excel(output_file)
    
    # Save simulation statistics
    stats_data = {
        'Metric': ['Target valid samples', 'Actual valid samples', 'Total attempts', 
                   'Success rate (%)', 'Failed simulations'],
        'Value': [N_target, len(valid_results), total_attempts,
                 f"{len(valid_results)/total_attempts*100:.2f}",
                 total_attempts - len(valid_results)]
    }
    stats_df = pd.DataFrame(stats_data)
    
    with pd.ExcelWriter(output_file, mode='a', engine='openpyxl') as writer:
        stats_df.to_excel(writer, sheet_name='Simulation Statistics', index=False)
    
    # =============================================================================
    # Spearman's Rank Correlation Analysis
    # =============================================================================
    
    print("\n" + "="*80)
    print("SPEARMAN'S RANK CORRELATION ANALYSIS")
    print("="*80)
    
    # Calculate Spearman's rank correlation
    print("\nCalculating Spearman's rank correlation...")
    
    # No need to clean the table since all samples are already valid
    param_data = results_table[param_indices]
    metric_data = results_table[metric_indices]
    
    # Calculate Spearman correlation manually
    from scipy.stats import spearmanr

    # Initialize with proper MultiIndex
    df_rho = pd.DataFrame(index=pd.MultiIndex.from_tuples(param_indices), 
                          columns=pd.MultiIndex.from_tuples(metric_indices))
    df_p = pd.DataFrame(index=pd.MultiIndex.from_tuples(param_indices), 
                        columns=pd.MultiIndex.from_tuples(metric_indices))

    for param_idx in param_indices:
        for metric_idx in metric_indices:
            rho, p_val = spearmanr(param_data[param_idx], metric_data[metric_idx])
            df_rho.at[param_idx, metric_idx] = rho
            df_p.at[param_idx, metric_idx] = p_val

    print("\nSpearman's Correlation Coefficients (ρ):")
    print(df_rho)
    print("\nP-values:")
    print(df_p)
    
    # Get metric indices
    msp_index = ('Biorefinery', 'MSP [$/kg]')
    tci_index = ('Biorefinery', 'TCI [10^6 $]')
    aoc_index = ('Biorefinery', 'AOC [10^6 $/yr]')
    production_index = ('Biorefinery', 'Actual production [kg/hr]')
    
    # Extract correlation values for each metric
    msp_rho = df_rho[msp_index]
    tci_rho = df_rho[tci_index]
    aoc_rho = df_rho[aoc_index]
    
    # IMPORTANT: Verify that correlations contain both positive and negative values
    print(f"\nMSP Spearman correlation values:")
    print(f"  Min: {msp_rho.min():.4f}")
    print(f"  Max: {msp_rho.max():.4f}")
    print(f"  Number of negative values: {(msp_rho < 0).sum()}")
    print(f"  Number of positive values: {(msp_rho > 0).sum()}")
    
    # Create parameter descriptions for plots
    parameter_descriptions = [p.element_name + ': ' + p.name for p in model.parameters]
    
    # Plot 1D Spearman correlation for MSP
    print("\nGenerating Spearman correlation plots...")
    
    fig_msp, ax_msp = bst.plots.plot_spearman_1d(
        msp_rho,
        index=parameter_descriptions,
        name='MSP [$/kg]',
        color='#A97802',
        sort=True,
    )
    plt.tight_layout()
    plt.savefig(f'legH_spearman_MSP_{timestamp}.png', dpi=300, bbox_inches='tight')
    print(f"Saved MSP Spearman plot to 'legH_spearman_MSP_{timestamp}.png'")
    plt.close()
    
    # Plot 1D Spearman correlation for TCI
    fig_tci, ax_tci = bst.plots.plot_spearman_1d(
        tci_rho,
        index=parameter_descriptions,
        name='TCI [10^6 $]',
        color='#607429',
        sort=True,
    )
    plt.tight_layout()
    plt.savefig(f'legH_spearman_TCI_{timestamp}.png', dpi=300, bbox_inches='tight')
    print(f"Saved TCI Spearman plot to 'legH_spearman_TCI_{timestamp}.png'")
    plt.close()
    
    # Plot 1D Spearman correlation for AOC
    fig_aoc, ax_aoc = bst.plots.plot_spearman_1d(
        aoc_rho,
        index=parameter_descriptions,
        name='AOC [10^6 $/yr]',
        color='#A100A1',
        sort=True,
    )
    plt.tight_layout()
    plt.savefig(f'legH_spearman_AOC_{timestamp}.png', dpi=300, bbox_inches='tight')
    print(f"Saved AOC Spearman plot to 'legH_spearman_AOC_{timestamp}.png'")
    plt.close()
    
    # Plot 2D Spearman correlation (all three metrics)
    print("\nGenerating 2D Spearman correlation plot...")
    
    fig_2d, ax_2d = bst.plots.plot_spearman_2d(
        [msp_rho, tci_rho, aoc_rho],
        index=parameter_descriptions,
        color_wheel=(
            mcolors.to_rgb('#A97802'),  # MSP
            mcolors.to_rgb('#607429'),  # TCI
            mcolors.to_rgb('#A100A1'),  # AOC
        ),
        sort=True,
    )
    fig_2d.set_figwidth(8)
    fig_2d.set_figheight(10)
    plt.tight_layout()
    plt.savefig(f'legH_spearman_2D_{timestamp}.png', dpi=300, bbox_inches='tight')
    print(f"Saved 2D Spearman plot to 'legH_spearman_2D_{timestamp}.png'")
    plt.close()
    
    # =============================================================================
    # Monte Carlo Box Plots
    # =============================================================================
    
    print("\n" + "="*80)
    print("MONTE CARLO UNCERTAINTY BOX PLOTS")
    print("="*80)
    
    # Extract metric values (no need to clean since all data is valid)
    msp_values = results_table[msp_index]
    tci_values = results_table[tci_index]
    aoc_values = results_table[aoc_index]
    
    print(f"\nAll {len(msp_values)} samples are valid (no cleaning required)")
    
    # Define standard colors (hex strings)
    msp_light_hex = "#ffd580"
    msp_dark_hex = "#a97802"
    tci_light_hex = "#b4d7a8"
    tci_dark_hex = "#607429"
    aoc_light_hex = "#e6b3e6"
    aoc_dark_hex = "#a100a1"
    
    # Convert to RGB tuples for plotting
    msp_light_color = mcolors.to_rgb(msp_light_hex)
    msp_dark_color = mcolors.to_rgb(msp_dark_hex)
    tci_light_color = mcolors.to_rgb(tci_light_hex)
    tci_dark_color = mcolors.to_rgb(tci_dark_hex)
    aoc_light_color = mcolors.to_rgb(aoc_light_hex)
    aoc_dark_color = mcolors.to_rgb(aoc_dark_hex)
    
    # Monte Carlo plots
    fig_box, axes_box = plt.subplots(1, 3, figsize=(15, 5))
    plt.sca(axes_box[0])
    bst.plots.plot_montecarlo(
        msp_values.values,
        light_color=msp_light_color,
        dark_color=msp_dark_color,
        positions=[1],
    )
    plt.sca(axes_box[1])
    bst.plots.plot_montecarlo(
        tci_values.values,
        light_color=tci_light_color,
        dark_color=tci_dark_color,
        positions=[1],
    )
    plt.sca(axes_box[2])
    bst.plots.plot_montecarlo(
        aoc_values.values,
        light_color=aoc_light_color,
        dark_color=aoc_dark_color,
        positions=[1],
    )

    axes_box[0].set_ylabel('MSP [$/kg]')
    axes_box[0].set_xticks([])
    axes_box[0].set_title('Minimum Selling Price')

    axes_box[1].set_ylabel('TCI [10^6 $]')
    axes_box[1].set_xticks([])
    axes_box[1].set_title('Total Capital Investment')

    axes_box[2].set_ylabel('AOC [10^6 $/yr]')
    axes_box[2].set_xticks([])
    axes_box[2].set_title('Annual Operating Cost')

    plt.tight_layout()
    plt.savefig(f'legH_montecarlo_boxes_{timestamp}.png', dpi=300, bbox_inches='tight')
    print(f"Saved Monte Carlo box plots to 'legH_montecarlo_boxes_{timestamp}.png'")
    plt.close()
    
    # =============================================================================
    # KDE (Kernel Density Estimation) Plots
    # =============================================================================
    
    print("\n" + "="*80)
    print("KERNEL DENSITY ESTIMATION (KDE) PLOTS")
    print("="*80)
    
    # 1D KDE plot for MSP
    fig_kde1d, ax_kde1d = plt.subplots(figsize=(10, 6))
    ax_kde1d.hist(msp_values, bins=50, density=True, alpha=0.7, 
                  edgecolor='black', color='#FFD580')
    ax_kde1d.set_xlabel('MSP [$/kg]')
    ax_kde1d.set_ylabel('Probability Density')
    ax_kde1d.set_title('Distribution of Minimum Selling Price')
    
    # Add baseline line
    ax_kde1d.axvline(baseline[msp_index], color='r', linestyle='--', linewidth=2, 
                     label=f'Baseline: ${baseline[msp_index]:.4f}/kg')
    ax_kde1d.legend()
    
    plt.tight_layout()
    plt.savefig(f'legH_msp_distribution_{timestamp}.png', dpi=300, bbox_inches='tight')
    print(f"Saved MSP distribution plot to 'legH_msp_distribution_{timestamp}.png'")
    plt.close()
    
    # 2D KDE plot (MSP vs TCI)
    print("\nGenerating 2D KDE plot (MSP vs TCI)...")
    
    try:
        format_units = tmo.units_of_measure.format_units
        ylabel = f"MSP [{format_units('$/kg')}]"
        xlabel = f"TCI [{format_units('10^6 $')}]"
        
        fig_kde2d, ax_kde2d, axes_kde2d = bst.plots.plot_kde(
            y=msp_values.values,
            x=tci_values.values,
            ylabel=ylabel,
            xlabel=xlabel,
            aspect_ratio=1.1,
        )
        plt.savefig(f'legH_kde_MSP_vs_TCI_{timestamp}.png', dpi=300, bbox_inches='tight')
        print(f"Saved 2D KDE plot to 'legH_kde_MSP_vs_TCI_{timestamp}.png'")
        plt.close()
    except Exception as e:
        print(f"WARNING: 2D KDE plot failed: {e}")
        print("  Creating scatter plot as fallback...")
        
        # Fallback: simple scatter plot without KDE
        fig_scatter, ax_scatter = plt.subplots(figsize=(10, 8))
        ax_scatter.scatter(tci_values, msp_values, alpha=0.5, c='#A97802',
                          edgecolors='black', linewidth=0.5)
        ax_scatter.set_xlabel(xlabel)
        ax_scatter.set_ylabel(ylabel)
        ax_scatter.set_title('MSP vs TCI (Scatter Plot)')
        ax_scatter.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'legH_scatter_MSP_vs_TCI_{timestamp}.png', dpi=300, bbox_inches='tight')
        print(f"Saved scatter plot to 'legH_scatter_MSP_vs_TCI_{timestamp}.png'")
        plt.close()
    
    # =============================================================================
    # CRITICAL FIX: Reset Model to Baseline Before Single-Point Sensitivity
    # =============================================================================
    
    print("\n" + "="*80)
    print("RESETTING MODEL TO BASELINE STATE")
    print("="*80)
    
    print("\n⚠️  After Monte Carlo evaluation, the system is in a random state.")
    print("Resetting model to ensure clean baseline for single-point sensitivity analysis...")
    
    # Re-create the model from scratch to get a clean baseline state
    model = create_model(baseline_production_kg_hr=baseline_production_kg_hr)
    
    # Verify baseline is correct
    baseline_verification = model.metrics_at_baseline()
    print("\n✓ Model reset complete. Verified baseline metrics:")
    for name, value in baseline_verification.items():
        print(f"  {name}: {value:.4f}")
    
    # =============================================================================
    # Single-Point Sensitivity Analysis (WITH ROBUST ERROR HANDLING)
    # =============================================================================
    
    print("\n" + "="*80)
    print("SINGLE-POINT SENSITIVITY ANALYSIS")
    print("="*80)
    
    print("\nPerforming single-point sensitivity analysis...")
    print("Note: This may take a few minutes as each parameter bound is evaluated...")
    
    try:
        baseline_sp, lower_sp, upper_sp = model.single_point_sensitivity()
        
        print("\nBASELINE VALUES:")
        print(baseline_sp)
        print("\nLOWER BOUND VALUES:")
        print(lower_sp)
        print("\nUPPER BOUND VALUES:")
        print(upper_sp)
        
        # Check for NaN values in single-point sensitivity results
        has_nan_baseline = baseline_sp.isnull().any()
        has_nan_lower = lower_sp.isnull().any().any() if isinstance(lower_sp, pd.DataFrame) else lower_sp.isnull().any()
        has_nan_upper = upper_sp.isnull().any().any() if isinstance(upper_sp, pd.DataFrame) else upper_sp.isnull().any()
        
        if has_nan_baseline or has_nan_lower or has_nan_upper:
            print("\n⚠️  WARNING: Some single-point sensitivity simulations failed (NaN detected)")
            print("This typically occurs when parameter bounds cause the system to fail.")
            print("Plotting only parameters with valid results...")
            
            # Filter out parameters that caused failures
            valid_params_mask = ~(lower_sp[msp_index].isnull() | upper_sp[msp_index].isnull())
            
            if valid_params_mask.sum() == 0:
                print("\n❌ ERROR: All single-point sensitivity simulations failed!")
                print("Skipping single-point sensitivity plots.")
                print("\nPossible causes:")
                print("  1. Parameter bounds may be too extreme")
                print("  2. System constraints may be violated at bounds")
                print("  3. Model may have convergence issues at boundary conditions")
                baseline_sp = None
                lower_sp = None
                upper_sp = None
            else:
                print(f"\nValid parameters for sensitivity analysis: {valid_params_mask.sum()}/{len(valid_params_mask)}")
                
                # Create filtered data
                baseline_sp_filtered = baseline_sp[valid_params_mask]
                lower_sp_filtered = lower_sp.loc[valid_params_mask, :]
                upper_sp_filtered = upper_sp.loc[valid_params_mask, :]
                
                # Create filtered index
                index_sp = [p.describe(distribution=False) for p in model.parameters]
                index_sp_filtered = [idx for idx, valid in zip(index_sp, valid_params_mask) if valid]
                
                # Plot single-point sensitivity for MSP (filtered)
                print("\nGenerating single-point sensitivity plots for valid parameters...")
                
                fig_sp_msp, ax_sp_msp = bst.plots.plot_single_point_sensitivity(
                    baseline_sp_filtered[msp_index],
                    lower_sp_filtered[msp_index],
                    upper_sp_filtered[msp_index],
                    name='MSP [$/kg]',
                    index=index_sp_filtered,
                    sort=True,
                )
                plt.tight_layout()
                plt.savefig(f'legH_single_point_sensitivity_MSP_{timestamp}.png', 
                            dpi=300, bbox_inches='tight')
                print(f"Saved single-point sensitivity plot for MSP (filtered)")
                plt.close()
                
                # Plot single-point sensitivity for TCI (filtered)
                valid_tci_mask = ~(lower_sp_filtered[tci_index].isnull() | upper_sp_filtered[tci_index].isnull())
                if valid_tci_mask.sum() > 0:
                    fig_sp_tci, ax_sp_tci = bst.plots.plot_single_point_sensitivity(
                        baseline_sp_filtered[tci_index][valid_tci_mask],
                        lower_sp_filtered[tci_index][valid_tci_mask],
                        upper_sp_filtered[tci_index][valid_tci_mask],
                        name='TCI [10^6 $]',
                        index=[idx for idx, valid in zip(index_sp_filtered, valid_tci_mask) if valid],
                        sort=True,
                    )
                    plt.tight_layout()
                    plt.savefig(f'legH_single_point_sensitivity_TCI_{timestamp}.png', 
                                dpi=300, bbox_inches='tight')
                    print(f"Saved single-point sensitivity plot for TCI (filtered)")
                    plt.close()
                else:
                    print("⚠️  Skipping TCI single-point plot: no valid data")
        else:
            # All results are valid - proceed normally
            print("  ✓ All single-point sensitivity simulations successful")
            
            # Create index for plotting (without distribution info)
            index_sp = [p.describe(distribution=False) for p in model.parameters]
            
            # Plot single-point sensitivity for MSP
            fig_sp_msp, ax_sp_msp = bst.plots.plot_single_point_sensitivity(
                baseline_sp[msp_index],
                lower_sp[msp_index],
                upper_sp[msp_index],
                name='MSP [$/kg]',
                index=index_sp,
                sort=True,
            )
            plt.tight_layout()
            plt.savefig(f'legH_single_point_sensitivity_MSP_{timestamp}.png', 
                        dpi=300, bbox_inches='tight')
            print("Saved single-point sensitivity plot for MSP")
            plt.close()
            
            # Plot single-point sensitivity for TCI
            fig_sp_tci, ax_sp_tci = bst.plots.plot_single_point_sensitivity(
                baseline_sp[tci_index],
                lower_sp[tci_index],
                upper_sp[tci_index],
                name='TCI [10^6 $]',
                index=index_sp,
                sort=True,
            )
            plt.tight_layout()
            plt.savefig(f'legH_single_point_sensitivity_TCI_{timestamp}.png', 
                        dpi=300, bbox_inches='tight')
            print("Saved single-point sensitivity plot for TCI")
            plt.close()
    
    except Exception as e:
        print(f"\n❌ ERROR in single-point sensitivity analysis: {e}")
        print("Skipping single-point sensitivity plots.")
        import traceback
        print("\nFull error traceback:")
        traceback.print_exc()
        baseline_sp = None
        lower_sp = None
        upper_sp = None
    
    # =============================================================================
    # Summary Statistics and Reports (CONTINUES AS BEFORE)
    # =============================================================================
    
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    print(f"\nNumber of valid simulations: {len(results_table)}")
    
    print(f"\nMSP Statistics:")
    print(f"  Mean:   ${msp_values.mean():.4f}/kg")
    print(f"  Median: ${msp_values.median():.4f}/kg")
    print(f"  Std:    ${msp_values.std():.4f}/kg")
    print(f"  Min:    ${msp_values.min():.4f}/kg")
    print(f"  Max:    ${msp_values.max():.4f}/kg")
    print(f"  5th percentile:  ${np.percentile(msp_values, 5):.4f}/kg")
    print(f"  95th percentile: ${np.percentile(msp_values, 95):.4f}/kg")
    
    print(f"\nTCI Statistics:")
    print(f"  Mean:   ${tci_values.mean():.2f} million")
    print(f"  Median: ${tci_values.median():.2f} million")
    print(f"  Std:    ${tci_values.std():.2f} million")
    print(f"  Min:    ${tci_values.min():.2f} million")
    print(f"  Max:    ${tci_values.max():.2f} million")
    print(f"  5th percentile:  ${np.percentile(tci_values, 5):.2f} million")
    print(f"  95th percentile: ${np.percentile(tci_values, 95):.2f} million")
    
    print(f"\nAOC Statistics:")
    print(f"  Mean:   ${aoc_values.mean():.2f} million/yr")
    print(f"  Median: ${aoc_values.median():.2f} million/yr")
    print(f"  Std:    ${aoc_values.std():.2f} million/yr")
    print(f"  Min:    ${aoc_values.min():.2f} million/yr")
    print(f"  Max:    ${aoc_values.max():.2f} million/yr")
    print(f"  5th percentile:  ${np.percentile(aoc_values, 5):.2f} million/yr")
    print(f"  95th percentile: ${np.percentile(aoc_values, 95):.2f} million/yr")
    
    print("\nMost Influential Parameters for MSP (by absolute Spearman ρ):")
    sorted_params = msp_rho.abs().sort_values(ascending=False)
    for i, (param_name, abs_rho) in enumerate(sorted_params.items(), 1):
        actual_rho = msp_rho[param_name]
        direction = "↑ increases" if actual_rho > 0 else "↓ decreases"
        print(f"  {i}. {param_name}")
        print(f"      ρ = {actual_rho:+.4f} ({direction} MSP)")
    
    # =============================================================================
    # Save Comprehensive Excel Report (WITH NULL CHECK)
    # =============================================================================
    
    print("\n" + "="*80)
    print("SAVING COMPREHENSIVE EXCEL REPORT")
    print("="*80)
    
    comprehensive_output = f'legH_comprehensive_analysis_{timestamp}_{N_target}sims_robust.xlsx'
    print(f"\nSaving comprehensive results to {comprehensive_output}...")
    
    with pd.ExcelWriter(comprehensive_output) as writer:
        # Raw data
        results_table.to_excel(writer, sheet_name='Raw Monte Carlo Data')
        
        # Spearman correlation
        df_rho.to_excel(writer, sheet_name='Spearman rho')
        df_p.to_excel(writer, sheet_name='Spearman p-values')
        
        # Single-point sensitivity (only if available)
        if baseline_sp is not None:
            try:
                baseline_sp.to_excel(writer, sheet_name='SP Baseline')
                lower_sp.to_excel(writer, sheet_name='SP Lower Bound')
                upper_sp.to_excel(writer, sheet_name='SP Upper Bound')
            except Exception as e:
                print(f"Warning: Could not save single-point sensitivity data: {e}")
        
        # Summary statistics
        summary_stats = pd.DataFrame({
            'MSP [$/kg]': [
                msp_values.mean(), msp_values.median(), msp_values.std(),
                msp_values.min(), msp_values.max(),
                np.percentile(msp_values, 5), np.percentile(msp_values, 95)
            ],
            'TCI [10^6 $]': [
                tci_values.mean(), tci_values.median(), tci_values.std(),
                tci_values.min(), tci_values.max(),
                np.percentile(tci_values, 5), np.percentile(tci_values, 95)
            ],
            'AOC [10^6 $/yr]': [
                aoc_values.mean(), aoc_values.median(), aoc_values.std(),
                aoc_values.min(), aoc_values.max(),
                np.percentile(aoc_values, 5), np.percentile(aoc_values, 95)
            ]
        }, index=['Mean', 'Median', 'Std Dev', 'Min', 'Max', '5th %ile', '95th %ile'])
        summary_stats.to_excel(writer, sheet_name='Summary Statistics')
        
        # Simulation statistics
        stats_df.to_excel(writer, sheet_name='Simulation Statistics', index=False)
        
        # Spearman interpretation for MSP
        spearman_summary = []
        for param_name in msp_rho.index:
            rho = msp_rho[param_name]
            p_val = df_p[msp_index][param_name]
            significance = "Significant (p<0.05)" if p_val < 0.05 else "Not significant"
            direction = "Positive" if rho > 0 else "Negative"
            spearman_summary.append({
                'Parameter': param_name,
                'Spearman ρ (MSP)': rho,
                'P-value': p_val,
                'Significance': significance,
                'Direction': direction,
                'Interpretation': f"{'Increasing' if rho > 0 else 'Decreasing'} this parameter {'increases' if rho > 0 else 'decreases'} MSP"
            })
        
        spearman_df = pd.DataFrame(spearman_summary)
        spearman_df = spearman_df.sort_values('Spearman ρ (MSP)', key=abs, ascending=False)
        spearman_df.to_excel(writer, sheet_name='Spearman Interpretation', index=False)
    
    print(f"Saved comprehensive analysis to: {comprehensive_output}")
    
    # =============================================================================
    # Final Summary
    # =============================================================================
    
    print("\n" + "="*80)
    print("ROBUST PARALLEL ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nSimulation Statistics:")
    print(f"  Target valid samples:     {N_target}")
    print(f"  Actual valid samples:     {len(valid_results)}")
    print(f"  Total simulation attempts: {total_attempts}")
    print(f"  Success rate:             {len(valid_results)/total_attempts*100:.2f}%")
    print(f"  Failed simulations:       {total_attempts - len(valid_results)}")
    print(f"\nFiles generated:")
    print(f"  1. Raw data: {output_file}")
    print(f"  2. Comprehensive report: {comprehensive_output}")
    print(f"  3. MSP Spearman plot: legH_spearman_MSP_{timestamp}.png")
    print(f"  4. TCI Spearman plot: legH_spearman_TCI_{timestamp}.png")
    print(f"  5. AOC Spearman plot: legH_spearman_AOC_{timestamp}.png")
    print(f"  6. 2D Spearman plot: legH_spearman_2D_{timestamp}.png")
    print(f"  7. Monte Carlo box plots: legH_montecarlo_boxes_{timestamp}.png")
    print(f"\nParallel execution used {n_workers} CPU cores")
    print("="*80)