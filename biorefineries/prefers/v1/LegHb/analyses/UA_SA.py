# -*- coding: utf-8 -*-
"""
Uncertainty Analysis & Sensitivity Analysis (UA/SA) for LegHemoglobin Production
================================================================================

This script bridges the BioSTEAM model with PreFerS-QSD style visualization.

Strategy:
1. Load the BioSTEAM system and create Model object
2. Define Parameters (Inputs) with distributions
3. Define Metrics (Outputs)
4. Run Monte Carlo simulation
5. Generate PreFerS-style plots:
   - Plot 1: Spearman Rank Correlation (parameter importance)
   - Plot 2: Joint Marginal Plot (Cost vs. Carbon trade-off)
   - Plot 3: Tornado Diagram (sensitivity analysis)

Output Structure:
    analyses/
    ├── data/      # Excel, CSV data files
    └── figure/    # PNG, PDF figures

Created on 2026-01-21

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from warnings import filterwarnings
filterwarnings('ignore')

import os
import argparse
import numpy as np
import pandas as pd
import biosteam as bst
from scipy.stats import spearmanr
from datetime import datetime

# Import the LegHb model
from biorefineries.prefers.v1.LegHb._models import create_model
from biorefineries.prefers.v1.LegHb.system import get_available_configs

# Import PreFerS plotting suite
from biorefineries.prefers.v1.plot import style, plots, utils


# =============================================================================
# Analysis Configuration
# =============================================================================

def parse_arguments():
    """Parse command line arguments for configuration."""
    parser = argparse.ArgumentParser(description='LegHemoglobin UA/SA Analysis')
    parser.add_argument('--config', type=str, default='config1',
                        choices=get_available_configs(),
                        help='Process configuration (default: config1)')
    parser.add_argument('--production', type=float, default=275,
                        help='Baseline production rate in kg/hr (default: 275)')
    parser.add_argument('--samples', type=int, default=500,
                        help='Number of Monte Carlo samples (default: 500)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')
    args, _ = parser.parse_known_args()
    return args


# =============================================================================
# Monte Carlo Simulation
# =============================================================================

def run_monte_carlo(model, n_samples, seed=42, exclude_production_scale=True):
    """
    Run Monte Carlo simulation with the given model.
    
    Parameters
    ----------
    model : bst.Model
        BioSTEAM model with parameters and metrics defined
    n_samples : int
        Number of Monte Carlo samples to generate
    seed : int
        Random seed for reproducibility
    exclude_production_scale : bool
        If True, fix production scale to baseline value
        
    Returns
    -------
    results_df : pd.DataFrame
        DataFrame with parameter values and metric results
    """
    print(f"\n{'='*60}")
    print(f"Running Monte Carlo Simulation")
    print(f"{'='*60}")
    print(f"  Target samples: {n_samples}")
    print(f"  Exclude production scale: {exclude_production_scale}")
    
    np.random.seed(seed)
    
    # Generate samples using model's sample method
    samples = model.sample(N=n_samples, rule='L')  # Latin Hypercube sampling
    
    # If excluding production scale, fix it to baseline
    if exclude_production_scale:
        for i, param in enumerate(model.parameters):
            if 'production scale' in param.name.lower():
                # Fix all samples for this parameter to baseline
                samples[:, i] = param.baseline
                print(f"  Fixed '{param.name}' to baseline: {param.baseline}")
    
    # Load samples into model (disable sorting to avoid divide-by-zero with fixed params)
    model.load_samples(samples, sort=False)
    
    # Run evaluation
    print(f"\n  Evaluating {n_samples} samples...")
    model.evaluate(notify=10)
    
    # Get results table
    results_df = model.table.copy()
    
    # Count valid samples
    metric_indices = [m.index for m in model.metrics]
    valid_mask = results_df[metric_indices].notna().all(axis=1)
    n_valid = valid_mask.sum()
    
    print(f"\n  Valid samples: {n_valid} / {n_samples} ({100*n_valid/n_samples:.1f}%)")
    
    return results_df[valid_mask]


# =============================================================================
# Spearman Rank Correlation Analysis
# =============================================================================

def compute_spearman_correlations(results_df, model, exclude_fixed=True):
    """
    Compute Spearman rank correlations between parameters and metrics.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Monte Carlo results
    model : bst.Model
        BioSTEAM model
    exclude_fixed : bool
        Exclude parameters with zero variance (fixed values)
        
    Returns
    -------
    rho_df : pd.DataFrame
        Correlation matrix (parameters × metrics)
    """
    param_indices = [p.index for p in model.parameters]
    metric_indices = [m.index for m in model.metrics]
    
    # Filter out fixed parameters
    if exclude_fixed:
        variable_params = []
        for idx in param_indices:
            if results_df[idx].std() > 1e-10:
                variable_params.append(idx)
            else:
                print(f"  Excluding fixed parameter: {idx}")
        param_indices = variable_params
    
    # Compute correlations
    rho_data = {}
    for metric_idx in metric_indices:
        rho_row = {}
        for param_idx in param_indices:
            try:
                rho, _ = spearmanr(results_df[param_idx], results_df[metric_idx])
                rho_row[param_idx] = rho
            except:
                rho_row[param_idx] = np.nan
        rho_data[metric_idx] = rho_row
    
    rho_df = pd.DataFrame(rho_data).T
    return rho_df


# =============================================================================
# Single-Point Sensitivity for Tornado
# =============================================================================

def compute_single_point_sensitivity(model, metric_name='MSP [$/kg]', 
                                      perturbation=0.1):
    """
    Compute single-point sensitivity for tornado diagram.
    
    Parameters
    ----------
    model : bst.Model
        BioSTEAM model
    metric_name : str
        Name of metric to analyze
    perturbation : float
        Fractional perturbation (±10% by default)
        
    Returns
    -------
    sensitivity_df : pd.DataFrame
        DataFrame with columns ['Parameter', 'Low', 'High']
    baseline_value : float
        Baseline metric value
    """
    print(f"\n{'='*60}")
    print(f"Computing Single-Point Sensitivity for {metric_name}")
    print(f"{'='*60}")
    
    # Get baseline
    baseline_metrics = model.metrics_at_baseline()
    
    # Find the metric
    baseline_value = None
    for key, val in baseline_metrics.items():
        if metric_name in str(key):
            baseline_value = val
            break
    
    if baseline_value is None:
        raise ValueError(f"Metric '{metric_name}' not found")
    
    print(f"  Baseline {metric_name}: {baseline_value:.4f}")
    
    sensitivity_data = []
    
    for param in model.parameters:
        param_baseline = param.baseline
        lb = param_baseline * (1 - perturbation)
        ub = param_baseline * (1 + perturbation)
        
        # Evaluate at low bound
        param.setter(lb)
        try:
            model.system.simulate()
            metrics_low = {m.index: m() for m in model.metrics}
            value_low = None
            for key, val in metrics_low.items():
                if metric_name in str(key):
                    value_low = val
                    break
        except:
            value_low = np.nan
        
        # Evaluate at high bound
        param.setter(ub)
        try:
            model.system.simulate()
            metrics_high = {m.index: m() for m in model.metrics}
            value_high = None
            for key, val in metrics_high.items():
                if metric_name in str(key):
                    value_high = val
                    break
        except:
            value_high = np.nan
        
        # Reset to baseline
        param.setter(param_baseline)
        
        # Extract parameter name (handle tuple indices)
        if isinstance(param.index, tuple):
            param_name = f"{param.index[0]}: {param.index[1]}"
        else:
            param_name = str(param.index)
        
        sensitivity_data.append({
            'Parameter': param_name,
            'Low': value_low if value_low is not None else np.nan,
            'High': value_high if value_high is not None else np.nan
        })
        
        print(f"  {param_name}: Low={value_low:.4f}, High={value_high:.4f}")
    
    # Reset model to baseline
    model.metrics_at_baseline()
    
    return pd.DataFrame(sensitivity_data), baseline_value


# =============================================================================
# Main Execution
# =============================================================================

def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Setup
    CONFIG = args.config
    BASELINE_PRODUCTION = args.production
    N_SAMPLES = args.samples
    SEED = args.seed
    
    print("="*70)
    print("LEGHEMOGLOBIN PRODUCTION - UNCERTAINTY & SENSITIVITY ANALYSIS (UA/SA)")
    print("PreFerS-QSD Style Analysis with New Visualization Suite")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  Process config: {CONFIG}")
    print(f"  Baseline production: {BASELINE_PRODUCTION} kg/hr")
    print(f"  Monte Carlo samples: {N_SAMPLES}")
    print(f"  Random seed: {SEED}")
    
    # Create output directories using new utility
    dirs = utils.get_analysis_dirs(__file__, config=CONFIG)
    print(f"\nOutput directories:")
    print(f"  Data:   {dirs['data']}")
    print(f"  Figure: {dirs['figure']}")
    
    # =========================================================================
    # Step 1: Create Model
    # =========================================================================
    print("\n" + "="*70)
    print("STEP 1: Creating BioSTEAM Model")
    print("="*70)
    
    model = create_model(
        baseline_production_kg_hr=BASELINE_PRODUCTION,
        config=CONFIG,
        verbose=True
    )
    
    print("\nModel Parameters:")
    for param in model.parameters:
        print(f"  - {param.name} [{param.units}]: baseline={param.baseline:.4f}")
    
    print("\nModel Metrics:")
    for metric in model.metrics:
        print(f"  - {metric.name} [{metric.units}]")
    
    # =========================================================================
    # Step 2: Evaluate Baseline
    # =========================================================================
    print("\n" + "="*70)
    print("STEP 2: Evaluating Baseline")
    print("="*70)
    
    baseline_metrics = model.metrics_at_baseline()
    print("\nBaseline Metrics:")
    for name, value in baseline_metrics.items():
        print(f"  {name}: {value:.4f}")
    
    # Get baseline MSP and GWP for later use
    msp_idx = ('PreFerS', 'MSP [$/kg]')
    gwp_idx = ('PreFerS', 'GWP [kg CO2-eq/kg]')
    baseline_msp = baseline_metrics.get(msp_idx, 0)
    baseline_gwp = baseline_metrics.get(gwp_idx, 0)
    
    # =========================================================================
    # Step 3: Run Monte Carlo Simulation
    # =========================================================================
    print("\n" + "="*70)
    print("STEP 3: Monte Carlo Simulation")
    print("="*70)
    
    mc_results = run_monte_carlo(
        model=model,
        n_samples=N_SAMPLES,
        seed=SEED,
        exclude_production_scale=True  # Fix production scale for fair comparison
    )
    
    # Save MC results
    mc_file = os.path.join(dirs['data'], f'monte_carlo_results.xlsx')
    # Flatten MultiIndex columns before saving
    mc_results_flat = mc_results.copy()
    mc_results_flat.columns = [f"{col[0]}_{col[1]}" if isinstance(col, tuple) else str(col) 
                               for col in mc_results.columns]
    mc_results_flat.to_excel(mc_file, index=True)
    print(f"\n  ✓ Saved MC results to: {mc_file}")
    
    # =========================================================================
    # Step 4: Plot 1 - Spearman Rank Correlation Heatmap
    # =========================================================================
    print("\n" + "="*70)
    print("STEP 4: Spearman Rank Correlation Analysis")
    print("="*70)
    
    rho_df = compute_spearman_correlations(mc_results, model, exclude_fixed=True)
    
    # Format index/columns for display
    rho_display = rho_df.copy()
    rho_display.index = [f"{idx[0]}: {idx[1]}" if isinstance(idx, tuple) else str(idx) 
                         for idx in rho_display.index]
    rho_display.columns = [f"{col[1]}" if isinstance(col, tuple) else str(col) 
                           for col in rho_display.columns]
    
    # Generate heatmap
    fig1, ax1 = plots.plot_spearman_heatmap(rho_display.T, figsize=(12, 6))
    ax1.set_title('Spearman Rank Correlation: Parameters vs Metrics', fontweight='bold')
    
    utils.save_figure(fig1, 'spearman_heatmap', dirs['figure'], formats=('png', 'pdf'))
    
    # Save correlation table
    corr_file = os.path.join(dirs['data'], 'spearman_correlations.xlsx')
    rho_display.T.to_excel(corr_file)
    print(f"  ✓ Saved correlation table to: {corr_file}")
    
    # =========================================================================
    # Step 5: Plot 2 - Joint Marginal Plot (MSP vs GWP)
    # =========================================================================
    print("\n" + "="*70)
    print("STEP 5: Joint Marginal Plot (Cost vs. Carbon Trade-off)")
    print("="*70)
    
    # Check if both metrics exist in results
    if msp_idx in mc_results.columns and gwp_idx in mc_results.columns:
        g = plots.plot_joint_marginal(
            data=mc_results,
            x_col=gwp_idx,
            y_col=msp_idx,
            x_label='GWP [kg CO₂-eq/kg]',
            y_label='MSP [$/kg]'
        )
        
        # Add title
        g.figure.suptitle('Cost vs. Carbon Trade-off Analysis', 
                          fontweight='bold', y=1.02)
        
        utils.save_figure(g.figure, 'joint_msp_gwp', dirs['figure'], 
                          formats=('png', 'pdf'))
    else:
        print("  ⚠️  MSP or GWP metric not found in results, skipping joint plot")
    
    # =========================================================================
    # Step 6: Plot 3 - Tornado Diagram for MSP
    # =========================================================================
    print("\n" + "="*70)
    print("STEP 6: Tornado Diagram (MSP Sensitivity)")
    print("="*70)
    
    sensitivity_df, msp_baseline = compute_single_point_sensitivity(
        model=model,
        metric_name='MSP [$/kg]',
        perturbation=0.10  # ±10%
    )
    
    # Filter out rows with NaN
    sensitivity_df = sensitivity_df.dropna()
    
    # Generate tornado plot
    fig3, ax3 = plots.plot_tornado(
        sensitivity_df=sensitivity_df,
        baseline=msp_baseline,
        metric_name='MSP [$/kg]'
    )
    
    utils.save_figure(fig3, 'tornado_msp', dirs['figure'], formats=('png', 'pdf'))
    
    # Save sensitivity table
    sens_file = os.path.join(dirs['data'], 'tornado_sensitivity.xlsx')
    sensitivity_df.to_excel(sens_file, index=False)
    print(f"  ✓ Saved sensitivity table to: {sens_file}")
    
    # =========================================================================
    # Step 7: Additional Distribution Plots
    # =========================================================================
    print("\n" + "="*70)
    print("STEP 7: Distribution Plots")
    print("="*70)
    
    # MSP Distribution
    if msp_idx in mc_results.columns:
        fig4, ax4 = plots.plot_monte_carlo_dist(
            results=mc_results,
            metric_col=msp_idx,
            units='$/kg',
            baseline=baseline_msp,
            color=style.get_color('gold')
        )
        ax4.set_title('MSP Uncertainty Distribution', fontweight='bold')
        utils.save_figure(fig4, 'distribution_msp', dirs['figure'], formats=('png',))
    
    # GWP Distribution
    if gwp_idx in mc_results.columns:
        fig5, ax5 = plots.plot_monte_carlo_dist(
            results=mc_results,
            metric_col=gwp_idx,
            units='kg CO₂-eq/kg',
            baseline=baseline_gwp,
            color=style.get_color('orange')
        )
        ax5.set_title('GWP Uncertainty Distribution', fontweight='bold')
        utils.save_figure(fig5, 'distribution_gwp', dirs['figure'], formats=('png',))
    
    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nResults saved to:")
    print(f"  Data:    {dirs['data']}")
    print(f"  Figures: {dirs['figure']}")
    
    print(f"\nGenerated figures:")
    print(f"  1. spearman_heatmap.png/pdf - Parameter importance ranking")
    print(f"  2. joint_msp_gwp.png/pdf    - Cost vs. Carbon trade-off")
    print(f"  3. tornado_msp.png/pdf      - MSP sensitivity analysis")
    print(f"  4. distribution_msp.png     - MSP uncertainty distribution")
    print(f"  5. distribution_gwp.png     - GWP uncertainty distribution")
    
    print(f"\nGenerated data files:")
    print(f"  1. monte_carlo_results.xlsx   - Raw MC simulation data")
    print(f"  2. spearman_correlations.xlsx - Correlation coefficients")
    print(f"  3. tornado_sensitivity.xlsx   - Sensitivity analysis data")
    
    return mc_results, rho_df, sensitivity_df


if __name__ == '__main__':
    results = main()
