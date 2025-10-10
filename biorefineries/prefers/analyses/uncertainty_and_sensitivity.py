# -*- coding: utf-8 -*-
"""
Created on 2025-01-XX

Uncertainty and Sensitivity Analysis for LegH Biorefinery

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

# Set random seed for reproducibility
np.random.seed(1234)

# Get timestamp for file naming
now = datetime.now()
timestamp = now.strftime('%Y.%m.%d-%H.%M')

# Create model
print("Creating model...")
model = create_model()

# Display model structure
print("\nModel Parameters:")
model.show()

# Set sample size
N = 1000
print(f"\nGenerating {N} samples using Latin Hypercube sampling...")

# Sample parameters using Latin Hypercube sampling
samples = model.sample(N, rule='L')
model.load_samples(samples, sort=True)

# Set exception handling
model.exception_hook = 'warn'  # Warn on errors but continue

# Evaluate baseline
print("\nEvaluating baseline scenario...")
baseline = model.metrics_at_baseline()
print("\nBaseline Metrics:")
for name, value in baseline.items():
    print(f"  {name}: {value:.4f}")

# Run Monte Carlo simulation
print(f"\nRunning Monte Carlo simulation with {N} samples...")
print("This may take several minutes...")

model.evaluate(
    notify=100  # Print progress every 100 simulations
)

print("\nMonte Carlo simulation complete!")

# Get results table
results_table = model.table.copy()

# Save results to Excel
output_file = f'legH_uncertainty_analysis_{timestamp}_{N}sims.xlsx'
print(f"\nSaving results to {output_file}...")
results_table.to_excel(output_file)

# =============================================================================
# Spearman's Rank Correlation Analysis
# =============================================================================

print("\n" + "="*60)
print("SPEARMAN'S RANK CORRELATION ANALYSIS")
print("="*60)

# Calculate Spearman's rank correlation
print("\nCalculating Spearman's rank correlation...")

# Clean the table before calculating Spearman correlation
table_backup = model.table
model.table = model.table.dropna()

# Calculate Spearman correlation
df_rho, df_p = model.spearman_r()

# Restore original table
model.table = table_backup

print("\nSpearman's Correlation Coefficients (ρ):")
print(df_rho)
print("\nP-values:")
print(df_p)

# Get metric indices
msp_index = ('Biorefinery', 'MSP [$/kg]')
tci_index = ('Biorefinery', 'TCI [10^6 $]')
aoc_index = ('Biorefinery', 'AOC [10^6 $/yr]')

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
print("Saved MSP Spearman plot to 'legH_spearman_MSP_{}.png'".format(timestamp))
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
print("Saved TCI Spearman plot to 'legH_spearman_TCI_{}.png'".format(timestamp))
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
print("Saved AOC Spearman plot to 'legH_spearman_AOC_{}.png'".format(timestamp))
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
print("Saved 2D Spearman plot to 'legH_spearman_2D_{}.png'".format(timestamp))
plt.close()

# =============================================================================
# Monte Carlo Box Plots
# =============================================================================

print("\n" + "="*60)
print("MONTE CARLO UNCERTAINTY BOX PLOTS")
print("="*60)

# Extract metric values
msp_values = results_table[msp_index]
tci_values = results_table[tci_index]
aoc_values = results_table[aoc_index]

# Create box plots
fig_box, axes_box = plt.subplots(1, 3, figsize=(15, 5))

# Convert hex colors to RGB normalized to 1
# MSP colors
msp_light_color = mcolors.to_rgb('#FFD580')
msp_dark_color = mcolors.to_rgb('#A97802')

# TCI colors
tci_light_color = mcolors.to_rgb('#B4D7A8')
tci_dark_color = mcolors.to_rgb('#607429')

# AOC colors
aoc_light_color = mcolors.to_rgb('#E6B3E6')
aoc_dark_color = mcolors.to_rgb('#A100A1')

# MSP box plot
plt.sca(axes_box[0])
bst.plots.plot_montecarlo(
    msp_values,
    light_color=msp_light_color,
    dark_color=msp_dark_color,
    positions=[1],
)
axes_box[0].set_ylabel('MSP [$/kg]')
axes_box[0].set_xticks([])
axes_box[0].set_title('Minimum Selling Price')

# TCI box plot
plt.sca(axes_box[1])
bst.plots.plot_montecarlo(
    tci_values,
    light_color=tci_light_color,
    dark_color=tci_dark_color,
    positions=[1],
)
axes_box[1].set_ylabel('TCI [10^6 $]')
axes_box[1].set_xticks([])
axes_box[1].set_title('Total Capital Investment')

# AOC box plot
plt.sca(axes_box[2])
bst.plots.plot_montecarlo(
    aoc_values,
    light_color=aoc_light_color,
    dark_color=aoc_dark_color,
    positions=[1],
)
axes_box[2].set_ylabel('AOC [10^6 $/yr]')
axes_box[2].set_xticks([])
axes_box[2].set_title('Annual Operating Cost')

plt.tight_layout()
plt.savefig(f'legH_montecarlo_boxes_{timestamp}.png', dpi=300, bbox_inches='tight')
print("Saved Monte Carlo box plots to 'legH_montecarlo_boxes_{}.png'".format(timestamp))
plt.close()

# =============================================================================
# KDE (Kernel Density Estimation) Plots
# =============================================================================

print("\n" + "="*60)
print("KERNEL DENSITY ESTIMATION (KDE) PLOTS")
print("="*60)

# 1D KDE plot for MSP
fig_kde1d, ax_kde1d = plt.subplots(figsize=(10, 6))
ax_kde1d.hist(msp_values, bins=50, density=True, alpha=0.7, 
              edgecolor='black', color='#FFD580')
ax_kde1d.set_xlabel('MSP [$/kg]')
ax_kde1d.set_ylabel('Probability Density')
ax_kde1d.set_title('Distribution of Minimum Selling Price')
ax_kde1d.axvline(baseline[msp_index], color='r', linestyle='--', linewidth=2, 
                 label=f'Baseline: ${baseline[msp_index]:.4f}/kg')
ax_kde1d.legend()
plt.tight_layout()
plt.savefig(f'legH_msp_distribution_{timestamp}.png', dpi=300, bbox_inches='tight')
print("Saved MSP distribution plot to 'legH_msp_distribution_{}.png'".format(timestamp))
plt.close()

# 2D KDE plot (MSP vs TCI)
print("\nGenerating 2D KDE plot (MSP vs TCI)...")
format_units = tmo.units_of_measure.format_units
ylabel = f"MSP [{format_units('$/kg')}]"
xlabel = f"TCI [{format_units('10^6 $')}]"

fig_kde2d, ax_kde2d, axes_kde2d = bst.plots.plot_kde(
    y=msp_values,
    x=tci_values,
    ylabel=ylabel,
    xlabel=xlabel,
    aspect_ratio=1.1,
)
plt.savefig(f'legH_kde_MSP_vs_TCI_{timestamp}.png', dpi=300, bbox_inches='tight')
print("Saved 2D KDE plot to 'legH_kde_MSP_vs_TCI_{}.png'".format(timestamp))
plt.close()

# =============================================================================
# Single-Point Sensitivity Analysis
# =============================================================================

print("\n" + "="*60)
print("SINGLE-POINT SENSITIVITY ANALYSIS")
print("="*60)

print("\nPerforming single-point sensitivity analysis...")
baseline_sp, lower_sp, upper_sp = model.single_point_sensitivity()

print("\nBASELINE VALUES:")
print(baseline_sp)
print("\nLOWER BOUND VALUES:")
print(lower_sp)
print("\nUPPER BOUND VALUES:")
print(upper_sp)

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

# =============================================================================
# Summary Statistics and Reports
# =============================================================================

print("\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)

print(f"\nNumber of successful simulations: {len(results_table)}")

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
# Save Comprehensive Excel Report
# =============================================================================

print("\n" + "="*60)
print("SAVING COMPREHENSIVE EXCEL REPORT")
print("="*60)

spearman_output = f'legH_comprehensive_analysis_{timestamp}_{N}sims.xlsx'
print(f"\nSaving comprehensive results to {spearman_output}...")

with pd.ExcelWriter(spearman_output) as writer:
    # Raw data
    results_table.to_excel(writer, sheet_name='Raw Monte Carlo Data')
    
    # Spearman correlation
    df_rho.to_excel(writer, sheet_name='Spearman rho')
    df_p.to_excel(writer, sheet_name='Spearman p-values')
    
    # Single-point sensitivity
    baseline_sp.to_excel(writer, sheet_name='SP Baseline')
    lower_sp.to_excel(writer, sheet_name='SP Lower Bound')
    upper_sp.to_excel(writer, sheet_name='SP Upper Bound')
    
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

print(f"Saved comprehensive analysis to: {spearman_output}")

# =============================================================================
# Final Summary
# =============================================================================

print("\n" + "="*60)
print("ANALYSIS COMPLETE!")
print("="*60)
print(f"\nFiles generated:")
print(f"  1. Raw data: {output_file}")
print(f"  2. Comprehensive analysis: {spearman_output}")
print(f"  3. MSP Spearman plot: legH_spearman_MSP_{timestamp}.png")
print(f"  4. TCI Spearman plot: legH_spearman_TCI_{timestamp}.png")
print(f"  5. AOC Spearman plot: legH_spearman_AOC_{timestamp}.png")
print(f"  6. 2D Spearman plot: legH_spearman_2D_{timestamp}.png")
print(f"  7. Monte Carlo box plots: legH_montecarlo_boxes_{timestamp}.png")
print(f"  8. MSP distribution: legH_msp_distribution_{timestamp}.png")
print(f"  9. 2D KDE plot: legH_kde_MSP_vs_TCI_{timestamp}.png")
print(f" 10. MSP single-point sensitivity: legH_single_point_sensitivity_MSP_{timestamp}.png")
print(f" 11. TCI single-point sensitivity: legH_single_point_sensitivity_TCI_{timestamp}.png")
print("="*60)