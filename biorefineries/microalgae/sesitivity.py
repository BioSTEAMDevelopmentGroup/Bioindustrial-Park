#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BioSTEAM Sensitivity Analysis for Microalgae MCCA Biorefinery
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import biosteam as bst
from biosteam.evaluation import Model, Metric
from chaospy import distributions as shape
import os
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy.stats import spearmanr
import warnings
from microalgae.system import compute_labor_cost
from microalgae._chemicals import chems
from biosteam import main_flowsheet
from microalgae import system
from microalgae.system import create_microalgae_MCCA_production_sys
from microalgae.tea import create_tea

# Suppress hot utility load warnings
warnings.filterwarnings("ignore", message="Hot utility load is negative.*", category=RuntimeWarning)

# Set up plotting configuration
def setup_plotting_config():
    """Configure plotting environment"""
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

# Initialize system
def initialize_system():
    """Initialize the biorefinery system and TEA"""
    # Set up thermodynamics
    bst.settings.set_thermo(chems)
    main_flowsheet.clear()
    
    # Create system
    from microalgae.system import create_microalgae_MCCA_production_sys
    system = create_microalgae_MCCA_production_sys()
    system.simulate()
    
    # Create TEA object
    from microalgae.tea import create_tea
    u = system.flowsheet.unit
    dry_tpd = u.U101.ins[0].F_mass * 24 / 1000  # kg/h -> t/d
    tea = create_tea(
        system=system, 
        IRR=0.10, 
        duration=(2024, 2045),
        depreciation='MACRS7', 
        income_tax=0.21, 
        operating_days=330,
        lang_factor=None, 
        construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, 
        startup_FOCfrac=1, 
        startup_salesfrac=0.5,
        startup_VOCfrac=0.75, 
        WC_over_FCI=0.05,
        finance_interest=0.08, 
        finance_years=10, 
        finance_fraction=0.4,
        OSBL_units=(u.CT, u.CWP, u.ADP, u.PWC),
        warehouse=0.04, 
        site_development=0.09, 
        additional_piping=0.045,
        proratable_costs=0.10, 
        field_expenses=0.10, 
        construction=0.20,
        contingency=0.10, 
        other_indirect_costs=0.10, 
        labor_cost=max(0.5e6, compute_labor_cost(dry_tpd)),
        labor_burden=0.90, 
        property_insurance=0.007, 
        maintenance=0.03, 
        boiler_turbogenerator=None,
        steam_power_depreciation='MACRS20'
    )
    
    return system, tea

# Create evaluation model
def create_evaluation_model(microalgae_sys, microalgae_tea):
    from biosteam.evaluation import Model, Metric
    import chaospy as cp
    
    u = microalgae_sys.flowsheet.unit  # 修正：使用flowsheet.unit而不是unit
    s = microalgae_sys.flowsheet.stream  # 修正：使用flowsheet.stream而不是stream
    tea = microalgae_tea
    
    # 首先定义所有setter函数
    def set_microalgae_price(price):
        u.U101.ins[0].price = price
        
    def set_caproic_acid_price(price):
        s.caproic_acid_product.price = price
        
    def set_caproic_acid_yield_factor(factor):
        u.R301.caproic_acid_yield_factor = factor
    
    def set_operating_days(days):
        tea.operating_days = days
        
    def set_maintenance_factor(factor):
        tea.maintenance = factor
        
    # 创建模型
    model = Model(microalgae_sys)
    
    # 定义指标
    model.metrics = [
        Metric('ROI', lambda: tea.ROI * 100, '%'),
        Metric('MSP', lambda: tea.solve_price(s.caproic_acid_product), '$/kg'),
        Metric('TCI', lambda: tea.TCI, '$'),
        Metric('VOC', lambda: tea.VOC, '$'),
        Metric('FOC', lambda: tea.FOC, '$')
    ]
    
    # 之后再注册参数 - 修改为使用命名参数
    model.parameter(setter=set_microalgae_price, element=u.U101.ins[0], 
                  name='Microalgae Price', kind='isolated', 
                  units='$/kg', distribution=shape.Triangle(-0.2, -0.1, 0.1))
    
    model.parameter(setter=set_caproic_acid_price, element=s.caproic_acid_product,
                  name='Caproic Acid Price', kind='isolated',
                  units='$/kg', distribution=shape.Triangle(2.0, 2.89, 4.0))
    
    model.parameter(setter=set_caproic_acid_yield_factor, element=u.R301,
                  name='Caproic Acid Yield Factor', kind='isolated',
                  distribution=shape.Triangle(0.7, 1.0, 1.3))
    
    model.parameter(setter=set_operating_days, element=tea,
                  name='Operating Days', kind='isolated',
                  distribution=shape.Triangle(300, 330, 350))
    
    model.parameter(setter=set_maintenance_factor, element=tea,
                  name='Maintenance Factor', kind='isolated',
                  distribution=shape.Triangle(0.02, 0.03, 0.05))
    
    # Print registered parameters to debug
    print("Registered parameters:")
    for i, param in enumerate(model.parameters):
        print(f"{i+1}. {param.name} (index: {param.index})")
    
    return model

# Univariate sensitivity analysis
def run_univariate_analysis(model, parameter_name, values=None, metrics=None):
    """Run univariate sensitivity analysis"""
    # Find parameter
    parameter = None
    for param in model.parameters:
        if param.name == parameter_name:
            parameter = param
            break
            
    if parameter is None:
        raise ValueError(f"Parameter '{parameter_name}' not found")
    
    # Generate value range
    if values is None:
            dist = parameter.distribution
            if hasattr(dist, 'lower') and hasattr(dist, 'upper'):
                values = np.linspace(dist.lower, dist.upper, 11)
            else:
                raise ValueError(f"Parameter '{parameter_name}' has no defined distribution range")
    
    print(f"Running univariate sensitivity analysis: {parameter_name}")
    print(f"Parameter range: {min(values)} - {max(values)}")
    
    # Create sample points
    baseline = model.get_baseline_sample()
    results = []
    
    # Evaluate each value
    for value in values:
        sample = baseline.copy()
        sample[parameter.index] = value
        # Directly call model for results
        result = model(sample)
        result[parameter.name] = value
        results.append(result)
    
    # Combine results
    df = pd.DataFrame(results)
    
    return df

# Bivariate sensitivity analysis
def run_bivariate_analysis(model, param_name1, param_name2, metric_name='ROI', 
                          values1=None, values2=None, save=True):
    """Run bivariate sensitivity analysis"""
    # Find matching parameters
    param1 = param2 = None
    for param in model.parameters:
        if param.name == param_name1:
            param1 = param
        elif param.name == param_name2:
            param2 = param
        if param1 and param2:
            break
    
    if not param1:
        raise ValueError(f"Parameter '{param_name1}' not found")
    if not param2:
        raise ValueError(f"Parameter '{param_name2}' not found")
    
    # Set default value ranges
    if values1 is None:
        if hasattr(param1, 'distribution'):
            dist = param1.distribution
            if hasattr(dist, 'lower') and hasattr(dist, 'upper'):
                values1 = np.linspace(dist.lower, dist.upper, 6)
    
    if values2 is None:
        if hasattr(param2, 'distribution'):
            dist = param2.distribution
            if hasattr(dist, 'lower') and hasattr(dist, 'upper'):
                values2 = np.linspace(dist.lower, dist.upper, 5)
    
    # Ensure values1 and values2 are 1D arrays
    values1 = np.ravel(values1)
    values2 = np.ravel(values2)
    
    print(f"Running bivariate sensitivity analysis: {param_name1} vs {param_name2}")
    print(f"{param_name1} range: {min(values1)} - {max(values1)}")
    print(f"{param_name2} range: {min(values2)} - {max(values2)}")
    
    # Get baseline sample
    baseline = model.get_baseline_sample()
    
    # Create result grid
    grid_data = np.zeros((len(values1), len(values2)))
    
    # Find correct metric index
    metric_index = None
    for i, metric in enumerate(model.metrics):
        if metric.name == metric_name:
            metric_index = i
            break
    
    if metric_index is None:
        raise ValueError(f"Metric '{metric_name}' not found. Available metrics: {[m.name for m in model.metrics]}")
    
    # Evaluate each parameter combination
    for i, val1 in enumerate(values1):
        for j, val2 in enumerate(values2):
            sample = baseline.copy()
            sample[param1.index] = val1
            sample[param2.index] = val2
            
            # Directly evaluate model to get series result
            result_series = model(sample)
            # Use iloc instead of direct indexing
            grid_data[i, j] = result_series.iloc[metric_index] if metric_index is not None else result_series[metric_name]
    
    # Create DataFrame for saving and plotting
    result_df = pd.DataFrame(grid_data, index=values1, columns=values2)
    
    # Save results
    if save:
        os.makedirs('results', exist_ok=True)
        safe_name1 = 'param1_' + str(hash(param_name1) % 10000)
        safe_name2 = 'param2_' + str(hash(param_name2) % 10000)
        safe_metric = 'metric_' + str(hash(metric_name) % 10000)
        result_df.to_csv(f'results/bivariate_{safe_name1}_{safe_name2}_{safe_metric}.csv')
    
    # Create colormap
    if metric_name == 'ROI':
        cmap = 'RdYlGn'
        center = 0
    else:
        cmap = 'viridis'
        center = None
    
    # Plot heatmap
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(result_df, annot=True, cmap=cmap, center=center, 
                fmt=".2f",  # Two decimal places
                annot_kws={"size": 9},  # Adjust annotation size
                    cbar_kws={'label': f'{metric_name}'})
    
    ax.set_title(f'Two-Variable Sensitivity Analysis: {metric_name}')
    ax.set_xlabel(param_name2)
    ax.set_ylabel(param_name1)
    plt.tight_layout()
    
    # Save figure
    if save:
        os.makedirs('figures', exist_ok=True)
        plt.savefig(f'figures/heatmap_{safe_name1}_{safe_name2}_{safe_metric}.png', dpi=300, bbox_inches='tight')
    
    return result_df

# Monte Carlo analysis
def run_monte_carlo(model, n_samples=1000, sample_type='LHS', save=True):
    """Run Monte Carlo uncertainty analysis using manual sample generation"""
    print(f"Running Monte Carlo simulation ({n_samples} samples)...")
    
    # Get parameter distributions
    param_distributions = {}
    for param in model.parameters:
        if hasattr(param, 'distribution'):
            param_distributions[param.index] = param.distribution
    
    # Get metric name mapping
    metric_names = {}
    for i, metric in enumerate(model.metrics):
        metric_names[i] = metric.name  # Store position and metric name mapping
    
    # Get baseline sample as template
    baseline = model.get_baseline_sample()
    results = []
    
    # Generate samples and evaluate each one
    for i in range(n_samples):
        if i % 100 == 0 and i > 0:
            print(f"Completed {i}/{n_samples} samples")
            
        # Copy baseline sample
        sample = baseline.copy()
        
        # Generate random values for each parameter
        for param_name, distribution in param_distributions.items():
            if sample_type == 'LHS':
                # Simplified Latin Hypercube Sampling
                u = (i + 0.5) / n_samples  # Ensure uniform coverage
                if hasattr(distribution, 'ppf'):
                    sample[param_name] = distribution.ppf(u)
                else:
                    # Fall back to regular random sampling
                    sample[param_name] = distribution.sample()
            else:
                # Regular random sampling
                sample[param_name] = distribution.sample()
        
        # Evaluate sample
        try:
            # Get metric values
            result_series = model(sample)
            
            # Create dictionary with parameter and metric values
            result_dict = {}
            
            # Add parameter values
            for param_name in param_distributions.keys():
                result_dict[param_name] = sample[param_name]
            
            # Add metric values using our defined names
            for i, name in metric_names.items():
                result_dict[name] = result_series.iloc[i]
            
            results.append(result_dict)
        except Exception as e:
            print(f"Sample {i} evaluation failed: {str(e)}")
    
    # Combine results
    df_results = pd.DataFrame(results)
    
    # Save results
    if save:
        os.makedirs('results', exist_ok=True)
        df_results.to_csv('results/monte_carlo_results.csv', index=False)
    
    # Plot histograms - use correct column names
    for metric_name in metric_names.values():
        if metric_name in df_results.columns:  # Confirm column exists
            plt.figure(figsize=(10, 6))
        
            # Plot histogram and kernel density estimate
        sns.histplot(df_results[metric_name], kde=True)
        
            # Add vertical lines for mean and median
        mean_val = df_results[metric_name].mean()
        median_val = df_results[metric_name].median()
        plt.axvline(mean_val, color='r', linestyle='-', 
                   label=f'Mean: {mean_val:.3f}')
        plt.axvline(median_val, color='g', linestyle='--', 
                   label=f'Median: {median_val:.3f}')
        
        plt.xlabel(metric_name)
        plt.ylabel('Frequency')
        plt.title(f'Monte Carlo Simulation: {metric_name} Distribution')
        plt.legend()
        
        if save:
            os.makedirs('figures', exist_ok=True)
            plt.savefig(f'figures/monte_carlo_{metric_name}.png', dpi=300)
    
    return df_results

# Run scenario analysis
def run_scenario_analysis(model, scenarios, save=True):
    """Run predefined scenario analysis"""
    print("Running scenario analysis...")
    
    # Get metric name mapping
    metric_names = {}
    for i, metric in enumerate(model.metrics):
        metric_names[i] = metric.name
        
    # Create result containers
    metrics = list(metric_names.values())
    results = {metric: [] for metric in metrics}
    results['Scenario'] = []
    
    # Get baseline sample as template
    baseline = model.get_baseline_sample()
    
    for scenario_name, params in scenarios.items():
        results['Scenario'].append(scenario_name)
        print(f"\nApplying scenario: {scenario_name}")
        
        # Display parameters
        for param_name, value in params.items():
            print(f"  {param_name} = {value}")
        
        # Create scenario sample
        sample = baseline.copy()
        
        # Set scenario parameter values
        for param in model.parameters:
            if param.name in params:
                sample[param.index] = params[param.name]
        
        # Directly evaluate model
        try:
            result_series = model(sample)
            
            # Add metric results
            for i, metric_name in metric_names.items():
                value = result_series.iloc[i]
                results[metric_name].append(value)
                print(f"  {metric_name}: {value}")
        except Exception as e:
            print(f"Scenario evaluation failed: {str(e)}")
            # Add missing values
            for metric_name in metric_names.values():
                results[metric_name].append(None)
    
    # Create result DataFrame
    df_results = pd.DataFrame(results)
    
    # Save results
    if save:
        os.makedirs('results', exist_ok=True)
        df_results.to_csv('results/scenario_analysis.csv', index=False)
    
    # Plot bar charts
    for metric_name in metric_names.values():
        plt.figure(figsize=(10, 6))
        plt.bar(df_results['Scenario'], df_results[metric_name])
        plt.title(f'{metric_name} Across Scenarios')
        plt.ylabel(f'{metric_name}')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        if save:
            os.makedirs('figures', exist_ok=True)
            plt.savefig(f'figures/scenario_{metric_name}.png', dpi=300, bbox_inches='tight')
    
    return df_results

def plot_contour_from_heatmap(result_df, param_name1, param_name2, metric_name, 
                             cmap='viridis', levels=15, save=True):
    """Convert heatmap data to contour plot, similar to cane project style"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Get coordinate grid
    x = result_df.columns.astype(float)
    y = result_df.index.astype(float)
    X, Y = np.meshgrid(x, y)
    Z = result_df.values
    
    # Create custom colormap
    if metric_name == 'ROI':
        # ROI uses red-to-green gradient
        cmap = LinearSegmentedColormap.from_list("ROI", ["#FF5555", "#FFFFFF", "#55FF55"])
        # Find zero point as center
        vmin, vmax = Z.min(), Z.max()
        if vmin < 0 < vmax:
            vabs = max(abs(vmin), abs(vmax))
            norm = plt.Normalize(-vabs, vabs)
        else:
            norm = plt.Normalize(vmin, vmax)
    else:
        norm = None
    
    # Plot contour
    cs = ax.contourf(X, Y, Z, levels=levels, cmap=cmap, norm=norm)
    cs2 = ax.contour(X, Y, Z, levels=levels, colors='k', linewidths=0.5, alpha=0.3)
    
    # Add labels
    ax.clabel(cs2, inline=True, fontsize=8, fmt='%.2f')
    
    # Add colorbar
    cbar = plt.colorbar(cs, ax=ax)
    cbar.set_label(metric_name)
    
    # Set labels
    ax.set_xlabel(param_name2)
    ax.set_ylabel(param_name1)
    ax.set_title(f'Contour Map: {metric_name} as function of {param_name1} and {param_name2}')
    
    # Save figure
    if save:
        os.makedirs('figures', exist_ok=True)
        safe_name1 = 'param1_' + str(hash(param_name1) % 10000)
        safe_name2 = 'param2_' + str(hash(param_name2) % 10000)
        safe_metric = 'metric_' + str(hash(metric_name) % 10000)
        plt.savefig(f'figures/contour_{safe_name1}_{safe_name2}_{safe_metric}.png', 
                   dpi=300, bbox_inches='tight')
    
    return fig, ax

def plot_improved_sensitivity(df, parameter_name, metrics=None, save=True):
    """Plot improved sensitivity curves similar to cane project style"""
    if metrics is None:
        metrics = [col for col in df.columns if col != parameter_name]
    
    # Create separate chart for each metric
    for metric in metrics:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Get data
        x = df[parameter_name]
        y = df[metric]
        
        # Plot curve with markers
        ax.plot(x, y, '-o', lw=2.5, markersize=8, color='#3366CC', 
                markerfacecolor='white', markeredgecolor='#3366CC', markeredgewidth=2)
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.6)
        
        # Set title and labels
        ax.set_title(f'Sensitivity: {metric} vs {parameter_name}', fontsize=14)
        ax.set_xlabel(parameter_name, fontsize=12)
        ax.set_ylabel(metric, fontsize=12)
        
        # Beautify axes
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Save chart
        if save:
            os.makedirs('figures', exist_ok=True)
            safe_param = 'param_' + str(hash(parameter_name) % 10000)
            safe_metric = 'metric_' + str(hash(metric) % 10000)
            plt.savefig(f'figures/sensitivity_{safe_param}_{safe_metric}_improved.png', 
                       dpi=300, bbox_inches='tight')
            
    return fig, ax

def plot_improved_monte_carlo(df, metric_name, save=True):
    """Plot enhanced Monte Carlo distribution plot similar to cane project style"""
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Separate histogram and KDE curve plotting
    # First plot histogram without KDE
    sns.histplot(df[metric_name], kde=False, bins=25, ax=ax,
                 color='#3366CC', edgecolor='white', alpha=0.7)
    
    # Then plot KDE curve separately
    sns.kdeplot(df[metric_name], ax=ax, color='red', lw=2)
    
    # Calculate and add statistics
    mean_val = df[metric_name].mean()
    median_val = df[metric_name].median()
    std_val = df[metric_name].std()
    
    # Add vertical lines
    ax.axvline(mean_val, color='red', linestyle='-', lw=2,
              label=f'Mean: {mean_val:.3f}')
    ax.axvline(median_val, color='green', linestyle='--', lw=2,
              label=f'Median: {median_val:.3f}')
    
    # Add distribution area
    p10 = df[metric_name].quantile(0.1)
    p90 = df[metric_name].quantile(0.9)
    ax.axvspan(p10, p90, alpha=0.2, color='grey', 
               label=f'10-90% CI: [{p10:.3f}, {p90:.3f}]')
    
    # Set title and labels
    ax.set_title(f'Monte Carlo Distribution: {metric_name}', fontsize=14)
    ax.set_xlabel(metric_name, fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    
    # Add legend
    ax.legend()
    
    # Save chart
    if save:
        os.makedirs('figures', exist_ok=True)
        safe_metric = 'metric_' + str(hash(metric_name) % 10000)
        plt.savefig(f'figures/monte_carlo_{safe_metric}_improved.png', 
                   dpi=300, bbox_inches='tight')
            
    return fig, ax

def plot_tornado(df, metric_name, parameters=None, n_parameters=10, save=True):
    """Plot tornado diagram showing parameter influence on specific metric"""
    # Calculate Spearman correlation coefficients
    if parameters is None:
        # Exclude metric columns, keep only parameter columns
        parameters = [col for col in df.columns if col != metric_name]
    
    # Calculate correlation for each parameter with the metric
    correlations = []
    for param in parameters:
        corr, _ = spearmanr(df[param], df[metric_name])
        correlations.append((param, corr))
    
    # Sort by correlation coefficient absolute value
    correlations.sort(key=lambda x: abs(x[1]), reverse=True)
    
    # Keep only top n_parameters
    if n_parameters < len(correlations):
        correlations = correlations[:n_parameters]
    
    # Extract data for plotting
    param_names, corr_values = zip(*correlations)
    
    # Plot tornado diagram
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Determine colors
    colors = ['#3366CC' if c > 0 else '#CC3366' for c in corr_values]
    
    # Plot horizontal bars
    bars = ax.barh(range(len(param_names)), corr_values, color=colors, 
                  height=0.6, edgecolor='white', linewidth=1)
    
    # Add parameter name labels
    ax.set_yticks(range(len(param_names)))
    ax.set_yticklabels(param_names)
    
    # Set axes
    ax.set_xlabel('Spearman Correlation Coefficient', fontsize=12)
    ax.set_title(f'Parameter Importance for {metric_name}', fontsize=14)
    
    # Add vertical line at zero
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5)
    
    # Add grid
    ax.grid(axis='x', linestyle='--', alpha=0.3)
    
    # Beautify
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add value labels for each bar
    for i, bar in enumerate(bars):
        width = bar.get_width()
        label_x_pos = width + 0.01 if width >= 0 else width - 0.08
        ax.text(label_x_pos, bar.get_y() + bar.get_height()/2, 
               f'{width:.3f}', va='center')
    
    # Save chart
    if save:
        os.makedirs('figures', exist_ok=True)
        safe_metric = 'metric_' + str(hash(metric_name) % 10000)
        plt.savefig(f'figures/tornado_{safe_metric}.png', dpi=300, bbox_inches='tight')
    
    return fig, ax

# Main program
if __name__ == "__main__":
    print("Running BioSTEAM Sensitivity Analysis...")
    
    # Initialize system and model
    bst.settings.set_thermo(chems)
    main_flowsheet.clear()
    microalgae_sys = create_microalgae_MCCA_production_sys()
    u = microalgae_sys.flowsheet.unit
    s = microalgae_sys.flowsheet.stream
    
    # Create TEA object
    dry_tpd = u.U101.ins[0].F_mass * 24 / 1000  # kg/h -> t/d
    microalgae_tea = create_tea(
        system=microalgae_sys, 
        IRR=0.10, 
        duration=(2024, 2045),
        depreciation='MACRS7', 
        income_tax=0.21, 
        operating_days=330,
        lang_factor=None, 
        construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, 
        startup_FOCfrac=1, 
        startup_salesfrac=0.5,
        startup_VOCfrac=0.75, 
        WC_over_FCI=0.05,
        finance_interest=0.08, 
        finance_years=10, 
        finance_fraction=0.4,
        OSBL_units=(u.CT, u.CWP, u.ADP, u.PWC),
        warehouse=0.04, 
        site_development=0.09, 
        additional_piping=0.045,
        proratable_costs=0.10, 
        field_expenses=0.10, 
        construction=0.20,
        contingency=0.10, 
        other_indirect_costs=0.10, 
        labor_cost=max(0.5e6, compute_labor_cost(dry_tpd)),
        labor_burden=0.90, 
        property_insurance=0.007, 
        maintenance=0.03, 
        boiler_turbogenerator=None,
        steam_power_depreciation='MACRS20'
    )

    # Create evaluation model
    model = create_evaluation_model(microalgae_sys, microalgae_tea)
    
    # Define scenarios
    scenarios = {
    'Optimistic': {
        'Microalgae price': -0.2,  # Higher tipping fee income
        'Caproic acid price': 3.5,   # Higher product price
        'Caproic acid yield factor': 1.2,  # Higher than baseline yield
        'Operating days': 350,   # More operating days
        'Maintenance factor': 0.02,  # Lower maintenance cost
    },
    'Baseline': {
        'Microalgae price': -0.1,  # Medium tipping fee income
        'Caproic acid price': 2.89,
        'Caproic acid yield factor': 1.0,
        'Operating days': 330,
        'Maintenance factor': 0.03,
    },
    'Pessimistic': {
        'Microalgae price': 0.0,   # No tipping fee income
        'Caproic acid price': 2.0,   # Lower product price
        'Caproic acid yield factor': 0.8, # Lower than baseline yield
        'Operating days': 300,   # Fewer operating days
        'Maintenance factor': 0.05, # Higher maintenance cost
    },
    'Purchase': {
        'Microalgae price': 0.1,   # Need to purchase feedstock
        'Caproic acid price': 3.2,   # Higher product price
        'Caproic acid yield factor': 1.15, # Higher yield
        'Operating days': 330,
        'Maintenance factor': 0.03,
    }
    }
    
    # Single variable sensitivity analysis
    print("\n1. Analyzing the effect of microalgae price on ROI and MSP...")
    microalgae_price_results = run_univariate_analysis(model, 'Microalgae price')
    plot_improved_sensitivity(microalgae_price_results, 'Microalgae price')
    
    # Bivariate analysis
    print("\n2. Analyzing the effect of microalgae price and caproic acid price on ROI...")
    heatmap_results = run_bivariate_analysis(model, 'Microalgae price', 'Caproic acid price', 'ROI')
    plot_contour_from_heatmap(heatmap_results, 'Microalgae price', 'Caproic acid price', 'ROI')
    
    # Monte Carlo analysis
    print("\n3. Running Monte Carlo simulation...")
    monte_carlo_results = run_monte_carlo(model, n_samples=10000) 
    plot_improved_monte_carlo(monte_carlo_results, 'ROI')
    
    # Tornado diagram
    plot_tornado(monte_carlo_results, 'ROI')
    
    # Scenario analysis
    print("\n4. Running scenario analysis...")
    scenario_results = run_scenario_analysis(model, scenarios)
    
    print("\nSensitivity analysis complete! Results saved in 'results' and 'figures' folders")