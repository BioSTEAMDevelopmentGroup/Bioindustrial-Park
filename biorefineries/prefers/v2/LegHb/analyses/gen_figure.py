# -*- coding: utf-8 -*-
"""
Figure Generation for LegHb Monte Carlo Results
===============================================

Loads Monte Carlo results and generates publication-ready figures including
breakdown stacked column plots for TEA/LCA.
"""

from warnings import filterwarnings
filterwarnings('ignore')

import sys
# Fix for Windows console encoding (cp1252) when printing unicode characters (e.g. checkmarks)
try:
    if sys.stdout.encoding != 'utf-8':
        sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import biosteam as bst

try:
    import plotly.graph_objects as go
except Exception:
    go = None

from biorefineries.prefers.v2.LegHb.system import get_available_configs
from biorefineries.prefers.v2.utils import plots, style, utils, plot_utils, diagram_utils


# =============================================================================
# Category Colors (Consistent Across All Plots)
# =============================================================================
# =============================================================================
# Category Colors (PreFerS Palette)
# =============================================================================
CATEGORY_COLORS = {
    'Conversion': style.PREFERS_COLORS[0],      # Navy
    'Concentration': style.PREFERS_COLORS[1],   # Indigo
    'Purification': style.PREFERS_COLORS[2],    # Cerulean
    'Formulation': style.PREFERS_COLORS[3],     # Teal
    'WasteTreatment': style.PREFERS_COLORS[4],  # Emerald
    'Facilities': style.PREFERS_COLORS[5],      # Leaf
}

MSP_COLORS = {
    'Conversion': style.PREFERS_COLORS[0],
    'Concentration': style.PREFERS_COLORS[1],
    'Purification': style.PREFERS_COLORS[2],
    'Formulation': style.PREFERS_COLORS[3],
    'WasteTreatment': style.PREFERS_COLORS[4],
    'Facilities': style.PREFERS_COLORS[5],
    'Labor salary': style.PREFERS_COLORS[6],    # Lime
    'Labor burden': style.PREFERS_COLORS[8],    # Orange
    'Maintenance': style.PREFERS_COLORS[7],     # Gold
    'Administration': style.adjust_lightness(style.PREFERS_COLORS[0], 1.2),
    'Property tax': style.adjust_lightness(style.PREFERS_COLORS[2], 1.2),
    'Property insurance': style.adjust_lightness(style.PREFERS_COLORS[4], 1.2),
    'Electricity': style.adjust_lightness(style.PREFERS_COLORS[1], 0.8),
    'Feedstock': '#8CD17D', # Keep specific if needed, or map to Leaf
    'Other Utilities': '#FABFD2',
}


# =============================================================================
# CLI
# =============================================================================

def parse_arguments():
    parser = argparse.ArgumentParser(description='LegHb figure generation')
    parser.add_argument('--config', type=str, default='config1',
                        choices=get_available_configs(),
                        help='Process configuration (default: config1)')
    parser.add_argument('--timestamp', type=str, default=None,
                        help='Custom timestamp (YYYYMMDD_HHMM) to locate results')
    parser.add_argument('--results-dir', type=str, default=None,
                        help='Existing results directory to load data from')
    args, _ = parser.parse_known_args()
    return args


# =============================================================================
# Utilities
# =============================================================================

def find_latest_results_dir(script_path, config):
    analyses_dir = os.path.dirname(os.path.abspath(script_path))
    prefix = f"results_{config}_" if config else "results_"
    if not os.path.isdir(analyses_dir):
        return None

    candidates = [
        os.path.join(analyses_dir, name)
        for name in os.listdir(analyses_dir)
        if name.startswith(prefix) and os.path.isdir(os.path.join(analyses_dir, name))
    ]

    if not candidates:
        return None

    candidates.sort()
    return candidates[-1]


def resolve_output_dirs(script_path, config, timestamp=None, results_dir=None):
    if results_dir:
        base_dir = os.path.abspath(results_dir)
    elif timestamp:
        dirs = utils.get_analysis_dirs(script_path, config=config, timestamp=timestamp, create=False)
        base_dir = dirs['base']
    else:
        base_dir = find_latest_results_dir(script_path, config)

    if not base_dir:
        raise FileNotFoundError('No results directory found. Provide --results-dir or --timestamp.')

    data_dir = os.path.join(base_dir, 'data')
    figure_dir = os.path.join(base_dir, 'figure')
    os.makedirs(figure_dir, exist_ok=True)
    return {
        'base': base_dir,
        'data': data_dir,
        'figure': figure_dir,
        'config': config,
    }


def load_breakdown_summary(data_dir):
    filepath = os.path.join(data_dir, 'Breakdown_Summary.xlsx')
    if not os.path.isfile(filepath):
        # Fallback: Check sibling directories (often gen_data_base runs in a different timestamp folder)
        try:
            parent_dir = os.path.dirname(os.path.dirname(data_dir)) # analyses/
            base_dir_name = os.path.basename(os.path.dirname(data_dir)) # results_...
            
            # Find all result directories
            candidates = [
                os.path.join(parent_dir, d, 'data', 'Breakdown_Summary.xlsx')
                for d in os.listdir(parent_dir)
                if d.startswith('results_') and os.path.isdir(os.path.join(parent_dir, d))
            ]
            # Sort by modification time (reverse) to find most recent
            candidates.sort(key=lambda x: os.path.getmtime(os.path.dirname(x)) if os.path.exists(os.path.dirname(x)) else 0, reverse=True)
            
            for candidate in candidates:
                if os.path.isfile(candidate):
                    print(f"  [Note] Found breakdown summary in sibling dir: {candidate}")
                    return pd.ExcelFile(candidate)
        except Exception:
            pass
        return None
    return pd.ExcelFile(filepath)


def deserialize_index(value):
    if isinstance(value, list) and len(value) == 2:
        return tuple(value)
    return value


def load_metadata(data_dir):
    metadata_file = os.path.join(data_dir, 'analysis_metadata.json')
    if not os.path.isfile(metadata_file):
        return {}
    with open(metadata_file, 'r', encoding='utf-8') as f:
        return json.load(f)


def load_column_map(data_dir):
    column_map_file = os.path.join(data_dir, 'column_map.json')
    if not os.path.isfile(column_map_file):
        return None
    with open(column_map_file, 'r', encoding='utf-8') as f:
        return json.load(f)


def load_mc_results(data_dir):
    """Load both Fixed Scale and Variable Scale MC results."""
    results = {}
    
    # helper to load a specific type
    def load_type(suffix):
        # Try pickles first
        names = [f'monte_carlo_{suffix}.pkl', f'monte_carlo_{suffix}.xlsx', f'monte_carlo_{suffix}.csv']
        for name in names:
            path = os.path.join(data_dir, name)
            if os.path.isfile(path):
                print(f"  [+] Loading {suffix} results from: {name}")
                if name.endswith('.pkl'):
                    return pd.read_pickle(path)
                elif name.endswith('.xlsx'):
                    return pd.read_excel(path, index_col=0)
                elif name.endswith('.csv'):
                    return pd.read_csv(path, index_col=0)
        return None

    results['no_scale'] = load_type('no_scale')
    results['with_scale'] = load_type('with_scale')
    results['ferm_only'] = load_type('ferm_only')
    results['dsp_only'] = load_type('dsp_only')
    results['econ_only'] = load_type('econ_only')
    
    # Fallback: try generic 'monte_carlo_results' and assign to 'no_scale' if others missing
    if results['no_scale'] is None and results['with_scale'] is None:
         results['no_scale'] = load_type('results') # generic
    
    # Apply column mapping if available
    column_map = load_column_map(data_dir)
    if column_map:
        mapping = {item['flat']: deserialize_index(item['original'])
                   for item in column_map.get('columns', [])}
        
        for key in results:
            if results[key] is not None:
                results[key].columns = [mapping.get(col, col) for col in results[key].columns]

    return results


def load_baseline_metrics(data_dir):
    baseline_json = os.path.join(data_dir, 'baseline_metrics.json')
    if not os.path.isfile(baseline_json):
        return {}
    with open(baseline_json, 'r', encoding='utf-8') as f:
        raw = json.load(f)

    baseline = {}
    for key, val in raw.items():
        if '|' in key:
            group, metric = key.split('|', 1)
            baseline[(group, metric)] = val
        else:
            baseline[key] = val
    return baseline


def load_sensitivity_table(data_dir):

    # Load sensitivity table for Tornado plots
    xlsx_path = os.path.join(data_dir, 'tornado_sensitivity.xlsx')
    if os.path.isfile(xlsx_path):
        try:
             # Try to read 'Summary' sheet first which we created
             return pd.read_excel(xlsx_path, sheet_name='Summary', index_col=0)
        except:
             # Fallback to default if sheet not found
             return pd.read_excel(xlsx_path, index_col=0)
    
    csv_path = os.path.join(data_dir, 'tornado_sensitivity.csv')
    if os.path.isfile(csv_path):
        return pd.read_csv(csv_path, index_col=0)
        
    return None



# =============================================================================
# Parameter Categorization
# =============================================================================

def get_parameter_category(param_name):
    """
    Categorize key parameters into 5-6 groups for coloring.
    """
    name = param_name.lower()
    
    if any(x in name for x in ['yield', 'titer', 'productivity', 'rate', 'tau', 'days']):
        return 'Process Performance'
    elif any(x in name for x in ['price', 'cost', 'catalyst', 'enzyme', 'feedstock', 'tax', 'irr']):
        return 'Material Prices'
    elif any(x in name for x in ['electricity', 'utility', 'steam', 'cooling']):
        return 'Utility Prices'
    elif any(x in name for x in ['capacity', 'investment', 'lang', 'capex', 'base']):
        return 'CapEx Factors'
    elif any(x in name for x in ['efficiency', 'recovery', 'split', 'capture', 'retention']):
        return 'Efficiency/Recovery'
    else:
        return 'Other'

# =============================================================================
# Spearman Rank Correlation
# =============================================================================

def plot_spearman_categorized(rho_series, metric_name, figure_dir):
    """Plot Spearman correlation with categories."""
    # Prepare data
    df = rho_series.to_frame(name='rho')
    df['Parameter'] = df.index
    df['Category'] = df['Parameter'].apply(get_parameter_category)
    df = df.sort_values('rho', key=abs, ascending=True)
    
    # Colors
    unique_cats = sorted(list(set(df['Category'])))
    cat_colors = dict(zip(unique_cats, style.PREFERS_COLORS[:len(unique_cats)]))
    if 'Other' in cat_colors and 'Other' not in unique_cats:
        cat_colors['Other'] = style.PREFERS_COLORS[-1]

    df['Color'] = df['Category'].map(cat_colors)
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, max(5, len(df) * 0.4 + 1.5)))
    
    ax.barh(np.arange(len(df)), df['rho'], color=df['Color'], alpha=0.85, edgecolor='black', linewidth=0.5)
    
    # Text annotations
    for i, (rho, name) in enumerate(zip(df['rho'], df['Parameter'])):
        ax.text(rho, i, f' {rho:.2f} ', va='center', ha='left' if rho > 0 else 'right', 
                fontsize=8, fontweight='bold', color='black')

    # Legend
    legend_elements = [plt.Rectangle((0,0),1,1, color=c, label=l) for l,c in cat_colors.items()]
    ax.legend(handles=legend_elements, loc='lower right', framealpha=0.95, title='Category')
    
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df['Parameter'])
    ax.set_xlabel(f"Spearman's ρ ({metric_name})", fontweight='bold')
    ax.set_title(f'Sensitivity Analysis: {metric_name}', fontweight='bold', pad=15)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlim(-1.1, 1.1)
    
    plt.tight_layout()
    utils.save_figure(fig, f'spearman_tornado_{metric_name}', figure_dir, formats=('png',))
    plt.close(fig)


def compute_spearman_correlations(results_df, parameter_indices, metric_indices, exclude_fixed=True):
    """Compute Spearman rank correlations between parameters and metrics."""
    # Filter out fixed parameters
    if exclude_fixed:
        variable_params = []
        for col in parameter_indices:
            # Check if column exists and has variance
            if col in results_df.columns and results_df[col].nunique() > 1:
                variable_params.append(col)
        parameter_indices = variable_params

    rho_data = {}
    
    for metric in metric_indices:
        if metric not in results_df.columns:
            continue
            
        rhos = []
        for param in parameter_indices:
            try:
                # Use scipy.stats.spearmanr
                s, p = spearmanr(results_df[param], results_df[metric], nan_policy='omit')
                if np.isnan(s): s = 0
                rhos.append(s)
            except:
                rhos.append(0)
        
        rho_data[metric] = rhos

    # DataFrame: Rows=Parameters, Cols=Metrics (initially)
    df = pd.DataFrame(rho_data, index=parameter_indices)
    
    # Return Transposed: Rows=Metrics, Cols=Parameters
    return df.T



# =============================================================================
# Breakdown Plots
# =============================================================================

def plot_stacked_bar_sheets_1to5(excel_file, figure_dir):
    """
    Plot 1: Stacked column for Sheets 1-5, X=Sheet, Y=%, stacked by Category.
    Allows negative values.
    """
    sheets = ['Total Capital', 'Annual Material', 'Heating Duty', 'Cooling Duty', 'Electricity Duty']
    
    data = {}
    for sheet in sheets:
        try:
            df = pd.read_excel(excel_file, sheet_name=sheet)
            # Find the % column and Category column
            cat_col = 'Category' if 'Category' in df.columns else 'Name'
            pct_col = '%'
            if pct_col in df.columns and cat_col in df.columns:
                data[sheet] = dict(zip(df[cat_col], df[pct_col]))
        except Exception:
            pass
    
    if not data:
        print("  [!] No data for stacked bar sheets 1-5")
        return
    
    # Build matrix
    # Collect all unique categories from data
    all_categories = set()
    for sheet_data in data.values():
        all_categories.update(sheet_data.keys())
    categories = sorted(list(all_categories)) # Sort for consistency
    
    sheet_names = list(data.keys())
    
    # Generate dynamic palette
    colors = style.get_palette(len(categories))
    cat_colors = dict(zip(categories, colors))

    matrix = np.zeros((len(categories), len(sheet_names)))
    for j, sheet in enumerate(sheet_names):
        sheet_data = data[sheet]
        for i, cat in enumerate(categories):
            matrix[i, j] = sheet_data.get(cat, 0.0)
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(sheet_names))
    width = 0.6
    
    bottom_pos = np.zeros(len(sheet_names))
    bottom_neg = np.zeros(len(sheet_names))
    
    for i, cat in enumerate(categories):
        vals = matrix[i, :]
        pos_vals = np.where(vals >= 0, vals, 0)
        neg_vals = np.where(vals < 0, vals, 0)
        
        c = cat_colors[cat]
        ax.bar(x, pos_vals, width, bottom=bottom_pos, label=cat, color=c)
        ax.bar(x, neg_vals, width, bottom=bottom_neg, color=c)
        
        bottom_pos += pos_vals
        bottom_neg += neg_vals
    
    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace(' ', '\n') for s in sheet_names], fontsize=9)
    ax.set_ylabel('Contribution (%)')
    ax.set_title('TEA Breakdown by Process Area')
    
    # Shrink current axis by 20% to fit legend outside
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8, title='Category')
    
    plt.tight_layout()
    
    filepath = os.path.join(figure_dir, 'breakdown_sheets_1to5.png')
    fig.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  [+] Saved: {filepath}")


def plot_stacked_bar_msp(excel_file, figure_dir):
    """
    Plot 2: MSP Breakdown stacked column (single bar with categories).
    """
    try:
        df = pd.read_excel(excel_file, sheet_name='MSP Breakdown')
    except Exception:
        print("  [!] MSP Breakdown sheet not found")
        return
    
    if df.empty:
        return
    
    fig, ax = plt.subplots(figsize=(6, 8))
    
    names = df['Name'].tolist()
    values = df['Cost (USD/kg product)'].tolist()
    
    bottom_pos = 0
    bottom_neg = 0
    
    # Generate dynamic palette for MSP items
    # Interpolate if many items to create spectrum effect
    colors = style.get_palette(len(names))
    
    for i, (name, val) in enumerate(zip(names, values)):
        color = colors[i]
        if val >= 0:
            ax.bar('MSP', val, bottom=bottom_pos, label=name, color=color)
            bottom_pos += val
        else:
            ax.bar('MSP', val, bottom=bottom_neg, label=name, color=color)
            bottom_neg += val
    
    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_ylabel('Cost (USD/kg product)')
    ax.set_title('MSP Cost Breakdown')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)
    
    plt.tight_layout()
    
    filepath = os.path.join(figure_dir, 'breakdown_msp.png')
    fig.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  [+] Saved: {filepath}")


def plot_stacked_bar_lca(excel_file, figure_dir):
    """
    Plot 3: LCA Breakdown stacked column.
    """
    try:
        df = pd.read_excel(excel_file, sheet_name='LCA Allocation')
    except Exception:
        print("  [!] LCA Allocation sheet not found")
        return
    
    if df.empty:
        return
    
    # Filter out Total row for stacking
    df_plot = df[~df['Name'].str.contains('Total', case=False, na=False)]
    
    fig, ax = plt.subplots(figsize=(6, 8))
    
    names = df_plot['Name'].tolist()
    values = df_plot['GWP (kg CO2e/kg product)'].tolist()
    
    bottom_pos = 0
    bottom_neg = 0
    
    bottom_pos = 0
    bottom_neg = 0
    
    # Color palette (PreFerS)
    colors = style.get_palette(len(names))
    
    for i, (name, val) in enumerate(zip(names, values)):
        color = colors[i % len(colors)]
        if val >= 0:
            ax.bar('GWP', val, bottom=bottom_pos, label=name[:30], color=color)
            bottom_pos += val
        else:
            ax.bar('GWP', val, bottom=bottom_neg, label=name[:30], color=color)
            bottom_neg += val
    
    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_ylabel('GWP (kg CO2e/kg product)')
    ax.set_title('LCA Breakdown by Source')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=7)
    
    plt.tight_layout()
    
    filepath = os.path.join(figure_dir, 'breakdown_lca.png')
    fig.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  [+] Saved: {filepath}")


# =============================================================================
# Scale Effects Plot
# =============================================================================

def plot_scale_effects(scale_values, metric_values, baseline_val,
                       xlabel="Production Scale [kg/hr]", ylabel="Metric",
                       title="Scale Effects", color=None, figsize=(9, 6), n_bins=20):
    """
    Generate percentile band plot showing metric vs scale.
    """
    if color is None:
        color = '#EDC948' # Gold

    # Clean data
    scale_arr = np.asarray(scale_values)
    metric_arr = np.asarray(metric_values)
    mask = np.isfinite(scale_arr) & np.isfinite(metric_arr)
    scale_clean = scale_arr[mask]
    metric_clean = metric_arr[mask]
    
    if len(scale_clean) < 10:
        return None, None

    # Bin the data
    bins = np.linspace(scale_clean.min(), scale_clean.max(), n_bins + 1)
    bin_indices = np.digitize(scale_clean, bins) - 1
    
    centers, p05, p25, p50, p75, p95 = [], [], [], [], [], []
    
    for i in range(n_bins):
        in_bin = metric_clean[bin_indices == i]
        if len(in_bin) >= 3:
            centers.append((bins[i] + bins[i+1]) / 2)
            p05.append(np.percentile(in_bin, 5))
            p25.append(np.percentile(in_bin, 25))
            p50.append(np.percentile(in_bin, 50))
            p75.append(np.percentile(in_bin, 75))
            p95.append(np.percentile(in_bin, 95))
    
    if not centers:
        return None, None

    centers = np.array(centers)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # 25-75 percentile band
    ax.fill_between(centers, p25, p75, color=color, alpha=0.3, 
                    label='25th-75th percentile')
    
    # 5-95 boundaries
    ax.plot(centers, p05, '--', color=color, linewidth=1.5, alpha=0.7,
            label='5th percentile')
    ax.plot(centers, p95, '--', color=color, linewidth=1.5, alpha=0.7,
            label='95th percentile')
    
    # Median
    ax.plot(centers, p50, '-', color=color, linewidth=2.5, alpha=0.9,
            label='Median')
    
    # Baseline
    if baseline_val is not None:
        ax.axhline(baseline_val, color='red', 
                   linestyle='-.', linewidth=2,
                   label=f'Baseline ({baseline_val:.4g})')
    
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    ax.set_title(title, fontweight='bold', pad=15)
    
    ax.legend(loc='best', framealpha=0.95)
    plt.tight_layout()
    return fig, ax


# =============================================================================
# Scatter with Color-coded Metric
# =============================================================================

def plot_colored_scatter(x_data, y_data, z_data, xlabel, ylabel, zlabel,
                         title="", cmap='PreFerS', figsize=(10, 8), reverse_cmap=False):
    """
    Generate scatter plot with color-coded third variable.
    """
    # Clean data
    x_arr = np.asarray(x_data)
    y_arr = np.asarray(y_data)
    z_arr = np.asarray(z_data)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr) & np.isfinite(z_arr)
    
    x_clean = x_arr[mask]
    y_clean = y_arr[mask]
    z_clean = z_arr[mask]

    if len(x_clean) < 10:
        return None, None
    
    cmap_val = cmap + '_r' if reverse_cmap else cmap

    fig, ax = plt.subplots(figsize=figsize)
    
    scatter = ax.scatter(x_clean, y_clean, c=z_clean,
                         s=30, cmap=cmap_val, alpha=0.7,
                         edgecolors='black', linewidths=0.3)
    
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label(zlabel, rotation=270, labelpad=20, fontweight='bold')
    
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    ax.set_title(title, fontweight='bold', pad=15)
    
    # Stats box
    stats_text = (
        f'Statistics (n={len(x_clean):,}):\n'
        f'{zlabel.split("[")[0].strip()}:\n'
        f'  Mean: {z_clean.mean():.4f}\n'
        f'  Median: {np.median(z_clean):.4f}\n'
        f'  Range: [{z_clean.min():.4f}, {z_clean.max():.4f}]'
    )
    ax.text(
        0.02, 0.98, stats_text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'),
        family='monospace'
    )
    
    plt.tight_layout()
    return fig, ax


# =============================================================================
# Composition Box Plots
# =============================================================================

def plot_composition_boxplots(results_df, figure_dir, product_type='LegHb'):
    """
    Create combined box plots for product composition metrics.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Monte Carlo results containing composition metric columns
    figure_dir : str
        Directory to save figure
    product_type : str
        'LegHb' or 'HemDx' to select appropriate specs
    """
    # Define specifications for LegHb
    leghb_specs = {
        'Fat Content': {'target': (0, 2), 'unit': 'wt%'},
        'Carbohydrates': {'target': (0, 4), 'unit': 'wt%'},
        'Product Content': {'target': (6, 9), 'unit': 'wt%'},
        'Total Solids': {'target': (0, 24), 'unit': 'wt%'},
        'Protein Purity': {'target': (65, 100), 'unit': '%'},
    }
    
    # Find composition metric columns in results
    comp_cols = []
    specs_to_plot = leghb_specs if product_type == 'LegHb' else {}
    
    for col in results_df.columns:
        col_str = str(col)
        for spec_name in specs_to_plot.keys():
            if spec_name.lower() in col_str.lower():
                comp_cols.append((col, spec_name))
                break
    
    if len(comp_cols) == 0:
        print("  [!] No composition columns found in results")
        return None
    
    n_plots = min(5, len(comp_cols))
    fig, axes = plt.subplots(1, n_plots, figsize=(4 * n_plots, 5))
    if n_plots == 1:
        axes = [axes]
    
    colors = style.get_palette(n_plots)
    
    for i, (col, spec_name) in enumerate(comp_cols[:n_plots]):
        ax = axes[i]
        data = results_df[col].dropna()
        
        if len(data) == 0:
            continue
        
        # Box plot
        bp = ax.boxplot(data, patch_artist=True, widths=0.6)
        bp['boxes'][0].set_facecolor(colors[i])
        bp['boxes'][0].set_alpha(0.7)
        
        # Target range as grey band
        spec = specs_to_plot.get(spec_name, {})
        target = spec.get('target')
        if target:
            ax.axhspan(target[0], target[1], alpha=0.2, color='grey', 
                      label='Target Range', zorder=0)
            # Target center line
            center = (target[0] + target[1]) / 2
            ax.axhline(center, color='red', linestyle='--', linewidth=1.5, 
                      alpha=0.8, label='Target Center')
        
        ax.set_title(spec_name, fontweight='bold', fontsize=10)
        ax.set_ylabel(spec.get('unit', '%'))
        ax.set_xticklabels([''])
        
        # Stats annotation
        stats_text = f'μ={data.mean():.2f}\nσ={data.std():.2f}'
        ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
                fontsize=8, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.suptitle(f'{product_type} Product Composition Distribution', 
                fontweight='bold', fontsize=12)
    plt.tight_layout()
    
    filepath = os.path.join(figure_dir, 'composition_boxplots.png')
    fig.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  [+] Saved: {filepath}")
    return fig


# =============================================================================
# Main Generation
# =============================================================================

def generate_figures(config='config1', timestamp=None, results_dir=None):
    style.set_style() # Ensure style is set
    dirs = resolve_output_dirs(__file__, config, timestamp=timestamp, results_dir=results_dir)

    print("=" * 80)
    print("LEGHEMOGLOBIN - FIGURE GENERATION")
    print("=" * 80)
    print(f"\nData:   {dirs['data']}")
    print(f"Figure: {dirs['figure']}")

    # Load data
    metadata = load_metadata(dirs['data'])
    parameter_indices = [deserialize_index(i) for i in metadata.get('parameter_indices', [])]
    metric_indices = [deserialize_index(i) for i in metadata.get('metric_indices', [])]

    # Load MC results (returns dict with 'no_scale' and 'with_scale')
    all_results = load_mc_results(dirs['data'])
    results_no_scale = all_results.get('no_scale')
    results_with_scale = all_results.get('with_scale')

    # If only one type found and assigned to wrong key (fallback logic in load_mc_results might need check)
    if results_no_scale is None and results_with_scale is not None:
        print("  [Note] Using 'with_scale' data for 'no_scale' analyses as fallback.")
        results_no_scale = results_with_scale

    baseline_metrics = load_baseline_metrics(dirs['data'])
    sensitivity_df = load_sensitivity_table(dirs['data'])

    # Helper indices
    msp_idx = ('PreFerS', 'MSP [$/kg]')
    gwp_idx = ('PreFerS', 'GWP [kg CO2-eq/kg]')
    baseline_msp = metadata.get('baseline_msp', baseline_metrics.get(msp_idx, np.nan))
    baseline_gwp = metadata.get('baseline_gwp', baseline_metrics.get(gwp_idx, np.nan))

    # Load breakdown summary
    excel_file = load_breakdown_summary(dirs['data'])
    
    # -------------------------------------------------------------------------
    # STEP 1: Spearman Rank Correlation Analysis (1D and Heatmap)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 1: Spearman Rank Correlation Analysis")
    print("=" * 80)

    if results_no_scale is not None:
        # Calculate correlations
        if parameter_indices and metric_indices:
            rho_df = compute_spearman_correlations(results_no_scale, parameter_indices, metric_indices, exclude_fixed=True)
        else:
             # Fallback
            rho_df = compute_spearman_correlations(
                results_no_scale,
                parameter_indices=[c for c in results_no_scale.columns if c not in metric_indices],
                metric_indices=metric_indices or [c for c in results_no_scale.columns if isinstance(c, tuple)],
                exclude_fixed=True
            )
        
        # Heatmap
        rho_display = rho_df.copy()
        rho_display.index = [f"{idx[0]}: {idx[1]}" if isinstance(idx, tuple) else str(idx) for idx in rho_display.index]
        rho_display.columns = [f"{col[1]}" if isinstance(col, tuple) else str(col) for col in rho_display.columns]

        fig1, ax1 = plots.plot_spearman_heatmap(
            rho_display.T, 
            figsize=(12, 6),
            cmap='PreFerS' # Match scatter plots (Spectrum)
        )
        ax1.set_title('Spearman Rank Correlation: Parameters vs Metrics', fontweight='bold')
        utils.save_figure(fig1, 'spearman_heatmap', dirs['figure'], formats=('png', 'pdf'))

        corr_file = os.path.join(dirs['data'], 'spearman_correlations.xlsx')
        rho_display.T.to_excel(corr_file)
        print(f"  [+] Saved correlation table to: {corr_file}")

        # 1D Spearman Plots (Tornado style) for separate metrics
        print("  Generating 1D Spearman plots...")
        for metric_name, metric_idx, color in [
            ('MSP', msp_idx, style.get_color('msp')),
            ('GWP', gwp_idx, style.get_color('gwp')),
        ]:
            if metric_idx in rho_df.index:
                series = rho_df.loc[metric_idx].dropna()
                if not series.empty:
                    # Map indices to readable names
                    series.index = [f"{idx[0]}: {idx[1]}" if isinstance(idx, tuple) else str(idx) for idx in series.index]
                    
                    try:
                        # Use custom categorized plot instead of bst default
                        plot_spearman_categorized(series, metric_name, dirs['figure'])

                    except Exception as e:
                        print(f"    [!] Failed to plot 1D Spearman for {metric_name}: {e}")
            else:
                 print(f"    [!] Metric {metric_idx} not found in correlation matrix")
    else:
        print("  [!] MC Results not found, skipping Spearman Analysis")

    # -------------------------------------------------------------------------
    # STEP 2: Joint KDE Plots (5 Scenarios)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 2: Joint KDE Plots (5 Scenarios)")
    print("=" * 80)
    
    scenarios = [
        ('no_scale', 'Fixed Scale', 'PreFerS_blue'),
        ('with_scale', 'Variable Scale', 'PreFerS_orange'),
        ('ferm_only', 'Fermentation Params Only', 'PreFerS_purple'),
        ('dsp_only', 'DSP Params Only', 'PreFerS_green'),
        ('econ_only', 'Econ & Env Params Only', 'PreFerS_red'),
    ]
    
    comparative_data = []

    for key, title, color_name in scenarios:
        data = all_results.get(key)
        if data is not None and msp_idx in data.columns and gwp_idx in data.columns:
             print(f"  Generating Joint Plot ({title})...")
             
             g = plots.plot_joint_marginal(
                data=data, x_col=gwp_idx, y_col=msp_idx,
                x_label='GWP [kg CO2-eq/kg]', y_label='MSP [$/kg]',
                kde_color=color_name,
                x_box_color=style.get_color(color_name.split('_')[-1]), # approximate mapping
                y_box_color=style.get_color(color_name.split('_')[-1]),
                baseline_point=(baseline_gwp, baseline_msp),
                median_point=(data[gwp_idx].median(), data[msp_idx].median())
             )
             g.figure.suptitle(f'GWP vs MSP ({title})', fontweight='bold', y=1.02)
             utils.save_figure(g.figure, f'joint_kde_{key}', dirs['figure'])
             plt.close(g.figure)
             
             # Collect for comparative (if needed)
             comparative_data.append({
                 'name': title,
                 'msp_mean': data[msp_idx].mean(),
                 'msp_std': data[msp_idx].std(),
                 'gwp_mean': data[gwp_idx].mean(),
                 'gwp_std': data[gwp_idx].std()
             })

    # -------------------------------------------------------------------------
    # STEP 3: Tornado Diagram (Single Point Sensitivity)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 3: Tornado Diagram (MSP Sensitivity)")
    print("=" * 80)

    if sensitivity_df is not None:
        # sensitivity_df = sensitivity_df.dropna() # Don't dropna blindly, it might drop rows with some NaNs in unrelated columns
        if not sensitivity_df.empty:
            # Create categories map
            categories = {p: get_parameter_category(p) for p in sensitivity_df['Parameter']}
            
            # ---------------------------
            # Plot MSP
            # ---------------------------
            # Map columns to Low/High for plot_tornado if they exist
            msplow = 'MSP [$/kg] Low'
            msphigh = 'MSP [$/kg] High'
            
            if msplow in sensitivity_df.columns:
                sens_msp = sensitivity_df.copy()
                sens_msp['Low'] = sensitivity_df[msplow]
                sens_msp['High'] = sensitivity_df[msphigh]
                sens_msp = sens_msp.dropna(subset=['Low', 'High'])
                
                if not sens_msp.empty:
                    fig3, ax3 = plots.plot_tornado(
                        sensitivity_df=sens_msp,
                        baseline=baseline_msp,
                        metric_name='MSP [$/kg]',
                        categories=categories
                    )
                    utils.save_figure(fig3, 'tornado_msp', dirs['figure'], formats=('png', 'pdf'))
                    plt.close(fig3)
                    print(f"  [+] Saved MSP tornado diagram")

            # ---------------------------
            # Plot GWP
            # ---------------------------
            gwplow = 'GWP [kg CO2-eq/kg] Low'
            gwphigh = 'GWP [kg CO2-eq/kg] High'
            
            if gwplow in sensitivity_df.columns:
                 sens_gwp = sensitivity_df.copy()
                 sens_gwp['Low'] = sensitivity_df[gwplow]
                 sens_gwp['High'] = sensitivity_df[gwphigh]
                 sens_gwp = sens_gwp.dropna(subset=['Low', 'High'])
                 
                 if not sens_gwp.empty:
                     categories_gwp = {p: get_parameter_category(p) for p in sens_gwp['Parameter']}
                     
                     fig_gwp, ax_gwp = plots.plot_tornado(
                         sensitivity_df=sens_gwp,
                         baseline=baseline_gwp,
                         metric_name='GWP [kg CO2-eq/kg]',
                         categories=categories_gwp
                     )
                     utils.save_figure(fig_gwp, 'tornado_gwp', dirs['figure'], formats=('png', 'pdf'))
                     plt.close(fig_gwp)
                     print(f"  [+] Saved GWP tornado diagram")
        else:
             print("  [!] Sensitivity df empty after dropna")
    else:
        print("  [!] No sensitivity table found")

    # -------------------------------------------------------------------------
    # STEP 4: Distribution Plots
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 4: Distribution Plots")
    print("=" * 80)

    data_for_dist = results_with_scale if results_with_scale is not None else results_no_scale
    if data_for_dist is not None:
        if msp_idx in data_for_dist.columns:
            fig4, ax4 = plots.plot_monte_carlo_dist(
                results=data_for_dist,
                metric_col=msp_idx,
                units='$/kg',
                baseline=baseline_msp,
                color=style.get_color('msp')
            )
            ax4.set_title('MSP Uncertainty Distribution', fontweight='bold')
            utils.save_figure(fig4, 'distribution_msp', dirs['figure'], formats=('png',))

        if gwp_idx in data_for_dist.columns:
            fig5, ax5 = plots.plot_monte_carlo_dist(
                results=data_for_dist,
                metric_col=gwp_idx,
                units='kg CO2-eq/kg',
                baseline=baseline_gwp,
                color=style.get_color('gwp')
            )
            ax5.set_title('GWP Uncertainty Distribution', fontweight='bold')
            utils.save_figure(fig5, 'distribution_gwp', dirs['figure'], formats=('png',))

    # -------------------------------------------------------------------------
    # STEP 5: Comprehensive Parameter Sensitivity Bands (Including Scale)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 5: Comprehensive Parameter Sensitivity Bands")
    print("=" * 80)
    
    # We want to plot Metric vs Parameter for ALL variable parameters.
    # We will use 'with_scale' if available, as authorized.
    
    data_for_sen = results_with_scale if results_with_scale is not None else results_no_scale
    
    if data_for_sen is not None:
        # Define all target metrics (MSP, GWP, and Composition)
        leghb_composition_metrics = [
            ('Composition', 'Fat Content [wt%]'), ('Composition', 'Fat content [wt%]'),
            ('Composition', 'Carbohydrates [wt%]'),
            ('Composition', 'Product Content [wt%]'), ('Composition', 'Product content [wt%]'),
            ('Composition', 'Total Solids [wt%]'), ('Composition', 'Total solids [wt%]'),
            ('Composition', 'Protein Purity [%]'), ('Composition', 'Protein purity [%]'),
        ]
        
        # Find actual composition columns in data
        comp_cols_in_data = {}
        for col in data_for_sen.columns:
            col_str = str(col).lower()
            if 'fat' in col_str: comp_cols_in_data['Fat'] = col
            elif 'carbo' in col_str: comp_cols_in_data['Carbohydrates'] = col
            elif 'product content' in col_str or 'product_content' in col_str: comp_cols_in_data['Product'] = col
            elif 'total solid' in col_str or 'total_solid' in col_str: comp_cols_in_data['Total Solids'] = col
            elif 'protein purity' in col_str or 'protein_purity' in col_str: comp_cols_in_data['Protein Purity'] = col
        
        # All metrics to plot (MSP, GWP, + found composition)
        all_metrics = [(msp_idx, 'MSP', '$/kg', baseline_msp), (gwp_idx, 'GWP', 'kg CO2-eq/kg', baseline_gwp)]
        for name, col in comp_cols_in_data.items():
            all_metrics.append((col, name, 'wt%', None))
        
        print(f"  Metrics to analyze: {[m[1] for m in all_metrics]}")
        
        # Identify parameter columns (exclude metrics)
        metric_col_set = {msp_idx, gwp_idx} | set(comp_cols_in_data.values())
        
        param_cols = []
        for col in data_for_sen.columns:
            if col in metric_col_set:
                continue
            # Check variance (ensure it actually varies)
            if data_for_sen[col].nunique() > 5:
                param_cols.append(col)
        
        print(f"  Found {len(param_cols)} varying parameters to plot.")
        
        for param_col in param_cols:
            param_name = param_col[1] if isinstance(param_col, tuple) else str(param_col)
            safe_name = "".join([c if c.isalnum() else "_" for c in param_name])
            
            category = get_parameter_category(param_name)
            band_color = style.PREFERS_COLORS[0]
            if 'Process' in category: band_color = style.PREFERS_COLORS[1]
            elif 'Price' in category: band_color = style.PREFERS_COLORS[4]
            elif 'CapEx' in category or 'Scale' in param_name: band_color = style.PREFERS_COLORS[7]
            
            for metric_idx, metric_name, metric_unit, baseline_val in all_metrics:
                if metric_idx not in data_for_sen.columns:
                    continue
                try:
                    fig_sen, ax_sen = plot_scale_effects(
                        data_for_sen[param_col], data_for_sen[metric_idx], baseline_val,
                        xlabel=param_name, ylabel=f'{metric_name} [{metric_unit}]',
                        title=f'{metric_name} Sensitivity to {param_name}',
                        color=band_color
                    )
                    if fig_sen:
                        utils.save_figure(fig_sen, f'sensitivity_{metric_name.lower()}_{safe_name}', dirs['figure'], formats=('png',))
                        plt.close(fig_sen)
                except Exception as e:
                    print(f"    [!] Failed plot {metric_name} vs {param_name}: {e}")

    else:
        print("  [!] No MC data available for sensitivity bands")

    # -------------------------------------------------------------------------
    # STEP 6: Contour Scatter Plots (All Fermentation Params vs MSP/GWP)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 6: Contour Scatter Plots (Fermentation Parameters)")
    print("=" * 80)

    # Use no_scale data to isolate biological effects
    data_for_scatter = results_no_scale
    if data_for_scatter is not None:
        # Find ALL fermentation parameter columns
        ferm_params = {}
        for col in data_for_scatter.columns:
            col_str = str(col).lower()
            if 'titer' in col_str: ferm_params['Titer'] = col
            elif 'product yield' in col_str or ('yield' in col_str and 'biomass' not in col_str): 
                if 'Yield_P' not in ferm_params: ferm_params['Yield_P'] = col
            elif 'biomass yield' in col_str or 'biomass_yield' in col_str: ferm_params['Yield_B'] = col
            elif 'tau' in col_str: ferm_params['Tau'] = col
            elif 'productivity' in col_str: ferm_params['Productivity'] = col
        
        print(f"  Found fermentation parameters: {list(ferm_params.keys())}")
        
        # Generate all pairwise combinations
        from itertools import combinations
        param_pairs = list(combinations(ferm_params.keys(), 2))
        
        # Target metrics for coloring
        target_metrics = [
            (msp_idx, 'MSP', '$/kg'),
            (gwp_idx, 'GWP', 'kg CO2-eq/kg'),
        ]
        
        for metric_idx, metric_name, metric_unit in target_metrics:
            if metric_idx not in data_for_scatter.columns:
                continue
            
            for p1_name, p2_name in param_pairs:
                p1_col = ferm_params[p1_name]
                p2_col = ferm_params[p2_name]
                
                if p1_col not in data_for_scatter.columns or p2_col not in data_for_scatter.columns:
                    continue
                
                try:
                    xlabel = p1_col[1] if isinstance(p1_col, tuple) else str(p1_col)
                    ylabel = p2_col[1] if isinstance(p2_col, tuple) else str(p2_col)
                    
                    fig_s, ax_s = plot_colored_scatter(
                        data_for_scatter[p1_col], data_for_scatter[p2_col], data_for_scatter[metric_idx],
                        xlabel=xlabel, ylabel=ylabel, zlabel=f'{metric_name} [{metric_unit}]',
                        title=f'{metric_name} vs {p1_name} & {p2_name}',
                        cmap='PreFerS', figsize=(9, 7)
                    )
                    if fig_s:
                        safe_title = f"scatter_{metric_name}_{p1_name}_{p2_name}"
                        utils.save_figure(fig_s, safe_title, dirs['figure'], formats=('png',))
                        plt.close(fig_s)
                except Exception as e:
                    print(f"    [!] Failed scatter {metric_name} vs {p1_name}/{p2_name}: {e}")
    else:
        print("  [!] No MC data available for scatter plots")
    
    # -------------------------------------------------------------------------
    # STEP 7: Breakdown Stacked Bar Charts (Existing)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 7: Breakdown Stacked Bar Charts")
    print("=" * 80)
    
    if excel_file:
        plot_stacked_bar_sheets_1to5(excel_file, dirs['figure'])
        plot_stacked_bar_msp(excel_file, dirs['figure'])
        plot_stacked_bar_lca(excel_file, dirs['figure'])
    else:
        print("  [!] Breakdown_Summary.xlsx not found")

    # -------------------------------------------------------------------------
    # STEP 8: Composition Box Plots
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 8: Composition Box Plots")
    print("=" * 80)
    
    # Use the no_scale results for composition analysis
    comp_data = results_no_scale if results_no_scale is not None else results_with_scale
    if comp_data is not None:
        try:
            plot_composition_boxplots(comp_data, dirs['figure'], product_type='LegHb')
        except Exception as e:
            print(f"  [!] Failed to generate composition box plots: {e}")
    else:
        print("  [!] No MC data available for composition plots")

    print("\n" + "=" * 80)
    print("FIGURE GENERATION COMPLETE")
    print("=" * 80)

    return dirs


def main():
    args = parse_arguments()
    generate_figures(
        config=args.config,
        timestamp=args.timestamp,
        results_dir=args.results_dir,
    )


if __name__ == '__main__':
    main()
