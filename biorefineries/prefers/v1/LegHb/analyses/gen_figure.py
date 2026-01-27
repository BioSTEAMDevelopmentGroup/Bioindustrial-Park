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

from biorefineries.prefers.v1.LegHb.system import get_available_configs
from biorefineries.prefers.v1.utils import plots, style, utils, plot_utils, diagram_utils


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
    xlsx_path = os.path.join(data_dir, 'tornado_sensitivity.xlsx')
    csv_path = os.path.join(data_dir, 'tornado_sensitivity.csv')
    if os.path.isfile(xlsx_path):
        return pd.read_excel(xlsx_path)
    if os.path.isfile(csv_path):
        return pd.read_csv(csv_path)
    return None



# =============================================================================
# Parameter Categorization
# =============================================================================

def get_parameter_category(param_name):
    """
    Categorize key parameters into 5-6 groups for coloring.
    """
    name = param_name.lower()
    
    if any(x in name for x in ['yield', 'titer', 'productivity', 'rate']):
        return 'Process Performance'
    elif any(x in name for x in ['price', 'cost', 'catalyst', 'enzyme', 'feedstock']):
        return 'Material Prices'
    elif any(x in name for x in ['electricity', 'utility', 'steam', 'cooling']):
        return 'Utility Prices'
    elif any(x in name for x in ['capacity', 'investment', 'lang', 'capex', 'base']):
        return 'CapEx Factors'
    elif any(x in name for x in ['efficiency', 'recovery', 'split']):
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
                         title="", cmap='viridis', figsize=(10, 8), reverse_cmap=False):
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
            cmap='viridis' # Match scatter plots (Spectrum)
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
    # STEP 2: Joint Marginal Plots (Fixed, Variable, Overlap)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 2: Joint KDE Plots (Cost vs. Carbon)")
    print("=" * 80)

    # 1. Fixed Scale
    if results_no_scale is not None and msp_idx in results_no_scale.columns and gwp_idx in results_no_scale.columns:
        print("  Generating Joint Plot (Fixed Scale)...")
        g = plots.plot_joint_marginal(
            data=results_no_scale, x_col=gwp_idx, y_col=msp_idx,
            x_label='GWP [kg CO2-eq/kg]', y_label='MSP [$/kg]',
            kde_color='PreFerS_blue',
            x_box_color=style.get_color('cerulean'),
            y_box_color=style.get_color('navy'),
            baseline_point=(baseline_gwp, baseline_msp),
            median_point=(results_no_scale[gwp_idx].median(), results_no_scale[msp_idx].median())
        )
        g.figure.suptitle('GWP vs MSP (Fixed Scale)', fontweight='bold', y=1.02)
        utils.save_figure(g.figure, 'kde_joint_fixed_scale', dirs['figure'], formats=('png',))
        plt.close('all')

    # 2. Variable Scale
    if results_with_scale is not None and msp_idx in results_with_scale.columns and gwp_idx in results_with_scale.columns:
        print("  Generating Joint Plot (Variable Scale)...")
        g = plots.plot_joint_marginal(
            data=results_with_scale, x_col=gwp_idx, y_col=msp_idx,
            x_label='GWP [kg CO2-eq/kg]', y_label='MSP [$/kg]',
            kde_color='PreFerS_orange',
            x_box_color=style.get_color('orange'),
            y_box_color=style.get_color('gold'),
            baseline_point=(baseline_gwp, baseline_msp),
            median_point=(results_with_scale[gwp_idx].median(), results_with_scale[msp_idx].median())
        )
        g.figure.suptitle('GWP vs MSP (Variable Scale)', fontweight='bold', y=1.02)
        utils.save_figure(g.figure, 'kde_joint_variable_scale', dirs['figure'], formats=('png',))
        plt.close('all')

    # 3. Overlap / Comparison
    if (results_no_scale is not None and results_with_scale is not None and 
        msp_idx in results_no_scale.columns and gwp_idx in results_no_scale.columns and
        msp_idx in results_with_scale.columns and gwp_idx in results_with_scale.columns):
        
        print("  Generating Comparison KDE Plot...")
        try:
            fig_comp, ax_comp = plt.subplots(figsize=(10, 8))
            
            # No Scale (Blue theme)
            sns.kdeplot(
                x=results_no_scale[gwp_idx], y=results_no_scale[msp_idx],
                ax=ax_comp, fill=False, levels=5, color=style.get_color('cerulean'), 
                linewidths=2, alpha=0.8, label='Fixed Scale'
            )
            
            # With Scale (Orange theme)
            sns.kdeplot(
                x=results_with_scale[gwp_idx], y=results_with_scale[msp_idx],
                ax=ax_comp, fill=False, levels=5, color=style.get_color('orange'), 
                linewidths=2, alpha=0.8, linestyles='--', label='Variable Scale'
            )
            
            ax_comp.set_xlabel('GWP [kg CO2-eq/kg]', fontweight='bold')
            ax_comp.set_ylabel('MSP [$/kg]', fontweight='bold')
            ax_comp.set_title('Impact of Scale on GWP-MSP Trade-off', fontweight='bold', pad=15)
            
            # Add annotations manually (Baseline + 2 Medians)
            # Baseline
            ax_comp.scatter(baseline_gwp, baseline_msp,
                           c=style.get_color('baseline'), s=250, marker='*',
                           edgecolors='black', linewidths=1.5,
                           label='Baseline', zorder=10)
            
            # Median No Scale (Blue diamond)
            ax_comp.scatter(results_no_scale[gwp_idx].median(), results_no_scale[msp_idx].median(),
                           c=style.get_color('navy'), s=150, marker='D',
                           edgecolors='black', linewidths=1.5,
                           label='Median (Fixed)', zorder=10)

            # Median With Scale (Orange X)
            ax_comp.scatter(results_with_scale[gwp_idx].median(), results_with_scale[msp_idx].median(),
                           c=style.get_color('gold'), s=150, marker='X',
                           edgecolors='black', linewidths=1.5,
                           label='Median (Varied)', zorder=10)

            ax_comp.legend()
            ax_comp.grid(True, alpha=0.3)
            
            utils.save_figure(fig_comp, 'kde_comparison_overlap', dirs['figure'], formats=('png',))
            plt.close(fig_comp)
        except Exception as e:
            print(f"    [!] Failed to generate comparison KDE: {e}")

    # -------------------------------------------------------------------------
    # STEP 3: Tornado Diagram (Single Point Sensitivity)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 3: Tornado Diagram (MSP Sensitivity)")
    print("=" * 80)

    if sensitivity_df is not None:
        sensitivity_df = sensitivity_df.dropna()
        if not sensitivity_df.empty:
            # Create categories map
            categories = {p: get_parameter_category(p) for p in sensitivity_df['Parameter']}
            
            # Plot MSP
            fig3, ax3 = plots.plot_tornado(
                sensitivity_df=sensitivity_df,
                baseline=baseline_msp,
                metric_name='MSP [$/kg]',
                categories=categories
            )
            utils.save_figure(fig3, 'tornado_msp', dirs['figure'], formats=('png', 'pdf'))
            plt.close(fig3)

            # Check for GWP sensitivity data
            gwp_file_xlsx = os.path.join(dirs['data'], 'tornado_sensitivity_gwp.xlsx')
            gwp_file_csv = os.path.join(dirs['data'], 'tornado_sensitivity_gwp.csv')
            
            sensitivity_gwp = None
            if os.path.isfile(gwp_file_xlsx):
                sensitivity_gwp = pd.read_excel(gwp_file_xlsx)
            elif os.path.isfile(gwp_file_csv):
                sensitivity_gwp = pd.read_csv(gwp_file_csv)
            
            if sensitivity_gwp is not None:
                sensitivity_gwp = sensitivity_gwp.dropna()
                if not sensitivity_gwp.empty:
                    categories_gwp = {p: get_parameter_category(p) for p in sensitivity_gwp['Parameter']}
                    
                    fig_gwp, ax_gwp = plots.plot_tornado(
                        sensitivity_df=sensitivity_gwp,
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
    # STEP 5: Scale Effects
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 5: Scale Effects")
    print("=" * 80)
    
    scale_idx = ('Design', 'Production scale [kg/hr]')
    
    if results_with_scale is not None:
        # Check if scale column exists (handle tuple index)
        scale_col = None
        if scale_idx in results_with_scale.columns:
            scale_col = results_with_scale[scale_idx]
        else:
             # Try Partial match or string match
             for c in results_with_scale.columns:
                 if 'Production scale' in str(c):
                     scale_col = results_with_scale[c]
                     break
        
        if scale_col is not None:
            # MSP vs Scale
            if msp_idx in results_with_scale.columns:
                fig_sc1, ax_sc1 = plot_scale_effects(
                    scale_col, results_with_scale[msp_idx], baseline_msp,
                    ylabel='MSP [$/kg]', title='MSP vs Production Scale', 
                    color=style.get_color('msp') # Update color
                )
                if fig_sc1:
                    utils.save_figure(fig_sc1, 'scale_effect_msp', dirs['figure'], formats=('png',))
            
            # GWP vs Scale
            if gwp_idx in results_with_scale.columns:
                fig_sc2, ax_sc2 = plot_scale_effects(
                    scale_col, results_with_scale[gwp_idx], baseline_gwp,
                    ylabel='GWP [kg CO2-eq/kg]', title='GWP vs Production Scale', 
                    color=style.get_color('gwp') # Update color
                )
                if fig_sc2:
                    utils.save_figure(fig_sc2, 'scale_effect_gwp', dirs['figure'], formats=('png',))
        else:
            print("  [!] Production scale column not found in variable scale results")
    else:
        print("  [!] 'with_scale' results not available for scale effects analysis")

    # -------------------------------------------------------------------------
    # STEP 5b: Parameter Sensitivity Bands (Extended)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 5b: Extended Parameter Sensitivity Bands")
    print("=" * 80)
    
    # Parameters to analyze (Metric vs Parameter)
    # Using 'no_scale' results to isolate effects from scale variation
    data_for_sen = results_no_scale
    if data_for_sen is not None:
        # Define search terms for target parameters
        target_params = {
            'Productivity': ['productivity', 'rate'],
            'Titer': ['titer', 'concentration'],
            'Yield': ['yield', 'conversion'],
            'Electricity Price': ['electricity', 'power'],
            'Glucose Price': ['glucose', 'feedstock', 'sugar']
        }
        
        # Find column indices
        param_cols = {}
        for friendly_name, keywords in target_params.items():
            for col in data_for_sen.columns:
                col_str = str(col).lower()
                if any(k in col_str for k in keywords) and 'scale' not in col_str:
                    param_cols[friendly_name] = col
                    break
        
        # Plot for both MSP and GWP
        metrics = [
            (msp_idx, 'MSP', '$/kg', style.get_color('msp')),
            (gwp_idx, 'GWP', 'kg CO2-eq/kg', style.get_color('gwp'))
        ]

        for p_name, p_col in param_cols.items():
            if p_col is None:
                continue
                
            print(f"  Generating sensitivity bands for: {p_name}")
            
            # Extract readable unit if possible (usually in tuple)
            x_label = f"{p_name}"
            if isinstance(p_col, tuple) and len(p_col) > 1:
                x_label += f" [{p_col[1].split('[')[-1] if '[' in p_col[1] else ''}"
                # Clean up dirty label extraction
                if not x_label.endswith(']'): x_label = f"{p_name} [-]"
            
            for m_idx, m_name, m_unit, m_color in metrics:
                if m_idx in data_for_sen.columns:
                    try:
                        fig_sen, ax_sen = plot_scale_effects(
                            data_for_sen[p_col], 
                            data_for_sen[m_idx], 
                            baseline_metrics.get(m_idx, np.nan),
                            xlabel=x_label,
                            ylabel=f"{m_name} [{m_unit}]", 
                            title=f"{m_name} vs {p_name} (Uncertainty)", 
                            color=m_color
                        )
                        fname = f"sensitivity_{m_name}_vs_{p_name.replace(' ', '_')}"
                        utils.save_figure(fig_sen, fname, dirs['figure'], formats=('png',))
                        plt.close(fig_sen)
                    except Exception as e:
                        print(f"    [!] Failed to plot {m_name} vs {p_name}: {e}")

    # -------------------------------------------------------------------------
    # STEP 6: Contour Scatter Plots (Titer/Yield/Prod vs MSP/GWP)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("STEP 6: Contour Scatter Plots")
    print("=" * 80)

    # Use no_scale data to isolate biological effects
    data_for_scatter = results_no_scale
    if data_for_scatter is not None:
        # Find parameter indices for Titer, Yield, Productivity
        titer_idx = next((c for c in data_for_scatter.columns if 'titer' in str(c).lower()), None)
        yield_idx = next((c for c in data_for_scatter.columns if 'yield' in str(c).lower()), None)
        prod_idx = next((c for c in data_for_scatter.columns if 'productivity' in str(c).lower()), None)
        
        plot_configs = [
             # MSP Plots
             {'x': titer_idx, 'y': yield_idx, 'z': msp_idx, 'title': 'MSP_vs_Titer_Yield'},
             {'x': titer_idx, 'y': prod_idx, 'z': msp_idx, 'title': 'MSP_vs_Titer_Prod'},
             {'x': yield_idx, 'y': prod_idx, 'z': msp_idx, 'title': 'MSP_vs_Yield_Prod'},
             # GWP Plots (New)
             {'x': titer_idx, 'y': yield_idx, 'z': gwp_idx, 'title': 'GWP_vs_Titer_Yield'},
             {'x': titer_idx, 'y': prod_idx, 'z': gwp_idx, 'title': 'GWP_vs_Titer_Prod'},
             {'x': yield_idx, 'y': prod_idx, 'z': gwp_idx, 'title': 'GWP_vs_Yield_Prod'},
        ]
        
        if titer_idx and yield_idx and prod_idx and msp_idx in data_for_scatter.columns:
            for cfg in plot_configs:
                try:
                    xlabel = str(cfg['x']).split("'")[3] if "'" in str(cfg['x']) else str(cfg['x'])
                    ylabel = str(cfg['y']).split("'")[3] if "'" in str(cfg['y']) else str(cfg['y'])
                    zlabel_raw = str(cfg['z'])
                    z_name = 'MSP [$/kg]' if 'MSP' in zlabel_raw else 'GWP [kg CO2-eq/kg]'
                    
                    # Use PreFerS spectrum (reversed if cost, standard if metric)
                    # Actually MSP is usually lower=better (so reversed might be good if Red=bad, Blue=Good?)
                    # PreFerS cmap is Navy->...->Gold.
                    fig_s, ax_s = plot_colored_scatter(
                        data_for_scatter[cfg['x']], data_for_scatter[cfg['y']], data_for_scatter[cfg['z']],
                        xlabel=xlabel, ylabel=ylabel, zlabel=z_name,
                        title=cfg['title'].replace('_', ' '),
                        cmap='PreFerS', figsize=(9, 7)
                    )
                    if fig_s:
                        utils.save_figure(fig_s, f"scatter_{cfg['title']}", dirs['figure'], formats=('png',))
                        plt.close(fig_s)
                except Exception as e:
                    print(f"    [!] Failed scatter {cfg['title']}: {e}")
        else:
            print("  [!] Could not identify Titer, Yield, or Productivity columns")
    
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

    print("\n" + "=" * 80)
    print("FIGURE GENERATION COMPLETE")
    print("=" * 80)

    return dirs

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
