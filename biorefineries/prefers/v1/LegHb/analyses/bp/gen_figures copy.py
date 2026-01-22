# -*- coding: utf-8 -*-
"""
Figure Generation for LegHemoglobin UA/SA
=========================================

This script loads precomputed "cold data" from the analyses/data folder
and generates all plots and derived datasets without re-running simulations.

Created on 2026-01-21
"""

from warnings import filterwarnings
filterwarnings('ignore')

import os
import json
import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from biorefineries.prefers.v1.LegHb.system import get_available_configs
from biorefineries.prefers.v1.utils import style, plots, utils


# =============================================================================
# CLI
# =============================================================================

def parse_arguments():
    """Parse command line arguments for configuration."""
    parser = argparse.ArgumentParser(description='LegHemoglobin UA/SA Figure Generation')
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

def deserialize_index(value):
    if isinstance(value, list) and len(value) == 2:
        return tuple(value)
    return value


def resolve_output_dirs(script_path, config, timestamp=None, results_dir=None):
    if results_dir:
        base_dir = os.path.abspath(results_dir)
        data_dir = os.path.join(base_dir, 'data')
        figure_dir = os.path.join(base_dir, 'figure')
        return {
            'base': base_dir,
            'data': data_dir,
            'figure': figure_dir,
            'timestamp': os.path.basename(base_dir).split('_')[-1],
            'config': config
        }
    return utils.get_analysis_dirs(script_path, config=config, timestamp=timestamp, create=False)


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
    pickle_path = os.path.join(data_dir, 'monte_carlo_results.pkl')
    if os.path.isfile(pickle_path):
        return pd.read_pickle(pickle_path)

    xlsx_path = os.path.join(data_dir, 'monte_carlo_results.xlsx')
    csv_path = os.path.join(data_dir, 'monte_carlo_results.csv')
    column_map = load_column_map(data_dir)

    if os.path.isfile(xlsx_path):
        df = pd.read_excel(xlsx_path, index_col=0)
    elif os.path.isfile(csv_path):
        df = pd.read_csv(csv_path, index_col=0)
    else:
        raise FileNotFoundError('No Monte Carlo results found in data directory.')

    if column_map:
        mapping = {item['flat']: deserialize_index(item['original'])
                   for item in column_map.get('columns', [])}
        df.columns = [mapping.get(col, col) for col in df.columns]

    return df


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
    raise FileNotFoundError('No tornado sensitivity table found in data directory.')


# =============================================================================
# Spearman Rank Correlation
# =============================================================================

def compute_spearman_correlations(results_df, parameter_indices, metric_indices, exclude_fixed=True):
    """Compute Spearman rank correlations between parameters and metrics."""
    if exclude_fixed:
        variable_params = []
        for idx in parameter_indices:
            if results_df[idx].std() > 1e-10:
                variable_params.append(idx)
            else:
                print(f"  Excluding fixed parameter: {idx}")
        parameter_indices = variable_params

    rho_data = {}
    for metric_idx in metric_indices:
        rho_row = {}
        for param_idx in parameter_indices:
            try:
                rho, _ = spearmanr(results_df[param_idx], results_df[metric_idx])
                rho_row[param_idx] = rho
            except Exception:
                rho_row[param_idx] = np.nan
        rho_data[metric_idx] = rho_row

    return pd.DataFrame(rho_data).T


# =============================================================================
# Figure Generation
# =============================================================================

def generate_figures(config='config1', timestamp=None, results_dir=None):
    dirs = resolve_output_dirs(__file__, config, timestamp=timestamp, results_dir=results_dir)

    print("="*70)
    print("LEGHEMOGLOBIN - UA/SA FIGURE GENERATION")
    print("="*70)
    print(f"\nLoading data from:")
    print(f"  Data:   {dirs['data']}")
    print(f"  Figure: {dirs['figure']}")

    metadata = load_metadata(dirs['data'])
    parameter_indices = [deserialize_index(i) for i in metadata.get('parameter_indices', [])]
    metric_indices = [deserialize_index(i) for i in metadata.get('metric_indices', [])]

    mc_results = load_mc_results(dirs['data'])
    baseline_metrics = load_baseline_metrics(dirs['data'])
    sensitivity_df = load_sensitivity_table(dirs['data'])

    msp_idx = ('PreFerS', 'MSP [$/kg]')
    gwp_idx = ('PreFerS', 'GWP [kg CO2-eq/kg]')
    baseline_msp = metadata.get('baseline_msp', baseline_metrics.get(msp_idx, np.nan))
    baseline_gwp = metadata.get('baseline_gwp', baseline_metrics.get(gwp_idx, np.nan))

    print("\n" + "="*70)
    print("STEP 1: Spearman Rank Correlation Analysis")
    print("="*70)

    if parameter_indices and metric_indices:
        rho_df = compute_spearman_correlations(mc_results, parameter_indices, metric_indices, exclude_fixed=True)
    else:
        rho_df = compute_spearman_correlations(
            mc_results,
            parameter_indices=[c for c in mc_results.columns if c not in metric_indices],
            metric_indices=metric_indices or [c for c in mc_results.columns if isinstance(c, tuple)],
            exclude_fixed=True
        )

    rho_display = rho_df.copy()
    rho_display.index = [f"{idx[0]}: {idx[1]}" if isinstance(idx, tuple) else str(idx)
                         for idx in rho_display.index]
    rho_display.columns = [f"{col[1]}" if isinstance(col, tuple) else str(col)
                           for col in rho_display.columns]

    fig1, ax1 = plots.plot_spearman_heatmap(rho_display.T, figsize=(12, 6))
    ax1.set_title('Spearman Rank Correlation: Parameters vs Metrics', fontweight='bold')

    utils.save_figure(fig1, 'spearman_heatmap', dirs['figure'], formats=('png', 'pdf'))

    corr_file = os.path.join(dirs['data'], 'spearman_correlations.xlsx')
    rho_display.T.to_excel(corr_file)
    print(f"  [+] Saved correlation table to: {corr_file}")

    print("\n" + "="*70)
    print("STEP 2: Joint Marginal Plot (Cost vs. Carbon Trade-off)")
    print("="*70)

    if msp_idx in mc_results.columns and gwp_idx in mc_results.columns:
        g = plots.plot_joint_marginal(
            data=mc_results,
            x_col=gwp_idx,
            y_col=msp_idx,
            x_label='GWP [kg CO₂-eq/kg]',
            y_label='MSP [$/kg]'
        )

        g.figure.suptitle('Cost vs. Carbon Trade-off Analysis',
                          fontweight='bold', y=1.02)

        utils.save_figure(g.figure, 'joint_msp_gwp', dirs['figure'],
                          formats=('png', 'pdf'))
    else:
        print("  [!]️  MSP or GWP metric not found in results, skipping joint plot")

    print("\n" + "="*70)
    print("STEP 3: Tornado Diagram (MSP Sensitivity)")
    print("="*70)

    sensitivity_df = sensitivity_df.dropna()

    fig3, ax3 = plots.plot_tornado(
        sensitivity_df=sensitivity_df,
        baseline=baseline_msp,
        metric_name='MSP [$/kg]'
    )

    utils.save_figure(fig3, 'tornado_msp', dirs['figure'], formats=('png', 'pdf'))

    print("\n" + "="*70)
    print("STEP 4: Distribution Plots")
    print("="*70)

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

    print("\n" + "="*70)
    print("FIGURE GENERATION COMPLETE")
    print("="*70)
    print(f"\nResults saved to:")
    print(f"  Figures: {dirs['figure']}")
    print(f"  Data:    {dirs['data']}")

    return {
        'dirs': dirs,
        'spearman': rho_df,
        'sensitivity_df': sensitivity_df,
        'baseline_msp': baseline_msp,
        'baseline_gwp': baseline_gwp
    }


def main():
    args = parse_arguments()
    generate_figures(
        config=args.config,
        timestamp=args.timestamp,
        results_dir=args.results_dir
    )


if __name__ == '__main__':
    main()

