# -*- coding: utf-8 -*-
"""
Figure Generation for LegHb Monte Carlo Results
===============================================

Loads Monte Carlo results and generates publication-ready figures using
prefers.v1.utils.plots.
"""

from warnings import filterwarnings
filterwarnings('ignore')

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    import plotly.graph_objects as go
except Exception:  # pragma: no cover
    go = None

from biorefineries.prefers.v1.LegHb.system import get_available_configs
from biorefineries.prefers.v1.utils import plots, style, utils, plot_utils, diagram_utils


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

def deserialize_index(value):
    if isinstance(value, list) and len(value) == 2:
        return tuple(value)
    return value


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


def load_results(data_dir, filename):
    file_path = os.path.join(data_dir, filename)
    if not os.path.isfile(file_path):
        return None

    if filename.endswith('.pkl'):
        return pd.read_pickle(file_path)
    if filename.endswith('.xlsx'):
        return pd.read_excel(file_path, index_col=0)
    if filename.endswith('.csv'):
        return pd.read_csv(file_path, index_col=0)
    return None


def first_available(*items):
    for item in items:
        if item is not None:
            return item
    return None


def restore_columns(df, column_map):
    if df is None or column_map is None:
        return df

    mapping = {
        item['flat']: deserialize_index(item['original'])
        for item in column_map.get('columns', [])
    }
    df = df.copy()
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


def load_sankey_data(data_dir, filename):
    file_path = os.path.join(data_dir, filename)
    if not os.path.isfile(file_path):
        return None
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)


def save_sankey_figure(sankey_data, figure_dir, filename='sankey_baseline'):
    if go is None:
        print("  [!] Plotly not available; skipping Sankey figure.")
        return None

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color='black', width=0.5),
            label=sankey_data.get('nodes', []),
        ),
                cost_breakdown_file = metadata.get('cost_breakdown_steps_csv', 'cost_breakdown_steps.csv')
                cost_breakdown_path = os.path.join(dirs['data'], cost_breakdown_file)
                breakdown_summary = metadata.get('breakdown_summary_xlsx', 'Breakdown_Summary.xlsx')
                breakdown_summary_path = os.path.join(dirs['data'], breakdown_summary)
                step_costs = None

                if os.path.isfile(cost_breakdown_path):
                    cost_df = pd.read_csv(cost_breakdown_path)
                    step_costs = {row['Step']: row['Installed cost [$]'] for _, row in cost_df.iterrows()}
                elif os.path.isfile(breakdown_summary_path):
                    try:
                        capital_df = pd.read_excel(breakdown_summary_path, sheet_name='Total Capital')
                        step_costs = {
                            row['Name']: row['Cost (USD)']
                            for _, row in capital_df.iterrows()
                            if row.get('Name')
                        }
                    except Exception:
                        step_costs = None
        ),
    )])

    fig.update_layout(title_text=sankey_data.get('title', 'Sankey Diagram'), font_size=10)

    html_path = os.path.join(figure_dir, f"{filename}.html")
    fig.write_html(html_path)
    print(f"  ✓ Saved: {html_path}")

    try:
        png_path = os.path.join(figure_dir, f"{filename}.png")
        fig.write_image(png_path)
                    print("  [!] Cost breakdown data not found; skipping economic summary.")
    except Exception:
        pass

    return fig


def compute_spearman(results_df, parameter_indices, metric_indices, exclude_fixed=True):
    if results_df is None:
        return None

    if exclude_fixed:
        variable_params = []
        for idx in parameter_indices:
            if idx in results_df.columns and results_df[idx].std() > 1e-10:
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

    rho_df = pd.DataFrame(rho_data).T
    return rho_df


def find_param_indices(parameter_indices):
    titer_idx = next((idx for idx in parameter_indices if 'titer' in str(idx).lower()), None)
    yield_idx = next((idx for idx in parameter_indices if 'yield' in str(idx).lower()), None)
    productivity_idx = next((idx for idx in parameter_indices if 'productivity' in str(idx).lower()), None)
    return titer_idx, yield_idx, productivity_idx


# =============================================================================
# Figure Generation
# =============================================================================

def generate_figures(config='config1', timestamp=None, results_dir=None):
    style.set_style()
    dirs = resolve_output_dirs(__file__, config, timestamp=timestamp, results_dir=results_dir)

    print("=" * 80)
    print("LEGHEMOGLOBIN - FIGURE GENERATION")
    print("=" * 80)
    print(f"\nLoading data from:")
    print(f"  Data:   {dirs['data']}")
    print(f"  Figure: {dirs['figure']}")

    metadata = load_metadata(dirs['data'])
    column_map = load_column_map(dirs['data'])
    baseline_metrics = load_baseline_metrics(dirs['data'])
    sankey_carbon = load_sankey_data(dirs['data'], metadata.get('sankey_carbon_json', 'sankey_carbon.json'))
    sankey_energy = load_sankey_data(dirs['data'], metadata.get('sankey_energy_json', 'sankey_energy.json'))

    parameter_indices = [deserialize_index(i) for i in metadata.get('parameter_indices', [])]
    metric_indices = [deserialize_index(i) for i in metadata.get('metric_indices', [])]

    scenario_files = metadata.get('scenario_files', {})

    no_scale_file = scenario_files.get('no_scale', {})
    with_scale_file = scenario_files.get('with_scale', {})

    results_no_scale = first_available(
        load_results(dirs['data'], no_scale_file.get('pickle', '')),
        load_results(dirs['data'], no_scale_file.get('xlsx', '')),
        load_results(dirs['data'], no_scale_file.get('csv', '')),
    )
    results_with_scale = first_available(
        load_results(dirs['data'], with_scale_file.get('pickle', '')),
        load_results(dirs['data'], with_scale_file.get('xlsx', '')),
        load_results(dirs['data'], with_scale_file.get('csv', '')),
    )

    results_no_scale = restore_columns(results_no_scale, column_map)
    results_with_scale = restore_columns(results_with_scale, column_map)

    msp_idx = ('PreFerS', 'MSP [$/kg]')
    gwp_idx = ('PreFerS', 'GWP [kg CO2-eq/kg]')
    baseline_msp = metadata.get('baseline_msp', baseline_metrics.get(msp_idx, np.nan))
    baseline_gwp = metadata.get('baseline_gwp', baseline_metrics.get(gwp_idx, np.nan))

    print("\n" + "=" * 80)
    print("STEP 0: Sankey Diagrams (Step-wise Carbon & Energy)")
    print("=" * 80)

    if sankey_carbon:
        save_sankey_figure(sankey_carbon, dirs['figure'], filename='sankey_carbon')
    else:
        print("  [!] Sankey carbon data not found; skipping.")

    if sankey_energy:
        save_sankey_figure(sankey_energy, dirs['figure'], filename='sankey_energy')
    else:
        print("  [!] Sankey energy data not found; skipping.")

    print("\n" + "=" * 80)
    print("STEP 1: Spearman Rank Correlation (No-Scale Scenario)")
    print("=" * 80)

    if results_no_scale is not None and parameter_indices and metric_indices:
        rho_df = compute_spearman(results_no_scale, parameter_indices, metric_indices, exclude_fixed=True)
    else:
        rho_df = None

    if rho_df is not None and not rho_df.empty:
        rho_display = rho_df.copy()
        rho_display.index = [f"{idx[0]}: {idx[1]}" if isinstance(idx, tuple) else str(idx)
                             for idx in rho_display.index]
        rho_display.columns = [f"{col[1]}" if isinstance(col, tuple) else str(col)
                               for col in rho_display.columns]

        fig, ax = plots.plot_spearman_heatmap(rho_display.T, figsize=(12, 6))
        ax.set_title('Spearman Rank Correlation: Parameters vs Metrics', fontweight='bold')
        utils.save_figure(fig, 'spearman_heatmap', dirs['figure'], formats=('png', 'pdf'))

        corr_file = os.path.join(dirs['data'], 'spearman_correlations.xlsx')
        rho_display.T.to_excel(corr_file)
        print(f"  [+] Saved correlation table: {corr_file}")
    else:
        print("  [!] No valid data for Spearman correlation.")

    print("\n" + "=" * 80)
    print("STEP 2: Distribution Plots")
    print("=" * 80)

    if results_no_scale is not None and msp_idx in results_no_scale.columns:
        fig, _ = plots.plot_monte_carlo_dist(
            results=results_no_scale,
            metric_col=msp_idx,
            units='$/kg',
            baseline=baseline_msp,
            color=style.get_color('gold'),
        )
        utils.save_figure(fig, 'distribution_msp_no_scale', dirs['figure'], formats=('png',))

    if results_no_scale is not None and gwp_idx in results_no_scale.columns:
        fig, _ = plots.plot_monte_carlo_dist(
            results=results_no_scale,
            metric_col=gwp_idx,
            units='kg CO2-eq/kg',
            baseline=baseline_gwp,
            color=style.get_color('orange'),
        )
        utils.save_figure(fig, 'distribution_gwp_no_scale', dirs['figure'], formats=('png',))

    if results_with_scale is not None and msp_idx in results_with_scale.columns:
        fig, _ = plots.plot_monte_carlo_dist(
            results=results_with_scale,
            metric_col=msp_idx,
            units='$/kg',
            baseline=baseline_msp,
            color=style.get_color('gold'),
        )
        utils.save_figure(fig, 'distribution_msp_with_scale', dirs['figure'], formats=('png',))

    if results_with_scale is not None and gwp_idx in results_with_scale.columns:
        fig, _ = plots.plot_monte_carlo_dist(
            results=results_with_scale,
            metric_col=gwp_idx,
            units='kg CO2-eq/kg',
            baseline=baseline_gwp,
            color=style.get_color('orange'),
        )
        utils.save_figure(fig, 'distribution_gwp_with_scale', dirs['figure'], formats=('png',))

    print("\n" + "=" * 80)
    print("STEP 3: 2D KDE (MSP vs GWP)")
    print("=" * 80)

    if results_no_scale is not None and msp_idx in results_no_scale.columns and gwp_idx in results_no_scale.columns:
        fig, _ = plots.plot_2d_kde(
            x_data=results_no_scale[gwp_idx],
            y_data=results_no_scale[msp_idx],
            xlabel='GWP [kg CO2-eq/kg]',
            ylabel='MSP [$/kg]',
            title='Joint Distribution: GWP vs MSP (Fixed Scale)',
            cmap='PreFerS_positive',
        )
        utils.save_figure(fig, 'kde_2d_gwp_msp_no_scale', dirs['figure'], formats=('png',))

    if results_with_scale is not None and msp_idx in results_with_scale.columns and gwp_idx in results_with_scale.columns:
        fig, _ = plots.plot_2d_kde(
            x_data=results_with_scale[gwp_idx],
            y_data=results_with_scale[msp_idx],
            xlabel='GWP [kg CO2-eq/kg]',
            ylabel='MSP [$/kg]',
            title='Joint Distribution: GWP vs MSP (Variable Scale)',
            cmap='PreFerS_positive',
        )
        utils.save_figure(fig, 'kde_2d_gwp_msp_with_scale', dirs['figure'], formats=('png',))

    print("\n" + "=" * 80)
    print("STEP 3B: Joint KDE with Marginal Box Plots")
    print("=" * 80)

    if results_no_scale is not None and msp_idx in results_no_scale.columns and gwp_idx in results_no_scale.columns:
        fig, _ = plot_utils.plot_joint_kde_with_marginals(
            x=results_no_scale[gwp_idx],
            y=results_no_scale[msp_idx],
            xlabel='GWP [kg CO2-eq/kg]',
            ylabel='MSP [$/kg]',
            title='Joint KDE: GWP vs MSP (Fixed Scale)',
            cmap='PreFerS_positive',
        )
        utils.save_figure(fig, 'kde_joint_boxplot_no_scale', dirs['figure'], formats=('png',))

    print("\n" + "=" * 80)
    print("STEP 4: Colored Scatter (Key Fermentation Parameters)")
    print("=" * 80)

    if results_no_scale is not None and parameter_indices:
        titer_idx, yield_idx, productivity_idx = find_param_indices(parameter_indices)
        if all([titer_idx, yield_idx, productivity_idx]):
            scatter_configs = [
                (titer_idx, yield_idx, msp_idx, 'Titer [g/L]', 'Yield [%]', 'MSP [$/kg]', 'scatter_msp_titer_yield', 'PreFerS'),
                (titer_idx, productivity_idx, msp_idx, 'Titer [g/L]', 'Productivity [g/L/hr]', 'MSP [$/kg]', 'scatter_msp_titer_productivity', 'PreFerS'),
                (yield_idx, productivity_idx, msp_idx, 'Yield [%]', 'Productivity [g/L/hr]', 'MSP [$/kg]', 'scatter_msp_yield_productivity', 'PreFerS'),
                (titer_idx, yield_idx, gwp_idx, 'Titer [g/L]', 'Yield [%]', 'GWP [kg CO2-eq/kg]', 'scatter_gwp_titer_yield', 'PreFerS_positive'),
                (titer_idx, productivity_idx, gwp_idx, 'Titer [g/L]', 'Productivity [g/L/hr]', 'GWP [kg CO2-eq/kg]', 'scatter_gwp_titer_productivity', 'PreFerS_positive'),
                (yield_idx, productivity_idx, gwp_idx, 'Yield [%]', 'Productivity [g/L/hr]', 'GWP [kg CO2-eq/kg]', 'scatter_gwp_yield_productivity', 'PreFerS_positive'),
            ]

            for x_idx, y_idx, z_idx, xlabel, ylabel, zlabel, fname, cmap in scatter_configs:
                if z_idx not in results_no_scale.columns:
                    continue
                fig, _ = plots.plot_colored_scatter(
                    x_data=results_no_scale[x_idx],
                    y_data=results_no_scale[y_idx],
                    z_data=results_no_scale[z_idx],
                    xlabel=xlabel,
                    ylabel=ylabel,
                    zlabel=zlabel,
                    title=f"{zlabel} vs {xlabel} and {ylabel}",
                    cmap=cmap,
                )
                utils.save_figure(fig, fname, dirs['figure'], formats=('png',))
        else:
            print("  [!] Key fermentation parameters not found for scatter plots.")

    print("\n" + "=" * 80)
    print("STEP 5: Scale Effects")
    print("=" * 80)

    production_scale_idx = ('Design', 'Production scale [kg/hr]')
    if results_with_scale is not None and production_scale_idx in results_with_scale.columns:
        if msp_idx in results_with_scale.columns:
            fig, _ = plots.plot_scale_effects(
                scale_values=results_with_scale[production_scale_idx],
                metric_values=results_with_scale[msp_idx],
                baseline_val=baseline_msp,
                xlabel='Production Scale [kg/hr]',
                ylabel='MSP [$/kg]',
                title='Production Scale Effects on MSP',
            )
            utils.save_figure(fig, 'scale_effects_msp', dirs['figure'], formats=('png',))

        if gwp_idx in results_with_scale.columns:
            fig, _ = plots.plot_scale_effects(
                scale_values=results_with_scale[production_scale_idx],
                metric_values=results_with_scale[gwp_idx],
                baseline_val=baseline_gwp,
                xlabel='Production Scale [kg/hr]',
                ylabel='GWP [kg CO2-eq/kg]',
                title='Production Scale Effects on GWP',
            )
            utils.save_figure(fig, 'scale_effects_gwp', dirs['figure'], formats=('png',))
    else:
        print("  [!] Production scale data not found, skipping scale effects plots.")

    print("\n" + "=" * 80)
    print("STEP 6: Economics Summary (MSP + Cost Breakdown)")
    print("=" * 80)

    cost_breakdown_file = metadata.get('cost_breakdown_steps_csv', 'cost_breakdown_steps.csv')
    cost_breakdown_path = os.path.join(dirs['data'], cost_breakdown_file)
    if results_no_scale is not None and os.path.isfile(cost_breakdown_path):
        cost_df = pd.read_csv(cost_breakdown_path)
        step_costs = {row['Step']: row['Installed cost [$]'] for _, row in cost_df.iterrows()}
        msp_by_config = {config: results_no_scale[msp_idx].dropna().values}
        cost_by_config = {config: step_costs}

        fig, _ = plot_utils.plot_economic_panel(
            msp_by_config=msp_by_config,
            cost_by_config=cost_by_config,
            config_labels=[config],
        )
        utils.save_figure(fig, 'economic_summary', dirs['figure'], formats=('png',))
    else:
        print("  [!] Cost breakdown data not found; skipping economic summary.")

    print("\n" + "=" * 80)
    print("STEP 7: Block Flow Diagram")
    print("=" * 80)

    try:
        labels = diagram_utils.get_default_block_labels()
        fig, _ = diagram_utils.plot_block_flow_diagram(labels, title='LegHb Block Flow Diagram')
        utils.save_figure(fig, 'block_flow_diagram', dirs['figure'], formats=('png',))
    except Exception as exc:
        print(f"  [!] Failed to generate block flow diagram: {exc}")

    print("\n" + "=" * 80)
    print("FIGURE GENERATION COMPLETE")
    print("=" * 80)
    print(f"\nResults saved to:")
    print(f"  Figures: {dirs['figure']}")
    print(f"  Data:    {dirs['data']}")

    return {
        'dirs': dirs,
        'baseline_msp': baseline_msp,
        'baseline_gwp': baseline_gwp,
        'rho_df': rho_df,
    }


def main():
    args = parse_arguments()
    generate_figures(
        config=args.config,
        timestamp=args.timestamp,
        results_dir=args.results_dir,
    )


if __name__ == '__main__':
    main()
