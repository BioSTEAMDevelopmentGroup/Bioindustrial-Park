# -*- coding: utf-8 -*-
"""
Cold Data Generation for LegHemoglobin UA/SA
===========================================

This script generates baseline, Monte Carlo, and single-point sensitivity
results and saves them to the analyses/data folder for downstream plotting.

Created on 2026-01-21
"""

from warnings import filterwarnings
filterwarnings('ignore')

import os
import json
import argparse
import numpy as np
import pandas as pd
import biosteam as bst

from biorefineries.prefers.v1.LegHb._models import create_model
from biorefineries.prefers.v1.LegHb.system import get_available_configs
from biorefineries.prefers.v1.utils import utils


# =============================================================================
# CLI
# =============================================================================

def parse_arguments():
    """Parse command line arguments for configuration."""
    parser = argparse.ArgumentParser(description='LegHemoglobin UA/SA Data Generation')
    parser.add_argument('--config', type=str, default='config1',
                        choices=get_available_configs(),
                        help='Process configuration (default: config1)')
    parser.add_argument('--production', type=float, default=275,
                        help='Baseline production rate in kg/hr (default: 275)')
    parser.add_argument('--samples', type=int, default=500,
                        help='Number of Monte Carlo samples (default: 500)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--include-production-scale', action='store_true',
                        help='Include production scale variation in Monte Carlo')
    parser.add_argument('--timestamp', type=str, default=None,
                        help='Custom timestamp (YYYYMMDD_HHMM). Defaults to current time.')
    parser.add_argument('--results-dir', type=str, default=None,
                        help='Existing results directory to write data into')
    args, _ = parser.parse_known_args()
    return args


# =============================================================================
# Utilities
# =============================================================================

def serialize_index(index):
    """Serialize a tuple index for JSON storage."""
    return list(index) if isinstance(index, tuple) else index


def deserialize_index(value):
    """Deserialize a JSON-stored index back into tuple if applicable."""
    if isinstance(value, list) and len(value) == 2:
        return tuple(value)
    return value


def metric_key_to_str(key):
    if isinstance(key, tuple):
        return f"{key[0]}|{key[1]}"
    return str(key)


def flatten_column(col):
    if isinstance(col, tuple):
        return f"{col[0]}_{col[1]}"
    return str(col)


def resolve_output_dirs(script_path, config, timestamp=None, results_dir=None):
    if results_dir:
        base_dir = os.path.abspath(results_dir)
        data_dir = os.path.join(base_dir, 'data')
        figure_dir = os.path.join(base_dir, 'figure')
        os.makedirs(data_dir, exist_ok=True)
        os.makedirs(figure_dir, exist_ok=True)
        return {
            'base': base_dir,
            'data': data_dir,
            'figure': figure_dir,
            'timestamp': os.path.basename(base_dir).split('_')[-1],
            'config': config
        }
    return utils.get_analysis_dirs(script_path, config=config, timestamp=timestamp)


# =============================================================================
# Monte Carlo Simulation
# =============================================================================

def run_monte_carlo(model, n_samples, seed=42, exclude_production_scale=True):
    """Run Monte Carlo simulation with the given model."""
    print(f"\n{'='*60}")
    print("Running Monte Carlo Simulation")
    print(f"{'='*60}")
    print(f"  Target samples: {n_samples}")
    print(f"  Exclude production scale: {exclude_production_scale}")

    np.random.seed(seed)

    samples = model.sample(N=n_samples, rule='L')  # Latin Hypercube sampling

    if exclude_production_scale:
        for i, param in enumerate(model.parameters):
            if 'production scale' in param.name.lower():
                samples[:, i] = param.baseline
                print(f"  Fixed '{param.name}' to baseline: {param.baseline}")

    model.load_samples(samples, sort=False)

    print(f"\n  Evaluating {n_samples} samples...")
    model.evaluate(notify=10)

    results_df = model.table.copy()

    metric_indices = [m.index for m in model.metrics]
    valid_mask = results_df[metric_indices].notna().all(axis=1)
    n_valid = valid_mask.sum()

    print(f"\n  Valid samples: {n_valid} / {n_samples} ({100*n_valid/n_samples:.1f}%)")

    return results_df[valid_mask]


# =============================================================================
# Single-Point Sensitivity for Tornado
# =============================================================================

def compute_single_point_sensitivity(model, metric_name='MSP [$/kg]', perturbation=0.1):
    """Compute single-point sensitivity for tornado diagram."""
    print(f"\n{'='*60}")
    print(f"Computing Single-Point Sensitivity for {metric_name}")
    print(f"{'='*60}")

    baseline_metrics = model.metrics_at_baseline()

    baseline_value = None
    for key, val in baseline_metrics.items():
        if metric_name in str(key):
            baseline_value = val
            break

    if baseline_value is None:
        raise ValueError(f"Metric '{metric_name}' not found")

    print(f"  Baseline {metric_name}: {baseline_value:.4f}")

    def _fmt(value):
        if value is None or pd.isna(value):
            return "nan"
        return f"{value:.4f}"

    sensitivity_data = []

    for param in model.parameters:
        param_baseline = param.baseline
        lb = param_baseline * (1 - perturbation)
        ub = param_baseline * (1 + perturbation)

        param.setter(lb)
        try:
            model.system.simulate()
            metrics_low = {m.index: m() for m in model.metrics}
            value_low = None
            for key, val in metrics_low.items():
                if metric_name in str(key):
                    value_low = val
                    break
        except Exception:
            value_low = np.nan

        param.setter(ub)
        try:
            model.system.simulate()
            metrics_high = {m.index: m() for m in model.metrics}
            value_high = None
            for key, val in metrics_high.items():
                if metric_name in str(key):
                    value_high = val
                    break
        except Exception:
            value_high = np.nan

        param.setter(param_baseline)

        if isinstance(param.index, tuple):
            param_name = f"{param.index[0]}: {param.index[1]}"
        else:
            param_name = str(param.index)

        sensitivity_data.append({
            'Parameter': param_name,
            'Low': value_low if value_low is not None else np.nan,
            'High': value_high if value_high is not None else np.nan
        })

        print(f"  {param_name}: Low={_fmt(value_low)}, High={_fmt(value_high)}")

    model.metrics_at_baseline()

    return pd.DataFrame(sensitivity_data), baseline_value


# =============================================================================
# Data Generation
# =============================================================================

def generate_data(config='config1', baseline_production_kg_hr=275, n_samples=500,
                  seed=42, exclude_production_scale=True, timestamp=None,
                  results_dir=None):
    """Generate baseline, Monte Carlo, and sensitivity datasets."""
    dirs = resolve_output_dirs(__file__, config, timestamp=timestamp, results_dir=results_dir)

    print("="*70)
    print("LEGHEMOGLOBIN - UA/SA COLD DATA GENERATION")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  Process config: {config}")
    print(f"  Baseline production: {baseline_production_kg_hr} kg/hr")
    print(f"  Monte Carlo samples: {n_samples}")
    print(f"  Random seed: {seed}")
    print(f"  Output data dir: {dirs['data']}")

    model = create_model(
        baseline_production_kg_hr=baseline_production_kg_hr,
        config=config,
        verbose=True
    )

    print("\nModel Parameters:")
    for param in model.parameters:
        print(f"  - {param.name} [{param.units}]: baseline={param.baseline:.4f}")

    print("\nModel Metrics:")
    for metric in model.metrics:
        print(f"  - {metric.name} [{metric.units}]")

    print("\n" + "="*70)
    print("STEP 1: Evaluating Baseline")
    print("="*70)

    baseline_metrics = model.metrics_at_baseline()
    baseline_rows = []
    metrics_map = {m.index: m for m in model.metrics}

    for idx, value in baseline_metrics.items():
        metric = metrics_map.get(idx, None)
        if isinstance(idx, tuple):
            group, name = idx
        else:
            group, name = "", str(idx)
        baseline_rows.append({
            'Group': group,
            'Metric': name,
            'Units': metric.units if metric else "",
            'Value': value
        })

    baseline_table = pd.DataFrame(baseline_rows)

    baseline_csv = os.path.join(dirs['data'], 'baseline_metrics.csv')
    baseline_xlsx = os.path.join(dirs['data'], 'baseline_metrics.xlsx')
    baseline_json = os.path.join(dirs['data'], 'baseline_metrics.json')

    baseline_table.to_csv(baseline_csv, index=False)
    baseline_table.to_excel(baseline_xlsx, index=False)

    baseline_serialized = {metric_key_to_str(k): v for k, v in baseline_metrics.items()}
    with open(baseline_json, 'w', encoding='utf-8') as f:
        json.dump(baseline_serialized, f, indent=2)

    print(f"  ✓ Saved baseline metrics: {baseline_csv}")
    print(f"  ✓ Saved baseline metrics: {baseline_xlsx}")
    print(f"  ✓ Saved baseline metrics: {baseline_json}")

    print("\n" + "="*70)
    print("STEP 2: Monte Carlo Simulation")
    print("="*70)

    mc_results = run_monte_carlo(
        model=model,
        n_samples=n_samples,
        seed=seed,
        exclude_production_scale=exclude_production_scale
    )

    mc_pickle = os.path.join(dirs['data'], 'monte_carlo_results.pkl')
    mc_results.to_pickle(mc_pickle)

    mc_results_flat = mc_results.copy()
    mc_results_flat.columns = [flatten_column(col) for col in mc_results.columns]

    mc_csv = os.path.join(dirs['data'], 'monte_carlo_results.csv')
    mc_xlsx = os.path.join(dirs['data'], 'monte_carlo_results.xlsx')

    mc_results_flat.to_csv(mc_csv, index=True)
    mc_results_flat.to_excel(mc_xlsx, index=True)

    print(f"  ✓ Saved MC results: {mc_pickle}")
    print(f"  ✓ Saved MC results: {mc_csv}")
    print(f"  ✓ Saved MC results: {mc_xlsx}")

    column_map = {
        'columns': [
            {
                'flat': flatten_column(col),
                'original': serialize_index(col)
            }
            for col in mc_results.columns
        ]
    }

    column_map_file = os.path.join(dirs['data'], 'column_map.json')
    with open(column_map_file, 'w', encoding='utf-8') as f:
        json.dump(column_map, f, indent=2)

    print(f"  ✓ Saved column map: {column_map_file}")

    print("\n" + "="*70)
    print("STEP 3: Tornado Sensitivity (Single-Point)")
    print("="*70)

    sensitivity_df, msp_baseline = compute_single_point_sensitivity(
        model=model,
        metric_name='MSP [$/kg]',
        perturbation=0.10
    )

    sensitivity_df = sensitivity_df.dropna()

    sens_csv = os.path.join(dirs['data'], 'tornado_sensitivity.csv')
    sens_xlsx = os.path.join(dirs['data'], 'tornado_sensitivity.xlsx')
    sensitivity_df.to_csv(sens_csv, index=False)
    sensitivity_df.to_excel(sens_xlsx, index=False)

    print(f"  ✓ Saved sensitivity table: {sens_csv}")
    print(f"  ✓ Saved sensitivity table: {sens_xlsx}")

    msp_idx = ('PreFerS', 'MSP [$/kg]')
    gwp_idx = ('PreFerS', 'GWP [kg CO2-eq/kg]')
    baseline_msp = baseline_metrics.get(msp_idx, np.nan)
    baseline_gwp = baseline_metrics.get(gwp_idx, np.nan)

    metadata = {
        'config': config,
        'timestamp': dirs.get('timestamp'),
        'baseline_production_kg_hr': baseline_production_kg_hr,
        'n_samples': n_samples,
        'seed': seed,
        'exclude_production_scale': exclude_production_scale,
        'parameter_indices': [serialize_index(p.index) for p in model.parameters],
        'metric_indices': [serialize_index(m.index) for m in model.metrics],
        'baseline_msp': baseline_msp,
        'baseline_gwp': baseline_gwp,
        'baseline_metric_key_format': 'group|metric'
    }

    metadata_file = os.path.join(dirs['data'], 'analysis_metadata.json')
    with open(metadata_file, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2)

    print(f"  ✓ Saved metadata: {metadata_file}")

    manifest = {
        'baseline_metrics_csv': os.path.basename(baseline_csv),
        'baseline_metrics_xlsx': os.path.basename(baseline_xlsx),
        'baseline_metrics_json': os.path.basename(baseline_json),
        'monte_carlo_results_pickle': os.path.basename(mc_pickle),
        'monte_carlo_results_csv': os.path.basename(mc_csv),
        'monte_carlo_results_xlsx': os.path.basename(mc_xlsx),
        'column_map_json': os.path.basename(column_map_file),
        'tornado_sensitivity_csv': os.path.basename(sens_csv),
        'tornado_sensitivity_xlsx': os.path.basename(sens_xlsx),
        'analysis_metadata_json': os.path.basename(metadata_file)
    }

    manifest_file = os.path.join(dirs['data'], 'analysis_manifest.json')
    with open(manifest_file, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2)

    print(f"  ✓ Saved manifest: {manifest_file}")

    print("\n" + "="*70)
    print("DATA GENERATION COMPLETE")
    print("="*70)
    print(f"\nResults saved to:")
    print(f"  Data:    {dirs['data']}")

    return {
        'dirs': dirs,
        'baseline_metrics': baseline_metrics,
        'mc_results': mc_results,
        'sensitivity_df': sensitivity_df,
        'metadata': metadata
    }


def main():
    args = parse_arguments()
    exclude_production_scale = not args.include_production_scale

    generate_data(
        config=args.config,
        baseline_production_kg_hr=args.production,
        n_samples=args.samples,
        seed=args.seed,
        exclude_production_scale=exclude_production_scale,
        timestamp=args.timestamp,
        results_dir=args.results_dir
    )


if __name__ == '__main__':
    main()
