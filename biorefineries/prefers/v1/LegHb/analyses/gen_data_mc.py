# -*- coding: utf-8 -*-
"""
Monte Carlo Data Generation for LegHb
====================================

Runs robust Monte Carlo simulations with/without production scale variation.
"""

from warnings import filterwarnings
filterwarnings('ignore')

import argparse
import json
import multiprocessing as mp
import os
from functools import partial

import numpy as np
import pandas as pd

from biorefineries.prefers.v1.LegHb._models import create_model
from biorefineries.prefers.v1.LegHb.system import get_available_configs
from biorefineries.prefers.v1.utils import utils, sankey_utils

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    tqdm = None


# =============================================================================
# CLI
# =============================================================================

def parse_arguments():
    parser = argparse.ArgumentParser(description='LegHb Monte Carlo data generation')
    parser.add_argument('--config', type=str, default='config1',
                        choices=get_available_configs(),
                        help='Process configuration (default: config1)')
    parser.add_argument('--production', type=float, default=275,
                        help='Baseline production rate in kg/hr (default: 275)')
    parser.add_argument('--samples', type=int, default=120000,
                        help='Target number of valid samples per scenario')
    parser.add_argument('--batch-size', type=int, default=30000,
                        help='Samples per batch (default: 30000)')
    parser.add_argument('--cores', type=int, default=None,
                        help='Number of worker processes (default: max-2)')
    parser.add_argument('--no-multiprocessing', action='store_true',
                        help='Disable multiprocessing and run sequentially')
    parser.add_argument('--seed', type=int, default=1234,
                        help='Random seed for reproducibility (default: 1234)')
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
    return list(index) if isinstance(index, tuple) else index


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
            'config': config,
        }
    return utils.get_analysis_dirs(script_path, config=config, timestamp=timestamp)


# =============================================================================
# Worker Function for Parallel Execution
# =============================================================================

def evaluate_single_sample(sample_index_and_data, baseline_production_kg_hr, exclude_production_scale=False,
                          config='config1', param_indices=None, metric_indices=None):
    """
    Worker function to evaluate a single Monte Carlo sample.
    """
    index, sample = sample_index_and_data

    model = None
    try:
        model = create_model(baseline_production_kg_hr=baseline_production_kg_hr, config=config, verbose=False)

        param_values = {}
        for i, param in enumerate(model.parameters):
            is_production_scale = 'production scale' in param.name.lower()

            if exclude_production_scale and is_production_scale:
                param.setter(baseline_production_kg_hr)
                param_values[param.index] = baseline_production_kg_hr
            else:
                param.setter(sample[i])
                param_values[param.index] = sample[i]

        model.system.simulate()
        metric_values = {metric.index: metric() for metric in model.metrics}

        is_valid = all(np.isfinite(v) for v in metric_values.values())
        return (index, param_values, metric_values, is_valid)

    except Exception:
        param_values = {}
        if model is not None:
            for i, param in enumerate(model.parameters):
                is_production_scale = 'production scale' in param.name.lower()
                if exclude_production_scale and is_production_scale:
                    param_values[param.index] = baseline_production_kg_hr
                else:
                    param_values[param.index] = sample[i]
            metric_values = {metric.index: np.nan for metric in model.metrics}
        else:
            if param_indices is None:
                param_indices = []
            if metric_indices is None:
                metric_indices = []
            for i, param_idx in enumerate(param_indices):
                if 'production scale' in str(param_idx).lower() and exclude_production_scale:
                    param_values[param_idx] = baseline_production_kg_hr
                else:
                    param_values[param_idx] = sample[i]
            metric_values = {metric_idx: np.nan for metric_idx in metric_indices}
        return (index, param_values, metric_values, False)


# =============================================================================
# Robust Parallel Monte Carlo
# =============================================================================

def run_monte_carlo(model, n_target, baseline_production_kg_hr, exclude_production_scale=False,
                   batch_size=1000, n_workers=None, scenario_name="", config='config1',
                   use_multiprocessing=True):
    """
    Run robust parallel Monte Carlo simulation until n_target valid samples.
    """
    print(f"\n{'='*80}")
    print(f"RUNNING MONTE CARLO: {scenario_name}")
    print(f"{'='*80}")

    if exclude_production_scale:
        print(f"Production scale FIXED at {baseline_production_kg_hr} kg/hr")
    else:
        print(f"Production scale VARYING (baseline: {baseline_production_kg_hr} kg/hr)")

    if n_target < 1000:
        use_multiprocessing = False

    n_cores = mp.cpu_count()
    if n_workers is None:
        n_workers = max(1, n_cores - 2)

    if not use_multiprocessing or n_workers <= 1:
        print(f"\nUsing sequential evaluation (workers={n_workers})")
        return run_monte_carlo_sequential(
            model=model,
            n_target=n_target,
            baseline_production_kg_hr=baseline_production_kg_hr,
            exclude_production_scale=exclude_production_scale,
            batch_size=batch_size,
            scenario_name=scenario_name,
        )

    print(f"\nUsing {n_workers} parallel workers (of {n_cores} cores)")

    valid_results = []
    total_attempts = 0

    param_indices = [param.index for param in model.parameters]
    metric_indices = [metric.index for metric in model.metrics]

    worker_func = partial(
        evaluate_single_sample,
        baseline_production_kg_hr=baseline_production_kg_hr,
        exclude_production_scale=exclude_production_scale,
        config=config,
        param_indices=param_indices,
        metric_indices=metric_indices,
    )

    while len(valid_results) < n_target:
        remaining = n_target - len(valid_results)
        batch_to_generate = min(batch_size, int(remaining * 1.2) + 10)

        print(f"\n--- Batch {len(valid_results)//batch_size + 1} ---")
        print(f"Generating {batch_to_generate} samples...")

        batch_samples = model.sample(batch_to_generate, rule='L')
        sample_data = [(total_attempts + i, batch_samples[i]) for i in range(batch_to_generate)]

        try:
            with mp.Pool(processes=n_workers) as pool:
                iterator = pool.imap_unordered(worker_func, sample_data)
                if tqdm is not None:
                    iterator = tqdm(iterator, total=len(sample_data), desc="Evaluating", ncols=80)
                batch_results = [result for result in iterator]
                pool.close()
                pool.join()
        except Exception as exc:
            print(f"\n[!] Multiprocessing failed: {exc}")
            print("[!] Falling back to sequential evaluation for this scenario...")
            return run_monte_carlo_sequential(
                model=model,
                n_target=n_target,
                baseline_production_kg_hr=baseline_production_kg_hr,
                exclude_production_scale=exclude_production_scale,
                batch_size=batch_size,
                scenario_name=scenario_name,
            )

        batch_valid_count = 0
        batch_invalid_count = 0

        for _, param_vals, metric_vals, is_valid in batch_results:
            total_attempts += 1
            if is_valid:
                valid_results.append((len(valid_results), param_vals, metric_vals))
                batch_valid_count += 1
            else:
                batch_invalid_count += 1

        success_rate = (len(valid_results) / total_attempts * 100) if total_attempts > 0 else 0
        print(f"  Batch results: {batch_valid_count} valid, {batch_invalid_count} invalid")
        print(f"  Progress: {len(valid_results)}/{n_target} valid samples ({success_rate:.1f}% success rate)")

        if total_attempts > 100 and success_rate < 50:
            print(f"\n⚠️  WARNING: Low success rate ({success_rate:.1f}%)")

    print(f"\n✓ Monte Carlo simulation complete!")
    print(f"  Valid samples: {len(valid_results)}")
    print(f"  Total attempts: {total_attempts}")
    print(f"  Success rate: {len(valid_results)/total_attempts*100:.1f}%")

    param_indices = [param.index for param in model.parameters]
    metric_indices = [metric.index for metric in model.metrics]
    all_indices = param_indices + metric_indices

    data = []
    for _, param_vals, metric_vals in valid_results:
        row = {}
        row.update(param_vals)
        row.update(metric_vals)
        data.append(row)

    results_table = pd.DataFrame(data, columns=all_indices)

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
        'success_rate': len(valid_results) / total_attempts * 100,
        'exclude_production_scale': exclude_production_scale,
        'batch_size': batch_size,
        'n_workers': n_workers,
    }

    return results_table, simulation_stats


def run_monte_carlo_sequential(model, n_target, baseline_production_kg_hr, exclude_production_scale=False,
                              batch_size=1000, scenario_name=""):
    """
    Run Monte Carlo sequentially using a single model instance.
    """
    print(f"\n{'='*80}")
    print(f"RUNNING MONTE CARLO (SEQUENTIAL): {scenario_name}")
    print(f"{'='*80}")

    if exclude_production_scale:
        print(f"Production scale FIXED at {baseline_production_kg_hr} kg/hr")
    else:
        print(f"Production scale VARYING (baseline: {baseline_production_kg_hr} kg/hr)")

    param_indices = [param.index for param in model.parameters]
    metric_indices = [metric.index for metric in model.metrics]

    valid_chunks = []
    total_attempts = 0

    while sum(len(chunk) for chunk in valid_chunks) < n_target:
        remaining = n_target - sum(len(chunk) for chunk in valid_chunks)
        batch_to_generate = min(batch_size, int(remaining * 1.2) + 10)

        print(f"\n--- Batch {sum(len(chunk) for chunk in valid_chunks)//batch_size + 1} ---")
        print(f"Generating {batch_to_generate} samples...")

        batch_samples = model.sample(batch_to_generate, rule='L')

        if exclude_production_scale:
            for i, param in enumerate(model.parameters):
                if 'production scale' in param.name.lower():
                    batch_samples[:, i] = baseline_production_kg_hr

        model.load_samples(batch_samples, sort=False)
        model.evaluate(notify=max(1, batch_to_generate // 10))

        results_df = model.table
        valid_mask = results_df[metric_indices].apply(np.isfinite).all(axis=1)
        valid_rows = results_df.loc[valid_mask, param_indices + metric_indices]

        total_attempts += len(results_df)
        valid_chunks.append(valid_rows)

        n_valid = sum(len(chunk) for chunk in valid_chunks)
        success_rate = (n_valid / total_attempts * 100) if total_attempts else 0
        print(f"  Batch results: {len(valid_rows)} valid, {len(results_df) - len(valid_rows)} invalid")
        print(f"  Progress: {n_valid}/{n_target} valid samples ({success_rate:.1f}% success rate)")

    results_table = pd.concat(valid_chunks, ignore_index=True).iloc[:n_target]

    simulation_stats = {
        'n_valid': len(results_table),
        'n_attempts': total_attempts,
        'success_rate': len(results_table) / total_attempts * 100,
        'exclude_production_scale': exclude_production_scale,
        'batch_size': batch_size,
        'n_workers': 1,
    }

    return results_table, simulation_stats


# =============================================================================
# Data Generation
# =============================================================================

def generate_monte_carlo(config='config1', baseline_production_kg_hr=275, n_target=120000,
                          batch_size=30000, n_workers=None, seed=1234,
                          timestamp=None, results_dir=None, use_multiprocessing=True):
    dirs = resolve_output_dirs(__file__, config, timestamp=timestamp, results_dir=results_dir)

    print("=" * 80)
    print("LEGHEMOGLOBIN - MONTE CARLO DATA GENERATION")
    print("=" * 80)
    print("\nConfiguration:")
    print(f"  Process config: {config}")
    print(f"  Baseline production: {baseline_production_kg_hr} kg/hr")
    print(f"  Target samples per scenario: {n_target}")
    print(f"  Batch size: {batch_size}")
    print(f"  Random seed: {seed}")
    if n_target < 1000:
        use_multiprocessing = False
        print("  Multiprocessing disabled (n_target < 1000)")
    print(f"  Output data dir: {dirs['data']}")

    np.random.seed(seed)

    model = create_model(
        baseline_production_kg_hr=baseline_production_kg_hr,
        config=config,
        verbose=True,
    )

    baseline_metrics = model.metrics_at_baseline()

    print("\nComputing installed equipment cost breakdown...")
    step_costs = sankey_utils.compute_step_costs(model.system)
    step_totals = sankey_utils.compute_step_totals(step_costs)
    cost_breakdown_df = pd.DataFrame([
        {
            'Step': step,
            'Installed cost [$]': cost,
            'Category': 'OSBL' if step.lower() == 'facilities' else 'ISBL',
        }
        for step, cost in step_costs.items()
    ])
    cost_breakdown_summary_df = pd.DataFrame([
        {'Category': key, 'Installed cost [$]': val}
        for key, val in step_totals.items()
    ])

    cost_breakdown_csv = os.path.join(dirs['data'], 'cost_breakdown_steps.csv')
    cost_breakdown_xlsx = os.path.join(dirs['data'], 'cost_breakdown_steps.xlsx')
    cost_summary_csv = os.path.join(dirs['data'], 'cost_breakdown_summary.csv')

    cost_breakdown_df.to_csv(cost_breakdown_csv, index=False)
    cost_breakdown_df.to_excel(cost_breakdown_xlsx, index=False)
    cost_breakdown_summary_df.to_csv(cost_summary_csv, index=False)
    print(f"  [+] Saved cost breakdown: {cost_breakdown_csv}")

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
            'Value': value,
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

    results_no_scale, stats_no_scale = run_monte_carlo(
        model=model,
        n_target=n_target,
        baseline_production_kg_hr=baseline_production_kg_hr,
        exclude_production_scale=True,
        batch_size=batch_size,
        n_workers=n_workers,
        scenario_name="Fermentation Parameters Only (Fixed Scale)",
        config=config,
        use_multiprocessing=use_multiprocessing,
    )

    results_with_scale, stats_with_scale = run_monte_carlo(
        model=model,
        n_target=n_target,
        baseline_production_kg_hr=baseline_production_kg_hr,
        exclude_production_scale=False,
        batch_size=batch_size,
        n_workers=n_workers,
        scenario_name="All Parameters Including Production Scale",
        config=config,
        use_multiprocessing=use_multiprocessing,
    )

    def save_results(results_df, tag):
        pickle_path = os.path.join(dirs['data'], f'monte_carlo_{tag}.pkl')
        results_df.to_pickle(pickle_path)

        flat_df = results_df.copy()
        flat_df.columns = [flatten_column(col) for col in flat_df.columns]

        csv_path = os.path.join(dirs['data'], f'monte_carlo_{tag}.csv')
        xlsx_path = os.path.join(dirs['data'], f'monte_carlo_{tag}.xlsx')
        flat_df.to_csv(csv_path, index=True)
        flat_df.to_excel(xlsx_path, index=True)

        print(f"  [+] Saved MC results ({tag}): {csv_path}")
        print(f"  [+] Saved MC results ({tag}): {xlsx_path}")

        return {
            'pickle': os.path.basename(pickle_path),
            'csv': os.path.basename(csv_path),
            'xlsx': os.path.basename(xlsx_path),
        }

    file_no_scale = save_results(results_no_scale, 'no_scale')
    file_with_scale = save_results(results_with_scale, 'with_scale')

    column_map = {
        'columns': [
            {
                'flat': flatten_column(col),
                'original': serialize_index(col),
            }
            for col in results_no_scale.columns
        ]
    }

    column_map_file = os.path.join(dirs['data'], 'column_map.json')
    with open(column_map_file, 'w', encoding='utf-8') as f:
        json.dump(column_map, f, indent=2)

    msp_idx = ('PreFerS', 'MSP [$/kg]')
    gwp_idx = ('PreFerS', 'GWP [kg CO2-eq/kg]')

    metadata = {
        'config': config,
        'timestamp': dirs.get('timestamp'),
        'baseline_production_kg_hr': baseline_production_kg_hr,
        'n_target': n_target,
        'batch_size': batch_size,
        'seed': seed,
        'parameter_indices': [serialize_index(p.index) for p in model.parameters],
        'metric_indices': [serialize_index(m.index) for m in model.metrics],
        'baseline_msp': baseline_metrics.get(msp_idx, np.nan),
        'baseline_gwp': baseline_metrics.get(gwp_idx, np.nan),
        'baseline_metric_key_format': 'group|metric',
        'scenario_files': {
            'no_scale': file_no_scale,
            'with_scale': file_with_scale,
        },
        'scenario_stats': {
            'no_scale': stats_no_scale,
            'with_scale': stats_with_scale,
        },
        'cost_breakdown_steps_csv': os.path.basename(cost_breakdown_csv),
        'cost_breakdown_steps_xlsx': os.path.basename(cost_breakdown_xlsx),
        'cost_breakdown_summary_csv': os.path.basename(cost_summary_csv),
    }

    metadata_file = os.path.join(dirs['data'], 'analysis_metadata.json')
    with open(metadata_file, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2)

    manifest = {
        'baseline_metrics_csv': os.path.basename(baseline_csv),
        'baseline_metrics_xlsx': os.path.basename(baseline_xlsx),
        'baseline_metrics_json': os.path.basename(baseline_json),
        'column_map_json': os.path.basename(column_map_file),
        'analysis_metadata_json': os.path.basename(metadata_file),
        'monte_carlo_no_scale_pickle': file_no_scale['pickle'],
        'monte_carlo_no_scale_csv': file_no_scale['csv'],
        'monte_carlo_no_scale_xlsx': file_no_scale['xlsx'],
        'monte_carlo_with_scale_pickle': file_with_scale['pickle'],
        'monte_carlo_with_scale_csv': file_with_scale['csv'],
        'monte_carlo_with_scale_xlsx': file_with_scale['xlsx'],
        'cost_breakdown_steps_csv': os.path.basename(cost_breakdown_csv),
        'cost_breakdown_steps_xlsx': os.path.basename(cost_breakdown_xlsx),
        'cost_breakdown_summary_csv': os.path.basename(cost_summary_csv),
    }

    manifest_file = os.path.join(dirs['data'], 'analysis_manifest.json')
    with open(manifest_file, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2)

    print("\n" + "=" * 80)
    print("MONTE CARLO DATA GENERATION COMPLETE")
    print("=" * 80)
    print(f"\nResults saved to:")
    print(f"  Data: {dirs['data']}")

    return {
        'dirs': dirs,
        'baseline_metrics': baseline_metrics,
        'results_no_scale': results_no_scale,
        'results_with_scale': results_with_scale,
        'metadata': metadata,
    }


def main():
    mp.freeze_support()
    args = parse_arguments()

    generate_monte_carlo(
        config=args.config,
        baseline_production_kg_hr=args.production,
        n_target=args.samples,
        batch_size=args.batch_size,
        n_workers=args.cores,
        seed=args.seed,
        timestamp=args.timestamp,
        results_dir=args.results_dir,
        use_multiprocessing=not args.no_multiprocessing,
    )


if __name__ == '__main__':
    main()
