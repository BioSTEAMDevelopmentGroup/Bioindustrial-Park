# -*- coding: utf-8 -*-
"""
Cold Data Generation for LegHemoglobin UA/SA
===========================================

This script generates baseline, Monte Carlo, and single-point sensitivity
results and saves them to the analyses/data folder for downstream plotting.

Features:
- Baseline metrics evaluation
- Monte Carlo uncertainty analysis with Latin Hypercube sampling
- Single-point tornado sensitivity analysis
- Dual-mode processing:
  * Sequential (default): Safe, cross-platform, ~11 samples/sec
  * Parallel (experimental): Optional multiprocessing with --use-multiprocessing flag
    WARNING: Multiprocessing may not work on Windows due to BioSTEAM/Fortran 
    compatibility issues. Use sequential mode for reliable production runs.

Usage:
  python gen_data.py --samples 1000                      # Sequential (default)
  python gen_data.py --samples 1000 --use-multiprocessing  # Parallel (may fail on Windows)

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
import multiprocessing as mp
from functools import partial

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
    parser.add_argument('--use-multiprocessing', action='store_true',
                        help='Use multiprocessing (experimental, may not work on Windows)')
    parser.add_argument('--mp-workers', type=int, default=None,
                        help='Number of multiprocessing workers to use (overrides auto).')
    parser.add_argument('--mp-batch-size', type=int, default=200,
                        help='Batch size for multiprocessing (default: 200).')
    parser.add_argument('--mp-report-every', type=int, default=200,
                        help='Progress report frequency in samples (default: 200).')
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
# Worker Function for Parallel Execution
# =============================================================================

def evaluate_single_sample(sample_index_and_data, baseline_production_kg_hr, config, param_indices):
    """
    Worker function to evaluate a single Monte Carlo sample in parallel.
    
    WARNING: This function is designed for multiprocessing.Pool workers. On Windows,
    BioSTEAM model creation in worker processes may cause Fortran runtime errors due to
    object serialization limitations. Use sequential mode (default) for reliable operation.
    
    Parameters
    ----------
    sample_index_and_data : tuple
        (index, sample_array) where sample_array contains parameter values
    baseline_production_kg_hr : float
        Baseline production rate
    config : str
        Process configuration
    param_indices : list
        List of parameter indices
    
    Returns
    -------
    tuple
        (index, param_values_dict, metric_values_dict, is_valid)
    """
    index, sample = sample_index_and_data
    
    try:
        # Create a fresh model for this worker
        model = create_model(baseline_production_kg_hr=baseline_production_kg_hr, config=config, verbose=False)
        
        # Dictionary to store parameter values used in simulation
        param_values = {}
        
        # Apply the sample to the model parameters
        for i, param in enumerate(model.parameters):
            value = sample[i]
            param.setter(value)
            param_values[param.index] = value
        
        # Run simulation
        model.system.simulate()
        
        # Collect metric values
        metric_values = {metric.index: metric() for metric in model.metrics}
        
        # Check validity (all finite values)
        is_valid = all(np.isfinite(v) for v in metric_values.values())
        
        return (index, param_values, metric_values, is_valid)
        
    except Exception as e:
        # On error, report NaN metrics
        param_values = {}
        for i, param_idx in enumerate(param_indices):
            param_values[param_idx] = sample[i]
        
        # Create dummy model to get metric indices
        try:
            model = create_model(baseline_production_kg_hr=baseline_production_kg_hr, config=config, verbose=False)
            metric_values = {metric.index: np.nan for metric in model.metrics}
        except:
            metric_values = {}
        
        return (index, param_values, metric_values, False)


# =============================================================================
# Monte Carlo Simulation - Sequential (Safe) and Parallel (Experimental) Options
# =============================================================================

def run_monte_carlo_sequential(model, n_samples, seed=42, exclude_production_scale=True):
    """
    Run Monte Carlo simulation sequentially (safe, cross-platform).
    
    This is the default and recommended method for all production runs.
    Performance: ~11 samples/second on typical workstation.
    
    Parameters
    ----------
    model : biosteam.Model
        BioSTEAM model with defined parameters and metrics
    n_samples : int
        Number of Monte Carlo samples to generate
    seed : int
        Random seed for reproducibility
    exclude_production_scale : bool
        If True, fix production scale to baseline value
        
    Returns
    -------
    pd.DataFrame
        Results dataframe with parameter values and metric outputs
    """
    print(f"\n{'='*60}")
    print("Running Monte Carlo Simulation (SEQUENTIAL)")
    print(f"{'='*60}")
    print(f"  Target samples: {n_samples}")
    print(f"  Exclude production scale: {exclude_production_scale}")

    np.random.seed(seed)

    samples = model.sample(N=n_samples, rule='L')

    if exclude_production_scale:
        for i, param in enumerate(model.parameters):
            if 'production scale' in param.name.lower():
                samples[:, i] = param.baseline
                print(f"  Fixed '{param.name}' to baseline: {param.baseline}")

    model.load_samples(samples, sort=False)

    print(f"\n  Evaluating {n_samples} samples...")
    import time
    start_time = time.time()
    
    model.evaluate(notify=max(1, n_samples // 20))  # Show progress every 5%

    elapsed = time.time() - start_time
    print(f"  Completed in {elapsed:.1f} seconds (~{n_samples/elapsed:.1f} samples/sec)")

    results_df = model.table.copy()

    metric_indices = [m.index for m in model.metrics]
    valid_mask = results_df[metric_indices].notna().all(axis=1)
    n_valid = valid_mask.sum()

    print(f"  Valid samples: {n_valid} / {n_samples} ({100*n_valid/n_samples:.1f}%)")

    return results_df[valid_mask]


def run_monte_carlo_parallel(model, n_samples, seed=42, exclude_production_scale=True, config='config1', baseline_production_kg_hr=275, n_workers=None):
    """
    Run Monte Carlo simulation with batched multiprocessing (experimental).
    
    WARNING: This function may not work on Windows due to BioSTEAM's Fortran dependencies
    and multiprocessing serialization issues. It is primarily intended for Linux/Mac systems.
    
    The function processes samples in batches using multiprocessing.Pool to improve throughput.
    If any batch fails, the entire operation will raise an exception and fall back to sequential.
    
    Known Issues:
    - Windows: "forrtl: error (200)" when creating BioSTEAM models in worker processes
    - Pickle serialization limitations with complex BioSTEAM objects
    - Worker process startup overhead may negate performance gains for small sample sizes
    
    Parameters
    ----------
    model : biosteam.Model
        BioSTEAM model (used only for sampling; workers create their own models)
    n_samples : int
        Number of Monte Carlo samples
    seed : int
        Random seed
    exclude_production_scale : bool
        If True, fix production scale to baseline
    config : str
        Process configuration name
    baseline_production_kg_hr : float
        Baseline production rate
        
    Returns
    -------
    pd.DataFrame
        Results dataframe with valid samples only
    """
    print(f"\n{'='*60}")
    print("Running Monte Carlo Simulation (PARALLEL - EXPERIMENTAL)")
    print(f"{'='*60}")
    print(f"  Target samples: {n_samples}")
    print(f"  Exclude production scale: {exclude_production_scale}")
    print(f"  WARNING: Multiprocessing may not work reliably on Windows with BioSTEAM")
    
    # Determine batch size
    if n_samples <= 100:
        batch_size = 100
        n_workers = n_workers or 2
    elif n_samples <= 1000:
        batch_size = 500
        n_workers = n_workers or max(2, mp.cpu_count() - 2)
    else:
        batch_size = 1000
        n_workers = n_workers or max(2, mp.cpu_count() - 2)
    
    print(f"  Using {n_workers} parallel workers")
    print(f"  Batch size: {batch_size} samples per batch")

    np.random.seed(seed)
    samples = model.sample(N=n_samples, rule='L')

    if exclude_production_scale:
        for i, param in enumerate(model.parameters):
            if 'production scale' in param.name.lower():
                samples[:, i] = param.baseline
                print(f"  Fixed '{param.name}' to baseline: {param.baseline}")

    print(f"\n  Evaluating {n_samples} samples...")
    
    param_indices = [param.index for param in model.parameters]
    metric_indices = [metric.index for metric in model.metrics]
    
    all_results = []
    valid_count = 0
    
    import time
    start_time = time.time()
    
    for batch_start in range(0, n_samples, batch_size):
        batch_end = min(batch_start + batch_size, n_samples)
        batch_samples = samples[batch_start:batch_end]
        
        sample_data = [(batch_start + i, batch_samples[i]) for i in range(len(batch_samples))]
        worker_func = partial(evaluate_single_sample, 
                             baseline_production_kg_hr=baseline_production_kg_hr,
                             config=config,
                             param_indices=param_indices)
        
        try:
            with mp.Pool(processes=n_workers) as pool:
                batch_results = pool.map(worker_func, sample_data, chunksize=5)
            
            for idx, param_vals, metric_vals, is_valid in batch_results:
                if is_valid:
                    row = {}
                    row.update(param_vals)
                    row.update(metric_vals)
                    all_results.append(row)
                    valid_count += 1
                    
        except Exception as e:
            print(f"  [!] Batch failed: {e}")
            raise
    
    elapsed = time.time() - start_time
    print(f"  Completed in {elapsed:.1f} seconds (~{n_samples/elapsed:.1f} samples/sec)")
    print(f"  Valid samples: {valid_count} / {n_samples} ({100*valid_count/n_samples:.1f}%)")
    
    if valid_count == 0:
        raise ValueError("No valid samples generated!")

    all_indices = param_indices + metric_indices
    results_df = pd.DataFrame(all_results, columns=all_indices)
    
    return results_df


def run_monte_carlo(model, n_samples, seed=42, exclude_production_scale=True, 
                   config='config1', baseline_production_kg_hr=275, use_multiprocessing=False, n_workers=None):
    """
    Run Monte Carlo simulation.
    
    Parameters
    ----------
    use_multiprocessing : bool
        If True, use parallel processing (experimental, may fail on Windows).
        If False, use sequential processing (safe, slower but reliable).
    """
    if use_multiprocessing:
        try:
            return run_monte_carlo_parallel(model, n_samples, seed, exclude_production_scale,
                                           config, baseline_production_kg_hr, n_workers=n_workers)
        except Exception as e:
            print(f"\n[!] Multiprocessing failed: {e}")
            print("[!] Falling back to sequential processing...")
            return run_monte_carlo_sequential(model, n_samples, seed, exclude_production_scale)
    else:
        return run_monte_carlo_sequential(model, n_samples, seed, exclude_production_scale)


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
                  results_dir=None, use_multiprocessing=False, mp_workers=None):
    """Generate baseline, Monte Carlo, and sensitivity datasets."""
    dirs = resolve_output_dirs(__file__, config, timestamp=timestamp, results_dir=results_dir)

    print("="*70)
    print("LEGHEMOGLOBIN - UA/SA DATA GENERATION")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  Process config: {config}")
    print(f"  Baseline production: {baseline_production_kg_hr} kg/hr")
    print(f"  Monte Carlo samples: {n_samples}")
    print(f"  Random seed: {seed}")
    print(f"  Multiprocessing: {'ENABLED (experimental)' if use_multiprocessing else 'DISABLED (sequential)'}")
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

    print(f"  [+] Saved baseline metrics: {baseline_csv}")
    print(f"  [+] Saved baseline metrics: {baseline_xlsx}")
    print(f"  [+] Saved baseline metrics: {baseline_json}")

    print("\n" + "="*70)
    print("STEP 2: Monte Carlo Simulation")
    print("="*70)

    mc_results = run_monte_carlo(
        model=model,
        n_samples=n_samples,
        seed=seed,
        exclude_production_scale=exclude_production_scale,
        config=config,
        baseline_production_kg_hr=baseline_production_kg_hr,
        use_multiprocessing=use_multiprocessing,
        n_workers=mp_workers
    )

    mc_pickle = os.path.join(dirs['data'], 'monte_carlo_results.pkl')
    mc_results.to_pickle(mc_pickle)

    mc_results_flat = mc_results.copy()
    mc_results_flat.columns = [flatten_column(col) for col in mc_results.columns]

    mc_csv = os.path.join(dirs['data'], 'monte_carlo_results.csv')
    mc_xlsx = os.path.join(dirs['data'], 'monte_carlo_results.xlsx')

    mc_results_flat.to_csv(mc_csv, index=True)
    mc_results_flat.to_excel(mc_xlsx, index=True)

    print(f"  [+] Saved MC results: {mc_pickle}")
    print(f"  [+] Saved MC results: {mc_csv}")
    print(f"  [+] Saved MC results: {mc_xlsx}")

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

    print(f"  [+] Saved column map: {column_map_file}")

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

    print(f"  [+] Saved sensitivity table: {sens_csv}")
    print(f"  [+] Saved sensitivity table: {sens_xlsx}")

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

    print(f"  [+] Saved metadata: {metadata_file}")

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

    print(f"  [+] Saved manifest: {manifest_file}")

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
    mp.freeze_support()  # Required for Windows multiprocessing
    args = parse_arguments()
    exclude_production_scale = not args.include_production_scale
    
    # Platform warning for multiprocessing
    if args.use_multiprocessing:
        import platform
        if platform.system() == 'Windows':
            print("\n" + "="*70)
            print("[!] WINDOWS MULTIPROCESSING WARNING")
            print("="*70)
            print("  Multiprocessing on Windows may cause runtime errors with BioSTEAM.")
            print("  If you experience crashes, re-run WITHOUT --use-multiprocessing flag.")
            print("  Proceeding with multiprocessing (will auto-fallback if it fails)...")
            print("="*70 + "\n")

    generate_data(
        config=args.config,
        baseline_production_kg_hr=args.production,
        n_samples=args.samples,
        seed=args.seed,
        exclude_production_scale=exclude_production_scale,
        timestamp=args.timestamp,
        results_dir=args.results_dir,
        use_multiprocessing=args.use_multiprocessing,
        mp_workers=args.mp_workers
    )


if __name__ == '__main__':
    main()

