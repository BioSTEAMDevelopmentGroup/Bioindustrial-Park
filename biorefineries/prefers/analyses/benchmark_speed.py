# -*- coding: utf-8 -*-
"""
Benchmark Script: Serial vs Parallel Monte Carlo Simulation

Compares execution time between single-core and multi-core Monte Carlo simulations
for LegHemoglobin production uncertainty analysis.

Author: Dr. Ouwen Peng
Institute: Illinois ARCS
"""

import time
import os
import contextlib
import numpy as np
import multiprocessing as mp
from functools import partial
from warnings import filterwarnings

# Suppress warnings
filterwarnings('ignore')

# Import model creation function
from biorefineries.prefers._models import create_model


# =============================================================================
# Worker Function (Same as in uncertainty_and_sensitivity.py)
# =============================================================================

def evaluate_single_sample_quiet(sample_index_and_data, baseline_production_kg_hr):
    """
    Evaluate a single Monte Carlo sample with output suppression.
    
    Parameters
    ----------
    sample_index_and_data : tuple
        (index, sample_array) where sample_array contains parameter values
    baseline_production_kg_hr : float
        Baseline production rate [kg/hr]
        
    Returns
    -------
    tuple
        (index, param_values, metric_values, is_valid)
    """
    index, sample = sample_index_and_data
    
    try:
        # Redirect stdout to suppress simulation output
        with contextlib.redirect_stdout(open(os.devnull, 'w')):
            model = create_model(baseline_production_kg_hr=baseline_production_kg_hr)
            
            # Apply sample to parameters
            param_values = {}
            for i, param in enumerate(model.parameters):
                param.setter(sample[i])
                param_values[param.index] = sample[i]
            
            # Simulate
            model.system.simulate()
            
            # Collect metrics
            metric_values = {metric.index: metric() for metric in model.metrics}
            
            # Check validity
            is_valid = all(np.isfinite(v) for v in metric_values.values())
            
        return (index, param_values, metric_values, is_valid)
        
    except Exception as e:
        # On error, return NaN metrics
        param_values = {param.index: sample[i] for i, param in enumerate(model.parameters)}
        metric_values = {metric.index: np.nan for metric in model.metrics}
        return (index, param_values, metric_values, False)


# =============================================================================
# Serial Benchmark
# =============================================================================

def benchmark_serial(N_samples, baseline_production_kg_hr=275):
    """
    Run Monte Carlo simulation using single-core serial processing.
    
    Parameters
    ----------
    N_samples : int
        Number of Monte Carlo samples to evaluate
    baseline_production_kg_hr : float
        Production scale baseline value
        
    Returns
    -------
    float
        Elapsed time in seconds
    """
    print(f"\n[Serial Benchmark] Starting with {N_samples} samples...")
    
    # Create model
    model = create_model(baseline_production_kg_hr=baseline_production_kg_hr)
    
    # Generate samples
    np.random.seed(1234)
    samples = model.sample(N_samples, rule='L')
    
    # Prepare sample data
    sample_data = [(i, samples[i]) for i in range(N_samples)]
    
    # Start timer
    start_time = time.perf_counter()
    
    # Serial execution with output suppression
    valid_count = 0
    with contextlib.redirect_stdout(open(os.devnull, 'w')):
        for data in sample_data:
            idx, param_vals, metric_vals, is_valid = evaluate_single_sample_quiet(
                data, baseline_production_kg_hr
            )
            if is_valid:
                valid_count += 1
    
    # Stop timer
    elapsed_time = time.perf_counter() - start_time
    
    print(f"[Serial Benchmark] Completed in {elapsed_time:.2f} seconds")
    print(f"  Valid samples: {valid_count}/{N_samples} ({valid_count/N_samples*100:.1f}%)")
    
    return elapsed_time


# =============================================================================
# Parallel Benchmark
# =============================================================================

def benchmark_parallel(N_samples, baseline_production_kg_hr=275):
    """
    Run Monte Carlo simulation using multi-core parallel processing.
    
    Parameters
    ----------
    N_samples : int
        Number of Monte Carlo samples to evaluate
    baseline_production_kg_hr : float
        Production scale baseline value
        
    Returns
    -------
    float
        Elapsed time in seconds
    """
    # Determine number of workers
    n_cores = mp.cpu_count()
    n_workers = max(1, n_cores - 2)
    
    print(f"\n[Parallel Benchmark] Starting with {N_samples} samples ({n_workers} workers)...")
    
    # Create model
    model = create_model(baseline_production_kg_hr=baseline_production_kg_hr)
    
    # Generate samples
    np.random.seed(1234)
    samples = model.sample(N_samples, rule='L')
    
    # Prepare sample data
    sample_data = [(i, samples[i]) for i in range(N_samples)]
    
    # Create worker function with fixed baseline
    worker_func = partial(evaluate_single_sample_quiet, 
                         baseline_production_kg_hr=baseline_production_kg_hr)
    
    # Start timer
    start_time = time.perf_counter()
    
    # Parallel execution with output suppression
    valid_count = 0
    with contextlib.redirect_stdout(open(os.devnull, 'w')):
        with mp.Pool(processes=n_workers) as pool:
            results = []
            for result in pool.imap_unordered(worker_func, sample_data):
                results.append(result)
                if result[3]:  # is_valid
                    valid_count += 1
            pool.close()
            pool.join()
    
    # Stop timer
    elapsed_time = time.perf_counter() - start_time
    
    print(f"[Parallel Benchmark] Completed in {elapsed_time:.2f} seconds")
    print(f"  Valid samples: {valid_count}/{N_samples} ({valid_count/N_samples*100:.1f}%)")
    
    return elapsed_time


# =============================================================================
# Main Execution
# =============================================================================

if __name__ == '__main__':
    mp.freeze_support()
    
    # Configuration
    N_SAMPLES = 1000  # Number of samples for benchmark
    BASELINE_PRODUCTION = 275  # kg/hr
    
    print("="*80)
    print("MONTE CARLO SIMULATION BENCHMARK: SERIAL vs PARALLEL")
    print("="*80)
    print(f"Configuration:")
    print(f"  Samples: {N_SAMPLES}")
    print(f"  Production Scale: {BASELINE_PRODUCTION} kg/hr")
    print(f"  CPU Cores: {mp.cpu_count()}")
    print("="*80)
    
    # Run benchmarks
    serial_time = benchmark_serial(N_SAMPLES, BASELINE_PRODUCTION)
    parallel_time = benchmark_parallel(N_SAMPLES, BASELINE_PRODUCTION)
    
    # Calculate speedup
    speedup = serial_time / parallel_time if parallel_time > 0 else 0
    
    # Print summary
    print("\n" + "="*80)
    print(f"BENCHMARK RESULTS ({N_SAMPLES} samples)")
    print("="*80)
    print(f"Serial (1 Core):      {serial_time:>8.2f} seconds")
    print(f"Parallel ({mp.cpu_count() - 2} Workers):   {parallel_time:>8.2f} seconds")
    print("-"*80)
    print(f"Speedup Factor:       {speedup:>8.2f}x")
    print("="*80)
    
    # Performance insights
    if speedup > 1:
        efficiency = (speedup / (mp.cpu_count() - 2)) * 100
        print(f"\n✓ Parallel execution is {speedup:.2f}x faster")
        print(f"  Parallel efficiency: {efficiency:.1f}%")
    elif speedup < 1:
        print(f"\n⚠️  Serial execution is faster (overhead dominates)")
        print(f"  Parallel overhead: {(1/speedup - 1)*100:.1f}%")
    else:
        print(f"\n  No significant difference detected")
    
    print("\n" + "="*80)