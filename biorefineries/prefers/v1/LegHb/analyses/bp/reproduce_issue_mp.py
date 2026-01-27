# -*- coding: utf-8 -*-
import multiprocessing as mp
import biosteam as bst
from biorefineries.prefers.v1.LegHb._models import create_model
import numpy as np

def worker_func(idx_sample):
    idx, sample = idx_sample
    try:
        model = create_model(baseline_production_kg_hr=150, config='config1', verbose=False)
        for i, param in enumerate(model.parameters):
            param.setter(sample[i])
        model.system.simulate()
        
        # Check metrics
        metrics = [m() for m in model.metrics]
        is_finite = all(np.isfinite(metrics))
        return (idx, is_finite, None)
    except Exception as e:
        return (idx, False, str(e))

def run_mp_test():
    print(f"Testing Multiprocessing with {mp.cpu_count()} cores...")
    
    model = create_model(baseline_production_kg_hr=150, config='config1', verbose=False)
    samples_count = 1000
    np.random.seed(1234)
    samples = model.sample(samples_count, rule='L')
    
    data = [(i, samples[i]) for i in range(samples_count)]
    
    n_workers = max(1, mp.cpu_count() - 2)
    print(f"Using {n_workers} workers.")
    
    try:
        with mp.Pool(processes=n_workers) as pool:
            results = pool.map(worker_func, data)
            
        valid = [r for r in results if r[1]]
        failures = [r for r in results if not r[1]]
        
        print(f"Valid: {len(valid)}/{samples_count}")
        print(f"Failures: {len(failures)}/{samples_count}")
        
        if failures:
            print("Sample failure reasons:")
            reasons = [r[2] for r in failures if r[2] is not None]
            unique_reasons = set(reasons)
            for r in list(unique_reasons)[:5]:
                print(f" - {r}")
            if not reasons and failures:
                 print(" - Non-finite metrics (NaN/Inf)")

    except Exception as e:
        print(f"Multiprocessing execution failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    mp.freeze_support()
    run_mp_test()
