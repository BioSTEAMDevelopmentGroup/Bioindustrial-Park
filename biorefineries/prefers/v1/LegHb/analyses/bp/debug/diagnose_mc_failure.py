# -*- coding: utf-8 -*-
import biosteam as bst
import numpy as np
import pandas as pd
from biorefineries.prefers.v1.LegHb._models import create_model
import traceback
import warnings
warnings.filterwarnings('ignore')

def diagnose():
    print("="*80)
    print("DIAGNOSTIC RUN: Stopping on First Failure")
    print("="*80)
    
    # 1. Create Model
    print("Creating model...")
    model = create_model(baseline_production_kg_hr=150, config='config1', verbose=False)
    
    # 2. Sample (User requested 5000)
    samples_count = 20000
    seed = 1234
    print(f"Generating {samples_count} samples (Seed: {seed})...")
    np.random.seed(seed)
    parameter_samples = model.sample(samples_count, rule='L')
    
    # 3. Iterate and Run
    print("Starting evaluation...")
    
    valid_count = 0
    
    for i, sample in enumerate(parameter_samples):
        if i % 10000 == 0:
            print(f"Processing sample {i}/{samples_count}...")
            
        try:
            # Set parameters manually
            for param, value in zip(model.parameters, sample):
                param.setter(value)
            
            # Simulate System
            model.system.simulate()
            
            # Check Metrics (Finite Check)
            metrics = [m() for m in model.metrics]
            if not all(np.isfinite(metrics)):
                raise ValueError(f"Non-finite metrics detected: {metrics}")
            
            valid_count += 1
            
        except Exception as e:
            print(f"\n[!] FAILURE DETECTED at Sample {i}")
            print(f"Error Type: {type(e).__name__}")
            print(f"Error Message: {e}")
            
            print("\nParameter Values:")
            param_dict = {}
            for param, value in zip(model.parameters, sample):
                print(f"  {param.name}: {value} {param.units}")
                param_dict[param.name] = value
            
            # Save to CSV
            try:
                df = pd.DataFrame([param_dict])
                df['Error'] = str(e)
                df.to_csv('failed_params.csv', index=False)
                print(f"\nSaved failed parameters to 'failed_params.csv'")
            except:
                print("Could not save CSV.")

            print("\nTraceback:")
            traceback.print_exc()
            
            print("\nStopping diagnosis due to failure.")
            return

    print(f"\nDiagnosis Complete: All {samples_count} samples were VALID locally.")
    print("If you are seeing failures, it implies an environment difference or specific parameter sensitivity not hit here.")

if __name__ == '__main__':
    diagnose()
