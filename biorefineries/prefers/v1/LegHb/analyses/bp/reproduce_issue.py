# -*- coding: utf-8 -*-
import biosteam as bst
import numpy as np
import pandas as pd
from biorefineries.prefers.v1.LegHb._models import create_model

def reproduce():
    print("Creating model...")
    # Use config1 as default per user request context
    model = create_model(baseline_production_kg_hr=150, config='config1', verbose=True)
    
    samples = 1000
    seed = 1234
    
    print(f"Generating {samples} samples with seed {seed}...")
    np.random.seed(seed)
    parameter_samples = model.sample(samples, rule='L')
    model.load_samples(parameter_samples)
    
    print("Evaluating...")
    # Run evaluation
    # We catch exceptions to see what's happening
    
    print("Evaluating manually to catch errors...")
    
    valid_count = 0
    errors = []
    
    # Run a few manually to see errors
    for i in range(min(50, len(parameter_samples))):
        sample = parameter_samples[i]
        try:
            # model.parameters includes setters that update the system
            # We need to manually set them if we aren't using evaluate
            # But model.parameters is a list of Parameter objects.
            # We can use model.load_samples but that prepares for evaluate.
            # To run manually:
            for param, value in zip(model.parameters, sample):
                param.setter(value)
            
            model.system.simulate()
            valid_count += 1
        except Exception as e:
            errors.append(str(e))
            if len(errors) <= 5:
                print(f"Sample {i} failed: {e}")
            
    print(f"\nManual check of first 50: {valid_count} valid")
    if errors:
        print(f"Common errors (first 5 unique): {list(set(errors))[:5]}")

    # Now run full evaluate to get the table
    print("\nRunning model.evaluate()...")
    model.evaluate(notify=100)
    table = model.table
    
    print("\nTable Columns:", table.columns.tolist())
    
    # Check validity again
    # Usually success is metrics not being NaN
    if not table.empty:
        # Check for NaNs
        nan_rows = table[table.isna().any(axis=1)]
        print(f"Total rows: {len(table)}")
        print(f"Rows with NaNs: {len(nan_rows)}")
        
        if len(nan_rows) > 0:
             print("First row with NaNs:")
             print(nan_rows.iloc[0])
    else:
        print("Table is empty!")

if __name__ == '__main__':
    reproduce()
