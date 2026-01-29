
import sys
import os
# Adjust path to include the project root
sys.path.append('c:/Programming/PREFERS/Bioindustrial-Park')

from biorefineries.prefers.v1.LegHb.system._config1 import create_LegHb_system, set_production_rate, optimize_NH3_loading
import biosteam as bst

def test_scaling(nn):
    print(f"\nTesting scaling factor nn={nn}...")
    bst.preferences.N = 50
    LegHb_sys = create_LegHb_system()
    
    # Run baseline
    print("Running baseline simulation...")
    LegHb_sys.simulate()
    
    target_production = 150 * nn
    print(f"Target production: {target_production}")
    
    # Try to reset streams to fix the issue
    for s in LegHb_sys.streams:
        if s not in LegHb_sys.ins:
            s.empty()

    try:
        optimize_NH3_loading(LegHb_sys, verbose=False)
        set_production_rate(LegHb_sys, target_production)
    except Exception as e:
        print(f"FAILED at nn={nn}. Error: {e}")

if __name__ == "__main__":
    test_scaling(0.1)
    test_scaling(10)
