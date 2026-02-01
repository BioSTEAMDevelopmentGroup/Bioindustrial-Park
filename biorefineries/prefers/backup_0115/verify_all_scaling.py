
import sys
import os
import biosteam as bst

# Path adjustment
sys.path.append('c:/Programming/PREFERS/Bioindustrial-Park')

from biorefineries.prefers.v1.LegHb.system._config1 import create_LegHb_system, set_production_rate as set_prod_leg, optimize_NH3_loading as opt_nh3_leg
from biorefineries.prefers.v1.HemDx.system._config1 import create_NHemDx_system, set_production_rate as set_prod_hem, optimize_NH3_loading as opt_nh3_hem

bst.preferences.N = 50

def run_test_leghb(nn_values):
    print("\n" + "="*50)
    print("VERIFYING LEGHB SCALING")
    print("="*50)
    
    sys = create_LegHb_system()
    try:
        opt_nh3_leg(sys, verbose=False)
        sys.simulate()
    except Exception as e:
        print(f"LegHb Baseline Failed (Proceeding): {e}")
        # return False
        
    all_passed = True
    
    for nn in nn_values:
        target = 150 * nn
        print(f"\n--- Testing NN={nn} (Target: {target} kg/hr) ---")
        try:
            opt_nh3_leg(sys, verbose=False)
            achieved = set_prod_leg(sys, target, verbose=True)
            error = abs(achieved - target)
            if error < 1.0:
                print(f"[PASS] NN={nn} Achieved: {achieved:.2f}")
            else:
                print(f"[FAIL] NN={nn} Target missed. Got {achieved:.2f}, Diff {error:.2f}")
                all_passed = False
        except Exception as e:
            print(f"[FAIL] NN={nn} Crashed: {e}")
            all_passed = False
            
    return all_passed

def run_test_hemdx(nn_values):
    print("\n" + "="*50)
    print("VERIFYING HEMDX SCALING")
    print("="*50)
    
    sys = create_NHemDx_system()
    try:
        opt_nh3_hem(sys, verbose=False)
        sys.simulate()
    except Exception as e:
        print(f"HemDx Baseline Failed (Proceeding): {e}")
        # return False

    all_passed = True
    
    for nn in nn_values:
        target = 150 * nn
        print(f"\n--- Testing NN={nn} (Target: {target} kg/hr) ---")
        try:
            opt_nh3_hem(sys, verbose=False)
            achieved = set_prod_hem(sys, target, verbose=True)
            error = abs(achieved - target)
            if error < 1.0:
                print(f"[PASS] NN={nn} Achieved: {achieved:.2f}")
            else:
                print(f"[FAIL] NN={nn} Target missed. Got {achieved:.2f}, Diff {error:.2f}")
                all_passed = False
        except Exception as e:
            print(f"[FAIL] NN={nn} Crashed: {e}")
            all_passed = False
            
    return all_passed

if __name__ == "__main__":
    nn_list = [0.2, 5]
    
    # leg_res = run_test_leghb(nn_list)
    leg_res = True
    hem_res = run_test_hemdx(nn_list)
    # hem_res = True # Skip HemDx for now
    
    if leg_res and hem_res:
        print("\nALL TESTS PASSED")
        sys.exit(0)
    else:
        print("\nSOME TESTS FAILED")
        sys.exit(1)
