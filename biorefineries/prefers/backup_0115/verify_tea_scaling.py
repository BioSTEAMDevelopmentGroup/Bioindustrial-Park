
import sys
import os
import biosteam as bst
import warnings

# Path adjustment
sys.path.append('c:/Programming/PREFERS/Bioindustrial-Park')

from biorefineries.prefers.v1.LegHb.system._config1 import create_LegHb_system
from biorefineries.prefers.v1.LegHb._tea_config1 import PreFerSTEA as LegHbTEA
from biorefineries.prefers.v1.LegHb._chemicals import create_chemicals_LegHb
from biorefineries.prefers.v1._process_settings import load_process_settings

def run_test_leghb_tea(scales):
    print("\n" + "="*50)
    print("VERIFYING LEGHB TEA SCALING")
    print("="*50)
    
    bst.settings.set_thermo(create_chemicals_LegHb(), skip_checks=True)
    load_process_settings()
    
    sys = create_LegHb_system()
    sys.simulate()
    
    # Create TEA
    tea = LegHbTEA(
        system=sys, 
        IRR=0.18, 
        duration=(2024, 2044), 
        depreciation='IRAS6',
        income_tax=0.17, 
        operating_days=333, 
        lang_factor=None,
        construction_schedule=(0.15, 0.60, 0.25), 
        WC_over_FCI=0.15,
        labor_cost=10*6e4, 
        fringe_benefits=0.17+0.07, 
        property_tax=0.005,
        property_insurance=0.005, 
        supplies=0.02, 
        maintenance=0.03,
        administration=0.05,
        target_production_kg_hr=150
    )
    
    all_passed = True
    
    for scale in scales:
        target = 150 * scale
        print(f"\n--- Testing TEA Scale={scale}x (Target: {target:.2f} kg/hr) ---")
        try:
            achieved = tea.set_production_rate(target)
            error = abs(achieved - target)
            if error < 1.0:
                print(f"[PASS] Achieved: {achieved:.2f} kg/hr")
            else:
                print(f"[FAIL] Target missed. Got {achieved:.2f}, Diff {error:.2f}")
                all_passed = False
        except Exception as e:
            print(f"[FAIL] Crashed: {e}")
            all_passed = False
            
    return all_passed


from biorefineries.prefers.v1.HemDx.system._config1 import create_NHemDx_system
from biorefineries.prefers.v1.HemDx._tea_config1 import PreFerSTEA as HemDxTEA
from biorefineries.prefers.v1.HemDx._chemicals import create_chemicals_Hemodextrin

def run_test_hemdx_tea(scales):
    print("\n" + "="*50)
    print("VERIFYING HEMDX TEA SCALING")
    print("="*50)
    
    bst.settings.set_thermo(create_chemicals_Hemodextrin(), skip_checks=True)
    load_process_settings()
    
    sys = create_NHemDx_system()
    sys.simulate()
    
    # Create TEA
    tea = HemDxTEA(
        system=sys, 
        IRR=0.18, 
        duration=(2024, 2044), 
        depreciation='IRAS6',
        income_tax=0.17, 
        operating_days=333, 
        lang_factor=None,
        construction_schedule=(0.15, 0.60, 0.25), 
        WC_over_FCI=0.15,
        labor_cost=10*6e4, 
        fringe_benefits=0.17+0.07, 
        property_tax=0.005,
        property_insurance=0.005, 
        supplies=0.02, 
        maintenance=0.03,
        administration=0.05,
        target_production_kg_hr=25
    )
    
    all_passed = True
    
    for scale in scales:
        target = 25 * scale # HemDx base is smaller
        print(f"\n--- Testing TEA Scale={scale}x (Target: {target:.2f} kg/hr) ---")
        try:
            achieved = tea.set_production_rate(target)
            error = abs(achieved - target)
            if error < 1.0:
                print(f"[PASS] Achieved: {achieved:.2f} kg/hr")
            else:
                print(f"[FAIL] Target missed. Got {achieved:.2f}, Diff {error:.2f}")
                all_passed = False
        except Exception as e:
            print(f"[FAIL] Crashed: {e}")
            all_passed = False
            
    return all_passed

if __name__ == "__main__":
    scales = [0.2, 5]
    
    leg_res = run_test_leghb_tea(scales)
    hem_res = run_test_hemdx_tea(scales)
    
    if leg_res and hem_res:
        print("\nALL TEA SCALING TESTS PASSED")
        sys.exit(0)
    else:
        print("\nSOME TEA SCALING TESTS FAILED")
        sys.exit(1)
