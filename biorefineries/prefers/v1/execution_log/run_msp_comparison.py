import time
import sys
import os
import biosteam as bst
from biosteam import settings

# Inject Heartbeat
print("[HEARTBEAT] Script started", flush=True)

# Add project root to path
sys.path.append(r"c:\Programming\PreFerS\Bioindustrial-Park")

try:
    from biorefineries.prefers.v1.LegHb.system._config2 import create_LegHb_system, set_production_rate
    from biorefineries.prefers.v1.LegHb._chemicals import create_chemicals_LegHb
    from biorefineries.prefers.v1._process_settings import load_process_settings
    from biorefineries.prefers.v1.LegHb._tea import create_LegHb_tea
    
    # Initialize basic settings
    bst.preferences.classic_mode()
    bst.settings.set_thermo(create_chemicals_LegHb(), skip_checks=True)
    load_process_settings()
    
    multipliers = [0.5, 1, 2]
    results = {}
    
    print(f"\n{'='*60}")
    print("RUNNING MSP COMPARISON FOR LegHb")
    print(f"{'='*60}")
    
    for n in multipliers:
        print(f"\n--- Running Case: n = {n} ---")
        
        target_production = 275 * n
        print(f"Target Production: {target_production} kg/hr")
        
        LegHb_sys = create_LegHb_system()
        LegHb_sys.operating_hours = 8000
        
        achieved = set_production_rate(LegHb_sys, target_production)
        print(f"Achieved Production: {achieved:.2f} kg/hr")
        
        tea = create_LegHb_tea(LegHb_sys)
        tea.IRR = 0.10
        
        try:
            msp = tea.solve_price(LegHb_sys.flowsheet.stream.LegHb_3)
            print(f"MSP: ${msp:.4f}/kg")
            results[n] = msp
        except Exception as e:
            print(f"Error calculating MSP for n={n}: {e}")
            results[n] = None
        print("-" * 30)

    print(f"\n{'='*60}")
    print("COMPARISON RESULTS")
    print(f"{'='*60}")
    print(f"{'Multiplier (n)':<15} | {'Production (kg/hr)':<20} | {'MSP ($/kg)':<15}")
    print("-" * 56)
    
    for n in multipliers:
        prod = 275 * n
        msp = results.get(n)
        msp_str = f"{msp:.4f}" if msp is not None else "ERROR"
        print(f"{n:<15} | {prod:<20.2f} | {msp_str:<15}")
        
    print(f"{'='*60}")

except Exception as e:
    print(f"\nCRITICAL ERROR: {e}")
    import traceback
    traceback.print_exc()

print("DONE")
