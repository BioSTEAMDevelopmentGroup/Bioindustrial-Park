
import warnings
warnings.filterwarnings('ignore')
import biosteam as bst
import sys
import os

# Ensure local path is in sys.path
sys.path.append('c:/Programming/PreFerS/Bioindustrial-Park')

from biorefineries.prefers.v1.LegHb._system import create_LegHb_system
from biorefineries.prefers.v1.LegHb import _streams as s
from biorefineries.prefers.v1.LegHb import _chemicals as c

def verify():
    print("Initializing LegHb System...")
    
    # Debug: Direct creation
    print("Testing create_chemicals_LegHb()...")
    chems = c.create_chemicals_LegHb()
    print("Directly created chemicals count:", len(chems))
    if 'Globin' in chems.IDs:
        print("DEBUG: Globin IS in created chemicals.")
    else:
        print("DEBUG: Globin IS NOT in created chemicals.")
        if 'YeastExtract' in chems.IDs:
             print("DEBUG: YeastExtract found in created chemicals (unexpected?).")
    
    # Set it explicitly
    bst.settings.set_thermo(chems)
    
    # Setup placeholder outputs
    outs = [bst.Stream(f'out{i}') for i in range(6)]
    
    # Create streams from kwargs
    input_stream_kwargs = [
        s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt,
        s.DfUltraBuffer, s.IXEquilibriumBuffer, s.IXElutionBuffer,
        s.IXRegenerationSolution, s.DfNanoBuffer
    ]
    
    # Diagnostic: Check chemicals
    print("Checking chemicals...")
    if 'Globin' in bst.settings.chemicals.IDs:
        print("Globin FOUND in settings.")
    else:
        print("Globin NOT FOUND in settings.")
        print("Available:", bst.settings.chemicals.IDs)

    input_streams = [bst.Stream(**kwargs) for kwargs in input_stream_kwargs]
    
    # Create system
    results = create_LegHb_system(
        ins=input_streams,
        outs=outs
    )
    
    sys_obj = bst.main_flowsheet.system.LegHb_sys
    
    print("Running Simulation...")
    try:
        sys_obj.simulate()
        print("Simulation Converged Successfully.")
    except Exception as e:
        print(f"Simulation Failed: {e}")
        # Continue to inspect whatever state exists
        
    u = bst.main_flowsheet.unit
    
    # Verification checks
    print("\n--- Filtration Verification (S402) ---")
    if hasattr(u, 'S402'):
        s402 = u.S402
        print(f"Class: {type(s402).__name__}")
        print(f"Preset: {getattr(s402, 'preset', 'None')}")
        print(f"Solids Loading: {s402.solids_loading} kg/m2/hr")
        area = s402.design_results.get('Filter Area', 0)
        power = s402.design_results.get('Power', 0)
        cost = s402.baseline_purchase_costs.get('Rotary Vacuum Filter', 0)
        print(f"Calculated Area: {area:.2f} m2")
        print(f"Calculated Power: {power:.2f} kW")
        print(f"Calculated Cost: ${cost:,.2f}")
        
        # Check mass balance
        m_in = s402.ins[0].F_mass
        m_out = sum(s.F_mass for s in s402.outs)
        print(f"Mass Balance Error: {abs(m_in - m_out):.6f} kg/hr")
        
        if area > 0 and cost > 0:
             print("PASS: S402 is sized and costed.")
        else:
             print("FAIL: S402 sizing/costing suspicious.")
             
    else:
        print("FAIL: S402 not found.")

    print("\n--- Diafiltration Verification (U401) ---")
    if hasattr(u, 'U401'):
        u401 = u.U401
        print(f"Class: {type(u401).__name__}")
        print(f"Preset: {getattr(u401, 'preset', 'None')}")
        print(f"Flux: {u401.membrane_flux_LMH} LMH")
        area = u401.design_results.get('Membrane Area', 0)
        cost = u401.baseline_purchase_costs.get('Membrane System', 0) + u401.baseline_purchase_costs.get('Membrane replacement', 0)
        print(f"Calculated Area: {area:.2f} m2")
        print(f"Total Membrane Cost: ${cost:,.2f}")
        
        if area > 0:
             print("PASS: U401 is sized.")
        else:
             print("FAIL: U401 size is zero.")

if __name__ == '__main__':
    verify()
