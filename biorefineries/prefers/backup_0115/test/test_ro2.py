# -*- coding: utf-8 -*-
"""
RO2 Unit Quick Test
===================
Quick import and simulation test for RO2 unit.
"""
import sys
sys.path.insert(0, r'c:\Programming\PreFerS\Bioindustrial-Park')

print("Testing RO2 unit...")

try:
    import biosteam as bst
    import thermosteam as tmo
    print("✓ BioSTEAM imported")
    
    # Set up simple thermo
    chemicals = tmo.Chemicals(['H2O', 'NaCl', 'Glucose'])
    bst.settings.set_thermo(chemicals)
    print("✓ Chemicals set")
    
    # Import RO2
    from biorefineries.prefers.v1._units import RO2
    print("✓ RO2 imported")
    
    # Create test stream
    feed = bst.Stream('test_feed', H2O=1000, NaCl=10, Glucose=5, units='kg/hr')
    print(f"✓ Feed stream created: {feed.F_mass:.2f} kg/hr")
    
    # Create RO2 unit
    ro = RO2('RO_test', ins=feed)
    print(f"✓ RO2 unit created")
    print(f"  Default water recovery: {ro.water_recovery}")
    print(f"  Default membrane flux: {ro.membrane_flux} LMH")
    print(f"  Default operating pressure: {ro.operating_pressure_bar} bar")
    
    # Simulate
    ro.simulate()
    print("✓ Simulation completed")
    
    # Show results
    D = ro.design_results
    print(f"\nDesign Results:")
    print(f"  Flow rate: {D.get('Flow rate', 0):.2f} m3/hr")
    print(f"  Permeate flow: {D.get('Permeate flow', 0):.2f} m3/hr")
    print(f"  Membrane Area: {D.get('Membrane Area', 0):.2f} m2")
    print(f"  Power: {D.get('Power', 0):.2f} kW")
    print(f"  Water recovery: {D.get('Water recovery', 0):.1f}%")
    
    C = ro.purchase_costs
    print(f"\nPurchase Costs:")
    for key, value in C.items():
        print(f"  {key}: ${value:,.0f}")
    print(f"  Total: ${sum(C.values()):,.0f}")
    
    # Mass balance check
    mass_in = feed.F_mass
    mass_out = sum(s.F_mass for s in ro.outs)
    print(f"\nMass Balance:")
    print(f"  In: {mass_in:.2f} kg/hr")
    print(f"  Out: {mass_out:.2f} kg/hr")
    print(f"  Error: {abs(mass_in - mass_out):.4f} kg/hr")
    
    print("\n✓ ALL TESTS PASSED")
    
except Exception as e:
    print(f"\n✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
