# -*- coding: utf-8 -*-
"""
VacuumPSA Unit Verification Test
================================
Created: 2026-01-18
Author: PreFerS Development Team

Verification script for VacuumPSA unit operation.
Tests mass balance, separation performance, and costing.
"""

import sys
import biosteam as bst
import thermosteam as tmo
import numpy as np

# Change: Import VacuumPSA at module level first to handle thermo initialization issues
try:
    from biorefineries.prefers.v1._units import VacuumPSA
except ImportError:
    sys.path.append(r'c:\Programming\PreFerS\Bioindustrial-Park')
    from biorefineries.prefers.v1._units import VacuumPSA


# === Configuration ===
TOLERANCE_MASS = 1e-6  # kg/hr tolerance for mass balance
TOLERANCE_RELATIVE = 0.01  # 1% relative tolerance for comparisons


def setup_chemicals():
    """Set up the chemical property package for gas separation."""
    # Define chemicals AFTER import
    chemicals = tmo.Chemicals([
        'H2',
        'CO',
        'C2H4',  # Ethylene
        'CO2',
        'N2',
        'CH4',
    ])
    bst.settings.set_thermo(chemicals)
    print("✓ Chemicals loaded: H2, CO, C2H4, CO2, N2, CH4")
    return chemicals


def create_test_feed():
    """Create a representative syngas feed stream for testing."""
    feed = bst.Stream(
        'syngas_feed',
        H2=50,      # kmol/hr
        CO=30,      # kmol/hr
        C2H4=15,    # kmol/hr
        CO2=3,      # kmol/hr
        N2=2,       # kmol/hr
        units='kmol/hr',
        phase='g',
        T=298.15,   # K (25°C)
        P=6e5,      # Pa (6 bar)
    )
    return feed


def test_instantiation():
    """Test that the unit can be instantiated without errors."""
    print("\n=== Test 1: Instantiation ===")
    
    try:
        # Pass factory chemicals explicitly via thermo to ensure compatibility
        # even if module imports messed with defaults
        current_thermo = bst.settings.get_thermo()
        
        feed = create_test_feed()
        unit = VacuumPSA('VPSA1', ins=feed, thermo=current_thermo)
        print(f"✓ Unit instantiated: {unit}")
        print(f"  Split factors: {unit.split}")
        return True, unit
    except Exception as e:
        print(f"✗ Instantiation failed: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def test_simulation(unit):
    """Test that the unit simulates without errors."""
    print("\n=== Test 2: Simulation ===")
    
    try:
        unit.simulate()
        print(f"✓ Simulation completed")
        print(f"  Design results: {dict(unit.design_results)}")
        print(f"  Purchase costs: {dict(unit.purchase_costs)}")
        return True
    except Exception as e:
        print(f"✗ Simulation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_mass_balance(unit):
    """Test that mass is conserved."""
    print("\n=== Test 3: Mass Balance ===")
    
    mass_in = sum(s.F_mass for s in unit.ins if s)
    mass_out = sum(s.F_mass for s in unit.outs if s)
    
    error = abs(mass_in - mass_out)
    
    print(f"  Mass in:  {mass_in:.6f} kg/hr")
    print(f"  Mass out: {mass_out:.6f} kg/hr")
    print(f"  Error:    {error:.6e} kg/hr")
    
    if error < TOLERANCE_MASS:
        print(f"✓ Mass balance within tolerance ({TOLERANCE_MASS})")
        return True
    else:
        print(f"✗ Mass balance error exceeds tolerance")
        return False


def test_molar_balance(unit):
    """Test that moles are conserved."""
    print("\n=== Test 4: Molar Balance ===")
    
    mol_in = sum(s.F_mol for s in unit.ins if s)
    mol_out = sum(s.F_mol for s in unit.outs if s)
    
    error = abs(mol_in - mol_out)
    
    print(f"  Moles in:  {mol_in:.6f} kmol/hr")
    print(f"  Moles out: {mol_out:.6f} kmol/hr")
    print(f"  Error:     {error:.6e} kmol/hr")
    
    if error < 1e-9:
        print(f"✓ Molar balance within tolerance")
        return True
    else:
        print(f"✗ Molar balance error exceeds tolerance")
        return False


def test_separation_performance(unit):
    """Test separation performance against expected splits."""
    print("\n=== Test 5: Separation Performance ===")
    
    feed = unit.ins[0]
    product = unit.outs[0]
    purge = unit.outs[1]
    
    expected_splits = {
        'H2': 0.95,
        'CO': 0.30,
        'C2H4': 0.20,
        'CO2': 0.10,
        'N2': 0.85,
    }
    
    all_correct = True
    print("  Component | Feed (kmol/hr) | Product | Purge | Actual Split | Expected")
    print("  " + "-" * 75)
    
    for comp, expected_split in expected_splits.items():
        feed_mol = feed.imol[comp]
        product_mol = product.imol[comp]
        purge_mol = purge.imol[comp]
        
        if feed_mol > 0:
            actual_split = product_mol / feed_mol
        else:
            actual_split = 0
        
        match = abs(actual_split - expected_split) < 0.001
        status = "✓" if match else "✗"
        
        print(f"  {comp:9} | {feed_mol:14.2f} | {product_mol:7.2f} | {purge_mol:5.2f} | {actual_split:12.3f} | {expected_split:.3f} {status}")
        
        if not match:
            all_correct = False
    
    if all_correct:
        print("✓ All separation factors match expected values")
    else:
        print("✗ Some separation factors do not match")
    
    return all_correct


def test_design_results(unit):
    """Test that expected design results are present."""
    print("\n=== Test 6: Design Results ===")
    
    D = unit.design_results
    
    expected_keys = [
        'Feed gas flow',
        'Adsorbent mass',
        'Bed volume',
        'Vacuum power',
        'H2 recovery',
    ]
    
    all_present = True
    for key in expected_keys:
        if key in D:
            print(f"  ✓ {key}: {D[key]:.4g} {unit._units.get(key, '')}")
        else:
            print(f"  ✗ {key}: MISSING")
            all_present = False
    
    if all_present:
        print(f"✓ All expected design results present")
    else:
        print(f"✗ Some design results missing")
    
    return all_present


def test_purchase_costs(unit):
    """Test that expected purchase costs are present and reasonable."""
    print("\n=== Test 7: Purchase Costs ===")
    
    C = unit.purchase_costs
    
    expected_keys = [
        'Pressure vessels',
        'Vacuum pump',
        'Adsorbent',
        'Piping and valves',
    ]
    
    all_valid = True
    total_cost = 0
    
    for key in expected_keys:
        if key in C:
            cost = C[key]
            total_cost += cost
            if cost >= 0:
                print(f"  ✓ {key}: ${cost:,.0f}")
            else:
                print(f"  ✗ {key}: ${cost:,.0f} (negative!)")
                all_valid = False
        else:
            print(f"  ✗ {key}: MISSING")
            all_valid = False
    
    print(f"\n  Total purchase cost: ${total_cost:,.0f}")
    print(f"  Installed cost: ${unit.installed_cost:,.0f}")
    
    if all_valid:
        print("✓ All purchase costs valid")
    else:
        print("✗ Some purchase costs invalid")
    
    return all_valid


def test_stream_properties(unit):
    """Test that output stream properties are correctly set."""
    print("\n=== Test 8: Stream Properties ===")
    
    feed = unit.ins[0]
    product = unit.outs[0]
    purge = unit.outs[1]
    
    all_correct = True
    
    # Check phases
    if product.phase == 'g':
        print(f"  ✓ Product phase: gas")
    else:
        print(f"  ✗ Product phase: {product.phase} (expected 'g')")
        all_correct = False
    
    if purge.phase == 'g':
        print(f"  ✓ Purge phase: gas")
    else:
        print(f"  ✗ Purge phase: {purge.phase} (expected 'g')")
        all_correct = False
    
    # Check temperatures
    if abs(product.T - feed.T) < 1e-6:
        print(f"  ✓ Product temperature: {product.T:.2f} K")
    else:
        print(f"  ✗ Product temperature mismatch: {product.T} vs feed {feed.T}")
        all_correct = False
    
    # Check pressures
    if product.P == unit.P_ads:
        print(f"  ✓ Product pressure: {product.P/1e5:.2f} bar (adsorption)")
    else:
        print(f"  ✗ Product pressure: {product.P/1e5} bar (expected {unit.P_ads/1e5})")
        all_correct = False
    
    if purge.P == unit.P_des:
        print(f"  ✓ Purge pressure: {purge.P/1e5:.3f} bar (desorption)")
    else:
        print(f"  ✗ Purge pressure: {purge.P/1e5} bar (expected {unit.P_des/1e5})")
        all_correct = False
    
    if all_correct:
        print("✓ All stream properties correct")
    else:
        print("✗ Some stream properties incorrect")
    
    return all_correct


def run_all_tests():
    """Run the complete verification suite."""
    print("=" * 70)
    print("VacuumPSA Unit Verification")
    print("=" * 70)
    
    # Setup
    setup_chemicals()
    
    # Run tests
    results = {}
    
    success, unit = test_instantiation()
    results['instantiation'] = success
    
    if success:
        results['simulation'] = test_simulation(unit)
        
        if results['simulation']:
            results['mass_balance'] = test_mass_balance(unit)
            results['molar_balance'] = test_molar_balance(unit)
            results['separation'] = test_separation_performance(unit)
            results['design_results'] = test_design_results(unit)
            results['purchase_costs'] = test_purchase_costs(unit)
            results['stream_properties'] = test_stream_properties(unit)
            
            # Show unit summary
            print("\n" + "=" * 70)
            print("UNIT SUMMARY")
            print("=" * 70)
            unit.show()
    
    # Summary
    print("\n" + "=" * 70)
    print("VERIFICATION SUMMARY")
    print("=" * 70)
    
    passed = sum(results.values())
    total = len(results)
    
    for test_name, result in results.items():
        status = "PASS" if result else "FAIL"
        print(f"  {test_name}: {status}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n✓ VERIFICATION SUCCESSFUL")
        return True
    else:
        print("\n✗ VERIFICATION FAILED")
        return False


if __name__ == '__main__':
    run_all_tests()
