# -*- coding: utf-8 -*-
"""
LegHb Production System - Config 2 (Titer = 3.5 g/L)

Same flowsheet and settings as Config 1, but with reduced titer:
    - titer_LegHb: 5.0 → 3.5 g/L
    - productivity_LegHb: 5/72 → 3.5/72 g/L/h
    - Y_p: 5*4/650 → 3.5*4/650

All other parameters, process areas, and specifications are identical to Config 1.

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biorefineries.prefers.v2.LegHb.system import _config1 as _base

# Save original function reference BEFORE any monkey-patching
_original_get_fermentation_parameters = _base.get_fermentation_parameters

# Re-export shared objects
seed_targets = _base.seed_targets

__all__ = (
    'create_LegHb_system',
    'set_production_rate',
    'check_LegHb_specifications',
    'get_fermentation_parameters',
    'seed_targets',
)


def get_fermentation_parameters():
    """
    Return fermentation parameters for Config 2 (titer = 3.5 g/L).
    All other parameters identical to Config 1.
    """
    params = _original_get_fermentation_parameters()
    params['titer_LegHb'] = 3.5
    params['productivity_LegHb'] = 3.5 / 72
    params['Y_p'] = 3.5 * 4 / 650
    params['yield_LegHb'] = params['Y_p']
    return params


def create_LegHb_system(**kwargs):
    """
    Create LegHb system with Config 2 parameters (titer = 3.5 g/L).

    Temporarily patches _config1.get_fermentation_parameters so that
    the SystemFactory-decorated function uses Config 2 titer values.
    """
    original_fn = _base.get_fermentation_parameters
    _base.get_fermentation_parameters = get_fermentation_parameters
    try:
        sys = _base.create_LegHb_system(**kwargs)
    finally:
        _base.get_fermentation_parameters = original_fn
    return sys


# Re-export unchanged functions
set_production_rate = _base.set_production_rate
check_LegHb_specifications = _base.check_LegHb_specifications


# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == '__main__':
    import biosteam as bst
    
    bst.preferences.N = 50
    nn = 1
    TARGET_PRODUCTION = 150 * nn
    
    print("="*85)
    print("LEGHEMOGLOBIN PRODUCTION SYSTEM - CONFIG2 (Titer = 3.5 g/L)")
    print("="*85)
    
    print("\n1. Creating system...")
    LegHb_sys = create_LegHb_system()
    sys = LegHb_sys
    f = sys.flowsheet
    u = f.unit
    ss = f.stream
    sys.operating_hours = 8000
    
    print("\n2. Running baseline simulation...")
    try:
        LegHb_sys.simulate()
        baseline_production = ss.LegHb_3.F_mass
        print(f"   Baseline production rate: {baseline_production:.2f} kg/hr")
    except Exception as e:
        print(f"   Baseline simulation failed: {e}")
        raise

    print(f"\n3. Applying design specification: TARGET_PRODUCTION = {TARGET_PRODUCTION} kg/hr")
    try:
        target_titer = u.R302.titer if u.R302.titer is not None else u.R302.target_titer
        if target_titer is not None:
            u.R302.titer = float(target_titer)
            u.R302.target_titer = float(target_titer)
            
        print(f"   Target Titer: {target_titer} g/L")
        
        set_production_rate(sys, TARGET_PRODUCTION, verbose=True)
        
        final_production = ss.LegHb_3.F_mass
        
        if abs(final_production - TARGET_PRODUCTION) > 1.0:
            print(f"\n   WARNING: Production rate drifted after final simulation!")
            print(f"   Target:  {TARGET_PRODUCTION:.2f} kg/hr")
            print(f"   Actual:  {final_production:.2f} kg/hr")
            
    except Exception as e:
        print(f"\n   Could not achieve target production: {e}")
        print("   Continuing with baseline production rate...")
        import traceback; traceback.print_exc()
    
    print(f"\n4. Verifying product specifications...")
    try:
        check_LegHb_specifications(ss.LegHb_3)
    except ValueError as e:
        print(f"\n   SPECIFICATION CHECK FAILED: {e}")
    
    print(f"\n5. System Summary")
    print("="*85)
    LegHb_sys.show()
    
    legh_purity = ss.LegHb_3.imass['Leghemoglobin'] / ss.LegHb_3.F_mass * 100
    
    # Calculate Heme Equivalent
    legh_mol = ss.LegHb_3.imol['Leghemoglobin']
    heme_equiv_mass = 0
    try:
        heme_mw = _base.LEGHB_THERMO.Heme_b.MW
        heme_equiv_mass = (legh_mol / 763.0) * heme_mw
    except:
        pass

    print(f"\n{'='*85}")
    print("KEY PERFORMANCE INDICATORS")
    print(f"{'='*85}")
    print(f"  Product Stream:           {ss.LegHb_3.ID}")
    print(f"  Production Rate:          {ss.LegHb_3.F_mass:.2f} kg/hr")
    print(f"  Leghemoglobin Content:    {legh_purity:.2f}%")
    print(f"  Annual Production:        {ss.LegHb_3.F_mass * sys.operating_hours / 1000:.2f} metric tons/year")
    print(f"  Yield (Heme eq):          {heme_equiv_mass:.4f} kg/hr")
    print(f"{'='*85}\n")
    
    print(f"\n6. Broth Analysis")
    print(f"   target titer: {u.R302.titer:.2f} g/L")
    print(f"   Achieved titer: {u.R302.actual_titer:.2f} g/L")
    print(f"   Achieved titer2: {(ss.Broth.imass['Leghemoglobin']+ss.Broth.imass['Leghemoglobin_In']) / ss.Broth.F_vol:.2f} g/L")
    print(ss.Broth)
    
    print(f"\n7. Performing LCA analysis...")
    try:
        r1 = bst.report.lca_inventory_table(
            systems=[sys],
            keys='GWP',
            items=[ss.LegHb_3],
        )
        print("   LCA Inventory Table generated")
        print(r1)
        r2 = bst.report.lca_displacement_allocation_table(
            systems=[sys],
            key='GWP',
            items=[ss.LegHb_3],
        )
        print("   LCA Displacement Allocation Table generated")
        print(r2)
    except Exception as e:
        print(f"   LCA analysis failed: {e}")
    
    print(f"\n{'='*85}")
    print("SIMULATION COMPLETE")
    print(f"Target Production:   {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"Achieved Production: {ss.LegHb_3.F_mass:.2f} kg/hr")
    print(f"{'='*85}\n")
