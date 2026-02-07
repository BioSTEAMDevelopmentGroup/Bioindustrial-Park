# -*- coding: utf-8 -*-
"""
TEA Configuration for HemDx system using _config2_new (Internalized Specifications).

Refactored from _tea_config2.py to work with the new config where:
    - adjust_glucose_for_titer → embedded in R302 add_specification
    - optimize_NH3_loading → embedded in R302 add_specification (empirical)
    - set_production_rate → simplified (no internal optimize_NH3_loading calls)
    - No run_titer_convergence needed; system.simulate() handles everything

Config 2: Intracellular production (SF=0.1).

Key changes vs _tea_config2.py:
    - Imports from _config2_new instead of _config2
    - optimize_NH3_loading method REMOVED (handled by R302 spec)
    - set_production_rate simplified (no adjust_glucose_for_titer call)
    - check_product_specifications imports from _config2_new

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from warnings import filterwarnings
filterwarnings('ignore')

import biosteam as bst
import pandas as pd
import flexsolve as flx
import warnings
from biorefineries.prefers.v1._tea import PreFerSTEA as PreFerSTEA_Base
# Import seed_targets from the NEW config module
from biorefineries.prefers.v1.HemDx.system._config2_new import seed_targets

class PreFerSTEA(PreFerSTEA_Base):
    
    # NOTE: optimize_NH3_loading is NOT needed anymore.
    # NH3 optimization is handled internally by R302's add_specification
    # in _config2_new.py (empirical measurement per titer iteration).
    
    def set_production_rate(self, target_production_kg_hr, verbose=True):
        """
        Set the target production rate and adjust system inputs accordingly.
        Delegates to _config2_new.set_production_rate which:
          - Scales all input streams proportionally
          - Relies on R302 spec for titer/NH3 convergence (no external loops)
        """
        from biorefineries.prefers.v1.HemDx.system._config2_new import set_production_rate
        
        self._target_production_kg_hr = target_production_kg_hr
        # No adjust_glucose_for_titer call needed — R302 spec handles it
        return set_production_rate(self.system, target_production_kg_hr, verbose=verbose)
    
    def check_product_specifications(self):
        """
        Check if the product stream meets all specifications.
        """
        from biorefineries.prefers.v1.HemDx.system._config2_new import check_HemDx_specifications
        
        product_stream = self.system.flowsheet.stream.NHemDx_Product
        return check_HemDx_specifications(product_stream)
    
    @property
    def target_production_kg_hr(self):
        """Target production rate in kg/hr"""
        return self._target_production_kg_hr
    
    @target_production_kg_hr.setter
    def target_production_kg_hr(self, value):
        """Set target production rate and adjust system"""
        if value is not None:
            self.set_production_rate(value)
        else:
            self._target_production_kg_hr = None

if __name__ == '__main__':
    import argparse
    from biorefineries.prefers.v1.HemDx.system._config2_new import create_NHemDx_system
    from biorefineries.prefers.v1.HemDx._chemicals import create_chemicals_Hemodextrin
    from biorefineries.prefers.v1._process_settings import load_process_settings
    
    parser = argparse.ArgumentParser(description='HemDx TEA Analysis (Config 2 NEW)')
    parser.add_argument('--production', type=float, default=150,
                        help='Target production rate in kg/hr (default: 150)')
    args, _ = parser.parse_known_args()
    
    TARGET_PRODUCTION = args.production
    
    print("="*85)
    print(f"HEMDX TEA (Config 2 NEW) - WITH INTERNALIZED DESIGN SPECIFICATION")
    print("="*85)
    
    print("\n0. Setting HemDx thermodynamics...")
    bst.settings.set_thermo(create_chemicals_Hemodextrin(), skip_checks=True)
    load_process_settings()
    
    print(f"\n1. Creating HemDx system (Config 2 NEW)...")
    HemDx_sys = create_NHemDx_system()
    HemDx_sys.operating_hours = 8000
    
    print(f"\n2. Creating TEA with target production = {TARGET_PRODUCTION} kg/hr...")
    HemDx_tea = PreFerSTEA(
        system=HemDx_sys, 
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
        target_production_kg_hr=TARGET_PRODUCTION
    )
    
    print(f"\n3. Setting production rate...")
    achieved_production = HemDx_tea.set_production_rate(TARGET_PRODUCTION)
    
    print(f"\n4. Checking product specifications...")
    try:
        HemDx_tea.check_product_specifications()
        print("   [OK] All specifications met!")
    except Exception as e:
        print(f"   [FAIL] Specification check failed: {e}")
    
    products = HemDx_sys.flowsheet.stream.NHemDx_Product
    
    print(f"\n{'='*85}")
    print("TEA RESULTS")
    print(f"{'='*85}")
    
    annual_production_kg = products.F_mass * HemDx_sys.operating_hours
    annual_production_MT = annual_production_kg / 1000
    
    print(f"\nProduction Summary:")
    print(f"  Hourly production:  {products.F_mass:.2f} kg/hr")
    print(f"  Annual production:  {annual_production_MT:.2f} metric tons/year")
    print(f"  Operating hours:    {HemDx_sys.operating_hours} hr/year")
    
    print(f"\n{'='*85}")
    HemDx_tea.show()
    
    print(f"\n{'='*85}")
    print("MINIMUM SELLING PRICE")
    print(f"{'='*85}")
    price = HemDx_tea.solve_price(products)
    print(f"  HemDx MPSP: ${price:.4f}/kg")
    print(f"  HemDx MPSP: ${price*1000:.2f}/metric ton")
    annual_revenue = price * annual_production_kg
    print(f"  Annual revenue (at MPSP): ${annual_revenue/1e6:.2f} million")
    
    print(f"\n{'='*85}")
    print("COST BREAKDOWN TABLES")
    print(f"{'='*85}")
    
    print("\nVariable Operating Costs:")
    df2 = bst.report.voc_table(HemDx_sys, 'NHemDx_Product')
    print(df2)
    
    print("\nCAPEX Breakdown:")
    df8 = HemDx_tea.CAPEX_table()
    print(df8)
    
    print("\nFixed Operating Costs:")
    df9 = HemDx_tea.FOC_table()
    print(df9)
    
    print(f"\n{'='*85}")
    print("TEA ANALYSIS COMPLETE (Config 2 NEW)")
    print(f"{'='*85}")
    print(f"Target Production:   {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"Achieved Production: {products.F_mass:.2f} kg/hr")
    print(f"HemDx MPSP:          ${price:.4f}/kg")
    print(f"{'='*85}\n")
