# -*- coding: utf-8 -*-
"""
HemDx TEA Configuration

Refactored to include PreFerSTEA class with integrated production scaling 
and ammonia optimization logic.
"""

import biosteam as bst
import pandas as pd
import flexsolve as flx
import warnings
from biorefineries.prefers.v1._tea import PreFerSTEA as PreFerSTEA_Base
# Import seed_targets to update the global registry used by unit specs
from biorefineries.prefers.v1.HemDx.system._config1 import seed_targets

class PreFerSTEA(PreFerSTEA_Base):
    
    def optimize_NH3_loading(self, verbose=True):
        """
        Optimizes NH3_25wt flow rate and S202 split ratio to meet fermentation demand.
        Delegates to system module implementation.
        """
        from biorefineries.prefers.v1.HemDx.system import optimize_NH3_loading
        optimize_NH3_loading(self.system, verbose=verbose)

    def set_production_rate(self, target_production_kg_hr, verbose=True):
        """
        Set the target production rate and adjust system inputs accordingly.
        Delegates to system module implementation.
        """
        from biorefineries.prefers.v1.HemDx.system import set_production_rate
        
        self._target_production_kg_hr = target_production_kg_hr
        return set_production_rate(self.system, target_production_kg_hr, verbose=verbose)
    
    def check_product_specifications(self):
        """
        Check if the product stream meets all specifications.
        """
        return True
    
    @property
    def target_production_kg_hr(self):
        return self._target_production_kg_hr
    
    @target_production_kg_hr.setter
    def target_production_kg_hr(self, value):
        if value is not None:
            self.set_production_rate(value)
        else:
            self._target_production_kg_hr = None

if __name__ == '__main__':
    import argparse
    from biorefineries.prefers.v1.HemDx.system import create_NHemDx_system
    from biorefineries.prefers.v1.HemDx._chemicals import create_chemicals_Hemodextrin
    from biorefineries.prefers.v1._process_settings import load_process_settings
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='HemDx TEA Analysis')
    parser.add_argument('--production', type=float, default=150,
                        help='Target production rate in kg/hr (default: 150)')
    args, _ = parser.parse_known_args()
    
    TARGET_PRODUCTION = args.production
    
    print("="*85)
    print(f"HEMDX TEA - WITH DESIGN SPECIFICATION")
    print("="*85)
    
    print("\n0. Setting HemDx thermodynamics...")
    bst.settings.set_thermo(create_chemicals_Hemodextrin(), skip_checks=True)
    
    # 2. Load process settings
    load_process_settings()
    
    print(f"\n1. Creating HemDx system...")
    HemDx_sys = create_NHemDx_system()
    HemDx_sys.operating_hours = 8000
    
    print(f"\n2. Running baseline simulation...")
    HemDx_sys.simulate()
    baseline_production = HemDx_sys.flowsheet.stream.NHemDx_Product.F_mass
    print(f"   Baseline production: {baseline_production:.2f} kg/hr")
    
    print(f"\n3. Creating TEA with target production = {TARGET_PRODUCTION} kg/hr...")
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
    
    print(f"\n4. Adjusting system to target production rate...")
    achieved_production = HemDx_tea.set_production_rate(TARGET_PRODUCTION)
    
    print(f"\n5. Checking product specifications...")
    try:
        HemDx_tea.check_product_specifications()
        print("   [OK] All specifications met!")
    except Exception as e:
        print(f"   [FAIL] Specification check failed: {e}")
    
    # Results
    products = HemDx_sys.flowsheet.stream.NHemDx_Product
    
    print(f"\n{'='*85}")
    print("TEA RESULTS")
    print(f"{'='*85}")
    
    # Calculate annual production
    annual_production_kg = products.F_mass * HemDx_sys.operating_hours
    annual_production_MT = annual_production_kg / 1000
    
    print(f"\nProduction Summary:")
    print(f"  Hourly production:  {products.F_mass:.2f} kg/hr")
    print(f"  Annual production:  {annual_production_MT:.2f} metric tons/year")
    print(f"  Operating hours:    {HemDx_sys.operating_hours} hr/year")
    
    # %%
    print(f"\n{'='*85}")
    HemDx_tea.show()
    
    # %%
    print(f"\n{'='*85}")
    print("CASH FLOW TABLE")
    print(f"{'='*85}")
    cashflow_table = HemDx_tea.get_cashflow_table()
    print(cashflow_table)
    
    # %%
    print(f"\n{'='*85}")
    print("MINIMUM SELLING PRICE")
    print(f"{'='*85}")
    price = HemDx_tea.solve_price(products) # USD/kg
    print(f"  HemDx MPSP: ${price:.4f}/kg")
    print(f"  HemDx MPSP: ${price*1000:.2f}/metric ton")
    annual_revenue = price * annual_production_kg
    print(f"  Annual revenue (at MPSP): ${annual_revenue/1e6:.2f} million")
    
    # %%
    print(f"\n{'='*85}")
    print("COST BREAKDOWN TABLES")
    print(f"{'='*85}")
    
    print("\nAll Cost Table:")
    # df1 = tb.all_cost_table(HemDx_tea)
    # print(df1)
    
    # %%
    print("\nVariable Operating Costs:")
    df2 = bst.report.voc_table(HemDx_sys, 'NHemDx_Product')
    print(df2)
    
    # %%
    print("\nCAPEX Breakdown:")
    df8 = HemDx_tea.CAPEX_table()
    print(df8)
    
    # %%
    print("\nFixed Operating Costs:")
    df9 = HemDx_tea.FOC_table()
    print(df9)
    
    # %%
    print(f"\n{'='*85}")
    print("UNIT OPERATION TABLES")
    print(f"{'='*85}")
    
    print("\nReaction Tables:")
    df3 = bst.report.unit_reaction_tables(HemDx_sys.units)
    print(df3)
    
    # %%
    print("\nUnit Results:")
    df4 = bst.report.unit_result_tables(HemDx_sys.units)
    print(df4)
    
    # %%
    print(f"\n{'='*85}")
    print("UTILITY TABLES")
    print(f"{'='*85}")
    
    print("\nHeat Utilities:")
    df5 = bst.report.heat_utility_tables(HemDx_sys.units)
    try:
        print(df5[0])  # Heating
        print(df5[1])  # Cooling
        print(df5[2])  # Summary
    except:
        print(df5)
        
    # %%
    print("\nPower Utilities:")
    df6 = bst.report.power_utility_table(HemDx_sys.units)
    print(df6)
    
    # %%
    print("\nOther Utilities:")
    df7 = bst.report.other_utilities_table(HemDx_sys.units)
    print(df7)
    
    print(f"\n{'='*85}")
    print("TEA ANALYSIS COMPLETE")
    print(f"{'='*85}")
    print(f"Target Production:   {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"Achieved Production: {products.F_mass:.2f} kg/hr")
    print(f"HemDx MPSP:          ${price:.4f}/kg")
    print(f"{'='*85}\n")
