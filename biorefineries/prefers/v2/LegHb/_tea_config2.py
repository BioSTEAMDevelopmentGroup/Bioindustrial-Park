# -*- coding: utf-8 -*-
"""
TEA Configuration for LegHb system using _config2 (Titer = 3.5 g/L).

Same as _tea_config1.py but imports from _config2 instead of _config1.
Config 2 uses identical flowsheet with reduced titer (3.5 g/L vs 5.0 g/L).

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from warnings import filterwarnings
filterwarnings('ignore')

import biosteam as bst
import numpy as np
import pandas as pd
from biorefineries.tea.conventional_ethanol_tea import *
from numba import njit
from numpy import ndarray as NDArray
from biorefineries.prefers.v2._tea import PreFerSTEA as PreFerSTEA_Base
import flexsolve as flx
import warnings
# Import seed_targets from the Config 2 module
from biorefineries.prefers.v2.LegHb.system._config2 import seed_targets

class PreFerSTEA(PreFerSTEA_Base):
    
    def set_production_rate(self, target_production_kg_hr, verbose=True):
        """
        Set the target production rate and adjust system inputs accordingly.
        Delegates to _config2.set_production_rate.
        """
        from biorefineries.prefers.v2.LegHb.system._config2 import set_production_rate
        
        self._target_production_kg_hr = target_production_kg_hr
        return set_production_rate(self.system, target_production_kg_hr, verbose=verbose)
    
    def check_product_specifications(self):
        """
        Check if the product stream meets all specifications.
        """
        from biorefineries.prefers.v2.LegHb.system._config2 import check_LegHb_specifications
        
        # Get the product stream (LegHb_3)
        product_stream = self.system.flowsheet.stream.LegHb_3
        
        return check_LegHb_specifications(product_stream)
    
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

# %%
if __name__ == '__main__':

    import argparse
    import biosteam as bst
    from biorefineries.prefers.v2.LegHb.system._config2 import create_LegHb_system
    from biorefineries.prefers.v2.LegHb._chemicals import create_chemicals_LegHb
    from biorefineries.prefers.v2._process_settings import load_process_settings
    from biosteam import settings
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='LegHemoglobin TEA Analysis (Config 2 - Titer 3.5 g/L)')
    parser.add_argument('--production', type=float, default=150,
                        help='Target production rate in kg/hr (default: 150)')
    args, _ = parser.parse_known_args()
    
    TARGET_PRODUCTION = args.production
    
    print("="*85)
    print(f"LEGHEMOGLOBIN TEA - CONFIG2 (Titer = 3.5 g/L)")
    print("="*85)
    
    # 1. FORCE the correct thermodynamics FIRST
    print("\n0. Setting LegHb thermodynamics...")
    bst.settings.set_thermo(create_chemicals_LegHb(), skip_checks=True)
    
    # 2. Load process settings
    load_process_settings()
    
    # 3. Create system using config2 (R302 spec handles titer + NH3 internally)
    print(f"\n1. Creating LegHb system (config2 - titer 3.5 g/L)...")
    LegHb_sys = create_LegHb_system()
    LegHb_sys.operating_hours = 8000
    
    # Multiple production rate multiplier
    n = 1
    TARGET_PRODUCTION = TARGET_PRODUCTION * n
    
    print(f"\n2. Running baseline simulation...")
    LegHb_sys.simulate()
    baseline_production = LegHb_sys.flowsheet.stream.LegHb_3.F_mass
    print(f"   Baseline production: {baseline_production:.2f} kg/hr")
    
    # Create TEA object
    print(f"\n3. Creating TEA with target production = {TARGET_PRODUCTION} kg/hr...")
    LegHb_tea = PreFerSTEA(
        system=LegHb_sys, 
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
    
    # Set production rate using the TEA method
    print(f"\n4. Adjusting system to target production rate...")
    achieved_production = LegHb_tea.set_production_rate(TARGET_PRODUCTION)
    
    # Verify specifications
    print(f"\n5. Checking product specifications...")
    try:
        LegHb_tea.check_product_specifications()
        print("   [OK] All specifications met!")
    except ValueError as e:
        print(f"   [FAIL] Specification check failed: {e}")
    
    # Get product stream
    products = LegHb_sys.flowsheet.stream.LegHb_3
    
    print(f"\n{'='*85}")
    print("TEA RESULTS")
    print(f"{'='*85}")
    
    # Calculate annual production
    annual_production_kg = products.F_mass * LegHb_sys.operating_hours
    annual_production_MT = annual_production_kg / 1000
    
    print(f"\nProduction Summary:")
    print(f"  Hourly production:  {products.F_mass:.2f} kg/hr")
    print(f"  Annual production:  {annual_production_MT:.2f} metric tons/year")
    print(f"  Operating hours:    {LegHb_sys.operating_hours} hr/year")
    
    # %%
    print(f"\n{'='*85}")
    LegHb_tea.show()
    
    # %%
    print(f"\n{'='*85}")
    print("CASH FLOW TABLE")
    print(f"{'='*85}")
    cashflow_table = LegHb_tea.get_cashflow_table()
    print(cashflow_table)
    
    # %%
    print(f"\n{'='*85}")
    print("MINIMUM SELLING PRICE")
    print(f"{'='*85}")
    price = LegHb_tea.solve_price(products)  # USD/kg
    print(f"  LegHb MPSP: ${price:.4f}/kg")
    print(f"  LegHb MPSP: ${price*1000:.2f}/metric ton")
    annual_revenue = price * annual_production_kg
    print(f"  Annual revenue (at MPSP): ${annual_revenue/1e6:.2f} million")
    
    # %%
    print(f"\n{'='*85}")
    print("COST BREAKDOWN TABLES")
    print(f"{'='*85}")
    
    print("\nVariable Operating Costs:")
    df2 = bst.report.voc_table(LegHb_sys, 'LegHb_3')
    print(df2)
    
    print("\nCAPEX Breakdown:")
    df8 = LegHb_tea.CAPEX_table()
    print(df8)
    
    print("\nFixed Operating Costs:")
    df9 = LegHb_tea.FOC_table()
    print(df9)
    
    # %%
    print(f"\n{'='*85}")
    print("UNIT OPERATION TABLES")
    print(f"{'='*85}")
    
    print("\nReaction Tables:")
    df3 = bst.report.unit_reaction_tables(LegHb_sys.units)
    print(df3)
    
    print("\nUnit Results:")
    df4 = bst.report.unit_result_tables(LegHb_sys.units)
    print(df4)
    
    # %%
    print(f"\n{'='*85}")
    print("UTILITY TABLES")
    print(f"{'='*85}")
    
    print("\nHeat Utilities:")
    df5 = bst.report.heat_utility_tables(LegHb_sys.units)
    print(df5[0])
    print(df5[1])
    print(df5[2])
    df55 = pd.concat([df5[0], df5[1], df5[2]], axis=0, ignore_index=True)
    
    print("\nPower Utilities:")
    df6 = bst.report.power_utility_table(LegHb_sys.units)
    print(df6)
    
    print("\nOther Utilities:")
    df7 = bst.report.other_utilities_table(LegHb_sys.units)
    print(df7)
    
    print(f"\n{'='*85}")
    print("TEA ANALYSIS COMPLETE (Config 2 - Titer 3.5 g/L)")
    print(f"{'='*85}")
    print(f"Target Production:   {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"Achieved Production: {products.F_mass:.2f} kg/hr")
    print(f"LegHb MPSP:          ${price:.4f}/kg")
    print(f"{'='*85}\n")
