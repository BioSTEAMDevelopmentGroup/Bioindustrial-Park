# -*- coding: utf-8 -*-
"""
Created on 2025-07-02 15:20:03

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

# NOTE: Removed conflicting import that overwrites thermodynamics
# from biorefineries.corn import tea
import biosteam as bst
import numpy as np
import pandas as pd
from biorefineries.tea.conventional_ethanol_tea import *
# bst.nbtutorial()
from numba import njit
from numpy import ndarray as NDArray
from biorefineries.prefers.v1._tea import PreFerSTEA as PreFerSTEA_Base

class PreFerSTEA(PreFerSTEA_Base):
    
    def set_production_rate(self, target_production_kg_hr):
        """
        Set the target production rate and adjust system inputs accordingly.
        
        Parameters
        ----------
        target_production_kg_hr : float
            Target production rate in kg/hr for LegHb_3 stream
            
        Returns
        -------
        float
            Achieved production rate [kg/hr]
        """
        # Import the function from LegHb module
        from biorefineries.prefers.v1.LegHb.system import set_production_rate
        
        self._target_production_kg_hr = target_production_kg_hr
        
        # Call the design specification function
        achieved_production = set_production_rate(self.system, target_production_kg_hr)
        
        return achieved_production
    
    def check_product_specifications(self):
        """
        Check if the product stream meets all specifications.
        
        Returns
        -------
        bool
            True if all specifications are met
            
        Raises
        ------
        ValueError
            If any specification is not met
        """
        # Import the function from LegHb module
        from biorefineries.prefers.v1.LegHb.system import check_LegHb_specifications
        
        # Get the product stream (LegHb_3)
        product_stream = self.system.flowsheet.stream.LegHb_3
        
        # Check specifications
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
    
    # S-curve 0.10, 0.6, 0.3
    # Food tech 0.15，0.60， 0.25
    # Pharmaceutical 0.25, 0.50, 0.25
    # Biochemical 0.10, 0.65, 0.25
    # agri food 0.1 0.35 0.4 0.15

# %%
if __name__ == '__main__':

    import argparse
    import biosteam as bst
    from biorefineries.prefers.v1.LegHb.system import create_LegHb_system, get_available_configs
    from biorefineries.prefers.v1.LegHb._chemicals import create_chemicals_LegHb
    from biorefineries.prefers.v1._process_settings import load_process_settings
    from biosteam import settings
    
    # Parse command line arguments for config selection
    parser = argparse.ArgumentParser(description='LegHemoglobin TEA Analysis')
    parser.add_argument('--config', type=str, default='config1',
                        choices=get_available_configs(),
                        help='Process configuration to use (default: config1)')
    parser.add_argument('--production', type=float, default=150,
                        help='Target production rate in kg/hr (default: 150)')
    args, _ = parser.parse_known_args()
    
    CONFIG = args.config
    TARGET_PRODUCTION = args.production
    
    print("="*85)
    print(f"LEGHEMOGLOBIN TEA - WITH DESIGN SPECIFICATION (CONFIG: {CONFIG.upper()})")
    print("="*85)
    
    # 1. FORCE the correct thermodynamics FIRST (before any system creation)
    print("\n0. Setting LegHb thermodynamics...")
    bst.settings.set_thermo(create_chemicals_LegHb(), skip_checks=True)
    
    # 2. Load process settings
    load_process_settings()
    
    # 3. Create system with selected configuration
    print(f"\n1. Creating LegHb system (config={CONFIG})...")
    LegHb_sys = create_LegHb_system(config=CONFIG)
    LegHb_sys.operating_hours = 8000
    # multiple production rate
    n=1
    # Define target production rate
    TARGET_PRODUCTION = TARGET_PRODUCTION * n  # kg/hr
    
    print(f"\n2. Running baseline simulation...")
    LegHb_sys.simulate()
    baseline_production = LegHb_sys.flowsheet.stream.LegHb_3.F_mass
    print(f"   Baseline production: {baseline_production:.2f} kg/hr")
    
    # Create TEA object WITH target production rate
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
        target_production_kg_hr=TARGET_PRODUCTION  # SET PRODUCTION RATE HERE
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
    LegHb_tea.show()  # Display the TEA summary
    
    # %%
    print(f"\n{'='*85}")
    print("CASH FLOW TABLE")
    print(f"{'='*85}")
    cashflow_table = LegHb_tea.get_cashflow_table()
    print(cashflow_table)
    #cashflow_table.to_excel('LegHb_cashflow_table.xlsx',index=True)
    
    # %%
    print(f"\n{'='*85}")
    print("MINIMUM SELLING PRICE")
    print(f"{'='*85}")
    price = LegHb_tea.solve_price(products) # USD/kg
    print(f"  LegHb MPSP: ${price:.4f}/kg")
    print(f"  LegHb MPSP: ${price*1000:.2f}/metric ton")
    annual_revenue = price * annual_production_kg
    print(f"  Annual revenue (at MPSP): ${annual_revenue/1e6:.2f} million")
    
    # %%
    print(f"\n{'='*85}")
    print("COST BREAKDOWN TABLES")
    print(f"{'='*85}")
    
    # from biorefineries.prefers import _table as tb
    
    print("\nAll Cost Table:")
    # df1 = tb.all_cost_table(LegHb_tea)
    # print(df1)
    #df1.to_excel('LegHb_cost_table.xlsx',index=True)
    
    # %%
    print("\nVariable Operating Costs:")
    df2 = bst.report.voc_table(LegHb_sys, 'LegHb_3')
    print(df2)
    #df2.to_excel('LegHb_voc_table.xlsx',index=True)
    
    # %%
    print("\nCAPEX Breakdown:")
    df8 = LegHb_tea.CAPEX_table()
    print(df8)
    #df8.to_excel('LegHb_CAPEX_table.xlsx',index=True)
    
    # %%
    print("\nFixed Operating Costs:")
    df9 = LegHb_tea.FOC_table()
    print(df9)
    #df9.to_excel('LegHb_FOC_table.xlsx',index=True)
    
    # %%
    print(f"\n{'='*85}")
    print("UNIT OPERATION TABLES")
    print(f"{'='*85}")
    
    print("\nReaction Tables:")
    df3 = bst.report.unit_reaction_tables(LegHb_sys.units)
    print(df3)
    #df3.to_excel('LegHb_reaction_table.xlsx',index=True)
    
    # %%
    print("\nUnit Results:")
    df4 = bst.report.unit_result_tables(LegHb_sys.units)
    print(df4)
    #df4.to_excel('LegHb_result_table.xlsx',index=True)
    
    # %%
    print(f"\n{'='*85}")
    print("UTILITY TABLES")
    print(f"{'='*85}")
    
    print("\nHeat Utilities:")
    df5 = bst.report.heat_utility_tables(LegHb_sys.units)
    print(df5[0])  # Heating
    print(df5[1])  # Cooling
    print(df5[2])  # Summary
    df55 = pd.concat([df5[0], df5[1], df5[2]], axis=0, ignore_index=True)
    #df55.to_excel('LegHb_heat_utility_table_combined.xlsx',index=True)
    
    # %%
    print("\nPower Utilities:")
    df6 = bst.report.power_utility_table(LegHb_sys.units)
    print(df6)
    #df6.to_excel('LegHb_power_utility_table.xlsx',index=True)
    
    # %%
    print("\nOther Utilities:")
    df7 = bst.report.other_utilities_table(LegHb_sys.units)
    print(df7)
    #df7.to_excel('LegHb_other_utility_table.xlsx',index=True)
    
    print(f"\n{'='*85}")
    print("TEA ANALYSIS COMPLETE")
    print(f"{'='*85}")
    print(f"Target Production:   {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"Achieved Production: {products.F_mass:.2f} kg/hr")
    print(f"LegHb MPSP:          ${price:.4f}/kg")
    print(f"{'='*85}\n")
