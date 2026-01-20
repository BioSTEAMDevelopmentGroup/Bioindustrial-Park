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

@njit(cache=True)
def generate_DDB_schedule(years):
    val = 1.
    arr = np.ones(years)
    factor = 2. / years
    for i in range(years):
        depreciation = val * factor
        arr[i] = depreciation
        val -= depreciation
    return arr
    
@njit(cache=True)
def generate_SYD_schedule(years):
    digit_sum = years * (years + 1.) * 0.5
    arr = np.ones(years)
    for i in range(years):
        arr[i] = (years - i) / digit_sum
    return arr

class PreFerSTEA(bst.TEA):

    def __init__(self, system, IRR, duration, depreciation, income_tax,
                operating_days, lang_factor, construction_schedule, WC_over_FCI,
                labor_cost, fringe_benefits, property_tax,
                property_insurance, supplies, maintenance, administration,
                target_production_kg_hr=None):
        # Huang et. al. does not take into account financing or startup
        # so these parameters are 0 by default
        super().__init__(system, IRR, duration, depreciation, income_tax,
                        operating_days, lang_factor, construction_schedule,
                        startup_months=6, startup_FOCfrac=0.9, startup_VOCfrac=0.4,
                        startup_salesfrac=0.4, finance_interest=0.05, finance_years=10,
                        finance_fraction=0.7, WC_over_FCI=WC_over_FCI)
        self.labor_cost = labor_cost
        self.fringe_benefits = fringe_benefits
        self.property_tax = property_tax
        self.property_insurance = property_insurance
        self.supplies= supplies
        self.maintenance = maintenance
        self.administration = administration
        
        # Store target production rate
        self._target_production_kg_hr = target_production_kg_hr

    depreciation_schedules: dict[tuple[str, int], NDArray[float]] = {

        ('MACRS', 3): np.array([.3333, .4445, .1481, .0741]),

        ('MACRS', 5): np.array([.2000, .3200, .1920,
                                .1152, .1152, .0576]),

        ('MACRS', 7):  np.array([.1429, .2449, .1749,
                                 .1249, .0893, .0892,
                                 .0893, .0446]),
        
        ('MACRS', 10): np.array([.1000, .1800, .1440,
                                 .1152, .0922, .0737,
                                 .0655, .0655, .0656,
                                 .0655, .0328]),
      
        ('MACRS', 15): np.array([.0500, .0950, .0855,
                                 .0770, .0693, .0623,
                                 .0590, .0590, .0591,
                                 .0590, .0591, .0590,
                                 .0591, .0590, .0591,
                                 .0295]),
      
        ('MACRS', 20): np.array([0.03750, 0.07219, 0.06677,
                                 0.06177, 0.05713, 0.05285,
                                 0.04888, 0.04522, 0.04462,
                                 0.04461, 0.04462, 0.04461,
                                 0.04462, 0.04461, 0.04462,
                                 0.04461, 0.04462, 0.04461,
                                 0.04462, 0.04461, 0.02231]),
                            
        ('IRAS', 1): np.array([1.0]),

        ('IRAS', 2): np.array([0.75, 0.25]),

        ('IRAS', 3): np.array([0.3333, 0.3333, 0.3334]),

        # This method involves a 20% initial allowance in the first year, with the remaining 80% written off over the asset's life.
        # The prescribed working lives are simplified into 6, 12, or 16 years.
        ('IRAS', 6): np.array([0.3333, 0.1333, 0.1333, 0.1333, 0.1334, 0.1334]),

        ('IRAS', 12): np.array([0.2667, 0.0667, 0.0667, 0.0667, 0.0667, 0.0667,
                            0.0667, 0.0667, 0.0667, 0.0667, 0.0667, 0.0666]),

        ('IRAS', 16): np.array([0.25, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                            0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
    }

    @property
    def depreciation(self) -> str|NDArray[float]:
        """
        Depreciation schedule array or a string with format '{schedule}{years}', 
        where years is the number of years until the property value is zero 
        and schedule is one of the following: 'MACRS' (Modified Accelerated Cost Recovery System), 
        'SL' (straight line), 'DDB' (double declining balance), or 
        'SYD' (sum-of-the-years' digits). If years is not given, it defaults 
        to the venture years at run time.
        """
        return self._depreciation
    @depreciation.setter
    def depreciation(self, depreciation):
        if isinstance(depreciation, str):
            self._depreciation_key = self._depreciation_key_from_name(depreciation)
            self._depreciation = depreciation
        else:
            try:
                self._depreciation = np.array(depreciation, dtype=float)
            except:
                raise TypeError(
                   f"invalid depreciation type '{type(depreciation).__name__}'; "
                    "depreciation must be either an array or a string"
                ) from None
            else:
                self._depreciation_key = None
    
    @classmethod
    def _depreciation_key_from_name(cls, name):
        for prefix in ('MACRS', 'SL', 'DDB', 'SYD', 'IRAS'):
            if name.startswith(prefix):
                years = name[len(prefix):]
                key = (prefix, int(years) if years else None)
                if prefix == 'MACRS' and key not in cls.depreciation_schedules:
                    raise ValueError(
                        f"depreciation name {repr(name)} has a valid "
                         "format, but is not yet implemented in BioSTEAM"
                    )
                return key
        raise ValueError(
               f"invalid depreciation name {repr(name)}; "
                "name must have format '{schedule}{years}', "
                "where years is the number of years until the property value is zero "
                "and schedule is one of the following: 'MACRS' (Modified Accelerated Cost Recovery System), "
                "'SL' (straight line), 'DDB' (double declining balance), or "
                "'SYD' (sum-of-the-years' digits)"
            )

    @classmethod
    def _depreciation_array_from_key(cls, key):
        depreciation_schedules = cls.depreciation_schedules
        if key in depreciation_schedules:
            return depreciation_schedules[key]
        else:
            schedule, years = key
            if schedule == 'SL':
                arr = np.full(years, 1./years)
            elif schedule == 'DDB':
                arr = generate_DDB_schedule(years)
            elif schedule == 'SYD':
                arr = generate_SYD_schedule(years)
            else: # pragma: no cover
                raise RuntimeError(f'unknown depreciation schedule {repr(schedule)}')
            depreciation_schedules[key] = arr
            return arr

    # The abstract _FOC method should take fixed capital investment
    # and return the fixed operating cost.
    def _FOC(self, FCI):
        return (self.purchase_cost*(self.property_tax + self.property_insurance
                    + self.maintenance + self.administration)
                + self.labor_cost*(1+self.fringe_benefits+self.supplies))
    
    def CAPEX_table(self):
        purchase_cost = self.purchase_cost /1e6
        lang_factor = self.lang_factor
        if lang_factor is None:
            FCI = purchase_cost
        else:
            FCI = lang_factor * purchase_cost
        working_capital = FCI * self.WC_over_FCI
        TCI = FCI + working_capital
        accounting = self.Accounting(units='MM$')
        accounting.entry('Purchase cost', purchase_cost)
        accounting.entry('Fixed capital investment (FCI)', FCI, f"Purchase cost x {lang_factor}")
        accounting.entry('Working capital',  working_capital, f"{self.WC_over_FCI:.1%} of FCI")
        accounting.entry('Total capital investment (TCI)',  TCI, "Working capital + FCI")
        return accounting.table()
    
    def FOC_table(self):
        accounting = self.Accounting(units='MM$ / yr')
        FCI = self.purchase_cost / 1e6
        labor_cost = self.labor_cost / 1e6
        labor_burden = self.fringe_benefits + self.supplies
        accounting.entry('Labor salary', np.array(labor_cost))
        accounting.entry('Labor burden', np.array(labor_cost * labor_burden), f'{labor_burden:.0%} of labor salary')
        accounting.entry('Maintenance', np.array(self.maintenance * FCI), f'{self.maintenance:.0%} of FCI')
        accounting.entry('Administration', np.array(self.administration * FCI), f'{self.administration:.1%} of FCI')
        accounting.entry('Property tax', np.array(self.property_tax * FCI), f'{self.property_tax:.1%} of FCI')
        accounting.entry('Property insurance', np.array(self.property_insurance * FCI), f'{self.property_insurance:.1%} of FCI')
        return accounting.table()
    
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

    import biosteam as bst
    from biorefineries.prefers.v1.LegHb.system._config2 import create_LegHb_system
    from biorefineries.prefers.v1.LegHb._chemicals import create_chemicals_LegHb
    from biorefineries.prefers.v1._process_settings import load_process_settings
    from biosteam import settings
    
    print("="*85)
    print("LEGHEMOGLOBIN TEA - WITH DESIGN SPECIFICATION")
    print("="*85)
    
    # 1. FORCE the correct thermodynamics FIRST (before any system creation)
    print("\n0. Setting LegHb thermodynamics...")
    bst.settings.set_thermo(create_chemicals_LegHb(), skip_checks=True)
    
    # 2. Load process settings
    load_process_settings()
    
    # 3. Create system
    print("\n1. Creating LegHb system...")
    LegHb_sys = create_LegHb_system()
    LegHb_sys.operating_hours = 8000
    # multiple production rate
    n=1
    # Define target production rate
    TARGET_PRODUCTION = 275 * n  # kg/hr
    
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
        print("   ✓ All specifications met!")
    except ValueError as e:
        print(f"   ✗ Specification check failed: {e}")
    
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

# %%
