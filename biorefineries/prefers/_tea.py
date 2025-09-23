# -*- coding: utf-8 -*-
"""
Created on 2025-07-02 15:20:03

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biorefineries.corn import tea
import biosteam as bst
import numpy as np
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
                property_insurance, supplies, maintenance, administration):
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
        return (FCI*(self.property_tax + self.property_insurance
                    + self.maintenance + self.administration)
                + self.labor_cost*(1+self.fringe_benefits+self.supplies))
    
    def CAPEX_table(self):
        purchase_cost = self.purchase_cost /1e6
        lang_factor = self.lang_factor
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
        FCI = self.FCI / 1e6
        labor_cost = self.labor_cost / 1e6
        labor_burden = self.fringe_benefits + self.supplies
        accounting.entry('Labor salary', np.array(labor_cost))
        accounting.entry('Labor burden', np.array(labor_cost * labor_burden), f'{labor_burden:.0%} of labor salary')
        accounting.entry('Maintenance', np.array(self.maintenance * FCI), f'{self.maintenance:.0%} of FCI')
        accounting.entry('Administration', np.array(self.administration * FCI), f'{self.administration:.1%} of FCI')
        accounting.entry('Property tax', np.array(self.property_tax * FCI), f'{self.property_tax:.1%} of FCI')
        accounting.entry('Property insurance', np.array(self.property_insurance * FCI), f'{self.property_insurance:.1%} of FCI')
        return accounting.table()
    
    # S-curve 0.10, 0.6, 0.3
    # Food tech 0.15，0.60， 0.25
    # Pharmaceutical 0.25, 0.50, 0.25
    # Biochemical 0.10, 0.65, 0.25
    # agri food 0.1 0.35 0.4 0.15

# %%
if __name__ == '__main__':
    import biosteam as bst
    from biorefineries.prefers.systems.LegH.LegH import create_LegH_system
    from biorefineries.prefers._process_settings import load_process_settings
    from biosteam import settings
    load_process_settings()
    legH_sys = create_LegH_system()
    legH_sys.simulate()  # Simulate the system to ensure all units are ready
    legH_tea = PreFerSTEA(
        system=legH_sys, IRR=0.18, duration=(2025, 2045), depreciation='IRAS6',
        income_tax=0.17, operating_days=333, lang_factor=None,
        construction_schedule=(0.15, 0.60, 0.25), WC_over_FCI=0.15,
        labor_cost=3e6, fringe_benefits=0.17, property_tax=0.005,
        property_insurance=0.005, supplies=0.02, maintenance=0.03,
        administration=0.05
    )
    products = legH_sys.flowsheet('LegH_3')
    # %%
    legH_tea.show()  # Display the TEA summary
    # %%
    legH_tea.get_cashflow_table()
    # %%
    cost=legH_tea.solve_price(products) # USD/kg
    print(f"\nLegH price: ${round(cost, 5)}/kg")   
    # %%
    from biorefineries.prefers import _table as tb
    df1 = tb.all_cost_table(legH_tea)
    print(df1)
    # %%
    df2=bst.report.voc_table(legH_sys,'LegH_3')
    print(df2)
    # # %%
    # df3 = bst.report.unit_reaction_tables(legH_sys.units)
    # print(df3)
    # # %%
    # df4 = bst.report.unit_result_tables(legH_sys.units)
    # print(df4)
    # # %%
    # df5 = bst.report.heat_utility_tables(legH_sys.units)
    # print(df5)     # %%
    # df6 = bst.report.power_utility_table(legH_sys.units)
    # print(df6)
    # %%
    # df7 = bst.report.other_utilities_table(legH_sys.units)
    # print(df7)

# %% Purchased equipment costs
    # print('Purchased Costs:')
    # for i in legH_tea.units:
    #     #if i in legH_tea.units: continue
    #     val = i.purchase_cost
    #     print(i.ID, ', ', val)
    # print('Installed Costs:')
    # for i in legH_tea.units:
    #     #if i in legH_tea.units: continue
    #     val = i.installed_cost
    #     print(i.ID, ', ', val)
    # for i in legH_tea.units:
    #     #if i in legH_tea.units: continue
    #     val = i.installed_cost
    #     print(i.ID, ', ', val)
    #     # if val or i.ID=='M':
    #     #     print(i.ID, ', ', val)
# %% Production costs
    
    # # %%
    # costs = legH_tea.production_costs(products)# USD/yr
    # np.round(costs / 1e6) # million USD / yr

    # # %%
    # # Debug: Print system feeds and products to understand structure
    # print("System feeds:")
    # for i, feed in enumerate(legH_sys.feeds):
    #     print(f"  {i}: {feed.ID} ({feed})")
    # print("\nSystem products:")
    # for i, product in enumerate(legH_sys.products):
    #     print(f"  {i}: {product.ID} ({product})")
    
    # # Alternative approach: Find the Glucose feed by name
    # glucose_feed = None
    # for feed in legH_sys.feeds:
    #     if feed.ID == 'Glucose':
    #         glucose_feed = feed
    #         break
    
    # if glucose_feed is None:
    #     # If not found by ID, try to access from system inputs by index
    #     # The Glucose stream is the 3rd input (index 2) in the system inputs
    #     glucose_feed = legH_sys.feeds[2]  # Index 2 corresponds to s.Glucose
    
    # price = legH_tea.solve_price(glucose_feed) # USD/kg
    # print(f"\nGlucose price: ${round(price, 5)}/kg")
    # # %%
    # legH_tea.IRR = legH_tea.solve_IRR()
    # legH_tea.show()
# %%
