# -*- coding: utf-8 -*-
"""
Created on 2026-01-28 09:10:00

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import biosteam as bst
import numpy as np
from biosteam import TEA
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
