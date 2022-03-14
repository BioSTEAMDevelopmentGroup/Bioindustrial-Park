# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 16:59:06 2022

@author: lavanyakudli
"""
import biosteam as bst
import numpy as np

class OzonolysisTEA(bst.TEA):
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule, WC_over_FCI,
                 labor_cost, fringe_benefits, property_tax,
                 property_insurance, supplies, maintenance, administration):
        
        # Huang et. al. does not take into account financing or startup
        # so these parameters are 0 by default
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months=0, startup_FOCfrac=0, startup_VOCfrac=0,
                         startup_salesfrac=0, finance_interest=0, finance_years=0,
                         finance_fraction=0, WC_over_FCI=WC_over_FCI)
        self.labor_cost = labor_cost
        self.fringe_benefits = fringe_benefits
        self.property_tax = property_tax
        self.property_insurance = property_insurance
        self.supplies= supplies
        self.maintenance = maintenance
        self.administration = administration

    # The abstract _DPI method should take installed equipment cost
    # and return the direct permanent investment. Huang et. al. assume
    # these values are equal
    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost

    # The abstract _TDC method should take direct permanent investment
    # and return the total depreciable capital. Huang et. al. assume
    # these values are equal
    def _TDC(self, DPI):
        return DPI

    # The abstract _FCI method should take total depreciable capital
    # and return the fixed capital investment. Again, Huang et. al.
    # assume these values are equal.
    def _FCI(self, TDC):
        return TDC

    # The abstract _FOC method should take fixed capital investment
    # and return the fixed operating cost.
    def _FOC(self, FCI):
        return (FCI*(self.property_tax + self.property_insurance
                     + self.maintenance + self.administration)
                + self.labor_cost*(1+self.fringe_benefits+self.supplies))