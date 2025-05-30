# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
import pandas as pd

__all__ = (
    'ConventionalEthanolTEA', 
    'create_conventional_ethanol_tea', 
)

class ConventionalEthanolTEA(bst.TEA):
    """
    Create a ConventionalEthanolTEA object for techno-economic analysis of a
    biorefinery [1]_.
    
    Parameters
    ----------    
    system : System
        Should contain feed and product streams.
    IRR : float
        Internal rate of return (fraction).
    duration : tuple[int, int]
        Start and end year of venture (e.g. (2018, 2038)).
    depreciation : str
        'MACRS' + number of years (e.g. 'MACRS7').
    operating_days : float
        Number of operating days per year.
    income_tax : float
        Combined federal and state income tax rate (fraction).
    lang_factor : float
        Lang factor for getting fixed capital investment from
        total purchase cost. If no lang factor, estimate capital investment
        using bare module factors.
    startup_schedule : tuple[float]
        Startup investment fractions per year 
        (e.g. (0.5, 0.5) for 50% capital investment in the first year and 50%
        investment in the second).
    WC_over_FCI : float
        Working capital as a fraction of fixed capital investment.
    labor_cost : float
        Total labor cost (USD/yr).
    fringe_benefits : float
        Cost of fringe benefits as a fraction of labor cost.
    property_tax : float
        Fee as a fraction of fixed capital investment.
    property_insurance : float
        Fee as a fraction of fixed capital investment.    
    supplies : float
        Yearly fee as a fraction of labor cost.
    maintenance : float
        Yearly fee as a fraction of fixed capital investment.
    administration : float
        Yearly fee as a fraction of fixed capital investment.

    References
    ----------
    .. [1] Huang, H., Long, S., & Singh, V. (2016). Techno-economic analysis of biodiesel
        and ethanol co-production from lipid-producing sugarcane. Biofuels, Bioproducts
        and Biorefining, 10(3), 299–315. https://doi.org/10.1002/bbb.1640
    
    """
    __slots__ = ('labor_cost', 'fringe_benefits', 'maintenance',
                 'property_tax', 'property_insurance', '_FCI_cached',
                 'supplies', 'maintenance', 'administration')
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule, WC_over_FCI,
                 labor_cost, fringe_benefits, property_tax,
                 property_insurance, supplies, maintenance, administration):
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
    
    def _FCI(self, TDC):
        self._FCI_cached = TDC
        return TDC
    
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
        accounting.entry('Labor salary', labor_cost)
        accounting.entry('Labor burden', labor_cost * labor_burden, f'{labor_burden:.0%} of labor salary')
        accounting.entry('Maintenance', self.maintenance * FCI, f'{self.maintenance:.0%} of FCI')
        accounting.entry('Administration', self.administration * FCI, f'{self.administration:.1%} of FCI')
        accounting.entry('Property tax', self.property_tax * FCI, f'{self.property_tax:.1%} of FCI')
        accounting.entry('Property insurance', self.property_insurance * FCI, f'{self.property_insurance:.1%} of FCI')
        return accounting.table()
        

def create_conventional_ethanol_tea(
        system, cls=ConventionalEthanolTEA, IRR=0.15,
        duration=(2018, 2038), depreciation='MACRS7', income_tax=0.35,
        operating_days=200, lang_factor=3, construction_schedule=(0.4, 0.6), 
        WC_over_FCI=0.05, labor_cost=2.5e6, fringe_benefits=0.4,
        property_tax=0.001, property_insurance=0.005, supplies=0.20, 
        maintenance=0.01, administration=0.005
    ): # Values default to sugarcane biorefinery
    return cls(system, 
        IRR=IRR,
        duration=duration,
        depreciation=depreciation, income_tax=income_tax,
        operating_days=operating_days, lang_factor=lang_factor,
        construction_schedule=construction_schedule, WC_over_FCI=WC_over_FCI,
        labor_cost=labor_cost, fringe_benefits=fringe_benefits,
        property_tax=property_tax, property_insurance=property_insurance,
        supplies=supplies, maintenance=maintenance, administration=administration
    )
