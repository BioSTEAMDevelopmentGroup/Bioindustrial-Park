# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 21:48:12 2019

@author: yoelr
"""

from biosteam import TEA
import thermosteam as tmo
import biosteam as bst

__all__ = ('CellulosicEthanolTEA', 'create_tea')

class CellulosicEthanolTEA(TEA):
    
    __slots__ = ('OSBL_units', 'warehouse', 'site_development',
                 'additional_piping', 'proratable_costs', 'field_expenses',
                 'construction', 'contingency', 'other_indirect_costs', 
                 'labor_cost', 'labor_burden', 'property_insurance',
                 'maintenance', '_ISBL_DPI_cached', '_FCI_cached',
                 '_utility_cost_cached')
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI,  finance_interest,
                 finance_years, finance_fraction, OSBL_units, warehouse,
                 site_development, additional_piping, proratable_costs,
                 field_expenses, construction, contingency,
                 other_indirect_costs, labor_cost, labor_burden,
                 property_insurance, maintenance):
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months, startup_FOCfrac, startup_VOCfrac,
                         startup_salesfrac, WC_over_FCI,  finance_interest,
                         finance_years, finance_fraction)
        self.OSBL_units = OSBL_units
        self.warehouse = warehouse
        self.site_development = site_development
        self.additional_piping = additional_piping
        self.proratable_costs = proratable_costs
        self.field_expenses = field_expenses
        self.construction = construction
        self.contingency = contingency
        self.other_indirect_costs = other_indirect_costs
        self.labor_cost = labor_cost
        self.labor_burden = labor_burden
        self.property_insurance = property_insurance
        self.maintenance = maintenance
    
    @property
    def utility_cost(self):
        self._utility_cost_cached = utility_cost = super().utility_cost
        return utility_cost
    
    def _ISBL_DPI(self, DPI):
        """Direct permanent investment of units inside battery limits."""
        if self.lang_factor:
            self._ISBL_DPI_cached = DPI - sum([i.purchase_cost for i in self.OSBL_units]) * self.lang_factor
        else:
            self._ISBL_DPI_cached = DPI - sum([i.installed_cost for i in self.OSBL_units])
        return self._ISBL_DPI_cached
        
    def _TDC(self, DPI):
        return DPI + self._ISBL_DPI(DPI) * (self.warehouse
                                            + self.site_development
                                            + self.additional_piping)
    
    def _indirect_costs(self, TDC):
        return TDC*(self.proratable_costs + self.field_expenses
                    + self.construction + self.contingency
                    + self.other_indirect_costs)
    
    def _FCI(self, TDC):
        self._FCI_cached = FCI = TDC + self._indirect_costs(TDC)
        return FCI
    
    def _FOC(self, FCI):
        return (FCI * self.property_insurance
                + self._ISBL_DPI_cached * self.maintenance
                + self.labor_cost * (1 + self.labor_burden))
    
def create_tea(sys, OSBL_units, ignored_units=()):
    BT = tmo.utils.get_instance(OSBL_units, bst.facilities.BoilerTurbogenerator)
    OSBL_units = list(OSBL_units)
    OSBL_units.remove(BT)
    ethanol_tea = CellulosicEthanolTEA(
        system=sys, 
        IRR=0.10, 
        duration=(2007, 2037),
        depreciation='MACRS7', 
        income_tax=0.35,
        operating_days=350.4,
        lang_factor=None, 
        construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, 
        startup_FOCfrac=1,
        startup_salesfrac=0.5,
        startup_VOCfrac=0.75,
        WC_over_FCI=0.05,
        finance_interest=0.08,
        finance_years=10,
        finance_fraction=0.4,
        OSBL_units=OSBL_units, # BT not included
        warehouse=0.04, 
        site_development=0.09, 
        additional_piping=0.045,
        proratable_costs=0.10,
        field_expenses=0.10,
        construction=0.20,
        contingency=0.10,
        other_indirect_costs=0.10, 
        labor_cost=2.5e6,
        labor_burden=0.90,
        property_insurance=0.007, 
        maintenance=0.03)
    for i in ignored_units: ethanol_tea.units.remove(i)
    ethanol_tea.units.remove(BT)
    Area700 = bst.TEA.like(bst.System('Area700', (BT,)), ethanol_tea)
    Area700.labor_cost = 0
    Area700.depreciation = 'MACRS20'
    Area700.OSBL_units = (BT,)
    return bst.CombinedTEA([ethanol_tea, Area700], IRR=0.10)