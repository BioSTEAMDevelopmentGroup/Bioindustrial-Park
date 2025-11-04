# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 15:18:58 2025

@author: IGB
"""


import numpy as np
from biosteam import TEA
import biosteam as bst
import thermosteam as tmo

#%%
class CellulosicEthanolTEA(TEA):
    
    __slots__ = ('OSBL_units', 'warehouse', 'site_development',
                 'additional_piping', 'proratable_costs', 'field_expenses',
                 'construction', 'contingency', 'other_indirect_costs', 
                 'labor_cost_hourly', 'labor_burden', 'property_insurance',
                 'maintenance', '_ISBL_DPI_cached', '_FCI_cached',
                 '_utility_cost_cached', '_steam_power_depreciation',
                 '_steam_power_depreciation_array',
                 'boiler_turbogenerator', '_DPI_cached')
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI,  finance_interest,
                 finance_years, finance_fraction, OSBL_units, warehouse,
                 site_development, additional_piping, proratable_costs,
                 field_expenses, construction, contingency,
                 other_indirect_costs, labor_cost_hourly, labor_burden,
                 property_insurance, maintenance, steam_power_depreciation,
                 boiler_turbogenerator):
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
        self.labor_cost_hourly = labor_cost_hourly
        self.labor_burden = labor_burden
        self.property_insurance = property_insurance
        self.maintenance = maintenance
        self.steam_power_depreciation = steam_power_depreciation
        self.boiler_turbogenerator = boiler_turbogenerator
        
    @property
    def labor_cost(self):
        return 2 * self.labor_cost_hourly * self.operating_days * 24
    
    @property
    def steam_power_depreciation(self):
        """[str] 'MACRS' + number of years (e.g. 'MACRS7')."""
        return self._steam_power_depreciation
    @steam_power_depreciation.setter
    def steam_power_depreciation(self, depreciation):
        self._steam_power_depreciation_array = self._depreciation_array_from_key(
            self._depreciation_key_from_name(depreciation)
        )
        self._steam_power_depreciation = depreciation
    
    @property
    def ISBL_installed_equipment_cost(self):
        return self._ISBL_DPI(self.installed_equipment_cost)
    
    @property
    def OSBL_installed_equipment_cost(self):
        if self.lang_factor:
            raise NotImplementedError('lang factor cannot yet be used')
        elif isinstance(self.system, bst.AgileSystem):
            unit_capital_costs = self.system.unit_capital_costs
            OSBL_units = self.OSBL_units
            return sum([unit_capital_costs[i].installed_cost for i in OSBL_units])
        else:
            return sum([i.installed_cost for i in self.OSBL_units])
    
    def _fill_depreciation_array(self, D, start, years, TDC):
        depreciation_array = self._get_depreciation_array()
        N_depreciation_years = depreciation_array.size
        if N_depreciation_years > years:
            dummy = depreciation_array[:years]
            dummy[-1] = depreciation_array[years-1:].sum()
            depreciation_array = dummy
        system = self.system
        BT = self.boiler_turbogenerator
        if BT is None:
            D[start:start + N_depreciation_years] = TDC * depreciation_array
        else:
            if isinstance(system, bst.AgileSystem): BT = system.unit_capital_costs[BT]
            BT_TDC = BT.installed_cost 
            D[start:start + N_depreciation_years] = (TDC - BT_TDC) * depreciation_array
            depreciation_array = self._steam_power_depreciation_array
            N_depreciation_years = depreciation_array.size
            if N_depreciation_years > years:
                dummy = depreciation_array[:years]
                dummy[-1] = depreciation_array[years-1:].sum()
                depreciation_array = dummy
            D[start:start + N_depreciation_years] += BT_TDC * depreciation_array
    
    def _ISBL_DPI(self, installed_equipment_cost):
        """Direct permanent investment of units inside battery limits."""
        if self.lang_factor:
            raise NotImplementedError('lang factor cannot yet be used')
        else:
            self._ISBL_DPI_cached = installed_equipment_cost - self.OSBL_installed_equipment_cost
        return self._ISBL_DPI_cached
        
    def _DPI(self, installed_equipment_cost): # Direct Permanent Investment
        factors = self.warehouse + self.site_development + self.additional_piping
        self._DPI_cached = DPI = installed_equipment_cost + self._ISBL_DPI(installed_equipment_cost) * factors
        return DPI
    
    def _TDC(self, DPI): # Total Depreciable Capital
        return DPI + self._depreciable_indirect_costs(DPI)
    
    def _nondepreciable_indirect_costs(self, DPI):
        return DPI * self.other_indirect_costs
    
    def _depreciable_indirect_costs(self, DPI):
        return DPI * (self.proratable_costs + self.field_expenses
                      + self.construction + self.contingency)
    
    def _FCI(self, TDC): # Fixed Capital Investment
        self._FCI_cached = FCI = TDC + self._nondepreciable_indirect_costs(self._DPI_cached)
        return FCI
    
    def _FOC(self, FCI): # Fixed Operating Costs
        return (FCI * self.property_insurance
                + self._ISBL_DPI_cached * self.maintenance
                + self.labor_cost * (1 + self.labor_burden))
    
    def CAPEX_table(self):
        return capex_table(self)
    
    def FOC_table(self):
        return foc_table(self)

def capex_table(teas, names=None, dataframe=True):
    if isinstance(teas, bst.TEA): teas = [teas]
    if names is None: 
        if len(teas) == 1:
            names = None
        else:
            names = [i.system.ID for i in teas]
    tea, *_ = teas
    capex = tea.Accounting('MM$', names=names)
    ISBL_installed_equipment_costs = np.array([i.ISBL_installed_equipment_cost / 1e6 for i in teas])
    OSBL_installed_equipment_costs = np.array([i.OSBL_installed_equipment_cost / 1e6 for i in teas])
    capex.entry(('Direct costs', 'ISBL installed equipment cost'), ISBL_installed_equipment_costs)
    capex.entry(('Direct costs', 'OSBL installed equipment cost'), OSBL_installed_equipment_costs)
    ISBL_factor_entry = lambda name, value: capex.entry(name, ISBL_installed_equipment_costs * value, f"{value:.1%} of ISBL")
    ISBL_factor_entry(('Direct costs', 'Warehouse'), tea.warehouse)
    ISBL_factor_entry(('Direct costs', 'Site development'), tea.site_development)
    ISBL_factor_entry(('Direct costs', 'Additional piping'), tea.additional_piping)
    TDC = np.array(capex.total_costs)
    capex.entry(('Total direct cost (TDC)', ''), TDC)
    TDC_factor_entry = lambda name, value: capex.entry(name, TDC * value, f"{value:.1%} of TDC")
    TDC_factor_entry(('Indirect costs', 'Proratable costs'), tea.proratable_costs)
    TDC_factor_entry(('Indirect costs', 'Field expenses'), tea.field_expenses)
    TDC_factor_entry(('Indirect costs', 'Construction'), tea.construction)
    TDC_factor_entry(('Indirect costs', 'Contingency'), tea.contingency)
    TDC_factor_entry(('Indirect costs', 'Other (start-up, permits, etc.)'), tea.other_indirect_costs)
    TIC = np.array(capex.total_costs) - 2 * TDC
    capex.entry(('Total indirect cost (TIC)', ''), TIC)
    FCI = TDC + TIC
    capex.entry(('Fixed capital investment (FCI)', ''), FCI, 'TDC + TIC')
    working_capital = FCI * tea.WC_over_FCI
    capex.entry(('Working capital (WC)', ''), working_capital, f"{tea.WC_over_FCI:.1%} of FCI")
    TCI = FCI + working_capital
    capex.entry(('Total capital investment (TCI)', ''), TCI, 'FCI + WC')
    return capex.table() if dataframe else capex

voc_table = bst.report.voc_table

def foc_table(teas, names=None, dataframe=True):
    if isinstance(teas, bst.TEA): teas = [teas]
    tea, *_ = teas
    if names is None: 
        if len(teas) == 1: 
            names = None
        else:
            names = [i.system.ID for i in teas]
    accounting = tea.Accounting('MM$ / yr', names=names)
    ISBL = np.array([i.ISBL_installed_equipment_cost / 1e6 for i in teas])
    labor_cost = np.array([i.labor_cost / 1e6 for i in teas])
    FCI = np.array([i.FCI / 1e6 for i in teas])
    accounting.entry('Labor salary', labor_cost)
    accounting.entry('Labor burden', tea.labor_burden * labor_cost, '90% of labor salary')
    accounting.entry('Maintenance', tea.maintenance * ISBL, f'{tea.maintenance:.1%} of ISBL')
    accounting.entry('Property insurance', tea.property_insurance * FCI, f'{tea.property_insurance:.1%} of FCI')
    return accounting.table() if dataframe else accounting


#%%

class Plasma_TEA(CellulosicEthanolTEA):
    _TCI_ratio_cached = 1


def create_plasma_tea(sys, IRR, OSBL_units=None):
    if OSBL_units is None: OSBL_units = bst.get_OSBL(sys.cost_units)
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
    except:
        BT = None
    tea_plasma = Plasma_TEA(
        system=sys,
        IRR=IRR,
        duration=(2023, 2043),
        depreciation='MACRS7', 
        income_tax=0,
        operating_days=0.9*365,
        lang_factor=None, 
        construction_schedule=(1,0),
        startup_months=0, 
        startup_FOCfrac=1, 
        startup_salesfrac=1,
        startup_VOCfrac=1,
        WC_over_FCI=0.0,
        finance_interest=0.0,
        finance_years=0,
        finance_fraction=0.,
        OSBL_units=OSBL_units,
        warehouse=0,
        site_development=0.0,
        additional_piping=0.0,
        proratable_costs=0.0, 
        field_expenses=0.0, 
        construction=0.0,
        contingency=0., 
        other_indirect_costs=0.0,
        labor_cost_hourly=0,
        labor_burden=0,
        property_insurance=0.00, 
        maintenance=0.0,
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=BT,
        )
    return tea_plasma
        
