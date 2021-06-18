# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 21:48:12 2019

@author: yoelr
"""

from biosteam import TEA, AgileTEA
import thermosteam as tmo
import biosteam as bst
import pandas as pd
import numpy as np

__all__ = ('CellulosicEthanolTEA', 'AgileCellulosicEthanolTEA', 
           'create_tea', 'create_agile_tea',
           'capex_table', 'voc_table', 'foc_table')

class CAPEXTableBuilder:
    __slots__ = ('index', 'cost', 'notes')
    
    def __init__(self):
        self.index = []
        self.cost = []
        self.notes = []
    
    def entry(self, index: str, cost: float, notes: str = '-'):
        self.index.append(index)
        self.cost.append(cost)
        self.notes.append(notes)
    
    def table(self):
        data = tuple(zip(self.notes, self.cost))
        return pd.DataFrame(data, 
                            index=self.index,
                            columns=('Notes', 'Cost [MM$]')
        )


class VOCTableBuilder:
    __slots__ = ('index', 'price', 'cost', 'operating_days')
    
    def __init__(self, operating_days):
        self.operating_days = operating_days
        self.index = []
        self.price = []
        self.cost = []
        
    def entry(self, stream, name=None, price=None):
        if not name:
            name = stream.ID.replace('_', ' ')
        if name.islower(): name= name.capitalize()
        if price: 
            cost = price * self.operating_days * stream.get_total_flow('ton/day') / 1e6
        else:
            price = stream.price * 907.185
            cost = stream.cost * self.operating_days * 24. / 1e6
        if stream.isfeed():
            index = ('Raw materials', name)
        elif stream.isproduct():
            index = ('By-products and credits', name)
        else:
            raise ValueError('stream must be either a feed or a product')
        self.index.append(index)
        self.price.append(price)
        self.cost.append(cost)
        
    def table(self):
        VOC = 0.
        for index, cost in zip(self.index, self.cost):
            kind = index[0]
            if kind == 'Raw materials':
                VOC += cost
            elif kind == 'By-products and credits':
                VOC += -cost
            else:
                raise RuntimeError(f"invalid index '{kind}'")
        data = (*zip(self.price, self.cost), ('', VOC))
        index = pd.MultiIndex.from_tuples(self.index + [('Total variable operating cost', '')])
        return pd.DataFrame(data, 
                            index=index,
                            columns=('Price [$/ton]', 'Cost [MM$ / yr]'))


class FOCTableBuilder:
    __slots__ = ('operating_days', 'index', 'price', 'notes')
    
    def __init__(self, operating_days):
        self.operating_days = operating_days
        self.index = []
        self.price = []
        self.notes = []
        
    def entry(self, index, cost, notes='-'):
        self.index.append(index)
        self.price.append(cost)
        self.notes.append(notes)
    
    def table(self):
        yearly_cost = np.array(self.price) * self.operating_days / 365.
        data = tuple(zip(self.notes, self.price, yearly_cost))
        return pd.DataFrame(data, 
                            index=self.index,
                            columns=(
                                'Notes',
                                'Price [MM$ / yr]',
                                'Cost [MM$ / yr]',
                            )
        )
    
        
class CellulosicEthanolTEA(TEA):
    
    __slots__ = ('OSBL_units', 'warehouse', 'site_development',
                 'additional_piping', 'proratable_costs', 'field_expenses',
                 'construction', 'contingency', 'other_indirect_costs', 
                 'labor_cost', 'labor_burden', 'property_insurance',
                 'maintenance', '_ISBL_DPI_cached', '_FCI_cached',
                 '_utility_cost_cached', '_steam_power_depreciation',
                 '_steam_power_depreciation_array',
                 'boiler_turbogenerator', '_TCI_ratio_cached')
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI,  finance_interest,
                 finance_years, finance_fraction, OSBL_units, warehouse,
                 site_development, additional_piping, proratable_costs,
                 field_expenses, construction, contingency,
                 other_indirect_costs, labor_cost, labor_burden,
                 property_insurance, maintenance, steam_power_depreciation,
                 boiler_turbogenerator, _TCI_ratio_cached=1.):
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
        self.steam_power_depreciation = steam_power_depreciation
        self.boiler_turbogenerator = boiler_turbogenerator
        self._TCI_ratio_cached = _TCI_ratio_cached
        
    @property
    def steam_power_depreciation(self):
        """[str] 'MACRS' + number of years (e.g. 'MACRS7')."""
        return self._steam_power_depreciation
    @steam_power_depreciation.setter
    def steam_power_depreciation(self, depreciation):
        try:
            self._steam_power_depreciation_array = self.depreciation_schedules[depreciation]
        except KeyError:
            raise ValueError(f"depreciation must be either 'MACRS5', 'MACRS7', 'MACRS10' or 'MACRS15 (not {repr(depreciation)})")
        self._steam_power_depreciation = depreciation
    
    @property
    def ISBL_installed_equipment_cost(self):
        return self._ISBL_DPI(self.DPI)
    
    @property
    def OSBL_installed_equipment_cost(self):
        if self.lang_factor:
            raise NotImplementedError('lang factor cannot yet be used')
        else:
            return sum([i.installed_cost for i in self.OSBL_units])
    
    def _fill_depreciation_array(self, D, start, years, TDC):
        depreciation_array = self._depreciation_array
        N_depreciation_years = depreciation_array.size
        if N_depreciation_years > years:
            raise RuntimeError('depreciation schedule is longer than plant lifetime')
        BT = self.boiler_turbogenerator
        if BT is None:
            try:
                D[start:start + N_depreciation_years] = TDC * depreciation_array
            except:
                print(D)
                print(D[start:start + N_depreciation_years].shape)
                print((TDC * depreciation_array).shape)
                breakpoint()
        else:
            BT_TDC = BT.installed_cost 
            D[start:start + N_depreciation_years] = (TDC - BT_TDC) * depreciation_array
            
            depreciation_array = self._steam_power_depreciation_array
            N_depreciation_years = depreciation_array.size
            if N_depreciation_years > years:
                raise RuntimeError('steam power depreciation schedule is longer than plant lifetime')
            D[start:start + N_depreciation_years] += BT_TDC * depreciation_array
    
    def _ISBL_DPI(self, installed_equipment_cost):
        """Direct permanent investment of units inside battery limits."""
        if self.lang_factor:
            raise NotImplementedError('lang factor cannot yet be used')
        else:
            self._ISBL_DPI_cached = installed_equipment_cost - self.OSBL_installed_equipment_cost
        return self._ISBL_DPI_cached
        
    def _DPI(self, installed_equipment_cost):
        factors = self.warehouse + self.site_development + self.additional_piping
        return installed_equipment_cost + self._ISBL_DPI(installed_equipment_cost) * factors
    
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
                + self.labor_cost * (1 + self.labor_burden)) * self.operating_hours / 8760
    
class AgileCellulosicEthanolTEA(AgileTEA):
    
    __slots__ = CellulosicEthanolTEA.__slots__
    steam_power_depreciation = CellulosicEthanolTEA.steam_power_depreciation
    ISBL_installed_equipment_cost = CellulosicEthanolTEA.ISBL_installed_equipment_cost
    _ISBL_DPI = CellulosicEthanolTEA._ISBL_DPI
    _DPI = CellulosicEthanolTEA._DPI
    _indirect_costs = CellulosicEthanolTEA._indirect_costs
    _FCI = CellulosicEthanolTEA._FCI
    _FOC = CellulosicEthanolTEA._FOC
    
    def __init__(self, IRR, duration, depreciation, income_tax,
                 lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI,  finance_interest,
                 finance_years, finance_fraction, OSBL_units, warehouse,
                 site_development, additional_piping, proratable_costs,
                 field_expenses, construction, contingency,
                 other_indirect_costs, labor_cost, labor_burden,
                 property_insurance, maintenance, steam_power_depreciation,
                 boiler_turbogenerator):
        super().__init__(IRR, duration, depreciation, income_tax,
                         lang_factor, construction_schedule,
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
        self.steam_power_depreciation = steam_power_depreciation
        self.boiler_turbogenerator = boiler_turbogenerator
    
    @property
    def OSBL_installed_equipment_cost(self):
        if self.lang_factor:
            raise NotImplementedError('lang factor cannot yet be used')
        else:
            units = self.units
            OSBL_units = self.OSBL_units
            return sum([i.installed_cost for i in units if i.unit in OSBL_units])
    
    def _fill_depreciation_array(self, D, start, years, TDC):
        depreciation_array = self._depreciation_array
        N_depreciation_years = depreciation_array.size
        if N_depreciation_years > years:
            raise RuntimeError('depreciation schedule is longer than plant lifetime')
        for BT in self.units:
            if BT.unit is self.boiler_turbogenerator: break
        BT_TDC = BT.installed_cost 
        D[start:start + N_depreciation_years] = (TDC - BT_TDC) * depreciation_array
        depreciation_array = self._steam_power_depreciation_array
        N_depreciation_years = depreciation_array.size
        if N_depreciation_years > years:
            raise RuntimeError('steam power depreciation schedule is longer than plant lifetime')
        D[start:start + N_depreciation_years] += BT_TDC * depreciation_array
    
    
def create_agile_tea(units, OSBL_units=None, ignored_units=()):
    for i in ignored_units: units.remove(i)
    if OSBL_units is None: OSBL_units = bst.get_OSBL(units)
    BT = tmo.utils.get_instance(OSBL_units, bst.BoilerTurbogenerator)
    tea = AgileCellulosicEthanolTEA( 
        IRR=0.10, 
        duration=(2007, 2037),
        depreciation='MACRS7', 
        income_tax=0.35,
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
        OSBL_units=OSBL_units, 
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
        maintenance=0.03,
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=BT)
    return tea

def create_tea(sys, OSBL_units=None, ignored_units=(), TEA=None):
    if OSBL_units is None: OSBL_units = bst.get_OSBL(sys.units)
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
    except:
        BT = None
    if TEA is None: TEA = CellulosicEthanolTEA
    tea = TEA(
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
        OSBL_units=OSBL_units,
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
        maintenance=0.03,
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=BT)
    for i in ignored_units: tea.units.remove(i)
    return tea

def capex_table(tea):
    capex = CAPEXTableBuilder()
    ISBL_installed_equipment_cost = tea.ISBL_installed_equipment_cost / 1e6
    OSBL_installed_equipment_cost = tea.OSBL_installed_equipment_cost / 1e6
    capex.entry('ISBL installed equipment cost', ISBL_installed_equipment_cost)
    capex.entry('OSBL installed equipment cost', OSBL_installed_equipment_cost)
    ISBL_factor_entry = lambda name, value: capex.entry(name, ISBL_installed_equipment_cost * value, f"{value:.1%} of ISBL")
    ISBL_factor_entry('Warehouse', tea.warehouse)
    ISBL_factor_entry('Site development', tea.site_development)
    ISBL_factor_entry('Additional piping', tea.additional_piping)
    TDC = sum(capex.cost)
    capex.entry('Total direct cost (TDC)', TDC)
    TDC_factor_entry = lambda name, value: capex.entry(name, TDC * value, f"{value:.1%} of TDC")
    TDC_factor_entry('Proratable costs', tea.proratable_costs)
    TDC_factor_entry('Field expenses', tea.field_expenses)
    TDC_factor_entry('Construction', tea.construction)
    TDC_factor_entry('Contingency', tea.contingency)
    TDC_factor_entry('Other indirect costs (start-up, permits, etc.)', tea.other_indirect_costs)
    TIC = sum(capex.cost) - 2 * TDC
    capex.entry('Total indirect cost', TIC)
    FCI = TDC + TIC
    capex.entry('Fixed capital investment (FCI)', FCI)
    working_capital = FCI * tea.WC_over_FCI
    capex.entry('Working capital', working_capital, f"{tea.WC_over_FCI:.1%} of FCI")
    TCI = FCI + working_capital
    capex.entry('Total capital investment (TCI)', TCI)
    return capex.table()

def voc_table(system, tea, main_products):
    voc = VOCTableBuilder(tea.operating_days)
    for i in system.feeds + system.products: 
        if i in main_products: continue
        if i.price and not i.isempty(): voc.entry(i)
    isa = isinstance
    for i in set(system.facilities):
        if isa(i, bst.BoilerTurbogenerator):
            natural_gas = i.ins[3]
            if natural_gas.isempty(): continue
            voc.entry(natural_gas, price=i.natural_gas_price * 907.185)
    return voc.table()

def foc_table(tea):
    foc = FOCTableBuilder(tea.operating_days)
    tea, _ = tea.TEAs
    ISBL = tea.ISBL_installed_equipment_cost / 1e6
    labor_cost = tea.labor_cost / 1e6
    foc.entry('Labor salary', labor_cost)
    foc.entry('Labor burden', tea.labor_burden * labor_cost, '90% of labor salary')
    foc.entry('Maintenance', tea.maintenance * ISBL, f'{tea.maintenance:.1%} of ISBL')
    foc.entry('Property insurance', tea.property_insurance * ISBL, f'{tea.maintenance:.1%} of ISBL')
    return foc.table()