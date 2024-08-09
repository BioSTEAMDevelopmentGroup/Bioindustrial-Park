#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 16:38:37 2024

@author: wenjun
"""

from biosteam import TEA
import biosteam as bst
import thermosteam as tmo
from biorefineries.cornstover import CellulosicEthanolTEA



__all__ = ('SAF_TEA', 'create_SAF_tea', 'ATJ_TEA', 'create_ATJ_tea', 'SAF_coprocessing_TEA', 'create_SAF_coprocessing_tea')


class SAF_TEA(CellulosicEthanolTEA):
    _TCI_ratio_cached = 1




def create_SAF_tea(sys, OSBL_units=None):
    if OSBL_units is None: OSBL_units = bst.get_OSBL(sys.cost_units)
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
    except:
        BT = None
    SAF_tea = SAF_TEA(
        system=sys,
        IRR=0.10,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        duration=(2023, 2053),  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        depreciation='MACRS7', 
        income_tax=0.35, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        operating_days=0.9*365,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        lang_factor=None, 
        construction_schedule=(0.08, 0.60, 0.32),  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        startup_months=3, 
        startup_FOCfrac=1, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        startup_salesfrac=0.5, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        startup_VOCfrac=0.75,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        WC_over_FCI=0.05,
        finance_interest=0.08, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        finance_years=10,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        finance_fraction=0.4,
        OSBL_units=OSBL_units,
        warehouse=0.04, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        site_development=0.09,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        additional_piping=0.045, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        proratable_costs=0.10, 
        field_expenses=0.10, 
        construction=0.20,
        contingency=0.10, 
        other_indirect_costs=0.10,
        labor_cost=2.4e6, # =90*0.9*365*24*33.64 more workers than 60 workers in 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks'(Ling Tao); 33.64 is employment cost of average 2023 from https://data.bls.gov/cgi-bin/srgate
        labor_burden=0.90,
        property_insurance=0.007, 
        maintenance=0.03,
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=BT,
        )
    return SAF_tea

#%%

class SAF_Coprocessing_TEA(TEA):
    
    __slots__ = ('OSBL_units', 'warehouse', 'site_development',
                 'additional_piping', 'proratable_costs', 'field_expenses',
                 'construction', 'contingency', 'other_indirect_costs', 
                 'labor_cost', 'labor_burden', 'property_insurance',
                 'maintenance', '_ISBL_DPI_cached', '_FCI_cached',
                 '_utility_cost_cached', '_steam_power_depreciation',
                 '_steam_power_depreciation_array',
                 'boiler_turbogenerator', '_DPI_cached',
                 'steam_distribution', 'water_supply_cooling_pumping', 
                 'water_distribution', 'electric_substation_and_distribution',
                 'gas_supply_and_distribution', 'comminication', 
                 'safety_installation', 'building', 'yard_works', 'land')
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI,  finance_interest,
                 finance_years, finance_fraction, OSBL_units, warehouse,
                 site_development, additional_piping, labor_cost, labor_burden,
                 property_insurance, maintenance, steam_power_depreciation,
                 boiler_turbogenerator, 
                 steam_distribution, 
                 water_supply_cooling_pumping, water_distribution, 
                 electric_substation_and_distribution,
                 gas_supply_and_distribution,
                 comminication, safety_installation,
                 building, yard_works, land, contingency_new, sanitary_waste_disposal,ESCCL):
                
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months, startup_FOCfrac, startup_VOCfrac,
                         startup_salesfrac, WC_over_FCI,  finance_interest,
                         finance_years, finance_fraction)
        self.OSBL_units = OSBL_units
        self.warehouse = warehouse
        self.site_development = site_development
        self.additional_piping = additional_piping
        self.labor_cost = labor_cost
        self.labor_burden = labor_burden
        self.property_insurance = property_insurance
        self.maintenance = maintenance
        self.steam_power_depreciation = steam_power_depreciation
        self.boiler_turbogenerator = boiler_turbogenerator
        self.steam_distribution = steam_distribution
        self.water_supply_cooling_pumping = water_supply_cooling_pumping
        self.water_distribution = water_distribution
        self.electric_substation_and_distribution = electric_substation_and_distribution
        self.gas_supply_and_distribution = gas_supply_and_distribution
        self.comminication = comminication
        self.safety_installation = safety_installation
        self.building = building
        self.yard_works = yard_works
        self.land = land
        self.contingency_new = contingency_new
        self.sanitary_waste_disposal = sanitary_waste_disposal
        self.ESCCL = ESCCL
    
    
    
    
    @property
    def DPI(self) -> float:
        """Direct permanent investment [USD]."""
        return self._DPI(self.installed_equipment_cost, self.purchase_cost)
    @property
    def TDC(self) -> float:
        """Total depreciable capital [USD]."""
        return self._TDC(self.DPI, self.installed_equipment_cost)
    
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
    
    @property
    def OSBL_purchase_equipment_cost(self):
        return sum([i.purchase_cost for i in self.OSBL_units])
            
                                                                                                                                  
    def _fill_depreciation_array(self, D, start, years, TDC):
        depreciation_array = self._get_depreciation_array()
        N_depreciation_years = depreciation_array.size
        if N_depreciation_years > years:
            raise RuntimeError('depreciation schedule is longer than plant lifetime')
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
                raise RuntimeError('steam power depreciation schedule is longer than plant lifetime')
            D[start:start + N_depreciation_years] += BT_TDC * depreciation_array
    
    def _ISBL_DPI(self, installed_equipment_cost):
        """Direct permanent investment of units inside battery limits."""
        if self.lang_factor:
            raise NotImplementedError('lang factor cannot yet be used')
        else:
            self._ISBL_DPI_cached = installed_equipment_cost - self.OSBL_installed_equipment_cost
        return self._ISBL_DPI_cached
    
    
    def _ISBL_purchase_equipment_cost(self, purchase_cost):
        """Purchase cost of units inside battery limits."""
        return purchase_cost - self.OSBL_purchase_equipment_cost
    
    
    def _service_facilites_installed_cost(self, installed_equipment_cost):
        """Service facilities installed costs"""
        return self._ISBL_DPI(installed_equipment_cost) * (self.steam_distribution + self.water_supply_cooling_pumping + self.water_distribution +\
                                                           self.electric_substation_and_distribution + self.gas_supply_and_distribution + \
                                                           self.comminication + self.safety_installation + self.sanitary_waste_disposal)
    


    def _building_cost(self, purchase_cost):
        """Building cost replaces warehouse and construction cost """
        return self.building * self._ISBL_purchase_equipment_cost(purchase_cost)
    
    
    
    def _yard_works(self, installed_equipment_cost):
        return self.yard_works * self._ISBL_DPI(installed_equipment_cost)
    
    
    def _land(self, installed_equipment_cost):
        return self.land * self._ISBL_DPI(installed_equipment_cost)
    
    
    def _DPI(self, installed_equipment_cost, purchase_cost): # Direct Permanent Investment
        factors = self.warehouse + self.site_development + self.additional_piping
        self._DPI_cached = DPI = installed_equipment_cost + self._ISBL_DPI(installed_equipment_cost) * factors +\
            self._service_facilites_installed_cost(installed_equipment_cost) + self._building_cost(purchase_cost) + self._yard_works(installed_equipment_cost)+\
                self._land(installed_equipment_cost)
        return DPI
    
    # Indirect cost
    def _engineering_supervision_construction_contractor_legal_costs(self, installed_equipment_cost):
        return self.ESCCL * self._ISBL_DPI(installed_equipment_cost)
        
    def _contingency(self, installed_equipment_cost):
        return self.contingency_new * self._ISBL_DPI(installed_equipment_cost)
        
    def _TDC(self, DPI, installed_equipment_cost): # Total Depreciable Capital
        return DPI + self._depreciable_indirect_costs(installed_equipment_cost)
    
    # def _nondepreciable_indirect_costs(self, DPI):
    #     return DPI * self.other_indirect_costs
    
    def _depreciable_indirect_costs(self, installed_equipment_cost):
        return self._engineering_supervision_construction_contractor_legal_costs(installed_equipment_cost) + self._contingency(installed_equipment_cost)
    
    def _FCI(self, TDC): # Fixed Capital Investment
        self._FCI_cached = FCI = TDC 
        # + self._nondepreciable_indirect_costs(self._DPI_cached)
        return FCI
    
    def _FOC(self, FCI): # Fixed Operating Costs
        return (FCI * self.property_insurance
                + self._ISBL_DPI_cached * self.maintenance
                + self.labor_cost * (1 + self.labor_burden))
    






class SAF_coprocessing_TEA(SAF_Coprocessing_TEA):
    _TCI_ratio_cached = 1




def create_SAF_coprocessing_tea(sys,
                                steam_distribution, water_supply_cooling_pumping, water_distribution, 
                                electric_substation_and_distribution,
                                gas_supply_and_distribution,
                                comminication, 
                                safety_installation,
                                building, yard_works, contingency_new, land, labor_cost, sanitary_waste_disposal,
                                OSBL_units=None,):
    if OSBL_units is None: OSBL_units = bst.get_OSBL(sys.cost_units)            
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
    except:
        BT = None
    SAF_tea = SAF_coprocessing_TEA(
        system=sys,
        IRR=0.10,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        duration=(2023, 2053),  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        depreciation='MACRS7', 
        income_tax=0.35, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        operating_days=0.9*365,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        lang_factor=None, 
        construction_schedule=(0.08, 0.60, 0.32),  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        startup_months=3, 
        startup_FOCfrac=1, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        startup_salesfrac=0.5, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        startup_VOCfrac=0.75,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        WC_over_FCI=0.05,
        finance_interest=0.08, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        finance_years=10,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        finance_fraction=0.4,
        OSBL_units=OSBL_units,
        warehouse=0.0, 
        site_development=0.0,
        additional_piping=0.045, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        # proratable_costs=0.0, 
        # field_expenses=0.0, 
        # construction=0.0,
        # contingency=0.0, 
        # other_indirect_costs=0.00,
        labor_burden=0.90,
        property_insurance=0.007, 
        maintenance=0.03,
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=BT,
        steam_distribution=steam_distribution,
        water_supply_cooling_pumping=water_supply_cooling_pumping, 
        water_distribution=water_distribution,
        electric_substation_and_distribution=electric_substation_and_distribution,
        gas_supply_and_distribution=gas_supply_and_distribution,
        comminication=comminication, 
        safety_installation=safety_installation,
        building=building,
        yard_works=yard_works,
        contingency_new=contingency_new,
        land=land,
        labor_cost=labor_cost,
        sanitary_waste_disposal = sanitary_waste_disposal,
        ESCCL=0.89, # 0.32+0.34+0.19+0.04
        )
    return SAF_tea


#%%

# class ATJ_TEA(SAF_Coprocessing_TEA):
#     def __init__(self, system, IRR, duration, depreciation, income_tax,
#                  operating_days, lang_factor, construction_schedule,
#                  startup_months, startup_FOCfrac, startup_VOCfrac,
#                  startup_salesfrac, WC_over_FCI,  finance_interest,
#                  finance_years, finance_fraction, OSBL_units, warehouse,
#                  site_development, additional_piping, proratable_costs,
#                  field_expenses, construction, contingency,
#                  other_indirect_costs, labor_cost, labor_burden,
#                  property_insurance, maintenance, steam_power_depreciation,
#                  boiler_turbogenerator, unit_transportation_cost):
#         super().__init__(system, IRR, duration, depreciation, income_tax,
#                          operating_days, lang_factor, construction_schedule,
#                          startup_months, startup_FOCfrac, startup_VOCfrac,
#                          startup_salesfrac, WC_over_FCI,  finance_interest,
#                          finance_years, finance_fraction)
#         self.OSBL_units = OSBL_units
#         self.warehouse = warehouse
#         self.site_development = site_development
#         self.additional_piping = additional_piping
#         self.proratable_costs = proratable_costs
#         self.field_expenses = field_expenses
#         self.construction = construction
#         self.contingency = contingency
#         self.other_indirect_costs = other_indirect_costs
#         self.labor_cost = labor_cost
#         self.labor_burden = labor_burden
#         self.property_insurance = property_insurance
#         self.maintenance = maintenance
#         self.steam_power_depreciation = steam_power_depreciation
#         self.boiler_turbogenerator = boiler_turbogenerator
#         self.unit_transportation_cost = unit_transportation_cost # $/kg ethanol 
        
#     @property
#     def steam_power_depreciation(self):
#         """[str] 'MACRS' + number of years (e.g. 'MACRS7')."""
#         return self._steam_power_depreciation
#     @steam_power_depreciation.setter
#     def steam_power_depreciation(self, depreciation):
#         self._steam_power_depreciation_array = self._depreciation_array_from_key(
#             self._depreciation_key_from_name(depreciation)
#         )
#         self._steam_power_depreciation = depreciation
    
#     @property
#     def ISBL_installed_equipment_cost(self):
#         return self._ISBL_DPI(self.DPI)
    
#     @property
#     def OSBL_installed_equipment_cost(self):
#         if self.lang_factor:
#             raise NotImplementedError('lang factor cannot yet be used')
#         elif isinstance(self.system, bst.AgileSystem):
#             unit_capital_costs = self.system.unit_capital_costs
#             OSBL_units = self.OSBL_units
#             return sum([unit_capital_costs[i].installed_cost for i in OSBL_units])
#         else:
#             return sum([i.installed_cost for i in self.OSBL_units])
    
#     @property
#     def eth_trans_cost(self):
#         ethanol = self.system.flowsheet.ethanol
#         return ethanol.F_mass * self.unit_transportation_cost * self.operating_hours
        
#     @property
#     def VOC(self) -> float:
#         """Variable operating costs [USD/yr]."""
#         return self.material_cost + self.utility_cost + self._eth_trans_cost
    
#     def _fill_depreciation_array(self, D, start, years, TDC):
#         depreciation_array = self._get_depreciation_array()
#         N_depreciation_years = depreciation_array.size
#         if N_depreciation_years > years:
#             raise RuntimeError('depreciation schedule is longer than plant lifetime')
#         system = self.system
#         BT = self.boiler_turbogenerator
#         if BT is None:
#             D[start:start + N_depreciation_years] = TDC * depreciation_array
#         else:
#             if isinstance(system, bst.AgileSystem): BT = system.unit_capital_costs[BT]
#             BT_TDC = BT.installed_cost 
#             D[start:start + N_depreciation_years] = (TDC - BT_TDC) * depreciation_array
            
#             depreciation_array = self._steam_power_depreciation_array
#             N_depreciation_years = depreciation_array.size
#             if N_depreciation_years > years:
#                 raise RuntimeError('steam power depreciation schedule is longer than plant lifetime')
#             D[start:start + N_depreciation_years] += BT_TDC * depreciation_array
    
#     def _ISBL_DPI(self, installed_equipment_cost):
#         """Direct permanent investment of units inside battery limits."""
#         if self.lang_factor:
#             raise NotImplementedError('lang factor cannot yet be used')
#         else:
#             self._ISBL_DPI_cached = installed_equipment_cost - self.OSBL_installed_equipment_cost
#         return self._ISBL_DPI_cached
        
#     def _DPI(self, installed_equipment_cost): # Direct Permanent Investment
#         factors = self.warehouse + self.site_development + self.additional_piping
#         self._DPI_cached = DPI = installed_equipment_cost + self._ISBL_DPI(installed_equipment_cost) * factors
#         return DPI
    
#     def _TDC(self, DPI): # Total Depreciable Capital
#         return DPI + self._depreciable_indirect_costs(DPI)
    
#     def _nondepreciable_indirect_costs(self, DPI):
#         return DPI * self.other_indirect_costs
    
#     def _depreciable_indirect_costs(self, DPI):
#         return DPI * (self.proratable_costs + self.field_expenses
#                       + self.construction + self.contingency)
    
#     def _FCI(self, TDC): # Fixed Capital Investment
#         self._FCI_cached = FCI = TDC + self._nondepreciable_indirect_costs(self._DPI_cached)
#         return FCI
    
#     def _FOC(self, FCI): # Fixed Operating Costs
#         return (FCI * self.property_insurance
#                 + self._ISBL_DPI_cached * self.maintenance
#                 + self.labor_cost * (1 + self.labor_burden))
    
    
    
    
    
# class ATJ_TEA(ATJ_TEA):
#     _TCI_ratio_cached = 1





# def create_ATJ_tea(sys, unit_transportation_cost, OSBL_units=None):
#     if OSBL_units is None: OSBL_units = bst.get_OSBL(sys.cost_units)
#     try:
#         BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
#     except:
#         BT = None
#     SAF_tea = ATJ_TEA(
#         system=sys,
#         IRR=0.10,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         duration=(2023, 2053),  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         depreciation='MACRS7', 
#         income_tax=0.35, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         operating_days=0.9*365,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         lang_factor=None, 
#         construction_schedule=(0.08, 0.60, 0.32),  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         startup_months=3, 
#         startup_FOCfrac=1, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         startup_salesfrac=0.5, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         startup_VOCfrac=0.75,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         WC_over_FCI=0.05,
#         finance_interest=0.08, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         finance_years=10,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         finance_fraction=0.4,
#         OSBL_units=OSBL_units,
#         warehouse=0.04, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         site_development=0.09,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         additional_piping=0.045, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
#         proratable_costs=0.10, 
#         field_expenses=0.10, 
#         construction=0.20,
#         contingency=0.10, 
#         other_indirect_costs=0.10,
#         labor_cost=2.4e6, # =90*0.9*365*24*33.64 more workers than 60 workers in 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks'(Ling Tao); 33.64 is employment cost of average 2023 from https://data.bls.gov/cgi-bin/srgate
#         labor_burden=0.90,
#         property_insurance=0.007, 
#         maintenance=0.03,
#         steam_power_depreciation='MACRS20',
#         boiler_turbogenerator=BT,
#         unit_transportation_cost = unit_transportation_cost
#         )
#     return SAF_tea
