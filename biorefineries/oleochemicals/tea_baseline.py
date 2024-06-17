# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 16:17:35 2022 
@author: Lavanya
"""

import biosteam as bst

#Modelled based on assumptions for a nth degree plant
class TEA_baseline(bst.TEA):
    def __init__(self, system, IRR, duration, 
                 depreciation,startup_months,
                 startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac,
                 finance_interest,
                 finance_years, finance_fraction,
                 income_tax,operating_days, lang_factor,
                 construction_schedule,
                 WC_over_FCI,
                 operating_labor_cost,
                 direct_supervisory_clerical_labor,
                 factor_maintenance_and_repairs,
                 factor_for_contingency_fees,
                 location_factor,
                 factor_operating_supplies,          
                 factor_laboratory_charges,
                 factor_local_taxes_and_insurance,
                 OSBL_units,
                 factor_for_site_prep,
                 factor_for_building_construction,
                 factor_for_land_purchase,
                 factor_for_royalties,
                 factor_for_plant_startup,
                 factor_plant_overhead_fci,
                 factor_plant_overhead_lc,
                 factor_administration_costs_fci,
                 factor_administration_costs_lc):
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months, startup_FOCfrac, startup_VOCfrac,
                         startup_salesfrac, WC_over_FCI,  finance_interest,
                         finance_years, finance_fraction)        
        self.IRR = IRR 
        self.duration = duration 
        self.depreciation = depreciation
        self.income_tax = income_tax 
        self.operating_days 
#Direct manufacturing costs include operating_labor_cost,direct_supervisory_clerical_labor,maintenance_and_repairs,operating_supplies,laboratory_charges [1]  
        self.location_factor = location_factor
        self.direct_supervisory_clerical_labor = direct_supervisory_clerical_labor 
        self.factor_maintenance_and_repairs = factor_maintenance_and_repairs
        self.factor_for_contingency_fees = factor_for_contingency_fees 
        self.factor_operating_supplies = factor_operating_supplies 
        self.factor_laboratory_charges = factor_laboratory_charges 
        self.operating_labor_cost = operating_labor_cost 
#Fixed manufacturing costs(depreciation,local taxes and insurance and plant overhead)     
        self.factor_local_taxes_and_insurance = factor_local_taxes_and_insurance  
#General manufacturing expenses (admin costs, distribution and selling costs and RnD) 
        self.OSBL_units = OSBL_units 
        self.factor_for_site_prep =factor_for_site_prep
        self.factor_for_building_construction =factor_for_building_construction
        self.factor_for_land_purchase = factor_for_land_purchase
        self.factor_for_royalties = factor_for_royalties
        self.factor_for_plant_startup = factor_for_plant_startup
        self.factor_plant_overhead_fci =factor_plant_overhead_fci
        self.factor_plant_overhead_lc =factor_plant_overhead_lc
        self.factor_administration_costs_fci = factor_administration_costs_fci
        self.factor_administration_costs_lc = factor_administration_costs_lc
        
        
#The total installed equipment cost is also called the total bare module cost
#All assumptions for DPI are based on Warren Sieder's chapter 16 [1]
    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost*(1 + self.factor_for_site_prep + self.factor_for_building_construction)
#All assumptions for TDC are based on Warren Sieder's chapter 16     
    def _TDC(self, DPI):
        self.total_depreciable_capital = DPI*(1 + self.factor_for_contingency_fees)
        return self.total_depreciable_capital
#ALl assumps for FCI are based on Warren Sieder's chp 16[1]
    def _FCI(self, TDC):
        total_depreciable_capital = TDC
        Sum = total_depreciable_capital*(1 + self.factor_for_land_purchase+ self.factor_for_royalties+ self.factor_for_plant_startup)
        FCI = self.location_factor*Sum
        return FCI 

    def _FOC(self, FCI):
#All assumptions for FOC based on [2]      
        return (FCI*(self.factor_local_taxes_and_insurance +
                     self.factor_maintenance_and_repairs + 
                     self.factor_operating_supplies +
                     self.factor_plant_overhead_fci +
                     self.factor_administration_costs_fci)
               + self.operating_labor_cost*(1 + self.direct_supervisory_clerical_labor + 
                                       self.factor_laboratory_charges +
                                       self.factor_administration_costs_lc +
                                       self.factor_plant_overhead_lc))
    
#[1] Warren Sieder
#[2] (Turton et al., 2013)