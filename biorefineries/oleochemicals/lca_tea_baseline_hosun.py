# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 16:17:35 2022 

@author: Lavanya
"""

import biosteam as bst
import chaospy
from chaospy import distributions as shape
from biorefineries.oleochemicals.systems_baseline_hosun import aa_baseline_sys,F

import numpy as np

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
                 maintenance_and_repairs,
                 operating_supplies,          
                 laboratory_charges,
                 local_taxes_and_insurance,
                 plant_overhead,
                 administration_costs,
                 # distribution_and_selling_costs,
                 OSBL_units,
                 ):
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months, startup_FOCfrac, startup_VOCfrac,
                         startup_salesfrac, WC_over_FCI,  finance_interest,
                         finance_years, finance_fraction)        
        self.IRR = IRR #10%
        self.duration = duration #TODO: not sure
        self.depreciation = depreciation
        self.income_tax = income_tax #35%
        
#Direct manufacturing costs (operating_labor_cost,direct_supervisory_clerical_labor,maintenance_and_repairs,operating_supplies,laboratory_charges)        
#Ref: (Turton et al., 2013), Chapter 18 -  Table 8.2. Multiplication Factors for Estimating Manufacturing Cost        
        self.operating_labor_cost = operating_labor_cost #https://www.bls.gov/iag/tgs/iag325.htm#prices
        #Avg value of labour cost for a chemicals manufacturing plant
        #Mean annual wage is 54,874$
        self.direct_supervisory_clerical_labor = direct_supervisory_clerical_labor #Ref(Turton et al., 2013), 0.18*C_operating_labor #range (0.1-0.25)*C_operating_labor
        self.maintenance_and_repairs = maintenance_and_repairs #Ref(Turton et al., 2013) 0.06*FCI, range:(0.02-0.1)*FCI
        self.operating_supplies = operating_supplies #Ref(Turton et al., 2013) 0.009*FCI range:(0.1-0.2)*FCI
        self.laboratory_charges = laboratory_charges #Ref(Turton et al., 2013) 0.15*C_operating_labor, range:(0.1-0.2)*C_operating_labor
        #TODO: depreciation 0.1*FCI?
#Fixed manufacturing costs(depreciation,local taxes and insurance and plant overhead)
#Ref: (Turton et al., 2013), Chapter 18 -  Table 8.2. Multiplication Factors for Estimating Manufacturing Cost        
        self.local_taxes_and_insurance = local_taxes_and_insurance  #0.1*FCI  
        self.plant_overhead = plant_overhead #0.708*C_operating_labor + 0.036*FCI
#General manufacturing expenses (admin costs, distribution and selling costs and RnD)                
        self.administration_costs = administration_costs #0.177*C_operating_labour + 0.009*FCI
        #self.distribution_and_selling_costs = distribution_and_selling_costs#0.11*COM, function of manufacturing cost
        #Same with research and develop       
        self.OSBL_units = OSBL_units #Based on the 
        
#The total installed equipment cost is also called the total bare module cost
#All assumptions for DPI are based on Warren Sieder's chapter 16
    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost + 0.04*installed_equipment_cost + 0.1*installed_equipment_cost 
        # installed_equipment_cost +
        # + "C_site is 10 -20% of TBM" (for grass roots)  and C_site is 4-6% for nth plant
        # + "C_buildings is 10% of TBM" 
        # + "C_spare" = ?
        # + "C_storage" = covered in the process simulation
        # + "C_catalyst" = covered in the process simulation
        # + "C_comp" = ?
        # + "C_serv" = ?
        # + "C_alloc" : to include or upgrade utilities
        # + "C_offsite" facilities ? 
        
#All assumptions for TDC are based on Warren Sieder's chapter 16     
    def _TDC(self, DPI):
        #"Contigency fees": 0.18*DPI + 0.03*DPI + 0.15*DPI 
        self.total_depreciable_capital = 1.18* DPI
        return self.total_depreciable_capital
#ALl assumps for FCI are based on Warren Sieder's chp 16
    def _FCI(self, TDC):
        total_depreciable_capital = TDC
        F_site = 1.15 #For US mid
        Sum = total_depreciable_capital + 0.02*total_depreciable_capital + 0.02* total_depreciable_capital + 0.1*total_depreciable_capital
        FCI = F_site*Sum
        # FCI = F_site*(total_depreciable_capital 
        #               + "C_land" =  is 2% TDC= " 0.02*total_depreciable_capital 
        #               + "C_royalties" = initial C_royalties is 2% of TDC" =  0.02*total_depreciable_capital 
        #               + "C_plant_startup" = typicaly 10% = 0.1*total_depreciable_capital )
    #The fixed capital investment is also called the  total permanent investment    
        return FCI 

    def _FOC(self, FCI):
        return (FCI*(self.local_taxes_and_insurance +
                     self.maintenance_and_repairs + 
                     self.operating_supplies +
                     self.plant_overhead*0.0508 +
                     self.administration_costs*0.0508
                     )
               + self.operating_labor_cost*(1 + self.direct_supervisory_clerical_labor + 
                                       self.laboratory_charges +
                                       self.administration_costs +
                                       self.plant_overhead                                  
                                   ))
    
#TODO: what to do about prices  
# #Given std deviation: 0.34, mean = 6.8
# HOSO_price_USD_per_ton = shape.normal(mu = 6.8, sigma = 0.34)

#