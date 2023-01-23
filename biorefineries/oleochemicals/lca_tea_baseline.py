# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 16:17:35 2022 

@author: Lavanya
"""

import biosteam as bst
from chaospy import distributions as shape
from biorefineries.oleochemicals import systems_baseline 
#TODO: what to do about prices
##Prices from the literature
##Ref: https://doi.org/10.1016/j.indcrop.2020.112411
#Raw material prices
# #Given std deviation: 0.34, mean = 6.8
# HOSO_price_USD_per_ton = shape.normal(mu = 6.8, sigma = 0.34)
# ##Given std deviation: 5.8, mean = 0.302
# Methanol_price_USD_per_ton = shape.normal(mu = 5.8, sigma = 0.302)
# ##Given std deviation: 92, mean = 975
# H2O2_USD_per_ton = shape.normal(mu = 975, sigma = 92)
##Given std deviation: 1000, mean = 33000
# Tungstic_acid_USD_per_ton = shape.normal(mu = 3000, sigma = 1000)
##Given std deviation: 500, mean = 10918
# Cobalt_acetate_USD_per_ton = shape.normal(mu = 10918, sigma = 500)
#Given triangular distribution a:100, b: 500, c:200
# Crude_glycerol_USD_per_ton = shape.Triangle(100, 500 ,200)
#Product_prices
# ##Given std deviation: 0.11, mean = 8.82
# Azelaic_acid_USD_per_ton =  shape.normal(mu = 8.82, sigma = 0.11)
# ##Given std deviation: 0.2, mean = 8.2
# Pelargonic_acid_USD_per_ton = shape.normal(mu = 8.2, sigma = 0.2)

#This plant is a grass roots plant

class TEA_baseline(bst.TEA):
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule, WC_over_FCI,
                 labor_cost,  property_tax,
                 property_insurance, supplies, maintenance, 
                 administration, operating_supervision, laboratory_charges,
                 plant_overhead,OSBL_units
                 ):
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months=0, startup_FOCfrac=0, startup_VOCfrac=0,
                         startup_salesfrac=0, finance_interest=0, finance_years=0,
                         finance_fraction=0, WC_over_FCI=WC_over_FCI)
        
        self.labor_cost = labor_cost #The average labor cost of one operator is 60,000 $/y,
                                     #giving a total of 3.6 m $/y for the labor costs (LC).
        self.property_tax = property_tax #2% of TCC
        self.property_insurance = property_insurance #2% of TCC        
        self.laboratory_charges = laboratory_charges
        self.operating_supervision = operating_supervision
        self.supplies= supplies # 1% of TCC
        self.maintenance = maintenance #2% of TCC
        self.administration = administration #0.2*self.labor_cost
        self.plant_overhead = plant_overhead #0.6*self.labor_cost
        self.IRR = IRR #10%
        self.income_tax = income_tax #35%
        self.duration = duration #TODO: not sure
        self.OSBL_units = OSBL_units
        

    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost + 0.2*installed_equipment_cost + 0.1*installed_equipment_cost 
        # installed_equipment_cost + "C_site is 20% of TBM"  0.2*installed_equipment_cost + "C_buildings is 10% of TBM" 0.1*installed_equipment_cost "C_serv is NA" #C_alloc #C_comp #C_offsite facilities #TODO: ask YOEL
    

    def _TDC(self, DPI):
        #1.18% = DPI + 0.03*DPI + 0.15*DPI
        self.total_depreciable_capital = 1.18* DPI
        return self.total_depreciable_capital

    def _FCI(self, TDC):
        total_depreciable_capital = TDC
        F_site = 1.15 #For US mid
        FCI = F_site* total_depreciable_capital + 0.02*total_depreciable_capital + 0.02* total_depreciable_capital + 0.2*total_depreciable_capital
        # FCI = F_site*(total_depreciable_capital + "C_land is 2% TDC=" 0.02*total_depreciable_capital + "initial C_royalties is 2% of TDC" 0.02*total_depreciable_capital + "C_plant_startup is 10% to 30%" 0.2*total_depreciable_capital )
    #The fixed capital investment is also called the  total permanent investment    
        return FCI 

    def _FOC(self, FCI):
        return (FCI*(self.property_tax + self.property_insurance
                     + self.maintenance + self.supplies)
                + self.labor_cost*(1+self.operating_supervision +self.laboratory_charges +self.plant_overhead +
                                   self.administration +self.plant_overhead                                  
                                   ))
    
azelaic_baseline = TEA_baseline(system = systems_baseline.azelaic_acid_sys,
                                operating_days = 200,
                                IRR = 0.1,duration=(2013,2023),
                                depreciation = 'MACRS7',
                                lang_factor = 3,income_tax = 0.35,
                                construction_schedule = (0.4,0.6),
                                WC_over_FCI=0.05,
                                labor_cost = 3600000,
                                property_tax = 0.02,
                                property_insurance = 0.02,
                                supplies = 0.01,
                                maintenance = 0.02,
                                administration = 0.02,
                                laboratory_charges = 0.18,
                                operating_supervision = 0.18,
                                plant_overhead = 0.6,
                                OSBL_units = [systems_baseline.azelaic_acid_sys.facilities])    
    
    