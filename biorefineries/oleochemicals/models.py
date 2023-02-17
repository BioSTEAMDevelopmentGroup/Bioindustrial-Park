# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:45:33 2023
@author: Lavanya
"""
import biosteam as bst
import chaospy
from chaospy import distributions as shape
from biorefineries.oleochemicals import systems_baseline 
from systems_baseline import aa_baseline_sys
from systems_baseline  import F_baseline
import numpy as np
from lca_tea_baseline import tea_azelaic_baseline
from biosteam.evaluation import Model, Metric

#TODO: which reactor to use titer on and why
aa_sys_op_hours = aa_baseline_sys.operating_hours = tea_azelaic_baseline.operating_days * 24
#functions to solve for different indicators
# Minimum product selling price of lactic_acid stream
azelaic_acid = F_baseline.stream.azelaic_acid_product_stream
# funcs = create_funcs(lactic_tea=lactic_sys.TEA, flowsheet=flowsheet)
azelaic_acid_tea = aa_baseline_sys.TEA
def get_MPSP():
        azelaic_acid.price = 0
        for i in range(3):
            MPSP = azelaic_acid.price = azelaic_acid_tea.solve_price(azelaic_acid)
        return MPSP

crude_HOSO = F_baseline.stream.crude_vegetable_oil
# Yield in 10^6 kg/yr
get_annual_factor = lambda: azelaic_acid_tea.operating_days*24
get_total_yield = lambda: azelaic_acid.F_mass/1e6
# Purity (%) of LacticAcid in the final product
get_purity = lambda: azelaic_acid.imass['Azelaic_acid']/azelaic_acid.F_mass
get_overall_TCI = lambda: azelaic_acid_tea.TCI/1e6
get_azelaic_acid_sale = lambda: get_total_yield*azelaic_acid.price
get_material_cost = lambda: azelaic_acid_tea.material_cost/1e6 #TODO: why should ash dsposal be subtracted?
get_operating_cost = lambda: azelaic_acid_tea.AOC/1e6 #TODO: why to sub? -get_gypsum_sale()-get_ash_sale()
##### Utilities related #####
get_system_heating_demand = lambda: aa_baseline_sys.get_heating_duty()/1e9
get_system_cooling_water_duty = lambda: aa_baseline_sys.get_cooling_duty()/1e9
# To see if TEA converges well for each simulation #TODO: ask why?
get_NPV = lambda: azelaic_acid_tea.NPV
metrics = [    
               # Metric('MPSP', get_MPSP, '$/kg'),
               # Metric('Total product yield', get_total_yield, '10^6 kg/yr'),
               # Metric('Product purity', get_purity, '%'),
               # Metric('Total capital investment', get_overall_TCI, '10^6 $'),
               # Metric('Annual operating cost', get_operating_cost, '10^6 $/yr'),
               # Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
               # Metric('Annual product sale', get_azelaic_acid_sale, '10^6 $/yr'),
                Metric('Total', get_system_heating_demand, '10^6 MJ/yr', ' Heating demand'),
                Metric('Total', get_system_cooling_water_duty, '10^6 MJ/yr', 'Cooling demand'),
               # Metric('NPV', get_NPV, '$', 'TEA')               
               ]
model = Model(aa_baseline_sys, metrics)
#TODO: ask what check functions are in lactic acid models

#parameters
#Process related parameters
#1 Transesterification reaction conversion
#TODO: or you could acess this F_baseline.X_tes, what would the element be?
# lb1 = 0.90
# ub1 = 0.99
# @model.parameter(element=F_baseline.U1001,
#                  kind='coupled',
#                  distribution=shape.Uniform(lb1,ub1))
# def set_transesterification_conversion_1(X_tes):
#     for i in range (0,11):
#         F_baseline.U1001.transesterification.X[i] = X_tes
        
# @model.parameter(element=F_baseline.U1002,
#                  kind='coupled',
#                  distribution=shape.Uniform(lb1,ub1))
# def set_transesterification_conversion_2(X_tes):
#     for i in range (0,11):
#         F_baseline.U1002.transesterification.X[i] = X_tes         

#2 Dihydroxylation reaction conversion
# F_baseline.DihydroxylationReactor.R101
# @model.parameter(element=F_baseline.U1002,
#                  kind='coupled',
#                  distribution=shape.Uniform(lb1,ub1))
# def set_transesterification_conversion_2(X_tes):
#     for i in range (0,11):
#         F_baseline.U1002.transesterification.X[i] = X_tes
  
# #3 Oxidative cleavage reaction conversion
# F_baseline.OxidativeCleavageReactor.X_oxidativecleavage = 0.8
# F_baseline.OxidativeCleavageReactor.X_decarboxylation = 0.2
# F_baseline.OxidativeCleavageReactor.X_side_rxn = 0.8
#############################################################################3
#4 Pressure for degassing_the_oily_phase
lb4 = 10000
ub4 = 20000
@model.parameter(element=F_baseline.F2001,
                 kind='coupled',
                 distribution=shape.Uniform(lb4,ub4))
def set_degassingcolumn_P(P_d):
    F_baseline.unit.F2001.P = P_d
#5 nonanoic_acid_fraction_separation
#Setting the lighter key recovery
lb_Lr_D501 = 0.9
ub_Lr_D501 = 0.99
@model.parameter(element=F_baseline.D501,
                 kind='coupled',
                 distribution=shape.Uniform(lb_Lr_D501,ub_Lr_D501))
def set_Lr_D501(Lr):
    F_baseline.D501.Lr = Lr
    
##Setting the heavier key recovery     
lb_Hr_D501 = 0.9
ub_Hr_D501 = 0.99
@model.parameter(element=F_baseline.D501,
                 kind='coupled',
                 distribution=shape.Uniform(lb_Hr_D501,ub_Hr_D501))
def set_Hr_D501(Hr):
    F_baseline.D501.Hr = Hr
    
#TODO: F_baseline.unit.D501.P maybe
#6 separation of azelaic acid rich fraction and diols
#Setting the lighter key recovery
lb_Lr_D604 = 0.9
ub_Lr_D604 = 0.99
@model.parameter(element=F_baseline.D604,
                 kind='coupled',
                 distribution=shape.Uniform(lb_Lr_D604,ub_Lr_D604))
def set_Lr_D604(Lr):
    F_baseline.D604.Lr = Lr
    
##Setting the heavier key recovery     
lb_Hr_D604 = 0.9
ub_Hr_D604 = 0.99
@model.parameter(element=F_baseline.D604,
                 kind='coupled',
                 distribution=shape.Uniform(lb_Hr_D604,ub_Hr_D604))
def set_Hr_D604(Hr):
    F_baseline.D604.Hr = Hr
    
#7 separation of ligher boiling impurities from the azelaic acid product stream
#Setting the lighter key recovery
lb_Lr_D605 = 0.9
ub_Lr_D605 = 0.99
@model.parameter(element=F_baseline.D605,
                 kind='coupled',
                 distribution=shape.Uniform(lb_Lr_D605,ub_Lr_D605))
def set_Lr_D605(Lr):
    F_baseline.D605.Lr = Lr
    
##Setting the heavier key recovery     
lb_Hr_D605 = 0.9
ub_Hr_D605 = 0.99
@model.parameter(element=F_baseline.D605,
                 kind='coupled',
                 distribution=shape.Uniform(lb_Hr_D605,ub_Hr_D605))
def set_Hr_D605(Hr):
    F_baseline.D501.Hr = Hr

############################################################################
#Isolated parameters
crude_oil_feedstock = F_baseline.stream.crude_vegetable_oil# The feedstock stream
lb_o = crude_oil_feedstock.price * 0.9 # Minimum price #TODO: how else can i change this
ub_o = crude_oil_feedstock.price * 1.1 # Maximum price
@model.parameter(element=crude_oil_feedstock, kind='isolated', units='USD/kg',
                 distribution=shape.Uniform(lb_o, ub_o))
def set_feedstock_price(feedstock_price):
    crude_oil_feedstock.price = feedstock_price
    
fresh_tungsten_catalyst = F_baseline.stream.fresh_tungsten_catalyst
lb_t = fresh_tungsten_catalyst.price*0.9
ub_t = fresh_tungsten_catalyst.price*1.1
@model.parameter(element=fresh_tungsten_catalyst, kind='isolated', units='USD/kg',
                 distribution=shape.Uniform(lb_t, ub_t))
def set_tungstencat_price(tungstencat_price):
   fresh_tungsten_catalyst.price = tungstencat_price

fresh_cobalt_catalyst = F_baseline.stream.fresh_cobalt_catalyst_stream
lb_c = fresh_cobalt_catalyst.price*0.9
ub_c = fresh_cobalt_catalyst.price*1.1
@model.parameter(element=fresh_cobalt_catalyst, kind='isolated', units='USD/kg',
                 distribution=shape.Uniform(lb_c, ub_c))
def set_tungstencat_price(fresh_cobalt_catalyst_price):
   fresh_cobalt_catalyst.price = fresh_cobalt_catalyst_price   

model.show()

