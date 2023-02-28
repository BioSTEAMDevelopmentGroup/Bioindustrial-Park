# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:45:33 2023
@author: Lavanya
"""
import biosteam as bst
import numpy as np
import chaospy
from chaospy import distributions as shape
from biorefineries.oleochemicals.systems_baseline import F_baseline
from biorefineries.oleochemicals.systems_baseline import aa_baseline_sys
import numpy as np
from lca_tea_baseline import TEA_baseline
from biosteam.evaluation import Model, Metric
from biorefineries.lipidcane._process_settings import price #TODO: were these prices adjusted to 2013 prices?
from biorefineries.cane.data.lca_characterization_factors import GWP_characterization_factors 
from biorefineries.cornstover import CellulosicEthanolTEA

#Settings to set GWP100 as the main characterisation factor
GWP = 'GWP100'
bst.settings.define_impact_indicator(key=GWP, units='kg*CO2e')
# NG-Fired Simple-Cycle Gas Turbine CHP Plant, no transmission included
bst.settings.set_electricity_CF(GWP, 0.36, basis='kWhr', units='kg*CO2e')
aa_baseline_sys = aa_baseline_sys()
# aa_baseline_sys.prioritize_unit(F_baseline.unit.S702)
aa_baseline_sys.set_tolerance(mol=0.0,
                              rmol=0.001,
                              method = 'fixedpoint',
                              subsystems = True)
aa_baseline_sys.simulate()  
#############################################################################################################
# renaming the first system factory for biodiesel prep as the 1000 series
biodiesel_prep_units = (F_baseline.unit.S402,F_baseline.unit.T401,
                        F_baseline.unit.P401,F_baseline.unit.T402,
                        F_baseline.unit.P402,F_baseline.unit.T403,
                        F_baseline.unit.P403,F_baseline.unit.T404,
                        F_baseline.unit.P404,F_baseline.unit.S401,
                        F_baseline.unit.R401,F_baseline.unit.C401,
                        F_baseline.unit.P405,F_baseline.unit.R402,
                        F_baseline.unit.C402,F_baseline.unit.T405,
                        F_baseline.unit.P406,F_baseline.unit.C403,
                        F_baseline.unit.F401,F_baseline.unit.P407,
                        F_baseline.unit.H401,F_baseline.unit.P408,
                        F_baseline.unit.T406,F_baseline.unit.P409,
                        F_baseline.unit.C404,F_baseline.unit.T407,
                        F_baseline.unit.P410,F_baseline.unit.D401,
                        F_baseline.unit.H402,F_baseline.unit.D402,
                        F_baseline.unit.P413,F_baseline.unit.H403,
                        F_baseline.unit.P411,F_baseline.unit.H404,
                        F_baseline.unit.P412,F_baseline.unit.T408,
                        F_baseline.unit.T409,F_baseline.unit.B401)
bst.rename_units(units = biodiesel_prep_units, area = 1000)
aa_baseline_groups = bst.UnitGroup.group_by_area(aa_baseline_sys.units)
groups_by_area = bst.UnitGroup.df_from_groups(aa_baseline_groups, fraction = True) #populates the dataframe as percentages
plot_by_area = groups_by_area.plot.bar(stacked = True) #plot for all the areas stacked by metrics
plot_by_metric =  groups_by_area.T
plot_by_metric.plot.bar(stacked = True) #plot for all the metrics stacked by area
#TODO: decide how to name areas
#What to create boxplots for?
# boxplot = plot_by_metric.boxplot(column=['Heating duty']) #for a specific column name, .boxplot() for the entire dataframe


# #########################################################################################################
#TODO: add price later insiide the system
F_baseline.pelargonic_acid_rich_fraction.price = 30
# # Streams specs belonging to the cane biorefinery used for biodisel prep
# #Methanol
F_baseline.stream.methanol.price = 0.792*401.693/275.700 #Based on Catbio costs adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals
F_baseline.stream.methanol.characterization_factors = {'GWP100': GWP_characterization_factors['methanol']}

#Catalyst
F_baseline.stream.catalyst.price = 0.25*(price['NaOCH3']*401.693/259.900) + 0.75*( 0.792*401.693/275.700) #Adjusted from 2019 to 2022, Fred's PPI for industrial chemicals
F_baseline.stream.catalyst.characterization_factors = {'GWP100': GWP_characterization_factors['methanol catalyst mixture']}

#Biodiesel wash water
F_baseline.stream.biodiesel_wash_water.price = 3.945 *(217.9/132.9)/(3.78541*1000)#Ref: DOE Annual water rates pdf,adjusted using FRED's PPI> Industry based> Utilities.(1kgal = 1000gal, 1gal = 3.78541 Kg)
F_baseline.stream.biodiesel_wash_water.characterization_factors={'GWP100': 0.00035559}#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)

#HCl 
F_baseline.stream.HCl.price = 0.88*401.693/275.700 #Based on Catbio price ($/Kg) for Hydrogen peroxide, 35%, tech., tankcars, works, frt. equald. adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals 
F_baseline.stream.HCl.characterization_factors = {'GWP100': GWP_characterization_factors['HCl']}

#NaOH
F_baseline.stream.NaOH.price = 0.93*401.693/275.700 #Based on Catbio price ($/Kg) for Caustic soda (sodium hydroxide), liq., dst spot barge f.o.b. USG adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals 
F_baseline.stream.NaOH.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}


#ask Yoel if this should be considered or the prices from economic assessment paper should be taken instead
#crude_glycerol
F_baseline.stream.crude_glycerol.price = price['Crude glycerol']*401.693/275.700 #Adjusted from 2021 Jan to 2022 Dec based on FRED's PPI
F_baseline.stream.crude_glycerol.characterization_factors = {'GWP100': GWP_characterization_factors['crude-glycerol']}
 
# #######################################################################################################################33
# #####################################################################################################33
tea_azelaic_baseline = TEA_baseline(
                                    system = aa_baseline_sys,
                                    lang_factor = None,                                    
                                    IRR =  0.1,#TODO:?
                                    duration=(2023,2053), #TODO:?
                                    depreciation = 'MACRS7',#Ref: Cellulosic ethanol 
                                    income_tax = 0.35,#TODO:?
                                    operating_days = 300,#TODO:?                                    
                                    construction_schedule=(0.08, 0.60, 0.32),#Ref: Cellulosic ethanol 
                                    startup_months=3, #Ref: Cellulosic ethanol 
                                    startup_FOCfrac=1,#Ref: Cellulosic ethanol 
                                    startup_salesfrac=0.5,#Ref: Cellulosic ethanol 
                                    startup_VOCfrac=0.75,#Ref: Cellulosic ethanol 
                                    WC_over_FCI=0.05,#Ref: Cellulosic ethanol 
                                    finance_interest=0.08,#Ref: Cellulosic ethanol 
                                    finance_years=10,#Ref: Cellulosic ethanol 
                                    finance_fraction=0.4,#Ref: Cellulosic ethanol                 
                                    operating_labor_cost = 54874*5*5*300,#US beaurue of stats and Ref: economic risk asessment paper 
                                    #A total of 5 operators/shift and 5 shifts per day are needed
                                    direct_supervisory_clerical_labor = 0.18,
                                    maintenance_and_repairs = 0.06,#Ref(Turton et al., 2013)
                                    operating_supplies = 0.009,    #Ref(Turton et al., 2013)      
                                    laboratory_charges = 0.15,#Ref(Turton et al., 2013)
                                    local_taxes_and_insurance = 0.1, #Ref(Turton et al., 2013)
                                    plant_overhead = 0.708,#Ref(Turton et al., 2013)
                                    administration_costs = 0.177,#Ref(Turton et al., 2013)
                                    OSBL_units = [  F_baseline.CW901,
                                                    F_baseline.CT901,
                                                    F_baseline.BT901,
                                                    F_baseline.PWT901, 
                                                    F_baseline.W901,
                                                    F_baseline.ADP901,
                                                  ])                                                         
                                   
aa_sys_op_hours = aa_baseline_sys.operating_hours = tea_azelaic_baseline.operating_days * 24

#functions to solve for different indicators
# Minimum product selling price of lactic_acid stream
azelaic_acid = F_baseline.stream.azelaic_acid_product_stream
# funcs = create_funcs(lactic_tea=lactic_sys.TEA, flowsheet=flowsheet)
azelaic_acid_tea = aa_baseline_sys.TEA
def get_MPSP():
        azelaic_acid.price = 0
        #TODO: ASK WHY? for i in range(3):
        MPSP = azelaic_acid.price = azelaic_acid_tea.solve_price(azelaic_acid)
        return MPSP

crude_HOSO = F_baseline.stream.crude_vegetable_oil
# Yield in 10^6 kg/yr
get_annual_factor = lambda: azelaic_acid_tea.operating_days*24
get_total_yield = lambda: azelaic_acid.F_mass/1e6
# Purity (%) of LacticAcid in the final product
get_purity = lambda: azelaic_acid.imass['Azelaic_acid']/azelaic_acid.F_mass
get_overall_TCI = lambda: azelaic_acid_tea.TCI/1e6
# get_azelaic_acid_sale = lambda: get_total_yield*azelaic_acid.price
get_material_cost = lambda: azelaic_acid_tea.material_cost/1e6 #TODO: why should ash dsposal be subtracted?
get_operating_cost = lambda: azelaic_acid_tea.AOC/1e6 #TODO: why to sub? -get_gypsum_sale()-get_ash_sale()
##### Utilities related #####
get_system_heating_demand = lambda: aa_baseline_sys.get_heating_duty()/1e9
get_system_cooling_water_duty = lambda: aa_baseline_sys.get_cooling_duty()/1e9
# To see if TEA converges well for each simulation #TODO: ask why?
get_NPV = lambda: azelaic_acid_tea.NPV
#Environmental impact
# get_GWP = lambda: get_process_impact(key)

metrics = [    
            Metric('MPSP', get_MPSP, '$/kg'),
            Metric('Total product yield', get_total_yield, '10^6 kg/yr'),
            Metric('Product purity', get_purity, '%'),
            Metric('Total capital investment', get_overall_TCI, '10^6 $'),
            Metric('Annual operating cost', get_operating_cost, '10^6 $/yr'),
            Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
            # Metric('Annual product sale', get_azelaic_acid_sale, '10^6 $/yr'),
            Metric('Total', get_system_heating_demand, '10^6 MJ/yr', 'Heating demand'),
            Metric('Total', get_system_cooling_water_duty, '10^6 MJ/yr', 'Cooling demand'),
            Metric('NPV', get_NPV, '$', 'TEA'),
            # Metric('Total GWP', get_GWP, 'kg CO2-eq/kg', 'LCA'),
          ]
model = Model(aa_baseline_sys, metrics, exception_hook='raise')
#TODO: ask what check functions are in lactic acid models
#################################################################################################
# #Isolated parameters
crude_oil_feedstock = F_baseline.stream.crude_vegetable_oil# The feedstock stream
lb_o = crude_oil_feedstock.price * 0.9 # Minimum price #TODO: how else can i change this
ub_o = crude_oil_feedstock.price * 1.1 # Maximum price
@model.parameter(name = 'HOSO price',
                 element=crude_oil_feedstock, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_o, ub_o))
def set_feedstock_price(feedstock_price):
    crude_oil_feedstock.price = feedstock_price
    
fresh_tungsten_catalyst = F_baseline.stream.fresh_tungsten_catalyst
lb_t = fresh_tungsten_catalyst.price*0.9
ub_t = fresh_tungsten_catalyst.price*1.1
@model.parameter(name = 'Tungsten catalyst price',
                 element=fresh_tungsten_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_t, ub_t))
def set_tungstencat_price(tungstencat_price):
    fresh_tungsten_catalyst.price = tungstencat_price

fresh_cobalt_catalyst = F_baseline.stream.fresh_cobalt_catalyst_stream
lb_c = fresh_cobalt_catalyst.price*0.9
ub_c = fresh_cobalt_catalyst.price*1.1
@model.parameter(name = 'Cobalt acetate catalyst price',
                 element=fresh_cobalt_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_c, ub_c))
def set_cobaltcat_price(fresh_cobalt_catalyst_price):
    fresh_cobalt_catalyst.price = fresh_cobalt_catalyst_price   

fresh_solvent = F_baseline.stream.solvent_for_extraction
fresh_solvent.price = 566*2.8535/55
lb_s = fresh_solvent.price*0.9
ub_s = fresh_solvent.price*1.1
@model.parameter(name = 'VM&P Solvent price',
                 element=fresh_solvent, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_s, ub_s))
def set_solvent_price(fresh_solvent_price):
        fresh_solvent.price = fresh_solvent_price  

fresh_HP = F_baseline.stream.fresh_HP
fresh_HP.price = 1.21*152.650/77.500
lb_hp = fresh_HP.price*0.9
ub_hp = fresh_HP.price*1.1
@model.parameter(name = 'Hydrogen peroxide price',
                 element=fresh_HP, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_hp, ub_hp))
def set_fresh_HP_price(fresh_HP_price):
        fresh_HP.price = fresh_HP_price 

F_baseline.BT901.natural_gas_price = 0.253 #REF: lactic acid SI
@model.parameter(name = 'Natural gas price',
                 element = F_baseline.BT901,
                 kind = 'isolated',
                 units='USD/kg',
                 distribution=shape.Triangle(0.198, 0.253, 0.304) #REF: lactic acid SI
                  )
def set_natural_gas_price(natural_gas_price):
        F_baseline.BT901.natural_gas_price = natural_gas_price
                 

#############################################################################################################33
#coupled parameters
#Process related parameters
#TODO: ask YOEL how to access this?s
# # 1 Transesterification reaction conversion
# lb1 = 0.90
# ub1 = 0.99
# @model.parameter(element=F_baseline.U1001,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb1,ub1))
# def set_transesterification_conversion_1(X_tes):
#     for i in range (0,11):
#         F_baseline.U1001.X_tes = X_tes
        
# @model.parameter(element=F_baseline.U1002,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb1,ub1))
# def set_transesterification_conversion_2(X_tes):
#     for i in range (0,11):
#         F_baseline.U1002.transesterification.X[i] = X_tes         

#2 Dihydroxylation reaction conversion
lb3 = 0.90
ub3 = 0.99
@model.parameter( name = 'Dihydroxylation reaction conversion',
                  element=F_baseline.R101,
                  kind='coupled',
                  distribution=shape.Uniform(lb3,ub3))
def set_dihydroxylation_conversion(X_dih):
        F_baseline.unit.R101.X_dih = X_dih
  
#3 Oxidative cleavage reaction conversion
# F_baseline.OxidativeCleavageReactor.X_oxidativecleavage = 0.8
lb4 = 0.90
ub4 = 0.99
@model.parameter( name = 'Oxidative cleavage reaction conversion', 
                  element=F_baseline.R202,
                  kind='coupled',
                  distribution=shape.Uniform(lb4,ub4))
def set_X_oxidativecleavage_conversion(X_oxidativecleavage):
        F_baseline.unit.R202.X_oxidativecleavage = X_oxidativecleavage        


#3 Oxidative cleavage reaction conversion
# F_baseline.OxidativeCleavageReactor.X_decarboxylation = 0.2
lb5 = 0.90
ub5 = 0.99
@model.parameter( name = 'Decarboxylation reaction conversion' , 
                  element=F_baseline.R202,
                  kind='coupled',
                  distribution=shape.Uniform(lb5,ub5))
def set_X_decarboxylation_conversion(X_decarboxylation):
        F_baseline.unit.R202.X_decarboxylation =  X_decarboxylation

#3 Oxidative cleavage reaction conversion
# F_baseline.OxidativeCleavageReactor.X_side_rxn = 0.8        
lb6 = 0.90
ub6 = 0.99
@model.parameter( name = 'Side reaction conversion' ,
                  element=F_baseline.R202,
                  kind='coupled',
                  distribution=shape.Uniform(lb6,ub6))
def set_X_side_rxn_conversion(X_side_rxn):
        F_baseline.unit.R202.X_side_rxn =  X_side_rxn

lb7 = 10000 
ub7 = 1000000
@model.parameter(name = 'Feedstock_capacity',
                 element = F_baseline.stream.crude_vegetable_oil,
                 kind = 'coupled',
                 distribution=shape.Uniform(lb7,ub7)
                 )
def set_feedstock_input_flow(Cap):
    F_baseline.stream.crude_vegetable_oil.F_mass = Cap    
    
N_samples = 1000
rule = 'L' # For Latin-Hypercube sampling
np.random.seed(1234) # For consistent results
samples = model.sample(N_samples, rule)
model.load_samples(samples)
model.evaluate(
    notify=100 # Also print elapsed time after 50 simulations
              )
model.show()

 ############################################################################3
# #4 Pressure for degassing_the_oily_phase
# lb4 = 10000
# ub4 = 20000
# @model.parameter(element=F_baseline.F2001,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb4,ub4))
# def set_degassingcolumn_P(P_d):
#     F_baseline.unit.F2001.P = P_d
# #5 nonanoic_acid_fraction_separation
# #Setting the lighter key recovery
# lb_Lr_D501 = 0.9
# ub_Lr_D501 = 0.99
# @model.parameter(element=F_baseline.D501,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Lr_D501,ub_Lr_D501))
# def set_Lr_D501(Lr):
#     F_baseline.D501.Lr = Lr
    
# ##Setting the heavier key recovery     
# lb_Hr_D501 = 0.9
# ub_Hr_D501 = 0.99
# @model.parameter(element=F_baseline.D501,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Hr_D501,ub_Hr_D501))
# def set_Hr_D501(Hr):
#     F_baseline.D501.Hr = Hr
    
# #TODO: F_baseline.unit.D501.P maybe
# #6 separation of azelaic acid rich fraction and diols
# #Setting the lighter key recovery
# lb_Lr_D604 = 0.9
# ub_Lr_D604 = 0.99
# @model.parameter(element=F_baseline.D604,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Lr_D604,ub_Lr_D604))
# def set_Lr_D604(Lr):
#     F_baseline.D604.Lr = Lr
    
# ##Setting the heavier key recovery     
# lb_Hr_D604 = 0.9
# ub_Hr_D604 = 0.99
# @model.parameter(element=F_baseline.D604,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Hr_D604,ub_Hr_D604))
# def set_Hr_D604(Hr):
#     F_baseline.D604.Hr = Hr
    
# #7 separation of ligher boiling impurities from the azelaic acid product stream
# #Setting the lighter key recovery
# lb_Lr_D605 = 0.9
# ub_Lr_D605 = 0.99
# @model.parameter(element=F_baseline.D605,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Lr_D605,
#                                              ub_Lr_D605))
# def set_Lr_D605(Lr):
#     F_baseline.D605.Lr = Lr
    
# ##Setting the heavier key recovery     
# lb_Hr_D605 = 0.9
# ub_Hr_D605 = 0.99
# @model.parameter(element=F_baseline.D605,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Hr_D605,
#                                              ub_Hr_D605))
# def set_Hr_D605(Hr):
#     F_baseline.D501.Hr = Hr
     
# ############################################################################