# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:45:33 2023
@author: Lavanya
"""
import biosteam as bst
import numpy as np
import chaospy
from chaospy import distributions as shape
from biorefineries.oleochemicals.systems_baseline_hoysoy import F_baseline
from biorefineries.oleochemicals.systems_baseline_hoysoy import aa_baseline_sys
import numpy as np
from lca_tea_baseline_hoysoy import TEA_baseline
from biosteam.evaluation import Model, Metric
from biorefineries.lipidcane._process_settings import price #TODO: were these prices adjusted to 2013 prices?
from biorefineries.cane.data.lca_characterization_factors import GWP_characterization_factors 
from biorefineries.tea.cellulosic_ethanol_tea import CellulosicEthanolTEA,create_cellulosic_ethanol_tea
from thermosteam.utils import GG_colors,GG_light_colors
from units_baseline import HydrolysisReactor

#Settings to set GWP100 as the main characterisation factor
GWP = 'GWP100'
bst.settings.CEPCI = 708 #for 2021, based on CEPCI_by_year, biosteam #TODO: change
bst.settings.define_impact_indicator(key=GWP, units='kg*CO2e')
# NG-Fired Simple-Cycle Gas Turbine CHP Plant, no transmission included
bst.settings.set_electricity_CF(GWP, 0.36, basis='kWhr', units='kg*CO2e')
bst.PowerUtility.price = 0.070 #Lactic acid
aa_baseline_sys = aa_baseline_sys()
aa_baseline_sys.set_tolerance(mol=0.0,
                              rmol=0.3,
                              subsystems = True)
aa_baseline_sys.simulate()
    
#############################################################################################################
# renaming the first system factory for biodiesel prep as the 1000 series
biodiesel_prep_units = (F_baseline.unit.S402,F_baseline.unit.T401,F_baseline.unit.P401,F_baseline.unit.T402,F_baseline.unit.P402,F_baseline.unit.T403,F_baseline.unit.P403,F_baseline.unit.T404,
                        F_baseline.unit.P404,F_baseline.unit.S401,F_baseline.unit.R401,F_baseline.unit.C401,F_baseline.unit.P405,F_baseline.unit.R402,F_baseline.unit.C402,F_baseline.unit.T405,
                        F_baseline.unit.P406,F_baseline.unit.C403,F_baseline.unit.F401,F_baseline.unit.P407,F_baseline.unit.H401,F_baseline.unit.P408,F_baseline.unit.T406,F_baseline.unit.P409,
                        F_baseline.unit.C404,F_baseline.unit.T407,F_baseline.unit.P410,F_baseline.unit.D401,F_baseline.unit.H402,F_baseline.unit.D402,F_baseline.unit.P413,F_baseline.unit.H403,
                        F_baseline.unit.P411,F_baseline.unit.H404,F_baseline.unit.P412,F_baseline.unit.T408,F_baseline.unit.T409,F_baseline.unit.B401)
bst.rename_units(units = biodiesel_prep_units, area = 1000)

aa_baseline_groups = bst.UnitGroup.group_by_area(aa_baseline_sys.units)
groups_by_area = bst.UnitGroup.df_from_groups(aa_baseline_groups, fraction = True) #populates the dataframe as percentages
plot_by_area = groups_by_area.plot.bar(stacked = True) #plot for all the areas stacked by metrics
plot_by_metric =  groups_by_area.T
#plot for all the metrics stacked by area
groups_by_area.T.plot.bar(stacked = True,color =  [ GG_colors.blue.RGBn,GG_light_colors.yellow.RGBn,GG_colors.green.RGBn,GG_colors.red.RGBn,GG_colors.purple.RGBn,GG_light_colors.orange.RGBn,GG_light_colors.green.RGBn,]).set_yticks((10,20,30,40,50,60,70,80,90,100),minor=True)
# legend(loc = 'upper right',bbox_to_anchor=(1.25,0.75))
# legend(['Dihydroxylation','Oxidative cleavage','Catalyst Recovery','Pelargonic acid recovery and C5_C8_fraction  recovery','Azelaic acid recovery','CT,PWT,BT and CWP','Biodiesel production'],loc = 'upper right',bbox_to_anchor=(1.958,0.8))
# plot.bar(stacked = True,color =  [ GG_colors.blue.RGBn,GG_light_colors.yellow.RGBn,GG_colors.green.RGBn,GG_colors.red.RGBn,GG_colors.purple.RGBn,GG_light_colors.orange.RGBn,GG_light_colors.green.RGBn],yticks = (10,20,30,40,50,60,70,80,90,100)).legend(['Dihydroxylation','Oxidative cleavage','Catalyst Recovery','Pelargonic acid recovery and C5_C8_fraction  recovery','Azelaic acid recovery','CT,PWT,BT and CWP','Biodiesel production'],loc = 'upper right',bbox_to_anchor=(1.958,0.8))
#Code to set ticks
# boxplot = plot_by_metric.boxplot(column=['Heating duty']) #for a specific column name, .boxplot() for the entire dataframe
# plot_by_area = groups_by_area.plot.bar(stacked = True,color = GG_colors.blue.RGBn)
# #########################################################################################################
# # Streams specs belonging to the cane biorefinery used for biodisel prep
# #Methanol
F_baseline.stream.methanol.price = 0.792*401.693/275.700 #Based on Catbio costs adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals
F_baseline.stream.methanol.characterization_factors = {'GWP100': GWP_characterization_factors['methanol']}

#Catalyst (25% methoxide and 75% methanol)
F_baseline.stream.catalyst.price = 0.25*(price['NaOCH3']*401.693/259.900) + 0.75*( 0.792*401.693/275.700) #Adjusted from 2019 to 2022, Fred's PPI for industrial chemicals
F_baseline.stream.catalyst.characterization_factors = {'GWP100': GWP_characterization_factors['methanol catalyst mixture']}

#Biodiesel wash water
F_baseline.stream.biodiesel_wash_water.price = 3.945/(3.78541*1000)#Ref: DOE Annual water rates pdf,adjusted using FRED's PPI> Industry based> Utilities.(1kgal = 1000gal, 1gal = 3.78541 Kg)
# F_baseline.stream.biodiesel_wash_water.characterization_factors={'GWP100': 0.00035559}#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)

#HCl 
F_baseline.stream.HCl.price = 0.14*401.693/275.700 #Based on Catbio price ($/Kg) for 22 deg baume, US gulf dom. adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals 
F_baseline.stream.HCl.characterization_factors = {'GWP100': GWP_characterization_factors['HCl']*0.35 + 0.00035559*0.65}

#NaOH
F_baseline.stream.NaOH.price = 0.93*401.693/275.700 #Based on Catbio price ($/Kg) for Caustic soda (sodium hydroxide), liq., dst spot barge f.o.b. USG adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals 
F_baseline.stream.NaOH.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}


#ask Yoel if this should be considered or the prices from economic assessment paper should be taken instead
#crude_glycerol
F_baseline.stream.crude_glycerol.characterization_factors = {'GWP100': GWP_characterization_factors['crude-glycerol']}
 
#Crude methanol
#TODO:??? how is the price same methanol and crude methanol
F_baseline.stream.crude_methanol.characterization_factors = {'GWP100': GWP_characterization_factors['methanol']}

#https://www.teknorapex.com/why-esters-di-esters-polyol-esters-trimelliate-esters#:~:text=Synthetic%20esters%2C%20with%20their%20polarity%20and%20uniform%20structure%2C,metal%20surface%2C%20and%20improve%20cleanliness%20and%20sludge%20control.
F_baseline.stream.recovered_C5_to_C8_MCA_fraction.price = 4.47*401.693/275.7 #Based on Adipic acid resin, resin grade bulk, hopper cars, frt. equald. #This is based on the fact that adipates are used in the making of lubricants 
F_baseline.stream.pelargonic_acid_rich_fraction.price = 5.56*357.592/255.400 #thesis check #EXW price of glyphosate 95% AI was around 5.56 USD/Kg
F_baseline.stream.fatty_acid_blend.price = 1.23 #Stearic acid, rubber grade, bulk
#TODO: change above stream to oleic acid
F_baseline.stream.azelaic_acid_product_stream.price = 6.88*401.693/275.7 #Based on sebacic acid bulk price
F_baseline.stream.crude_glycerol.price = price['Crude glycerol']*401.693/259.9 #check thesis
F_baseline.crude_methanol.price =  0.792*401.693/275.700 #Based on Catbio costs adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals

azelaic_acid = F_baseline.stream.azelaic_acid_product_stream
recovered_C5_to_C8_MCA_fraction = F_baseline.stream.recovered_C5_to_C8_MCA_fraction
pelargonic_acid_rich_fraction = F_baseline.stream.pelargonic_acid_rich_fraction
fatty_acid_blend = F_baseline.stream.fatty_acid_blend

#######################################################################################################################
#####################################################################################################

tea_azelaic_baseline = TEA_baseline(
                                    system = aa_baseline_sys,
                                    lang_factor = None,                                    
                                    IRR =  0.10,#Ref: economic risk asessment paper 
                                    duration=(2022,2032), #TODO:?
                                    depreciation = 'MACRS7',#Ref: Cellulosic ethanol 
                                    income_tax = 0.35,#Ref: Cellulosic ethanol
                                    operating_days = 300,#TODO:?                                    
                                    construction_schedule=(2/3, 1/3),#Ref: economic risk asessment paper  
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
                                                    F_baseline.ADP901,
                                                  ])         
                                
aa_sys_op_hours = aa_baseline_sys.operating_hours = tea_azelaic_baseline.operating_days * 24

#functions to solve for different indicators
azelaic_acid_tea = aa_baseline_sys.TEA
def get_MPSP():
        azelaic_acid.price = 0
        MPSP = azelaic_acid.price = azelaic_acid_tea.solve_price(azelaic_acid)
        return MPSP

crude_HOSO = F_baseline.stream.crude_vegetable_oil
# Yield in 10^6 kg/yr
get_annual_factor = lambda: azelaic_acid_tea.operating_days*24
get_total_yield = lambda: azelaic_acid.F_mass/1e6
# Purity (%) of LacticAcid in the final product
get_purity = lambda: azelaic_acid.imass['Azelaic_acid']/azelaic_acid.F_mass
get_overall_TCI = lambda: azelaic_acid_tea.TCI/1e6
# get_azelaic_acid_sale = lambda: get_total_yield*azelaic_acid.price/1e6
get_material_cost = lambda: azelaic_acid_tea.material_cost/1e6 #TODO: why should ash dsposal be subtracted?
get_operating_cost = lambda: azelaic_acid_tea.AOC/1e6 #TODO: why to sub? -get_gypsum_sale()-get_ash_sale()
##### Utilities related #####
get_system_heating_demand = lambda: aa_baseline_sys.get_heating_duty()/1e9
get_system_cooling_water_duty = lambda: aa_baseline_sys.get_cooling_duty()/1e9
# To see if TEA converges well for each simulation #TODO: ask why?
get_NPV = lambda: azelaic_acid_tea.NPV
solve_IRR = lambda: azelaic_acid_tea.solve_IRR()
#Environmental impact
#TODO: net impact works but get_process_impact does not
# get_net_GWP = lambda: aa_baseline_sys.get_net_impact(key = GWP)
#Mass allocation method: TODO: check
# GWP_per_Kg_of_crude_oil = aa_baseline_sys.get_net_impact(key = GWP)/aa_baseline_sys.get_mass_flow(stream = F_baseline.crude_vegetable_oil)
# GWP_ofAA_per_Kg_of_oil = GWP_per_Kg_of_crude_oil*aa_baseline_sys.get_mass_flow(stream = azelaic_acid)/aa_baseline_sys.get_mass_flow(stream = F_baseline.crude_vegetable_oil)
#Revenue allocation
# GWP_per_USD = aa_baseline_sys.get_property_allocated_impact(
    # key=GWP, name='revenue', units='USD/hr',
# )
#  GWP_AA_revenue = (
    # GWP_per_USD * azelaic_acid.price
# ) # kg

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
            Metric('Internal rate of return', solve_IRR, '%'),
            # Metric('Total Net GWP', get_net_GWP, 'kg CO2-eq', 'LCA'),
          ]
model = Model(aa_baseline_sys,
              metrics,
              )
#TODO: ask what check functions are in lactic acid models
#################################################################################################
# #Isolated parameters
crude_oil_feedstock = F_baseline.stream.crude_vegetable_oil# The feedstock stream
lb_o = crude_oil_feedstock.price * 0.75 # Minimum price #TODO: how else can i change this
ub_o = crude_oil_feedstock.price * 1.25 # Maximum price
@model.parameter(name = 'HoSun oil price',
                 element=crude_oil_feedstock, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_o, ub_o))
def set_feedstock_price(feedstock_price):
    crude_oil_feedstock.price = feedstock_price
    
fresh_tungsten_catalyst = F_baseline.stream.fresh_tungsten_catalyst
lb_t = fresh_tungsten_catalyst.price*0.75
ub_t = fresh_tungsten_catalyst.price*1.25
@model.parameter(name = 'Tungstic acid catalyst price',
                 element=fresh_tungsten_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_t, ub_t))
def set_tungstencat_price(tungstencat_price):
    fresh_tungsten_catalyst.price = tungstencat_price

fresh_cobalt_catalyst = F_baseline.stream.fresh_cobalt_catalyst_stream
lb_c = fresh_cobalt_catalyst.price*0.75
ub_c = fresh_cobalt_catalyst.price*1.25
@model.parameter(name = 'Cobalt acetate catalyst price',
                 element=fresh_cobalt_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_c, ub_c))
def set_cobaltcat_price(fresh_cobalt_catalyst_price):
    fresh_cobalt_catalyst.price = fresh_cobalt_catalyst_price   

fresh_solvent = F_baseline.stream.solvent_for_extraction
fresh_solvent.price = 0.996
lb_s = fresh_solvent.price*0.75
ub_s = fresh_solvent.price*1.25
@model.parameter(name = 'Heptane solvent price',
                 element=fresh_solvent, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_s, ub_s))
def set_solvent_price(fresh_solvent_price):
        fresh_solvent.price = fresh_solvent_price  

fresh_HP = F_baseline.stream.fresh_HP
lb_hp = fresh_HP.price*0.75
ub_hp = fresh_HP.price*1.25
@model.parameter(name = 'Hydrogen peroxide price',
                 element=fresh_HP, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_hp, ub_hp))
def set_fresh_HP_price(fresh_HP_price):
        fresh_HP.price = fresh_HP_price 
#TODO: ask the units in cane.data.price for mean_natural??

F_baseline.BT901.natural_gas_price = 0.253 #REF: lactic acid SI
@model.parameter(name = 'Natural gas price',
                 element = F_baseline.BT901,
                 kind = 'isolated',
                 units='USD/kg',
                 distribution=shape.Triangle(0.198,
                                             0.253, 
                                             0.304) #REF: lactic acid SI
                  )
def set_natural_gas_price(natural_gas_price):
        F_baseline.BT901.natural_gas_price = natural_gas_price

citric_acid = F_baseline.stream.citricacid_for_degumming
lb_ca = citric_acid.price*0.75
ub_ca = citric_acid.price*1.25
@model.parameter(name = 'Citric acid price',
                 element = citric_acid,
                 kind = 'isolated',
                 units='USD/kg',
                 distribution=shape.Uniform(lb_ca,ub_ca)
                  )
def set_citric_acid_price(citric_acid_price):
       citric_acid.price  = citric_acid_price

conc_hydrochloric_acid = F_baseline.stream.conc_hydrochloric_acid
lb_hcl = conc_hydrochloric_acid.price*0.75
ub_hcl = conc_hydrochloric_acid.price*1.25
@model.parameter(name = 'Conc HCl price',
                 element = conc_hydrochloric_acid,
                 kind = 'isolated',
                 units='USD/kg',
                 distribution=shape.Uniform(lb_hcl,ub_hcl)
                  )
def set_HCl_price(conc_hydrochloric_acid_price):
       conc_hydrochloric_acid.price  = conc_hydrochloric_acid_price


calcium_chloride = F_baseline.stream.calcium_chloride
calcium_chloride.price = 0.493
lb_cc = calcium_chloride.price*0.75
ub_cc = calcium_chloride.price*1.25      
@model.parameter(name = 'Calcium chloride price',
                 element = calcium_chloride,
                 kind = 'isolated',
                 units='USD/kg',
                 distribution=shape.Uniform(lb_cc,ub_cc)
                  )
def set_calcium_chloride_price(calcium_chloride_price):
       calcium_chloride.price  = calcium_chloride_price

polystyrene_based_catalyst = F_baseline.stream.polystyrene_based_catalyst
polystyrene_based_catalyst.price = 5.55*280.446/240.300       
lb_pbc =  1.84*280.446/240.300
ub_pbc = 9.199*280.446/240.300
@model.parameter(name = 'Hydrolysis resin price',
                 element = polystyrene_based_catalyst,
                 kind = 'isolated',
                 units='USD/kg',
                 distribution=shape.Uniform(lb_pbc,ub_pbc)
                  )
def set_polystyrene_based_catalyst_price(polystyrene_based_catalyst_price):
       polystyrene_based_catalyst.price  = polystyrene_based_catalyst_price       

     
#############################################################################################################33
#coupled parameters
#Process related parameters
#1 Amount ratios varied
#Tungstic acid mole ratio
lb1 = 0.006
ub1 = 0.15
@model.parameter(name = 'Tungstic acid moles',
                 element = F_baseline.M101,
                 kind = 'coupled',
                 distribution=shape.Uniform(lb1,ub1)
                 )
def set_tungstic_acid_moles(X_tam):
        F_baseline.unit.M101.specifications[0].args[0] = X_tam
        
#Cobalt acetate catalyst ratio
lb2 = 0.003
ub2 = 0.015
@model.parameter(name = 'Cobalt acetate moles',
                  element = F_baseline.R202,
                  kind = 'coupled',
                  distribution=shape.Uniform(lb2,ub2)
                  )
def set_cobalt_acetate_moles(X_cam):
        F_baseline.unit.R202.specifications[0].args[0] = X_cam  

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
lb4 = 0.80
ub4 = 0.99
@model.parameter( name = 'Oxidative cleavage reaction conversion', 
                  element=F_baseline.R202,
                  kind='coupled',
                  distribution=shape.Uniform(lb4,ub4))
def set_X_oxidativecleavage_conversion(X_oxidativecleavage):
        F_baseline.unit.R202.X_oxidativecleavage = X_oxidativecleavage        


#3 Oxidative cleavage reaction conversion
# F_baseline.OxidativeCleavageReactor.X_decarboxylation = 0.2
lb5 = 0.1
ub5 = 0.5
@model.parameter( name = 'Decarboxylation reaction conversion' , 
                  element=F_baseline.R202,
                  kind='coupled',
                  distribution=shape.Uniform(lb5,ub5))
def set_X_decarboxylation_conversion(X_decarboxylation):
        F_baseline.unit.R202.X_decarboxylation =  X_decarboxylation

#3 Oxidative cleavage reaction conversion
# F_baseline.OxidativeCleavageReactor.X_side_rxn = 0.8        
lb6 = 0.1
ub6 = 0.5
@model.parameter( name = 'Side reaction conversion' ,
                  element=F_baseline.R202,
                  kind='coupled',
                  distribution=shape.Uniform(lb6,ub6))
def set_X_side_rxn_conversion(X_side_rxn):
        F_baseline.unit.R202.X_side_rxn =  X_side_rxn

lb7 = 35000*0.95
ub7 = 35000*1.05
@model.parameter(name = 'Feedstock_capacity',
                  element = F_baseline.stream.crude_vegetable_oil,
                  kind = 'coupled',
                  distribution=shape.Uniform(lb7,ub7)
                  )
def set_feedstock_input_flow(Cap):
    F_baseline.stream.crude_vegetable_oil.F_mass = Cap  

lb8 = 1.052
ub8 = 1.315
@model.parameter(name = 'Moles of air',
                  element = F_baseline.R202,
                  kind = 'coupled',
                  distribution=shape.Uniform(lb8,ub8)
                  )
def set_air_moles(X_air):
    F_baseline.unit.R202.specifications[0].args[1] = X_air

lb9 = 0.3
ub9 = 2.5
# @model.parameter(name = 'Extraction solvent mass',
#                   element = F_baseline.M602,
#                   kind = 'coupled',
#                   distribution=shape.Uniform(lb9,ub9)
#                   )
# def set_solvent_mass(X_sm):
#     F_baseline.unit.M602.specifications[0].args[0] = X_sm
    
    
lb10 = 0.08
ub10 = 0.15
@model.parameter(name = 'System IRR',
                  element = azelaic_acid_tea.IRR,
                  kind = 'cost',
                  distribution=shape.Uniform(lb10,ub10)
                  )
def set_irr(X_irr):
    azelaic_acid_tea.IRR = X_irr
    

N_samples = 1000
rule = 'L' # For Latin-Hypercube sampling
np.random.seed(1234) # For consistent results
samples = model.sample(N_samples, rule)
model.load_samples(samples)
model.evaluate(
    notify=20 # Also print elapsed time after 100 simulations
              )
model.show()

