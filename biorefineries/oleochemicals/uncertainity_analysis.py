# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:58:21 2023
@author: lavan
"""

#TODO: how would we consider electricity use if BT901 purchases electricity from the grid

from biorefineries.oleochemicals.systems_baseline import F
from biorefineries.oleochemicals.chemicals_baseline import chems
from biorefineries.oleochemicals.system_simulate import azelaic_acid_tea,aa_baseline
from tea_baseline import TEA_baseline
from biosteam.evaluation import Model, Metric
from prices_and_GWP_factors import prices_per_stream,utility_prices,GWP_per_stream,Utility_GWP_factors
from warnings import filterwarnings; filterwarnings('ignore')
from chaospy import distributions as shape
import numpy as np

#unit conversions
KJpersec_to_KJhr = 3600 

#shortform stream names for ease
azelaic_acid = F.stream.azelaic_acid_product_stream
recovered_C5_to_C8_MCA_fraction = F.stream.recovered_C5_to_C8_MCA_fraction
pelargonic_acid_rich_fraction = F.stream.pelargonic_acid_rich_fraction
fatty_acid_blend = F.stream.fatty_acid_blend
crude_methanol = F.stream.crude_methanol
crude_glycerol = F.stream.crude_glycerol
sodium_oleate = F.stream.sodium_oleate_product
crude_vegetable_oil = F.stream.crude_vegetable_oil 
fresh_tungsten_catalyst = F.stream.fresh_tungsten_catalyst
fresh_cobalt_catalyst = F.stream.fresh_cobalt_catalyst_stream
fresh_solvent = F.stream.solvent_for_extraction
fresh_HP = F.stream.fresh_HP
citric_acid = F.stream.citricacid_for_degumming
polystyrene_based_catalyst = F.stream.polystyrene_based_catalyst
conc_hydrochloric_acid = F.stream.conc_hydrochloric_acid
calcium_chloride = F.stream.calcium_chloride_for_cat_sep
Resin = F.stream.polystyrene_based_catalyst
Liquid_HCl = F.stream.Liquid_HCl
natural_gas = F.stream.natural_gas


#functions to calculate all the metrics
def get_MPSP():
    azelaic_acid.price  = 0
    MPSP = azelaic_acid.price = azelaic_acid_tea.solve_price(azelaic_acid)
    return MPSP
def get_annual_factor(): return azelaic_acid_tea.operating_days*24
def get_total_yield(): return azelaic_acid.F_mass*get_annual_factor()/1e6 #To get total yield in .10^6 Kg
def get_purity(): return azelaic_acid.imass['Azelaic_acid']/azelaic_acid.F_mass
def get_overall_TCI(): return azelaic_acid_tea.TCI/1e6        
def get_material_cost(): return azelaic_acid_tea.material_cost /1e6 #includes the money spent on ash disposal 
def get_operating_cost(): return azelaic_acid_tea.AOC/1e6
def get_system_heating_demand(): 
    #The boilerturbogenerator is designed to satisfy all the system heating demand
    heat_duty_handled_by_BT901 = -aa_baseline.operating_hours * sum([i.duty for i in F.BT901.heat_utilities if i.flow * i.duty > 0])
    return  heat_duty_handled_by_BT901 #KJ/yr
def get_system_cooling_water_duty(): 
    CT901_duty = -aa_baseline.operating_hours * sum([i.duty for i in F.CT901.heat_utilities if i.flow * i.duty < 0])
    CW901_duty = -aa_baseline.operating_hours * sum([i.duty for i in F.CW901.heat_utilities if i.flow * i.duty < 0])
    total_cooling_duty = -(CT901_duty+CW901_duty)
    return total_cooling_duty #KJ/yr
def get_NPV(): return azelaic_acid_tea.NPV
def solve_IRR(): return azelaic_acid_tea.solve_IRR()

# #######################################################################################################################3

# Environmental impact calculation
products = [azelaic_acid,
            recovered_C5_to_C8_MCA_fraction,
            pelargonic_acid_rich_fraction,
            fatty_acid_blend,
            crude_glycerol,
            crude_methanol,
            sodium_oleate]

#Input and recycle streams in the system do not count as emissions
extra_streams = [F.resin_for_hydrolysis_1,
                 F.resin_for_hydrolysis_2,
                 F.resin_for_hydrolysis_3,
                 F.regeneration_acid_to_HydrolysisSystem]

#The emissions originate from R300, D403,D404,D601,HX608,ADP901,U902,R902,BT901,CT901,PWT901
#Emissions were assumed to be carbon emissions
emissions = [i for i in F.stream if i.source and not i.sink and not i in products and not i in extra_streams]
#Total emissions (TE)
get_total_emissions_GWP = TE_GWP = lambda: sum([stream.get_atomic_flow('C') for stream in emissions]) * chems.CO2.MW / azelaic_acid.F_mass
#End of life emissions (EOL) emissions are indirect emissions
get_EOL_GWP = EOL_GWP = lambda: sum([stream.get_atomic_flow('C')*chems.CO2.MW for stream in products])/azelaic_acid.F_mass
#Direct emissions are emissions without the EOL emissions
Direct_emmisions = get_total_direct_emissions_GWP = lambda: get_total_emissions_GWP() - get_EOL_GWP()

#Boilerturbogenerator also contributes to direct emissions
#The ratio of BT direct emissions and total direct emissions is the same as
#ratio of BT total emissions and overall total emissions
Direct_emmisions_from_BT = get_total_direct_BT_emissions_GWP = lambda: (((sum([i.get_atomic_flow('C') for i in F.BT901.outs])*chems.CO2.MW/ azelaic_acid.F_mass))/get_total_emissions_GWP())*get_total_direct_emissions_GWP()
Direct_emmisions_without_BT = get_total_non_BT_direct_emissions_GWP = lambda: get_total_direct_emissions_GWP() - get_total_direct_BT_emissions_GWP()


#Material related impacts
#TODO: check whether these impacts are direct or indirect
get_feedstock_GWP = lambda: aa_baseline.get_material_impact(crude_vegetable_oil,key = 'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
get_internal_materials_impact = lambda:(aa_baseline.get_material_impact(Resin,key ='GWP100')+aa_baseline.get_material_impact(Liquid_HCl,key ='GWP100'))/aa_baseline.get_mass_flow(azelaic_acid)
#TODO: check again if the materials are right
get_other_materials_impact = lambda: (aa_baseline.get_total_feeds_impact('GWP100')/aa_baseline.get_mass_flow(azelaic_acid))-get_feedstock_GWP() +get_internal_materials_impact()

#Electrcity related impacts
#Electricity is produced using steam and natural gas/or is directly purchased from the grid

#Net electricity use includes net electricity requirement of all the units including the facilities
#Includes production and consumption both
get_net_electricity_use = lambda: sum(i.power_utility.rate for i in aa_baseline.units)
get_net_electricity_use_GWP = lambda: get_net_electricity_use()*Utility_GWP_factors['Electricity']/azelaic_acid.F_mass


#Since in our system, BT901 is satisfying the system electricity demand, 
#all the electricity in the system is supplied by BT901
#natural gas is bought if there is not enough heat generated from the feeds to generate enough electricity
#In this case no electricity from the grid is purchased
#Therefore, total electricity demand can be calculated from the BT901 power utility rate
get_total_electricity_demand = get_electricity_use = lambda: -F.BT901.power_utility.rate 
#Total electricity demand required for cooling 
#Can be calulated from the power utility rates of the cooling tower and the chilled water package
get_cooling_electricity_demand = lambda: F.CT901.power_utility.rate + F.CW901.power_utility.rate 
#Total non cooling electricity demands for pumps,agitators etc for all the units and facilities both
get_non_cooling_electricity_demand = lambda: get_total_electricity_demand() - get_cooling_electricity_demand()
#Fractions of electricity use for both cooling and non cooling
get_elec_frac_cooling_demand = lambda: get_cooling_electricity_demand()/get_total_electricity_demand()
get_elec_frac_non_cooling_demand = lambda:get_non_cooling_electricity_demand()/get_total_electricity_demand()


#Impact of the total natural gas consumed in the process via the boiler turbogenerator
get_ng_GWP = lambda: aa_baseline.get_material_impact(natural_gas,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
#The steam is generated using heat obtained from burning the feed and some natural gas
#The generation of steam causes direct emissions (resulting from burning the feeds)
#Therefore,we attribute the impact of steam to NG use and emissions resulting from generation of steam
get_total_steam_GWP = lambda: get_ng_GWP() + get_total_direct_BT_emissions_GWP()

#The steam utilities (heat and electricity)
#steam produced satisfies all of the heating demand
#steam produced can also be used to satisfy some portion of the heating and cooling demand

#steam used for heating
get_BT_steam_kJph_heating = lambda: sum([i.duty for i in F.BT901.steam_utilities]) 
#steam used for generating electricity
#This electricity can be for both cooling and non-cooling purposes)
#since the efficiency of the turbogen is not 100%, not all the steam gets used in the production of electricity
get_BT_steam_kJph_turbogen = lambda: KJpersec_to_KJhr*F.BT901.electricity_demand/F.BT901.turbogenerator_efficiency 
get_BT_steam_kJph_total = lambda: get_BT_steam_kJph_heating() + get_BT_steam_kJph_turbogen()

#fractions of total steam used for both heating and electricity gen
get_steam_frac_heating = lambda: get_BT_steam_kJph_heating()/get_BT_steam_kJph_total()
get_steam_frac_turbogen = lambda: get_BT_steam_kJph_turbogen()/get_BT_steam_kJph_total()

#Out of the steam used for electricity generation in the turbogen
#Some is used for cooling and other is used for non-cooling purposes
get_steam_frac_electricity_cooling = lambda: get_steam_frac_turbogen() * get_elec_frac_cooling_demand()
get_steam_frac_electricity_non_cooling = lambda: get_steam_frac_turbogen() * (1-get_elec_frac_cooling_demand())

get_heating_demand_GWP = lambda: get_steam_frac_heating()*get_total_steam_GWP()
get_cooling_demand_GWP = lambda: get_steam_frac_electricity_cooling()*get_total_steam_GWP()

get_electricity_demand_cooling_GWP = lambda:get_steam_frac_electricity_cooling()*get_total_steam_GWP() + get_elec_frac_cooling_demand()*get_net_electricity_use_GWP()
get_electricity_demand_non_cooling_GWP = lambda: get_steam_frac_electricity_non_cooling() * get_total_steam_GWP() + get_elec_frac_non_cooling_demand()*get_net_electricity_use_GWP() 

def get_net_GWP():
    GWP_breakdown_sum = get_feedstock_GWP() + get_other_materials_impact() + \
                        get_heating_demand_GWP() +\
                        get_electricity_demand_cooling_GWP()+\
                        get_electricity_demand_non_cooling_GWP()+\
                        get_total_non_BT_direct_emissions_GWP()
    PA_GWP = aa_baseline.get_material_impact(F.pelargonic_acid_rich_fraction,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
    C5_C8_fraction_GWP = aa_baseline.get_material_impact(recovered_C5_to_C8_MCA_fraction,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
    fatty_acid_blend_GWP = aa_baseline.get_material_impact(fatty_acid_blend,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
    glycerol_GWP = aa_baseline.get_material_impact(crude_glycerol,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
    methanol_GWP = aa_baseline.get_material_impact(crude_methanol,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
    sodium_oleate_GWP = aa_baseline.get_material_impact(sodium_oleate,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
    products_GWP_sum = PA_GWP +C5_C8_fraction_GWP +fatty_acid_blend_GWP+fatty_acid_blend_GWP +glycerol_GWP +methanol_GWP+sodium_oleate_GWP
    return GWP_breakdown_sum - products_GWP_sum


# ######################################################################################################################3
metrics = [
    Metric('MPSP', get_MPSP, '$/kg'),
    Metric('Total product yield', get_total_yield, '10^6 kg/yr'),
    Metric('Product purity', get_purity, '%'),
    Metric('Total capital investment', get_overall_TCI, '10^6 $'),
    Metric('Annual operating cost', get_operating_cost, '10^6 $/yr'),
    Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
    Metric('Total', get_system_heating_demand, 'KJ/yr', 'Heating demand'),
    Metric('Total', get_system_cooling_water_duty,'KJ/yr', 'Cooling demand'),
    Metric('NPV', get_NPV, '$', 'TEA'),
    Metric('Internal rate of return', solve_IRR, '%'),
    Metric('Net GWP', get_net_GWP, 'kg CO2-eq/kg', 'LCA'),
]
model = Model(aa_baseline,metrics,  
              exception_hook='warn' #raise
              )

#ranges can be changed using lb_fac and ub_fac
lb_fac = 0.5
ub_fac = 1.5

#price distributions
#TODO: fit historical prices for crude veg oil into different distributions
lb_o = crude_vegetable_oil.price * lb_fac  # Minimum price
ub_o = crude_vegetable_oil.price * ub_fac  # Maximum price
@model.parameter(name='Crude oil price',
                  element=crude_vegetable_oil, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_o, ub_o))
def set_feedstock_price(feedstock_price):
    crude_vegetable_oil.price = feedstock_price


#find multiple lab prices and use those
lb_t = fresh_tungsten_catalyst.price*lb_fac
ub_t = fresh_tungsten_catalyst.price*ub_fac
@model.parameter(name='Tungstic acid catalyst price',
                  element=fresh_tungsten_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_t, ub_t))
def set_tungstencat_price(tungstencat_price):
    fresh_tungsten_catalyst.price = tungstencat_price


#find multiple lab prices and use those
lb_c = fresh_cobalt_catalyst.price*lb_fac
ub_c = fresh_cobalt_catalyst.price*ub_fac
@model.parameter(name='Cobalt acetate catalyst price',
                  element=fresh_cobalt_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_c, ub_c))
def set_cobaltcat_price(fresh_cobalt_catalyst_price):
    fresh_cobalt_catalyst.price = fresh_cobalt_catalyst_price


lb_s = fresh_solvent.price*lb_fac
ub_s = fresh_solvent.price*ub_fac
@model.parameter(name='Heptane solvent price',
                  element=fresh_solvent, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_s, ub_s))
def set_solvent_price(fresh_solvent_price):
    fresh_solvent.price = fresh_solvent_price

#TODO: find multiple lab prices for different concentrations of HP
lb_hp = fresh_HP.price*lb_fac
ub_hp = fresh_HP.price*ub_fac
@model.parameter(name='Hydrogen peroxide price',
                  element=fresh_HP, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_hp, ub_hp))
def set_fresh_HP_price(fresh_HP_price):
    fresh_HP.price = fresh_HP_price

#TODO: find historical prices and distributions
lb_ngp = natural_gas.price*lb_fac
ub_ngp = natural_gas.price*ub_fac
@model.parameter(name='Natural gas price',
                  element=F.BT901,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_ngp,ub_ngp))
def set_natural_gas_price(natural_gas_price):
    F.BT901.natural_gas_price = natural_gas_price


lb_ca = citric_acid.price*lb_fac
ub_ca = citric_acid.price*ub_fac
@model.parameter(name='Citric acid price',
                  element=citric_acid,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_ca, ub_ca))
def set_citric_acid_price(citric_acid_price):
    citric_acid.price = citric_acid_price
    
    
lb_hcl = conc_hydrochloric_acid.price*lb_fac
ub_hcl = conc_hydrochloric_acid.price*ub_fac
@model.parameter(name='Conc HCl price',
                  element=conc_hydrochloric_acid,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_hcl, ub_hcl))
def set_HCl_price(conc_hydrochloric_acid_price):
    conc_hydrochloric_acid.price = conc_hydrochloric_acid_price
    
lb_cc = calcium_chloride.price*lb_fac
ub_cc = calcium_chloride.price*ub_fac
@model.parameter(name='Calcium chloride price',
                  element=calcium_chloride,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_cc, ub_cc))
def set_calcium_chloride_price(calcium_chloride_price):
    calcium_chloride.price = calcium_chloride_price


lb_pbc =  Resin.price*lb_fac
ub_pbc =  Resin.price*ub_fac
@model.parameter(name='Hydrolysis resin price',
                  element=polystyrene_based_catalyst,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_pbc, ub_pbc))
def set_polystyrene_based_catalyst_price(polystyrene_based_catalyst_price):
    polystyrene_based_catalyst.price = polystyrene_based_catalyst_price


# coupled parameters
# Process related parameters
# #Tungstic acid mole ratio
#TODO: tungstic acid m
lb1 = 0.001
ub1 = 0.15
@model.parameter(name='Tungstic acid moles',
                  element=F.M200,
                  kind='coupled',
                  distribution=shape.Uniform(lb1, ub1))
def set_tungstic_acid_moles(X_tam):
    F.unit.M200.specifications[0].args[0] = X_tam

# # #Cobalt acetate catalyst ratio
lb2 = 0.001 
ub2 = 0.015
@model.parameter(name='Cobalt acetate moles',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb2, ub2))
def set_cobalt_acetate_moles(X_cam):
    F.unit.R300.specifications[0].args[1] = X_cam


# Dihydroxylation reaction conversion
lb3 = 0.70
ub3 = 0.99
@model.parameter(name='Dihydroxylation reaction conversion',
                  element=F.R200,
                  kind='coupled',
                  distribution=shape.Uniform(lb3, ub3))
def set_dihydroxylation_conversion(X_dih):
    F.unit.R200.X_dih = X_dih
    


## Oxidative cleavage reaction conversion (F.R300.X_ox_rxn_1)
lb4 = 0.70
ub4 = 0.95
@model.parameter(name='Oxidative cleavage reaction conversion primary',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb4, ub4))
def set_X_oxidativecleavage_conversion1(X_ox_rxn_1):
    F.unit.R300.X_ox_rxn_1 = X_ox_rxn_1
    
    

lb9 = 0.90
ub9 = 0.95
@model.parameter(name='Oxidative cleavage reaction conversion',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb9, ub9))
def set_X_oxidativecleavage_conversion(X_oxidativecleavage):
    F.unit.R300.X_oxidativecleavage = X_oxidativecleavage


##Oxidative cleavage section reactions
lb5 = 0.05
ub5 = 0.30
@model.parameter(name='Decarboxylation reaction conversion',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb5, ub5))
def set_X_decarboxylation_conversion(X_decarboxylation):
    F.unit.R300.X_decarboxylation = X_decarboxylation


##Oxidative cleavage reaction conversion
lb6 = 0.05
ub6 = 0.30
@model.parameter(name='Side reaction conversion',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb6, ub6))
def set_X_side_rxn_conversion(X_side_rxn):
    F.unit.R300.X_side_rxn = X_side_rxn

lb7 = 2
ub7 = 10
@model.parameter(name='Dihydroxylation reaction time',
                  element=F.R200,
                  kind='coupled',
                  distribution=shape.Uniform(lb7, ub7))
def set_dihydroxylation_process_time(tau):
    F.unit.R200.tau = tau
    
lb8 = 3
ub8 = 10
@model.parameter(name='Oxidative reaction time',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb8, ub8))
def set_oxidative_rxn_conversion_time(tau):
     F.unit.R300.tau = tau
    

lb10 = 0.08
ub10 = 0.15
@model.parameter(name='System IRR',
                  element=azelaic_acid_tea.IRR,
                  kind='cost',
                  distribution=shape.Uniform(lb10, ub10)
                  )
def set_irr(X_irr):
    azelaic_acid_tea.IRR = X_irr

lb11 = 0.70
ub11 = 0.90
@model.parameter(name='Turbogen efficiency',
                  element=F.BT901,
                  kind='coupled',
                  distribution=shape.Uniform(lb11, ub11)
                  )
def set_tubeff(X_beff):
    F.BT901.turbogenerator_efficiency = X_beff


lb14 = 3 #current state
ub14 = 20 #optimistic 
@model.parameter(name='Tungstic acid reusability',
                  element=F.M200,
                  kind='coupled',
                  distribution=shape.Uniform(lb14, ub14)
                  )
def set_cta(cycles_of_reuse):
    F.unit.M200.specifications[0].args[1] = cycles_of_reuse
    

# N_samples = 1000
# rule = 'L' # For Latin-Hypercube sampling
# np.random.seed(1234) # For consistent results
# samples = model.sample(N_samples, rule)
# model.load_samples(samples)
# model.evaluate(notify=2)
# model.show()