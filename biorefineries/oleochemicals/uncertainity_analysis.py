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
from biosteam.evaluation import Model, Metric, Parameter
from prices_and_GWP_factors import prices_per_stream,utility_prices,GWP_per_stream,Utility_GWP_factors,isolated_para_dists,coupled_para_dist,environmental_facs_dist,GWP_factors_per_Kg
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
crude_vegetable_oil = F.stream.crude_vegetable_oil 
fresh_tungsten_catalyst = F.stream.fresh_tungsten_catalyst
fresh_cobalt_catalyst = F.stream.fresh_cobalt_catalyst_stream
fresh_solvent = F.stream.solvent_for_extraction
fresh_HP = F.stream.fresh_HP
citric_acid = F.stream.citricacid_for_degumming
conc_hydrochloric_acid = F.stream.conc_hydrochloric_acid
calcium_chloride = F.stream.calcium_chloride_for_cat_sep
resin = F.stream.polystyrene_based_catalyst
biodiesel_catalyst = F.stream.catalyst
liquid_HCl = F.stream.Liquid_HCl
natural_gas = F.stream.natural_gas



#functions to calculate all the metrics
def get_MPSP():
    azelaic_acid.price  = 0
    MPSP = azelaic_acid.price = azelaic_acid_tea.solve_price(azelaic_acid)
    return MPSP
def get_annual_factor(): return azelaic_acid_tea.operating_days*24
def get_total_yield(): return azelaic_acid.F_mass*get_annual_factor()/1e6 #To get total yield in .10^6 Kg
def get_purity(): return azelaic_acid.imass['Azelaic_acid']*100/azelaic_acid.F_mass
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
main_product = [azelaic_acid]
products = [recovered_C5_to_C8_MCA_fraction,
            pelargonic_acid_rich_fraction,
            fatty_acid_blend,
            crude_glycerol,
            crude_methanol,
            ]

#Input and recycle streams in the system do not count as emissions
extra_streams = [F.resin_for_hydrolysis_1,
                 F.resin_for_hydrolysis_2,
                 F.resin_for_hydrolysis_3,
                 F.regeneration_acid_to_HydrolysisSystem]

#The emissions originate from R300, D403,D404,D601,HX608,U902,R902,BT901,CT901,PWT901
#Emissions were assumed to be carbon emissions
#TODO: confirm if the ADP outs should be an emission or not, also outs from PWT901 should be there or npt
emissions = [i for i in F.stream if i.source and not i.sink and i not in products and i not in extra_streams and i not in main_product and i not in F.ADP901.outs]

#Total emissions (TE)
get_total_emissions_GWP = TE_GWP = lambda: sum([stream.get_atomic_flow('C') for stream in emissions]) * chems.CO2.MW / azelaic_acid.F_mass
#End of life emissions (EOL) emissions are indirect emissions
#TODO: check if in end of life GWP all products are to be included or not

get_EOL_GWP = EOL_GWP = lambda: sum([stream.get_atomic_flow('C')*chems.CO2.MW for stream in main_product])/azelaic_acid.F_mass
#Direct emissions are emissions without the EOL emissions
Direct_emmisions = get_total_direct_emissions_GWP = lambda: get_total_emissions_GWP() - get_EOL_GWP()

#Boilerturbogenerator also contributes to direct emissions
Direct_emmisions_from_BT = get_total_direct_BT_emissions_GWP = lambda: (sum([i.get_atomic_flow('C') for i in F.BT901.outs])*chems.CO2.MW/ azelaic_acid.F_mass)*get_total_direct_emissions_GWP()/get_total_emissions_GWP()
Direct_emmisions_without_BT = get_total_non_BT_direct_emissions_GWP = lambda: get_total_direct_emissions_GWP() - get_total_direct_BT_emissions_GWP()
#Material related impacts

#TODO: check whether these impacts are direct or indirect
#TODO: your feedstock includes biogenic emissions
get_feedstock_GWP = lambda: aa_baseline.get_material_impact(crude_vegetable_oil,key = 'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
get_internal_materials_impact = lambda:(aa_baseline.get_material_impact(resin,key ='GWP100')+aa_baseline.get_material_impact(liquid_HCl,key ='GWP100'))/aa_baseline.get_mass_flow(azelaic_acid)
get_other_materials_impact = lambda: (aa_baseline.get_total_feeds_impact('GWP100')/aa_baseline.get_mass_flow(azelaic_acid))-get_feedstock_GWP() +get_internal_materials_impact()
get_other_products_impact = lambda: sum([aa_baseline.get_material_impact(i,key = 'GWP100')/aa_baseline.get_mass_flow(azelaic_acid) for i in products])

#Electrcity related impacts
#Electricity is produced using steam and natural gas/or is directly purchased from the grid
#Net electricity use includes net electricity requirement of all the units including the facilities
#Includes production and consumption both
#TODO: add annual factor to account for time
get_net_electricity_use = lambda: sum(i.power_utility.rate for i in aa_baseline.units)
get_net_electricity_use_GWP = lambda: get_net_electricity_use()*Utility_GWP_factors['Electricity']/aa_baseline.get_mass_flow(azelaic_acid)


#Since in our system, BT901 is satisfying the system electricity demand, 
#all the electricity in the system is supplied by BT901
#natural gas is bought if there is not enough heat generated from the feeds to generate enough electricity
#and no electricity from the grid is purchased
#Therefore, total electricity demand can be calculated from the BT901 power utility rate
get_total_electricity_demand = get_electricity_use = lambda: -F.BT901.power_utility.rate 
#Total electricity demand required for cooling 
#Can be calulated from the power utility rates of the cooling tower and the chilled water package
get_cooling_electricity_demand = lambda: F.CT901.power_utility.rate + F.CW901.power_utility.rate 
#Total non cooling electricity demands for pumps,agitators etc for all the units and facilities both
get_non_cooling_electricity_demand = lambda: get_total_electricity_demand() - get_cooling_electricity_demand()
#Fractions of electricity use for both cooling and non cooling

#TODO: are these both purchased electricity (power units) that is used to satisfy the SCD,
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
get_BT_steam_kJph_total = lambda: get_BT_steam_kJph_heating() + get_BT_steam_kJph_turbogen() #total steam generation

#fractions of total steam used for both heating and electricity gen
get_steam_frac_heating = lambda: get_BT_steam_kJph_heating()/get_BT_steam_kJph_total()
get_steam_frac_turbogen = lambda: get_BT_steam_kJph_turbogen()/get_BT_steam_kJph_total()

#Out of the steam used for electricity generation in the turbogen
#Some is used for cooling and other is used for non-cooling purposes (like pumps, agitator etc)
get_steam_frac_electricity_cooling = lambda: get_steam_frac_turbogen() * get_elec_frac_cooling_demand()
get_steam_frac_electricity_non_cooling = lambda: get_steam_frac_turbogen() * (1-get_elec_frac_cooling_demand())

get_heating_demand_GWP = lambda: get_steam_frac_heating()*get_total_steam_GWP()
get_steam_electricity_demand_cooling_GWP = lambda: get_steam_frac_electricity_cooling()*get_total_steam_GWP()
get_cooling_demand_GWP = lambda: get_steam_electricity_demand_cooling_GWP() + get_elec_frac_cooling_demand()*get_net_electricity_use_GWP()
get_electricity_demand_non_cooling_GWP = lambda: get_steam_frac_electricity_non_cooling() * get_total_steam_GWP() + get_elec_frac_non_cooling_demand()*get_net_electricity_use_GWP() 

#TODO: check again if the below is right
def get_net_GWP():
    GWP_breakdown_sum = sum([get_feedstock_GWP(),-get_other_products_impact(),
                             get_other_materials_impact(),
                             get_heating_demand_GWP(),
                             get_cooling_demand_GWP(),#includes steam used for electricity production 
                             # and other electricity used for cooling purposes 
                            get_electricity_demand_non_cooling_GWP(),#TODO; not very significant
                            get_total_non_BT_direct_emissions_GWP() #TODO; not very significant
                            ])
    return GWP_breakdown_sum 

def get_net_GWP_2():
    GWP_breakdown_2 = sum([get_feedstock_GWP(),-get_other_products_impact(),
                             get_other_materials_impact(),
                             get_heating_demand_GWP(),
                             get_cooling_demand_GWP(),#includes steam used for electricity production 
                             # and other electricity used for cooling purposes 
                            get_electricity_demand_non_cooling_GWP(),
                            get_total_emissions_GWP(),
                            -get_EOL_GWP(),-Direct_emmisions_from_BT()])
    return GWP_breakdown_2

def get_net_GWP_3():
    return(sum([get_feedstock_GWP(),
                get_other_materials_impact(),
                get_ng_GWP(),
                get_net_electricity_use_GWP(),
                Direct_emmisions(),
                -get_other_products_impact()]))


# ######################################################################################################################3
all_metrics = [
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

model = Model(aa_baseline,metrics = all_metrics,exception_hook= 'warn')    

#TODO; remove the number of runs argument, return the model
def conduct_uncertainity_analysis(system = aa_baseline, #name of the system
                                  metrics = all_metrics, #list of Metric objects
                                  number_of_runs = 3, #number of 
                                  exception_hook = 'warn',#raise/warn
                                  rule = 'L', #for For Latin-Hypercube sampling
                                  notify_runs = 50, #runs after which you are notified 
                                  indicator = 'GWP100',
                                  feedstock_type = 'HoySoy_oil' #can also provide 'HoSun_oil
                                 ):
    model = Model(aa_baseline,metrics = metrics,exception_hook= 'warn')     
    ## price distributions
    @model.parameter(name='Crude oil price',
                      element=crude_vegetable_oil, 
                      kind='isolated',
                      units='USD/kg',
                      distribution= isolated_para_dists['Crude oil price'][feedstock_type])
    def set_feedstock_price(feedstock_price):
            crude_vegetable_oil.price = feedstock_price
    
    @model.parameter(name='Tungstic acid price',
                      element=fresh_tungsten_catalyst, 
                      kind='isolated', 
                      units='USD/kg',
                      distribution= isolated_para_dists['Tungstic acid price'])
    def set_tungstencat_price(tungstencat_price):
        fresh_tungsten_catalyst.price = tungstencat_price
    
    
    @model.parameter(name='Cobalt acetate catalyst price',
                      element=fresh_cobalt_catalyst, 
                      kind='isolated', 
                      units='USD/kg',
                      distribution= isolated_para_dists['Cobalt acetate price'])
    def set_cobaltcat_price(fresh_cobalt_catalyst_price):
        fresh_cobalt_catalyst.price = fresh_cobalt_catalyst_price
    
    @model.parameter(name='Heptane solvent price',
                      element=fresh_solvent, 
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Heptane solvent price'])
    def set_solvent_price(fresh_solvent_price):
        fresh_solvent.price = fresh_solvent_price
    
    @model.parameter(name='Hydrogen peroxide price',
                      element=fresh_HP, 
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Hydrogen peroxide price'])
    def set_fresh_HP_price(fresh_HP_price):
        fresh_HP.price = fresh_HP_price
    
    @model.parameter(name='Natural gas price',
                      element=F.BT901,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Natural gas price'])
    def set_natural_gas_price(natural_gas_price):
        F.BT901.natural_gas_price = natural_gas_price
    
    @model.parameter(name='Citric acid price',
                      element=citric_acid,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Citric acid price'])
    def set_citric_acid_price(citric_acid_price):
        citric_acid.price = citric_acid_price
    
    
    @model.parameter(name='Conc HCl price',
                      element=conc_hydrochloric_acid,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Conc HCl price'])
    def set_HCl_price(conc_hydrochloric_acid_price):
        conc_hydrochloric_acid.price = conc_hydrochloric_acid_price
    
    
    @model.parameter(name='Calcium chloride price',
                  element=calcium_chloride,
                  kind='isolated',
                  units='USD/kg',
                  distribution=isolated_para_dists['Calcium chloride price'])
    def set_calcium_chloride_price(calcium_chloride_price):
        calcium_chloride.price = calcium_chloride_price
    
    
    @model.parameter(name='Hydrolysis resin price',
                      element=resin,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Hydrolysis resin price'])
    def set_polystyrene_based_catalyst_price(resin_price):
        resin.price = resin_price 
    
    #Product prices    
    @model.parameter(name='Crude glycerol price',
                      element=crude_glycerol,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Crude glycerol price'])
    def set_crude_glycerol_price(crude_glycerol_price):
        crude_glycerol.price = crude_glycerol_price
        
    @model.parameter(name='Pelargonic acid price',
                      element=pelargonic_acid_rich_fraction,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Pelargonic acid price'])
    def set_pelargonic_acid_price(pelargonic_acid_price):
        pelargonic_acid_rich_fraction.price = pelargonic_acid_price            
        
    @model.parameter(name='Fatty acid blend price',
                      element=fatty_acid_blend,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Fatty acid blend price'])
    def set_fatty_acid_blend_price(fatty_acid_blend_price):
        fatty_acid_blend.price = fatty_acid_blend_price        
    
    # # coupled parameters
    # # Process related parameters
    # # #Tungstic acid mole ratio
    @model.parameter(name='Tungstic acid moles',
                      element=F.M200,
                      kind='coupled',
                      distribution= coupled_para_dist['Tungstic acid moles'])
    def set_tungstic_acid_moles(X_tam):
        F.unit.M200.specifications[0].args[0] = X_tam
    
    @model.parameter(name='Cobalt acetate moles',
                      element=F.R300,
                      kind='coupled',
                      distribution= coupled_para_dist['Cobalt acetate moles'])
    def set_cobalt_acetate_moles(X_cam):
        F.unit.R300.specifications[0].args[1] = X_cam
    
    
    ## Oxidative cleavage reaction conversion (F.R300.X_ox_rxn_1)
    @model.parameter(name='Oxidative cleavage reaction conversion primary',
                      element=F.R300,
                      kind='coupled',
                      distribution=coupled_para_dist['Oxidative cleavage reaction conversion primary'])
    def set_X_oxidativecleavage_conversion1(X_ox_rxn_1):
        F.unit.R300.X_ox_rxn_1 = X_ox_rxn_1
        
    
    @model.parameter(name='Oxidative cleavage reaction conversion',
                      element=F.R300,
                      kind='coupled',
                      distribution=coupled_para_dist['Oxidative cleavage reaction conversion'])
    def set_X_oxidativecleavage_conversion(X_oxidativecleavage):
        F.unit.R300.X_oxidativecleavage = X_oxidativecleavage
    
    
    ##Oxidative cleavage section reactions
    @model.parameter(name='Decarboxylation reaction conversion',
                      element=F.R300,
                      kind='coupled',
                      distribution=coupled_para_dist['Decarboxylation reaction conversion'])
    def set_X_decarboxylation_conversion(X_decarboxylation):
        F.unit.R300.X_decarboxylation = X_decarboxylation
    
    
    ##Oxidative cleavage reaction conversion
    @model.parameter(name='Side reaction conversion',
                      element=F.R300,
                      kind='coupled',
                      distribution=coupled_para_dist['Side reaction conversion'])
    def set_X_side_rxn_conversion(X_side_rxn):
        F.unit.R300.X_side_rxn = X_side_rxn
        
    @model.parameter(name='Dihydroxylation reaction time',
                      element=F.R200,
                      kind='coupled',
                      distribution=coupled_para_dist['Dihydroxylation reaction time'])
    def set_dihydroxylation_process_time(tau):
        F.unit.R200.tau = tau
        
    @model.parameter(name='Oxidative reaction time',
                      element=F.R300,
                      kind='coupled',
                      distribution=coupled_para_dist['Oxidative reaction time'])
    def set_oxidative_rxn_conversion_time(tau):
          F.unit.R300.tau = tau
          
    @model.parameter(name='Turbogen efficiency',
                      element=F.BT901,
                      kind='coupled',
                      distribution=coupled_para_dist['Turbogen efficiency'])
    def set_tubeff(X_tubeff):
        F.BT901.turbogenerator_efficiency = X_tubeff
        
    @model.parameter(name='Boiler efficiency',
                          element=F.BT901,
                          kind='coupled',
                          distribution=coupled_para_dist['Boiler efficiency'])
    def set_beff(X_beff):
            F.BT901.boiler_efficiency = X_beff
        
    @model.parameter(name='Tungstic acid reusability',
                      element=F.M200,
                      kind='coupled',
                      distribution=coupled_para_dist['Tungstic acid reusability']
                      )
    def set_cta(cycles_of_reuse):
        F.unit.M200.specifications[0].args[1] = cycles_of_reuse
                                    
    #environmental_facs_dist
    @model.parameter(name='crude oil GWP',
                      element=crude_vegetable_oil,
                      kind='isolated',
                      distribution=environmental_facs_dist['Oil_GWP'][feedstock_type]
                      )
    def crude_vegetable_oil_gwp(crude_vegetable_oil_gwp):
        crude_vegetable_oil.characterization_factors[indicator] = crude_vegetable_oil_gwp
        
    @model.parameter(name='tungstic acid GWP',
                      element=fresh_tungsten_catalyst,
                      kind='isolated',
                      distribution=environmental_facs_dist['Tungstic_acid']
                      )
    def tungstic_acid_gwp(tungstic_acid_gwp):
        fresh_tungsten_catalyst.characterization_factors[indicator] = tungstic_acid_gwp 
        
    @model.parameter(name='cobalt acetate GWP',
                      element=fresh_cobalt_catalyst,
                      kind='isolated',
                      distribution=environmental_facs_dist['Cobalt_acetate']
                      )
    def cobalt_acetate_gwp(cobalt_acetate_gwp):
        fresh_cobalt_catalyst.characterization_factors[indicator] = cobalt_acetate_gwp 
    
    @model.parameter(name='hcl GWP',
                      element=liquid_HCl,
                      kind='isolated',
                      distribution=environmental_facs_dist['HCl']
                      )
    def hcl_gwp(hcl_gwp):
        liquid_HCl.characterization_factors[indicator]  =  hcl_gwp
    
    @model.parameter(name='resin GWP',
                      element=resin,
                      kind='isolated',
                      distribution=environmental_facs_dist['Resin']
                      )
    def Resin_gwp(Resin_gwp):
        resin.characterization_factors[indicator]  = Resin_gwp     
    
    @model.parameter(name='Sodium methoxide GWP',
                      element=biodiesel_catalyst,
                      kind='isolated',
                      distribution=environmental_facs_dist['Sodium methoxide']
                      )
    def biodiesel_catalyst_gwp(biodiesel_catalyst_gwp):
        biodiesel_catalyst.characterization_factors[indicator]  = 0.25*biodiesel_catalyst_gwp + 0.75*GWP_factors_per_Kg['Methanol']
        
    @model.parameter(name='heptane GWP',
                      element=fresh_solvent,
                      kind='isolated',
                      distribution=environmental_facs_dist['Heptane']
                      )
    def heptane_solvent_gwp(heptane_solvent_gwp):
        fresh_solvent.characterization_factors[indicator]  = heptane_solvent_gwp      
        
    @model.parameter(name='fatty acid blend GWP',
                      element=fatty_acid_blend,
                      kind='isolated',
                      distribution=environmental_facs_dist['Fatty_acid_blend']
                      )
    def fatty_acid_blend_gwp(fatty_acid_blend_gwp):
        fatty_acid_blend.characterization_factors[indicator]  = fatty_acid_blend_gwp     
        
     
    @model.parameter(name='C5_C8 fraction GWP',
                      element=recovered_C5_to_C8_MCA_fraction,
                      kind='isolated',
                      distribution=environmental_facs_dist['C5_C8_fraction']
                      )
    def recovered_C5_to_C8_MCA_fraction_gwp(recovered_C5_to_C8_MCA_fraction_gwp):
        recovered_C5_to_C8_MCA_fraction.characterization_factors[indicator]  = recovered_C5_to_C8_MCA_fraction_gwp    
    
    @model.parameter(name='pelargonic acid GWP',
                      element=pelargonic_acid_rich_fraction,
                      kind='isolated',
                      distribution=environmental_facs_dist['Pelargonic_acid']
                      )
    def pelargonic_acid_rich_fraction_gwp(pelargonic_acid_rich_fraction_gwp):
        pelargonic_acid_rich_fraction.characterization_factors[indicator]  = pelargonic_acid_rich_fraction_gwp    
        
    @model.parameter(name='cooling tower chems GWP',
                      element=F.cooling_tower_chemicals,
                      kind='isolated',
                      distribution=environmental_facs_dist['Cooling_tower_chemicals']
                      )
    def cooling_tower_chemicals_gwp(cooling_tower_chemicals_gwp):
        F.cooling_tower_chemicals.characterization_factors[indicator]  = cooling_tower_chemicals_gwp      
        
    @model.parameter(name='boiler chems GWP',
                      element=F.boiler_chems,
                      kind='isolated',
                      distribution=environmental_facs_dist['Boiler_chems']
                      )
    def boiler_chems_gwp(boiler_chems_gwp):
        F.boiler_chems.characterization_factors[indicator]  = boiler_chems_gwp      
        
    @model.parameter(name='lime boiler GWP',
                      element=F.lime_boiler,
                      kind='isolated',
                      distribution=environmental_facs_dist['Lime_boiler']
                      )
    def lime_boiler_gwp(lime_boiler_gwp):
        F.lime_boiler.characterization_factors[indicator]  = lime_boiler_gwp      
        

    if number_of_runs > 0:
        rule = rule
        np.random.seed(1234) # For consistent results
        samples = model.sample(number_of_runs,rule)
        model.load_samples(samples)
        model.evaluate(notify=50) #notify each time after 2 runs
        model.show()
        model.table.to_excel('model_table.xlsx')
        df_rho,df_p = model.spearman_r()
        df_rho.to_excel('df_rho.xlsx')
    else:
        model.show()
        