# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:58:21 2023
@author: lavanya
"""
from biorefineries.oleochemicals.systems_baseline import F
import pandas as pd
from biorefineries.oleochemicals.chemicals_baseline import chems
from biorefineries.oleochemicals.system_simulate import aa_baseline,azelaic_acid_tea,aa_baseline_groups
from biosteam.evaluation import Model, Metric
from prices_and_GWP_factors import Utility_GWP_factors,isolated_para_dists,coupled_para_dist,environmental_facs_dist,GWP_factors_per_Kg,correlation_based_bulk_prices
from warnings import filterwarnings; filterwarnings('ignore')
import chaospy
from chaospy import distributions as shape
import numpy as np
from biorefineries.oleochemicals.system_simulate import azelaic_acid_tea,aa_baseline,aa_baseline_groups
import biosteam as bst
from biosteam import report
from biosteam import settings
from pathlib import Path

#unit conversions
KJpersec_to_KJhr = 3600 

#shortform stream names for ease
azelaic_acid = F.stream.azelaic_acid_product_stream
recovered_C5_to_C9_MCA_fraction = F.stream.recovered_C5_to_C9_MCA_fraction
pelargonic_acid_rich_fraction = F.stream.pelargonic_acid_rich_fraction
fatty_acid_blend = F.stream.fatty_acid_blend
crude_methanol = F.stream.crude_methanol
crude_glycerol = F.stream.crude_glycerol
crude_vegetable_oil = F.stream.crude_vegetable_oil 
fresh_tungsten_catalyst = F.stream.fresh_tungsten_catalyst
fresh_cobalt_catalyst = F.stream.fresh_cobalt_catalyst_stream
fresh_solvent = F.stream.solvent_for_extraction
fresh_HP = F.stream.fresh_HP
resin = F.stream.polystyrene_based_catalyst
biodiesel_catalyst = F.stream.catalyst
conc_HCl_for_resin_wash = F.stream.conc_HCl_for_resin_wash
natural_gas = F.stream.natural_gas


#functions to calculate all the metrics
def get_MPSP(): # $/kg of AA
    MPSP = azelaic_acid_tea.solve_price(azelaic_acid)
    return MPSP

def get_MFPP():    
    MFPP = azelaic_acid_tea.solve_price(crude_vegetable_oil)
    return MFPP

def get_total_yield(): return aa_baseline.get_mass_flow(azelaic_acid)/1e6 #To get total yield in .10^6 Kg
def get_purity(): return azelaic_acid.imass['Azelaic_acid']*100/azelaic_acid.F_mass
def get_material_cost(): return azelaic_acid_tea.material_cost/1e6 #includes the money spent on ash disposal and brine disposal, [USD/yr]
def get_utility_cost(): return aa_baseline.utility_cost/1000
def get_AOC(): return azelaic_acid_tea.AOC/1e6 #[mil USD/yr]
def get_TCI(): return azelaic_acid_tea.TCI/1e6 #[mil USD]  
def get_FCI(): return azelaic_acid_tea.FCI/1e6 #[mil USD]      
def get_ROI(): return azelaic_acid_tea.ROI #[USD/yr]
def get_IRR(): return azelaic_acid_tea.solve_IRR()
def get_DPI(): return azelaic_acid_tea.DPI/1e6 #[mil USD]        
def get_TDC(): return azelaic_acid_tea.TDC/1e6 #[mil USD]        
def get_total_production_cost():return azelaic_acid_tea.total_production_cost([F.azelaic_acid_product_stream])/1e6

def get_system_heating_demand(): 
    #The boilerturbogenerator is designed to satisfy all the system heating demand
    total_heating_demand = sum([sum([i.duty for i in unit.heat_utilities
                                     if i.flow > 0 and i.duty > 0])
                                for unit in aa_baseline.units])*aa_baseline.operating_hours/1e9 #MJ
    return total_heating_demand #10^6*MJ/yr    
    
def get_system_cooling_demand(): 
    total_heating_demand = -1*sum([sum([i.duty for i in unit.heat_utilities
                                     if i.flow > 0 and i.duty < 0])
                                for unit in aa_baseline.units])*aa_baseline.operating_hours/1e9 #MJ
    return total_heating_demand #10^6*MJ/yr    
def solve_NPV(): return azelaic_acid_tea.NPV/1e6
def get_amt_of_system_makeup_water(): return F.system_makeup_water.F_mass
def electricity_consumption(): return aa_baseline.get_electricity_consumption()

#metrics to support facts for the system performance under uncertainty section
#function to estimate the fraction splits for each of the metrics
def df_related_fraction(): 
    df_unit_groups = bst.UnitGroup.df_from_groups(unit_groups = aa_baseline_groups,fraction = True)
    
    df_unit_groups.index = [100,200,300,400,500,600,900,700,800]

    df_unit_groups = df_unit_groups.sort_index()    
    df_unit_groups.index = [
                            'Transesterification',
                            'Dihydroxylation',
                            'Oxidative cleavage',
                            'Catalyst recovery',
                            'Pelargonic acid and C5-C8 fraction recovery',
                            'Azelaic acid and heavy tails recovery',
                            'Boilerturbogenerator',
                            'Wastewater treatment',
                            'Other Facilities (cooling tower,air distribution and cleaning auxiliaries)'
                            ]
    return df_unit_groups

# #######################################################################################################################3
# Environmental impact calculation
main_product = [azelaic_acid]
products = [recovered_C5_to_C9_MCA_fraction,
            pelargonic_acid_rich_fraction,
            fatty_acid_blend,
            crude_glycerol,
            crude_methanol,
            ]

#Input and recycle streams in the system do not count as emissions
extra_streams = [F.resin_for_hydrolysis_1,
                 F.resin_for_hydrolysis_2,
                 F.resin_for_hydrolysis_3,
                 F.regeneration_acid_to_HydrolysisSystem,]

#The emissions originate from R300, D403,D404,D601,HX608,U902,R902,BT801,CT901,PWT901
#Emissions were assumed to be carbon emissions
#assumed that carbon related compounds in brine get passively oxidized    

emissions = [i for i in aa_baseline.products if i not in extra_streams and i not in main_product and i not in products and i not in F.ADP701.outs]

#Total emissions (TE) 
emissions_GWP = TE_GWP = lambda: sum([stream.get_atomic_flow('C') for stream in emissions]) * chems.CO2.MW / azelaic_acid.F_mass

#End of life emissions (EOL) emissions are indirect emissions
# --- 1. EOL GWP from azelaic acid (main product) ---
EOL_GWP_main_product = lambda: azelaic_acid.get_atomic_flow('C') * chems.CO2.MW / azelaic_acid.F_mass

# --- 2. EOL GWP from coproducts (e.g., pelargonic acid, C5–C9 fraction, etc.) ---
EOL_GWP_coproducts = lambda: sum([stream.get_atomic_flow('C') for stream in products]) * chems.CO2.MW / azelaic_acid.F_mass

# --- 3. Combined EOL emissions (total cradle-to-grave end-of-life burden) ---
EOL_emissions = lambda: EOL_GWP_main_product() + EOL_GWP_coproducts()
direct_emmisions = get_total_direct_emissions_GWP = lambda: emissions_GWP() + EOL_emissions()

#Boilerturbogenerator also contributes to direct emissions #(unitless)
direct_emmisions_from_BT = get_total_direct_BT_emissions_GWP = lambda: (sum([i.get_atomic_flow('C') for i in F.BT801.outs])*chems.CO2.MW/ azelaic_acid.F_mass)*get_total_direct_emissions_GWP()/emissions_GWP()
direct_emmisions_without_BT = get_total_non_BT_direct_emissions_GWP = lambda: get_total_direct_emissions_GWP() - get_total_direct_BT_emissions_GWP()

#Impact of the total natural gas consumed in the process via the boiler turbogenerator
get_ng_GWP = lambda: aa_baseline.get_material_impact(natural_gas,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid) #annual GWP of ng/annual flowrate of AA

#total steam production
KW_to_KJ = 3600
net_electricity_consumption = lambda: sum([u.power_utility.consumption for u in aa_baseline.units])/ azelaic_acid.F_mass
total_electricity_produced_by_BT_using_steam = lambda:F.BT801.power_utility.production
net_electricity_purchased = lambda: net_electricity_consumption()-total_electricity_produced_by_BT_using_steam()
net_electricity_purchased_GWP = lambda: net_electricity_purchased()*Utility_GWP_factors['Electricity']

total_steam_GWP = lambda: get_ng_GWP()+ get_total_direct_BT_emissions_GWP()
cooling_electricity_demand = lambda: (F.CT701.power_utility.rate + F.CW701.power_utility.rate)/ azelaic_acid.F_mass
electricity_frac_cooling = lambda: cooling_electricity_demand()/net_electricity_consumption()
electricity_frac_non_cooling = lambda: 1 - electricity_frac_cooling()

BT_steam_kJph_heating = lambda: sum([i.duty for i in F.BT801.steam_utilities])/azelaic_acid.F_mass #kJ/kg of AA
BT_steam_kJph_turbogen = lambda: KW_to_KJ*total_electricity_produced_by_BT_using_steam()/F.BT801.turbogenerator_efficiency/azelaic_acid.F_mass #kJ/kg of AA
BT_steam_kJph_total = lambda: BT_steam_kJph_heating() + BT_steam_kJph_turbogen()
steam_frac_heating = lambda: BT_steam_kJph_heating()/BT_steam_kJph_total()
steam_frac_turbogen= lambda: BT_steam_kJph_turbogen()/ BT_steam_kJph_total()
steam_frac_cooling= lambda: steam_frac_turbogen()*electricity_frac_cooling()
steam_frac_electricity_non_cooling = lambda: steam_frac_turbogen() * (1-electricity_frac_cooling())

heating_demand_GWP = lambda: steam_frac_heating() * total_steam_GWP() 
cooling_demand_GWP = lambda: steam_frac_cooling() * total_steam_GWP() + electricity_frac_cooling()*net_electricity_purchased_GWP()
electricity_demand_non_cooling_GWP = lambda: steam_frac_electricity_non_cooling()* total_steam_GWP() + electricity_frac_non_cooling()*net_electricity_purchased_GWP()

#Material related impacts
get_feedstock_GWP = lambda: aa_baseline.get_material_impact(crude_vegetable_oil,key = 'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
get_internal_materials_impact = lambda:(aa_baseline.get_material_impact(resin,key ='GWP100')+aa_baseline.get_material_impact(conc_HCl_for_resin_wash,key ='GWP100'))/aa_baseline.get_mass_flow(azelaic_acid)
get_other_materials_impact = lambda: (aa_baseline.get_total_feeds_impact('GWP100')/aa_baseline.get_mass_flow(azelaic_acid))-get_feedstock_GWP() -get_internal_materials_impact() - get_ng_GWP()

#%% functions to estimate gwp

def get_system_GWP():
    return(sum([get_feedstock_GWP(),
                get_other_materials_impact(),
                get_ng_GWP(),
                net_electricity_purchased_GWP(),
                get_total_direct_emissions_GWP()
                ]))

def get_system_GWP1():
    return(sum([get_feedstock_GWP(),
                get_other_materials_impact(),
                heating_demand_GWP(),
                cooling_demand_GWP(),
                electricity_demand_non_cooling_GWP(),
                direct_emmisions_without_BT() 
                ]))

def get_total_product_market_value():
    return(sum([aa_baseline.get_market_value(azelaic_acid),
                aa_baseline.get_market_value(recovered_C5_to_C9_MCA_fraction),
                aa_baseline.get_market_value(pelargonic_acid_rich_fraction),
                aa_baseline.get_market_value(crude_glycerol),                                
                aa_baseline.get_market_value(fatty_acid_blend),    
                aa_baseline.get_market_value(crude_methanol),
                ]))

#Economic allocation
economic_C5_C9_fraction = aa_baseline.get_market_value(recovered_C5_to_C9_MCA_fraction)*get_system_GWP()/get_total_product_market_value()
economic_pa_fraction = aa_baseline.get_market_value(pelargonic_acid_rich_fraction)*get_system_GWP()/get_total_product_market_value()
economic_glycerol_fraction = aa_baseline.get_market_value(crude_glycerol)*get_system_GWP()/get_total_product_market_value() 
economic_fa_fraction = aa_baseline.get_market_value(fatty_acid_blend)*get_system_GWP()/get_total_product_market_value()          
economic_methanol_fraction = aa_baseline.get_market_value(crude_methanol)*get_system_GWP()/get_total_product_market_value()          
sum_economic = economic_C5_C9_fraction+economic_pa_fraction+economic_glycerol_fraction+economic_fa_fraction+economic_methanol_fraction 

def get_economic_based_AA_GWP():
    def get_total_product_market_value():
        return(sum([aa_baseline.get_market_value(azelaic_acid),
                    aa_baseline.get_market_value(recovered_C5_to_C9_MCA_fraction),
                    aa_baseline.get_market_value(pelargonic_acid_rich_fraction),
                    aa_baseline.get_market_value(crude_glycerol),                                
                    aa_baseline.get_market_value(fatty_acid_blend),    
                    aa_baseline.get_market_value(crude_methanol),
                    ]))

    #Economic allocation
    economic_C5_C9_fraction = aa_baseline.get_market_value(recovered_C5_to_C9_MCA_fraction)*get_system_GWP()/get_total_product_market_value()
    economic_pa_fraction = aa_baseline.get_market_value(pelargonic_acid_rich_fraction)*get_system_GWP()/get_total_product_market_value()
    economic_glycerol_fraction = aa_baseline.get_market_value(crude_glycerol)*get_system_GWP()/get_total_product_market_value() 
    economic_fa_fraction = aa_baseline.get_market_value(fatty_acid_blend)*get_system_GWP()/get_total_product_market_value()          
    economic_methanol_fraction = aa_baseline.get_market_value(crude_methanol)*get_system_GWP()/get_total_product_market_value()          
    sum_economic = economic_C5_C9_fraction+economic_pa_fraction+economic_glycerol_fraction+economic_fa_fraction+economic_methanol_fraction 
#this function estimates GWP via economic allocation per kg of azelaic acid produced
    economic_aa_fraction = aa_baseline.get_market_value(azelaic_acid)*get_system_GWP()/get_total_product_market_value()
    return economic_aa_fraction #Kg CO2 per kg of AA


def get_total_product_mass():
    return(sum([aa_baseline.get_mass_flow(azelaic_acid),
                    aa_baseline.get_mass_flow(recovered_C5_to_C9_MCA_fraction),
                    aa_baseline.get_mass_flow(pelargonic_acid_rich_fraction),
                    aa_baseline.get_mass_flow(crude_glycerol),       
                    aa_baseline.get_mass_flow(fatty_acid_blend),
                    aa_baseline.get_mass_flow(crude_methanol)
                    ]))
        
#Mass allocation
mass_C5_C9_fraction = lambda: aa_baseline.get_mass_flow(recovered_C5_to_C9_MCA_fraction)*get_system_GWP()/get_total_product_mass()
mass_pa_fraction = lambda: aa_baseline.get_mass_flow(pelargonic_acid_rich_fraction)*get_system_GWP()/get_total_product_mass()
mass_glycerol_fraction = lambda: aa_baseline.get_mass_flow(crude_glycerol)*get_system_GWP()/get_total_product_mass() 
mass_fa_fraction = lambda: aa_baseline.get_mass_flow(fatty_acid_blend)*get_system_GWP()/get_total_product_mass()
mass_methanol_fraction = lambda: aa_baseline.get_mass_flow(crude_methanol)*get_system_GWP()/get_total_product_mass() 
#this function estimates GWP via mass allocation per kg of azelaic acid produced
mass_aa_fraction = lambda: aa_baseline.get_mass_flow(azelaic_acid)*get_system_GWP()/get_total_product_mass()
    
#GWP related functions
def sum_displacement():
    # Calculate GWP credits from co-products (per kg azelaic acid basis)
    PA = aa_baseline.get_material_impact(F.pelargonic_acid_rich_fraction, 'GWP100') / aa_baseline.get_mass_flow(azelaic_acid)
    FA = aa_baseline.get_material_impact(F.fatty_acid_blend, 'GWP100') / aa_baseline.get_mass_flow(azelaic_acid)    
    MCA = aa_baseline.get_material_impact(F.recovered_C5_to_C9_MCA_fraction, 'GWP100') / aa_baseline.get_mass_flow(azelaic_acid)
    CG = aa_baseline.get_material_impact(F.crude_glycerol, 'GWP100') / aa_baseline.get_mass_flow(azelaic_acid)
    M = aa_baseline.get_material_impact(F.crude_methanol, 'GWP100') / aa_baseline.get_mass_flow(azelaic_acid)
    sum_displacement = PA + FA + MCA + CG + M
    return sum_displacement

#this function estimates GWP via displacement per kg of azelaic acid produced
def get_net_GWP_displacement():
    return get_system_GWP() - sum_displacement()

#%%
def disp_feedstock(): return get_feedstock_GWP() * 100 / abs(get_system_GWP())
def disp_materials(): return get_other_materials_impact() * 100 / abs(get_system_GWP())
def disp_ng(): return get_ng_GWP() * 100 / abs(get_system_GWP())
def disp_electricity(): return net_electricity_purchased_GWP() * 100 / abs(get_system_GWP())
def disp_direct(): return get_total_direct_emissions_GWP() * 100 / abs(get_system_GWP())

def disp_other_co_products():
    FA = aa_baseline.get_material_impact(fatty_acid_blend, 'GWP100') / aa_baseline.get_mass_flow(azelaic_acid)
    MCA = aa_baseline.get_material_impact(recovered_C5_to_C9_MCA_fraction, 'GWP100') / aa_baseline.get_mass_flow(azelaic_acid)
    CG = aa_baseline.get_material_impact(crude_glycerol, 'GWP100') /aa_baseline.get_mass_flow(azelaic_acid)
    M = aa_baseline.get_material_impact(crude_methanol, 'GWP100') /aa_baseline.get_mass_flow(azelaic_acid)
    return abs(FA + MCA + CG + M) * 100 / abs(get_system_GWP())

def disp_PA():
    PA = aa_baseline.get_material_impact(pelargonic_acid_rich_fraction, 'GWP100') /aa_baseline.get_mass_flow(azelaic_acid)
    return abs(PA) * 100 / abs(get_system_GWP())

#mass allocation fractions based on stacked bar plot
def mass_other_co_products():
    mass_C5_C9_fraction = aa_baseline.get_mass_flow(recovered_C5_to_C9_MCA_fraction)*get_system_GWP()/get_total_product_mass()
    mass_glycerol_fraction = aa_baseline.get_mass_flow(crude_glycerol)*get_system_GWP()/get_total_product_mass() 
    mass_fa_fraction = aa_baseline.get_mass_flow(fatty_acid_blend)*get_system_GWP()/get_total_product_mass()
    mass_methanol_fraction = aa_baseline.get_mass_flow(crude_methanol)*get_system_GWP()/get_total_product_mass()
    sum_mass_others =  mass_C5_C9_fraction+mass_glycerol_fraction+mass_fa_fraction+mass_methanol_fraction
    return  abs(sum_mass_others)*100/abs(get_system_GWP())

def mass_PA():
    mass_pa_fraction = aa_baseline.get_mass_flow(pelargonic_acid_rich_fraction)*get_system_GWP()/get_total_product_mass()
    return abs(mass_pa_fraction) *100/abs(get_system_GWP())

#economic allocation fractions based on stacked bar plot

def econ_other_co_products():
    economic_C5_C9_fraction = aa_baseline.get_market_value(recovered_C5_to_C9_MCA_fraction)*get_system_GWP()/get_total_product_market_value()
    economic_glycerol_fraction = aa_baseline.get_market_value(crude_glycerol)*get_system_GWP()/get_total_product_market_value() 
    economic_fa_fraction = aa_baseline.get_market_value(fatty_acid_blend)*get_system_GWP()/get_total_product_market_value()          
    economic_methanol_fraction = aa_baseline.get_market_value(crude_methanol)*get_system_GWP()/get_total_product_market_value()          
    sum_economic_others = economic_C5_C9_fraction+economic_glycerol_fraction+economic_fa_fraction+economic_methanol_fraction 
    return abs(sum_economic_others)*100/abs(get_system_GWP())
def econ_PA():
    economic_pa_fraction = aa_baseline.get_market_value(pelargonic_acid_rich_fraction)*get_system_GWP()/get_total_product_market_value()
    return abs(economic_pa_fraction)*100/abs(get_system_GWP())

#%% Keeping track of each section
#all fractions
def transesterification_IEC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Installed equipment cost']['Transesterification']

def dih_IEC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Installed equipment cost']['Dihydroxylation'] 

def ox_cl_IEC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Installed equipment cost']['Oxidative cleavage'] 

def cat_recovery_IEC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Installed equipment cost']['Catalyst recovery'] 

def pa_recovery_IEC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Installed equipment cost']['Pelargonic acid and C5-C8 fraction recovery']

def aa_recovery_IEC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Installed equipment cost']['Azelaic acid and heavy tails recovery'] 

def boiler_IEC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Installed equipment cost']['Boilerturbogenerator']

def wastewater_IEC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Installed equipment cost']['Wastewater treatment']
def other_facility_IEC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Installed equipment cost']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)']

#All MC's
def transesterification_MC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Material cost']['Transesterification']

def dih_MC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Material cost']['Dihydroxylation'] 

def ox_cl_MC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Material cost']['Oxidative cleavage'] 

def cat_recovery_MC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Material cost']['Catalyst recovery'] 

def aa_recovery_MC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Material cost']['Azelaic acid and heavy tails recovery'] 

def boiler_MC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Material cost']['Boilerturbogenerator']
def wastewater_MC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Material cost']['Wastewater treatment']
def other_facility_MC():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Material cost']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)']

#All Heating demand
def transesterification_heating_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Heating duty']['Transesterification']

def dih_heating_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Heating duty']['Dihydroxylation'] 

def ox_cl_heating_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Heating duty']['Oxidative cleavage'] 

def cat_recovery_heating_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Heating duty']['Catalyst recovery'] 

def pa_recovery_heating_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Heating duty']['Pelargonic acid and C5-C8 fraction recovery']

def aa_recovery_heating_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Heating duty']['Azelaic acid and heavy tails recovery'] 

#All cooling demand
def transesterification_cooling_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Cooling duty']['Transesterification']

def dih_cooling_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Cooling duty']['Dihydroxylation'] 

def ox_cl_cooling_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Cooling duty']['Oxidative cleavage'] 

def cat_recovery_cooling_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Cooling duty']['Catalyst recovery'] 

def pa_recovery_cooling_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Cooling duty']['Pelargonic acid and C5-C8 fraction recovery']

def aa_recovery_cooling_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Cooling duty']['Azelaic acid and heavy tails recovery'] 

def boiler_cooling_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Cooling duty']['Boilerturbogenerator']
def other_facility_cooling_demand():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Cooling duty']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)']

#All electricity consumption
def transesterification_electricity_consumption():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Electricity consumption']['Transesterification']

def dih_electricity_consumption():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Electricity consumption']['Dihydroxylation'] 

def ox_cl_electricity_consumption():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Electricity consumption']['Oxidative cleavage'] 

def cat_recovery_electricity_consumption():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Electricity consumption']['Catalyst recovery'] 


def pa_recovery_electricity_consumption():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Electricity consumption']['Pelargonic acid and C5-C8 fraction recovery']

def aa_recovery_electricity_consumption():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Electricity consumption']['Azelaic acid and heavy tails recovery'] 

def boiler_electricity_consumption():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Electricity consumption']['Boilerturbogenerator']
def wastewater_electricity_consumption():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Electricity consumption']['Wastewater treatment']
def other_facility_electricity_consumption():
    df_unit_groups = df_related_fraction()
    return df_unit_groups['Electricity consumption']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)']

#Additional specific functions
def natural_gas_material_cost():
    total_annual_NG_cost = aa_baseline.get_mass_flow(F.natural_gas)*F.natural_gas.price
    return total_annual_NG_cost*100/aa_baseline.material_cost
def feedstock_material_cost():
    feedstock_cost = aa_baseline.get_mass_flow(F.crude_vegetable_oil)*F.crude_vegetable_oil.price
    return feedstock_cost*100/aa_baseline.material_cost

def hp_material_cost():
    hp_cost = aa_baseline.get_mass_flow(F.stream.fresh_HP)*F.stream.fresh_HP.price
    return hp_cost*100/aa_baseline.material_cost

def hydrolysis_aa_IEC_frac():
    sum_hydrolysis_IEC = []
    for i in [F.T602, F.S606, F.R605, F.T606, F.R606, F.T607, F.R607, F.T601, 
    F.P601, F.M608, F.S607, F.T603]:
        sum_hydrolysis_IEC.append(i.installed_cost)
    return sum(sum_hydrolysis_IEC)*100/aa_baseline.installed_cost

def HD_azelaic_acid_purification():
    sum_1 = sum([j.duty for j in F.HX602.heat_utilities if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9
    sum_2 = sum([j.duty for j in F.HX603.heat_utilities if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9
    sum_3 = sum([j.duty for j in F.H608.heat_utilities if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9    
    sum_5 = sum([j.duty for j in F.H612.heat_utilities if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9   
    sum_6 = sum([j.duty for j in F.MMS601.heat_utilities if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9
    sum_7 = sum([j.duty for j in F.E601.heat_utilities if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9  
    sum_8 = sum([j.duty for j in F.HX606.heat_utilities if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9        
    fraction = (sum_1+sum_2+sum_3+sum_5+sum_6+sum_7+sum_8)*100/get_system_heating_demand()
    return fraction
def CD_azelaic_acid_purification():
    sum_1 = sum([j.duty for j in F.HX602.heat_utilities if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9
    sum_2 = sum([j.duty for j in F.HX603.heat_utilities if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9
    sum_3 = sum([j.duty for j in F.H608.heat_utilities if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9    
    sum_5 = sum([j.duty for j in F.H612.heat_utilities if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9   
    sum_6 = sum([j.duty for j in F.MMS601.heat_utilities if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9
    sum_7 = sum([j.duty for j in F.E601.heat_utilities if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9  
    sum_8 = sum([j.duty for j in F.HX606.heat_utilities if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9        
    fraction = -(sum_1+sum_2+sum_3+sum_5+sum_6+sum_7+sum_8)*100/get_system_cooling_demand()
    return fraction

def CD_during_hydrolysis_all():
    # List of hydrolysis units
    hydrolysis_units = [
        F.T602, F.S606, F.R605, F.T606, F.R606, F.T607,
        F.R607, F.T601, F.P601, F.M608, F.S607, F.T603
    ]

    # Sum up cooling duties (negative values) for each unit
    total_CD = 0
    for unit in hydrolysis_units:
        unit_CD = sum([j.duty for j in unit.heat_utilities if j.flow > 0 and j.duty < 0])
        total_CD += unit_CD

    # Convert to GJ/year
    total_CD_GJ_per_year = total_CD * aa_baseline.operating_hours / 1e9

    # Calculate as fraction of total plant cooling demand
    fraction = -(total_CD_GJ_per_year / get_system_cooling_demand()) * 100
    return fraction

def HD_during_hydrolysis_all():
    # List of hydrolysis units
    hydrolysis_units = [
        F.T602, F.S606, F.R605, F.T606, F.R606, F.T607,
        F.R607, F.T601, F.P601, F.M608, F.S607, F.T603
    ]

    # Sum up cooling duties (negative values) for each unit
    total_CD = 0
    for unit in hydrolysis_units:
        unit_CD = sum([j.duty for j in unit.heat_utilities if j.flow > 0 and j.duty > 0])
        total_CD += unit_CD

    # Convert to GJ/year
    total_CD_GJ_per_year = total_CD * aa_baseline.operating_hours / 1e9

    # Calculate as fraction of total plant cooling demand
    fraction = -(total_CD_GJ_per_year / get_system_cooling_demand()) * 100
    return fraction

def HD_mma_azelaic_acid_sep():
    sum_1 = sum([j.duty for j in F.D608.heat_utilities if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9
    return sum_1*100/get_system_heating_demand()

def CD_mma_azelaic_acid_sep():
    sum_1 = sum([j.duty for j in F.D608.heat_utilities if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9
    return -sum_1*100/get_system_cooling_demand()

def D502_elec_fraction():
    return  F.D502.power_utility.consumption*100/aa_baseline.power_utility.consumption

def F608_elec_fraction():
    return F.F608.power_utility.consumption*100/aa_baseline.power_utility.consumption

def CW701_elec_fraction():
    return F.CW701.power_utility.consumption*100/aa_baseline.power_utility.consumption

##############################################################################################3
all_metrics = [
        # Economic performance and TEA metrics
    Metric('MPSP', get_MPSP, '$/kg', 'TEA'),
    Metric('Maximum feedstock purchase price', get_MFPP, '$/kg', 'TEA'),
    Metric('Total product yield', get_total_yield, '10^6 kg/yr', 'TEA'),
    Metric('Product purity', get_purity, '%', 'TEA'),
    Metric('Material cost', get_material_cost, 'MM$/yr', 'TEA'),
    Metric('Utility cost', get_utility_cost, 'MM$/yr', 'TEA'),
    Metric('Annual operating cost', get_AOC, 'MM$/yr', 'TEA'),
    Metric('Total production cost', get_total_production_cost, 'MM$/yr', 'TEA'),
    Metric('Net present value', solve_NPV, 'MM$', 'TEA'),
    Metric('Internal rate of return (IRR)', get_IRR, '%', 'TEA'),
    Metric('Return on investment (ROI)', get_ROI, '$/yr', 'TEA'),
    
    # Capital investment metrics
    Metric('Total capital investment', get_TCI, 'MM$', 'TEA'),
    Metric('Fixed capital investment (FCI)', get_FCI, 'MM$', 'TEA'),
    Metric('Direct permanent investment (DPI)', get_DPI, 'MM$', 'TEA'),
    Metric('Total depreciable capital (TDC)', get_TDC, 'MM$', 'TEA'),
    
    # Resource consumption metrics
    Metric('Quantity of system makeup water', get_amt_of_system_makeup_water, 'kg/h', 'TEA'),
    Metric('Total electricity consumption', electricity_consumption, 'kWh/yr', 'TEA'),
    Metric('Total heating demand', get_system_heating_demand, 'MMJ/yr', 'Heating demand'),
    Metric('Total cooling demand', get_system_cooling_demand, 'MMJ/yr', 'Cooling demand'),
    
    # LCA metrics – mass allocation
    Metric('Mass AA LCA', mass_aa_fraction, 'kg CO2-eq/kg', 'LCA'),
    # LCA metrics – economic allocation and displacement
    Metric('Economic allocation LCA', get_economic_based_AA_GWP, 'kg CO2-eq/kg', 'LCA'),
    Metric('Net CI displacement LCA', get_net_GWP_displacement, 'kg CO2-eq/kg', 'LCA'),  
    Metric('Total emissions', get_system_GWP, 'kg CO2-eq/kg', 'LCA'),  
    
    Metric('Feedstock GWP – Displacement', disp_feedstock, '%', 'LCA'),
    Metric('Other materials GWP – Displacement', disp_materials, '%', 'LCA'),
    Metric('Natural gas GWP – Displacement', disp_ng, '%', 'LCA'),
    Metric('Electricity GWP – Displacement', disp_electricity, '%', 'LCA'),
    Metric('Direct emissions GWP – Displacement', disp_direct, '%', 'LCA'),
    
    Metric('Other coproducts GWP – Displacement', disp_other_co_products, '%', 'LCA'),
    Metric('Pelargonic acid GWP – Displacement', disp_PA, '%', 'LCA'),
    
   
    Metric('Other coproducts GWP – Mass allocation', mass_other_co_products, '%', 'LCA'),
    Metric('Pelargonic acid GWP – Mass allocation', mass_PA, '%', 'LCA'),
    
    Metric('Other coproducts GWP – Economic allocation', econ_other_co_products, '%', 'LCA'),
    Metric('Pelargonic acid GWP – Economic allocation', econ_PA, '%', 'LCA'),
    
    
    # Installed equipment cost (%)
    Metric('Installed equipment cost – Transesterification', transesterification_IEC, '%', 'TEA'),
    Metric('Installed equipment cost – Dihydroxylation', dih_IEC, '%', 'TEA'),
    Metric('Installed equipment cost – Oxidative cleavage', ox_cl_IEC, '%', 'TEA'),
    Metric('Installed equipment cost – Catalyst recovery', cat_recovery_IEC, '%', 'TEA'),
    Metric('Installed equipment cost – PA + C5–C8 recovery', pa_recovery_IEC, '%', 'TEA'),
    Metric('Installed equipment cost – AA + heavy tails recovery', aa_recovery_IEC, '%', 'TEA'),
    Metric('Installed equipment cost – Boilerturbogenerator', boiler_IEC, '%', 'TEA'),
    Metric('Installed equipment cost – Wastewater treatment', wastewater_IEC, '%', 'TEA'),
    Metric('Installed equipment cost – Other facilities', other_facility_IEC, '%', 'TEA'),
    
    # Material cost (%)
    Metric('Material cost – Transesterification', transesterification_MC, '%', 'TEA'),
    Metric('Material cost – Dihydroxylation', dih_MC, '%', 'TEA'),
    Metric('Material cost – Oxidative cleavage', ox_cl_MC, '%', 'TEA'),
    Metric('Material cost – Catalyst recovery', cat_recovery_MC, '%', 'TEA'),
    Metric('Material cost – AA + heavy tails recovery', aa_recovery_MC, '%', 'TEA'),
    Metric('Material cost – Boilerturbogenerator', boiler_MC, '%', 'TEA'),
    Metric('Material cost – Wastewater treatment', wastewater_MC, '%', 'TEA'),
    Metric('Material cost – Other facilities', other_facility_MC, '%', 'TEA'),
    
    # Heating demand (%)
    Metric('Heating duty – Transesterification', transesterification_heating_demand, '%', 'Utility'),
    Metric('Heating duty – Dihydroxylation', dih_heating_demand, '%', 'Utility'),
    Metric('Heating duty – Oxidative cleavage', ox_cl_heating_demand, '%', 'Utility'),
    Metric('Heating duty – Catalyst recovery', cat_recovery_heating_demand, '%', 'Utility'),
    Metric('Heating duty – PA + C5–C8 recovery', pa_recovery_heating_demand, '%', 'Utility'),
    Metric('Heating duty – AA + heavy tails recovery', aa_recovery_heating_demand, '%', 'Utility'),
    
    # Cooling demand (%)
    Metric('Cooling duty – Transesterification', transesterification_cooling_demand, '%', 'Utility'),
    Metric('Cooling duty – Dihydroxylation', dih_cooling_demand, '%', 'Utility'),
    Metric('Cooling duty – Oxidative cleavage', ox_cl_cooling_demand, '%', 'Utility'),
    Metric('Cooling duty – Catalyst recovery', cat_recovery_cooling_demand, '%', 'Utility'),
    Metric('Cooling duty – PA + C5–C8 recovery', pa_recovery_cooling_demand, '%', 'Utility'),
    Metric('Cooling duty – AA + heavy tails recovery', aa_recovery_cooling_demand, '%', 'Utility'),
    Metric('Cooling duty – Boilerturbogenerator', boiler_cooling_demand, '%', 'Utility'),
    Metric('Cooling duty – Other facilities', other_facility_cooling_demand, '%', 'Utility'),
    
    # Electricity consumption (%)
    Metric('Electricity – Transesterification', transesterification_electricity_consumption, '%', 'Utility'),
    Metric('Electricity – Dihydroxylation', dih_electricity_consumption, '%', 'Utility'),
    Metric('Electricity – Oxidative cleavage', ox_cl_electricity_consumption, '%', 'Utility'),
    Metric('Electricity – Catalyst recovery', cat_recovery_electricity_consumption, '%', 'Utility'),
    Metric('Electricity – PA + C5–C8 recovery', pa_recovery_electricity_consumption, '%', 'Utility'),
    Metric('Electricity – AA + heavy tails recovery', aa_recovery_electricity_consumption, '%', 'Utility'),
    Metric('Electricity – Boilerturbogenerator', boiler_electricity_consumption, '%', 'Utility'),
    Metric('Electricity – Wastewater treatment', wastewater_electricity_consumption, '%', 'Utility'),
    Metric('Electricity – Other facilities', other_facility_electricity_consumption, '%', 'Utility'),

    Metric('Material cost – Natural gas', natural_gas_material_cost, '%', 'TEA'),
    Metric('Material cost – Feedstock', feedstock_material_cost, '%', 'TEA'),
    Metric('Material cost – Hydrogen peroxide', hp_material_cost, '%', 'TEA'),
    
    Metric('IEC – Hydrolysis section', hydrolysis_aa_IEC_frac, '%', 'TEA'),
    Metric('HD – Hydrolysis section', HD_during_hydrolysis_all, '%', 'Utility'),
    Metric('CD – Hydrolysis section', CD_during_hydrolysis_all, '%', 'Utility'),
    
    Metric('HD – Azelaic acid purification', HD_azelaic_acid_purification, '%', 'Utility'),
    Metric('CD – Azelaic acid purification', CD_azelaic_acid_purification, '%', 'Utility'),
    
    Metric('HD – MMA/AA separation', HD_mma_azelaic_acid_sep, '%', 'Utility'),
    Metric('CD – MMA/AA separation', CD_mma_azelaic_acid_sep, '%', 'Utility'),
    
    Metric('Electricity – D502 fraction', D502_elec_fraction, '%', 'Utility'),
    Metric('Electricity – F608 fraction', F608_elec_fraction, '%', 'Utility'),
    Metric('Electricity – CW701 fraction', CW701_elec_fraction, '%', 'Utility'),
    ]

def _slice_rows_from_excel(xlsx_path, sheet_name, row_indices, header_rows=3):
    """
    Return a NumPy array of samples for the given row(s).
    
    Parameters
    ----------
    row_indices : int or list[int]
        - If int → run a single Excel row (e.g. 77).
        - If list[int] → run multiple Excel rows (e.g. [4, 77, 100]).
    
    header_rows : int
        Number of non-data rows at the top (default=3 → data starts at Excel row 4).
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=None)
    
    # Normalize row_indices into a list
    if isinstance(row_indices, int):
        row_indices = [row_indices]
    
    # Convert Excel row numbers to pandas iloc indices
    adjusted = [r - 1 for r in row_indices]  # Excel rows are 1-based
    return df.iloc[adjusted].to_numpy()


def conduct_uncertainity_analysis(system = aa_baseline, #name of the system
                                  metrics = all_metrics, #list of Metric objects
                                  number_of_runs = 100, #number of 
                                  exception_hook = 'warn',#raise/warn
                                  rule = 'L', #for For Latin-Hypercube sampling
                                  notify_runs = 10, #runs after which you are notified 
                                  indicator = 'GWP100',
                                  feedstock_type = 'HoySoy_oil', #can also provide 'HoSun_oil
                                  xlsx_path = r"C:/Users/lavan/gglk_code/Bioindustrial-Park/biorefineries/oleochemicals/model_table_10000_w_aa_price_frac_lca.xlsx",
                                  sheet_name = "Sheet1",
                                  batch_id = 1,
                                  batch_size = 1000,
                                  row_indices = None,
                                  header_rows = 3,          # <-- set to 3 to make batch 1 = Excel rows 4–1003 (inclusive)
                                  outdir = "C:/Users/lavan/gglk_code/Bioindustrial-Park/biorefineries/oleochemicals/batch_wise_uncertainty_analysis"
                                                              
                                 ):
    model = Model(aa_baseline,metrics = all_metrics,
                  exception_hook= exception_hook,
                  retry_evaluation = True,                  
                  )
    ## price distributions    
    @model.parameter(name='Crude vegetable oil cost',
                      element=crude_vegetable_oil, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 2.13,
                      distribution= isolated_para_dists['Crude oil cost'][feedstock_type])
    def set_feedstock_cost(feedstock_cost):
            crude_vegetable_oil.price = feedstock_cost
            
            
    @model.parameter(name='Methanol cost',
                      element=F.methanol, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 0.79,
                      distribution= isolated_para_dists['Methanol'])
    def set_methanol_cost(methanol_cost):
            F.methanol.price = methanol_cost
            
    @model.parameter(name='Catalyst biodiesel cost',
                      element=F.catalyst, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 1.53,
                      distribution= isolated_para_dists['Catalyst_for_trans'])
    def set_trans_catalyst_cost(trans_catalyst_cost):
            F.catalyst.price = trans_catalyst_cost                 
        
    @model.parameter(name='sodium hydroxide biodiesel cost',
                      element=F.NaOH, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 0.91,
                      distribution= isolated_para_dists['NaOH'])
    def set_naoh_cost(naoh_cost):
            F.NaOH.price = naoh_cost

    @model.parameter(name='sodium hydroxide catalyst recovery cost',
                      element=F.sodium_hydroxide_stream, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 0.91,
                      distribution= isolated_para_dists['NaOH for cat'])
    def set_naoh_cat_cost(naoh_cost):
            F.sodium_hydroxide_stream.price = naoh_cost
            
    @model.parameter(name='HCl biodiesel cost',
                      element=F.HCl, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 0.13,
                      distribution= isolated_para_dists['Conc HCl'])
    def set_hcl_cost(hcl_cost):
            F.HCl.price = hcl_cost

    @model.parameter(name='HCl catalyst recovery cost',
                      element=F.conc_hydrochloric_acid, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 0.13,
                      distribution= isolated_para_dists['Conc HCl cat'])
    def set_hcl_cat_cost(hcl_cost):
            F.conc_hydrochloric_acid.price = hcl_cost    
            
    @model.parameter(name='HCl resin wash cost',
                      element=F.conc_HCl_for_resin_wash, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 0.13,
                      distribution= isolated_para_dists['Conc HCl resin'])
    def set_hcl_resin_cost(hcl_cost):
            F.conc_HCl_for_resin_wash.price = hcl_cost    
            
    @model.parameter(name='Tungstic acid cost',
                      element=F.fresh_tungsten_catalyst, 
                      kind='isolated', 
                      units='USD/kg',
                      baseline = 35.25,
                      distribution= isolated_para_dists['Tungstic acid cost'])
    def set_tungstencat_cost(tungstencat_cost):
        F.fresh_tungsten_catalyst.price = tungstencat_cost

    @model.parameter(name='Hydrogen peroxide cost',
                      element=F.fresh_HP, 
                      kind='isolated',
                      baseline = 1.5,
                      units='USD/kg',
                      distribution=isolated_para_dists['Hydrogen peroxide cost'])
    def set_fresh_HP_cost(fresh_HP_cost):
        fresh_HP.price = fresh_HP_cost

    @model.parameter(name='Cobalt acetate catalyst cost',
                      element=fresh_cobalt_catalyst, 
                      kind='isolated', 
                      units='USD/kg',
                      baseline=20.76,
                      distribution= isolated_para_dists['Cobalt acetate cost'])
    def set_cobaltcat_cost(fresh_cobalt_catalyst_cost):
        fresh_cobalt_catalyst.price = fresh_cobalt_catalyst_cost
        
    @model.parameter(name='Calcium chloride price',
                  element=F.calcium_chloride_for_cat_sep,
                  kind='isolated',
                  units='USD/kg',
                  baseline = 0.82,
                  distribution=isolated_para_dists['Calcium chloride cost'])
    def set_calcium_chloride_cost(calcium_chloride_cost):
        F.calcium_chloride_for_cat_sep.price = calcium_chloride_cost
        
    @model.parameter(name='Hydrolysis resin cost',
                      element=F.polystyrene_based_catalyst,
                      kind='isolated',
                      units='USD/kg',
                      baseline = 106,
                      distribution=isolated_para_dists['Hydrolysis resin cost'])
    def set_polystyrene_based_catalyst_cost(resin_cost):
        F.polystyrene_based_catalyst.price = resin_cost
    
    @model.parameter(name='Heptane solvent cost',
                      element=fresh_solvent, 
                      kind='isolated',
                      units='USD/kg',
                      baseline = 0.84,
                      distribution=isolated_para_dists['Heptane solvent cost'])
    def set_solvent_cost(fresh_solvent_cost):
        fresh_solvent.price = fresh_solvent_cost    
                        
#Product prices    
    @model.parameter(name='Crude glycerol price',
                       element=F.crude_glycerol,
                       kind='isolated',
                       baseline = 0.16,
                       units='USD/kg',
                       distribution=isolated_para_dists['Crude glycerol price'])
    def set_crude_glycerol_price(crude_glycerol_price):
         F.crude_glycerol.price = crude_glycerol_price
         
    @model.parameter(name='C5_C9 price',
                       element=F.recovered_C5_to_C9_MCA_fraction,
                       kind='isolated',
                       baseline = 2.37,
                       units='USD/kg',
                       distribution=isolated_para_dists['C5_C9 fraction price'])
    def set_recovered_C5_to_C9_MCA_fraction_price(recovered_C5_to_C9_MCA_fraction_price):
         F.recovered_C5_to_C9_MCA_fraction.price = recovered_C5_to_C9_MCA_fraction_price
                    
    @model.parameter(name='Pelargonic acid price',
                      element=pelargonic_acid_rich_fraction,
                      kind='isolated',
                      baseline=5.37,
                      units='USD/kg',
                      distribution=isolated_para_dists['Pelargonic acid price'])
    def set_pelargonic_acid_price(pelargonic_acid_price):
        pelargonic_acid_rich_fraction.price = pelargonic_acid_price 
            
    @model.parameter(name='Methanol price',
                      element=F.crude_methanol, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 0.79,
                      distribution= isolated_para_dists['Methanol product'])
    def set_methanol_product_price(methanol_price):
            F.crude_methanol.price = methanol_price            

    @model.parameter(name='Fatty acid blend price',
                      element= F.fatty_acid_blend,
                      kind='isolated',
                      baseline= 1.91,
                      units='USD/kg',
                      distribution= isolated_para_dists['fatty_acid_blend_price'])
    def set_fatty_acid_blend_price(fa_blend_price):
        F.fatty_acid_blend.price = fa_blend_price

    @model.parameter(name='Azelaic acid price',
                      element= F.azelaic_acid_product_stream,
                      kind='isolated',
                      baseline=11.02,
                      units='USD/kg',
                      distribution= isolated_para_dists['Azelaic acid price'])
    def set_aa_price(aa_price):
        F.azelaic_acid_product_stream.price = aa_price
   
        
    @model.parameter(name='Natural gas cost',
                      element=F.BT801,
                      baseline = 0.235,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Natural gas cost'])
    def set_natural_gas_cost(natural_gas_cost):
        F.BT801.natural_gas_cost = natural_gas_cost    
    
    @model.parameter(name='Electricty cost',
                      element=bst.settings.electricity_price,
                      baseline =  0.09,
                      kind='isolated',
                      units='USD/kWh',
                      distribution=isolated_para_dists['Electricity cost'])
    def set_electricity_cost(electricity_cost):
        bst.settings.electricity_price = electricity_cost
    
    @model.parameter(name='Ash disposal cost',
                      element=F.ash_disposal,
                      baseline = -0.032,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Ash_disposal_cost'])
    def set_ash_disposal_cost(ash_disposal_cost):
        F.ash_disposal.price = ash_disposal_cost

    @model.parameter(name='Cooling tower chemicals cost',
                      element=F.cooling_tower_chemicals,
                      baseline = 2.41,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Cooling_tower_chemicals'])
    def set_cooling_tower_chems_cost(cooling_tower_chemicals_cost):
        F.cooling_tower_chemicals.price = cooling_tower_chemicals_cost
        
        
    @model.parameter(name='Lime used in boiler cost',
                      element=F.lime_boiler,
                      baseline = 0.11,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Lime_boiler'])
    def set_lime_cost(lime_cost):
        F.lime_boiler.price = lime_cost
        
    #environmental_facs_dist
    @model.parameter(name='crude oil GWP',
                      element=crude_vegetable_oil,
                      kind='isolated',
                      baseline = 0.48,
                      distribution=environmental_facs_dist['Oil_GWP'][feedstock_type]
                      )
    def crude_vegetable_oil_gwp(crude_vegetable_oil_gwp):
        crude_vegetable_oil.characterization_factors[indicator] = crude_vegetable_oil_gwp

         
    @model.parameter(name='Catalyst biodiesel GWP',
                     element=F.catalyst, 
                     kind='isolated',
                     baseline= 1.32,
                     distribution= environmental_facs_dist['Catalyst_for_trans'])
    def set_trans_catalyst_gwp(trans_catalyst_gwp):
           F.catalyst.characterization_factors[indicator] = trans_catalyst_gwp                 
       
    @model.parameter(name='sodium hydroxide biodiesel GWP',
                     element=F.NaOH, 
                     kind='isolated',
                     baseline= 2.05,
                     distribution= environmental_facs_dist['Sodium_hydroxide'])
    def set_naoh_gwp(naoh_gwp):
           F.NaOH.characterization_factors[indicator] = naoh_gwp

    @model.parameter(name='sodium hydroxide catalyst recovery gwp',
                     element=F.sodium_hydroxide_stream, 
                     kind='isolated',
                     baseline= 2.05,
                     distribution= environmental_facs_dist['Sodium_hydroxide cat'])
    def set_naoh_cat_gwp_cat(naoh_gwp):
           F.sodium_hydroxide_stream.characterization_factors[indicator] = naoh_gwp

    @model.parameter(name='HCl biodiesel gwp',
                      element=F.HCl, 
                      kind='isolated',
                      baseline= 0.67,
                      distribution= environmental_facs_dist['HCl'])
    def set_hcl_biodiesel(hcl_gwp):
            F.HCl.characterization_factors[indicator] = hcl_gwp

    @model.parameter(name='HCl catalyst recovery gwp',
                      element=F.conc_hydrochloric_acid, 
                      kind='isolated',
                      baseline= 0.67,
                      distribution= environmental_facs_dist['HCl cat'])
    def set_hcl_cat_gwp(hcl_gwp):
            F.conc_hydrochloric_acid.characterization_factors[indicator]= hcl_gwp
            
    @model.parameter(name='HCl resin wash gwp',
                      element=F.conc_HCl_for_resin_wash, 
                      kind='isolated',
                      baseline= 0.67,
                      distribution= environmental_facs_dist['HCl resin'])
    def set_hcl_resin_gwp(hcl_gwp):
            F.conc_HCl_for_resin_wash.characterization_factors[indicator] = hcl_gwp
            
    @model.parameter(name='Crude glycerol gwp',
                       element=F.crude_glycerol,
                       kind='isolated',
                       baseline = 4.29,
                       distribution= environmental_facs_dist['Glycerol'])
    def set_crude_glycerol_gwp(crude_glycerol_gwp):            
        F.crude_glycerol.characterization_factors[indicator]= crude_glycerol_gwp
                       
    @model.parameter(name='tungstic acid GWP',
                      element=fresh_tungsten_catalyst,
                      kind='isolated',
                      baseline = 68.5,
                      distribution=environmental_facs_dist['Tungstic_acid']
                      )
    def tungstic_acid_gwp(tungstic_acid_gwp):
        fresh_tungsten_catalyst.characterization_factors[indicator] = tungstic_acid_gwp 
        
    @model.parameter(name='HP_GWP',
                      element=fresh_HP,
                      kind='isolated',
                      baseline = 0.54,
                      distribution=environmental_facs_dist['Hydrogen peroxide']
                      )
    def HP_GWP(HP_GWP):
        fresh_HP.characterization_factors[indicator]  = HP_GWP               
        
    @model.parameter(name='cobalt acetate GWP',
                      element=fresh_cobalt_catalyst,
                      kind='isolated',
                      baseline = 8.59,
                      distribution=environmental_facs_dist['Cobalt_acetate']
                      )
    def cobalt_acetate_gwp(cobalt_acetate_gwp):
        fresh_cobalt_catalyst.characterization_factors[indicator] = cobalt_acetate_gwp 
    
    @model.parameter(name='Calcium chloride gwp',
                  element=F.calcium_chloride_for_cat_sep,
                  kind='isolated',
                  baseline = 1.61,
                  distribution=environmental_facs_dist['Calcium_chloride'])
    def set_calcium_chloride_gwp(calcium_chloride_gwp):
        F.calcium_chloride_for_cat_sep.characterization_factors[indicator] = calcium_chloride_gwp
    
    @model.parameter(name='pelargonic acid GWP',
                      element=pelargonic_acid_rich_fraction,
                      kind='isolated',
                      baseline = 11.24,
                      distribution=environmental_facs_dist['Pelargonic_acid']
                      )
    def pelargonic_acid_rich_fraction_gwp(pelargonic_acid_rich_fraction_gwp):
        pelargonic_acid_rich_fraction.characterization_factors[indicator]  = pelargonic_acid_rich_fraction_gwp          
    
    @model.parameter(name='C5_C9 fraction GWP',
                      element=recovered_C5_to_C9_MCA_fraction,
                      kind='isolated',
                      baseline = 2.53,
                      distribution=environmental_facs_dist['C5_C9_MCA']
                      )
    def recovered_C5_to_C9_MCA_fraction_gwp(recovered_C5_to_C9_MCA_fraction_gwp):
        recovered_C5_to_C9_MCA_fraction.characterization_factors[indicator]  = recovered_C5_to_C9_MCA_fraction_gwp    
     
    @model.parameter(name='resin GWP',
                      element=resin,
                      kind='isolated',
                      baseline = 3.06,
                      distribution=environmental_facs_dist['Resin']
                      )
    def resin_gwp(resin_gwp):
        resin.characterization_factors[indicator]  = resin_gwp

    @model.parameter(name='Heptane solvent gwp',
                      element=fresh_solvent, 
                      kind='isolated',
                      baseline = 0.8,
                      distribution=environmental_facs_dist['Heptane'])
    def set_solvent_gwp(fresh_solvent_gwp):
        fresh_solvent.characterization_factors[indicator] = fresh_solvent_gwp        
    
    @model.parameter(name='fatty acid blend GWP',
                      element=fatty_acid_blend,
                      kind='isolated',
                      baseline = 0.57,
                      distribution=environmental_facs_dist['fatty_acid_blend']
                      )
    def fatty_acid_blend_gwp(fatty_acid_blend_gwp):
        fatty_acid_blend.characterization_factors[indicator]  = fatty_acid_blend_gwp     
        
    @model.parameter(name='Methanol product gwp',
                      element=F.crude_methanol, 
                      kind='isolated',
                      baseline= 1.14,
                      distribution= environmental_facs_dist['Methanol product'])
    def set_methanol_product_gwp(methanol_gwp):
            F.crude_methanol.characterization_factors[indicator] = methanol_gwp

    @model.parameter(name='Natural gas gwp',
                      element=F.BT801,
                      baseline = 0.38,
                      kind='isolated',
                      distribution=environmental_facs_dist['Natural gas'])
    def set_natural_gas_gwp(natural_gas_gwp):
        F.BT801.natural_gas.characterization_factors[indicator] = natural_gas_gwp
    
    @model.parameter(name='Electricty gwp',
                      element=bst.settings.set_electricity_CF(indicator,0.39,0.39),
                      baseline =  0.39,
                      kind='isolated',
                      distribution=environmental_facs_dist['Electricity'])
    def set_electricity_gwp(electricity_gwp):
        bst.settings.set_electricity_CF(indicator,electricity_gwp, electricity_gwp,
                                        basis='kWhr', units='kg*CO2e')
   
    @model.parameter(name='Cooling tower chemicals gwp',
                      element=F.cooling_tower_chemicals,
                      baseline = 0.56,
                      kind='isolated',
                      distribution=environmental_facs_dist['Cooling_tower_chemicals'])
    def set_cooling_tower_chems_gwp(cooling_tower_chemicals_gwp):
        F.cooling_tower_chemicals.characterization_factors[indicator]  = cooling_tower_chemicals_gwp
        
    @model.parameter(name='Boiler chems cost',
                     element=F.boiler_chems,
                     baseline =2.05,
                     kind='isolated',
                     distribution=environmental_facs_dist['Boiler_chems'])
    def set_boiler_chems_gwp(boiler_chems_gwp):
       F.boiler_chems.characterization_factors[indicator]   = boiler_chems_gwp
       
    @model.parameter(name='Lime used in boiler gwp',
                      element=F.lime_boiler,
                      baseline = 1.28,
                      kind='isolated',
                      distribution=environmental_facs_dist['Lime_boiler'])
    def set_lime_gwp(lime_gwp):
        F.lime_boiler.characterization_factors[indicator] = lime_gwp

    # # coupled parameters
    # # Process related parameters
    # # #Tungstic acid mole ratio
    @model.parameter(name='Tungstic acid moles',
                      element=F.M200,
                      kind='coupled',
                      baseline = 0.0078,
                      distribution= coupled_para_dist['Tungstic acid moles'])
    def set_tungstic_acid_moles(X_tam):
        F.unit.M200.specifications[0].args[0] = X_tam
    
    @model.parameter(name='Cobalt acetate moles',
                      element=F.R300,
                      kind='coupled',
                      baseline= 0.015,
                      distribution= coupled_para_dist['Cobalt acetate moles'])
    def set_cobalt_acetate_moles(X_cam):
        F.unit.R300.specifications[0].args[1] = X_cam
        
    @model.parameter(name='Dihydroxylation reaction conversion',
                      element=F.R200,
                      kind='coupled',
                      baseline = 0.86,
                      distribution=coupled_para_dist['Dihydroxylation reaction conversion'])
    def set_X_dih_conversion(X_dih):
        F.unit.R200.X_dih = X_dih
          
    ## Oxidative cleavage reaction conversion (F.R300.X_ox_rxn_1)
    @model.parameter(name='Oxidative cleavage reaction',
                      element=F.R300,
                      kind='coupled',
                      baseline = 0.93,
                      distribution=coupled_para_dist['Oxidative cleavage reaction conversion primary'])
    def set_X_oxidativecleavage_conversion1(X_ox_rxn_1):
        F.unit.R300.X_ox_rxn_1 = X_ox_rxn_1
        
    
    @model.parameter(name='Oxidation of intermediates reaction',
                      element=F.R300,
                      kind='coupled',
                      baseline = 0.96,
                      distribution=coupled_para_dist['Oxidative cleavage reaction conversion'])
    def set_X_oxidativecleavage_conversion(X_oxidativecleavage):
        F.unit.R300.X_oxidativecleavage = X_oxidativecleavage
        
    @model.parameter(name='Oxidative decarboxylation reaction',
                         element=F.R300,
                         kind='coupled',
                         baseline = 0.464,
                         distribution=coupled_para_dist['Decarboxylation reaction'])
    def set_decarboxylation_ratio(decarboxylation_ratio):
           F.unit.R300.decarboxylation_ratio = decarboxylation_ratio
      
    @model.parameter(name='Diol recycle fraction',
                      element=F.S611,
                      kind='coupled',
                      baseline = 0.5,
                      distribution=chaospy.Uniform(0.1,0.9))
    def set_diol_recycle(diol_recycle):
        F.unit.S611.specifications[0].args[0]  = diol_recycle
        
    @model.parameter(name='water factor for extraction',
                      element=F.M605,
                      kind='coupled',
                      baseline = 5.5,
                      distribution= chaospy.Uniform(3,9))
    def set_water_factor(X_wat):
        F.M605.specifications[0].args[0]  = X_wat  
          
    @model.parameter(name='Turbogen efficiency',
                      element=F.BT801,
                      kind='coupled',
                      baseline = 0.85,
                      distribution=coupled_para_dist['Turbogen efficiency'])
    def set_tubeff(X_tubeff):
        F.BT801.turbogenerator_efficiency = X_tubeff
        
    @model.parameter(name='Boiler efficiency',
                          element=F.BT801,
                          kind='coupled',
                          baseline = 0.8,
                          distribution=coupled_para_dist['Boiler efficiency'])
    def set_beff(X_beff):
            F.BT801.boiler_efficiency = X_beff
        
    @model.parameter(name='Tungstic acid reusability',
                      element=F.M200,
                      kind='coupled',
                      baseline = 6,
                      distribution=coupled_para_dist['Tungstic acid reusability']
                      )
    def set_cta(cycles_of_reuse):
        F.unit.M200.specifications[0].args[1] = cycles_of_reuse
        
    # --- Use pre-generated samples from Excel for this batch ---
    samples = _slice_rows_from_excel(xlsx_path, "Sheet1", row_indices=row_indices, header_rows=3)

    # Sanity check: columns in Excel must match the number of model parameters
    n_params = len(model.parameters)
    if samples.shape[1] != n_params:
        raise ValueError(f"Excel has {samples.shape[1]} columns but model expects {n_params} parameters "
                         f"(order must match model.parameters).")    

    rule = rule
    np.random.seed(1234) # For consistent results
    model.load_samples(samples)                 # shape: (1000, n_params)
    model.evaluate(notify=notify_runs)
    # Save with batch tag
    if type(row_indices) == int:
        row_tag = f"rows{row_indices}"
    else:
        row_tag = f"rows{row_indices[0]}-{row_indices[-1]}"

    btag = f"batch_{batch_id:02d}"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    table = model.table.copy()

    if isinstance(table.columns, pd.MultiIndex):
        # Convert MultiIndex to single-level string names
        table.columns = (
            table.columns
            .to_flat_index()
            .map(lambda tup: " | ".join(str(x) for x in tup if x not in (None, "", " ")))
        )
    
    table.to_excel(outdir / f"model_table_{btag}_{row_tag}.xlsx", index=False)

    # df_rho, df_p = model.spearman_r()
    # df_rho.to_excel(outdir / f"df_rho_{btag}_{row_tag}.xlsx", index=True)
    # df_p.to_excel(outdir / f"df_p_{btag}_{row_tag}.xlsx", index=True)
    return model

# #%% spearmans rho
# import pandas as pd
# import numpy as np
# from scipy.stats import spearmanr

# # 1. Load Excel
# df = pd.read_excel('C:/Users/lavan/gglk_code/Bioindustrial-Park/biorefineries/oleochemicals/uncertainty_and_sensitivity_analysis_2000/uncertainty_analysis_results_2000.xlsx', sheet_name=0, header=1)

# # 2. Select parameter and metric columns
# param_df = df.iloc[:, 1:60].copy()
# metric_df = df.iloc[:, 60:157].copy()

# # 3. Force numeric types
# param_df = param_df.apply(pd.to_numeric, errors='coerce')
# metric_df = metric_df.apply(pd.to_numeric, errors='coerce')

# # # 4. Drop rows with NaNs
# combined = pd.concat([param_df, metric_df], axis=1).dropna()


# # Count rows before dropping NaNs
# initial_row_count = df.shape[0]
# final_row_count = combined.shape[0]
# rows_dropped = initial_row_count - final_row_count

# print(f"🧹 Dropped {rows_dropped} rows with missing values (from {initial_row_count} to {final_row_count})")

# # 5. Split again
# param_data = combined.iloc[:, :param_df.shape[1]]
# metric_data = combined.iloc[:, param_df.shape[1]:]

# # 6. Identify and print constant columns
# param_const = param_data.columns[param_data.nunique() <= 1].tolist()
# metric_const = metric_data.columns[metric_data.nunique() <= 1].tolist()

# print(f"Dropping {len(param_const)} constant parameter columns:", param_const)
# print(f"Dropping {len(metric_const)} constant metric columns:", metric_const)

# # 7. Drop constant columns
# param_data = param_data.drop(columns=param_const)
# metric_data = metric_data.drop(columns=metric_const)

# # 8. Compute Spearman correlation
# try:
#     rho_matrix, pval_matrix = spearmanr(param_data, metric_data, axis=0)

#     # Extract dimensions
#     n_param = param_data.shape[1]
#     n_metric = metric_data.shape[1]

#     # Extract cross-correlation matrices
#     rho_cross = rho_matrix[:n_param, n_param:]
#     pval_cross = pval_matrix[:n_param, n_param:]

#     # Wrap in DataFrames
#     rho_df = pd.DataFrame(rho_cross, index=param_data.columns, columns=metric_data.columns)
#     pval_df = pd.DataFrame(pval_cross, index=param_data.columns, columns=metric_data.columns)

#     # Export to Excel
#     output_file = 'spearman_results.xlsx'
#     with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
#         rho_df.to_excel(writer, sheet_name='Spearman_rho')
#         pval_df.to_excel(writer, sheet_name='Spearman_pval')

#     print(f"\n✅ Spearman correlation exported to '{output_file}'")

# except FloatingPointError:
#     print("❌ Still encountered a floating-point error during correlation.")

# # #%%
# # #Explaination on LCA methods below
# # #EOL is always included in cradle to grave
# # # If the feedstock CF does not include that biogenic CO₂ was reabsorbed by the crop, then:
# # # You must subtract the amount of fixed carbon from your process to correct for this.
# # # Otherwise, you'd overestimate the GWP, because it would include carbon that was never fossil-derived.
# # # You would base this subtraction on atomic carbon flow (e.g., how much C is in the oil).
# # # If the feedstock CF does not include:
# # # iLUC = indirect land use change
# # # SOC = soil organic carbon losses/gains
# # # Then you should:
# # # Add uncertainty (e.g., ±X%) to the CF to reflect the possible unaccounted emissions or sequestration.
# # # This is commonly done when relying on partial or optimistic LCA models that ignore landscape-level effects.
# # #cradle to gate will be a higher number since EOL is not included
# # #gate to gate will not be able let you do a fair comparision between fossil fuel and bioproducts
