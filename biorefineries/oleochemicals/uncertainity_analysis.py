# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:58:21 2023
@author: lavanya
"""
from biorefineries.oleochemicals.systems_baseline import F
from biorefineries.oleochemicals.chemicals_baseline import chems
from biorefineries.oleochemicals.system_simulate import aa_baseline,azelaic_acid_tea,aa_baseline_groups
from biosteam.evaluation import Model, Metric
from prices_and_GWP_factors import Utility_GWP_factors,isolated_para_dists,coupled_para_dist,environmental_facs_dist,GWP_factors_per_Kg
from warnings import filterwarnings; filterwarnings('ignore')
import chaospy
from chaospy import distributions as shape
import numpy as np
from biorefineries.oleochemicals.system_simulate import azelaic_acid_tea,aa_baseline,aa_baseline_groups
import biosteam as bst


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
conc_hydrochloric_acid = F.stream.conc_hydrochloric_acid
calcium_chloride = F.stream.calcium_chloride_for_cat_sep
resin = F.stream.polystyrene_based_catalyst
biodiesel_catalyst = F.stream.catalyst
liquid_HCl = F.stream.Liquid_HCl
natural_gas = F.stream.natural_gas



#functions to calculate all the metrics
def get_MPSP(): # $/kg of AA
    MPSP = azelaic_acid_tea.solve_price(azelaic_acid)
    return MPSP

def get_MFPP():
    MFPP = azelaic_acid_tea.solve_price(crude_vegetable_oil)
    return MFPP

def get_annual_factor(): return aa_baseline.operating_hours 
def get_total_yield(): return aa_baseline.get_mass_flow(azelaic_acid)/1e6 #To get total yield in .10^6 Kg
def get_purity(): return azelaic_acid.imass['Azelaic_acid']*100/azelaic_acid.F_mass
def get_overall_TCI(): return azelaic_acid_tea.TCI/1e6 #[USD]        
def get_material_cost(): return azelaic_acid_tea.material_cost/1e6 #includes the money spent on ash disposal and brine disposal, [USD/yr]
def get_operating_cost(): return azelaic_acid_tea.AOC/1e6 #[USD/yr]
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
def fatty_acid_price():return fatty_acid_blend.price

#metrics to support facts for the system performance under uncertainty section
def df_related():
    df_unit_groups = bst.UnitGroup.df_from_groups(unit_groups = aa_baseline_groups, fraction=True, 
                                                  scale_fractions_to_positive_values=True)
    
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

def df_related_non_fraction():
    df_unit_groups = bst.UnitGroup.df_from_groups(unit_groups = aa_baseline_groups, fraction=False, 
                                                  scale_fractions_to_positive_values=False)
    
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

#IEC related
def boiler_IEC():
    df_unit_groups = df_related()
    return df_unit_groups['Installed equipment cost']['Boilerturbogenerator']

def Transesterification_IEC():
    df_unit_groups = df_related()
    return df_unit_groups['Installed equipment cost']['Transesterification']

def AA_recovery_IEC():
    df_unit_groups = df_related()
    return df_unit_groups['Installed equipment cost']['Azelaic acid and heavy tails recovery'] 

def Dih_IEC():
    df_unit_groups = df_related()
    return df_unit_groups['Installed equipment cost']['Dihydroxylation'] 

def OX_cl_IEC():
    df_unit_groups = df_related()
    return df_unit_groups['Installed equipment cost']['Oxidative cleavage'] 

def Wastewater_IEC():
    df_unit_groups = df_related()
    return df_unit_groups['Installed equipment cost']['Wastewater treatment']

def Other_facility_IEC():
    df_unit_groups = df_related()
    return df_unit_groups['Installed equipment cost']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)']


#Heating demand related
def Dih_heating_demand():
    df_unit_groups = df_related()
    return df_unit_groups['Heating duty']['Dihydroxylation']

def OX_cl_heating_demand():
    df_unit_groups = df_related()
    return df_unit_groups['Heating duty']['Oxidative cleavage']


def AA_recovery_heating_demand():
    df_unit_groups = df_related()
    return df_unit_groups['Heating duty']['Azelaic acid and heavy tails recovery']

def PA_recovery_heating_demand():
    df_unit_groups = df_related()
    return df_unit_groups['Heating duty']['Pelargonic acid and C5-C8 fraction recovery']

def Transes_heating_demand():
    df_unit_groups = df_related()
    return df_unit_groups['Heating duty']['Transesterification']

#Cooling demand related
def AA_recovery_cooling_demand():
    df_unit_groups = df_related()
    return df_unit_groups['Cooling duty']['Azelaic acid and heavy tails recovery']

def other_facilities_cooling_demand():
    df_unit_groups = df_related()
    return df_unit_groups['Cooling duty']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)']

def Dihydroxylation_cooling_demand():
    df_unit_groups = df_related()
    return df_unit_groups['Cooling duty']['Dihydroxylation']

def Oxidative_cleavage_cooling_demand():
    df_unit_groups = df_related()
    return df_unit_groups['Cooling duty']['Oxidative cleavage']

#Electricity consumption
def PA_recovery_electricity_consumption():
    df_unit_groups = df_related()
    return df_unit_groups['Electricity consumption']['Pelargonic acid and C5-C8 fraction recovery']

def other_facilities_electricity_consumption():
    df_unit_groups = df_related()
    return df_unit_groups['Electricity consumption']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)']

def AA_recovery_electricity_consumption():
    df_unit_groups = df_related()
    return df_unit_groups['Electricity consumption']['Azelaic acid and heavy tails recovery']

def Dih_electricity_consumption():
    df_unit_groups = df_related()
    return df_unit_groups['Electricity consumption']['Dihydroxylation']

def OX_cl_electricity_consumption():
    df_unit_groups = df_related()
    return df_unit_groups['Electricity consumption']['Oxidative cleavage']

def Transesterification_MC():
    df_unit_groups = df_related()
    return df_unit_groups['Material cost']['Transesterification']

def Dih_MC():
    df_unit_groups = df_related()
    return df_unit_groups['Material cost']['Dihydroxylation']

def OX_cl_MC():
    df_unit_groups = df_related()
    return df_unit_groups['Material cost']['Oxidative cleavage'] 
    
def Cat_rec_MC():
    df_unit_groups = df_related()
    return df_unit_groups['Material cost']['Catalyst recovery']  
    
def AA_rec_MC():
    df_unit_groups = df_related()
    return df_unit_groups['Material cost']['Azelaic acid and heavy tails recovery'] 
    
def BT_MC():
    df_unit_groups = df_related()
    return df_unit_groups['Material cost']['Boilerturbogenerator'] 

def WWT_MC():
    df_unit_groups = df_related()
    return df_unit_groups['Material cost']['Wastewater treatment'] 

def OF_MC():
    df_unit_groups = df_related()
    return df_unit_groups['Material cost']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)'] 

def Pel_rec_MC():
    df_unit_groups = df_related()
    return df_unit_groups['Material cost']['Pelargonic acid and C5-C8 fraction recovery']
    
#IEC related
def boiler_IEC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Installed equipment cost [MM$]']['Boilerturbogenerator']

def Transesterification_IEC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Installed equipment cost [MM$]']['Transesterification']

def AA_recovery_IEC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Installed equipment cost [MM$]']['Azelaic acid and heavy tails recovery'] 

def Dih_IEC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Installed equipment cost [MM$]']['Dihydroxylation'] 

def OX_cl_IEC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Installed equipment cost [MM$]']['Oxidative cleavage'] 

#Heating demand related
def Dih_heating_demand_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Heating duty [GJ/hr]']['Dihydroxylation']

def OX_cl_heating_demand_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Heating duty [GJ/hr]']['Oxidative cleavage']


def AA_recovery_heating_demand_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Heating duty [GJ/hr]']['Azelaic acid and heavy tails recovery']

def PA_recovery_heating_demand_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Heating duty [GJ/hr]']['Pelargonic acid and C5-C8 fraction recovery']

def Transes_heating_demand_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Heating duty [GJ/hr]']['Transesterification']

#Cooling demand related
def AA_recovery_cooling_demand_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Cooling duty [GJ/hr]']['Azelaic acid and heavy tails recovery']

def other_facilities_cooling_demand_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Cooling duty [GJ/hr]']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)']

def Dihydroxylation_cooling_demand_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Cooling duty [GJ/hr]']['Dihydroxylation']

def Oxidative_cleavage_cooling_demand_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Cooling duty [GJ/hr]']['Oxidative cleavage']

#Electricity consumption
def PA_recovery_electricity_consumption_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Electricity consumption [MW]']['Pelargonic acid and C5-C8 fraction recovery']

def other_facilities_electricity_consumption_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Electricity consumption [MW]']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)']

def AA_recovery_electricity_consumption_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Electricity consumption [MW]']['Azelaic acid and heavy tails recovery']

def Dih_electricity_consumption_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Electricity consumption [MW]']['Dihydroxylation']

def OX_cl_electricity_consumption_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Electricity consumption [MW]']['Oxidative cleavage']

def Transesterification_MC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Material cost [USD/hr]']['Transesterification']

def Dih_MC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Material cost [USD/hr]']['Dihydroxylation']

def OX_cl_MC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Material cost [USD/hr]']['Oxidative cleavage'] 
    
def Cat_rec_MC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Material cost [USD/hr]']['Catalyst recovery']  
    
def AA_rec_MC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Material cost [USD/hr]']['Azelaic acid and heavy tails recovery'] 
    
def BT_MC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Material cost [USD/hr]']['Boilerturbogenerator'] 

def WWT_MC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Material cost [USD/hr]']['Wastewater treatment'] 

def OF_MC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Material cost [USD/hr]']['Other Facilities (cooling tower,air distribution and cleaning auxiliaries)'] 

def Pel_rec_MC_1():
    df_unit_groups = df_related_non_fraction()
    return df_unit_groups['Material cost [USD/hr]']['Pelargonic acid and C5-C8 fraction recovery']
    
def feedstock_frac_LCA():
    return (get_feedstock_GWP()*100/get_net_GWP_PA_1())/5

def other_materials_frac_LCA():
    return (get_other_materials_impact()*100/get_net_GWP_PA_1())/5

def natural_gas_frac_LCA():
    return (get_ng_GWP()*100/get_net_GWP_PA_1())/5

def electricity_cons_frac_LCA():
    return  (get_electricity_consumption_GWP()*100/get_net_GWP_PA_1())/5

def direct_emission_frac_LCA():
    return (get_total_direct_emissions_GWP()*100/get_net_GWP_PA_1())/5

def products_displacement_total_frac_LCA():
    return (-sum_displacement*100/get_net_GWP_PA_1())/5

    

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
                 F.regeneration_acid_to_HydrolysisSystem,
                 F.crude_methanol]

#The emissions originate from R300, D403,D404,D601,HX608,U902,R902,BT801,CT901,PWT901
#Emissions were assumed to be carbon emissions
emissions = [i for i in aa_baseline.products if i not in extra_streams and i not in main_product and i not in products and i not in F.ADP701.outs]
#Total emissions (TE) (unitless)
get_total_emissions_GWP = TE_GWP = lambda: sum([stream.get_atomic_flow('C') for stream in emissions]) * chems.CO2.MW / azelaic_acid.F_mass
#End of life emissions (EOL) emissions are indirect emissions
#All products were included in the EOL calculations
get_EOL_GWP = EOL_GWP = lambda: (sum([stream.get_atomic_flow('C')*chems.CO2.MW for stream in products]) + azelaic_acid.get_atomic_flow('C')*chems.CO2.MW )/azelaic_acid.F_mass #(unitless)
direct_emmisions = get_total_direct_emissions_GWP = lambda: get_total_emissions_GWP() + get_EOL_GWP() #(unitless)

#Boilerturbogenerator also contributes to direct emissions #(unitless)
direct_emmisions_from_BT = get_total_direct_BT_emissions_GWP = lambda: (sum([i.get_atomic_flow('C') for i in F.BT801.outs])*chems.CO2.MW/ azelaic_acid.F_mass)*get_total_direct_emissions_GWP()/get_total_emissions_GWP()
direct_emmisions_without_BT = get_total_non_BT_direct_emissions_GWP = lambda: get_total_direct_emissions_GWP() - get_total_direct_BT_emissions_GWP()

#Material related impacts
get_feedstock_GWP = lambda: aa_baseline.get_material_impact(crude_vegetable_oil,key = 'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
get_internal_materials_impact = lambda:(aa_baseline.get_material_impact(resin,key ='GWP100')+aa_baseline.get_material_impact(liquid_HCl,key ='GWP100'))/aa_baseline.get_mass_flow(azelaic_acid)
get_other_materials_impact = lambda: (aa_baseline.get_total_feeds_impact('GWP100')/aa_baseline.get_mass_flow(azelaic_acid))-get_feedstock_GWP() +get_internal_materials_impact()

#Electrcity related impacts
#Electricity is produced using steam and natural gas 
get_electricity_consumption = lambda: aa_baseline.get_electricity_consumption() #KWH per year
get_electricity_consumption_GWP = lambda: get_electricity_consumption()*Utility_GWP_factors['Electricity']/aa_baseline.get_mass_flow(azelaic_acid)
#Impact of the total natural gas consumed in the process via the boiler turbogenerator
get_ng_GWP = lambda: aa_baseline.get_material_impact(natural_gas,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid) #annual GWP of ng/annual flowrate of AA


def get_system_GWP():
    return(sum([get_feedstock_GWP(),
                get_other_materials_impact(),
                get_ng_GWP(),
                get_electricity_consumption_GWP(),#TODO: check if summation is right
                get_total_emissions_GWP(),
                get_EOL_GWP()]))

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
#TODO: should this be done or not
economic_fa_fraction = aa_baseline.get_market_value(fatty_acid_blend)*get_system_GWP()/get_total_product_market_value()          
economic_methanol_fraction = aa_baseline.get_market_value(crude_methanol)*get_system_GWP()/get_total_product_market_value()          

sum_economic = economic_C5_C9_fraction+economic_pa_fraction+economic_glycerol_fraction+economic_fa_fraction+economic_methanol_fraction

def get_economic_based_AA_GWP():
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
mass_C5_C9_fraction = aa_baseline.get_mass_flow(recovered_C5_to_C9_MCA_fraction)*get_system_GWP()/get_total_product_mass()
mass_pa_fraction = aa_baseline.get_mass_flow(pelargonic_acid_rich_fraction)*get_system_GWP()/get_total_product_mass()
mass_glycerol_fraction = aa_baseline.get_mass_flow(crude_glycerol)*get_system_GWP()/get_total_product_mass() 
mass_fa_fraction = aa_baseline.get_mass_flow(fatty_acid_blend)*get_system_GWP()/get_total_product_mass()
mass_methanol_fraction = aa_baseline.get_mass_flow(crude_methanol)*get_system_GWP()/get_total_product_mass()

sum_mass =  mass_C5_C9_fraction+mass_pa_fraction+mass_glycerol_fraction+mass_fa_fraction+mass_methanol_fraction

def get_mass_based_AA_GWP():
    mass_aa_fraction = aa_baseline.get_mass_flow(azelaic_acid)*get_system_GWP()/get_total_product_mass()
    return mass_aa_fraction #Kg CO2 per kg of AA
    
#GWP related functions
#Displacement scenarios
#Azelaic acid - Adipic acid
#Pelargonic acid - Glyohosate, or acetic acid
#C5-C9 - This fraction contains octanoic acid, which can be used to produce octanol-1
#TODO: discuss fatty acid production from coconut oil
#Fatty acid blend - biodiesel

#Values per kg of product

#herbicide
#Glyphosate production (ecoinvent)
PA_1 = 11.239* aa_baseline.get_mass_flow(pelargonic_acid_rich_fraction)/aa_baseline.get_mass_flow(azelaic_acid)

#replaces conventionally produced acetic acid which is used to make peroxyacetic acid
#pelargonic acid based peroxy ca are used for the same purposes
PA_2 = 1.77*aa_baseline.get_mass_flow(pelargonic_acid_rich_fraction)/aa_baseline.get_mass_flow(azelaic_acid)

#biodiesel
#Based on conventional diesel from US Crude refineries
FA = 0.67*aa_baseline.get_mass_flow(fatty_acid_blend)/aa_baseline.get_mass_flow(azelaic_acid)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5480234/
#novel production method for caproic acid production
MCA = 0.23*aa_baseline.get_mass_flow(recovered_C5_to_C9_MCA_fraction)/aa_baseline.get_mass_flow(azelaic_acid)

#glycerine production, from epichlorohydrin, RoW (Ecoinvent)
CG = 4.2879e+0*aa_baseline.get_mass_flow(crude_glycerol)/aa_baseline.get_mass_flow(azelaic_acid)

#methanol production conventional from natural gas (GREET)
CM = 0.4955*aa_baseline.get_mass_flow(crude_methanol)/aa_baseline.get_mass_flow(azelaic_acid)


sum_displacement = PA_1 + FA + MCA + + CG + CM
def get_net_GWP_PA_1():
    return(sum([get_feedstock_GWP(),
                get_other_materials_impact(),
                get_ng_GWP(),
                get_electricity_consumption_GWP(),#TODO: check if summation is right
                get_total_emissions_GWP(),
                get_EOL_GWP(), -sum_displacement
                ]))

#######################################################################################################################3
all_metrics = [
    Metric('MPSP', get_MPSP, '$/kg'),
    Metric('Total product yield', get_total_yield, '10^6 kg/yr'),
    Metric('Product purity', get_purity, '%'),
    Metric('Total capital investment', get_overall_TCI, '10^6 $'),
    Metric('Annual operating cost', get_operating_cost, '10^6 $/yr'),
    Metric('Net present value',solve_NPV,'10^6 $'),
    Metric('Maximum feedstock purchase price',get_MFPP,'$/kg'),

#TODO: units of these metrics are wrong
    Metric('Transesterification_MC', Transesterification_MC, '10^6 $/yr'),
    Metric('Dih_MC', Dih_MC, '10^6 $/yr'),
    Metric('OX_cl_MC', OX_cl_MC, '10^6 $/yr'),
    Metric('Cat_rec_MC', Cat_rec_MC, '10^6 $/yr'),
    Metric('Pel_rec_MC', Pel_rec_MC, '10^6 $/yr'),
    Metric('AA_rec_MC', AA_rec_MC, '10^6 $/yr'),
    Metric('BT_MC', BT_MC, '10^6 $/yr'),
    Metric('WWT_MC', WWT_MC, '10^6 $/yr'),
    Metric('OF_MC', OF_MC, '10^6 $/yr'),
    
    Metric('Transesterification_MC_1', Transesterification_MC_1, '10^6 $/yr'),
    Metric('Dih_MC_1', Dih_MC_1, '10^6 $/yr'),
    Metric('OX_cl_MC_1', OX_cl_MC_1, '10^6 $/yr'),
    Metric('Cat_rec_MC_1', Cat_rec_MC_1, '10^6 $/yr'),
    Metric('Pel_rec_MC_1', Pel_rec_MC_1, '10^6 $/yr'),
    Metric('AA_rec_MC_1', AA_rec_MC_1, '10^6 $/yr'),
    Metric('BT_MC_1', BT_MC_1, '10^6 $/yr'),
    Metric('WWT_MC_1', WWT_MC_1, '10^6 $/yr'),
    Metric('OF_MC_1', OF_MC_1, '10^6 $/yr'),
    
    Metric('Total', get_system_heating_demand, '10^6 MJ/yr', 'Heating demand'),
    Metric('Total', get_system_cooling_demand,'10^6 MJ/yr', 'Cooling demand'),
    Metric('Price of fatty acid blend',fatty_acid_price,'$/kg'),

    
    Metric('boiler_IEC',boiler_IEC,'10^6'),
    Metric('Transesterification_IEC',Transesterification_IEC,'10^6'),
    Metric('AA_recovery_IEC',AA_recovery_IEC,'10^6'),
    Metric('Dih_IEC',Dih_IEC,'10^6'),
    Metric('OX_cl_IEC',OX_cl_IEC,'10^6'),
    Metric('Wastewater_IEC',Wastewater_IEC,'10^6'),
    Metric('Other_facility_IEC',Other_facility_IEC,'10^6'),
    
    
    Metric('boiler_IEC_1',boiler_IEC_1,'10^6'),
    Metric('Transesterification_IEC_1',Transesterification_IEC_1,'10^6'),
    Metric('AA_recovery_IEC_1',AA_recovery_IEC_1,'10^6'),
    Metric('Dih_IEC_1',Dih_IEC_1,'10^6'),
    Metric('OX_cl_IEC_1',OX_cl_IEC_1,'10^6'),
        
    
    Metric('Dih_heating_demand',Dih_heating_demand,'10^6 MJ/yr', 'Heating demand'),
    Metric('OX_cl_heating_demand',OX_cl_heating_demand,'10^6 MJ/yr', 'Heating demand'),    
    Metric('AA_recovery_heating_demand',AA_recovery_heating_demand,'10^6 MJ/yr', 'Heating demand'),
    Metric('PA_recovery_heating_demand',PA_recovery_heating_demand,'10^6 MJ/yr', 'Heating demand'),
    Metric('Transes_heating_demand',Transes_heating_demand,'10^6 MJ/yr', 'Heating demand'),
    Metric('AA_recovery_cooling_demand',AA_recovery_cooling_demand,'10^6 MJ/yr', 'Cooling demand'),
    Metric('other_facilities_cooling_demand',other_facilities_cooling_demand,'10^6 MJ/yr', 'Cooling demand'),
    Metric('Dihydroxylation_cooling_demand',Dihydroxylation_cooling_demand,'10^6 MJ/yr', 'Cooling demand'),
    Metric('Oxidative_cleavage_cooling_demand',Oxidative_cleavage_cooling_demand,'10^6 MJ/yr', 'Cooling demand'),
    Metric('PA_recovery_electricity_consumption',PA_recovery_electricity_consumption,'KW', 'Electricity consumption'),
    Metric('other_facilities_electricity_consumption',other_facilities_electricity_consumption,'KW', 'Electricity consumption'),
    Metric('AA_recovery_electricity_consumption',AA_recovery_electricity_consumption,'KW', 'Electricity consumption'),
    Metric('Dih_electricity_consumption',Dih_electricity_consumption,'KW', 'Electricity consumption'),
    Metric('OX_cl_electricity_consumption',OX_cl_electricity_consumption,'KW', 'Electricity consumption'),

    Metric('Dih_heating_demand_1',Dih_heating_demand_1,'10^6 MJ/yr', 'Heating demand'),
    Metric('OX_cl_heating_demand_1',OX_cl_heating_demand_1,'10^6 MJ/yr', 'Heating demand'),    
    Metric('AA_recovery_heating_demand_1',AA_recovery_heating_demand_1,'10^6 MJ/yr', 'Heating demand'),
    Metric('PA_recovery_heating_demand_1',PA_recovery_heating_demand_1,'10^6 MJ/yr', 'Heating demand'),
    Metric('Transes_heating_demand_1',Transes_heating_demand_1,'10^6 MJ/yr', 'Heating demand'),
    Metric('AA_recovery_cooling_demand_1',AA_recovery_cooling_demand_1,'10^6 MJ/yr', 'Cooling demand'),
    Metric('other_facilities_cooling_demand_1',other_facilities_cooling_demand_1,'10^6 MJ/yr', 'Cooling demand'),
    Metric('Dihydroxylation_cooling_demand_1',Dihydroxylation_cooling_demand_1,'10^6 MJ/yr', 'Cooling demand'),
    Metric('Oxidative_cleavage_cooling_demand_1',Oxidative_cleavage_cooling_demand_1,'10^6 MJ/yr', 'Cooling demand'),
    Metric('PA_recovery_electricity_consumption_1',PA_recovery_electricity_consumption_1,'KW', 'Electricity consumption'),
    Metric('other_facilities_electricity_consumption_1',other_facilities_electricity_consumption_1,'KW', 'Electricity consumption'),
    Metric('AA_recovery_electricity_consumption_1',AA_recovery_electricity_consumption_1,'KW', 'Electricity consumption'),
    Metric('Dih_electricity_consumption_1',Dih_electricity_consumption_1,'KW', 'Electricity consumption'),
    Metric('OX_cl_electricity_consumption_1',OX_cl_electricity_consumption_1,'KW', 'Electricity consumption'),
        
    Metric('Mass allocation LCA', get_mass_based_AA_GWP, 'kg CO2-eq/kg', 'LCA'),
    Metric('Economic allocation LCA', get_economic_based_AA_GWP, 'kg CO2-eq/kg', 'LCA'),
    Metric('Net GWP displacement 1', get_net_GWP_PA_1, 'kg CO2-eq/kg', 'LCA'),
    
    Metric('feedstock_frac_LCA', feedstock_frac_LCA, 'kg CO2-eq/kg', 'LCA'),
    Metric('other_materials_frac_LCA', other_materials_frac_LCA, 'kg CO2-eq/kg', 'LCA'),
    Metric('natural_gas_frac_LCA', natural_gas_frac_LCA, 'kg CO2-eq/kg', 'LCA'),
    Metric('electricity_cons_frac_LCA', electricity_cons_frac_LCA, 'kg CO2-eq/kg', 'LCA'),
    Metric('direct_emission_frac_LCA', direct_emission_frac_LCA, 'kg CO2-eq/kg', 'LCA'),
    Metric('products_displacement_total_frac_LCA', products_displacement_total_frac_LCA, 'kg CO2-eq/kg', 'LCA'),
    
    ]

#TODO: try changing feedstock flow rate from 5% of HOSOY production to 30% of HOSun

def conduct_uncertainity_analysis(system = aa_baseline, #name of the system
                                  metrics = all_metrics, #list of Metric objects
                                  number_of_runs = 3, #number of 
                                  exception_hook = 'raise',#raise/warn
                                  rule = 'L', #for For Latin-Hypercube sampling
                                  notify_runs = 2, #runs after which you are notified 
                                  indicator = 'GWP100',
                                  feedstock_type = 'HoySoy_oil' #can also provide 'HoSun_oil
                                 ):
    model = Model(aa_baseline,metrics = all_metrics,exception_hook= exception_hook)     
    ## price distributions
    @model.parameter(name='Crude oil price',
                      element=crude_vegetable_oil, 
                      kind='isolated',
                      units='USD/kg',
                      baseline= 1.70,
                      distribution= isolated_para_dists['Crude oil price'][feedstock_type])
    def set_feedstock_price(feedstock_price):
            crude_vegetable_oil.price = feedstock_price
    
    @model.parameter(name='Tungstic acid price',
                      element=fresh_tungsten_catalyst, 
                      kind='isolated', 
                      units='USD/kg',
                      baseline = 35.2,
                      distribution= isolated_para_dists['Tungstic acid price'])
    def set_tungstencat_price(tungstencat_price):
        fresh_tungsten_catalyst.price = tungstencat_price
    
    
    @model.parameter(name='Cobalt acetate catalyst price',
                      element=fresh_cobalt_catalyst, 
                      kind='isolated', 
                      units='USD/kg',
                      baseline=26,
                      distribution= isolated_para_dists['Cobalt acetate price'])
    def set_cobaltcat_price(fresh_cobalt_catalyst_price):
        F.fresh_cobalt_catalyst_stream.price = fresh_cobalt_catalyst_price
    
    @model.parameter(name='Heptane solvent price',
                      element=fresh_solvent, 
                      kind='isolated',
                      units='USD/kg',
                      baseline = 0.68,
                      distribution=isolated_para_dists['Heptane solvent price'])
    def set_solvent_price(fresh_solvent_price):
        fresh_solvent.price = fresh_solvent_price
    
    @model.parameter(name='Hydrogen peroxide price',
                      element=fresh_HP, 
                      kind='isolated',
                      baseline = 1.5,
                      units='USD/kg',
                      distribution=isolated_para_dists['Hydrogen peroxide price'])
    def set_fresh_HP_price(fresh_HP_price):
        fresh_HP.price = fresh_HP_price
    
    @model.parameter(name='Natural gas price',
                      element=F.BT801,
                      baseline = 0.234,
                      kind='isolated',
                      units='USD/kg',
                      distribution=isolated_para_dists['Natural gas price'])
    def set_natural_gas_price(natural_gas_price):
        F.BT801.natural_gas_price = natural_gas_price    
    
    
    @model.parameter(name='Conc HCl price',
                      element=conc_hydrochloric_acid,
                      kind='isolated',
                      baseline = 0.129,
                      units='USD/kg',
                      distribution=isolated_para_dists['Conc HCl price'])
    def set_HCl_price(conc_hydrochloric_acid_price):
        conc_hydrochloric_acid.price = conc_hydrochloric_acid_price
    
    
    @model.parameter(name='Calcium chloride price',
                  element=calcium_chloride,
                  kind='isolated',
                  units='USD/kg',
                  baseline = 0.823,
                  distribution=isolated_para_dists['Calcium chloride price'])
    def set_calcium_chloride_price(calcium_chloride_price):
        calcium_chloride.price = calcium_chloride_price
    
    
    @model.parameter(name='Hydrolysis resin price',
                      element=resin,
                      kind='isolated',
                      units='USD/kg',
                      baseline = 106,
                      distribution=isolated_para_dists['Hydrolysis resin price'])
    def set_polystyrene_based_catalyst_price(resin_price):
        resin.price = resin_price 
    
    #Product prices    
    @model.parameter(name='Crude glycerol price',
                      element=crude_glycerol,
                      kind='isolated',
                      baseline = 0.155,
                      units='USD/kg',
                      distribution=isolated_para_dists['Crude glycerol price'])
    def set_crude_glycerol_price(crude_glycerol_price):
        crude_glycerol.price = crude_glycerol_price
        
    @model.parameter(name='Pelargonic acid price',
                      element=pelargonic_acid_rich_fraction,
                      kind='isolated',
                      baseline=5.36,
                      units='USD/kg',
                      distribution=isolated_para_dists['Pelargonic acid price'])
    def set_pelargonic_acid_price(pelargonic_acid_price):
        pelargonic_acid_rich_fraction.price = pelargonic_acid_price 
        
    @model.parameter(name='Fatty acid blend price',
                      element= fatty_acid_blend,
                      kind='isolated',
                      baseline= 0.0231,
                      units='USD/kg',
                      distribution=chaospy.Uniform(0.88*0.9,0.88*1.1))
    def set_fatty_acid_blend_price(fa_blend_price):
        joule_to_megajoule = 1.0E-6
        ratio_2023from2021 = 0.992
        LHV_of_high_oleic_soybean_oil_biodiesel = 37.7
        fatty_acid_blend.price = fatty_acid_blend.LHV*joule_to_megajoule*ratio_2023from2021*fa_blend_price/LHV_of_high_oleic_soybean_oil_biodiesel

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
        
    @model.parameter(name='Hydrogen peroxide composition',
                      element=F.stream.fresh_HP,
                      kind='coupled',
                      baseline = 0.5,
                      distribution= coupled_para_dist['HP_concentration'])
    def set_hydrogen_peroxide_comp(X_hp):
        F.stream.fresh_HP.imass['Hydrogen_peroxide'] = X_hp
        F.stream.fresh_HP.imass['Water'] = 1-X_hp
      
    ## Oxidative cleavage reaction conversion (F.R300.X_ox_rxn_1)
    @model.parameter(name='Oxidative cleavage reaction conversion primary',
                      element=F.R300,
                      kind='coupled',
                      baseline = 0.95,
                      distribution=coupled_para_dist['Oxidative cleavage reaction conversion primary'])
    def set_X_oxidativecleavage_conversion1(X_ox_rxn_1):
        F.unit.R300.X_ox_rxn_1 = X_ox_rxn_1
        
    
    @model.parameter(name='Oxidative cleavage reaction conversion',
                      element=F.R300,
                      kind='coupled',
                      baseline = 0.96,
                      distribution=coupled_para_dist['Oxidative cleavage reaction conversion'])
    def set_X_oxidativecleavage_conversion(X_oxidativecleavage):
        F.unit.R300.X_oxidativecleavage = X_oxidativecleavage
    
    @model.parameter(name='Dihydroxylation reaction conversion',
                      element=F.R200,
                      kind='coupled',
                      baseline = 0.86,
                      distribution=coupled_para_dist['Dihydroxylation reaction conversion'])
    def set_X_dih_conversion(X_dih):
        F.unit.R200.X_dih = X_dih
    
    @model.parameter(name='Diol recycle fraction',
                      element=F.S611,
                      kind='coupled',
                      baseline = 0.5,
                      distribution=chaospy.Uniform(0.1,0.9))
    def set_diol_recycle(diol_recycle):
        F.unit.S611.specifications[0].args[0]  = diol_recycle
        
          
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
                                    
    #environmental_facs_dist
    @model.parameter(name='crude oil GWP',
                      element=crude_vegetable_oil,
                      kind='isolated',
                      baseline = 0.493,
                      distribution=environmental_facs_dist['Oil_GWP'][feedstock_type]
                      )
    def crude_vegetable_oil_gwp(crude_vegetable_oil_gwp):
        crude_vegetable_oil.characterization_factors[indicator] = crude_vegetable_oil_gwp
        
    @model.parameter(name='tungstic acid GWP',
                      element=fresh_tungsten_catalyst,
                      kind='isolated',
                      baseline = 68.5,
                      distribution=environmental_facs_dist['Tungstic_acid']
                      )
    def tungstic_acid_gwp(tungstic_acid_gwp):
        fresh_tungsten_catalyst.characterization_factors[indicator] = tungstic_acid_gwp 
        
    @model.parameter(name='cobalt acetate GWP',
                      element=fresh_cobalt_catalyst,
                      kind='isolated',
                      baseline = 8.89,
                      distribution=environmental_facs_dist['Cobalt_acetate']
                      )
    def cobalt_acetate_gwp(cobalt_acetate_gwp):
        fresh_cobalt_catalyst.characterization_factors[indicator] = cobalt_acetate_gwp 
    
    @model.parameter(name='resin GWP',
                      element=resin,
                      kind='isolated',
                      baseline = 0.1324,
                      distribution=environmental_facs_dist['Resin']
                      )
    def Resin_gwp(Resin_gwp):
        resin.characterization_factors[indicator]  = Resin_gwp    
        
    @model.parameter(name='HP_GWP',
                      element=fresh_HP,
                      kind='isolated',
                      baseline = 1.08,
                      distribution=environmental_facs_dist['Hydrogen peroxide']
                      )
    def HP_GWP(HP_GWP):
        fresh_HP.characterization_factors[indicator]  = HP_GWP       
         
    
    @model.parameter(name='fatty acid blend GWP',
                      element=fatty_acid_blend,
                      kind='isolated',
                      baseline = 0.67,
                      distribution=chaospy.Uniform(0.67*0.5,0.67*1.5),
                      )
    def fatty_acid_blend_gwp(fatty_acid_blend_gwp):
        fatty_acid_blend.characterization_factors[indicator]  = fatty_acid_blend_gwp     
        
     
    @model.parameter(name='C5_C9 fraction GWP',
                      element=recovered_C5_to_C9_MCA_fraction,
                      kind='isolated',
                      baseline = 1.05,
                      distribution=chaospy.Uniform(1.05*0.5,1.05*1.5),
                      )
    def recovered_C5_to_C9_MCA_fraction_gwp(recovered_C5_to_C9_MCA_fraction_gwp):
        recovered_C5_to_C9_MCA_fraction.characterization_factors[indicator]  = recovered_C5_to_C9_MCA_fraction_gwp    
    
    @model.parameter(name='pelargonic acid GWP 1',
                      element=pelargonic_acid_rich_fraction,
                      kind='isolated',
                      baseline = 11.239,
                      distribution=chaospy.Uniform(11.239*0.5,11.239*1.5),
                      )
    def pelargonic_acid_rich_fraction_gwp_1(pelargonic_acid_rich_fraction_gwp):
        pelargonic_acid_rich_fraction.characterization_factors[indicator]  = pelargonic_acid_rich_fraction_gwp    

    @model.parameter(name='pelargonic acid GWP 2',
                      element=pelargonic_acid_rich_fraction,
                      kind='isolated',
                      baseline = 1.77,
                      distribution=chaospy.Uniform(1.77*0.5,1.77*1.5),
                      )
    def pelargonic_acid_rich_fraction_gwp_2(pelargonic_acid_rich_fraction_gwp):
        pelargonic_acid_rich_fraction.characterization_factors[indicator]  = pelargonic_acid_rich_fraction_gwp    
        

        
    if number_of_runs > 0:
        rule = rule
        np.random.seed(1234) # For consistent results
        samples = model.sample(number_of_runs,rule)
        model.load_samples(samples)
        model.evaluate(notify=notify_runs) #notify each time after 2 runs
        model.show()
        model.table.to_excel('model_table_.xlsx')
        df_rho,df_p = model.spearman_r()
        df_rho.to_excel('df_rho_.xlsx')
        df_p.to_excel('df_p_.xlsx')
    else:
        model.show()
    return model


#TODO: discuss composition variation
#TODO: dont vary things that do not affect the prodyct
    # @model.parameter(name='Crude oil composition',
    #                   element=F.stream.fresh_HP,
    #                   kind='coupled',
    #                   distribution= coupled_para_dist['Crude_oil_composition'])
    # def set_crude_oil_comp(X_co):
    #     vsg= high_oleic_vari_adjusted['Vistive gold']
    #     ratios_wrt_ooo = {'ppp/rem':vsg['PPP']/(98-vsg['OOO']),
    #                       'lll/rem':vsg['LLL']/(98-vsg['OOO']),
    #                       'sss/rem':vsg['SSS']/(98-vsg['OOO']),
    #                       'lnlnln/rem':vsg['LnLnLn']/(98-vsg['OOO'])} 
    #     F.crude_vegetable_oil.imass['OOO'] = X_co    
    #     total_fatty_acid_mass = 98
    #     remaining_fraction = total_fatty_acid_mass - X_co
    #     F.crude_vegetable_oil.imass['PPP'] = round(ratios_wrt_ooo['ppp/rem']*remaining_fraction)
    #     F.crude_vegetable_oil.imass['SSS'] = round(ratios_wrt_ooo['sss/rem']*remaining_fraction)
    #     F.crude_vegetable_oil.imass['LLL'] = round(ratios_wrt_ooo['lll/rem']*remaining_fraction)
    #     F.crude_vegetable_oil.imass['LnLnLn'] =round(ratios_wrt_ooo['lnlnln/rem']*remaining_fraction)
    #     F.crude_vegetable_oil.imass['PL']= 1
    #     F.crude_vegetable_oil.imass['MAG']=0
    #     F.crude_vegetable_oil.imass['DAG']=0
    #     F.crude_vegetable_oil.imass['Water']=0.05
    #     F.crude_vegetable_oil.imass['Oleic_acid']=0.95 
    #     F.crude_vegetable_oil.F_mass = 34635.99842

        