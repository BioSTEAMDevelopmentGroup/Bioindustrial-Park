#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:52:43 2024

@author: bianco3

# results from this model are used to create the files:
    - Jet fuel price for python
    - Jet_fuel_GHG_for_python
    
These are used to combine ethanol with SAF results in the next model file 05_combine_ethanol_with_SAF


"""

#%% 
# To run the model

# First: change ethanol flow in line 685 to the jet capacity simulated. 
# Just in case, also change it in console doing-> ethanol.F_mass = 1960.7 (eg)
# These capacities come from the 01_feedstock_transport_model (according to the blending target)
# I am using
    ## 80 MM Gal ethanol/yr = 28712.2 kg/hr (it's actually around 76 MM Gal)
    ## 71 MM Gal ethanol/yr = 26751.5 kg/hr
    ## 62 MM Gal ethanol/yr = 23407.5 kg/hr
    ## 45 MM Gal ethanol/yr = 17008.4 kg/hr
    ## 44 MM Gal ethanol/yr = 16719.7 kg/hr
    ## 32 MM Gal ethanol/yr = 11992.5 kg/hr
    ## 31 MM Gal ethanol/yr = 11703.79 kg/hr
    ## 14 MM Gal ethanol/yr = 5304.6 kg/hr
    ## 5 MM Gal ethanol/yr = 1960.7 kg/hr
   
# For each size, run for 2 prices and 2 GWP of ethanol (Only for first size run for 3 to check linearity)
# change ethanol price-> ethanol.price = 0
# I used 0, 2.5 and 1  (last one only for first size)

# change ethanol CI -> F.ethanol.characterization_factors['GWP100'] = -0.19
# I used 0, -1.08 and 0.809 (last one only for first size)

# do-> model = create_states_model()
# run doing -> evaluate_SAF(N=1000, notify_runs=10, model = model)
# If there is an error, the model will keep running but the capital investment would be 0, so we have to discard that sample for all runs
# If that happens, run more than 1000 times each to be able to have extra samples.
# Always use same samples for all runs.

#%%
import numpy as np
from chaospy import distributions as shape
import biosteam as bst
import blocs as blc
import thermosteam as tmo
import pandas as pd
from biosteam import TEA
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools.parameter import Setter
from biorefineries.SAF._chemicals import SAF_chemicals
# from biorefineries.SAF.ATJ.separate.system_ethanol import MESP
from biorefineries.SAF.ATJ.separate.system_upgrading import F, SAF_sys
from biorefineries.SAF.ATJ._process_settings import price, GWP_CFs, load_preferences_and_process_settings
from tea_loc_model import create_SAF_copr_tea_with_blocs
# from biorefineries.SAF._tea import create_SAF_coprocessing_tea, SAF_coprocessing_TEA

import os
# import pandas as pd
import time

from warnings import warn
from warnings import filterwarnings; filterwarnings('ignore')



bst.settings.set_thermo(SAF_chemicals, cache=True)

_gal_per_m3 = 1000/3.785412
_L_per_gal = 3.785412
_GGE_per_J = 1 / 116090 / 1055.06 # 1 BTU = 1055.06 J, gasoline LHV is 116 090 BTU per gal from https://afdc.energy.gov/fuels/properties
_kJpersec_to_kJhr = 3600

__all__ = ('create_model')

#%% 

sys = SAF_sys(product_storage=False, WWTC=False, BoilerTurbo=False, hydrogenation_distillation=False,h2_purchase=False,opportunity_cost=True)
df_unit_groups = bst.UnitGroup.group_by_area(sys.units)
# '400' is the group for SAF production

# F (Flowsheet) defined in system_upgrading.py

# To change ethanol flow> the default value is set in system_upgrading.py
F.ethanol.set_total_flow(30223.17,'kg/hr') # 80 MM Gal ethanol per year 
# If I want to check the value in the console I can do:
    # F.ethanol.get_total_flow('kg/hr')


load_preferences_and_process_settings(T='K',
                                      flow_units='kg/hr',
                                      N=100,
                                      P_units='Pa',
                                      CE=798, # Average 2023 https://toweringskills.com/financial-analysis/cost-indices/
                                      indicator='GWP100',
                                      electricity_EI=GWP_CFs['electricity'],
                                      electricity_price=price['electricity'])
sys.set_tolerance(rmol=1e-6, mol=1e-5, maxiter=400)

feedstock = F.ethanol

#%%
tea_SAF = create_SAF_copr_tea_with_blocs(sys=sys,steam_distribution=0.0, water_supply_cooling_pumping=0.0, 
                                      water_distribution=0.0, electric_substation_and_distribution=0.0,
                                      gas_supply_and_distribution=0.009, comminication=0.0, safety_installation=0.013,
                                      building=0.07, yard_works=0.0, contingency_new=0.1, land=0.0, labor_cost=3763935,
                                      sanitary_waste_disposal=0.0, ethanol_feed = feedstock, jet_fuel_product = F.jet_fuel,
                                      state = 'Delaware')
sys.operating_hours = tea_SAF.operating_days * 24
tea_SAF.federal_income_tax = 0.21 # corporate federal income tax
tea_SAF.ethanol_feed = F.ethanol
tea_SAF.jet_fuel_product = F.jet_fuel
tea_SAF.jet_fuel_group = bst.UnitGroup('Jet fuel group', sys.units) # Assuming all units qualify
# tea_SAF.biodiesel_product = F.diesel

# PROBAR BORRANDO LO OTRO

def create_SAF_tea_withBlocs():
    tea_SAF = create_SAF_copr_tea_with_blocs(sys=sys,steam_distribution=0.0, water_supply_cooling_pumping=0.0, 
                                          water_distribution=0.0, electric_substation_and_distribution=0.0,
                                          gas_supply_and_distribution=0.009, comminication=0.0, safety_installation=0.013,
                                          building=0.07, yard_works=0.0, contingency_new=0.1, land=0.0, labor_cost=3763935,
                                          sanitary_waste_disposal=0.0, ethanol_feed = feedstock, jet_fuel_product = F.jet_fuel,
                                          state = 'Delaware')
    sys.operating_hours = tea_SAF.operating_days * 24
    tea_SAF.federal_income_tax = 0.21 # corporate federal income tax
    tea_SAF.ethanol_feed = F.ethanol
    tea_SAF.jet_fuel_product = F.jet_fuel
    tea_SAF.jet_fuel_group = bst.UnitGroup('Jet fuel group', sys.units) # Assuming all units qualify
    
    return tea_SAF

tea = create_SAF_tea_withBlocs()

#%%

MESP_with_transport = 2.5 # price in USD/kg ethanol # price includes transportation

def set_price_of_streams():
    for i in sys.streams:
        if i.ID in price.keys():
            i.price = price[i.ID]
    F.ethanol.price = MESP_with_transport # price in USD/kg ethanol


    
    
def set_GWP_of_streams(indicator):
    F.caustic.set_CF(key='GWP100', value=GWP_CFs['caustic']) # caustic in WWT
    for i in sys.streams:
        if i.ID in GWP_CFs.keys():                
            i.characterization_factors[indicator]= GWP_CFs[i.ID] 
    
    
set_prices = set_price_of_streams()
set_GWP = set_GWP_of_streams(indicator='GWP100')


F.natural_gas_price = bst.stream_utility_prices['Natural gas']


# For simplification 
ethanol = F.ethanol
jet_fuel = F.jet_fuel
diesel = F.diesel
gasoline = F.gasoline


natural_gas = F.natural_gas
BT = F.BT
HXN = F.HXN




#%%
# To try the tea functions with state data

folder = os.path.dirname(__file__)
st_data_file = os.path.join(folder, 'state_scenarios_for_import_2023.xlsx') # Electricity and natural gas prices actualized to 2023
st_data = pd.read_excel(st_data_file, index_col=[0])

state = 'Alabama' # example

tea_SAF.state_income_tax = st_data.loc[state]['Income Tax Rate (decimal)']
tea_SAF.property_tax = st_data.loc[state]['Property Tax Rate (decimal)']
tea_SAF.fuel_tax = st_data.loc[state]['State Motor Fuel Tax (decimal)']
# tea_SAF.sales_tax = st_data.loc[state]['State Sales Tax Rate (decimal)']
bst.PowerUtility.price = st_data.loc[state]['Electricity Price (USD/kWh)']
tea_SAF.F_investment = st_data.loc[state]['Location Capital Cost Factor (dimensionless)']
price['natural gas'] = st_data.loc[state]['Natural Gas Price (USD/kg)']

GWP_CFs['electricity'] = st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)']






# Actualize for later
price['electricity'] = bst.PowerUtility.price
F.natural_gas.price = price['natural gas']

bst.stream_utility_prices['Natural gas'] = F.natural_gas.price = round(st_data.loc[state]['Natural Gas Price (USD/kg)'],3)


#%%
##### Functions to calculate all the metrics #####

get_annual_factor = lambda: tea_SAF.operating_days * 24
# 1. Product characteristics

get_ethanol_yield = lambda: ethanol.F_vol * _gal_per_m3 * get_annual_factor() / 1e6 # in MMGal (million gallon)
get_jet_yield = lambda:  jet_fuel.F_vol * _gal_per_m3 * get_annual_factor() / 1e6 # in MMGal (million gallon)
get_diesel_yield = lambda:  diesel.F_vol * _gal_per_m3 * get_annual_factor() / 1e6
get_gasoline_yield = lambda:  gasoline.F_vol * _gal_per_m3 * get_annual_factor() / 1e6
get_total_yield = lambda:  get_jet_yield() + get_diesel_yield() + get_gasoline_yield()
get_jet_vol_ratio = lambda:  get_jet_yield() / get_total_yield() * 100
get_jet_to_eth_ratio = lambda:  get_jet_yield() / get_ethanol_yield() * 100


get_ethanol_energy = lambda:  ethanol.LHV * _GGE_per_J * get_annual_factor() / 1000 # in MM GGE
get_jet_energy = lambda:  jet_fuel.LHV * _GGE_per_J * get_annual_factor() / 1000
get_diesel_energy = lambda:  diesel.LHV * _GGE_per_J * get_annual_factor() / 1000
get_gasoline_energy = lambda:  gasoline.LHV * _GGE_per_J * get_annual_factor() / 1000
get_total_energy = lambda:  get_jet_energy() + get_diesel_energy() + get_gasoline_energy()

get_jet_energy_ratio = lambda: get_jet_energy() / get_total_energy() * 100
get_diesel_energy_ratio = lambda: get_diesel_energy() / get_total_energy() * 100
get_gasoline_energy_ratio = lambda: get_gasoline_energy() / get_total_energy() * 100


# 2. TEA

# get_ethanol_transportation_cost = lambda: tea_SAF.eth_trans_cost / 1e6 # in a year

get_installed_cost = lambda: tea_SAF.installed_equipment_cost / 1e6
get_installed_cost_OSBL = lambda: sum(i.installed_cost for i in tea_SAF.OSBL_units) / 1e6
get_installed_cost_ISBL = lambda: get_installed_cost() - get_installed_cost_OSBL()
get_DPI = lambda: tea_SAF.DPI / 1e6
get_TDC = lambda: tea_SAF.TDC / 1e6
get_FCI = lambda: tea_SAF.FCI / 1e6
get_TCI = lambda: tea_SAF.TCI / 1e6
get_FOC = lambda: tea_SAF.FOC / 1e6
get_VOC = lambda: tea_SAF.VOC / 1e6
get_AOC = lambda: tea_SAF.AOC / 1e6 # Excludes electricity credit

get_MPSP = lambda: tea_SAF.solve_price(F.jet_fuel)  # $/kg

_jet_price_conversion_index_vol = lambda: _L_per_gal * jet_fuel.rho / 1000 # from $/kg to $/gal
_jet_price_conversion_index_energy = lambda: jet_fuel.F_mass / jet_fuel.LHV * 1000000 # from $/kg to $/GJ

get_MPSP_per_gallon = lambda: get_MPSP() * _jet_price_conversion_index_vol()

get_jet_price_per_GJ = lambda: get_MPSP() * _jet_price_conversion_index_energy()

kerosene_price_per_GJ = 2.95/0.1899  # price is $2.95/gal (2022-202405 average) from https://www.eia.gov/dnav/pet/hist/LeafHandler.ashx?n=PET&s=EER_EPJK_PF4_RGC_DPG&f=M; LHV of jet-kerosene is 5.670 million Btu per barrel = 0.189910053 GJ/gallon from https://www.eia.gov/totalenergy/data/monthly/pdf/sec12.pdf
    
get_NPV = lambda: tea_SAF.NPV
solve_IRR = lambda: tea_SAF.solve_IRR()

get_material_cost = lambda: tea_SAF.material_cost / 1e6 #  Excludes electricity credit but includes the money spent on ash disposal 
get_annual_sale = lambda: tea_SAF.sales / 1e6



# Electricity 
# BT.power_utility.production - BT.power_utility.consumption = -BT.power_utility.rate
# i.power_utility.rate = i.power_utility.power; meaning requirement
get_excess_power = lambda: -sum(i.power_utility.rate for i in sys.units) * sys.operating_hours
electricity_price = price['electricity']
get_electricity_credit = lambda: get_excess_power() * electricity_price / 1e6

jet_fuel_gal_per_hr = lambda: jet_fuel.F_vol * _gal_per_m3
jet_fuel_gal_per_year = lambda: jet_fuel_gal_per_hr() * sys.operating_hours

# In per gal of jet
get_cost_electricity_credit = lambda: get_excess_power() * electricity_price / jet_fuel_gal_per_year()


metrics = [Metric('Installed cost', get_installed_cost, '10^6 $'),
            Metric('Installed cost OSBL', get_installed_cost_OSBL, '10^6 $'),
            Metric('Installed cost ISBL', get_installed_cost_ISBL, '10^6 $'),
            Metric('DPI', get_DPI, '10^6 $'),
            Metric('TDC', get_TDC, '10^6 $'),
            Metric('FCI', get_FCI, '10^6 $'),
            Metric('TCI', get_TCI, '10^6 $'),
            Metric('FOC', get_FOC, '10^6 $/yr'),
            Metric('VOC', get_VOC, '10^6 $/yr'),
            Metric('AOC', get_AOC, '10^6 $/yr'),
           
           Metric('Minimum selling price', get_MPSP_per_gallon, '$/gal'),
           Metric('Jet volume yield', get_jet_yield, '10^6 Gal/yr'),
            Metric('Total volume yield', get_total_yield, '10^6 Gal/yr'),
            Metric('Jet volume ratio', get_jet_vol_ratio, '%'),
            Metric('Jet energy ratio', get_jet_energy_ratio, '%'),
           Metric('Jet volume ratio to ethanol', get_jet_to_eth_ratio, '%'),
           
            Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
            Metric('Annual product sale', get_annual_sale, '10^6 $/yr'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr'),
           Metric('Electricity credit to jet', get_cost_electricity_credit, '$/gal'),]


# 3. LCA

# in g CO2 eq / MJ blend fuel
_total_energy_per_hr = lambda: (jet_fuel.LHV + diesel.LHV + gasoline.LHV) / 1000 # from kJ to MJ
_total_energy_per_year = lambda: _total_energy_per_hr() * sys.operating_hours 

main_product = [jet_fuel]
coproducts = [diesel, gasoline]
#impurities = [CH4_C2H6] # not counted here

emissions = [i for i in F.stream if i.source and not i.sink and i not in main_product and i not in coproducts]


# Carbon balance
# total_C_in = sum([feed.get_atomic_flow('C') for feed in sys.feeds])
# total_C_out = sum([i.get_atomic_flow('C') for i in emissions]) + sum([i.get_atomic_flow('C') for i in main_product]) +\
#               sum([i.get_atomic_flow('C') for i in coproducts]) 
# C_bal_error = (total_C_out - total_C_in)/total_C_in


# Displacement GWP (gasoline, diesel)
get_GWP_displacement = lambda: (sys.get_material_impact(F.gasoline, key='GWP100') + sys.get_material_impact(F.diesel, key='GWP100')) * 1000 / _total_energy_per_year() 
                                  
get_GWP_emissions_non_BT = 0

# Non-biogenic missions (BT)
get_GWP_emissions_BT = lambda: (F.natural_gas.get_atomic_flow('C') + F.natural_gas.get_atomic_flow('C')) * SAF_chemicals.CO2.MW * 1000 / _total_energy_per_hr()

# Material contribution
get_GWP_material_total = lambda: sys.get_total_feeds_impact('GWP100') * 1000 / _total_energy_per_year()

get_GWP_NG = lambda: (sys.get_material_impact(F.natural_gas, key='GWP100') + sys.get_material_impact(F.natural_gas_for_h2, key='GWP100')) * 1000 / _total_energy_per_year()


get_GWP_other_materials = lambda: get_GWP_material_total() - get_GWP_NG()

get_GWP_ethanol_tranporation = lambda: sys.get_material_impact(F.ethanol, key='GWP100') * 1000 / _total_energy_per_year() # GHG of ethanol transportation

# Total = emission + material
#get_GWP_total = lambda: get_GWP_emissions_total() + get_GWP_material_total()
get_GWP_total = lambda: get_GWP_material_total() + get_GWP_emissions_BT() + get_GWP_displacement() #+ get_GWP_ethanol_tranporation() # ethanol is already in feeds

# Electricity (BT satisfies all electricity in system by buying natural gas if needed, no buying electricity)
get_electricity_use_offset_total = lambda: sum(i.power_utility.rate for i in sys.units) # .= 0 per hour

get_GWP_electricity_use_total = lambda: get_electricity_use_offset_total() * GWP_CFs['electricity'] * 1000 / _total_energy_per_hr() # . = 0

get_electricity_demand_total = get_electricity_use = lambda: -BT.power_utility.rate # Outside BT + BT consumption 

get_electricity_demand_cooling = lambda: F.CT.power_utility.rate

get_electricity_demand_non_cooling = lambda: get_electricity_demand_total() - get_electricity_demand_cooling()

get_electricity_demand_cooling_frac = lambda: get_electricity_demand_cooling() / get_electricity_demand_total()

get_electricity_demand_non_cooling_frac = lambda: get_electricity_demand_non_cooling() / get_electricity_demand_total()


# Steam contribution
get_steam_heating = lambda: sum([i.duty for i in BT.steam_utilities]) + sum([i.duty for i in F.BT.steam_utilities])# in kJ/hr

get_steam_electricity = lambda: 3600. * BT.electricity_demand / BT.turbogenerator_efficiency + 3600. * F.BT.electricity_demand / F.BT.turbogenerator_efficiency# in kJ/hr

get_steam_total = lambda: get_steam_heating() + get_steam_electricity()

get_steam_frac_heating = lambda: get_steam_heating() / get_steam_total()

get_steam_frac_electricity = lambda: get_steam_electricity() / get_steam_total()

get_steam_frac_cooling = lambda: get_steam_frac_electricity() * get_electricity_demand_cooling_frac()

get_steam_frac_non_cooling = lambda: get_steam_frac_electricity() * get_electricity_demand_non_cooling_frac()

# Steam GWP allocation

get_GWP_steam_total = lambda: get_GWP_NG() + get_GWP_emissions_BT()

get_GWP_heating_demand = lambda: get_steam_frac_heating() * get_GWP_steam_total()

get_GWP_cooling_demand = lambda: get_steam_frac_cooling() * get_GWP_steam_total() + \
                                  get_electricity_demand_cooling_frac() * get_GWP_electricity_use_total()
                                 
get_GWP_electricity_non_cooling = lambda: get_steam_frac_non_cooling() * get_GWP_steam_total() + \
                                          get_electricity_demand_non_cooling_frac() * get_GWP_electricity_use_total()

# GWP allocation (displacement method)
get_GWP_electricity_credit = lambda: -get_excess_power() * GWP_CFs['electricity'] * 1000 / _total_energy_per_year()

get_GWP_total_with_eletricity_credit = lambda: get_GWP_total() + get_GWP_electricity_credit()

get_GWP_jet = lambda: get_GWP_total_with_eletricity_credit()

get_GWP_jet = lambda: get_GWP_total_with_eletricity_credit() * get_jet_energy_ratio() / 100 # in g CO2-eq/MJ jet

get_GWP_diesel = lambda: get_GWP_total_with_eletricity_credit() * get_diesel_energy_ratio() / 100

get_GWP_gasoline = lambda: get_GWP_total_with_eletricity_credit() * get_gasoline_energy_ratio() / 100


# metrics.extend((Metric('GWP - ethanol transportation', get_GWP_ethanol_tranporation, 'g CO2-eq/MJ blend fuel', 'LCA'),)) 
metrics.extend((Metric('GWP - total', get_GWP_total, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - electricity credit', get_GWP_electricity_credit, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - total-electricity credit', get_GWP_total_with_eletricity_credit, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - NG', get_GWP_NG, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - other materials', get_GWP_other_materials, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - electricity', get_GWP_electricity_use_total, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - non biogenic emissions', get_GWP_emissions_BT, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - other non biogenic emissions', 0, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - heating_demand', get_GWP_heating_demand, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - cooling_demand', get_GWP_cooling_demand, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - electricity non cooling', get_GWP_electricity_non_cooling, 'g CO2-eq/MJ blend fuel', 'LCA'),))
metrics.extend((Metric('GWP - jet', get_GWP_jet, 'g CO2-eq/MJ jet fuel', 'LCA'),))
metrics.extend((Metric('GWP - displacement', get_GWP_displacement, 'g CO2-eq/MJ', 'LCA'),))


#%% CREATE MODEL WITH STATE SPECIFIC COSTS 


#%%
folder = os.path.dirname(__file__)
st_data_file = os.path.join(folder, 'state_scenarios_for_import_2023.xlsx') # Electricity and natural gas prices actualized to 2023
st_data = pd.read_excel(st_data_file, index_col=[0])
#%%
# Model for state specific analysis ===========================================
def create_states_model():
    
    tea = create_SAF_tea_withBlocs()
    sys.operating_hours = tea.operating_days * 24
    tea.federal_income_tax = 0.21 # corporate federal income tax
    tea.ethanol_feed = F.ethanol
    tea.jet_fuel_product = F.jet_fuel
    tea.jet_fuel_group = bst.UnitGroup('Jet fuel group', sys.units)
    if tea.jet_fuel_product == None:
        print('print: ', tea_SAF.jet_fuel_product)
    
    all_states = [
                'Alabama',
                'Arkansas',
                'California', # not rainfed
                'Colorado', # not rainfed but included for jet producers
                'Connecticut', # included now
                'Delaware',
                'Florida',
                'Georgia',
                'Illinois',
                'Indiana',
                'Iowa',
                'Kansas',
                'Kentucky',
                'Louisiana',
                'Maine', # included now
                'Maryland',
                'Massachusetts', # included now
                'Michigan',
                'Minnesota',
                'Mississippi',
                'Missouri',
                'Montana', # not rainfed but included for jet producers
                'New Hampshire', # included now
                'Nebraska',
                'New Jersey',
                'New Mexico', # not rainfed but included for jet producers
                'New York',
                'North Carolina',
                'North Dakota',
                'Ohio',
                'Oklahoma',
                'Oregon', # not rainfed
                'Pennsylvania',
                'Rhode Island', # included now
                'South Carolina',
                'South Dakota',
                'Tennessee',
                'Texas',
                'Utah', # not rainfed
                'Vermont', # included now
                'Virginia',
                'West Virginia',
                'Wisconsin',
                'Wyoming', # not rainfed but included for jet producers 
                ]
    states_w_inc = [
                    # 'Alabama', #
                    # 'Colorado',#
                    # 'Iowa', #
                    # 'Kansas',#
                    # 'Kentucky', #
                    # 'Louisiana', #
                    # 'Montana', #
                    # 'Nebraska', #
                    # 'New Mexico', #
                    # 'Oregon', #
                    # 'South Carolina', #
                    # 'Utah', #
                    # 'Virginia', #
                    ]
    model = bst.Model(sys, exception_hook='raise')
    #model = bst.Model(tea.system, exception_hook='raise')
    param = model.parameter
    bst.stream_utility_prices['Natural gas'] = F.natural_gas.price = 0.894
    bst.PowerUtility.price = 0.0779 # electricity price Delaware (using as baseline example)
    bst.CE = 798 # WENJUN'S VALUE AVERAGE 2023
    tea.fuel_tax = 0.
    tea.sales_tax = 0 # Delaware (using as baseline example) # 0.05785
    tea.federal_income_tax = 0.21
    tea.state_income_tax = 0.087 # Delaware (using as baseline example) # 0.065
    tea.property_tax = 0.003655 # Delaware (using as baseline example) # 0.0136
    tea.F_investment = 1.04588235294117 #Delaware (using as baseline example) # 1.02
    
    # tea.feedstock.price = feed_price
    
    # def get_state_incentives(state):
    # NOT USED
        
    # def solve_price():
    #     MFSP = float(tea.solve_price(F.jet_fuel) * _jet_price_conversion_index_vol())
    #     # try:
    #     #     MFSP = tea.solve_price(tea.ethanol_product)
    #     # except:
    #     #     MFSP = tea.solve_price([tea.ethanol_product], [4])
            
    #     # return 2.98668849 * MFSP # IN $/GAL
    #     return MFSP #$/gal
    
    # def solve_price(set_price=False): 
    #     MFSP = tea.solve_price(tea.jet_fuel_product)
    #     if set_price: tea.jet_fuel_product.price = MFSP
    #     return 2.98668849 * MFSP
    
    def MFSP_getter(state):
        def MFSP():
            names = (
                  'state_income_tax',
                  'property_tax',
                  'fuel_tax',
                  'sales_tax',
                  'F_investment',
                  'state_tax_by_gross_receipts',
                  'deduct_federal_income_tax_to_state_taxable_earnings',
                  'deduct_half_federal_income_tax_to_state_taxable_earnings',
                  'incentive_numbers',
              )
            # # original_feedstock_price = tea.feedstock.price
            original_electricity_price = bst.PowerUtility.price
            dct = {i: getattr(tea, i) for i in names}
            tea.state_income_tax = st_data.loc[state]['Income Tax Rate (decimal)']
            tea.property_tax = st_data.loc[state]['Property Tax Rate (decimal)']
            tea.fuel_tax = st_data.loc[state]['State Motor Fuel Tax (decimal)']
            # tea.sales_tax = st_data.loc[state]['State Sales Tax Rate (decimal)']
            bst.PowerUtility.price = st_data.loc[state]['Electricity Price (USD/kWh)']
            bst.stream_utility_prices['Natural gas'] = F.natural_gas.price = round(st_data.loc[state]['Natural Gas Price (USD/kg)'],3)
            tea.F_investment = st_data.loc[state]['Location Capital Cost Factor (dimensionless)']
            F.gasoline.price = price['gasoline'] = st_data.loc[state]['Gasoline Price (USD/kg)']
            F.diesel.price = price['diesel'] = st_data.loc[state]['Diesel Price (USD/kg)']
            
            # GWP_CFs['electricity'] = st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)']
            # bst.PowerUtility.characterization_factors['GWP100'] = (st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)'], st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)'])
            # tea.incentive_numbers = ()
            
            if state == 'Ohio' or state == 'Texas':
                tea.state_tax_by_gross_receipts = True
            else:
                tea.state_tax_by_gross_receipts = False

            if state == 'Alabama' or state == 'Louisiana':
                tea.deduct_federal_income_tax_to_state_taxable_earnings = True
            else:
                tea.deduct_federal_income_tax_to_state_taxable_earnings = False

            if state == 'Iowa' or state == 'Missouri':
                tea.deduct_half_federal_income_tax_to_state_taxable_earnings = True
            else:
                tea.deduct_half_federal_income_tax_to_state_taxable_earnings = False

            # MFSP = float(tea.solve_price(F.jet_fuel) * _jet_price_conversion_index_vol())
            MFSP = tea.solve_price(tea.jet_fuel_product) * _jet_price_conversion_index_vol()

            # MFSP = solve_price()
            # tea.feedstock.price = original_feedstock_price
            bst.PowerUtility.price = original_electricity_price
            for i in names: setattr(tea, i, dct[i])
            return MFSP
        return MFSP
    
    def GWP_getter(state):
        def GWP():
            
            bst.PowerUtility.characterization_factors['GWP100'] = (st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)'], st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)'])
            GWP_CFs['electricity'] = st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)']
            GWP = get_GWP_jet()
            
            return GWP
        return GWP
    
    # def MFSP_w_inc_getter(state): # NOT INCLUDED
    
    
    @model.metric(name='Jet fuel production', units='MMgal/yr')
    def get_Jet_production():
        return get_jet_yield()
        #return tea.ethanol_product.F_mass / 2.98668849 * tea.operating_hours # gal /yr
    
    @model.metric(name='Ethanol input', units='MMgal/yr')
    def get_ethanol_production():
        return get_ethanol_yield()
    
    @model.metric(name='Jet to Ethanol ratio', units='%')
    def get_jet_conv():
        return get_jet_to_eth_ratio()
    
    @model.metric(name='Annual factor', units='hr/yr')
    def annual_fac():
        return get_annual_factor()
    
    @model.metric(name='Ethanol price', units='$/kg')
    def eth_price():
        return F.ethanol.price
    
    @model.metric(name='Total capital investment', units='USD')
    def get_TCI():
        return tea.TCI
    
    @model.metric(name='Ethanol GWP', units='kgCO2e/kg')
    def get_ethanol_GWP():
        return F.ethanol.characterization_factors['GWP100']
    
    # @model.metric(name='Natural gas price', units='USD/kg')
    # def get_NGasPrice():
    #     return bst.stream_utility_prices['Natural gas']
    # @model.metric(name='Natural gas price2', units='USD/kg')
    # def get_NGasPrice2():
    #     return F.natural_gas.price


    
    @model.metric(name="Baseline MFSP", units='USD/gal') #within this function, set whatever parameter values you want to use as the baseline
    def MFSP_baseline():
        names = (
            'state_income_tax', 'property_tax', 'fuel_tax', 'sales_tax', 'F_investment',
            #'incentive_numbers',
        )
        dct_old = {i: getattr(tea, i) for i in names}
        tea.state_income_tax = 0.087 # Delaware example
        tea.property_tax = 0.003655 # Delaware example
        tea.fuel_tax = 0
        tea.sales_tax = 0
        tea.F_investment = 1.04588235294117 # Delaware example
        # tea.incentive_numbers = ()
        old_price = bst.PowerUtility.price
        bst.PowerUtility.price = 0.0779 # Delaware example
        bst.stream_utility_prices['Natural gas'] = F.natural_gas.price = 0.894 # Delaware
        F.gasoline.price = price['gasoline'] = 1.205 # Delaware
        F.diesel.price = price['diesel'] = 1.142 # Delaware
        # tea.feedstock.price = tea.feedstock.price
        MFSP = tea.solve_price(F.jet_fuel) * _jet_price_conversion_index_vol()
        # MFSP = 2.98668849 * tea.solve_price(tea.ethanol_product) 
        for i, j in dct_old.items(): setattr(tea, i, j)
        bst.PowerUtility.price = old_price
        return MFSP
    
    
    get_inc_value = lambda: tea.exemptions.sum() + tea.deductions.sum() + tea.credits.sum()+ tea.refunds.sum()

    for state in all_states:
        model.metric(MFSP_getter(state), 'MFSP', 'USD/gal', state)
        model.metric(GWP_getter(state), 'GWP Jet', 'g CO2-eq/MJ jet fuel', state)


    
    # ============================================================================
    # TEA parameters
    # ============================================================================  
    ##### Financial parameters #####
    
    baseline_flow =  1960.7 # kg/hr
    D = shape.Triangle(baseline_flow*0.9, baseline_flow,baseline_flow*1.1)
    @param(name='Ethanol flow', element='ethanol', kind='coupled', units='kg/hr',
            baseline=baseline_flow, distribution=D)
    def set_ethanol_flow(flow):
        ethanol.F_mass = flow
        
        
    # For now, keep constant 
    # D = shape.Uniform(0.8,1.6)
    # @param(name='Ethanol price', element='ethanol', kind='isolated', units='$/kg',
    #         baseline=1.25, distribution=D)
    # def set_ethanol_price(price):
    #     ethanol.price = price
        
    D = shape.Triangle(0.84, 0.9, 0.96)
    @param(name='Plant uptime', element='TEA', kind='isolated', units='%',
           baseline=0.9, distribution=D)
    def set_operating_days(uptime):
        tea_SAF.operating_days = 365. * uptime



    D = shape.Triangle(0.75, 1, 1.25)
    @param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
            baseline=1, distribution=D)
    def set_TCI_ratio(new_ratio):
        old_ratio = tea_SAF._TCI_ratio_cached
        for unit in sys.units:
            if hasattr(unit, 'cost_items'):
                for item in unit.cost_items:
                    unit.cost_items[item].cost /= old_ratio
                    unit.cost_items[item].cost *= new_ratio
        tea_SAF._TCI_ratio_cached = new_ratio


    # state specific
    ##### Material price #####
    #D = shape.Triangle(0.0583, 0.0637, 0.069) # 
    baseline_price=bst.PowerUtility.price
    D = shape.Triangle(baseline_price*0.9, baseline_price, baseline_price*1.1)
    @param(name='Electricity price', element='Electricity', kind='isolated', units='$/kWh',
           baseline=baseline_price, distribution=D)
    def set_electricity_price(price): 
        bst.PowerUtility.price = price


    # state specific
    # F.natural_gas.price = st_data.loc[state]['Natural Gas Price (USD/kg)']
    natural_gas_price = F.natural_gas.price
    D = shape.Triangle(natural_gas_price*0.9, natural_gas_price, natural_gas_price*1.1)
    @param(name='Natural gas price', element='Natural gas', kind='isolated', units='$/kg',
            baseline=natural_gas_price, distribution=D)
    def set_natural_gas_price(price): 
        F.natural_gas.price = price
        F.natural_gas_for_h2.price = price
        bst.stream_utility_prices['Natural gas'] = price
               


    ash_disposal_price = price['ash disposal']
    D = shape.Triangle(ash_disposal_price*1.5, ash_disposal_price, ash_disposal_price*0.5)
    @param(name='Ash disposal price', element='Ash disposal', kind='isolated', units='$/kg',
            baseline=ash_disposal_price, distribution=D)
    def set_ash_disposal_price(price):
        F.ash.price = price



    Syndol_catalyst_price = price['Syndol catalyst']
    D = shape.Triangle(Syndol_catalyst_price*0.7, Syndol_catalyst_price, Syndol_catalyst_price*1.3)
    @param(name='Syndol catalyst price', element='Syndol catalyst', kind='isolated', units='$/kg',
            baseline=Syndol_catalyst_price, distribution=D)
    def set_Syndol_catalyst_price(price):
        F.R401.catalyst_price = price



    Oligomerization1_catalyst_price = price['Ni-loaded aluminosilicate catalyst']
    D = shape.Triangle(Oligomerization1_catalyst_price*0.7, Oligomerization1_catalyst_price, Oligomerization1_catalyst_price*1.3)
    @param(name='Oligomerization1 catalyst price', element='Oligomerization1 catalyst', kind='isolated', units='$/kg',
            baseline=Oligomerization1_catalyst_price, distribution=D)
    def set_oligomerization1_catalyst_price(price):
        F.R402.catalyst_price = price



    Oligomerization2_catalyst_price = price['Aluminosilicate catalyst']
    D = shape.Triangle(Oligomerization2_catalyst_price*0.7, Oligomerization2_catalyst_price, Oligomerization2_catalyst_price*1.3)
    @param(name='Oligomerization2 catalyst price', element='Oligomerization2 catalyst', kind='isolated', units='$/kg',
            baseline=Oligomerization2_catalyst_price, distribution=D)
    def set_oligomerization2_catalyst_price(price):
        F.R403.catalyst_price = price



    Hydrogenation_catalyst_price = price['Como catalyst']
    D = shape.Triangle(Hydrogenation_catalyst_price*0.5, Hydrogenation_catalyst_price, Hydrogenation_catalyst_price*1.5)
    @param(name='Hydrogenation catalyst price', element='Hydrogenation catalyst', kind='isolated', units='$/kg',
            baseline=Hydrogenation_catalyst_price, distribution=D)
    def set_hydrogenation_catalyst_price(price):
        F.R404.catalyst_price = price



    H2_price = price['h2']
    D = shape.Triangle(H2_price*0.7, H2_price, H2_price*1.3)
    @param(name='H2 price', element='H2', kind='isolated', units='$/kg',
            baseline=H2_price, distribution=D)
    def set_hydrogen_price(price):
        F.hydrogen.price = price

     
     
    diesel_price = price['diesel']
    D = shape.Triangle(diesel_price*0.7, diesel_price, diesel_price*1.5) # Fluctuation maximum rate=0.4
    @param(name='Diesel price', element='Diesel', kind='isolated', units='$/kg',
            baseline=diesel_price, distribution=D)
    def set_diesel_price(price):
        diesel.price = price
        
        

    gasoline_price = price['gasoline']
    D = shape.Triangle(gasoline_price*0.8, gasoline_price, gasoline_price*1.2) # Fluctuation maximum rate=0.2 (Last 10 year average variation is +-15%)
    @param(name='Gasoline price', element='Gasoline', kind='isolated', units='$/kg',
            baseline=gasoline_price, distribution=D)
    def set_gasoline_price(price):
        gasoline.price = price



    ##### Upgrading parameters #####
    R401 = F.R401
    R402 = F.R402
    R403 = F.R403
    R404 = F.R404

    D = shape.Triangle(0.9, .988, .988)# (0.92, .988)# (0.995*0.988*0.9, 0.995*0.988) 
    @param(name='Dehydration ethanol-to-ethylene', element=R401, kind='coupled', units='%',
            baseline=0.988, distribution=D)#.995*0.988
    def set_R401_ethanol_conversion(X):
        R401.dehydration_rxns[0].X = X
        # print(R401.dehydration_rxns[0].X)



    # D = shape.Triangle(0.25, 0.38, 1.8)#(0.3, 0.43, 2)
    # @param(name='Dehydration WHSV', element=R401, kind='coupled', units='h^-1',
    #         baseline=0.43, distribution=D)
    # def set_R401_WHSV(X):
    #     R401.WHSV = X
    #     # print('WHSV: ', X)



    # D = shape.Uniform(3.5/3600, 4.5/3600)#(3.14/3600, 4.05/3600)
    # @param(name='Dehydration residence time', element=R401, kind='coupled', units='hr',
    #         baseline=3.14/3600, distribution=D)
    # def set_R401_residence_time(X):
    #     R401.tau = X
    #     # print('tau: ', X)


    # D = shape.Uniform(7884*0.85, 7884)#(7884*0.8,7884)
    # @param(name='Dehydration catalyst longevity', element=R401, kind='coupled', units='hr',
    #         baseline=7884, distribution=D)
    # def set_R401_catalyst_longevity(t):
    #     R401.catalyst_lifetime = t
    #     # print('catalyst_lifetime: ', t)



    D = shape.Uniform(0.988*0.97, 0.988)
    @param(name='1st oligomerization ethylene-to-C4H8', element=R402, kind='coupled', units='%',
            baseline=0.988, distribution=D)
    def set_R402_ethylene_conversion(X):
        R402.oligomerization_rxns[0].X = X
    #     # print(X)
        
        
    D = shape.Uniform(0.5, 5) 
    @param(name='1st oligomerization WHSV', element=R402, kind='coupled', units='h^-1',
            baseline=5, distribution=D)
    def set_R402_WHSV(X):
        R402.WHSX = X



    D = shape.Triangle(48*0.8, 48, 48*1.2) 
    @param(name='1st oligomerization residence time', element=R402, kind='coupled', units='hr',
            baseline=48, distribution=D)
    def set_R402_residence_time(X):
        R402.tau = X


    D = shape.Uniform(7884*0.8,7884)
    @param(name='1st oligomerization catalyst longevity', element=R402, kind='coupled', units='hr',
            baseline=7884, distribution=D)
    def set_R402_catalyst_longevity(t):
        R402.catalyst_lifetime = t



    D = shape.Uniform(0.5, 10) 
    @param(name='2nd oligomerization WHSV', element=R403, kind='coupled', units='h^-1',
            baseline=10, distribution=D)
    def set_R403_WHSV(X):
        R403.WHSX = X



    D = shape.Triangle(48*0.8, 48, 48*1.2) 
    @param(name='2nd oligomerization residence time', element=R403, kind='coupled', units='hr',
            baseline=48, distribution=D)
    def set_R403_residence_time(X):
        R403.tau = X



    D = shape.Uniform(7884*0.8,7884)
    @param(name='2nd oligomerization catalyst longevity', element=R403, kind='coupled', units='hr',
            baseline=7884, distribution=D)
    def set_R403_catalyst_longevity(t):
        R403.catalyst_lifetime = t



    D = shape.Uniform(1, 3)
    @param(name='Hydrogenation WHSV', element=R404, kind='coupled', units='h^-1',
            baseline=3, distribution=D)
    def set_R404_WHSV(X):
        R404.WHSX = X



    D = shape.Triangle(1*0.8, 1, 1*1.2) 
    @param(name='Hydrogenation residence time', element=R404, kind='coupled', units='hr',
            baseline=1, distribution=D)
    def set_R404_residence_time(X):
        R404.tau = X



    D = shape.Uniform(7884*0.8,7884)
    @param(name='Hydrogenation catalyst longevity', element=R404, kind='coupled', units='hr',
            baseline=7884, distribution=D)
    def set_R404_catalyst_longevity(t):
        R404.catalyst_lifetime = t



    ##### Facilities parameter #####
    D = shape.Uniform(0.8*(1-0.1), 0.8*(1+0.1))
    @param(name='Boiler efficiency', element=BT, kind='coupled', units='',
            baseline=0.8, distribution=D)
    def set_boiler_efficiency(efficiency):
        F.BT.boiler_efficiency = efficiency



    # =============================================================================
    # LCA parameters
    # =============================================================================
    D = shape.Uniform(2.11*(1-0.5), 2.11*(1+0.5))
    @param(name='NaOH GWP', element='NaOH', kind='isolated', units='kg CO2-eq/kg',
            baseline=2.11, distribution=D)
    def set_NaOH_GWP(X):
        F.NaOH.characterization_factors['GWP100'] = X



    D = shape.Uniform(0.4*(1-0.5), 0.4*(1+0.5))
    @param(name='Natural gas GWP', element='Natural gas', kind='isolated', units='kg CO2-eq/kWh',
            baseline=0.4, distribution=D)
    def set_natural_gas_GWP(X):
        F.natural_gas.characterization_factors['GWP100'] = X
        F.natural_gas_for_h2.characterization_factors['GWP100'] = X

    
    D = shape.Uniform(2.01*(1-0.5), 2.01*(1+0.5))
    @model.parameter(name='caustic', element='caustic', kind='isolated', units='kg CO2-eq/kg',
                     baseline=2.01, distribution=D)
    def set_caustic_GWP(X):
        F.caustic.characterization_factors['GWP100'] = X 
        
    
    D = shape.Uniform(0.01*(1-0.5), 0.01*(1+0.5))
    @model.parameter(name='ash', element='ash', kind='isolated', units='kg CO2-eq/kg',
                     baseline=0.01, distribution=D)
    def set_ash_GWP(X):
        F.ash.characterization_factors['GWP100'] = X
        
    D = shape.Uniform(1.5568*(1-0.5), 1.5568*(1+0.5))
    @model.parameter(name='boiler_chemicals', element='boiler_chemicals', kind='isolated', units='kg CO2-eq/kg',
                     baseline=1.5568, distribution=D)
    def set_boiler_chemicals_GWP(X):
        F.boiler_chems.characterization_factors['GWP100'] = X
        
    base_value = GWP_CFs['electricity']
    D = shape.Uniform(base_value*(1-0.2), base_value*(1+0.2))
    @model.parameter(name='electricity', element='electricity', kind='isolated', units='kg CO2-eq/kWh',
                     baseline=base_value, distribution=D)
    def set_electricity_GWP(X):
        bst.PowerUtility.characterization_factors['GWP100'] = (X, X)
        GWP_CFs['electricity'] = X
        
    D = shape.Uniform(1.5624*(1-0.4), 1.5624*(1+0.4))
    @model.parameter(name='hydrogen', element='hydrogen', kind='isolated', units='kg CO2-eq/kg',
                     baseline=1.5624, distribution=D)
    def set_hydrogen_GWP(X):
        F.hydrogen.characterization_factors['GWP100'] = X
        
    D = shape.Uniform(11*(1-0.5), 11*(1+0.5))
    @model.parameter(name='Syndol_catalyst', element='Syndol_catalyst', kind='isolated', units='kg CO2-eq/kg',
                     baseline=11, distribution=D)
    def set_Syndol_catalyst_GWP(X):
        F.Syndol_catalyst.characterization_factors['GWP100'] = X
    
    D = shape.Uniform(11*(1-0.5), 11*(1+0.5))
    @model.parameter(name='first_catalyst', element='first_catalyst', kind='isolated', units='kg CO2-eq/kg',
                     baseline=11, distribution=D)
    def set_first_catalyst_GWP(X):
        F.first_catalyst.characterization_factors['GWP100'] = X
    
    D = shape.Uniform(11*(1-0.5), 11*(1+0.5))
    @model.parameter(name='second_catalyst', element='second_catalyst', kind='isolated', units='kg CO2-eq/kg',
                     baseline=11, distribution=D)
    def set_second_catalyst_GWP(X):
        F.second_catalyst.characterization_factors['GWP100'] = X
        
    D = shape.Uniform(10.52*(1-0.5), 10.52*(1+0.5))
    @model.parameter(name='Como_catalyst', element='Como_catalyst', kind='isolated', units='kg CO2-eq/kg',
                     baseline=10.52, distribution=D)
    def set_Como_catalyst_GWP(X):
        F.Como_catalyst.characterization_factors['GWP100'] = X
        
    D = shape.Uniform(-0.8415*(1+0.5), -0.8415*(1-0.5))
    @model.parameter(name='gasoline', element='gasoline', kind='isolated', units='kg CO2-eq/kg',
                     baseline=-0.8415, distribution=D)
    def set_gasoline_GWP(X):
        F.gasoline.characterization_factors['GWP100'] = X
     
    D = shape.Uniform(-0.6535*(1+0.5), -0.6535*(1-0.5))
    @model.parameter(name='diesel', element='diesel', kind='isolated', units='kg CO2-eq/kg',
                     baseline=-0.6535, distribution=D)
    def set_diesel_GWP(X):
        F.diesel.characterization_factors['GWP100'] = X
        
    D = shape.Uniform(1.29*(1-0.5), 1.29*(1+0.5))
    @model.parameter(name='lime_boiler', element='lime_boiler', kind='isolated', units='kg CO2-eq/kg',
                     baseline=1.29, distribution=D)
    def set_lime_boiler_GWP(X):
        F.lime_boiler.characterization_factors['GWP100'] = X
    
   
    
    
    # # =============================================================================
    # # Ethanol tranportation parameters
    # # =============================================================================
    # D = shape.Uniform(0.009*0.8, 0.009*1.2)
    # @param(name='Ethanol transportation price', element='Ethanol', kind='isolated', units='$/kg',
    #        baseline=0.009, distribution=D)
    # def set_ethanol_transportation_price(X):
    #     tea_SAF.eth_trans_cost = X
    
    
    
    # D = shape.Uniform(0.0036*0.8, 0.0036*1.2)
    # @param(name='Ethanol transportation GWP', element='Ethanol', kind='isolated', units='kg CO2e/kg ethanol',
    #        baseline=0.0036, distribution=D)
    # def set_ethanol_transportation_GWP(X):
    #     F.ethanol.characterization_factors['GWP100'] = X
    
    return model

model = create_states_model()
system=model._system
def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
def run_bugfix_barrage():
    try:
        reset_and_reload()
    except Exception as e:
        print(str(e))
        try:
            reset_and_switch_solver('fixedpoint')
        except Exception as e:
            print(str(e))
            try:
                reset_and_switch_solver('aitken')
            except Exception as e:
                print(str(e))
                # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                print("Bugfix barrage failed.\n")
                # breakpoint()
                raise e
def model_specification():
    try:
        system.simulate()
    except Exception as e:
        str_e = str(e).lower()
        print('Error in model spec: %s'%str_e)
        run_bugfix_barrage()
        
model._specification = model_specification

def evaluate_SAF(N=1000, notify_runs=0, model = model): 
    start_time = time.time()
    np.random.seed(1237) 
    rule = 'L' # For Latin-Hypercube sampling
    samples = model.sample(N, rule)
    # samples = np.load('samples_jet.npy')
    # print(samples)

    # filtered_samples_45 = np.delete(samples, 398, axis=0)
    model.load_samples(samples)
    model.evaluate(notify = notify_runs) #**evaluate_args('SS'),notify=notify_runs
    model.show()
    model.get_baseline_sample()
    model.table.to_excel('model_table_Jet.xlsx')
    df_rho,df_p = model.spearman_r()
    df_rho.to_excel('df_rho_Jet.xlsx')
    df_p.to_excel('df_p_Jet.xlsx')
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")
    return model

# results_folder = os.path.join(folder, 'results')

# def get_file_name(name):
#     return os.path.join(results_folder, name)

# def evaluate_args(name, nbox=None):
#     if nbox is None:
#         return {'autoload': True,
#                 'autosave': 20,
#                 'file': get_file_name(name)}
#     else:
#         nbox[0] += 1
#         return {'autoload': True,
#                 'autosave': 20,
#                 'file': get_file_name(os.path.join(name, str(nbox[0])))}

