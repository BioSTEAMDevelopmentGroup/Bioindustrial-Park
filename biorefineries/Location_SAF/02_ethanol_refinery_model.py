#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:06:22 2024

@author: bianco3

Feedstock must be defined in line 26

Second model to run, from the results of this model we created the files:
    - ethanol_price_for_python
    - Ethanol_GWP_for_python
These files will be used to couple the uncertainty in feedstock transport model with uncertainty in ethanol results.

Next file to run: 03_combine_feedstock_and_ethanol
"""
#%% RUN Model

# change feedstock price -> tea.feedstock.price = 0.2 # I used 0, 1, 0.098 in different runs to test linearity
# change feedstock CI -> F.switchgrass.set_CF('GWP100', -0.19) # I used 0, 0.13 and  -0.19 in different runs to test linearity

# evaluate_SS('cellulosic', N=1000, notify_runs = 0) # refinery name is 'cellulosic' but feedstock should be changed below:
# N is the number of samples for uncertainty
    
# define feedstock 
crop = 'switchgrass'

#%%

import biosteam as bst
import thermosteam
from biorefineries.cellulosic.chemicals import create_cellulosic_ethanol_chemicals
from biorefineries.cellulosic.systems import create_cellulosic_ethanol_system
#from biorefineries.SAF.ATJ._process_settings import price, GWP_CFs, load_preferences_and_process_settings
from _process_settings import price, GWP_CFs, load_preferences_and_process_settings
import blocs as blc
import os
import pandas as pd
import numpy as np
from chaospy import distributions as shape
import time

_gal_per_m3 = 1000/3.785412
_L_per_gal = 3.785412

F = bst.Flowsheet('switchgrass_flor')
bst.main_flowsheet.set_flowsheet(F)

# #create chemicals
chem = create_cellulosic_ethanol_chemicals() 

chem.set_synonym('Extract','Extractives') 
bst.settings.set_thermo(chem, cache= True)


#create streams

if crop == 'switchgrass':
# # SWITCHGRASS
    swg = bst.Stream('switchgrass',
    Arabinan=0.02789023841655421,
    Galactan=0.010436347278452543,
    Glucan=0.2717049032838507,
    Xylan=0.21214574898785432,
    Mannan=0.005937921727395412,
    Lignin=0.17112010796221322,
    Ash=0.016194331983805668,
    Extractives=0.08457040035987407,
    Water=0.2,
    total_flow=127538.3908, # from feedstock transport model 80MMGal ethanol refinery
    units='kg/hr',
    price=0.08572288) # placeholder for price $/kg
    
elif crop == 'miscanthus':
    #MISCANTHUS
    swg = bst.Stream('switchgrass', #named switchgrass but using miscanthus composition
    Water=0.2,
    Ash=0.0301,
    Glucan=0.3472,
    Xylan=0.2048,
    Galactan=0.0112,
    Arabinan=0.0136,
    Lignin=0.1616,
    Extract=0.0223,
    Acetate=0.0034,
    Protein=0.0058,
    units='kg/hr',
    total_flow= 136545.25, #109236.2081, # change value if you want different refinery size
    price=0.08)
    


#create system
sys_ethanol = create_cellulosic_ethanol_system('sys_switchgrass',ins = swg)

bst.rename_units([i for i in F.unit if i.ID[1] == '6'], 700)

# F.switchgrass # to access stream
# F.ethanol # to access stream

# to simplify
ethanol = F.ethanol
feedstock = F.switchgrass

BT = F.BT

# 
# get_ethanol_conversion = lambda: ethanol.F_vol * _gal_per_m3 / (F.switchgrass.F_mass * 0.8/1000) #Gal per dry ton of biomass


#%%
load_preferences_and_process_settings(T='K',
                                      flow_units='kg/hr',
                                      N=100,
                                      P_units='Pa',
                                      CE=798, # Average 2023 https://toweringskills.com/financial-analysis/cost-indices/
                                      indicator='GWP100',
                                      electricity_EI=GWP_CFs['electricity'],
                                      electricity_price=price['electricity'])
sys= sys_ethanol
sys.set_tolerance(rmol=1e-6, mol=1e-5, maxiter=400)


#%%
from biorefineries.tea import cellulosic_ethanol_tea

__all__ = [*cellulosic_ethanol_tea.__all__, 'create_tea']

from biorefineries.tea.cellulosic_ethanol_tea import *

create_tea = create_cellulosic_ethanol_tea
#%%
def create_tea_Ethanol():
    tea = create_tea(sys=sys, cls=blc.incentives.incentives_tea.CellulosicIncentivesTEA) 
    tea.incentive_numbers = () # Empty for now
    tea.fuel_tax = 0.
    tea.sales_tax = 0.
    tea.federal_income_tax = 0.21
    tea.state_income_tax = 0.065
    tea.property_tax = 0.013
    tea.utility_tax = 0.
    tea.ethanol_product = F.ethanol
    tea.ethanol_group = bst.UnitGroup('Ethanol group', F.unit) # Assume all unit operations qualify
    tea.BT = F.BT
    tea.feedstock = F.switchgrass
    tea.cellulase = F.cellulase
    tea.cellulase.price = 0.258
    tea.R201 = F.R201
    tea.R303 = F.R303
    # tea.property_tax = 0.001 # repeated
    tea.sulfuric_acid = F.sulfuric_acid
    tea.sulfuric_acid.price = 0.1070
    tea.FGD_lime = F.FGD_lime
    tea.FGD_lime.price = 0.2958
    tea.denaturant = F.denaturant
    tea.denaturant.price = 0.496
    tea.cooling_tower_chemicals = F.cooling_tower_chemicals
    tea.cooling_tower_chemicals.price = 2.1447049186042744 # Wenjun's number # 4.4385 Dalton's number
    tea.makeup_water = F.makeup_water
    tea.makeup_water.price = 0.0005
    tea.CSL = F.CSL
    tea.CSL.price = 0.0843
    tea.DAP = F.DAP
    tea.DAP.price = 0.4092
    tea.ammonia = F.ammonia
    tea.ammonia.price = 0.4727
    tea.boiler_chemicals = F.boiler_chemicals
    tea.boiler_chemicals.price = 3.5635187930642878 # Wenjun's number # 7.4062 Dalton's number
    tea.caustic = F.caustic
    tea.caustic.price = 0.5931
    tea.natural_gas = F.natural_gas
    
    return tea

tea = create_tea_Ethanol()

#%%
sys.operating_hours = tea.operating_days * 24


def set_price_of_streams():
    for i in sys.streams:
        if i.ID in price.keys():
            i.price = price[i.ID]
    #F.natural_gas.price = price['natural gas']
            
def set_GWP_of_streams(indicator):
    F.caustic.set_CF(key='GWP100', value=GWP_CFs['caustic']) # caustic in WWT
    for i in sys.streams:
        if i.ID in GWP_CFs.keys():                
            i.characterization_factors[indicator]= GWP_CFs[i.ID]  
            
# GWP in kgCO2e/kg

set_prices = set_price_of_streams()
set_GWP = set_GWP_of_streams(indicator='GWP100')

#VERIFICAR SI FUNCIONA ######

#%%
#CORRRER FUNCION PARA DEFINIR POR ESTADO LOS SIGUIENTES VALORES
folder = os.path.dirname(__file__)
st_data_file = os.path.join(folder, 'state_scenarios_for_import_2023.xlsx') # Electricity and natural gas prices actualized to 2023
st_data = pd.read_excel(st_data_file, index_col=[0])

state = 'Alabama' # example

tea.state_income_tax = st_data.loc[state]['Income Tax Rate (decimal)']
tea.property_tax = st_data.loc[state]['Property Tax Rate (decimal)']
tea.fuel_tax = st_data.loc[state]['State Motor Fuel Tax (decimal)']
# tea.sales_tax = st_data.loc[state]['State Sales Tax Rate (decimal)']
bst.PowerUtility.price = st_data.loc[state]['Electricity Price (USD/kWh)']
tea.F_investment = st_data.loc[state]['Location Capital Cost Factor (dimensionless)']
price['natural gas'] = st_data.loc[state]['Natural Gas Price (USD/kg)']

GWP_CFs['electricity'] = st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)']




# Actualize for later
price['electricity'] = bst.PowerUtility.price
F.natural_gas.price = price['natural gas']

#%%
from biosteam.evaluation import Model, Metric

##### Functions to calculate all the metrics #####

get_annual_factor = lambda: tea.operating_days * 24 # same as tea.operating_hours

# 1. Product characteristics

get_ethanol_conversion = lambda: ethanol.F_vol * _gal_per_m3 / (F.switchgrass.F_mass * 0.8/1000) # in Gal ethanol per dry ton of feedstock (20% MC)
get_ethanol_yield = lambda: ethanol.F_vol * _gal_per_m3 * get_annual_factor() / 1e6 # in MMGal (million gallon)


# 2. TEA
get_MPSP = lambda: tea.solve_price(F.ethanol)
get_NPV = lambda: tea.NPV
solve_IRR = lambda: tea.solve_IRR()

_ethanol_price_conversion_index_vol = lambda: _L_per_gal * ethanol.rho / 1000 # from $/kg to $/gal
# _ethanol_price_conversion_index_energy = lambda: ethanol.F_mass / ethanol.LHV * 1000000 # from $/kg to $/GJ

get_MPSP_per_gallon = lambda: get_MPSP() * _ethanol_price_conversion_index_vol()


get_installed_cost = lambda: tea.installed_equipment_cost / 1e6
get_installed_cost_OSBL = lambda: sum(i.installed_cost for i in tea.OSBL_units) / 1e6
get_installed_cost_ISBL = lambda: get_installed_cost() - get_installed_cost_OSBL()
get_DPI = lambda: tea.DPI / 1e6
get_TDC = lambda: tea.TDC / 1e6
get_FCI = lambda: tea.FCI / 1e6
get_TCI = lambda: tea.TCI / 1e6
get_FOC = lambda: tea.FOC / 1e6
get_VOC = lambda: tea.VOC / 1e6
get_AOC = lambda: tea.AOC / 1e6 

get_material_cost = lambda: tea.material_cost / 1e6 #  Excludes electricity credit but includes the money spent on ash disposal 
get_annual_sale = lambda: tea.sales / 1e6


# Electricity 
# BT.power_utility.production - BT.power_utility.consumption = -BT.power_utility.rate
# i.power_utility.rate = i.power_utility.power; meaning requirement
get_excess_power = lambda: -sum(i.power_utility.rate for i in sys.units) * sys.operating_hours
electricity_price = price['electricity'] 
get_electricity_credit = lambda: get_excess_power() * electricity_price / 1e6

ethanol_gal_per_hr = lambda: ethanol.F_vol * _gal_per_m3
ethanol_gal_per_year = lambda: ethanol_gal_per_hr() * sys.operating_hours

# In per gal of ethanol
get_cost_electricity_credit = lambda: get_excess_power() * electricity_price / ethanol_gal_per_year()



metrics = [
            Metric('Installed cost', get_installed_cost, '10^6 $'),
            Metric('Installed cost OSBL', get_installed_cost_OSBL, '10^6 $'),
            Metric('Installed cost ISBL', get_installed_cost_ISBL, '10^6 $'),
            Metric('DPI', get_DPI, '10^6 $'),
            Metric('TDC', get_TDC, '10^6 $'),
            Metric('FCI', get_FCI, '10^6 $'),
            Metric('TCI', get_TCI, '10^6 $'),
            Metric('FOC', get_FOC, '10^6 $/yr'),
            Metric('VOC', get_VOC, '10^6 $/yr'),
            Metric('AOC', get_AOC, '10^6 $/yr'),
           Metric('NPV', get_NPV, ''), # change unit???
           Metric('Minimum selling price', get_MPSP_per_gallon, '$/gal'),
           Metric('Ethanol yield', get_ethanol_yield, '10^6 Gal/yr'),
           Metric('Ethanol conversion', get_ethanol_conversion, 'Gal/dry ton'),
           
            Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
            Metric('Annual product sale', get_annual_sale, '10^6 $/yr'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr'),
           Metric('Electricity credit to ethanol', get_cost_electricity_credit, '$/gal')
           ]


# 3. LCA

# in g CO2 eq / MJ fuel
_total_energy_per_hr = lambda: (ethanol.LHV) / 1000 # from kJ to MJ
_total_energy_per_year = lambda: _total_energy_per_hr() * sys.operating_hours 

main_product = [ethanol]
coproducts = []
# coproducts = [diesel, gasoline] # No coproducts
#impurities = [CH4_C2H6] # not counted here

emissions = [i for i in F.stream if i.source and not i.sink and i not in main_product and i not in coproducts]

# Carbon balance
# total_C_in = sum([feed.get_atomic_flow('C') for feed in sys.feeds])
# total_C_out = sum([i.get_atomic_flow('C') for i in emissions]) + sum([i.get_atomic_flow('C') for i in main_product]) +\
#               sum([i.get_atomic_flow('C') for i in coproducts]) 
# C_bal_error = (total_C_out - total_C_in)/total_C_in

# Feedstock contribution
get_GWP_feedstock_input = lambda: sys.get_material_impact(feedstock, key='GWP100') * 1000 / _total_energy_per_year()

# Non-biogenic emissions (ash)
get_GWP_missions_waste = lambda: (sys.get_material_impact(F.ash, key='GWP100') 
                                  #+ sys.get_material_impact(F.ash_1, key='GWP100')
                                  ) * 1000 / _total_energy_per_year() 

# Non-biogenic emissions (enzyme+CSL)
get_GWP_emissions_C_source = lambda: (F.CSL.get_atomic_flow('C') + F.cellulase.get_atomic_flow('C'))\
                                  * chem.CO2.MW * 1000 / _total_energy_per_hr()         
                                  
# Displacement GWP (ethnol_to_sold displacement, gasoline, diesel) # WE DON'T HAVE THIS HERE
# get_GWP_displacement = lambda: (sys.get_material_impact(F.gasoline, key='GWP100') + sys.get_material_impact(F.diesel, key='GWP100')) * 1000 / _total_energy_per_year() 
                                  
get_GWP_emissions_non_BT = lambda: get_GWP_missions_waste() + get_GWP_emissions_C_source()

# Non-biogenic emissions (BT)
get_GWP_emissions_BT = lambda: (F.natural_gas.get_atomic_flow('C') 
                                #+ F.natural_gas_1.get_atomic_flow('C')
                                ) * chem.CO2.MW * 1000 / _total_energy_per_hr()

# Material contribution
get_GWP_material_total = lambda: sys.get_total_feeds_impact('GWP100') * 1000 / _total_energy_per_year()

get_GWP_NG = lambda: (sys.get_material_impact(F.natural_gas, key='GWP100') + sys.get_material_impact(F.natural_gas_1, key='GWP100')) * 1000 / _total_energy_per_year()


get_GWP_other_materials = lambda: get_GWP_material_total()  - get_GWP_feedstock_input() - get_GWP_NG()

# Total = emission + material
#get_GWP_total = lambda: get_GWP_emissions_total() + get_GWP_material_total()
get_GWP_total = lambda: get_GWP_material_total() + get_GWP_emissions_BT() + get_GWP_emissions_non_BT() #+ get_GWP_displacement() # No coproducts

# Electricity (BT satisfies all electricity in system by buying natural gas if needed, no buying electricity)
get_electricity_use_offset_total = lambda: sum(i.power_utility.rate for i in sys.units) # .= 0 per hour
get_GWP_electricity_use_total = lambda: get_electricity_use_offset_total() * GWP_CFs['electricity'] * 1000 / _total_energy_per_hr() # . = 0
get_electricity_demand_total = get_electricity_use = lambda: -BT.power_utility.rate # Outside BT + BT consumption 
get_electricity_demand_cooling = lambda: F.CT.power_utility.rate + F.CT_1.power_utility.rate+ F.CWP.power_utility.rate
get_electricity_demand_non_cooling = lambda: get_electricity_demand_total() - get_electricity_demand_cooling()
get_electricity_demand_cooling_frac = lambda: get_electricity_demand_cooling() / get_electricity_demand_total()
get_electricity_demand_non_cooling_frac = lambda: get_electricity_demand_non_cooling() / get_electricity_demand_total()

# Steam contribution
get_steam_heating = lambda: sum([i.duty for i in BT.steam_utilities]) + sum([i.duty for i in F.BT_1.steam_utilities])# in kJ/hr
get_steam_electricity = lambda: 3600. * BT.electricity_demand / BT.turbogenerator_efficiency + 3600. * F.BT_1.electricity_demand / F.BT_1.turbogenerator_efficiency# in kJ/hr
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
get_GWP_ethanol = lambda: get_GWP_total_with_eletricity_credit()

# get_GWP_jet = lambda: get_GWP_total_with_eletricity_credit() * get_jet_energy_ratio() / 100 # in g CO2-eq/MJ jet
# get_GWP_diesel = lambda: get_GWP_total_with_eletricity_credit() * get_diesel_energy_ratio() / 100
# get_GWP_gasoline = lambda: get_GWP_total_with_eletricity_credit() * get_gasoline_energy_ratio() / 100
# get_SAF_abatement_cost = lambda: (get_jet_price_per_GJ() - kerosene_price_per_GJ) / (89-get_GWP_jet()) * 1000 # in $/tonne CO2


metrics.extend((Metric('GWP - total', get_GWP_total, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - electricity credit', get_GWP_electricity_credit, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - total-electricity credit', get_GWP_total_with_eletricity_credit, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - feedstock', get_GWP_feedstock_input, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - NG', get_GWP_NG, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - other materials', get_GWP_other_materials, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - electricity', get_GWP_electricity_use_total, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - non biogenic emissions', get_GWP_emissions_BT, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - other non biogenic emissions', get_GWP_emissions_non_BT, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - heating_demand', get_GWP_heating_demand, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - cooling_demand', get_GWP_cooling_demand, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - electricity non cooling', get_GWP_electricity_non_cooling, 'g CO2-eq/MJ fuel', 'LCA'),))
metrics.extend((Metric('GWP - ethanol', get_GWP_ethanol, 'g CO2-eq/MJ ethanol', 'LCA'),))

# metrics.extend((Metric('GWP - displacement', get_GWP_displacement, 'g CO2-eq/MJ', 'LCA'),)) # No displacement


#%% CREATE MODEL from evaluations.py - Blocs Dalton


# Model for state specific analysis ===========================================
def create_states_model(biorefinery):#, feed_price):
    biorefinery = biorefinery.lower()
    
    if biorefinery == 'cellulosic':
        name = 'CS'
        tea = create_tea_Ethanol() # changed to use funcion from this file
        all_states = [
                    'Alabama',
                    'Arkansas',
                    # 'California', # not rainfed
                    # 'Colorado', # not rainfed
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
                    # 'Montana', # not rainfed
                    'New Hampshire', # included now
                    'Nebraska',
                    'New Jersey',
                    # 'New Mexico', # not rainfed
                    'New York',
                    'North Carolina',
                    'North Dakota',
                    'Ohio',
                    'Oklahoma',
                    # 'Oregon', # not rainfed
                    'Pennsylvania',
                    'Rhode Island', # included now
                    'South Carolina',
                    'South Dakota',
                    'Tennessee',
                    'Texas',
                    # 'Utah', # not rainfed
                    'Vermont', # included now
                    'Virginia',
                    'West Virginia',
                    'Wisconsin',
                    ]
        states_w_inc = [
                        'Alabama', #
                        'Colorado',#
                        'Iowa', #
                        'Kansas',#
                        'Kentucky', #
                        'Louisiana', #
                        'Montana', #
                        'Nebraska', #
                        'New Mexico', #
                        'Oregon', #
                        'South Carolina', #
                        'Utah', #
                        'Virginia', #
                        ]

    model = bst.Model(tea.system, exception_hook='warn') # changed from raise to warn
    bst.PowerUtility.price = 0.0779 # electricity price Delaware (using as baseline example)
    bst.CE = 798 #  AVERAGE 2023
    tea.fuel_tax = 0.
    tea.sales_tax = 0 # Delaware (using as baseline example) 
    tea.federal_income_tax = 0.21
    tea.state_income_tax = 0.087 # Delaware (using as baseline example) # 0.065
    tea.property_tax = 0.003655 # Delaware (using as baseline example) # 0.0136
    tea.F_investment = 1.04588235294117 #Delaware (using as baseline example) # 1.02
    tea.feedstock.price = F.switchgrass.price
        
    
    def get_state_incentives(state):
            avail_incentives = st_data.loc[state]['Incentives Available']
            avail_incentives = None if pd.isna(avail_incentives) else avail_incentives # no incentives
            if avail_incentives is not None:
                try: # multiple incentives
                    avail_incentives = [int(i) for i in avail_incentives if i.isnumeric()]
                except TypeError: # only one incentive
                    avail_incentives = [int(avail_incentives)]
            return avail_incentives
        
    def solve_price():
        MFSP = tea.solve_price(F.ethanol)* _ethanol_price_conversion_index_vol()
        # try:
        #     MFSP = tea.solve_price(tea.ethanol_product)
        # except:
        #     MFSP = tea.solve_price([tea.ethanol_product], [4])
            
        # return 2.98668849 * MFSP # IN $/GAL
        return MFSP #$/gal
    
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
            # original_feedstock_price = tea.feedstock.price
            original_electricity_price = bst.PowerUtility.price
            dct = {i: getattr(tea, i) for i in names}
            tea.state_income_tax = st_data.loc[state]['Income Tax Rate (decimal)']
            tea.property_tax = st_data.loc[state]['Property Tax Rate (decimal)']
            tea.fuel_tax = st_data.loc[state]['State Motor Fuel Tax (decimal)']
            # tea.sales_tax = st_data.loc[state]['State Sales Tax Rate (decimal)'] # Don't consider sales tax because this one is only for final consumer
            bst.PowerUtility.price = st_data.loc[state]['Electricity Price (USD/kWh)']
            tea.F_investment = st_data.loc[state]['Location Capital Cost Factor (dimensionless)']
            tea.incentive_numbers = ()
            
            # tea.feedstock.price = st_data.loc[state][f'{name} Price (USD/kg)'] 
            tea.feedstock.price = tea.feedstock.price
            
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

            MFSP = solve_price()
            # tea.feedstock.price = original_feedstock_price
            bst.PowerUtility.price = original_electricity_price
            for i in names: setattr(tea, i, dct[i])
            return MFSP
        return MFSP
    
    def GWP_getter(state):
        def GWP():
            
            bst.PowerUtility.characterization_factors['GWP100'] = (st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)'], st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)'])
            GWP_CFs['electricity'] = st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)']
            GWP = get_GWP_ethanol()
            
            return GWP
        return GWP
    
    def MFSP_w_inc_getter(state):
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
            # original_feedstock_price = tea.feedstock.price
            original_electricity_price = bst.PowerUtility.price
            dct = {i: getattr(tea, i) for i in names}
            tea.state_income_tax = st_data.loc[state]['Income Tax Rate (decimal)']
            tea.property_tax = st_data.loc[state]['Property Tax Rate (decimal)']
            tea.fuel_tax = st_data.loc[state]['State Motor Fuel Tax (decimal)']
            # tea.sales_tax = st_data.loc[state]['State Sales Tax Rate (decimal)'] # Don't consider sales tax because this one is only for final consumer
            bst.PowerUtility.price = st_data.loc[state]['Electricity Price (USD/kWh)']
            tea.F_investment = st_data.loc[state]['Location Capital Cost Factor (dimensionless)']
            price['natural gas'] = st_data.loc[state]['Natural Gas Price (USD/kg)']
            # tea.incentive_numbers = get_state_incentives(state)
            # tea.feedstock.price = st_data.loc[state][f'{name} Price (USD/kg)']
            tea.feedstock.price = tea.feedstock.price
            
            
            GWP_CFs['electricity'] = st_data.loc[state]['Electricity GWP-100 (kg CO2-eq/kWh)']
            
            # Actualize for later
            price['electricity'] = bst.PowerUtility.price
            F.natural_gas.price = price['natural gas']

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

            if state == 'Alabama':
                tea.incentive_numbers = (7,)
            # elif state == 'Colorado':
            #     tea.incentive_numbers = (8,) # NO LONGER IN EFFECT
            # elif state == 'Hawaii':
            #     tea.incentive_numbers = (9,) # NO LONGER IN EFFECT
            # elif state == 'Iowa':
            #     tea.incentive_numbers = (1,10,18) # NO LONGER IN EFFECT
            elif state == 'Kansas':
                tea.incentive_numbers = (2,) 
            elif state == 'Kentucky':
                if biorefinery == 'sugarcane':
                    tea.incentive_numbers = (11,19)
                else:
                    tea.incentive_numbers = (11,12,19)
            elif state == 'Louisiana':
                tea.incentive_numbers = (13,)
            elif state == 'Montana':
                if biorefinery == 'corn':
                    tea.incentive_numbers = (3,20)
                else:
                    tea.incentive_numbers = (3,)
            elif state == 'Nebraska':
                tea.incentive_numbers = (4,)
            elif state == 'New Mexico':
                tea.incentive_numbers = (6,)
            elif state == 'Oregon':
                tea.incentive_numbers = (5,)
            # elif state == 'South Carolina':
            #     tea.incentive_numbers = (14,15)  # NO LONGER IN EFFECT
            elif state == 'Utah':
                tea.incentive_numbers = (16,)
            elif state == 'Virginia':
                tea.incentive_numbers = (17,)

            MFSP = solve_price()
            # tea.feedstock.price = original_feedstock_price
            bst.PowerUtility.price = original_electricity_price
            for i in names: setattr(tea, i, dct[i])
            return MFSP
        return MFSP
    
    @model.metric(name='Ethanol production', units='MMgal/yr')
    def get_ethanol_production():
        return get_ethanol_yield()
        #return tea.ethanol_product.F_mass / 2.98668849 * tea.operating_hours # gal /yr
    
    @model.metric(name='Ethanol conversion', units='Gal/dry ton feedstock')
    def get_ethanol_conv():
        return get_ethanol_conversion()
    
    @model.metric(name='Annual factor', units='hr/yr')
    def annual_fac():
        return get_annual_factor()
    
    @model.metric(name='Feedstock price', units='$/kg')
    def crop_price():
        return tea.feedstock.price
    
    @model.metric(name='Total capital investment', units='USD')
    def get_TCI():
        return tea.TCI
    
    @model.metric(name='Feedstock GWP', units='kgCO2e/wet kg')
    def get_feedstock_GWP():
        return F.switchgrass.characterization_factors['GWP100']
    
    # @model.metric(name='Electricity credit', units='$/Gal')
    # def get_electricity_cred():
    #     return get_cost_electricity_credit()
    
    # @model.metric(name='Electricity GWP credit', units='kgCO2/Gal')
    # def get_electricity_GWP_cred():
    #     return get_GWP_electricity_credit()
    
    # Add feedstock contribution
    
    
    @model.metric(name="Baseline MFSP", units='USD/gal') #within this function, set whatever parameter values you want to use as the baseline
    def MFSP_baseline():
        names = (
            'state_income_tax', 'property_tax', 'fuel_tax', 'sales_tax', 'F_investment',
            'incentive_numbers',
        )
        dct_old = {i: getattr(tea, i) for i in names}
        tea.state_income_tax = 0.087 # Delaware example
        tea.property_tax = 0.003655 # Delaware example
        tea.fuel_tax = 0
        tea.sales_tax = 0
        tea.F_investment = 1.04588235294117 # Delaware example
        tea.incentive_numbers = ()
        old_price = bst.PowerUtility.price
        bst.PowerUtility.price = 0.0779 # Delaware example
        tea.feedstock.price = tea.feedstock.price
        MFSP = tea.solve_price(F.ethanol)* _ethanol_price_conversion_index_vol()
        # MFSP = 2.98668849 * tea.solve_price(tea.ethanol_product) 
        for i, j in dct_old.items(): setattr(tea, i, j)
        bst.PowerUtility.price = old_price
        return MFSP
    
    
    get_inc_value = lambda: tea.exemptions.sum() + tea.deductions.sum() + tea.credits.sum()+ tea.refunds.sum()

    for state in all_states:
        model.metric(MFSP_getter(state), 'MFSP', 'USD/gal', state)
        model.metric(GWP_getter(state), 'GWP', 'g CO2-eq/MJ ethanol', state)

    for state in states_w_inc:
        model.metric(MFSP_w_inc_getter(state), 'Inc MFSP', 'USD/gal', state)
        model.metric(get_inc_value, 'Inc value', 'USD', state)
   
    
    ### Add Parameters =============================================================
    feedstock = tea.feedstock
    kg_per_ton = 1000 # using metric tonnes #907.185 if using ton
    def param(name, baseline, bounds=None, **kwargs):
        lb = 0.9 * baseline
        ub = 1.1 * baseline
        if bounds is not None:
            if lb < bounds[0]:
                lb = bounds[0]
            if ub > bounds[1]:
                ub = bounds[1]
        distribution = shape.Uniform(lb, ub)
        return model.parameter(name=name, bounds=bounds, distribution=distribution, baseline=baseline, **kwargs)
    
    if biorefinery == 'cellulosic':
        pretreatment_conversions = tea.R201.reactions.X
        cofermentation_conversions = tea.R303.cofermentation.X
        saccharification_conversions = tea.R303.saccharification.X

        # ============================================================================
        # TEA parameters
        # =============================================================================
        
        # I'M WORKING WITH A LINEAR RELATION WITH FEEDSTOCK PRICE
        # @model.parameter(element=cornstover, kind='isolated', units='USD/ton',
        #                 distribution=shape.Triangle(0.8*0.0994, 0.0994, 1.2*0.0994))
        # def set_cornstover_price(price):
        #     cornstover.price = price
        
        @param(name='Enzyme price', element=tea.cellulase, kind='isolated',
                description='price of cellulase enzyme mixture containing 50 g of cellulase per 1000g of cocktail',
                units='$USD/ton', baseline=tea.cellulase.price * kg_per_ton)
        def set_cellulase_price(price):
            tea.cellulase.price = price / kg_per_ton

        # @param(name='Plant capacity', element=cornstover, kind='coupled', units='dry US ton/yr',
        #        baseline=(cornstover.F_mass - cornstover.imass['H2O']) * tea.operating_hours / kg_per_ton,
        #        description="annual feestock processing capacity")
        # def set_plant_size(flow_rate):
        #     dry_content = 1 - cornstover.imass['H2O'] / cornstover.F_mass
        #     cornstover.F_mass = flow_rate / tea.operating_hours / dry_content * kg_per_ton
        
        # I'M KEEPING THE CAPACITY CONSTANT FOR NOW
        # D = shape.Triangle(127538.3908*0.8, 127538.3908, 127538.3908*1.2) 
        # @model.parameter(name='Feedstock flow', element='feedstock', kind='coupled', units='kg/hr',
        #                   baseline=127538.3908, distribution=D)
        # def set_feedstock_flow(flow):
        #     feedstock.F_mass = flow
            
        # D = shape.Triangle(0.078, 0.101, 0.13) 
        # @model.parameter(name='Feedstock price', element='feedstock', kind='coupled', units='$/wet kg',
        #                   baseline=0.101, distribution=D)
        # def set_feedstock_price(price):
        #     feedstock.price = price
        
        
            
        D = shape.Triangle(0.84, 0.9, 0.96)
        @model.parameter(name='Plant uptime', element='TEA', kind='isolated', units='%',
                          baseline=0.9, distribution=D)
        def set_operating_days(uptime):
            tea.operating_days = 365. * uptime
        
        @param(name='PT glucan-to-glucose', element=tea.R201, kind='coupled', units='% theoretical',
                description='extent of reaction, glucan + water -> glucose, in pretreatment reactor',
                baseline=pretreatment_conversions[0] * 100,
                bounds=(0, 100))
        def set_PT_glucan_to_glucose(X):
            X /= 100.
            pretreatment_conversions[0] = X
            corxns = pretreatment_conversions[1:3]
            corxns[:] = 0.003
            if pretreatment_conversions[:3].sum() > 1.:
                f = corxns / corxns.sum()
                corxns[:] = f * (1. - X) * 0.9999999
        
        @param(name='PT xylan-to-xylose', element=tea.R201, kind='coupled', units='% theoretical',
                description='extent of reaction, xylan + water -> xylose, in pretreatment reactor',
                baseline=pretreatment_conversions[8] * 100,
                bounds=(0, 100))
        def set_PT_xylan_to_xylose(X):
            X /= 100.
            pretreatment_conversions[8] = X
            corxns = pretreatment_conversions[9:11]
            corxns[:] = [0.024, 0.05]
            if pretreatment_conversions[8:11].sum() > 1.:
                f = corxns / corxns.sum()
                corxns[:] = f * (1. - X) * 0.9999999

        @param(name='PT xylan-to-furfural', element=tea.R201, kind='coupled', units='% theoretical',
                description='extent of reaction, xylan -> furfural + 2 water, in pretreatment reactor',
                baseline=pretreatment_conversions[10] * 100,
                bounds=(0, 100))
        def set_PT_xylan_to_furfural(X):
            # To make sure the overall xylan conversion doesn't exceed 100%
            lb = 1. - pretreatment_conversions[8] - pretreatment_conversions[9]
            pretreatment_conversions[10] = min(lb, X / 100.) * 0.9999999

        @param(name='EH cellulose-to-glucose', element=tea.R303, kind='coupled', units='% theoretical',
                description='extent of reaction, gluan + water -> glulose, in enzyme hydrolysis',
                baseline=saccharification_conversions[2] * 100,
                bounds=(0, 100))
        def set_EH_glucan_to_glucose(X):
            X /= 100.
            saccharification_conversions[2] = X
            corxns = saccharification_conversions[:2]
            corxns[:] = [0.04, 0.0012]
            if saccharification_conversions[:3].sum() > 1.:
                f = corxns / corxns.sum()
                corxns[:] = f * (1. - X) * 0.9999999
        
        @param(name='FERM glucose-to-ethanol', element=tea.R303, kind='coupled', units='% theoretical',
                description='extent of reaction, glucose -> 2 ethanol + 2 CO2, in enzyme hydrolysis',
                baseline=cofermentation_conversions[0] * 100,
                bounds=(0, 100))
        def set_FERM_glucose_to_ethanol(X):
            X /= 100.
            cofermentation_conversions[0] = X
            corxns = cofermentation_conversions[1:4]
            corxns[:] = [0.02, 0.0004, 0.006]
            if cofermentation_conversions[:4].sum() > 1.:
                f = corxns / corxns.sum()
                corxns[:] = f * (1. - X) * 0.9999999

        @param(name='Boiler efficiency', element=tea.BT, kind='coupled', units='%',
                description='efficiency of burning fuel to produce steam',
                baseline=tea.BT.boiler_efficiency * 100,
                bounds=(0, 100))
        def set_boiler_efficiency(X):
            tea.BT.boiler_efficiency = X / 100.

        @param(name='Turbogenerator efficiency', element=tea.BT, kind='coupled', units='%',
                description='efficiency of converting steam to power',
                baseline=tea.BT.turbogenerator_efficiency * 100,
                bounds=(0, 100))
        def set_turbogenerator_efficiency(X):
            tea.BT.turbogenerator_efficiency = X / 100.
        
        @param(name='Electricity price', element='TEA', kind='isolated', units='USD/kWh',
                baseline=bst.PowerUtility.price)
        def set_electricity_price(price):
            bst.PowerUtility.price = price
        
        # moisture = feedstock.imass['Water'] / feedstock.F_mass
        # D = shape.Triangle(75, 87.5, 100)
        # @model.parameter(name='Feedstock unit price', element='Feedstock', kind='isolated', units='$/dry-ton',
        #                   baseline=87.5, distribution=D)
        # def set_feedstock_price(price):
        #     feedstock.price = price * (1-moisture) / 1000 # in $/kg

        ash_disposal_price = price['ash disposal']
        D = shape.Triangle(ash_disposal_price*1.5, ash_disposal_price, ash_disposal_price*0.5)
        @model.parameter(name='Ash disposal price', element='Ash disposal', kind='isolated', units='$/kg',
                          baseline=ash_disposal_price, distribution=D)
        def set_ash_disposal_price(price):
            F.ash.price = price
        
        enzymeM301_price = F.cellulase.price
        D = shape.Triangle(enzymeM301_price*0.8, enzymeM301_price, enzymeM301_price*1.2)
        @model.parameter(name='Enzyme_price', element='Enzyme', kind='isolated', units='$/kg',
                          baseline=enzymeM301_price, distribution=D)
        def set_enzymeM301_price(price):
            F.cellulase.price = price
          
        CSL_price = price['CSL']
        D = shape.Triangle(CSL_price*0.8, CSL_price, CSL_price*1.2)
        @model.parameter(name='CSL price', element='CSL', kind='isolated', units='$/kg',
                          baseline=CSL_price, distribution=D)
        def set_CSL_price(price):
            F.CSL.price = price


        
        DAP_price = price['DAP']
        D = shape.Triangle(DAP_price*0.8, DAP_price, DAP_price*1.2)
        @model.parameter(name='DAP price', element='DAP', kind='isolated', units='$/kg',
                          baseline=DAP_price, distribution=D)
        def set_DAP_price(price):
            F.DAP.price = price

            
        ###### Fermentation parameters ######
        M301 = F.M301
        # R303 = F.R303
        
        D = shape.Triangle(0.175, 0.2, 0.25)
        @model.parameter(name='Enzymatic hydrolysis solids loading', element=M301, kind='coupled', units='%',
                          baseline=0.2, distribution=D)
        def set_M301_solids_loading(loading):
            M301.solids_loading = loading


        
        D = shape.Triangle(0.01, 0.02, 0.03)
        @model.parameter(name='Enzyme loading', element=M301, kind='coupled', units='g/g',
                          baseline=0.02, distribution=D)
        def set_M301_enzyme_loading(loading):
            M301.enzyme_loading = loading


        # D = shape.Triangle(0.75, 0.9, 0.948-1e-6)
        # @param(name='Enzymatic hydrolysis glucan-to-glucose', element=R303, kind='coupled', units='%',
        #         baseline=0.9, distribution=D)
        # def set_R303_glucan_to_glucose_conversion(X):
        #     R303.saccharification[2].X = X
               

        # D = shape.Triangle(0.9, 0.95, 0.97)
        # @param(name='Fermentation glucose-to-ethanol', element=R303, kind='coupled', units='%',
        #         baseline=0.95, distribution=D)
        # def set_R301_glucose_to_ethanol_conversion(X):
        #     R303.cofermentation[0].X = X
            
        # =============================================================================
        # LCA parameters
        # =============================================================================

        # D = shape.Triangle(-0.122, -0.0497, 0.018) #
        # @model.parameter(name='Feedstock CI', element='feedstock', kind='coupled', units='kgCO2e/wet kg',
        #                   baseline= -0.0497, distribution=D)
        # def set_feedstock_CI(CI):
        #     feedstock.characterization_factors['GWP100'] = CI
        
        D = shape.Uniform(2.24*0.05*(1-0.5), 2.24*0.05*(1+0.5))
        @model.parameter(name='Enzyme GWP', element='Enzyme', kind='isolated', units='kg CO2-eq/kg',
                         baseline=2.24*0.05, distribution=D)
        def set_enzyme_GWP(X):
            F.cellulase.characterization_factors['GWP100'] = X



        D = shape.Uniform(1.55*(1-0.5), 1.55*(1+0.5))
        @model.parameter(name='CSL GWP', element='CSL', kind='isolated', units='kg CO2-eq/kg',
                         baseline=1.55, distribution=D)
        def set_CSL_GWP(X):
            F.CSL.characterization_factors['GWP100'] = X



        D = shape.Uniform(1.6445*(1-0.5), 1.6445*(1+0.5))
        @model.parameter(name='DAP GWP', element='DAP', kind='isolated', units='kg CO2-eq/kg',
                         baseline=1.6445, distribution=D)
        def set_DAP_GWP(X):
            F.DAP.characterization_factors['GWP100'] = X
            
        
        D = shape.Uniform(2.01*(1-0.5), 2.01*(1+0.5))
        @model.parameter(name='caustic', element='caustic', kind='isolated', units='kg CO2-eq/kg',
                         baseline=2.01, distribution=D)
        def set_caustic_GWP(X):
            F.caustic.characterization_factors['GWP100'] = X  
        
        
        D = shape.Uniform(0.04509*(1-0.5), 0.04509*(1+0.5))
        @model.parameter(name='sulfuric_acid', element='sulfuric_acid', kind='isolated', units='kg CO2-eq/kg',
                         baseline=0.04509, distribution=D)
        def set_sulfuric_GWP(X):
            F.sulfuric_acid.characterization_factors['GWP100'] = X  
            
        
        D = shape.Uniform(3.1996*(1-0.5), 3.1996*(1+0.5))
        @model.parameter(name='ammonia', element='ammonia', kind='isolated', units='kg CO2-eq/kg',
                         baseline=3.1996, distribution=D)
        def set_ammonia_GWP(X):
            F.ammonia.characterization_factors['GWP100'] = X 
         
        D = shape.Uniform(0.01*(1-0.5), 0.01*(1+0.5))
        @model.parameter(name='ash', element='ash', kind='isolated', units='kg CO2-eq/kg',
                         baseline=0.01, distribution=D)
        def set_ash_GWP(X):
            F.ash.characterization_factors['GWP100'] = X
        
        
        D = shape.Uniform(0.645*(1-0.5), 0.645*(1+0.5))
        @model.parameter(name='FGD_lime', element='FGD_lime', kind='isolated', units='kg CO2-eq/kg',
                         baseline=0.645, distribution=D)
        def set_FGD_lime_GWP(X):
            F.FGD_lime.characterization_factors['GWP100'] = X
            
        D = shape.Uniform(0.4*(1-0.5), 0.4*(1+0.5))
        @model.parameter(name='natural_gas', element='natural_gas', kind='isolated', units='kg CO2-eq/kg',
                         baseline=0.4, distribution=D)
        def set_natural_gas_GWP(X):
            F.natural_gas.characterization_factors['GWP100'] = X

        D = shape.Uniform(0.8415*(1-0.5), 0.8415*(1+0.5))
        @model.parameter(name='denaturant', element='denaturant', kind='isolated', units='kg CO2-eq/kg',
                         baseline=0.8415, distribution=D)
        def set_denaturant_GWP(X):
            F.denaturant.characterization_factors['GWP100'] = X
            
        D = shape.Uniform(1.5568*(1-0.5), 1.5568*(1+0.5))
        @model.parameter(name='boiler_chemicals', element='boiler_chemicals', kind='isolated', units='kg CO2-eq/kg',
                         baseline=1.5568, distribution=D)
        def set_boiler_chemicals_GWP(X):
            F.boiler_chemicals.characterization_factors['GWP100'] = X
            
           
        base_value = GWP_CFs['electricity']
        D = shape.Uniform(base_value*(1-0.2), base_value*(1+0.2))
        @model.parameter(name='electricity', element='electricity', kind='isolated', units='kg CO2-eq/kWh',
                         baseline=base_value, distribution=D)
        def set_electricity_GWP(X):
            bst.PowerUtility.characterization_factors['GWP100'] = (X, X)
            GWP_CFs['electricity'] = X
            
    
    return model


def evaluate_SS(biorefinery, N=3000, notify_runs=0):#, feed_price = tea.feedstock.price): 
    start_time = time.time()
    model = create_states_model(biorefinery)#, feed_price)
    np.random.seed(1688)
    rule = 'L' # For Latin-Hypercube sampling
    samples = model.sample(N, rule)
    model.load_samples(samples)
    model.evaluate( notify=notify_runs) # **evaluate_args('SS'),
    model.show()
    model.get_baseline_sample()
    model.table.to_excel('model_table_Ethanol.xlsx')
    df_rho,df_p = model.spearman_r()
    df_rho.to_excel('df_rho_Eth.xlsx')
    df_p.to_excel('df_p_Eth.xlsx')
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")
    return model


#%%
# if N > 0:
#     rule=rule
#     np.random.seed(1234)
#     samples=model.sample(N, rule)
#     model.load_samples(samples)
#     model.evaluate(notify=notify_runs)
#     model.show()
#     model.get_baseline_sample()
#     model.table.to_excel('model_table_ATJ.xlsx')
#     df_rho,df_p = model.spearman_r()
#     df_rho.to_excel('df_rho_ATJ.xlsx')
#     df_p.to_excel('df_p_ATJ.xlsx')
# else:
#     model.show()
# return model

def evaluate_correlation(biorefinery,N=10000):
    model = create_IPs_model(biorefinery)
    np.random.seed(1688)
    rule = 'L' # For Latin-Hypercube sampling
    samples = model.sample(N, rule)
    model.load_samples(samples)
    model.evaluate(**evaluate_args('correlation'))
    sp_rho_table, sp_p_table = model.spearman_r()
    sp_rho_table.to_excel(get_file_name('correlation.xlsx'))
    return sp_rho_table

    
results_folder = os.path.join(folder, 'results')

def get_file_name(name):
    return os.path.join(results_folder, name)

def evaluate_args(name, nbox=None):
    if nbox is None:
        return {'autoload': True,
                'autosave': 20,
                'file': get_file_name(name)}
    else:
        nbox[0] += 1
        return {'autoload': True,
                'autosave': 20,
                'file': get_file_name(os.path.join(name, str(nbox[0])))}
    
    
