# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 09:14:21 2023

@author: lavan
"""

import biosteam as bst
import numpy as np
import chaospy
from chaospy import distributions as shape
from biorefineries.oleochemicals.systems_baseline_hosun import F
from biorefineries.oleochemicals.systems_baseline_hosun import aa_baseline_sys,crude_HO_oil_to_biodiesel,dihydroxylation_system,oxidative_cleavage_system,organic_phase_separation_and_catalyst_recovery,nonanoic_acid_fraction_separation,azelaic_acid_production
from biorefineries.oleochemicals.chemicals_baseline import chems
import numpy as np
from tea_baseline import TEA_baseline
from biosteam.evaluation import Model, Metric
from biorefineries.lipidcane._process_settings import price
from biorefineries.cane.data.lca_characterization_factors import GWP_characterization_factors
from biorefineries.tea.cellulosic_ethanol_tea import CellulosicEthanolTEA, create_cellulosic_ethanol_tea
from units_baseline import HydrolysisReactor
from biorefineries.oleochemicals import prices_and_GWP_factors
from prices_and_GWP_factors import prices_per_stream,GWP_per_stream
from biosteam import preferences
from biosteam import report



def load_preferences_and_process_settings(T,flow_units,N,P_units,CE,
                                          indicator,electricity_price,
                                          electricity_EI,heat_transfer_efficiency,
                                          steam_regeneration_price,
                                          steam_utility_T,
                                          power_utility_price,
                                          ):    
    bst.preferences.T = T
    bst.preferences.flow = flow_units
    bst.preferences.N = N
    bst.preferences.P = P_units
    bst.preferences.composition = True
    bst.preferences.light_mode()
    bst.preferences.save()
    bst.settings.CEPCI = CE  
    bst.settings.define_impact_indicator(key=indicator, units='kg*CO2e')
    bst.settings.set_electricity_CF(indicator,electricity_EI, basis='kWhr', units='kg*CO2e')
    bst.process_tools.default_utilities()
    bst.settings.electricity_price = electricity_price
    bst.steam_utility = bst.settings.get_agent('high_pressure_steam')
    bst.steam_utility.heat_transfer_efficiency = heat_transfer_efficiency
    bst.steam_utility.regeneration_price = steam_regeneration_price
    bst.steam_utility.T = steam_utility_T
    bst.steam_utility.P = chems.Water.Psat(steam_utility_T)
    bst.settings.get_agent('cooling_water').regeneration_price = 0
    bst.settings.get_agent('chilled_water').heat_transfer_price = 0
    bst.PowerUtility.price = power_utility_price 


#Feedstock type can be either HoySoy or HoSun
def set_price_of_all_streams(feedstock_type): 
    #Feedstock stream                    
    if feedstock_type == 'HoSun':
          F.crude_vegetable_oil.price = prices_per_stream['crude_vegetable_Hosun']
    else: F.crude_vegetable_oil.price = prices_per_stream['crude_vegetable_Hoysoy']
    for i in [ #streams defined inside the system
              F.stream.polystyrene_based_catalyst,
              F.stream.Liquid_HCl,
              #streams which are a part of the imported lipidcane transesterification unit
              F.stream.NaOH,
              F.stream.methanol,
              F.stream.HCl]:
        i.price = prices_per_stream[i.ID]
    #Other streams
    for i in [F.crude_HO_oil_to_biodiesel,
              F.dihydroxylation_system,
              F.oxidative_cleavage_system,
              F.organic_phase_separation_and_catalyst_recovery,
              F.nonanoic_acid_fraction_separation,
              F.azelaic_acid_production,
              F.BT901,
              F.CT901,
              F.PWT901
              ]:
            for j in i.ins:
                if j.ID in prices_per_stream.keys():
                    j.price = prices_per_stream[j.ID]
            for k in i.outs:
                if k.ID in prices_per_stream.keys():
                    k.price = prices_per_stream[k.ID]
    
                    

def set_environmental_impact_of_all_streams(indicator,feedstock_type):
    #Feedstock stream                    
    if feedstock_type == 'HoSun':
          F.crude_vegetable_oil.characterization_factors = {indicator : GWP_per_stream['crude_vegetable_HoSun']}
    else: F.crude_vegetable_oil.characterization_factors = {indicator : GWP_per_stream['crude_vegetable_Hoysoy']}
    for i in [F.crude_HO_oil_to_biodiesel,
              F.dihydroxylation_system,
              F.oxidative_cleavage_system,
              F.organic_phase_separation_and_catalyst_recovery,
              F.nonanoic_acid_fraction_separation,
              F.azelaic_acid_production]:
        for j in i.ins:
            if j.ID in GWP_per_stream.keys():                
                j.characterization_factors = {indicator : GWP_per_stream[j.ID]}               
        for k in i.outs:
            if k.ID in GWP_per_stream.keys():                
                k.characterization_factors = {indicator : GWP_per_stream[k.ID]}    
    for i in [ #streams defined inside the system
              F.stream.polystyrene_based_catalyst,
              F.stream.Liquid_HCl,
              #streams which are a part of the imported lipidcane transesterification unit
              F.stream.NaOH,
              F.stream.methanol,
              F.stream.HCl]:
        i.characterization_factors = {indicator : GWP_per_stream[k.ID]}       
       


#Techno-economic analysis
def tea_azelaic_baseline(system,WC_over_FCI,operating_days,payrate,IRR):
    tea_azelaic_baseline = TEA_baseline(lang_factor=None,
    system = system,    
    IRR=IRR,  
    duration=(2022, 2032),
    depreciation='MACRS7',#Since MARCS was not available in Biosteam  
    income_tax=0.29+0.5,  #First value denotes federal corporate income tax for Illinois, second value denotes 
    #state income tax for Illinois [2],[3]
    operating_days=operating_days,  # Uncertain, can be varied
    construction_schedule=(0.4, 0.6),  #[1]
    startup_months=0,  # Based on conventional ethanol TEA assumptions in the tea module in bioindustrial park.
    #The assumptions are derived from [4]
    startup_FOCfrac=0,  # Based on ConventionalEthanolTEA class
    startup_salesfrac=0,  # Based on ConventionalEthanolTEA class
    startup_VOCfrac=0,  #Based on ConventionalEthanolTEA class
    WC_over_FCI=WC_over_FCI,  
    finance_interest=0,  #Based on ConventionalEthanolTEA class
    finance_years=0,  #Based on ConventionalEthanolTEA class
    finance_fraction=0,  #Based on ConventionalEthanolTEA class
    operating_labor_cost=5*5*payrate*operating_days/4,#Weekly pay rate based on [5].
    #Number of shifts based on [6]
    direct_supervisory_clerical_labor=0.18,#[6]
    OSBL_units=[ #Based on the designed system
                F.CW901,
                F.CT901,
                F.BT901,
                F.PWT901,
                F.ADP901,
                F.WWT901,
                # F.HXN901
                ],
    #DPI = IEC(1+ Factor for site preparation+ Factor for building construction)
    factor_for_site_prep=0.04,#[7]
    factor_for_building_construction=0.1,#[7]
    #TDC = DPI (1+ Factor for contingency fees)
    factor_for_contingency_fees =0.18,#[7]
    #Fixed capital investment (FCI) =  TDC (1+ Factor for land purchase+ 
    #Factor for royalties+ Factor for plant startup)
    factor_for_land_purchase = 0.02,#[7]
    factor_for_plant_startup = 0.10,#[7]
    factor_for_royalties = 0.02,#[7]
    #FCI_corrected = location_factor*FCI
    location_factor = 1.15,#[7]
    # FOC = FCI* (Factor for local taxes and insurance + 
    # Factor for maintenance and repairs + 
    # Factor for operating supplies + Factor for plant overhead + 
    # Factor for administration costs ) + OLC *(Factor for direct supervisory and clerical cost + 
    # Factor for laboratory charges + Factor for administrative costs + Factor for plant overhead costs)
    factor_local_taxes_and_insurance=0.1,#0.1*FCI  [8],
    factor_maintenance_and_repairs=0.06, #0.06*FCI,range:(0.02-0.1)*FCI[8] 
    factor_operating_supplies=0.009, #[7] 0.009*FCI,range:(0.1-0.2)*FCI[8]
    factor_laboratory_charges=0.15,  #0.15*operating_labor_cost, range:(0.1-0.2)*operating_labor_cost[8]
    factor_administration_costs_fci= 0.009,
    factor_administration_costs_lc=0.177,#0.177*operating_labor_cost and 0.09*FCI[8]
    factor_plant_overhead_fci = 0.708,#0.708*operating_labor_cost + 0.036*FCI[8]
    factor_plant_overhead_lc = 0.306)#[8]  
    
    return tea_azelaic_baseline
#References    
#[1]https://doi.org/10.1016/j.fuproc.2009.04.017
#[2]https://taxfoundation.org/data/all/state/combined-federal-state-corporate-tax-rates-2022/
#[3]https://www.tax-rates.org/taxtables/income-tax-by-state
#[4]https://doi.org/10.1002/bbb.1640 
#[5]https://www.bls.gov/iag/tgs/iag325.htm#earnings
#[6]https://doi.org/10.1016/j.indcrop.2020.112411
#[7]Warren Sider
#[7]# Ref(Turton et al., 2013)
 

