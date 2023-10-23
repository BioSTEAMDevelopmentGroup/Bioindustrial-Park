"""
Created on Wed Feb 15 13:45:33 2023
@author: Lavanya
"""
#TODO: check why M602 dissapoears and can not be accessed
import biosteam as bst
import chaospy
from chaospy import distributions as shape
from biorefineries.oleochemicals.systems_baseline_hosun import F
from biorefineries.oleochemicals.systems_baseline_hosun import aa_baseline_sys
from biorefineries.oleochemicals.chemicals_baseline import chems
from biorefineries.oleochemicals._process_settings import load_preferences_and_process_settings,set_price_of_all_streams,set_environmental_impact_of_all_streams,tea_azelaic_baseline
from biorefineries.tea.cellulosic_ethanol_tea import CellulosicEthanolTEA, create_cellulosic_ethanol_tea
from units_baseline import HydrolysisReactor
from biorefineries.oleochemicals import prices_and_GWP_factors
from prices_and_GWP_factors import prices_per_stream,utility_prices,GWP_per_stream,Utility_GWP_factors
from biosteam import preferences
from biosteam import report
from biorefineries.tea.conventional_ethanol_tea import *



load_preferences_and_process_settings(T = 'K',flow_units = 'kg/hr',
                                      N = 100,P_units = 'Pa',
                                      CE = 708,indicator = 'GWP100',
                                      electricity_EI=Utility_GWP_factors['Electricity'],
                                      electricity_price = utility_prices['Electricity'],
                                      heat_transfer_efficiency = 0.9,
                                      steam_regeneration_price = utility_prices['Steam_regeneration_price'],
                                      steam_utility_T = 620,
                                      power_utility_price = 0.065)

aa_baseline = aa_baseline_sys()     
aa_baseline.simulate()
bst.rename_units(units=F.WWT901.units, area=900)
bst.rename_units(units=F.crude_HO_oil_to_biodiesel.units, area=100)
aa_baseline_groups = bst.UnitGroup.group_by_area(aa_baseline.units)

#setting prices for all inputs,products, and utilities
setting_prices = set_price_of_all_streams(feedstock_type='HoySoy')
setting_EI_factors = set_environmental_impact_of_all_streams(indicator = 'GWP100',
                                                             feedstock_type='HoySoy')
# tea_azelaic_baseline = create_conventional_ethanol_tea
tea_azelaic_baseline= tea_azelaic_baseline(system = aa_baseline,
                                           operating_days=300,
                                           WC_over_FCI = 0.05,# WC_over_FCI=0.05,  # Ref: Cellulosic ethanol
                                           payrate = 41,
                                           IRR = 0.10)#https://www.bls.gov/iag/tgs/iag325.htm#earnings

aa_sys_op_hours = aa_baseline.operating_hours = tea_azelaic_baseline.operating_days * 24
azelaic_acid_tea = aa_baseline.TEA

  


