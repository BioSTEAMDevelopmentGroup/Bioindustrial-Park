"""
Created on Wed Feb 15 13:45:33 2023
@author: Lavanya
"""
import biosteam as bst
from biorefineries.oleochemicals.systems_baseline import F
from biorefineries.oleochemicals.systems_baseline import aa_baseline_sys
from biorefineries.oleochemicals._process_settings import load_preferences_and_process_settings,set_price_of_all_streams,set_environmental_impact_of_all_streams,tea_azelaic_baseline
from biorefineries.oleochemicals.prices_and_GWP_factors import utility_prices,Utility_GWP_factors
from tag_compositions import high_oleic_vari_adjusted
from biorefineries.tea import ConventionalEthanolTEA,create_conventional_ethanol_tea



load_preferences_and_process_settings(T = 'K',flow_units = 'kg/hr',
                                      N = 100,P_units = 'Pa',
                                      CE = 708,indicator = 'GWP100',
                                      electricity_EI=Utility_GWP_factors['Electricity'],
                                      electricity_price = utility_prices['Electricity'],
                                      heat_transfer_efficiency = 0.9,
                                      power_utility_price = 0.065)

aa_baseline = aa_baseline_sys(tag_compositions=high_oleic_vari_adjusted['Vistive gold'])  
aa_baseline.simulate()
bst.rename_units(units=F.crude_HO_oil_to_biodiesel.units, area=100)
bst.rename_units(units=F.WWT901.units, area=900)

aa_baseline_groups = bst.UnitGroup.group_by_area(aa_baseline.units)
#setting prices for all inputs,products, and utilities
setting_prices = set_price_of_all_streams(feedstock_type='Vistive gold')
setting_EI_factors = set_environmental_impact_of_all_streams(indicator = 'GWP100',
                                                              feedstock_type='HoySoy')

tea_azelaic_baseline= tea_azelaic_baseline(system = aa_baseline,
                                            operating_days=300,
                                            WC_over_FCI = 0.05,#[1]
                                            payrate = 41,#[2]
                                            IRR = 0.15) 
aa_sys_op_hours = aa_baseline.operating_hours = tea_azelaic_baseline.operating_days * 24
azelaic_acid_tea = aa_baseline.TEA
  
# Refs
#[1] Cellulosic ethanol
#[2] https://www.bls.gov/iag/tgs/iag325.htm#earnings


