"""
Created on Wed Feb 15 13:45:33 2023
@author: Lavanya
"""
import biosteam as bst
from biorefineries.oleochemicals.systems_baseline import F
from biorefineries.oleochemicals.systems_baseline import aa_baseline_sys
from biorefineries.oleochemicals._process_settings import load_preferences_and_process_settings,set_price_of_all_streams,set_environmental_impact_of_all_streams,tea_azelaic_baseline
from biorefineries.oleochemicals.prices_and_GWP_factors import utility_prices,Utility_GWP_factors,prices_per_Kg
from tag_compositions import high_oleic_vari_adjusted
import pandas as pd


load_preferences_and_process_settings(T = 'K',flow_units = 'kg/hr',
                                      N = 100,P_units = 'Pa',
                                      CE = 778,#2023 CE [3]
                                      indicator = 'GWP100',
                                      electricity_EI=Utility_GWP_factors['Electricity'],
                                      electricity_price = utility_prices['Electricity'])

aa_baseline = aa_baseline_sys(tag_compositions=high_oleic_vari_adjusted['Plenish'], 
                              )  
aa_baseline.set_tolerance(mol=1e-2, rmol=1e-4)

bst.rename_units(units=F.crude_HO_oil_to_biodiesel.units, area=100)
bst.rename_units(units=F.WWT901.units, area=900)

aa_baseline_groups = bst.UnitGroup.group_by_area(aa_baseline.units)
#setting prices for all inputs,products, and utilities
setting_prices = set_price_of_all_streams(feedstock_type='Plenish')
setting_EI_factors = set_environmental_impact_of_all_streams(indicator = 'GWP100',
                                                              feedstock_type='Plenish')

tea_azelaic_baseline= tea_azelaic_baseline(system = aa_baseline,
                                            operating_days=300,
                                            WC_over_FCI = 0.05,#[1]
                                            payrate = 41,#[2]
                                            IRR = 0.15) 
aa_sys_op_hours = aa_baseline.operating_hours = tea_azelaic_baseline.operating_days * 24
azelaic_acid_tea = aa_baseline.TEA
F.T101.tau = 504

aa_baseline.simulate()

# #%% To generate an excel file with total installed equipment costs of each equipment in each area
# area_1 = []
# for i in F.crude_HO_oil_to_biodiesel.units:
#     area_1.append([str(i.ID),i.installed_cost/1e6])
# df1 = pd.DataFrame(area_1, columns=['IDs', 'IEC'])

# area_2 = []
# for i in F.dihydroxylation_system.units:
#     area_2.append([str(i.ID),i.installed_cost/1e6])
# df2 = pd.DataFrame(area_2, columns=['IDs', 'IEC'])

# area_3 = []
# for i in F.oxidative_cleavage_system.units:
#     area_3.append([str(i.ID),i.installed_cost/1e6])
# df3 = pd.DataFrame(area_3, columns=['IDs', 'IEC'])

# area_4 = []
# for i in F.organic_phase_separation_and_catalyst_recovery.units:
#     area_4.append([str(i.ID),i.installed_cost/1e6])
# df4 = pd.DataFrame(area_4, columns=['IDs', 'IEC'])

# area_5 = []
# for i in F.nonanoic_acid_fraction_separation.units:
#     area_5.append([str(i.ID),i.installed_cost/1e6])
# df5 = pd.DataFrame(area_5, columns=['IDs', 'IEC'])

# area_6 = []
# for i in F.azelaic_acid_production.units:
#     area_6.append([str(i.ID),i.installed_cost/1e6])
# df6 = pd.DataFrame(area_6, columns=['IDs', 'IEC'])    

# area_7 = []
# for i in [  F.BT801,
#             F.CT701,
#             F.ADP701,
#             F.CW701,
#             F.PWT701,
#             F.WWT901,
#             F.FP701,
#             F.CIP701,
#             F.ADP701
#             ]:
      
#     area_7.append([str(i.ID),i.installed_cost/1e6])
# df7 = pd.DataFrame(area_7, columns=['IDs', 'IEC'])    
# import pandas as pd
# # Combine all the DataFrames vertically
# df_all = pd.concat([df1, df2, df3, df4, df5, df6, df7], ignore_index=True)
# # Save to Excel
# df_all.to_excel('all_metrics_across_areas.xlsx', index=False)

# #Heat duty
# # area_6_hd = []

# # for i in F.azelaic_acid_production.units:
# #     area_6_hd.append([str(i.ID),sum([j.duty for j in i.heat_utilities
# #                                      if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9 
# #                      ])
# # df6_hd = pd.DataFrame(area_6_hd, columns=['IDs', 'Heating duty [MJ]'])  

# # #Cooling duty
# # area_6_cd = []
# # for i in F.azelaic_acid_production.units:
# #     area_6_cd.append([str(i.ID),sum([j.duty for j in i.heat_utilities
# #                                      if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9 
# #                      ])
# # df6_cd = pd.DataFrame(area_6_cd, columns=['IDs', 'Cooling duty [MJ]'])  

# # All units throughtout
# list_of_all_sys = [F.crude_HO_oil_to_biodiesel,
#                     F.dihydroxylation_system,
#                     F.oxidative_cleavage_system,
#                     F.organic_phase_separation_and_catalyst_recovery,
#                     F.nonanoic_acid_fraction_separation,
#                     F.azelaic_acid_production]
# other_sys = [F.BT801,
#             F.CT701,
#             F.ADP701,
#             F.CW701,
#             F.PWT701,
#             F.WWT901,
#             F.FP701,
#             F.CIP701,
#             F.ADP701
#             ]
      
# area_all_hd = []
# for j in list_of_all_sys:
#     for i in j.units:
#         area_all_hd.append([str(i.ID),sum([j.duty for j in i.heat_utilities
#                                           if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9 
#                           ])
# for i in other_sys:
#     area_all_hd.append([str(i.ID),sum([j.duty for j in i.heat_utilities
#                                       if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9 ])        

# area_all_hd = pd.DataFrame(area_all_hd, columns=['IDs', 'Heating duty [MJ]'])  
    
# area_all_cd = []
# for j in list_of_all_sys:
#     for i in j.units:
#         area_all_cd.append([str(i.ID),sum([j.duty for j in i.heat_utilities
#                                           if j.flow > 0 and j.duty < 0])*aa_baseline.operating_hours/1e9 
#                           ])
# for i in other_sys:
#     area_all_cd.append([str(i.ID),sum([j.duty for j in i.heat_utilities
#                                       if j.flow > 0 and j.duty > 0])*aa_baseline.operating_hours/1e9 ])        
            
# area_all_cd = pd.DataFrame(area_all_cd, columns=['IDs', 'Cooling duty [MJ]'])  
# # 

# First sheet already saved
# with pd.ExcelWriter('all_metrics_across_areas.xlsx', engine='openpyxl', mode='a') as writer:
#     # Your new DataFrame, for example `area_all_hd`
#     area_all_hd.to_excel(writer, sheet_name='Cooling_Duty', index=False)

# #%%
  
# Refs
#[1] Cellulosic ethanol
#[2] https://www.bls.gov/iag/tgs/iag325.htm#earnings
#[3] https://toweringskills.com/financial-analysis/cost-indices/

