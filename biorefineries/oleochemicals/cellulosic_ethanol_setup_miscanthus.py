# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:02:04 2024

@author: lavan
"""
import biosteam as bst
import thermosteam
from biorefineries.cellulosic.chemicals import create_cellulosic_ethanol_chemicals
from biorefineries.cellulosic.systems import create_cellulosic_ethanol_system
from biorefineries.tea.cellulosic_ethanol_tea import create_cellulosic_ethanol_tea
F = bst.Flowsheet('switchgrass_flor')
bst.main_flowsheet.set_flowsheet(F)

#create chemicals
chems = create_cellulosic_ethanol_chemicals()
#TODO: change extractives to extract
chems.set_synonym('Extract','Extractives')
bst.settings.set_thermo(chems, cache= True)

# change below as needed
# bst.preferences.T = T
# bst.preferences.flow = flow_units
# bst.preferences.N = N
# bst.preferences.P = P_units
# bst.preferences.composition = True
# bst.preferences.light_mode()
# bst.preferences.save()
# bst.settings.CEPCI = CE  
# bst.settings.define_impact_indicator(key=indicator, units='kg*CO2e')
# bst.settings.set_electricity_CF(indicator,electricity_EI, electricity_EI,
#                                 basis='kWhr', units='kg*CO2e')

#create streams

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
total_flow=104229.16,
units='kg/hr',
price=0.08)

#create system
sys = create_cellulosic_ethanol_system('sys_switchgrass',ins = swg)
sys.simulate()

#saving report
# sys.save_report('swg.xlsx')

#usefull commands

#F.unit to see all units
#F.system
#F.unit_name.results()
#F.stream #can directly impact
#%% How to get heating demand, cooling demand and installed equipment cost for all units
#TODO: change example system names
#list_of_all_sys = [F.crude_HO_oil_to_biodiesel, 
#                     F.dihydroxylation_system,
#                     F.oxidative_cleavage_system,
#                     F.organic_phase_separation_and_catalyst_recovery,
#                     F.nonanoic_acid_fraction_separation,
#                     F.azelaic_acid_production]
#TODO: add facilities separately 
# other_sys = [F.BT801,F.CT701,
#               F.ADP701,F.CW701,
#               F.PWT701,F.WWT901]
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
#use the dataframe to see everything    

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
#use the dataframe to see everything cooling demand

# area_all_iec = []
# for j in list_of_all_sys:
#     for i in j.units:
#         area_all_iec.append([str(i.ID),i.installed_cost/1e6])
# for i in other_sys:
#    area_all_iec.append([str(i.ID),i.installed_cost/1e6])    
# area_all_iec = pd.DataFrame(area_all_iec, columns=['IDs', 'MMS$'])
#%% TEA code
# tea = create_cellulosic_ethanol_tea(sys = sys, OSBL_units= [F.BT, F.PWC])


#MPSP
# tea.solve_price(F.stream_name_change)
# def solve_NPV(): return tea.NPV/1e6
# def get_overall_TCI(): return tea.TCI/1e6 #[USD]  


#%% system breakdown

# df_unit_groups = bst.UnitGroup.group_by_area(sys.units)
#
#%% LCA

#add character fac for all input streams
#example
# F.ethanol.characterization_factors['GWP100'] = 1

#look at them in the report

# from biosteam import report
# report.lca_inventory_table(
#     systems=[a],
#     key=GWP,
#     items=[F.ethanol], # For including products without characterization factors
# )

#TODO: make sure to account for direct emissions
#example F.D401
