# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 11:31:23 2023

@author: lavan
"""



# F.stream.crude_vegetable_oil.characterization_factors = {'GWP100':GWP_factors['HoySoy_oil']}
# F.stream.base_for_saponification_of_FFA.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}
# F.stream.crude_glycerol.characterization_factors = {'GWP100': GWP_factors['Crude_glycerol']}
# F.stream.citricacid_for_degumming.characterization_factors = {'GWP100':GWP_factors['Citric_acid']}                                                      
# F.stream.fresh_HP.characterization_factors= {'GWP100':GWP_factors['50%_HP_mix']}
# F.stream.fresh_tungsten_catalyst.characterization_factors= {'GWP100':GWP_factors['Tungstic_acid']}
# F.stream.fresh_cobalt_catalyst_stream.characterization_factors = {'GWP100':GWP_factors['Cobalt_acetate']}
# F.stream.sodium_hydroxide_for_cat_sep.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}
# F.stream.calcium_chloride_for_cat_sep.characterization_factors =  {'GWP100':GWP_factors['Cobalt_nitrate']}
# F.stream.conc_hydrochloric_acid.characterization_factors = {'GWP100':GWP_factors['Conc_HCl']}
# F.stream.solvent_for_MCA_extraction.characterization_factors = {'GWP100':GWP_factors['Solvent']}
# F.stream.natural_gas.characterization_factors={'GWP100': GWP_characterization_factors['CH4']}
# F.stream.lime_boiler.characterization_factors={'GWP100': GWP_characterization_factors['lime']}
# F.stream.boiler_chems.characterization_factors={'GWP100': GWP_characterization_factors['lime']}
# F.stream.system_makeup_water.characterization_factors={'GWP100': 0.00035559} 
# F.stream.NaOH.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}
# F.stream.conc_hydrochloric_acid.price = prices_per_Kg['HCl']
# F.stream.conc_hydrochloric_acid.characterization_factors = {'GWP100':GWP_factors['Conc_HCl']}

# F.stream.calcium_chloride_for_cat_sep.price = prices_per_Kg['Calcium_chloride']
# F.stream.calcium_chloride_for_cat_sep.characterization_factors =  {'GWP100':GWP_factors['Cobalt_nitrate']}

# # Streams specs belonging to create_transesterification_and_biodiesel_separation_system()
# F.stream.methanol.price = prices_per_Kg['Methanol']
# F.stream.methanol.characterization_factors = {'GWP100': GWP_characterization_factors['methanol']}
# F.stream.catalyst.price = transesterification_catalyst_price
# F.stream.catalyst.characterization_factors = {'GWP100': GWP_characterization_factors['methanol catalyst mixture']}
# F.stream.HCl.price = prices_per_Kg['HCl']
# F.stream.HCl.characterization_factors = {'GWP100': GWP_factors['Conc_HCl']}
# F.stream.NaOH.price = prices_per_Kg['Sodium_hydroxide']
# F.stream.crude_vegetable_oil.price = prices_per_Kg['Crude_HoySoy_oil']
# F.stream.citricacid_for_degumming.characterization_factors = {'GWP100':GWP_factors['Citric_acid']}
# F.stream.citricacid_for_degumming.price = prices_per_Kg['Citric_acid']
# F.stream.fresh_HP.price = prices_per_Kg['Hydrogen_peroxide']
# F.stream.fresh_HP.characterization_factors= {'GWP100':GWP_factors['50%_HP_mix']},
# F.stream.fresh_cobalt_catalyst_stream.characterization_factors = {'GWP100':GWP_factors['Cobalt_acetate']},
# F.stream.fresh_cobalt_catalyst_stream.price = prices_per_Kg['Cobalt_acetate']
# F.stream.base_for_saponification_of_FFA.price = prices_per_Kg['Sodium_hydroxide']
# F.stream.base_for_saponification_of_FFA.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']} #TODO: check values again
# F.stream.sodium_hydroxide_for_cat_sep.price = prices_per_Kg['Sodium_hydroxide']
# F.stream.sodium_hydroxide_for_cat_sep.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}#TODO: add latest
# F.stream.solvent_for_extraction.price = prices_per_Kg['Heptane']
# F.stream.solvent_for_extraction.characterization_factors = {'GWP100':GWP_factors['Solvent']}

# # #########################################################################################################
# F.stream.crude_vegetable_oil.characterization_factors = {'GWP100':GWP_factors['HoySoy_oil']}
# F.stream.base_for_saponification_of_FFA.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}
# F.stream.crude_glycerol.characterization_factors = {'GWP100': GWP_factors['Crude_glycerol']}
# F.stream.citricacid_for_degumming.characterization_factors = {'GWP100':GWP_factors['Citric_acid']}                                                      
# F.stream.fresh_HP.characterization_factors= {'GWP100':GWP_factors['50%_HP_mix']}
# F.stream.fresh_tungsten_catalyst.characterization_factors= {'GWP100':GWP_factors['Tungstic_acid']}
# F.stream.fresh_cobalt_catalyst_stream.characterization_factors = {'GWP100':GWP_factors['Cobalt_acetate']}
# F.stream.sodium_hydroxide_for_cat_sep.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}
# F.stream.calcium_chloride_for_cat_sep.characterization_factors =  {'GWP100':GWP_factors['Cobalt_nitrate']}
# F.stream.conc_hydrochloric_acid.characterization_factors = {'GWP100':GWP_factors['Conc_HCl']}
# F.stream.solvent_for_MCA_extraction.characterization_factors = {'GWP100':GWP_factors['Solvent']}
# F.stream.natural_gas.characterization_factors={'GWP100': GWP_characterization_factors['CH4']}
# F.stream.lime_boiler.characterization_factors={'GWP100': GWP_characterization_factors['lime']}
# F.stream.boiler_chems.characterization_factors={'GWP100': GWP_characterization_factors['lime']}
# F.stream.system_makeup_water.characterization_factors={'GWP100': 0.00035559} 
# F.stream.NaOH.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}

# F.stream.conc_hydrochloric_acid.price = prices_per_Kg['HCl']
# F.stream.conc_hydrochloric_acid.characterization_factors = {'GWP100':GWP_factors['Conc_HCl']}

# F.stream.calcium_chloride_for_cat_sep.price = prices_per_Kg['Calcium_chloride']
# F.stream.calcium_chloride_for_cat_sep.characterization_factors =  {'GWP100':GWP_factors['Cobalt_nitrate']}

# # Streams specs belonging to create_transesterification_and_biodiesel_separation_system()
# F.stream.methanol.price = prices_per_Kg['Methanol']
# F.stream.methanol.characterization_factors = {'GWP100': GWP_characterization_factors['methanol']}
# F.stream.catalyst.price = transesterification_catalyst_price
# F.stream.catalyst.characterization_factors = {'GWP100': GWP_characterization_factors['methanol catalyst mixture']}
# F.stream.HCl.price = prices_per_Kg['HCl']
# F.stream.HCl.characterization_factors = {'GWP100': GWP_factors['Conc_HCl']}
# F.stream.NaOH.price = prices_per_Kg['Sodium_hydroxide']
# F.stream.crude_vegetable_oil.price = prices_per_Kg['Crude_HoySoy_oil']
# F.stream.citricacid_for_degumming.characterization_factors = {'GWP100':GWP_factors['Citric_acid']}
# F.stream.citricacid_for_degumming.price = prices_per_Kg['Citric_acid']
# F.stream.fresh_HP.price = prices_per_Kg['Hydrogen_peroxide']
# F.stream.fresh_HP.characterization_factors= {'GWP100':GWP_factors['50%_HP_mix']},
# F.stream.fresh_cobalt_catalyst_stream.characterization_factors = {'GWP100':GWP_factors['Cobalt_acetate']},
# F.stream.fresh_cobalt_catalyst_stream.price = prices_per_Kg['Cobalt_acetate']
# F.stream.base_for_saponification_of_FFA.price = prices_per_Kg['Sodium_hydroxide']
# F.stream.base_for_saponification_of_FFA.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']} #TODO: check values again
# F.stream.sodium_hydroxide_for_cat_sep.price = prices_per_Kg['Sodium_hydroxide']
# F.stream.sodium_hydroxide_for_cat_sep.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}#TODO: add latest
# F.stream.solvent_for_extraction.price = prices_per_Kg['Heptane']
# F.stream.solvent_for_extraction.characterization_factors = {'GWP100':GWP_factors['Solvent']}
# F.BT901.natural_gas_price = 0.253  # REF: lactic acid SI
# #############################################################################################################
# #Since the units were imported from the cane biorefinery, the IDs of the units need to be updated in order to group them
# #under biodiesel production (Area 100)
# aa_baseline_groups = bst.UnitGroup.group_by_area(aa_baseline_sys.units)

# ############################################################################################################################

# #BASELINE PRICES OF THE PRODUCT
# #Price of C5-C9 fraction is based on adipic acid, common ester for lubricant manufacture
# F.stream.recovered_C5_to_C8_MCA_fraction.price = prices_per_Kg['C5_C9_fraction']
# F.stream.recovered_C5_to_C8_MCA_fraction.characterization_factors = {'GWP100': GWP_factors['C5_C9_fraction']}
# #Price of pelargonic acid is based on the market price of glyphosate (conventional herbicide chemical)
# F.stream.pelargonic_acid_rich_fraction.price = prices_per_Kg['Pelargonic_acid_rich_fraction']
# F.stream.pelargonic_acid_rich_fraction.characterization_factors = {'GWP100': GWP_factors['Pelargonic_acid']}
# #Price of fatty acid blend is based on stearic acid which is a significant portion of Matricia's fatty acid blend intermeddiate
# F.stream.fatty_acid_blend.price = prices_per_Kg['Fatty_acid_blend']
# F.stream.fatty_acid_blend.characterization_factors = {'GWP100': GWP_factors['Fatty_acid_blend']}
# #Price of azelaic acid is based on sebacic acid (competitor chemical)
# F.stream.azelaic_acid_product_stream.price =  prices_per_Kg['Azelaic_acid']
# F.stream.crude_glycerol.price = prices_per_Kg['Crude_glycerol']
# F.stream.crude_glycerol.characterization_factors = {'GWP100': GWP_factors['Crude_glycerol']}
# F.crude_methanol.price = prices_per_Kg['Crude_methanol']
# F.stream.crude_methanol.characterization_factors = {'GWP100': GWP_characterization_factors['methanol']}
# F.stream.sodium_oleate_product.price = prices_per_Kg['Sodium_oleate']


# prices_per_stream_ID = [{}]
# for i in range(len(a)):
#     for j in a[i]:
#         for k in a[i][j].keys():
#             if k in prices_per_Kg.keys():
#                 print(k,a[i][j][k])
#             else:
#                 print('nf')

# def obtain_all_streams_with_mass_comps():
#     mass_fraction_with_IDs = [{}]
#     for i in [F.crude_HO_oil_to_biodiesel,
#               F.dihydroxylation_system,
#               F.oxidative_cleavage_system,
#               F.organic_phase_separation_and_catalyst_recovery,
#               F.nonanoic_acid_fraction_separation,
#               F.azelaic_acid_production]:
#             for j in i.ins:
#                 mass_fraction_with_IDs[-1][j.ID] = {}
#                 for i in j.available_chemicals:
#                     mass_fraction_with_IDs[-1][j.ID][i.ID] = j.get_mass_fraction(i.ID)
#                 mass_fraction_with_IDs.append({})      
#     return(mass_fraction_with_IDs)                
                
# def set_price_of_all_product_streams():
#     for i in [F.crude_HO_oil_to_biodiesel,
#               F.dihydroxylation_system,
#               F.oxidative_cleavage_system,
#               F.organic_phase_separation_and_catalyst_recovery,
#               F.nonanoic_acid_fraction_separation,
#               F.azelaic_acid_production]:
#             for j in i.ins:
#                 j.price = prices_per_Kg[j.ID]         

# azelaic_acid = F.stream.azelaic_acid_product_stream
# recovered_C5_to_C8_MCA_fraction = F.stream.recovered_C5_to_C8_MCA_fraction
# pelargonic_acid_rich_fraction = F.stream.pelargonic_acid_rich_fraction
# fatty_acid_blend = F.stream.fatty_acid_blend
# crude_methanol = F.stream.crude_methanol
# crude_glycerol = F.stream.crude_glycerol
# crude_HOSO = F.stream.crude_vegetable_oil
# crude_oil_feedstock = F.stream.crude_vegetable_oil 
# fresh_tungsten_catalyst = F.stream.fresh_tungsten_catalyst
# fresh_cobalt_catalyst = F.stream.fresh_cobalt_catalyst_stream
# fresh_solvent = F.stream.solvent_for_extraction
# fresh_HP = F.stream.fresh_HP
# citric_acid = F.stream.citricacid_for_degumming
# polystyrene_based_catalyst = F.stream.polystyrene_based_catalyst
# conc_hydrochloric_acid = F.stream.conc_hydrochloric_acid
# calcium_chloride = F.stream.calcium_chloride

# #######################################################################################################################
# #####################################################################################################
# def tea_azelaic_baseline(system, IRR, duration, 
#                          depreciation,startup_months,
#                          startup_FOCfrac, startup_VOCfrac,
#                          startup_salesfrac,
#                          finance_interest,
#                          finance_years, finance_fraction,
#                          income_tax,operating_days, lang_factor,
#                          construction_schedule,
#                          WC_over_FCI,
#                          operating_labor_cost,
#                          direct_supervisory_clerical_labor,
#                          maintenance_and_repairs,
#                          operating_supplies,          
#                          laboratory_charges,
#                          local_taxes_and_insurance,
#                          plant_overhead,
#                          administration_costs,
#                          OSBL_units):
    
#     tea_azelaic_baseline = TEA_baseline(
#     system=system,lang_factor=None,
#     IRR = IRR,duration = duration,
#     depreciation = depreciation,
#     income_tax=income_tax,
#     operating_days =operating_days,
#     construction_schedule=construction_schedule,
#     startup_months = startup_months,
#     startup_FOCfrac = startup_FOCfrac,
#     startup_VOCfrac = startup_VOCfrac,
#     WC_over_FCI = WC_over_FCI,
#     finance_interest = finance_interest,
#     finance_years=finance_years,
#     finance_fraction= finance_fraction,
#     operating_labor_cost=operating_labor_cost,
#     direct_supervisory_clerical_labor=direct_supervisory_clerical_labor,
#     maintenance_and_repairs=maintenance_and_repairs,
#     operating_supplies=operating_supplies,
#     laboratory_charges=laboratory_charges,
#     local_taxes_and_insurance=local_taxes_and_insurance,
#     plant_overhead=plant_overhead,
#     administration_costs=administration_costs)
    
#     return tea_azelaic_baseline

       
   
# print(
#     report.lca_inventory_table(
#         systems=[aa_baseline_sys],
#         key=GWP,
#         items=[azelaic_acid], # For including products without characterization factors
#     ))
# print(
#     report.lca_displacement_allocation_table(
#         systems=[aa_baseline_sys],
#         key=GWP,
#         items=[azelaic_acid], # For dividing yearly impact by ethanol production
#     ))       
         
        
    #14: METHOD FOR PURIFYING AZELAIC ACID , patent number : US 2003/0032825 A1   
    #15: 72 - 1080 Kg/hour.m2 #Ref: Rule of thumb,#Power: 0.9–1.1 kW/m2#Ref: Rule of thumb for Grooved surface drums for liquids that will not wet the surface..
    #16:: https://www.usgs.gov/mission-areas/water-resources/science/industrial-water-use#science
    #17:#Composition based on latest report by NSA 2022
    #Ref book for tungstic acid waste disposal: Chemical Engineering Design Principles, Practice and Economics of Plant and Process Design By Gavin Towler, Ray Sinnott   
    #18:Storage MOC:
#Water for industrial use comes from public water supply [16]
#Composition by [17]    

 

#TAG compositions
# ooo = np.linspace(75,85,10,dtype = int)
# lll = np.linspace(10,1,9,dtype = int)
# lnlnln = np.linspace(9,1,8,dtype = int)
# ppp = np.linspace(8,4,4,dtype = int)
# sss = np.linspace(6,2,4, dtype = int)


# tag = []
# for o in ooo:
#     tag.append([])
#     for p in ppp:
#         for l in lll:
#             for s in sss:
#                 for ln in lnlnln:
#                     sumi = o+p+l+ln+s
#                     if int(sumi) == 98:
#                         tag[-1].append({'PPP':p,'SSS':s,'OOO':o,'LLL':l,'LnLnLn':ln})
                        
#hoysoy commercial varieties


# hoysoy_commercial_vari =  {'Commercial': {'PPP': 11, 'SSS': 4, 'OOO': 22, 'LLL': 55,'LnLnLn':8},
#                            'Vistive gold': {'PPP': 3,  'SSS': 4, 'OOO': 72, 'LLL': 16,'LnLnLn':3},
#                            'Plenish': {'PPP': 6,  'SSS': 4, 'OOO': 76, 'LLL': 7, 'LnLnLn':2},
#                            'Calyno': {'PPP': 7,  'SSS': 3, 'OOO': 78, 'LLL': 3, 'LnLnLn':4},
#                            'Soyoleic':{'PPP': 7,  'SSS': 4, 'OOO': 81, 'LLL': 6, 'LnLnLn':2},
#                            'Veri':{'PPP': 7,  'SSS': 4, 'OOO': 77, 'LLL': 9, 'LnLnLn':1}}

                        
# hoysoy_fcc_tag_combos = [
  # {'PPP': 8, 'SSS': 4, 'OOO': 75, 'LLL': 10, 'LnLnLn': 1},
  # {'PPP': 8, 'SSS': 3, 'OOO': 75, 'LLL': 10, 'LnLnLn': 2},
  # {'PPP': 8, 'SSS': 2, 'OOO': 75, 'LLL': 10, 'LnLnLn': 3},
   # {'PPP': 8, 'SSS': 6, 'OOO': 75, 'LLL': 8, 'LnLnLn': 1},
   # {'PPP': 8, 'SSS': 4, 'OOO': 75, 'LLL': 8, 'LnLnLn': 3},
   # {'PPP': 8, 'SSS': 3, 'OOO': 75, 'LLL': 8, 'LnLnLn': 4},
   # {'PPP': 8, 'SSS': 2, 'OOO': 75, 'LLL': 8, 'LnLnLn': 5}]
  # {'PPP': 8, 'SSS': 6, 'OOO': 75, 'LLL': 7, 'LnLnLn': 2},
  # {'PPP': 8, 'SSS': 4, 'OOO': 75, 'LLL': 7, 'LnLnLn': 4},
  # {'PPP': 8, 'SSS': 3, 'OOO': 75, 'LLL': 7, 'LnLnLn': 5},
  # {'PPP': 8, 'SSS': 2, 'OOO': 75, 'LLL': 7, 'LnLnLn': 6},
  # {'PPP': 8, 'SSS': 6, 'OOO': 75, 'LLL': 6, 'LnLnLn': 3},
  # {'PPP': 8, 'SSS': 4, 'OOO': 75, 'LLL': 6, 'LnLnLn': 5},
  # {'PPP': 8, 'SSS': 3, 'OOO': 75, 'LLL': 6, 'LnLnLn': 6},
  # {'PPP': 8, 'SSS': 2, 'OOO': 75, 'LLL': 6, 'LnLnLn': 7},
  # {'PPP': 8, 'SSS': 6, 'OOO': 75, 'LLL': 5, 'LnLnLn': 4},
  # {'PPP': 8, 'SSS': 4, 'OOO': 75, 'LLL': 5, 'LnLnLn': 6},
  # {'PPP': 8, 'SSS': 3, 'OOO': 75, 'LLL': 5, 'LnLnLn': 7},
  # {'PPP': 8, 'SSS': 6, 'OOO': 75, 'LLL': 4, 'LnLnLn': 5},
  # {'PPP': 8, 'SSS': 4, 'OOO': 75, 'LLL': 4, 'LnLnLn': 7},
  # {'PPP': 8, 'SSS': 2, 'OOO': 75, 'LLL': 4, 'LnLnLn': 9},
  # {'PPP': 8, 'SSS': 6, 'OOO': 75, 'LLL': 3, 'LnLnLn': 6},
  # {'PPP': 8, 'SSS': 3, 'OOO': 75, 'LLL': 3, 'LnLnLn': 9},
  # {'PPP': 8, 'SSS': 6, 'OOO': 75, 'LLL': 2, 'LnLnLn': 7},
  # {'PPP': 8, 'SSS': 4, 'OOO': 75, 'LLL': 2, 'LnLnLn': 9}]
   # {'PPP': 6, 'SSS': 6, 'OOO': 75, 'LLL': 10, 'LnLnLn': 1},
   # {'PPP': 6, 'SSS': 4, 'OOO': 75, 'LLL': 10, 'LnLnLn': 3},
   # {'PPP': 6, 'SSS': 3, 'OOO': 75, 'LLL': 10, 'LnLnLn': 4},
   # {'PPP': 6, 'SSS': 2, 'OOO': 75, 'LLL': 10, 'LnLnLn': 5},
   # {'PPP': 6, 'SSS': 6, 'OOO': 75, 'LLL': 8, 'LnLnLn': 3},
   # {'PPP': 6, 'SSS': 4, 'OOO': 75, 'LLL': 8, 'LnLnLn': 5},
   # {'PPP': 6, 'SSS': 3, 'OOO': 75, 'LLL': 8, 'LnLnLn': 6},
   # {'PPP': 6, 'SSS': 2, 'OOO': 75, 'LLL': 8, 'LnLnLn': 7},
   # {'PPP': 6, 'SSS': 6, 'OOO': 75, 'LLL': 7, 'LnLnLn': 4},
   # {'PPP': 6, 'SSS': 4, 'OOO': 75, 'LLL': 7, 'LnLnLn': 6},
   # {'PPP': 6, 'SSS': 3, 'OOO': 75, 'LLL': 7, 'LnLnLn': 7},
   # {'PPP': 6, 'SSS': 6, 'OOO': 75, 'LLL': 6, 'LnLnLn': 5},
   # {'PPP': 6, 'SSS': 4, 'OOO': 75, 'LLL': 6, 'LnLnLn': 7},
   # {'PPP': 6, 'SSS': 2, 'OOO': 75, 'LLL': 6, 'LnLnLn': 9},
   # {'PPP': 6, 'SSS': 6, 'OOO': 75, 'LLL': 5, 'LnLnLn': 6},
   # {'PPP': 6, 'SSS': 3, 'OOO': 75, 'LLL': 5, 'LnLnLn': 9},
   # {'PPP': 6, 'SSS': 6, 'OOO': 75, 'LLL': 4, 'LnLnLn': 7},
   # {'PPP': 6, 'SSS': 4, 'OOO': 75, 'LLL': 4, 'LnLnLn': 9},
   # {'PPP': 6, 'SSS': 6, 'OOO': 75, 'LLL': 2, 'LnLnLn': 9},
   # {'PPP': 5, 'SSS': 6, 'OOO': 75, 'LLL': 10, 'LnLnLn': 2},
   # {'PPP': 5, 'SSS': 4, 'OOO': 75, 'LLL': 10, 'LnLnLn': 4},
   # {'PPP': 5, 'SSS': 3, 'OOO': 75, 'LLL': 10, 'LnLnLn': 5},
   # {'PPP': 5, 'SSS': 2, 'OOO': 75, 'LLL': 10, 'LnLnLn': 6},
   # {'PPP': 5, 'SSS': 6, 'OOO': 75, 'LLL': 8, 'LnLnLn': 4},
   # {'PPP': 5, 'SSS': 4, 'OOO': 75, 'LLL': 8, 'LnLnLn': 6},
   # {'PPP': 5, 'SSS': 3, 'OOO': 75, 'LLL': 8, 'LnLnLn': 7},
   # {'PPP': 5, 'SSS': 6, 'OOO': 75, 'LLL': 7, 'LnLnLn': 5},
   # {'PPP': 5, 'SSS': 4, 'OOO': 75, 'LLL': 7, 'LnLnLn': 7},
   # {'PPP': 5, 'SSS': 2, 'OOO': 75, 'LLL': 7, 'LnLnLn': 9},
   # {'PPP': 5, 'SSS': 6, 'OOO': 75, 'LLL': 6, 'LnLnLn': 6},
   # {'PPP': 5, 'SSS': 3, 'OOO': 75, 'LLL': 6, 'LnLnLn': 9},
   # {'PPP': 5, 'SSS': 6, 'OOO': 75, 'LLL': 5, 'LnLnLn': 7},
   # {'PPP': 5, 'SSS': 4, 'OOO': 75, 'LLL': 5, 'LnLnLn': 9},
   # {'PPP': 5, 'SSS': 6, 'OOO': 75, 'LLL': 3, 'LnLnLn': 9},
   # {'PPP': 4, 'SSS': 6, 'OOO': 75, 'LLL': 10, 'LnLnLn': 3},
   # {'PPP': 4, 'SSS': 4, 'OOO': 75, 'LLL': 10, 'LnLnLn': 5},
   # {'PPP': 4, 'SSS': 3, 'OOO': 75, 'LLL': 10, 'LnLnLn': 6},
   # {'PPP': 4, 'SSS': 2, 'OOO': 75, 'LLL': 10, 'LnLnLn': 7},
   # {'PPP': 4, 'SSS': 6, 'OOO': 75, 'LLL': 8, 'LnLnLn': 5},
   # {'PPP': 4, 'SSS': 4, 'OOO': 75, 'LLL': 8, 'LnLnLn': 7},
   # {'PPP': 4, 'SSS': 2, 'OOO': 75, 'LLL': 8, 'LnLnLn': 9},
   # {'PPP': 4, 'SSS': 6, 'OOO': 75, 'LLL': 7, 'LnLnLn': 6},
   # {'PPP': 4, 'SSS': 3, 'OOO': 75, 'LLL': 7, 'LnLnLn': 9},
   # {'PPP': 4, 'SSS': 6, 'OOO': 75, 'LLL': 6, 'LnLnLn': 7},
   # {'PPP': 4, 'SSS': 4, 'OOO': 75, 'LLL': 6, 'LnLnLn': 9},
   # {'PPP': 4, 'SSS': 6, 'OOO': 75, 'LLL': 4, 'LnLnLn': 9}]
   
  #  {'PPP': 8, 'SSS': 3, 'OOO': 76, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 76, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 76, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 76, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 76, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 76, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 76, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 76, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 76, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 76, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 76, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 76, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 76, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 76, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 76, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 76, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 76, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 76, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 76, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 76, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 76, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 76, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 76, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 76, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 76, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 76, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 76, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 76, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 76, 'LLL': 10, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 76, 'LLL': 10, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 76, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 76, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 76, 'LLL': 8, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 76, 'LLL': 8, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 76, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 76, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 76, 'LLL': 7, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 76, 'LLL': 7, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 76, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 76, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 76, 'LLL': 6, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 76, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 76, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 76, 'LLL': 5, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 76, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 76, 'LLL': 4, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 76, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 76, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 76, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 76, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 76, 'LLL': 10, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 76, 'LLL': 10, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 76, 'LLL': 10, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 76, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 76, 'LLL': 8, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 76, 'LLL': 8, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 76, 'LLL': 8, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 76, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 76, 'LLL': 7, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 76, 'LLL': 7, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 76, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 76, 'LLL': 6, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 76, 'LLL': 6, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 76, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 76, 'LLL': 5, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 76, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 76, 'LLL': 4, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 76, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 76, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 76, 'LLL': 10, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 76, 'LLL': 10, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 76, 'LLL': 10, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 76, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 76, 'LLL': 8, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 76, 'LLL': 8, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 76, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 76, 'LLL': 7, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 76, 'LLL': 7, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 76, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 76, 'LLL': 6, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 76, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 76, 'LLL': 5, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 76, 'LLL': 3, 'LnLnLn': 9},
  # {'PPP': 8, 'SSS': 2, 'OOO': 77, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 77, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 77, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 77, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 77, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 77, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 77, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 77, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 77, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 77, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 77, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 77, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 77, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 77, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 77, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 77, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 77, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 77, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 77, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 77, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 77, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 77, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 77, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 77, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 77, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 77, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 77, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 77, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 77, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 77, 'LLL': 10, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 77, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 77, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 77, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 77, 'LLL': 8, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 77, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 77, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 77, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 77, 'LLL': 7, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 77, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 77, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 77, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 77, 'LLL': 6, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 77, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 77, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 77, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 77, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 77, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 77, 'LLL': 4, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 77, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 77, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 77, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 77, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 77, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 77, 'LLL': 10, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 77, 'LLL': 10, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 77, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 77, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 77, 'LLL': 8, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 77, 'LLL': 8, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 77, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 77, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 77, 'LLL': 7, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 77, 'LLL': 7, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 77, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 77, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 77, 'LLL': 6, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 77, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 77, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 77, 'LLL': 5, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 77, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 77, 'LLL': 4, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 77, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 77, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 77, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 77, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 77, 'LLL': 10, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 77, 'LLL': 10, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 77, 'LLL': 10, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 77, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 77, 'LLL': 8, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 77, 'LLL': 8, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 77, 'LLL': 8, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 77, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 77, 'LLL': 7, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 77, 'LLL': 7, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 77, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 77, 'LLL': 6, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 77, 'LLL': 6, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 77, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 77, 'LLL': 5, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 77, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 77, 'LLL': 4, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 77, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 78, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 78, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 78, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 78, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 78, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 78, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 78, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 78, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 78, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 78, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 78, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 78, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 78, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 78, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 78, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 78, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 78, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 78, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 78, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 78, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 78, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 78, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 78, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 78, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 78, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 78, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 78, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 78, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 78, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 78, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 78, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 78, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 78, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 78, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 78, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 78, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 78, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 78, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 78, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 78, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 78, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 78, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 78, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 78, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 78, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 78, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 78, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 78, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 78, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 78, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 78, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 78, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 78, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 78, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 78, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 78, 'LLL': 10, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 78, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 78, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 78, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 78, 'LLL': 8, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 78, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 78, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 78, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 78, 'LLL': 7, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 78, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 78, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 78, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 78, 'LLL': 6, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 78, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 78, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 78, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 78, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 78, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 78, 'LLL': 4, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 78, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 78, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 78, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 78, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 78, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 78, 'LLL': 10, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 78, 'LLL': 10, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 78, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 78, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 78, 'LLL': 8, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 78, 'LLL': 8, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 78, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 78, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 78, 'LLL': 7, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 78, 'LLL': 7, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 78, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 78, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 78, 'LLL': 6, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 78, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 78, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 78, 'LLL': 5, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 78, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 78, 'LLL': 4, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 78, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 78, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 78, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 79, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 79, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 79, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 79, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 79, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 79, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 79, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 79, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 79, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 79, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 79, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 79, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 79, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 79, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 79, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 79, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 79, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 79, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 79, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 79, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 79, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 79, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 79, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 79, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 79, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 79, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 79, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 79, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 79, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 79, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 79, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 79, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 79, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 79, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 79, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 79, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 79, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 79, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 79, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 79, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 79, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 79, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 79, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 79, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 79, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 79, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 79, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 79, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 79, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 79, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 79, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 79, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 79, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 79, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 79, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 79, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 79, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 79, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 79, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 79, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 79, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 79, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 79, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 79, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 79, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 79, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 79, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 79, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 79, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 79, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 79, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 79, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 79, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 79, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 79, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 79, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 79, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 79, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 79, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 79, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 79, 'LLL': 10, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 79, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 79, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 79, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 79, 'LLL': 8, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 79, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 79, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 79, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 79, 'LLL': 7, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 79, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 79, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 79, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 79, 'LLL': 6, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 79, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 79, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 79, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 79, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 79, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 79, 'LLL': 4, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 79, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 79, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 79, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 79, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 80, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 80, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 80, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 80, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 80, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 80, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 80, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 80, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 80, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 80, 'LLL': 3, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 80, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 80, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 80, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 80, 'LLL': 2, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 80, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 80, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 80, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 80, 'LLL': 1, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 80, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 80, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 80, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 80, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 80, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 80, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 80, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 80, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 80, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 80, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 80, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 80, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 80, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 80, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 80, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 80, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 80, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 80, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 80, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 80, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 80, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 80, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 80, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 80, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 80, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 80, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 80, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 80, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 80, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 80, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 80, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 80, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 80, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 80, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 80, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 80, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 80, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 80, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 80, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 80, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 80, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 80, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 80, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 80, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 80, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 80, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 80, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 80, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 80, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 80, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 80, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 80, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 80, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 80, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 80, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 80, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 80, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 80, 'LLL': 10, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 80, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 80, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 80, 'LLL': 8, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 80, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 80, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 80, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 80, 'LLL': 7, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 80, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 80, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 80, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 80, 'LLL': 6, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 80, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 80, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 80, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 80, 'LLL': 5, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 80, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 80, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 80, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 80, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 80, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 80, 'LLL': 3, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 80, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 80, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 80, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 80, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 81, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 81, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 81, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 81, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 81, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 81, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 81, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 81, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 81, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 81, 'LLL': 2, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 81, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 81, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 81, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 81, 'LLL': 1, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 81, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 81, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 81, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 81, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 81, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 81, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 81, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 81, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 81, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 81, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 81, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 81, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 81, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 81, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 81, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 81, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 81, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 81, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 81, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 81, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 81, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 81, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 81, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 81, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 81, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 81, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 81, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 81, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 81, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 81, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 81, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 81, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 81, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 81, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 81, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 81, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 81, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 81, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 81, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 81, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 81, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 81, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 81, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 81, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 81, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 81, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 81, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 81, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 81, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 81, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 81, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 81, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 81, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 81, 'LLL': 10, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 81, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 81, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 81, 'LLL': 8, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 81, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 81, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 81, 'LLL': 7, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 81, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 81, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 81, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 81, 'LLL': 6, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 81, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 81, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 81, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 81, 'LLL': 5, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 81, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 81, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 81, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 81, 'LLL': 4, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 81, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 81, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 81, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 81, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 81, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 81, 'LLL': 2, 'LnLnLn': 9},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 81, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 81, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 82, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 82, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 82, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 82, 'LLL': 3, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 82, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 82, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 82, 'LLL': 2, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 82, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 82, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 6, 'OOO': 82, 'LLL': 1, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 82, 'LLL': 1, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 82, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 82, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 82, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 82, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 82, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 82, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 82, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 82, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 82, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 82, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 82, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 82, 'LLL': 3, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 82, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 82, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 82, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 82, 'LLL': 2, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 82, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 82, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 82, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 82, 'LLL': 1, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 82, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 82, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 82, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 82, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 82, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 82, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 82, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 82, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 82, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 82, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 82, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 82, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 82, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 82, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 82, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 82, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 82, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 82, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 82, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 82, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 82, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 82, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 82, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 82, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 82, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 82, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 82, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 82, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 82, 'LLL': 8, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 82, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 82, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 82, 'LLL': 7, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 82, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 82, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 82, 'LLL': 6, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 82, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 82, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 82, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 82, 'LLL': 5, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 82, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 82, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 82, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 82, 'LLL': 4, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 82, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 82, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 82, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 82, 'LLL': 3, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 82, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 82, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 82, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 82, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 82, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 82, 'LLL': 1, 'LnLnLn': 9},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 83, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 83, 'LLL': 3, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 83, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 83, 'LLL': 2, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 83, 'LLL': 2, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 83, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 4, 'OOO': 83, 'LLL': 1, 'LnLnLn': 2},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 83, 'LLL': 1, 'LnLnLn': 3},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 83, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 83, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 83, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 83, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 83, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 83, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 83, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 83, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 83, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 83, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 83, 'LLL': 2, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 83, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 83, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 83, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 6, 'OOO': 83, 'LLL': 1, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 83, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 83, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 83, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 83, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 83, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 83, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 83, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 83, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 83, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 83, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 83, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 83, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 83, 'LLL': 3, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 83, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 83, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 83, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 83, 'LLL': 2, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 83, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 83, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 83, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 83, 'LLL': 1, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 83, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 83, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 83, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 83, 'LLL': 8, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 83, 'LLL': 7, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 83, 'LLL': 7, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 83, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 83, 'LLL': 6, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 83, 'LLL': 6, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 83, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 83, 'LLL': 5, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 83, 'LLL': 5, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 83, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 83, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 83, 'LLL': 4, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 83, 'LLL': 4, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 83, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 83, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 83, 'LLL': 3, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 83, 'LLL': 3, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 83, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 83, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 83, 'LLL': 2, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 83, 'LLL': 2, 'LnLnLn': 7},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 83, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 83, 'LLL': 1, 'LnLnLn': 6},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 83, 'LLL': 1, 'LnLnLn': 7},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 85, 'LLL': 2, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 3, 'OOO': 85, 'LLL': 1, 'LnLnLn': 1},
  #  {'PPP': 8, 'SSS': 2, 'OOO': 85, 'LLL': 1, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 85, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 85, 'LLL': 3, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 85, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 85, 'LLL': 2, 'LnLnLn': 1},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 85, 'LLL': 2, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 85, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 4, 'OOO': 85, 'LLL': 1, 'LnLnLn': 2},
  #  {'PPP': 6, 'SSS': 3, 'OOO': 85, 'LLL': 1, 'LnLnLn': 3},
  #  {'PPP': 6, 'SSS': 2, 'OOO': 85, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 85, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 85, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 85, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 85, 'LLL': 3, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 85, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 85, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 85, 'LLL': 2, 'LnLnLn': 2},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 85, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 85, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 6, 'OOO': 85, 'LLL': 1, 'LnLnLn': 1},
  #  {'PPP': 5, 'SSS': 4, 'OOO': 85, 'LLL': 1, 'LnLnLn': 3},
  #  {'PPP': 5, 'SSS': 3, 'OOO': 85, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 5, 'SSS': 2, 'OOO': 85, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 85, 'LLL': 6, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 85, 'LLL': 5, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 85, 'LLL': 5, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 85, 'LLL': 4, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 85, 'LLL': 4, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 85, 'LLL': 4, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 85, 'LLL': 3, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 85, 'LLL': 3, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 85, 'LLL': 3, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 85, 'LLL': 2, 'LnLnLn': 1},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 85, 'LLL': 2, 'LnLnLn': 3},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 85, 'LLL': 2, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 85, 'LLL': 2, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 6, 'OOO': 85, 'LLL': 1, 'LnLnLn': 2},
  #  {'PPP': 4, 'SSS': 4, 'OOO': 85, 'LLL': 1, 'LnLnLn': 4},
  #  {'PPP': 4, 'SSS': 3, 'OOO': 85, 'LLL': 1, 'LnLnLn': 5},
  #  {'PPP': 4, 'SSS': 2, 'OOO': 85, 'LLL': 1, 'LnLnLn': 6}]


# #plot for different combinations of high oleic soybean oil
# MPSPs_at_diff_hoysoy_varieties = []
# for i in hoysoy_vari_adjusted:
#           # MPSPs_at_diff_hoysoy_varieties.append([])
#           feed =  F.T101.ins[0]
# #changing compositions for commercial varities
#           feed.imass['OOO'] = hoysoy_vari_adjusted[i]['OOO']
#           feed.imass['PPP'] =  hoysoy_vari_adjusted[i]['PPP']
#           feed.imass['SSS'] = hoysoy_vari_adjusted[i]['SSS']
#           feed.imass['LLL'] = hoysoy_vari_adjusted[i]['LLL']
#           feed.imass['LnLnLn'] = hoysoy_vari_adjusted[i]['LnLnLn']
#           feed.imass['PL'] = 1
#           feed.imass['MAG'] = 0
#           feed.imass['DAG'] = 0
#           feed.imass['Water'] =0.05 
#           feed.imass['Oleic_acid']= 0.95
#           feed.characterization_factors = {'GWP100':GWP_factors['HoySoy_oil']}
#           feed.set_total_flow(34635.99842,'kg/hr')
#           feed.phase = 'l'
#           feed.price = prices_per_Kg['Crude_HoySoy_oil']
#           try:
#               aa_baseline_sys.simulate()
#           except:
#               pass
#           MPSP_at_i  = tea_azelaic_baseline.solve_price(azelaic_acid)
#           MPSPs_at_diff_hoysoy_varieties.append(MPSP_at_i)
#for all varieties 
MPSPs_at_diff_hoysoy_varieties = []          
for i in range(len(hoysoy_fcc_tag_combos)):
          feed =  F.T101.ins[0]
#changing compositions
          feed.imass['OOO'] = hoysoy_fcc_tag_combos[i]['OOO']
          feed.imass['PPP'] =  hoysoy_fcc_tag_combos[i]['PPP']
          feed.imass['SSS'] = hoysoy_fcc_tag_combos[i]['SSS']
          feed.imass['LLL'] = hoysoy_fcc_tag_combos[i]['LLL']
          feed.imass['LnLnLn'] = hoysoy_fcc_tag_combos[i]['LnLnLn']
          feed.imass['PL'] = 1
          feed.imass['MAG'] = 0
          feed.imass['DAG'] = 0
          feed.imass['Water'] =0.05 
          feed.imass['Oleic_acid']= 0.95
          feed.characterization_factors = {'GWP100':GWP_factors['HoySoy_oil']}
          feed.set_total_flow(34635.99842,'kg/hr')
          feed.phase = 'l'
          feed.price = prices_per_Kg['Crude_HoySoy_oil']
          try:
              aa_baseline_sys.simulate()
          except:
              pass
          MPSP_at_i  = tea_azelaic_baseline.solve_price(azelaic_acid)
          MPSPs_at_diff_hoysoy_varieties.append(MPSP_at_i)
