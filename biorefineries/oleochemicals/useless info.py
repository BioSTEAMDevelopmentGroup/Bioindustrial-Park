# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 11:31:23 2023

@author: lavan
"""

#The steam is generated using heat obtained from burning the feed and some natural gas
#The generation of steam causes direct emissions (resulting from burning the feeds)
#Therefore,we attribute the impact of steam to NG use and emissions resulting from generation of steam
get_total_steam_GWP = lambda: get_ng_GWP() + get_total_direct_BT_emissions_GWP() #per Kg of azelaic acid

#The steam utilities (heat and electricity)
#steam produced satisfies all of the heating demand
#steam produced can also be used to satisfy some portion of the heating and cooling demand
#steam used for heating
get_BT_steam_kJph_heating = lambda: sum([i.duty for i in F.BT901.steam_utilities])*aa_baseline.operating_hours
#steam used for generating electricity
#This electricity can be for both cooling and non-cooling purposes)
#since the efficiency of the turbogen is not 100%, not all the steam gets used in the production of electricity
get_BT_steam_kJph_turbogen = lambda: KJpersec_to_KJhr*F.BT901.electricity_demand*aa_baseline.operating_hours/F.BT901.turbogenerator_efficiency 
get_BT_steam_kJph_total = lambda: get_BT_steam_kJph_heating() + get_BT_steam_kJph_turbogen() #total steam generation

#fractions of total steam used for both heating and electricity gen
get_steam_frac_heating = lambda: get_BT_steam_kJph_heating()/get_BT_steam_kJph_total()
get_steam_frac_turbogen = lambda: get_BT_steam_kJph_turbogen()/get_BT_steam_kJph_total()

#Out of the steam used for electricity generation in the turbogen
#Some is used for cooling and other is used for non-cooling purposes (like pumps, agitator etc)
get_steam_frac_electricity_cooling = lambda: get_steam_frac_turbogen() * get_elec_frac_cooling_demand()
get_steam_frac_electricity_non_cooling = lambda: get_steam_frac_turbogen() * (1-get_elec_frac_cooling_demand())

get_heating_demand_GWP = lambda: get_steam_frac_heating()*get_total_steam_GWP()
get_steam_electricity_demand_cooling_GWP = lambda: get_steam_frac_electricity_cooling()*get_total_steam_GWP()
get_cooling_demand_GWP = lambda: get_steam_electricity_demand_cooling_GWP() + get_elec_frac_cooling_demand()*get_electricity_consumption_GWP()
get_electricity_demand_non_cooling_GWP = lambda: get_steam_frac_electricity_non_cooling() * get_total_steam_GWP() + get_elec_frac_non_cooling_demand()*get_electricity_consumption_GWP() 

def get_net_GWP():
    GWP_breakdown_sum = sum([get_feedstock_GWP(),-get_other_products_impact(),
                              get_other_materials_impact(),
                              get_heating_demand_GWP(),
                              get_cooling_demand_GWP(),#includes steam used for electricity production 
                              # and other electricity used for cooling purposes 
                            get_electricity_demand_non_cooling_GWP(),#TODO; not very significant
                            get_total_non_BT_direct_emissions_GWP() #TODO; not very significant
                            ])
    return GWP_breakdown_sum 

def get_net_GWP_2():
    GWP_breakdown_2 = sum([get_feedstock_GWP(),-get_other_products_impact(),
                              get_other_materials_impact(),
                              get_heating_demand_GWP(),
                              get_cooling_demand_GWP(),#includes steam used for electricity production 
                              # and other electricity used for cooling purposes 
                            get_electricity_demand_non_cooling_GWP(),
                            get_total_emissions_GWP(),
                            -get_EOL_GWP(),-direct_emmisions_from_BT()])
    return GWP_breakdown_2

#Since in our system, BT901 is satisfying the system electricity demand, 
#all the electricity in the system is supplied by BT901
#natural gas is bought if there is not enough heat generated from the feeds to generate enough electricity
#and no electricity from the grid is purchased
#Therefore, total electricity demand can be calculated from the BT901 power utility rate
#TODO: #only consumption, replace sys.power(consum)

get_total_electricity_demand = get_electricity_use = lambda: -F.BT901.power_utility.rate*aa_baseline.operating_hours 
#Total electricity demand required for cooling 
#Can be calulated from the power utility rates of the cooling tower and the chilled water package
get_cooling_electricity_demand = lambda: (F.CT901.power_utility.rate + F.CW901.power_utility.rate )*aa_baseline.operating_hours #only consumption
#Total non cooling electricity demands for pumps,agitators etc for all the units and facilities both
get_non_cooling_electricity_demand = lambda: get_total_electricity_demand() - get_cooling_electricity_demand()
#Fractions of electricity use for both cooling and non cooling

get_elec_frac_cooling_demand = lambda: get_cooling_electricity_demand()/get_total_electricity_demand()
get_elec_frac_non_cooling_demand = lambda:get_non_cooling_electricity_demand()/get_total_electricity_demand()


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


# #######################################################################################################################
# #####################################################################################################
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
 

#TAG compositions
ooo = np.linspace(75,85,11,dtype = int)
lll = np.linspace(10,1,10,dtype = int)
lnlnln = np.linspace(9,1,9,dtype = int)
ppp = np.linspace(8,4,5,dtype = int)
sss = np.linspace(6,2,5, dtype = int)


tag = []
for o in ooo:
    tag.append([])
    for p in ppp:
        for l in lll:
            for s in sss:
                for ln in lnlnln:
                    sumi = o+p+l+ln+s
                    if int(sumi) == 98:
                        tag[-1].append({'PPP':p,'SSS':s,'OOO':o,'LLL':l,'LnLnLn':ln})
                        
#hoysoy commercial varieties


# hoysoy_commercial_vari =  {'Commercial': {'PPP': 11, 'SSS': 4, 'OOO': 22, 'LLL': 55,'LnLnLn':8},
#                            'Vistive gold': {'PPP': 3,  'SSS': 4, 'OOO': 72, 'LLL': 16,'LnLnLn':3},
#                            'Plenish': {'PPP': 6,  'SSS': 4, 'OOO': 76, 'LLL': 7, 'LnLnLn':2},
#                            'Calyno': {'PPP': 7,  'SSS': 3, 'OOO': 78, 'LLL': 3, 'LnLnLn':4},
#                            'Soyoleic':{'PPP': 7,  'SSS': 4, 'OOO': 81, 'LLL': 6, 'LnLnLn':2},
#                            'Veri':{'PPP': 7,  'SSS': 4, 'OOO': 77, 'LLL': 9, 'LnLnLn':1}}


# Convert the list of dictionaries to a DataFrame
df = pd.DataFrame(hoysoy_fcc_tag_combos)

# Compute the angle of each axis
angles = np.linspace(0, 2 * np.pi, len(df.columns), endpoint=False).tolist()

# Make the plot cyclic
angles += angles[:1]

# Initialize the spider plot
fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

# Draw one axe per variable + add labels
plt.xticks(angles[:-1], df.columns, color='grey', size=8)

# Draw filled area
for idx, row in df.iterrows():
    values = row.values.flatten().tolist()
    values += values[:1] # repeat the first value to close the circular graph
    ax.fill(angles, values, color='b', alpha=0.1)

plt.show()
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
# MPSPs_at_diff_hoysoy_varieties = []          
# for i in range(len(hoysoy_fcc_tag_combos)):
#           feed =  F.T101.ins[0]
# #changing compositions
#           feed.imass['OOO'] = hoysoy_fcc_tag_combos[i]['OOO']
#           feed.imass['PPP'] =  hoysoy_fcc_tag_combos[i]['PPP']
#           feed.imass['SSS'] = hoysoy_fcc_tag_combos[i]['SSS']
#           feed.imass['LLL'] = hoysoy_fcc_tag_combos[i]['LLL']
#           feed.imass['LnLnLn'] = hoysoy_fcc_tag_combos[i]['LnLnLn']
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

#TODO: discuss with someone about the below
# def set_feedstock_price(feedstock_price):
#         crude_vegetable_oil.price = feedstock_price
# def get_feedstock_price(sample):
#     return(sample)
# p1 = Parameter(name = 'feed price',setter = set_feedstock_price,
#                element = crude_vegetable_oil,distribution = chaospy.Uniform(),
#                kind = 'isolated',baseline = 2, bounds = (1,3),scale = 1,units = 'kg/hr',
#                system = aa_baseline, hook = get_feedstock_price, description = 'lol')


# PA_GWP = aa_baseline.get_material_impact(F.pelargonic_acid_rich_fraction,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
# C5_C8_fraction_GWP = aa_baseline.get_material_impact(recovered_C5_to_C8_MCA_fraction,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
# fatty_acid_blend_GWP = aa_baseline.get_material_impact(fatty_acid_blend,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
# glycerol_GWP = aa_baseline.get_material_impact(crude_glycerol,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
# methanol_GWP = aa_baseline.get_material_impact(crude_methanol,'GWP100')/aa_baseline.get_mass_flow(azelaic_acid)
# products_GWP_sum = PA_GWP +C5_C8_fraction_GWP +fatty_acid_blend_GWP+fatty_acid_blend_GWP +glycerol_GWP +methanol_GWP


# class Methanolseparationsystem(bst.Unit,isabstract = True):
#       _N_ins = 1
#       _N_outs = 2 
      
#       auxiliary_unit_names = ('distillation_column1',
#                                 'heat_exchanger1')
#       def __init__(self, ID='', ins=(), outs=(),
#                     thermo=None):
#           Unit.__init__(self, ID, ins, outs, thermo)
#           self.distillation_column1 = distillation_column1 = bst.BinaryDistillation(ID ='D601',
#                                                                                 ins = 'methanol_water_mixture',
                                                                                # LHK = ('Methanol',
                                                                                      # 'Water'),
                                                                                # Lr = 0.999, Hr = 0.999,
                                                                                # k = 2,
                                                                                # P = 101325)
#           self.heat_exchanger1 = heat_exchanger1 = bst.HXutility(ID = 'HX608',
#                                                             ins = 'purified_methanol_stream',
#                                                             T = 25+273,
#                                                             cool_only = True,
#                                                             rigorous= True
#                                                             )

# #This unit only recovers methanol if it is greater than 32 wt.% in the incoming feed 
# #Only cooling takes place if the stream has less than 32 wt.% of methanol in the incoming feed
# #In that case the cooled water is sent out as wastewater

#       def _run(self):
#         feed = self.ins[0]
#         top,bottom, = self.outs
#         if  feed.imass['Methanol']/feed.F_mass < 0.32:
#             self.heat_exchanger1.ins[0].copy_like(feed)
#             # self.heat_exchanger1._run()
#             self.heat_exchanger1.simulate()
#             bottom.copy_like(self.heat_exchanger1.outs[0])
#         else:
#             self.distillation_column1.ins[0].copy_like(feed)
#             # self.distillation_column1._run()
#             self.distillation_column1.simulate()
#             self.heat_exchanger1.ins[0].copy_like(self.distillation_column1.outs[0])
#             # self.heat_exchanger1._run()
#             self.heat_exchanger1.simulate()
#             top.copy_like(self.heat_exchanger1.outs[0])
#             bottom.copy_like(self.distillation_column1.outs[1])

# #TODO: figure out why cost is not working
        # def _design(self):
        #     self.distillation_column1._design()
        #     self.heat_exchanger1._design()    
          
        # def _cost(self):
        #     self.distillation_column1._cost()
        #     self.heat_exchanger1._cost()
# fraction_of_impacts_1 = [get_feedstock_GWP()*100/get_net_GWP(),
#                         get_other_materials_impact()*100/get_net_GWP(),
#                         get_heating_demand_GWP()*100/get_net_GWP(),\
#                         get_cooling_demand_GWP()*100/get_net_GWP(),\
#                         get_electricity_demand_non_cooling_GWP()*100/get_net_GWP(),\
#                         get_total_non_BT_direct_emissions_GWP()*100/get_net_GWP(),\
#                         -get_other_products_impact()*100/get_net_GWP()
#                         ]

# impacts_2 = [get_feedstock_GWP()*100/get_net_GWP(),
#                          get_other_materials_impact()*100/get_net_GWP(),
#                          get_heating_demand_GWP()*100/get_net_GWP(),
#                          get_cooling_demand_GWP()*100/get_net_GWP(),
#                         get_electricity_demand_non_cooling_GWP()*100/get_net_GWP(),
#                         get_total_emissions_GWP()*100/get_net_GWP(),
#                          -get_other_products_impact()*100/get_net_GWP(),
#                         -get_EOL_GWP()*100/get_net_GWP(),
#                         -Direct_emmisions_from_BT()*100/get_net_GWP()]
# i1 = ['Feedstock',
    # 'Other input materials',
    # 'Heating demand',
    # 'Cooling demand',
    # 'Co-products']
# i2 = ['Feedstock','Other input materials','Heating demand',
      # 'Cooling demand','Electricity demand (non cooling)',
      # 'Total emissions','Co-products','EOL emissions','Direct emissions from BT']

# # #%%HoySoy composition charts

# species = (
#     "Vistive gold",
#     "Plenish",
#     "Calyno",
#     "Soyoleic",
#     "Veri")

# weight_counts = {
#     "LnLnLn": [3,2,4,2,1],
#     "LLL": [16,7,3,6,9],
#     "OOO": [72,76,78,81,77],
#     "SSS": [4,4,3,4,4],
#     "PPP": [3,6,7,7,7]}

# width = 0.4
# fig, ax = plt.subplots()
# bottom = np.zeros(5)
# #%%
# for boolean, weight_count in weight_counts.items():
#     p = ax.bar(species, weight_count, width, label=boolean, bottom=bottom)
#     bottom += weight_count
#     ax.bar_label(p, label_type='center')
# ax.plot([10.44510,11.896538696039581,12.126500510940005,12.172363122939533,12.076725393460933],'ro')
# ax.set_title("MPSP across different compositions of Hoysoy")

# ax.legend(loc="upper right")

# plt.show()

##########################################################################################################################################################################################################################################################################################################################################    
    MT902 = bst.Mixer(ID = 'MT902',
                        ins = (
                              F.wastewater_from_cat_sep_2,
                              F.wastewater_from_cat_sep_1,
                              F.moisture_nonanal_stream_at_1atm_to_BT901,
                              F.wastewater,
                              #TODO; check how to recycle below directly instead of treating them
                              # F.R200.outs[0],
                              # F.wastewater3,
                              ),
                      outs = ('Combined_wastewater_stream_to_WWT'))

    total_air_required = sum([F.air_1.F_mass,
                              F.air_2.F_mass,
                              F.air_for_oxidative_cleavage.F_mass,
                              ])
    
    create_facilities(feedstock=None,
                      plant_air_over_feedstock=None,
                      fire_water_over_feedstock=None,
                      recycle_process_water_streams= (F.stream.water_for_emulsification1,#Water used for hydrolysis and emulsification
                                                      F.stream.water_for_emulsification2,#Water used for hydrolysis
                                                      F.stream.water_for_emulsification3,#Water used for hydrolysis
                                                      F.stream.water_for_degumming,#Water used for degumming the oils from the polar lipids
                                                      F.stream.water_for_degumming_2,#Second Water stream used for degumming the oils from the polar lipids
                                                      F.stream.biodiesel_wash_water, #Wash water for biodiesel
                                                      F.stream.water_for_azelaicacid_extraction,                                                                        
                                                      # F.CT901.outs[1],#Cooling_tower_blowdown_from_cooling tower
                                                      # F.BT901.outs[1],#rejected_water_and_blowdown from boilerturbogen   
                                                      ),
                      treated_water_streams=None, # WWT901.outs[2] #RO water from WWT901
                      CT=True,
                      CWP=True,
                      CIP=True,
                      FWT=True,
                      ADP=True,
                      WWT=True,
                      CHP=True,
                      PWC=True,
                      WWT_kwargs=dict(
                        ins = MT902-0,
                        autopopulate = False),
                      CHP_kwargs=dict(
                              combustible_slurries = (
                                                        F.polar_lipids_to_boilerturbogenerator, 
                                                        F.free_fatty_acid_waste_mixture,
                                                        F.waste_to_boilerturbogenerator,
                                                        # F.WWT901.outs[1],#sludge from WWT,
                                                      ),
                              # combustible_gases = 
                              # F.WWT901.outs[0],#biogas from WWT
                                turbogenerator_efficiency = 0.85,
                                boiler_efficiency = 0.8,
                                satisfy_system_electricity_demand =True),
                      area=900,
                      blowdown_recycle=False,
                      ADP_kwargs=dict(ins = bst.Stream('plant_air',
                                                        Oxygen = total_air_required*0.22,
                                                        Nitrogen = total_air_required*0.78)))
    
#     H100 = bst.HXutility('H100',
#                           ins = P100-0,
#                           outs = ('heated_crude_oil'),
#                           T = 273.15 + 80)#Based on acid degumming detailed in [3]
 
# #Mixing 30% of citric acid and then adding the solution 2% by vol to oil in a mixtank [3]
# #Degumming ensures that the refined oil is free from gums that are highly prone to transesterification [4]. 
#     T101 = bst.MixTank(ID = 'T101',
#                         ins = (H100-0,
#                               water_for_degumming,
#                               citricacid_for_degumming),
#                         outs = ('acid_water_mixture'),
#                         vessel_material='Stainless steel',
#                         tau = 20/60)
#     def adjust_degumming_components(): 
#       citricacid_for_degumming.imass['Citric_acid'] = 0.3 * water_for_degumming.F_mass
#       (citricacid_for_degumming+ water_for_degumming).F_vol =  0.02 * H100-0
#       T101.add_specification(adjust_degumming_components, run=True)   
# #Cooling the mixture       
#     H101 = bst.HXutility('H101',
#                           ins = T101-0,
#                           outs = ('cooled_crude_oil'),
#                           T = 273.15 + 25)
# #Adding 1% water solution to the mix [2]
#     M100 = bst.Mixer(ID = 'M100',
#                            ins = (H101-0,
#                                   water_for_degumming_2))
    
# #Holding tank to hold the mixture                           
#     T102 = bst.units.StorageTank(ID = 'T102',
#                        ins = M100-0,
#                        outs = ('water_oil_mixture'), 
#                        tau = 1)
                                                                                                   
#     def adjust_degumming_components_2():
#         water_for_degumming_2.F_vol = 0.01*H101-0
#         M100.run()            
#         M100.add_specification(adjust_degumming_components_2) 
        
# #Centrifuging the degummed oil out, 97% removal of PL using acid degumming[3]
#     C100 = bst.LiquidsSplitCentrifuge(ID = 'C100',
#                          ins = T102-0,
#                          outs = ('degummed_oil', 
#                                   polar_lipids_to_boilerturbogenerator),
#                          split = dict(PL = 0.4,
#                                       TAG = 1,
#                                       Water = 0,
#                                       Citric_acid = 0,
#                                       Oleic_acid = 1))
    
# # Neutralization is conducted to control the oil flavor and rancidity of the oil.
# # But it is also essential for preventing the saponification of free fatty acids during the transesterification [5]    
# #Assuming complete neutralisation in the Mixtank (99%)
#     M101 = units_baseline.FFA_neutralisation_tank('M101',
#                               ins = (C100-0,
#                                      base_for_saponification_of_FFA),
#                               outs = 'saponified_mix_for_separation',
#                              tau = 0.5)
#     def adjust_neutralisation_components(): 
#           base_for_saponification_of_FFA.imol['Sodium_hydroxide'] = 1.5*M101.ins[0].imol['Oleic_acid']
#     M101.add_specification(adjust_neutralisation_components, run=True)     
    
#     S100 = bst.units.Splitter('S100',
#                               ins = M101-0,
#                               outs = ('neutralised_oil_for_transesterification',
#                                        'Recovered_sodium_oleate_mixture'),
#                               split = dict(
#                                        Sodium_hydroxide = 0,
#                                        Phosphatidylinositol = 1,
#                                        OOO = 1,LnLnLn = 1,
#                                        LLL = 1,PPP = 1,
#                                        SSS = 1,Water = 1,
#                                        Oleic_acid  = 0,
#                                        Sodium_oleate =0,
#                                        ))
                          
# #Bleaching process is for bleaching colour related impurities
# #Deodorisation process is for removing odors
# #Dewaxing processes are intended for waxes present in the oil crystallizes, which give hazy appearance to the oil, 
# #These processes were ignored for this process as the produced biodiesel is an intermediate product in the process 
# # and is not sold as a co-product of the biorefinery [18]
# # Even though High oleic oils generally have an increased stability over conventional oils, 
# # Due to limited availability and high cost of the feedstock, the biodiesel derived from it might not be market competitive [7]    

kg_of_oil_per_bushel = 11.7*pounds_to_kg #Approximately 11.7 pounds of oil is present in 1 buschel of soybean oilseeds[5]
#yearly commodity soybean oil prices (2017-2023) obtained from USDA database [6]
premiums_2017_2018 = [0.2,0.6] #$perbuschel of premium over conventional soybean oilseed [7]
premiums_2018_2019 = [0.4,0.5] #$perbuschel of premium over conventional soybean oilseed, harvest delivery and onfarm storage options [8]
premiums_2019_2020 = [0.5]#$perbuschel of premium over conventional soybean oilseed [9]
premiums_2020_2021 = [0.5]#$perbuschel of premium over conventional soybean oilseed [9]
premiums_2021_2022 = [1.50,0.75,0.80,1.15,1.25,0.85,1.00,0.85,1.30,0.55,0.80] #includes on farm storage and harvest delivery options for major manufacturers[10]
premiums_2022_2023 = [1.00,2.20,1.15,2.05]#[5],[11]

soyoil_and_premium = {'2017_2018':{'soycommodityoil':0.66088,'premium':np.average(premiums_2017_2018)/kg_of_oil_per_bushel},
                      '2018_2019':{'soycommodityoil':0.62172,'premium':np.average(premiums_2018_2019)/kg_of_oil_per_bushel},
                      '2019_2020':{'soycommodityoil':0.6523,'premium':np.average(premiums_2019_2020)/kg_of_oil_per_bushel},
                      '2020_2021':{'soycommodityoil':1.25114,'premium':np.average(premiums_2020_2021)/kg_of_oil_per_bushel},
                      '2021_2022':{'soycommodityoil':1.60556,'premium':np.average(premiums_2021_2022)/kg_of_oil_per_bushel},
                      '2022_2023':{'soycommodityoil':1.452,'premium':np.average(premiums_2022_2023)/kg_of_oil_per_bushel}}

(soyoil_and_premium['2022_2023']['soycommodityoil'] + soyoil_and_premium['2022_2023']['premium'])*ratio_2023from2022
#[5] US High Oleic Soybeans & High Oleic Soybean Oil Sourcing Guide for International CustomersThird Edition – May 2023Prepared for US Soybean Export Council
#[6] https://www.ers.usda.gov/data-products/oil-crops-yearbook/oil-crops-yearbook/#All%20Tables
#[7] https://www.agriculture.com/crops/soybeans/high-oleic-soybeans-have-promise
#[8] https://corporate.perduefarms.com/news/press-releases/expansion-of-high-oleic-soybean-production/
#[9] https://www.hoosieragtoday.com/2019/09/24/consider-high-oleic-soybeans-in-2020-to-add-premium-to-your-price/



#[3]Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008)
#[4]U. N. Annal, A. Natarajan, B. Gurunathan, V. Mani, and R. Sahadevan, “Challenges and opportunities in large-scale production of biodiesel,” in Biofuels and Bioenergy, Elsevier, 2022, 
# pp. 385–407. doi: 10.1016/B978-0-323-90040-9.00033-3.
#[5] L. C. Meher, D. Vidya Sagar, and S. N. Naik, “Technical aspects of biodiesel production by 
# transesterification—a review,” Renew. Sustain. Energy Rev., vol. 10, no. 3, pp. 248–268, 
# Jun. 2006, doi: 10.1016/j.rser.2004.09.002.
#[6] “Sunflower Oil Refining Process | Sunflower Oil Refinery Plant - Goyum.” 
# https://www.oilexpeller.com/sunflower-oil-refining-process/ (accessed Mar. 31, 2023)
#[7] N. T. Dunford, E. Martínez-Force, and J. J. Salas, “High-oleic sunflower seed oil,” High 
#Oleic Oils Dev. Prop. Uses, pp. 109–124, Jan. 2022, doi: 10.1016/B978-0-12-822912-5.00004-6.

CT901_duty = -aa_baseline.operating_hours * sum([i.duty for i in F.CT901.heat_utilities if i.flow * i.duty < 0])
CW901_duty = -aa_baseline.operating_hours * sum([i.duty for i in F.CW901.heat_utilities if i.flow * i.duty < 0])
total_cooling_duty = -(CT901_duty+CW901_duty)
heat_duty_handled_by_BT901 = -aa_baseline.operating_hours * sum([i.duty for i in F.BT901.heat_utilities if i.flow * i.duty > 0])
return  heat_duty_handled_by_BT901 #KJ/yr

# %% plotting feedstock variation


# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np

# # Load the data from the Excel file
# file_path = 'just tag composition data.xlsx'
# sheet_name = '85'
# df = pd.read_excel(file_path, sheet_name=sheet_name)

# # Remove rows with  values in 'P' to avoid plotting issues
# df = df.dropna(subset=['P'])

# # Get unique values of P excluding 
# unique_p_values = df['P'].unique()

# # Create subplots side by side for each unique value of P
# fig, axs = plt.subplots(1, len(unique_p_values), figsize=(32, 7))
# plt.rcParams.update({'font.size': 14})  # Increase the font size

# for i, p_value in enumerate(unique_p_values):
#     df_filtered = df[df['P'] == p_value]
#     unique_l_values = df_filtered['L'].unique()

#     for l_value in unique_l_values:
#         df_filtered_l = df_filtered[df_filtered['L'] == l_value]
#         axs[i].scatter(df_filtered_l['Ln'], df_filtered_l['MPSP'], label=f'L = {l_value}')
#     axs[i].set_title(f'%P = {p_value}')
#     axs[i].set_xlabel('%Ln')    
#     if i == 0:
#         axs[i].set_ylabel('MPSP[$/kg]')

#     axs[i].set_xlim(0, 9.5)  # Adjust X-axis to ensure dots lie inside the frame
#     axs[i].set_xticks(np.arange(0, 10, 1))  # Ensure all integers from 0 to 10 are shown
#     axs[i].set_yticks(np.arange(9.7, 10.7, 0.1))  # Ensure all integers from 0 to 10 are shown
#     axs[i].legend()
#     axs[i].legend(loc='right', bbox_to_anchor=(0.95, 0.7))
# plt.suptitle('Impact of change in %L and %Ln on MPSP[$/kg] at constant 85% O at different %P')
# plt.tight_layout()
# plt.show()