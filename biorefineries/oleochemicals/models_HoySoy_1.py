# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:45:33 2023
@author: Lavanya
"""
import biosteam as bst
import numpy as np
import chaospy
from chaospy import distributions as shape
from biorefineries.oleochemicals.Systems_baseline_HoySoy_1 import F
# biodiesel_prep_units
from biorefineries.oleochemicals.Systems_baseline_HoySoy_1 import aa_baseline_sys
import numpy as np
from lca_tea_baseline_hosun import TEA_baseline
from biosteam.evaluation import Model, Metric
# TODO: were these prices adjusted to 2013 prices?
from biorefineries.lipidcane._process_settings import price
from biorefineries.cane.data.lca_characterization_factors import GWP_characterization_factors
from biorefineries.tea.cellulosic_ethanol_tea import CellulosicEthanolTEA, create_cellulosic_ethanol_tea
from thermosteam.utils import GG_colors, GG_light_colors
from units_baseline import HydrolysisReactor
from biorefineries.oleochemicals import prices_and_GWP_factors
from prices_and_GWP_factors import prices_per_Kg, GWP_factors,transesterification_catalyst_price
from biosteam import preferences
from biosteam import report
import contourplots
from contourplots import stacked_bar_plot

bst.preferences.T = 'K'
# bst.preferences.flow = 'kg/hr'
bst.preferences.N = 100
bst.preferences.P = 'Pa'
bst.preferences.composition = True
bst.preferences.light_mode()
bst.preferences.save()

# Settings to set GWP100 as the main characterisation factor
GWP = 'GWP100'
bst.settings.CEPCI = 708  # for 2021. All the analyses are based for Dec,2021
bst.settings.define_impact_indicator(key=GWP, units='kg*CO2e')
bst.settings.set_electricity_CF(GWP,0.36, basis='kWhr', units='kg*CO2e')# From cane biorefinery #TODO:change
bst.PowerUtility.price = 0.065  # Electricity price [USD/hr], based on lipidcane biorefinery #TODO: change
aa_baseline_sys = aa_baseline_sys() #Includes the processes and facilities
aa_baseline_sys.set_tolerance(mol=1e-5, 
                              rmol=1e-2,
                              subsystems=True)
aa_baseline_sys.simulate()


#########################################################################################################
F.stream.crude_hoysoy_vegetable_oil.characterization_factors = {'GWP100':GWP_factors['HoySoy_oil']}
F.stream.base_for_saponification_of_FFA.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}
F.stream.citricacid_for_degumming.characterization_factors = {'GWP100':GWP_factors['Citric_acid']}                                                      
F.stream.fresh_HP.characterization_factors= {'GWP100':GWP_factors['50%_HP_mix']}
F.stream.fresh_tungsten_catalyst.characterization_factors= {'GWP100':GWP_factors['Tungstic_acid']}
F.stream.fresh_cobalt_catalyst_stream.characterization_factors = {'GWP100':GWP_factors['Cobalt_acetate']}
F.stream.sodium_hydroxide_for_cat_sep.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}
F.stream.calcium_chloride_for_cat_sep.characterization_factors =  {'GWP100':GWP_factors['Cobalt_nitrate']}
F.stream.conc_hydrochloric_acid.characterization_factors = {'GWP100':GWP_factors['Conc_HCl']}
F.stream.solvent_for_MCA_extraction.characterization_factors = {'GWP100':GWP_factors['Solvent']}
F.stream.natural_gas.characterization_factors={'GWP100': GWP_characterization_factors['CH4']}
F.stream.lime_boiler.characterization_factors={'GWP100': GWP_characterization_factors['lime']}
F.stream.boiler_chems.characterization_factors={'GWP100': GWP_characterization_factors['lime']}
F.stream.system_makeup_water.characterization_factors={'GWP100': 0.00035559} 

# Streams specs belonging to create_transesterification_and_biodiesel_separation_system()
# Methanol
F.stream.methanol.price = prices_per_Kg['Methanol']
F.stream.methanol.characterization_factors = {'GWP100': GWP_characterization_factors['methanol']}
#Catalyst (25% methoxide and 75% methanol)
F.stream.catalyst.price = transesterification_catalyst_price
F.stream.catalyst.characterization_factors = {'GWP100': GWP_characterization_factors['methanol catalyst mixture']}
# HCl
F.stream.HCl.price = prices_per_Kg['HCl']
F.stream.HCl.characterization_factors = {'GWP100': GWP_factors['Conc_HCl']}
# NaOH
F.stream.NaOH.price = prices_per_Kg['Sodium_hydroxide']
F.stream.NaOH.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}
#############################################################################################################
#Since the units were imported from the cane biorefinery, the IDs of the units need to be updated in order to group them
#under biodiesel production (Area 100)
bst.rename_units(units=F.crude_HOSO_oil_to_biodiesel.units, area=100)
aa_baseline_groups = bst.UnitGroup.group_by_area(aa_baseline_sys.units)
############################################################################################################################

#BASELINE PRICES OF THE PRODUCT
#Price of C5-C9 fraction is based on adipic acid, common ester for lubricant manufacture
F.stream.recovered_C5_to_C8_MCA_fraction.price = prices_per_Kg['C5_C9_fraction']
F.stream.recovered_C5_to_C8_MCA_fraction.characterization_factors = {'GWP100': GWP_factors['C5_C9_fraction']}
#Price of pelargonic acid is based on the market price of glyphosate (conventional herbicide chemical)
F.stream.pelargonic_acid_rich_fraction.price = prices_per_Kg['Pelargonic_acid_rich_fraction']
F.stream.pelargonic_acid_rich_fraction.characterization_factors = {'GWP100': GWP_factors['Pelargonic_acid']}
#Price of fatty acid blend is based on stearic acid which is a significant portion of Matricia's fatty acid blend intermeddiate
#TODO: change description
F.stream.fatty_acid_blend.price = prices_per_Kg['Fatty_acid_blend']
F.stream.fatty_acid_blend.characterization_factors = {'GWP100': GWP_factors['Fatty_acid_blend']}
#Price of azelaic acid is based on sebacic acid (competitor chemical)
F.stream.azelaic_acid_product_stream.price =  prices_per_Kg['Azelaic_acid']
# Crude glycerol price (obtained from lipidcane.process_settings), adjusted based on Fred's producer price index to obtain price for Dec 2021
F.stream.crude_glycerol.price = prices_per_Kg['Crude_glycerol']
F.stream.crude_glycerol.characterization_factors = {'GWP100': GWP_factors['Crude_glycerol']}
# Based on Catbio costs adjusted t
F.crude_methanol.price = prices_per_Kg['Crude_methanol']
F.stream.crude_methanol.characterization_factors = {'GWP100': GWP_characterization_factors['methanol']}

azelaic_acid = F.stream.azelaic_acid_product_stream
recovered_C5_to_C8_MCA_fraction = F.stream.recovered_C5_to_C8_MCA_fraction
pelargonic_acid_rich_fraction = F.stream.pelargonic_acid_rich_fraction
fatty_acid_blend = F.stream.fatty_acid_blend

#######################################################################################################################
#####################################################################################################

tea_azelaic_baseline = TEA_baseline(
    system=aa_baseline_sys,
    lang_factor=None,
    IRR=0.10,  # Ref: Economic risk asessment paper
    duration=(2022, 2032),  # TODO:?
    depreciation='MACRS7',  # Ref: Cellulosic ethanol
    income_tax=0.35,  # Ref: Cellulosic ethanol
    operating_days=300,  # TODO:?
    construction_schedule=(2/3, 1/3),  # Ref: economic risk asessment paper
    startup_months=3,  # Ref: Cellulosic ethanol
    startup_FOCfrac=1,  # Ref: Cellulosic ethanol
    startup_salesfrac=0.5,  # Ref: Cellulosic ethanol
    startup_VOCfrac=0.75,  # Ref: Cellulosic ethanol
    WC_over_FCI=0.05,  # Ref: Cellulosic ethanol
    finance_interest=0.08,  # Ref: Cellulosic ethanol
    finance_years=10,  # Ref: Cellulosic ethanol
    finance_fraction=0.4,  # Ref: Cellulosic ethanol
    # US beaurue of stats and Ref: economic risk asessment paper
    operating_labor_cost=54874*5*5*300,#Based on US beurau of stats
    # A total of 5 operators/shift and 5 shifts per day are needed
    direct_supervisory_clerical_labor=0.18,#TODO: check economic risk assessment paper
    maintenance_and_repairs=0.06,  # Ref(Turton et al., 2013)
    operating_supplies=0.009,  # Ref(Turton et al., 2013)
    laboratory_charges=0.15,  # Ref(Turton et al., 2013)
    local_taxes_and_insurance=0.1,  # Ref(Turton et al., 2013)
    plant_overhead=0.708,  # Ref(Turton et al., 2013)
    administration_costs=0.177,  # Ref(Turton et al., 2013)
    OSBL_units=[
                F.CW901,
                F.CT901,
                F.BT901,
                F.PWT901,
                F.ADP901,
                ])

aa_sys_op_hours = aa_baseline_sys.operating_hours = tea_azelaic_baseline.operating_days * 24
azelaic_acid_tea = aa_baseline_sys.TEA
  

print(
    report.lca_inventory_table(
        systems=[aa_baseline_sys],
        key=GWP,
        items=[azelaic_acid], # For including products without characterization factors
    )
)

def get_MPSP():
    azelaic_acid.price = 0
    MPSP = azelaic_acid.price = azelaic_acid_tea.solve_price(azelaic_acid)
    return MPSP

crude_HOSO = F.stream.crude_vegetable_oil
def get_annual_factor(): return azelaic_acid_tea.operating_days*24
def get_total_yield(): return azelaic_acid.F_mass*get_annual_factor()/1e6 #To get total yield in 10^6 Kg
def get_purity(): return azelaic_acid.imass['Azelaic_acid']/azelaic_acid.F_mass
def get_overall_TCI(): return azelaic_acid_tea.TCI/1e6        
def get_material_cost(): return azelaic_acid_tea.material_cost /1e6 #includes the money spent on ash disposal and gypsum sale
def get_operating_cost(): return azelaic_acid_tea.AOC/1e6
#The boilerturbogenerator is designed to satisfy all the system heating demand
def get_system_heating_demand(): 
    heat_duty_handled_by_BT901 = -aa_baseline_sys.operating_hours * sum([i.duty for i in F.BT901.heat_utilities if i.flow * i.duty > 0])
    return  heat_duty_handled_by_BT901 #KJ/yr

def get_system_cooling_water_duty(): 
    CT901_duty = -aa_baseline_sys.operating_hours * sum([i.duty for i in F.CT901.heat_utilities if i.flow * i.duty < 0])
    CW901_duty = -aa_baseline_sys.operating_hours * sum([i.duty for i in F.CW901.heat_utilities if i.flow * i.duty < 0])
    total_cooling_duty = -(CT901_duty+CW901_duty)
    return total_cooling_duty #KJ/yr
def get_NPV(): return azelaic_acid_tea.NPV
def solve_IRR(): return azelaic_acid_tea.solve_IRR()


# Environmental impact
# TODO: net impact works but get_process_impact does not
# get_net_GWP = lambda: aa_baseline_sys.get_net_impact(key = GWP)
# Mass allocation method: TODO: check
# GWP_per_Kg_of_crude_oil = aa_baseline_sys.get_net_impact(key = GWP)/aa_baseline_sys.get_mass_flow(stream = F.crude_vegetable_oil)
# GWP_ofAA_per_Kg_of_oil = GWP_per_Kg_of_crude_oil*aa_baseline_sys.get_mass_flow(stream = azelaic_acid)/aa_baseline_sys.get_mass_flow(stream = F.crude_vegetable_oil)
# Revenue allocation
# GWP_per_USD = aa_baseline_sys.get_property_allocated_impact(
# key=GWP, name='revenue', units='USD/hr',
# )
#  GWP_AA_revenue = (
# GWP_per_USD * azelaic_acid.price
# ) # kg


metrics = [
    Metric('MPSP', get_MPSP, '$/kg'),
    Metric('Total product yield', get_total_yield, '10^6 kg/yr'),
    Metric('Product purity', get_purity, '%'),
    Metric('Total capital investment', get_overall_TCI, '10^6 $'),
    Metric('Annual operating cost', get_operating_cost, '10^6 $/yr'),
    Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
    # Metric('Annual product sale', get_azelaic_acid_sale, '10^6 $/yr'),
    Metric('Total', get_system_heating_demand, 'KJ/yr', 'Heating demand'),
    Metric('Total', get_system_cooling_water_duty,'KJ/yr', 'Cooling demand'),
    Metric('NPV', get_NPV, '$', 'TEA'),
    Metric('Internal rate of return', solve_IRR, '%'),
    # Metric('Total Net GWP', get_net_GWP, 'kg CO2-eq', 'LCA'),
]
model = Model(aa_baseline_sys,
              metrics,
              #   exception_hook='raise'
              )
# # #################################################################################################
# # #Isolated parameters
crude_oil_feedstock = F.stream.crude_hoysoy_vegetable_oil  # The feedstock stream
# Minimum price #TODO: how else can i change this
lb_o = crude_oil_feedstock.price * 0.75
ub_o = crude_oil_feedstock.price * 1.25  # Maximum price
@model.parameter(name='HoSun oil price',
                  element=crude_oil_feedstock, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_o, ub_o))
def set_feedstock_price(feedstock_price):
    crude_oil_feedstock.price = feedstock_price


fresh_tungsten_catalyst = F.stream.fresh_tungsten_catalyst
lb_t = fresh_tungsten_catalyst.price*0.75
ub_t = fresh_tungsten_catalyst.price*1.25
@model.parameter(name='Tungstic acid catalyst price',
                  element=fresh_tungsten_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_t, ub_t))
def set_tungstencat_price(tungstencat_price):
    fresh_tungsten_catalyst.price = tungstencat_price


fresh_cobalt_catalyst = F.stream.fresh_cobalt_catalyst_stream
lb_c = fresh_cobalt_catalyst.price*0.75
ub_c = fresh_cobalt_catalyst.price*1.25
@model.parameter(name='Cobalt acetate catalyst price',
                  element=fresh_cobalt_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_c, ub_c))
def set_cobaltcat_price(fresh_cobalt_catalyst_price):
    fresh_cobalt_catalyst.price = fresh_cobalt_catalyst_price


fresh_solvent = F.stream.solvent_for_extraction
fresh_solvent.price = prices_per_Kg['Heptane']
lb_s = fresh_solvent.price*0.75
ub_s = fresh_solvent.price*1.25


@model.parameter(name='Heptane solvent price',
                  element=fresh_solvent, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_s, ub_s))
def set_solvent_price(fresh_solvent_price):
    fresh_solvent.price = fresh_solvent_price


fresh_HP = F.stream.fresh_HP
lb_hp = fresh_HP.price*0.75
ub_hp = fresh_HP.price*1.25
@model.parameter(name='Hydrogen peroxide price',
                  element=fresh_HP, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_hp, ub_hp))
def set_fresh_HP_price(fresh_HP_price):
    fresh_HP.price = fresh_HP_price

F.BT901.natural_gas_price = 0.253  # REF: lactic acid SI
@model.parameter(name='Natural gas price',
                  element=F.BT901,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Triangle(0.198,
                                              0.253,
                                              0.304)  # REF: lactic acid SI
                  )
def set_natural_gas_price(natural_gas_price):
    F.BT901.natural_gas_price = natural_gas_price


citric_acid = F.stream.citricacid_for_degumming
lb_ca = citric_acid.price*0.75
ub_ca = citric_acid.price*1.25
@model.parameter(name='Citric acid price',
                  element=citric_acid,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_ca, ub_ca)
                  )
def set_citric_acid_price(citric_acid_price):
    citric_acid.price = citric_acid_price
conc_hydrochloric_acid = F.stream.conc_hydrochloric_acid
lb_hcl = conc_hydrochloric_acid.price*0.75
ub_hcl = conc_hydrochloric_acid.price*1.25


@model.parameter(name='Conc HCl price',
                  element=conc_hydrochloric_acid,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_hcl, ub_hcl)
                  )
def set_HCl_price(conc_hydrochloric_acid_price):
    conc_hydrochloric_acid.price = conc_hydrochloric_acid_price
calcium_chloride = F.stream.calcium_chloride
calcium_chloride.price = prices_per_Kg['Calcium_chloride']
lb_cc = calcium_chloride.price*0.75
ub_cc = calcium_chloride.price*1.25


@model.parameter(name='Calcium chloride price',
                  element=calcium_chloride,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_cc, ub_cc)
                  )
def set_calcium_chloride_price(calcium_chloride_price):
    calcium_chloride.price = calcium_chloride_price


polystyrene_based_catalyst = F.stream.polystyrene_based_catalyst
polystyrene_based_catalyst.price = prices_per_Kg['Resin']
lb_pbc =  prices_per_Kg['Resin']*0.75
ub_pbc =   prices_per_Kg['Resin']*1.25


@model.parameter(name='Hydrolysis resin price',
                  element=polystyrene_based_catalyst,
                  kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_pbc, ub_pbc)
                  )
def set_polystyrene_based_catalyst_price(polystyrene_based_catalyst_price):
    polystyrene_based_catalyst.price = polystyrene_based_catalyst_price


# 33
# coupled parameters
# Process related parameters
# 1 Amount ratios varied

# # Tungstic acid mole ratio
# lb1 = 0.006
# ub1 = 0.15
# @model.parameter(name='Tungstic acid moles',
#                  element=F.M200,
#                  kind='coupled',
#                  distribution=shape.Uniform(lb1, ub1)
#                  )
# def set_tungstic_acid_moles(X_tam):
#     F.unit.M200.specifications[0].args[0] = X_tam


# # Cobalt acetate catalyst ratio
# lb2 = 0.003
# ub2 = 0.015


# @model.parameter(name='Cobalt acetate moles',
#                  element=F.R300,
#                  kind='coupled',
#                  distribution=shape.Uniform(lb2, ub2)
#                  )
# def set_cobalt_acetate_moles(X_cam):
#     F.unit.R300.specifications[0].args[0] = X_cam


# # 2 Dihydroxylation reaction conversion
# lb3 = 0.90
# ub3 = 0.99


# @model.parameter(name='Dihydroxylation reaction conversion',
#                   element=F.R200,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb3, ub3))
# def set_dihydroxylation_conversion(X_dih):
#     F.unit.R200.X_dih = X_dih


# # 3 Oxidative cleavage reaction conversion
# # F.OxidativeCleavageReactor.X_oxidativecleavage = 0.8
# lb4 = 0.80
# ub4 = 0.90


# @model.parameter(name='Oxidative cleavage reaction conversion',
#                   element=F.R300,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb4, ub4))
# def set_X_oxidativecleavage_conversion(X_oxidativecleavage):
#     F.unit.R300.X_oxidativecleavage = X_oxidativecleavage


# # 3 Oxidative cleavage reaction conversion
# # F.OxidativeCleavageReactor.X_decarboxylation = 0.2
# lb5 = 0.1
# ub5 = 0.3


# @model.parameter(name='Decarboxylation reaction conversion',
#                   element=F.R300,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb5, ub5))
# def set_X_decarboxylation_conversion(X_decarboxylation):
#     F.unit.R300.X_decarboxylation = X_decarboxylation


# # 3 Oxidative cleavage reaction conversion
# # F.OxidativeCleavageReactor.X_side_rxn = 0.8
# lb6 = 0.1
# ub6 = 0.3


# @model.parameter(name='Side reaction conversion',
#                   element=F.R300,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb6, ub6))
# def set_X_side_rxn_conversion(X_side_rxn):
#     F.unit.R300.X_side_rxn = X_side_rxn


# lb7 = 35000*0.95
# ub7 = 35000*1.05


# @model.parameter(name='Feedstock_capacity',
#                   element=F.stream.crude_vegetable_oil,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb7, ub7)
#                   )
# def set_feedstock_input_flow(Cap):
#     F.stream.crude_vegetable_oil.F_mass = Cap


# lb8 = 1.052
# ub8 = 1.315


# # @model.parameter(name='Moles of air',
# #                   element=F.R300,
# #                   kind='coupled',
# #                   distribution=shape.Uniform(lb8, ub8)
# #                   )
# # def set_air_moles(X_air):
# #     F.unit.R300.specifications[0].args[1] = X_air


# # lb9 = 0.3
# # ub9 = 2.5


# # @model.parameter(name='Extraction solvent mass',
# #                   element=F.M602,
# #                   kind='coupled',
# #                   distribution=shape.Uniform(lb9, ub9)
# #                   )
# # def set_solvent_mass(X_sm):
# #     F.unit.M602.specifications[0].args[0] = X_sm


# lb10 = 0.08
# ub10 = 0.15


# @model.parameter(name='System IRR',
#                  element=azelaic_acid_tea.IRR,
#                  kind='cost',
#                  distribution=shape.Uniform(lb10, ub10)
#                  )
# def set_irr(X_irr):
#     azelaic_acid_tea.IRR = X_irr


# N_samples = 1000
# rule = 'L' # For Latin-Hypercube sampling
# np.random.seed(1234) # For consistent results
# samples = model.sample(N_samples, rule)
# model.load_samples(samples)
# model.evaluate(
#     notify=20 # Also print elapsed Time after 100 simulations
#               )
# model.show()


# plot for all the metrics stacked by area# legend(loc = 'upper right',bbox_to_anchor=(1.25,0.75))
# legend(['Dihydroxylation','Oxidative cleavage','Catalyst Recovery','Pelargonic acid recovery and C5_C8_fraction  recovery','Azelaic acid recovery','CT,PWT,BT and CWP','Biodiesel production'],loc = 'upper right',bbox_to_anchor=(1.958,0.8))
# plot.bar(stacked = True,color =  [ GG_colors.blue.RGBn,GG_light_colors.yellow.RGBn,GG_colors.green.RGBn,GG_colors.red.RGBn,GG_colors.purple.RGBn,GG_light_colors.orange.RGBn,GG_light_colors.green.RGBn],yticks = (10,20,30,40,50,60,70,80,90,100)).legend(['Dihydroxylation','Oxidative cleavage','Catalyst Recovery','Pelargonic acid recovery and C5_C8_fraction  recovery','Azelaic acid recovery','CT,PWT,BT and CWP','Biodiesel production'],loc = 'upper right',bbox_to_anchor=(1.958,0.8))
# Code to set ticks
# boxplot = plot_by_metric.boxplot(column=['Heating duty']) #for a specific column name, .boxplot() for the entire dataframe
# plot_by_area = groups_by_area.plot.bar(stacked = True,color = GG_colors.blue.RGBn)

# groups_by_area = bst.UnitGroup.df_from_groups(
#     aa_baseline_groups, fraction=True)  # populates the dataframe as percentages
# # plot for all the areas stacked by metrics
# plot_by_area = groups_by_area.plot.bar(stacked=True)
# plot_by_metric = groups_by_area.T
# groups_by_area.T.plot.bar(stacked=True, color=[GG_colors.blue.RGBn, GG_light_colors.yellow.RGBn, GG_colors.green.RGBn, GG_colors.red.RGBn,
#                           GG_colors.purple.RGBn, GG_light_colors.orange.RGBn, GG_light_colors.green.RGBn,]).set_yticks((10, 20, 30, 40, 50, 60, 70, 80, 90, 100), minor=True)

## https://www.teknorapex.com/why-esters-di-esters-polyol-esters-trimelliate-esters#:~:text=Synthetic%20esters%2C%20with%20their%20polarity%20and%20uniform%20structure%2C,metal%20surface%2C%20and%20improve%20cleanliness%20and%20sludge%20control.

# Based on Adipic acid resin, resin grade bulk, hopper cars, frt. equald. #

# bst.plots.plot_spearman_1d(rhos = (-0.0237,0.0201,0.0378,0.0131),
# index = ['Dihydroxylation reaction conversion (Uniform; 0.90 to 0.99)',
# 'Oxidative cleavage reaction conversion (Uniform; 0.80 to 0.90)',
# 'Decarboxylation reaction conversion (Uniform; 0.1 to 0.3)',
# 'Side reaction conversion (Uniform; 0.1 to 0.3)'],
# name = "NPV [$]",color = [GG_colors.blue.RGBn])

#model.table.drop([80,162,244,326,408,490,572,654,736,818,900,982],inplace = True)


# bst.plots.plot_spearman_1d(rhos = (-0.959,-0.00497,-0.0645,0.00686,-0.16,-0.0296,0.0266,-0.0345,0.0309,-0.077),
# index = ['Crude HoSun oil price (Uniform; 2.19 to 3.65)',
# 'Tungtic acid catalyst price (Uniform; 187.53 to 312.55)',
# 'Cobalt acetate catalyst price (Uniform; 36.40 to 60.675)',
# 'MCA extraction solvent price (Uniform; 0.747 to 1.245)',
# 'Hydrogen peroxide price (Uniform; 1.787 to 2.97)',
# 'Natural gas price (Triangular; 0.198,0.253,0.304)',
# 'Citric acid price (Uniform; 1.86 to 3.10)',
# 'Conc. HCl price (Uniform; 0.152 to 0.253)',
# 'Calcium chloride price (Uniform; 0.369 to 0.616)',
# 'Hydrolysis resin price (Uniform; 4.85 to 8.096)'],
# name = "NPV [$]",color = [GG_colors.blue.RGBn])

#code to get stacked bar plot
#df_breakdown = bst.UnitGroup.df_from_groups(aa_baseline_groups,
#fraction = True, scale_fractions_to_positive_values = True)
#stacked_bar_plot(dataframe = df_breakdown)
# df_breakdown.rename(columns={'Electricity consumption': 'Net electricity demand'})
 # df_breakdown.rename(index = {'100' : 'Biodiesel production','200' : 'Dihydroxylation process', '300' : 'Oxidative cleavage process', '400' : 'Catalyst recovery process', '500' : 'Pelargonic acid recovery','600' : 'Azelaic acid recovery','900':'Facilities (boiler, cooling utilities, and air distribution utlity)'},inplace = True)


# a = stacked_bar_plot(dataframe = df_breakdown,y_ticks = (0,10,20,30,40,50,60,70,80,90,100),y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
#                 y_units = "%",colors=contourplots.utils.defaults_dict['colors']['Guest_Group_TEA_Breakdown'],n_minor_ticks=4,
#                  fig_height=5.5*1.1777*0.94*1.0975)

# df_breakdown.rename(index = {'Facilities (boiler, cooling utilities, and air distribution utlity)':'Facilities (boilerturbogenerator, cooling utilities and, air distribution utility)'},inplace = True)

# a = stacked_bar_plot(dataframe = df_breakdown,y_ticks = (0,10,20,30,40,50,60,70,80,90,100),y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
#                 y_units = "%",colors=contourplots.utils.defaults_dict['colors']['Guest_Group_TEA_Breakdown'],n_minor_ticks=4,
#                  fig_height=5.5*1.1777*0.94*1.0975)


