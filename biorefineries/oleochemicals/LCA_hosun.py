#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 08:34:22 2023

@author: lavanyakudli
"""
import biosteam as bst
import numpy as np
import chaospy
from chaospy import distributions as shape
from biorefineries.oleochemicals.systems_baseline_hosun import F
from biorefineries.oleochemicals.systems_baseline_hosun import aa_baseline_sys
import numpy as np
from lca_tea_baseline_hosun import TEA_baseline
from biosteam.evaluation import Model, Metric
# TODO: were these prices adjusted to 2013 prices?
from biorefineries.lipidcane._process_settings import price
from biorefineries.cane.data.lca_characterization_factors import GWP_characterization_factors
from biorefineries.tea.cellulosic_ethanol_tea import CellulosicEthanolTEA, create_cellulosic_ethanol_tea
from units_baseline import HydrolysisReactor
from biorefineries.oleochemicals import prices_and_GWP_factors
from prices_and_GWP_factors import prices_per_Kg, GWP_factors,transesterification_catalyst_price
from biosteam import preferences
from biosteam import report
import contourplots
from contourplots import stacked_bar_plot,box_and_whiskers_plot
import models_hosun
from models_hosun import *
import chemicals_baseline
from chemicals_baseline import *
import systems_baseline_hosun 
from systems_baseline_hosun import *
import units_baseline
from units_baseline import *

# # #################################################################################################
# # #Isolated parameters
# # # # Minimum price #TODO: how else can i change this

lb_o = crude_oil_feedstock.price * 0.75
ub_o = crude_oil_feedstock.price * 1.25  # Maximum price
@model.parameter(name='HoSun oil price',
                  element=crude_oil_feedstock, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_o, ub_o))
def set_feedstock_price(feedstock_price):
    crude_oil_feedstock.price = feedstock_price


lb_t = fresh_tungsten_catalyst.price*0.75
ub_t = fresh_tungsten_catalyst.price*1.25
@model.parameter(name='Tungstic acid catalyst price',
                  element=fresh_tungsten_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_t, ub_t))
def set_tungstencat_price(tungstencat_price):
    fresh_tungsten_catalyst.price = tungstencat_price


lb_c = fresh_cobalt_catalyst.price*0.75
ub_c = fresh_cobalt_catalyst.price*1.25
@model.parameter(name='Cobalt acetate catalyst price',
                  element=fresh_cobalt_catalyst, kind='isolated', units='USD/kg',
                  distribution=shape.Uniform(lb_c, ub_c))
def set_cobaltcat_price(fresh_cobalt_catalyst_price):
    fresh_cobalt_catalyst.price = fresh_cobalt_catalyst_price


fresh_solvent.price = prices_per_Kg['Heptane']
lb_s = fresh_solvent.price*0.75
ub_s = fresh_solvent.price*1.25


@model.parameter(name='Heptane solvent price',
                  element=fresh_solvent, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_s, ub_s))
def set_solvent_price(fresh_solvent_price):
    fresh_solvent.price = fresh_solvent_price


lb_hp = fresh_HP.price*0.75
ub_hp = fresh_HP.price*1.25
@model.parameter(name='Hydrogen peroxide price',
                  element=fresh_HP, kind='isolated',
                  units='USD/kg',
                  distribution=shape.Uniform(lb_hp, ub_hp))
def set_fresh_HP_price(fresh_HP_price):
    fresh_HP.price = fresh_HP_price


# @model.parameter(name='Natural gas price',
#                   element=F.BT901,
#                   kind='isolated',
#                   units='USD/kg',
#                   distribution=shape.Triangle(0.198,
#                                               0.253,
#                                               0.304)  # REF: lactic acid SI
#                   )
# def set_natural_gas_price(natural_gas_price):
#     F.BT901.natural_gas_price = natural_gas_price


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


# coupled parameters
# Process related parameters
##Tungstic acid mole ratio
lb1 = 0.006
ub1 = 0.15
@model.parameter(name='Tungstic acid moles',
                  element=F.M200,
                  kind='coupled',
                  distribution=shape.Uniform(lb1, ub1)
                  )
def set_tungstic_acid_moles(X_tam):
    F.unit.M200.specifications[0].args[0] = X_tam

##Cobalt acetate catalyst ratio
lb2 = 0.003
ub2 = 0.015

 
@model.parameter(name='Cobalt acetate moles',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb2, ub2)
                  )
def set_cobalt_acetate_moles(X_cam):
    F.unit.R300.specifications[0].args[0] = X_cam


## Dihydroxylation reaction conversion
lb3 = 0.90
ub3 = 0.99
@model.parameter(name='Dihydroxylation reaction conversion',
                  element=F.R200,
                  kind='coupled',
                  distribution=shape.Uniform(lb3, ub3))
def set_dihydroxylation_conversion(X_dih):
    F.unit.R200.X_dih = X_dih


## Oxidative cleavage reaction conversion
F.R300.X_oxidativecleavage = 0.8
lb4 = 0.80
ub4 = 0.90
@model.parameter(name='Oxidative cleavage reaction conversion',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb4, ub4))
def set_X_oxidativecleavage_conversion(X_oxidativecleavage):
    F.unit.R300.X_oxidativecleavage = X_oxidativecleavage


##Oxidative cleavage reaction conversion
F.R300.X_decarboxylation = 0.2
lb5 = 0.1
ub5 = 0.3
@model.parameter(name='Decarboxylation reaction conversion',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb5, ub5))
def set_X_decarboxylation_conversion(X_decarboxylation):
    F.unit.R300.X_decarboxylation = X_decarboxylation


##Oxidative cleavage reaction conversion
F.R300.X_side_rxn = 0.8
lb6 = 0.1
ub6 = 0.3
@model.parameter(name='Side reaction conversion',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb6, ub6))
def set_X_side_rxn_conversion(X_side_rxn):
    F.unit.R300.X_side_rxn = X_side_rxn


lb7 = 35000*0.95
ub7 = 35000*1.05
@model.parameter(name='Feedstock_capacity',
                  element=F.stream.crude_vegetable_oil,
                  kind='coupled',
                  distribution=shape.Uniform(lb7, ub7)
                  )
def set_feedstock_input_flow(Cap):
    F.stream.crude_vegetable_oil.F_mass = Cap


lb9 = 0.3
ub9 = 2.5
@model.parameter(name='Extraction solvent mass',
                  element=F.M602,
                  kind='coupled',
                  distribution=shape.Uniform(lb9, ub9)
                  )
def set_solvent_mass(X_sm):
    F.unit.M602.specifications[0].args[0] = X_sm

#amaze
lb10 = 0.08
ub10 = 0.15


@model.parameter(name='System IRR',
                  element=azelaic_acid_tea.IRR,
                  kind='cost',
                  distribution=shape.Uniform(lb10, ub10)
                  )
def set_irr(X_irr):
    azelaic_acid_tea.IRR = X_irr

lb11 = 0.80
ub11 = 0.90
@model.parameter(name='Turbogen efficiency',
                  element=F.BT901,
                  kind='coupled',
                  distribution=shape.Uniform(lb11, ub11)
                  )
def set_tubeff(X_beff):
    F.BT901.turbogenerator_efficiency = X_beff


lb14 = 3 #current state
ub14 = 10 #optimistic 
@model.parameter(name='Tungstic acid reusability',
                  element=F.M200,
                  kind='coupled',
                  distribution=shape.Uniform(lb14, ub14)
                  )
def set_cta(cycles_of_reuse):
    F.unit.M200.specifications[0].args[1] = cycles_of_reuse

lb15 = 3 #current state assumed to be same as TA
ub15 = 10 #optimistic state randomly assigned
@model.parameter(name='Cobalt acetate reusability',
                  element=F.R300,
                  kind='coupled',
                  distribution=shape.Uniform(lb15, ub15)
                  )
def set_cca(cycles_of_reuse):
    F.unit.R300.specifications[0].args[0] = cycles_of_reuse



N_samples = 10
rule = 'L' # For Latin-Hypercube sampling
np.random.seed(1234) # For consistent results
samples = model.sample(N_samples, rule)
model.load_samples(samples)
model.evaluate(
    notify=20 # Also print elapsed Time after 100 simulations
              )
model.show()


##PLOTS


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

# b = sum([get_feedstock_GWP(),get_other_materials_impact(),
# get_total_non_BT_direct_emissions_GWP(),   get_heating_demand_GWP(),get_cooling_demand_GWP(),get_electricity_demand_non_cooling_GWP()])
# a = sum ([get_feedstock_GWP(),get_other_materials_impact(), get_ng_GWP(),get_electricity_use_GWP(),get_total_direct_emissions_GWP()])
# #code to get stacked bar plot

# df_breakdown = bst.UnitGroup.df_from_groups(aa_baseline_groups,
# fraction = True, scale_fractions_to_positive_values = True)
# stacked_bar_plot(dataframe = df_breakdown)
# df_breakdown.rename(columns={'Electricity consumption': 'Net electricity demand'})
# df_breakdown.rename(index = {'100' : 'Biodiesel production','200' : 'Dihydroxylation process', '300' : 'Oxidative cleavage process', '400' : 'Catalyst recovery process', '500' : 'Pelargonic acid recovery','600' : 'Azelaic acid recovery','900':'Facilities (boiler, cooling utilities, and air distribution utlity)'},inplace = True)


# a = stacked_bar_plot(dataframe = df_breakdown,y_ticks = (0,10,20,30,40,50,60,70,80,90,100),y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
#                 y_units = "%",colors=contourplots.utils.defaults_dict['colors']['Guest_Group_TEA_Breakdown'],n_minor_ticks=4,
#                   fig_height=5.5*1.1777*0.94*1.0975)

# df_breakdown.rename(index = {'Facilities (boiler, cooling utilities, and air distribution utlity)':'Facilities (boilerturbogenerator, cooling utilities and, air distribution utility)'},inplace = True)

# a = stacked_bar_plot(dataframe = df_breakdown,y_ticks = (0,10,20,30,40,50,60,70,80,90,100),y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
#                 y_units = "%",colors=contourplots.utils.defaults_dict['colors']['Guest_Group_TEA_Breakdown'],n_minor_ticks=4,
#                   fig_height=5.5*1.1777*0.94*1.0975)

# b = sum([get_feedstock_GWP(),get_other_materials_impact(),
# get_total_non_BT_direct_emissions_GWP(),   get_heating_demand_GWP(),get_cooling_demand_GWP(),get_electricity_demand_non_cooling_GWP()])
# c = [get_feedstock_GWP()/b,get_other_materials_impact()/b,
# get_total_non_BT_direct_emissions_GWP()/b,   get_heating_demand_GWP()/b,get_cooling_demand_GWP()/b,get_electricity_demand_non_cooling_GWP()/b]

# data = [[0.110], [0.115], [0.14]]
# df = pd.DataFrame(data, columns=['Age'])

# bst.plots.plot_spearman_1d(rhos = (0.0235,-0.0202,-0.0373,-0.0129,-0.0425,0.0141),
# index = ['Dihydroxylation reaction conversion (Uniform; 0.90 to 0.99)',
# 'Oxidative cleavage reaction conversion (Uniform; 0.80 to 0.90)',
# 'Decarboxylation reaction conversion (Uniform; 0.1 to 0.3)',
# 'Side reaction conversion (Uniform; 0.1 to 0.3)',
# 'Boilerturbogenerator efficiency (Uniform; 0.80 to 0.90)',
# 'Amount of tungstic acid (Uniform; 0.006 to 0.15)'],
# name = "Total capital investment",
# color = [GG_colors.blue.RGBn])

# bst.plots.plot_spearman_1d(rhos = (0.0236,-0.0202,-0.0373,-0.0128,-0.0406,-0.113),
# index = ['Dihydroxylation reaction conversion (Uniform; 0.90 to 0.99)',
# 'Oxidative cleavage reaction conversion (Uniform; 0.80 to 0.90)',
# 'Decarboxylation reaction conversion (Uniform; 0.1 to 0.3)',
# 'Side reaction conversion (Uniform; 0.1 to 0.3)',
# 'Boilerturbogenerator efficiency (Uniform; 0.80 to 0.90)',
# 'Amount of tungstic acid (Uniform; 0.006 to 0.15)'],
# name = "Internal rate of return",
# color = [GG_colors.blue.RGBn])


# bst.plots.plot_spearman_1d(rhos = (0.0243,-0.0194,-0.0372,-0.0121,-0.0269,-0.0335),
# index = ['Dihydroxylation reaction conversion (Uniform; 0.90 to 0.99)',
# 'Oxidative cleavage reaction conversion (Uniform; 0.80 to 0.90)',
# 'Decarboxylation reaction conversion (Uniform; 0.1 to 0.3)',
# 'Side reaction conversion (Uniform; 0.1 to 0.3)',
# 'Boilerturbogenerator efficiency (Uniform; 0.80 to 0.90)',
# 'Amount of tungstic acid (Uniform; 0.006 to 0.15)'],
# name = "Total product yield [10^6 kg/yr]",
# color = [GG_colors.blue.RGBn])

# bst.plots.plot_spearman_1d(rhos = (0.0243,-0.0194,-0.0372,-0.0121,-0.0338,-0.0269),
# index = ['Dihydroxylation reaction conversion (Uniform; 0.90 to 0.99)',
# 'Oxidative cleavage reaction conversion (Uniform; 0.80 to 0.90)',
# 'Decarboxylation reaction conversion (Uniform; 0.1 to 0.3)',
# 'Side reaction conversion (Uniform; 0.1 to 0.3)',
# 'Boilerturbogenerator efficiency (Uniform; 0.80 to 0.90)',
# 'Amount of tungstic acid (Uniform; 0.006 to 0.15)'],
# name = "Purity [%]",
# color = [GG_colors.blue.RGBn])

# bst.plots.plot_spearman_1d(rhos = (0.0232,-0.0194,-0.0376,-0.0129,-0.0383,0.447),
# index = ['Dihydroxylation reaction conversion (Uniform; 0.90 to 0.99)',
# 'Oxidative cleavage reaction conversion (Uniform; 0.80 to 0.90)',
# 'Decarboxylation reaction conversion (Uniform; 0.1 to 0.3)',
# 'Side reaction conversion (Uniform; 0.1 to 0.3)',
# 'Boilerturbogenerator efficiency (Uniform; 0.80 to 0.90)',
# 'Amount of tungstic acid (Uniform; 0.006 to 0.15)'],
# name = "Annual material cost [10^6 $/yr]",
# color = [GG_colors.blue.RGBn])

# bst.plots.plot_spearman_1d(rhos = (0.0228,-0.0183,-0.0368,-0.0116,-0.0465,0.24),
# index = ['Dihydroxylation reaction conversion (Uniform; 0.90 to 0.99)',
# 'Oxidative cleavage reaction conversion (Uniform; 0.80 to 0.90)',
# 'Decarboxylation reaction conversion (Uniform; 0.1 to 0.3)',
# 'Side reaction conversion (Uniform; 0.1 to 0.3)',
# 'Boilerturbogenerator efficiency (Uniform; 0.80 to 0.90)',
# 'Amount of tungstic acid (Uniform; 0.006 to 0.15)'],
# name = "Annual operating cost [10^6 $/yr]",
# color = [GG_colors.blue.RGBn])

# bst.plots.plot_spearman_1d(rhos = (0.0224,-0.0185,-0.038,-0.0122,-0.0266,-0.0207),
# index = ['Dihydroxylation reaction conversion (Uniform; 0.90 to 0.99)',
# 'Oxidative cleavage reaction conversion (Uniform; 0.80 to 0.90)',
# 'Decarboxylation reaction conversion (Uniform; 0.1 to 0.3)',
# 'Side reaction conversion (Uniform; 0.1 to 0.3)',
# 'Boilerturbogenerator efficiency (Uniform; 0.80 to 0.90)',
# 'Amount of tungstic acid (Uniform; 0.006 to 0.15)'],
# name = "Heating demand [KJ/yr]",
# color = [GG_colors.blue.RGBn])

# bst.plots.plot_spearman_1d(rhos = (0.0232,-0.0194,-0.0375,-0.0129,-0.337,-0.0207),
# index = ['Dihydroxylation reaction conversion (Uniform; 0.90 to 0.99)',
# 'Oxidative cleavage reaction conversion (Uniform; 0.80 to 0.90)',
# 'Decarboxylation reaction conversion (Uniform; 0.1 to 0.3)',
# 'Side reaction conversion (Uniform; 0.1 to 0.3)',
# 'Boilerturbogenerator efficiency (Uniform; 0.80 to 0.90)',
# 'Amount of tungstic acid (Uniform; 0.006 to 0.15)'],
# name = "Cooling demand [KJ/yr]",
# color = [GG_colors.blue.RGBn])

#GWP uncertainity
# -10.1013715391653,-10.1013715391653,-9.2194383948433,-10.0992571308872,-10.0982008385692,-9.90186159253726,-9.4548428868478,-10.0950355050746,-9.42430561268723,-8.33361692388758,-8.31372653621201,-10.0908230821497,-10.0897714418251,-10.0887203522706,-10.0876698047761,-10.0866197915065,-10.1067083182161,-9.84009484194162,-10.1045898705105,-9.47812570724711,-9.88217752134524,-10.0323101119251,-10.1003594950231,-9.68667960690314,-9.93553136755859,-9.23125591180045,-8.02104220210868,-9.37684640998147,-10.0940313318843,-8.23718114409215,-8.75282217013377,-7.29193486618047,-7.91590047987156,-9.44722783698728,-10.0469311342924,-9.69502046600994,-10.1067068534771,-7.87964827991652,-8.37969349582928,-9.95887303530911,-10.1024718300382,-9.85613400908745,-10.1003576854123,-9.61717065072607,-7.14749857290316,-9.661466697058,-9.44246155951915,-9.83833794700436,-10.0940292443226,-8.60938351610883,-9.65695864596739,-10.0908727329402,-9.50574023541691,-10.0887714349952,-10.0670164516835,-7.38443753283331,-10.1066792875445,-9.87090032696554,-10.1045596471412,-9.65664127358544,-10.1024421550502,-10.1013842586725,-9.48769456071162,-10.099270196704,-9.19945226985314,-8.49180672611547,-10.0211630459577,-10.0950491461355,-9.46588236409236,-10.0929421820126,-10.09188957704,-10.0908373968771,-10.0897860566266,-8.21165803150189,-10.0876846416924,-10.0866347827689,-10.1066949695785,-10.1056351591784,-10.1045758901031,-9.00838175345125,-9.80149725456785,-10.0230681514018,-10.1003442378284,-10.0992877809441,-9.51928730579778,-8.78526207195747,-10.0458468598402,-8.2795439863606,-10.0940142480005,-10.0929613261427,-10.0362091963592,-9.55493397323389,-9.31606612753456,-10.0887552206163,-10.0877050775955,-10.08665546629,-10.1085212465847,-9.69986445374034,-10.1064022850433,-10.1034755311526,-10.1024167474912,-10.1013585330618,-8.98846809933749,-10.0992438121833,-10.0981873716628,-10.0534437517036,-10.0960762851052,-9.97244396777282,-8.91083365655507,-10.0929140282936,-10.0918611232773,-10.0908086448998,-10.0897570081094,-10.0887057817079,-10.087655098407,-10.0866049500592,-10.1084724175319,-10.1074125329221,-10.1063531911362,-10.1052943818313,-10.1042361210697,-10.1031784272942,-10.1021213153454,-10.1010647946832,-10.1000088696347,-10.0989535419332,-8.83893166352619,-10.0968445071972,-10.0957909519627,-9.50771475439876,-10.0936855001503,-10.092633522642,-10.0915823846749,-10.0905316551282,-10.0894814667329,-9.78210964342588,-10.1067016539007,-10.1056419894841,-9.45593186860443,-10.1035242522871,-10.1024662036572,-9.99790820955933,-10.1003518119296,-9.68622868606781,-10.0982397601285,-9.81091585423554,-9.57842538773895,-8.31211536073528,-10.0940226955943,-9.42643641241228,-8.59384622525658,-9.60407205197404,-10.0898148422959,-10.0887642913297,-9.53046923428395,-10.08666477529,-10.1067119584365,-8.80779825039237,-9.70614325824136,-7.44796002491679,-8.65323185217118,-8.73242377939026,-9.31954049391593,-9.48851665078555,-10.0734897608802,-9.61038457270073,-9.59648784447914,-9.68091340096497,-9.88911732488145,-9.5519845493472,-10.0220335991904,-7.11487201521683,-8.54241617410398,-9.63841455501212,-8.53237518145583,-9.46302161519456,-10.1083361474529,-10.1072758952252,-10.1062161878311,-10.1051570149521,-10.104098392487,-10.1030403390682,-10.1019828693642,-10.100925992875,-10.0998697139499,-10.098814034251,-10.0977589517946,-10.0967044633527,-9.66088235728334,-10.0945971856523,-10.093544441973,-9.71652961373123,-10.0914405972143,-10.0903895310149,-10.0893390077208,-10.0882890191555,-10.106711959225,-10.047243414349,-10.1045937446872,-10.1035354353829,-8.09100826301391,-9.10500328323839,-10.0299035700308,-9.38559254805554,-9.97459744735152,-8.5822642675017,-8.93857355173676,-9.93330311534331,-9.44002711143028,-10.0828680662542,-9.53986023251582,-9.78694162636132,-9.95051738582172,-10.0530496257721,-8.62974739089554,-9.82993428622228,-10.1085113668546,-10.1074516065036,-10.1063923882612,-10.1053337018121,-10.1042755630989,-10.1032179907521,-10.1021609994672,-8.7954811615957,-10.1000486296161,-9.86141772577289,-10.0979387745008,-10.0950167966048,-9.67391690854324,-10.0929091416979,-9.99461137172568,-9.39616225915678,-10.0897519434041,-10.0887006713666,-10.08764994271,-10.0865997492496,-10.1084153090709,-10.1073552590016,-10.106295752644,-9.48706340403406,-10.1041782766894,-10.1031204205268,-10.1020631470474,-10.1010064658197,-10.0999503810887,-9.78162001055098,-10.0978399630671,-10.0967856655452,-9.399818677219,-10.094678733927,-10.0936261780141,-10.0925740475272,-10.0915227574226,-10.0904718765606,-10.089421537664,-10.0883717325402,-10.1067051506836,-9.93768153604116,-10.1045865309424,-8.90857702063655,-10.1008386302955,-9.60281147149444,-10.100355788139,-10.0992995548212,-10.098243914241,-10.0971888681246,-10.0067561087303,-9.20382043609539,-9.42159236504598,-10.0212197213012,-10.0919223702417,-9.3524155436089,-10.0898197538716,-10.0887692831836,-8.93785243973378,-8.69306851260585,-10.1066938040774,-7.62945367361295,-10.1045745938194,-10.1035158302153,-10.1024576137425,-8.7910631453093,-10.1003428454075,-10.099286363887,-10.0982304767659,-10.0971751857189,-10.0961204889399,-9.5387686589409,-9.49045223166861,-10.0929598734822,-10.0919075023736,-9.64534451667581,-10.0873294575775,-10.0887537239751,-8.14463137996709,-9.68631788422778,-10.1067084814421,-9.93114970184612,-9.78936082454083,-10.0581159982637,-9.29566786742395,-9.88132831084859,-10.1003596750896,-9.48794282470274,-9.03224448258948,-9.86974604030765,-8.63399871302036,-7.25762288703698,-10.0420417916614,-9.76874949789443,-9.94547885243165,-9.797396020573,-8.60122561858394,-10.0078625499769,-9.63035840843117,-9.97944466761917,-10.1067004371855,-8.96267701202536,-9.95178521063257,-10.10352293028,-7.95235274010922,-10.1014072699281,-9.6443582709614,-9.79274261479871,-10.0982382117105,-10.0971830540844,-10.0961284897628,-9.90937641920048,-9.11635198193249,-9.76702025171744,-10.0919160150519,-9.30685951252538,-10.0898131819277,-10.0887626055463,-10.0572914160656,-10.0866630610061,-10.1066762219857,-10.1056160876335,-10.1045564966877,-10.1034974388075,-9.28057943044121,-10.1013809552467,-9.06441711870125,-10.0992667833122,-10.0982106102655,-10.0971550350766,-10.0961000558084,-10.0950456692366,-10.0939918692681,-10.0929386503414,-8.69372945171929,-10.0908337170264,-9.198225906687,-10.0887313187285,-10.0876808864584,-9.83739377950547,-10.1083994378095,-10.1073393448063,-10.1062797958033,-9.77900040959522,-9.38373345717851,-10.1031042810983,-10.1020469654945,-10.1009902423192,-9.75027193665981,-10.0988785422399,-10.0978236104546,-10.0967692719373,-10.0957155205775,-10.0946623508876,-10.0936097552892,-10.0925575853186,-10.0915062559328,-10.0904553359684,-9.56334916981733,-10.0883550435491,-10.1083410874127,-10.1072808469703,-10.1062211513322,-10.1051619901304,-10.104103379303,-10.1030453374231,-10.1019878791967,-10.1009310141542,-10.0998747466037,-10.0988190782472,-10.0977640070831,-10.0967095297968,-10.0956556404753,-10.0946023334421,-10.0935496011996,-10.0924972952809,-10.0914458306048,-10.0903947760366,-9.47299817319375,-10.0882942022085,-10.1067113284028,-9.58469096872724,-9.99351035445146,-10.0387219955239,-9.89116083619445,-9.68623654472846,-9.34448426991256,-9.60396242394826,-10.098251490427,-9.16733386921167,-8.39831235686112,-7.03964076236946,-8.51325532550932,-9.72241142890889,-9.97313204172971,-8.57307841341697,-9.4172739499261,-8.16039100929687,-9.8170979047498,-9.78128812034179,-10.1084439430168,-10.1073839735876,-10.1063245474627,-10.1052656543404,-10.1042073101496,-10.1031495334769,-10.1020923390602,-10.1010357364095,-9.57840007855023,-10.0989242539786,-8.7638210788091,-10.0968150453666,-10.0957614099126,-10.0947083554744,-10.0936558745373,-10.0926038186083,-10.091552602701,-10.0905017956226,-10.0894515301177,-9.8241125635556,-10.1084704596275,-10.1074105690579,-10.1063512212929,-9.54307371171521,-10.1042340691437,-10.103176369026,-9.97078055901327,-10.101062707283,-10.1000067759222,-10.0989514419117,-10.0978967033373,-10.0968425569239,-10.0957889965945,-9.78533572323712,-10.093683570566,-10.0926315882428,-10.0915804454874,-10.090529711162,-9.39634719466098,-10.0884297669358,-10.1066985178536,-8.62981226127921,-10.1045795409718,-9.18276177744284,-10.1024627260025,-10.059401436211,-7.936584579311,-10.099291731348,-7.32396103415071,-9.74350940081292,-10.096126013553,-10.0950719956326,-10.0940185620935,-10.0929657073326,-9.54460031472344,-8.91619602168291,-9.79051763232529,-8.84783307470962,-10.0877097076941,-10.0866601591785,-10.1067043901164,-9.31826881340713,-9.72042916607917,-9.90031063999667,-8.63920239335042,-9.13134185540263,-10.0777700294052,-9.76803226542033,-10.0424736269491,-10.097187853732,-9.43868667595746,-10.0194109916117,-9.87954177803217,-7.59470895202831,-10.0919211862627,-9.90745491235736,-9.23996994588226,-10.0887680231146,-10.0877180694992,-8.69044179847352,-10.106707601452,-8.00195692495965,-9.34876816482353,-9.40706371194132,-8.04587643535313,-10.1014152816702,-8.7785018793937,-9.65836216958238,-8.87491185333263,-9.41679041198097,-10.0961372664663,-8.25140810962097,-7.52402971749675,-8.1854365761341,-10.0919252754819,-8.97531430844349,-8.94755086578418,-8.55223207552086,-8.70170907096305,-9.15927672419556,-10.1067118791971,-9.67979144619107,-9.35636336868988,-9.64814455773682,-9.83678421757592,-9.12838718878588,-9.87998678928768,-9.44124577222083,-9.01554718654439,-9.01938767893353,-8.1796116851078,-9.67647296366382,-9.97443920772357,-8.00522890358712,-9.96448990166881,-9.39747375741478,-9.7333434465068,-10.0329759550093,-9.6556892608019,-9.85331093339863,-10.1066598060067,-10.105599448562,-10.1045396358415,-10.1034803576139,-10.102421629723,-10.1013634707606,-9.73230442767326,-9.94620053976751,-10.0981925043851,-10.0971367181439,-10.0960815290431,-10.0950269339047,-10.0939729265919,-10.092919501557,-10.0918666512445,-9.63013086598198,-10.089762620762,-10.0887114474123,-10.0876608168878,-10.0866107209287,-10.1083831131808,-10.107322977402,-10.1062633858183,-10.0552801335509,-10.1041458138169,-9.63934489295993,-10.1020304572795,-10.100973693136,-10.0999175259449,-10.0988619574028,-10.0978069856355,-10.0967526072362,-10.0956988162633,-10.0946456070964,-10.0935929722119,-10.0925407631805,-10.0914893948883,-10.0904384362617,-10.0893880199677,-10.0883381378303,-10.1067116171283,-10.100320236598,-9.88441082073349,-10.0935086999914,-10.0095465520441,-9.94983591657836,-10.0237584720855,-7.05475468595305,-8.82068844795433,-9.72028394700412,-9.20716868171098,-9.77727959645839,-9.64748608645761,-7.81158096740534,-7.48241454198151,-9.7315007226836,-9.851308581376,-9.01125423713557,-9.98873753480142,-10.0866794838852,-10.1084017463831,-10.1073416595576,-10.1062821166388,-10.10522310745,-10.1041646477905,-10.1031067562824,-10.1020494476585,-10.1009927314237,-10.0999366118816,-10.0988810907546,-10.097826166072,-10.0967718345566,-9.26723472410015,-10.0946648171077,-10.0936122275884,-10.0925600636544,-10.0915087402579,-10.0904578262777,-10.0894074543948,-10.0883576164946,-10.1084057652903,-10.1073456892419,-10.1062861570965,-10.1052271585337,-10.1041687094833,-10.1031108285269,-10.1020535303981,-10.1009968245494,-9.38227496301139,-10.0988851104788,-10.0978301953855,-10.0967758734471,-10.095722138598,-10.094668985254,-9.63147136533709,-10.0925641907948,-10.091512877093,-10.0904619727358,-10.0894116104331,-10.0883617820668,-10.108484673786,-10.1074248272414,-9.95661312902122,-10.105306733183,-10.1042485095602,-10.1031908528462,-10.1021337777073,-10.1010772936294,-10.1000214049972,-10.0989661134896,-10.097911417148,-10.0968573127572,-10.0958037942179,-10.0947508560339,-10.093698490696,-10.0926465496951,-9.39592085307714,-10.0905446644388,-10.0894945112698,-8.50837418613987,-10.1067080439485,-7.4946802130914,-10.104589517097,-7.40763404569427,-10.1024731204182,-9.55431563740985,-7.97140743499863,-9.64806798580413,-9.03515051127575,-9.96255834859906,-8.31020524060118,-10.0950839631863,-9.02270998646968,-7.35651987532309,-10.0919259029568,-9.2271979622208,-9.44931011481172,-9.86396649917473,-10.0405823326348,-9.94086447196181,-10.1083424617662,-10.1072822246308,-10.1062225322608,-10.105163374347,-9.59751793218019,-10.103046660734,-10.1019892052485,-10.1009323429164,-10.0998760780864,-10.0988204124067,-10.0977653439357,-10.0967108694439,-10.0956569827968,-10.0946036784648,-10.0935509489208,-10.0924986456782,-10.0914471836689,-10.0903961317852,-10.0893456227418,-10.0882956483081,-10.1067122027065,-9.99111357094993,-7.7194781730707,-7.01876359013512,-10.0032641362761,-9.986613064899,-8.74321875001099,-6.81140037670863,-10.0213187951918,-9.96791429139757,-9.51583618882525,-9.60987980715764,-7.89805019145659,-9.88123743846721,-9.16691002752804,-9.53928824106156,-9.77578498632759,-8.05616279592113,-9.67599383011282,-9.66221532076022,-10.106707567263,-9.81781538592796,-10.1045890757054,-10.1035306320253,-9.54188353481925,-9.92389806377183,-10.100358620425,-9.0695221528672,-9.90919119512314,-10.0971918523557,-10.0324165265978,-8.8502955919518,-8.15715979442423,-9.56974528940284,-9.63424797542729,-10.0887984457418,-7.69062111628673,-10.0878244059796,-9.92297081470712,-9.79224667967573,-10.1067112143077,-10.0830414944165,-9.88498762939694,-9.91039951308215,-9.18911552390902,-9.18083882567148,-8.29954092859631,-9.46627709144713,-9.70665672982529,-8.60601908694239,-10.0250787169309,-9.88543013078826,-10.0303916935514,-9.33672243421223,-9.63364681106933,-9.47472028318356,-9.84086398806365,-9.95801485578523,-10.0138099993531,-8.52266379239452,-10.1084937572781,-10.1074339394629,-8.53041181822731,-10.1053157244364,-10.1042575280157,-9.92282378923505,-8.47247127269113,-10.1010861641631,-10.100030300576,-10.0989750340519,-10.0979203625945,-10.0968662830224,-10.0958127891978,-10.0947598756132,-10.0937075347378,-9.21874264490906,-10.0916044281374,-8.84967176365308,-10.0895034672866,-10.0884538697269,-10.1067026218964,-10.1056429806627,-9.41935144964126,-9.8306326059393,-9.94523990046721,-10.1014097882141,-8.76307465233301,-10.0992965675025,-10.0982408655207,-10.0971857583812,-10.0961312443013,-10.088242171945,-8.67880690403901,-7.98714882137875,-9.89275637170712,-9.52938047493363,-7.85868147403967,-10.0887655285159,-10.0877155351155,-9.91251496618722,-10.1066657194495,-10.1056054376929,-10.104545700192,-8.97463746753811,-10.0195608962103,-10.1013697041169,-10.1003122007776,-10.0992552902112,-10.0981989766562,-10.0971432617968,-10.0960881437877,-8.91132004024011,-10.0939796242102,-10.0929262676432,-10.0918734854036,-10.090821129054,-10.0897696135457,-10.036047569897,-10.0876679416879,-10.0866179124264,-10.1066727561288,-10.1056125710541,-10.1045529297366,-10.1034938217417,-10.1024352631448,-10.1013772724122,-8.71934838242953,-10.0992629875374,-10.0459122643963,-10.0971511413061,-10.0960961145778,-10.057565989411,-10.0939878321566,-10.0929345665153,-9.02834280179974,-10.090829556333,-8.3085214392523,-9.54867816135293,-10.0876765156184,-9.97068676532392,-10.1066993057059,-9.96554908604545,-10.0808763561392,-10.1035217576783,-10.102463656731,-9.92890142998472,-8.32644474607707,-10.0992927302427,-9.01060081649798,-10.0971817360654,-10.0961271485123,-9.41850775128292,-9.44124877057274,-10.0929668524184,-10.0919145876971,-8.97258912724618,-10.0898116955139,-9.62928471702504,-9.75463333057918,-10.0522171292151,-10.1067057698344,-9.37306956149885,-9.68262355104698,-10.1035286630435,-7.7532593331154,-9.08326432272125,-9.88292157757004,-10.042083998446,-9.98586376324786,-7.32252031691292,-10.0961349911164,-10.0950811368564,-8.85220189638011,-9.93324006109088,-9.44444699710576,-9.60191751203723,-9.525199500392,-10.0471086079441,-9.42290614530476,-10.0656814490491,-10.108411136003,-10.1073510744996,-10.1062915568074,-10.1052325725814,-10.1041741378756,-10.1031162711414,-10.1020589871413,-10.1010022954345,-10.099946200269,-10.0988907033893,-10.0978358028078,-10.0967814952999,-10.0957277747761,-10.0946746356771,-10.0936220705384,-10.0925699308247,-10.0915186315138,-10.0904677414449,-10.0894173933796,-10.0883675791288,-10.1066758555279,-10.1056157157294,-10.1045561193316,-10.1034970560756,-10.1024385418549,-10.1013805953021,-9.71412555853947,-9.8643681127517,-10.0982102564557,-10.0971546773885,-10.0960996943124,-10.095045303852,-10.0939915000625,-10.0929382773226,-10.0918856281113,-10.0908334040144,-7.7918972558776,-10.0887309288447,-10.0876804929182,-10.0866305904409,-10.1066872046304,-10.1056272479601,-10.1045678336358,-10.1035089512412,-9.21994740185997,-9.10827682888575,-10.1003355902091,-8.99334626840867,-7.60225238681905,-10.0971674336142,-10.0961126186809,-10.0950583954115,-10.0940047577926,-9.93113998064764,-9.34493969408984,-8.1621960786335,-8.54366454236248,-9.87505061071799,-9.44235268856542,-10.0866448519368,-10.1066765611606,-10.0498105238568,-10.0368835249889,-10.1034977885439,-10.1024392845596,-9.85324857686061,-10.1003239840759,-10.0992672214408,-10.0982110549195,-10.0971554861961,-9.09992778431009,-10.095046088417,-10.093992293627,-10.0929390798456,-10.0918864396164,-10.0062821370639,-9.62045559661795
#GWP breakdown
#Major feed impacts
#Feeds
#Natural gas
#aa_baseline_sys.get_material_impact(F.natural_gas,GWP)/aa_baseline_sys.get_mass_flow(azelaic_acid)
#Pelargonic acid 
# aa_baseline_sys.get_material_impact(F.pelargonic_acid_rich_fraction,GWP)/aa_baseline_sys.get_mass_flow(azelaic_acid)
#Fatty acid blend
# aa_baseline_sys.get_material_impact(F.fatty_acid_blend,GWP)/aa_baseline_sys.get_mass_flow(azelaic_acid)
#Hydrogen peroxide
# aa_baseline_sys.get_material_impact(F.fresh_HP,GWP)/aa_baseline_sys.get_mass_flow(azelaic_acid)
#Other feeds
# get_other_materials_impact() - 0.4404211826499017-0.2453338195419956-0.2971806321427208
#other products

#

