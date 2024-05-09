#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 16:43:50 2023
@author: lavanyakudli
"""
import numpy as np
import chaospy
from scipy import stats
from chaospy import distributions as shape

__all__ = ('prices_per_Kg',
           'product_prices_per_kg',
           'utility_prices',
           'prices_per_stream',
           'GWP_per_stream'
           )

#TODO: ask about heat transfer regen price

#%% Data for prices of materials and utility

#Freds Producers price index was used to update the prices of chemicals to 2023 prices[1]
#Index was chosen for the general category of all chemicals and allied products
PPI_1984 = 100 # Dec 1984 (Chemical manufacturing industry)[39]
PPI_2016 = 270.400 #Dec 2016
PPI_2021 = 349.906 #Dec 2021
PPI_2022 = 357.404 #Dec 2022
PPI_2023 = 347.123 #Sept 2023, latest PPI
PPI_1997 = 143.500 #Dec 1997
PPI_2018 = 293.000 #Dec 2018
PPI_utility_2016 = 136.200 #Dec 2016
PPI_utility_2017 = 139.60#Dec 2017
PPI_utility_2018 = 151.10#Dec 2018
PPI_utility_2019 = 141.10#Dec 2019
PPI_utility_2020 = 146.50#Dec 2020
PPI_utility_2021 = 171.89#Dec 2021
PPI_utility_2022 = 222.043#Dec 2022
PPI_utility_2017_2022 = [PPI_utility_2017,
                         PPI_utility_2018,
                         PPI_utility_2019,
                         PPI_utility_2020,
                         PPI_utility_2021,
                         PPI_utility_2022]
PPI_utility_2023 = 197.546 #Sept 2023 

ratio_2023from2021 = PPI_2023/PPI_2021
ratio_2023from2022 = PPI_2023/PPI_2022
ratio_2023from1997 = PPI_2023/PPI_1997
ratio_2023from2016 = PPI_2023/PPI_2016
ratio_2023from2018 = PPI_2023/PPI_2018
ratio_2023from1984 = PPI_2023/PPI_1984
ratio_utility_2023from2016 = PPI_utility_2023/PPI_utility_2016
ratio_utility_2023from2022 = PPI_utility_2023/PPI_utility_2022 
ratio_utility_2023from2017_2022_period_mean = PPI_utility_2023/np.mean(PPI_utility_2017_2022)

prices_per_Kg = {
        #2021 prices to 2023
        'HCl' : 0.13*ratio_2023from2021, # Catcost materials library [2]
        'Sodium_hydroxide': 0.915*ratio_2023from2021,#2021 price, Catcost materials library [2]
        'Citric_acid': 3.26*ratio_2023from2021,#Ranges [$/kg] between 1.60 - 4.93, 2021 price [2]
        'Hydrogen_peroxide': 1.505*ratio_2023from2021,#Based on 50% of hydrogen peroxide, range [$/kg] 1.21 - 1.80 [2]
        'Calcium_chloride': 0.83*ratio_2023from2021,#Based on 94-97% pure flakes,ranges [$/kg] 0.47 to 1.19 [2]             
        'Heptane':0.684*ratio_2023from2021,# [2]
        'Methanol':0.792*ratio_2023from2021,#[2]
        'Resin': 106, #2023 price, [3]      
        'Sodium_methoxide': 2.93*ratio_2023from2016,# [4]
        'biodiesel_sold_as_product': 1.35, #[37],[38] based on general biodiesel prices and mass of biodiesel
        'Oleic acid_1': 1.62*ratio_2023from2021,
        'Oleic acid_2': 1.60*ratio_2023from2021,
        'Price of ozone1': 1.03*ratio_2023from1984 #[39]
        }

#feedstock prices
#high oleic soybean oil
pounds_to_kg = 0.453
rbd_soybean_oil_pricing = {'2020_2021':0,#TODO: #[5]
                           '2022_2023':85.87/100,#[6]
                           '2023_2024':66.40/100} #[7]
premium_value_23 = 0.3 #[8]

#high oleic sunflower oil
price_differential = [0.1,0.18] #Typically, high oleic sunflower oil is priced 10%-18% higher than sunflower oil #[12]
commodity_sunflower_oil_pricing = {'2017_2018':1.20054,#$/kg #[6]
                                   '2018_2019':1.17216,#$/kg #[6]
                                   '2019_2020':1.43066,#$/kg #[6]
                                   '2020_2021':1.738,#$/kg #[6]
                                   '2021_2022':2.45058,#$/kg #[6]
                                   '2022_2023':1.936}#$/kg #[6]


feedstock_prices_per_Kg = {'HoSun_oil': commodity_sunflower_oil_pricing['2022_2023']*(1+np.average(price_differential))*ratio_2023from2022,
                           'HoySoy_oil': (rbd_soybean_oil_pricing['2023_2024']+premium_value_23)/pounds_to_kg
                                }

#lower heating value of veg oil derived fatty acid blends
LHV_of_high_oleic_soybean_oil_biodiesel = 37.7 #MJ/kg #[36]
LHV_of_sunflower_oil_biodiesel = 38.49 #MJ/kg

product_prices_per_kg = {        
        'C5_C9_fraction': 2.00*ratio_2023from2018,#2018 price for caproic acid (major component),1500–2500 USD/ton [13]
        'Crude_glycerol': 0.16*ratio_2023from2022,#Prices based on [14]
        'Crude_methanol':0.792*ratio_2023from2021,#[2]
        }



#Bulk prices for several chemicals were not available
#Laboratory prices for lab quantities were upgraded to bulk prices using correlation provided by [15]
#The correlation uses a bulk quantity of 60lbs (30 kgs) for estimating the bulk price [15]

def convert_lab_price_to_bulk_price(
                                    lab_quantity_in_kg = np.array([]), 
                                    lab_price_per_kg = np.array([]),
                                    bulk_quantity_in_kg = 30, #used by the paper [15]
                                    bulk_coeff = -0.75 #as mentioned in the paper [15]
                                    ):
    bulk_price = np.array([(lab_price_per_kg[i]*(bulk_quantity_in_kg/lab_quantity_in_kg[i])**(bulk_coeff))for i in range(len(lab_quantity_in_kg))])
    return np.average(bulk_price)

correlation_based_bulk_prices = {
        #2023 prices
        #catalysts
        'Cobalt_acetate': convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,
                                                                                        0.5,
                                                                                       ]),
                                                          lab_price_per_kg = np.array([379.58,#[16] GFS chemicals
                                                                                       735.42,#[17] BEantown chemicals
                                                                                       ]),
                                                          bulk_quantity_in_kg = 30, bulk_coeff = -0.75),
        
        
        'Tungstic_acid': convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,
                                                                                        0.5,
                                                                                        ]),
                                                          lab_price_per_kg = np.array([682.00,#[18] Sigma aldrich
                                                                                       838.00]),#[18] Sigma aldrich
                                                          bulk_quantity_in_kg = 30, bulk_coeff = -0.75),
        #products
        'Pelargonic_acid':convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.906]),#[18] Sigma aldrich
                                                          lab_price_per_kg = np.array([74.06 #[18] Sigma aldrich
                                                                                       ]),
                                                          bulk_quantity_in_kg = 30, bulk_coeff = -0.75),
        #products
        'Azelaic_acid': convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([2,1,1]),#[18] Sigma aldrich
                                                          lab_price_per_kg = np.array([77.00,224.00,112.00 #[18] Sigma aldrich
                                                                                       ]),
                                                          bulk_quantity_in_kg = 30, bulk_coeff = -0.75)}

    
    
#unit conversions
lb_to_kg = 0.453
cooling_tower_chem_density = 8.65 #lbs/gal [22]
boiler_chem_density = 9.10 #lb/ga[23]
cents_to_dollar = 100 #cents to dollar
tcf_to_cf = 1000 #thousand cubic feet to cubic feet
cf_to_m3 = 0.029 #cubic feet to cubic meter
kgal_to_L = 3.785412 *1000

natural_gas_prices_2017_to_2022 = [4.08,4.19,3.90,3.32,5.44,7.66] #Provided by EIA, Dollars per Thousand Cubic Feet
average_electricity_prices_2017_to_2022 = [6.88,6.92,6.81,6.67,7.18,8.45]#Provided by EIA, cents/Kwh
average_water_rates_2008_2012_2016 = [2.44,3.02,3.38]#$/kGal Provided by USDOE

utility_prices = {   
        'System_makeup_water':np.mean(average_water_rates_2008_2012_2016)*ratio_utility_2023from2016/kgal_to_L ,#[19],[20]
        'Natural_gas_price': 7.66*ratio_utility_2023from2022/(tcf_to_cf*cf_to_m3),#Industrial natural gas price provided by EIA for 2017-2022 period[21]
        'Cooling_tower_chemicals':519/(55*cooling_tower_chem_density*lb_to_kg),#2023 price,Based on [22], 55 gal drum costs 519$
        'Boiler_chems':1.216/(55*boiler_chem_density*lb_to_kg), #2023 price,Based on [23],55 gal drum costs 609$
        'Lime_boiler':0.11*ratio_2023from2021,#2021 price based on average of all lime values provided in [2]
        'Electricity':np.mean(average_electricity_prices_2017_to_2022)*ratio_utility_2023from2017_2022_period_mean/cents_to_dollar,#Industrial retail electricity price provided by EIA for 2017-2022 period[24]
        'Ash_disposal_price':-0.0318,#Biosteam default price
        'Heat_transfer_price': 0,
        'brine_disposal': -round(550/160.3,2) #density of the stream at baseline,
        #In practise the disposal cost varies from 10 to 1,000 $/m3 [32]
        }


prices_per_stream = { #input prices
                      'crude_vegetable_Hosun': feedstock_prices_per_Kg['HoSun_oil'],
                      'crude_vegetable_Hoysoy':feedstock_prices_per_Kg['HoySoy_oil'],
                      'base_for_saponification_of_FFA':prices_per_Kg['Sodium_hydroxide'],
                      'citricacid_for_degumming':prices_per_Kg['Citric_acid'],
                      'fresh_HP':prices_per_Kg['Hydrogen_peroxide'],
                      'fresh_tungsten_catalyst':correlation_based_bulk_prices['Tungstic_acid'],
                      'fresh_cobalt_catalyst_stream':correlation_based_bulk_prices['Cobalt_acetate'],
                      'sodium_hydroxide_for_cat_sep':prices_per_Kg['Sodium_hydroxide'],
                      'calcium_chloride_for_cat_sep':prices_per_Kg['Calcium_chloride'],
                      'conc_hydrochloric_acid': prices_per_Kg['HCl'],
                      'solvent_for_extraction':prices_per_Kg['Heptane'],
                      'polystyrene_based_catalyst':prices_per_Kg['Resin'],
                      'Liquid_HCl':prices_per_Kg['HCl'],
                      'NaOH':prices_per_Kg['Sodium_hydroxide'],
                      'methanol':prices_per_Kg['Methanol'],
                      'catalyst':0.25*prices_per_Kg['Sodium_methoxide'] + 0.75*prices_per_Kg['Methanol'],
                      'HCl': prices_per_Kg['HCl'],
                      'biodiesel_sold_as_product': prices_per_Kg['biodiesel_sold_as_product'],
                      
                      #product prices
                      'crude_glycerol':product_prices_per_kg['Crude_glycerol'],
                      'recovered_C5_to_C9_MCA_fraction':product_prices_per_kg['C5_C9_fraction'],
                      'pelargonic_acid_rich_fraction':correlation_based_bulk_prices['Pelargonic_acid'],
                      'fatty_acid_blend':{'HoySoy':0.88*ratio_2023from2021/LHV_of_high_oleic_soybean_oil_biodiesel,#2021 price for biodiesel as Matricia uses the complex heavy tail mixture for production of biofuels[2]
                                          'HoSun':0.88*ratio_2023from2021/LHV_of_sunflower_oil_biodiesel,#2021 price for biodiesel as Matricia uses the complex heavy tail mixture for production of biofuels[2]
                                          },
                      'azelaic_acid_product_stream':correlation_based_bulk_prices['Azelaic_acid'],
                      'crude_methanol':product_prices_per_kg['Crude_methanol'],
                      # 'malonic_acid_stream':correlation_based_bulk_prices['Malonic_acid'],
                      
                      #utility prices
                      'natural_gas': utility_prices['Natural_gas_price'],
                      'lime_boiler':utility_prices['Lime_boiler'],
                      'boiler_chems':utility_prices['Boiler_chems'],
                      'system_makeup_water':utility_prices['System_makeup_water'],
                      'ash_disposal':utility_prices['Ash_disposal_price'],
                      'cooling_tower_chemicals':utility_prices['Cooling_tower_chemicals'],
                      'brine_from_wastewater_treatment': utility_prices['brine_disposal']
                      }

#%% Global warming potential factors (kg CO2-eq/Kg)

lbsperMWh_to_KgsperKWh = 4.53592e-4  

GWP_factors_per_Kg = {  #inputs
                        'HoSun_oil':0.76, #Based on GWP of sunflower oil (incl. iLUC and biogenic CO2 uptake) [25]                                               
                        'HoySoy_oil': 0.485,#Based on GWP of commodity soybean oil as it is same for HO variety (excluding biogenic emissions), obtained value for GREET [35],[26]
                        'Citric_acid': 1.4821,#GHG-100 GREET 2022 [26]
                        'Hydrogen_peroxide':1.0813, #GHG-100 GREET 2022 [26]
                        'Tungstic_acid':68.5,#Value based on midpoint results available for tungsten carbide[27]
                        'Cobalt_acetate':8.8913, #Value based on cobalt nitrate [26]
                        'HCl':2.0486, #Value based on conc HCl production in the U.S. GREET[26]
                        'Resin': 0.1324, #Value for GHG-100 high impact polystyrene GREET [26]
                        'Sodium_hydroxide':0.3602, #GHG-100 value for sodium hydroxide, GREET 2022 [26]
                        'Sodium methoxide': 1.85, #Ecoinvent, production of sodium methoxide, RoW [28]
                        'Calcium chloride':0.017,#Based on calcium nitrate, GREET 2022 [26]
                        'Heptane':6.4337e-1,#Ecoinvent,heptane to generic market for solvent, organic [28]
                        'Methanol':2.5154, #Methanol from coal,GREET 2022[28],
                        }
Utility_GWP_factors =   {                      
                        #utilities
                        'Electricity': round(852.3*lbsperMWh_to_KgsperKWh,2), #Based on NAtional US average for electricty CO2 emissions [30]
                        'System_makeup_water':0.00035559,#Tap water production through conventional treatment [28]
                        'Natural_gas': 0.3983, #Based on GREET 2022, NA NG from Shale and Conventional Recovery [26]
                        'Cooling_tower_chemicals':0.2574, #Based on sodium bicarbonate which is a common corrosion inhibitor [26].
                        # Data for other cooling tower chemicals such as pH controllers and anti-scalants not found.
                        'Boiler_chems':1.5568, #Based on the production of sulfite which is commonly used as an oxygen scavenger in the boiler [26]
                        # Data for other boiler chemicals (neutralising amines) not found.
                        'Lime_boiler':1.2844, #Based on lime production from lime stone [26]
                        }                                     
                                                                                                                     

GWP_per_stream = { #input factors
                      'crude_vegetable_HoSun': GWP_factors_per_Kg['HoSun_oil'],
                      'crude_vegetable_Hoysoy':GWP_factors_per_Kg['HoySoy_oil'],
                      'base_for_saponification_of_FFA':GWP_factors_per_Kg['Sodium_hydroxide'],
                      'citricacid_for_degumming':GWP_factors_per_Kg['Citric_acid'],
                      'fresh_HP':GWP_factors_per_Kg['Hydrogen_peroxide'],
                      'fresh_tungsten_catalyst':GWP_factors_per_Kg['Tungstic_acid'],
                      'fresh_cobalt_catalyst_stream':GWP_factors_per_Kg['Cobalt_acetate'],
                      'sodium_hydroxide_for_cat_sep':GWP_factors_per_Kg['Sodium_hydroxide'],
                      'calcium_chloride_for_cat_sep':GWP_factors_per_Kg['Calcium chloride'],
                      'conc_hydrochloric_acid': GWP_factors_per_Kg['HCl'],
                      'solvent_for_extraction':GWP_factors_per_Kg['Heptane'],
                      'polystyrene_based_catalyst':GWP_factors_per_Kg['Resin'],
                      'Liquid_HCl':GWP_factors_per_Kg['HCl'],
                      
                      'NaOH':GWP_factors_per_Kg['Sodium_hydroxide'],
                      'methanol':GWP_factors_per_Kg['Methanol'],
                      'catalyst':0.25*GWP_factors_per_Kg['Sodium methoxide'] + 0.75*GWP_factors_per_Kg['Methanol'],
                      'HCl': GWP_factors_per_Kg['HCl'],                                           
                      #utility factors
                      'ash_disposal':0,
                          #TODO: 3.5141e+0,[28]
                      'cooling_tower_evaporation':0,
                      'natural_gas': Utility_GWP_factors['Natural_gas'],
                      'lime_boiler':Utility_GWP_factors['Lime_boiler'],
                      'boiler_chems':Utility_GWP_factors['Boiler_chems'],
                      'system_makeup_water':Utility_GWP_factors['System_makeup_water'],
                      'cooling_tower_chemicals':Utility_GWP_factors['Cooling_tower_chemicals'],
                      'brine_from_wastewater_treatment':0,
                      #TODO:0.0013,#sodium brine in U.S [26]
                      }


#%% Uncertainity analysis distributions
#Distributions for prices obtained from bulk correlations
#Catalyst prices
bulk_coeff_lower = -0.80
bulk_coeff_upper = -0.70
Tungstic_acid_upper_bound =  convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,0.5,]),#[18] Sigma aldrich
                                                  lab_price_per_kg = np.array([682.00,838.00]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff =bulk_coeff_upper)
Tungstic_acid_lower_bound =  convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,0.5,]),#[18] Sigma aldrich
                                                  lab_price_per_kg = np.array([682.00,838.00]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff = bulk_coeff_lower)

Cobalt_acetate_upper_bound = convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,0.5]),
                                                  lab_price_per_kg = np.array([379.58,#[16] GFS chemicals
                                                                               735.42,#[17] BEantown chemicals
                                                                               ]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff = bulk_coeff_upper)
Cobalt_acetate_lower_bound = convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,0.5]),
                                                  lab_price_per_kg = np.array([379.58,#[16] GFS chemicals
                                                                               735.42,#[17] BEantown chemicals
                                                                               ]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff = bulk_coeff_lower)

#for uniform distributions 10%, 25% and 50% ranges were used
#10% was used for distributions where the sources for prices/environmental factors were reliable (ICIS chemical database,US soybean export council, GREET, ecoinvent)
#50% was used for distributions where the prices/ environmental factors for the commodity was unavailable and had to be assumed based on the competitor commodity or in case of LCA via displacement of fossil fuel based routes
Lower_bound_factor_10_per = 0.90
Upper_bound_factor_10_per = 1.1
Lower_bound_factor_25_per = 0.75
Upper_bound_factor_25_per = 1.25
Lower_bound_factor_50_per = 0.5
Upper_bound_factor_50_per = 1.5


isolated_para_dists = {#Crude vegetable oil was priced based on the ratio of 10-year historical
                       #feedstock and oil price, along with an assumption of 2$ premium recieved by high oleic farmers
                       #over regular soybean feedstock
                       'Crude oil price' : {'HoySoy_oil':chaospy.Uniform(feedstock_prices_per_Kg['HoySoy_oil']*Lower_bound_factor_25_per, feedstock_prices_per_Kg['HoySoy_oil']*Upper_bound_factor_25_per)},
                       'Tungstic acid price': chaospy.Uniform(Tungstic_acid_lower_bound,Tungstic_acid_upper_bound),
                       'Cobalt acetate price':chaospy.Uniform(Cobalt_acetate_lower_bound,Cobalt_acetate_upper_bound),
                       'Heptane solvent price':chaospy.Uniform(prices_per_Kg['Heptane']*Lower_bound_factor_10_per,prices_per_Kg['Heptane']*Upper_bound_factor_10_per),
                       
                       'Hydrogen peroxide price':chaospy.Uniform(prices_per_Kg['Hydrogen_peroxide']*Lower_bound_factor_25_per,prices_per_Kg['Hydrogen_peroxide']*Upper_bound_factor_25_per),
                       
                       'Citric acid price':chaospy.Uniform(prices_per_Kg['Citric_acid']*Lower_bound_factor_10_per,prices_per_Kg['Citric_acid']*Upper_bound_factor_10_per),
                       
                       'Conc HCl price': chaospy.Uniform(prices_per_Kg['HCl']*Lower_bound_factor_10_per,prices_per_Kg['HCl']*Upper_bound_factor_10_per),
                       
                       'Calcium chloride price':chaospy.Uniform(prices_per_Kg['Calcium_chloride']*Lower_bound_factor_10_per,prices_per_Kg['Calcium_chloride']*Upper_bound_factor_10_per),
                       
                       'Hydrolysis resin price':chaospy.Uniform(prices_per_Kg['Resin']*Lower_bound_factor_50_per,prices_per_Kg['Resin']*Upper_bound_factor_50_per),
                       'Crude glycerol price': chaospy.Triangle(100*ratio_2023from2022/1000,160*ratio_2023from2022/1000,220*ratio_2023from2022/1000),#Based on historical crude glycerol prices provided in [14]
                       'Pelargonic acid price':chaospy.Uniform(correlation_based_bulk_prices['Pelargonic_acid']*Lower_bound_factor_25_per,correlation_based_bulk_prices['Pelargonic_acid']*Upper_bound_factor_25_per),
                       'C5_C9 fraction price':chaospy.Uniform(product_prices_per_kg['C5_C9_fraction']*Lower_bound_factor_25_per,product_prices_per_kg['C5_C9_fraction']*Upper_bound_factor_25_per),
                       
                       #Utility prices
                       'Natural gas price':chaospy.Uniform(utility_prices['Natural_gas_price']*Lower_bound_factor_10_per,utility_prices['Natural_gas_price']*Upper_bound_factor_10_per),
                       'Electricity price':chaospy.Triangle( lower = 0.0583,midpoint = 0.0637,upper = 0.069),#Based on historical prices provided in [14]
                       'Ash_disposal_price':chaospy.Uniform(utility_prices['Ash_disposal_price']*Lower_bound_factor_50_per,utility_prices['Ash_disposal_price']*Upper_bound_factor_50_per)}

#Limits obtained from literature
#Tungstic acid mole fraction wrt the unsaturations can be 
#between 0.03% to 3% [31]
TA_moles_lower_bound = 0.03/100 
TA_moles_upper_bound = 3/100
#tunggstic acid moles are typically used in the range 0.06% to 1.5% [31]
TA_moles_mid_point = 0.0078

#Cobalt acetate mole fraction wrt the moles of diol can be
#between 0.05% to 3% [31] 
CA_moles_lower_bound = 0.05/100
CA_moles_upper_bound = 3/100
#Cobalt acetate moles are typicallt used in the range of 0.3% to 1.5%
CA_moles_midpoint = 0.9/100

#Dihydroxylation reaction and Oxidative cleavage time can vary between 2 to 10 hours [31]
Dih_time_lower = 2
Dih_time_upper = 10
Ox_time_lower =  2
Ox_time_upper = 10

#Oxidative cleavage reactions
#Oxidative cleavage reaction primary
Oxcp_lower_bound = 0.80#TODO: random
Oxcp_upper_bound = 0.98#TODO: random

#Oxidative cleavage reaction of intermediates to form azelaic acid and pelargonic acid
Ox_lower_bound = 0.80 #TODO: random
Ox_upper_bound = 0.98 #TODO: random

#Dihydroxylation reaction conversion
Dih_lower_bound = 0.70 #Lowest reported
Dih_upper_bound = 0.99 #Highest reported

#Ranges for mimimum and maximum catalyst reuse were not available and were assumed
#TODO: find ranges 
TA_reuse_lower = 1
TA_reuse_upper = 10
CA_reuse_lower = 1
CA_reuse_upper = 10

#Turbogen efficiency and boiler efficiencies were varied uniformly based on ranges available in literature
Turbogen_eff_lower = 0.80 #[33] #TODO:random 
Turbogen_eff_upper = 0.90 #[33]
Boiler_eff_lower = 0.75 #[34]
Boiler_eff_upper = 0.85 #[34]

#Concentrations of hydrogen peroxide used in the process
Hp_conc_lower = 0.40#[31]
Hp_conc_upper = 0.70#[31]

Oil_compo_lower = 72
Oil_compo_upper = 85


coupled_para_dist = {'Tungstic acid moles': chaospy.Triangle(TA_moles_lower_bound,TA_moles_mid_point,TA_moles_upper_bound),
                     'Cobalt acetate moles': chaospy.Triangle(CA_moles_lower_bound,CA_moles_midpoint,CA_moles_upper_bound),
                     'Oxidative cleavage reaction conversion primary':chaospy.Uniform(Oxcp_lower_bound,Oxcp_upper_bound),
                     'Oxidative cleavage reaction conversion':chaospy.Uniform(Ox_lower_bound,Ox_upper_bound),
                     'Dihydroxylation reaction conversion':chaospy.Uniform(Dih_lower_bound,Dih_upper_bound),
                     'Dihydroxylation reaction time':chaospy.Uniform(Dih_time_lower ,Dih_time_upper),
                     'Oxidative reaction time':chaospy.Uniform(Ox_time_lower,Ox_time_upper),
                     'Turbogen efficiency':chaospy.Uniform(Turbogen_eff_lower,Turbogen_eff_upper),
                     'Boiler efficiency':chaospy.Uniform(Boiler_eff_lower,Boiler_eff_upper),
                     'Tungstic acid reusability':chaospy.Uniform(TA_reuse_lower,TA_reuse_upper),
                     'Cobalt acetate reusability':chaospy.Uniform(CA_reuse_lower,CA_reuse_upper),
                     'HP_concentration':chaospy.Uniform(Hp_conc_lower,Hp_conc_upper),
                     'Crude_oil_composition':chaospy.Uniform(Oil_compo_lower,Oil_compo_upper)
                     }

environmental_facs_dist = {'Oil_GWP':{'HoSun_oil':chaospy.Uniform(Lower_bound_factor_25_per*GWP_factors_per_Kg['HoSun_oil'],Upper_bound_factor_25_per*GWP_factors_per_Kg['HoSun_oil']),
                                      'HoySoy_oil':chaospy.Uniform(Lower_bound_factor_25_per*GWP_factors_per_Kg['HoySoy_oil'],Upper_bound_factor_25_per*GWP_factors_per_Kg['HoySoy_oil'])},
                           'Tungstic_acid':chaospy.Uniform(Lower_bound_factor_25_per*GWP_factors_per_Kg['Tungstic_acid'],Upper_bound_factor_25_per*GWP_factors_per_Kg['Tungstic_acid']),
                           'Cobalt_acetate':chaospy.Uniform(Lower_bound_factor_25_per*GWP_factors_per_Kg['Cobalt_acetate'],Upper_bound_factor_25_per*GWP_factors_per_Kg['Cobalt_acetate']),
                           'HCl':chaospy.Uniform(Lower_bound_factor_25_per*GWP_factors_per_Kg['HCl'],Upper_bound_factor_25_per*GWP_factors_per_Kg['HCl']),
                           'Resin':chaospy.Uniform(Lower_bound_factor_25_per*GWP_factors_per_Kg['Resin'],Upper_bound_factor_25_per*GWP_factors_per_Kg['Resin']),
                           'Sodium methoxide':chaospy.Uniform(Lower_bound_factor_25_per*GWP_factors_per_Kg['Sodium methoxide'],Upper_bound_factor_25_per*GWP_factors_per_Kg['Sodium methoxide']),
                           'Heptane':chaospy.Uniform(Lower_bound_factor_25_per*GWP_factors_per_Kg['Heptane'],Upper_bound_factor_25_per*GWP_factors_per_Kg['Heptane']),
                           'Cooling_tower_chemicals': chaospy.Uniform(Lower_bound_factor_25_per*Utility_GWP_factors['Cooling_tower_chemicals'],Upper_bound_factor_25_per*Utility_GWP_factors['Cooling_tower_chemicals']),
                           'Boiler_chems':chaospy.Uniform(Lower_bound_factor_25_per*Utility_GWP_factors['Boiler_chems'],Upper_bound_factor_25_per*Utility_GWP_factors['Boiler_chems']),
                           'Lime_boiler':chaospy.Uniform(Lower_bound_factor_25_per*Utility_GWP_factors['Lime_boiler'],Upper_bound_factor_25_per*Utility_GWP_factors['Lime_boiler']),
                           'Electricity':chaospy.Uniform(Lower_bound_factor_25_per*Utility_GWP_factors['Electricity'],Upper_bound_factor_25_per*Utility_GWP_factors['Electricity']),
                           'Hydrogen peroxide':chaospy.Uniform(Lower_bound_factor_25_per*GWP_factors_per_Kg['Hydrogen_peroxide'],Upper_bound_factor_25_per*GWP_factors_per_Kg['Hydrogen_peroxide']),
}
                            
#%%
#References:
#[1] https://fred.stlouisfed.org/release/tables?rid=46&eid=142872#snid=142874
#[2] https://catcost.chemcatbio.org/materials-library
#[3] https://www.chemicalbook.com/SupplierPriceList_EN.aspx?cbn=CB6292238&c=5kg#price
#[4] 10.1002/bbb.1640; Biofuels, Bioprod. Bioref. 10:299–315 (2016)
#[5] 
#[6]https://www.fastmarkets.com/insights/margins-rise-for-biodiesel-producers-using-soybean-oil-feedstock/#:~:text=Central%20IL%20RBD%20soybean%20oil%20prices%20averaged%2085.87,per%20pound%2C%20up%20from%2085.46%20cents%20per%20pound.
#[7]https://www.fastmarkets.com/insights/rbd-soybean-oil-trading/#:~:text=RBD%20currently%20carries%20a%2010.1%20cents%20per%20pound,remains%20in%20a%20downward%20trend%20since%20June%202022.
#[8]US High Oleic Soybeans & High Oleic Soybean Oil Sourcing Guide for International CustomersThird Edition – May 2023Prepared for US Soybean Export Council
#[10]https://www.pioneer.com/us/news-and-events/news/media-release/plenish-soybean-premiums.html
#[11]https://ussec.org/wp-content/uploads/2023/05/USSEC-High-Oleic-Sourcing-Guide_2023-Update-Third-Edition-Final-003.pdf
#[12]How High Oleic Sunflower Oil (HOSO) is Making a Statement in the QSR Industry by Suresh Thamotharan
#[13]https://doi.org/10.3390/en11061551
#[14]Yoel oilcane paper
#[15]https://www.researchgate.net/publication/271702798
#[16] GFS chemicals
#[17] Beantown chemicals
#[18] Sigma aldrich
#[19]https://www.energy.gov/sites/prod/files/2017/10/f38/water_wastewater_escalation_rate_study.pdf
#[20]https://fred.stlouisfed.org/series/PCU221221
#[21]https://www.eia.gov/dnav/ng/ng_pri_sum_dcu_nus_a.htm
#[22]https://www.coolingtowerchemicals.com/ProductDetails.asp?ProductCode=CTC1334NM
#[23]https://www.chemworld.com/Boiler-Chemical-p/chemworld-1394.htm?msclkid=59bcee3dabfc15d1fc428b13737412bb&utm_term=4585513243667217&utm_medium=cpc&utm_content=ChemWorld%203-29-17&utm_campaign=1234%20Products&utm_source=bing
#[24]https://www.eia.gov/electricity/data/browser/#/topic/7?agg=0,1&geo=g&endsec=vg&linechart=ELEC.PRICE.US-ALL.A~ELEC.PRICE.US-RES.A~ELEC.PRICE.US-COM.A~ELEC.PRICE.US-IND.A&columnchart=ELEC.PRICE.US-ALL.A~ELEC.PRICE.US-RES.A~ELEC.PRICE.US-COM.A~ELEC.PRICE.US-IND.A&map=ELEC.PRICE.US-ALL.A&freq=A&start=2011&end=2022&ctype=linechart&ltype=pin&rtype=s&maptype=0&rse=0&pin=
#[25]http://dx.doi.org/10.1016/j.jclepro.2014.10.011                           
#[26]GREET 2022
#[27]Ref: http://dx.doi.org/10.1016/j.jclepro.2017.02.184
#[28]Ecoinvent
#[29]https://doi.org/10.1021/acssuschemeng.2c05764
#[30]https://www.epa.gov/egrid/power-profiler#/
#[31] CONTINUOUS PROCESS FOR THE PRODUCTION OF DERVATIVES OF SATURATED CARBOXYLIC ACDS,US 9.272,975 B2
#[32]https://www.lenntech.com/processes/desalination/general/brine-faq.htm
#[33]https://www.energy.gov/sites/prod/files/2016/09/f33/CHP-Steam%20Turbine.pdf
#[34]https://iea-etsap.org/E-TechDS/PDF/I01-ind_boilers-GS-AD-gct.pdf#:~:text=The%20amount%20of%20input%20fuel%20depends%20on%20the,be%20improved%20by%20preventing%20and%2For%20recovering%20heat%20loss.

#[35]https://onlinelibrary-wiley-com.proxy2.library.illinois.edu/doi/full/10.1002/bbb.2462
#[36]https://doi-org.proxy2.library.illinois.edu/10.1111/j.1467-7652.2009.00408.x
#[37] https://www.aqua-calc.com/calculate/volume-to-weight
#[38]https://afdc.energy.gov/fuels/prices.html
#[39]https://fred.stlouisfed.org/series/PCU325325
#[40]Large-Scale Production and Application of Highly Concentrated Ozone HANS-PETER KLEIN, BBC Brown Boveri & Cie., CH-8050, Z0rich, Switzerland