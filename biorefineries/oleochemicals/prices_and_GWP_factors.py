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
ratio_utility_2023from2016 = PPI_utility_2023/PPI_utility_2016
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
        }

#feedstock prices
#high oleic soybean oil
pounds_to_kg = 0.453
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

#high oleic sunflower oil
price_differential = [0.1,0.18] #Typically, high oleic sunflower oil is priced 10%-18% higher than sunflower oil #[12]
commodity_sunflower_oil_pricing = {'2017_2018':1.20054,#$/kg #[6]
                                   '2018_2019':1.17216,#$/kg #[6]
                                   '2019_2020':1.43066,#$/kg #[6]
                                   '2020_2021':1.738,#$/kg #[6]
                                   '2021_2022':2.45058,#$/kg #[6]
                                   '2022_2023':1.936}#$/kg #[6]

feedstock_prices_per_Kg = {'HoSun_oil':  commodity_sunflower_oil_pricing['2022_2023']*(1+np.average(price_differential))*ratio_2023from2022,
                           'HoySoy_oil': (soyoil_and_premium['2022_2023']['soycommodityoil'] + soyoil_and_premium['2022_2023']['premium'])*ratio_2023from2022
                            }


product_prices_per_kg = {        
        'C5_C8_fraction': 2.00*ratio_2023from2018,#2018 price for caproic acid (major component),1500–2500 USD/ton [13]
        'Fatty_acid_blend':1.37*ratio_2023from2021,#2021 price for Stearic acid[2]
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
cents_to_dollar = 100
tcf_to_cf = 1000 #thousand cubic feet to cubic feet
cf_to_m3 = 0.029 #cubic feet to cubic meter
kgal_to_L = 3.785412 *1000

natural_gas_prices_2017_to_2022 = [4.08,4.19,3.90,3.32,5.44,7.66] #Provided by EIA, Dollars per Thousand Cubic Feet
average_electricity_prices_2017_to_2022 = [6.88,6.92,6.81,6.67,7.18,8.45]#Provided by EIA, cents/Kwh
average_water_rates_2008_2012_2016 = [2.44,3.02,3.38]#$/kGal Provided by USDOE

#TODO; should you also include bag house bag price
utility_prices = {   
        'System_makeup_water':np.mean(average_water_rates_2008_2012_2016)*ratio_utility_2023from2016/kgal_to_L ,#[19],[20]
        'Natural_gas_price': np.mean(natural_gas_prices_2017_to_2022)*ratio_utility_2023from2017_2022_period_mean/(tcf_to_cf*cf_to_m3),#Industrial natural gas price provided by EIA for 2017-2022 period[21]
        'Cooling_tower_chemicals':519/(55*cooling_tower_chem_density*lb_to_kg),#2023 price,Based on [22], 55 gal drum costs 519$
        'Boiler_chems':1.216/(55*boiler_chem_density*lb_to_kg), #2023 price,Based on [23],55 gal drum costs 609$
        'Lime_boiler':0.11*ratio_2023from2021,#2021 price based on average of all lime values provided in [2]
        'Electricity':np.mean(average_electricity_prices_2017_to_2022)*ratio_utility_2023from2017_2022_period_mean/cents_to_dollar,#Industrial retail electricity price provided by EIA for 2017-2022 period[24]
        'Ash_disposal_price':-0.0318,#Biosteam default price
        'Steam_regeneration_price': 0,#Steam, cooling water, and chilled water are regenerated by on-site facilities. 
        #The regeneration and heat transfer prices given are accounted for by the capital cost.
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
                      
                      #product prices
                      'crude_glycerol':product_prices_per_kg['Crude_glycerol'],
                      'recovered_C5_to_C8_MCA_fraction':product_prices_per_kg['C5_C8_fraction'],
                      'pelargonic_acid_rich_fraction':correlation_based_bulk_prices['Pelargonic_acid'],
                      'fatty_acid_blend':product_prices_per_kg['Fatty_acid_blend'],
                      'azelaic_acid_product_stream':correlation_based_bulk_prices['Azelaic_acid'],
                      'crude_methanol':product_prices_per_kg['Crude_methanol'],
                      
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
                        'HoySoy_oil':0.7050,#Based on GWP of soybean oil includes soybean meal displacement[26]   
                        'Citric_acid': 1.4821,#GHG-100 GREET 2022 [26]
                        'Hydrogen_peroxide':1.0813, #GHG-100 GREET 2022 [26]
                        'Tungstic_acid':68.5,#Value based on midpoint results available for tungsten carbide[27]
                        'Cobalt_acetate':8.8913, #Value based on cobalt nitrate [26]
                        'HCl':2.0486, #Value based on conc HCl production in the U.S.[26]
                        'Resin': 0.1324, #Value for GHG-100 high impact polystyrene [26]
                        'Sodium_hydroxide':0.3602, #GHG-100 value for sodium hydroxide, GREET 2022 [26]
                        'Sodium methoxide': 1.85, #Ecoinvent, production of sodium methoxide, RoW [28]
                        'Calcium chloride':0.017,#Based on calcium nitrate, GREET 2022 [26]
                        'Heptane':6.4337e-1,#Ecoinvent,heptane to generic market for solvent, organic [28]
                        'Methanol':2.5154, #Methanol from coal,GREET 2022[28],
                        #products
                        'Fatty_acid_blend': 0.00105,#Fatty acid blend can be used as a lubricant. 
                        #GWP based on solvent (lubricant pathway) available in GREET [26]
                        'C5_C8_fraction': 10, #Based on adipic acid, C5_C9 fraction can displace fossil based adipic acid[29]
                        'Crude_glycerol': 2,#Based on biobased pathway for production of glycerol[29]
                        'Pelargonic_acid': 11.239, #Based on glyphosate production, RoW, Ecoinvent [28]
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
                      'crude_vegetable_Hosun': GWP_factors_per_Kg['HoSun_oil'],
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
                      
                      #product factors
                      'crude_glycerol':GWP_factors_per_Kg['Crude_glycerol'],
                      'recovered_C5_to_C8_MCA_fraction':GWP_factors_per_Kg['C5_C8_fraction'],
                      'pelargonic_acid_rich_fraction':GWP_factors_per_Kg['Pelargonic_acid'],
                      'fatty_acid_blend':GWP_factors_per_Kg['Fatty_acid_blend'],
                      'crude_methanol':GWP_factors_per_Kg['Methanol'],
                      
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
Tungstic_acid_upper_bound =  convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,0.5,]),#[18] Sigma aldrich
                                                  lab_price_per_kg = np.array([682.00,838.00]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff = -0.65)
Tungstic_acid_lower_bound =  convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,0.5,]),#[18] Sigma aldrich
                                                  lab_price_per_kg = np.array([682.00,838.00]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff = -0.85)

Cobalt_acetate_upper_bound = convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,0.5]),
                                                  lab_price_per_kg = np.array([379.58,#[16] GFS chemicals
                                                                               735.42,#[17] BEantown chemicals
                                                                               ]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff = -0.65)
Cobalt_acetate_lower_bound = convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.5,0.5]),
                                                  lab_price_per_kg = np.array([379.58,#[16] GFS chemicals
                                                                               735.42,#[17] BEantown chemicals
                                                                               ]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff = -0.85)
#Product prices
Pelargonic_acid_upper_bound = convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.906]),#[18] Sigma aldrich
                                                  lab_price_per_kg = np.array([74.06 #[18] Sigma aldrich
                                                                               ]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff = -0.65)

Pelargonic_acid_lower_bound = convert_lab_price_to_bulk_price(lab_quantity_in_kg= np.array([0.906]),#[18] Sigma aldrich
                                                  lab_price_per_kg = np.array([74.06 #[18] Sigma aldrich
                                                                               ]),
                                                  bulk_quantity_in_kg = 30, bulk_coeff = -0.85)
#for uniform distributions a +-25% range was assumed 
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
                       'Heptane solvent price':chaospy.Uniform(prices_per_Kg['Heptane']*Lower_bound_factor_25_per,prices_per_Kg['Heptane']*Upper_bound_factor_25_per),
                       'Hydrogen peroxide price':chaospy.Uniform(prices_per_Kg['Hydrogen_peroxide']*Lower_bound_factor_25_per,prices_per_Kg['Hydrogen_peroxide']*Upper_bound_factor_25_per),
                       'Citric acid price':chaospy.Uniform(prices_per_Kg['Citric_acid']*Lower_bound_factor_25_per,prices_per_Kg['Citric_acid']*Upper_bound_factor_25_per),
                       'Conc HCl price': chaospy.Uniform(prices_per_Kg['HCl']*Lower_bound_factor_25_per,prices_per_Kg['HCl']*Upper_bound_factor_25_per),
                       'Calcium chloride price':chaospy.Uniform(prices_per_Kg['Calcium_chloride']*Lower_bound_factor_25_per,prices_per_Kg['Calcium_chloride']*Upper_bound_factor_25_per),
                       'Hydrolysis resin price':chaospy.Uniform(prices_per_Kg['Resin']*Lower_bound_factor_50_per,prices_per_Kg['Resin']*Upper_bound_factor_50_per),
                       'Crude glycerol price': chaospy.Triangle(100*ratio_2023from2022/1000,160*ratio_2023from2022/1000,220*ratio_2023from2022/1000),#Based on historical crude glycerol prices provided in [14]
                       'Fatty acid blend price':chaospy.Uniform(product_prices_per_kg['Fatty_acid_blend']*Lower_bound_factor_50_per,product_prices_per_kg['Fatty_acid_blend']*Upper_bound_factor_50_per),
                       'Pelargonic acid price':chaospy.Uniform(Pelargonic_acid_lower_bound,Pelargonic_acid_upper_bound),
                       'C5_C8 fraction price':chaospy.Uniform(product_prices_per_kg['C5_C8_fraction']*Lower_bound_factor_50_per,product_prices_per_kg['C5_C8_fraction']*Upper_bound_factor_50_per),
                       #Utility prices
                       'Natural gas price':chaospy.Uniform(utility_prices['Natural_gas_price']*Lower_bound_factor_25_per,utility_prices['Natural_gas_price']*Upper_bound_factor_25_per),
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
Oxcp_lower_bound = 0.7#current state of technology[31]
Oxcp_upper_bound = 0.9#potential improvement

#Oxidative cleavage reaction of intermediates to form azelaic acid and pelargonic acid
Ox_lower_bound = 0.9#current state of technology [31]
Ox_upper_bound = 0.99#potential improvement

#Products of the reaction participate in side reactions like
#esterification to produce long chain esters
#Baseline value not available in literature
#Ranges were assumed based on #TODO: find a citation
Side_lower_bound = 0.01
Side_upper_bound = 0.10

#Decarboxylation reaction of the Methyl_oxo_nonanoic acid and Nonanal
#Specific value not available in literature
#Ranges were assumed based on #TODO: find a citation
Dc_lower_bound = 0.1
Dc_upper_bound = 0.5

#Ranges for mimimum and maximum catalyst reuse were not available and were assumed
#TODO: find ranges 
TA_reuse_lower = 3
TA_reuse_upper = 9
CA_reuse_lower = 3
CA_reuse_upper = 9 

#Turbogen efficiency and boiler efficiencies were varied uniformly in the range +-50% efficiency range from baseline
Turbogen_eff_lower = Lower_bound_factor_50_per*0.85
Turbogen_eff_upper = Upper_bound_factor_50_per*0.85
Boiler_eff_lower = Lower_bound_factor_50_per*0.80
Boiler_eff_upper = Upper_bound_factor_50_per*0.80




coupled_para_dist = {'Tungstic acid moles': chaospy.Triangle(TA_moles_lower_bound,TA_moles_mid_point,TA_moles_upper_bound),
                     'Cobalt acetate moles': chaospy.Triangle(CA_moles_lower_bound,CA_moles_midpoint,CA_moles_upper_bound),
                     'Oxidative cleavage reaction conversion primary':chaospy.Uniform(Oxcp_lower_bound,Oxcp_upper_bound),
                     'Oxidative cleavage reaction conversion':chaospy.Uniform(Ox_lower_bound,Ox_upper_bound),
                     'Decarboxylation reaction conversion':chaospy.Uniform(Dc_lower_bound,Dc_upper_bound),
                     'Side reaction conversion':chaospy.Uniform(Side_lower_bound,Side_upper_bound),
                     'Dihydroxylation reaction time':chaospy.Uniform(Dih_time_lower ,Dih_time_upper),
                     'Oxidative reaction time':chaospy.Uniform(Ox_time_lower,Ox_time_upper),
                     'Turbogen efficiency':chaospy.Uniform(Turbogen_eff_lower,Turbogen_eff_upper),
                     'Boiler efficiency':chaospy.Uniform(Boiler_eff_lower,Boiler_eff_upper),
                     'Tungstic acid reusability':chaospy.Uniform(TA_reuse_lower,TA_reuse_upper),
                     'Cobalt acetate reusability':chaospy.Uniform(CA_reuse_lower,CA_reuse_upper)
                     }

environmental_facs_dist = {'Oil_GWP':{'HoSun_oil':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['HoSun_oil'],Upper_bound_factor_50_per*GWP_factors_per_Kg['HoSun_oil']),
                                      'HoySoy_oil':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['HoySoy_oil'],Upper_bound_factor_50_per*GWP_factors_per_Kg['HoySoy_oil'])},
                           'Tungstic_acid':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['Tungstic_acid'],Upper_bound_factor_50_per*GWP_factors_per_Kg['Tungstic_acid']),
                           'Cobalt_acetate':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['Cobalt_acetate'],Upper_bound_factor_50_per*GWP_factors_per_Kg['Cobalt_acetate']),
                           'HCl':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['HCl'],Upper_bound_factor_50_per*GWP_factors_per_Kg['HCl']),
                           'Resin':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['Resin'],Upper_bound_factor_50_per*GWP_factors_per_Kg['Resin']),
                           'Sodium methoxide':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['Sodium methoxide'],Upper_bound_factor_50_per*GWP_factors_per_Kg['Sodium methoxide']),
                           'Heptane':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['Heptane'],Upper_bound_factor_50_per*GWP_factors_per_Kg['Heptane']),
                           'Fatty_acid_blend':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['Fatty_acid_blend'],Upper_bound_factor_50_per*GWP_factors_per_Kg['Fatty_acid_blend']),
                           'C5_C8_fraction':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['C5_C8_fraction'],Upper_bound_factor_50_per*GWP_factors_per_Kg['C5_C8_fraction']),
                           'Pelargonic_acid':chaospy.Uniform(Lower_bound_factor_50_per*GWP_factors_per_Kg['Pelargonic_acid'],Upper_bound_factor_50_per*GWP_factors_per_Kg['Pelargonic_acid']),
                           'Cooling_tower_chemicals': chaospy.Uniform(Lower_bound_factor_50_per*Utility_GWP_factors['Cooling_tower_chemicals'],Upper_bound_factor_50_per*Utility_GWP_factors['Cooling_tower_chemicals']),
                           'Boiler_chems':chaospy.Uniform(Lower_bound_factor_50_per*Utility_GWP_factors['Boiler_chems'],Upper_bound_factor_50_per*Utility_GWP_factors['Boiler_chems']),
                           'Lime_boiler':chaospy.Uniform(Lower_bound_factor_50_per*Utility_GWP_factors['Lime_boiler'],Upper_bound_factor_50_per*Utility_GWP_factors['Lime_boiler']),
                           }
                            
#%%
#References:
#[1] https://fred.stlouisfed.org/release/tables?rid=46&eid=142872#snid=142874
#[2] https://catcost.chemcatbio.org/materials-library
#[3] https://www.chemicalbook.com/SupplierPriceList_EN.aspx?cbn=CB6292238&c=5kg#price
#[4] 10.1002/bbb.1640; Biofuels, Bioprod. Bioref. 10:299–315 (2016)
#[5] US High Oleic Soybeans & High Oleic Soybean Oil Sourcing Guide for International CustomersThird Edition – May 2023Prepared for US Soybean Export Council
#[6] https://www.ers.usda.gov/data-products/oil-crops-yearbook/oil-crops-yearbook/#All%20Tables
#[7] https://www.agriculture.com/crops/soybeans/high-oleic-soybeans-have-promise
#[8] https://corporate.perduefarms.com/news/press-releases/expansion-of-high-oleic-soybean-production/
#[9] https://www.hoosieragtoday.com/2019/09/24/consider-high-oleic-soybeans-in-2020-to-add-premium-to-your-price/
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


