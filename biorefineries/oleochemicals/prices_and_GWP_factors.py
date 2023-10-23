#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 16:43:50 2023
@author: lavanyakudli
"""
import numpy as np
__all__ = ('prices_per_Kg',
           'product_prices_per_kg',
           'utility_prices',
           'prices_per_stream',
           'GWP_per_stream'
           )

#%% Data for prices of materials and utility

#Freds Producers price index was used to update the prices of chemicals to 2023 prices[1]
#Index was chosen for the general category of all chemicals and allied products
PPI_2016 = 270.400 #Dec 2016
PPI_2021 = 349.906 #Dec 2021
PPI_2022 = 357.404 #Dec 2022
PPI_2023 = 347.123 #Sept 2023
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
        #feedstock prices, 
        'Crude_HoSun_oil': 1.21, #2023 price, Price based on high oleic sunflower seed price [4] converted to
        # high oleic sunflower oil based on price differential provided in [5],
        
        #2022 price to 2023
        'Crude_HoySoy_oil': 0.82*ratio_2023from2022, #Typically, farmers recieve 2$/buchel higher for high oleic soybeans[6].
        #Price evaluated by adding $2/buchel to historic average of soybean price and
        #then multiplying historic average of price differential between soybean and soybean oil to get price of oil[5].     
        #catalysts
        'Cobalt_acetate': 0.67,#2023 price,Using correlations provided in [7]. Correlation uses laboratory scale values provided in [8]
        # Bulk price = laboratory scale price*(bulk quantity/lab scale quantity)^(-0.75) 
        # Based on [8], 50g of cobalt acetate tetrahydrate costs 1140$
        'Tungstic_acid': 0.835 ,#2023 price,Using correlations provided in [7]. Correlation uses laboratory scale values provided in [9]
        # Bulk price = laboratory scale price*(bulk quantity/lab scale quantity)^(-0.75) 
        # Based on [9], 250 g of tungstic acid costs 500 $
        
        #2016 price to 2023
        'Sodium_methoxide': 2.93*ratio_2023from2016,# [10]
        }


product_prices_per_kg = {
        #products
        'Pelargonic_acid': 2.42*ratio_2023from1997,#1997 price, [8]
        'C5_C8_fraction': 2.00*ratio_2023from2018,#2018 price for caproic acid (major component),1500–2500 USD/ton [11]
        'Fatty_acid_blend':1.37*ratio_2023from2021,#2021 price for Stearic acid[2]
        'Crude_glycerol': 2.74*ratio_2023from2021,#2021 price for vegtable glycerine [2]
        'Azelaic_acid': 6.67*ratio_2023from2021,#2021 price based on the price of sebacic acid [2]
        'Sodium_oleate':1.04,#2023 price, using correlations provided in [7]. Correlation uses laboratory scale values provided in [12]
        # Bulk price = laboratory scale price*(bulk quantity/lab scale quantity)^(-0.75) 
        # Based on [12], 500 g of sodium oleate costs 313 $
        'Crude_methanol':0.792*ratio_2023from2021,#[2]
        }

#unit conversions
lb_to_kg = 0.453
cooling_tower_chem_density = 8.65 #lbs/gal [16]
boiler_chem_density = 9.10 #lb/ga[17]
cents_to_dollar = 100
tcf_to_cf = 1000 #thousand cubic feet to cubic feet
cf_to_m3 = 0.029 #cubic feet to cubic meter
kgal_to_L = 3.785412 *1000

natural_gas_prices_2017_to_2022 = [4.08,4.19,3.90,3.32,5.44,7.66] #Provided by EIA, Dollars per Thousand Cubic Feet
average_electricity_prices_2017_to_2022 = [6.88,6.92,6.81,6.67,7.18,8.45]#Provided by EIA, cents/Kwh
average_water_rates_2008_2012_2016 = [2.44,3.02,3.38]#$/kGal Provided by USDOE

#TODO: add cooling tower chems
#TODO: ask about heat transfer regen price
utility_prices = {   
        'System_makeup_water':np.mean(average_water_rates_2008_2012_2016)*ratio_utility_2023from2016/kgal_to_L ,#[13],[14]
        'Natural_gas_price': np.mean(natural_gas_prices_2017_to_2022)*ratio_utility_2023from2017_2022_period_mean/(tcf_to_cf*cf_to_m3),#Industrial natural gas price provided by EIA for 2017-2022 period[15]
        'Cooling_tower_chemicals':519/(55*cooling_tower_chem_density*lb_to_kg),#2023 price,Based on [16], 55 gal drum costs 519$
        'Boiler_chems':1.216/(55*boiler_chem_density*lb_to_kg), #2023 price,Based on [17],55 gal drum costs 609$
        'Lime_boiler':0.11*ratio_2023from2021,#2021 price based on average of all lime values provided in [2]
        'Electricity':np.mean(average_electricity_prices_2017_to_2022)*ratio_utility_2023from2017_2022_period_mean/cents_to_dollar,#Industrial retail electricity price provided by EIA for 2017-2022 period[18]
        'Ash_disposal_price':-0.0318,#Biosteam default price
        'Steam_regeneration_price': 0.03,
        'Heat_transfer_price': 0
        }


prices_per_stream = { #input prices
                      'crude_vegetable_Hosun': prices_per_Kg['Crude_HoSun_oil'],
                      'crude_vegetable_Hoysoy':prices_per_Kg['Crude_HoySoy_oil'],
                      'base_for_saponification_of_FFA':prices_per_Kg['Sodium_hydroxide'],
                      'citricacid_for_degumming':prices_per_Kg['Citric_acid'],
                      'fresh_HP':prices_per_Kg['Hydrogen_peroxide'],
                      'fresh_tungsten_catalyst':prices_per_Kg['Tungstic_acid'],
                      'fresh_cobalt_catalyst_stream':prices_per_Kg['Cobalt_acetate'],
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
                      'pelargonic_acid_rich_fraction':product_prices_per_kg['Pelargonic_acid'],
                      'fatty_acid_blend':product_prices_per_kg['Fatty_acid_blend'],
                      'azelaic_acid_product_stream':product_prices_per_kg['Azelaic_acid'],
                      'crude_methanol':product_prices_per_kg['Crude_methanol'],
                      'sodium_oleate_product':product_prices_per_kg['Sodium_oleate'],
                      
                      #utility prices
                      'natural_gas': utility_prices['Natural_gas_price'],
                      'lime_boiler':utility_prices['Lime_boiler'],
                      'boiler_chems':utility_prices['Boiler_chems'],
                      'system_makeup_water':utility_prices['System_makeup_water'],
                      'ash_disposal':utility_prices['Ash_disposal_price'],
                      'cooling_tower_chemicals':utility_prices['Cooling_tower_chemicals']
                      }

#%% Global warming potential factors (kg CO2-eq/Kg)

lbsperMWh_to_KgsperKWh = 4.53592e-4                 
GWP_factors_per_Kg = {  #inputs
                        'HoSun_oil':0.76, #Based on GWP of sunflower oil (incl. iLUC and biogenic CO2 uptake) [19]                                               
                        'HoySoy_oil':0.7050,#Based on GWP of soybean oil includes soybean meal displacement[20]   
                        'Citric_acid': 1.4821,#GHG-100 GREET 2022 [20]
                        'Hydrogen_peroxide':1.0813, #GHG-100 GREET 2022 [20]
                        'Tungstic_acid':68.5,#Value based on midpoint results available for tungsten carbide[21]
                        'Cobalt_acetate':8.8913, #Value based on cobalt nitrate [20]
                        'HCl':2.0486, #Value based on conc HCl production in the U.S.[20]
                        'Resin': 0.1324, #Value for GHG-100 high impact polystyrene [20]
                        'Sodium_hydroxide':0.3602, #GHG-100 value for sodium hydroxide, GREET 2022 [20]
                        'Sodium methoxide': 1.85, #Ecoinvent, production of sodium methoxide, RoW [22]
                        'Calcium chloride':0.017,#Based on calcium nitrate, GREET 2022 [20]
                        'Heptane':6.4337e-1,#Ecoinvent,heptane to generic market for solvent, organic [22]
                        'Methanol':2.5154, #Methanol from coal,GREET 2022[22],
                        #products
                        'Fatty_acid_blend': 0.00105,#Fatty acid blend can be used as a lubricant. 
                        #GWP based on solvent (lubricant pathway) available in GREET [20]
                        'C5_C8_fraction': 10, #Based on adipic acid, C5_C9 fraction can displace fossil based adipic acid[23]
                        'Crude_glycerol': 2,#Based on biobased pathway for production of glycerol[23]
                        'Pelargonic_acid': 11.239, #Based on glyphosate production, RoW, Ecoinvent [22]
                        'Sodium oleate': 2.4787} #Based on soap production, RoW, Ecoinvent [22]}
Utility_GWP_factors =   {                      
                        #utilities
                        'Electricity': 852.3*lbsperMWh_to_KgsperKWh, #Based on NAtional US average for electricty CO2 emissions [24]
                        'System_makeup_water':0.00035559,#Tap water production through conventional treatment [22]
                        'Natural_gas': 0.3983, #Based on GREET 2022, NA NG from Shale and Conventional Recovery [20]
                        'Cooling_tower_chemicals':0.2574, #Based on sodium bicarbonate which is a common corrosion inhibitor [20].
                        # Data for other cooling tower chemicals such as pH controllers and anti-scalants not found.
                        'Boiler_chems':1.5568, #Based on the production of sulfite which is commonly used as an oxygen scavenger in the boiler [20]
                        # Data for other boiler chemicals (neutralising amines) not found.
                        'Lime_boiler':1.2844, #Based on lime production from lime stone [20]
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
                      'sodium_oleate_product':GWP_factors_per_Kg['Sodium oleate'],
                      
                      #utility factors
                      'natural_gas': Utility_GWP_factors['Natural_gas'],
                      'lime_boiler':Utility_GWP_factors['Lime_boiler'],
                      'boiler_chems':Utility_GWP_factors['Boiler_chems'],
                      'system_makeup_water':Utility_GWP_factors['System_makeup_water'],
                      'cooling_tower_chemicals':Utility_GWP_factors['Cooling_tower_chemicals']
                      }





#%%
#References:
#[1] https://fred.stlouisfed.org/release/tables?rid=46&eid=142872#snid=142874
#[2] https://catcost.chemcatbio.org/materials-library
#[3] https://www.chemicalbook.com/SupplierPriceList_EN.aspx?cbn=CB6292238&c=5kg#price
#[4] https://www.sunflowernsa.com/growers/Marketing/daily-oilseed-sunflower-price/
#[5] https://www.ers.usda.gov/data-products/oil-crops-yearbook/oil-crops-yearbook/#All%20Tables
#[6] https://www.unitedsoybean.org/issue-briefs/high-oleic-soybeans/
#[7] https://www.researchgate.net/publication/271702798
#[8] https://www.thermofisher.com/order/catalog/product/044345.18
#[9] https://www.fishersci.com/shop/products/tungstic-acid-99-thermo-scientific/AC207615000
#[10]10.1002/bbb.1640; Biofuels, Bioprod. Bioref. 10:299–315 (2016)
#[11]https://doi.org/10.3390/en11061551
#[12]https://www.fishersci.com/shop/products/sodium-oleate-tci-america-3/O0057500G
#[13]https://www.energy.gov/sites/prod/files/2017/10/f38/water_wastewater_escalation_rate_study.pdf
#[14]https://fred.stlouisfed.org/series/PCU221221
#[15]https://www.eia.gov/dnav/ng/ng_pri_sum_dcu_nus_a.htm
#[16]https://www.coolingtowerchemicals.com/ProductDetails.asp?ProductCode=CTC1334NM
#[17]https://www.chemworld.com/Boiler-Chemical-p/chemworld-1394.htm?msclkid=59bcee3dabfc15d1fc428b13737412bb&utm_term=4585513243667217&utm_medium=cpc&utm_content=ChemWorld%203-29-17&utm_campaign=1234%20Products&utm_source=bing
#[18]https://www.eia.gov/electricity/data/browser/#/topic/7?agg=0,1&geo=g&endsec=vg&linechart=ELEC.PRICE.US-ALL.A~ELEC.PRICE.US-RES.A~ELEC.PRICE.US-COM.A~ELEC.PRICE.US-IND.A&columnchart=ELEC.PRICE.US-ALL.A~ELEC.PRICE.US-RES.A~ELEC.PRICE.US-COM.A~ELEC.PRICE.US-IND.A&map=ELEC.PRICE.US-ALL.A&freq=A&start=2011&end=2022&ctype=linechart&ltype=pin&rtype=s&maptype=0&rse=0&pin=
#[19]http://dx.doi.org/10.1016/j.jclepro.2014.10.011                           
#[20]GREET 2022
#[21]Ref: http://dx.doi.org/10.1016/j.jclepro.2017.02.184
#[22]Ecoinvent
#[23]https://doi.org/10.1021/acssuschemeng.2c05764
#[24]https://www.epa.gov/egrid/power-profiler#/









