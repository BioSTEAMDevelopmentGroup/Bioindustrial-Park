# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 16:17:35 2022

@author: Lavanya
"""
import math
import numpy as np
from scipy.stats import lognorm
import random
import biosteam as bst

##Prices from the literature
##Ref: https://doi.org/10.1016/j.indcrop.2020.112411

#Raw material prices
#Given std deviation: 0.34, mean = 6.8
HOSO_price_USD_per_ton = lognorm.rvs(s=0.34, scale=math.exp(6.8), size=2000)
##Given std deviation: 5.8, mean = 0.302
Methanol_price_USD_per_ton = lognorm.rvs(s=5.8, scale=math.exp(0.302), size=2000)
##Given std deviation: 92, mean = 975
H2O2_USD_per_ton = lognorm.rvs(s=5.8, scale=math.exp(0.302), size=2000)
##Given std deviation: 1000, mean = 33000
Tungstic_acid_USD_per_ton = lognorm.rvs(s=1000, scale=math.exp(33000), size=2000)
##Given std deviation: 500, mean = 10918
Cobalt_acetate_USD_per_ton = lognorm.rvs(s=500, scale=math.exp(10918), size=2000)
#Given triangular distribution a:100, b: 500, c:200
Crude_glycerol_USD_per_ton = random.triangular(100, 500, 200)


#Product_prices
##Given std deviation: 0.11, mean = 8.82
Azelaic_acid_USD_per_ton = lognorm.rvs(s=0.11, scale=math.exp(8.8), size=2000)
##Given std deviation: 0.2, mean = 8.2
Pelargonic_acid_USD_per_ton = lognorm.rvs(s=0.2, scale=math.exp(8.2), size=2000)

#Characterisation factors
GWP = 'GWP 100yr'
bst.settings.define_impact_indicator(key=GWP, units='kg*CO2e')

#TODO: how to set Cf for HOSO
#Global warming (incl. iLUC and biogenic CO2 uptake) in kg CO2-eq
#HOSO: 760 Kg CO2-eq/ton
#http://dx.doi.org/10.1016/j.jclepro.2014.10.011 

#H2O2: 899.2177 kg of CO2/ton
#Ref: Greet

# Cobalt Nitrate instead of acetate because it can be used as a catalyst as well
# TODO: Try finding cobalt acetate as well
# Ref for cobalt acetate: Greet
# Cobalt_acetate: 7.2691 Kg of CO2/ Kg of cobalt nitrate

# Calcium nitrate
# TODO: Find values for CaOH2
# Calcium_hydroxide: 555.42 Kg of CO2/ ton of calcium nitrate
# Ref: Greet

# HCl
# Ref: Greet
# HCl: 1.816 CO2/Kg
