# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 16:17:35 2022 

@author: Lavanya
"""
from chaospy import distributions as shape
from biorefineries.oleochemicals.systems_baseline import F
import biosteam as bst
#Prices and GWP
#All GWP's are in [kg*CO2*eq / kg]

##Prices from the literature
##Ref: https://doi.org/10.1016/j.indcrop.2020.112411
#Raw material prices
#Given std deviation: 0.34, mean = 6.8
HOSO_price_USD_per_ton = shape.normal(mu = 6.8, sigma = 0.34)
##Given std deviation: 5.8, mean = 0.302
Methanol_price_USD_per_ton = shape.normal(mu = 5.8, sigma = 0.302)
##Given std deviation: 92, mean = 975
H2O2_USD_per_ton = shape.normal(mu = 975, sigma = 92)
##Given std deviation: 1000, mean = 33000
Tungstic_acid_USD_per_ton = shape.normal(mu = 3000, sigma = 1000)
##Given std deviation: 500, mean = 10918
Cobalt_acetate_USD_per_ton = shape.normal(mu = 10918, sigma = 500)
#Given triangular distribution a:100, b: 500, c:200
Crude_glycerol_USD_per_ton = shape.Triangle(100, 500 ,200)
#Calcium hydroxide price, found on alibaba - $130.00 - $300.00/ ton Ref: https://www.alibaba.com/product-detail/Slaked-Lime-Ca-OH-2-min_1600553478160.html
# Calcium_hydroxide_USD_per_ton = shape?
#Product_prices
##Given std deviation: 0.11, mean = 8.82
Azelaic_acid_USD_per_ton =  shape.normal(mu = 8.82, sigma = 0.11)
##Given std deviation: 0.2, mean = 8.2
Pelargonic_acid_USD_per_ton = shape.normal(mu = 8.2, sigma = 0.2)

# Tungstic acid : 6.85*10^4 kg CO2 eq /ton of of WC (tungsten carbide) # Ref: http://dx.doi.org/10.1016/j.jclepro.2017.02.184



# HCl
# Ref: Greet
# HCl: 1.816 CO2/Kg
