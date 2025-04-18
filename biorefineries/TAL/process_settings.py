#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module is a modified implementation of modules from the following:
[i]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[ii]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[iii]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040
"""


'''
References
----------
[1] Argonne National Laboratory. The Greenhouse gases, Regulated Emissions,
    and Energy use in Transportation (GREET) Model https://greet.es.anl.gov/
    (accessed Aug 25, 2020).
[2] ecoinvent 3.8 https://www.ecoinvent.org/ (accessed May 24, 2024).

'''

import biosteam as bst
import thermosteam as tmo
from biorefineries.TAL.chemicals_data import chems

bst.CE = 607.5 # year 2019

_kg_per_ton = 907.18474
_lb_per_kg = 2.20462
_liter_per_gallon = 3.78541
_ft3_per_m3 = 35.3147

#### Method for converting to 2019$ ####
# Where available, chemical indices are used (2010 - 2022 indices are available).
# To convert from Y1 to Y2: price_Y2 = price_Y1 * chem_index[Y2]/chem_index[Y1]
# For years prior to 2010, GDP indices are first used to convert from the given year to 2010,
# following which chemical indices are used.
# E.g., to convert from Y0 to Y2: price_Y2 = price_Y0 * GDP_index[Y1]/GDP_index[Y0] * chem_index[Y2]/chem_index[Y1]

_GDP_2007_to_2016 = 1.114 / 0.961
_GDP_2008_to_2016 = 1.114 / 0.990
_GDP_2008_to_2010 = 1.012 / 0.990
_GDP_2007_to_2010 = 1.012 / 0.961

_chemical_2011to2016 = 102.5 / 91.7
_chemical_2013to2016 = 102.5 / 101.3
_chemical_2014to2016 = 102.5 / 105.3

_chemical_2017to2016 = 102.5 / 106.9
_chemical_2020to2016 = 102.5 / 113.8 # average of Jan and Feb
_chemical_2022to2016 = 102.5 / 145.3

chem_index = { # Dictionary of chemical indices
                    2010: 82.2,
                    2011: 79.5,
                    2012: 83.7,
                    2013: 87.9,
                    2014: 91.3,
                    2015: 93.1,
                    2016: 88.8,
                    2017: 92.7,
                    2018: 93.3,
                    2019: 97.0, # the TEA year for the 2023 TAL production study
                    2020: 100.2,
                    2021: 112.0,
                    2022: 125.829,
                    }

_corn_bushel_to_kg = 25.402 # https://www.ers.usda.gov/webdocs/publications/41880/33132_ah697_002.pdf

# From USD/dry-ton to USD/kg in 2016$, 20% moisture content
# changed from Humbird et al., 2011 to Davis et al., 2018
corn_stover_price = 71.26 / _kg_per_ton * 0.8 

# 2.2 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal
# based on Annual Energy Outlook (AEO) from Energy Information Adiministration (EIA)
# (https://www.eia.gov/outlooks/aeo/), which is $0.7402/kg and similar to the
# 2.2/(2988/1e3) = $0.736/kg based on a density of 2988 g/gal from H2 Tools
# Lower and upper bounds are $1.37/gal and $2.79/gal, or $0.460/kg and $0.939/kg
# _ethanol_kg_2_gal =  2.9721617599077566
# price range in 2016$/kg: 
_ethanol_V = chems.Ethanol.V('l', 298.15, 101325) # molar volume in m3/mol
_ethanol_MW = chems.Ethanol.MW
_ethanol_kg_2_gal = _liter_per_gallon/_ethanol_V*_ethanol_MW/1e6
ethanol_price = 2.2 / _ethanol_kg_2_gal
	

# Dipotassium hydrogen phosphate (DPHP)
# https://www.alibaba.com/product-detail/Food-Grade-Dipotassium-Dydrogen-Phosphate-Trihydrate_60842047866.html?spm=a2700.7724857.normalList.2.4ef2457e3gPbfv&s=p
# DISREGARD: https://www.sigmaaldrich.com/catalog/product/mm/105104?lang=en&region=US
DPHP_price = 1.15

# 2.86 is the average motor gasoline price between 2010-2019 in 2016 $/gal	
# based on AEO from EIA, density of gasoline is 2.819 kg/gal	
# based on Lower and Higher Heating Values of Hydrogen and Other Fuels	
# from H2 Tools maintained by Pacific Northwest National Laboratory	
# (https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels)	
denaturant_price = 2.86 / 2.819

# 1.41e6 is $/yr and 4279 in kg/hr from Table 33 of Davis et al., 2018 (TAL scenario)
# 7880 is operating hours/yr on Page 10 of Davis et al., 2018,
# cost is negative because it's a product stream
ash_disposal_price = -1.41e6 / (4279*7880)

# Assums no cost/credit for baseline, the same as ash disposal for the lower bound,	
# for the upper bound (i.e., positive selling price indicating profit), use 	
# USGS 2015-2019 average  free on bard price in $/metric ton for crude gypsum. 
# National Minerals Information Center. Mineral Commodity Summaries 2020; 
# U.S. Geological Survey, 2020.
# Assuming all prices were in their nominal year (e.g., 2015 price in 2015$)
# and adjusted to 2016$
# 2015: 7.80 * 1.114 / 1.100 = 7.90
# 2016: 8.00
# 2017: 7.50 * 1.114 / 1.134 = 7.37
# 2018: 8.30 * 1.114 / 1.157 = 7.99
# 2019: 8.00 * 1.114 / 1.185 = 7.52
# (7.90+8.00+7.37+7.99+7.52) / 5 = 7.76 (in metric tonne)
# For the lower bound (i.e., negative selling price indicating cost), use price from
# Aden et al., 2002: $0.0094/lb in 2000$ = 0.0094*1.114/0.802*2.20462 = $0.0288/kg
# in 2016$
gypsum_price = 0

# Baseline from Davis et al., 2018, lower bound is 2015-2019 average of 	
# hydrate lime in $/ton at plant from Mineral Commodity Summaries 2020.	
# 2015: 146.40 * (1.114/1.100) / 907.18474 = 0.163	
# 2016: 145.50 / 907.18474 = 0.160	
# 2017: 147.10 * (1.114/1.134) / 907.18474 = 0.159	
# 2018: 151.50 * (1.114/1.157) / 907.18474 = 0.161	
# 2019: 151.00 * (1.114/1.185) / 907.18474 = 0.156	
# (0.163+0.160+0.159+0.161+0.156) / 5 = 0.160	
# Upper bound is +10% from baseline = 0.1189 * _lb_per_kg * 1.1 = 0.288
lime_price = 0.1189 * _lb_per_kg

# Mentioned in P53 of Humbird et al., not into any units, but a cashflow
# The original cost is $466,183 every 5 years, converted to per hour assuming 96% uptime
baghouse_bag_price = 466833 / 5 / (24*365*0.96)

# 4.70 is the average natural gas price in 2016$/Mcf based on AEO from EIA,
# which is $0.231/kg at 273.15 K or $0.253/kg at 298.15 K using BioSTEAM,   
# similar to the 4.7/1000/22*1000 = $0.214/kg at 273.15 K using 22 g/ft3 from H2 Tools
# Using the same conversion, lower and upper bounds should be 
# $3.68/Mcf and $5.65/Mcf, or $0.198/kg and $0.304/kg
CH4_V = chems.CH4.V(298.15, 101325) # molar volume in m3/mol
CH4_MW = chems.CH4.MW
natural_gas_price = 4.70/1e3*_ft3_per_m3*CH4_V * (1e3/CH4_MW) *\
    chem_index[2019]/chem_index[2016]

# https://www.alibaba.com/product-detail/Tricalcium-Phosphate-Tricalcium-Phosphate-TCP-Tricalcium_60744013678.html?spm=a2700.galleryofferlist.0.0.42f16684C9iJhz&s=p
TCP_price = 850 / _kg_per_ton # tricalcium (di)phosphate


# TAL_price = 1.88 # initial value
SA_price = 1.88 # initial value

# Currently not in use
# Methanol price from Goellner et al., production from natural gas (Case 3),
# average of two load structures in 2011$,
# crude methanol with ~1% CO2 and 1% H2O
methanol_price = (311.17+345.39)/2/_kg_per_ton * _chemical_2011to2016


# Acetoin product selling price
# https://www.alibaba.com/product-detail/Acetoin-CAS-NO-513-86-0_60673118759.html?spm=a2700.galleryofferlist.0.0.4d906f82dIhSkn
# acetoin_price = 5.925
acetoin_price = 3. # assumed
# Isobutyraldehyde product selling price
# https://www.alibaba.com/product-detail/China-manufacture-Isobutyraldehyde-price_60837108075.html?spm=a2700.galleryofferlist.0.0.753369fcbZcNhe
# IBA_price = 1.2
IBA_price = 0. # assumed

TAL_price = 4. # assumed; when solving for MPSP, this is merely the initial value and has no effect on results

#%% Feedstocks

# # Glucose / D-glucose / dextrose - from USDA
# # $/lb # USDA 2015-2019 mean
# # https://www.ers.usda.gov/data-products/sugar-and-sweeteners-yearbook-tables/sugar-and-sweeteners-yearbook-tables/#World,%20U.S.,%20and%20Mexican%20Sugar%20and%20Corn%20Sweetener%20Prices
# glucose_price = 0.3798 * _lb_per_kg
# # in $/kg:
# # 2015-2019 mean: 37.98	
# # 2015-2019 5th percentile: 36.00 
# # 2015-2019 95th percentile: 39.5

# Glucose / D-glucose / dextrose - from review paper Cheng et al. 2019 # https://doi.org/10.1002/bbb.1976
# $/kg in 2017$: 
# 0.22 # https://doi.org/10.1016/j.indcrop.2005.08.004  
# 0.26 # https://doi.org/10.22004/ag.econ.28658
# 0.33 # https://doi.org/10.1002/bbb.1475
# $/kg converted to 2019$:
# 0.23 # https://doi.org/10.1016/j.indcrop.2005.08.004  
# 0.27 # https://doi.org/10.22004/ag.econ.28658
# 0.34 # https://doi.org/10.1002/bbb.1475
# mean of 3 values in 2019$: 0.28
# mean of 3 values in 2016$:
glucose_price = ((0.22 + 0.26 + 0.33)/3) * chem_index[2016]/chem_index[2017]

# Corn stover
from biorefineries.lactic._process_settings import feedstock_price # in 2016$
cornstover_price = feedstock_price

# Corn
# $/bushel # USDA 2015-2019 mean
# https://www.nass.usda.gov/Charts_and_Maps/Agricultural_Prices/pricecn.php
corn_price = 3.543 / _corn_bushel_to_kg
# in $/kg:
# 2015-2019 mean: 0.139
# 2015-2019 5th percentile: 0.127
# 2015-2019 95th percentile: 0.150

# Sugarcane
from biorefineries.cane.streams import sugarcane # in 2018$
sugarcane_price = sugarcane['price']

#%%
# https://www.alibaba.com/product-detail/Manufacturer-of-Hexyl-Alcohol-Hexanol-n_60403061175.html?spm=a2700.galleryofferlist.0.0.1021992cx1VYY8
hexanol_price = 6.

# https://www.alibaba.com/product-detail/N-heptane-Heptane-Heptane-Heptane-Supply_62451341262.html?spm=a2700.galleryofferlist.0.0.2b407553zLGZvf&s=p
heptane_price = 1900./_kg_per_ton

# Viswanathan et al. 2020
toluene_price = 1.1

# search pages (all with 'trade assurance' and 'verified supplier' filters):
    # https://www.alibaba.com/trade/search?spm=a2700.galleryofferlist.leftFilter.d_filter.6fbb7c9clLc8dv&fsb=y&IndexArea=product_en&assessmentCompany=true&keywords=hydrogen+gas&ta=y&tab=all&
    # https://www.alibaba.com/trade/search?spm=a2700.galleryofferlist.leftFilter.d_filter.76512d31TN44FR&fsb=y&IndexArea=product_en&assessmentCompany=true&keywords=cas+1333-74-0&ta=y&tab=all&
    # https://www.alibaba.com/trade/search?spm=a2700.galleryofferlist.leftFilter.d_filter.542c64bepOk1ZA&fsb=y&IndexArea=product_en&assessmentCompany=true&keywords=1333-74-0&ta=y&tab=all&

# vendor listing 1: https://www.alibaba.com/product-detail/Industrial-Grade-H2-Hydrogen-Gas-Cylinder_1600624231996.html?spm=a2700.galleryofferlist.normal_offer.d_title.577364beUs9l40
# vendor listing 2 (same vendor): https://www.alibaba.com/product-detail/Promotion-Compressed-Hydrogen-Gas-With-Cylinder_1600948500210.html?spm=a2700.galleryofferlist.p_offer.d_title.253d7c9c5yZj8z&s=p
hydrogen_price = 1.

# https://www.energy.gov/eere/fuelcells/hydrogen-shot
hydrogen_renewable_price = 5.

# Ni-Al2O3-SiO2 catalyst
### --- 2024.06.03 search
    ## Alibaba listings
        # search pages: 
            # Raney Nickel (Nickel-Alumina alloy): https://www.alibaba.com/trade/search?spm=a2700.galleryofferlist.leftFilter.d_filter.33fa26503Y6HGm&fsb=y&IndexArea=product_en&assessmentCompany=true&keywords=raney+nickel&ta=y&tab=all&
            # Nickel: https://www.alibaba.com/trade/search?spm=a2700.galleryofferlist.leftFilter.d_filter.37af154aqb1xS5&fsb=y&IndexArea=product_en&assessmentCompany=true&keywords=7440-02-0&ta=y&tab=all&
        # vendor listings:
            # 1. https://www.alibaba.com/product-detail/Pure-99-99-99-Ni-Powdery_1600796014023.html?spm=a2700.galleryofferlist.normal_offer.d_title.b7d0154ab6qIRO
                # Nickel: $35/kg
            # 2. https://www.alibaba.com/product-detail/High-purity-nickel-powder-CAS-7440_1601043582452.html?spm=a2700.galleryofferlist.normal_offer.d_title.b7d0154ab6qIRO
                # Nickel: $40/kg
            # 3. https://www.alibaba.com/product-detail/Nickel-Alumina-modium-catalyst-cost-Price_62039092901.html?spm=a2700.galleryofferlist.normal_offer.d_title.31592650y6YlHm
                # Raney Nickel: $5/kg
            # 4. https://www.alibaba.com/product-detail/High-quality-customized-NiAl20-Raney-nickel_1601019811920.html?spm=a2700.galleryofferlist.normal_offer.d_title.31592650ztR5t4
                # $50/kg
            # 5. https://www.alibaba.com/product-detail/Customized-design-increase-specific-surface-area_1601019920484.html?spm=a2700.galleryofferlist.normal_offer.d_title.31592650ztR5t4
                # $35/kg
    ##
### ---
# Lowest = $5/kg
# Highest = $50/kg
# Mean = $33/kg
NiSiO2_price = RaneyNi_price = 33.

# https://www.alibaba.com/product-detail/Reagent-Grade-90-caustic-potash-potassium_62118969650.html?spm=a2700.galleryofferlist.0.0.28555ed4pKlEVC&s=p
KOH_price = 1.6

# https://www.alibaba.com/product-detail/Hydrochloric-acid-HCl-7647-01-0_60439085052.html?spm=a2700.galleryofferlist.0.0.4fa42c515nP2GH
HCl_price = 0.3

activated_carbon_price = 41. # $/ft^3 # Seader et al.


PdC_price = 0.075*(2045./0.0311035) \
      + (1-0.075)*0.45 # Pd : 2045 EIB (USD/troy-ounce) # https://apps.catalysts.basf.com/apps/eibprices/mp/ (accessed 6/3/2022)
                       # activated carbon: $0.45/kg # https://www-sciencedirect-com.proxy2.library.illinois.edu/science/article/pii/S2590174522000411?via%3Dihub
spent_PdC_price = 1. # assumed

acetone_price = 0.63 * _GDP_2008_to_2016 * _lb_per_kg # average of range ($0.44 - $0.82 /lb) from https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/

# Q1 2022 - Q4 2023 range from https://www.chemanalyst.com/Pricing-data/isopropyl-alcohol-31
# as reported:
# min: $1.225 /kg; max: 1.662/kg; mean: $1.387/kg
# converted to 2019$ (reported * chem_index[2019]/chem_index[2022]): 
# min: $0.944 /kg; max: 1.281/kg; mean: $1.069/kg
isopropanol_price = 1.387 * chem_index[2019]/chem_index[2022]

acetic_acid_price = 0.38 *_lb_per_kg * _GDP_2008_to_2010 * chem_index[2019]/chem_index[2010] # average of 2008$ range ($ 0.35 - 0.41 /lb) from # https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/
sodium_acetate_price = acetic_acid_price # unused

CSL_price = 0.0339 * _lb_per_kg * chem_index[2019]/chem_index[2016] # from lactic acid paper

# Amberlyst 70 catalyst
### --- 2024.06.03 search
    ## Chemicalbook listings
        # 1. $100/100kg Amberlyst 15 from bulk vendor ("VIP, 6-year, Enterprise Certified") listing https://www.chemicalbook.com/ProductDetail_EN_451808.htm
        # 2. $1/kg Amberlyst 15 from bulk vendor ("VIP, 5-year, Enterprise Certified") listing https://www.chemicalbook.com/ProductDetail_EN_916657.htm
    ##
    ## Alibaba listings # search page: https://www.alibaba.com/trade/search?spm=a2700.galleryofferlist.leftFilter.d_filter.1b624164QCx1xc&fsb=y&IndexArea=product_en&assessmentCompany=true&keywords=amberlyst&ta=y&tab=all&
        # 3. Amberlyst-15 https://www.alibaba.com/product-detail/Ion-exchange-resin-for-MTBE-equal_62529933735.html?spm=a2700.galleryofferlist.normal_offer.d_title.22a24164j8EA7y
            # $1.5 /L
            # 20 kg / 25 L
            # => $1.875 / kg
        # 4. Amberlyst DT https://www.alibaba.com/product-detail/High-Temperature-Resistance-Catalyst-Resin-equal_62530999330.html?spm=a2700.galleryofferlist.normal_offer.d_title.22a24164j8EA7y
            # $2000/metric ton        
            # => $2/kg
        # 5. Amberlyst-15: https://www.alibaba.com/product-detail/Chinese-manufacturers-Ion-exchange-resin-for_1600848354148.html?spm=a2700.galleryofferlist.normal_offer.d_title.22a24164j8EA7y
            # $1/kg
        # 6. Amberlyst(R) 15: https://www.alibaba.com/product-detail/purity-99-AMBERLYST-R-15-with_1600342344644.html?spm=a2700.galleryofferlist.normal_offer.d_title.22a24164j8EA7y
            # $1/kg
        # 7. Amberlyst(R) 15: https://www.alibaba.com/product-detail/AMBERLYST-R-15-cas-9037-24_1600993084228.html?spm=a2700.galleryofferlist.normal_offer.d_title.22a24164j8EA7y
            # $2100/metric ton
            # => $2.1/kg
        # 8. Amberlyst A45: https://www.alibaba.com/product-detail/High-Temperature-Resistance-Catalyst-Resin-equal_62530910597.html?spm=a2700.galleryofferlist.normal_offer.d_title.22a24164j8EA7y
            # $1/kg
        # 9. Amberlyst(R) 15: https://www.alibaba.com/product-detail/99-AMBERLYST-R-15-wet-type_1600267306406.html?spm=a2700.galleryofferlist.normal_offer.d_title.22a24164j8EA7y
            # $2/kg
    ##
### ---
# Lowest = $1/kg
# Highest = $2.1/kg
# Mean = $1.442/kg
amberlyst70_price = 1.442

# 2.0/kg from https://www.alibaba.com/product-detail/Desiccant-for-paints-and-varnishes-Acetylacetone_10000004008185.html?spm=a2700.galleryofferlist.normal_offer.d_price.e82458eazJeqeC
acetylacetone_price = 2.0 # 2,4-pentanedione or acetylacetone


# 225 - 230 $/tonne in 2007 # https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/
# 895.32 $/US ton in 2011 #  Humbird 2011 https://doi.org/10.2172/1013269
DAP_price = 0.5*\
            (0.2275 * _GDP_2007_to_2010 * chem_index[2019]/chem_index[2010]
             +
             0.89532 * chem_index[2019]/chem_index[2011])

# 390.09-789.25 USD / short ton in 2008; https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/            
NaOH_price = ((390.09+789.25)/2.)/_kg_per_ton * _GDP_2008_to_2010 * chem_index[2019]/chem_index[2010]

# 1.55 - 1.7/lb # https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/
# => $3.42 - 3.75 /kg in 2007$, mean of $3.58/kg in 2007$
# => $4.25 - 4.66 /kg in 2019$, mean of $4.45/kg in 2019$
THF_price = 3.58 * _GDP_2007_to_2010 * chem_index[2019]/chem_index[2010]

price = {'SA': SA_price,
         'PD': acetylacetone_price, # 2,4-pentanedione or acetylacetone
         'TCP': TCP_price,
         'IBA': IBA_price,
         'TAL': TAL_price,
         'KOH': KOH_price,
         'HCl': HCl_price,
         'Hydrogen': hydrogen_price,
         'Renewable hydrogen': hydrogen_renewable_price,
         'Acetoin': acetoin_price,
         'RaneyNi': RaneyNi_price,
         'Ni-SiO2': NiSiO2_price,
         'Amberlyst-70': amberlyst70_price,
         'PdC': PdC_price,
         'Spent PdC': spent_PdC_price,
         'Hexanol': hexanol_price,
         'Heptane': heptane_price,
         'Toluene': toluene_price,
         'Isopropanol': isopropanol_price,
         'Acetone': acetone_price,
         'Sulfuric acid': 0.0430 * _lb_per_kg,	
         # 0.1900 is for NH3	
         'AmmoniumHydroxide': 0.1900 * _lb_per_kg * 17.031/35.046,	
         'CSL': CSL_price,
         'Caustics': NaOH_price * 0.5, # 50 wt% NaOH/water mixture	
         'Sodium hydroxide': NaOH_price,
         'Boiler chems': 2.9772 * _lb_per_kg,	
         'Lime': lime_price,
         'Cooling tower chems': 1.7842 * _lb_per_kg,	
         'Makeup water': 0.0002 * _lb_per_kg,	
         # Cost of ash is negative because it's a product stream	
         'Ash disposal': ash_disposal_price,
         'Electricity': 0.070, # AEO from EIA, 2010-2019 average (0.067-0.074 range)	
         # $6.16/kg protein in 2016$, P25 of Davis et al., 2018	
         'Enzyme': 6.16,
         'DPHP': DPHP_price,
         'Baghouse bag': baghouse_bag_price,	
         'Natural gas': natural_gas_price,
         'Methanol': methanol_price,
         'Ethanol': ethanol_price,
         # Below currently not in use
         'Gypsum': gypsum_price,
         'Denaturant': denaturant_price,
         'DAP': DAP_price,
         'Activated carbon': activated_carbon_price,
         'Sodium acetate': sodium_acetate_price,
         'Acetic acid': acetic_acid_price,
         'Tetrahydrofuran': THF_price,
         'Glucose': glucose_price,
         'Corn stover': cornstover_price,
         'Corn': corn_price,
         'Sugarcane': sugarcane_price,
         }
    
#!!! Round all prices to 4 *decimal places*
for k in price.keys():
    price[k] = round(price[k], 4)

bst.PowerUtility.price = price['Electricity']

_lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
_mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
_hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
_mps.T = 233 + 273.15
_hps.T = 266 + 273.15

_cooling = bst.HeatUtility.get_cooling_agent('cooling_water')
_chilled = bst.HeatUtility.get_cooling_agent('chilled_water')
_chilled_brine = bst.HeatUtility.get_cooling_agent('chilled_brine')

_cooling.T = 28 + 273.15
_cooling.T_limit = _cooling.T + 9
# Side steam in CHP not a heat utility, thus will cause problem in TEA utility
# cost calculation if price not set to 0 here, costs for regeneration of heating
# and cooling utilities will be considered as CAPEX and OPEX of CHP and CT, respectively
for i in (_lps, _mps, _hps, _cooling, _chilled, _chilled_brine):
    i.heat_transfer_price = i.regeneration_price = 0
    # if i == _cooling: continue
    # i.heat_transfer_efficiency = 0.85


# %%

# =============================================================================
# Characterization factors (CFs) for life cycle analysis (LCA), all from ref [1] if not noted otherwise
# =============================================================================

CFs = {}

# =============================================================================
# 100-year global warming potential (GWP) in kg CO2-eq/kg
# =============================================================================
GWP_CFs = {
    'CH4': 0.40, # NA NG from shale and conventional recovery
    'CSL': 1.55,
    'DAP': 1.6354, # ecoinvent 3.8 [2] diammonium phosphate production, RoW
    
    # 'Enzyme': 2.24, 
    'Ethanol': 1.44,
    'Acetone': 2.5435, #  ecoinvent 3.8 [2] market for acetone, liquid, RoW
    
    'H2SO4': 44.47/1e3,   
    'Lime': 1.29 * 56.0774/74.093, # CaO to Ca(OH)2
    'CaO': 1.29,
    'NaOH': 2.11,
    # 'NH4OH': 2.64 * 0.4860, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,   
    # 'MEA': 3.4062, # ecoinvent 3.8 [2] ethanolamine production, RoW [monoethanolamine]
    'H3PO4': 1.3598, # ecoinvent 3.8 [2] purification of wet-process phosphoric acid to industrial grade, product in 85% solution state, RoW # cradle-to-gate
    'CO2': 0.87104, # ecoinvent 3.8 [2] carbon dioxide production, liquid, RoW
    'H2': 2.3716, # ecoinvent 3.8 [2] market for hydrogen, liquid, RoW
    
    'AceticAcid': 1.6198, # market for acetic acid, without water, in 98% solution state, GLO
    'SodiumAcetate': 1.6198*60.05196/82.033789, # adjusted acetic acid CF 
    
    'PD': (2*(58.080)*3.5917 + 102.089*2.5435)/(2*100.117), # Acetylacetone; based on GLO/RoW IPCC 2013 CFs of markets for precursors acetone and acetyl anhydride
    # 'DiammoniumSulfate': 1.2901, # ecoinvent 3.8 [2] market for ammonium sulfate, RoW
    # 'MagnesiumSulfate': 1.0411, # ecoinvent 3.8 [2] market for magnesium sulfate, GLO
    
    'NiSiO2':18.711, # ecoinvent 3.8 [2] market for nickel, class 1
    'Amberlyst70_':1.3803, # ecoinvent 3.8 [2] market for naphthalene sulfonic acid, GLO # Amberlyst-15 is 1,2-diethenylbenzene; 2-ethenylbenzene-1-sulfonic acid
    'Isopropanol':2.3219, # ecoinvent 3.8 [2] market for isopropanol, RoW
    'THF': 6.0475, # ecoinvent 3.8 [2] market for tetrahydrofuran, GLO
    'KOH': 2.63, # ecoinvent 3.8 [2] market for potassium hydroxide, GLO
    }



GWP_CF_array = chems.kwarray(GWP_CFs)
# In kg CO2-eq/kg of material
GWP_CF_stream = tmo.Stream('GWP_CF_stream', GWP_CF_array, units='kg/hr')
# CFs['GWP_CF_stream'] = GWP_CF_stream

GWP_CFs['Electricity'] = 0.4490 # kg CO2-eq/kWh GREET 2022 US Mix  # assume production==consumption, both in kg CO2-eq/kWh


GWP_CFs['Corn stover'] = 0.10945 # see Table S4 of the SI of Bhagwat et al. 2021
GWP_CFs['Sugarcane'] = 0.12158 # ecoinvent 3.6, sugarcane production, RoW, IPCC 2013 GWP-100a
# GWP_CFs['Sugarcane'] = 0.044535 # GREET 2023, Sugarcane Production for Brazil Plant
GWP_CFs['Corn'] = 0.2610 # GREET 2023, Corn Production for Biofuel Refinery
# GWP_CFs['Glucose'] = 1.2127 # ecoinvent 3.8 glucose production, GLO
GWP_CFs['Glucose'] = 0.7539 * 0.909 # GREET 2023, Glucose (from corn; based on Fuel-Cycle Fossil Energy Use and Greenhouse Gas Emissions of Fuel Ethanol Produced from U.S. Midwest Corn)
                                    # multiplied by 0.909 as feedstock dextrose monohydrate stream is 90.9 wt% glucose
CFs['GWP_100'] = GWP_CFs

# =============================================================================
# Fossil energy consumption (FEC), in MJ/kg of material
# =============================================================================


FEC_CFs = {
    'CH4': 50, # NA NG from shale and conventional recovery
    'CSL': 12,
    'DAP': 22.028, # ecoinvent 3.8 [2] diammonium phosphate production, RoW
    
    'Ethanol': 16,
    'Acetone': 66.852, #  ecoinvent 3.8 [2] market for acetone, liquid, RoW
    # 'Enzyme': 26,
    
    'H2SO4': 568.98/1e3,
    'Lime': 4.896 * 56.0774/74.093, # CaO to Ca(OH)2
    'CaO': 4.896, 
    'NaOH': 29,
    # 'NH4OH': 42 * 0.4860, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,
    # 'MEA': 67.898, # ecoinvent 3.8 [2] ethanolamine production, RoW [monoethanolamine]
    'H3PO4': 16.538, # ecoinvent 3.8 [2] purification of wet-process phosphoric acid to industrial grade, product in 85% solution state, RoW # cradle-to-gate
    'CO2': 7.4243, # ecoinvent 3.8 [2] carbon dioxide production, liquid, RoW
    'H2': 75.747, # ecoinvent 3.8 [2] market for hydrogen, liquid, RoW

    'AceticAcid': 45.611, # market for acetic acid, without water, in 98% solution state, GLO
    'SodiumAcetate': 45.611*60.05196/82.033789, # adjusted acetic acid CF 
    
    'PD': (2*(58.080)*66.852 + 102.089*70.817)/(2*100.117), # Acetylacetone; based on GLO/RoW cumulative energy demand CFs of markets for precursors acetone and acetyl anhydride
    
    # 'DiammoniumSulfate': 15.166, # ecoinvent 3.8 [2] market for ammonium sulfate, RoW
    # 'MagnesiumSulfate': 13.805, # ecoinvent 3.8 [2] market for magnesium sulfate, GLO
    
    'NiSiO2':229.44	, # ecoinvent 3.8 [2] market for nickel, class 1
    'Amberlyst70_':31.273, # ecoinvent 3.8 [2] market for naphthalene sulfonic acid, GLO # Amberlyst-15 is 1,2-diethenylbenzene; 2-ethenylbenzene-1-sulfonic acid
    'Isopropanol':63.878, # ecoinvent 3.8 [2] market for isopropanol, RoW
    'THF': 101.77, # ecoinvent 3.8 [2] market for tetrahydrofuran, GLO
    'KOH': 33.33, # ecoinvent 3.8 [2] market for potassium hydroxide, GLO
    }

FEC_CF_array = chems.kwarray(FEC_CFs)
# In MJ/kg of material
FEC_CF_stream = tmo.Stream('FEC_CF_stream', FEC_CF_array, units='kg/hr')

# CFs['FEC_CF_stream'] = FEC_CF_stream

FEC_CFs['Electricity'] = 5.724 # MJ/kWh # GREET 2022 US Mix #assume production==consumption, both in MJ/kWh

FEC_CFs['Corn stover'] = 1.68000 # see Table S4 in the SI of Bhagwat et al. 2021
FEC_CFs['Sugarcane'] = 0.37338 # ecoinvent 3.6, sugarcane production, RoW, IPCC 2013 GWP-100a
# FEC_CFs['Sugarcane'] = 0.28832 # GREET 2023, Sugarcane Production for Brazil Plant
FEC_CFs['Corn'] = 1.684 # GREET 2023, Corn Production for Biofuel Refinery
# FEC_CFs['Glucose'] = 14.507 # ecoinvent 3.8 glucose production, GLO
FEC_CFs['Glucose'] = 7.74 * 0.909 # GREET 2023, Glucose (from corn; based on Fuel-Cycle Fossil Energy Use and Greenhouse Gas Emissions of Fuel Ethanol Produced from U.S. Midwest Corn)
                                  # multiplied by 0.909 as feedstock dextrose monohydrate stream is 90.9 wt% glucose

CFs['FEC'] = FEC_CFs
