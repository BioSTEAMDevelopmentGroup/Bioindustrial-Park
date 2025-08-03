#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat July 07 13:15:00 2025

Microalgae biorefinery to produce medium chain fatty acids 
by anaerobic fermentation without external electron donor addition

References
----------
[1] BioSTEAM Documentation: 
    https://biosteam.readthedocs.io/en/latest/API/thermosteam/Chemicals.html
[2] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
[3] 3-Hydroxypropionic acid biorefineries project:
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/HP

@author: Xingdong Shi
@version: 0.0.1
"""

import numpy as np
import pandas as pd
import thermosteam as tmo
from ._chemicals import chems

_kg_per_ton = 907.18474
# Baseline feedstock flow rate
def get_feedstock_flow(dry_composition, moisture_content, dry_flow):
    dry_array = chems.kwarray(dry_composition)
    wet_flow = dry_flow / (1-moisture_content)
    moisture_array = chems.kwarray(dict(Water=moisture_content))
    feedstock_flow = wet_flow * (dry_array*(1-moisture_content)+moisture_array)
    return feedstock_flow

# !!! This is the dry composition of corn stover; update to the composition you need
dry_composition = dict(
    Carbohydrate=0.234, Lipid=0.059, Protein=0.5369, Ash=0.0512, Extract=0.1189)

moisture_content = 0.05 #!!! This is the moisture content of preprocessed corn stover; update to the moisture content you need
dry_feedstock_flow = 2205. * _kg_per_ton / 24. # !!! Changing this will alter the total flow without altering the composition
 
baseline_feedflow = get_feedstock_flow(dry_composition, moisture_content, 
                                       dry_feedstock_flow)

# Price dictionary for all relevant chemicals and products
price = {
    'GlucoAmylase': 6.16, # $/kg (NREL 2018)
    'AlphaAmylase': 6.16, # $/kg (NREL 2018)
    'Yeast': 2,         # $/kg Alibaba
    'Ethanol': 1.78,     # $/kg https://catcost.chemcatbio.org/materials-library
    'Octanol': 1.00,       # $/lb https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/
    'AceticAcid': 1.31, # $/kg https://catcost.chemcatbio.org/materials-library
    'ButyricAcid': 1.72, # $/kg  https://www.imarcgroup.com/butyric-acid-pricing-report?
    'CaproicAcid': 2.89, # $/kg Increasing the economic value of lignocellulosic stillage through medium-chain fatty acid production
    'Butanol': 2.27,     # $/kg https://catcost.chemcatbio.org/materials-library
    #'Butanol': 0.7, # $/lb https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/
    'HeptanoicAcid': 3.17, # $/kg https://www.expertmarketresearch.com/prefeasibility-reports/heptanoic-acid-manufacturing-plant-project-report
    'CaprylicAcid': 5.07,  # $/kg Increasing the economic value of lignocellulosic stillage through medium-chain fatty acid production
    'SulfuricAcid': 0.11, # $/kg https://catcost.chemcatbio.org/materials-library
    'Lime': 0.12, # https://catcost.chemcatbio.org/materials-library
    'AmmoniumHydroxide': 0.204, # BDO, succinic, HP program
    'NaOH': 1.01, #https://catcost.chemcatbio.org/materials-library
}

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
CEPCI = {1997: 386.5,
         1998: 389.5,
         2007: 525.4,
         2009: 521.9,
         2010: 550.8,
         2011: 585.7,
         2012: 584.6,
         2013: 567.3,
         2014: 576.1,
         2016: 541.7,
         2017: 567.5,
         2018: 603.1,
         2019: 607.5,
         2020: 596.2,
         2021: 708.8,
         2022: 816.0,
         2023: 797.9}


 


