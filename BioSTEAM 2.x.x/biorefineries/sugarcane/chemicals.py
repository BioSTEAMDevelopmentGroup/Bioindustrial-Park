# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 06:42:02 2020

@author: yoelr
"""
from biorefineries.lipidcane.chemicals import lipidcane_chemicals

sugarcane_chemicals = lipidcane_chemicals.subgroup(
    ['Ash', 'Cellulose', 'Hemicellulose', 'Lignin',
     'Glucose', 'Sucrose', 'Solids', 'Water', 'Ethanol',
     'Octane', 'DryYeast', 'H3PO4', 'CaO', 'Flocculant', 'CO2']
)