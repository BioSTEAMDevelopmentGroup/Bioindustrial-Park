"""
Global Constants for Biorefinery Optimization
"""
import numpy as np

# --- Conversion Factors ---
ethanol_conversion_factor = 90.4 # gal of ethanol per wet metric ton of miscanthus
gal_to_MMgal = 1000000 # gallons to million gallons
ton_to_kg = 1000 # metric ton to kilograms

# --- Refinery Cost Parameters ---
a_ethanol = 11.07

# --- Product Configurations ---

PRODUCT_CONFIG = {
    'EtOH': {
        'market_price': np.mean((2.68, 3)),         # USD/gal
        'a_coeff': 11.07,             
        'conv_factor': ethanol_conversion_factor / gal_to_MMgal,          
        'target_demand': 320,      # Total system target (MMgal)
        # 'bc_func': b_c_vectorized 
    },
    'AA': {
        'market_price': np.mean((1.40,1.65)), # USD/kg
        'a_coeff': 4.15, 
        'conv_factor': 0.000245528, # to transform from metric ton of feedstock to 10^6kg per yr of product
        'target_demand': 426.0,
        # 'bc_func': b_c_AA_vectorized
    },
    'SA': {
        'market_price': np.mean((3.23,3.29 )),# USD/kg
        'a_coeff': 3.32,
        'conv_factor': 0.000216572, # to transform from metric ton of feedstock to 10^6kg per yr of product
        'target_demand': 26.5,
        # 'bc_func': b_c_SA_vectorized
    },
    'KS': {
        'market_price': 7.12,# USD/kg
        'a_coeff': 15.7,
        'conv_factor': 6.57157E-05, # to transform from metric ton of feedstock to 10^6kg per yr of product
        'target_demand': 20.0,
        # 'bc_func': b_c_KS_vectorized
    },
    'LA': {
        'market_price': np.mean((1.13,5.10)),# USD/kg
        'a_coeff': 3.3,
        'conv_factor': 0.00032046, # to transform from metric ton of feedstock to 10^6kg per yr of product
        'target_demand': 250.0,
        # 'bc_func': b_c_LA_vectorized
    }
}

__all__ = [a_ethanol, ethanol_conversion_factor, gal_to_MMgal, ton_to_kg, PRODUCT_CONFIG]