# -*- coding: utf-8 -*-
"""
"""

__all__ = ()

# https://www.sgh2energy.com/economics
green_hydrogen_price_range = (10, 15) # USD / kg
blue_hydrogen_price_range = (5, 7) # USD / kg
brown_hydrogen_price_range = (2, 3) # USD / kg

baseline_dodecanol_production = dict(
    titer = 1.9, # g / L    
    yield_ = 0.03, # g-product / g-acetate
    productivity = 0.05 / 24 # g / L / h
)

