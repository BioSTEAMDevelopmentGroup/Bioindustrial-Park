# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 21:14:00 2020

@author: yrc2
"""
from biorefineries import LAOs

__all__ = (
    'get_product_purity',
    'get_product_flow',
    'get_alcohol_loss',
    'get_LAOs_MPSP',
    'set_LAOs_MPSP',
)

fatty_alcohols = ('Hexanol','Octanol', 'Decanol', 'Dodecanol', 'Tetradecanol')
products = ('Hexene','Octene', 'Decene', 'Dodecene', 'Tetradecene')

def get_product_purity(product):
    """Get product purity in olefin wt. %."""
    return product.imass[products].sum() / product.F_mass
        
def get_product_flow():
    """Get product flow in ton / yr."""
    return sum([i.get_total_flow('ton/day') for i in LAOs.products]) * LAOs.LAOs_tea.operating_days
    
def get_alcohol_loss():
    """Get alcohol loss in ton / yr."""
    F_ton_per_yr = LAOs.wastewater.get_flow('ton/day', fatty_alcohols)
    return F_ton_per_yr * LAOs.LAOs_tea.operating_days

def get_LAOs_MPSP():
    """Get product price in USD/ton."""
    products = LAOs.products
    sales = sum([i.cost for i in products])
    F_mass_LAOs = sum([i.F_mass for i in products])
    LAOs_tea = LAOs.LAOs_tea
    return (LAOs_tea.solve_sales() / LAOs_tea.operating_days / 24  + sales) / F_mass_LAOs * 907.185 # To USD / ton

def set_LAOs_MPSP(MPSP):
    """Set product price in USD/ton."""
    price = MPSP / 907.185
    for i in LAOs.products:
        i.price = price