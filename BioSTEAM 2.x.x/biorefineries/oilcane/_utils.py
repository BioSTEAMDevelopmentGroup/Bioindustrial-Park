# -*- coding: utf-8 -*-
"""
Created on Thu May  6 07:29:54 2021

@author: yrc2
"""
from biorefineries import lipidcane2g as lc
from thermosteam.units_of_measure import convert

def get_results_summary(name):
    lc.load(name)
    MFPP = lc.lipidcane_tea.solve_price(lc.lipidcane)
    feedstock = lc.lipidcane.get_total_flow('ton/hr')
    GGE = (lc.ethanol.F_mass * 2.98668849 / 1.5
           + lc.biodiesel.get_total_flow('gal/hr') / 0.9536
           - lc.natural_gas.get_total_flow('ft3/hr') / 126.67)
    return {
        'MFPP [USD/ton]': convert(MFPP, 'USD/kg', 'USD/ton'),
        'Productivity [GGE/ton]': GGE / feedstock,
        'TCI [10^6 USD]': lc.lipidcane_tea.TCI,
    }
            
    