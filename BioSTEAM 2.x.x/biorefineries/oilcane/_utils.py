# -*- coding: utf-8 -*-
"""
Created on Thu May  6 07:29:54 2021

@author: yrc2
"""
from biorefineries import oilcane as oc
from thermosteam.units_of_measure import convert

def get_results_summary(name):
    oc.load(name)
    MFPP = oc.tea.solve_price(oc.oilcane)
    feedstock = oc.oilcane.get_total_flow('ton/hr')
    GGE = (oc.ethanol.F_mass * 2.98668849 / 1.5
           + oc.biodiesel.get_total_flow('gal/hr') / 0.9536
           - oc.natural_gas.get_total_flow('ft3/hr') / 126.67)
    return {
        'MFPP [USD/ton]': convert(MFPP, 'USD/kg', 'USD/ton'),
        'Productivity [GGE/ton]': GGE / feedstock,
        'TCI [10^6 USD]': oc.oilcane_tea.TCI,
    }
            
    