# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:33:19 2025

@author: yoelr
"""
import biosteam as bst
bst.settings.set_thermo(['CO2', 'H2', 'AceticAcid', 'O2', 'Dodecanol', 'H2O', 'Glucose'])

price_H2 = 3
rx_AcOH = bst.Reaction(
    'CO2 + H2 + H2O -> AceticAcid', reactant='H2', basis='wt', 
    correct_atomic_balance=True
)
rx_Dodecanol = bst.Reaction(
    'AceticAcid -> Dodecanol + H2O + CO2', reactant='AceticAcid', basis='wt', 
    correct_atomic_balance=True
)
theoretical_yield_AcOH = rx_AcOH.istoichiometry['AceticAcid']
theoretical_yield_Dodecanol = rx_Dodecanol.istoichiometry['Dodecanol']

theoretical_yield_net = theoretical_yield_AcOH * theoretical_yield_Dodecanol
H2_per_Dodecanol_wt = 1 / theoretical_yield_net
Dodecanol_production_cost_H2 = price_H2 * H2_per_Dodecanol_wt

price_glucose = 0.215
rx_Dodecanol = bst.Reaction(
    'Glucose -> Dodecanol + H2O + CO2', reactant='Glucose', basis='wt', 
    correct_atomic_balance=True
)
theoretical_yield_net = rx_Dodecanol.istoichiometry['Dodecanol']
Glucose_per_Dodecanol_wt = 1 / theoretical_yield_net
Dodecanol_production_cost_Glucose = price_glucose * Glucose_per_Dodecanol_wt