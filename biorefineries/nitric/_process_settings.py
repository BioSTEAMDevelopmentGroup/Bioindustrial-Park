# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 15:40:58 2025

@author: IGB
"""
import biosteam as bst

#%%
def load_preferences_and_process_settings(T,flow_units,N,P_units,CE,
                                          indicator,
                                          electricity_EI):
    bst.preferences.T = T
    bst.preferences.flow = flow_units
    bst.preferences.N = N
    bst.preferences.P = P_units
    bst.preferences.composition = True
    bst.preferences.light_mode()
    bst.preferences.save()
    bst.settings.CEPCI = CE  
    bst.settings.define_impact_indicator(key=indicator, units='kg*CO2e')
    # bst.settings.electricity_price = electricity_price
    bst.settings.set_electricity_CF(indicator,electricity_EI, basis='kWhr', units='kg*CO2e')
    
    # bst.settings.get_agent('chilled_brine').heat_transfer_price = 0