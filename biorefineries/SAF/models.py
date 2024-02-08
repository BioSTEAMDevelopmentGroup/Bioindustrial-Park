#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 13:57:09 2024

@author: wenjun
"""

import biosteam as bst
from biosteam.evaluation import Model, Metric
from chaospy import distributions as shape
from biorefineries.SAF.system import create_SAF_sys
from biorefineries.SAF import load_process_settings
from biorefineries.SAF._tea import SAF_TEA, create_energycane_SAF_tea

__all__ = {'ceate_model'}

funcs = create_funcs(SAF_tea=SAF_tea, flowsheet=flowsheet)

# %%

# =============================================================================
# Models for uncertainty and sensitivity analyses
# =============================================================================

def create_model(flowsheet=None):
    if not flowsheet:
        load_process_settings()
        SAF_sys = create_SAF_sys()
        flowsheet = SAF_sys.flowsheet
    else: SAF_sys = flowsheet.system.SAF_sys
    SAF_sys.simulate() # need this to initialize some settings
    s = flowsheet.stream
    u = flowsheet.unit
    
    # =============================================================================
    # Overall biorefinery metrics
    # =============================================================================

    # Minimum product selling price of jet fuel stream
    
    jet_fuel = s.jet_fuel

    SAF_tea = SAF_sys.TEA
    def get_MPSP():
        jet_fuel.price = 0
        for i in range(3):
            MPSP = jet_fuel.price = jet_fuel.solve_price(jet_fuel)
        return MPSP
    
    energycane = s.energycane
    
    # Yield in 10^6 kg/yr
    get_annual_factor = lambda: SAF_tea.operating_days*24
    get_total_yield = lambda: funcs['get_jet_fuel_flow']()/1e6
    
    # Yield in % of dry feedstock
    get_mass_yield = lambda: jet_fuel.F_mass/(energycane.F_mass-energycane.imass['Water'])
    
    get_overall_TCI = lambda:SAF_tea.TCI/1e6
    get_jet_fuel_sale = lambda:get_total_yield()*jet_fuel.price
    
    # Including negative product sales (ash disposal/diesel/gasoline) but excluding electricity credit
    ash = s.ash_disposal
    diesel = s.diesel
    gasoline = s.gasoline
    
    get_ash_sale = lambda: ash.F_mass*ash.price*get_annual_factor()/1e6
    get_diesel_sale = lambda: diesel.F_mass*diesel.price*get_annual_factor()/1e6
    get_gasoline_sale = lambda: gasoline.F_mass*gasoline.price*get_annual_factor()/1e6
    
    get_operating_cost = lambda: SAF_tea.AOC/1e6-get_ash_sale()-get_diesel_sale()-get_gasoline_sale()
    
    # Including negative product sales (ash/gypsum disposal) but excluding electricity credit
    get_material_cost = lambda: SAF_tea.material_cost/1e6-get_ash_sale()-get_diesel_sale()-get_gasoline_sale()
    
    metrics = [Metric('MPSP', get_MPSP, '$/gal'),
               Metric('Total product yield', get_total_yield, '10^6 kg/yr'),
               Metric('Product mass yield', get_mass_yield, '%'),
               Metric('Fermentation titer', get_titer, 'g/L'),
               Metric('Product purity', get_purity, '%'),
               Metric('Product recovery', get_recovery, '%'),
               Metric('Total capital investment', get_overall_TCI, '10^6 $'),
               Metric('Annual operating cost', get_operating_cost, '10^6 $/yr'),
               Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
               Metric('Annual product sale', get_lactic_sale, '10^6 $/yr'),
               ]

    
