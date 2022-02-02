# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 01:57:46 2021

@author: yrc2
"""
import os
import pandas as pd
import biosteam as bst
import biorefineries.cornstover as cs
import biorefineries.oilcane as oc

__all__ = (
    'save_detailed_expenditure_tables',
    'save_detailed_life_cycle_tables',   
)

def save_detailed_expenditure_tables():
    folder = os.path.dirname(__file__)
    folder = os.path.join(folder, 'results')
    filename = 'expenditures.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    
    def get_sys(name):
        oc.load(name)
        return oc.sys
    
    def get_tea(name):
        oc.load(name)
        return oc.tea
    IDs = ('S1', 'S2', 'O1', 'O2')
    names = ['Conventional Sugarcane Biorefinery',
             'Cellulosic Sugarcane Biorefinery',
             'Conventional Oilcane Biorefinery',
             'Cellulosic Oilcane Biorefinery',]
    syss = [get_sys(i) for i in IDs]
    teas = [get_tea(i) for i in IDs]
    product_IDs = [oc.ethanol.ID, oc.biodiesel.ID]
    bst.report.voc_table(syss, product_IDs, names).to_excel(writer, 'VOC')
    cs.foc_table(teas, names).to_excel(writer, 'FOC')
    cs.capex_table(teas, names).to_excel(writer, 'CAPEX')
    writer.save()
    
def save_detailed_life_cycle_tables():
    try:
        # Energy allocation by gasoline gallon equivalent (GGE)
        bst.PowerUtility.define_property(
            name='energy', units='kW',
            fget=lambda power_utility: -power_utility.rate,
        )
        bst.Stream.define_property(
            name='energy', units='kJ/hr',
            # Ignore streams that are not ethanol (e.g. process water)
            fget=lambda stream: stream.LHV if stream.price else 0.,
        )
        
        # Economic/revenue allocation
        bst.PowerUtility.define_property(
            name='revenue', units='USD/hr',
            fget=lambda power_utility: -power_utility.cost,
        )
        bst.Stream.define_property(
            name='revenue', units='USD/hr',
            # Ignore streams that are not ethanol (e.g. ash disposal)
            fget=lambda stream: stream.cost if stream.price > 0. else 0.,
        )
    except:
        pass
    def get(configuration, name):
        oc.load(configuration)
        return getattr(oc, name)
    bst.settings.define_impact_indicator(oc.GWP, 'kg*CO2e')
    folder = os.path.dirname(__file__)
    folder = os.path.join(folder, 'results')
    filename = 'life_cycle.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    IDs = ('S1', 'S2', 'O1', 'O2')
    names = ['Conventional Sugarcane Biorefinery',
             'Cellulosic Sugarcane Biorefinery',
             'Conventional Oilcane Biorefinery',
             'Cellulosic Oilcane Biorefinery',]
    syss = [get(i, 'sys') for i in IDs]
    ethanol_streams = [get(i, 'ethanol') for i in IDs]
    bst.report.lca_inventory_table(
        syss, oc.GWP, ethanol_streams, names
    ).to_excel(writer, 'Inventory')
    bst.report.lca_displacement_allocation_table(
        syss, oc.GWP, ethanol_streams, 'ethanol', names
    ).to_excel(writer, 'Displacement allocation')
    bst.report.lca_property_allocation_factor_table(
        syss, property='energy', units='GGE/hr', system_names=names,
    ).to_excel(writer, 'Energy allocation factors')
    bst.report.lca_property_allocation_factor_table(
        syss, property='revenue', units='USD/hr', system_names=names,
    ).to_excel(writer, 'Economic allocation factors')
    bst.report.lca_displacement_allocation_factor_table(
        syss, ethanol_streams, oc.GWP, names
    ).to_excel(writer, 'Displacement allocation factors')
    writer.save()