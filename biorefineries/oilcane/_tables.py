# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 01:57:46 2021

@author: yrc2
"""
import os
import numpy as np
import pandas as pd
import biosteam as bst
import biorefineries.cornstover as cs
import biorefineries.oilcane as oc
from thermosteam.utils import array_roundsigfigs

__all__ = (
    'save_detailed_expenditure_tables',
    'save_detailed_life_cycle_tables',   
)

def save_detailed_expenditure_tables(sigfigs=3):
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
    IDs = (
        # 'S1',
        # 'S2', 
        'O1', 
        'O2'
    )
    names = [
        # 'Conventional Sugarcane Biorefinery',
        # 'Cellulosic Sugarcane Biorefinery',
        'Conventional Oilcane Biorefinery',
        'Cellulosic Oilcane Biorefinery',
    ]
    syss = [get_sys(i) for i in IDs]
    teas = [get_tea(i) for i in IDs]
    product_IDs = [oc.advanced_ethanol.ID, oc.cellulosic_ethanol.ID, oc.biodiesel.ID]
    tables = {
        'VOC': bst.report.voc_table(syss, product_IDs, names, with_products=True),
        'FOC': cs.foc_table(teas, names),
        'CAPEX': cs.capex_table(teas, names),
    }
    for key, table in tables.items(): 
        values = array_roundsigfigs(table.values, sigfigs=3, inplace=True)
        if key == 'CAPEX': # Bug in pandas
            for i, col in enumerate(table):
                table[col] = values[:, i]
        # if key == 'VOC':
        #     new_index = sorted(
        #         table.index, 
        #         key=lambda x: -abs(table.loc[x][1:].sum()) if x[0] == 'Raw materials' else 0,
        #     )
        #     tables[key] = table.reindex(new_index)
        table.to_excel(writer, key)
    writer.save()
    return tables
    
def save_detailed_life_cycle_tables(sigfigs=3):
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
    def get(configuration, name, RIN=False):
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
    ethanol_streams = [(get(i, 'advanced_ethanol'), get(i, 'cellulosic_ethanol')) for i in IDs]
    index = ['Energy allocation',
             'Economic allocation',
             'Displacement allocation']
    columns = [
        'Sugarcane DC [kg∙CO2e∙kg-1]',
        'Sugarcane ICF [kg∙CO2e∙kg-1]',
        'Oilcane DC [kg∙CO2e∙kg-1]',
        'Oilcane ICF [kg∙CO2e∙kg-1]',
    ]
    methods = ('GWP_ethanol_allocation', 'GWP_ethanol', 'GWP_ethanol_displacement')
    values = np.zeros([len(methods), len(IDs)])
    for i, method in enumerate(methods):
        for j, name in enumerate(IDs):
            oc.load(name)
            values[i, j] = getattr(oc, method)() 
    df_gwp = pd.DataFrame(values, index=index, columns=columns)
    tables = {
        'Inventory': bst.report.lca_inventory_table(
            syss, oc.GWP, ethanol_streams, names
        ),
        'Displacement allocation': bst.report.lca_displacement_allocation_table(
            syss, oc.GWP, ethanol_streams, 'ethanol', names
        ),
        'Energy allocation factors': bst.report.lca_property_allocation_factor_table(
            syss, property='energy', units='GGE/hr', system_names=names, groups=('ethanol',),
        ),
        'Economic allocation factors': bst.report.lca_property_allocation_factor_table(
            syss, property='revenue', units='USD/hr', system_names=names, groups=('ethanol',),
        ),
        'Displacement allocation factors': bst.report.lca_displacement_allocation_factor_table(
            syss, ethanol_streams, oc.GWP, names, groups=('ethanol',),
        ),
        'GWP ethanol': df_gwp,
    }
    for key, table in tables.items(): 
        array_roundsigfigs(table.values, sigfigs=3, inplace=True)
        table.to_excel(writer, key) 
    writer.save()
    return tables