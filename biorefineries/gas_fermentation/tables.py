# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 01:57:46 2021

@author: yrc2
"""
import os
import numpy as np
import pandas as pd
import biosteam as bst
from biorefineries.gas_fermentation  import Biorefinery
from biorefineries.tea.cellulosic_ethanol_tea import foc_table, capex_table
from thermosteam.utils import array_roundsigfigs

__all__ = (
    'save_system_reports',
    'save_detailed_expenditure_tables',
    'save_detailed_life_cycle_tables',  
    'save_distribution_table',
)
names = key_scenarios = ('acetate', 'acetate/glucose-seed', 'glucose')

def save_system_reports():
    folder = os.path.dirname(__file__)
    folder = os.path.join(folder, 'results')
    for scenario in key_scenarios:
        pm = Biorefinery(scenario=scenario)
        scenario = scenario.replace('/', '_')
        filename = f'{scenario}_detailed_report.xlsx'
        file = os.path.join(folder, filename)
        pm.system.save_report(file)

def save_detailed_expenditure_tables(sigfigs=3):
    folder = os.path.dirname(__file__)
    folder = os.path.join(folder, 'results')
    filename = 'expenditures.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    product = 'product'
    process_models = [
        Biorefinery(scenario=i)
        for i in key_scenarios
    ]
    for pm in process_models: pm.product.price = pm.tea.solve_price(pm.product)
    systems = [i.system for i in process_models]
    teas = [i.tea for i in process_models]
    tables = {
        'VOC': bst.report.voc_table(systems, product, names, with_products=True),
        'FOC': foc_table(teas, names),
        'CAPEX': capex_table(teas, names),
    }
    for key, table in tables.items(): 
        values = array_roundsigfigs(table.values, sigfigs=3, inplace=True)
        if key == 'CAPEX': # Bug in pandas
            for i, col in enumerate(table):
                table[col] = values[:, i]
        table.to_excel(writer, sheet_name=key)
    writer.close()
    return tables
    
def save_detailed_life_cycle_tables(sigfigs=3):
    process_models = [Biorefinery(scenario=i) for i in key_scenarios]
    systems = [i.system for i in process_models]
    folder = os.path.dirname(__file__)
    folder = os.path.join(folder, 'results')
    filename = 'life_cycle.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    streams = [i.product for i in process_models]
    tables = {
        'LCA breakdown': bst.report.lca_displacement_allocation_table(
            systems, 'GWP', streams, 'dodecanol', names
        ),
    }
    index = ['Dodecanol [kg∙CO2e∙kg-1]']
    values = np.zeros([len(index), len(names)])
    for j, pm in enumerate(process_models):
        GWPs = [pm.carbon_intensity]
        for i, GWP in enumerate(GWPs):
            values[i, j] = GWP()
    df_gwp = pd.DataFrame(values, index=index, columns=names)
    tables['Estimated environmental impact'] = df_gwp
    for key, table in tables.items(): 
        array_roundsigfigs(table.values, sigfigs=3, inplace=True)
        table.to_excel(writer, sheet_name=key) 
    writer.close()
    return tables

def save_distribution_table():
    folder = os.path.dirname(__file__)
    folder = os.path.join(folder, 'results')
    filename = 'parameters.xlsx'
    file = os.path.join(folder, filename)
    br = Biorefinery(scenario=key_scenarios[0])
    df_dct = br.model.get_distribution_summary()
    table = df_dct['Uniform']
    table['Baseline'] = [i.baseline for i in br.parameters]
    assert len(df_dct) == 1
    table.to_excel(file)
    return table

if __name__ == '__main__':
    save_system_reports()
    save_detailed_expenditure_tables()
    save_detailed_life_cycle_tables()
    save_distribution_table()