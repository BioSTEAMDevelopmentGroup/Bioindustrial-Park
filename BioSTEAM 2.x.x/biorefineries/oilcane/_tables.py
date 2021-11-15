# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 01:57:46 2021

@author: yrc2
"""
import os
import pandas as pd
import biorefineries.cornstover as cs
import biorefineries.oilcane as oc
from ._parse_configuration import (
    parse
)


def save_detailed_expenditure_tables(name):
    number, agile = parse(name)
    folder = os.path.dirname(__file__)
    folder = os.path.join(folder, 'results')
    filename = f'expenditures_{number}'
    if agile: filename += '_agile'
    filename += '.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    oc.load(name)
    cs.voc_table(oc.sys, oc.tea, [oc.ethanol, oc.biodiesel]).to_excel(writer, 'VOC', startrow=1)
    cs.foc_table(oc.tea).to_excel(writer, 'FOC', startrow=1)
    cs.capex_table(oc.tea).to_excel(writer, 'CAPEX', startrow=1)
    writer.save()