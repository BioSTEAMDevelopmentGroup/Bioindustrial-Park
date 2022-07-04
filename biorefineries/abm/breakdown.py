# -*- coding: utf-8 -*-
"""
Example module for calculating ethanol yield

"""
import os
import pandas as pd 
import biosteam as bst
from biorefineries import (
    cornstover as cs, 
    corn as cn, 
    sugarcane as sc
)

cn.load()
bst.CE = 607.5

folder = os.path.dirname(__file__)
file = f'corn_tables.xlsx'
file = os.path.join(folder, file)
writer = pd.ExcelWriter(file)
cn.corn_sys.simulate()
bst.report.voc_table(cn.corn_sys, ['ethanol'], ['Corn Ethanol Biorefinery']).to_excel(writer, 'VOC')
sc.foc_table(cn.corn_tea).to_excel(writer, 'FOC')
sc.capex_table(cn.corn_tea).to_excel(writer, 'CAPEX')
writer.save()

cs.load()
bst.CE = 607.5

file = f'cornstover_tables.xlsx'
file = os.path.join(folder, file)
writer = pd.ExcelWriter(file)
cs.cornstover_sys.simulate()
bst.report.voc_table(cs.cornstover_sys, ['ethanol'], ['Corn Stover Ethanol Biorefinery']).to_excel(writer, 'VOC')
cs.foc_table(cs.cornstover_tea).to_excel(writer, 'FOC')
cs.capex_table(cs.cornstover_tea).to_excel(writer, 'CAPEX')
writer.save()