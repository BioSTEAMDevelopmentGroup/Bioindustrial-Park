# -*- coding: utf-8 -*-
"""
Created on 2025-08-05 19:22:05

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import pubchempy as pcp

c = pcp.Compound.from_cid(5090)

# print(c.molecular_formula)
# print(c.molecular_weight)
# print(c.smiles)
# print(c.xlogp)
# print(c.iupac_name)
# print(c.synonyms)

# %%
Glucose = pcp.get_compounds('Glucose', 'name')
print(Glucose)
print(Glucose[0].molecular_formula)
print(Glucose[0].molecular_weight)
print(Glucose[0].smiles)  # %%

# %%
from chemprice import PriceCollector
import pandas as pd
import biosteam as bst
pc = PriceCollector()
pc.setMolportApiKey('57fde93d-e1d0-4618-9edb-48a03ff1234a')
pc.setMolportUsername('ouwenp@illinois.edu')
pc.setMolportPassword('03Mar25%np')

pc.status()

pc.check()
# # %%
# Glucose_all_prices = pc.collect(Glucose[0].smiles)
# Best_Glucose_all_prices = pc.selectBest(Glucose_all_prices)

# # %%
# Best_Glucose_all_prices.to_excel('Glucose_best_prices.xlsx', index=False)
# Glucose_all_prices.to_excel('Glucose_all_prices.xlsx', index=False)
# glucosedf = pd.read_excel('Glucose_all_prices.xlsx')
# glucosedf.head()
# %%
# KOH = pcp.get_compounds('KOH', 'name')
# print(KOH[0].smiles)
KOH2 = bst.Chemical('KOH')
print(KOH2.smiles)


# %%
# KOH_all_prices = pc.collect(KOH[0].smiles)
KOH2_all_prices = pc.collect([KOH2.smiles])
# %%
KOH2_all_prices.to_excel('KOH_all_prices.xlsx', index=False)

KOHdf = pd.read_excel('KOH_all_prices.xlsx')
KOHdf.head()

