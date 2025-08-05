# coding: utf-8
# %%
from chemprice import PriceCollector

smiles_list = ["CC(=O)NC1=CC=C(C=C1)O", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "O=C(C)Oc1ccccc1C(=O)O"]

pc = PriceCollector()
pc.setMolportApiKey('57fde93d-e1d0-4618-9edb-48a03ff1234a')
pc.setMolportUsername('ouwenp@illinois.edu')
pc.setMolportPassword('03Mar25%np')

pc.status()

pc.check()
# %%
all_prices = pc.collect(smiles_list)
# %%
Best_all_prices = pc.selectBest(all_prices)
Best_all_prices
