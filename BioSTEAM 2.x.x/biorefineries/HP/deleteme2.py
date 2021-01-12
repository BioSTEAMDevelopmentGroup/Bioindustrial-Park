# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 16:51:14 2021

@author: yrc2
"""
import thermosteam as tmo
from biorefineries.HP.chemicals_data import HP_chemicals

ethyl3HP = '623-72-3'
methyl3HP = '6149-41-3'
solute = tmo.Chemical(methyl3HP)
solvents = tmo.Chemicals(['Butanol', 'Hexanol', 'Octanol', 'Decanol', 'Undecanol', 'Dodecanol'])
chemicals = tmo.Chemicals([solute, 'Water', *solvents])
Tbs = {str(i): i.Tb for i in chemicals} 
print(Tbs)


chemicals = tmo.Chemicals({*chemicals, *HP_chemicals})
chemicals.compile(skip_checks=True)
tmo.settings.set_thermo(chemicals)
chemicals.set_synonym('HP', '3-Hydroxypropionic acid')

s = tmo.Stream()
partition_coefficients = {}
for i in solvents:
    if Tbs[i.ID] < solute.Tb: continue
    s.empty()
    s.phase = 'l'
    key = ('Water', i.ID)
    s.set_flow([1, 1], 'L/hr', key)
    s.imass[solute.ID] = 1e-3 
    s.lle(T=300, top_chemical=i.ID)
    IDs = tuple([i.ID for i in s.available_chemicals])
    Ks = tmo.separations.partition_coefficients(IDs, s['l'], s['L'])
    partition_coefficients[i.ID] = Ks[IDs.index(solute.ID)]

print(partition_coefficients)