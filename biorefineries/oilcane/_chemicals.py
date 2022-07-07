# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 06:42:02 2020

@author: yoelr
"""
import thermosteam as tmo

__all__ = ('create_chemicals',)

def create_chemicals():
    from biorefineries import lipidcane as lc
    from biorefineries import cornstover as cs
    removed = {'SuccinicAcid', 'H2SO4', 'Z_mobilis', 'Oil'}
    chemicals = tmo.Chemicals([
        i for i in (lc.chemicals.tuple + cs.chemicals.tuple) if i.ID not in removed
    ])
    chemicals.compile()
    chemicals.define_group('Lipid', ['PL', 'FFA', 'MAG', 'DAG', 'TAG'])
    chemicals.define_group('Oil', ['PL', 'FFA', 'MAG', 'DAG', 'TAG'])
    chemicals.set_synonym('DryYeast', 'Cellmass')
    return chemicals
    
