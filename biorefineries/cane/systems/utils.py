# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import thermosteam as tmo

def set_sugarcane_composition(stream, water, fiber, sugar):
    chemicals = stream.chemicals
    if 'Sugar' not in chemicals:
        IDs = ('Sucrose', 'Glucose')
        chemicals.define_group('Sugar', IDs, composition=stream.imass[IDs], wt=True)
    if 'Fiber' not in chemicals:
        IDs = ('Cellulose', 'Hemicellulose', 'Lignin')
        chemicals.define_group('Fiber', IDs, composition=stream.imass[IDs], wt=True)
    if 'Other' not in chemicals:
        IDs = ('Ash', 'Solids')
        chemicals.define_group('Other', IDs, composition=stream.imass[IDs], wt=True)
    F_mass = stream.F_mass
    other = 1 - water - fiber - sugar
    assert other > 0, (water, fiber, sugar, other)
    stream.imass['Water', 'Fiber', 'Sugar', 'Other'] = F_mass * np.array([water, fiber, sugar, other])

def convert_fiber_to_lignocelluosic_components(stream, ignore_acetate=False):
    chemicals = stream.chemicals
    if chemicals is convert_fiber_to_lignocelluosic_components.last_chemicals:
        prxn = convert_fiber_to_lignocelluosic_components.last_reaction
    else:
        cellulose_rxn = tmo.Reaction('Cellulose -> Glucan', 'Cellulose', 1.0,
                                     basis='wt', chemicals=chemicals)
        cellulose_rxn.basis = 'mol'
        # Bagasse composition https://www.sciencedirect.com/science/article/pii/S0144861710005072
        # South american; by HPLC
        # Glucan: 41.3%
        # Xylan: 24.9%
        # Galactan: 0.6%
        # Arabinan: 1.7%
        # Lignin: 23.2%
        # Acetyl: 3.0%
        if ignore_acetate:
            hemicellulose_rxn = tmo.Reaction(
                '27.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan', 'Hemicellulose', 
                1.0, basis='wt', chemicals=chemicals
            )
        else:
            hemicellulose_rxn = tmo.Reaction(
                '30.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan + 3 Acetate', 'Hemicellulose', 
                1.0, basis='wt', chemicals=chemicals
            )
        hemicellulose_rxn.basis = 'mol'
        convert_fiber_to_lignocelluosic_components.last_chemicals = chemicals
        convert_fiber_to_lignocelluosic_components.last_reaction = prxn = tmo.ParallelReaction(
            [cellulose_rxn, hemicellulose_rxn]
        )
    prxn(stream)

convert_fiber_to_lignocelluosic_components.last_chemicals = None