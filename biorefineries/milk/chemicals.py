# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import numpy as np
from biorefineries import cellulosic

__all__ = (
    'create_chemicals',
)

# PROTEIN: DO NOT DELETE
# from thermosteam.chemicals.elements import array_to_atoms
# AAs = ['Threonine', 'Methionine',
#        'Phenylalanine', 'Histidine',
#        'Lysine', 'Valine', 'Isoleucine',
#        'Leucine', 'Serine', 'Glycine',
#        'Glutamic acid', 'Proline',
#        'Cysteine', 'Alanine', 'Tyrosine',
#        'Arginine']
# whey_composition_wt = np.array([5.4, 1.8, 2.5, 1.4, 7.1, 3.5, 3.8, 8.6,
#                                 4, 1.5, 15.5, 4.8, 0.8, 4.2, 2.4, 1.7])
# chemicals = bst.Chemicals(AAs)
# bst.settings.set_thermo(chemicals)
# formula_array = chemicals.formula_array
# MWs = np.array([i.MW for i in chemicals])
# whey_composition_by_mol = whey_composition_wt / MWs
# whey_composition_by_mol = whey_composition_by_mol / whey_composition_by_mol.sum()
# whey_composition_by_mol = whey_composition_by_mol.reshape([1, 16])
# formula_array = formula_array * whey_composition_by_mol
# protein_formula = array_to_atoms(formula_array.sum(axis=1))

def create_chemicals(yeast_includes_nitrogen=False):
    Yeast = bst.Chemical(
        'Yeast', 
        phase='l',
        formula='CH1.61O0.56N0.16' if yeast_includes_nitrogen else 'CH1.61O0.56',
        rho=1540,
        Cp=1.54,
        default=True,
        search_db=False,
        aliases=['Cellmass'],
    )
    
    # Sample work
    Water = bst.Chemical('Water')
    
    # Composition of Acid Whey Protein
    Solids = bst.Chemical(
        'Solids',
        phase = 'l',
        default=True,
        search_db=False,
    )
    Lactose = bst.Chemical(
        'Lactose',
        phase = 'l'
    )
    Galactose = bst.Chemical(
        'Galactose',
        phase = 'l'
    )
    LacticAcid = bst.Chemical(
        'LacticAcid',
        phase = 'l'
    )
    CitricAcid = bst.Chemical(
        'CitricAcid',
        phase = 'l'
    )
    GlutaricAcid = bst.Chemical(
        'GlutaricAcid',
        phase = 'l'
    )
    protein_formula = {
        'H': 10.019566174405634,
        'C': 4.974848530454082,
        'N': 1.1787790508313507,
        'O': 2.573004354214718,
        'S': 0.03480161876387264
    }
    formula = {i: round(j, 2) for i, j in protein_formula.items()}
    Protein = bst.Chemical(
        'Protein',
        search_db = False,
        default = True,
        atoms = formula,
        phase='l'
    )
    mineral_names = [
        'Ca', 'Na', 'P', 'K', 'Mg'
    ]
    mineral_weights = np.array([
        121, 37.9, 66.8, 164, 10.6
    ])
    minerals = [bst.Chemical(i, search_db=False, default=True, phase='l', formula=i) for i in mineral_names]
    Ash = bst.Chemical('Ash', default=True, phase='s', MW=1, search_db=False)
    gases = [
        bst.Chemical(i, phase='g') for i in 
        ('N2', 'CO', 'H2', 'O2', 'CO2', 'CH4')
    ]
    chemicals = bst.Chemicals([
        Water, Lactose, Galactose, LacticAcid, Protein, CitricAcid, GlutaricAcid,
        Yeast, Ash,
        Solids,
        'Dodecanol',
        'DodecylAcetate',        
        'DodecanoicAcid',
        'Tetradecanol',
        'Ca(OH)2',
        'Hexane',
        *gases,
        *minerals,
        
    ])
    IDs = set([i.ID for i in chemicals])
    CASs = set([i.CAS for i in chemicals])
    for i in cellulosic.create_cellulosic_ethanol_chemicals():
        if i.ID in IDs or i.CAS in CASs: continue
        chemicals.append(i)
    bst.settings.set_thermo(chemicals)
    chemicals.define_group(
        'Minerals', mineral_names, mineral_weights, wt=True
    )
    # Make acid whey a feed/stream
    acid_whey = dict(
        Lactose=3.33e1, Galactose=0.59e1, LacticAcid=0.65e1,
        CitricAcid=0.18e1, GlutaricAcid=0.06e1, Protein = 3.71, Minerals = 400.3e-2
    )
    acid_whey['Water'] = (
        1000 - acid_whey['Lactose']- acid_whey['Galactose']
        - acid_whey['LacticAcid'] - acid_whey['CitricAcid']
        - acid_whey['Protein'] - acid_whey['Minerals']
    )
    chemicals.define_group(
        'AcidWhey',
        [*acid_whey.keys()],
        [*acid_whey.values()],
        wt=True,
    )
    percent_solids = 0.30
    solids = sum([0.2, 4.4, 0.4])
    water_to_solids = (1 - percent_solids) / percent_solids
    chemicals.define_group(
        'UFPermeate',
        ['Protein', 'Lactose', 'Ash', 'LacticAcid', 'Water'],
        [0.2, 4.4, 0.4, 0.6, water_to_solids * solids - 0.6],
        wt=True,
    )
    chemicals.define_group(
        'SaccharifiedCornSlurry',
        ['Protein', 'Glucose', 'Ash', 'Water', 'Cellulose'],
        [2.77, 20.7, 0.425, 71.2, 3.65],
        wt=True,
    )
    chemicals.define_group(
        'GlucoseMedia',
        ['Protein', 'Glucose', 'Ash', 'Water'],
        [2.77, 20.7, 0.425, 71.2],
        wt=True,
    )
    chemicals.set_alias('Yeast', 'Cellmass')
    return chemicals

