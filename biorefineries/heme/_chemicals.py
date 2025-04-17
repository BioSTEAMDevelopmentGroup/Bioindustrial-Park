# -*- coding: utf-8 -*-
"""
Created on 2025-04-16 15:19:55

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import biosteam as bst

bst.Chemical('water')

def create_chemicals(set_thermo=True):
    chems = Chemicals([])
    def add_chemical(ID, ref=None, **data):
        chemical = Chemical(ID, **data) if ref is None else ref.copy(ID, **data)
        chems.append(chemical)
        return chemical










def create_leghemoglobin_chemicals():
    # Sample work

    # Culture Media
    ####
    Water = bst.Chemical('Water')
    Glucose = bst.Chemical('Glucose', phase='l', default=True)
    # FeSO4 //// 7 H2O
    FeSO4 = bst.Chemical('FeSO4', phase='s', default=True)
    KH2PO4 = bst.Chemical('KH2PO4', phase='s', default=True)
    NH4_2HPO4 = bst.Chemical('(NH4)2HPO4', phase='s', default=True)
    KOH = bst.Chemical('KOH', phase='s', default=True)
    # MgSO4 //// 7 H2O
    MgSO4 = bst.Chemical('MgSO4', phase='s', default=True)
    CitricAcid = bst.Chemical('CitricAcid', phase='s', default=True)
    ####

    # feed
    ####
    # Glucose
    # MgSO4
    IPTG = bst.Chemical('IPTG', phase='s', default=True)
    NH4_2SO4 = bst.Chemical('(NH4)2SO4', phase='s', default=True)
    # FeSO4
    # ZnSO4
    ####

    # Main Product
    ####
    # heme molecule formula: C34H32FeN4O4
    # heme_b = bst.Chemical('heme')
    heme_b = bst.Chemical('PubChem=26945')

    # Only protein: 144 amino acids Formula:C729H1166N200O219S2​
    # https://www.uniprot.org/uniprotkb/P02236/entry#sequences
    protein_formula = {
        'H': 1166,
        'C': 729,
        'N': 200,
        'O': 219,
        'S': 2
    }
    formula = {i: round(j, 2) for i, j in protein_formula.items()}
    Protein = bst.Chemical(
        'Protein',
        search_db = False,
        default = True,
        atoms = formula,
        phase='l'
    )



    Ash = bst.Chemical('Ash', default=True, phase='s', MW=1, search_db=False)
    gases =[
        bst.Chemical(i, phase='g') for i in 
        ('N2', 'CO', 'H2', 'O2', 'CO2', 'CH4')
    ]
    return chemicals

def create_trace_metal_solution():
    Water = bst.Chemical('Water')
    # trace metal solution
    ####
    HCl = bst.Chemical('HCl', phase='l', default=True)
    CaCl2 = bst.Chemical('CaCl2', phase='s', default=True)
    # ZnSO4 //// 7 H2O
    ZnSO4 = bst.Chemical('ZnSO4',search_ID='PubChem=24424', phase='s', default=True)
    # MnSO4 //// 7 H2O
    MnSO4 = bst.Chemical('MnSO4',search_ID='PubChem=24580', phase='s', default=True)
    # CoCl2 //// 6 H2O
    CoCl2 = bst.Chemical('CoCl2',search_ID='PubChem=24288', phase='s', default=True)
    # CuSO4 //// 5H2O
    CuSO4 = bst.Chemical('CuSO4',search_ID='PubChem=24462', phase='s', default=True)
    # NH4_6Mo7O24 //// 4 H2O
    NH4_6Mo7O24 = bst.Chemical('(NH4)6Mo7O24',search_ID='PubChem=61578', phase='s', default=True)
    # Na2B4O7 //// 10 H2O
    Na2B4O7 = bst.Chemical('Na2B4O7',search_ID='PubChem=10219853', phase='s', default=True)
    ####

    # chemicals = bst.Chemicals([
    #     Water,
    #     HCl, CaCl2, ZnSO4, MnSO4, CoCl2, CuSO4, NH4_6Mo7O24, Na2B4O7
    # ])
    # bst.settings.set_thermo(chemicals)
    # chemicals.define_group( # Source: https://www.thinkusadairy.org/products/permeate-(dairy-product-solids)/permeate-categories/whey-permeate#:~:text=Whey%20permeate%20(also%20called%20dairy%20product%20solids%2C,pleasant%20dairy%20flavor%20make%20whey%20permeate%20formulator%2Dfriendly.
    #     'minerals',
    #     #['Protein', 'Lactose', 'ButterFat', 'Ash', 'Water'],
    #     [100, 5*0.3, 2, 2.3, 2, 0.25, 1, 0.3, 0.5],
    #     wt=True,
    chemicals = bst.Chemicals([
        Water,
        HCl, CaCl2, ZnSO4, MnSO4, CoCl2, CuSO4, NH4_6Mo7O24, Na2B4O7
    ])
    bst.settings.set_thermo(chemicals)

    chemicals.define_group( # Source: https://www.thinkusadairy.org/products/permeate-(dairy-product-solids)/permeate-categories/whey-permeate#:~:text=Whey%20permeate%20(also%20called%20dairy%20product%20solids%2C,pleasant%20dairy%20flavor%20make%20whey%20permeate%20formulator%2Dfriendly.
        'minerals',
        ['Water', 'HCl', 'CaCl2', 'ZnSO4', 'MnSO4', 'CoCl2', 'CuSO4', '(NH4)6Mo7O24', 'Na2B4O7'],
        [1000, 5*0.3, 2, 2.3, 2, 0.25, 1, 0.3, 0.5],
        wt=True,
    )
    
    return chemicals

def create_antibiotics():
    # antibiotics
    ####
    ampicillin = bst.Chemical('Ampicillin', phase='s', default=True)
    kanamycin = bst.Chemical('Kanamycin', phase='s', default=True)
    streptomycin = bst.Chemical('Streptomycin', phase='s', default=True)
    chloramphenicol = bst.Chemical('Chloramphenicol', phase='s', default=True)
    ####
    chemicals = bst.Chemicals([
        ampicillin, kanamycin, streptomycin, chloramphenicol
    ])
    bst.settings.set_thermo(chemicals)

    return chemicals