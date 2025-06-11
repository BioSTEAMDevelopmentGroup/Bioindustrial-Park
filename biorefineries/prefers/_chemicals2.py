# -*- coding: utf-8 -*-
"""
Created on 2025-04-16 15:19:55

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""
# %%
import biosteam as bst
import thermosteam as tmo
import numpy as np
from biosteam.utils import chemical_cache

__all__ = (
    # 'create_culture_media',
    # 'create_trace_metal_solution',
    # 'create_antibiotics',
    # 'create_feed1',
    # 'create_feed2',
    # 'create_leghemoglobin_chemicals',
    'create_chemical_LegH',
)
_cal2joule = 4.184 # auom('cal').conversion_factor('J')

chemical_groups = dict(
    Culture_Media = ('Yeast', 'H2O', 'Glucose', 'FeSO4', 'KH2PO4', '(NH4)2HPO4', 'KOH', 'MgSO4', 'CitricAcid'),
    Trace_Metal_Solution = ('H2O', 'HCl', 'CaCl2', 'ZnSO4', 'MnSO4', 'CoCl2', 'CuSO4', '(NH4)6Mo7O24', 'Na2B4O7'),
    Feed1 = ('H2O', 'Glucose', 'MgSO4'),
    Feed2 = ('H2O', 'Glucose', 'MgSO4', 'IPTG', '(NH4)2SO4', 'FeSO4', 'ZnSO4'),
    Leghemoglobin = ('heme_b', 'Protein'),

    Other_InsolubleSolids = ()
)

# %%
#
def create_chemical_LegH(set_thermo=True,yeast_includes_nitrogen=False):
    
    # create new chemical
    chems = bst.Chemicals([])
    def add_chemical(ID, ref=None, **data):
        chemical = bst.Chemical(ID, **data) if ref is None else ref.copy(ID, **data)
        chems.append(chemical)
        return chemical
    # Culture Media
    ####
    add_chemical('H2O')

    add_chemical('H2SO4', phase='l')
    ##### Gases #####
    O2 = add_chemical('O2', phase='g', Hf=0)
    N2 = add_chemical('N2', phase='g', Hf=0)
    CH4 = add_chemical('CH4', phase='g')
    CO = add_chemical('CO', search_ID='CarbonMonoxide', phase='g', Hf=-26400*_cal2joule)
    CO2 = add_chemical('CO2', phase='g')
    NH3 = add_chemical('NH3', phase='g', Hf=-10963*_cal2joule)
    NO = add_chemical('NO', search_ID='NitricOxide', phase='g')
    NO2 = add_chemical('NO2', phase='g')
    H2S = add_chemical('H2S', phase='g', Hf=-4927*_cal2joule)
    SO2 = add_chemical('SO2', phase='g')


    ##### Soluble inorganics #####
    # culture_media
    FeSO4 = add_chemical('FeSO4', phase='l', default=True)
    KH2PO4 = add_chemical('KH2PO4', phase='l', default=True)
    NH4_2HPO4 = add_chemical('(NH4)2HPO4', phase='l', default=True)
    MgSO4 = add_chemical('MgSO4', phase='l', default=True)
    KOH = add_chemical('KOH', phase='l', default=True)
    NaCl = add_chemical('NaCl', phase='l', default=True)
    # trace_metal_solution
    HCl = add_chemical('HCl', phase='l', default=True)
    CaCl2 = add_chemical('CaCl2', phase='s', default=True)
    ZnSO4 = add_chemical('ZnSO4', phase='s', default=True)
    MnSO4 = add_chemical('MnSO4', phase='s', default=True)
    CoCl2 = add_chemical('CoCl2', phase='s', default=True)
    CuSO4 = add_chemical('CuSO4', phase='s', default=True)
    NH4_6Mo7O24 = add_chemical('(NH4)6Mo7O24', search_ID='PubChem=61578', phase='s', default=True)
    Na2B4O7 = add_chemical('Na2B4O7', phase='s', default=True)

    # feed2
    NH4_2SO4 = add_chemical('(NH4)2SO4', phase='s', default=True, Hf=-288994*_cal2joule)
    # FeSO4
    # ZnSO4

    ##### Soluble organics #####
    # culture_media & feed
    CitricAcid = add_chemical('CitricAcid', phase='l', default=True)
    Glucose = add_chemical('Glucose', phase='l', default=True)
    IPTG = add_chemical('IPTG', phase='l', default=True)
    #add_chemical('Tryptone', phase='l', default=True)
    # antibiotics   
    Ampicillin = add_chemical('Ampicillin', phase='l', default=True)
    Kanamycin = add_chemical('Kanamycin', phase='l', default=True)
    Streptomycin = add_chemical('Streptomycin', phase='l', default=True)
    Chloramphenicol = add_chemical('Chloramphenicol', phase='l', default=True)

    Yeast = add_chemical(
        'Yeast',
        phase='l',
        formula='CH1.61O0.56N0.16' if yeast_includes_nitrogen else 'CH1.61O0.56',
        rho=1540,
        Cp=Glucose.Cp(298.15), # 1.54
        default=True,
        search_db=False,
        aliases=['Cellmass'],
    )
    Yeast.Hf = Glucose.Hf / Glucose.MW * Yeast.MW 
    # Same as glucose to ignore heats related to growth 

    ####
    # heme molecule formula: C34H32FeN4O4
    # heme_b = add_chemical('heme')
    Heme_b = add_chemical('Heme_b', search_ID='PubChem=26945', phase='s', default=True)

    # Only protein: 144 amino acids Formula:C729H1166N200O219S2​
    # https://www.uniprot.org/uniprotkb/P02236/entry#sequences
    protein_formula = {
        'H': 1166 / 729,
        'C': 729 / 729,
        'N': 200 / 729,
        'O': 219 / 729,
        'S': 2 / 729
    }
    formula = {i: round(j, 2) for i, j in protein_formula.items()}
    Protein = add_chemical(
        'Protein',
        search_db=False,
        default=True,
        atoms=formula,
        phase='s'
    )
    Leghemoglobin_formula = {
        'H': (1166+32) / 729,
        'C': (729+34) / 729,
        'N': (200 +4) / 729,
        'O': (219+4) / 729,
        'S': 2 / 729,
        'Fe': 1 / 729
    }
    formula2 = {i: round(j, 2) for i, j in Leghemoglobin_formula.items()}
    Leghemoglobin = add_chemical(
        'Leghemoglobin',
        search_db=False,
        default=True,
        atoms=formula2,
        phase='s'
    )

    # Default missing properties of chemicals to those of water
    for chemical in chems: chemical.default()

    ##### Group #####
    chems.compile()

    chems.define_group(
        'air',
        ['O2', 'N2'],
        [28, 72],
        wt=True
    )
    chems.define_group(
        'Culture_Media',
        ['Yeast', 'H2O', 'NaCl', 'Glucose', 'FeSO4', 'KH2PO4', '(NH4)2HPO4', 'KOH', 'MgSO4', 'CitricAcid'],
        [5, 1000, 10, 20, 0.05, 10, 4, 5, 1, 1.3],
        wt=True
    )

    chems.define_group(
        'Trace_Metal_Solution',
        ['H2O', 'HCl', 'CaCl2', 'ZnSO4', 'MnSO4', 
        'CoCl2', 'CuSO4', '(NH4)6Mo7O24', 'Na2B4O7'],
        [1000, 5*0.3, 2, 2.3, 2, 0.25, 1, 0.3, 0.5],
        wt=True,
    )

    chems.define_group(
        'Antibiotics',
        ['H2O', 'Ampicillin', 'Kanamycin', 'Streptomycin', 
        'Chloramphenicol', 'IPTG'],
        [1000, 0.05, 0.025, 0.05, 0.017, 0.01],
        wt=True,
    )

    chems.define_group(
        'Feed1',
        ['H2O', 'Glucose', 'MgSO4'],
        [1000, 800, 14],
        wt=True,
    )

    chems.define_group(
        'Feed2',
        ['H2O', 'Glucose', 'MgSO4', 'IPTG', '(NH4)2SO4', 
        'FeSO4', 'ZnSO4'],
        [1000, 800, 14, 0.01, 5, 4, 4],
        wt=True,
    )
    
    # chems.define_group(
    #     'Leghemoglobin',
    #     ['Heme_b', 'Protein'],
    #     np.array([1, 729]) / 729,
    #     wt=True,
    # )

    
    if set_thermo: bst.settings.set_thermo(chems)
    return chems
