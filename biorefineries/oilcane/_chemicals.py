# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 06:42:02 2020

@author: yoelr
"""
import thermosteam as tmo
from chemicals import atoms_to_Hill
from thermosteam import functional as fn

__all__ = ('create_chemicals',)

def create_acetyl_diolein():
    # Assume properties are similar between oleic and palmitic chains
    Hf_oleic_acid = -764.80e3
    Hf_triolein = -1776e3
    Hf_acetic_acid = -483.58e3
    ID_model = 'Palmitin'
    ID = 'AcetylDiOlein'
    # Ac + TAG -> acTAG + FFA
    # Hrxn = acTAG + FFA - (TAG + Ac)  # Assume Hrxn is negligible
    # acTAG = Hrxn - FFA + TAG + Ac
    Hf = Hf_triolein - Hf_oleic_acid + Hf_acetic_acid
    model = tmo.Chemical(ID_model).at_state(phase='l', copy=True)
    formula = atoms_to_Hill({'C': 41, 'H':75, 'O': 6})
    chemical = tmo.Chemical(ID, search_db=False, phase='l', phase_ref='l',
                            formula=formula, Hf=Hf)
    chemical._Dortmund = group_counts = model.Dortmund.copy()
    group_counts.set_group_counts_by_name({'CH=CH': 2, 'CH':1, 'CH2': 28, 'CH3': 3, 'CH2COO': 2, 'COOH': 1}, reset=True)
    chemical.V.add_model(fn.rho_to_V(rho=900, MW=chemical.MW))
    chemical.mu.add_method(900 * 20.3e-6)
    chemical.Cn.add_method(model.Cp(298.15) * chemical.MW)
    chemical.copy_models_from(model, ['sigma', 'kappa'])
    return chemical

def create_chemicals():
    from biorefineries import lipidcane as lc
    from biorefineries import cornstover as cs
    removed = {'SuccinicAcid', 'H2SO4', 'Z_mobilis', 'Oil'}
    chemicals = tmo.Chemicals([
        i for i in (lc.chemicals.tuple + cs.chemicals.tuple) if i.ID not in removed
    ])
    chemicals.extend([
        create_acetyl_diolein(),
        tmo.Chemical('Urea', default=True, phase='l'),
    ])
    chemicals.compile()
    chemicals.set_synonym('AcetylDiOlein', 'AcTAG')
    chemicals.define_group('Lipid', ['PL', 'FFA', 'MAG', 'DAG', 'TAG', 'AcTAG'])
    chemicals.define_group('lipid', ['PL', 'FFA', 'MAG', 'DAG', 'TAG', 'AcTAG'])
    chemicals.define_group('Oil', ['PL', 'FFA', 'MAG', 'DAG', 'TAG'])
    chemicals.set_synonym('DryYeast', 'Cellmass')
    chemicals.AcetylDiOlein.V.method_P = chemicals.TriOlein.V.method_P = None
    chemicals.set_synonym('Cellmass', 'Cells')
    chemicals.set_synonym('Cellmass', 'cellmass')
    chemicals.set_synonym('TriOlein', 'TAG')
    return chemicals
    
