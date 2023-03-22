# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autofunction:: biorefineries.cane.chemicals.create_sugarcane_chemicals
.. autofunction:: biorefineries.cane.chemicals.create_oilcane_chemicals
.. autofunction:: biorefineries.cane.chemicals.create_cellulosic_oilcane_chemicals
.. autofunction:: biorefineries.cane.chemicals.create_acyl_olein
.. autofunction:: biorefineries.cane.chemicals.create_acetyl_diolein

"""
import thermosteam as tmo
from thermosteam import functional as fn
from chemicals import atoms_to_Hill
from thermosteam.utils import chemical_cache
from biorefineries import cellulosic

__all__ = (
    'create_sugarcane_chemicals',
    'create_oilcane_chemicals',
    'create_cellulosic_oilcane_chemicals',
    'create_acyl_olein',
    'create_acetyl_diolein',
)

@chemical_cache
def create_sugarcane_chemicals():
    (Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, N2, CH4) = chemicals = tmo.Chemicals(
        ['Water', 'Ethanol', 'Glucose', 'Sucrose', 'H3PO4', 'P4O10',
         'CO2', 'Octane', 'O2', 'N2', 'CH4']
    )
    O2.at_state(phase='g')
    N2.at_state(phase='g')
    CH4.at_state(phase='g')
    CO2.at_state(phase='g')
    H3PO4.at_state(phase='s')
    P4O10.at_state(phase='s')
    Glucose.at_state(phase='s')
    Sucrose.at_state(phase='s')
    Glucose.N_solutes = 1
    Sucrose.N_solutes = 2
    
    def create_new_chemical(ID, phase='s', **constants):
        chemical = tmo.Chemical(ID, phase=phase, phase_ref=phase, search_db=False, **constants)
        chemicals.append(chemical)
        return chemical
    
    Ash = create_new_chemical('Ash', MW=1.)
    Cellulose = create_new_chemical('Cellulose',
                                    formula="C6H10O5", # Glucose monomer minus water
                                    Hf=-975708.8)
    Hemicellulose = create_new_chemical('Hemicellulose',
                                        formula="C5H8O4", # Xylose monomer minus water
                                        Hf=-761906.4)
    Flocculant = create_new_chemical('Flocculant',
                                     MW=1.)
    Lignin = create_new_chemical('Lignin',
                                 formula='C8H8O3', # Vainillin
                                 Hf=-452909.632)
    Solids = create_new_chemical('Solids', MW=1.)
    Yeast = create_new_chemical(
        'Yeast', 
        formula='CH1.61O0.56N0.16',
        rho=1540,
        default=True,
        Hf=-130412.73,
    )
    CaO = create_new_chemical('CaO', formula='CaO')

    
    ### Fill missing properties ###
    
    # Insolubles occupy a significant volume
    insoluble_solids = (Ash, Cellulose, Hemicellulose,
                        Flocculant, Lignin, Solids, Yeast, P4O10)
    
    # Solubles don't occupy much volume
    soluble_solids = (CaO, H3PO4, Glucose, Sucrose) 
    
    for chemical in insoluble_solids:
        V = fn.rho_to_V(rho=1540, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
    
    for chemical in soluble_solids:
        V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
    
    # Add constant models for molar heat capacity of solids
    # Lignocellulosic heat capacities:
    # https://link.springer.com/article/10.1007/s10853-013-7815-6
    # https://www.sciencedirect.com/science/article/pii/0032386182901252
    Ash.Cn.add_model(0.09 * 4.184 * Ash.MW) 
    CaO.Cn.add_model(1.02388 * CaO.MW) 
    Cellulose.Cn.add_model(1.364 * Cellulose.MW) 
    Hemicellulose.Cn.add_model(1.364 * Hemicellulose.MW)
    Flocculant.Cn.add_model(4.184 * Flocculant.MW)
    Lignin.Cn.add_model(1.364 * Lignin.MW)
    Solids.Cn.add_model(1.100 * Solids.MW)
    
    for chemical in chemicals: chemical.default()
    
    chemicals.compile()
    chemicals.set_synonym('Water', 'H2O')
    chemicals.set_synonym('Yeast', 'DryYeast')
    
    return chemicals

@chemical_cache
def create_acyl_olein(N_acyl):
    # Assume properties are similar between oleic and palmitic chains
    Hf_oleic_acid = -764.80e3
    Hf_triolein = -1776e3
    Hf_water = -285.825e3
    Hf_glycerol = -668.6e3
    # TAG + 3H2O -> GLY + 3FFA
    # Hrxn = (GLY + 3FFA) - (TAG + 3H2O)
    Hf_tag_to_ffa = (Hf_glycerol + 3 * Hf_oleic_acid) - (Hf_triolein + 3 * Hf_water)
    Hf_minus_one_acyl = Hf_tag_to_ffa / 3. 
    if N_acyl == 0:
        ID_model = 'Palmitic acid'
        ID = 'OleicAcid'
        Hf = Hf_oleic_acid
    elif N_acyl == 1:
        ID_model = 'Monopalmitin'
        ID = 'MonoOlein'
        # TAG + 2H2O -> MAG + 2FFA
        # Hrxn = MAG + 2FFA - (TAG + 2H2O) # Assume to be Hf_tag_to_ffa * 2 / 3
        # MAG = Hrxn - 2FFA + TAG + 2H2O
        Hf = 2 * Hf_minus_one_acyl - 2 * Hf_oleic_acid + Hf_triolein + 2 * Hf_water
    elif N_acyl == 2:
        ID_model = 'Dipalmitin'
        ID = 'DiOlein'
        # TAG + H2O -> DAG + FFA
        # Hrxn = DAG + FFA - (TAG + H2O) # Assume to be Hf_tag_to_ffa / 3
        # DAG = Hrxn - FFA + TAG + H2O
        Hf = Hf_minus_one_acyl - Hf_oleic_acid + Hf_triolein + Hf_water
    elif N_acyl == 3:
        ID_model = 'Tripalmitin'
        ID = 'TriOlein'
        Hf = Hf_triolein
    model = tmo.Chemical(ID_model).at_state(phase='l', copy=True)
    atoms = model.atoms
    atoms['C'] += N_acyl * 2
    atoms['H'] += N_acyl * 2
    formula = atoms_to_Hill(atoms)
    chemical = tmo.Chemical(ID, search_db=False, phase='l', phase_ref='l',
                            formula=formula, Hf=Hf)
    chemical._Dortmund = group_counts = model.Dortmund.copy()
    group_counts.set_group_counts_by_name({'CH=CH': N_acyl}, reset=False)
    chemical.V.add_model(fn.rho_to_V(rho=900, MW=chemical.MW))
    chemical.Cn.add_model(model.Cp(298.15) * chemical.MW)
    chemical.copy_models_from(model, ['mu', 'sigma', 'kappa'])
    if ID == 'MonoOlein': chemical.mu.add_model(lambda T: 0.0001)
    return chemical

@chemical_cache
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

@chemical_cache
def create_oilcane_chemicals():
    chemicals = create_sugarcane_chemicals().copy()
    (Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, N2, CH4, 
     Ash, Cellulose, Hemicellulose, Flocculant, Lignin, Solids, DryYeast, CaO) = chemicals
    
    ### Define common chemicals ###
    Biodiesel = tmo.Chemical('Biodiesel',
                             search_ID='Methyl oleate')
    Methanol = tmo.Chemical('Methanol')
    Glycerol = tmo.Chemical('Glycerol')
    chemicals.extend([Biodiesel, Methanol, Glycerol])
    
    ### Define new chemicals ###
    
    def create_new_chemical(ID, phase='s', **constants):
        solid = tmo.Chemical.blank(ID, phase=phase, **constants)
        chemicals.append(solid)
        return solid
    
    HCl = create_new_chemical('HCl', formula='HCl')
    NaOH = create_new_chemical('NaOH', formula='NaOH')
    NaOCH3 = create_new_chemical('NaOCH3', formula='NaOCH3')
    ### Fill missing properties ###
    
    # Solubles don't occupy much volume
    for chemical in (HCl, NaOH):
        V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
    
    # Assume sodium methoxide has some of the same properties as methanol
    LiquidMethanol = Methanol.at_state(phase='l', copy=True)
    NaOCH3.copy_models_from(LiquidMethanol, ['V', 'sigma', 'kappa', 'Cn'])
    chemicals.extend([
        tmo.Chemical('Phosphatidylinositol', formula='C47H83O13P',
                     search_db=False, CAS='383907-36-6', default=True,
                     Hf=-1.779e6, # Assume same as TAG on a weight basis
                     aliases={'PL', 'PolarLipid'}, phase='l'),
        # tmo.Chemical('BetaSitosterol', search_ID='83-46-5',
        #               phase='l', Hf=-1000. * (533.79 + 93.90),
        #               aliases=['Sterol']),
        create_acyl_olein(0),
        create_acyl_olein(1),
        create_acyl_olein(2),
        create_acyl_olein(3),
        tmo.Chemical('Acetone')
    ])
    for chemical in chemicals: chemical.default()
    chemicals.compile()
    chemicals.set_synonym('OleicAcid', 'FFA')
    chemicals.set_synonym('MonoOlein', 'MAG')
    chemicals.set_synonym('DiOlein', 'DAG')
    chemicals.set_synonym('TriOlein', 'TAG')
    chemicals.define_group('Lipid', ('PL', 'FFA', 'MAG', 'DAG', 'TAG'))
    chemicals.define_group('Oil', ('PL', 'FFA', 'MAG', 'DAG', 'TAG'))
    chemicals.set_synonym('Water', 'H2O')
    chemicals.set_synonym('Yeast', 'DryYeast')
    return chemicals

@chemical_cache
def create_cellulosic_oilcane_chemicals():
    oilcane_chemicals = create_oilcane_chemicals()
    cellulosic_chemicals = cellulosic.create_cellulosic_ethanol_chemicals()
    removed = {'SuccinicAcid', 'H2SO4', 'Z_mobilis', 'Oil', 'Yeast'}
    chemicals = tmo.Chemicals([
        i for i in (oilcane_chemicals.tuple + cellulosic_chemicals.tuple) if i.ID not in removed
    ])
    chemicals.extend([
        create_acetyl_diolein(),
        tmo.Chemical('Urea', default=True, phase='l'),
        chemicals.Glucose.copy('Yeast'),
    ])
    chemicals.compile()
    chemicals.set_synonym('AcetylDiOlein', 'AcTAG')
    chemicals.define_group('Lipid', ['PL', 'FFA', 'MAG', 'DAG', 'TAG', 'AcTAG'])
    chemicals.define_group('lipid', ['PL', 'FFA', 'MAG', 'DAG', 'TAG', 'AcTAG'])
    chemicals.define_group('Oil', ['PL', 'FFA', 'MAG', 'DAG', 'TAG'])
    chemicals.set_synonym('Yeast', 'DryYeast')
    chemicals.set_synonym('Yeast', 'Cellmass')
    chemicals.AcetylDiOlein.V.method_P = chemicals.TriOlein.V.method_P = None
    chemicals.set_synonym('Cellmass', 'Cells')
    chemicals.set_synonym('Cellmass', 'cellmass')
    chemicals.set_synonym('TriOlein', 'TAG')
    return chemicals