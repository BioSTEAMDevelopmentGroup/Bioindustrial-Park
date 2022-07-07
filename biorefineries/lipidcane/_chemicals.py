# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 06:42:02 2020

@author: yoelr
"""
from thermosteam import functional as fn
from chemicals import atoms_to_Hill
import thermosteam as tmo

__all__ = ('create_chemicals',)

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

def create_chemicals():
    from biorefineries import sugarcane as sc
    (Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, N2, CH4, 
     Ash, Cellulose, Hemicellulose, Flocculant, Lignin, Solids, DryYeast, CaO) = sc.chemicals
    
    ### Define common chemicals ###
    Biodiesel = tmo.Chemical('Biodiesel',
                             search_ID='Methyl oleate')
    Methanol = tmo.Chemical('Methanol')
    Glycerol = tmo.Chemical('Glycerol')
    lipidcane_chemicals = sc.chemicals.copy()
    lipidcane_chemicals.extend([Biodiesel, Methanol, Glycerol])
    
    ### Define new chemicals ###
    
    def create_new_chemical(ID, phase='s', **constants):
        solid = tmo.Chemical.blank(ID, phase=phase, **constants)
        lipidcane_chemicals.append(solid)
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
    lipidcane_chemicals.extend([
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
    for chemical in lipidcane_chemicals: chemical.default()
    lipidcane_chemicals.compile()
    lipidcane_chemicals.set_synonym('OleicAcid', 'FFA')
    lipidcane_chemicals.set_synonym('MonoOlein', 'MAG')
    lipidcane_chemicals.set_synonym('DiOlein', 'DAG')
    lipidcane_chemicals.set_synonym('TriOlein', 'TAG')
    lipidcane_chemicals.define_group('Lipid', ('PL', 'FFA', 'MAG', 'DAG', 'TAG'))
    lipidcane_chemicals.define_group('Oil', ('PL', 'FFA', 'MAG', 'DAG', 'TAG'))
    lipidcane_chemicals.set_synonym('Water', 'H2O')
    lipidcane_chemicals.set_synonym('Yeast', 'DryYeast')
    return lipidcane_chemicals

