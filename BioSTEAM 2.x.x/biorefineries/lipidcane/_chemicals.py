# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 06:42:02 2020

@author: yoelr
"""
from thermosteam import functional as fn
import thermosteam as tmo

__all__ = ('create_chemicals',)

def create_chemicals():
    from biorefineries import sugarcane as sc
    (Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, CH4, 
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
    Lipid = create_new_chemical(
        'Lipid',
        phase = 'l',
        formula = 'C57H104O6',
        Hf = -2193.7e3
    )
    Lipid.Dortmund.set_group_counts_by_name({'CH3':3, 'CH2':41, 'CH':1, 'CH=CH':3, 'CH2COO':3})
    
    ### Fill missing properties ###
    
    # Assume properties are similar for trioleate and tripalmitin
    Tripalmitin = tmo.Chemical('Tripalmitin').at_state(phase='l', copy=True)
    Lipid.copy_models_from(Tripalmitin, ['V', 'sigma', 'kappa', 'Cn'])
    
    # Assume a constant volume for lipid
    lipid_molar_volume = fn.rho_to_V(rho=900, MW=Lipid.MW)
    Lipid.V.add_model(lipid_molar_volume)
    
    # Solubles don't occupy much volume
    for chemical in (HCl, NaOH):
        V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
    
    # Assume sodium methoxide has some of the same properities as methanol
    LiquidMethanol = Methanol.at_state(phase='l', copy=True)
    NaOCH3.copy_models_from(LiquidMethanol, ['V', 'sigma', 'kappa', 'Cn'])
    
    for chemical in lipidcane_chemicals: chemical.default()
    
    lipidcane_chemicals.compile()
    lipidcane_chemicals.set_synonym('Water', 'H2O')
    lipidcane_chemicals.set_synonym('Yeast', 'DryYeast')
    return lipidcane_chemicals

