# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 06:42:02 2020

@author: yoelr
"""
from thermosteam import functional as fn
import thermosteam as tmo

__all__ = ('lipidcane_chemicals', 'pretreatment_chemicals',
           'ethanol_chemicals', 'biodiesel_chemicals')


# %% Define common chemicals
Biodiesel = tmo.Chemical('Biodiesel',
                         search_ID='Methyl oleate')
lipidcane_chemicals = tmo.Chemicals(
    ['Water', 'Methanol', 'Ethanol', 'Glycerol',
     'Glucose', 'Sucrose', 'H3PO4', 'P4O10', 'CO2',
     'Octane', 'O2', Biodiesel])

(Water, Methanol, Ethanol,
 Glycerol, Glucose, Sucrose,
 H3PO4, P4O10, CO2, Octane, O2, Biodiesel) = lipidcane_chemicals

O2.at_state(phase='g')
CO2.at_state(phase='g')
H3PO4.at_state(phase='s')
P4O10.at_state(phase='s')
Glucose.at_state(phase='s')
Sucrose.at_state(phase='s')

# %% Define new chemicals

def create_new_chemical(ID, phase='s', **constants):
    solid = tmo.Chemical.blank(ID, phase=phase, **constants)
    lipidcane_chemicals.append(solid)
    return solid

Ash = create_new_chemical('Ash', MW=1.)
Cellulose = create_new_chemical('Cellulose',
                                formula="C6H10O5", # Glucose monomer
                                MW=162.14,
                                Hf=-975708.8)
Hemicellulose = create_new_chemical('Hemicellulose',
                                    formula="C5H8O5", # Xylose monomer
                                    MW=132.12,
                                    Hf=-761906.4)
Flocculant = create_new_chemical('Flocculant',
                                 MW=1.)
Lignin = create_new_chemical('Lignin',
                             formula='C8H8O3',
                             MW=152.15,
                             Hf=-452909.632)
Solids = create_new_chemical('Solids', MW=1.)
DryYeast = create_new_chemical('DryYeast', MW=1., CAS='Yeast')
CaO = create_new_chemical('CaO', MW=56.0774)
HCl = create_new_chemical('HCl', MW=36.46094)
NaOH = create_new_chemical('NaOH', MW=39.997109)
NaOCH3 = create_new_chemical('NaOCH3', MW=54.023689)

Lipid = create_new_chemical(
    'Lipid',
    phase = 'l',
    MW = 885.432,
    formula = 'C57H104O6',
    Hf = -2193.7e3
)
Lipid.load_combustion_data()

# %% Fill missing properties

# Assume properties are similar for trioleate and tripalmitin
Tripalmitin = tmo.Chemical('Tripalmitin').at_state(phase='l', copy=True)
Lipid.copy_missing_slots_from(Tripalmitin, slots=['V', 'sigma',
                                                  'kappa', 'Cn'])

# Assume a constant volume for lipid
lipid_molar_volume = fn.rho_to_V(rho=900, MW=Lipid.MW)
Lipid.V.add_model(lipid_molar_volume)

# Insolubles occupy a significant volume
insoluble_solids = (Ash, Cellulose, Hemicellulose, Sucrose,
                    Flocculant, Lignin, Solids, DryYeast, P4O10)

# Solubles don't occupy much volume
soluble_solids = (CaO, HCl, NaOH, H3PO4, Glucose) 

for chemical in insoluble_solids:
    V = fn.rho_to_V(rho=1540, MW=chemical.MW)
    chemical.V.add_model(V, top_priority=True)

for chemical in soluble_solids:
    V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
    chemical.V.add_model(V, top_priority=True)

# Assume sodium methoxide has some of the same properities as methanol
LiquidMethanol = Methanol.at_state(phase='l', copy=True)
NaOCH3.copy_missing_slots_from(LiquidMethanol, slots=['V', 'sigma',
                                                      'kappa', 'Cn',
                                                      'H', 'S'])

# Add constant models for molar heat capacity of solids
Ash.Cn.add_model(0.09 * 4.184 * Ash.MW) 
CaO.Cn.add_model(1.02388 * CaO.MW) 
Cellulose.Cn.add_model(1.364 * Cellulose.MW) 
Hemicellulose.Cn.add_model(1.364 * Hemicellulose.MW)
Flocculant.Cn.add_model(4.184 * Flocculant.MW)
Lignin.Cn.add_model(1.364 * Lignin.MW)
Solids.Cn.add_model(1.100 * Solids.MW)

for chemical in lipidcane_chemicals:
    chemical.default()

lipidcane_chemicals.compile()
lipidcane_chemicals.set_synonym('Water', 'H2O')

pretreatment_chemicals = lipidcane_chemicals.subgroup(
    ['Lipid', 'Ash', 'Cellulose',
     'Flocculant', 'Hemicellulose',
     'Lignin', 'Solids', 'CaO',
     'Water', 'Ethanol',
     'Glucose', 'Sucrose',
     'DryYeast', 'H3PO4', 'O2', 'CO2', 'P4O10']
)
ethanol_chemicals = lipidcane_chemicals.subgroup(
    ['Water', 'CO2', 'Ethanol',
     'Glucose', 'Sucrose', 'DryYeast',
     'H3PO4', 'Octane']
)
biodiesel_chemicals = lipidcane_chemicals.subgroup(
    ['Lipid', 'Biodiesel', 'Methanol', 'Glycerol', 'Water',
     'HCl', 'NaOH', 'NaOCH3']
)