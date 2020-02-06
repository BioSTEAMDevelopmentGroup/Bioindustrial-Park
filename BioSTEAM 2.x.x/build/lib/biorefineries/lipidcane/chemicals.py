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
     'Glucose', 'Sucrose', 'H3PO4', 'CO2',
     'Octane', Biodiesel])

(Water, Methanol, Ethanol,
 Glycerol, Glucose, Sucrose,
 H3PO4, CO2, Octane, Biodiesel) = lipidcane_chemicals

CO2.at_state(phase='g')
H3PO4.at_state(phase='s')
Glucose.at_state(phase='s')
Sucrose.at_state(phase='s')
Glucose.phase_ref = 's'
Glucose.load_free_energies()

# %% Define new chemicals

def create_new_chemical(ID, phase='s', **constants):
    solid = tmo.Chemical.blank(ID, phase=phase, **constants)
    lipidcane_chemicals.append(solid)
    return solid

Ash = create_new_chemical('Ash')
Cellulose = create_new_chemical('Cellulose')
Hemicellulose = create_new_chemical('Hemicellulose')
Flocculant = create_new_chemical('Flocculant')
Lignin = create_new_chemical('Lignin')
Solids = create_new_chemical('Solids')
DryYeast = create_new_chemical('DryYeast', CAS='Yeast')
CaO = create_new_chemical('CaO', MW=56.0774)
HCl = create_new_chemical('HCl', MW=36.46094)
NaOH = create_new_chemical('NaOH', MW=39.997109)
NaOCH3 = create_new_chemical('NaOCH3', MW=54.023689)

Lipid = create_new_chemical(
    'Lipid',
    phase='l',
    MW=885.432,
)


# %% Fill missing properties

# Assume properties are similar for trioleate and tripalmitin
Tripalmitin = tmo.Chemical('Tripalmitin').at_state(phase='l', copy=True)
Lipid.copy_missing_slots_from(Tripalmitin, slots=['V', 'sigma',
                                                  'kappa', 'Cn'])

# Assume a constant volume for lipid
lipid_molar_volume = fn.rho_to_V(rho=900, MW=Lipid.MW)
Lipid.V.add_model(lipid_molar_volume)

# Insolubles occupy a significant volume
insoluble_solid_molar_volume = fn.rho_to_V(rho=1540, MW=1.)
insoluble_solids = (Ash, Cellulose, Hemicellulose, Sucrose,
                    Flocculant, Lignin, Solids, DryYeast)

# Solubles don't occupy much volume
soluble_solids = (CaO, HCl, NaOH, Glucose) 

for chemical in insoluble_solids:
    chemical.V.add_model(insoluble_solid_molar_volume)

for chemical in soluble_solids:
    V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
    chemical.V.add_model(V, top_priority=True)

# Assume sodium methoxide has some of the same properities as methanol
LiquidMethanol = Methanol.at_state(phase='l', copy=True)
NaOCH3.copy_missing_slots_from(LiquidMethanol, slots=['V', 'sigma',
                                                      'kappa', 'Cn',
                                                      'H', 'S'])

# Add constant models for molar heat capacity of solids
Ash.Cn.add_model(0.09 * 4.184) 
CaO.Cn.add_model(1.02388 * 56.0774) 
Cellulose.Cn.add_model(1.364) 
Hemicellulose.Cn.add_model(1.364)
Flocculant.Cn.add_model(4.184)
Lignin.Cn.add_model(1.364)
Solids.Cn.add_model(1.100)

# Set heats of formation and combustion (for those that need it)
Hemicellulose.Hc = 17e3
Cellulose.Hc = 17e3
Lignin.Hc = 21e3
Glucose.Hc = 2.8e6
Sucrose.Hc = 5.7e6
Lipid.Hc = 35e6

Ethanol.Hf = -277.690e3
Water.Hf = -241.820e3
CO2.Hf = -393.520e3 # in gas
Glucose.Hf = -1274e3
Sucrose.Hf = -2221.2e3
Octane.Hf = -250e3
H3PO4.Hf = -1271.66e3
Lipid.Hf = -2.1937e3
Biodiesel.Hf = -0.72764e3
Glycerol.Hf = -0.6696e3

for chemical in lipidcane_chemicals:
    chemical.default()

lipidcane_chemicals.compile()

pretreatment_chemicals = lipidcane_chemicals.subgroup(
    ['Lipid', 'Ash', 'Cellulose',
     'Flocculant', 'Hemicellulose',
     'Lignin', 'Solids', 'CaO',
     'Water', 'Ethanol',
     'Glucose', 'Sucrose',
     'DryYeast', 'H3PO4']
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