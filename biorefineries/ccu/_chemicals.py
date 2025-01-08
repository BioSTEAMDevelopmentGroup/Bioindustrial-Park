# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 08:37:07 2024

@author: IGB
"""

import thermosteam as tmo



#%% 
# Constants

_cal2joule=4.184

#%% 

chems=tmo.Chemicals([])

database_chemicals_dict={}
copied_chemicals_dict={}
defined_chemicals_dict={}

def chemical_database(ID, search_ID=None, phase=None, **kwargs):
    chemical = tmo.Chemical(ID,search_ID=search_ID, **kwargs)
    if phase:
        chemical.at_state(phase)
        chemical.phase_ref = phase
    chems.append(chemical)
    database_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

def chemical_copied(ID,ref_chemical,**data):
    chemical=ref_chemical.copy(ID)
    chems.append(chemical)
    for i,j in data.items():setattr(chemical,i,j)
    copied_chemicals_dict[ID]=f'{ID}:{chemical.formula}/{chemical.MW}'
    return chemical

def chemical_defined(ID,**kwargs):
    chemical=tmo.Chemical.blank(ID,**kwargs)
    chems.append(chemical)
    defined_chemicals_dict[ID]=f'{ID}:{chemical.formula}/{chemical.MW}'
    return chemical

# %% 

# =============================================================================
# Create chemical objects available in database
# =============================================================================

H2O = chemical_database('H2O')

O2 = chemical_database('O2', phase='g', Hf=0)
N2 = chemical_database('N2', phase='g', Hf=0)
H2 = chemical_database('H2', Hf=0)
CH4 = chemical_database('Methane')
CO = chemical_database('CarbonMonoxide', Hf=-26400*_cal2joule)
CO2 = chemical_database('CO2')

CH3OH = chemical_database('Methanol')

HCOOH = chemical_database('HCOOH')
C18H39N = chemical_database('C18H39N') # amine
C19H41NO2 = chemical_database('C19H41NO2')

Triphenylphosphine = chemical_database('Triphenylphosphine', search_ID='603-35-0', phase='s') # surrogate model for DCPE

DCPE = chemical_database('DCPE', search_ID='23743-26-2', phase='s') # FA catalyst; second component
for i in DCPE.get_missing_properties():
    if i == 'V':
        try:
            DCPE.copy_models_from(Triphenylphosphine, [i])
        except:
            pass

phase_change_chemicals = ['H2O', 'H2', 'CH4', 'CO', 'CO2', 'CH3OH', 'HCOOH'
                          'C18H39N', 'C19H41NO2']

for chem in chems:
    if chem.ID in phase_change_chemicals: pass
    elif chem.locked_state: pass
    else: 
        # Set phase_ref to avoid missing model errors
        if chem.phase_ref == 'g':
            chem.at_state('g')
        if chem.phase_ref == 'l':
            chem.at_state('l')

chems.compile()
tmo.settings.set_thermo(chems, cache=True)

chems.set_synonym('H2O', 'Water')
chems.set_synonym('CO2', 'CarbonDioxide')
chems.set_synonym('CO', 'CarbonMonoxide')
chems.set_synonym('Methanol', 'CH3OH')