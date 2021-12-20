# -*- coding: utf-8 -*-
"""
"""
import thermosteam as tmo
import biosteam as bst
from biosteam import Unit, Stream, settings, main_flowsheet

#%% Chemicals - Definition
# Create chemicals here
ozo_chemicals = tmo.Chemicals(
    ['Water','Hydrogen_peroxide','Oleic_acid',
     'Nonanal','Nonanoic_acid','Azelaic_acid'])

(Water,Hydrogen_peroxide,Oleic_acid,
Nonanal,Nonanoic_acid,Azelaic_acid) = ozo_chemicals


def create_new_chemical(ID, phase='s', **constants):
    # Create a new solid chemical without any data
    solid = tmo.Chemical(ID, search_db=False, phase=phase, **constants)

    # Add chemical to the Chemicals object
    ozo_chemicals.append(solid)

    return solid
#Writing chemicals not defined in the database and other additional chemicals

Catalyst = create_new_chemical(
    'Phosphotungstic_acid',
    formula="H3PW12O40", # Chemical Formula
    MW=2880.2,
    CAS='1343-93-7'
    )

Oxononanoic_acid = create_new_chemical('Oxononanoic_acid',
    phase='l',
    Hf=-579480,
    formula = 'C9H16O3',
    MW = 172.22,
    CAS = '2553-17-5'
)

Epoxide = create_new_chemical('Epoxy_stearic_acid',
    phase='l',
    #Hf,
    formula = 'C18H34O3',
    MW = 298.5,
    CAS = '2443-39-2'
 )
for chemical in ozo_chemicals: chemical.default()
 

 
     





