# -*- coding: utf-8 -*-
"""
"""
import thermosteam as tmo
from thermosteam import Chemical
import biosteam as bst
from biosteam import Unit, Stream, settings, main_flowsheet

# Chemicals - Definition

#These chemicals are available in the database
ozo_chemicals = tmo.Chemicals(
    ['Water','Hydrogen_peroxide','Oleic_acid',
     'Nonanal','Nonanoic_acid','Azelaic_acid'])

#For unavailable solids
def create_new_solid_chemical(ID, phase='s', **constants):
    # Create a new solid chemical without any data
    solid = tmo.Chemical(ID, search_db=False, phase=phase, **constants)

    # Add chemical to the Chemicals object
    ozo_chemicals.append(solid)

    return solid

#Solid Catalyst not available in the database
#TODO.xxx Add ref for the below
Catalyst = create_new_solid_chemical(
    'Phosphotungstic_acid',
    formula="H3PW12O40", # Chemical Formula
    MW=2880.2,
    CAS='1343-93-7'
    )
#For unavailable liquids
def create_new_liquid_chemical(ID, phase='l', **constants):
    # Create a new liquid chemical without any data
    liquid = tmo.Chemical(ID, search_db=False, phase=phase, **constants)

    # Add chemical to the Chemicals object
    ozo_chemicals.append(liquid)

    return liquid

#Liquid chemical not available in the database
#TODO.xxx Add ref for the below
Oxononanoic_acid = create_new_liquid_chemical('Oxononanoic_acid',
    phase='l',
    Hf=-579480,
    formula = 'C9H16O3',
    MW = 172.22,
    CAS = '2553-17-5'
)
# =============================================================================
#Tried below alternatives for Oxononanoic acid
#Recognised CAS but does not have Tb, P or Hf
#Methyl_oxo_nonanoate = Chemical('1931-63-1')
#Methyl_oxo_nonanoate.show()
# =============================================================================

# =============================================================================
# Adding data to the existing chemicals 
def updating_existing_chemical(ID, phase='l', **constants):
    # Update an existing chemical without any data
    existing_chemical = tmo.Chemical(ID, search_db=True, phase=phase, **constants)
    
    # Add underdefined chemical to the Chemicals object
    ozo_chemicals.append(existing_chemical)

    return existing_chemical

#TODO.xxx Find Hf for the below
#Reference for boiling point:https://www.guidechem.com/encyclopedia/trans-9-10-epoxyoctadecanoic-a-dic125128.html
#Reference for the compound: https://pubchem.ncbi.nlm.nih.gov/compound/15868
Epoxy_stearic_acid = updating_existing_chemical(
    '9,10-Epoxyoctadecanoic_acid',
     phase='l',
     #Hf,
    formula = 'C18H34O3',
    CAS = '2443-39-2',
    Tb = 696.05,
   )
Epoxy_stearic_acid.show()

# =============================================================================
#Prospective alternative chemicals for epoxide 
#Has tb, lacks Pt and Hf
#cis_epoxy_stearic_acid = Chemical('24560-98-3')
#cis_epoxy_stearic_acid.show()

#Tried the below did't work
#Oxiraneoctanoic acid, 3-octyl-, trans= '13980-07-9'
#Does not P, Hf or Tb
#Oxiran = Chemical('13980-07-9')
#Oxiran.show()
#a valid CAS number was recognized, but its not in the database
#Methyl_epoxystearate =Chemical('2500-59-6')
#Methyl_epoxystearate.show()
#epoxy_stearic_acid_methyl_ester = Chemical('6084-76-0')
#=============================================================================

for chemical in ozo_chemicals: chemical.default()
ozo_chemicals.compile()
#TODO.XXX Change name for Epoxy_stearic_acid 
ozo_chemicals.set_synonym('8-(3-octyloxiran-2-yl)octanoic acid','Epoxy_stearic_acid')
ozo_chemicals.show()

# =============================================================================
#IS THIS NEEDED? 
#(Water,Hydrogen_peroxide,Oleic_acid,
# Nonanal,Nonanoic_acid,Azelaic_acid) = ozo_chemicals
# =============================================================================
     





