# -*- coding: utf-8 -*-
"""
"""
import thermosteam as tmo
from thermosteam import Chemical
import biosteam as bst
from biosteam import Unit, Stream, settings, main_flowsheet

# Chemicals - Definition

#For unavailable chemicals
def create_new_chemical(ID, phase = 's', **constants):
       na_chemical = tmo.Chemical(ID, search_db=False,
                                  phase=phase, **constants)
       ozo_chemicals.append(na_chemical)
       return na_chemical

#Chemicals that are available in the database
ozo_chemicals = tmo.Chemicals(
    ['Water','Hydrogen_peroxide','Oleic_acid',
     'Nonanal','Nonanoic_acid','Azelaic_acid',
     'oxiraneoctanoic_acid,_3-octyl-'])

#Solid Catalyst not available in the database
#TODO.xxx Add ref for the below
Catalyst = create_new_chemical(
    'Phosphotungstic_acid',
    formula="H3PW12O40",
    MW=2880.2,
    CAS='1343-93-7'
    )

#Liquid chemical not available in the database
#TODO.xxx Add ref for the below
Oxononanoic_acid = create_new_chemical(
   'Oxononanoic_acid',
    phase='l',
    Hf=-579480,
    formula = 'C9H16O3',
    MW = 172.22,
    CAS = '2553-17-5'
)

for chemical in ozo_chemicals: chemical.default()
ozo_chemicals.compile()
ozo_chemicals.set_synonym('oxiraneoctanoic_acid,_3-octyl-' ,'Epoxy_stearic_acid') 
ozo_chemicals.show()

# =============================================================================
#Tried below alternatives for Oxononanoic acid
#Recognised CAS but does not have Tb, P or Hf
#Methyl_oxo_nonanoate = Chemical('1931-63-1')
#Methyl_oxo_nonanoate.show()
# =============================================================================






