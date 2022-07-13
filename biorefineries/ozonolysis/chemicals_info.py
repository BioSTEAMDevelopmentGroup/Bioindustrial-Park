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
       na_chemical = tmo.Chemical(ID,
                                  search_db=False,
                                  phase=phase,
                                  **constants)
       ozo_chemicals.append(na_chemical)
       return na_chemical

#Chemicals that are available in the database
DHSA = tmo.Chemical('DHSA', search_ID = '120-87-6')
ozo_chemicals = tmo.Chemicals(
    ['Water','Hydrogen_peroxide','Oleic_acid',
     'Nonanal','Nonanoic_acid','Azelaic_acid',
     'oxiraneoctanoic_acid,_3-octyl-','Hexane',
     'Ethyl_acetate','Octane',
     'pentene',
     'Methylcyclohexane',
     'Cyclopentane',
     '540-84-1',
     'heptane',
     'pentane',
     'Sulphuric_acid',
      DHSA])

#Solid Catalyst not available in the database
#TODO.xxx Add ref for the below
Catalyst = create_new_chemical(
        'Phosphotungstic_acid',
        # formula="H3PW12O40",
        # MW=2880.2,
        #CAS = '7783-03-1',
        CAS='1343-93-7',
       phase = 'l',
       # rho = 2850 
    )

#https://www.fishersci.com/shop/products/phosphotungstic-acid-hydrate-thermo-scientific/AA4011614
#rho = 2.852g/cm3
#conversion = 

#####I think the chemmical is actually this###
# https://www.sigmaaldrich.com/US/en/product/aldrich/455970

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

Dortmund = ozo_chemicals.Oxononanoic_acid.Dortmund
Dortmund.set_group_counts_by_name(dict(COOH=1, CH2=7, CHO=1))
ozo_chemicals.compile()
ozo_chemicals.set_synonym('oxiraneoctanoic_acid,_3-octyl-' ,'Epoxy_stearic_acid') 
ozo_chemicals.show()

for i in ozo_chemicals:
    if not i.locked_state: i.V.g.method_P = 'IDEAL'







