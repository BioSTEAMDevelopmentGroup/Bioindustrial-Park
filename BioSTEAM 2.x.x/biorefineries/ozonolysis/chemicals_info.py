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
       chemicals.append(na_chemical)
       return na_chemical

chemicals = tmo.Chemicals(
    ['Water','Hydrogen_peroxide','Oleic_acid',
     'Nonanal','Nonanoic_acid','Azelaic_acid',
     tmo.Chemical('Epoxy_stearic_acid', search_ID='oxiraneoctanoic_acid,_3-octyl-'),
     Chemical('Naphtol', search_ID='90-15-3'),
     'Ethyl_acetate', 
     # 'Octane',
     # 'pentene',
     # 'Methylcyclohexane',
     # 'Cyclopentane',
     # tmo.Chemical('Trimethylpentene', search_ID='540-84-1'),
     # 'heptane',
     # 'pentane'
])

# Solid Catalyst not available in the database
# TODO.xxx Add ref for the below
Catalyst = create_new_chemical(
    'Phosphotungstic_acid',
    formula="H3PW12O40",
    MW=2880.2,
    CAS='1343-93-7',
    phase = 'l',
    rho = 2850 
    )
#https://www.fishersci.com/shop/products/phosphotungstic-acid-hydrate-thermo-scientific/AA4011614
#rho = 2.852g/cm3
#conversion = 

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

for chemical in chemicals: chemical.default()
chemicals.compile()
chemicals.set_synonym('Epoxy_stearic_acid', 'oxiraneoctanoic_acid,_3-octyl-') 
chemicals.show()

for i in chemicals:
    if not i.locked_state: i.V.g.method_P = 'IDEAL'







