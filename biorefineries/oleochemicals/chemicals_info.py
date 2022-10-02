# -*- coding: utf-8 -*-
"""
"""
import thermosteam as tmo
from thermosteam import Chemical
import biosteam as bst
from biosteam import Unit, Stream, settings, main_flowsheet

# Chemicals - Definition

def create_chemicals():
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
    Methyl_oct_3 = tmo.Chemical('Methyl_oct_3', search_ID ='3-Methyloctane')
    Methyl_oct_4 = tmo.Chemical('Methyl_oct_4', search_ID = '4-Methyloctane')
    Dimethyl_heptane_3_3 = tmo.Chemical('Dimethyl_heptane_3_3', search_ID ='3,3-Dimethylheptane')
    Ethylheptane_3 = tmo.Chemical('Ethylheptane_3',search_ID ='3-Ethylheptane')
    Ethylheptane_4 = tmo.Chemical('Ethylheptane_4', search_ID ='4-Ethylheptane') 
    bicyclo_octane = tmo.Chemical('bicyclo_octane',search_ID = '6221-55-2')
    Linoleic_acid = tmo.Chemical('Linoleic_acid',search_ID = '60-33-3')
    Palmitic_acid = tmo.Chemical('Palmitic_acid',search_ID = '57-10-3')
    Stearic_acid =  tmo.Chemical('Stearic Acid',search_ID = '57-11-4')
    Other_FAs = tmo.Chemical('Other_FAs', search_ID = '67701-02-4')
    EPS = tmo.Chemical('EPS', search_ID = '24560-98-3')   
   
    ozo_chemicals = tmo.Chemicals(
        ['Water','Hydrogen_peroxide','Oleic_acid',
         'Nonanal','Nonanoic_acid','Azelaic_acid',
         'Hexane',
         'Ethyl_acetate','Octane',
         'pentene',
         'Methylcyclohexane',
         'Cyclopentane',
         '540-84-1',
         'heptane',
         'pentane',
         'Sulphuric_acid',
         'Nonane',
         'Decane',
         'Undecane',
         'Cyclohexane',
         'Cycloheptane',
         'Benzene',
         'toluene',
          bicyclo_octane,
          Methyl_oct_3,
          Methyl_oct_4,
          Dimethyl_heptane_3_3,
          Ethylheptane_3,
          Ethylheptane_4,
          Palmitic_acid,
          Stearic_acid,
          Other_FAs,
          EPS
          ])      

    Catalyst = create_new_chemical(
            'Phosphotungstic_acid',
            # formula="H3PW12O40",
            # MW=2880.2,
            #CAS = '7783-03-1',
            CAS='1343-93-7',
            phase = 'l',
            # rho = 2850 
        )
## Impurities in oleic acid

    #https://www.fishersci.com/shop/products/phosphotungstic-acid-hydrate-thermo-scientific/AA4011614
    #rho = 2.852g/cm3
    
    #####I think the chemmical is actually this###
    # https://www.sigmaaldrich.com/US/en/product/aldrich/455970
    
    #Liquid chemical not available in the database
    #TODO.xxx Add ref for the below
    Oxononanoic_acid = create_new_chemical(
       'Oxononanoic_acid',
        search_ID = '2553-17-5',
        phase='l',
        Tb = 304.8 + 25.0 + 273.15,
        Hf=-579480,
        formula = 'C9H16O3',
        MW = 172.22,
        CAS = '2553-17-5'
    )
    
#TODO.xxx use this catalyst later for heterogeneous system
   
    # bea_zeolite = create_new_chemical(
    #               'bea_zeolite',
    #               search_ID = '1318-02-1',
    #               phase = 's',
    #               Tb = 10000+273,
    #               CAS = '1318-02-1')
    # aluminium_oxide = tmo.Chemical('1344-28-1')
    # bea_zeolite.copy_models_from(aluminium_oxide, ['V', 'Hvap', 'Psat', 'mu', 'kappa'])

    for chemical in ozo_chemicals: chemical.default()
    
    Dortmund = ozo_chemicals.Oxononanoic_acid.Dortmund
    Dortmund.set_group_counts_by_name(dict(COOH=1, CH2=7, CHO=1))
    ozo_chemicals.compile()
    ozo_chemicals.set_synonym('Phosphotungstic_acid', 'fresh_Cat')
    
    for i in ozo_chemicals:
        if not i.locked_state: i.V.g.method_P = 'IDEAL'
    
    return ozo_chemicals
    






