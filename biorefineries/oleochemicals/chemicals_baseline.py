"""
Created on Sat Aug 20 21:47:53 2022

@author: Lavanya_Kudli
"""
import biosteam as bst
import biosteam.units
import thermosteam as tmo
from biorefineries import lipidcane

from biorefineries import sugarcane as sc
from thermosteam import functional as fn
from chemicals import atoms_to_Hill
import thermosteam as tmo


##Chemicals that are already in the data base
# Biodiesel = chemical_database('Methyl_oleate', phase = 'l')
# # available methods are {'ROWLINSON_POLING', 'DADGOSTAR_SHAW', 'ROWLINSON_BONDI'}
# Biodiesel.Cn.method = 'ROWLINSON_BONDI'

##All the TAGs are based on Ruiz-Gutiérrez et. al (1998)
# chems is the object containing all chemicals used in this biorefinery
Cobalt_chloride = tmo.Chemical('Cobalt_chloride',
                               search_ID = '7646-79-9',
                               phase = 's')
tungsten = tmo.Chemical('tungsten')
chems = tmo.Chemicals([
    tmo.Chemical('Hydrogen_peroxide', phase='l'),
    tmo.Chemical('Water'),
    #look into phase of the below
    tmo.Chemical('MDHSA', search_ID = '1115-01-1'),
    tmo.Chemical('Pelargonic_acid'),
    tmo.Chemical('Azelaic_acid'),
    tmo.Chemical('Monomethyl_azelate'),
    tmo.Chemical('Suberic_acid'),
    tmo.Chemical('Caprylic_acid'),
    tmo.Chemical('Nitrogen'),
    tmo.Chemical('Oxygen'),
    tmo.Chemical('Methanol'),
    tmo.Chemical('Glycerol'),
    tmo.Chemical('Sodium_methoxide',formula ='NaOCH3'),
    tmo.Chemical('HCl'),
    tmo.Chemical('NaOH'),
    
###All the chemicals related to TAGs in HOSO
    tmo.Chemical('OOO', search_ID = '122-32-7'),
    tmo.Chemical('LLL', search_ID = '537-40-6'),
    
    tmo.Chemical('OOL',
                 search_ID = '104485-08-7',
                 phase ='l',
                 search_db = False),
    
    tmo.Chemical('LLO', 
                 search_ID = '2190-22-9',
                 phase = 'l',
                 search_db = False),
    
    tmo.Chemical('SOO',
                 search_ID = '79517-06-9',
                 phase ='l',
                 search_db = False),
    
    tmo.Chemical('PLO', 
                 search_ID = '2680-59-3',
                 phase = 'l',
                 search_db = False),
      
    tmo.Chemical('PoOO',
                 search_ID = 'PubChem= 9544083',
                 phase ='l',
                 search_db = False),

    tmo.Chemical('POO',
                 search_ID = '',
                 phase ='l',
                 search_db = False),
    
    tmo.Chemical('POS', 
                 search_ID = 'PubChem = 129723993',
                 phase = 'l',
                 search_db = False),
    
    tmo.Chemical('POP',
                 search_ID = '2190-25-2',
                 phase ='l',
                 search_db = False),
    
    tmo.Chemical('PLS',
                 search_ID = 'PubChem = 102332193',
                 phase ='l',
                 search_db = False),
    
##chemical for representing phospholipids, taken directly from lipidcane biorefinery
    tmo.Chemical('Phosphatidylinositol', 
                 formula='C47H83O13P',
                     search_db=False, 
                     CAS='383907-36-6',
                     default=True,
                     Hf=-1.779e6, # Assume same as TAG on a weight basis
                     phase='l'),
  
##All the chemicals that go in the Biodiesel

    tmo.Chemical('Methyl_oleate'),
    tmo.Chemical('Methyl_palmitate'),
    tmo.Chemical('Methyl_stearate'),
    tmo.Chemical('Methyl_linoleate'), 
    tmo.Chemical('Methyl_palmitoleate',search_ID ='1120-25-8'),
    
## Extra chems for lipidcane compatibility
    tmo.Chemical('MonoOlein', search_ID = '111-03-5'),
    tmo.Chemical('DiOlein', search_ID = '2465-32-9'),
   
    
    ### Chemicals that were missing some properties
    ### TODO.xxx check if this is a good assumption with Yoel
    
    #Tungstic_acid boiling point: https://en.wikipedia.org/wiki/Tungstic_acid
    tmo.Chemical('tungstic_acid', 
                            Tb = 1746, 
                            phase = 's', 
                            Hvap=tungsten.Hvap,
                            Psat=tungsten.Psat,
                            default = True),
    
    ###Using cobalt chloride instead of acetate as allowed by the patent
    ###cobalt acetate has a few missing properties, further GWP data not available
    
    # cobalt_acetate_Tb: https://www.chemsrc.com/en/cas/71-48-7_34110.html
    Cobalt_chloride,
    tmo.Chemical('Cobalt_acetate',Tb = 117.1+273.15, 
                phase = 's',
                Hvap=Cobalt_chloride.Hvap,
                Psat=Cobalt_chloride.Psat,
                default=True),
    
    ###Modelling amberlyte catalyst like a solid catalyst
    ##Using sunfonated_polystyrene boiling point
    ##https://www.chemsrc.com/en/cas/39389-20-3_843683.html#:~:text=amberlyst%28r%29%2015%20CAS%20Number%3A%2039389-20-3%3A%20Molecular%20Weight%3A%20314.39900%3A,Point%3A%20266.3%C2%BAC%3A%20Symbol%3A%20GHS07%3A%20Signal%20Word%3A%20Warning%20%C3%97
    
    ##Chemicals not in the database
    tmo.Chemical('polystyrene_based_catalyst',
                 search_db=False,
                 Tb = 516.7+273.15,
                 phase = 's',
                 default=True),
    
    tmo.Chemical('OOL',
                 search_ID = '104485-08-7',
                 phase ='l',
                 search_db = False),
    
    tmo.Chemical('LLO', 
                 search_ID = '2190-22-9',
                 phase = 'l',
                 search_db = False),
    ##Below lacks Hvap etc models
    # Sulfonated_polystyrene = chemical_database('Sulfonated_polystyrene',
    #                                            search_ID = '98-70-4',
    #                                            )
    ##Hence using polystyrene
    tmo.Chemical('Polystyrene'),
    ##Hence using polystyrene
    tmo.Chemical('Tripalmitin')
  ])
chems.polystyrene_based_catalyst.copy_models_from(chems.Polystyrene,
                                        ['Hvap','Psat','sigma', 
                                         'epsilon', 'kappa', 'V',
                                         'Cn', 'mu'])
PPP = tmo.Chemical('Tripalmitin')
TAGs = chems['OOO','LLL','OOL','LLO','SOO','PLO',
             'PoOO','POO','POS','POP','PLS']
for i in TAGs:
    i.copy_models_from(PPP, ['mu'])
    
for chemical in chems: chemical.default()

chems.compile()
chems.define_group('TAG', ('OOO','LLL','OOL',
                           'LLO','SOO','PLO',
                           'PoOO','POO','POS',
                           'POP','PLS'))

chems.define_group('Lipid', ('OOO','LLL','OOL',
                           'LLO','SOO','PLO',
                           'PoOO','POO','POS',
                           'POP','PLS'))
chems.define_group('Oil', ('OOO','LLL','OOL',
                           'LLO','SOO','PLO',
                           'PoOO','POO','POS',
                           'POP','PLS'))

chems.define_group('Biodiesel', ('Methyl_oleate',
                                 'Methyl_palmitate',
                                 'Methyl_stearate',
                                 'Methyl_linoleate',
                                 'Methyl_palmitoleate'))

chems.set_synonym('Water', 'H2O')
chems.set_synonym('Phosphatidylinositol','PL')
chems.set_synonym('MonoOlein', 'MAG')
chems.set_synonym('DiOlein', 'DAG')
bst.settings.set_thermo(chems)
chems.show()