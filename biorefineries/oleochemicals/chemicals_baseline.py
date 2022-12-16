"""
Created on Sat Aug 20 21:47:53 2022

@author: Lavanya_Kudli
"""
import biosteam as bst
import biosteam.units
import thermosteam as tmo
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
Tungsten = tmo.Chemical('Tungsten')
Cobalt = tmo.Chemical('Cobalt')
Hydrogen = tmo.Chemical('Hydrogen',phase = 'l')

chems = tmo.Chemicals([
    #Dihydroxylation chemicals
    tmo.Chemical('Hydrogen_peroxide', phase='l'),
    tmo.Chemical('Water'),
    # #Chemical for acid degumming 
    tmo.Chemical('Citric_acid'),
    
#TODO: look into phase of the below
    tmo.Chemical('MDHSA', 
                 search_ID = '1115-01-1', 
                 phase = 'l'),
#Ref for dihydroxy_palmitic_acid: https://www.chemsrc.com/en/cas/29242-09-9_803992.html    
    tmo.Chemical('Dihydroxy_palmitic_acid',
                 search_ID = '29242-09-9',
                 search_db = False,
                 Tb = 458 + 273.15,
                 formula = 'C16H32O4',
                 phase = 'l'
                 ),
#Ref for tetrahydroxy_octadecanoic_acid: https://www.chemsrc.com/en/cas/541-82-2_148112.html
    tmo.Chemical('Tetrahydroxy_octadecanoic_acid',
                 search_ID = '541-82-2', 
                 search_db = False, 
                 Tb = 583.1+273.15,
                 formula = 'C18H36O6',
                 phase = 'l'),
    
# Products of oxidative_cleavage
    tmo.Chemical('Monomethyl_azelate'),
    tmo.Chemical('Suberic_acid'),
    tmo.Chemical('Caprylic_acid'),
    tmo.Chemical('Hexanoic_acid'),
    tmo.Chemical('Heptanoic_acid'),
    tmo.Chemical('Hexanoic_acid'),
    tmo.Chemical('Malonic_acid'),
    tmo.Chemical('Pelargonic_acid'),
    
# Products of hydrolysis
    tmo.Chemical('Palmitic_acid'),
    tmo.Chemical('Stearic_acid'),
    tmo.Chemical('Oleic_acid'),
    tmo.Chemical('Linoleic_acid', search_ID = '60-33-3'),
    tmo.Chemical('Palmitoleic_acid', search_ID = '373-49-9'),
#TODO:should I set the phase or not?    
    tmo.Chemical('Azelaic_acid', phase = 's'),

# Oxidants used and other gaseous products
    tmo.Chemical('Nitrogen'),
    tmo.Chemical('Oxygen'),
    tmo.Chemical('Carbon_dioxide'),
    tmo.Chemical('Methanol'),
    tmo.Chemical('Glycerol', rho = 1261.3, phase = 'l'),
    tmo.Chemical('Sodium_methoxide',formula ='NaOCH3'),
   
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
    tmo.Chemical('Methyl_oleate', phase = 'l'),
    tmo.Chemical('Methyl_palmitate', phase = 'l'),
    tmo.Chemical('Methyl_stearate', phase = 'l'),
    tmo.Chemical('Methyl_linoleate', phase = 'l'), 
    tmo.Chemical('Methyl_palmitoleate',search_ID ='1120-25-8',
                                       phase = 'l'),
     
## Catalyst data    
## Tungstic_acid boiling point: https://en.wikipedia.org/wiki/Tungstic_acid
## TODO: Ask Yoel about the phase of tungstic acid
    tmo.Chemical('Tungstic_acid', 
                  Tb = 1746, 
                  phase = 's', 
                  Hvap=Tungsten.Hvap,
                  Psat=Tungsten.Psat,
                  default = True),
    tmo.Chemical('Tungstate_ion',
                 search_db = False,
                 CAS = '12737-86-9',
                 Hvap=Tungsten.Hvap,
                 Psat=Tungsten.Psat,
                 formula = 'O4W-2',
                 MW =  tmo.Chemical('Tungstic_acid').MW,
                 phase = 'l',
                 Tb = Tungsten.Tb
                 ),
    tmo.Chemical('Hydrogen_ion',
                 Tb = Hydrogen.Tb,
                 phase = 'l',
                 ),    
## TODO: GWP data for cobalt acetate missing, figure it out 
## cobalt_acetate_tetrahydrate Tb: 
    Cobalt_chloride,
    tmo.Chemical('Cobalt_acetate_tetrahydrate',
                 search_ID = '6147-53-1',
                 Tb = 117.1+273.15, 
                 phase = 's',
                 Hvap=Cobalt_chloride.Hvap,
                 Psat=Cobalt_chloride.Psat,
                 default=True),
    tmo.Chemical('Acetate_ion', phase = 'l', search_ID = '71-50-1'),
## TODO: defaulting the below to water for now, maybe look for a better assumption
    tmo.Chemical('Cobalt(2+)',
                 Hvap=Cobalt.Hvap,
                 Psat=Cobalt.Psat,
                 phase = 'l',
                 Tb=Cobalt.Tb,
                 default = True),
    
    ###Modelling amberlyte catalyst like a solid catalyst
    ##Using sunfonated_polystyrene boiling point
    ##https://www.chemsrc.com/en/cas/39389-20-3_843683.html#:~:text=amberlyst%28r%29%2015%20CAS%20Number%3A%2039389-20-3%3A%20Molecular%20Weight%3A%20314.39900%3A,Point%3A%20266.3%C2%BAC%3A%20Symbol%3A%20GHS07%3A%20Signal%20Word%3A%20Warning%20%C3%97
    
    ##Resin for hydrolysis
    tmo.Chemical('polystyrene_based_catalyst',
                 search_db=False,
                 Tb = 516.7+273.15,
                 phase = 's',
                 default=True),
    ##For catalyst recovery chemicals required
    tmo.Chemical('Calcium_hydroxide',
                 phase = 's',
                 default = True),
    tmo.Chemical('Calcium_chloride'),    
    tmo.Chemical.blank('Calcium_tungstate',
                       CAS = '7790-75-2',
                       MW = 287.92,
                       Tb = Tungsten.Tb,   
                       phase = 's',
                       formula = 'CaWO4'),
    
    tmo.Chemical('Calcium_acetate',                      
                  phase = 'l'),
    
    tmo.Chemical('Cobalt_hydroxide',                      
                  phase = 's'),
    
    
    ##Below lacks Hvap etc models
    # Sulfonated_polystyrene = chemical_database('Sulfonated_polystyrene',
    #                                            search_ID = '98-70-4',
    #                                            )
    ##Hence using polystyrene
    tmo.Chemical('Polystyrene', phase = 's'),
    
##For imported lipidcane module compatibility
    tmo.Chemical('MonoOlein',search_ID = '111-03-5'),
    # tmo.Chemical('DiOlein',search_ID = 'PubChem = 6505653'),
    tmo.Chemical('Dipalmitin'),
    tmo.Chemical('Monopalmitin'),
    
    #Below TAG's not a part of the database    
    tmo.Chemical('OOL',
                 search_ID = '104485-08-7',
                 phase ='l',
                 search_db = False),
    
    tmo.Chemical('LLO', 
                 search_ID = '2190-22-9',
                 phase = 'l',
                 search_db = False),
#Natural gas for heating purposes 
    tmo.Chemical('Natural_gas',
                 search_ID = 'CH4'),
#Solvent for countercurrent extraction of azelaic acid
    tmo.Chemical('Octane'),
    tmo.Chemical('Cycloheptane'),
    tmo.Chemical('Bicyclo_octane',search_ID = '6221-55-2'),
    tmo.Chemical('Toluene')
    
                 ])

##Modelling the properties of resin used for hydrolysis based on polystyrene
chems.polystyrene_based_catalyst.copy_models_from(chems.Polystyrene,
                                        ['Hvap','Psat','sigma', 
                                         'epsilon', 'kappa', 'V',
                                         'Cn', 'mu'])
LiquidMethanol = chems['Methanol'].at_state(phase='l', copy=True)
chems['Sodium_methoxide'].copy_models_from (LiquidMethanol,
                                            ['V', 'sigma',
                                             'kappa', 'Cn',
                                             'Psat'])
## The oxidative cleavage catalyst properties are based on cobalt chloride
chems['Cobalt_acetate_tetrahydrate'].copy_models_from(Cobalt_chloride,
                                            ['sigma',
                                             'kappa', 'Cn',
                                             ])
## The density of cobalt acetate tetrahydrate is available in the literature, rho is 1730 Kg/m3
## https://scifinder-n.cas.org/searchDetail/substance/634d77a73c1f076117e95f61/substanceDetails
## Cobalt acetate dissolution reaction data
V_of_cobalt_acetate_tetrahydrate = fn.rho_to_V(1730,chems['Cobalt_acetate_tetrahydrate'].MW)
chems['Cobalt_acetate_tetrahydrate'].V.add_model(V_of_cobalt_acetate_tetrahydrate,
                                                 top_priority = True)
## Tungstate ion default to properties
chems['Tungstate_ion'].copy_models_from(chems['Tungstic_acid'],
                                        ['Hvap',
                                         'Psat',
                                         'Cn',
                                         'V',
                                         'mu'
                                         ])
## Hydrogen ion defaults to properties of hydrogen
chems['Hydrogen_ion'].copy_models_from(Hydrogen,
                                       ['Hvap',
                                        'Psat',
                                        'mu',
                                        'V',
                                        'Cn'])
## Products of the precipitation reaction
## https://scifinder-n.cas.org/searchDetail/substance/634db2343c1f076117eb77ce/substanceDetails
chems['Calcium_tungstate'].copy_models_from(Tungsten,
                                            ['Hvap',
                                             'Psat'])
chems['Calcium_tungstate'].V.add_model(fn.rho_to_V(5800,
                                                   chems['Calcium_tungstate'].MW),
                                                   top_priority = True)
chems['Acetate_ion'].copy_models_from(tmo.Chemical('Acetate'),
                                      ['mu',
                                       'Psat',
                                       'Hvap',
                                       'sigma',
                                       ])



## Changing this for Azelaic acid as it probably didnt work during a distillation process
chems['Azelaic_acid'].Cn.method = 'LASTOVKA_S'

## Modelling properties of dihydroxylated compounds as MDHSA
chems['Dihydroxy_palmitic_acid'].copy_models_from(chems['MDHSA'],
                                                  ['Hvap',
                                                   'Psat',
                                                   'Cn',
                                                   'V',
                                                   'mu'
                                                   ])
# chems['Dihydroxy_palmitic_acid'].Cn.method = 'LASTOVKA_S'

                                                  
chems['Tetrahydroxy_octadecanoic_acid'].copy_models_from(chems['MDHSA'],
                                                         ['Hvap',
                                                          'Psat',
                                                          'Cn',
                                                          'V',
                                                          'mu'
                                                          ])
# chems['Tetrahydroxy_octadecanoic_acid'].Cn.method = 'LASTOVKA_S'


#TODO.xxx check if Psat from liquid methanol is a good idea
def create_new_chemical(ID, phase='s', **constants):
        solid = tmo.Chemical.blank(ID, phase=phase, **constants)
        chems.append(solid)
        return solid
    
NaOH = create_new_chemical('NaOH', formula='NaOH')
HCl = create_new_chemical('HCl', formula='HCl')

# Solubles don't occupy much volume
for chemical in (HCl, NaOH,):
        V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
 

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

#composition of VM_Naphtha based on https://www.cdc.gov/niosh/npg/npgd0664.html#:~:text=None%20reported%20%5BNote%3A%20VM%26P%20Naphtha%20is%20a%20refined,Exposure%20Routes%20inhalation%2C%20ingestion%2C%20skin%20and%2For%20eye%20contact
chems.define_group('VM_Naphtha',['Octane',
                                 'Cycloheptane',
                                 'Bicyclo_octane',
                                 'Toluene'], 
                                 composition = [0.55,
                                                0.30,
                                                0.02,
                                                0.12])

chems.define_group('Air', ['Oxygen', 'Nitrogen'],composition=[0.21,0.79])
chems.set_synonym('Water', 'H2O')
chems.set_synonym('Carbon_dioxide','CO2')
chems.set_synonym('Phosphatidylinositol','PL')
chems.set_synonym('MonoOlein', 'MAG')
chems.set_synonym('Dipalmitin', 'DAG')
chems.set_synonym('Cobalt(2+)','Cobalt_ion')
chems.set_synonym('Hydrogen_ion', 'H+')
chems.set_synonym('Pelargonic_acid','Nonanoic_acid')
bst.settings.set_thermo(chems)
# chems.show()