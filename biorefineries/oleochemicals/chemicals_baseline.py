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
from thermosteam import Chemicals
from thermo import TDependentProperty
from biorefineries.oleochemicals import TAG_properties
from TAG_properties import *
#chems is a list of all the chemicals used in this biorefinery

chems = tmo.Chemicals([
#Biodiesel preparation section chemicals
        tmo.Chemical('Sodium_methoxide',formula ='NaOCH3',phase = 'l',default = True), #Catalyst
        tmo.Chemical('Methanol'), 
        tmo.Chemical('Citric_acid'),#Chemical for acid degumming 
        tmo.Chemical('Phosphatidylinositol', ##chemical for representing phospholipids
                      formula='C47H83O13P',
                      search_db=False, 
                      CAS='383907-36-6',
                      MW = 	886.56,
                      default=True,
                      Hf=-1.779e6,#taken from cane biorefinery
                      phase = 'l',#taken from cane biorefinery
                      ),
        tmo.Chemical('Glycerol'),#By-product
#TAGs of high oleic oils
        tmo.Chemical('OOO', search_ID = '122-32-7',Hf = 1000*(-76.494*18 - 815.18),phase_ref = 'l' ),#from cane>chemicals
        tmo.Chemical('LnLnLn', search_ID = '537-40-6',Hf = 1000*(-316.67*3 - 2466.4)),#This is Trilinolein
#Below TAG's not a part of the database    
         tmo.Chemical('LLL',
                      search_ID = '7049-66-3',#Based on reference search in CAS Finder
                      search_db = False,
                      formula = 'C54H96O6',#Based on reference search in CAS Finder
                      phase = 'l',
                      Hf = (((2/3)*(-76.494*18 - 815.18))+((1/3)*(-316.67*2 - 2466.4)))*1000 #TODO: make a dictionary for these values
                      ),  
        tmo.Chemical('OOL',
                     search_ID = '28880-78-6',#Based on reference search in CAS Finder
                     search_db = False,
                     formula = 'C57H102O6',#Based on reference search in CAS Finder
                     phase = 'l',
                     Hf = (((2/3)*(-76.494*18 - 815.18))+((1/3)*(-316.67*2 - 2466.4)))*1000 
                     ),           
        tmo.Chemical('LLO', 
                         search_ID = '28409-91-8',#Based on reference search in CAS Finder
                         search_db = False,
                         formula = 'C57H100O6',#Based on reference search in CAS Finder
                         phase = 'l',
                         Hf = (((2/3)*(-316.67*2 - 2466.4))+((1/3)*(-76.494*18 - 815.18)))*1000
                         ),  
        tmo.Chemical('SOO',
                         search_ID = '29590-02-1',#Based on reference search in CAS Finder
                         search_db = False,
                         formula = 'C57H106O6',#Based on reference search in CAS Finder
                         phase = 'l',
                         Hf = 1000*(((1/3)*( -59.571*18 - 1358.7))+((2/3)*(-76.494*18 - 815.18)))
                         ),      
        tmo.Chemical('PLO', 
                         search_ID = '26836-35-1',#Based on reference search in CAS Finder
                         search_db = False,
                         formula = 'C55H100O6',#Based on reference search in CAS Finder
                         phase = 'l',
                         Hf = 1000*(((1/3)*(-76.494*18 - 815.18)) + ((1/3)*(-316.67*2 - 2466.4)) + ((1/3)*(-59.571*16 - 1358.7)))
                         ),         
        
        tmo.Chemical('PoOO',
                     search_ID = '38703-17-2',#Based on reference search in CAS Finder
                     search_db = False,
                     formula = 'C55H100O6',#Based on reference search in CAS Finder
                     phase = 'l',
                     Hf = 1000*(((1/3)*(-76.494*16 - 815.18))+ ((2/3)*(-76.494*18 - 815.18)))),
        tmo.Chemical('POO',
                 search_ID = '27071-84-7',#Based on reference search in CAS Finder
                 search_db = False,
                 formula = 'C55H102O6',#Based on reference search in CAS Finder
                 phase = 'l',
                 Hf = 1000*(((2/3)*(-76.494*18 - 815.18)) + ((1/3)*(-59.571*16 - 1358.7)))
                 ),
    
        tmo.Chemical('POS', 
                 search_ID = '26836-31-7',#Based on reference search in CAS Finder
                 search_db = False,
                 formula = 'C55H104O6', #Based on reference search in CAS Finder
                 phase = 'l',
                 Hf = 1000*(((1/3)*(-76.494*18 - 815.18)) + ((1/3)*(-59.571*16 - 1358.7)) + ((1/3)*( -59.571*18 - 1358.7)))
                 ),
    
        tmo.Chemical('POP',
                 search_ID = '28409-94-1',#Based on reference search in CAS Finder
                 search_db = False,
                 formula = 'C53H100O6',#Based on reference search in CAS Finder
                 phase = 'l',
                 Hf = (((2/3)*(-59.571*16 - 1358.7)) + ((1/3)*(-76.494*18 - 815.18)))*1000
                 ),
    
        tmo.Chemical('PLS',
                 search_ID = '26836-32-8',#Based on reference search in CAS Finder 
                 search_db = False,
                 formula = 'C55H102O6',#Based on reference search in CAS Finder
                 phase = 'l',
                 Hf = 1000*(((1/3)*(-59.571*16 - 1358.7)) +((1/3)*(-59.571*18 - 1358.7))+ ((1/3)*( -316.67*2 - 2466.4)))
                 ),
        tmo.Chemical('PPP',
                     search_ID = 'Tripalmitin',
                     Hf = 1000*(-59.571*16 - 1358.7)
                     ),
        tmo.Chemical('SSS',
                     search_ID = 'Tristearin',
                     Hf = 1000*(-59.571*18 - 1358.7),
                     phase = 'l'),
        
##Chemicals part of Biodiesel
        tmo.Chemical('Methyl_oleate'),
        tmo.Chemical('Methyl_palmitate'),
        tmo.Chemical('Methyl_stearate'),
        tmo.Chemical('Methyl_linoleate'), 
        tmo.Chemical('Methyl_linolenate',search_ID = '301-00-8'),
        tmo.Chemical('Methyl_palmitoleate',search_ID ='1120-25-8'),
     
#Dihydroxylation chemicals
        tmo.Chemical('Hydrogen_peroxide',phase ='l'),
        tmo.Chemical('Water',default = True),
#Dihydroxylation products
#Product 1       
        tmo.Chemical('MDHSA', 
                     search_ID = '1115-01-1', 
                     Hf = -947.7*1000, #Wikepedia value for https://en.wikipedia.org/wiki/Stearic_acid
                    phase ='l'
                    ),
        #Diketo deriviative: 7108-68-1
        #methyl 9,10-dioxooctadecanoate
        #Monoketo derivative: 7297-29-2
        
#Product 2
        tmo.Chemical('Methyl_dihydroxy_palmitate',
                 search_ID = '908329-09-9',#Reaxys CAS number for methyl 9,10-dihydroxy-palmitate
                 search_db = False,
                 Tb = 458 + 273.15,#Based on 9,10 dihydroxy palmitic acid CAS finder (29242-09-9): 458.0±25.0 °C,Press: 760 Torr
                 formula = 'C17H34O4',
                 Hf = -892*1000, #Based on palmitic acid, Ref: https://en.wikipedia.org/wiki/Palmitic_acid
                 phase = 'l'
                 ),
#Product 3
#Below properties based on DWSIM
        tmo.Chemical('Tetrahydroxy_octadecanoate',
                 search_ID = '61117-79-1', #Reaxys
                 search_db = False, 
                 MW = 362.51,
                 Tb = 1077.57,#DWSIM                
                 Hf = -3614.63*1000, #DWSIM
                 phase = 'l'
                 ),
# Product 4 
#Below properties based on DWSIM   
        tmo.Chemical.blank('Hexahydroxy_octadecanoate',#Structure obtained from Reaxys
                           MW = 394.51,
                           Tb = 1261.05,
                           formula = 'C19H38O8',
                           Hf = -4119.97*1000, #Assumed to be same as Linolenic acid
                           phase = 'l'
                           ),
#Possible monoester of MDHSA with Pelargonic acid        
        tmo.Chemical.blank('Monoester_MDHSA_PA',#Structure obtained from Reaxys                     
                     CAS = '950910-15-3',
                     MW = 470.73,
                     Tb = 1084.12,#DWSIM
                     formula = 'C28H54O5',#CAS generated compound
                     Hf = -2705.64*1000, #DWSIM
                     phase = 'l'
                     ),
#Possible diester of MDHSA with Monomethyl azelate acid        
        tmo.Chemical.blank('Diester_MDHSA_PA',#Structure obtained from Reaxys
                     
                     MW = 610.76,
                     Tb = 1274.15,
                     formula = 'C37H70O6',
                     Hf = -2540.21*1000,
                     phase = 'l'),      
        tmo.Chemical.blank('Monoester_MDHSA_MMA',#Structure obtained from Reaxys
                     
                     MW = 514.74,
                     Tb = 1183.29,
                     formula = 'C29H54O7' ,
                     Hf =-2989.99*1000,
                     phase = 'l'
                     ),
        tmo.Chemical.blank('Diester_MDHSA_MMA',#Structure obtained from Reaxys
                     
                     MW = 698.98,
                     Tb =1472.49,
                     formula = 'C39H70O10' ,
                     Hf =-2979.85*1000,
                     phase = 'l'),
        
    
# Products of oxidative_cleavage
        tmo.Chemical('Monomethyl_azelate'),
        tmo.Chemical('Monomethyl_suberate',search_ID = '3946-32-5'),
        tmo.Chemical('Caprylic_acid'),
        tmo.Chemical('Hexanoic_acid'),
        tmo.Chemical('Heptanoic_acid'),
        tmo.Chemical('Hexanoic_acid'),
        tmo.Chemical('Malonic_acid'),
        tmo.Chemical('Pelargonic_acid'),
        tmo.Chemical('Propanoic_acid'),
        tmo.Chemical('Methyl_caprylate',search_ID = '111-11-5'),
        tmo.Chemical('Nonanal'),
        tmo.Chemical('Methyl_oxo_nonanoicacid',
                     search_ID = '1931-63-1'
                     ),
        
# Products of hydrolysis
        tmo.Chemical('Palmitic_acid'),
        tmo.Chemical('Stearic_acid'),
        tmo.Chemical('Oleic_acid',phase = 'l'),
        tmo.Chemical('Sodium_oleate',Tb = 633 ,Hf = -764800.0,phase ='l'),
        tmo.Chemical('Linoleic_acid', search_ID = '60-33-3'),
        tmo.Chemical('Linolenic_acid'),
        tmo.Chemical('Palmitoleic_acid', search_ID = '373-49-9'),
        tmo.Chemical('Azelaic_acid'),

# Oxidants used and other gaseous products
#Liquid phase reaction?   
        tmo.Chemical('Nitrogen',phase = 'g'),
        tmo.Chemical('Oxygen',phase = 'g'),
        tmo.Chemical('Carbon_dioxide',phase = 'g'),
        
#Solvent for countercurrent extraction of monocarboxylic acids
        tmo.Chemical('Octane'),
        tmo.Chemical('Cycloheptane'),
        tmo.Chemical('Bicyclo_octane',search_ID = '6221-55-2'),
        tmo.Chemical('Toluene'), 
        tmo.Chemical('Heptane'),
     
## Catalysts used in the process   
## Tungstic_acid boiling point: https://en.wikipedia.org/wiki/Tungstic_acid
       tmo.Chemical('Tungsten',default = True),
       tmo.Chemical('Tungstic_acid', 
                  Tb = 1746,
                  default = True,
                  formula = 'H2WO4',
                  phase = 'l'),
       tmo.Chemical('Sodium_tungstate',
                    phase = 'l',
                    Tb = 1746, #same as assumed above
                    ),
       tmo.Chemical('Sodium_chloride'),
        tmo.Chemical('Tungstate_ion',
                 search_db = False,
                 CAS = '12737-86-9',
                 formula = 'O4W-2',
                 MW =  tmo.Chemical('Tungstic_acid').MW,#TODO: check this assumption
                 phase = 'l',
                 default = True,
                 Tb = 5828.15 #Based on tungsten's BP
                 ),        
        tmo.Chemical('Cobalt',default = True),
        tmo.Chemical('Hydrogen',default = True),
        tmo.Chemical('Hydrogen_ion',
                 Tb = 20.271,#Based on hydrogen's BP
                 phase = 'l',
                 formula = 'H+',
                 default = True
                 ),    
        tmo.Chemical('Cobalt_chloride',
                     search_ID = '7646-79-9', 
                     phase = 'l'
                     ),
        tmo.Chemical('Cobalt_acetate_tetrahydrate',
                      phase = 'l',#TODO: explore phase
                      search_ID = '6147-53-1',
                      Tb = 117.1+273.15,
                      Hf = -2.7151e+05 #Assumed to be the same as cobalt chloride for now
                      ),
        tmo.Chemical('Acetate',phase = 'l', default = True),
        tmo.Chemical('Acetate_ion', phase = 'l', search_ID = '71-50-1',default = True),
        tmo.Chemical('Cobalt_ion',search_ID ='Cobalt(2+)',  phase = 'l', 
                      Tb=3200.15,#Based on cobalt's BP
                      default = True),
##Chemicals required for catalyst recovery
        tmo.Chemical('Calcium_hydroxide',default = True,phase = 'l'),#assumed this phase because it is soluble in the reaction media it is used in
        tmo.Chemical('Calcium_chloride',phase = 'l'),    
        tmo.Chemical.blank('Calcium_tungstate',
                       CAS = '7790-75-2',
                       MW = 287.92,
                       Tb = 5828.15, #Based on tungsten's BP
                       formula = 'CaWO4',
                       phase = 'l',
                       Hf = -1002.82*1000),#Ref: Value for calcium hydroxide https://www.chemeurope.com/en/encyclopedia/Standard_enthalpy_change_of_formation_%28data_table%29.html
        tmo.Chemical('Calcium_acetate',phase = 'l', default = True),
        tmo.Chemical('Cobalt_hydroxide',phase = 'l', default = True),
        tmo.Chemical('HCl2',phase = 'l',search_ID = 'HCl'),
        tmo.Chemical('Sodium_hydroxide_liquid',phase = 'l', search_ID ='NaOH'),
        tmo.Chemical('Sodium_acetate', phase = 'l'),
##Resin for hydrolysis
##Sulfonated_polystyrene with a search_ID = '98-70-4', lacks Hvap data, hence polystyrene's properties were used
##Boiling point based on amberlyte
# Ref:##https://www.chemsrc.com/en/cas/39389-20-3_843683.html#:~:text=amberlyst%28r%29%2015%20CAS%20Number%3A%2039389-20-3%3A%20Molecular%20Weight%3A%20314.39900%3A,Point%3A%20266.3%C2%BAC%3A%20Symbol%3A%20GHS07%3A%20Signal%20Word%3A%20Warning%20%C3%97
        tmo.Chemical('polystyrene_based_catalyst',
                     search_ID='Polystyrene',
                     Tb = 516.7+273.15,
                     phase = 's',
                     default=True),   
##For imported lipidcane module compatibility
        tmo.Chemical('MonoOlein',search_ID = '111-03-5'),
        tmo.Chemical('Dipalmitin'),
        tmo.Chemical('Monopalmitin'),
#Natural gas for heating purposes 
        tmo.Chemical('Natural_gas',
                     search_ID = 'CH4'), 
#Chemicals required by the boilerturbogenerator        
        tmo.Chemical('Ash',
                     MW = 1,
                     phase = 's',
                     search_db=False,
                     default = True,
                     Hf = 0),
        tmo.Chemical('P4O10',
                     phase = 's',
                     default = True),
        tmo.Chemical('SO2',
                     phase = 'g',
                     ),
    tmo.Chemical('CO', search_ID='CarbonMonoxide',
                 phase='g', Hf=-26400* 4.184),#Ref: lactic acid chemicals 
    tmo.Chemical('NH3', phase='g', Hf=-10963*4.184),
    tmo.Chemical('NO', search_ID='NitricOxide', phase='g'),
    tmo.Chemical('NO2', phase='g'),
    tmo.Chemical('H2S', phase='g', Hf=-4927*4.184),
    tmo.Chemical('SO2', phase='g'),
    tmo.Chemical('CaSO4',search_ID = 'Calcium_sulfate'),
    tmo.Chemical('Calcium_dihydroxide', phase='s', Hf=-235522*4.184),        
#WWTs sludge based on cellulosic.chemicals     
        tmo.Chemical('WWTsludge',
                     search_db = False,
                     formula="CH1.64O0.39N0.23S0.0035",
                     Hf=-23200.01*4.184,
                     phase = 'l',
                     default = True,
                     )
                 ])


chems['Monomethyl_azelate'].Pc = 2.39587E+06#Chemical compound generator DWSIM
chems['Monomethyl_azelate'].Tc = 837.971#Chemical compound generator DWSIM
chems['Monomethyl_azelate'].omega = 1.09913#Chemical compound generator DWSIM
chems['Monomethyl_azelate'].Tb = 650.2#Chemical compound generator DWSIM

chems.Methyl_palmitate.Cn.l.method = 'ROWLINSON_BONDI'
# chems.Methyl_oleate.Cn.g.method = 'LASTOVKA_SHAW'
chems.Palmitic_acid.Cn.l.method = 'ROWLINSON_BONDI'
chems.Stearic_acid.Cn.l.method = 'ROWLINSON_BONDI'
chems.Oleic_acid.Cn.method = 'ROWLINSON_BONDI'
chems.Octane.Cn.l.method = 'DADGOSTAR_SHAW' #This method works for hydrocarbons
chems.Natural_gas.Cn.l.method = 'DADGOSTAR_SHAW' #This method works for hydrocarbons
chems.OOO.copy_models_from(chems.PPP,['mu',
                                      ])
chems.LnLnLn.copy_models_from(chems.PPP,['mu'])
TAGs_with_unknown_props = [ 'LLL',
                            'OOL','LLO','SOO',
                            'PLO','PoOO','POO',
                            'POS','POP','PLS'
                          ]
for i in TAGs_with_unknown_props:
    chems[i].Tb = chems.OOO.Tb
    chems[i].copy_models_from(chems.PPP,['mu'])  

chems.Sodium_oleate.copy_models_from(chems.Oleic_acid,
                                     ['mu','Hvap','Psat','Cn','V'])

## The oxidative cleavage catalyst properties are based on cobalt chloride
chems.Cobalt_acetate_tetrahydrate.copy_models_from(chems.Cobalt_chloride,
                                                      ['sigma',                                                      
                                                       #'Hvap',
                                                       #'Psat',
                                                       'kappa',
                                                       #'Cn',
                                                       ])
## The density of cobalt acetate tetrahydrate is available in the literature, rho is 1730 Kg/m3
## https://scifinder-n.cas.org/searchDetail/substance/634d77a73c1f076117e95f61/substanceDetails
## Cobalt acetate dissolution reaction data
# V_of_cobalt_acetate_tetrahydrate = fn.rho_to_V(1730,chems['Cobalt_acetate_tetrahydrate'].MW)
# chems.Cobalt_acetate_tetrahydrate.V.add_model(V_of_cobalt_acetate_tetrahydrate,
                                                 # top_priority = True)
chems.Tungstic_acid.copy_models_from(chems.Tungsten,
                                     ['Hvap','Psat'])
chems.Sodium_tungstate.copy_models_from(chems.Tungsten,
                                     ['Hvap','Psat'])
## Tungstate ion default to properties
chems.Tungstate_ion.copy_models_from(chems.Tungsten,
                                     ['Hvap','Psat',
                                      'Cn','mu'
                                     ])
## Hydrogen ion defaults to properties of hydrogen
chems.Hydrogen_ion.copy_models_from(chems.Hydrogen,
                                       ['Hvap',
                                        'Psat',
                                        'mu',
                                        'Cn'])
for i in ['Hydrogen_ion','Tungstate_ion','Cobalt_ion',
          'Acetate_ion','Cobalt_acetate_tetrahydrate','Sodium_acetate',
          'Cobalt_hydroxide','Tungstic_acid','Hydrogen_peroxide','Sodium_tungstate']:
    V = fn.rho_to_V(rho=1e5, MW=chems[i].MW)
    chems[i].V.add_model(V, top_priority=True)
## Products of the precipitation reaction
## https://scifinder-n.cas.org/searchDetail/substance/634db2343c1f076117eb77ce/substanceDetails
chems.Calcium_hydroxide.V.add_model(fn.rho_to_V(2.34,
                                                chems.Calcium_hydroxide.MW),
                                                top_priority = True)
    
chems.Calcium_tungstate.copy_models_from(chems.Tungsten,
                                            ['Hvap',
                                             'Psat','Cn'])
chems.Calcium_tungstate.V.add_model(fn.rho_to_V(5800,
                                                chems.Calcium_tungstate.MW),
                                                top_priority = True)
chems.Calcium_acetate.V.add_method(chems.Calcium_acetate.MW*0.000001/1.5)#density from: https://pubchem.ncbi.nlm.nih.gov/compound/Calcium-acetate#section=Density
#TODO: change the below
chems.Calcium_hydroxide.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
chems.Calcium_chloride.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
chems.Calcium_acetate.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
## Modelling properties of dihydroxylated compounds as MDHSA
##TODO: find better assumptions for these below
chems.Methyl_dihydroxy_palmitate.copy_models_from(chems.MDHSA,
                                                  ['Hvap',
                                                   'Psat',
                                                   'Cn',
                                                   'V',
                                                   'mu'
                                                   ])
                                              
chems.Tetrahydroxy_octadecanoate.copy_models_from(chems.MDHSA,['Hvap','Psat',
                                                               'Cn','V','mu',
                                                 ])
chems.Tetrahydroxy_octadecanoate.Tc = 1379.78
chems.Tetrahydroxy_octadecanoate.Pc = 1483852.13
chems.Tetrahydroxy_octadecanoate.omega  = 0.84
chems.Hexahydroxy_octadecanoate.copy_models_from(chems.MDHSA,['Hvap','Psat',
                                                              'Cn','V','mu'])
chems.Hexahydroxy_octadecanoate.Pc = 1741912.65
chems.Hexahydroxy_octadecanoate.Tc = 1785.56
chems.Hexahydroxy_octadecanoate.omega = 0.28

chems.Monoester_MDHSA_PA.Pc = 728096.88
chems.Monoester_MDHSA_PA.Tc = 1377.87
chems.Monoester_MDHSA_PA.omega = 0.38
chems.Monoester_MDHSA_PA.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu'])

chems.Diester_MDHSA_PA.Pc = 469319.28
chems.Diester_MDHSA_PA.Tc = 1753.52
chems.Diester_MDHSA_PA.omega = -0.24
chems.Diester_MDHSA_PA.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu'])

chems.Monoester_MDHSA_MMA.Pc = 708841.01
chems.Monoester_MDHSA_MMA.Tc = 1547.14
chems.Monoester_MDHSA_MMA.omega = 0.19
chems.Monoester_MDHSA_MMA.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu'])

chems.Diester_MDHSA_MMA.Pc = 449627.03
chems.Diester_MDHSA_MMA.Tc = 2269.6
chems.Diester_MDHSA_MMA.omega = -0.46
chems.Diester_MDHSA_MMA.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu'])

chems.Cobalt_hydroxide.copy_models_from(chems.Cobalt_chloride,
                                        ['Hvap','Psat'])
chems.Cobalt_hydroxide.Tb = 100+273.15 #https://www.chembk.com/en/chem/Cobalt%20hydroxide
##Adding missing dortmund properties for the following chemicals
chems.Methyl_palmitoleate.Dortmund.set_group_counts_by_name({'CH3': 2,
                                                                'CH2': 11,
                                                                'CH2COO': 1,
                                                                'CH=CH' : 1})

chems.Nitrogen.PSRK.set_group_counts_by_name({'N2':1})

chems.LLL.Dortmund.set_group_counts_by_name({'CH3': 3,
                                             'CH2': 11+11+11+2,
                                             'CH=CH': 2+2+2,
                                             'CH2COO':3,
                                             'CH':1
                                             })

chems.OOL.Dortmund.set_group_counts_by_name({'CH3': 3,
                                             'CH2': 13+13+11+2,
                                             'CH=CH': 4,
                                             'CH2COO':3,
                                             'CH':1
                                             })
chems.LLO.Dortmund.set_group_counts_by_name({'CH3': 3,
                                             'CH2': 11+13 +11+2,
                                             'CH=CH': 5,
                                             'CH2COO': 3,                                             
                                             'CH':1
                                             })

chems.SOO.Dortmund.set_group_counts_by_name({'CH3':3,
                                             'CH2':15+13+13+2,
                                             'CH=CH':2,
                                             'CH2COO':3,
                                             'CH':1
                                             })
chems.PLO.Dortmund.set_group_counts_by_name({'CH3':3,
                                             'CH2':11+13+13+2,
                                             'CH=CH':2+1,
                                             'CH2COO': 1+1+1,
                                             'CH':1
                                             })
chems.PoOO.Dortmund.set_group_counts_by_name({'CH3':3,
                                              'CH2':11+13+13+2,
                                              'CH=CH':1+1+1,
                                              'CH2COO': 1+1+1,
                                              'CH':1                                                   
                                                    })
chems.POO.Dortmund.set_group_counts_by_name({'CH3':3,
                                             'CH2':13+13+13+2,
                                             'CH=CH':1,
                                             'CH2COO': 1+1+1,
                                             'CH':1  
                                             })
chems.POS.Dortmund.set_group_counts_by_name({'CH3':3,
                                             'CH2':13+13+15+2,
                                             'CH=CH':1,
                                             'CH2COO': 1,
                                             'CH':1  
                                             })
chems.POO.Dortmund.set_group_counts_by_name({'CH3':3,
                                             'CH2':13+13+13+2,
                                             'CH=CH':1+1,
                                             'CH2COO': 1+1+1,
                                             'CH':1  
                                             })

chems.POP.Dortmund.set_group_counts_by_name({'CH3':3,
                                    'CH2':13+13+13+2,
                                    'CH=CH':1,
                                    'CH2COO': 1+1+1,
                                    'CH':1  
                                              })
chems.PLS.Dortmund.set_group_counts_by_name({'CH3':3,
                                             'CH2':13+11+15+2,
                                             'CH=CH':2,
                                             'CH2COO': 1+1+1,
                                             'CH':1 })

chems.Methyl_dihydroxy_palmitate.Dortmund.set_group_counts_by_name({'CH3':2,
                                                        'CH2COO':1,
                                                        'CH2':13,
                                                        'OH(S)':2})
chems.Tetrahydroxy_octadecanoate.Dortmund.set_group_counts_by_name({'CH3':2,
                                                                    'CH2':15,
                                                                    'CH2COO':1,
                                                                    'OH(S)':4})
chems.Hexahydroxy_octadecanoate.Dortmund.set_group_counts_by_name({'CH3':2,
                                                                    'CH2':9,
                                                                    'CH2COO':1,
                                                                    'OH(S)':6})

chems.Monoester_MDHSA_PA.Dortmund.set_group_counts_by_name({'CH3':3,
                                                                 'CH2':13,
                                                                 'CH2COO':2,
                                                                 'OH(S)':1,
                                                                 'CH': 2})                                                                    
chems.Diester_MDHSA_PA.Dortmund.set_group_counts_by_name({'CH3':4,
                                                                 'CH2':25,
                                                                 'CH2COO':3,
                                                                 'CH': 2})   
   
chems.Diester_MDHSA_MMA.Dortmund.set_group_counts_by_name({'CH3':4,
                                                                 'CH2':23,
                                                                 'CH2COO':5,
                                                                 'CH': 2})   
chems.Monoester_MDHSA_MMA.Dortmund.set_group_counts_by_name({'CH3':3,
                                                                 'CH2':18,
                                                                 'CH2COO':3,
                                                                 'CH': 2,
                                                                 'OH(S)': 1,
                                                                 }) 
chems.Dipalmitin.Dortmund.set_group_counts_by_name({'CH3':2,
                                                    'CH2COO':2,
                                                    'CH2':13+13+2,
                                                    'CH':1,
                                                    'OH(P)':1
                                                    })
############################################################################################################### 
#The below Psat models return values in Pa (N/m)
# chems.OOO.Psat.add_method(f=OOO_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.LLL.Psat.add_method(f=LLL_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.OOL.Psat.add_method(f=OOL_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.LLO.Psat.add_method(f=LLO_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.SOO.Psat.add_method(f=SOO_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.PLO.Psat.add_method(f=PLO_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.PoOO.Psat.add_method(f=PoOO_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.POO.Psat.add_method(f=POO_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.POS.Psat.add_method(f=POS_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.POP.Psat.add_method(f=POP_CCPsat_model, Tmin=323.15, Tmax=573.15)
chems.PLS.Psat.add_method(f=PLS_CCPsat_model, Tmin=323.15, Tmax=573.15)
###############################################################################################33
# chems.OOO.Cn.l.add_method(f=OOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.LLL.Cn.add_method(f=OOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.LLO.Cn.add_method(f=OOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.OOL.Cn.add_method(f=OOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.SOO.Cn.add_method(f=SOO_Cnl_model,Tmin= 298.15, Tmax=453.15)
chems.PoOO.Cn.add_method(f=PoOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.PLO.Cn.add_method(f=PoOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.POO.Cn.add_method(f=PoOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.POP.Cn.add_method(f=POP_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.PLS.Cn.add_method(f=POP_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.POS.Cn.add_method(f=POP_Cnl_model, Tmin= 298.15, Tmax=453.15)
###############################################################################################33
#Adding the molar volumes
# chems.OOO.V.l.add_method(f=OOO_Vl_model, Tmin= 258.15, Tmax=516.15)
# TODO: check what methodP means ...in the V.l definition
chems.LLL.V.add_method(f=LLL_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.OOL.V.add_method(f=OOL_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.LLO.V.add_method(f=LLO_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.SOO.V.add_method(f=SOO_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.PLO.V.add_method(f=PLO_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.POO.V.add_method(f=POO_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.POS.V.add_method(f=POS_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.POP.V.add_method(f=POP_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.PLS.V.add_method(f=PLS_Vl_model, Tmin= 258.15, Tmax=516.15)
# chems.SSS.V.add_method(f=SSS_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.PoOO.V.add_method(f=PoOO_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.SSS.copy_models_from(chems.OOO, ['V','mu'])
chems.PPP.copy_models_from(chems.OOO, ['V'])
########Heat of formation ##########################################################################
#LiquidMethanol and Sodium_methoxide were added for cane biorefinery compatibility
LiquidMethanol = tmo.Chemical('Methanol').at_state(phase='l', copy=True)
chems.Sodium_methoxide.copy_models_from (LiquidMethanol,
                                            ['V', 'sigma',
                                              'kappa', 'Cn',
                                                'Psat','mu'])
#Properties of chemicals required by the boilerturbogenerator based on cane.chemicals.py
chems.Ash.Cn.add_model(0.09 * 4.184 * chems.Ash.MW) 
for chemical in [chems.Ash,chems.P4O10]:
    V = fn.rho_to_V(rho=1540, MW=chemical.MW)
    chemical.V.add_model(V, top_priority=True)
    
#The below NaOH and HCl are based on the cane refinery chemicals.py module
def create_new_chemical(ID, phase='s', **constants):
        solid = tmo.Chemical.blank(ID, phase=phase, **constants)
        chems.append(solid)
        return solid
    
NaOH = create_new_chemical('NaOH', formula='NaOH')
HCl = create_new_chemical('HCl', formula='HCl')
# Assuming that solubles don't occupy much volume
for chemical in (HCl, NaOH,):
        V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
        chemical.default()
NaOH.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
HCl.copy_models_from(tmo.Chemical('HCl'),['Psat'])        
        
chems.Ash.copy_models_from(tmo.Chemical('Water'),['Psat'])
chems.P4O10.copy_models_from(tmo.Chemical('Water'),['Psat'])
chems.WWTsludge.copy_models_from(tmo.Chemical('Water'),['Psat'])
    
for chemical in chems: chemical.default()
        
chems.compile()
chems.define_group('TAG', ('OOO','LLL','OOL',
                           'LLO','SOO','PLO',
                           'PoOO','POO','POS',
                           'POP','PLS','PPP',
                           'LnLnLn','SSS'))

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
                                 'Methyl_linolenate'))

#composition of VM_Naphtha based on https://www.cdc.gov/niosh/npg/npgd0664.html#:~:text=None%20reported%20%5BNote%3A%20VM%26P%20Naphtha%20is%20a%20refined,Exposure%20Routes%20inhalation%2C%20ingestion%2C%20skin%20and%2For%20eye%20contact
chems.define_group('VM_Naphtha',['Octane',
                                 'Cycloheptane',
                                 'Bicyclo_octane',
                                 'Toluene'], 
                                 composition = [0.55,
                                                0.30,
                                                0.02,
                                                0.12])
chems.define_group('COSOxNOxH2S',['NO', 'NO2',
                                  'SO2','CO', 'H2S'])

chems.define_group('Air', ['Oxygen', 'Nitrogen'],composition=[0.21,0.79])
chems.set_synonym('Water', 'H2O')
chems.set_synonym('Carbon_dioxide','CO2')
chems.set_synonym('Phosphatidylinositol','PL')
chems.set_synonym('MonoOlein', 'MAG')
chems.set_synonym('Dipalmitin', 'DAG')
chems.set_synonym('Hydrogen_ion', 'H+')
chems.set_synonym('Pelargonic_acid','Nonanoic_acid')
chems.set_synonym('HCl2','Liquid_HCl')
chems.set_synonym('CaSO4', 'Gypsum')

