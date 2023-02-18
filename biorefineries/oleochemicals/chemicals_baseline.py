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
#TODO: how to tackle? RuntimeError: Failed to extrapolate integral of liquid heat capacity method 'ZABRANSKY_SPLINE_C' between T=298.15 to nan K for component with CASRN '122-32-7'
#TODO: ask if every inorganic chemical neeeds Dortmund groups
#TODO: the important chemicals with no PSAT are PL,CaOH2, CaCl2, 62-54-4, Ash,HCl,NaOH, WWTs_sludge, P4O10
#chems is a list of all the chemicals used in this biorefinery
#Hf_formation_(OOO):-76.494*18 - 815.18
#Hf_formation_(PPP):- -59.571*16 - 1358.7
#Hf_formation_(PoPoPo):- -76.494*16 - 815.18
#Hf_formation_(SSS):- -59.571*18 - 1358.7
#Hf_formation_(LLL):-  -316.67*2 - 2466.4
#Hf_formation_(LnLnLn):-  -316.67*3 - 2466.4
chems = tmo.Chemicals([
#Chemicals used in biodiesel production
        tmo.Chemical('Sodium_methoxide',formula ='NaOCH3',phase = 'l',default = True),   
        tmo.Chemical('Methanol'),
#Chemical for acid degumming 
        tmo.Chemical('Citric_acid'),

##chemical for representing phospholipids, taken directly from lipidcane biorefinery
        tmo.Chemical('Phosphatidylinositol', 
                      formula='C47H83O13P',
                      search_db=False, 
                      CAS='383907-36-6',
                      MW = 	886.56,
                      default=True,
                      Hf=1000*(-316.67*4 - 2466.4),#TODO: find a better assumption,assuming LLL with 4 bonds HF
                      phase = 'l',
                      ),
#Crude glycerol is a product of transesterification
        tmo.Chemical('Glycerol'),
#TAGs of High oleic sunflower oil
        tmo.Chemical('OOO', search_ID = '122-32-7',Hf = 1000*(-76.494*18 - 815.18) ),#from cane>chemicals
        tmo.Chemical('LnLnLn', search_ID = '537-40-6',Hf = 1000*(-316.67*3 - 2466.4)),#This is Trilinolein
#Below TAG's not a part of the database    
         tmo.Chemical('LLL',
                      search_ID = '7049-66-3',#Based on reference search in CAS Finder
                      search_db = False,
                      formula = 'C54H96O6',#Based on reference search in CAS Finder
                      phase = 'l',
                      Hf = (((2/3)*(-76.494*18 - 815.18))+((1/3)*(-316.67*2 - 2466.4)))*1000 
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
        
##Chemicals part of Biodiesel
        tmo.Chemical('Methyl_oleate'),
        tmo.Chemical('Methyl_palmitate'),
        tmo.Chemical('Methyl_stearate'),
        tmo.Chemical('Methyl_linoleate'), 
        tmo.Chemical('Methyl_palmitoleate',search_ID ='1120-25-8'),

#Dihydroxylation chemicals
            tmo.Chemical('Hydrogen_peroxide'),
            tmo.Chemical('Water'),
#Dihydroxylation products
#Product 1       
        tmo.Chemical('MDHSA', 
                     search_ID = '1115-01-1', 
                     Hf = -947.7*1000 #Wikepedia value for https://en.wikipedia.org/wiki/Stearic_acid
                    ),
#Product 2
#Ref for dihydroxy_palmitic_acid: https://www.chemsrc.com/en/cas/29242-09-9_803992.html    
        tmo.Chemical('Dihydroxy_palmitic_acid',
                 search_ID = '29242-09-9',
                 search_db = False,
                 Tb = 458 + 273.15,#CAS finder: 458.0±25.0 °C,Press: 760 Torr
                 formula = 'C16H32O4',
                 MW = 288.42+273.15,
                 Hf = -892*1000 #Wikepedia value for https://en.wikipedia.org/wiki/Palmitic_acid
                 ),
#Product 3
        tmo.Chemical('Methyl_9_10_dihydroxylinoleate',
                 search_ID = '62071-03-8', 
                 search_db = False, 
                 MW = 328.492,
                 Tb = 449.4+273.15,#CAS finder: Tb =449.4±35.0 °C,Press: 760 Torr
                 formula = 'C19H36O4',
                 Hf = -634.7*1000 #TODO: Ref:https://webbook.nist.gov/cgi/cbook.cgi?ID=C60333&Mask=2
                 ),
    
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
        tmo.Chemical('Azelaic_acid'),

# Oxidants used and other gaseous products
        tmo.Chemical('Nitrogen',phase = 'g'),
        tmo.Chemical('Oxygen',phase = 'g'),
        tmo.Chemical('Carbon_dioxide',phase = 'g'),
        
#Solvent for countercurrent extraction of monocarboxylic acids
        tmo.Chemical('Octane'),
        tmo.Chemical('Cycloheptane'),
        tmo.Chemical('Bicyclo_octane',search_ID = '6221-55-2'),
        tmo.Chemical('Toluene'),        
     
## Catalysts used in the process   
## Tungstic_acid boiling point: https://en.wikipedia.org/wiki/Tungstic_acid
       tmo.Chemical('Tungsten',default = True),
       tmo.Chemical('Tungstic_acid', 
                  Tb = 1746,
                  default = True,
                  formula = 'H2WO4',
                  phase = 'l'),
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
                      Hf = 0.45 #TODO: change
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
                       Hf = 0.45),#TODO: change    
        tmo.Chemical('Calcium_acetate',phase = 'l', default = True),
        tmo.Chemical('Cobalt_hydroxide',phase = 'l', default = True),
        tmo.Chemical('HCl2',phase = 'l',search_ID = 'HCl'),
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
#WWTs sludge based on cellulosic.chemicals     
        tmo.Chemical('WWTsludge',
                     search_db = False,
                     formula="CH1.64O0.39N0.23S0.0035",
                     Hf=-23200.01*4.184,
                     phase = 'l',
                     default = True,
                     )
                 ])

#Fitting data for MMA based on ChemSep
# Ts = [i + 273.15 for i in  (148, 159, 120, 185.5, )]
# Psats = [i / 760 * 101325 for i in (1, 3, 0.03, 11, )]
# res, stats = TDependentProperty.fit_data_to_model(Ts=Ts, data=Psats, model='Antoine', do_statistics=True, multiple_tries=True, model_kwargs={'base': 10.0})
# method = 'ANTOINE_POLING'
# chems['Monomethyl_azelate'].Psat.ANTOINE_POLING_coefs = res['A'], res['B'], res['C']
# chems['Monomethyl_azelate'].Psat.all_methods.add(method)
# chems['Monomethyl_azelate'].Psat.method = method
# chems['Monomethyl_azelate'].Psat.T_limits[method] = (100, chems['Monomethyl_azelate'].Psat.Tc)
chems['Monomethyl_azelate'].Pc = 2.39587E+06#Chemical compound generator DWSIM
chems['Monomethyl_azelate'].Tc = 837.971#Chemical compound generator DWSIM
chems['Monomethyl_azelate'].omega = 1.09913#Chemical compound generator DWSIM
chems['Monomethyl_azelate'].Tb = 650.2#Chemical compound generator DWSIM

#properties of TAGs that are missing
# chems.OOO.Cn.l.method = 'DADGOSTAR_SHAW'
# chems.PPP.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Methyl_palmitate.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Palmitic_acid.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Stearic_acid.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Oleic_acid.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Octane.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Natural_gas.Cn.l.method = 'DADGOSTAR_SHAW'
# chems.Cobalt_chloride.Cn.l.method = 'ROWLINSON_BONDI'


chems.OOO.copy_models_from(chems.PPP,['mu',
                                      # 'Cn'
                                      ])
chems.LnLnLn.copy_models_from(chems.PPP,['mu'])
TAGs_with_unknown_props = [ 'LLL',
                            'OOL','LLO','SOO',
                            'PLO','PoOO','POO',
                            'POS','POP','PLS'
                          ]
for i in TAGs_with_unknown_props:
    # chems[i].copy_models_from(chems.PPP,[
    #                                     #'Hvap',
    #                                     'Psat',
    #                                      #'Cn'
    #                                      ])
    chems[i].Tb = chems.OOO.Tb
    chems[i].copy_models_from(chems.PPP,['mu'])
    # chems[i].V.add_model(fn.rho_to_V(rho=chems.OOO.rho('l', 298.15, 101325),
    #                                  MW=chems[i].MW))
    
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
V_of_cobalt_acetate_tetrahydrate = fn.rho_to_V(1730,chems['Cobalt_acetate_tetrahydrate'].MW)
chems.Cobalt_acetate_tetrahydrate.V.add_model(V_of_cobalt_acetate_tetrahydrate,
                                                 top_priority = True)
chems.Tungstic_acid.copy_models_from(chems.Tungsten,
                                     ['Hvap','Psat'])
# chems.Tungstic_acid.V.add_model(fn.rho_to_V(rho=5.59,#https://en.wikipedia.org/wiki/Tungstic_acid
                                            # MW = chems.Tungstic_acid))
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
chems.Cobalt_ion.copy_models_from(chems.Cobalt,
                                       ['Hvap',
                                        'Psat'])

for i in ['Hydrogen_ion','Tungstate_ion','Cobalt_ion']:
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
chems.Acetate_ion.copy_models_from((chems.Acetate),
                                      ['mu',
                                       'Psat',
                                       'Hvap',
                                       'sigma',
                                       ])
Acetate_rho = chems.Acetate.rho(298.15,101325)
chems.Acetate_ion.V.add_method(chems.Acetate_ion.MW*10E-3/Acetate_rho)
chems.Calcium_acetate.V.add_method(chems.Calcium_acetate.MW*0.000001/1.5)#density from: https://pubchem.ncbi.nlm.nih.gov/compound/Calcium-acetate#section=Density
chems.Cobalt_hydroxide.V.add_method(chems.Cobalt_hydroxide.MW*1.0E-6/3.597)
# chems.Cobalt_hydroxide.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
#TODO: change the below
chems.Calcium_hydroxide.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
chems.Calcium_chloride.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
chems.Calcium_acetate.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
## Modelling properties of dihydroxylated compounds as MDHSA
##TODO: find better assumptions for these below
chems.Dihydroxy_palmitic_acid.copy_models_from(chems.MDHSA,
                                                  ['Hvap',
                                                   'Psat',
                                                   'Cn',
                                                   'V',
                                                   'mu'
                                                   ])
                                              
chems.Methyl_9_10_dihydroxylinoleate.copy_models_from(chems.MDHSA,
                                                         ['Hvap',
                                                          'Psat',
                                                          'Cn',
                                                          'V',#rho_methyldihydroxylinoleate = 0.980±0.06 g/cm3,Temp: 20 °C; Press: 760 Torr
                                                              #rho similar to MDHSA
                                                          'mu',
                                                          ])
chems.Methyl_9_10_dihydroxylinoleate.Pc = chems.MDHSA.Pc
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

chems.Dihydroxy_palmitic_acid.Dortmund.set_group_counts_by_name({'CH3':1,
                                                        'COOH':1,
                                                        'CH2':12,
                                                        'CH':2,
                                                        'OH(S)':2})
chems.Methyl_9_10_dihydroxylinoleate.Dortmund.set_group_counts_by_name({'CH3':2,
                                                               'CH2':11,
                                                               'CH=CH':1,
                                                               'CH':2,
                                                               'CH2COO':1,
                                                               'OH(S)':2})
chems.Dipalmitin.Dortmund.set_group_counts_by_name({'CH3':2,
                                                    'CH2COO':2,
                                                    'CH2':13+13+2,
                                                    'CH':1,
                                                    'OH(P)':1
                                                    })
###############################################################################################################33   
# def OOO_CCPsat_model(T):
#       R = 8.314
#       theta = 273.15
#       ln10 = 2.30258509
#       return (-91320000/(R*theta*ln10))+ ((169240000/R*ln10)((1/theta) - (1/T)))
# chems.OOO.Psat.add_method(f=OOO_CCPsat_model, Tmin=323.15, Tmax=573.15)

def LLL_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-91320000/(R*theta*ln10))+ ((169240000/R*ln10)((1/theta) - (1/T)))
chems.LLL.Psat.add_method(f=LLL_CCPsat_model, Tmin=323.15, Tmax=573.15)

def OOL_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-91320000/(R*theta*ln10))+ ((169240000/R*ln10)((1/theta) - (1/T)))
chems.OOL.Psat.add_method(f=OOL_CCPsat_model, Tmin=323.15, Tmax=573.15)

def LLO_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-91320000/(R*theta*ln10))+ ((169240000/R*ln10)((1/theta) - (1/T)))
chems.LLO.Psat.add_method(f=LLO_CCPsat_model, Tmin=323.15, Tmax=573.15)

def SOO_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-91320000/(R*theta*ln10))+ ((169240000/R*ln10)((1/theta) - (1/T)))
chems.SOO.Psat.add_method(f=SOO_CCPsat_model, Tmin=323.15, Tmax=573.15)

def PLO_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)((1/theta) - (1/T)))
chems.PLO.Psat.add_method(f=PLO_CCPsat_model, Tmin=323.15, Tmax=573.15)

def PoOO_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)((1/theta) - (1/T)))
chems.PoOO.Psat.add_method(f=PoOO_CCPsat_model, Tmin=323.15, Tmax=573.15)

def POO_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)((1/theta) - (1/T)))
chems.POO.Psat.add_method(f=POO_CCPsat_model, Tmin=323.15, Tmax=573.15)

def POS_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)((1/theta) - (1/T)))
chems.POS.Psat.add_method(f=POS_CCPsat_model, Tmin=323.15, Tmax=573.15)

def POP_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-87180000/(R*theta*ln10))+ ((164240000/R*ln10)((1/theta) - (1/T)))
chems.POP.Psat.add_method(f=POP_CCPsat_model, Tmin=323.15, Tmax=573.15)

def PLS_CCPsat_model(T):
      R = 8.314
      theta = 273.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)((1/theta) - (1/T)))
chems.PLS.Psat.add_method(f=PLS_CCPsat_model, Tmin=323.15, Tmax=573.15)
###############################################################################################33
#The model returns values in J/mol.K
# def OOO_Cnl_model(T):
#       return 3*(397600 + 540.89*T) + (61355 + 148.23*T)
# chems.OOO.Cn.add_method(f=OOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
# chems.LLL.Cn.add_method(f=OOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
# chems.LLO.Cn.add_method(f=OOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
# chems.OOL.Cn.add_method(f=OOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
#The below Psat models return values in Pa (N/m)
def SOO_Cnl_model(T):
      return (2*(397600 + 540.89*T) + (366930 + 685.76*T)+ (61355 + 148.23*T))*(1/1000)
chems.SOO.Cn.add_method(f=SOO_Cnl_model,Tmin= 298.15, Tmax=453.15)

def PoOO_Cnl_model(T):
      return (2*(397600 + 540.89*T) +(330360 + 616.35*T) + (61355 + 148.23*T))*(1/1000)
chems.PoOO.Cn.add_method(f=PoOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.PLO.Cn.add_method(f=PoOO_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.POO.Cn.add_method(f=PoOO_Cnl_model, Tmin= 298.15, Tmax=453.15)

def POP_Cnl_model(T):
      return (1*(397600 + 540.89*T) +2*(330360 + 616.35*T) + (61355 + 148.23*T))*(1/1000)
chems.POP.Cn.add_method(f=POP_Cnl_model, Tmin= 298.15, Tmax=453.15)

def PLS_Cnl_model(T):
      return (1*(397600 + 540.89*T) +1*(330360 + 616.35*T) + 1*(366930 + 685.76*T) + (61355 + 148.23*T))*(1/1000)
chems.PLS.Cn.add_method(f=POP_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.POS.Cn.add_method(f=POP_Cnl_model, Tmin= 298.15, Tmax=453.15)
###############################################################################################33
#Adding the molar volumes
# def OOO_Vl_model(T):
#       return (3*((1 + 0.0009865*T)/4.2924)) + ((1 + 20.048*T)/0.00076923)
# chems.OOO.V.add_method(f=OOO_Vl_model, Tmin= 258.15, Tmax=516.15)
def LLL_Vl_model(T):
      return ((3*((1 + 0.00074102*T)/4.1679)) + ((1 + 20.048*T)/0.00076923))*(1/1000)
chems.LLL.V.add_method(f=LLL_Vl_model, Tmin= 258.15, Tmax=516.15)
def OOL_Vl_model(T):
      return ((1*((1 + 0.00074102*T)/4.1679)) + (2*((1 + 0.0009865*T)/4.2924))+  ((1 + 20.048*T)/0.00076923))*1/1000
chems.OOL.V.add_method(f=OOL_Vl_model, Tmin= 258.15, Tmax=516.15)
def LLO_Vl_model(T):
      return ((2*((1 + 0.00074102*T)/4.1679)) + (1*((1 + 0.0009865*T)/4.2924))+  ((1 + 20.048*T)/0.00076923))*(1/1000)
chems.LLO.V.add_method(f=LLO_Vl_model, Tmin= 258.15, Tmax=516.15)
def SOO_Vl_model(T):
      return ((1*((1 + 0.0014091*T)/4.6326)) + (2*((1 + 0.0009865*T)/4.2924))+  ((1 + 20.048*T)/0.00076923))*(1/1000)
chems.SOO.V.add_method(f=SOO_Vl_model, Tmin= 258.15, Tmax=516.15)
def PLO_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (1*((1 + 0.00074102*T)/4.1679)) + (1*((1 + 0.0009865*T)/4.2924)) + ((1 + 20.048*T)/0.00076923))*(1/1000)
chems.PLO.V.add_method(f=PLO_Vl_model, Tmin= 258.15, Tmax=516.15)
def PoOO_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (2*((1 + 0.0009865*T)/4.2924)) + ((1 + 20.048*T)/0.00076923))*(1/1000)
chems.PoOO.V.add_method(f=PoOO_Vl_model, Tmin= 258.15, Tmax=516.15)
def POO_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (2*((1 + 0.0009865*T)/4.2924)) + ((1 + 20.048*T)/0.00076923))*(1/1000)
chems.POO.V.add_method(f=POO_Vl_model, Tmin= 258.15, Tmax=516.15)
def POS_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (1*((1 + 0.0009865*T)/4.2924)) +(1*((1 + 0.0014091*T)/4.6326))+ ((1 + 20.048*T)/0.00076923))*(1/1000)
chems.POS.V.add_method(f=POS_Vl_model, Tmin= 258.15, Tmax=516.15)
def POP_Vl_model(T):
      return ((2*((1 + 0.0013008*T)/5.0524))+ (1*((1 + 0.0009865*T)/4.2924)) + ((1 + 20.048*T)/0.00076923))*(1/1000)
chems.POP.V.add_method(f=POP_Vl_model, Tmin= 258.15, Tmax=516.15)
def PLS_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (1*((1 + 0.00074102*T)/4.1679)) +(1*((1 + 0.0014091*T)/4.6326))+  ((1 + 20.048*T)/0.00076923))*(1/1000)
chems.PLS.V.add_method(f=PLS_Vl_model, Tmin= 258.15, Tmax=516.15)

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
# Solubles don't occupy much volume
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
#TODO: add viscosity for all the catalyst related chems
        
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
chems.set_synonym('Hydrogen_ion', 'H+')
chems.set_synonym('Pelargonic_acid','Nonanoic_acid')
chems.set_synonym('HCl2','Liquid_HCl')
# chems.set_synonyn('WWTs_sludge','Biomass')
# chems.show()

#Dict for gibbs free energies
# {GOOO : 91320000,GLLL : 91320000,
#  GOOL : 91320000,GLLO : 91320000,
#  GSOO : 91320000,GPLO : 89250000,
#  GPoOO : 89250000,GPOO : 89250000,
#  GPOS : 89250000,GPOP : 87180000,
#  GPLS : 89250000}
#Dict for enthalpy of vapourisation
# {HOOO : 169240000,HLLL : 169240000,
#  HOOL : 169240000,HLLO : 169240000,
#  HSOO : 169240000,HPLO : 166740000,
#  HPoOO :166740000,HPOO : 166740000,
#  HPOS : 166740000,HPOP : 164240000,
#  HPLS : 166740000}