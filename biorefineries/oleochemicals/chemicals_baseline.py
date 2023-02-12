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
                      default=True,
                      Hf=-1.779e6, # Assume same as TAG on a weight basis
                      phase = 'l',
                      ),
#Crude glycerol is a product of transesterification
        tmo.Chemical('Glycerol'),
#TAGs of High oleic sunflower oil
        tmo.Chemical('OOO', search_ID = '122-32-7',Hf = -1776e3),#from cane>chemicals
        tmo.Chemical('LLL', search_ID = '537-40-6',Hf = -1776e3),#TODO: Don't Asumming to be the same as OOO  
#Below TAG's not a part of the database    
        tmo.Chemical('OOL',
                     search_ID = '28880-78-6',#Based on reference search in CAS Finder
                     search_db = False,
                     formula = 'C57H102O6',#Based on reference search in CAS Finder
                     phase = 'l',
                     Hf = -1776e3),           
        tmo.Chemical('LLO', 
                         search_ID = '28409-91-8',#Based on reference search in CAS Finder
                         search_db = False,
                         formula = 'C57H100O6',#Based on reference search in CAS Finder
                         phase = 'l',
                         Hf = -1776e3),  
        tmo.Chemical('SOO',
                         search_ID = '29590-02-1',#Based on reference search in CAS Finder
                         search_db = False,
                         formula = 'C57H106O6',#Based on reference search in CAS Finder
                         phase = 'l',
                         Hf = -1776e3),    
        tmo.Chemical('PLO', 
                         search_ID = '26836-35-1',#Based on reference search in CAS Finder
                         search_db = False,
                         formula = 'C55H100O6',#Based on reference search in CAS Finder
                         phase = 'l',
                         Hf = -1776e3),         
        tmo.Chemical('PoOO',
                     search_ID = '38703-17-2',#Based on reference search in CAS Finder
                     search_db = False,
                     formula = 'C55H100O6',#Based on reference search in CAS Finder
                     phase = 'l',
                     Hf = -1776e3),
        tmo.Chemical('POO',
                 search_ID = '27071-84-7',#Based on reference search in CAS Finder
                 search_db = False,
                 formula = 'C55H102O6',#Based on reference search in CAS Finder
                 phase = 'l',
                 Hf = -1776e3),
    
        tmo.Chemical('POS', 
                 search_ID = '26836-31-7',#Based on reference search in CAS Finder
                 search_db = False,
                 formula = 'C55H104O6', #Based on reference search in CAS Finder
                 phase = 'l',
                 Hf = -1776e3),
    
        tmo.Chemical('POP',
                 search_ID = '28409-94-1',#Based on reference search in CAS Finder
                 search_db = False,
                 formula = 'C53H100O6',#Based on reference search in CAS Finder
                 phase = 'l',
                 Hf = -1776e3),
    
        tmo.Chemical('PLS',
                 search_ID = '26836-32-8',#Based on reference search in CAS Finder 
                 search_db = False,
                 formula = 'C55H102O6',#Based on reference search in CAS Finder
                 phase = 'l',
                 Hf = -1776e3),
        tmo.Chemical('PPP',
                     search_ID = 'Tripalmitin',
                     Hf = -1776e3),
        
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
                    ),
#Product 2
#Ref for dihydroxy_palmitic_acid: https://www.chemsrc.com/en/cas/29242-09-9_803992.html    
        tmo.Chemical('Dihydroxy_palmitic_acid',
                 search_ID = '29242-09-9',
                 search_db = False,
                 Tb = 458 + 273.15,#CAS finder: 458.0±25.0 °C,Press: 760 Torr
                 formula = 'C16H32O4',
                 MW = 288.42+273.15,
                 Hf = -1776e3
                 ),
#Product 3
        tmo.Chemical('Methyl_9_10_dihydroxylinoleate',
                 search_ID = '62071-03-8', 
                 search_db = False, 
                 MW = 328.492,
                 Tb = 449.4+273.15,#CAS finder: Tb =449.4±35.0 °C,Press: 760 Torr
                 formula = 'C19H36O4',
                 Hf = -1776e3
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
chems.OOO.Cn.l.method = 'DADGOSTAR_SHAW'
chems.PPP.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Methyl_palmitate.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Palmitic_acid.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Stearic_acid.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Oleic_acid.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Octane.Cn.l.method = 'DADGOSTAR_SHAW'
chems.Natural_gas.Cn.l.method = 'DADGOSTAR_SHAW'
# chems.Cobalt_chloride.Cn.l.method = 'ROWLINSON_BONDI'


chems.OOO.copy_models_from(chems.PPP,['mu','Cn'])
chems.LLL.copy_models_from(chems.PPP,['mu'])
TAGs_with_unknown_props = ['OOL','LLO','SOO','PLO','PoOO',
                           'POO','POS','POP','PLS']
for i in TAGs_with_unknown_props:
    chems[i].copy_models_from(chems.PPP,['Hvap','Psat','Cn'])
    chems[i].Tb = chems.OOO.Tb
    chems[i].copy_models_from(chems.PPP,['mu'])
    chems[i].V.add_model(fn.rho_to_V(rho=chems.OOO.rho('l', 298.15, 101325),
                                     MW=chems[i].MW))
    
## The oxidative cleavage catalyst properties are based on cobalt chloride
chems.Cobalt_acetate_tetrahydrate.copy_models_from(chems.Cobalt_chloride,
                                                      ['sigma','Hvap','Psat',
                                                       'kappa', 'Cn',
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
