"""
Created on Sat Aug 20 21:47:53 2022
@author: Lavanya Kudli
"""

import thermosteam as tmo
from thermosteam import functional as fn
from biorefineries.oleochemicals.tag_properties import *

#chems is a list of all the chemicals used in the azelaic acid production process
chems = tmo.Chemicals([
#1.Biodiesel production
#Feedstock TAGs in the composition, properties modelled based on[1]
#TAGs of high oleic oils
        tmo.Chemical('OOO',#Triolein
                     search_ID = '122-32-7',
                     Hf = TAG_Hf(name = 'C18_unsaturated',double_bonds = 1),
                     phase_ref = 'l' ),
        tmo.Chemical('LnLnLn', #Trilinolenin 
                     search_ID = '14465-68-0',
                     Hf = TAG_Hf(name = 'C18_unsaturated',double_bonds = 3),
                     phase = 'l'),
        tmo.Chemical('LLL',#Trilinolein 
                      search_ID = ' 537-40-6',#Based on reference search in CAS Finder
                      search_db = False,
                      formula = 'C54H96O6',#Based on reference search in CAS Finder
                      Hf = TAG_Hf(name = 'C18_unsaturated',double_bonds = 2),
                      phase = 'l',
                      ),         
        tmo.Chemical('PPP',#Tripalmitin
                     search_ID = 'Tripalmitin',
                     Hf = TAG_Hf(name = 'Saturated',carbon_atoms = 16),
                     phase = 'l'),
        tmo.Chemical('SSS',#Tristearin
                     search_ID = 'Tristearin',
                     Hf = TAG_Hf(name = 'Saturated',carbon_atoms = 18),
                     phase = 'l'),
#For lipidcane module compatibility inported for using create_transesterification_and_biodiesel_separation_system()
        tmo.Chemical('MonoOlein',search_ID = '111-03-5'),
        tmo.Chemical('Dipalmitin'),
        tmo.Chemical('Monopalmitin'),
#Base used for neutralisation of free fatty acids in the crude oil mixture        
        tmo.Chemical('Sodium_hydroxide',
                     search_ID ='NaOH',  
                     ),     
#Chemical for acid degumming         
        tmo.Chemical('Citric_acid'),
#Chemical for representing the gumns (phospholipids)        
        tmo.Chemical('Phosphatidylinositol', 
                      formula='C47H83O13P',
                      search_db=False, 
                      CAS='383907-36-6',
                      MW = 	886.56,
                      default=True,
                      Hf=-1.779e6,#Hf value obtained from cane.chemicals.py from cane biorefinery   [2] 
                      phase = 'l',#phase based on cane.chemicals.py from cane biorefinery[2]
                      ), 
#Tranesterification section   
#Transesterification catalyst  
        tmo.Chemical('Sodium_methoxide',
                     formula ='NaOCH3',
                     phase = 's',#Phase based on cane.chemicals.py from cane biorefinery [2]  
                     default = True), 
#Tranesterification reagent        
        tmo.Chemical('Methanol'), 
#Tranesterification by-product        
        tmo.Chemical('Glycerol'),              
##Chemicals part of Biodiesel (both HoSun and HoySoy) product produced through transesterification
        tmo.Chemical('Methyl_oleate'),
        tmo.Chemical('Methyl_palmitate'),
        tmo.Chemical('Methyl_stearate'),
        tmo.Chemical('Methyl_linoleate'), 
        tmo.Chemical('Methyl_linolenate',search_ID = '301-00-8'),
        tmo.Chemical('Methyl_palmitoleate',search_ID ='1120-25-8'),
        tmo.Chemical('Methyl_myristate'),
        tmo.Chemical('Arachidic_acid_methyl_ester', search_ID = '2566-89-4'),
        tmo.Chemical('Methyl_behenate'),     
        
#Dihydroxylation section
#Reagents for dihydroxylation
        tmo.Chemical('Hydrogen_peroxide'),
        tmo.Chemical('Water',default = True),
#Product 1 (Methyl dihydroxy stearic acid)
#Dihydroxy product of methyl oleate       
        tmo.Chemical('MDHSA', 
                     search_ID = '1115-01-1', 
                     Hf = -947.7*1000, #Wikepedia value for stearic acid [3]
                     phase_ref ='l'
                    ),      
#Product 2 
#Tetrahydroxlated product of methyl linoleate (properties were obtained from DWSIM [4])
        tmo.Chemical('Tetrahydroxy_octadecanoate',
                     formula = 'C19H38O6',
                      search_ID = '61117-79-1', #Reaxys [5]
                      search_db = False, 
                      MW = 362.51,#Reaxys [5]
                      Tb = 1077.57,#DWSIM                
                      Hf = -3614.63*1000, #DWSIM
                      phase = 'l'
                     ),
#Product 3
#Hexahydroxlated product of methyl linolenate (properties were obtained from DWSIM [4]) 
        tmo.Chemical.blank('Hexahydroxy_octadecanoate',#Structure obtained from Reaxys
                           MW = 394.51,#Reaxys [5]
                            Tb = 1261.05,#DWSIM
                            formula = 'C19H38O8',
                            Hf = -508800.0,#Assumed to be same as Linolenic acid created above
                            phase = 'l'
                           ),
   
#Oxidative cleavage section
#Oxidants used 
        tmo.Chemical('Nitrogen',phase = 'g'),
        tmo.Chemical('Oxygen',phase = 'g'),        
#Products of oxidative cleavage
        tmo.Chemical('Monomethyl_azelate','l'),
        tmo.Chemical('Monomethyl_suberate',search_ID = '3946-32-5'),
        tmo.Chemical('Caprylic_acid'),
        tmo.Chemical('Hexanoic_acid'),
        tmo.Chemical('Heptanoic_acid'),
        tmo.Chemical('Hexanoic_acid'),
        tmo.Chemical('Malonic_acid'),
        tmo.Chemical('Adipic_acid'),
        tmo.Chemical('Pelargonic_acid'),
        tmo.Chemical('Propanoic_acid'),
        tmo.Chemical('Nonanal'),
        tmo.Chemical('Methyl_oxo_nonanoicacid',search_ID = '1931-63-1'),   
#Possible monoester of MDHSA with Pelargonic acid        
        tmo.Chemical.blank('Monoester_MDHSA_PA',#Structure obtained from Reaxys                     
                           CAS = '950910-15-3',
                           MW = 470.73,
                           Tb = 1084.12,#DWSIM [4]
                           formula = 'C28H54O5',#CAS generated compound
                           Hf = -2705.64*1000, #DWSIM [4]
                           phase = 'l'
                           ),
#Possible diester of MDHSA with Monomethyl azelate acid        
        tmo.Chemical.blank('Diester_MDHSA_PA',#Structure obtained from Reaxys
                            MW = 610.76,
                            Tb = 1274.15,#DWSIM [4]
                            formula = 'C37H70O6',#CAS generated compound
                            Hf = -2540.21*1000,#DWSIM [4]
                            phase = 'l'),   
#Possible monoester of MDHSA with Monomethyl azelate           
        tmo.Chemical.blank('Monoester_MDHSA_MMA',#Structure obtained from Reaxys                     
                            MW = 514.74,
                            Tb = 1183.29,#DWSIM [4]
                            formula = 'C29H54O7',#CAS generated compound
                            Hf =-2989.99*1000,#DWSIM [4]
                            phase = 'l'),
#Possible diester of MDHSA with Monomethyl azelate        
        tmo.Chemical.blank('Diester_MDHSA_MMA',#Structure obtained from Reaxys
                            MW = 698.98,
                            Tb =1472.49,#DWSIM [4]
                            formula = 'C39H70O10' ,
                            Hf =-2979.85*1000,#DWSIM [4]
                            phase = 'l'), 
#Gaseous product of oxidative cleavage from decarboxylation
        tmo.Chemical('Carbon_dioxide',phase = 'g'),     
#Catalysts used in both dihydroxylation and oxidative cleavage 
        tmo.Chemical('Tungstic_acid', 
                  Tb = 1746,#Tungstic_acid boiling point [6]
                  default = True,
                  formula = 'H2WO4',
                  phase = 'l' #because tungstic acid forms a hydrate in the presence of h2o2 [7]
                    ),
        tmo.Chemical('Cobalt_acetate_tetrahydrate',
                      search_ID = '6147-53-1',
                      Tb = 117.1+273.15,
                      Hf = -2.7151e+05, #Assumed to be same as cobalt chloride, declared above
                      phase = 'l',
                      ),
#Catalyst recovery section
#Reagents and products of the catalyst recovery section
        tmo.Chemical('Sodium_tungstate',
                    Tb = 1746,#same as Tungstic acid
                    phase = 'l',
                    ),
        tmo.Chemical('Sodium_chloride'), 
        tmo.Chemical('Hydrogen',default = True),  
        tmo.Chemical('Cobalt_chloride',
                     search_ID = '7646-79-9', 
                     phase = 's' #phase needed in order to copy models from this chemical
                     ),
        tmo.Chemical('Acetate',phase = 'l', default = True),
        tmo.Chemical('Calcium_hydroxide',default = True,phase = 'l'),
        tmo.Chemical('Calcium_chloride',
                     phase = 's',default = True),#phase needed in order to add molar volume
        tmo.Chemical.blank('Calcium_tungstate',
                       CAS = '7790-75-2',
                       MW = 287.92,
                       Tb = 5828.15, #Based on tungsten's BP
                       formula = 'CaWO4',
                       phase = 'l',
                       Hf = -1002.82*1000),#[8]
        tmo.Chemical('Calcium_acetate',phase = 'l', default = True),
        tmo.Chemical('Cobalt_hydroxide',phase = 's', default = True, Tb = 100+273.15), #Solid precipitate obtained after catalyst seperation which is used as such [7], Tb from [11]
        tmo.Chemical('HCl2',
                     search_ID = 'HCl'),
        tmo.Chemical('Sodium_acetate', phase = 'l'),
        tmo.Chemical('Chromic_acid',default = True,phase = 'l'),#Because properties in the same vertical column of a periodic table remain the same

#Hydrolysis section        
#Products of hydrolysis
        tmo.Chemical('Palmitic_acid'),
        tmo.Chemical('Stearic_acid'),
        tmo.Chemical('Oleic_acid'),
        tmo.Chemical('Sodium_oleate',Tb = 633 ,Hf = -764800.0,phase ='l'),#Hf based on Oleic acid declared above
        tmo.Chemical('Linoleic_acid', search_ID = '60-33-3'),
        tmo.Chemical('Linolenic_acid'),
        tmo.Chemical('Palmitoleic_acid', search_ID = '373-49-9'),
        tmo.Chemical('Azelaic_acid'),
#Resin for hydrolysis
##Sulfonated_polystyrene with a search_ID = '98-70-4', lacks Hvap data, hence polystyrene's properties were used
##Boiling point based on amberlyte
#Solid phase as the chemical beads are insoluble in water during hydrolysis
        tmo.Chemical('polystyrene_based_catalyst',
                     search_ID='Polystyrene',
                     Tb = 516.7+273.15,
                     phase = 's',
                     default=True),  
#Multistage countercurrent extraction        
#Solvent for countercurrent extraction of monocarboxylic acids
        tmo.Chemical('Heptane'),
#Facilities        
#Natural gas for heating purposes 
        # tmo.Chemical('Natural_gas',search_ID = 'CH4'), 
        tmo.Chemical('Methane',search_ID = 'CH4',phase = 'g'), 
        
#Chemicals required by the boilerturbogenerator (lactic acid chemicals from lactic acid bioerefinery) [10]
        tmo.Chemical('Ash',
                     MW = 1,
                     phase = 's',
                     search_db=False,
                     default = True,
                     Hf = 0),
        tmo.Chemical('P4O10',phase = 's',
                     default = True),
        tmo.Chemical('SO2',
                     phase = 'g',
                     ),
        tmo.Chemical('CO', search_ID='CarbonMonoxide',
                 phase='g', Hf=-26400* 4.184),
        tmo.Chemical('NH3', phase='g', Hf=-10963*4.184),
        tmo.Chemical('NO', search_ID='NitricOxide', phase='g'),
        tmo.Chemical('NO2', phase='g'),
        tmo.Chemical('H2S', phase='g', Hf=-4927*4.184),
        tmo.Chemical('SO2', phase='g'),
        tmo.Chemical('CaSO4',search_ID = 'Calcium_sulfate'),
        tmo.Chemical('Calcium_oxide',phase='s', Hf=-235522*4.184),     
        tmo.Chemical('BoilerChems', search_ID='DiammoniumPhosphate',
                                phase='l', Hf=0, HHV=0, LHV=0),
#WWTs sludge based on cellulosic.chemicals     
        tmo.Chemical('WWTsludge',
                     search_db = False,
                     formula="CH1.64O0.39N0.23S0.0035",
                     Hf=-23200.01*4.184,
                     phase = 'l',
                     default = True,
                     ),
                 ])

#Filling missing properties

#Properties of Monomethyl azelate obtained from DWSIM[4]
chems['Monomethyl_azelate'].Pc = 2.39587E+06#Chemical compound generator DWSIM
chems['Monomethyl_azelate'].Tc = 837.971#Chemical compound generator DWSIM
chems['Monomethyl_azelate'].omega = 1.09913#Chemical compound generator DWSIM
chems['Monomethyl_azelate'].Tb = 650.2#Chemical compound generator DWSIM

#Adding models for unknown TAG properties based on [1]
chems.LLL.Psat.add_method(f=TAG_Psat_model, Tmin=323.15, Tmax=573.15)
chems.LLL.Cn.add_method(f=LLL_Cnl_model, Tmin= 298.15, Tmax=453.15)
chems.LLL.V.add_method(f=LLL_Vl_model, Tmin= 258.15, Tmax=516.15)
chems.SSS.copy_models_from(chems.OOO, ['V'])
chems.PPP.copy_models_from(chems.OOO, ['V'])

#Changing liquid heat capacity methods to ensure compaibility with the process
chems.Methyl_palmitate.Cn.l.method = 'ROWLINSON_BONDI'
chems.Palmitic_acid.Cn.l.method = 'ROWLINSON_BONDI'
chems.Stearic_acid.Cn.l.method = 'ROWLINSON_BONDI'
chems.Oleic_acid.Cn.l.method = 'ROWLINSON_BONDI'
chems.Malonic_acid.Cn.g.method = 'LASTOVKA_SHAW'
chems.Malonic_acid.Cn.l.method = 'ROWLINSON_POLING'

# TODO:RuntimeError: <System: SYS19> <System: SYS17> <System: SYS16> <ShortcutColumn: D604> Failed to extrapolate integral of gas heat capacity method 'HEOS_FIT' between T=-998.6322733218122 to 175.61 K for component with CASRN '67-56-1'

#Making assumptions for missing viscosity of below TAGS 
TAGs_with_unknown_props = [ 'LLL','OOO','LnLnLn','SSS']
for i in TAGs_with_unknown_props:
    chems[i].copy_models_from(chems.PPP,['mu'])  

# Feedstock, reactants and products     
chems.Sodium_oleate.copy_models_from(chems.Oleic_acid,['mu','Hvap','Psat','Cn','V'])
chems.Methyl_linoleate.copy_models_from(chems.Linoleic_acid,['V'])
chems.Methyl_linolenate.copy_models_from(chems.Linolenic_acid,['V'])
chems.Malonic_acid.copy_models_from(chems.Adipic_acid,['mu'])
chems.Methyl_palmitate.copy_models_from(chems.Palmitic_acid,['mu'])
chems.Methyl_stearate.copy_models_from(chems.Stearic_acid,['mu'])
chems.Monomethyl_azelate.copy_models_from(chems.Azelaic_acid,['Cn','Hvap'])

## Modelling properties of dihydroxylated compounds, mono and diester side products as MDHSA                                    
chems.Tetrahydroxy_octadecanoate.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu',])
chems.Tetrahydroxy_octadecanoate.Tc = 1379.78 #DWSIM
chems.Tetrahydroxy_octadecanoate.Pc = 1483852.13#DWSIM
chems.Tetrahydroxy_octadecanoate.omega  = 0.84#DWSIM

chems.Hexahydroxy_octadecanoate.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu'])
chems.Hexahydroxy_octadecanoate.Pc = 1741912.65#DWSIM
chems.Hexahydroxy_octadecanoate.Tc = 1785.56#DWSIM
chems.Hexahydroxy_octadecanoate.omega = 0.28#DWSIM

chems.Monoester_MDHSA_PA.Pc = 728096.88#DWSIM
chems.Monoester_MDHSA_PA.Tc = 1377.87#DWSIM
chems.Monoester_MDHSA_PA.omega = 0.38#DWSIM
chems.Monoester_MDHSA_PA.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu'])

chems.Diester_MDHSA_PA.Pc = 469319.28#DWSIM
chems.Diester_MDHSA_PA.Tc = 1753.52#DWSIM
chems.Diester_MDHSA_PA.omega = -0.24#DWSIM
chems.Diester_MDHSA_PA.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu'])

chems.Monoester_MDHSA_MMA.Pc = 708841.01#DWSIM
chems.Monoester_MDHSA_MMA.Tc = 1547.14#DWSIM
chems.Monoester_MDHSA_MMA.omega = 0.19#DWSIM
chems.Monoester_MDHSA_MMA.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu'])

chems.Diester_MDHSA_MMA.Pc = 449627.03#DWSIM
chems.Diester_MDHSA_MMA.Tc = 2269.6#DWSIM
chems.Diester_MDHSA_MMA.omega = -0.46#DWSIM
chems.Diester_MDHSA_MMA.copy_models_from(chems.MDHSA,['Hvap','Psat','Cn','V','mu'])



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
#For lipidcane module compatibility inported for using create_transesterification_and_biodiesel_separation_system()
chems.Dipalmitin.Dortmund.set_group_counts_by_name({'CH3':2,
                                                    'CH2COO':2,
                                                    'CH2':13+13+2,
                                                    'CH':1,
                                                    'OH(P)':1
                                                    })

#Catalyst recovery chemicals
chems.Cobalt_acetate_tetrahydrate.copy_models_from(chems.Cobalt_chloride, ['sigma', 'kappa', 'mu'])
chems.Tungstic_acid.copy_models_from(chems.Chromic_acid,['Hvap','Psat','V'])
chems.Sodium_tungstate.copy_models_from(chems.Chromic_acid,['Hvap','Psat'])
chems.Calcium_tungstate.copy_models_from(chems.Chromic_acid,['Hvap','Psat','Cn'])


#set phase = s
chems.Calcium_hydroxide.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
chems.Calcium_chloride.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
chems.Calcium_acetate.copy_models_from(tmo.Chemical('NaOH'),['Psat'])
chems.Cobalt_hydroxide.copy_models_from(chems.Cobalt_chloride,['Hvap','Psat'])

#Adding molar volume 
chems.Calcium_tungstate.V.add_model(fn.rho_to_V(5800,
                                                chems.Calcium_tungstate.MW),
                                                top_priority = True)       
chems.Calcium_hydroxide.V.add_model(fn.rho_to_V(2.34*1000,
                                                chems.Calcium_hydroxide.MW),
                                                top_priority = True)
chems.Calcium_chloride.V.add_model(fn.rho_to_V(2.15*1000,#Wiki
                                                chems.Calcium_chloride.MW),
                                                top_priority = True)    

chems.Calcium_acetate.V.add_method(chems.Calcium_acetate.MW*0.000001/1.5)#density based on [12]


#LiquidMethanol and Sodium_methoxide were added for cane biorefinery compatibility[2]
#They are used as catalysts in the tranesterification section
LiquidMethanol = tmo.Chemical('Methanol').at_state(phase='l', copy=True)
chems.Sodium_methoxide.copy_models_from (LiquidMethanol,
                                         ['V', 'sigma',
                                          'kappa', 'Cn',
                                          'Psat','mu'])
chems.Methanol.sigma.method = 'BROCK_BIRD'
chems.Methanol.copy_models_from (LiquidMethanol,
                                         ['V',
                                          'kappa',
                                          'Psat','mu'])


#Properties of chemicals required by the boilerturbogenerator based on cane.chemicals.py [2]
chems.Ash.Cn.add_model(0.09 * 4.184 * chems.Ash.MW) 
for chemical in [chems.Ash,chems.P4O10]:
    V = fn.rho_to_V(rho=1540, MW=chemical.MW)
    chemical.V.add_model(V, top_priority=True)
    
#The below NaOH and HCl are based on the cane refinery chemicals.py module[2]
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
   
#Making assumptions for chemicals required/produced by the facilities   
chems.Ash.copy_models_from(tmo.Chemical('Water'),['Psat'])
chems.P4O10.copy_models_from(tmo.Chemical('Water'),['Psat'])
chems.WWTsludge.copy_models_from(tmo.Chemical('Water'),['Psat'])

#making assumption for Sfus
for i in ['Malonic_acid','Linoleic_acid',
          'Linolenic_acid','Azelaic_acid',
          'Methyl_oxo_nonanoicacid','Methanol',
          'Methyl_oleate',
          'Methyl_palmitate',
          'Methyl_stearate',
          'Methyl_linoleate',
          'Methyl_linolenate',
          'Water','MDHSA',
          'Monomethyl_azelate',
          'Palmitic_acid',
          'Stearic_acid',
          'Oleic_acid']:
    chems[i].Sfus = 22.0 #based on Sfus of water, no other value availabe
    #Sfus property is required when compressors are used [13]
for chemical in chems: chemical.default()
        
chems.compile()
chems.define_group('TAG', ('OOO','LLL',
                           'PPP',
                           'LnLnLn','SSS',
                           ))

chems.define_group('Lipid', ('OOO','LLL',
                           ))
chems.define_group('Oil', ('OOO','LLL',
                           ))

chems.define_group('Biodiesel', ('Methyl_oleate',
                                 'Methyl_palmitate',
                                 'Methyl_stearate',
                                 'Methyl_linoleate',
                                 'Methyl_linolenate',
                                 'Methyl_myristate',
                                 'Arachidic_acid_methyl_ester',
                                 'Methyl_behenate'
                                 ))


chems.define_group('COSOxNOxH2S',['NO', 'NO2',
                                  'SO2','CO', 'H2S'])

chems.define_group('Air', ['Oxygen', 'Nitrogen'],composition=[0.21,0.79])
chems.set_synonym('Water', 'H2O')
chems.set_synonym('Carbon_dioxide','CO2')
chems.set_synonym('Phosphatidylinositol','PL')
chems.set_synonym('MonoOlein', 'MAG')
chems.set_synonym('Dipalmitin', 'DAG')
chems.set_synonym('Pelargonic_acid','Nonanoic_acid')
chems.set_synonym('HCl2','Liquid_HCl')
chems.set_synonym('CaSO4', 'Gypsum')
chems.set_synonym('Calcium_oxide', 'Lime')
chems.set_synonym('Citric_acid', 'CitricAcid')
chems.set_synonym('Caprylic_acid','Octanoic_acid')

###References###
#[1] DOI: 10.1021/ie100160v
#[2] https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/cane/chemicals.py
#[3] https://en.wikipedia.org/wiki/Stearic_acid
#[4] DWSIM https://dwsim.org/
#[5] https://www.reaxys.com/#/search/quick
#[6] https://en.wikipedia.org/wiki/Tungstic_acid
#[7] Process for recovering cobalt and tungsten from reaction liquors, US5599514A
#[8] Value for calcium hydroxide https://www.chemeurope.com/en/encyclopedia/Standard_enthalpy_change_of_formation_%28data_table%29.html
#[9] https://www.chemsrc.com/en/cas/39389-20-3_843683.html#:~:text=amberlyst%28r%29%2015%20CAS%20Number%3A%2039389-20-3%3A%20Molecular%20Weight%3A%20314.39900%3A,Point%3A%20266.3%C2%BAC%3A%20Symbol%3A%20GHS07%3A%20Signal%20Word%3A%20Warning%20%C3%97
#[10] https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/lactic/_chemicals.py
#[11] #https://www.chembk.com/en/chem/Cobalt%20hydroxide),
#[12] density from: https://pubchem.ncbi.nlm.nih.gov/compound/Calcium-acetate#section=Density
#[13] #https://en.wikipedia.org/wiki/Water_%28data_page%29