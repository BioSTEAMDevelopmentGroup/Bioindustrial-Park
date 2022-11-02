import thermosteam as tmo
from thermosteam import Chemical
import biosteam as bst
# from biosteam import Unit, Stream, settings, main_flowsheet
chemicals_experimental= tmo.Chemicals([
#Raw materials to the process
#Feedstock: oleic acid   
                      tmo.Chemical('Oleic_acid',
                                   phase = 'l'),
#Impurities in the feedstock
#https://www.chemicalassociates.com/images/stories/virtuemart/product/sds/CA1331-High-Oleic-Vegetable-Oleic-Acid-SDS.pdf
                      tmo.Chemical('Linoleic_acid',search_ID = '60-33-3', phase = 'l'),#About 6%
                      tmo.Chemical('Palmitic_acid',search_ID = '57-10-3', phase = 'l'),#About 1%
                      tmo.Chemical('Stearic Acid',search_ID = '57-11-4',phase = 'l'),#About 1%
                      tmo.Chemical('Other_FAs', search_ID = '67701-02-4',phase = 'l'),#About 1%
                      
#Oxidative cleavage oxidants 
                      tmo.Chemical('Water'),
                      tmo.Chemical('Hydrogen_peroxide'),
#OXidative cleavage catalysts, it was not available in the database
#TODO: Essential properties like CAS, phase, volume were added
#TODO: ask yoel if phase of the catalyst needs to be added
#Ref for phosphotungstic acid: https://www.fishersci.com/shop/products/phosphotungstic-acid-hydrate-thermo-scientific/AA4011614#:~:text=Density%3A%202.852g%2Fcm%203%20at%2020%C2%B0C%3A%20Formula%20Weight%3A%205706.35,25g%3A%20Chemical%20Name%20or%20Material%3A%20Phosphotungstic%20acid%20hydrate
                      tmo.Chemical('Phosphotungstic_acid',
                                   search_db = False,
                                   formula="H3PW12O40",
                                   MW=2880.2,
                                   CAS='1343-93-7',
                                   phase = 'l',
                                   rho = 2.82       #From CAS finder                            
                                   ),  
#Intermediates of the process                      
                      tmo.Chemical('Nonanal'),
                      tmo.Chemical('Epoxy_stearic_acid', search_ID = '24560-98-3'),
                      tmo.Chemical('DHSA',phase = 's', search_ID = '120-87-6'), 
#Ref for oxononanoic acid: https://www.chemsrc.com/en/cas/2553-17-5_432314.html                      
                      tmo.Chemical('Oxononanoic_acid',
                                   search_db = False,
                                   MW = 2553-17-5,
                                   CAS = '2553-17-5',
                                   
                                   Tb = 552.15),                     
#Organic solvent for organic phase separation
                      tmo.Chemical('Ethyl_acetate'),  
#Products of the oxidative cleavage reaction
#Gaseous products
                      tmo.Chemical('Carbon_dioxide'),   

#Main products                    
                      tmo.Chemical('Nonanoic_acid'),
                      tmo.Chemical('Azelaic_acid'),
#Side products
                      tmo.Chemical('Octanal'),
                      tmo.Chemical('Octanoic_acid'),    
                      tmo.Chemical('Heptanal'),
                      tmo.Chemical('Heptanoic_acid'),
#Side products due to Linoleic acid impurity
                       tmo.Chemical('Malonic_acid'),
                       tmo.Chemical('Hexanoic_acid'),
                       
#Chemicals required for solvent extraction
                        tmo.Chemical('bicyclo_octane',search_ID = '6221-55-2'),
                        tmo.Chemical('Octane'),
                        tmo.Chemical('Cycloheptane'),
                        tmo.Chemical('toluene'),
#Heterogeneous catalyst chemicals
                        # tmo.Chemical('')              
                                  ]) 

#Adding properties of oxononaoic acid, important for distillation unit process
Dortmund = chemicals_experimental.Oxononanoic_acid.Dortmund
Dortmund.set_group_counts_by_name(dict(COOH=1, CH2=7, CHO=1))
chemicals_experimental['Oxononanoic_acid'].copy_models_from(chemicals_experimental['Azelaic_acid'],
                                            ['Hvap',
                                             'Psat',
                                             ])
Phosphoric_acid  = tmo.Chemical('Phosphoric_acid',phase = 'l')
chemicals_experimental['Phosphotungstic_acid'].copy_models_from(Phosphoric_acid,
                                                                ['Hvap',
                                                                 'Psat',
                                                                 'Cn',
                                                                 'V'])
#Setting all gas molar volume methods to IDEAL                                                                 
chemicals_experimental['Nonanoic_acid'].V.g.method_P = 'IDEAL'
chemicals_experimental['Azelaic_acid'].V.g.method_P = 'IDEAL'
chemicals_experimental['Malonic_acid'].V.g.method_P = 'IDEAL'





for chemical in chemicals_experimental: chemical.default()

chemicals_experimental.compile()
chemicals_experimental.set_synonym('Phosphotungstic_acid', 'fresh_Cat')
chemicals_experimental.set_synonym('Carbon_dioxide', 'CO2')
bst.settings.set_thermo(chemicals_experimental)

# for i in chemicals_process_1: 
#     if not chemicalsi.locked_state: i.V.g.method_P = 'IDEAL'   
# return chemicals_process_1

# for i in chemicals_experimental:
#     print(i,i.Psat)
    






