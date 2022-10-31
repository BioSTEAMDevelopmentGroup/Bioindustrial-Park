import thermosteam as tmo
from thermosteam import Chemical
import biosteam as bst
# from biosteam import Unit, Stream, settings, main_flowsheet
chemicals_experimental= tmo.Chemicals([
#Raw materials to the process
#Feedstock: oleic acid   
                      tmo.Chemical('Oleic_acid'),
#Impurities in the feedstock
#https://www.chemicalassociates.com/images/stories/virtuemart/product/sds/CA1331-High-Oleic-Vegetable-Oleic-Acid-SDS.pdf
                      tmo.Chemical('Linoleic_acid',search_ID = '60-33-3'),#About 6%
                      tmo.Chemical('Palmitic_acid',search_ID = '57-10-3'),#About 1%
                      tmo.Chemical('Stearic Acid',search_ID = '57-11-4'),#About 1%
                      tmo.Chemical('Other_FAs', search_ID = '67701-02-4'),#About 1%
                      
#Oxidative cleavage oxidants 
                      tmo.Chemical('Water'),
                      tmo.Chemical('Hydrogen_peroxide'),
#OXidative cleavage catalysts, it was not available in the database
#TODO: Essential properties like CAS, phase, volume were added
                      tmo.Chemical('Phosphotungstic_acid',
                                   search_db = False,
                                   # formula="H3PW12O40",
                                   # MW=2880.2,
                                   #CAS = '7783-03-1',
                                   CAS='1343-93-7',
                                   phase = 'l',
                                   ),  
                    
#Intermediates of the process
                      tmo.Chemical('Epoxy_stearic_acid', search_ID = '24560-98-3'),
                      tmo.Chemical('DHSA', search_ID = '120-87-6'), 
                      tmo.Chemical('Oxononanoic_acid',
                                   search_ID = '2553-17-5',
                                   search_db = False,
                                   phase='l',
                                   Tb = 304.8 + 25.0 + 273.15,
                                   Hf=-579480,
                                   formula = 'C9H16O3',
                                   MW = 172.22,
                                   CAS = '2553-17-5'),
                      
#Organic solvent for organic phase separation
                      tmo.Chemical('Ethyl_acetate'),  
#Products of the oxidative cleavage reaction
                      tmo.Chemical('Nonanal'),
                      tmo.Chemical('Nonanoic_acid'),
                      tmo.Chemical('Azelaic_acid'),
                      
                      tmo.Chemical('Methyl_oct_3', search_ID ='3-Methyloctane'),
                      tmo.Chemical('Methyl_oct_4', search_ID = '4-Methyloctane'),
                      tmo.Chemical('Dimethyl_heptane_3_3', search_ID ='3,3-Dimethylheptane'),
                      tmo.Chemical('Ethylheptane_3',search_ID ='3-Ethylheptane'),
                      tmo.Chemical('Ethylheptane_4', search_ID ='4-Ethylheptane'),
                      tmo.Chemical('bicyclo_octane',search_ID = '6221-55-2'),
                      tmo.Chemical('Hexane'),
                      tmo.Chemical('Octane'),
                      tmo.Chemical('pentene'),
                      tmo.Chemical('Methylcyclohexane'),
                      tmo.Chemical('Cyclopentane'),
                      tmo.Chemical('heptane'),
                      tmo.Chemical('pentane'),
                      tmo.Chemical('Sulphuric_acid'),
                      tmo.Chemical('Nonane'),
                      tmo.Chemical('Decane'),
                      tmo.Chemical('Undecane'),
                      tmo.Chemical('Cyclohexane'),
                      tmo.Chemical('Cycloheptane'),
                      tmo.Chemical('Benzene'),
                      tmo.Chemical('toluene'),
                                  ]) 

#Adding properties of oxononaoic acid, important for distillation unit process
Dortmund = chemicals_experimental.Oxononanoic_acid.Dortmund
Dortmund.set_group_counts_by_name(dict(COOH=1, CH2=7, CHO=1))
for chemical in chemicals_experimental: chemical.default()

chemicals_experimental.compile()
chemicals_experimental.set_synonym('Phosphotungstic_acid', 'fresh_Cat')
bst.settings.set_thermo(chemicals_experimental)

# for i in chemicals_process_1: 
#     if not i.locked_state: i.V.g.method_P = 'IDEAL'   
# return chemicals_process_1


    






