"""
Created on Sat Aug 20 21:47:53 2022

@author: Lavanya_Kudli
"""
import biosteam as bst
import biosteam.units
import thermosteam as tmo
# from biorefineries.lipidcane import create_lipid_pretreatment_system, create_transesterification_and_biodiesel_separation_system
# from biorefineries.lipidcane import chemicals
# bst.settings.set_thermo(chemicals)

#Figure the following out
# vegetable_oil = bst.Stream(ID='vegetable_oil',
#                            Water=0.0184,
#                            TAG=11.1)
# biodiesel = bst.Stream(ID = 'biodiesel')
# crude_glycerol = bst.Stream(ID = 'crude_glycerol')
# wastewater = bst.Stream(ID = 'wastewater')
# create_transesterification_and_biodiesel_separation_system(
#         ins= vegetable_oil, 
#         outs=[biodiesel, crude_glycerol,wastewater],)
# create_transesterification_and_biodiesel_separation_system().simulate()

# chems is the object containing all chemicals used in this biorefinery
chems = tmo.Chemicals([])

# To keep track of which chemicals are available in the database and which
# are created from scratch
database_chemicals_dict = {}
copied_chemicals_dict = {}
defined_chemicals_dict = {}

def chemical_database(ID, phase=None, **kwargs):
    chemical = tmo.Chemical(ID, **kwargs)
    if phase:
        chemical.at_state(phase)
        chemical.phase_ref = phase
    chems.append(chemical)
    database_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

def chemical_copied(ID, ref_chemical, **data):
    chemical = ref_chemical.copy(ID)
    chems.append(chemical)
    for i, j in data.items(): setattr(chemical, i, j)
    copied_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

def chemical_defined(ID, **kwargs):
    chemical = tmo.Chemical.blank(ID, **kwargs)
    chems.append(chemical)
    defined_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

##Chemicals that are already in the data base

Hydrogen_peroxide = chemical_database('Hydrogen_peroxide',phase = 'l')

Biodiesel = chemical_database('Methyl_oleate')
Water =  chemical_database('Water')
#look into phase of the below
MDHSA = chemical_database('MDHSA',search_ID = '1115-01-1')
Pelargonic_acid = chemical_database('Pelargonic_acid')
Azelaic_acid = chemical_database('Azelaic_acid')
Monomethyl_azelate = chemical_database('Monomethyl_azelate')
Suberic_acid = chemical_database('Suberic_acid')
Caprylic_acid = chemical_database('Caprylic_acid')
Nitrogen = chemical_database('Nitrogen')
Oxygen = chemical_database('Oxygen')
Methanol = chemical_database('Methanol')
### Chemicals that were missing some properties
### TODO.xxx check if this is a good assumption with Yoel
Tungsten = chemical_database('tungsten')
#Tungstic_acid boiling point: https://en.wikipedia.org/wiki/Tungstic_acid
Tungsten_catalyst = chemical_database('tungstic_acid', Tb = 1746, phase = 's')
Tungsten_catalyst.copy_models_from(Tungsten,['Hvap','Psat'])

###Using cobalt chloride instead of acetate as allowed by the patent
###cobalt acetate has a few missing properties, further GWP data not available
Cobalt_chloride = chemical_database('Cobalt_chloride',search_ID = '7646-79-9',
                                    phase = 's')
# cobalt_acetate_Tb: https://www.chemsrc.com/en/cas/71-48-7_34110.html
Cobalt_catalyst = chemical_database('Cobalt_acetate',Tb = 117.1+273.15, 
                                    phase = 's')
Cobalt_catalyst.copy_models_from(Cobalt_chloride,['Hvap','Psat'])
#defaulting the rest of the properties to that of water
Cobalt_catalyst.default()

###Modelling amberlyte catalyst like a solid catalyst
##Using sunfonated_polystyrene boiling point
##https://www.chemsrc.com/en/cas/39389-20-3_843683.html#:~:text=amberlyst%28r%29%2015%20CAS%20Number%3A%2039389-20-3%3A%20Molecular%20Weight%3A%20314.39900%3A,Point%3A%20266.3%C2%BAC%3A%20Symbol%3A%20GHS07%3A%20Signal%20Word%3A%20Warning%20%C3%97

##Chemicals not in the database
Amberlyst_catalyst = chemical_defined('polystyrene_based_catalyst',
                                       Tb = 516.7+273.15,
                                       phase = 's')
##Below lacks Hvap etc models
# Sulfonated_polystyrene = chemical_database('Sulfonated_polystyrene',
#                                            search_ID = '98-70-4',
#                                            )
##Hence using polystyrene
Polystyrene = chemical_database('Polystyrene')
##Hence using polystyrene
Amberlyst_catalyst.copy_models_from(Polystyrene,
                                    ['Hvap','Psat','sigma', 
                                     'epsilon', 'kappa', 'V',
                                     'Cn', 'mu'])
chems.compile()
tmo.settings.set_thermo(chems)
chems.show()