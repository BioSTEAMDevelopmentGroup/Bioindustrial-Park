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
    ### Chemicals that were missing some properties
    ### TODO.xxx check if this is a good assumption with Yoel
    
    #Tungstic_acid boiling point: https://en.wikipedia.org/wiki/Tungstic_acid
    tungsten,
    tmo.Chemical('tungstic_acid', Tb = 1746, phase = 's', 
                 Hvap=tungsten.Hvap, Psat=tungsten.Psat,
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
    ##Below lacks Hvap etc models
    # Sulfonated_polystyrene = chemical_database('Sulfonated_polystyrene',
    #                                            search_ID = '98-70-4',
    #                                            )
    ##Hence using polystyrene
    tmo.Chemical('Polystyrene'),
    ##Hence using polystyrene
    

])
chems.polystyrene_based_catalyst.copy_models_from(chems.Polystyrene,
                                        ['Hvap','Psat','sigma', 
                                         'epsilon', 'kappa', 'V',
                                         'Cn', 'mu'])

lipidcane_chems = lipidcane.chemicals.copy()
for chemical in lipidcane_chems: chemical.default()
chems.extend(lipidcane_chems)
chems.extend(sc.chemicals.copy())
chems.compile()
chems.set_synonym('OleicAcid', 'FFA')
chems.set_synonym('MonoOlein', 'MAG')
chems.set_synonym('DiOlein', 'DAG')
chems.set_synonym('TriOlein', 'TAG')
chems.define_group('Lipid', ('PL', 'FFA', 'MAG', 'DAG', 'TAG'))
chems.define_group('Oil', ('PL', 'FFA', 'MAG', 'DAG', 'TAG'))
chems.set_synonym('Water', 'H2O')
chems.set_synonym('Yeast', 'DryYeast')
chems.set_synonym('Biodiesel', 'Methyl_oleate')
bst.settings.set_thermo(chems)
chems.show()