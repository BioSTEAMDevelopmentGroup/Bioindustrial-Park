# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 15:20:08 2023

@author: sarangbhagwat
"""

data = {'lab scale':{'batch': {}, 'fed-batch': {}}, 'pilot scale':{'batch': {}, 'fed-batch': {}}}

data['lab scale']['batch'] = {
                                'yield': 0.473347, # g-succinic-acid-eq. / g-sugars-glucose-eq.
                                'titer': 63.1, # g-succinic-acid-eq. / L-effluent-broth
                                'productivity': 0.657292, # g-succinic-acid-eq. / L-effluent-broth / h-fermentation-time
                                'microbe_yield': 0.152281, # g-cell-mass / g-sugars-glucose-eq.
                                'inoculum_ratio': 0.04, # % of feed to fermentation tank
                                'CO2_safety_factor': 1.00153, # mol-CO2-req / mol-CO2-stoichiometry
                                'mol_air_per_mol_succinic_acid': 5.30096,
                                'mol_base_per_mol_succinic_acid': {'lime': 0.298866/2., 'NH4OH': 0.298866,},
                                'g_diammonium_sulfate_per_g_succinic_acid': 4.5/63.1,
                                'g_magnesium_sulfate_per_g_succinic_acid': 0.9/63.1,
                                'T': 30. + 273.15, # degrees Kelvin
                                'agitation': 173.4, # RPM
                                'working_volume_fraction': 0.4,
                             }
    
data['pilot scale']['batch'] = {
                                'yield': 0.473347, # g-succinic-acid-eq. / g-sugars-glucose-eq.
                                'titer': 63.1, # g-succinic-acid-eq. / L-effluent-broth
                                'productivity': 0.657292, # g-succinic-acid-eq. / L-effluent-broth / h-fermentation-time
                                'microbe_yield': 0.152281, # g-cell-mass / g-sugars-glucose-eq.
                                'inoculum_ratio': 0.04, # % of feed to fermentation tank
                                'CO2_safety_factor': 1.00153, # mol-CO2-req / mol-CO2-stoichiometry
                                'mol_air_per_mol_succinic_acid': 5.30096,
                                'mol_base_per_mol_succinic_acid': {'lime': 0.298866/2., 'NH4OH': 0.298866,},
                                'g_diammonium_sulfate_per_g_succinic_acid': 4.5/63.1,
                                'g_magnesium_sulfate_per_g_succinic_acid': 0.9/63.1,
                                'T': 30. + 273.15, # degrees Kelvin
                                'agitation': 173.4, # RPM
                                'working_volume_fraction': 0.4,
                             }
