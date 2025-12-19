# -*- coding: utf-8 -*-
# Copyright (C) 2025-, Sarang S. Bhagwat <sarangbhagwat.developer@gmail.com>
#
# References used in this module:
# [1] Smallbone et al., 2013. http://dx.doi.org/10.1016/j.febslet.2013.06.043
# [2] Palme et al., 1998. https://doi.org/10.1016/S0141-0229(96)00114-7
# [3] Schneidermann et al., 2014. https://dx.doi.org/10.1016/j.biortech.2014.11.087
# [4] Raj et al., 2014. https://doi.org/10.1021/bi5006442.
# [5] Cuny et al., 2022. https://doi.org/10.1038/s41467-022-30781-y

import nskinetics as nsk
import numpy as np

__all__ = ('rxn_sys')

MMC = nsk.MichaelisMentenCellular
Monod = nsk.Monod

#%%
chems = [
         'Yeast', 
         'Glucose', 
         'Ethanol',
         'Isobutanol',
         ]

MW_Yeast = ((10+100)/2) * 1e-12 * 6.02214e23 # [5]
MW_ADH = 145_000 # [4]
MW_Glucose = 180.156
MW_Ethanol = 46.06844

mol_ADH_per_mol_Yeast = 494_000 # [1]
MW_C5H7O2N_by_5 = 22.5931474
mol_ADH_per_mol_C5H7O2N_by_5 = mol_ADH_per_mol_Yeast * MW_C5H7O2N_by_5 / MW_Yeast

#%% Create a SpeciesSystem object with defined kinetic parameters
feed_concentrations = np.array([0.008946590564561485, 1.4697103500868698, 3.0390195636212934e-05, 0.0])
sp_sys = nsk.SpeciesSystem('sp_sys', 
                       chems,
                       concentrations=feed_concentrations)

# Describe reactions by writing chemical equations and kinetic parameter info
reactions = [
            Monod(C='Yeast', S='Glucose', 
                  stoich=1.0,
                  # umax=0.44/3600.0, # [2]
                  umax=0.3/3600.0, # [3]
                  KS=0.7/MW_Glucose, # 0.7 g/L from [3]
                  species_system=sp_sys,
                  inhibitors=['Ethanol', 'Isobutanol', 'Glucose'],
                  k_inhibs=[.058*MW_Ethanol, # [3]
                            .058*MW_Ethanol, 
                            .058*MW_Ethanol, # !!!
                            ],
                  
                  # k_inhibs=[3,
                  #           3, 
                  #           1,
                  #           ],
                  ),
            
            MMC(C='Yeast', S='Glucose', E_per_C=mol_ADH_per_mol_C5H7O2N_by_5, P='Ethanol',
                stoich=2.0,
                kcat=175, # [1]
                KM=0.462*1e-3, # [1]
                species_system=sp_sys),
            
            MMC(C='Yeast', S='Glucose', 
                E_per_C=mol_ADH_per_mol_C5H7O2N_by_5/5, # !!!
                P='Isobutanol',
                stoich=1.0,
                kcat=175,
                KM=0.462*1e-3,
                species_system=sp_sys),
            ]

# Generate a ReactionSystem from strings
rxn_sys = nsk.ReactionSystem(ID='rxn_sys', 
                                 reactions=reactions,
                                 species_system=sp_sys)

