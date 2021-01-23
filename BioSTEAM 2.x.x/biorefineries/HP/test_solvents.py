# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 19:09:17 2021

@author: saran
"""
# %% Imports and chemicals initialization

import numpy as np
import thermosteam as tmo
import biosteam as bst

from biorefineries.HP.chemicals_data import HP_chemicals

tmo.settings.set_thermo(['Water', 'octanol', 'hexanol', 'butyl acetate', HP_chemicals['Xylose'], HP_chemicals['Glucose'], HP_chemicals['HP'], 'isoamyl alcohol'])

Water = HP_chemicals['Water']
Glucose = HP_chemicals['Glucose']
HP = HP_chemicals['HP']
AQ336 = HP_chemicals['AQ336']
Octanol = HP_chemicals['Octanol']
Hexanol = tmo.Chemical('Hexanol')
Butyl_acetate = tmo.Chemical('Butyl acetate')
Propyl_acetate = tmo.Chemical('Propyl acetate')
Isoamyl_alcohol = tmo.Chemical('Isoamyl alcohol')
TOA = tmo.Chemical('Trioctylamine')
Dodecanol = tmo.Chemical('Dodecanol')
Nonanol = tmo.Chemical('Nonanol')
te_hexanol = tmo.Chemical('2-Ethyl hexanol')
Cyclohexanol = tmo.Chemical('Cyclohexanol')
Cyclohexanone = tmo.Chemical('Cyclohexanone')
Dioctyl_phthalate = tmo.Chemical('117-81-7')
Diethyl_sebacate = tmo.Chemical('Diethyl sebacate')
Diethyl_sebacate.copy_models_from(Water, ['Psat', 'Hvap'])
Octanediol = tmo.Chemical('1,8-Octanediol')
H2SO4 = tmo.Chemical('H2SO4')
tmo.settings.set_thermo(['Water', 'H2SO4', 'Cyclohexanol', Octanediol, Diethyl_sebacate, '2-Ethyl hexanol', '117-81-7', 'Cyclohexanone', 'Hexanol', 'Butyl acetate', 'Isoamyl alcohol', 'Dodecanol', 'Trioctylamine', HP_chemicals['Xylose'], HP_chemicals['Glucose'], HP_chemicals['HP']], HP_chemicals['AQ336'])
# %% Streams initialization

T = 320
process_stream = tmo.Stream('process_stream',
                            Water = 4000, HP = 250,
                            units = 'kmol/hr',
                            T = T)

#
solvent_chemical = Octanediol # set this
#

solvent_ID = solvent_chemical.ID
solvent_stream = tmo.Stream('solvent_stream',
                            T = T)
solvent_stream.imol[solvent_ID] = 2000
# solvent_stream.imol['H2SO4'] = 200
mixed_stream = tmo.Stream()

# %%% Functions

def get_K(chem_ID, stream, phase_1, phase_2):
    return (stream[phase_1].imol[chem_ID]/stream[phase_1].F_mol)/(stream[phase_2].imol[chem_ID]/stream[phase_2].F_mol)

# %% Stream-only tests
mixed_stream.empty()
mixed_stream.mix_from([process_stream, solvent_stream])
mixed_stream.lle(T=T, top_chemical = solvent_ID)
mixed_stream.show(N=100, composition=True)

extract_phase = 'l'
raffinate_phase = 'L'
K_HP_in_solvent = get_K('HP', mixed_stream, extract_phase, raffinate_phase)
# K_Glucose_in_solvent = get_K('Glucose', mixed_stream, extract_phase, raffinate_phase)
K_Water_in_solvent = get_K('Water', mixed_stream, extract_phase, raffinate_phase)
K_solvent_in_Water = get_K(solvent_ID, mixed_stream, raffinate_phase, extract_phase)
# K_Water_in_Water = get_K('Water', mixed_stream, raffinate_phase, extract_phase)

print(K_HP_in_solvent)
print(K_Water_in_solvent, K_solvent_in_Water)

# %% Unit initialization and tests

partition_data = dict(IDs=('HP', 'Water', 'Isoamyl alcohol'),
           K=np.array([1./K_HP_in_solvent, 1/K_Water_in_solvent, K_solvent_in_Water]),
           phi = 0.5)

MS = bst.units.MultiStageMixerSettlers('MS', ins = (process_stream, solvent_stream),
                                     outs = ('raffinate', 'extract'),
                                     N_stages = 5, partition_data = partition_data,) 

MS.simulate()
MS.show(N=100, composition=True)