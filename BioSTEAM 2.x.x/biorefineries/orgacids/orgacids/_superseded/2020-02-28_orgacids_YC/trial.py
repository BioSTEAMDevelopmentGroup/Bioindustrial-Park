#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""

# %% Set up

import numpy as np
import biosteam as bst
import thermosteam as tmo

from orgacids.chemicals import orgacids_chemicals, chemical_groups
tmo.settings.set_thermo(orgacids_chemicals)



#%%
from orgacids.system import *
feedstock_sys.simulate()
pretreatment_sys.simulate()
fermentation_sys.simulate()





update_stripping_water_loading()
U401.simulate()
M401.simulate()
S401.simulate()
R401.simulate()
update_separation_sulfuric_acid_loading()
T401.simulate()
M402.simulate()

S402.simulate()
print('S402.outs[0] is')
S402.outs[0].show(N=100)


S403.simulate()
print('S403.outs[0] is')
S403.outs[0].show(N=100)
print('S403.outs[1] is')
S403.outs[1].show(N=100)


#%%
R402.simulate()
print('R402.outs[0] is')
R402.outs[0].show(N=100)

S404.simulate()
print('S404.outs[0] is')
S404.outs[0].show(N=100)
print('S404.outs[1] is')
S404.outs[1].show(N=100)


#%%
H401.simulate()

R403.simulate()
print('R403.outs[0] is')
R403.outs[0].show(N=100)

S405.simulate()
print('S405.ins[0] is')
S405.ins[0].show(N=100)
print('S405.outs[0] is')
S405.outs[0].show(N=100)
print('S405.outs[1] is')
S405.outs[1].show(N=100)

M403.simulate()
H402.simulate()

#%%



separation_sys.simulate()
wastewater_sys.simulate()
facilities_sys.simulate()



#%%

import thermosteam as tmo
glycerol, water = tmo.Chemicals(['Glycerol', 'Water'])
water.show()
glycerol.show()


# %% Mass balance

from biosteam import System
from biosteam.units import Mixer, Splitter, StorageTank, Pump, Flash, MassBalance
from thermosteam import Chemicals, Stream, settings
chemicals = Chemicals(['Water', 'Ethanol'])
settings.set_thermo(chemicals)
water = Stream('water',
                Water=40,
                units='lb/s',
                T=350, P=101325)
ethanol = Stream('ethanol',
                 Ethanol=190, Water=30,
                 T=300, P=101325)
target = Stream('target', flow=[50, 50])
T1 = StorageTank('T1')
T2 = StorageTank('T2')
P1 = Pump('P1', P=101325)
P2 = Pump('P2', P=101325)
M1 = Mixer('M1', outs='s1')
S1 = Splitter('S1', outs=('s2', 's3'), split=0.5)
F1 = Flash('F1', outs=('s4', 's5'), V=0.5, P =101325)
# Connect units
water-T1-P1
ethanol-T2-P2
[P1-0, P2-0, S1-0]-M1

MB1 = MassBalance('MB1', streams=[0,1],
                  chemical_IDs=['Ethanol', 'Water'],
                  outs=target,
                  ins=(water, ethanol, S1-0))
mixSys = System('mixSys',
                 recycle=S1-0,
                 network=(MB1, T1, P1, T2, P2, M1, F1, S1))
# Make diagram to view system
mixSys.diagram()
mixSys.simulate()
MB1.show()

