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

from orgacids.system import *
feedstock_sys.simulate()
pretreatment_sys.simulate()
fermentation_sys.simulate()
separation_sys.simulate()
wastewater_sys.simulate()
facilities_sys.simulate()


# %%

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


# %% Cutout codes

'''
# Placeholder copying the AmmoniaAdditionTank unit from cornstover biorefinery
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=410369, ub=False, CE=521.9, cost=236000*10, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=410369, ub=False, CE=521.9, cost=21900*2, n=0.5, kW=7.457, BM=1.5)
class SeparationChemicalsAdditionTank(bst.Unit): pass

# Placeholder copying the BeerTank unit from cornstover biorefinery
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=425878, ub=False, CE=521.9, cost=636000*10, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=425878, ub=False, CE=521.9, cost=26800*10, n=0.8, kW=93.2125, BM=2.3)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=425878, ub=False, CE=521.9, cost=68300*10, n=0.5, kW=14.914, BM=1.5)
class SeparationSplitter(bst.units.Splitter): pass

# Placeholder copying the heat exchangers of the SaccharificationAndCoFermentation
# unit from cornstover biorefinery
@cost(basis='Duty', ID='Heat exchangers', units='Gcal/hr',
      S=5*_Gcal2kJ, ub=False, CE=521.9, cost=23900*10, n=0.7, BM=2.2)
class SeparationHX(bst.units.HXutility): pass
'''

'''
# Placeholder for separation process
# outs[0] is the water-rich phase and outs[1] is the organic acid
# Assuming 90% of target acid can be recovered,
# 1% of water is carried over with acid,
# and no other impurities exist in the acid product
water_carryover = 0.01
target_acid_recovery = 0.9
separation_split = np.ones(len(M401.outs[0].chemicals))
water_index = M401.outs[0].chemicals.index('Water')
target_acid_index = M401.outs[0].chemicals.index(target_acid)
separation_split[water_index] = 1 - water_carryover
separation_split[target_acid_index] = 1 - target_acid_recovery
S401 = units.SeparationSplitter('S401', ins=M401-0, outs=('', 'raw_acid'),
                                split=separation_split)
M402 = bst.Mixer('M402', ins=(S401-1, ''))
# Placeholder for heat exchagne during separation,
# temperature followed cornstover biorefienry design
H401 = units.SeparationHX('H401', ins=M402-0, T=115+298.15)
# Assuming 95% water is removed and 80% of acid is recovered in outs[0]
# Note that acid removed with water will to recycled
S402 = bst.MolecularSieve('S402', ins=H401-0, outs=('product_stream_for_storage', ''),
                          #??? check Yoel's number in the origianl script
                          split=(0.05, 0.8),
                          order=('Water', target_acid)
                          )
S402-1-1-M402

broth_for_separation = M401-0
separation_chemicals_over_target_acid = 1 # assuming 1:1 addition of separation chemicals
separation_chemicals_mass = separation_chemicals.imass['SeparationChemicals']
def update_separation_chemicals_loading():
    target_acid_mass = broth_for_separation.imass[target_acid]
    separation_chemicals_mass = target_acid_mass * separation_chemicals_over_target_acid
    separation_chemicals.imass['SeparationChemicals'] = separation_chemicals_mass

'''

'''
J501 = bst.Junction('J501', S401-0, 0-S501)

lignin_sys = System('lignin_sys',
                    network=(J501, S501)
                    )
Area500 = lignin_sys
'''