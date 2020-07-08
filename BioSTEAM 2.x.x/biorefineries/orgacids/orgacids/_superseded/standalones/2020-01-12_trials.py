#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 09:21:08 2019

@author: yalinli_cabbi
"""

# %% Setup for cornstover

import biosteam as bst
from biorefineries.cornstover import system


# %% Setup for orgacids

import biosteam as bst
from orgacids import system
from orgacids.system import *

#system.orgacids_sys.diagram(kind='surface')
#system.orgacids_sys.diagram(kind='thorough')
system.orgacids_sys.simulate()

#system.feedstock_sys.diagram()
system.feedstock_sys.simulate()
U101 = bst.find.unit.U101

#system.pretreatment_sys.diagram()
system.pretreatment_sys.simulate()
M202 = bst.find.unit.M202
M203 = bst.find.unit.M203
R201 = bst.find.unit.R201


#system.fermentation_sys.diagram()
system.fermentation_sys.simulate()
H301 = bst.find.unit.H301
M301 = bst.find.unit.M301
R301 = bst.find.unit.R301

#system.separation_sys.diagram()
system.separation_sys.simulate()
M401 = bst.find.unit.M401

#system.lignin_sys.diagram()
system.lignin_sys.simulate()
S501 = bst.find.unit.S501

#system.wastewater_sys.diagram()
system.wastewater_sys.simulate()
M601 = bst.find.unit.M601
R601 = bst.find.unit.R601

CT = bst.find.unit.CT
BT = bst.find.unit.BT


# %% Random ones

target_acid = 'Ethanol'
water_carryover = 0.01
target_acid_recovery = 0.9
separation_split = np.ones(len(M401.outs[0].species))
water_index = M401.outs[0].index('Water')
target_acid_index = M401.outs[0].index(target_acid)
separation_split[water_index] = 1 - water_carryover
separation_split[target_acid_index] = 1 - target_acid_recovery
S401 = bst.units.Splitter('S401', ins=M401-0, outs=('', 'recovered_acid'),
                                split=separation_split)
S401.simulate()
S401.show()


# %% 

import biosteam as bst
from biosteam import Species, Stream, units
Stream.species = Species('Water', 'Ethanol', 'Methanol')
s1 = Stream('s1', Water=10, Ethanol=10, Methanol=3, T=340)

Stream.species = Species('Water', 'Ethanol')
s2 = Stream('s2', Ethanol=20)
M1 = units.Mixer('M1', ins=(s1, s2))

s3 = bst.Stream.like(s1, 's3')
T1 = bst.units.StorageTank('T1', ins=s3, outs=s1)

J1 = bst.Junction(upstream=s1, downstream=s2)
J1.simulate()
J1

import biosteam as bst
from biosteam import Species, Stream, units
Stream.species = Species('Water', 'Ethanol', 'Methanol')
s1 = Stream('s1', Water=10, Ethanol=10, Methanol=3, T=340)
s2 = Stream('s2', Water=20, Ethanol=10, Methanol=3, T=340)
s3 = Stream('s3', Water=30, Ethanol=10, Methanol=3, T=340)
s4 = Stream('s4', Water=30, Ethanol=10, Methanol=3, T=340)


M2 = bst.MolecularSieve('M2', ins=M1-0, outs=(a, ''), split=(0.1, 0.1),
                        order=('Ethanol', a)
                        )

# %% 


centriifuge = bst.units.Centrifuge_LLE(ID='centrifuge', ins=s1, outs=())


from biosteam import Stream
from biosteam.units import Centrifuge_LLE
from biorefineries.lipidcane import biodiesel_species

# Set up stream
Stream.species = biodiesel_species
feed = Stream(Lipid=1, Methanol=51, Glycerol= 9, Biodiesel=27, T=333.15)

# Set up centrifuge
C1 = Centrifuge_LLE('C1',
                    ins = feed,
                    outs = ('light', 'heavy'),
                    species_IDs=('Lipid', 'Methanol', 'Biodiesel'),
                    split=(1, 0.5, 1),
                    solvents=('Glycerol',),
                    solvent_split=(0.05,))

# Run all methods
C1.simulate()
# See all results
C1.diagram()
C1.show(T='degC', P='atm', fraction=True)

# %%

J1 = upstream=BT-1 - bst.Junction('J1') - Stream()

J701 = bst.Junction('J701', BT.outs[1], Stream())
J701-0-2-M601


import biosteam as bst
from biosteam import Species, Stream, units
Stream.species = Species('Water')
s1 = Stream('s1', Water=10)
Stream.species = Species('Water', 'Ethanol', 'Methanol')

J1 = upstream=s1 - bst.Junction('J1') - Stream()
J1.simulate()
J1

from biosteam import Stream
from biosteam.units import Centrifuge_LLE
from biorefineries.lipidcane import biodiesel_species

# Set up stream
Stream.species = biodiesel_species
feed = Stream(Lipid=1, Methanol=51, Glycerol= 9, Biodiesel=27, T=333.15)

# Set up centrifuge
C1 = Centrifuge_LLE('C1',
                    ins = feed,
                    outs = ('light', 'heavy'),
                    species_IDs=('Lipid', 'Methanol', 'Biodiesel'),
                    split=(1, 0.5, 1),
                    solvents=('Glycerol',),
                    solvent_split=(0.05,))

# Run all methods
C1.simulate()

# See all results
C1.diagram()
C1.show(T='degC', P='atm', fraction=True)

for i, x in U101.cost_items.items():
    S = D[x._basis]
    q = S/x.S
    C[i] = N*bst.CE/x.CE*x.cost*q**x.n
    print(C[i])







