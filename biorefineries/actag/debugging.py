# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 14:31:22 2022

@author: yrc2
"""
import numpy as np
from biorefineries import actag as ac
import matplotlib.pyplot as plt
from warnings import filterwarnings
import biosteam as bst
filterwarnings('ignore')
ac.load(3)

titers = np.linspace(20, 40, 20)
ac.fermentation.product_yield = 0.6
MPSPs = []
TCIs = []
utility_costs = []
evaporators = []
dilution_water = []
MEE = ac.sys.flowsheet(bst.MultiEffectEvaporator)
# ac.HXN901.units = []
titer_err = []
for titer in titers:
    ac.fermentation.titer = titer
    ac.sys.simulate()
    MPSP = ac.MPSP()
    TCIs.append(ac.tea.TCI / 1e6)
    utility_costs.append(ac.sys.utility_cost / 1e6)
    evaporators.append(len(MEE.P))
    dilution_water.append(ac.sys.flowsheet.stream.dilution_water.F_mass / 1e3)
    MPSPs.append(MPSP)
    titer_err.append(ac.fermentation.titer - ac.fermentation.get_titer())

fig, axs = plt.subplots(2, 2)
plt.sca(axs[0, 0])
plt.ylabel('MPSP')
plt.plot(titers, MPSPs, label='MPSP', color='g')

plt.sca(axs[0, 1])
plt.ylabel('TCI')
plt.plot(titers, TCIs, label='TCIs', color='b')

plt.sca(axs[1, 0])
plt.ylabel('titer_err')
plt.plot(titers, titer_err, label='titer_err', color='r')

# plt.sca(axs[1, 1])
# plt.ylabel('Evaporators')
# plt.plot(titers, evaporators, label='evaporators', color='orange')

plt.sca(axs[1, 1])
plt.ylabel('dilution_water')
plt.plot(titers, dilution_water, label='dilution_water', color='k')



    