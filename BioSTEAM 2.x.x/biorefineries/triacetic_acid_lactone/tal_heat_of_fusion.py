# -*- coding: utf-8 -*-
"""
Estimation of heat of fusion of triacetic acid lactone (TAL)
by fitting solubility across temperature.

"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from thermosteam import (
    functional as fn,
    Thermo,
    Chemical,
    Chemicals,
    equilibrium as eq, 
    settings,
)
temperatures = np.array([0, 22, 42, 65], float) # C
temperatures += 273.15 # K
solubilities = np.array([3.52, 8.92, 16.46, 21.13], float)

Phenol = Chemical('Phenol').Cn.l[0]
TAL = Chemical('TAL',
               search_ID='Triacetic acid lactone')
TAL.Tm = 185 + 273.15
molar_volume_TAL = fn.rho_to_V(1.348e-3, TAL.MW)
TAL.V.l.add_model(molar_volume_TAL, top_priority=True)
TAL.V.s.add_model(molar_volume_TAL, top_priority=True)
TAL.Cn.l.add_model(2 * TAL.MW, top_priority=True)
TAL.Cn.s.add_model(2 * TAL.MW, top_priority=True)
Water, TAL = chemicals = Chemicals(['Water', TAL])
chemicals.compile(skip_checks=True)
thermo = Thermo(chemicals)
settings.set_thermo(thermo)
sle = eq.SLE()
sle.imol['l', 'Water'] = 100
sle.imol['l', 'TAL'] = 10

@np.vectorize
def solubility(T, Hfus): # g/L
    # TAL.Cn.s[0].value = Cns
    TAL.Hfus = Hfus
    sle('TAL', T=T)
    molar_volume = Water.V('l', T, 101325)
    return  TAL.MW * sle.imol['l', 'TAL'] / (molar_volume * sle.imol['l', 'Water'])

params, _ = curve_fit(solubility, temperatures, solubilities, p0=[1000])
plt.plot(temperatures, solubility(temperatures, *params), label='fit')
plt.scatter(temperatures, solubilities, label='actual')
plt.show()