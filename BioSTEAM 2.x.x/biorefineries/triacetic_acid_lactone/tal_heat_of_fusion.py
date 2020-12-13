# -*- coding: utf-8 -*-
"""
Estimation of heat of fusion of triacetic acid lactone (TAL)
by fitting solubility across temperature.

"""
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
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
solubilities = np.array([3.52, 8.92, 16.46, 21.13], float) # g/L

Phenol = Chemical('Phenol').Cn.l[0]
TAL = Chemical('TAL',
               search_ID='Triacetic acid lactone')
TAL.Tm = 185 + 273.15
rho_TAL = 1.348e-3
molar_volume_TAL = fn.rho_to_V(rho_TAL, TAL.MW)
TAL.V.l.add_model(molar_volume_TAL, top_priority=True)
TAL.V.s.add_model(molar_volume_TAL, top_priority=True)
Cn = 2 * TAL.MW
TAL.Cn.l.add_model(Cn, top_priority=True)
TAL.Cn.s.add_model(Cn, top_priority=True)
Water, TAL = chemicals = Chemicals(['Water', TAL])
chemicals.compile(skip_checks=True)
thermo = Thermo(chemicals, Gamma=eq.DortmundActivityCoefficients)
settings.set_thermo(thermo)
sle = eq.SLE()
sle.imol['l', 'Water'] = 100
sle.imol['l', 'TAL'] = 10

@np.vectorize
def solubility(T, Hfus, Cns): # g/L
    TAL.Cn.s[0].value = Cns
    TAL.Hfus = Hfus
    # sle.activity_coefficient = gamma
    sle('TAL', T=T)
    molar_volume = Water.V.l(T, 101325)
    tal_wt = TAL.MW * sle.imol['l', 'TAL']
    tal_vol = tal_wt / rho_TAL
    water_vol = molar_volume * sle.imol['l', 'Water']
    return  tal_wt / (1000 * water_vol)

Cn_lb = Cn/4
Cn_ub = 4*Cn
Hfus_lb = 0.
Hfus_ub = 2*95.15
gamma_lb = 0.000000001
gamma_ub = 20
params, _ = curve_fit(solubility, temperatures, solubilities,
                      p0=[40, Cn],
                      bounds=([Hfus_lb, Cn_lb], [Hfus_ub, Cn_ub])
)

def display_solubility_results(Hfus, gamma):
    plt.figure()
    plt.plot(temperatures, solubility(temperatures, Hfus, gamma), label='fit')
    plt.scatter(temperatures, solubilities, label='actual')
    plt.show()
    
    r2 = r2_score(solubilities, solubility(temperatures, *params))
    print('R2:', r2)
    
display_solubility_results(*params)