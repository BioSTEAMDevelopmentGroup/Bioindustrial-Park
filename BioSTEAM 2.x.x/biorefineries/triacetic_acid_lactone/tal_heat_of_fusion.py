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

# def solubility(T, dG, gamma):
#     RT = (R*T)
#     return np.log(-dG/RT)/gamma

# class SLE(SLE):
#     __slots__ = ('dG',)
#     def _solve_x(self, T):
#         solute_chemical = self.chemicals.tuple[self._solute_index]
#         gamma = 1.
#         x = solubility(T, dG, gamma) # Initial guess
#         args = (T, dG)
#         if isinstance(self._gamma, IdealActivityCoefficients):
#             return solubility(T, dG, self.activity_coefficient or 1.)
#         return flx.aitken(self._x_iter, x, xtol=1e-6, args=args, checkiter=False, maxiter=100)
        
#     def _x_iter(self, x, T, dG):
#         self._update_solubility(x)
#         liquid_mol = self._liquid_mol[self._index]
#         F_mol_liquid = liquid_mol.sum()
#         x_l = liquid_mol / F_mol_liquid
#         gamma = self._gamma(x_l, T)
#         return solubility(T, dG, gamma[self._solute_gamma_index])

temperatures = np.array([0, 22, 42, 65], float) # C
temperatures += 273.15 # K
solubilities = np.array([3.52, 8.92, 16.46, 21.13], float) # g/L

Phenol = Chemical('Phenol').Cn.l[0]
TAL = Chemical('TAL',
               search_ID='Triacetic acid lactone')
TAL.Tm = 185 + 273.15
TAL.Hfus = 15425
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
def solubility(T): # g/L
    sle('TAL', T=T)
    molar_volume = Water.V.l(T, 101325)
    tal_wt = TAL.MW * sle.imol['l', 'TAL']
    tal_vol = tal_wt / rho_TAL
    water_vol = molar_volume * sle.imol['l', 'Water']
    return  tal_wt / (1000 * water_vol)

# Cn_lb = -np.inf
# Cn_ub = np.inf
# Hfus_lb = 0
# Hfus_ub = np.inf
# params, _ = curve_fit(solubility, temperatures, solubilities,
#                       p0=[40, Cn], bounds=([Hfus_lb, Cn_lb], [Hfus_ub, Cn_ub]))


def display_solubility_results(*params):
    plt.figure()
    plt.plot(temperatures, solubility(temperatures), label='fit')
    plt.scatter(temperatures, solubilities, label='actual')
    plt.show()
    
    r2 = r2_score(solubilities, solubility(temperatures))
    print('R2:', r2)
    
    
display_solubility_results()
