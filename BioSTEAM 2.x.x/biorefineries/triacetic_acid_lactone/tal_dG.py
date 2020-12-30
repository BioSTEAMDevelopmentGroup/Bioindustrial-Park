# -*- coding: utf-8 -*-
"""
Estimation of heat of fusion of triacetic acid lactone (TAL)
by fitting solubility across temperature.

"""
import numpy as np
from thermosteam.constants import R
from scipy.optimize import curve_fit
import flexsolve as flx
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

# %% Solubility model

def _solubility(T, dG, gamma):
    RT = (R*T)
    return np.exp(-dG/RT)/gamma

class SLE(eq.SLE):
    __slots__ = ('dG',)
    def _solve_x(self, T):
        x = _solubility(T, self.dG, 1.) # Initial guess
        if isinstance(self._gamma, eq.IdealActivityCoefficients):
            return _solubility(T, self.dG, self.activity_coefficient or 1.)
        return flx.aitken(self._x_iter, x, xtol=1e-6, args=(T, self.dG), checkiter=False, maxiter=100)
        
    def _x_iter(self, x, T, dG):
        self._update_solubility(x)
        liquid_mol = self._liquid_mol[self._index]
        F_mol_liquid = liquid_mol.sum()
        x_l = liquid_mol / F_mol_liquid
        gamma = self._gamma(x_l, T)
        return _solubility(T, dG, gamma[self._solute_gamma_index])
# SLE = eq.SLE

# %% Data

temperatures = np.array([0, 22, 42, 65], float) # C
#temperatures = temperatures[:-1]
temperatures += 273.15 # K
TAL_g = np.array([3.52, 8.92, 16.46, 21.13], float) # g/L
#TAL_g = TAL_g[:-1]
TAL_mol = TAL_g / 126.11004
Water_mol = 55.5 # mol / L
solubilities = TAL_mol / (Water_mol + TAL_mol) # by wt.
RTlnx = R * temperatures * np.log(solubilities)
Phenol = Chemical('Phenol').Cn.l[0]
TAL = Chemical('TAL',
               search_ID='Triacetic acid lactone')
Water = Chemical('WAter')
Water.NIST.set_group_counts_by_name({'H2O':1})
TAL.NIST.set_group_counts_by_name({'CH3':1, 'c-CH=C':2, 'c-CO-O':1, 'OH tert':1})
TAL.Tm = 185 + 273.15
rho_TAL = 1.348e-3
molar_volume_TAL = fn.rho_to_V(rho_TAL, TAL.MW)
TAL.V.l.add_model(molar_volume_TAL, top_priority=True)
TAL.V.s.add_model(molar_volume_TAL, top_priority=True)
Cn = 2 * TAL.MW
TAL.Cn.l.add_model(Cn, top_priority=True)
TAL.Cn.s.add_model(Cn, top_priority=True)
Water, TAL = chemicals = Chemicals([Water, TAL])
chemicals.compile(skip_checks=True)
thermo = Thermo(chemicals, Gamma=eq.NISTActivityCoefficients)
settings.set_thermo(thermo)

# %% Activity coefficients for excel fitting

Gamma = eq.DortmundActivityCoefficients(chemicals)
gamma = [Gamma([1-x, x], T)[1] for x, T in zip(solubilities, temperatures)]


# %% Solubility fitting

# sle = SLE()
# sle.imol['l', 'Water'] = 100
# sle.imol['l', 'TAL'] = 10

# @np.vectorize
# def RTlnsolubility(T, dG): # by wt
#     sle.dG = dG
#     sle('TAL', T=T)
#     return  sle.imol['l', 'Water'] / sle.imol['l', ('Water', 'TAL')].sum()

# params, _ = curve_fit(RTlnsolubility, temperatures, RTlnx,
#                       p0=[-40], bounds=([-np.inf], [-0.1]))

# def display_solubility_results(dG):
#     plt.figure()
#     plt.plot(temperatures, RTlnsolubility(temperatures, dG), label='fit')
#     plt.scatter(temperatures, RTlnx, label='actual')
#     plt.show()
    
#     r2 = r2_score(RTlnx, RTlnsolubility(temperatures, *params))
#     print('R2:', r2)
    
# display_solubility_results(*params)
