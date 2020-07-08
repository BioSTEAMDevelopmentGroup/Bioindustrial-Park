#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 16:54:56 2019

Fermentation of lactic acid using environmental biotechnology (EB) model.

Structure of EB model based on:
    (1) Tchobanoglous, G.; Burton, F. L.; Stensel, H. D.; Inc, M. & E.; Burton, F. Wastewater Engineering: Treatment and Reuse; McGraw-Hill Education, 2003.
    (2) Rittmann, B. E.; McCarty, P. L. Environmental Biotechnology: Principles and Applications; McGraw-Hill, 2001.

    as:
        dX/dt = μ_max * X * (S/(K_S+S)) * (1-(P/P_max)^n)
        dS/dt = -dX/dt * (1/Y_X/S)
        dP/dt = -dX/dt * (Y_P/S/Y_X/S)
    
    where:
        X is biomasss
        S is substrate
        P is product, P_max is maximum product titer
        μ is specific cell growth rate, μ = dX/dt/X, μ_max is the maximum specific cell growth rate
        Y_X/S is apparent yield of biomass per substrate consumed, Y_X/S = -dX/dS
        Y_P/S is apparent yield of product per substrate consumed, Y_P/S = -dP/dS
        n is toxic power
        K_S is half-velocity constnat, the value of S when μ/μ_max = 0.5

@author: yalinli_cabbi
"""

#%% import modules

import numpy as np
#import pandas as pd
from scipy.integrate import odeint
from biosteam import Stream, Unit
from biosteam.reaction import Reaction
from biosteam.units._hx import HXutility
from biosteam.units._tank import MixTank
from biosteam.units.decorators import cost
from biosteam.units.designtools import size_batch


#%% set unit costs

@cost('Reactor volume', 'Clearning in place', CE=521.9, cost=421e3, S=3785,
      n=0.6, BM=1.8, N='N')

@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500, S=3785,
      n=0.5, kW=22.371, BM=1.5, N='N')

@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000, S=3785,
      n=0.5, BM=1.5, N='N')


#%% fermentation

class Ferm_EB_LA(Unit):
    
    """
    Create a Fermentation object which models large-scale batch fermentation of lactic acid, kinetic parameters were based on [REFS]    
    A compound with CAS 'Yeast' must be present in species.    
    Sucrose (will be hydrolyzed into glucose), glucose, and xylose are taken into account for conversion    
    Conversion is based on reaction time, `tau`    
    Cleaning and unloading time, `tau_0`, fraction of working volume, `V_wf`, and number of reactors, `N_reactors`, are attributes that can be changed    
    Cost of a reactor is based on the NREL batch fermentation tank cost assuming volumetric scaling with a 6/10th exponent 
    
    Parameters
    ----------
    ins
        [:] Inffluent streams
        
    outs
        [0] Effluent
        
        [1] CO2
        
    tau: reaction time, float    
    eff_g: % of glucose that can be converted to lactic, optional, default is 0.9    
    eff_x: % of xylose that can be converted to lactic, optional, default is 0.9
    iskinetic_g: bool, optional, default is False; if True, 'Fermentation.kinetic_mode' will be used
    iskinetic_x: bool, optional, default is False; if True, 'Fermentation.kinetic_mode' will be used        
    N: int, number of batch reactions
    
    References
    ----------
    TO BE ADDED
    """
       
    _units = {'Reactor volume': 'm3',
              'Cycle time': 'hr',
              'Loading time': 'hr',
              'Total dead time': 'hr'}
    
    _N_ins = _N_outs = 2
    _N_heat_utilities = 0
    _has_power_utility = True
    line = 'Fermentation'
    
    # If True, number of reactions (N) will be calculated to minimize installation cost, otherwise N is constant
    autoselect_N = False
    
    # Clearning and unloading time (hr)
    tau_0 = 3
    
    # Fraction of filled tank to total tank volume
    working_volume_fraction = MixTank.working_volume_fraction
    _V_wf = 0.9
    
    # Fitted kinetic parameters, g for glucose and x for xylose
    n = 3
    P_max = 100 # kg/m3
    Y_XS_g = 0.08
    Y_XS_x = 0.11
    Y_PS_g = 0.33
    Y_PS_x = 0.66
    mu_max_g = 0.21
    mu_max_x = 0.087
    K_S_g = 30
    K_S_x = 7.39
    kinetic_constants = (n, P_max, Y_XS_g, Y_XS_x, Y_PS_g, Y_PS_x,
                         mu_max_g, mu_max_x. K_S_g, K_S_x)
    

    def _more_design_specs(self):
        return(
            ('Cleaning and unloading time', self.tau_0, 'hr'),
            ('Working volume fraction', self.working_volume_fraction, ''),
            ('Number of reactors', self.N, '')
               )
    
    def __init__(self, ID='', ins=None, outs=(), *, tau, N, 
                 eff_g=0.9, eff_x=0.9, iskinetic_g=False, iskinetic_x=False):
        Unit.__init__(self, ID, ins, outs)
        self.hydrolysis = Reaction('Sucrose + Water -> 2Glucose', reactan='Sucrose', X=1.)
        self.ferm_g = Reaction('Glucose -> 2LacticAcid', 'Glucose', eff_g)
        self.ferm_x = Reaction('3Xylose -> 5LacticAcid', 'Xylose', eff_x)
        self.tau = tau
        self.N = N
        self.cooler = hx = HXutility(None)
        self._heat_utilities = hx._heat_utilities
        hx.ins =hx.outsvent, effluent = self.outs
        hx._outs[0].T = effluent.T = vent.T = 305.15
        vent.phase = 'g'


    def _cal_efficiency(self, feed, tau):
        # Get initial concentrations
        y, la, g, x, w = feed.indices(['Yeast',
                                   '50-21-5', # lactic acid
                                   '492-61-5', # beta-D-glucose
                                   '58-86-6', # D-xylose
                                   '7732-18-5' # water
                                   ])
        mass = feed.mass
        volnet = feed.volnet
        concentration_in = mass/volnet
        X0, P0, S0_g, S0_x = (concentration_in[i] for i in (y, la, g, x))
        S0 = S0_g + S0_x
        
        # Integrate to get final concentration
        t = np.linspace(0, tau, 1000)
        C_t = odeint(self.kinetic_model, (X0, P0, S0, S0_g, S0_x), t)
        
        # Cache data
        self._X = C_t[:, 0]
        self._P = C_t[:, 1]
        self._S = S = C_t[:, 2]
        self._S_g = S_g = C_t[:, 3]
        self._S_x = S_x = C_t[:, 4]
        
        # Calculate fermentation efficiency
        Sf_g = (S_g[-1], S_x[-1]
        Sf_g = Sf_g if Sf_g >0 else 0
        Y_PS_g = self.kinetic_constants['Y_PS'][0]
        eff_g = (S0 - Sf)/S0 * Y_PS/1.
        Sf_x = S_x[-1]
        Sf_x = Sf_g if Sf_g >0 else 0
        Y_PS_g = self.kinetic_constants['Y_PS'][0]
        eff_g = (S0 - Sf)/S0 * Y_PS/1.
        return eff

    def _cal_efficiency(self, feed, tau):
        # Get initial concentrations
        y, la, g, x, w = feed.indices(['Yeast',
                                   '50-21-5', # lactic acid
                                   '492-61-5', # beta-D-glucose
                                   '58-86-6', # D-xylose
                                   '7732-18-5' # water
                                   ])
        mass = feed.mass
        volnet = feed.volnet
        concentration_in = mass/volnet
        X0, P0, S_g0, S_x0 = (concentration_in[i] for i in (y, la, g, x))
    
        # Integrate to get final concentration
        t = np.linspace(0, tau, 1000)
        C_t = odeint(self.kinetic_model, (X0, P0, S_g0, S_x0), t,
                  args=self.kinetic_constants)
        
        # Cache data
        self._X = C_t[:, 0]
        self._P = C_t[:, 1]
        self._S_g = S_g = C_t[:, 2]
        self._S_x = S_x = C_t[:, 3]
        
        # Calculate fermentation efficiency
        S_gf = S_g[-1]
        S_xf = S_x[-1]
        S_gf = S_gf if S_gf >0 else 0
        S_xf = S_xf if S_xf >0 else 0
        eff_g = (S_g0 - S_gf)/S_g0 * self.Y_PS_g/1. # 100% g-to-la under ideal conditions
        eff_x = (S_x0 - S_xf)/S_x0 * self.Y_PS_x/1. # 100% x-to-la under ideal conditions
        return eff_g, eff_x
    
    @staticmethod
    def kinetic_model(z, t, *kinetic_constants) -> '(dXdt, dPdt, dS_gdt, dS_xdt)':
        """
        Return change of yeast, lactic acid, and substrate (glucose and xylose) concentrations in kg/m3
        
        Parameters
        ----------
        z: iterable with (X, P, S_g, S_x)
        t: float, time
        kinetic_constants: kinetic parameters as defined earlier
        """
        
        n, P_max, Y_XS_g, Y_XS_x, Y_PS_g, Y_PS_x, mu_max_g, mu_max_x, K_S_g, K_S_x = kinetic_constants
        
        # Current yeast, lactic acid, and glucose concentrations in kg/m3
        X, P, S_g, S_x = z
        
        # Compute coefficients
        toxic = 1 - (P/P_max)^n
        mu_g = mu_max_g * S_g/(K_S_g + K_S_x)
        mu_x = mu_g/mu_max_g * mu_max_x
    
        # Compute derivatives
        dXdt = X * toxic * (mu_g + mu_x)
        dS_gdt = -dXdt * (1/Y_XS_g)
        dS_xdt = -dXdt * (1/Y_XS_x)
        dPdt = dXdt * ((Y_PS_g/Y_XS_g) + (Y_PS_x/Y_XS_x))
        return (dXdt, dPdt, dS_gdt, dS_xdt)

    
    @property
    def N(self):
        """[int] Number of reactors"""
        return self._N
    @N.setter
    def N(self, N):
        if N <= 1:
            raise ValueError(f'number of reactors must be greater than 1, value {N} is infeasible')
        self._N = N
        
    @property
    def eff_g(self):
        return self.ferm_g.X
    @eff_g.setter
    def eff_g(self, eff_g):
        self.ferm_g.X = eff_g
        
    @property
    def eff_x(self):
        return self.ferm_x.X
    @eff_x.setter
    def eff_x(self, eff_x):
        self.ferm_x.X = eff_x

    @property
    def tau(self):
        return self.tau
    @tau.setter
    def tau(self, tau):
        self._tau = tau
        
    def _run(self):
        vent, effluent = self.outs
        Stream.sum(effluent, self.ins)
        effluent_mol = effluent.mol
        self.hydrolysis(effluent_mol)
        if self.iskinetic_g:
            self.ferm_g.X = self._cal_eff_g(effluent, self._tau)
        self.ferm_g(effluent_mol)
        if self.iskinetic_x:
            self.ferm_x.X = self._cal_eff_x(effluent, self._tau)
        self.ferm_x(effluent_mol)
        vent.copyflow(effluent, ('CO2',), remove=True)
        vent.recieve_vent(effluent)

    @property
    def N_at_minimum_capital_cost(self):
        cost_old = np.inf
        self._N, N = 2, self._N
        cost_new = self.purchase_cost
        self.summary()
        while cost_new < cost_old:
            self._N += 1
            self._summary()
            cost_old = cost_new
            cost_new = self.purchase_cost
        self._N, N = N, self._N
        return N - 1
    
    def _design(self):
        v_0 = self.outs[1].volnet
        tau = self._tau
        tau_0 = self.tau_0
        Design = self._Design
        if self.autoselect_N:
            self.autoselec_N = False
            self._N = self.N_at_minimum_capital_cost
            self.autoselec_N = True
        N = self._N
        Design.update(size_batch(v_0, tau, tau_0, N, self._V_wf))
        hx = self._cooler
        hx.outs[0]._mol[:] = self.outs[0].mol/N
        hu = hx._heat_utilities[0]
        hu(self.Hnet/N, self.outs[0].T)
        hx._design(hu.duty)
        hx._cost()
        hu.duty *= N
        hu.cost *= N
        hu.flow *= N
        self._Cost['Coolers'] = self.cooler._Cost['Heat exchanger'] * self._N











        
        
        