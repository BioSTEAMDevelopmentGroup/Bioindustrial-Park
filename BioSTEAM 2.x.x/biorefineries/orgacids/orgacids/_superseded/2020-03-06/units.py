#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 16:19:18 2019

Based on the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

All units are explicitly defined here for transparency and easy reference

@author: yalinli_cabbi
"""

'''
Considerations for different organic acids
    Fermentation kinetic parameters (update EB module)   
    Fermentation conditions
       Temperature
       Headspace (air/O2 or not?)
       Nutrient requirements
       Cell composition
'''

'''
TODO:
   Separation units
   Check equipment costs
'''


# %% Setup

import numpy as np
import biosteam as bst
import thermosteam as tmo
from flexsolve import aitken_secant
from scipy.integrate import odeint
from biosteam import Unit
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch
from biosteam.units._splitter import run_split_with_mixing
from thermosteam import Stream, MultiStream

# Used by Sarang's unit
#from scipy.optimize import fsolve
#import math

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
_lb2MT = 0.000453592 # MT is metric tonne
_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
_hp2kW = 0.7457
_Gcal2kJ = 4184e3


# %% Feedstock handling

# Feedstock handling system as a whole, capital and operating costs considered in feedstock cost
# CE = 521.9 is year 2009
@cost(basis='Flow rate', ID='System', units='kg/hr',
      S=94697, ub=False, CE=521.9, cost=3329690, n=0.6, kW=783, BM=1.7)
class FeedstockHandling(Unit): pass


# %% Pretreatment

# Sulfuric acid addition tank
# CE = 550.8 is year 2010
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1981, ub=False, CE=550.8, cost=6210, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=3720, ub=False, CE=521.9, cost=8000, n=0.8, kW=1, BM=2.3)
class SulfuricAcidAdditionTank(Unit): pass

# Sulfuric acid in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      S=136260, ub=False, CE=521.9, cost=6000, n=0.6, kW=0, BM=1)
class SulfuricAcidMixer(Unit):
    _N_ins = 2
    _N_outs = 1
    _graphics = bst.units.Mixer._graphics
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        
    def _run(self):
        acid, water = self.ins
        mixture = self.outs[0]
        # 0.05 is from 1842/36629 from streams 710 and 516 of Humbird et al.
        water.imass['Water'] = acid.imass['SulfuricAcid'] / 0.05
        mixture.mix_from([water, acid])
        
# Blowdown discharge pump
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=25365, n=0.8, kW=93.21, BM=2.3)
class BlowdownDischargePump(Unit): pass

# Pretreatment flash tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=252891, ub=False, CE=521.9, cost=511000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=252891, ub=False, CE=521.9, cost=30000, n=0.8, kW=55.9275, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=252891, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class PretreatmentFlash(bst.units.Flash): pass

# Oligomer conversion tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=203000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=17408, n=0.8, kW=55.9, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class OligomerConversionTank(Unit): pass

# Ammonia in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      S=157478, ub=False, CE=521.9, cost=5000, n=0.5, kW=0, BM=1)
class AmmoniaMixer(bst.units.Mixer): pass

# Ammonia addition tank
# Assumed all sulfuric acid has been reacted with biomass or neutralized by ammonia
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=410369, ub=False, CE=521.9, cost=236000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=410369, ub=False, CE=521.9, cost=21900, n=0.5, kW=7.457, BM=1.5)
class AmmoniaAdditionTank(Unit): 
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        
        self.reaction = Rxn('2 NH3 + H2SO4 -> NH4SO4 + 2 H2O', 'NH3', 1.)
    
    def _run(self):
        ins = self.ins[0]
        outs = self.outs[0]
        outs.copy_like(ins)
        
        self.reaction(outs.mol)
        outs.imol['H2SO4'] = 0
            
# Waste vapor condenser
@cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
      S=-2, ub=False, CE=521.9, cost=34000, n=0.7, kW=0, BM=2.2)
class WasteVaporCondenser(bst.units.HXutility):
    _graphics = bst.units.HXutility._graphics

# Steam mixer
class SteamMixer(Unit):
    """
    **ins**
    
        [0] Feed
        
        [1] Steam
    
    **outs**
    
        [0] Mixed
    
    """
    _N_ins = 2
    _N_outs = 1
    _N_heat_utilities = 1
    def __init__(self, ID='', ins=None, outs=(), *, P):
        super().__init__(ID, ins, outs)
        self.P = P
    
    @staticmethod
    def _P_at_flow(mol_water, P, steam, mixed, feed):
        steam.imol['Water'] = mol_water
        mixed.mol[:] = steam.mol + feed.mol
        mixed.H = feed.H + mol_water * 40798 #!!! not sure what 40798 is
        P_new = mixed.chemicals.Water.Psat(mixed.T)
        return P - P_new
    
    def _run(self):
        feed, steam = self.ins
        steam_mol = steam.F_mol
        mixed = self.outs[0]
        steam_mol = aitken_secant(self._P_at_flow,
                                  steam_mol, steam_mol+0.1, 
                                  1e-4, 1e-4,
                                  args=(self.P, steam, mixed, feed))
        mixed.P = self.P
        hu = self.heat_utilities[0]
        hu.ID = 'Low pressure steam'
        hu.flow = steam_mol
        hu.cost = steam_mol*bst.HeatUtility.heating_agents['Low pressure steam'].price_kmol
    
    @property
    def installation_cost(self): return 0
    @property
    def purchase_cost(self): return 0
    def _design(self): pass
    def _cost(self): pass
        
# Pretreatment reactor
@cost('Dry flow rate', units='kg/hr', S=83333, CE=521.9, ub=False,
      cost=19812400, n=0.6, kW=4578, BM=1.5)
class PretreatmentReactorSystem(Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = bst.units.Flash._graphics
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self._multistream = MultiStream(None)
        self.reactions = ParallelRxn([
    #            Reaction definition                 Reactant    Conversion
    Rxn('Glucan + H2O -> Glucose',                   'Glucan',   0.0990),
    Rxn('Glucan + H2O -> GlucoseOligomer',           'Glucan',   0.0030),
    Rxn('Glucan -> HMF + 2 H2O',                     'Glucan',   0.0030),
    Rxn('Galactan + H2O -> GalactoseOligomer',       'Galactan', 0.0030),
    Rxn('Galactan -> HMF + 2 H2O',                   'Galactan', 0.0030),
    Rxn('Mannan + H2O -> MannoseOligomer',           'Mannan',   0.0030),
    Rxn('Mannan -> HMF + 2 H2O',                     'Mannan',   0.0030),
    Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose',  1.0000),
    Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.9000),
    Rxn('Xylan + H2O -> XyloseOligomer',             'Xylan',    0.0024),
    Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.0050),
    Rxn('Arabinan + H2O -> Arabinose',               'Arabinan', 0.9000),
    Rxn('Arabinan + H2O -> ArabinoseOligomer',       'Arabinan', 0.0024),
    Rxn('Arabinan -> Furfural + 2 H2O',              'Arabinan', 0.0050),
    Rxn('Acetate -> AceticAcid',                     'Acetate',  1.0000),
    Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.0050)
    ])
        vapor, liquid = self.outs
        vapor.phase = 'g'
    
    def _run(self):
        ms = self._multistream
        feed = self.ins[0]
        vapor, liquid = self.outs
        liquid.copy_like(feed)
        self.reactions(liquid.mol) 
        ms.copy_like(liquid)
        H = ms.H + liquid.Hf - feed.Hf
        ms.vle(T=130+273.15, H=H)
        vapor.mol[:] = ms.imol['g']
        liquid.mol[:] = ms.imol['l']
        vapor.T = liquid.T = ms.T
        vapor.P = liquid.P = ms.P


# %% Saccharification and fermentation

# Seed hold tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=40414, ub=False, CE=521.9, cost=439000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=43149, ub=False, CE=521.9, cost=8200, n=0.8, kW=7.457, BM=2.3)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=40414, ub=False, CE=521.9, cost=31800, n=0.5, kW=11.1855, BM=1.5)
class SeedHoldTank(Unit): pass

# Hydrolysate cooler
@cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
      S=-8, ub=False, CE=521.9, cost=85000, n=0.7, kW=0, BM=2.2)
class HydrolysateCooler(bst.units.HXutility): pass

# Enzyme hydrolysate mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      S=380000, ub=False, CE=521.9, cost=109000, n=0.5, kW=74.57, BM=1.7)
class EnzymeHydrolysateMixer(bst.units.Mixer): pass

# Seed train (5 stages)
@cost('Flow rate', 'Pumps',
      S=43149, CE=521.9, cost=24800, n=0.8, kW=40, BM=2.3)
@cost('Stage #1 reactor volume', 'Stage #1 reactors',
      cost=37700, S=20*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #2 reactor volume', 'Stage #2 reactors',
      cost=58300, S=200*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #3 reactor volume', 'Stage #3 reactors',
      cost=78800, S=2e3*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 reactors',
      cost=176e3, S=20e3*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 agitators',
      cost=26e3/2, S=20e3*_gal2m3, kW=7.5, CE=521.9, n=0.5, BM=1.5)
@cost('Stage #5 reactor volume', 'Stage #5 reactors',
      cost=590e3, S=200e3*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #5 reactor volume', 'Stage #5 agitators',
      cost=43e3/2, S=200e3*_gal2m3, kW=10, CE=521.9, n=0.5, BM=1.5)
class SeedTrain(Unit):
    '''
    Fermentation of organic acids using environmental biotechnology (EB) model.

    Structure of EB model based on:
    (1) Tchobanoglous, G.; Burton, F. L.; Stensel, H. D.; Inc, M. & E.; Burton, F. Wastewater Engineering: Treatment and Reuse; McGraw-Hill Education, 2003.
    (2) Rittmann, B. E.; McCarty, P. L. Environmental Biotechnology: Principles and Applications; McGraw-Hill, 2001.

    as:
        dX/dt = mu_max * X * (S/(K_S+S)) * (1-(P/P_max)^n)
        dS/dt = -dX/dt * (1/Y_X/S)
        dP/dt = -dX/dt * (Y_P/S/Y_X/S)
    
    where:
        X is biomass
        S is substrate, g for glucose and x for xylose
        P is product, P_max is maximum product titer
        mu is specific cell growth rate, mu = dX/dt/X, mu_max is the maximum specific cell growth rate
        Y_X/S is apparent yield of biomass per substrate consumed, Y_X/S = -dX/dS
        Y_P/S is apparent yield of product per substrate consumed, Y_P/S = -dP/dS
        n is toxic power
        K_S is half-velocity constnat, the value of S when μ/μ_max = 0.5
    '''
    
    _N_heat_utilities = 1
    _units= {'Flow rate': 'kg/hr',
             'Stage #1 reactor volume': 'm3',
             'Stage #2 reactor volume': 'm3',
             'Stage #3 reactor volume': 'm3',
             'Stage #4 reactor volume': 'm3',
             'Stage #5 reactor volume': 'm3'}
    
    # Intial biomass loading to get things started, won't be needed later
    X_init = 0.5 # [g/L] (i.e., [kg/m3])
    
    # Number of parallel seed trains
    N_trains = 2
    
    # Cycle time for each batch (hr)
    tau_batch = 24

    # Operating temperature (K), updated for the target acid
    T = 32+273.15
    
    # Kinetic parameters
    n = 3 # [-]
    P_max = 61.5 # [kg/m3]
    Y_XS_g = 0.0783333333333333 # [kg/kg]
    Y_XS_x = 0.05728125 # [kg/kg]
    Y_PS_g = 0.273809523809524 # [kg/kg]
    Y_PS_x = 0.659722222222222 # [kg/kg]
    mu_max_g = 0.22773560652439 # [kg/kg]
    mu_max_x = 0.114883404331077 # [kg/kg]
    K_S_g = 0 # [kg/m3]
    K_S_x = 45 # [kg/m3]
    EtOH_over_LA = 1.04347826086956 #[g/g]
    AceA_over_LA = 0 #[g/g]
    
    @property
    def N_stages(self): 
        """Number of stages in series."""
        return 5

    @property
    def tau_turnover(self):
        """Turnover time (hr) calculated by batch time divided by number of trains."""
        return self.tau_batch/self.N_trains
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)

    def kinetic_model(self, C, t) -> '(dXdt, dPdt, dS_gdt, dS_xdt)':
        n = self.n
        P_max = self.P_max
        Y_XS_g = self.Y_XS_g
        Y_XS_x = self.Y_XS_x
        Y_PS_g = self.Y_PS_g
        Y_PS_x = self.Y_PS_x
        mu_max_g = self.mu_max_g
        mu_max_x = self.mu_max_x
        K_S_g = self.K_S_g
        K_S_x = self.K_S_x

        X, P, S_g, S_x = C
    
        # Compute coefficients
        toxic = 1 - (P/P_max)**n
        mu_g = max(0, mu_max_g * S_g/(K_S_g + S_g + S_x))
        mu_x = max(0, mu_max_x * S_x/(K_S_x + S_g + S_x))
    
        # Compute derivatives
        dX_gdt = X * toxic * mu_g
        dX_xdt = X * toxic * mu_x
        dXdt = dX_gdt +dX_xdt
        dP_gdt = dX_gdt * (Y_PS_g/Y_XS_g)
        dP_xdt = dX_xdt * (Y_PS_x/Y_XS_x)
        dPdt = dP_gdt + dP_xdt
        dS_gdt = -dX_gdt * (1/Y_XS_g) if S_g>0 else 0
        dS_xdt = -dX_xdt * (1/Y_XS_x) if S_x>0 else 0
        return (dXdt, dPdt, dS_gdt, dS_xdt)
    
    def EB_simulate(self, feed, tau):       
        # Get initial concentrations
        mass_in = feed.imass['FermentationMicrobe', 'LacticAcid', 'Glucose', 'Xylose']
        vol = feed.F_vol
        #concentration_in = mass_in/vol
        X0, P0, S_g0, S_x0 = mass_in/vol
        
        C0 = [X0, P0, S_g0, S_x0]
    
        # Integrate to get final concentration
        t = np.linspace(0, tau, 1000)
        Ct = odeint(self.kinetic_model, y0=C0, t=t)
        
        # Cache data
        Xt = Ct[:, 0]    
        Pt = Ct[:, 1]
        S_gt = Ct[:, 2]
        S_xt = Ct[:, 3]
     
        # Calculate final fermentation efficiency @ t = tau
        Xf = Xt[-1]
        Pf = Pt[-1]
        S_gf = S_gt[-1]
        S_gf = S_gf if S_gf >0 else 0
        S_xf = S_xt[-1]
        S_xf = S_xf if S_xf >0 else 0
        Cf = [Xf, Pf, S_gf, S_xf]
        EtOH = (S_g0 - S_gf)* self.Y_PS_g * self.EtOH_over_LA
        AceA = (S_x0 - S_xf)* self.Y_PS_x * self.AceA_over_LA
        loss = sum(C0) - sum(Cf) - EtOH - AceA
        Cf.append(EtOH)
        Cf.append(AceA)
        Cf.append(loss)
        mass_out = tuple([i * vol for i in Cf])
        return mass_out

    def _run(self):
        feed = self.ins[0]
        vent, effluent= self.outs
        effluent.copy_flow(feed)
        EtOH = 0
        AceA = 0
        loss = 0
        if effluent.imass['FermentationMicrobe'] == 0:
                effluent.imass['FermentationMicrobe'] = self.X_init * \
                    effluent.F_vol
        for i in range(0, 5):
            simulated_flow = self.EB_simulate(effluent, self.tau_batch)
            effluent.imass['FermentationMicrobe'] = simulated_flow[0]
            effluent.imass['LacticAcid'] = simulated_flow[1]
            effluent.imass['Glucose'] = simulated_flow[2]
            effluent.imass['Xylose'] = simulated_flow[3]
            EtOH += simulated_flow[4]
            AceA += simulated_flow[5]
            loss += simulated_flow[6]
        effluent.imass['Ethanol'] = feed.imass['Ethanol'] + EtOH
        effluent.imass['AceticAcid'] = feed.imass['AceticAcid'] + AceA
        # Assume all lost glucose and xylose changed to CO2
        effluent.imass['CO2'] = feed.imass['CO2'] + loss
        effluent.T = self.T
        vent.phase = 'g'
        # Did not use copy_flow as shown below as it will only remove those identified
        #vent.copy_flow(effluent, ('CO2', 'NH3', 'O2'), remove=True)
        vent.recieve_vent(effluent)


    def _design(self): 
        maxvol = self.outs[1].F_vol*self.tau_turnover
        vol = maxvol*10**-self.N_stages
        Design = self.design_results
        for i in range(1, self.N_stages+1):
            Design[f'Stage #{i} reactor volume'] = vol
            vol *= 10 
        Design['Flow rate'] = sum([i.F_mass for i in self.outs])
        self.heat_utilities[0](self.Hnet, self.T)

    def _cost(self):
        N = self.N_trains
        D = self.design_results
        C = self.purchase_costs
        kW = 0
        for i, x in self.cost_items.items():
            S = D[x._basis]
            q = S/x.S
            C[i] = N*bst.CE/x.CE*x.cost*q**x.n
            kW += N*x.kW*q
        self.power_utility(kW)

# Saccharification and co-fermentation reactor,
# "co-fermentation" means the microbe can simultaneously ferment glucose and xylose
# Did not include heat exchanger
@cost('Flow rate', 'Recirculation pumps', kW=30, S=340*_gpm2m3hr, ub=False,
      cost=47200, n=0.8, BM=2.3, CE=521.9, N='N_recirculation_pumps')
@cost('Reactor duty', 'Heat exchangers', CE=521.9, cost=23900, ub=False,
      S=5*_Gcal2kJ, n=0.7, BM=2.2, N='N_reactors') # Based on a similar heat exchanger
@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500, ub=False,
      S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_reactors')
@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000, ub=False,
      S=1e6*_gal2m3, n=0.5, BM=1.5, N='N_reactors')
@cost('Flow rate', 'Transfer pumps', kW=58, S=352*_gpm2m3hr, ub=False,
      cost=47200/5, CE=521.9, n=0.8, BM=2.3, N='N_transfer_pumps')
@cost('Tank volume', 'Tanks', cost=3840e3/8, S=250e3*_gal2m3, ub=False,
      CE=521.9, n=0.7, BM=2.0, N='N_tanks')
class SaccharificationAndCoFermentation(Unit):
    '''
    Fermentation of organic acids using environmental biotechnology (EB) model.

    Structure of EB model based on:
    (1) Tchobanoglous, G.; Burton, F. L.; Stensel, H. D.; Inc, M. & E.; Burton, F. Wastewater Engineering: Treatment and Reuse; McGraw-Hill Education, 2003.
    (2) Rittmann, B. E.; McCarty, P. L. Environmental Biotechnology: Principles and Applications; McGraw-Hill, 2001.

    as:
        dX/dt = mu_max * X * (S/(K_S+S)) * (1-(P/P_max)^n)
        dS/dt = -dX/dt * (1/Y_X/S)
        dP/dt = -dX/dt * (Y_P/S/Y_X/S)
    
    where:
        X is biomass
        S is substrate, g for glucose and x for xylose
        P is product, P_max is maximum product titer
        mu is specific cell growth rate, mu = dX/dt/X, mu_max is the maximum specific cell growth rate
        Y_X/S is apparent yield of biomass per substrate consumed, Y_X/S = -dX/dS
        Y_P/S is apparent yield of product per substrate consumed, Y_P/S = -dP/dS
        n is toxic power
        K_S is half-velocity constnat, the value of S when μ/μ_max = 0.5
    '''
    
    _N_ins = 4
    _N_outs = 3
    _N_heat_utilities = 2
    
    # Saccharification temperature (K)
    T_saccharification = 48+273.15
    
    # Fermentation temperature (K), updated for the target acid
    T_fermentation = 30+273.15
    
    # Residence time of countinuous saccharification tanks (hr)
    tau_tank = 24
    
    # Saccharification time (hr)
    tau_saccharification = 60
    
    # Co-Fermentation time (hr)
    tau_cofermentation = 36
    
    # Unload and clean up time (hr)
    tau_0 = 4
    
    #: Working volume fraction (filled tank to total tank volume)
    V_wf = 0.9
    
    # Number of reactors
    N_reactors = 12
    
    # Number of continuous saccharification tanks
    N_tanks = 8
    
    # Number of transfer pumps
    N_transfer_pumps = 5
    
    # Number of recirculation pumps
    N_recirculation_pumps = 5
    
    _units = {'Flow rate': 'm3/hr',
              'Tank volume': 'm3',
              'Reactor volume': 'm3',
              'Reactor duty': 'kJ/hr'}
    
    # Split to outs[2]
    saccharified_slurry_split = 0.1
    
    # Intial biomass loading to get things started, won't be needed later
    X_init = 0.5 # [g/L] (i.e., [kg/m3])
    
    # Kinetic parameters
    n = 3 # [-]
    P_max = 61.5 # [kg/m3]
    Y_XS_g = 0.08 # [kg/kg]
    Y_XS_x = 0.06 # [kg/kg]
    Y_PS_g = 0.27 # [kg/kg]
    Y_PS_x = 0.66 # [kg/kg]
    mu_max_g = 0.228 # [kg/kg]
    mu_max_x = 0.115 # [kg/kg]
    K_S_g = 0 # [kg/m3]
    K_S_x = 45 # [kg/m3]
    EtOH_over_LA = 1.04 #[g/g]
    AceA_over_LA = 0 #[g/g]
    
    def __init__(self, ID='', ins=None, outs=(), P=101325):
        Unit.__init__(self, ID, ins, outs)
        
        self.P = P
        
        # Enzymatic hydrolysis reactions including from downstream batch tank in co-fermentation.
        # Kept the same as Humbird et al.
        self.saccharification = ParallelRxn([
            #   Reaction definition                   Reactant   Conversion
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000),
            Rxn('Cellobiose + H2O -> Glucose',        'Cellobiose',  1.0000)
            ])
        
        self.saccharified_stream = Stream(None)
        
        # Assume all CSL and DAP are consumed
        self.CSL_DAP2constituents = ParallelRxn([
            #   Reaction definition                              Reactant    Conversion
            Rxn('CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL',    1.0000),
            # Assume biomass picks up phosphate
            Rxn('DAP -> 2 NH3',                                    'DAP',    1.0000)
            ])
        
        # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
        self.neutralization = ParallelRxn([
            #   Reaction definition                                             Reactant  Conversion
            Rxn('2 LacticAcid + CalciumDihydroxide -> CalciumLactate + 2 H2O', 'LacticAcid',   1.),
            Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O', 'AceticAcid',   1.)
            ])
    
    def kinetic_model(self, C, t) -> '(dXdt, dPdt, dS_gdt, dS_xdt)':
        n = self.n
        P_max = self.P_max
        Y_XS_g = self.Y_XS_g
        Y_XS_x = self.Y_XS_x
        Y_PS_g = self.Y_PS_g
        Y_PS_x = self.Y_PS_x
        mu_max_g = self.mu_max_g
        mu_max_x = self.mu_max_x
        K_S_g = self.K_S_g
        K_S_x = self.K_S_x

        X, P, S_g, S_x = C
    
        # Compute coefficients
        toxic = 1 - (P/P_max)**n
        mu_g = max(0, mu_max_g * S_g/(K_S_g + S_g + S_x))
        mu_x = max(0, mu_max_x * S_x/(K_S_x + S_g + S_x))
    
        # Compute derivatives
        dX_gdt = X * toxic * mu_g
        dX_xdt = X * toxic * mu_x
        dXdt = dX_gdt +dX_xdt
        dP_gdt = dX_gdt * (Y_PS_g/Y_XS_g)
        dP_xdt = dX_xdt * (Y_PS_x/Y_XS_x)
        dPdt = dP_gdt + dP_xdt
        dS_gdt = -dX_gdt * (1/Y_XS_g) if S_g>0 else 0
        dS_xdt = -dX_xdt * (1/Y_XS_x) if S_x>0 else 0
        return (dXdt, dPdt, dS_gdt, dS_xdt)
    
    def EB_simulate(self, feed, tau):       
        # Get initial concentrations
        mass_in = feed.imass['FermentationMicrobe', 'LacticAcid', 'Glucose', 'Xylose']
        vol = feed.F_vol
        #concentration_in = mass_in/vol
        X0, P0, S_g0, S_x0 = mass_in/vol
        
        C0 = [X0, P0, S_g0, S_x0]
    
        # Integrate to get final concentration
        t = np.linspace(0, tau, 1000)
        Ct = odeint(self.kinetic_model, y0=C0, t=t)
        
        # Cache data
        Xt = Ct[:, 0]    
        Pt = Ct[:, 1]
        S_gt = Ct[:, 2]
        S_xt = Ct[:, 3]
     
        # Calculate final fermentation efficiency @ t = tau
        Xf = Xt[-1]
        Pf = Pt[-1]
        S_gf = S_gt[-1]
        S_gf = S_gf if S_gf >0 else 0
        S_xf = S_xt[-1]
        S_xf = S_xf if S_xf >0 else 0
        Cf = [Xf, Pf, S_gf, S_xf]
        EtOH = (S_g0 - S_gf)* self.Y_PS_g * self.EtOH_over_LA
        AceA = (S_x0 - S_xf)* self.Y_PS_x * self.AceA_over_LA
        loss = sum(C0) - sum(Cf) - EtOH - AceA
        Cf.append(EtOH)
        Cf.append(AceA)
        Cf.append(loss)
        mass_out = tuple([i * vol for i in Cf])
        return mass_out
    
    def _run(self):
        feed, CSL, DAP, lime = self.ins
        vent, effluent, sidedraw = self.outs
        vent.P = effluent.P = sidedraw.P = self.P
        ss = self.saccharified_stream
        ss.T = sidedraw.T = self.T_saccharification
        vent.T = effluent.T = self.T_fermentation
        vent.phase = 'g'
        ss.copy_flow(feed)
        self.saccharification(ss.mol)
        # Sidedraw after saccharification but before co-fermentation
        sidedraw.mol[:] = ss.mol * self.saccharified_slurry_split
        # Now effluent is after saccharification but before co-fermentation
        effluent.mol[:] = ss.mol - sidedraw.mol + CSL.mol + DAP.mol + lime.mol
        if effluent.imass['FermentationMicrobe'] == 0:
            effluent.imass['FermentationMicrobe'] = self.X_init * \
                effluent.F_vol
        simulated_flow = self.EB_simulate(effluent, self.tau_cofermentation)
        effluent.imass['FermentationMicrobe'] = simulated_flow[0]
        effluent.imass['LacticAcid'] = simulated_flow[1]
        effluent.imass['Glucose'] = simulated_flow[2]
        effluent.imass['Xylose'] = simulated_flow[3]
        effluent.imass['Ethanol'] += simulated_flow[4]
        effluent.imass['AceticAcid'] += simulated_flow[5]
        # Assume all lost glucose and xylose changed to CO2
        effluent.imass['CO2'] += simulated_flow[6]
        # Changed from effluent.mass to effluent.mol
        self.CSL_DAP2constituents(effluent.mol)
        # Set feed lime mol to match rate of acids production, add 5% extra
        lime.imol['Lime'] = (effluent.imol['LacticAcid']/2/self.neutralization.X[0] \
                           + effluent.imol['AceticAcid']/2/self.neutralization.X[1]) \
                           * 1.05
        self.neutralization(effluent.mol)
        vent.recieve_vent(effluent)

    
    def _design(self):
        effluent = self.outs[1]
        v_0 = effluent.F_vol
        Design = self.design_results
        Design['Tank volume'] = v_0*self.tau_tank/self.V_wf/self.N_tanks
        Design['Flow rate'] = v_0/self.N_transfer_pumps
        tau = self.tau_saccharification + self.tau_cofermentation
        Design.update(size_batch(v_0, tau, self.tau_0, self.N_reactors, self.V_wf))
        hu_cooling, hu_fermentation = self.heat_utilities
        mixture = self.thermo.mixture
        ss = self.saccharified_stream
        mol = ss.mol
        hu_cooling(mixture.H('l', mol, self.T_fermentation, 101325) \
                   - mixture.H('l', mol, self.T_saccharification, 101325),
                   self.T_fermentation
                   )
        diff = sum([i.imol['CalciumLactate'] for i in self.outs]) \
             + sum([i.imol['CalciumAcetate'] for i in self.outs]) \
             - sum([i.imol['CalciumLactate'] for i in self.ins]) \
             - sum([i.imol['CalciumAcetate'] for i in self.ins])
        duty = diff * -5568 #??? Where does this come from?
        hu_fermentation(duty, effluent.T)
        Design['Reactor duty'] = -duty


# %% Organic acid separation

@cost('Flow rate', 'Flitrate tank agitator', ub=False,
      cost=26e3, CE=550.8, kW=7.5*_hp2kW, S=31815, n=0.5, BM=1.5)
@cost('Flow rate', 'Discharge pump', ub=False,
      cost=13040, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Filtrate tank', ub=False,
      cost=103e3, S=31815, CE=550.8, BM=2.0, n=0.7)
@cost('Flow rate', 'Feed pump', kW=74.57, ub=False,
      cost= 18173, S=31815, CE=550.8, n=0.8, BM=2.3)
@cost('Flow rate', 'Stillage tank 531', ub=False,
      cost=174800, CE=550.8, S=31815, n=0.7, BM=2.0)
@cost('Flow rate', 'Manifold flush pump', kW=74.57, ub=False,
      cost=17057, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Recycled water tank', ub=False,
      cost=1520, CE=550.8, S=31815, n=0.7, BM=3.0)
@cost('Flow rate', 'Lignin wet cake screw',  kW=15*_hp2kW, ub=False,
      cost=2e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Lignin wet cake conveyor', kW=10*_hp2kW, ub=False,
      cost=7e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressure filter', ub=False,
      cost=3294700, CE=550.8, S=31815, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressing air compressor reciever tank', ub=False,
      cost=8e3, CE=550.8, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Cloth wash pump', kW=150*_hp2kW, ub=False,
      cost=29154, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Dry air compressor reciever tank', ub=False,
      cost=17e3, CE=550.8, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Pressing air pressure filter', ub=False,
      cost=75200, CE=521.9, S=31815, n=0.6, kW=112, BM=1.6)
@cost('Flow rate', 'Dry air pressure filter (2)', ub=False,
      cost=405000, CE=521.9, S=31815, n=0.6, kW=1044, BM=1.6)
class CellMassFilter(bst.units.SolidsSeparator):
    _units = {'Flow rate': 'kg/hr'}
    
    def _run(self):
        cell_mass, recycled_water = self.ins
        # No information on recycled_water in flow table in Humbird et al.,
        # therefore 36538 is taken from the scaling basis of equipment T-532
        # 391501 is total flow of stream 508 in Humbird et al.
        recycled_water.imass['Water'] = 36538 * cell_mass.F_mass / 391501
        
        run_split_with_mixing(self)
        retentate, permeate = self.outs
        solids = retentate.F_mass
        mc = self.mositure_content
        retentate.imol['Water'] = water = (solids * mc/(1-mc))/18.01528
        permeate.imol['Water'] -= water
        if permeate.imol['Water'] < water:
            raise ValueError(f'not enough water for {repr(self)}')
        
    def _design(self):
        self.design_results['Flow rate'] = self.outs[0].F_mass

# Cost copied from OligomerConversionTank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=203000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=17408, n=0.8, kW=55.9, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class AcidulationReactor(Unit):
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        # Conversions are placeholders now
        self.acidulation = ParallelRxn([
            #   Reaction definition                                    Reactant   Conversion
            Rxn('CalciumLactate + H2SO4 -> 2 LacticAcid + CaSO4', 'CalciumLactate', 0.95),
            Rxn('CalciumAcetate + H2SO4 -> 2 AceticAcid + CaSO4', 'CalciumAcetate', 0.95)
            ])
        
    def _run(self):
        feed, acid = self.ins
        effluent = self.outs[0]
        
        effluent.copy_like(feed)
        # Set feed acid mol to match acidulation needs with 5% extra
        acid.imol['H2SO4'] = (feed.imol['CalciumLactate']/self.acidulation.X[0] \
                             +feed.imol['CalciumAcetate']/self.acidulation.X[1]) \
                             * 1.05
        effluent.mix_from([feed, acid])
        self.acidulation(effluent.mol)
        
# Now costs are the same as the CellMassFilter, but should consider changing to
# S-222, the hydroclone & rotary drum filter as in Aden et al.        
@cost('Flow rate', 'Flitrate tank agitator', ub=False,
      cost=26e3, CE=550.8, kW=7.5*_hp2kW, S=31815, n=0.5, BM=1.5)
@cost('Flow rate', 'Discharge pump', ub=False,
      cost=13040, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Filtrate tank', ub=False,
      cost=103e3, S=31815, CE=550.8, BM=2.0, n=0.7)
@cost('Flow rate', 'Feed pump', kW=74.57, ub=False,
      cost= 18173, S=31815, CE=550.8, n=0.8, BM=2.3)
@cost('Flow rate', 'Stillage tank 531', ub=False,
      cost=174800, CE=550.8, S=31815, n=0.7, BM=2.0)
@cost('Flow rate', 'Manifold flush pump', kW=74.57, ub=False,
      cost=17057, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Recycled water tank', ub=False,
      cost=1520, CE=550.8, S=31815, n=0.7, BM=3.0)
@cost('Flow rate', 'Gypsum wet cake screw',  kW=15*_hp2kW, ub=False,
      cost=2e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Gypsum wet cake conveyor', kW=10*_hp2kW, ub=False,
      cost=7e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressure filter', ub=False,
      cost=3294700, CE=550.8, S=31815, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressing air compressor reciever tank', ub=False,
      cost=8e3, CE=550.8, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Cloth wash pump', kW=150*_hp2kW, ub=False,
      cost=29154, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Dry air compressor reciever tank', ub=False,
      cost=17e3, CE=550.8, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Pressing air pressure filter', ub=False,
      cost=75200, CE=521.9, S=31815, n=0.6, kW=112, BM=1.6)
@cost('Flow rate', 'Dry air pressure filter (2)', ub=False,
      cost=405000, CE=521.9, S=31815, n=0.6, kW=1044, BM=1.6)
#!!! Need to look into Aden et al. and see how the filter was designed there
class GypsumFilter(bst.units.SolidsSeparator):
    _units = {'Flow rate': 'kg/hr'}
    
    def _design(self):
        self.design_results['Flow rate'] = self.outs[0].F_mass
        
class OrganicAcidsFlash(bst.units.Flash):
    vle_chemicals = ()
    V = None
    P = 101325
    
    def _run(self):
        feed = self.ins[0]
        vap, liq = self.outs
        
        # Vapor-liquid equilibrium
        ms = self._multistream
        ms.empty()
        vle_chemicals = self.vle_chemicals
        ms.imol['l', vle_chemicals] = feed.imol[vle_chemicals]
        ms.vle(P=self.P, V=self.V)

        # Set values
        vap.phase = 'g'
        liq.phase = 'l'
        liq.mol = feed.mol
        vap.imol[vle_chemicals] = ms.imol['g', vle_chemicals]
        liq.imol[vle_chemicals] = ms.imol['l', vle_chemicals]
        vap.T = liq.T = ms.T
        vap.P = liq.P = ms.P

# Cost copied from OligomerConversionTank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=203000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=17408, n=0.8, kW=55.9, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class EsterificationReactor(Unit):
    _N_ins = 2
    _N_outs = 1

    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        # Conversions are placeholders now
        self.esterification = ParallelRxn([
            #   Reaction definition                               Reactant  Conversion
            Rxn('LacticAcid + Methanol -> MethylLactate + H2O', 'LacticAcid', 0.98),
            Rxn('AceticAcid + Methanol -> MethylAcetate + H2O', 'AceticAcid', 0.98)
            ])
        
    # def _run(self):
    #     feed, recycled, methanol = self.ins
    #     effluent = self.outs[0]
        
    #     # Add 5% extra
    #     methanol_needed = (feed.imol['LacticAcid']/self.esterification.X[0] \
    #                         +feed.imol['AceticAcid']/self.esterification.X[1]) \
    #                       * 1.05
                          
    #     methanol.imol['Methanol'] = max(0, methanol_needed-recycled.imol['Methanol'])
    #     effluent.mix_from([feed, recycled, methanol])
    #     self.esterification(effluent.mol)
    
    def _run(self):
        feed, methanol = self.ins
        effluent = self.outs[0]
        
        # Add 5% extra
        methanol_needed = (feed.imol['LacticAcid']/self.esterification.X[0] \
                            +feed.imol['AceticAcid']/self.esterification.X[1]) \
                          * 1.05
                          
        methanol.imol['Methanol'] = methanol_needed
        effluent.mix_from([feed, methanol])
        self.esterification(effluent.mol)
        
# Cost copied from OligomerConversionTank
#!!! Need to add ethanol into consideration?
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=203000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=17408, n=0.8, kW=55.9, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class HydrolysisReactor(Unit):
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=(), outs=()):
        Unit.__init__(self, ID, ins, outs)
        # Conversions are placeholders now
        self.hydrolysis = ParallelRxn([
            #   Reaction definition                                Reactant    Conversion
            Rxn('MethylLactate + H2O -> LacticAcid + Methanol', 'MethylLactate', 0.98),
            Rxn('MethylAcetate + H2O -> AceticAcid + Methanol', 'MethylAcetate', 0.98)
            ])
    
    def _run(self):
        feed, water = self.ins
        effluent = self.outs[0]
        
        # Add 5% extra
        water_needed = (feed.imol['MethylLactate']/self.hydrolysis.X[0] \
                        +feed.imol['MethylAcetate']/self.hydrolysis.X[1]) \
                       * 1.05
        water.imol['Water'] = max(0, (water_needed-feed.imol['Water'])) 
        effluent.mix_from([feed, water])
        
        self.hydrolysis(effluent.mol)


# %% Wastewater treatment

# The total cost of wastewater treatment is combined into this placeholder
@cost('Flow rate', 'Wastewater system', units='kg/hr', CE=550.8, ub=False,
      cost=50280080., n=0.6, BM=1, kW=7139/1.05, S=393100)
class WastewaterSystemCost(Unit): pass

class AnaerobicDigestion(Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **reactions:** [ReactionSet] Anaerobic digestion reactions.
        
        **sludge_split:** [Array] Split between wastewater and sludge
        
    **ins**
    
        [0] Wastewater
        
        [1] Cool well water
        
    **outs**
    
        [0] Biogas
        
        [1] Wastewater
        
        [2] Sludge
        
        [3] Hot well water
    
    """
    purchase_cost = installation_cost = 0
    _N_ins = 2
    _N_outs = 4
    def __init__(self, ID='', ins=None, outs=(), *, reactions, sludge_split):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions
        self.sludge_split = sludge_split
        self.multi_stream = MultiStream()
    
    def _run(self):
        feed, cool_water = self.ins
        biogas, waste, sludge, hot_water = self.outs
        biogas.phase = 'g'
        hot_water.link_with(cool_water, TP=False)
        biogas.T = waste.T = sludge.T = T = 35+273.15
        # Hot water temperature is adjusted to maintain heat balance
        hot_water.T = feed.T - 5
        H_at_35C = feed.thermo.mixture.H(z=feed.mol, phase='l', T=T, P=101325)
        cool_water.mol[:] *= (feed.H - H_at_35C)/(hot_water.H - cool_water.H)
        while cool_water.F_mol < 0:
            hot_water.T -= 1
            cool_water.mol[:] *= (feed.H - H_at_35C)/(hot_water.H - cool_water.H)
        sludge.copy_flow(feed)
        self.reactions(sludge.mol)
        self.multi_stream.copy_flow(sludge)
        self.multi_stream.vle(P=101325, H=self.multi_stream.H)
        biogas.mol[:] = self.multi_stream.imol['g']
        liquid_mol = self.multi_stream.imol['l']
        sludge.mol[:] = liquid_mol * self.sludge_split
        waste.mol[:] = liquid_mol - sludge.mol
        biogas.recieve_vent(waste)     
    
class AerobicDigestion(Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **reactions:** [ReactionSet] Anaerobic digestion reactions.
        
        **sludge_split:** [Array] Split between wastewater and sludge
        
    **ins**
    
        [0] Wastewater
        
        [1] Air
        
        [2] Caustic
        
    **outs**
    
        [0] Vent
        
        [1] Treated wastewater

    """    
    _N_ins = 3
    _N_outs = 2
    purchase_cost = installation_cost = 0
    evaporation = 4/355
    
    def __init__(self, ID='', ins=None, outs=(), *, reactions):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions
    
    def _run(self):
        waste, air, caustic = self._ins
        vent, water = self.outs
        vent.phase = 'g'
        water.copy_like(waste)
        water.mol[:] += air.mol
        water.mol[:] += caustic.mol
        self.reactions(water.mol)
        vent.copy_flow(water, ('CO2', 'O2', 'N2'))
        vent.imol['Water'] = water.imol['Water'] * self.evaporation
        water.mol[:] -= vent.mol


# %% Facility units in bst.units

# Sulfuric acid storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1981, ub=False, CE=550.8, cost=96000, n=0.7, kW=0, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=1981, ub=False, CE=521.9, cost=7493, n=0.8, kW=0.37285, BM=2.3)
class SulfuricAcidStorageTank(bst.units.StorageTank): pass

# Ammonia storage tank, no pump
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1171, ub=False, CE=550.8, cost=196000, n=0.7, kW=0, BM=2)
class AmmoniaStorageTank(bst.units.StorageTank): pass

# CSL storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1393, ub=False, CE=521.9, cost=70000, n=0.7, kW=0, BM=2.6)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=1393, ub=False, CE=550.8, cost=3000, n=0.8, kW=0.37285, BM=3.1)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=1393, ub=False, CE=521.9, cost=21200, n=0.5, kW=7.457, BM=1.5)
class CSLStorageTank(Unit): pass

# DAP storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=163, ub=False, CE=521.9, cost=102000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=163, ub=False, CE=521.9, cost=3000, n=0.8, kW=0.37285, BM=3.1)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=163, ub=False, CE=521.9, cost=9800, n=0.5, kW=4.10135, BM=1.5)
@cost(basis='Flow rate', ID='Bag', units='kg/hr',
      S=163, ub=False, CE=521.9, cost=30000, n=0.6, kW=0, BM=1.7)
class DAPStorageTank(Unit): pass

# For storage of lime used in separation and waste treatment,
# design copied from cornstover biorefinery in Aden et al.
# CE = 386.5 and 389.5 is year 1997 and 1998, respectively
# Base flow from 227 of Aden et al.
@cost(basis='Flow rate', ID='Storage bin', units='kg/hr',
      S=227, ub=False, CE=386.5, cost=136370, n=0.46, kW=0, BM=1.3)
@cost(basis='Flow rate', ID='Dust vent baghouse', units='kg/hr',
      S=227, ub=False, CE=386.5, cost=140707, n=1, kW=0, BM=1.5)
# Power usage scaled based on M-106 of Humbird et al.
@cost(basis='Flow rate', ID='Unloading blower', units='kg/hr',
      S=227, ub=False, CE=389.5, cost=99594, n=0.5, kW=25*_hp2kW*7425/8500, BM=1.4)
# Cost not scaled, thus used n=0
# Power usage scaled based on C-101 of Humbird et al., 
# though C-101 is a belt and here is a rotary type
@cost(basis='Flow rate', ID='Feeder', units='kg/hr',
      S=227, ub=False, CE=386.5, cost=3900, n=0, kW=20*_hp2kW*7000*_lb2MT/160, BM=1.3)
class LimeStorageTank(Unit): pass

# Methanol storage tank, based on Ethanol storage tank from Humbird et al.
# scaled down from 7-d to 14-hr storage
# (similar to sulfuric acid/CSL/DAP storage in Humbird et al.)
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=21808/24/7*14, ub=False, CE=521.9, cost=1340000, n=0.7, kW=0, BM=1.7)
class MethanolStorageTank(Unit): pass

# Fire water tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=8343, ub=False, CE=521.9, cost=803000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=8343, ub=False, CE=521.9, cost=15000, n=0.8, kW=93.2125, BM=1.7)
class FireWaterTank(Unit): pass


# %% Sarang's work in progress
'''
class MultiComponentDistillation(Unit, isabstract=True):
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 P=101325, *, LHK, rLKD, rHKB, K_HK, K_specs):
        super().__init__(ID, ins, outs)
        self.outs[0].phase = 'g'
        
        self.P = P_in = ins.P
        self.T = T_in = ins.T
        self.LHK = LHK
        self.rLKD = rLKD
        self.rHKB = rHKB
        self.K_HK = K_HK
        self.K_specs = K_specs
        
        # equilibrium_chemicals should be an iterable of Chemical objects
        #??? Why use Water and Ethanol?
        equilibrium_chemicals = tmo.Chemicals(['Ethanol', 'Water'])
        F_l = tmo.equilibrium.LiquidFugacities(equilibrium_chemicals)
        liquid_molar_composition = np.array([0.9, 0.1])
        fugacities = F_l(x=liquid_molar_composition, T=T_in)
        
    def GetRelativeVolatilities(self):
        K_HK = self.K_HK
        K_specs = self.K_specs
        
        rel_vols = []
        for K_i in K_specs:
            rel_vols.append(K_i/K_HK)        
       
        self.rel_vols = rel_vols
        
        return rel_vols
    
    def EstimateMinStages(self):
        LHK = self.LHK
        rLKD = self.rLKD
        rHKB = self.rHKB
        
        alpha_LH = self.rel_vols[self.ins[0].chemicals.index(LHK[0])]
        Nm = math.log10((rLKD/(1-rHKB)) * (rHKB/(1-rLKD))) \
            /math.log10(alpha_LH)
        
        self.alpha_LH = alpha_LH
        self.Nm = Nm
        
        return Nm
    
    def EstimateDBCompositions(self):
        rHKB = self.rHKB
        Nm = self.Nm
        rel_vols = self.rel_vols
        
        A = math.log10((1-rHKB) / rHKB)*1**(-Nm)
        D_by_B_ratios = []
        for rv_i in rel_vols:
            D_by_B_ratios.append(10**(A + Nm*math.log10(rv_i)))        
        
        self.D_by_B_ratios = D_by_B_ratios
        
        DB = []
        mass = self.ins[0].F_mass
        for ind in range(len(D_by_B_ratios)):
            F = self.ins[0].mass[ind] / mass
            d = F / (1 + 1 / (D_by_B_ratios[ind]))
            b = F - d
            DB.append((d,b))
        
        self.DB = DB
        
        return DB
    
    def Underwood_1_single(self, ind, theta):
        return self.rel_vols[ind] * (self.DB[ind][0] + self.DB[ind][1]) \
               / (self.rel_vols[ind] - theta)

    """ Underwood Ean. # 1"""
    #!!! Ean. or Eqn.?
    def Underwood_1(self, theta):
        sum = 0.0000
        
        for ind in range(len(self.DB)):
            sum += self.Underwood_1_single(ind, theta)
            
        return sum

    def Underwood_2_single(self, ind):
        return self.rel_vols[ind] * (self.DB[ind][0]) / (self.rel_vols[ind] - self.theta)

    """Underwood Ean. #2"""
    def Underwood_2(self, Rm):
        sum = 0.0000
        
        for ind in range(len(self.DB)):
            sum += self.Underwood_2_single(ind)
        
        sum -= (Rm + 1)
            
        return sum

    """Uses Underwood Eqns."""
    def EstimateOperatingRefluxRatio(self):
        theta = fsolve(self.Underwood_1, [1, 1, 1])
        self.theta = theta
        
        Rm = fsolve(self.Underwood_2, [1, 1, 1])
        self.Rm = Rm
        
        self.R_heuristic = 1.2
        self.R_op = self.Rm * self.R_heuristic
        
        return self.R_op
    
    """Gilliland Eqn."""
    def Gilliland(self, N):
        N_ratio = ((N-self.Nm) / (N+1))
        fi = (self.R_op-self.Rm) / (self.R_op+1)
        exp = math.exp(((1+54.4*fi)/(11+117.2*fi)) * ((fi-1)/(fi**0.5)))
        
        return N_ratio - 1 + exp
        
    """Uses Gilliland Eqn."""
    def EstimateNumTheoreticalStages(self):
        N = fsolve(self.Gilliland, [4, 4, 4])
        
        return N

    """ WIP: Kirkbride Approximation for Feed Tray Location"""
    def Kirkbride(self): #WIP
        Nf = 1
        
        return Nf

    # """WIP"""
    # def _mass_balance(self):
    #     vap, liq = self.outs
        
    #     # Get all important flow rates (both light and heavy keys and non-keys)
    #     LHK_index = self.LHK_index
    #     LNK_index = self.LNK_index
    #     HNK_index = self.HNK_index
    #     mol = self._mol_in
        
    #     LHK_mol = mol[LHK_index]
    #     LNK_mol = mol[LNK_index]
    #     HNK_mol = mol[HNK_index]
        
    #     # Set light and heavy keys by lever rule
    #     light, heavy = LHK_mol
    #     LHK_mol = light + heavy
    #     zf = light / LHK_mol
    #     split_frac = (zf-self.x_bot) / (self.y_top-self.x_bot)
    #     top_net = LHK_mol * split_frac
        
    #     # Set output streams
    #     vap.mol[LHK_index] = top_net * self._y
    #     liq.mol[LHK_index] = LHK_mol - vap.mol[LHK_index]
    #     vap.mol[LNK_index] = LNK_mol
    #     liq.mol[HNK_index] = HNK_mol
        
    """WIP"""
    def _run(self):
        ins_0 = self.ins[0]
        vap, liq = self.outs
        LHK_0 = self.LHK[0]
        
        rel_vols = self.GetRelativeVolatilities()
        Nm = self.EstimateMinStages()
        DB = self.EstimateDBCompositions()
        R_op = self.EstimateOperatingRefluxRatio()
        
        vap.copy_flow(ins_0)
        liq.copy_flow(ins_0)
        
        FLK = ins_0.mass[ins_0.chemicals.index(LHK_0)]
        F = self.ins[0].F_mass
        
        dLK = DB[vap.chemicals.index(LHK_0)][0]
        bLK = DB[vap.chemicals.index(LHK_0)][1]
        
        D = (FLK-F*bLK) / (dLK-bLK)
        B = F - D
        print(D)
        print(F)
        print(B)
        
        for spec in vap.chemicals.IDs:
            vap.imol[spec] = DB[vap.chemicals.index(spec)][0] * D
            liq.imol[spec] = DB[vap.chemicals.index(spec)][1] * B
'''


