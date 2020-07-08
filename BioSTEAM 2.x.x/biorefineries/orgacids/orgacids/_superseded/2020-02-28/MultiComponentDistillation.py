#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 16:10:15 2020

@author: yalinli_cabbi

Based on Sarang's Distillation_Multicomponent moudule as of 2020-02-04,
    but compatible with biosteam v2.0.1
"""

import math
import numpy as np
import thermosteam as tmo
from scipy.optimize import fsolve
from biosteam import Unit

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
        equilibrium_chemicals = tmo.Chemicals(['Water', 'Ethanol'])
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
        
        for spec in vap.chemicals.IDs:
            vap.imol[spec] = DB[vap.chemicals.index(spec)][0] * D
            liq.imol[spec] = DB[vap.chemicals.index(spec)][1] * B

