#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import biosteam as bst
import flexsolve as flx
from .. import streams as s
from biosteam import SystemFactory

__all__ = (
    'create_sugar_crystallization_system',
)

@SystemFactory(
    ID='sugar_crystallization_sys',
    ins=[s.screened_juice, s.lime, s.H3PO4, s.polymer],
    outs=[s.sugar, s.molasses]
)
def create_sugar_crystallization_system(ins, outs):
    # TODO: Add conveyors, storage tanks, packing, sugar remelter
    # https://www.researchgate.net/profile/Maciej-Starzak/publication/311206128_Mass_and_Energy_Balance_Modelling_of_a_Sugar_Mill_A_comparison_of_MATLABR_and_SUGARS_simulations/links/583f240308ae2d217557dcd8/Mass-and-Energy-Balance-Modelling-of-a-Sugar-Mill-A-comparison-of-MATLABR-and-SUGARS-simulations.pdf?origin=publication_detail
    # https://www3.epa.gov/ttn/chief/ap42/ch09/final/c9s10-1a.pdf
    # http://sugartech.co.za/verticalcrystalliser/index.php
    screened_juice, lime, H3PO4, polymer = ins
    sugar, molasses, = outs
    
    if 'Sugar' not in sugar.chemicals:
        sugar.chemicals.define_group('Sugar', ('Glucose', 'Sucrose'))
    
    # Concentrate sugars
    P1 = bst.Pump('P1', ins=screened_juice, P=101325)
    
    MEE = bst.MultiEffectEvaporator('MEE', P1-0,
        P=(101325, 69682, 47057, 30953, 19781),
        V_definition='First-effect',
        V=0.3
    ) # fraction evaporated
    MEE.brix = 95
    
    def get_brix():
        effluent = MEE.outs[0]
        water = effluent.imass['Water']
        if water < 0.0001: water = 0.0001
        return 100 * effluent.imass['Sugar'] / water
    
    def brix_objective(V):
        MEE.V = V
        MEE._run()
        return MEE.brix - get_brix()
    
    @MEE.add_specification(run=False)
    def adjust_glucose_concentration():
        V_guess = MEE.V
        MEE.V = flx.IQ_interpolation(
            brix_objective, 0., 1., x=V_guess, ytol=1e-5
        )
    
    # Mix in flocculant
    T1 = bst.MixTank('T1', (MEE-0, lime, H3PO4, polymer))
    T1.tau = 0.10
    
    @T1.add_specification(run=True)
    def correct_flows():
        F_mass = T1.ins[0].F_mass
        # correct lime and phosphoric acid
        lime.imass['CaO', 'Water'] = 0.1 * F_mass * np.array([0.046, 0.954])
        H3PO4.imass['H3PO4', 'Water'] = 0.025 * F_mass
    
    # Separate residual solids
    C1 = bst.Clarifier('C1', T1-0, 
                           split=dict(Ash=0,
                                      Cellulose=0,
                                      Flocculant=1,
                                      Glucose=1,
                                      Hemicellulose=0,
                                      Lignin=0,
                                      CaO=1,
                                      H3PO4=1,
                                      Sucrose=1,
                                      Water=0.99))
    P2 = bst.Pump('P2', C1-0, P=101325)
    M1 = bst.Mixer('M1', (P2-0, '', ''))
    E1 = bst.Flash('E1', M1-0, V=0.5, P=15000)
    
    def get_purity(flash):
        effluent = flash.outs[1]
        return effluent.imass['Sugar'] / effluent.F_mass
    
    def purity_objective(V, flash):
        flash.V = V
        flash._run()
        return flash.purity - get_purity(flash)
    
    def adjust_purity(flash):
        V_guess = flash.V
        y0 = purity_objective(0., flash)
        if y0 < 0.: return
        y1 = purity_objective(1., flash)
        if y1 > 0.: return
        flash.V = flx.IQ_interpolation(
            purity_objective, 0., 1., y0, y1, x=V_guess, ytol=1e-5,
            args=(flash,),
        )
    
    E1.add_specification(adjust_purity, run=False, args=(E1,))
    E1.purity = 0.8623
    BC1 = bst.BatchCrystallizer('BC1', E1-1, tau=8, V=3785, T=55 + 273.15)
    
    def get_split(molasses_flow, molasses_purity, crystal_flow, crystal_purity):
        s_crystal = crystal_flow * crystal_purity
        s_molasses = molasses_flow * molasses_purity
        s_split = s_crystal / (s_crystal + s_molasses)
        o_crystal = crystal_flow * (100 - crystal_purity)
        o_molasses = molasses_flow * (100 - molasses_purity)
        o_split = o_crystal / (o_crystal + o_molasses)
        return dict(
            Water=o_split,
            H3PO4=o_split,
            CaO=o_split,
            Sugar=s_split,
        )
    
    C2 = bst.SolidsCentrifuge('C2', 
        BC1-0, 
        split=get_split(19.53, 62.91, 29.68, 98.61),
        moisture_content=None,
    )
    
    bst.StorageTank('S1', C2-0, sugar, tau=27 * 7)
    
    def correct_wash_water(mixer):
        mixer.ins[1].imass['Water'] = mixer.ins[0].imass['Sugar']
    
    M2 = bst.Mixer('M2', (C2-1, ''))
    M2.add_specification(correct_wash_water, run=True, args=(M2,))
    P3 = bst.Pump('P3', M2-0, P=101325)
    E2 = bst.Flash('E2', P3-0, V=0.5, P=15000)
    E2.add_specification(adjust_purity, run=False, args=(E2,))
    E2.purity = 0.6291
    BC2 = bst.BatchCrystallizer('BC2', E2-1, tau=24, V=3785, T=50 + 273.15)
    C3 = bst.SolidsCentrifuge('C3', 
        BC2-0, (2-M1, ''),
        split=get_split(4.34, 33.88, 3.15, 96.49),
        moisture_content=None,
    )
    M3 = bst.Mixer('M3', (C3-1, ''))
    M3.add_specification(correct_wash_water, run=True, args=(M3,))
    P4 = bst.Pump('P4', M3-0, P=101325)
    E3 = bst.Flash('E3', P4-0, V=0.5, P=15000)
    E3.add_specification(adjust_purity, run=False, args=(E3,))
    E3.purity = 0.5450
    BC3 = bst.BatchCrystallizer('BC3', E3-1, tau=40, V=3785, T=45 + 273.15)
    C4 = bst.SolidsCentrifuge('C4', 
        BC3-0, (1-M1, ''),
        split=get_split(9.04, 32.88, 4.48, 93.84),
        moisture_content=None,
    )
    bst.StorageTank('S2', C4-1, molasses, tau=24 * 7)
