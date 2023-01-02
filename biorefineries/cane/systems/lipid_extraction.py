# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import flexsolve as flx
import numpy as np
import biosteam as bst
from biosteam import SystemFactory
from .. import streams as s

from biorefineries.biodiesel import (
    create_lipid_wash_system,
)
__all__ = (
    'create_lipid_exctraction_system',
    'create_post_fermentation_oil_separation_system',
)

@SystemFactory(
    ID='lipid_exctraction_sys',
    ins=[s.fermentation_effluent],
    outs=[s.lipid, s.cellmass, s.wastewater],
)
def create_lipid_exctraction_system(ins, outs):
    fermentation_effluent, = ins
    lipid, cellmass, wastewater, = outs
    U401 = bst.SolidsCentrifuge('U401', fermentation_effluent, ['', ''],
        split=dict(cellmass=0.98, lipid=0.98),
        moisture_content=0.5,
        solids=['cellmass'],
    )
    U402 = bst.DrumDryer('U402', 
        (U401-0, 'dryer_air', 'dryer_natural_gas'), 
        ('', 'dryer_outlet_air', 'dryer_emissions'),
        moisture_content=0.18, split=0.,
    )
    # X401 = bst.ThermalOxidizer('X401', (U403-1, 'oxidizer_air'), 'oxidizer_emissions')
    U403 = bst.ScrewFeeder('U403', U402-0)
    U404 = bst.Splitter('U404', U403-0, 
        split=dict(cellmass=1, lipid=0.3, Water=0.8),
    )
    bst.ConveyingBelt('U405', U404-0, cellmass)
    lipid_wash_sys = create_lipid_wash_system(ins=U404-1, outs=lipid, mockup=True)
    washed_lipid, spent_wash_water = lipid_wash_sys.outs
    bst.Mixer(ins=[spent_wash_water, U401-1], outs=wastewater)

@SystemFactory(
    ID='post_fermentation_oil_separation_sys',
    ins=[s.stillage],
    outs=[s.lipid, s.wastewater, s.evaporator_condensate], 
)
def create_post_fermentation_oil_separation_system(ins, outs, wastewater_concentration=None,
                                                   target_oil_and_solids_content=60, 
                                                   separate_cellmass=False):
    lipid, wastewater, evaporator_condensate = outs
    if separate_cellmass:     
        cellmass = bst.Stream('cellmass')
        outs.insert(1, cellmass)
    V605 = bst.MixTank('V605', ins)
    P606 = bst.Pump('P606', V605-0)
    EvX = bst.MultiEffectEvaporator('Ev607',
        ins=P606-0,
        P=(101325, 69682, 47057, 30953),
        V=0.90, V_definition='First-effect',
        thermo=lipid.thermo.ideal(),
        flash=False,
    )
    EvX.target_oil_and_solids_content = target_oil_and_solids_content # kg / m3
    EvX.remove_evaporators = False
    P_original = tuple(EvX.P)
    Pstart = P_original[0]
    Plast = P_original[-1]
    N = len(P_original)
    def x_oil(V): # Objective function for specification
        EvX.V = V
        EvX.run()
        effluent = EvX.outs[0]
        moisture = effluent.imass['Water']
        total = effluent.F_mass
        return EvX.target_oil_and_solids_content - 1000. * (1. - moisture / total)
    
    @EvX.add_specification(run=False)
    def adjust_evaporation():
        V_last = EvX.V
        x0 = 0.
        x1 = 0.5
        EvX.P = P_original
        EvX._reload_components = True
        y0 = x_oil(x0)
        if y0 <= 0.:
            EvX.V = x0
            return
        elif EvX.remove_evaporators:
            EvX._load_components()
            for i in range(1, N):
                if x_oil(1e-6) < 0.:
                    EvX.P = np.linspace(Pstart, Plast, N - 1)
                    EvX._reload_components = True
                else:
                    break    
            y1 = x_oil(x1)
            EvX.V = flx.IQ_interpolation(x_oil, x0, x1, y0, y1, x=V_last, ytol=1e-5, xtol=1e-6)
        elif x_oil(1e-6) < 0.:
            EvX.V = 1e-6
        else:
            y1 = x_oil(x1)
            EvX.V = flx.IQ_interpolation(x_oil, 1e-6, x1, y0, y1, x=V_last, ytol=1e-5, xtol=1e-6)
        
    P607 = bst.Pump('P607', EvX-0, P=101325.)
    C603_2 = bst.LiquidsSplitCentrifuge('C603_2', P607-0, (lipid, ''), 
                                        split={'Oil': 0.99,
                                               'Water': 0.0001})
    if separate_cellmass:        
        C603_3 = bst.SolidsCentrifuge('C603_3', C603_2-1, (cellmass, ''), 
                                      split={'Cellmass': 0.99}, solids=('Cellmass',))
        stream = C603_3-1
    else:
        stream = C603_2-1
    S601 = bst.Splitter('S601', ins=EvX-1, outs=['', evaporator_condensate], split=0.5)
    M601 = bst.Mixer('M601', [S601-0, stream], wastewater)
    M601.target_wastewater_concentration = 60. # kg / m3
    @M601.add_specification(run=True)
    def adjust_wastewater_concentration():
        concentrated_wastewater = C603_2.outs[1]
        waste = concentrated_wastewater.F_mass - concentrated_wastewater.imass['Water'] 
        current_concentration = waste / concentrated_wastewater.F_vol
        required_water = (1./M601.target_wastewater_concentration - 1./current_concentration) * waste * 1000.
        F_mass = S601.ins[0].F_mass
        if F_mass:
            split = required_water / F_mass
            if split < 0:
                split = 0.
            elif split > 1.:
                split = 1.
            S601.split[:] = split
            for i in S601.path_until(M601): i.run()