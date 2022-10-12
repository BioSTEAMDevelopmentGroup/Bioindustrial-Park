# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import os
import biosteam as bst
import thermosteam as tmo
from biosteam import units
from biosteam import main_flowsheet as F
from biorefineries.LAOs import units as LAOs_units
from biorefineries import LAOs
from biorefineries.fattyalcohols import create_fattyalcohol_production_sys
from thermosteam.reaction import (Reaction as Rxn,
                                  ParallelReaction as ParallelRxn,
)

__all__ = ('create_system',)

# %% Reactions

def create_system(ID='LAOs_sys', stainless_steel=True):
    checks = {'check_atomic_balance': True}
    # LAOs production
    LAOs_production =  ParallelRxn([
        Rxn('Hexanol -> Hexene + H2O', 'Hexanol', 1, **checks),
        Rxn('Octanol -> Octene + H2O', 'Octanol', 1, **checks),
        Rxn('Decanol -> Decene + H2O', 'Decanol', 1, **checks),
        Rxn('Dodecanol -> Dodecene + H2O', 'Dodecanol', 1, **checks),
        Rxn('Tetradecanol -> Tetradecene + H2O', 'Tetradecanol', 1, **checks),
        Rxn('Hexadecanol -> Hexadecene + H2O', 'Hexadecanol', 1, **checks),
    ])
    
    stream_data_path = os.path.join(os.path.dirname(__file__), 'streams.yaml')
    stream_data = tmo.ThermoData.from_yaml(stream_data_path)
    hexene, octene, decene = stream_data.create_streams(['hexene', 'octene', 'decene'])
    
    fattyalcohol_production_sys = create_fattyalcohol_production_sys()
    F.unit.C101.liquids_isplit['Hexene', 'Octene', 'Decene', 'Dodecene', 'Tetradecene'] = 1.0
    H102 = units.HXprocess('H102', ins=F.unit.P107-0, dT=10, phase0='g', phase1='l')
    M103 = units.Mixer('M103', ins=(H102-0, None))
    H103 = units.HXutility('H103', ins=M103-0, V=1, T=623.15)
    R102 = LAOs_units.AdiabaticFixedbedGasReactor('R102',
        ins=H103-0,
        WHSV=1.,
        vessel_type='Vertical',
        P=5 * 101325, 
        reaction=LAOs_production,
        length_to_diameter=2,
        catalyst_price=9.69,
        catalyst_density=2650
    )
    R102-0-1-H102
    R102.dehydration_reactor_mass_fraction = 0.10 # A proces specification
    H105 = units.HXutility('H105', ins=H102-1, T=320)
    T107 = LAOs_units.SurgeTank('T107', H105-0, ('recycle_nitrogen', 'LAOs_to_separations'),
                                  tau=0.1)
    P105 = units.Pump('P105', T107-1)
    C102 = units.LiquidsSplitCentrifuge('C102', P105-0,
                               split=dict(Water=0.1,
                                          Tridecane=1.0,
                                          Hexene=1.0,
                                          Octene=1.0,
                                          Decene=1.0,
                                          Dodecene=1.0,
                                          Tetradecene=1.0))
    P106 = units.Pump('P106', C102-0)
    H104 = units.HXprocess('H104', ins=(P106-0, None), phase0='l', phase1='l')
    
    # Lab mass fraction is 0.0242 g / g
    # Baseline mass fraction is 0.1 g /g
    def adjust_nitrogen_flow():
        T107._run()
        recycle_nitrogen = T107.outs[0]
        feed = H102.outs[0]
        F_mass = feed.F_mass
        F_mass_alcohols = feed.imass[LAOs.fermentation_products].sum()
        x_alcohol = R102.dehydration_reactor_mass_fraction
        # Math to get amount of nitrogen flow:
        # x_alcohol = F_mass_alcohols / (F_mass + F_mass_N2)
        # F_mass_N2 * x_alcohol + x_alcohol * F_mass = F_mass_alcohols
        # F_mass_N2 * x_alcohol = F_mass_alcohols - x_alcohol * F_mass 
        F_mass_N2 = F_mass_alcohols / x_alcohol - F_mass 
        if F_mass_N2 < 0: F_mass_N2 = 0
        recycle_nitrogen.imass['Nitrogen'] = F_mass_N2
        recycle_nitrogen.phase = 'g'
    
    T107.add_specification(adjust_nitrogen_flow)
    T107 - 0 - 1 - M103
    D101 = units.ShortcutColumn('D101', H104-0, LHK=('Decene', 'Tridecane'), 
                                k=1.05, Lr=0.999, Hr=0.90, is_divided=True)
    
    def distillate_recoveries_hook(IDs, recoveries):
        light_keys = ('Water', 'Hexene', 'Octene')
        index = [n for n, i in enumerate(IDs) if i in light_keys]
        recoveries[index] = 1.0
    D102 = units.ShortcutColumn('D102', D101-0, LHK=('Decene', 'Tridecane'), 
                                k=1.05, y_top=0.98, x_bot=1e-9, is_divided=True)
    D102._distillate_recoveries_hook = distillate_recoveries_hook
    M105 = units.Mixer('M105', (D101-1, D102-1))
    P108 = units.Pump('P108', M105-0)
    P108-0-1-H104
    settler_data = {'split': dict(Water=0.000227,
                                  Tridecane=0.999999,
                                  Hexanol=0.342,
                                  Octanol=0.859,
                                  Decanol=0.932,
                                  Dodecanol=0.999,
                                  Hexene=1.0,
                                  Octene=1.0,
                                  Decene=1.0,
                                  Dodecene=1.0,
                                  Tetradecene=1.0,
                                  Tetradecanol=1.0)
    }
    M104 = units.MixerSettler('M104', (H104-1, F.unit.C101-1, C102-1, ''), ('', 'wastewater'),
                              settler_data=settler_data, model='split')
    
    def adjust_price_per_organics():
        M104._run()
        wastewater = M104.outs[1]
        wastewater.price = -0.33 * (1 - wastewater.imass['Water', 'NaCl'].sum() / wastewater.F_mass)
    
    M104.add_specification(adjust_price_per_organics)
    M104.outs[1].price = -0.33 # USD / kg organics
    P109 = units.Pump('P109', M104-0, 1**F.unit.M101)
    H109 = units.HXutility('H109', D102-0, T=320, V=0)
    C103 = units.LiquidsSplitCentrifuge('C103', H109-0, outs=('', 3**M104),
                               split=dict(Water=0.01,
                                          Hexene=0.995,
                                          Octene=0.999,
                                          Decene=1.0,
                                          Dodecene=1.0,
                                          Tetradecene=1.0))
    
    D103 = units.ShortcutColumn('D103', LHK=('Hexene', 'Octene'),
                                k=1.05, y_top=0.98, x_bot=1e-6)
    H106 = units.HXutility('H106', D103-0, T=320, V=0)
    D104 = units.BinaryDistillation('D104', D103-1, LHK=('Octene', 'Decene'),
                                k=1.05, y_top=0.98, x_bot=1e-6)
    H107 = units.HXutility('H107', D104-0, T=320, V=0)
    H108 = units.HXprocess('H108', (C103-0, D104-1), dT=5, phase1='l')
    H108-0-D103
    
    T108 = units.StorageTank('T108', H106-0, hexene,
                              vessel_type='Floating roof',
                              tau=7*24)
    T109 = units.StorageTank('T109', H107-0, octene,
                              vessel_type='Floating roof',
                              tau=7*24)
    T110 = units.StorageTank('T110', H108-1, decene,
                              vessel_type='Floating roof',
                              tau=7*24)
    
    ### Facilities ###
    
    *other_agents, high_pressure_steam = bst.HeatUtility.heating_agents
    BT = bst.BoilerTurbogenerator('BT', 
                                  ins=(None, None, 'boiler_makeup_water', 
                                       'natural_gas', 'lime', 'boiler_chemicals'),
                                  boiler_efficiency = 0.80,
                                  turbogenerator_efficiency=0.85,
                                  agent=high_pressure_steam,
                                  other_agents=other_agents)
    CT = bst.CoolingTower('CT')
    CWP = bst.ChilledWaterPackage('CWP')
    CCI = bst.ChemicalCapitalInvestment('CCI', 'Tridecane', 644.05)
    
    ### System ###
    
    sys = bst.System('LAOs_sys',
        [bst.System('tridecane_recycle',
            [fattyalcohol_production_sys,
             bst.System('reaction_heat_integration_sys',
                [H102,
                 H105,
                 T107,
                 M103,
                 H103,
                 R102],
                recycle=R102-0),
             P105,
             C102,
             P106,
             bst.System('distillation_heat_integration_sys',
                [H104,
                 D101,
                 D102,
                 C103,
                 H109,
                 M104,
                 M105,
                 P108],
                recycle=P108-0),
             P109],
            recycle=1**F.unit.M101),
         H108,
         D103,
         D104,
         H106,
         H107,
         H108,
         T108,
         T109,
         T110],
        facilities=[CWP, BT, CT, BT, CCI])

    if stainless_steel:
        isa = isinstance
        for i in sys.units:
            if isa(i, bst.HX):
                i.material = 'Stainless steel/stainless steel'
            elif isa(i, bst.BinaryDistillation):
                i.boiler.material = i.condenser.material = 'Stainless steel/stainless steel'
                i.vessel_material = 'Stainless steel 304'
                i.tray_material = 'Stainless steel 304'
            elif isa(i, bst.Pump):
                i.material = 'Stainless steel'
    return sys
