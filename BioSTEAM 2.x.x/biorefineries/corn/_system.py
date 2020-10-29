# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from biorefineries import corn as cn
from biorefineries.sugarcane import create_ethanol_separation_units, price

__all__ = ('create_system',)

def create_system(ID='corn_sys'):
    ### Streams ###
    chemicals = bst.settings.get_chemicals()
    
    z_mass_corn = chemicals.kwarray(
        dict(Starch=0.595,
             Water=0.15,
             Fiber=0.1237,
             SolubleProtein=0.034,
             InsolubleProtein=0.0493,
             Oil=0.034,
             Ash=0.014)
    )

    corn_slurry_water = bst.Stream('corn_slurry_water', Water=83e3, units='l/hr')
    solids_content = 0.311
    F_mass_corn = corn_slurry_water.F_mass * solids_content / (1 - solids_content)
    mass_corn = F_mass_corn * z_mass_corn
    corn = bst.Stream('corn', flow=mass_corn, units='kg/hr')
    ammonia = bst.Stream('ammonia', NH3=89.723, units='kg/hr') # TODO: Update flows
    lime = bst.Stream('lime', CaO=53.609, units='kg/hr') # TODO: Update flows
    alpha_amylase = bst.Stream('alpha_amylase',
                               SolubleProtein=0.00082,
                               Water=0.99918,
                               units='kg/hr')
    alpha_amylase.F_mass = 32.467
    backwater = bst.Stream('backwater', Water=79135.6528, units='kg/hr')
    gluco_amylase = bst.Stream('gluco_amylase', 
                               SolubleProtein=0.0011,
                               Water=0.9989,
                               units='kg/hr')
    gluco_amylase.F_mass = 46.895
    sulfuric_acid = bst.Stream('sulfuric_acid', 
                               H2SO4=0.00061,
                               Water=0.99939,
                               units='kg/hr')
    sulfuric_acid.F_mass = 92.59300
    yeast = bst.Stream('yeast', Yeast=3.045, Water=63.335, units='kg/hr')
    scrubber_water = bst.Stream('scrubber_water', Water=20139.4417, units='kg/hr')
    ethanol = bst.Stream('ethanol', price=price['Ethanol'])
    denaturant = bst.Stream('denaturant', price=price['Gasoline'])
    
    u = cn.units
    MH101 = u.GrainHandling('MH101', corn)
    V102 = u.CornStorage('V102', MH101-0)
    MH103 = u.CleaningSystem('MH103', V102-0, split=0.997)
    M104 = u.HammerMill('M104', MH103-0)
    V105 = u.MilledCornSurgeTank('V105', M104-0)
    W106 = u.MilledCornHopper('W106', V105-0)
    V107 = u.MilledCornWeighTank('V107', W106-0)
    V301 = u.AlphaAmylaseTank('V301', alpha_amylase)
    P302 = bst.Pump('P302', V301-0)
    V303 = u.AmmoniaTank('V303', ammonia)
    P304 = bst.Pump('P304', V303-0)
    V305 = u.LimeHopper('V305', lime)
    V307 = u.SlurryMixTank('V307', (V107-0, P302-0, P304-0, V305-0, backwater))
    P308 = bst.Pump('P308', V307-0)
    HX101 = bst.HXutility('HX101', P308-0, U=1.5, T=87 + 273.15, ft=1.0)
    V310 = u.Liquefaction('V310', HX101-0)
    P311 = bst.Pump('P311', V310-0)
    E312 = bst.HXprocess('E312', (P311-0, None), U=0.56783, ft=1.0, T_lim1=97.8 + 273.15)
    E313 = u.JetCooker('E313', E312-0)
    V314 = u.CookedSlurrySurgeTank('V314', E313-0)
    V314-0-1-E312
    E315 = bst.HXutility('E315', E312-1, U=0.9937, ft=1.0, T=97.8 + 273.15)
    E316 = bst.HXprocess('E316', (E315-0, None), U=0.85174, ft=1.0, dT=12)
    V317 = u.GlucoAmylaseTank('V317', gluco_amylase)
    P318 = bst.Pump('P318', V317-0)
    V319 = u.SulfuricAcidTank('V319', sulfuric_acid)
    P320 = bst.Pump('P320', V319-0)
    V321 = u.Saccharification('V321', (P318-0, P320-0, E316-0))
    P322 = bst.Pump('P322', V321-0)
    E401 = bst.HXprocess('E401', (P322-0, None), T_lim0=47.5 + 273.15, U=0.99370, ft=0.85)
    E402 = bst.HXutility('E402', E401-0, ft=0.95, U=0.9937, T=32.2 + 273.15)
    V403 = u.YeastTank('V403', yeast)
    P404 = bst.Pump('P404', V403-0)
    V405 = u.SSF('V405', (E402-0, P404-0), outs=('CO2', ''))
    P406 = bst.Pump('P406', V405-0)
    P406-0-1-E401
    P407 = bst.Pump('P407', E401-1)
    P407-0-1-E316
    V412 = bst.Flash('V412', V=0.005, P=101325)
    E408 = bst.HXutility('E408', V412-0, V=0.99, rigorous=True)
    E408_2 = bst.PhaseSplitter('E408_2', E408-0)
    MX1 = bst.Mixer('MX1', (V405-0, E408_2-0))
    V409 = bst.VentScrubber('V409', (MX1-0, scrubber_water), gas=('CO2', 'O2'))
    P410 = bst.Pump('P410', V409-1)
    E413 = bst.HXprocess('E413', (V412-1, ''), U=0.79496, ft=1.0)
    MX2 = bst.Mixer('MX2', (E408_2-1, E413-0))
    P411 = bst.Pump('P411', MX2-0)
    
    dct = create_ethanol_separation_units(
        ID='ethanol_separation_sys',
        degassed_beer=P411-0,
        ethanol_product=ethanol,
        denaturant=denaturant,
        IDs={
            'Beer column': 'T501',
            'Beer column bottoms product pump': 'P502',
            'Recycle mixer': 'MX3',
            'Distillation': 'T503_T507',
            'Heat exchanger to superheat vapor to molecular sieves': 'HX500',
            'Molecular sieves': 'X504',
            'Ethanol condenser': 'HX501',
            'Ethanol day tank': 'V511', 
            'Ethanol day tank pump': 'P304',
            'Denaturant storage': 'T303', 
            'Denaturant pump': 'P305', 
            'Ethanol-denaturant mixer': 'MX4',
            'Product tank': 'V513',
            'Distillation bottoms product pump': 'P508',
        }
    )
    P502 = dct['Beer column bottoms product pump']
    T507 = dct['Distillation bottoms product pump']
    P502-0-1-E413
    V601 = bst.MixTank('V601', P502-0)
    P602 = bst.Pump('P602', V601-0)
    C603 = u.DDGSCentrifuge('C603', P602-0,
        split=dict(
            Water=0.8285,
            Ethanol=0.8285,
            Yeast=0.2734,
            Lipid=0.331,
            H2SO4=0.8285,
            Starch=0.08,
            Fiber=0.5328,
            SolubleProtein=0.8285,
            InsolubleProtein=0.08)
    )
    MH604 = u.WetDDGSConveyor('MH604', C603-1)
    S1 = bst.Splitter('S1', C603-0, split=0.208)
    V605 = bst.MixTank('V605', S1-1)
    P606 = bst.Pump('P606', V605-0)
    Ev607 = bst.MultiEffectEvaporator('Ev607',
        ins=P606-0,
        P=(101325, 73581, 50892, 32777),
        V=0.65
    ) 
    Ev607.components['condenser'].U = 1.85
    MX5 = bst.Mixer('MX5', (Ev607-1, P410-0, T507-0))
    MX6 = bst.Mixer('MX6', Ev607-0, MH604-0)
    D610 = u.DDGSDryer('D610', (MX6-0, 'air', 'natural_gas'), moisture_content=0.10, split=dict(Ethanol=1.0))
    X611 = u.ThermalOxidizer('X611', D610-1)
    MH612 = u.DDGSHandling('MH612', X611-0)
    T608 = bst.facilities.ProcessWaterCenter('P608', (MX5-0, 'makeup_water'), (V307-0, 'wastewater'))
    