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
from biorefineries.sugarcane import create_ethanol_purification_system
from ._process_settings import price

__all__ = ('create_system',)

import warnings
import numpy as np
warnings.filterwarnings("error", category=np.ComplexWarning)

SLURRY_SOLIDS_CONTENT = 0.311   # g Corn / g Total
SLURRY_AMMONIA_LOADING = 0.002  # g Ammonia / g dry Corn
SLURRY_LIME_LOADING = 0.00012   # g Lime / g dry Corn
LIQUEFACTION_ALPHA_AMYLASE_LOADING = 0.0007    # g Enzyme / g dry Corn
SACCHARIFICATION_SULFURIC_ACID_LOADING = 0.001 # g Enzyme / g dry Corn
SACCHARIFICATION_GLUCO_AMYLASE_LOADING = 0.002 # g H2SO4 / g dry Corn
SCRUBBER_WASH_WATER_OVER_VENT = 1.21 # g Water / g Vent

def create_system(ID='corn_sys'):
    ### Streams ###
    
    chemicals = bst.settings.get_chemicals()
    
    z_mass_corn = chemicals.kwarray(
        dict(Starch=0.612,
             Water=0.15,
             Fiber=0.1067,
             SolubleProtein=0.034,
             InsolubleProtein=0.0493,
             Oil=0.034,
             Ash=0.014)
    )

    F_mass_corn = 46211.6723 # kg / hr
    mass_corn = F_mass_corn * z_mass_corn

    corn = bst.Stream('corn', flow=mass_corn, price=price['Corn'], units='kg/hr')
    crude_oil = bst.Stream('crude_oil', price=price['Crude oil'], units='kg/hr')
    ammonia = bst.Stream('ammonia', NH3=89.723, price=price['Ammonia'], units='kg/hr')
    lime = bst.Stream('lime', CaO=53.609, price=price['Lime'], units='kg/hr') 
    alpha_amylase = bst.Stream('alpha_amylase',
                               SolubleProtein=0.00082,
                               Water=0.99918,
                               price=price['Enzyme'],
                               units='kg/hr')
    alpha_amylase.F_mass = 32.467
    recycled_process_water = bst.Stream('recycled_process_water', 
                                        Water=1, 
                                        units='kg/hr')
    gluco_amylase = bst.Stream('gluco_amylase', 
                               SolubleProtein=0.0011,
                               Water=0.9989,
                               price=price['Enzyme'],
                               units='kg/hr')
    gluco_amylase.F_mass = 46.895
    sulfuric_acid = bst.Stream('sulfuric_acid', 
                               H2SO4=1.,
                               price=price['Sulfuric acid'],
                               units='kg/hr')
    sulfuric_acid.F_mass = 92.59300
    backwater = bst.Stream('backwater', Water=22484., units='kg/hr')
    yeast = bst.Stream('yeast', Yeast=3.045, Water=63.335, price=price['Yeast'], units='kg/hr')
    scrubber_water = bst.Stream('scrubber_water', Water=20139.4417, units='kg/hr')
    ethanol = bst.Stream('ethanol', price=price['Ethanol'])
    denaturant = bst.Stream('denaturant', price=price['Denaturant'])
    high_pressure_steam = bst.HeatUtility.get_agent('high_pressure_steam')
    DDGS = bst.Stream('DDGS', price=price['DDGS'])
    steam = bst.Stream('steam', Water=1, phase='g', 
                       price=price['Steam'],
                       T=high_pressure_steam.T,
                       P=high_pressure_steam.P)
    
    ### Process specifications ###
    
    def refresh_feed_specifications():
        MH101._F_mass_dry_corn = F_mass_dry_corn = corn.F_mass - corn.imass['Water']
        lime.F_mass = F_mass_lime = F_mass_dry_corn * SLURRY_LIME_LOADING
        ammonia.F_mass = F_mass_ammonia = F_mass_dry_corn * SLURRY_AMMONIA_LOADING
        alpha_amylase.F_mass = F_mass_aa = F_mass_dry_corn * LIQUEFACTION_ALPHA_AMYLASE_LOADING
        sulfuric_acid.F_mass = F_mass_sa = F_mass_dry_corn * SACCHARIFICATION_SULFURIC_ACID_LOADING
        gluco_amylase.F_mass = F_mass_ga = F_mass_dry_corn * SACCHARIFICATION_GLUCO_AMYLASE_LOADING
        MH101._F_mass_others = (F_mass_lime + F_mass_ammonia + F_mass_aa + F_mass_sa + F_mass_ga)
        MH101._run()
    
    def update_scrubber_wash_water():
        scrubber_water.F_mass =  V409.ins[1].F_mass * SCRUBBER_WASH_WATER_OVER_VENT
        V409._run()
        
    ### Units ###
    
    u = cn.units
    MH101 = u.GrainHandling('MH101', corn)
    MH101.specification = refresh_feed_specifications
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
    V307 = u.SlurryMixTank('V307', (V107-0, P302-0, P304-0, V305-0, recycled_process_water, backwater))
    @V307.add_specification(run=True)
    def correct_recycle_dilution_water():
        F_mass_others = MH101._F_mass_others + backwater.F_mass
        F_mass_dry_corn = MH101._F_mass_dry_corn
        recycled_process_water.F_mass = F_mass_dry_corn * (1. - SLURRY_SOLIDS_CONTENT) / SLURRY_SOLIDS_CONTENT - F_mass_others
    
    P311 = bst.Pump('P311', V307-0, P=1e6)
    E312 = bst.HXprocess('E312', (P311-0, None), U=0.56783, ft=1.0, T_lim0=410.)
    E313 = u.JetCooker('E313', (E312-0, steam))
    V314 = u.CookedSlurrySurgeTank('V314', E313-0)
    P308 = bst.Pump('P308', V314-0)
    P308-0-1-E312
    # E315 = bst.HXutility('E315', E312-1, U=0.9937, ft=1.0, T= + 273.15)
    
    HX101 = bst.HXutility('HX101', E312-1, U=1.5, T=87 + 273.15, ft=1.0)
    V310 = u.Liquefaction('V310', HX101-0)
    E316 = bst.HXprocess('E316', (V310-0, None), U=0.85174, ft=1.0, dT=12)
    
    V317 = u.GlucoAmylaseTank('V317', gluco_amylase)
    P318 = bst.Pump('P318', V317-0)
    V319 = u.SulfuricAcidTank('V319', sulfuric_acid)
    P320 = bst.Pump('P320', V319-0)
    V321 = u.Saccharification('V321', (P318-0, P320-0, E316-0))
    P322 = bst.Pump('P322', V321-0)
    E401 = bst.HXprocess('E401', (P322-0, None),
                         phase0='l', phase1='l',
                         U=0.99370, ft=0.85)
    E402 = bst.HXutility('E402', E401-0, ft=0.95, U=0.9937, T=32.2 + 273.15)
    V403 = u.YeastTank('V403', yeast)
    P404 = bst.Pump('P404', V403-0)
    V405 = u.SSF('V405', (E402-0, P404-0), outs=('CO2', ''), V=1.9e3)
    P406 = bst.Pump('P406', V405-1)
    P407 = bst.Pump('P407', E401-1)
    P407-0-1-E316
    V412 = bst.StorageTank('V412', E316-1, tau=4)
    PX = bst.Pump('PX', V412-0)
    V409 = bst.VentScrubber('V409', (scrubber_water, V405-0), gas=('CO2', 'O2'))
    V409.specification = update_scrubber_wash_water
    P410 = bst.Pump('P410', V409-1)
    MX = bst.Mixer('MX', [P410-0, P406-0])
    MX-0-1-E401
    E413 = bst.HXprocess('E413', (PX-0, None), U=0.79496, ft=1.0)
    P411 = bst.Pump('P411', E413-0)
    
    ethanol_purification_sys = create_ethanol_purification_system(
        ins=[P411-0, denaturant],
        outs=ethanol,
        beer_column_heat_integration=False,
        mockup=True,
        IDs={
            'Beer column': 'T501',
            'Beer column bottoms product pump': 'P502',
            'Recycle mixer': 'MX3',
            'Distillation': 'T503_T507',
            'Heat exchanger to superheat vapor to molecular sieves': 'HX500',
            'Molecular sieves': 'X504',
            'Ethanol condenser': 'HX501',
            'Ethanol day tank': 'V511', 
            'Ethanol day tank pump': 'P512',
            'Denaturant storage': 'V509', 
            'Denaturant pump': 'P510', 
            'Ethanol-denaturant mixer': 'MX4',
            'Product tank': 'V513',
            'Distillation bottoms product pump': 'P508',
        }
    )
    f = cn.flowsheet
    fu = f.unit
    fu.X504.approx_duty = False
    fu.T501.Rmin = 0.0001
    fu.T503_T507.k = 1.05
    fu.T501.P = 101325
    fu.P502-0-1-E413
    fu.MX4.denaturant_fraction = 0.04345
    V601 = bst.MixTank('V601', E413-1)
    P602 = bst.Pump('P602', V601-0)
    C603 = u.DDGSCentrifuge('C603', P602-0,
        split=dict(
            Water=0.8285,
            Ethanol=0.8285,
            Yeast=0.2734,
            Lipid=0.5, # 0.331 originally
            H2SO4=0.8285,
            Starch=0.08,
            Fiber=0.5328,
            SolubleProtein=0.8285,
            InsolubleProtein=0.08)
    )
    MH604 = u.WetDDGSConveyor('MH604', C603-1)
    S1 = bst.Splitter('S1', C603-0, (backwater, ''), split=0.208)
    V605 = bst.MixTank('V605', S1-1)
    P606 = bst.Pump('P606', V605-0)
    Ev607 = bst.MultiEffectEvaporator('Ev607',
        ins=P606-0,
        P=(101325, 69682, 47057, 30953),
        V=0.90
    ) 
    C603_2 = bst.LiquidsSplitCentrifuge('C603_2', Ev607-0, (crude_oil, ''), split={'Oil':0.99})
    
    MX5 = bst.Mixer('MX5', (Ev607-1, P410-0, fu.P508-0))
    MX6 = bst.Mixer('MX6', (C603_2-1, MH604-0))
    D610 = bst.DrumDryer('D610', (MX6-0, 'dryer_air', 'natural_gas'), moisture_content=0.10, split=dict(Ethanol=1.0))
    X611 = bst.ThermalOxidizer('X611', (D610-1, 'oxidizer_air'))
    MH612 = u.DDGSHandling('MH612', D610-0, DDGS)
    T608 = bst.facilities.ProcessWaterCenter(
        'T608', 
        (MX5-0, 'makeup_water'), 
        ('process_water', 'wastewater'),
        None,
        (),
        (recycled_process_water,)
    )
    other_facilities = u.PlantAir_CIP_WasteWater_Facilities('other_facilities', corn)
    HXN = bst.HeatExchangerNetwork('HXN', 
        units=lambda: [fu.Ev607.components['evaporators'][0], fu.T503_T507.condenser]
    )
    globals().update(f.unit.data)
    return f.create_system('corn_sys', feeds=[i for i in f.stream if i.isfeed()])
    
    # System = bst.System
    
    # return System('corn_sys',
    #     [MH101,
    #      V102,
    #      MH103,
    #      M104,
    #      V105,
    #      W106,
    #      V107,
    #      V303,
    #      P304,
    #      V305,
    #      V301,
    #      P302,
    #      V317,
    #      P318,
    #      V319,
    #      P320,
    #      V403,
    #      P404,
    #      System('SYS5',
    #         [V307,
    #          P311,
    #          System('SYS1',
    #             [E312,
    #              E313,
    #              V314,
    #              P308],
    #             recycle=P308-0),
    #          HX101,
    #          V310,
    #          System('SYS3',
    #             [E316,
    #              V321,
    #              P322,
    #              System('SYS2',
    #                 [E401,
    #                  E402,
    #                  V405,
    #                  P406],
    #                 recycle=P406-0),
    #              P407],
    #             recycle=P407-0),
    #          V412,
    #          System('SYS4',
    #             [E413,
    #              MX2,
    #              P411,
    #              P301,
    #              T501,
    #              P502],
    #             recycle=P502-0),
    #          V601,
    #          P602,
    #          C603,
    #          S1],
    #         recycle=S1-0),
    #      E408,
    #      E408_2,
    #      System('SYS6',
    #         [MX3,
    #          T503_T507,
    #          HX500,
    #          X504],
    #         recycle=X504-0),
    #      HX501,
    #      V511,
    #      P512,
    #      V509,
    #      P510,
    #      MX4,
    #      V513,
    #      V605,
    #      P606,
    #      Ev607,
    #      C603_2,
    #      MX6,
    #      D610,
    #      X611,
    #      MH612,
    #      P508,
    #      MX5,
    #      MX1,
    #      V409,
    #      P410,
    #      MH604],
    #     facilities=[other_facilities,
    #      T608,
    #      HXN])
                        