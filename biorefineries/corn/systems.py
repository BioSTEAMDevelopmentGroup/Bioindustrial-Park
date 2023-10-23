# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
#                     Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autofunction:: biorefineries.corn.systems.create_system

"""
import biosteam as bst
from biorefineries.ethanol import create_ethanol_purification_system
from .process_settings import BiorefinerySettings
from . import units

__all__ = ('create_system',)

import warnings
import numpy as np
warnings.filterwarnings("error", category=np.ComplexWarning)

def create_system(flowsheet=None, biorefinery_settings=None):
    settings = biorefinery_settings or BiorefinerySettings()
    system_ID = settings.system_ID
    feedstock_ID = 'corn' if 'corn' in system_ID else 'feedstock'
    parameters = settings.process_parameters
    prices = settings.stream_prices
    GWP_CFs = settings.stream_GWP_CFs
    
    get_stream_kwargs = lambda key: {
        'price': prices.get(key, 0.),
        'characterization_factors': {'GWP': GWP_CFs.get(key, 0.)},
        'units': 'kg/hr',
        }
    
    ### Streams ###
    
    chemicals = bst.settings.get_chemicals()
    
    z_mass_corn = chemicals.kwarray(settings.feedstock_composition)

    F_mass_corn = settings.feedstock_hourly_mass_flow
    mass_corn = F_mass_corn * z_mass_corn

    feedstock = bst.Stream(feedstock_ID, flow=mass_corn, **get_stream_kwargs(feedstock_ID.capitalize()))
    feedstock.register_alias('feedstock')
    crude_oil = bst.Stream('crude_oil', **get_stream_kwargs('Crude oil'))
    ammonia = bst.Stream('ammonia', NH3=89.723, **get_stream_kwargs('Ammonia'))
    lime = bst.Stream('lime', CaO=53.609, **get_stream_kwargs('Lime'))

    aa_price = prices.get('Alpha Amylase') or prices.get('Enzyme', 0.)
    aa_CF = GWP_CFs.get('Alpha Amylase') or prices.get('Enzyme', 0.)
    aa_kwargs = {
        'price': aa_price,
        'characterization_factors': {'GWP': aa_CF},
        'units': 'kg/hr',
        }

    alpha_amylase = bst.Stream('alpha_amylase',
                               SolubleProtein=0.00082,
                               Water=0.99918,
                               **aa_kwargs)
    # F_mass for this (and other feedstock-related streams)
    # will be adjusted through specification
    alpha_amylase.F_mass = 32.467
    recycled_process_water = bst.Stream('recycled_process_water', 
                                        Water=1, 
                                        units='kg/hr')
    ga_price = prices.get('Gluco Amylase') or prices.get('Enzyme', 0.)
    ga_CF = GWP_CFs.get('Gluco Amylase') or prices.get('Enzyme', 0.)
    ga_kwargs = {
        'price': ga_price,
        'characterization_factors': {'GWP': ga_CF},
        'units': 'kg/hr',
        }

    gluco_amylase = bst.Stream('gluco_amylase', 
                               SolubleProtein=0.0011,
                               Water=0.9989,
                               **ga_kwargs)
    gluco_amylase.F_mass = 46.895
    sulfuric_acid = bst.Stream('sulfuric_acid', 
                               H2SO4=1.,
                               **get_stream_kwargs('Sulfuric acid'))
    sulfuric_acid.F_mass = 92.59300
    backwater = bst.Stream('backwater', Water=22484., units='kg/hr')
    yeast = bst.Stream('yeast', Yeast=3.045, Water=63.335, **get_stream_kwargs('Yeast'))
    scrubber_water = bst.Stream('scrubber_water', Water=20139.4417, units='kg/hr')
    ethanol = bst.Stream('ethanol', **get_stream_kwargs('Ethanol'))
    denaturant = bst.Stream('denaturant', **get_stream_kwargs('Denaturant'))
    high_pressure_steam = bst.HeatUtility.get_agent('high_pressure_steam')
    DDGS = bst.Stream('DDGS', **get_stream_kwargs('DDGS'))
    steam = bst.Stream('steam', Water=1, phase='g',
                       T=high_pressure_steam.T,
                       P=high_pressure_steam.P,
                       **get_stream_kwargs('Steam'))
    
    ### Process specifications ###
    
    def refresh_feed_specifications():
        F_mass_dry_corn = feedstock.F_mass - feedstock.imass['Water']
        lime.F_mass = F_mass_dry_corn * parameters['slurry_lime_loading'] 
        ammonia.F_mass = F_mass_dry_corn * parameters['slurry_ammonia_loading']
        alpha_amylase.F_mass = F_mass_dry_corn * parameters['liquefaction_alpha_amylase_loading']
        sulfuric_acid.F_mass = F_mass_dry_corn * parameters['saccharification_sulfuric_acid_loading']
        MH101._run()
    
    def update_scrubber_wash_water():
        scrubber_water.F_mass =  V409.ins[1].F_mass * parameters['scrubber_wash_water_over_vent']
        V409._run()
        
    ### Units ###
    
    MH101 = units.GrainHandling('MH101', feedstock)
    MH101.add_specification(refresh_feed_specifications)
    V102 = units.CornStorage('V102', MH101-0)
    MH103 = units.CleaningSystem('MH103', V102-0, split=0.997)
    M104 = units.HammerMill('M104', MH103-0)
    V105 = units.MilledCornSurgeTank('V105', M104-0)
    W106 = units.MilledCornHopper('W106', V105-0)
    V107 = units.MilledCornWeighTank('V107', W106-0)
    V301 = units.AlphaAmylaseTank('V301', alpha_amylase)
    P302 = bst.Pump('P302', V301-0)
    V303 = units.AmmoniaTank('V303', ammonia)
    P304 = bst.Pump('P304', V303-0)
    V305 = units.LimeHopper('V305', lime)
    V307 = units.SlurryMixTank('V307', (V107-0, P302-0, P304-0, V305-0, recycled_process_water, backwater))
    slurry_solids_content = parameters['slurry_solids_content']
    @V307.add_specification(run=True)
    def correct_recycle_dilution_water():
        F_mass_dry_corn = feedstock.F_mass - feedstock.imass['Water']
        F_mass_others = lime.F_mass + ammonia.F_mass + alpha_amylase.F_mass + sulfuric_acid.F_mass + gluco_amylase.F_mass
        recycled_process_water.F_mass = F_mass_dry_corn * (1. - slurry_solids_content) / slurry_solids_content - F_mass_others
    
    P311 = bst.Pump('P311', V307-0, P=2e6)
    E312 = bst.HXprocess('E312', (P311-0, None), U=0.56783, ft=1.0)
    E313 = units.JetCooker('E313', (E312-0, steam))
    V314 = units.CookedSlurrySurgeTank('V314', E313-0)
    P308 = bst.Pump('P308', V314-0)
    P308-0-1-E312
    # E315 = bst.HXutility('E315', E312-1, U=0.9937, ft=1.0, T= + 273.15)
    
    HX101 = bst.HXutility('HX101', E312-1, U=1.5, T=87 + 273.15, ft=1.0)
    V310 = units.Liquefaction('V310', HX101-0)
    E316 = bst.HXprocess('E316', (V310-0, None), U=0.85174, ft=1.0, dT=5)
    
    V317 = units.GlucoAmylaseTank('V317', gluco_amylase)
    P318 = bst.Pump('P318', V317-0)
    V319 = units.SulfuricAcidTank('V319', sulfuric_acid)
    P320 = bst.Pump('P320', V319-0)
    V321 = units.Saccharification('V321', (P318-0, P320-0, E316-0))
    P322 = bst.Pump('P322', V321-0)
    E401 = bst.HXprocess('E401', (P322-0, None),
                         phase0='l', phase1='l',
                         U=0.99370, ft=0.85)
    E402 = bst.HXutility('E402', E401-0, ft=0.95, U=0.9937, T=32.2 + 273.15)
    V403 = units.YeastTank('V403', yeast)
    P404 = bst.Pump('P404', V403-0)
    V405 = units.SSF('V405', (E402-0, P404-0), outs=('CO2', ''), V=1.9e3)
    
    @V405.add_specification(run=True)
    def correct_saccharification_feed_flows():
        mash = V405.ins[0]
        mash_flow = mash.F_mass
        mash_dry_flow = mash_flow - mash.imass['Water']
        yeast.F_mass = parameters['yeast_loading'] * mash_flow
        gluco_amylase.F_mass = parameters['saccharification_gluco_amylase_loading'] * mash_dry_flow
    
    P406 = bst.Pump('P406', V405-1)
    P407 = bst.Pump('P407', E401-1)
    P407-0-1-E316
    V412 = bst.StorageTank('V412', E316-1, tau=4)
    PX = bst.Pump('PX', V412-0)
    V409 = bst.VentScrubber('V409', (scrubber_water, V405-0), gas=('CO2', 'O2'))
    V409.add_specification(update_scrubber_wash_water)
    P410 = bst.Pump('P410', V409-1)
    MX = bst.Mixer('MX', P406-0)
    MX-0-1-E401
    
    ethanol_purification_sys = create_ethanol_purification_system(
        ins=[PX-0, denaturant],
        outs=ethanol,
        beer_column_heat_integration=True,
        mockup=True,
        IDs={
            'Beer column': 'T501',
            'Beer column heat exchange': 'E413',
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
    f = flowsheet or bst.main_flowsheet
    u = f.unit
    u.X504.approx_duty = False
    u.T501.Rmin = 0.0001
    u.T501.k = 1.05
    u.T503_T507.k = 1.05
    u.T501.P = 101325
    u.T503_T507.P = 101325
    u.E413.U = 0.79496
    u.E413.ft = 1.0
    u.MX4.denaturant_fraction = 0.04345
    V601 = bst.MixTank('V601', u.E413-1)
    P602 = bst.Pump('P602', V601-0)
    C603 = units.DDGSCentrifuge('C603', P602-0,
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
    MH604 = units.WetDDGSConveyor('MH604', C603-1)
    S1 = bst.Splitter('S1', C603-0, (backwater, ''), split=0.208)
    V605 = bst.MixTank('V605', S1-1)
    P606 = bst.Pump('P606', V605-0)
    Ev607 = bst.MultiEffectEvaporator('Ev607',
        ins=P606-0,
        P=(101325, 69682, 47057, 30953, 19781),
        V=0.90
    ) 
    C603_2 = bst.LiquidsSplitCentrifuge('C603_2', Ev607-0, (crude_oil, ''), split={'Oil':0.99})
    
    MX5 = bst.Mixer('MX5', (Ev607-1, P410-0, u.P508-0))
    MX6 = bst.Mixer('MX6', (C603_2-1, MH604-0))
    D610 = bst.DrumDryer('D610', (MX6-0, 'dryer_air', 'natural_gas'), moisture_content=0.10, split=dict(Ethanol=1.0))
    X611 = bst.ThermalOxidizer('X611', (D610-1, 'oxidizer_air', ''))
    MH612 = units.DDGSHandling('MH612', D610-0, DDGS)
    T608 = bst.facilities.ProcessWaterCenter(
        'T608', 
        (MX5-0, 'makeup_water'), 
        ('process_water', 'wastewater'),
        None,
        (),
        (recycled_process_water,)
    )
    other_facilities = units.PlantAir_CIP_WasteWater_Facilities('other_facilities', feedstock)
    # HXN = bst.HeatExchangerNetwork('HXN', 
    #     units=lambda: [u.Ev607.evaporators[0], u.T503_T507.condenser]
    # )
    
    other_unit_parameters = settings.other_unit_parameters
    for ID, params in other_unit_parameters.items():
        attr, val = params.items()
        setattr(f.unit.search(ID), attr, val)
        
    globals().update(f.unit.data)
        
    return f.create_system(system_ID)
    
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
                        