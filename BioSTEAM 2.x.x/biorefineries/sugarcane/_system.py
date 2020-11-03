#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 11:05:24 2017

@author: Yoel
"""
import numpy as np
import biosteam as bst
from biosteam import units
from biosteam import main_flowsheet as F
from ._process_settings import price

__all__ = (
    'create_ethanol_separation_units',
    'create_ethanol_production_system',
    'create_system',)


# %% Ethanol separation

def create_ethanol_separation_units(degassed_beer=None,
                                    ethanol_product=None,
                                    denaturant=None,
                                    IDs={}):
    if not degassed_beer: # Feedstock
        degassed_beer = bst.Stream(
            ID='degassed_beer', 
            phase='l', T=348.34, P=101325, 
            Water=9.65e+04, Ethanol=2.255e+04, Glucose=4916, 
            H3PO4=83.33, Yeast=103, units='kg/hr'
        )
    if not ethanol_product:
        ethanol_product = bst.Stream('ethanol')
    if not denaturant:
        denaturant = bst.Stream('denaturant')
    
    # Beer column
    x_bot = 3.910570816782338e-06
    y_top = 0.2811210085806504
    D302 = units.BinaryDistillation(IDs.get('Beer column', 'D302'), P=101325,
                                y_top=y_top, x_bot=x_bot, k=1.25, Rmin=0.02,
                                LHK=('Ethanol', 'Water'))
    D302.tray_material = 'Stainless steel 304'
    D302.vessel_material = 'Stainless steel 304'
    D302.boiler.U = 1.85
    P302 = units.Pump(IDs.get('Beer column bottoms product pump', 'P302'))
    
    # Mix ethanol Recycle (Set-up)
    M303 = units.Mixer(IDs.get('Recycle mixer', 'M303'))
    
    y_top = 0.8080462715940735
    D303 = units.BinaryDistillation(IDs.get('Distillation', 'D303'), P=101325,
                                y_top=y_top, x_bot=x_bot, k=1.25,
                                LHK=('Ethanol', 'Water'),
                                tray_material='Stainless steel 304',
                                vessel_material='Stainless steel 304',
                                is_divided=True)
    D303.boiler.U = 1.85
    P303 = units.Pump(IDs.get('Distillation bottoms product pump', 'P303'))
    
    # Superheat vapor for mol sieve
    H303 = units.HXutility(IDs.get('Heat exchanger to superheat vapor to molecular sieves', 'H303'),
                           T=115+273.15, V=1)
    
    # Molecular sieve
    U301 = units.MolecularSieve(IDs.get('Molecular sieves', 'U301'),
                                split=(2165.14/13356.04, 1280.06/1383.85),
                                order=('Ethanol', 'Water'))
    
    # Condense ethanol product
    H304 = units.HXutility(IDs.get('Ethanol condenser', 'H304'), 'S149', V=0, T=340.)
    T302 = units.StorageTank(IDs.get('Ethanol day tank', 'T302'), tau=12,
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel')
    P304 = units.Pump(IDs.get('Ethanol day tank pump', 'P304'))
    
    # Storage for gasoline
    T303 = units.StorageTank(IDs.get('Denaturant storage', 'T303'), tau=7*24,
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel')
    P305 = units.Pump(IDs.get('Denaturant pump', 'P305'))
    
    # Denatured ethanol product
    M304 = units.Mixer(IDs.get('Ethanol-denaturant mixer', 'M304'))
    T304 = units.StorageTank(IDs.get('Product tank', 'T304'),
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel',
                             tau=6.5*24, outs=ethanol_product)
    
    ### Ethanol system set-up ###
    
    degassed_beer-D302-1-P302
    (D302-0, U301-0)-M303-0-D303-0-H303-U301
    D303-1-P303
    
    def adjust_denaturant():
        P304._run()
        pure_ethanol = P304.outs[0]
        denaturant.imol['Octane'] = 0.021*pure_ethanol.F_mass/114.232
    
    P304.specification = adjust_denaturant
    U301-1-H304-0-T302-0-P304
    denaturant-T303-P305
    (P305-0, P304-0)-M304-T304

    return {
        'Beer column': D302,
        'Beer column bottoms product pump': P302,
        'Ethanol dehydration system': 
            bst.System('ethanol_dehydration_sys',
                [M303,
                 D303,
                 H303,
                 U301],
                 recycle=U301-0
            ),
        'Recycle mixer': M303,
        'Distillation': D303,
        'Heat exchanger to superheat vapor to molecular sieves': H303,
        'Molecular sieves': U301,
        'Ethanol condenser': H304,
        'Ethanol day tank': T302, 
        'Ethanol day tank pump': P304,
        'Denaturant storage': T303, 
        'Denaturant pump': P305, 
        'Ethanol-denaturant mixer': M304,
        'Product tank': T304,
        'Distillation bottoms product pump': P303,
    }
        

# %% Ethanol production section (fermentation and separations)

def create_ethanol_production_system(ID='ethanol_production_sys',
                                     sugar_solution=None):
    ### Streams ###
    
    # Fresh water
    stripping_water = bst.Stream('stripping_water', Water=5000, units='kg/hr')
    
    # Gasoline
    denaturant = bst.Stream('denaturant', Octane=230.69,
                            units='kg/hr', price=price['Gasoline'])
    
    if not sugar_solution: # Feedstock
        sugar_solution = bst.Stream('sugar_solution',
            Glucose = 3802,
            Sucrose = 4.309e+04,
            Water   = 2.59e+05,
            H3PO4   = 83.33,
            units = 'kg/hr',
            T = 372,
        )
    
    # Yeast
    yeast = bst.Stream('yeast', Water=24700, DryYeast=10300, units='kg/hr')
    
    # Ethanol product
    ethanol = bst.Stream('ethanol', price=price['Ethanol'])
    
    ### Units ###
    
    # Split sugar solution
    S301 = units.Splitter('S301',
                        split=0.265)
    
    # Concentrate sugars
    F301 = units.MultiEffectEvaporator('F301',
                                       P=(101325, 73581, 50892, 32777),
                                       V=0.95) # fraction evaporated
    F301.components['condenser'].U = 1.85
    # Note: value of steam ~ 6.86 for the following 
    # (101325, 73580.467, 50891.17, 32777.406, 19999.925, 11331.5),
    
    # Mix sugar solutions
    M301 = units.Mixer('M301')
    
    # Cool for fermentation
    H301 = units.HXutility('H301', T=295.15)
    
    # Ethanol Production
    R301 = units.Fermentation('R301', outs=('CO2', ''), tau=9, efficiency=0.90, N=4) 
    T301 = units.StorageTank('T301', tau=4, vessel_material='Carbon steel')
    T301.line = 'Beer tank'
    
    D301 = units.VentScrubber('D301', ins=(stripping_water, R301-0), 
                              outs=('vent', ''),
                              gas=('CO2',))
    
    # Separate 99% of yeast
    C301 = units.SolidsCentrifuge('C301', outs=('', 'recycle_yeast'),
                                split=(1, 0.99999, 1, 0.96, 0.01),
                                order=('Ethanol', 'Glucose', 'H3PO4', 
                                       'Water', 'DryYeast'),
                                solids=('DryYeast',))
    
    # Mix in Water
    M302 = units.Mixer('M302')
    P301 = units.Pump('P301')
    
    # Heat up before beer column
    # Exchange heat with stillage
    H302 = units.HXprocess('H302', phase0='l', phase1='l', U=1.28)
    
    dct = create_ethanol_separation_units(
        degassed_beer=H302-0,
        ethanol_product=ethanol,
        denaturant=denaturant)
    D302 = dct['Beer column']
    P302 = dct['Beer column bottoms product pump']
    ethanol_dehydration_sys = dct['Ethanol dehydration system']
    H304 = dct['Ethanol condenser']
    T302 = dct['Ethanol day tank']
    P304 = dct['Ethanol day tank pump']
    T303 = dct['Denaturant storage']
    P305 = dct['Denaturant pump']
    M304 = dct['Ethanol-denaturant mixer']
    T304 = dct['Product tank']
    P303 = dct['Distillation bottoms product pump']
    
    # Waste water
    M305 = units.Mixer('M305', outs='wastewater')
    
    # Yeast mixing
    T305 = units.MixTank('T305')
    T305.tau = 0.1
    yeast-T305
    
    # Multi-effect evaporator pumps
    P306 = units.Pump('P306')
    
    
    ### Ethanol system set-up ###
    
    sugar_solution-S301-1-F301-0-P306
    (S301-0, P306-0)-M301-H301
    (H301-0, yeast-T305-0)-R301-1-T301-0-C301
    (C301-0, D301-1)-M302-P301
    (P303-0, F301-1, (P301-0, P302-0)-H302-1)-M305
    
    ### System ###
    
    return bst.System(ID, 
                [S301, 
                 F301, 
                 P306, 
                 M301, 
                 H301, 
                 T305,
                 R301,
                 T301, 
                 C301, 
                 M302, 
                 P301,
                 bst.System('beer_column_heat_integration',
                     [H302,
                      D302,
                      P302],
                     recycle=P302-0),
                 ethanol_dehydration_sys,
                 H304,
                 T302, 
                 P304,
                 T303, 
                 P305, 
                 M304,
                 T304, 
                 D301,
                 P303, 
                 M305])


# %% Complete system

def create_system(ID='sugarcane_sys'):
    ### Streams ###
    chemicals = bst.settings.get_chemicals()
    
    z_mass_sugarcane = chemicals.kwarray(
        dict(Glucose=0.0120811,
             Lignin=0.0327653,
             Solids=0.015,
             Sucrose=0.136919,
             Ash=0.006,
             Cellulose=0.0611531,
             Hemicellulose=0.036082,
             Water=0.7)
    )
    
    sugarcane = bst.Stream('sugarcane',
                           flow=333334 * z_mass_sugarcane,
                           units='kg/hr',
                           price=price['Sugar cane'])
    
    enzyme = bst.Stream('enzyme',
                        Cellulose=100, Water=900, units='kg/hr',
                        price=price['Protease'])
    
    imbibition_water = bst.Stream('imbibition_water',
                                  Water=87023.35, units='kg/hr',
                                  T = 338.15)
    
    H3PO4 = bst.Stream('H3PO4',
                       H3PO4=74.23, Water=13.10, units='kg/hr',
                       price=price['H3PO4'])  # to T203
    
    lime = bst.Stream('lime',
                      CaO=333.00, Water=2200.00, units='kg/hr',
                      price=price['Lime'])  # to P5
    
    polymer = bst.Stream('polymer',
                         Flocculant=0.83, units='kg/hr',
                         price=price['Polymer'])  # to T205
    
    rvf_wash_water = bst.Stream('rvf_wash_water',
                                Water=16770, units='kg/hr',
                                T=363.15)  # to C202
    
    ### Unit operations ###
    
    # Feed the shredder
    U101 = units.ConveyingBelt('U101', ins=sugarcane)
    
    # Separate metals
    U102 = units.MagneticSeparator('U102', ins=U101-0)
    
    # Shredded cane
    U103 = units.Shredder('U103', ins=U102-0)
    
    # Hydrolyze starch
    T201 = units.EnzymeTreatment('T201', T=323.15)  # T=50
    
    # Finely crush lipid cane
    U201 = units.CrushingMill('U201',
                              split=dict(Ash=0.92,
                                         Cellulose=0.92,
                                         Glucose=0.04,
                                         Hemicellulose=0.92,
                                         Lignin=0.92,
                                         Sucrose=0.04,
                                         Solids=1),
                              moisture_content=0.5)
    
    # Convey out bagasse
    U202 = units.ConveyingBelt('U202', ins=U201.outs[0], outs='Bagasse')
    
    # Mix in water
    M201 = units.Mixer('M201')
    # crushing_mill_recycle_sys = bst.System('crushing_mill_recycle_sys',
    #                                path=(U201, S201, M201),
    #                                recycle=M201-0)
    
    # Screen out fibers
    S201 = units.VibratingScreen('S201',
                                 split=dict(Ash=0.35,
                                            Cellulose=0.35,
                                            Glucose=0.88,
                                            Hemicellulose=0.35,
                                            Lignin=0.35,
                                            Solids=0,
                                            Sucrose=0.88,
                                            Water=0.88))
    
    # Store juice before treatment
    T202 = units.StorageTank('T202', tau=4, vessel_material='Carbon steel')
    
    # Heat up before adding acid
    H201 = units.HXutility('H201', T=343.15)
    
    # Mix in acid
    T203 = units.MixTank('T203')
    
    # Pump acid solution
    P201 = units.Pump('P201')
    
    # Mix lime solution
    T204 = units.MixTank('T204', tau=0.10)
    P202 = units.Pump('P202')
    
    # Blend acid lipid solution with lime
    T205 = units.MixTank('T205', tau=0.10)
    
    # Mix recycle
    M202 = units.Mixer('M202')
    
    # Heat before adding flocculant
    H202 = units.HXutility('H202', T=372.15)
    
    # Mix in flocculant
    T206 = units.MixTank('T206')
    T206.tau = 0.10
    
    # Separate residual solids
    C201 = units.Clarifier('C201',
                           split=dict(Ash=0,
                                      CaO=0,
                                      Cellulose=0,
                                      Flocculant=0.522,
                                      Glucose=0.522,
                                      Hemicellulose=0,
                                      Lignin=0,
                                      H3PO4=0.522,
                                      Sucrose=0.522,
                                      Water=0.522))
    
    # Remove solids as filter cake
    C202 = units.RVF('C202', 
                     outs=('filter_cake', ''),
                     moisture_content=0.80,
                     split=dict(Ash=0.85,
                                CaO=0.85,
                                Cellulose=0.85,
                                Glucose=0.01,
                                Hemicellulose=0.85,
                                Lignin=0.85,
                                Sucrose=0.01))
    P203 = units.Pump('P203')
    
    
    # Screen out small fibers from sugar stream
    S202 = units.VibratingScreen('S202', outs=('', 'fiber_fines'),
                                 split=dict(Ash=1.0,
                                            CaO=1.0,
                                            Cellulose=1.0,
                                            Flocculant=0.0,
                                            Glucose=0.998,
                                            Hemicellulose=1.0,
                                            Lignin=1.0,
                                            H3PO4=1.0,
                                            Sucrose=0.998,
                                            Water=0.998))
    S202.mesh_opening = 2
    
    ### Process specifications ###
    
    # Specifications dependent on lipid cane flow rate
    def correct_flows():
        U103._run()
        F_mass = sugarcane.F_mass
        # correct enzyme, lime, phosphoric acid, and imbibition water
        enzyme.imass['Cellulose', 'Water'] = 0.003 * F_mass * np.array([0.1, 0.9])
        lime.imass['CaO', 'Water'] = 0.001 * F_mass * np.array([0.046, 0.954])
        H3PO4.imass['H3PO4', 'Water'] = 0.00025 * F_mass
        imbibition_water.imass['Water'] = 0.25* F_mass
    
    U103.specification = correct_flows
    
    # Specifications within a system
    def correct_wash_water():
        P202._run()
        solids = P202.outs[0].imol['Ash', 'CaO', 'Cellulose',
                                   'Hemicellulose', 'Lignin'].sum()
        rvf_wash_water.imol['Water'] = 0.0574 * solids
    
    P202.specification = correct_wash_water
    
    ### System set-up ###
    
    (U103-0, enzyme)-T201
    (T201-0, M201-0)-U201-1-S201-0-T202
    (S201-1, imbibition_water)-M201
    
    T202-0-H201
    (H201-0, H3PO4)-T203-P201
    (P201-0, lime-T204-0)-T205-P202
    (P202-0, P203-0)-M202-H202
    (H202-0, polymer)-T206-C201
    (C201-1, rvf_wash_water)-C202-1-P203
    C201-0-S202
    
    ### Ethanol section ###
    
    ethanol_production_sys = create_ethanol_production_system(sugar_solution=S202-0)
    
    ### Facilities ###    
    
    s = F.stream
    BT = units.BoilerTurbogenerator('BT',
                                    (U202-0, '', 'boiler_makeup_water', 'natural_gas', '', ''),
                                    boiler_efficiency=0.80,
                                    turbogenerator_efficiency=0.85)
    
    CT = units.CoolingTower('CT')
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    process_water_streams = (s.imbibition_water,
                             rvf_wash_water,
                             s.stripping_water,
                             *makeup_water_streams)
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    
    CWP = units.ChilledWaterPackage('CWP')
    PWC = units.ProcessWaterCenter('PWC',
                                   (bst.Stream(), makeup_water),
                                   (),
                                   None,
                                   makeup_water_streams,
                                   process_water_streams)
    
    ### System ###
    
    return bst.System(ID,
        [U101,
         U102,
         U103,
         T201,
         bst.System("juice_extraction_sys",
            [U201,
             S201,
             M201],
            recycle=M201-0),
         T202,
         H201,
         T203,
         P201,
         T204,
         T205,
         P202,
         bst.System('juice_separation_sys',
            [M202,
             H202,
             T206,
             C201,
             C202,
             P203],
            recycle=P203-0),
         S202,
         ethanol_production_sys,
         U202],
        facilities=(CWP, BT, CT, PWC),
    )
