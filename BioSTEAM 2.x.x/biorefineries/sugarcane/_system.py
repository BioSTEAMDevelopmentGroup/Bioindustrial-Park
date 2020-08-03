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
from ._ethanol_production_system import create_ethanol_production_system

__all__ = ('create_system',)

# %% Helpful functions

def mass2molar_ethanol_fraction(x):
    """Return ethanol mol fraction in a ethanol water mixture"""
    return x/46.06844 / (x/46.06844 + (1-x)/18.01528)

# %% Pretreatment section

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
    
    (U103, enzyme)-T201
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
    
    return bst.main_flowsheet.create_system(ID)
