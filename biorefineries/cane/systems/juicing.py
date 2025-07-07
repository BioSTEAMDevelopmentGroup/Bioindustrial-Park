#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import biosteam as bst
from biosteam import SystemFactory
from biosteam import main_flowsheet as f
from .bagasse import (
    create_bagasse_drying_system,
    create_bagasse_pelleting_system,
)
from .. import streams as s

__all__ = (
    'create_feedstock_handling_system',
    'create_juicing_system_up_to_clarification',
    'create_juicing_system_without_treatment',
    'create_juicing_system',
)

@SystemFactory(
    ID='feedstock_handling_sys',
    ins=[s.sugarcane],
    outs=[s.shredded_cane]
)
def create_feedstock_handling_system(ins, outs):
    sugarcane, = ins
    shreaded_cane, = outs
    U101 = bst.ConveyingBelt('U101', sugarcane)
    U102 = bst.MagneticSeparator('U102', U101-0)
    U103 = bst.Shredder('U103', U102-0, shreaded_cane)

@SystemFactory(
    ID='juicing_sys',
    ins=[s.sugarcane],
    outs=[s.untreated_juice, s.bagasse]
)
def create_juicing_system_without_treatment(ins, outs, pellet_bagasse=None,
                                            dry_bagasse=None):
    if pellet_bagasse is None: pellet_bagasse = False
    if dry_bagasse is None: dry_bagasse = False
    
    ### Streams ###
    sugarcane, = ins
    untreated_juice, bagasse = outs
    
    imbibition_water = bst.Stream('imbibition_water',
                                  Water=87023.35, units='kg/hr',
                                  T = 350.15)
    
    # Finely crush lipid cane
    U201 = bst.CrushingMill('U201',
                              ins=(sugarcane, ''),
                              split=dict(Ash=0.92,
                                         Cellulose=0.92,
                                         Glucose=0.04,
                                         Hemicellulose=0.92,
                                         Lignin=0.92,
                                         Sucrose=0.04,
                                         Solids=1),
                              moisture_content=0.5)
    
    @U201.add_specification(run=True)
    def update_imbibition_water():
        feed = U201.ins[0]
        moisture = feed.imass['Water']
        dry = feed.F_mass - moisture
        imbibition_water.imass['Water'] = max(2.68333333 * dry - moisture, 0)
    
    U202 = bst.ConveyingBelt('U202', U201-0, [''] if (pellet_bagasse or dry_bagasse) else [bagasse])
    
    if pellet_bagasse:
        create_bagasse_pelleting_system(None, ins=U202-0, outs=bagasse, mockup=True)
    elif dry_bagasse:
        create_bagasse_drying_system(None, ins=U202-0, outs=bagasse, mockup=True)
    
    # Mix in water
    M201 = bst.Mixer('M201', ('', imbibition_water), 1-U201)
    
    # Screen out fibers
    S201 = bst.VibratingScreen('S201', U201-1, ('', 0-M201),
                                 split=dict(Ash=0.35,
                                            Cellulose=0.35,
                                            Glucose=0.88,
                                            Hemicellulose=0.35,
                                            Lignin=0.35,
                                            Solids=0,
                                            Sucrose=0.88,
                                            Water=0.88))
    
    # Store juice before treatment
    T202 = bst.StorageTank('T202', S201-0, untreated_juice,
                             tau=4, vessel_material='Carbon steel')
    
    if 'Lipid' in bst.settings.chemicals:
        u = f.unit
        u.U201.isplit['Lipid'] = 0.90 # Crushing mill
        u.S201.isplit['Lipid'] = 1.0 # Fiber screener #1


@SystemFactory(
    ID='juicing_sys',
    ins=[s.sugarcane, s.H3PO4, s.lime, s.polymer],
    outs=[s.clarified_juice, s.bagasse, s.filter_cake]
)
def create_juicing_system_up_to_clarification(ins, outs, pellet_bagasse=None,
                                              dry_bagasse=None):
    
    ### Streams ###
    sugarcane, H3PO4, lime, polymer = ins
    clarified_juice, bagasse, filter_cake = outs
    
    rvf_wash_water = bst.Stream('rvf_wash_water',
                                Water=16770, units='kg/hr',
                                T=363.15)  # to C202
    
    juicing_sys_without_treatment = create_juicing_system_without_treatment(
        ins=(sugarcane,),
        outs=('', bagasse),
        mockup=True,
        pellet_bagasse=pellet_bagasse,
        dry_bagasse=dry_bagasse,
    )
    U201 = f.unit.U201
    
    # Heat up before adding acid
    H201 = bst.HXutility('H201', juicing_sys_without_treatment-0, T=343.15)
    
    # Mix in acid
    T203 = bst.MixTank('T203', (H201-0, H3PO4))
    
    # Pump acid solution
    P201 = bst.Pump('P201', T203-0)
    
    # Mix lime solution
    T204 = bst.MixTank('T204', lime, tau=0.10)
    
    # Blend acid lipid solution with lime
    T205 = bst.MixTank('T205', (P201-0, T204-0), tau=0.10)
    P202 = bst.Pump('P202', T205-0)
    
    # Mix recycle
    M202 = bst.Mixer('M202', (P202-0, ''))
    
    # Heat before adding flocculant
    H202 = bst.HXutility('H202', M202-0, T=372.15)
    
    # Mix in flocculant
    T206 = bst.MixTank('T206', (H202-0, polymer))
    T206.tau = 0.10
    
    # Separate residual solids
    C201 = bst.Clarifier('C201', T206-0, (clarified_juice, ''),
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
    C202 = bst.RVF('C202', (C201-1, rvf_wash_water), (filter_cake, ''),
                     moisture_content=0.80,
                     split=dict(Ash=0.85,
                                CaO=0.85,
                                Cellulose=0.85,
                                Glucose=0.01,
                                Hemicellulose=0.85,
                                Lignin=0.85,
                                Sucrose=0.01))
    P203 = bst.Pump('P203', C202-1, 1-M202)
    
    ### Process specifications ###
    
    # Specifications dependent on lipid cane flow rate
    @U201.add_specification(run=True)
    def correct_flows():
        feed = U201.ins[0]
        other_liquids = feed.imass['Lipid'] if 'Lipid' in feed.chemicals else 0.
        F_mass = (feed.F_mass - feed.imass['Water'] - other_liquids) / 0.7
        # correct lime, phosphoric acid, and imbibition water
        lime.imass['CaO', 'Water'] = 0.001 * F_mass * np.array([0.046, 0.954])
        H3PO4.imass['H3PO4', 'Water'] = 0.00025 * F_mass
    
    # Specifications within a system
    @P202.add_specification
    def correct_wash_water():
        P202._run()
        solids = P202.outs[0].imol['Ash', 'CaO', 'Cellulose',
                                   'Hemicellulose', 'Lignin'].sum()
        rvf_wash_water.imol['Water'] = 0.0574 * solids
        
    if 'Lipid' in bst.settings.chemicals:
        u = f.unit
        u.C201.isplit['Lipid'] = 0.99 # Clarifier
    

@SystemFactory(
    ID='juicing_sys',
    ins=[s.sugarcane, s.H3PO4, s.lime, s.polymer],
    outs=[s.screened_juice, s.bagasse, s.fiber_fines]
)          
def create_juicing_system(ins, outs, pellet_bagasse=None, dry_bagasse=None):
    screened_juice, bagasse, fiber_fines = outs
    sys = create_juicing_system_up_to_clarification(
        None, ins, ['', bagasse],
        mockup=True,
        pellet_bagasse=pellet_bagasse,
        dry_bagasse=dry_bagasse,
    )
    # Screen out small fibers from sugar stream
    S202 = bst.VibratingScreen('S202', sys-0, (screened_juice, fiber_fines),
                                 split=dict(Ash=0.998,
                                            CaO=0.998,
                                            Cellulose=0.0,
                                            Flocculant=0.0,
                                            Glucose=0.998,
                                            Hemicellulose=0.0,
                                            Lignin=0.0,
                                            H3PO4=1.0,
                                            Sucrose=0.998,
                                            Water=0.998))
    S202.mesh_opening = 2
    if 'Lipid' in bst.settings.chemicals:
        u = f.unit
        u.S202.isplit['Lipid'] = 0.999  # Fiber screener #2
