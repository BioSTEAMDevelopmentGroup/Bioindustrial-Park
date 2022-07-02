# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
The complete lipid-cane biorefinery system is created here.

"""
import numpy as np
import biosteam as bst
from biosteam import main_flowsheet as f
from biosteam import units, SystemFactory
from ._process_settings import price
from .units import BlendingTankWithSkimming, GlycerolysisReactor
from ..sugarcane import (
    create_feedstock_handling_system,
    create_sucrose_to_ethanol_system, 
    create_juicing_system_up_to_clarification,
    create_juicing_system_with_fiber_screener,
)

__all__ = (
    'create_feedstock_handling_system',
    'create_juicing_system',
    'create_lipid_wash_system',
    'create_juicing_and_lipid_extraction_system',
    'create_lipid_pretreatment_system',
    'create_transesterification_and_biodiesel_separation_system',
    'create_lipidcane_to_biodiesel_and_conventional_ethanol_system',
)


# %% Pretreatment section

@SystemFactory(
    ID='juicing_and_lipid_extraction_sys',
    ins=[dict(ID='lipidcane',
              Ash=2000.042,
              Cellulose=26986.69,
              Glucose=2007.067,
              Hemicellulose=15922.734,
              Lignin=14459.241,
              TAG=10035.334,
              Solids=5017.667,
              Sucrose=22746.761,
              Water=234157.798,
              units='kg/hr',
              price=price['Lipid cane']),
         *create_juicing_system_up_to_clarification.ins[1:]],
    outs=[dict(ID='screened_juice'),
          dict(ID='lipid'),
          dict(ID='bagasse'),
          dict(ID='fiber_fines'),
          dict(ID='spent_oil_wash_water')]
)
def create_juicing_and_lipid_extraction_system(ins, outs, pellet_bagasse=None):
    lipidcane, H3PO4, lime, polymer = ins
    screened_juice, lipid, bagasse, fiber_fines, spent_oil_wash_water = outs
    
    oil_wash_water = bst.Stream('oil_wash_water',
                                Water=1350,
                                units='kg/hr',
                                T=358.15)  # to T207
    juicing_sys = create_juicing_system_up_to_clarification(
        pellet_bagasse=pellet_bagasse,
        ins=[lipidcane, H3PO4, lime, polymer], 
        outs=['', bagasse],
        mockup=True,
    )
    clarified_juice, bagasse = juicing_sys.outs
    u = f.unit
    u.U201.isplit['Lipid'] = 0.10 # Crushing mill
    u.S201.isplit['Lipid'] = 1.00 # Fiber screener 
    u.C201.isplit['Lipid'] = 0.98 # Clarifier
    
    # Separate oil and sugar
    T207 = units.MixTank('T207')
    T207_2 = units.Splitter('T207_2',
                            split=dict(Lipid=1,
                                       Water=1e-4))
    
    # Screen out small fibers from sugar stream
    S202 = units.VibratingScreen('S202', outs=(screened_juice, fiber_fines),
                                 split=dict(Ash=1.0,
                                            CaO=1.0,
                                            Cellulose=1.0,
                                            Flocculant=0.0,
                                            Glucose=0.998,
                                            Hemicellulose=1.0,
                                            Lignin=1.0,
                                            Lipid=1.0,
                                            H3PO4=1.0,
                                            Sucrose=0.998,
                                            Water=0.998))
    S202.mesh_opening = 2
    
    lipid_wash_sys = create_lipid_wash_system(
        ins=T207_2-0,
        outs=(lipid, spent_oil_wash_water),
        mockup=True
    )
    
    ### System set-up ###
    clarified_juice-T207-T207_2
    T207-T207_2-1-S202
 
@SystemFactory(
    ID='lipid_wash_sys',
    ins=[dict(ID='lipid'),
         dict(ID='lipid_wash_water')],
    outs=[dict(ID='washed_lipid'),
          dict(ID='spent_wash_water')]
)
def create_lipid_wash_system(ins, outs):  
    lipid, lipid_wash_water = ins
    washed_lipid, spent_wash_water = outs
    
    # Cool the oil
    H203 = units.HXutility('H203', lipid, T=343.15, V=0, cool_only=True)
    
    recycle = bst.Stream()
    
    # Add water to wash lipid
    T208 = units.MixTank('T208', (H203-0, lipid_wash_water, recycle))
    T208.tau = 0.10
    
    # Centrifuge out water
    C203 = units.LiquidsSplitCentrifuge('C203',
                                        T208-0,
                                        outs=('', spent_wash_water),
                                        split=dict(Lipid=0.99,
                                                   Water=0.01))
    
    # Vacume out water
    F201 = units.SplitFlash('F201', T=357.15, P=10000.,
                            ins=C203-0,
                            outs=('', washed_lipid),
                            split=dict(Lipid=0.0001,
                                       Water=0.999))
    
    H204 = units.HXutility('H204', ins=F201-0, T=320, V=0)
    P204 = units.Pump('P204', H204-0, recycle, P=101325)
    
    # Specifications within a system
    @T208.add_specification(run=True)
    def correct_lipid_wash_water():
        ins = T208.ins
        lipid, lipid_wash_water, recycle, *others = ins
        lipid_wash_water.imol['Water'] = 0.185 * sum([i.imass['Lipid'] for i in ins]) - recycle.imol['Water']

@SystemFactory(
    ID='juicing_sys',
    ins=create_juicing_and_lipid_extraction_system.ins,
    outs=[dict(ID='screened_juice'),
          dict(ID='bagasse'),
          dict(ID='fiber_fines')]
)
def create_juicing_system(ins, outs, pellet_bagasse=None):
    lipidcane, H3PO4, lime, polymer = ins
    screened_juice, bagasse, fiber_fines = outs
    
    oil_wash_water = bst.Stream('oil_wash_water',
                                Water=1350,
                                units='kg/hr',
                                T=358.15)  # to T207
    juicing_sys = create_juicing_system_with_fiber_screener(
        pellet_bagasse=pellet_bagasse,
        ins=[lipidcane, H3PO4, lime, polymer], 
        outs=[screened_juice, bagasse, fiber_fines],
        mockup=True,
    )
    u = f.unit
    u.U201.isplit['Lipid'] = 0.90 # Crushing mill
    u.S201.isplit['Lipid'] = 1.0 # Fiber screener #1
    u.C201.isplit['Lipid'] = 0.99 # Clarifier
    u.S202.isplit['Lipid'] = 0.999  # Fiber screener #2
    
@SystemFactory(
    ID="lipid_pretreatment_sys",
    ins=[dict(ID='crude_vegetable_oil',
              Water=0.0184,
              TAG=11.1,
              PL=0.1),
         dict(ID='acetone',
              price=0.80),
         dict(ID='pure_glycerine',
              price=0.65)],
    outs=[dict(ID='degummed_oil'),
          dict(ID='polar_lipids'),
          dict(ID='wastewater')]
)
def create_lipid_pretreatment_system(ins, outs):
    crude_vegetable_oil, acetone, pure_glycerine = ins
    degummed_oil, polar_lipids, wastewater = outs
    
    # Vacume out water
    F3 = units.SplitFlash('F3', T=357.15, P=10000.,
                          ins=crude_vegetable_oil,
                          split=dict(Lipid=0.0001,
                                     Water=1.0),)
    P10 = units.Pump('P10', F3-1, P=101325)
    H5 = units.HXutility('H5', ins=F3-0, T=320, V=0)
    P9 = units.Pump('P9', H5-0, wastewater, P=101325)
    T1 = bst.StorageTank('T1', acetone, tau=7*24)
    P1 = bst.Pump('P1', T1-0, P=101325.)
    T2 = bst.StorageTank('T2', pure_glycerine, tau=7*24)
    P2 = bst.Pump('P2', T2-0, P=101325.)
    N2 = bst.Stream('N2', phase='g')
    acetone_recycle = bst.Stream()
    T3 = BlendingTankWithSkimming('T3', 
        [P10-0, acetone_recycle], 
    )
    @T3.add_specification(run=True)
    def adjust_acetone_flow_rate():
        total_acetone_required = T3.ins[0].F_vol
        fresh_acetone = T1.ins[0]
        recycle_acetone = sum([i.outs[0].ivol['Acetone'] for i in (P5, P6)])
        acetone_required = total_acetone_required - recycle_acetone
        if acetone_required < 0.:
            fresh_acetone.ivol['Acetone'] = 0.
            P5.outs[0].ivol['Acetone'] = total_acetone_required
        else:
            fresh_acetone.ivol['Acetone'] = acetone_required
        for i in T1.path_until(M1): i._run()
        M1._run()
    
    P3 = bst.Pump('P3', T3-1, P=101325)
    F1 = bst.Flash('F1', P3-0, V=1., P=101325)
    H1 = bst.HXutility('H1', F1-0, V=0, rigorous=True)
    P5 = bst.Pump('P5', H1-0, P=101325)
    F2 = bst.Flash('F2', T3-0, ['', polar_lipids], V=1.0, P=101325)
    H2 = bst.HXutility('H2', F2-0, V=0, rigorous=True)
    P6 = bst.Pump('P6', H2-0, P=101325)
    M1 = bst.Mixer('M1', [P1-0, P5-0, P6-0], acetone_recycle)
    # bst.mark_disjunction(acetone_recycle)
    P4 = bst.Pump('P4', F1-1, P=101325)
    hx_stream = bst.Stream()
    H3 = bst.HXprocess('H3', [P4-0, hx_stream])
    glycerol_recycle = bst.Stream()
    M2 = bst.Mixer('M2', [P2-0, glycerol_recycle])
    R1 = GlycerolysisReactor('R1', [H3-0, M2-0, N2])
    P7 = bst.Pump('P7', R1-1, hx_stream, P=101325.)
    
    @R1.add_specification(run=True)
    def adjust_feed_flow_rates():
        lipid = R1.ins[0]
        required_glycerol = 1.5 * (
            + lipid.imol['FFA']
            + lipid.imol['TAG']
        ) - M2.ins[1].imol['Glycerol']
        if required_glycerol < 0.:
            T2.ins[0].imol['Glycerol'] = 0.
            M2.ins[1].imol['Glycerol'] += required_glycerol
        else:
            T2.ins[0].imol['Glycerol'] = required_glycerol
        R1.ins[2].ivol['N2'] = lipid.F_vol
        for i in T2.path_until(R1): i._run()
        
    H4 = bst.HXutility('H4', H3-1, T=333.15, V=0)
    C1 = bst.LiquidsSplitCentrifuge('C1', H4-0, ['', glycerol_recycle],
                                    split=dict(Lipid=1.))
    P8 = bst.Pump('P8', C1-0, degummed_oil)

@SystemFactory(
    ID="transesterification_and_biodiesel_separation_sys",
    ins=[dict(ID='vegetable_oil',
              Water=0.0184,
              TAG=11.1)],
    outs=[dict(ID='biodiesel'),
          dict(ID='crude_glycerol'),
          dict(ID='wastewater',
               price=price['Waste'])]
)
def create_transesterification_and_biodiesel_separation_system(ins, outs):
    ### Streams ###
    
    oil, = ins
    biodiesel, crude_glycerol, wastewater = outs
    
    # Fresh methanol
    methanol = bst.Stream('methanol',
                          Methanol=1,
                          price=price['Methanol'])
    
    # Catalyst
    catalyst = bst.Stream('catalyst', 
                          NaOCH3=0.25,
                          Methanol=0.75, 
                          units='kg/hr',
                          price=price['NaOCH3'])
    
    # Water to remove glycerol
    biodiesel_wash_water = bst.Stream('biodiesel_wash_water', Water=13.6, T=273.15+60, 
                                      price=price['Water'])
    
    HCl = bst.Stream('HCl', HCl=0.21, Water=0.79,
                      price=0.35 * price['HCl']) # 35% HCl by mass
    
    # Acid to neutralize catalyst after second centrifuge
    HCl1 = bst.Stream('HCl1', HCl=0.21, Water=0.79)
    
    # Acid to remove soaps after first centrifuge
    HCl2 = bst.Stream('HCl2', HCl=0.21, Water=0.79)
    
    S402 = bst.FakeSplitter('S402', ins=HCl, outs=(HCl1, HCl2))
    
    # Base to neutralize acid before distillation
    NaOH = bst.Stream('NaOH', NaOH=1, price=price['NaOH'])
    
    ### Units ###
    
    ### Biodiesel Transesterification Section ###
    
    # Aparently reactors are adiabatic, meoh coming in at 40C, lipid at 60C
    
    # From USDA Biodiesel model
    x_cat = 1.05e-05 # Catalyst molar fraction in methanol feed
    
    # Mix Recycle Methanol and Fresh Methanol
    T401 = units.StorageTank('T401')
    P401 = units.Pump('P401')
    
    # Storage Tank for Catalyst
    T402 = units.StorageTank('T402')
    P402 = units.Pump('P402')
    
    # Tank for oil
    T403 = units.StorageTank('T403')
    T403.tau = 4
    P403 = units.Pump('P403')
    
    # Mix Methanol and Catalyst stream
    T404 = units.MixTank('T404')
    P404 = units.Pump('P404')
    
    # Split Methanol/Catalyst to reactors (this is done through a process specification, so use a fake splitter)
    S401 = bst.FakeSplitter('S401')
    
    # First Reactor
    R401 = units.Transesterification('R401', efficiency=0.90, excess_methanol=1.,
                                     T=333.15, x_catalyst=x_cat)
    
    # Centrifuge to remove glycerol
    C401 = units.LiquidsSplitCentrifuge('C401',
                                 split=dict(Lipid=0.99,  
                                            Methanol=0.40,  
                                            Glycerol=0.06, 
                                            Biodiesel=0.999, 
                                            Water=0.40, 
                                            NaOH=0,
                                            HCl=0,
                                            NaOCH3=0.40)) 
    
    P405 = units.Pump('P405')
    
    # Second Reactor
    R402 = units.Transesterification('R402', efficiency=0.90, excess_methanol=1., 
                                     T=333.15, x_catalyst=x_cat)
    
    def adjust_feed_to_reactors():
        R402._run()
        for i in S401.outs:
            if isinstance(i, bst.Junction): i._run()
        S401.ins[0].mix_from(S401.outs)
    R402.specification = adjust_feed_to_reactors
    
    # Centrifuge to remove glycerol
    C402 = units.LiquidsSplitCentrifuge('C402',
                             split=dict(Lipid=0.90, 
                                        Methanol=0.10, 
                                        Glycerol=0.05, 
                                        Biodiesel=0.999, 
                                        Water=0.10, 
                                        NaOH=0,  
                                        HCl=0, 
                                        NaOCH3=0.10)) 
    
    # Acids and bases per catalyst by mol
    k1 = 0.323/1.5; k2 = 1.060/1.5; k3 = 0.04505/1.5
    def adjust_acid_and_base():
        T404._run()
        # Adjust according to USDA biodiesel model
        f = 0.79 / 0.21
        catalyst_mol = T404.outs[0].imol['NaOCH3']
        NaOH.imol['NaOH'] = k1 * catalyst_mol
        HCl1.imol['HCl'] = mol1 = k2 * catalyst_mol
        HCl2.imol['HCl'] = mol2 = k3 * catalyst_mol
        HCl.imol['HCl'] = mol12 = mol1 + mol2
        HCl.imol['Water'] = f * mol12
        HCl1.imol['Water'] = f * mol1
        HCl2.imol['Water'] = f * mol2
    
    T404.specification = adjust_acid_and_base
    
    ### Biodiesel Purification Section ###
    
    # Water wash
    T405 = units.MixTank('T405')
    P406 = units.Pump('P406')
    
    # Centrifuge out water
    C403 = units.LiquidsRatioCentrifuge('C403',
                             K_chemicals=('Methanol', 'Glycerol'),
                             Ks=np.array([0.382, 0.183]),
                             top_solvents=('Biodiesel',),
                             top_split=(0.999,),
                             bot_solvents=('Water', 'TAG', 'NaOH', 'HCl'),
                             bot_split=(0.999, 1, 1, 1))
    
    # Vacuum dry biodiesel
    # Consider doing this by split, and keeping this unit constant
    # 290 Pa, 324 K according to USDA Biodiesel Model
    F401 = units.SplitFlash('F401',
                    order=('Water', 'Methanol', 'Biodiesel'),
                    split=(0.9999, 0.9999, 0.00001),
                    P=10000., T=331.5, Q=0)
    F401.line = 'Vacuum dryer'
    F401.material = 'Stainless steel 304'
    P407 = units.Pump('P407', P=101325.)
    
    ### Glycerol Purification Section ###
    
    # Condense vacuumed methanol to recycle
    H401 = units.HXutility('H401', V=0, T=295)
    P408 = units.Pump('P408', P=101325.)
    
    # Mix recycled streams and HCl
    T406 = units.MixTank('T406')
    P409 = units.Pump('P409', P=101325.)
    
    # Centrifuge out waste fat
    # assume all the lipid, free_lipid and biodiesel is washed out
    C404 = units.LiquidsSplitCentrifuge('C404', outs=('', wastewater),
                                         order=('Methanol', 'Glycerol', 'Water'),
                                         split=(0.999, 0.999, 0.99))
    
    # Add and mix NaOH
    T407 = units.MixTank('T407')
    P410 = units.Pump('P410')
    
    # Methanol/Water distillation column
    D401 = units.BinaryDistillation('D401',
                      LHK=('Methanol', 'Water'), P=101325,
                      Lr=0.999, Hr=0.99999, k=2.5,
                      is_divided=True,
                      partial_condenser=False,
                      vessel_material='Stainless steel 304',
                      tray_material='Stainless steel 304')
    
    # Save column utilities
    H402 = units.HXprocess('H402')
    
    # Glycerol/Water flash (not a distillation column)
    chemicals = H402.chemicals
    w = 0.20/chemicals.Water.MW
    g = 0.80/chemicals.Glycerol.MW
    x_water = w/(w+g)
    # ideal_thermo = D401.thermo.ideal()
    D402 = units.BinaryDistillation('D402',
                        LHK=('Water', 'Glycerol'),
                        k=1.25,
                        P=101325,
                        y_top=0.999999,
                        x_bot=x_water,
                        tray_material='Stainless steel 304',
                        vessel_material='Stainless steel 304',
                        partial_condenser=False)
    P413 = units.Pump('P413', P=101325.,)
    
    def startup_water():
        imol = D402.ins[0].imol
        water, glycerol = imol['Water', 'Glycerol']
        minimum_water = 5 * x_water * glycerol
        if water < minimum_water:
            imol['Water'] = minimum_water
        D402._run()
        # Remove accumulation
        D402.outs[0].imol['Water'] = 800.*C402.outs[0].imol['Glycerol']
        D402.outs[0].T = D402.outs[0].bubble_point_at_P().T
            
    D402.specification = startup_water
    
    # Condense recycle methanol
    H403 = units.HXutility('H403', V=0, T=315)
    P411 = units.Pump('P411')
    
    # Condense recycle water
    H404 = units.HXutility('H404', V=0, T=315)
    P412 = units.Pump('P412')
    
    # Storage tank for glycerol
    T408 = units.StorageTank('T408', outs=crude_glycerol)
    
    # Storage tank for biodiesel
    T409 = units.StorageTank('T409', outs=biodiesel)
    F401-1-P407-0-T409
    
    ### Biodiesel system set-up ###
    
    # Biodiesel Transesterification Section
    oil-T403-P403
    (P403-0, S401-0)-R401-0-C401
    (C401-0, S401-1)-R402-0-C402-1
    
    # Specs for product https://www.afdc.energy.gov/fuels/biodiesel_specifications.html
    # minimum spec requirements:
    #  0.050 wt % water (lower to 0.045)
    #  0.2 wt % meoh (lower to 0.15)
    #  0.02 wt % glycerol (lower to 0.015)
    #  0.4 wt % free lipids (lower to 0.35)
    
    # Find Water Flow
    def adjust_biodiesel_wash_water():
        total_glycerol =  (C401.outs[1].imol['Glycerol'] + R402.outs[0].imol['Glycerol'])
        wash_water = (x_water / (1 - x_water) * total_glycerol
                      - HCl1.imol['Water']
                      - NaOH.imol['Water']
                      - oil.imol['Water']
                      - HCl2.imol['Water'])
        biodiesel_wash_water.imol['Water'] = wash_water if wash_water > 0 else 0.
        T405._run()
    
    T405.specification = adjust_biodiesel_wash_water
    D402-0-H404-0-P412
    
    # Biodiesel wash
    (C402-0, P412-0, biodiesel_wash_water, HCl1)-T405-P406-C403
    
    # Glycerol recycle and purification section
    C403-0-F401
    F401-0-H401-P408
    C401-1-P405
    (P405-0, C402-1, C403-1, P408-0, HCl2)-T406-P409-C404
    (C404-0, NaOH)-T407-P410
    P410-0-H402-0-D401-1-D402-1-P413-0-1-H402-1-T408
                 
    
    # Mass Balance for Methanol, Recycle Methanol, and Catalyst stream
    B401 = bst.MassBalance('B401',
                           description='Adjust catalyst and methanol feed to reactors',
                           variable_inlets=[catalyst, methanol],
                           constant_inlets=[D401-0],
                           constant_outlets=[1**R401, 1**R402],
                           chemical_IDs=('Methanol', 'NaOCH3'))
    
    # Find Fresh Methanol Flow
    D401-0-B401-H403-P411 # Recycle methanol
    methanol-T401-P401  # Mix fresh and recycled methanol
    catalyst-T402-P402  # Fresh catalyst
    (P411-0, P401-0, P402-0)-T404-P404-S401  # Mix Catalyst with Methanol
    
    # Initial guess
    D402.outs[0].imol['Methanol', 'Glycerol', 'Water'] = [0.00127, 3.59e-05, 0.0244]
    
    bst.mark_disjunction(P411-0)
    

@SystemFactory(
    ID='lipidcane_sys',
    ins=[*create_juicing_and_lipid_extraction_system.ins,
         create_sucrose_to_ethanol_system.ins[1]],
    outs=[dict(ID='ethanol', price=price['Ethanol']),
          dict(ID='biodiesel', price=price['Biodiesel']),
          dict(ID='crude_glycerol', price=price['Crude glycerol']),
          dict(ID='wastewater'),
          dict(ID='emissions'),
          dict(ID='ash_disposal')]
)
def create_lipidcane_to_biodiesel_and_conventional_ethanol_system(ins, outs):
    
    lipidcane, H3PO4, lime, polymer, denaturant = ins
    ethanol, biodiesel, crude_glycerol, wastewater, emissions, ash_disposal = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=lipidcane,
        mockup=True,
    )
    
    enzyme = bst.Stream(
        ID='enzyme',
        Cellulose=100,
        Water=900,
        units='kg/hr',
        price=0.5
    )
    
    # Hydrolyze lipid bodies
    T201 = units.EnzymeTreatment('T201', (feedstock_handling_sys-0, enzyme), T=323.15)  # T=50
    
    @T201.add_specification(run=True)
    def update_enzyme():
        sugarcane = T201.ins[0]
        F_mass = sugarcane.F_mass
        enzyme.imass['Cellulose', 'Water'] = 0.003 * F_mass * np.array([0.1, 0.9])
    
    ### Oil and juice separation ###
    
    juicing_and_lipid_extraction_sys = create_juicing_and_lipid_extraction_system(
        ins=[T201-0, H3PO4, lime, polymer],
        mockup=True,
    )
    
    ### Ethanol section ###
    
    ethanol_production_sys = create_sucrose_to_ethanol_system(
        ins=[juicing_and_lipid_extraction_sys-0, denaturant],
        outs=[ethanol],
        mockup=True,
    )
    
    ### Biodiesel section ###
    
    # Fresh degummed oil
    oil = juicing_and_lipid_extraction_sys-1
    create_transesterification_and_biodiesel_separation_system(
        ins=oil, 
        outs=[biodiesel, crude_glycerol],
        mockup=True,
    )

    ### Facilities ###
    
    s = f.stream
    u = f.unit
    
    M305 = bst.Mixer('M305', 
        ins=(juicing_and_lipid_extraction_sys-4, 
             juicing_and_lipid_extraction_sys-3,
             *ethanol_production_sys-[1, 2, 3]),
        outs=wastewater,
    )
    
    # Burn bagasse from conveyor belt
    BT = units.BoilerTurbogenerator('BT',
                                   (juicing_and_lipid_extraction_sys-2, '', 
                                    'boiler_makeup_water', 'natural_gas', '', ''),
                                   (emissions, 'rejected_water_and_blowdown', ash_disposal),
                                   boiler_efficiency=0.80,
                                   turbogenerator_efficiency=0.85)
    
    CT = units.CoolingTower('CT')
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    
    process_water_streams = (s.imbibition_water,
                             s.biodiesel_wash_water,
                             s.oil_wash_water,
                             s.rvf_wash_water,
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
    
    F301 = u.F301
    D303 = u.D303
    HXN = bst.HeatExchangerNetwork('HXN', units=[F301, D303])