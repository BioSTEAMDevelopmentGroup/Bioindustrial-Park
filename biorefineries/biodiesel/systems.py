# -*- coding: utf-8 -*-
"""
.. contents:: :local:
    
.. autofunction:: biorefineries.biodiesel.systems.create_lipid_wash_system
.. autofunction:: biorefineries.biodiesel.systems.create_lipid_pretreatment_system
.. autofunction:: biorefineries.biodiesel.systems.create_transesterification_and_biodiesel_separation_system

References
----------
.. [1] Cortés-Peña et al. BioSTEAM: A Fast and Flexible Platform for the Design,
 Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
 ACS Sustainable Chem. Eng. 2020. https://doi.org/10.1021/acssuschemeng.9b07040

.. [2] Cortés-Peña et al. Economic and Environmental Sustainability of Vegetative 
 Oil Extraction Strategies at Integrated Oilcane and Oil-sorghum Biorefineries. 
 ACS Sustainable Chem. Eng. 2022

"""
import biosteam as bst
from . import units
from . import streams as s
import numpy as np
from biosteam import SystemFactory

__all__ = (
    'create_lipid_wash_system',
    'create_lipid_pretreatment_system',
    'create_transesterification_and_biodiesel_separation_system',    
)

@SystemFactory(
    ID='lipid_wash_sys',
    ins=[s.lipid, s.lipid_wash_water],
    outs=[s.washed_lipid, s.spent_wash_water]
)
def create_lipid_wash_system(ins, outs):  
    lipid, lipid_wash_water = ins
    washed_lipid, spent_wash_water = outs
    
    # Cool the oil
    H203 = bst.HXutility('H203', lipid, T=343.15, V=0, cool_only=True)
    
    recycle = bst.Stream()
    
    # Add water to wash lipid
    T208 = bst.MixTank('T208', (H203-0, lipid_wash_water, recycle))
    T208.tau = 0.10
    
    # Centrifuge out water
    C203 = bst.LiquidsSplitCentrifuge('C203',
                                        T208-0,
                                        outs=('', spent_wash_water),
                                        split=dict(Lipid=0.99,
                                                   Water=0.01))
    
    # Vacume out water
    F201 = bst.SplitFlash('F201', T=357.15, P=10000.,
                            ins=C203-0,
                            outs=('', washed_lipid),
                            split=dict(Lipid=0.0001,
                                       Water=0.999))
    
    H204 = bst.HXutility('H204', ins=F201-0, T=320, V=0)
    P204 = bst.Pump('P204', H204-0, recycle, P=101325)
    
    # Specifications within a system
    @T208.add_specification(run=True)
    def correct_lipid_wash_water():
        ins = T208.ins
        lipid, lipid_wash_water, recycle, *others = ins
        lipid_wash_water.imol['Water'] = 0.185 * sum([i.imass['Lipid'] for i in ins]) - recycle.imol['Water']

@SystemFactory(
    ID="lipid_pretreatment_sys",
    ins=[s.crude_vegetable_oil, s.acetone, s.pure_glycerine],
    outs=[s.degummed_oil, s.polar_lipids, s.wastewater]
)
def create_lipid_pretreatment_system(ins, outs):
    crude_vegetable_oil, acetone, pure_glycerine = ins
    degummed_oil, polar_lipids, wastewater = outs
    
    # Vacume out water
    F3 = bst.SplitFlash('F3', T=357.15, P=10000.,
                          ins=crude_vegetable_oil,
                          split=dict(Lipid=0.00001,
                                     Water=1.0),)
    P10 = bst.Pump('P10', F3-1, P=101325)
    H5 = bst.HXutility('H5', ins=F3-0, T=320, V=0)
    P9 = bst.Pump('P9', H5-0, wastewater, P=101325)
    T1 = bst.StorageTank('T1', acetone, tau=7*24)
    P1 = bst.Pump('P1', T1-0, P=101325.)
    T2 = bst.StorageTank('T2', pure_glycerine, tau=7*24)
    P2 = bst.Pump('P2', T2-0, P=101325.)
    N2 = bst.Stream('N2', phase='g')
    acetone_recycle = bst.Stream()
    T3 = units.BlendingTankWithSkimming('T3', 
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
    R1 = units.GlycerolysisReactor('R1', [H3-0, M2-0, N2])
    P7 = bst.Pump('P7', R1-1, hx_stream, P=101325.)
    @R1.add_specification
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
        T2.run_until(P7)
        
    H4 = bst.HXutility('H4', H3-1, T=333.15, V=0)
    C1 = bst.LiquidsSplitCentrifuge('C1', H4-0, ['', glycerol_recycle],
                                    split=dict(Lipid=1.))
    P8 = bst.Pump('P8', C1-0, degummed_oil)

@SystemFactory(
    ID="transesterification_and_biodiesel_separation_sys",
    ins=[s.vegetable_oil],
    outs=[s.biodiesel, s.crude_glycerol, s.wastewater]
)
def create_transesterification_and_biodiesel_separation_system(ins, outs,
        transesterification_reactions=None,
    ):
    ### Streams ###
    
    oil, = ins
    biodiesel, crude_glycerol, wastewater = outs
    
    # Fresh methanol
    methanol = bst.Stream('methanol',
                          Methanol=1,
                          price=0.547)
    
    # Catalyst
    catalyst = bst.Stream('catalyst', 
                          NaOCH3=0.25,
                          Methanol=0.75, 
                          units='kg/hr',
                          price=2.93)
    
    # Water to remove glycerol
    biodiesel_wash_water = bst.Stream('biodiesel_wash_water', Water=13.6, T=273.15+60, 
                                      price=0.000353)
    
    HCl = bst.Stream('HCl', HCl=0.21, Water=0.79,
                      price=0.35 * 0.205) # 35% HCl by mass
    
    # Acid to neutralize catalyst after second centrifuge
    HCl1 = bst.Stream('HCl1', HCl=0.21, Water=0.79)
    
    # Acid to remove soaps after first centrifuge
    HCl2 = bst.Stream('HCl2', HCl=0.21, Water=0.79)
    
    S402 = bst.FakeSplitter('S402', ins=HCl, outs=(HCl1, HCl2))
    
    # Base to neutralize acid before distillation
    NaOH = bst.Stream('NaOH', NaOH=1, price=0.41)
    
    ### Units ###
    
    ### Biodiesel Transesterification Section ###
    
    # Aparently reactors are adiabatic, meoh coming in at 40C, lipid at 60C
    
    # From USDA Biodiesel model
    x_cat = 1.05e-05 # Catalyst molar fraction in methanol feed
    
    # Mix Recycle Methanol and Fresh Methanol
    T401 = bst.StorageTank('T401')
    P401 = bst.Pump('P401')
    
    # Storage Tank for Catalyst
    T402 = bst.StorageTank('T402')
    P402 = bst.Pump('P402')
    
    # Tank for oil
    T403 = bst.StorageTank('T403')
    T403.tau = 4
    P403 = bst.Pump('P403')
    
    # Mix Methanol and Catalyst stream
    T404 = bst.MixTank('T404')
    P404 = bst.Pump('P404')
    
    # Split Methanol/Catalyst to reactors (this is done through a process specification, so use a fake splitter)
    S401 = bst.FakeSplitter('S401')
    
    # First Reactor
    R401 = units.Transesterification('R401', 
        efficiency=0.90, excess_methanol=1., T=333.15, x_catalyst=x_cat,
        transesterification=transesterification_reactions,
    )
    
    # Centrifuge to remove glycerol
    C401 = bst.LiquidsSplitCentrifuge('C401',
                                 split=dict(Lipid=0.99,  
                                            Methanol=0.40,  
                                            Glycerol=0.06, 
                                            Biodiesel=0.999, 
                                            Water=0.40, 
                                            NaOH=0,
                                            HCl=0,
                                            NaOCH3=0.40)) 
    
    P405 = bst.Pump('P405')
    
    # Second Reactor
    R402 = units.Transesterification('R402',
        efficiency=0.90, excess_methanol=1., T=333.15, x_catalyst=x_cat,
        transesterification=transesterification_reactions,
    )
    @R402.add_specification
    def adjust_feed_to_reactors():
        R402._run()
        S401.ins[0].mix_from(S401.outs)
    
    # Centrifuge to remove glycerol
    C402 = bst.LiquidsSplitCentrifuge('C402',
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
    @T404.add_specification
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
    
    ### Biodiesel Purification Section ###
    
    # Water wash
    T405 = bst.MixTank('T405')
    P406 = bst.Pump('P406')
    
    # Centrifuge out water
    C403 = bst.LiquidsRatioCentrifuge('C403',
                             K_chemicals=('Methanol', 'Glycerol'),
                             Ks=np.array([0.382, 0.183]),
                             top_solvents=('Biodiesel',),
                             top_split=(0.999,),
                             bot_solvents=('Water', 'TAG', 'NaOH', 'HCl'),
                             bot_split=(0.999, 1, 1, 1))
    
    # Vacuum dry biodiesel
    # Consider doing this by split, and keeping this unit constant
    # 290 Pa, 324 K according to USDA Biodiesel Model
    F401 = bst.SplitFlash('F401',
                    order=('Water', 'Methanol', 'Biodiesel'),
                    split=(0.9999, 0.9999, 0.00001),
                    P=10000., T=331.5, Q=0)
    F401.line = 'Vacuum dryer'
    F401.material = 'Stainless steel 304'
    P407 = bst.Pump('P407', P=101325.)
    
    ### Glycerol Purification Section ###
    
    # Condense vacuumed methanol to recycle
    H401 = bst.HXutility('H401', V=0, T=295)
    P408 = bst.Pump('P408', P=101325.)
    
    # Mix recycled streams and HCl
    T406 = bst.MixTank('T406')
    P409 = bst.Pump('P409', P=101325.)
    
    # Centrifuge out waste fat
    # assume all the lipid, free_lipid and biodiesel is washed out
    C404 = bst.LiquidsSplitCentrifuge('C404', outs=('', wastewater),
                                         order=('Methanol', 'Glycerol', 'Water'),
                                         split=(0.999, 0.999, 0.99))
    
    # Add and mix NaOH
    T407 = bst.MixTank('T407')
    P410 = bst.Pump('P410')
    
    # Methanol/Water distillation column
    D401 = bst.BinaryDistillation('D401',
                      LHK=('Methanol', 'Water'), P=101325,
                      Lr=0.999, Hr=0.99999, k=2.5,
                      is_divided=True,
                      partial_condenser=False,
                      vessel_material='Stainless steel 304',
                      tray_material='Stainless steel 304')
    
    # Save column utilities
    H402 = bst.HXprocess('H402')
    
    # Glycerol/Water flash (not a distillation column)
    chemicals = H402.chemicals
    w = 0.20/chemicals.Water.MW
    g = 0.80/chemicals.Glycerol.MW
    x_water = w/(w+g)
    # ideal_thermo = D401.thermo.ideal()
    D402 = bst.BinaryDistillation('D402',
                        LHK=('Water', 'Glycerol'),
                        k=1.25,
                        P=101325,
                        y_top=0.999999,
                        x_bot=x_water,
                        tray_material='Stainless steel 304',
                        vessel_material='Stainless steel 304',
                        partial_condenser=False)
    P413 = bst.Pump('P413', P=101325.,)
    
    @D402.add_specification
    def startup_water():
        imol = D402.ins[0].imol
        glycerol = imol['Glycerol']
        imol['Water'] = water = (800. * C402.outs[0].imol['Glycerol'] + w / g * glycerol)
        z_distillate_LK = D402.y_top
        z_bottoms_LK = D402.x_bot
        z_feed_LK = water / (water + glycerol)
        if (z_bottoms_LK < z_feed_LK < z_distillate_LK):
            D402._run()
            # Remove accumulation
            D402.outs[0].imol['Water'] = 800. * C402.outs[0].imol['Glycerol']
            D402.outs[0].T = D402.outs[0].bubble_point_at_P().T
        else:
            z_feed_LK = 0.1 * z_bottoms_LK + 0.9 * z_distillate_LK
            imol['Water'] = water = glycerol / (1 - z_feed_LK)
            D402._run()
            # Remove accumulation
            D402.outs[0].imol['Water'] = water - w / g * glycerol
            D402.outs[0].T = D402.outs[0].bubble_point_at_P().T
            
    
    # Condense recycle methanol
    H403 = bst.HXutility('H403', V=0, T=315)
    P411 = bst.Pump('P411')
    
    # Condense recycle water
    H404 = bst.HXutility('H404', V=0, T=315)
    P412 = bst.Pump('P412')
    
    # Storage tank for glycerol
    T408 = bst.StorageTank('T408', outs=crude_glycerol)
    
    # Storage tank for biodiesel
    T409 = bst.StorageTank('T409', outs=biodiesel)
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
    @T405.add_specification(prioritize=True)
    def adjust_biodiesel_wash_water():
        total_glycerol =  (C401.outs[1].imol['Glycerol'] + R402.outs[0].imol['Glycerol'])
        wash_water = (x_water / (1 - x_water) * total_glycerol
                      - HCl1.imol['Water']
                      - NaOH.imol['Water']
                      - oil.imol['Water']
                      - HCl2.imol['Water'])
        biodiesel_wash_water.imol['Water'] = wash_water if wash_water > 0 else 0.
        T405._run()
    
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
    
    @B401.add_specification
    def vary_inlets():
        B401.run()
        catalyst.sink.run_until(T404)
        methanol.sink.run_until(T404)
        T404.run_until(S401)
        S401.ins[0].mix_from([1**R401, 1**R402])
    
    # Initial guess
    D402.outs[0].imol['Methanol', 'Glycerol', 'Water'] = [0.00127, 3.59e-05, 0.0244]
    
    bst.mark_disjunction(P411-0)
    