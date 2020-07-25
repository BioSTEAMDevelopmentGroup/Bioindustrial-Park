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
from biosteam import units
from ._process_settings import price
from ..sugarcane import create_ethanol_production_system

__all__ = ('create_system',)

# %% Pretreatment section

def create_system(ID='lipidcane_sys'):
    chemicals = bst.settings.get_chemicals()
    s = bst.main_flowsheet.stream
    u = bst.main_flowsheet.unit
    
    ### Streams ###
    
    lipidcane = lipidcane = bst.Stream('lipidcane',
                                       Ash=2000.042,
                                       Cellulose=26986.69,
                                       Glucose=2007.067,
                                       Hemicellulose=15922.734,
                                       Lignin=14459.241,
                                       Lipid=10035.334,
                                       Solids=5017.667,
                                       Sucrose=22746.761,
                                       Water=234157.798,
                                       units='kg/hr',
                                       price=price['Lipid cane'])
    
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
    
    oil_wash_water = bst.Stream('oil_wash_water',
                                Water=1350, units='kg/hr',
                                T=358.15)  # to T207
    
    ### Unit operations ###
    
    bst.Stream.ticket_name = 'd'
    bst.Stream.ticket_number = 100
    
    # Feed the shredder
    U101 = units.ConveyingBelt('U101', ins=lipidcane)
    U101.cost_items['Conveying belt'].ub = 5000
    
    # Separate metals
    U102 = units.MagneticSeparator('U102', ins=U101-0)
    
    # Shredded cane
    U103 = units.Shredder('U103', ins=U102-0)
    
    bst.Stream.ticket_number = 200
    
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
                                         Lipid=0.1,
                                         Solids=1),
                              moisture_content=0.5)
    
    # Convey out bagasse
    U202 = units.ConveyingBelt('U202', ins=U201.outs[0], outs='Bagasse')
    
    # Mix in water
    M201 = units.Mixer('M201')
    
    # Screen out fibers
    S201 = units.VibratingScreen('S201',
                                 split=dict(Ash=0.35,
                                            Cellulose=0.35,
                                            Glucose=0.88,
                                            Hemicellulose=0.35,
                                            Lignin=0.35,
                                            Lipid=0.88,
                                            Solids=0,
                                            Sucrose=0.88,
                                            Water=0.88))
    
    # Store juice before treatment
    T202 = units.StorageTank('T202', tau=4, vessel_material='Carbon steel')
    
    # Heat up before adding acid
    H201 = units.HXutility('H201', T=343.15, V=0)
    
    # Mix in acid
    T203 = units.MixTank('T203')
    T203.tau = 0.10
    
    # Pump acid solution
    P201 = units.Pump('P201')
    
    # Mix lime solution
    T204 = units.MixTank('T204')
    T204.tau = 0.10
    P202 = units.Pump('P202')
    
    # Blend acid lipid solution with lime
    T205 = units.MixTank('T205')
    
    # Mix recycle
    M202 = units.Mixer('M202')
    
    # Heat before adding flocculant
    H202 = units.HXutility('H202', T=372.15, V=0)
    
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
                                      Lipid=0.98,
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
    
    # Separate oil and sugar
    T207 = units.MixTank('T207')
    T207_2 = units.Splitter('T207_2',
                            split=dict(Lipid=1,
                                       Water=1e-4))
    
    # Cool the oil
    H203 = units.HXutility('H203', T=343.15, V=0)
    
    # Screen out small fibers from sugar stream
    S202 = units.VibratingScreen('S202', outs=('', 'fiber_fines'),
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
    sugar = S202-0
    S202.mesh_opening = 2
    
    # Add distilled water to wash lipid
    T208 = units.MixTank('T208')
    T208.tau = 0.10
    
    # Centrifuge out water
    C203 = units.LiquidsSplitCentrifuge('C203',
                                        split=dict(Lipid=0.99,
                                                   Water=0.01))
    
    # Vacume out water
    F201 = units.SplitFlash('F201', T=357.15, P=2026.5,
                            outs=('water_vapor', ''),
                            split=dict(Lipid=0.0001,
                                       Water=0.999))
    
    ### Process specifications ###
    def correct_flows():
        F_mass = lipidcane.F_mass
        # correct enzyme, lime, phosphoric acid, and imbibition water
        enzyme.imass['Cellulose', 'Water'] = 0.003 * F_mass * np.array([0.1, 0.9])
        lime.imass['CaO', 'Water'] = 0.001 * F_mass * np.array([0.046, 0.954])
        H3PO4.imass['H3PO4', 'Water'] = 0.00025 * F_mass
        imbibition_water.imass['Water'] = 0.25* F_mass
        T201._run()
    
    T201.specification = correct_flows
    
    # Specifications within a system
    def correct_lipid_wash_water():
        oil_wash_water.imol['Water'] = 100/11 * H202.outs[0].imol['Lipid']
        T208._run()
    
    T208.specification = correct_lipid_wash_water
    
    def correct_wash_water():
        P202._run()
        solids = P202.outs[0].imol['Ash', 'CaO', 'Cellulose', 'Hemicellulose', 'Lignin'].sum()
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
    
    C201-0-T207-T207_2-0-H203
    (H203-0, oil_wash_water)-T208-C203-0-F201
    T207-T207_2-1-S202
    
    lipid = F201-1
    
    ### Ethanol section ###
    
    ### Utilities ###
    
    MW_etoh = chemicals.Ethanol.MW
    MW_water = chemicals.Water.MW
    
    def mass2molar_ethanol_fraction(ethanol_mass_fraction):
        """Return ethanol mol fraction in a ethanol water mixture"""
        x = ethanol_mass_fraction
        return x/MW_etoh / (x/MW_etoh + (1-x)/MW_water)
    
    ### Input streams ###
    
    ethanol_production_sys = create_ethanol_production_system(sugar_solution=S202-0)
    u.M305.ins.append(C203-1)
    
    ### Biodiesel section ###
    
    ### Streams ###
    
    # Fresh degummed oil
    oil = lipid
    
    # Fresh methanol
    methanol = bst.Stream('methanol', Methanol=1,
                          price=price['Methanol'])
    
    # Catalyst
    catalyst = bst.Stream('catalyst', NaOCH3=0.25,
                          Methanol=0.75, units='kg/hr',
                          price=price['NaOCH3'])
                      # price=0.25*price['NaOCH3'] + 0.75*Methanol.price)
    
    # Water to remove glycerol
    biodiesel_wash_water = bst.Stream('biodiesel_wash_water', Water=13.6, T=273.15+60, 
                                      price=price['Water'])
    
    # Acid to neutralize catalyst after second centrifuge
    HCl1 = bst.Stream('HCl1', HCl=0.21, Water=0.79,
                      price=price['HCl'])
                  # price=0.21*price['HCl'] + 0.79*Water.price) # 35% HCl by mass
    
    # Acid to remove soaps after first centrifuge
    HCl2 = bst.Stream('HCl2', HCl=0.21, Water=0.79,
                      price=HCl1.price)
    
    # Base to neutralize acid before distillation
    NaOH = bst.Stream('NaOH', NaOH=1, price=price['NaOH'])
    
    # Products
    crude_glycerol = bst.Stream('crude_glycerol',
                                price=price['Crude glycerol'])
    biodiesel = bst.Stream('biodiesel',
                           price=price['Biodiesel'])
    
    # Waste
    waste = bst.Stream('waste', price=price['Waste'])
    
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
    R401 = units.Transesterification('R401', efficiency=0.90, methanol2lipid=6, T=333.15,
                             catalyst_molfrac=x_cat)
    
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
    R402 = units.Transesterification('R402', efficiency=0.90, methanol2lipid=6, T=333.15,
                             catalyst_molfrac=x_cat) 
    
    PS5 = S401.create_reversed_splitter_process_specification(
        'PS5', description='Adjust feed to reactors')
    
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
        catalyst_mol = T404.outs[0].imol['NaOCH3']
        NaOH.imol['NaOH'] = k1 * catalyst_mol
        HCl1.imol['HCl'] = k2 * catalyst_mol
        HCl2.imol['HCl'] = k3 * catalyst_mol
    
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
                             bot_solvents=('Water', 'Lipid', 'NaOH', 'HCl'),
                             bot_split=(0.999, 1, 1, 1))
    
    # Vacuum dry biodiesel
    # Consider doing this by split, and keeping this unit constant
    # 290 Pa, 324 K according to USDA Biodiesel Model
    F401 = units.SplitFlash('F401',
                    order=('Water', 'Methanol', 'Biodiesel'),
                    split=(0.9999, 0.9999, 0.00001),
                    P=2026.5, T=331.5, Q=0)
    F401.line = 'Vacuum dryer'
    F401.material = 'Stainless steel 304'
    P407 = units.Pump('P407')
    
    ### Glycerol Purification Section ###
    
    # Condense vacuumed methanol to recycle
    H401 = units.HXutility('H401', V=0, T=295)
    P408 = units.Pump('P408')
    
    # Mix recycled streams and HCl
    T406 = units.MixTank('T406')
    P409 = units.Pump('P409')
    
    # Centrifuge out waste fat
    # assume all the lipid, free_lipid and biodiesel is washed out
    C404 = units.LiquidsSplitCentrifuge('C404', outs=('', waste),
                                         order=('Methanol', 'Glycerol', 'Water'),
                                         split=(0.999, 0.999, 0.999))
    
    # Add and mix NaOH
    T407 = units.MixTank('T407')
    P410 = units.Pump('P410')
    
    # Methanol/Water distillation column
    D401 = units.BinaryDistillation('D401',
                      LHK=('Methanol', 'Water'), P=101325,
                      y_top=0.99999, x_bot=0.0001, k=2.5,
                      is_divided=True,
                      vessel_material='Stainless steel 304',
                      tray_material='Stainless steel 304')
    
    # Condense water to recycle
    H402 = units.HXutility('H402', T=353, V=0)
    
    # Glycerol/Water flash (not a distillation column)
    w = 0.20/chemicals.Water.MW
    g = 0.80/chemicals.Glycerol.MW
    x_water = w/(w+g)
    
    D402 = units.BinaryDistillation('D402',
                        LHK=('Water', 'Glycerol'),
                        k=1.25,
                        P=101325,
                        y_top=0.999,
                        x_bot=x_water,
                        tray_material='Stainless steel 304',
                        vessel_material='Stainless steel 304')
    
    def startup_water():
        imol = D402.ins[0].imol
        water, glycerol = imol['Water', 'Glycerol']
        minimum_water = 5 * (w / (w + g)) * glycerol
        if water < minimum_water:
            imol['Water'] = minimum_water
        D402._run()
        # Remove accumulation
        D402.outs[0].imol['Water'] = 1100*C402.outs[0].imol['Glycerol']
            
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
    (C401-0, S401-1)-R402-0-PS5-C402-1
    
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
    H402-0-D401-1-D402-1-T408
    P410-0-H402                
    
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
    
    ### Facilities ###
    
    # Burn bagasse from conveyor belt
    BT = units.BoilerTurbogenerator('BT',
                                   (U202-0, '', 'boiler_makeup_water', 'natural_gas', '', ''),
                                   boiler_efficiency=0.80,
                                   turbogenerator_efficiency=0.85)
    
    CT = units.CoolingTower('CT')
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    
    process_water_streams = (s.imbibition_water,
                             s.biodiesel_wash_water,
                             oil_wash_water,
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
    
    return bst.main_flowsheet.create_system(ID, ends=(*S401.outs, P408-0))
    
    
