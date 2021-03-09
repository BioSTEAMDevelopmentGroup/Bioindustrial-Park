#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 11:05:24 2017

@author: Yoel
"""
import numpy as np
import biosteam as bst
from biosteam import units, SystemFactory
from biosteam import main_flowsheet as f

__all__ = (
    'create_feedstock_handling_system',
    'create_juicing_system_up_to_clarification',
    'create_juicing_system_with_fiber_screener',
    'create_sucrose_fermentation_system',
    'create_beer_distillation_system',
    'create_ethanol_purification_system_after_beer_column',
    'create_sucrose_to_ethanol_system',
    'create_ethanol_purification_system',
    'create_sugarcane_to_ethanol_system',
)


# %% Juicing and evaporation

@SystemFactory(
    ID='feedstock_handling_sys',
    ins=[dict(ID='sugarcane',
              Water=0.7,
              Glucose=0.01208,
              Sucrose=0.1369,
              Ash=0.006,
              Cellulose=0.06115,
              Hemicellulose=0.03608,
              Lignin=0.03276,
              Solids=0.015,
              total_flow=333334.2,
              units='kg/hr',
              price=0.03455)],
    outs=[dict(ID='shreaded_cane')]
)
def create_feedstock_handling_system(ins, outs):
    sugarcane, = ins
    shreaded_cane, = outs
    U101 = units.ConveyingBelt('U101', sugarcane)
    U102 = units.MagneticSeparator('U102', U101-0)
    U103 = units.Shredder('U103', U102-0, shreaded_cane)

@SystemFactory(
    ID='juicing_sys',
    ins=[create_feedstock_handling_system.ins[0],
         dict(ID='enzyme',
              Cellulose=100,
              Water=900,
              units='kg/hr',
              price=0.5),
         dict(ID='H3PO4',
              H3PO4=74.23,
              Water=13.1,
              units='kg/hr',
              price=0),
         dict(ID='lime',
              CaO=333.0,
              Water=2200.0,
              units='kg/hr',
              price=0.077),
         dict(ID='polymer',
              Flocculant=0.83,
              units='kg/hr',
              price=0)],
    outs=[dict(ID='clarified_juice'),
          dict(ID='bagasse')]
)
def create_juicing_system_up_to_clarification(ins, outs):
    ### Streams ###
    sugarcane, enzyme, H3PO4, lime, polymer = ins
    clarified_juice, bagasse = outs
    
    imbibition_water = bst.Stream('imbibition_water',
                                  Water=87023.35, units='kg/hr',
                                  T = 338.15)
    
    rvf_wash_water = bst.Stream('rvf_wash_water',
                                Water=16770, units='kg/hr',
                                T=363.15)  # to C202
    
    ### Unit operations ###
    
    
    
    # Hydrolyze starch
    T201 = units.EnzymeTreatment('T201', (sugarcane, enzyme), T=323.15)  # T=50
    
    # Finely crush lipid cane
    U201 = units.CrushingMill('U201',
                              ins=(T201-0, ''),
                              split=dict(Ash=0.92,
                                         Cellulose=0.92,
                                         Glucose=0.04,
                                         Hemicellulose=0.92,
                                         Lignin=0.92,
                                         Sucrose=0.04,
                                         Solids=1),
                              moisture_content=0.5)
    
    # Convey out bagasse
    U202 = units.ConveyingBelt('U202', U201-0, bagasse)
    
    # Mix in water
    M201 = units.Mixer('M201', ('', imbibition_water), 1-U201)
    
    # Screen out fibers
    S201 = units.VibratingScreen('S201', U201-1, ('', 0-M201),
                                 split=dict(Ash=0.35,
                                            Cellulose=0.35,
                                            Glucose=0.88,
                                            Hemicellulose=0.35,
                                            Lignin=0.35,
                                            Solids=0,
                                            Sucrose=0.88,
                                            Water=0.88))
    
    # Store juice before treatment
    T202 = units.StorageTank('T202', S201-0,
                             tau=4, vessel_material='Carbon steel')
    
    # Heat up before adding acid
    H201 = units.HXutility('H201', T202-0, T=343.15)
    
    # Mix in acid
    T203 = units.MixTank('T203', (H201-0, H3PO4))
    
    # Pump acid solution
    P201 = units.Pump('P201', T203-0)
    
    # Mix lime solution
    T204 = units.MixTank('T204', lime, tau=0.10)
    
    # Blend acid lipid solution with lime
    T205 = units.MixTank('T205', (P201-0, T204-0), tau=0.10)
    P202 = units.Pump('P202', T205-0)
    
    # Mix recycle
    M202 = units.Mixer('M202', (P202-0, ''))
    
    # Heat before adding flocculant
    H202 = units.HXutility('H202', M202-0, T=372.15)
    
    # Mix in flocculant
    T206 = units.MixTank('T206', (H202-0, polymer))
    T206.tau = 0.10
    
    # Separate residual solids
    C201 = units.Clarifier('C201', T206-0, (clarified_juice, ''),
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
    C202 = units.RVF('C202', (C201-1, rvf_wash_water), ('filter_cake', ''),
                     moisture_content=0.80,
                     split=dict(Ash=0.85,
                                CaO=0.85,
                                Cellulose=0.85,
                                Glucose=0.01,
                                Hemicellulose=0.85,
                                Lignin=0.85,
                                Sucrose=0.01))
    P203 = units.Pump('P203', C202-1, 1-M202)
    
    ### Process specifications ###
    
    # Specifications dependent on lipid cane flow rate
    def correct_flows():
        F_mass = T201.ins[0].F_mass
        # correct enzyme, lime, phosphoric acid, and imbibition water
        enzyme.imass['Cellulose', 'Water'] = 0.003 * F_mass * np.array([0.1, 0.9])
        lime.imass['CaO', 'Water'] = 0.001 * F_mass * np.array([0.046, 0.954])
        H3PO4.imass['H3PO4', 'Water'] = 0.00025 * F_mass
        imbibition_water.imass['Water'] = 0.25* F_mass
        T201._run()
    
    T201.specification = correct_flows
    
    # Specifications within a system
    def correct_wash_water():
        P202._run()
        solids = P202.outs[0].imol['Ash', 'CaO', 'Cellulose',
                                   'Hemicellulose', 'Lignin'].sum()
        rvf_wash_water.imol['Water'] = 0.0574 * solids
    
    P202.specification = correct_wash_water

@SystemFactory(
    ID='juicing_sys',
    ins=create_juicing_system_up_to_clarification.ins,
    outs=[dict(ID='screened_juice'),
          dict(ID='bagasse'),
          dict(ID='fiber_fines')]
)          
def create_juicing_system_with_fiber_screener(ins, outs):
    screened_juice, bagasse, fiber_fines = outs
    sys = create_juicing_system_up_to_clarification(None, ins, ['', bagasse], mockup=True)

    # Screen out small fibers from sugar stream
    S202 = units.VibratingScreen('S202', sys-0, (screened_juice, fiber_fines),
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


# %% Ethanol separation

@SystemFactory(
    ID='beer_distillation_sys',
    ins=[dict(ID='beer',
              T=348.34,
              P=101325,
              Water=96500.0,
              Ethanol=22550.0,
              Glucose=4916,
              H3PO4=83.33,
              Yeast=103,
              units='kg/hr')],
    outs=[dict(ID='distilled_beer'),
          dict(ID='stillage')]
)
def create_beer_distillation_system(ins, outs,
                                    beer_column_heat_integration=True,
                                    IDs={}):
    beer, = ins
    distilled_beer, stillage = outs
    
    P301 = units.Pump(IDs.get('Beer pump', 'P301'))
    
    # Beer column
    x_bot = 3.910570816782338e-06
    y_top = 0.2811210085806504
    D302 = units.BinaryDistillation(IDs.get('Beer column', 'D302'), P=2. * 101325,
                                    outs=(distilled_beer, ''),
                                y_top=y_top, x_bot=x_bot, k=1.1, Rmin=0.001,
                                LHK=('Ethanol', 'Water'))
    D302.tray_material = 'Stainless steel 304'
    D302.vessel_material = 'Stainless steel 304'
    D302.boiler.U = 1.85
    P302 = units.Pump(IDs.get('Beer column bottoms product pump', 'P302'))
    
    # Heat up before beer column
    # Exchange heat with stillage    
    if beer_column_heat_integration:
        H302 = units.HXprocess(IDs.get('Beer column heat exchange', 'H302'), 
                               outs=('', stillage),
                               phase0='l', phase1='l', U=1.28)
        (beer-P301-0, P302-0)-H302-0-D302-1-P302
    else:
        beer-P301-0-D302-1-P302
        P302.outs[0] = stillage


@SystemFactory(
    ID='ethanol_purification_from_distilled_beer_sys',
    ins=[dict(ID='distilled_beer'),
         dict(ID='denaturant',
              Octane=230.69,
              units='kg/hr',
              price=0.756)],
    outs=[dict(ID='ethanol',
               price=0.789),
          dict(ID='stripper_bottoms_product')]
)
def create_ethanol_purification_system_after_beer_column(ins, outs, IDs={}):
    distilled_beer, denaturant = ins
    ethanol, stripper_bottoms_product = outs
  
    
    # Mix ethanol Recycle (Set-up)
    M303 = units.Mixer(IDs.get('Recycle mixer', 'M303'))
    
    D303 = units.BinaryDistillation(IDs.get('Distillation', 'D303'), P=101325,
                                x_bot=3.9106e-06, y_top=0.80805, k=1.25, Rmin=0.01,
                                LHK=('Ethanol', 'Water'),
                                tray_material='Stainless steel 304',
                                vessel_material='Stainless steel 304',
                                is_divided=True)
    D303.boiler.U = 1.85
    P303 = units.Pump(IDs.get('Distillation bottoms product pump', 'P303'), 
                      outs=stripper_bottoms_product)
    
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
                             tau=6.5*24, outs=ethanol)
    
    ### Ethanol system set-up ###
    (distilled_beer, U301-0)-M303-0-D303-0-H303-U301
    D303-1-P303
    
    def adjust_denaturant():
        P304._run()
        pure_ethanol = P304.outs[0]
        denaturant.imol['Octane'] = 0.022*pure_ethanol.F_mass/114.232
    
    P304.specification = adjust_denaturant
    (denaturant-T303-P305-0, U301-1-H304-0-T302-0-P304-0)-M304-T304

# @SystemFactory(
#     ID='ethanol_purification_from_distilled_beer_sys',
#     ins=[dict(ID='distilled_beer'),
#           dict(ID='denaturant', Octane=230.69, units='kg/hr', price=0.756)],
#     outs=[dict(ID='ethanol', price=0.789),
#           dict(ID='stripper_bottoms_product')]
# )
# def create_ethanol_purification_system_after_beer_column(ins, outs):
#     distilled_beer, denaturant = ins
#     ethanol, stripper_bottoms_product = outs
#     M303 = units.Mixer('M303')
#     D303 = units.BinaryDistillation('D303',
#         P=101325, x_bot=3.9106e-06, y_top=0.80805, k=1.25, Rmin=0.01,
#         LHK=('Ethanol', 'Water'), tray_material='Stainless steel 304',
#         vessel_material='Stainless steel 304', is_divided=True
#     )
#     D303.boiler.U = 1.85
#     P303 = units.Pump('P303', outs=stripper_bottoms_product)
#     H303 = units.HXutility('H303', T=115+273.15, V=1)
#     U301 = units.MolecularSieve('U301', split=dict(Ethanol=0.162, Water=0.925))
#     H304 = units.HXutility('H304', V=0, T=340.)
#     T302 = units.StorageTank('T302', 
#         tau=12, vessel_type='Floating roof', vessel_material='Carbon steel'
#     )
#     P304 = units.Pump('P304')
#     T303 = units.StorageTank('T303', 
#         tau=7*24, vessel_type='Floating roof', vessel_material='Carbon steel'
#     )
#     P305 = units.Pump('P305'); M304 = units.Mixer('M304')
#     T304 = units.StorageTank('T304',
#         vessel_type='Floating roof', vessel_material='Carbon steel',
#         tau=6.5*24, outs=ethanol
#     )
#     def adjust_denaturant():
#         P304._run()
#         denaturant.imol['Octane'] = 0.022 * P304.outs[0].F_mass / 114.232
    
#     P304.specification = adjust_denaturant
#     (distilled_beer, U301-0)-M303-0-D303-0-H303-U301
#     D303-1-P303
#     (denaturant-T303-P305, U301-1-H304-0-T302-0-P304-0)-M304-T304


@SystemFactory(
    ID='ethanol_purification_sys',
    ins=[dict(ID='beer',
              T=348.34,
              P=101325,
              Water=96500.0,
              Ethanol=22550.0,
              Glucose=4916,
              H3PO4=83.33,
              Yeast=103,
              units='kg/hr'),
         create_ethanol_purification_system_after_beer_column.ins[1]],
    outs=[create_ethanol_purification_system_after_beer_column.outs[0],
          dict(ID='stillage'),
          create_ethanol_purification_system_after_beer_column.outs[1]]
)
def create_ethanol_purification_system(ins, outs,
                                       beer_column_heat_integration=True,
                                       IDs={}):
    beer, denaturant = ins
    ethanol, stillage, stripper_bottoms_product = outs
    distilled_beer = bst.Stream('')
    create_beer_distillation_system(
        ins=beer,
        outs=[distilled_beer, stillage],
        beer_column_heat_integration=beer_column_heat_integration,
        IDs=IDs,
        mockup=True,
    )
    create_ethanol_purification_system_after_beer_column(
        ins=[distilled_beer, denaturant],
        outs=[ethanol, stripper_bottoms_product],
        IDs=IDs,
        mockup=True,
    )
    

# %% Ethanol production section (fermentation and separations)

@SystemFactory(
    ID='sucrose_fermentation_sys',
    ins=[dict(ID='screened_juice', 
              Glucose=3802,
              Sucrose=4.309e+04,
              Water=2.59e+05,
              H3PO4=83.33,
              units='kg/hr',
              T=372)],
    outs=[dict(ID='beer'),
          dict(ID='evaporator_condensate')],
)
def create_sucrose_fermentation_system(ins, outs):
    screened_juice, = ins
    beer, evaporator_condensate = outs
    
    ### Streams ###
    
    # Fresh water
    stripping_water = bst.Stream('stripping_water',
                                 Water=26836,
                                 units='kg/hr')
    
    ### Units ###
    
    # Split sugar solution
    S301 = units.Splitter('S301',
                          split=0.265)
    
    # Concentrate sugars
    F301 = units.MultiEffectEvaporator('F301',
                                       P=(101325, 69682, 47057, 30953, 19781),
                                       outs=('', evaporator_condensate),
                                       V=0.85) # fraction evaporated
    
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
    
    stripping_water_over_vent = stripping_water.mol / 21202.490455845436
    def update_stripping_water():
        stripping_water, vent = D301.ins
        stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
    
    D301 = units.VentScrubber('D301', ins=(stripping_water, R301-0), 
                              outs=('vent', ''),
                              gas=('CO2', 'O2'))
    
    # Separate 99% of yeast
    C301 = units.SolidsCentrifuge('C301', 
                                  split=(1-1e-6, 0.99, 1, 0.01),
                                  order=('Ethanol', 'Glucose', 'H3PO4', 'DryYeast'),
                                  solids=('DryYeast',))
    C301.split[:] = 1. - C301.split
    
    # Mix in Water
    M302 = units.Mixer('M302', outs=beer)
    
    # Yeast mixing
    T305 = units.MixTank('T305')
    T305.tau = 0.1
    
    # Multi-effect evaporator pumps
    P306 = units.Pump('P306')
    
    def adjust_yeast_recycle():
        recycle, beer = C301.outs
        feed, *_ = C301.ins
        yeast = 0.1 * feed.F_mass
        m = C301.moisture_content
        recycle.imass['Yeast', 'Water'] = [yeast, yeast * m / (1 - m)]
        beer.mol = feed.mol - recycle.mol
        beer.mol[beer.mol < 0.] = 0.
        beer.T = recycle.T = feed.T
    C301.specification = adjust_yeast_recycle
    
    ### Ethanol system set-up ###
    
    screened_juice-S301-1-F301-0-P306
    (S301-0, P306-0)-M301-H301
    (H301-0, C301-0-T305-0)-R301-1-T301-0-C301
    (C301-1, D301-1)-M302
    

@SystemFactory(
    ID='sucrose_to_ethanol_sys',
    ins=[create_sucrose_fermentation_system.ins[0],
         create_ethanol_purification_system.ins[1]], # denaturant
    outs=[*create_ethanol_purification_system.outs,
          create_sucrose_fermentation_system.outs[0]]
)
def create_sucrose_to_ethanol_system(ins, outs):
    screened_juice, denaturant = ins
    ethanol, stillage, stripper_bottoms_product, evaporator_condensate = outs
    
    beer = bst.Stream()
    
    create_sucrose_fermentation_system(
        ins=screened_juice,
        outs=[beer, evaporator_condensate],
        mockup=True
    )
    create_ethanol_purification_system(
        ins=[beer, denaturant], 
        outs=[ethanol, stillage, stripper_bottoms_product],
        mockup=True
    )
    
    

# %% Complete system

@SystemFactory(
    ID='sugarcane_sys', 
    ins=[*create_juicing_system_with_fiber_screener.ins,
         create_ethanol_purification_system.ins[1]], # denaturant
    outs=[create_ethanol_purification_system.outs[0], # ethanol
          dict(ID='wastewater'),
          dict(ID='emissions'),
          dict(ID='ash_disposal')]
)
def create_sugarcane_to_ethanol_system(ins, outs,
        evaporator_and_beer_column_heat_integration=True
    ):
    s = f.stream
    u = f.unit
    
    sugarcane, enzyme, H3PO4, lime, polymer, denaturant = ins
    ethanol, wastewater, emissions, ash_disposal = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=[sugarcane],
        outs=[''],
    )
    juicing_sys = create_juicing_system_with_fiber_screener(
        ins=[feedstock_handling_sys-0, enzyme, H3PO4, lime, polymer],
        mockup=True
    )
    ethanol_production_sys = create_sucrose_to_ethanol_system(
        ins=(juicing_sys-0, denaturant), outs=ethanol,
        mockup=True
    )
    M305 = units.Mixer('M305', 
        ins=(juicing_sys-2, *ethanol_production_sys-[1, 2, 3]),
        outs=wastewater
    )
    
    ### Facilities ###    
    
    BT = units.BoilerTurbogenerator('BT',
        (juicing_sys-1, '', 'boiler_makeup_water', 'natural_gas', '', ''),
        outs=(emissions, 'rejected_water_and_blowdown', ash_disposal),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85
    )
    CT = units.CoolingTower('CT')
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    process_water_streams = (s.imbibition_water,
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
    if evaporator_and_beer_column_heat_integration:
        def heat_integration():
            hu_mee = F301.heat_utilities[0]
            hu_dist = D303.heat_utilities[0]
            actual_duty = hu_mee.duty + hu_dist.duty
            if actual_duty > 0.:
                hu_mee(actual_duty, 373.15, 373.15)
                hu_dist.empty()
            else:
                hu_mee.empty()
                condenser = D303.condenser
                hu_dist(actual_duty, condenser.ins[0].T, condenser.outs[0].T)
        CWP.specification = heat_integration
        
        
