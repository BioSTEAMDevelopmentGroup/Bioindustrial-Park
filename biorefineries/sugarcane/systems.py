#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 11:05:24 2017
@author: Yoel
"""
import numpy as np
import biosteam as bst
import flexsolve as flx
from biosteam import units, SystemFactory, stream_kwargs as skw
from biosteam import main_flowsheet as f
import thermosteam as tmo

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
    'create_bagasse_pelleting_system',
    'create_sugar_crystallization_system',
    'create_sugarcane_to_sugar_and_molasses_system',
    'convert_fiber_to_lignocelluosic_components',
    'set_sugarcane_composition',
)


bagasse = skw('bagasse')
bagasse_pellets = skw('bagasse_pellets')
sugarcane = skw('sugarcane',
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
    price=0.03455
)
shredded_cane = skw('shredded_cane')
untreated_juice = skw('untreated_juice')
H3PO4 = skw('H3PO4',
    H3PO4=74.23,
    Water=13.1,
    units='kg/hr',
    price=0
)
lime = skw('lime',
    CaO=333.0,
    Water=2200.0,
    units='kg/hr',
    price=0.077
)
polymer = skw('polymer',
    Flocculant=0.83,
    units='kg/hr',
    price=0
)
clarified_juice = skw('clarified_juice')
screened_juice = skw('screened_juice', 
    Glucose=3802,
    Sucrose=4.309e+04,
    Water=2.59e+05,
    H3PO4=83.33,
    units='kg/hr',
    T=372
)
fiber_fines = skw('fiber_fines')
beer = skw('beer',
    T=348.34,
    P=101325,
    Water=96500.0,
    Ethanol=22550.0,
    Glucose=4916,
    H3PO4=83.33,
    Yeast=103,
    units='kg/hr'
)
distilled_beer = skw('distilled_beer')
stillage = skw('stillage')
denaturant = skw('denaturant',
    Octane=230.69,
    units='kg/hr',
    price=0.756
)
ethanol = skw('ethanol',
    price=0.789
)
stripper_bottoms_product = skw('stripper_bottoms_product')
evaporator_condensate = skw('evaporator_condensate')
vent = skw('vent')
vinasse = skw('vinasse')
wastewater = skw('wastewater')
emissions = skw('emissions')
ash_disposal = skw('ash_disposal')
molasses = skw('molasses')
sugar = skw('sugar', price=0.419) # https://markets.businessinsider.com/commodities/sugar-price?op=1 (3/18/2022)
      
def set_sugarcane_composition(stream, water, fiber, sugar):
    chemicals = stream.chemicals
    if 'Sugar' not in chemicals:
        IDs = ('Sucrose', 'Glucose')
        chemicals.define_group('Sugar', IDs, composition=stream.imass[IDs], wt=True)
    if 'Fiber' not in chemicals:
        IDs = ('Cellulose', 'Hemicellulose', 'Lignin')
        chemicals.define_group('Fiber', IDs, composition=stream.imass[IDs], wt=True)
    if 'Other' not in chemicals:
        IDs = ('Ash', 'Solids')
        chemicals.define_group('Other', IDs, composition=stream.imass[IDs], wt=True)
    F_mass = stream.F_mass
    other = 1 - water - fiber - sugar
    assert other > 0, (water, fiber, sugar, other)
    stream.imass['Water', 'Fiber', 'Sugar', 'Other'] = F_mass * np.array([water, fiber, sugar, other])

def convert_fiber_to_lignocelluosic_components(stream, ignore_acetate=False):
    chemicals = stream.chemicals
    if chemicals is convert_fiber_to_lignocelluosic_components.last_chemicals:
        prxn = convert_fiber_to_lignocelluosic_components.last_reaction
    else:
        cellulose_rxn = tmo.Reaction('Cellulose -> Glucan', 'Cellulose', 1.0,
                                     basis='wt', chemicals=chemicals)
        cellulose_rxn.basis = 'mol'
        # Bagasse composition https://www.sciencedirect.com/science/article/pii/S0144861710005072
        # South american; by HPLC
        # Glucan: 41.3%
        # Xylan: 24.9%
        # Galactan: 0.6%
        # Arabinan: 1.7%
        # Lignin: 23.2%
        # Acetyl: 3.0%
        if ignore_acetate:
            hemicellulose_rxn = tmo.Reaction(
                '27.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan', 'Hemicellulose', 
                1.0, basis='wt', chemicals=chemicals
            )
        else:
            hemicellulose_rxn = tmo.Reaction(
                '30.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan + 3 Acetate', 'Hemicellulose', 
                1.0, basis='wt', chemicals=chemicals
            )
        hemicellulose_rxn.basis = 'mol'
        convert_fiber_to_lignocelluosic_components.last_chemicals = chemicals
        convert_fiber_to_lignocelluosic_components.last_reaction = prxn = tmo.ParallelReaction(
            [cellulose_rxn, hemicellulose_rxn]
        )
    prxn(stream)

convert_fiber_to_lignocelluosic_components.last_chemicals = None

# %% Juicing and evaporation

@SystemFactory(
    ID='bagasse_pelleting_sys',
    ins=[bagasse],
    outs=[bagasse_pellets]
)
def create_bagasse_pelleting_system(ins, outs):
    bagasse, = ins
    bagasse_pellets, = outs
    U401 = units.HammerMill('U401', bagasse)
    U402 = units.DrumDryer('U402', 
        (U401-0, 'dryer_air', 'dryer_natural_gas'), 
        ('', 'dryer_outlet_air', 'dryer_emissions'),
        moisture_content=0.18, split=0.,
    )
    # X401 = bst.ThermalOxidizer('X401', (U403-1, 'oxidizer_air'), 'oxidizer_emissions')
    U403 = units.ScrewFeeder('U403', U402-0)
    U404 = units.BagassePelletMill('U404', U403-0)
    U405 = units.ConveyingBelt('U405', U404-0, bagasse_pellets)

@SystemFactory(
    ID='feedstock_handling_sys',
    ins=[sugarcane],
    outs=[shredded_cane]
)
def create_feedstock_handling_system(ins, outs):
    sugarcane, = ins
    shreaded_cane, = outs
    U101 = units.ConveyingBelt('U101', sugarcane)
    U102 = units.MagneticSeparator('U102', U101-0)
    U103 = units.Shredder('U103', U102-0, shreaded_cane)

@SystemFactory(
    ID='juicing_sys',
    ins=[sugarcane],
    outs=[untreated_juice, bagasse]
)
def create_juicing_system_without_treatment(ins, outs, pellet_bagasse=None):
    if pellet_bagasse is None: pellet_bagasse = False
    
    ### Streams ###
    sugarcane, = ins
    untreated_juice, bagasse = outs
    
    imbibition_water = bst.Stream('imbibition_water',
                                  Water=87023.35, units='kg/hr',
                                  T = 350.15)
    
    # Finely crush lipid cane
    U201 = units.CrushingMill('U201',
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
        imbibition_water.imass['Water'] = 0.245 * (feed.F_mass - feed.imass['Water']) / 0.7
    
    U202 = units.ConveyingBelt('U202', U201-0, [''] if pellet_bagasse else [bagasse])
    
    if pellet_bagasse:
        bagasse_pelleting_sys = create_bagasse_pelleting_system(None, ins=U202-0, outs=bagasse, mockup=True)
    
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
    T202 = units.StorageTank('T202', S201-0, untreated_juice,
                             tau=4, vessel_material='Carbon steel')


@SystemFactory(
    ID='juicing_sys',
    ins=[sugarcane, H3PO4, lime, polymer],
    outs=[clarified_juice, bagasse]
)
def create_juicing_system_up_to_clarification(ins, outs, pellet_bagasse=None):
    
    ### Streams ###
    sugarcane, H3PO4, lime, polymer = ins
    clarified_juice, bagasse = outs
    
    rvf_wash_water = bst.Stream('rvf_wash_water',
                                Water=16770, units='kg/hr',
                                T=363.15)  # to C202
    
    juicing_sys_without_treatment = create_juicing_system_without_treatment(
        ins=(sugarcane,),
        outs=('', bagasse),
        mockup=True,
        pellet_bagasse=pellet_bagasse,
    )
    U201 = f.unit.U201
    
    # Heat up before adding acid
    H201 = units.HXutility('H201', juicing_sys_without_treatment-0, T=343.15)
    
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
    @U201.add_specification(run=True)
    def correct_flows():
        feed = U201.ins[0]
        F_mass = (feed.F_mass - feed.imass['Water']) / 0.7
        # correct lime, phosphoric acid, and imbibition water
        lime.imass['CaO', 'Water'] = 0.001 * F_mass * np.array([0.046, 0.954])
        H3PO4.imass['H3PO4', 'Water'] = 0.00025 * F_mass
    
    # Specifications within a system
    def correct_wash_water():
        P202._run()
        solids = P202.outs[0].imol['Ash', 'CaO', 'Cellulose',
                                   'Hemicellulose', 'Lignin'].sum()
        rvf_wash_water.imol['Water'] = 0.0574 * solids
    
    P202.specification = correct_wash_water

@SystemFactory(
    ID='juicing_sys',
    ins=[sugarcane, H3PO4, lime, polymer],
    outs=[screened_juice, bagasse, fiber_fines]
)          
def create_juicing_system_with_fiber_screener(ins, outs, pellet_bagasse=None):
    screened_juice, bagasse, fiber_fines = outs
    sys = create_juicing_system_up_to_clarification(
        None, ins, ['', bagasse],
        mockup=True,
        pellet_bagasse=pellet_bagasse
    )
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
    ins=[beer],
    outs=[distilled_beer, stillage] # TODO: Find good selling price for stillage/vinasse and possibly yeast
)
def create_beer_distillation_system(ins, outs,
                                    beer_column_heat_integration=True,
                                    IDs={}):
    beer, = ins
    distilled_beer, stillage = outs
    
    P301 = units.Pump(IDs.get('Beer pump', 'P301'), P=2. * 101325)
    
    # Beer column
    x_bot = 3.91e-06
    y_top = 0.28
    D302 = units.BinaryDistillation(IDs.get('Beer column', 'D302'), P=2. * 101325,
                                    outs=(distilled_beer, ''),
                                y_top=y_top, x_bot=x_bot, k=1.1, Rmin=0.001,
                                LHK=('Ethanol', 'Water'))
    D302.tray_material = 'Stainless steel 304'
    D302.vessel_material = 'Stainless steel 304'
    D302.boiler.U = 1.85
    P302 = units.Pump(IDs.get('Beer column bottoms product pump', 'P302'), P=101325.)
    
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
    ins=[distilled_beer, denaturant],
    outs=[ethanol, stripper_bottoms_product]
)
def create_ethanol_purification_system_after_beer_column(ins, outs, IDs={}):
    distilled_beer, denaturant = ins
    ethanol, stripper_bottoms_product = outs
  
    
    # Mix ethanol Recycle (Set-up)
    M303 = units.Mixer(IDs.get('Recycle mixer', 'M303'))
    
    D303 = units.BinaryDistillation(IDs.get('Distillation', 'D303'), 
                                x_bot=3.9106e-06, y_top=0.80805, k=1.25, Rmin=0.01,
                                LHK=('Ethanol', 'Water'),
                                tray_material='Stainless steel 304',
                                vessel_material='Stainless steel 304',
                                P=1013250.,
                                is_divided=True)
    D303.boiler.U = 1.85
    P303 = units.Pump(IDs.get('Distillation bottoms product pump', 'P303'), 
                      outs=stripper_bottoms_product)
    
    # Superheat vapor for mol sieve
    H303 = units.HXutility(IDs.get('Heat exchanger to superheat vapor to molecular sieves', 'H303'),
                           T=115+273.15, V=1, heat_only=True)
    
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
    M304.denaturant_fraction = 0.022
    
    @M304.add_specification(run=True)
    def adjust_denaturant():
        pure_ethanol = M304.ins[1]
        denaturant.imol['Octane'] = M304.denaturant_fraction * pure_ethanol.F_mass / 114.232
        for i in T303.path_until(M304): i._run()
    
    (denaturant-T303-P305-0, U301-1-H304-0-T302-0-P304-0)-M304-T304

# @SystemFactory(
#     ID='ethanol_purification_from_distilled_beer_sys',
#     ins=[skw('distilled_beer'),
#           skw('denaturant', Octane=230.69, units='kg/hr', price=0.756)],
#     outs=[skw('ethanol', price=0.789),
#           skw('stripper_bottoms_product')]
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
    ins=[beer, denaturant],
    outs=[ethanol, stillage, stripper_bottoms_product]
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
    ID='sugar_crystallization_sys',
    ins=[screened_juice, lime, H3PO4, polymer],
    outs=[sugar, molasses]
)
def create_sugar_crystallization_system(ins, outs):
    # TODO: Add conveyors, storage tanks, packing, sugar remelter
    # https://www.researchgate.net/profile/Maciej-Starzak/publication/311206128_Mass_and_Energy_Balance_Modelling_of_a_Sugar_Mill_A_comparison_of_MATLABR_and_SUGARS_simulations/links/583f240308ae2d217557dcd8/Mass-and-Energy-Balance-Modelling-of-a-Sugar-Mill-A-comparison-of-MATLABR-and-SUGARS-simulations.pdf?origin=publication_detail
    # https://www3.epa.gov/ttn/chief/ap42/ch09/final/c9s10-1a.pdf
    # http://sugartech.co.za/verticalcrystalliser/index.php
    screened_juice, lime, H3PO4, polymer = ins
    sugar, molasses, = outs
    
    if 'Sugar' not in sugar.chemicals:
        sugar.chemicals.define_group('Sugar', ('Glucose', 'Sucrose'))
    
    # Concentrate sugars
    P1 = units.Pump('P1', ins=screened_juice, P=101325)
    
    MEE = units.MultiEffectEvaporator('MEE', P1-0,
        P=(101325, 69682, 47057, 30953, 19781),
        V_definition='First-effect',
        V=0.3
    ) # fraction evaporated
    MEE.brix = 95
    
    def get_brix():
        effluent = MEE.outs[0]
        water = effluent.imass['Water']
        if water < 0.0001: water = 0.0001
        return 100 * effluent.imass['Sugar'] / water
    
    def brix_objective(V):
        MEE.V = V
        MEE._run()
        return MEE.brix - get_brix()
    
    @MEE.add_specification(run=False)
    def adjust_glucose_concentration():
        V_guess = MEE.V
        MEE.V = flx.IQ_interpolation(
            brix_objective, 0., 1., x=V_guess, ytol=1e-5
        )
    
    # Mix in flocculant
    T1 = units.MixTank('T1', (MEE-0, lime, H3PO4, polymer))
    T1.tau = 0.10
    
    @T1.add_specification(run=True)
    def correct_flows():
        F_mass = T1.ins[0].F_mass
        # correct lime and phosphoric acid
        lime.imass['CaO', 'Water'] = 0.1 * F_mass * np.array([0.046, 0.954])
        H3PO4.imass['H3PO4', 'Water'] = 0.025 * F_mass
    
    # Separate residual solids
    C1 = units.Clarifier('C1', T1-0, 
                           split=dict(Ash=0,
                                      Cellulose=0,
                                      Flocculant=1,
                                      Glucose=1,
                                      Hemicellulose=0,
                                      Lignin=0,
                                      CaO=1,
                                      H3PO4=1,
                                      Sucrose=1,
                                      Water=0.99))
    P2 = units.Pump('P2', C1-0, P=101325)
    M1 = units.Mixer('M1', (P2-0, '', ''))
    E1 = units.Flash('E1', M1-0, V=0.5, P=15000)
    
    def get_purity(flash):
        effluent = flash.outs[1]
        return effluent.imass['Sugar'] / effluent.F_mass
    
    def purity_objective(V, flash):
        flash.V = V
        flash._run()
        return flash.purity - get_purity(flash)
    
    def adjust_purity(flash):
        V_guess = flash.V
        y0 = purity_objective(0., flash)
        if y0 < 0.: return
        y1 = purity_objective(1., flash)
        if y1 > 0.: return
        flash.V = flx.IQ_interpolation(
            purity_objective, 0., 1., y0, y1, x=V_guess, ytol=1e-5,
            args=(flash,),
        )
    
    E1.add_specification(adjust_purity, run=False, args=(E1,))
    E1.purity = 0.8623
    BC1 = units.BatchCrystallizer('BC1', E1-1, tau=8, V=3785, T=55 + 273.15)
    
    def get_split(molasses_flow, molasses_purity, crystal_flow, crystal_purity):
        s_crystal = crystal_flow * crystal_purity
        s_molasses = molasses_flow * molasses_purity
        s_split = s_crystal / (s_crystal + s_molasses)
        o_crystal = crystal_flow * (100 - crystal_purity)
        o_molasses = molasses_flow * (100 - molasses_purity)
        o_split = o_crystal / (o_crystal + o_molasses)
        return dict(
            Water=o_split,
            H3PO4=o_split,
            CaO=o_split,
            Sugar=s_split,
        )
    
    C2 = units.SolidsCentrifuge('C2', 
        BC1-0, 
        split=get_split(19.53, 62.91, 29.68, 98.61),
        moisture_content=None,
    )
    
    units.StorageTank('S1', C2-0, sugar, tau=27 * 7)
    
    def correct_wash_water(mixer):
        mixer.ins[1].imass['Water'] = mixer.ins[0].imass['Sugar']
    
    M2 = units.Mixer('M2', (C2-1, ''))
    M2.add_specification(correct_wash_water, run=True, args=(M2,))
    P3 = units.Pump('P3', M2-0, P=101325)
    E2 = units.Flash('E2', P3-0, V=0.5, P=15000)
    E2.add_specification(adjust_purity, run=False, args=(E2,))
    E2.purity = 0.6291
    BC2 = units.BatchCrystallizer('BC2', E2-1, tau=24, V=3785, T=50 + 273.15)
    C3 = units.SolidsCentrifuge('C3', 
        BC2-0, (2-M1, ''),
        split=get_split(4.34, 33.88, 3.15, 96.49),
        moisture_content=None,
    )
    M3 = units.Mixer('M3', (C3-1, ''))
    M3.add_specification(correct_wash_water, run=True, args=(M3,))
    P4 = units.Pump('P4', M3-0, P=101325)
    E3 = units.Flash('E3', P4-0, V=0.5, P=15000)
    E3.add_specification(adjust_purity, run=False, args=(E3,))
    E3.purity = 0.5450
    BC3 = units.BatchCrystallizer('BC3', E3-1, tau=40, V=3785, T=45 + 273.15)
    C4 = units.SolidsCentrifuge('C4', 
        BC3-0, (1-M1, ''),
        split=get_split(9.04, 32.88, 4.48, 93.84),
        moisture_content=None,
    )
    units.StorageTank('S2', C4-1, molasses, tau=24 * 7)

@SystemFactory(
    ID='sucrose_fermentation_sys',
    ins=[screened_juice],
    outs=[beer, evaporator_condensate, vent],
)
def create_sucrose_fermentation_system(ins, outs,
        scrubber=None, product_group=None, Fermentor=None, titer=None,
        productivity=None, ignored_volume=None, fermentation_reaction=None,
        fed_batch=None,
    ):
    screened_juice, = ins
    beer, evaporator_condensate, vent = outs
    if ignored_volume is None: ignored_volume = 'Lipid'
    if titer is None: titer = 117.0056 # g / L
    if productivity is None: productivity = 13
    if product_group is None: product_group = 'Ethanol'
    if Fermentor is None: Fermentor = units.Fermentation
    if fed_batch is None: fed_batch = False
    if Fermentor._N_outs == 2:
        if scrubber is None: scrubber = True
        fermentor_outs = [('CO2' if scrubber else vent), '']
    else:
        scrubber = False
        fermentor_outs = ['']
    dilution_water = bst.Stream('dilution_water')
    
    if fed_batch:
        if 'Sugar' not in dilution_water.chemicals:
            dilution_water.chemicals.define_group('Sugar', ('Glucose', 'Sucrose', 'Xylose'))
        
        SX0 = bst.Splitter(300, screened_juice, split=0.2)
        F301 = units.MultiEffectEvaporator('F301',
                                           SX0-1,
                                           P=(101325, 69682, 47057, 30953, 19781),
                                           V_definition='First-effect',
                                           thermo=dilution_water.thermo.ideal(),
                                           V=0.3) # fraction evaporated
        F301.brix = 95
        def get_brix():
            effluent = F301.outs[0]
            water = effluent.imass['Water']
            if water < 0.0001: water = 0.0001
            return 100 * effluent.imass['Sugar'] / water
        
        def brix_objective(V):
            F301.V = V
            F301._run()
            return F301.brix - get_brix()
        
        @F301.add_specification(run=False)
        def adjust_glucose_concentration():
            V_guess = F301.V
            F301.V = flx.IQ_interpolation(
                brix_objective, 0., 1., x=V_guess, ytol=1e-5
            )
        MT1 = bst.MixTank(300, F301-0)
        SX1 = bst.Splitter(300, ins=F301-1, outs=[evaporator_condensate, ''], split=0.9)
        P306 = units.Pump('P306', SX1-1, P=101325.)
        
        @SX1.add_specification(run=False)
        def sugar_concentration_adjustment():
            dilution_water = M301.ins[1]
            sugar_path = F301.path_until(R301, inclusive=False)[1:]
            for i in sugar_path: i.run()
            path = SX1.path_until(R301, inclusive=True)
            beer = R301.outs[1]
            target_titer = R301.titer
            def f(removed_water_split):
                SX1.split[:] = removed_water_split
                for unit in path: unit.run()
                return target_titer - get_titer()
            dilution_water.imass['Water'] = 0.
            x0 = 0
            x1 = 0.99
            y0 = f(x0)
            if y0 < 0.:
                product = float(beer.imass[product_group])
                current_titer = get_titer()
                ignored_product = P306.outs[0].imass[product_group]
                required_water = (1./target_titer - 1./current_titer) * (product - ignored_product) * 1000.
                dilution_water.imass['Water'] = max(required_water, 0)
            else:
                y1 = f(x1)
                if y1 > 0.:
                    long_path = [SX0, F301, *sugar_path]
                    for split in (0.20, 0.15, 0.10, 0.5, 0.):
                        SX0.split[:] = split
                        for i in long_path: i.run()
                        y1 = f(x1)
                        if y1 < 0.: break
                SX1.split[:] = flx.IQ_interpolation(f, x0, x1, y0, y1, x=SX1.split[0], ytol=1e-5, xtol=1e-6)
            R301.tau = target_titer / R301.productivity 
            SX0.split[:] = 0.2 # Restart
    else:
        F301 = units.MultiEffectEvaporator('F301',
                                           screened_juice,
                                           P=(101325, 69682, 47057, 30953, 19781),
                                           outs=('', evaporator_condensate),
                                           V_definition='First-effect',
                                           thermo=dilution_water.thermo.ideal(),
                                           V=0.3) # fraction evaporated
        P306 = units.Pump('P306', F301-0)
        # Note: value of steam ~ 6.86 for the following 
        # (101325, 73580.467, 50891.17, 32777.406, 19999.925, 11331.5),
        
        def titer_at_fraction_evaporated_objective(V, path):
            F301.V = V
            for i in path: i._run()
            return R301.titer - get_titer()
    
        F301.P_original = tuple(F301.P)
        @F301.add_specification(run=False)
        def evaporation():
            V_guess = F301.V
            s_dilution_water = M301.ins[-1]
            s_dilution_water.empty()
            path = F301.path_until(R301, inclusive=True)
            F301.V = 0
            for i in path: i._run()
            dilution_water = get_dilution_water()
            F301.P = F301.P_original
            F301._reload_components = True
            if dilution_water < 0.:
                x0 = 0.
                y0 = titer_at_fraction_evaporated_objective(x0, path)
                x1 = 0.95
                y1 = titer_at_fraction_evaporated_objective(x1, path)
                if y1 > 0.: raise RuntimeError('cannot evaporate to target sugar concentration')
                if y0 < 0.:
                    x0 = 0.
                    x1 = 0.1
                    F301.P = list(F301.P_original)
                    for i in range(F301._N_evap-1):
                        if f(1e-6) < 0.:
                            F301.P.pop()
                            F301._reload_components = True
                        else:
                            break  
                F301.V = flx.IQ_interpolation(
                    titer_at_fraction_evaporated_objective,
                    x0, x1, y0, y1, x=V_guess, ytol=1e-5,
                    args=(path,),
                )
            else:
                s_dilution_water.imass['Water'] = dilution_water
            R301.tau = R301.titer / R301.productivity
    
    # Mix sugar solutions
    M301 = units.Mixer('M301', ins=(P306-0, dilution_water))
    if fed_batch: M301.ins.append(SX0-0)
    
    # Cool for fermentation
    H301 = units.HXutility('H301', M301-0, T=295.15)
    
    # Ethanol Production
    R301 = Fermentor('R301', 
        ins=[H301-0, ''],
        outs=fermentor_outs, 
        tau=9, efficiency=0.90, N=4,
        fermentation_reaction=fermentation_reaction,
    ) 
    if fed_batch: R301.ins.append(MT1.outs[0])
    T301 = units.StorageTank('T301', R301-1, tau=4, vessel_material='Carbon steel')
    T301.line = 'Beer tank'
    
    # Separate 99% of yeast
    C301 = units.SolidsCentrifuge('C301', 
                                  ins=T301-0,
                                  outs=('', '' if scrubber else beer),
                                  split=(1-1e-6, 0.99, 1, 0.01),
                                  order=('Ethanol', 'Glucose', 'H3PO4', 'DryYeast'),
                                  solids=('DryYeast',))
    C301.split[:] = 1. - C301.split
    if 'Lipid' in C301.chemicals: C301.isplit['Lipid'] = 0.
    
    S302 = units.FakeSplitter('S302', C301-0, (1-R301, 'Yeast'))
    R301.titer = titer # g / L
    R301.productivity = productivity # g / L-h
    def adjust_yeast_recycle():
        recycle, beer = S302.outs
        feed, = S302.ins
        yeast = 0.1 * feed.F_mass
        m = C301.moisture_content
        recycle.imass['Yeast', 'Water'] = [yeast, yeast * m / (1 - m)]
        beer.copy_like(feed)
        beer.separate_out(recycle, energy_balance=False)
        beer.mol[beer.mol < 0.] = 0.
        beer.T = recycle.T = feed.T
    S302.specification = adjust_yeast_recycle
    
    def get_titer():
        s = R301.outs[1]
        ignored = s.ivol[ignored_volume] if ignored_volume in s.chemicals else 0.
        ignored_product = sum([i.imass[product_group] for i in R301.ins])
        return (s.imass[product_group] - ignored_product) / (s.F_vol - ignored)
    R301.get_titer = get_titer
    
    def get_dilution_water():
        target = R301.titer
        current = get_titer()
        ignored_product = sum([i.imass[product_group] for i in R301.ins])
        return (1./target - 1./current) * (R301.outs[1].imass[product_group] - ignored_product) * 1000.
    
    if scrubber:
        stripping_water = bst.Stream('stripping_water',
                                     Water=26836,
                                     units='kg/hr')
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D301.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
        
        D301 = units.VentScrubber('D301', ins=(stripping_water, R301-0), 
                                  outs=(vent, ''),
                                  gas=('CO2', 'O2'))
        units.Mixer('M302', ins=(C301-1, D301-1), outs=beer)
    

@SystemFactory(
    ID='sucrose_to_ethanol_sys',
    ins=[screened_juice, denaturant],
    outs=[ethanol, stillage, stripper_bottoms_product, evaporator_condensate]
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
    ins=[sugarcane, H3PO4, lime, polymer, denaturant], 
    outs=[ethanol, vinasse, wastewater, emissions, ash_disposal]
)
def create_sugarcane_to_ethanol_system(ins, outs, 
                                       use_area_convention=False,
                                       pellet_bagasse=None):
    s = f.stream
    u = f.unit
    
    sugarcane, H3PO4, lime, polymer, denaturant = ins
    ethanol, vinasse, wastewater, emissions, ash_disposal = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        area=100 if use_area_convention else None,
        ins=[sugarcane],
        outs=[''],
        mockup=True,
    )
    juicing_sys = create_juicing_system_with_fiber_screener(
        area=200 if use_area_convention else None,
        ins=[feedstock_handling_sys-0, H3PO4, lime, polymer],
        pellet_bagasse=pellet_bagasse,
        mockup=True
    )
    ethanol_production_sys, edct = create_sucrose_to_ethanol_system(
        area=300 if use_area_convention else None,
        udct=True,
        ins=(juicing_sys-0, denaturant), outs=(ethanol, vinasse),
        mockup=True
    )
    M305 = units.Mixer(400 if use_area_convention else 'M305', 
        ins=(juicing_sys-2, *ethanol_production_sys-[2, 3]),
        outs=wastewater
    )
    
    ### Facilities ###    
    
    BT = units.BoilerTurbogenerator(400 if use_area_convention else 'BT',
        (juicing_sys-1, '', 'boiler_makeup_water', 'natural_gas', '', ''),
        outs=(emissions, 'rejected_water_and_blowdown', ash_disposal),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85
    )
    CT = units.CoolingTower(500 if use_area_convention else 'CT')
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    process_water_streams = (s.imbibition_water,
                             s.rvf_wash_water,
                             s.stripping_water,
                             *makeup_water_streams)
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    CWP = units.ChilledWaterPackage(500 if use_area_convention else 'CWP')
    PWC = units.ProcessWaterCenter(500 if use_area_convention else 'PWC',
                                   (bst.Stream(), makeup_water),
                                   (),
                                   None,
                                   makeup_water_streams,
                                   process_water_streams)
    
    F301 = edct['F301']
    D303 = edct['D303']
    HXN = bst.HeatExchangerNetwork(600 if use_area_convention else 'HXN',
                                   units=[F301, D303.condenser])
    
    # if vinasse_to_wastewater:
    #     plant_air = bst.Stream('plant_air', N2=83333, units='kg/hr')
    #     ADP = bst.facilities.AirDistributionPackage('ADP', plant_air)
    #     @ADP.add_specification(run=True)
    #     def adjust_plant_air():
    #         plant_air.imass['N2'] = 0.8 * feedstock_handling_sys.ins[0].F_mass
            
    #     wastewater_treatment_sys = bst.create_wastewater_treatment_system(
    #         ins=[vinasse],
    #         mockup=True,
    #     )
    
@SystemFactory(
    ID='sugarcane_sys', 
    ins=[sugarcane, H3PO4, lime, polymer], 
    outs=[sugar, molasses, wastewater, emissions, ash_disposal]
)
def create_sugarcane_to_sugar_and_molasses_system(ins, outs, 
                                       use_area_convention=False,
                                       pellet_bagasse=None):
    s = f.stream
    u = f.unit
    
    sugarcane, H3PO4, lime, polymer, denaturant = ins
    sugar, molasses, wastewater, emissions, ash_disposal = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        area=100 if use_area_convention else None,
        ins=[sugarcane],
        outs=[''],
        mockup=True,
    )
    juicing_sys = create_juicing_system_with_fiber_screener(
        area=200 if use_area_convention else None,
        ins=[feedstock_handling_sys-0, H3PO4, lime, polymer],
        pellet_bagasse=pellet_bagasse,
        mockup=True
    )
    
    sugar_crystallization_sys, edct = create_sugar_crystallization_system(
        area=300 if use_area_convention else None,
        udct=True,
        ins=juicing_sys-0, outs=(sugar, molasses),
        mockup=True
    )
    M305 = units.Mixer(400 if use_area_convention else 'M305', 
        ins=(juicing_sys-2,),
        outs=wastewater
    )
    
    ### Facilities ###    
    
    BT = units.BoilerTurbogenerator(400 if use_area_convention else 'BT',
        (juicing_sys-1, '', 'boiler_makeup_water', 'natural_gas', '', ''),
        outs=(emissions, 'rejected_water_and_blowdown', ash_disposal),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85
    )
    CT = units.CoolingTower(500 if use_area_convention else 'CT')
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    process_water_streams = (s.imbibition_water,
                             s.rvf_wash_water,
                             *makeup_water_streams)
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    CWP = units.ChilledWaterPackage(500 if use_area_convention else 'CWP')
    PWC = units.ProcessWaterCenter(500 if use_area_convention else 'PWC',
                                   (bst.Stream(), makeup_water),
                                   (),
                                   None,
                                   makeup_water_streams,
                                   process_water_streams)
    HXN = bst.HeatExchangerNetwork(600 if use_area_convention else 'HXN')
    
    # if vinasse_to_wastewater:
    #     plant_air = bst.Stream('plant_air', N2=83333, units='kg/hr')
    #     ADP = bst.facilities.AirDistributionPackage('ADP', plant_air)
    #     @ADP.add_specification(run=True)
    #     def adjust_plant_air():
    #         plant_air.imass['N2'] = 0.8 * feedstock_handling_sys.ins[0].F_mass
            
    #     wastewater_treatment_sys = bst.create_wastewater_treatment_system(
    #         ins=[vinasse],
    #         mockup=True,
    #     )