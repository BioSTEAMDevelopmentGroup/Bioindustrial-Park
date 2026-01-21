# -*- coding: utf-8 -*-
"""
LegHb Production System - Config 1 (Modular Refactored Version)

Modular refactoring of LegHb production system into standard BioSTEAM Process Areas.

Process Areas:
    - Area 200: Media Preparation
    - Area 300: Conversion (Fermentation)
    - Area 400: Recovery (Harvest & Cell Disruption)
    - Area 500: Purification (UF/DF Membrane Separation)
    - Area 600: Formulation (Thermal Treatment & Product)
    - Area 900: Facilities & Utilities

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biorefineries.prefers.v1._process_settings import set_GWPCF, GWP_CFs, set_GWPCF_Multi, load_process_settings
from biorefineries.prefers.v1.LegHb import _streams as s
import biosteam as bst
from thermosteam import Stream
from biosteam import F
import thermosteam as tmo
import numpy as np
from biorefineries.prefers.v1.LegHb import _chemicals as c
from biorefineries.prefers.v1 import _units as u
from biorefineries.prefers.v1._process_settings import price

# =============================================================================
# MODULE INITIALIZATION
# =============================================================================
LEGHB_THERMO = c.create_chemicals_LegHb()
bst.settings.set_thermo(LEGHB_THERMO, skip_checks=True)
bst.preferences.classic_mode()
SHOW_RXN = False

__all__ = (
    'create_LegHb_system',
    'set_production_rate',
    'check_LegHb_specifications',
    # Area creation functions
    'create_area_200_media_prep',
    'create_area_300_conversion',
    'create_area_400_recovery',
    'create_area_500_purification',
    'create_area_600_formulation',
    'create_area_900_facilities',
    # Reaction setup
    'create_fermentation_reactions',
)

# =============================================================================
# FERMENTATION PARAMETERS (Shared across areas)
# =============================================================================
def get_fermentation_parameters():
    """
    Return fermentation parameters validated against literature.
    
    References:
        Yao et al. 2025, World J Microbiol Biotechnol 41:404
        Tian et al. 2024b - LegHb titer up to 7.27 g/L in K. marxianus
    
    Returns
    -------
    dict
        Dictionary of fermentation parameters
    """
    theta_O2 = 0.5  # Dissolved oxygen concentration [% saturation]
    agitation_power = 0.985  # [kW/m³] - Ref: Typical for STR bioreactors
    design = 'Stirred tank'  # Reactor type
    method = "Riet"  # kLa correlation method
    T_operation = 273.15 + 32  # [K] - Ref: 28-32°C optimal
    Q_O2_consumption = -110 * 4184  # [kJ/kmol] - Heat of aerobic metabolism
    dT_hx_loop = 8  # [°C] - Heat exchanger approach temperature
    cooler_pressure_drop = 20684  # [Pa]
    compressor_isentropic_efficiency = 0.85
    V_max = 500  # [m³] - Max vessel volume
    
    # LegHb production parameters
    titer_LegHb = 7.27  # [g/L] - Literature validated
    productivity_LegHb = 7.27 / 72  # [g/L/h] - 72h fermentation cycle
    Y_p = 7.27 * 4 / 1300  # [by wt] - ~2.2% product yield on glucose
    Y_b = 0.43  # [by wt] - Biomass yield, typical for yeast
    yield_LegHb = Y_p
    
    return {
        'theta_O2': theta_O2,
        'agitation_power': agitation_power,
        'design': design,
        'method': method,
        'T_operation': T_operation,
        'Q_O2_consumption': Q_O2_consumption,
        'dT_hx_loop': dT_hx_loop,
        'cooler_pressure_drop': cooler_pressure_drop,
        'compressor_isentropic_efficiency': compressor_isentropic_efficiency,
        'V_max': V_max,
        'titer_LegHb': titer_LegHb,
        'productivity_LegHb': productivity_LegHb,
        'Y_p': Y_p,
        'Y_b': Y_b,
        'yield_LegHb': yield_LegHb,
    }

# =============================================================================
# REACTION DEFINITIONS
# =============================================================================
def create_fermentation_reactions(params=None):
    """
    Create reaction systems for LegHb fermentation.
    
    Parameters
    ----------
    params : dict, optional
        Fermentation parameters. Uses defaults if not provided.
    
    Returns
    -------
    dict
        Dictionary containing:
        - fermentation_reaction: PRxn for LegHb production
        - neutralization_reaction: Rxn for pH control
        - cell_growth_reaction: Rxn for biomass formation
        - respiration_reaction1: Rxn for seed train respiration
        - respiration_reaction2: Rxn for main fermenter respiration
        - RXN: Complete ReactionSystem
    """
    if params is None:
        params = get_fermentation_parameters()
    
    yield_LegHb = params['yield_LegHb']
    Y_b = params['Y_b']
    
    # Fermentation reactions
    fermentation_reaction = bst.PRxn([
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                reactant='Glucose', X=yield_LegHb*0.01, check_atomic_balance=True),
        bst.Rxn('Glucose + 0.01646 (NH4)2SO4 + 1.61317 NH3 -> 6 Globin_In + 0.28807 O2 + 3.68724 H2O',
                reactant='Glucose', X=yield_LegHb*0.04, check_atomic_balance=True),
        bst.Rxn('Glucose + 0.00786 FeSO4 + 0.00786 (NH4)2SO4 + 1.58847 NH3 -> 6 Leghemoglobin_In + 0.30275 O2  + 3.70380 H2O',
                reactant='Glucose', X=yield_LegHb, check_atomic_balance=True),
    ])
    fermentation_reaction[2].product_yield('Leghemoglobin_In', basis='wt', product_yield=yield_LegHb)
    
    neutralization_reaction = bst.Rxn(
        'H2SO4 + 2 NH3 -> (NH4)2SO4', reactant='NH3', X=1,
        check_atomic_balance=True
    )
    
    cell_growth_reaction = bst.Rxn(
        'Glucose + 0.8364 NH3 + 0.0108 (NH4)2SO4 -> 2.01 H2O + 0.106 O2 + 6 Pichia_pastoris', 
        'Glucose', X=(1-yield_LegHb*1.1)*Y_b,
        correct_atomic_balance=True
    )
    cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=(1-yield_LegHb*1.1)*Y_b)
    
    respiration_reaction1 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', 1 - Y_b,
        check_atomic_balance=True
    )
    
    respiration_reaction2 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', 
        1 - cell_growth_reaction.X - fermentation_reaction[2].X * 1.1,
        check_atomic_balance=True
    )
    
    bst.settings.chemicals.set_alias('Pichia_pastoris', 'cellmass')
    
    RXN = bst.ReactionSystem(
        fermentation_reaction,
        bst.PRxn([cell_growth_reaction, respiration_reaction2])
    )
    
    if SHOW_RXN:
        RXN.show()
    
    return {
        'fermentation_reaction': fermentation_reaction,
        'neutralization_reaction': neutralization_reaction,
        'cell_growth_reaction': cell_growth_reaction,
        'respiration_reaction1': respiration_reaction1,
        'respiration_reaction2': respiration_reaction2,
        'RXN': RXN,
    }

# =============================================================================
# AREA 200: MEDIA PREPARATION
# =============================================================================
def create_area_200_media_prep(SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt):
    """
    Create Area 200: Media Preparation
    
    Purpose: Prepare all feed solutions for fermentation
    
    Units:
        - M301: Seed solution 1 preparation
        - M302: Seed solution 2 + culture media preparation
        - M303: Seed hold tank (combines both seed solutions)
        - M304: Glucose solution preparation (50% dilution)
        - T301: Glucose storage tank (96 hours residence)
        - T302: Ammonia storage tank (25 wt% aqueous ammonia)
    
    Parameters
    ----------
    SeedIn1 : Stream
        Seed input stream 1
    SeedIn2 : Stream
        Seed input stream 2
    CultureIn : Stream
        Culture media input stream
    Glucose : Stream
        Glucose feed stream
    NH3_25wt : Stream
        25 wt% aqueous ammonia stream
    
    Returns
    -------
    dict
        Dictionary of area units and output streams
    """
    # Seed solution 1 preparation
    M301 = bst.MixTank('M301', ins=[SeedIn1, 'Water1'], outs='M301Out', tau=16)
    
    @M301.add_specification(run=True)
    def update_seed1_inputs():
        target_stream = bst.Stream(**{**s.SeedSolution1, 'ID': None})
        SeedIn1.imass['Seed'] = target_stream.imass['Seed']
        M301.ins[1].imass['H2O'] = target_stream.imass['H2O']
        M301.ins[1].T = 25 + 273.15
    
    # Seed solution 2 + culture media preparation
    M302 = bst.MixTank('M302', ins=[SeedIn2, CultureIn, 'Water2'], outs='M302Out', tau=16)
    
    @M302.add_specification(run=True)
    def update_culture_inputs():
        target_stream = bst.Stream(**{**s.SeedSolution2, 'ID': None})
        SeedIn2.imass['Seed'] = target_stream.imass['Seed']
        M302.ins[2].imass['H2O'] = target_stream.imass['H2O']
        M302.ins[2].T = 25 + 273.15
        CultureIn.imass['Culture'] = target_stream.imass['SeedSolution'] * (0.1 + 60 + 0.15191) / 1000
    
    # Seed hold tank
    M303 = u.SeedHoldTank('M303', ins=[M301-0, M302-0], outs='M303Out')
    
    # Glucose solution preparation (50% dilution)
    M304 = bst.MixTank('M304', ins=[Glucose, 'Water3'], outs='M304Out', tau=16)
    
    @M304.add_specification(run=True)
    def update_water_content():
        M304.ins[1].imass['H2O'] = Glucose.imass['Glucose'] / 2
        M304.ins[1].T = 25 + 273.15
    
    # Glucose storage tank
    T301 = bst.StorageTank('T301', ins=M304-0, outs='T301Out', tau=16*4+72)
    
    # Ammonia storage tank
    T302 = u.AmmoniaStorageTank('T302', ins=NH3_25wt, outs='T302Out')
    
    return {
        'units': [M301, M302, M303, M304, T301, T302],
        'seed_out': M303.outs[0],
        'glucose_out': T301.outs[0],
        'ammonia_out': T302.outs[0],
        'M301': M301, 'M302': M302, 'M303': M303,
        'M304': M304, 'T301': T301, 'T302': T302,
    }

# =============================================================================
# AREA 300: CONVERSION (FERMENTATION)
# =============================================================================
def create_area_300_conversion(seed_in, glucose_in, ammonia_in, vent1, vent2, reactions, params=None):
    """
    Create Area 300: Conversion (Fermentation)
    
    Purpose: Convert glucose to LegHb through microbial fermentation
    
    Process Flow:
        1. Seed Train (R301): Grow inoculum culture
        2. Main Fermentation (R302): Batch aerobic fermentation (72h cycle)
    
    Parameters
    ----------
    seed_in : Stream
        Combined seed solution from Area 200
    glucose_in : Stream
        Diluted glucose solution from Area 200
    ammonia_in : Stream
        Ammonia solution from Area 200
    vent1 : Stream
        CO2 exhaust from seed train
    vent2 : Stream
        CO2 exhaust from main fermentation
    reactions : dict
        Reaction objects from create_fermentation_reactions()
    params : dict, optional
        Fermentation parameters
    
    Returns
    -------
    dict
        Dictionary of area units and output streams
    """
    if params is None:
        params = get_fermentation_parameters()
    
    fermentation_reaction = reactions['fermentation_reaction']
    cell_growth_reaction = reactions['cell_growth_reaction']
    respiration_reaction1 = reactions['respiration_reaction1']
    respiration_reaction2 = reactions['respiration_reaction2']
    neutralization_reaction = reactions['neutralization_reaction']
    RXN = reactions['RXN']
    
    # Seed train bioreactor
    R301 = u.SeedTrain(
        'R301',
        ins=[seed_in],
        outs=[vent1, 'R301Out'],
        reactions=bst.PRxn([cell_growth_reaction, respiration_reaction1]),
        saccharification=None,
        T=32 + 273.15,
    )
    R301.add_specification(run=True)
    
    # Main aerobic fermentation vessel
    R302 = u.AeratedFermentation(
        'R302',
        ins=[R301-1, glucose_in, ammonia_in, bst.Stream('FilteredAir', phase='g', P=2 * 101325)],
        outs=[vent2, 'Broth'],
        fermentation_reaction=fermentation_reaction,
        cell_growth_reaction=cell_growth_reaction,
        respiration_reaction=respiration_reaction2,
        neutralization_reaction=neutralization_reaction,
        design=params['design'],
        method=params['method'],
        theta_O2=params['theta_O2'],
        V_max=params['V_max'],
        Q_O2_consumption=params['Q_O2_consumption'],
        dT_hx_loop=params['dT_hx_loop'],
        T=params['T_operation'],
        batch=True,
        reactions=RXN,
        kW_per_m3=params['agitation_power'],
        tau=params['titer_LegHb'] / params['productivity_LegHb'],
        cooler_pressure_drop=params['cooler_pressure_drop'],
        compressor_isentropic_efficiency=params['compressor_isentropic_efficiency'],
        P=1 * 101325,
    )
    R302.target_titer = params['titer_LegHb']
    R302.target_productivity = params['productivity_LegHb']
    R302.target_yield = params['yield_LegHb']
    
    @R302.add_specification(run=True)
    def update_reaction_time_and_yield():
        R302.tau = R302.target_titer / R302.target_productivity
        fermentation_reaction[2].product_yield('Leghemoglobin_In', basis='wt', product_yield=R302.target_yield)
    
    return {
        'units': [R301, R302],
        'broth_out': R302.outs[1],
        'R301': R301, 'R302': R302,
    }

# =============================================================================
# AREA 400: RECOVERY (BIOMASS HARVEST & CELL DISRUPTION)
# =============================================================================
def create_area_400_recovery(broth_in, DfUltraBuffer_wash):
    """
    Create Area 400: Recovery (Biomass Harvest & Cell Disruption)
    
    Purpose: Release and clarify intracellular LegHb from yeast cells
    
    Process Flow:
        1. Biomass Harvest (C401): Centrifuge to concentrate cells
        2. Cell Washing (M401, H401_wash, M402_wash, C402): Remove spent media
        3. Cell Disruption (S401): High-pressure homogenization
        4. Debris Removal (H401, S402): Cool and centrifuge
        5. Lysate Clarification (S403, S404): Depth filtration + screw press
    
    Parameters
    ----------
    broth_in : Stream
        Fermentation broth from Area 300
    DfUltraBuffer_wash : Stream
        Buffer for cell washing
    
    Returns
    -------
    dict
        Dictionary of area units and output streams
    """
    # Primary harvest centrifuge
    C401 = bst.SolidsCentrifuge(
        'C401',
        ins=broth_in,
        outs=('CellCream', 'SpentMedia'),
        split={
            'cellmass': 0.98,
            'Leghemoglobin_In': 0.98,
            'Globin_In': 0.98,
            'Heme_b': 0.85,
            'Glucan': 0.95,
            'Mannoprotein': 0.95,
            'Chitin': 0.95,
            'OleicAcid': 0.90,
            'RNA': 0.90,
        },
        moisture_content=0.55,
    )
    
    # Wash buffer preparation
    M401 = bst.MixTank('M401', ins=(DfUltraBuffer_wash, 'WashWater'), 
                       outs='WashBufferOut', tau=0.5)
    
    @M401.add_specification(run=True)
    def update_wash_buffer():
        M401.ins[1].imass['H2O'] = C401.outs[0].F_mass * 1.5
        M401.ins[0].imol['DfUltraBuffer'] = M401.ins[1].imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000
    
    # Cool wash buffer
    H401_wash = bst.HXutility(
        'H401_wash',
        ins=M401-0,
        outs='ColdWashBuffer',
        T=10 + 273.15,
        cool_only=True,
    )
    
    # Cell wash mixer
    M402_wash = bst.MixTank('M402_wash', ins=(C401-0, H401_wash-0), 
                            outs='WashedCellSlurry', tau=0.25)
    
    # Second centrifuge for washed cells
    C402 = bst.SolidsCentrifuge(
        'C402',
        ins=M402_wash-0,
        outs=('WashedCellCream', 'WashEffluent'),
        split={
            'cellmass': 0.98,
            'Leghemoglobin_In': 0.98,
            'Globin_In': 0.98,
            'Heme_b': 0.85,
            'Glucan': 0.95,
            'Mannoprotein': 0.95,
        },
        moisture_content=0.55,
    )
    
    # Cell disruption (high-pressure homogenization)
    S401 = u.CellDisruption(
        'S401',
        ins=C402-0,
        outs='CrudeHomogenate',
        P_high=1000e5,  # 1000 bar
        cell_disruption_efficiency=0.90,
    )
    
    # Post-valve cooling
    H401 = bst.HXutility(
        'H401',
        ins=S401-0,
        outs='CooledHomogenate',
        T=15 + 273.15,
        cool_only=True,
    )
    
    # High-speed centrifuge for debris removal
    S402 = bst.SolidsCentrifuge(
        'S402',
        ins=H401-0,
        outs=('CellDebris', 'CrudeLysate'),
        split={
            'cellmass': 0.98,
            'Glucan': 0.98,
            'Chitin': 0.98,
            'OleicAcid': 0.98,
            'RNA': 0.85,
            'Leghemoglobin': 0.02,
            'Globin': 0.05,
            'Mannoprotein': 0.98,
        },
        moisture_content=0.20,
    )
    
    # Depth filtration
    S403 = u.Filtration.from_preset(
        'MF',
        'S403',
        ins=S402-1,
        outs=('FilterCake', 'ClarifiedLysate'),
        solid_capture_efficiency=0.95,
        cake_moisture_content=0.30,
    )
    S403.add_specification(run=True)
    
    # Screw press for debris dewatering
    S404 = bst.ScrewPress(
        'S404',
        ins=(S402-0, S403-0),
        outs=('DehydratedDebris', 'PressLiquor'),
        split=0.999,
        moisture_content=0.20
    )
    
    return {
        'units': [C401, M401, H401_wash, M402_wash, C402, S401, H401, S402, S403, S404],
        'clarified_lysate': S403.outs[1],
        'spent_media': C401.outs[1],
        'wash_effluent': C402.outs[1],
        'dehydrated_debris': S404.outs[0],
        'press_liquor': S404.outs[1],
        'C401': C401, 'M401': M401, 'H401_wash': H401_wash,
        'M402_wash': M402_wash, 'C402': C402, 'S401': S401,
        'H401': H401, 'S402': S402, 'S403': S403, 'S404': S404,
    }

# =============================================================================
# AREA 500: PURIFICATION (UF/DF MEMBRANE SEPARATION)
# =============================================================================
def create_area_500_purification(clarified_lysate, DfUltraBuffer):
    """
    Create Area 500: Purification (UF/DF Membrane Separation)
    
    Purpose: Concentrate LegHb and remove impurities via tangential flow filtration
    
    Process Flow:
        1. DF Buffer Preparation (M403, H402): Cold buffer for diafiltration
        2. UF/DF Stage 1 (U401): Concentration + buffer exchange (5-7 diavolumes)
        3. UF/DF Stage 2 (U404): Final concentration
    
    Parameters
    ----------
    clarified_lysate : Stream
        Clarified lysate from Area 400
    DfUltraBuffer : Stream
        DF buffer stream
    
    Returns
    -------
    dict
        Dictionary of area units and output streams
    """
    # DF buffer preparation
    M403 = bst.MixTank('M403', ins=(DfUltraBuffer, 'DFBufferWater'), 
                       outs='DFBufferOut', tau=0.5)
    
    # Reference to S403 for specification (will be set via system context)
    clarified_lysate_stream = clarified_lysate
    
    @M403.add_specification(run=True)
    def update_df_buffer():
        feed_water = clarified_lysate_stream.imass['H2O']
        M403.ins[1].imass['H2O'] = feed_water * 6
        M403.ins[0].imol['DfUltraBuffer'] = M403.ins[1].imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000
    
    # Cool DF buffer
    H402 = bst.HXutility(
        'H402',
        ins=M403-0,
        outs='ColdDFBuffer',
        T=5 + 273.15,
        cool_only=True,
    )
    
    # UF/DF unit - concentration + buffer exchange
    U401 = u.Diafiltration.from_preset(
        'UF',
        'U401',
        ins=(clarified_lysate, H402-0),
        outs=('UFConcentrate', 'UFPermeate'),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
        TargetProduct_Retention=0.99,
        Salt_Retention=0.05,
    )
    U401.add_specification(run=True)
    
    # Final concentration
    U404 = u.Diafiltration(
        'U404',
        ins=(U401-0, bst.Stream('U404_buffer', H2O=0.001)),
        outs=('ConcentratedLegH', 'ConcentrationPermeate'),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
        TMP_bar1=3,
        FeedWater_Recovery_to_Permeate=0.85,
    )
    
    @U404.add_specification(run=True)
    def U404_adjust_water_recovery():
        U404.Salt_Retention = 0.95
        U404._run()
    
    return {
        'units': [M403, H402, U401, U404],
        'concentrated_product': U404.outs[0],
        'uf_permeate': U401.outs[1],
        'concentration_permeate': U404.outs[1],
        'M403': M403, 'H402': H402, 'U401': U401, 'U404': U404,
    }

# =============================================================================
# AREA 600: FORMULATION (THERMAL TREATMENT & FINAL PRODUCT)
# =============================================================================
def create_area_600_formulation(concentrated_product, LegHb_3):
    """
    Create Area 600: Formulation (Thermal Treatment & Final Product)
    
    Purpose: Ensure microbiological safety and achieve final product specification
    
    Process Flow:
        1. HTST Pasteurization (H403, T401): Heat to 72°C for 30 seconds
        2. Formulation (M404): Add antioxidant + dilution water
        3. Final Cooling (H406): Rapid chill to 4°C
    
    Parameters
    ----------
    concentrated_product : Stream
        Concentrated LegHb from Area 500
    LegHb_3 : Stream
        Final product output stream
    
    Returns
    -------
    dict
        Dictionary of area units and output streams
    """
    # Pasteurization heater
    H403 = bst.HXutility(
        'H403',
        ins=concentrated_product,
        outs='HeatedConcentrate',
        T=72 + 273.15,
        heat_only=True,
    )
    
    # Hold time tank
    T401 = bst.MixTank(
        'T401',
        ins=H403-0,
        outs='PasteurizedConcentrate',
        tau=30/3600  # 30 second hold time
    )
    
    # Antioxidant stream
    AntioxidantStream = bst.Stream(
        'AntioxidantStream',
        SodiumAscorbate=0.1,
        H2O=1.0,
        units='kg/hr',
        price=price.get('SodiumAscorbate', 5.0)
    )
    
    # Dilution water stream
    DilutionWater = bst.Stream(
        'DilutionWater',
        H2O=100.0,
        units='kg/hr',
    )
    
    # Formulation mixer
    M404 = bst.MixTank(
        'M404',
        ins=(T401-0, AntioxidantStream, DilutionWater),
        outs='FormulatedProduct',
        tau=0.1
    )
    M404.target_total_solids_percent = 20.0
    M404.target_legh_percent = 7.5
    
    @M404.add_specification(run=True)
    def update_formulation():
        feed = T401.outs[0]
        feed_total_mass = feed.F_mass
        
        if feed_total_mass <= 0:
            DilutionWater.imass['H2O'] = 0
            AntioxidantStream.imass['SodiumAscorbate'] = 0
            AntioxidantStream.imass['H2O'] = 0
            return
        
        feed_water_mass = feed.imass['H2O']
        feed_solids_mass = feed_total_mass - feed_water_mass
        feed_legh_mass = feed.imass['Leghemoglobin']
        
        target_solids_pct = getattr(M404, 'target_total_solids_percent', 20.0)
        target_legh_pct = getattr(M404, 'target_legh_percent', 7.5)
        
        if target_solids_pct > 0:
            required_total_mass = feed_solids_mass * 100.0 / target_solids_pct
            required_dilution = required_total_mass - feed_total_mass
            
            if feed_legh_mass > 0 and target_legh_pct > 0:
                required_total_mass_legh = feed_legh_mass * 100.0 / target_legh_pct
                required_dilution_legh = required_total_mass_legh - feed_total_mass
                required_dilution = max(required_dilution, required_dilution_legh)
            
            required_dilution = max(0, required_dilution)
        else:
            required_dilution = 0
        
        DilutionWater.imass['H2O'] = required_dilution
        final_product_mass = feed_total_mass + required_dilution
        AntioxidantStream.imass['SodiumAscorbate'] = final_product_mass * 0.001
        AntioxidantStream.imass['H2O'] = AntioxidantStream.imass['SodiumAscorbate'] * 9
    
    # Final cooling
    H406 = bst.HXutility(
        'H406',
        ins=M404-0,
        outs=LegHb_3,
        T=4 + 273.15,
        cool_only=True,
    )
    
    return {
        'units': [H403, T401, M404, H406],
        'product': H406.outs[0],
        'H403': H403, 'T401': T401, 'M404': M404, 'H406': H406,
        'AntioxidantStream': AntioxidantStream,
        'DilutionWater': DilutionWater,
    }

# =============================================================================
# AREA 900: FACILITIES & UTILITIES
# =============================================================================
def create_area_900_facilities(area_400, area_500, effluent1, use_area_convention=False):
    """
    Create Area 900: Facilities & Utilities
    
    Purpose: Provide utilities and manage wastewater
    
    Subsystems:
        A. Wastewater Treatment & Water Recovery
        B. Utility Systems (CT, CWP, PWC)
    
    Parameters
    ----------
    area_400 : dict
        Area 400 outputs containing waste streams
    area_500 : dict
        Area 500 outputs containing permeate streams
    effluent1 : Stream
        RO reject output stream
    use_area_convention : bool
        Whether to use area ID convention for utility units
    
    Returns
    -------
    dict
        Dictionary of area units and output streams
    """
    # Collect permeate streams for RO treatment
    M501 = bst.MixTank(
        'M501',
        ins=(
            area_400['press_liquor'],
            area_400['spent_media'],
            area_400['wash_effluent'],
            area_500['uf_permeate'],
            area_500['concentration_permeate'],
        ),
        outs='CombinedWastewater',
        tau=1
    )
    M501.add_specification(run=True)
    
    # RO treatment for water recovery
    S501 = u.ReverseOsmosis('S501', ins=M501-0, outs=('RO_treated_water1', effluent1))
    
    # Utility systems
    CT = bst.CoolingTower(500 if use_area_convention else 'CT')
    CWP = bst.ChilledWaterPackage(500 if use_area_convention else 'CWP')
    
    # Route dehydrated debris to disposal
    debris_disposal = bst.Stream('debris_disposal')
    M503_debris = bst.Mixer('M503_debris', ins=area_400['dehydrated_debris'], outs=debris_disposal)
    
    makeup_water_streams = (
        F.cooling_tower_makeup_water,
        F.Water1, F.Water2,
        F.Water3, F.WashWater,
        F.DFBufferWater
    )
    process_water_streams = (
        F.S501.outs[1],
        *makeup_water_streams
    )
    
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    
    PWC = bst.ProcessWaterCenter(
        500 if use_area_convention else 'PWC',
        ins=('recycled_RO_water', makeup_water, 'recycled_process_water', 'makeup_process_water'),
        outs=('RO_water', 'process_water', 'excess_water'),
        makeup_water_streams=makeup_water_streams,
        process_water_streams=process_water_streams,
        reverse_osmosis_water_price=0.000254,
        process_water_price=0.000135,
    )
    
    return {
        'units': [M501, S501, CT, CWP, M503_debris, PWC],
        'effluent': effluent1,
        'debris_disposal': debris_disposal,
        'M501': M501, 'S501': S501, 'CT': CT, 'CWP': CWP,
        'M503_debris': M503_debris, 'PWC': PWC,
    }

# =============================================================================
# MAIN SYSTEM FACTORY
# =============================================================================
@bst.SystemFactory(
    ID='LegHb_sys',
    ins=[s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt, s.DfUltraBuffer, s.DfUltraBuffer_wash],
    outs=[s.LegHb_3, s.vent1, s.vent2, s.effluent1],
    fthermo=lambda chemicals=None: LEGHB_THERMO,
)
def create_LegHb_system(ins, outs, use_area_convention=False):
    """
    Creates the LegHb (Leghemoglobin) production system.
    
    Modular architecture using discrete Process Area functions:
        - Area 200: Media Preparation
        - Area 300: Conversion (Fermentation)
        - Area 400: Recovery (Harvest & Cell Disruption)
        - Area 500: Purification (UF/DF)
        - Area 600: Formulation (Thermal & Product)
        - Area 900: Facilities & Utilities
    
    Parameters
    ----------
    ins : tuple
        Input streams (SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, DfUltraBuffer_wash)
    outs : tuple
        Output streams (LegHb_3, vent1, vent2, effluent1)
    use_area_convention : bool
        Whether to use area ID convention for utility units
    
    Returns
    -------
    tuple
        Output streams (LegHb_3, vent1, vent2, effluent1, DfUltraBuffer, DfUltraBuffer_wash)
    """
    bst.preferences.N = 50
    
    # Unpack streams
    SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, DfUltraBuffer_wash = ins
    LegHb_3, vent1, vent2, effluent1 = outs
    
    # Set GWP characterization factors
    set_GWPCF(Glucose, 'Glucose')
    set_GWPCF_Multi(SeedIn1, ['AmmoniumSulfate', 'Glucose', 'MagnesiumSulfate', 'KH2PO4'],
                   [0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(SeedIn2, ['AmmoniumSulfate', 'Glucose', 'MagnesiumSulfate', 'KH2PO4'],
                   [0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(CultureIn, ['Glycine', 'Glucose', 'IronSulfate'],
                   [0.1/(0.1+60+0.15191), 60/(0.1+60+0.15191), 0.15191/(0.1+60+0.15191)])
    set_GWPCF(NH3_25wt, 'Ammonia_SEA', dilution=0.25)
    set_GWPCF_Multi(DfUltraBuffer, ['KH2PO4', 'NaCl'], [0.8472, 0.1455])
    set_GWPCF_Multi(DfUltraBuffer_wash, ['KH2PO4', 'NaCl'], [0.8472, 0.1455])
    
    load_process_settings()
    
    # Get fermentation parameters
    params = get_fermentation_parameters()
    
    # Create reaction systems
    reactions = create_fermentation_reactions(params)
    
    # =========================================================================
    # BUILD PROCESS AREAS
    # =========================================================================
    
    # Area 200: Media Preparation
    area_200 = create_area_200_media_prep(SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt)
    
    # Area 300: Conversion (Fermentation)
    area_300 = create_area_300_conversion(
        seed_in=area_200['seed_out'],
        glucose_in=area_200['glucose_out'],
        ammonia_in=area_200['ammonia_out'],
        vent1=vent1,
        vent2=vent2,
        reactions=reactions,
        params=params,
    )
    
    # Area 400: Recovery
    area_400 = create_area_400_recovery(
        broth_in=area_300['broth_out'],
        DfUltraBuffer_wash=DfUltraBuffer_wash,
    )
    
    # Area 500: Purification
    area_500 = create_area_500_purification(
        clarified_lysate=area_400['clarified_lysate'],
        DfUltraBuffer=DfUltraBuffer,
    )
    
    # Area 600: Formulation
    area_600 = create_area_600_formulation(
        concentrated_product=area_500['concentrated_product'],
        LegHb_3=LegHb_3,
    )
    
    # Area 900: Facilities
    area_900 = create_area_900_facilities(
        area_400=area_400,
        area_500=area_500,
        effluent1=effluent1,
        use_area_convention=use_area_convention,
    )
    
    # Update input stream prices
    s.update_all_input_stream_prices(
        streamlist=[SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, DfUltraBuffer_wash]
    )
    
    return LegHb_3, vent1, vent2, effluent1, DfUltraBuffer, DfUltraBuffer_wash


# =============================================================================
# DESIGN SPECIFICATION FUNCTIONS
# =============================================================================
def set_production_rate(system, target_production_rate_kg_hr, verbose=True):
    """
    Adjust system inputs to achieve target LegHb_3 production rate.
    
    Parameters
    ----------
    system : biosteam.System
        The LegHb production system
    target_production_rate_kg_hr : float
        Target mass flow rate for LegHb_3 stream [kg/hr]
    verbose : bool
        Whether to print progress information
    
    Returns
    -------
    float
        Achieved production rate [kg/hr]
    """
    import flexsolve as flx
    import warnings
    
    log = print if verbose else (lambda *args, **kwargs: None)
    
    product_stream = system.flowsheet.stream.LegHb_3
    
    baseline_flows = {}
    for stream in system.ins:
        if stream.F_mass > 0:
            baseline_flows[stream] = stream.F_mass
    
    if not baseline_flows:
        raise ValueError("No input streams with positive flow rates found")
    
    log(f"\n{'='*60}")
    log(f"Setting Production Rate to {target_production_rate_kg_hr:.2f} kg/hr")
    log(f"{'='*60}")
    log(f"Baseline input streams stored: {len(baseline_flows)} streams")
    
    log("Running initial simulation to establish baseline...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        system.simulate()
    
    initial_production = product_stream.F_mass
    
    if initial_production <= 0:
        raise ValueError("Initial production rate is zero. Cannot scale system.")
    
    initial_guess = target_production_rate_kg_hr / initial_production
    
    log(f"  Initial production: {initial_production:.2f} kg/hr")
    log(f"  Target production:  {target_production_rate_kg_hr:.2f} kg/hr")
    log(f"  Initial scaling guess: {initial_guess:.4f}x")
    
    iteration = [0]
    
    def objective_function(scaling_factor):
        iteration[0] += 1
        try:
            for stream, baseline_flow in baseline_flows.items():
                stream.F_mass = baseline_flow * scaling_factor
            system.simulate()
            achieved_rate = product_stream.F_mass
            error = achieved_rate - target_production_rate_kg_hr
            if iteration[0] % 5 == 0 or abs(error) < 1.0:
                log(f"    Iteration {iteration[0]}: scale={scaling_factor:.4f}x, "
                    f"production={achieved_rate:.2f} kg/hr, error={error:.4f} kg/hr")
            return error
        except Exception as e:
            log(f"    Iteration {iteration[0]} FAILED at scale={scaling_factor:.4f}x: {e}")
            return 1e6 if scaling_factor > initial_guess else -1e6
    
    try:
        log(f"\nSolving for optimal scaling factor...")
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            x0 = max(0.01, initial_guess * 0.1)
            x1 = max(5.0, initial_guess * 2.0)
            y0 = objective_function(x0)
            y1 = objective_function(x1)
            expand_iter = 0
            while y0 * y1 > 0 and expand_iter < 8:
                expand_iter += 1
                if abs(y1) <= abs(y0):
                    x0 = max(0.01, x0 / 2.0)
                    y0 = objective_function(x0)
                else:
                    x1 *= 2.0
                    y1 = objective_function(x1)
            log(f"  Search bounds: {x0:.4f}x to {x1:.4f}x (bracket {'OK' if y0*y1 <= 0 else 'NOT FOUND'})")
            
            scaling_factor = flx.IQ_interpolation(
                f=objective_function,
                x0=x0,
                x1=x1,
                x=initial_guess,
                xtol=0.0001,
                ytol=0.01,
                maxiter=100,
                checkbounds=(y0 * y1 <= 0),
                checkiter=True,
            )
        
        log(f"\nSolver converged. Running final validation simulation...")
        for stream, baseline_flow in baseline_flows.items():
            stream.F_mass = baseline_flow * scaling_factor
        
        system.simulate()
        achieved_rate = product_stream.F_mass
        
        log(f"\n✓ Successfully achieved production rate:")
        log(f"  Target:          {target_production_rate_kg_hr:.2f} kg/hr")
        log(f"  Achieved:        {achieved_rate:.2f} kg/hr")
        log(f"  Scaling factor:  {scaling_factor:.4f}x")
        log(f"  Error:           {abs(achieved_rate - target_production_rate_kg_hr):.4f} kg/hr")
        log(f"  Total iterations: {iteration[0]}")
        log(f"{'='*60}\n")
        
        return achieved_rate
        
    except Exception as e:
        log(f"\n✗ Failed to achieve target production rate: {e}")
        log(f"  Error type: {type(e).__name__}")
        log(f"  Total iterations attempted: {iteration[0]}")
        
        log(f"\n  Restoring baseline input flows...")
        for stream, baseline_flow in baseline_flows.items():
            stream.F_mass = baseline_flow
        
        try:
            system.simulate()
            log(f"  System restored to baseline: {product_stream.F_mass:.2f} kg/hr")
        except Exception as restore_error:
            log(f"  Warning: Could not restore baseline simulation: {restore_error}")
        
        raise ValueError(f"Could not achieve target production rate of {target_production_rate_kg_hr:.2f} kg/hr: {e}")


def check_LegHb_specifications(product_stream):
    """
    Verify that LegHb_3 product stream meets composition and purity specifications.
    
    Parameters
    ----------
    product_stream : thermosteam.Stream
        The LegHb_3 product stream to check
    
    Raises
    ------
    ValueError
        If any specification is not met
    """
    print(f"\n{'='*60}")
    print(f"Product Specification Check: {product_stream.ID}")
    print(f"{'='*60}")
    
    total_mass = product_stream.F_mass
    
    if total_mass <= 0:
        raise ValueError("Product stream has zero mass flow")
    
    def get_chemical_mass(chemical_id):
        try:
            return product_stream.imass[chemical_id]
        except:
            return 0
    
    def get_mass_percent(chemical_ids):
        if isinstance(chemical_ids, str):
            chemical_ids = [chemical_ids]
        total = sum(get_chemical_mass(chem) for chem in chemical_ids)
        return total / total_mass * 100
    
    specs = {
        'Fat (OleicAcid)': {
            'chemicals': ['OleicAcid'],
            'target': (0, 2),
            'actual': None
        },
        'Carbohydrates': {
            'chemicals': ['Glucan', 'Glucose', 'Chitin'],
            'target': (0, 4),
            'actual': None
        },
        'Product (Leghemoglobin)': {
            'chemicals': ['Leghemoglobin'],
            'target': (6, 9),
            'actual': None
        },
        'Total Solids': {
            'chemicals': [chem for chem in product_stream.chemicals.IDs if chem != 'H2O'],
            'target': (0, 24),
            'actual': None
        }
    }
    
    for spec_name, spec_data in specs.items():
        spec_data['actual'] = get_mass_percent(spec_data['chemicals'])
    
    legh_mass = get_chemical_mass('Leghemoglobin')
    globin_mass = get_chemical_mass('Globin')
    mannoprotein_mass = get_chemical_mass('Mannoprotein')
    total_protein_mass = legh_mass + globin_mass + mannoprotein_mass
    
    protein_purity = (legh_mass / total_protein_mass * 100) if total_protein_mass > 0 else 0
    
    all_passed = True
    
    print(f"\n{'Specification':<30} {'Target Range':<25} {'Actual':<15} {'Status'}")
    print(f"{'-'*85}")
    
    for spec_name, spec_data in specs.items():
        min_val, max_val = spec_data['target']
        actual = spec_data['actual']
        passed = min_val <= actual <= max_val
        
        status = '✓ PASS' if passed else '✗ FAIL'
        target_str = f"{min_val:.1f}% - {max_val:.1f}%"
        print(f"{spec_name:<30} {target_str:<25} {actual:>6.2f}%{'':<7} {status}")
        
        if not passed:
            all_passed = False
    
    purity_passed = protein_purity >= 65.0
    status = '✓ PASS' if purity_passed else '✗ FAIL'
    print(f"{'Protein Purity':<30} {'>= 65.0%':<25} {protein_purity:>6.2f}%{'':<7} {status}")
    
    if not purity_passed:
        all_passed = False
    
    print(f"{'='*85}")
    
    print(f"\nProduct Stream Summary:")
    print(f"  Total Flow Rate: {total_mass:.2f} kg/hr")
    print(f"  Water Content:   {get_mass_percent('H2O'):.2f}%")
    print(f"  Leghemoglobin:   {legh_mass:.4f} kg/hr ({get_mass_percent('Leghemoglobin'):.2f}%)")
    
    if not all_passed:
        print(f"\n{'='*85}")
        raise ValueError("Product does not meet one or more specifications. See details above.")
    
    print(f"\n✓ All specifications met!")
    print(f"{'='*85}\n")
    
    return True


# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == '__main__':
    bst.preferences.N = 50
    nn = 1
    TARGET_PRODUCTION = 275 * nn
    
    print("="*85)
    print("LEGHEMOGLOBIN PRODUCTION SYSTEM - MODULAR DESIGN")
    print("="*85)
    
    print("\n1. Creating system...")
    LegHb_sys = create_LegHb_system()
    sys = LegHb_sys
    f = sys.flowsheet
    u_units = f.unit
    ss = f.stream
    sys.operating_hours = 8000
    
    print("\n2. Running baseline simulation...")
    try:
        LegHb_sys.simulate()
        baseline_production = ss.LegHb_3.F_mass
        print(f"   Baseline production rate: {baseline_production:.2f} kg/hr")
    except Exception as e:
        print(f"   Baseline simulation failed: {e}")
        raise
    
    print(f"\n3. Applying design specification: TARGET_PRODUCTION = {TARGET_PRODUCTION} kg/hr")
    try:
        achieved_production = set_production_rate(LegHb_sys, TARGET_PRODUCTION)
        LegHb_sys.simulate()
        final_production = ss.LegHb_3.F_mass
        
        if abs(final_production - TARGET_PRODUCTION) > 1.0:
            print(f"\n   WARNING: Production rate drifted after final simulation!")
            print(f"   Target:  {TARGET_PRODUCTION:.2f} kg/hr")
            print(f"   Actual:  {final_production:.2f} kg/hr")
            
    except Exception as e:
        print(f"\n   Could not achieve target production: {e}")
        print("   Continuing with baseline production rate...")
        achieved_production = ss.LegHb_3.F_mass
    
    print(f"\n4. Verifying product specifications...")
    try:
        check_LegHb_specifications(ss.LegHb_3)
    except ValueError as e:
        print(f"\n   SPECIFICATION CHECK FAILED: {e}")
        print("   System may require process parameter adjustments to meet specifications.")
    
    print(f"\n5. System Summary")
    print("="*85)
    LegHb_sys.show()
    
    legh_purity = ss.LegHb_3.imass['Leghemoglobin'] / ss.LegHb_3.F_mass * 100
    print(f"\n{'='*85}")
    print("KEY PERFORMANCE INDICATORS")
    print(f"{'='*85}")
    print(f"  Product Stream:           {ss.LegHb_3.ID}")
    print(f"  Production Rate:          {ss.LegHb_3.F_mass:.2f} kg/hr")
    print(f"  Leghemoglobin Content:    {legh_purity:.2f}%")
    print(f"  Annual Production:        {ss.LegHb_3.F_mass * sys.operating_hours / 1000:.2f} metric tons/year")
    print(f"{'='*85}\n")
    
    print(f"\n6. Generating system diagram...")
    LegHb_sys.diagram(format='html', display=True)
    
    print(f"\n7. Performing LCA analysis...")
    try:
        r1 = bst.report.lca_inventory_table(
            systems=[sys],
            keys='GWP',
            items=[ss.LegHb_3],
        )
        print("   LCA Inventory Table generated")
        
        r2 = bst.report.lca_displacement_allocation_table(
            systems=[sys],
            key='GWP',
            items=[ss.LegHb_3],
        )
        print("   LCA Displacement Allocation Table generated")
    except Exception as e:
        print(f"   LCA analysis failed: {e}")
    
    print(f"\n{'='*85}")
    print("SIMULATION COMPLETE")
    print(f"Target Production:   {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"Achieved Production: {ss.LegHb_3.F_mass:.2f} kg/hr")
    print(f"{'='*85}\n")
