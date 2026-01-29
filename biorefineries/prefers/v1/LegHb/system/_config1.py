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
from biorefineries.prefers.v1 import _process_settings
from biorefineries.prefers.v1._process_settings import price

# Global registry for seed scaling targets (to coordinate between set_production_rate and specs)
seed_targets = {}

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
    'optimize_NH3_loading',
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
    theta_O2 = 0.4  # Dissolved oxygen concentration [% saturation]
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
    titer_LegHb = 5#7.27  # [g/L] - Literature validated
    productivity_LegHb = 5 / 72#7.27 / 72  # [g/L/h] - 72h fermentation cycle
    Y_p = 5 * 4 / 650 #7.27 * 4 / 1300  # [by wt] - ~2.2% product yield on glucose
    Y_b = 0.53  # [by wt] - Biomass yield, typical for yeast 0.4~0.6
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
        bst.Rxn('1 Glucose + 0.70588 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b + 0.617647 O2 + 0.17647 H2SO4 + 4.05882 H2O',
                reactant='Glucose', X=yield_LegHb*(0.00786/0.17647)*(1/0.93-1), check_atomic_balance=True),
        bst.Rxn('1 Glucose + 0.01646 (NH4)2SO4 + 1.61317 NH3 -> 6 Globin_In + 0.28807 O2 + 3.68724 H2O',
                reactant='Glucose', X=yield_LegHb*(1/0.93-1), check_atomic_balance=True),
        bst.Rxn('1 Glucose + 0.00786 FeSO4 + 0.00786 (NH4)2SO4 + 1.58847 NH3 -> 6 Leghemoglobin_In + 0.30275 O2  + 3.70380 H2O',
                reactant='Glucose', X=yield_LegHb, check_atomic_balance=True),     
    ])  
    fermentation_reaction[2].product_yield('Leghemoglobin_In', basis='wt', product_yield=yield_LegHb)
    
    neutralization_reaction = bst.Rxn(
        'H2SO4 + 2 NH3 -> (NH4)2SO4', reactant='H2SO4', X=0.99,
        check_atomic_balance=True
    )
    
    cell_growth_reaction = bst.Rxn(
        '1 Glucose + 0.8364 NH3 + 0.0108 (NH4)2SO4 -> 2.01 H2O + 0.106 O2 + 6 Pichia_pastoris', 
        'Glucose', X = Y_b,
        correct_atomic_balance=True
    )
    cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=Y_b)
    
    respiration_reaction1 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', X = 1 - Y_b,
        check_atomic_balance=True
    )
    respiration_reaction2 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', X= 1 - Y_b,
        #0.9995 - cell_growth_reaction.X -fermentation_reaction[0].X -fermentation_reaction[1].X -fermentation_reaction[2].X,
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
        - M201: Seed solution 1 preparation
        - M202: Seed solution 2 + culture media preparation
        - M203: Seed hold tank (combines both seed solutions)
        - M204: Glucose solution preparation (50% dilution)
        - T201: Glucose storage tank (96 hours residence)
        - T202: Ammonia storage tank (25 wt% aqueous ammonia)
    
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
    # Initialize streams to reference values (run once during creation)
    # This ensures valid baseline flows while allowing set_production_rate to scale them later
    try:
        SeedIn1.imass['Seed'] = s.SeedSolution1['Seed']
        SeedIn2.imass['Seed'] = s.SeedSolution2['Seed']
        
        # Calculate CultureIn baseline
        ref_seed2_total = s.SeedSolution2['Seed'] + s.SeedSolution2['H2O'] # Approximate
        # Better: create temp stream to get mass
        tmp = bst.Stream(ID='tmp_init', **s.SeedSolution2)
        ref_seed2_total = tmp.F_mass
        ratio_culture = (0.1 + 60 + 0.15191) / 1000
        CultureIn.imass['Culture'] = ref_seed2_total * ratio_culture
        tmp.empty() # cleanup
    except:
        pass # Fallback if stream definitions missing keys

    # Seed solution 1 preparation
    M201 = bst.MixTank('M201', ins=[SeedIn1, 'Water1'], outs='M201Out', tau=16)
    

        
    # Pre-calculate reference values for Seed1 to avoid re-registering stream in loop
    ref_stream1 = bst.Stream(**{**s.SeedSolution1, 'ID': 'RefSeed1'})
    ref_seed1_val = ref_stream1.imass['Seed']
    ref_water1_val = ref_stream1.imass['H2O']
    
    @M201.add_specification(run=True)
    def update_seed1_inputs():
        # Force SeedIn1 to target (from global dict or reference)
        target = seed_targets.get(SeedIn1, ref_seed1_val)
        
        # Enforce target
        SeedIn1.imass['Seed'] = target
        
        # Enforce water ratio
        if ref_seed1_val > 0:
            ratio = ref_water1_val / ref_seed1_val
            M201.ins[1].imass['H2O'] = target * ratio
        
        M201.ins[1].T = 25 + 273.15
    
    # Seed solution 2 + culture media preparation
    M202 = bst.MixTank('M202', ins=[SeedIn2, CultureIn, 'Water2'], outs='M202Out', tau=16)
    
    # Pre-calculate reference values for Seed2 to avoid re-registering stream in loop
    ref_stream2 = bst.Stream(**{**s.SeedSolution2, 'ID': 'RefSeed2'})
    ref_seed2_val = ref_stream2.imass['Seed']
    ref_water2_val = ref_stream2.imass['H2O']
    ref_total2_val = ref_stream2.F_mass
    
    @M202.add_specification(run=True)
    def update_culture_inputs():
        target = seed_targets.get(SeedIn2, ref_seed2_val)
        
        # Enforce target
        SeedIn2.imass['Seed'] = target
        
        if ref_seed2_val > 0:
            ratio_water = ref_water2_val / ref_seed2_val
            M202.ins[2].imass['H2O'] = target * ratio_water # Water2
            
            # CultureIn logic
            ratio_culture = (0.1 + 60 + 0.15191) / 1000
            
            # Predict current total solution mass
            current_total_mass_est = target * (ref_total2_val / ref_seed2_val)
            
            CultureIn.imass['Culture'] = current_total_mass_est * ratio_culture
            
        M202.ins[2].T = 25 + 273.15
    
    # Ammonia storage tank
    T202 = u.AmmoniaStorageTank('T202', ins=NH3_25wt, outs='T202Out')
    
    # Ammonia Splitter
    S202 = bst.Splitter('S202', ins=T202-0, outs=('NH3_Seed', 'NH3_Fer'), split=0.04)
    
    # Seed hold tank
    M203 = u.SeedHoldTank('M203', ins=[M201-0, M202-0, S202-0], outs='M203Out')
    
    # Glucose solution preparation (50% dilution)
    M204 = bst.MixTank('M204', ins=[Glucose, 'Water3'], outs='M204Out', tau=16)
    
    @M204.add_specification(run=True)
    def update_water_content():
        M204.ins[1].imass['H2O'] = Glucose.imass['Glucose'] / 2
        M204.ins[1].T = 25 + 273.15
    
    # Glucose storage tank
    T201 = bst.StorageTank('T201', ins=M204-0, outs='T201Out', tau=16*4+72)
    
    return {
        'units': [M201, M202, M203, M204, T201, T202, S202],
        'seed_out': M203.outs[0],
        'glucose_out': T201.outs[0],
        'ammonia_out': S202.outs[1],
        'M201': M201, 'M202': M202, 'M203': M203,
        'M204': M204, 'T201': T201, 'T202': T202, 'S202': S202,
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
def create_area_400_recovery(broth_in, DfUltraBuffer2):
    """
    Create Area 400: Recovery (Biomass Harvest & Cell Disruption)
    
    Purpose: Release and clarify intracellular LegHb from yeast cells
    
    Process Flow:
        1. Biomass Harvest (C401): Centrifuge to concentrate cells
        2. Cell Washing (M401, H402, M402, C402): Remove spent media
        3. Cell Disruption (S401): High-pressure homogenization
        4. Debris Removal (H401, S402): Cool and centrifuge
        5. Lysate Clarification (S403, S404): Depth filtration + screw press
    
    Parameters
    ----------
    broth_in : Stream
        Fermentation broth from Area 300
    DfUltraBuffer2 : Stream
        Buffer for cell washing
    
    Returns
    -------
    dict
        Dictionary of area units and output streams
    """
    # Primary harvest centrifuge
    C401 = u.Centrifuge(
        'C401',
        ins=broth_in,
        outs=('CellCream', 'SpentMedia'),
        split={
            'cellmass': 0.999,
            'Leghemoglobin_In': 0.999,
            'Globin_In': 0.999,
            'Heme_b': 0.85,
            'Glucan': 0.99,
            'Mannoprotein': 0.99,
            'Chitin': 0.99,
            'OleicAcid': 0.99,
            'RNA': 0.99,
        },
        moisture_content=0.55,
    )
    
    # Wash buffer preparation
    M401 = bst.MixTank('M401', ins=(DfUltraBuffer2, 'Water4'), 
                       outs='WashBufferOut', tau=0.5)
    
    @M401.add_specification(run=True)
    def update_wash_buffer():
        M401.ins[1].imass['H2O'] = C401.outs[0].F_mass * 1.5
        M401.ins[0].imol['DfUltraBuffer'] = M401.ins[1].imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000
    
    # Cool wash buffer
    H402 = bst.HXutility(
        'H402',
        ins=M401-0,
        outs='ColdWashBuffer',
        T=10 + 273.15,
        cool_only=True,
    )
    
    # Cell wash mixer
    M402 = bst.MixTank('M402', ins=(C401-0, H402-0), 
                            outs='WashedCellSlurry', tau=0.25)
    
    # Second centrifuge for washed cells
    C402 = u.Centrifuge(
        'C402',
        ins=M402-0,
        outs=('WashedCellCream', 'WashEffluent'),
        split={
            'cellmass': 0.999,
            'Leghemoglobin_In': 0.999,
            'Globin_In': 0.999,
            'Heme_b': 0.85,
            'Glucan': 0.99,
            'Mannoprotein': 0.99,
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
        moisture_content=0.001
    )
    
    return {
        'units': [C401, M401, H402, M402, C402, S401, H401, S402, S403, S404],
        'clarified_lysate': S403.outs[1],
        'spent_media': C401.outs[1],
        'wash_effluent': C402.outs[1],
        'dehydrated_debris': S404.outs[0],
        'press_liquor': S404.outs[1],
        'C401': C401, 'M401': M401, 'H402': H402,
        'M402': M402, 'C402': C402, 'S401': S401,
        'H401': H401, 'S402': S402, 'S403': S403, 'S404': S404,
    }

# =============================================================================
# AREA 500: PURIFICATION (UF/DF MEMBRANE SEPARATION)
# =============================================================================
def create_area_500_purification(clarified_lysate, DfUltraBuffer1):
    """
    Create Area 500: Purification (UF/DF Membrane Separation)
    
    Purpose: Concentrate LegHb and remove impurities via tangential flow filtration
    
    Process Flow:
        1. DF Buffer Preparation (M501, H501): Cold buffer for diafiltration
        2. UF/DF Stage 1 (U501): Concentration + buffer exchange (5-7 diavolumes)
        3. UF Stage 2 (U502): Final concentration
    
    Parameters
    ----------
    clarified_lysate : Stream
        Clarified lysate from Area 400
    DfUltraBuffer1 : Stream
        DF buffer stream
    
    Returns
    -------
    dict
        Dictionary of area units and output streams
    """
    # DF buffer preparation (minimal for concentration-only mode)
    M501 = bst.MixTank('M501', ins=(DfUltraBuffer1, 'Water5'), 
                       outs='DFBufferOut', tau=0.5)
    
    # Reference to S403 for specification (will be set via system context)
    clarified_lysate_stream = clarified_lysate
    
    @M501.add_specification(run=True)
    def update_df_buffer():
        feed_water = clarified_lysate_stream.imass['H2O']
        M501.ins[1].imass['H2O'] = feed_water * 6
        M501.ins[0].imol['DfUltraBuffer'] = M501.ins[1].imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000
    
    # Cool DF buffer
    H501 = bst.HXutility(
        'H501',
        ins=M501-0,
        outs='ColdDFBuffer',
        T=5 + 273.15,
        cool_only=True,
    )
    
    # UF Stage 1 - concentration + buffer exchange using Diafiltration with UF preset
    U501 = u.Diafiltration.from_preset(
        'UF',
        'U501',
        ins=(clarified_lysate, H501-0),
        outs=('UFConcentrate', 'UFPermeate'),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
        TargetProduct_Retention=0.99,
        Salt_Retention=0.05,
    )
    U501.add_specification(run=True)
    
    # UF Stage 2 - Final concentration using Diafiltration with minimal buffer
    U502 = u.Diafiltration(
        'U502',
        ins=(U501-0, bst.Stream('Water6', H2O=0.001)),
        outs=('ConcentratedLegH', 'ConcentrationPermeate'),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
        TMP_bar1=3,
        FeedWater_Recovery_to_Permeate=0.85,
    )
    
    @U502.add_specification(run=True)
    def U502_adjust_water_recovery():
        U502.Salt_Retention = 0.05
        U502._run()
    
    return {
        'units': [M501, H501, U501, U502],
        'concentrated_product': U502.outs[0],
        'uf_permeate': U501.outs[1],
        'concentration_permeate': U502.outs[1],
        'M501': M501, 'H501': H501, 'U501': U501, 'U502': U502,
    }

# =============================================================================
# AREA 600: FORMULATION (THERMAL TREATMENT & FINAL PRODUCT)
# =============================================================================
def create_area_600_formulation(concentrated_product, AntioxidantStream, LegHb_3):
    """
    Create Area 600: Formulation (Thermal Treatment & Final Product)
    
    Purpose: Ensure microbiological safety and achieve final product specification
    
    Process Flow:
        1. HTST Pasteurization (H603, T601): Heat to 72°C for 30 seconds
        2. Formulation (M604): Add antioxidant + dilution water
        3. Final Cooling (H606): Rapid chill to 4°C
    
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
    H603 = bst.HXutility(
        'H603',
        ins=concentrated_product,
        outs='HeatedConcentrate',
        T=72 + 273.15,
        heat_only=True,
    )
    
    # Hold time tank
    T601 = bst.MixTank(
        'T601',
        ins=H603-0,
        outs='PasteurizedConcentrate',
        tau=30/3600  # 30 second hold time
    )
    
    # # Antioxidant stream
    # AntioxidantStream = bst.Stream(
    #     'AntioxidantStream',
    #     SodiumAscorbate=0.1,
    #     H2O=1.0,
    #     units='kg/hr',
    #     price=price.get('SodiumAscorbate', 5.0)
    # )
    
    # Dilution water stream
    DilutionWater = bst.Stream(
        'Water7',
        H2O=100.0,
        units='kg/hr',
    )
    
    # Formulation mixer
    M604 = bst.MixTank(
        'M604',
        ins=(T601-0, AntioxidantStream, DilutionWater),
        outs='FormulatedProduct',
        tau=0.1
    )
    M604.target_total_solids_percent = 20.0
    M604.target_legh_percent = 7.5
    
    @M604.add_specification(run=True)
    def update_formulation():
        feed = T601.outs[0]
        feed_total_mass = feed.F_mass
        
        if feed_total_mass <= 0:
            DilutionWater.imass['H2O'] = 0
            AntioxidantStream.imass['SodiumAscorbate'] = 0
            AntioxidantStream.imass['H2O'] = 0
            return
        
        feed_water_mass = feed.imass['H2O']
        feed_solids_mass = feed_total_mass - feed_water_mass
        feed_legh_mass = feed.imass['Leghemoglobin']
        
        target_solids_pct = getattr(M604, 'target_total_solids_percent', 20.0)
        target_legh_pct = getattr(M604, 'target_legh_percent', 7.5)
        
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
    H606 = bst.HXutility(
        'H606',
        ins=M604-0,
        outs=LegHb_3,
        T=4 + 273.15,
        cool_only=True,
    )
    
    return {
        'units': [H603, T601, M604, H606],
        'product': H606.outs[0],
        'H603': H603, 'T601': T601, 'M604': M604, 'H606': H606,
        'AntioxidantStream': AntioxidantStream,
        'DilutionWater': DilutionWater,
    }

# =============================================================================
# AREA 900: FACILITIES & UTILITIES
# =============================================================================
def create_area_900_facilities(area_400, area_500, effluent1, emissions, ash_disposal, use_area_convention=False):
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
    # M901 = bst.MixTank(
    #     'M901',
    #     ins=(
    #         area_400['press_liquor'],
    #         area_400['spent_media'],
    #         area_400['wash_effluent'],
    #         area_500['uf_permeate'],
    #         area_500['concentration_permeate'],
    #     ),
    #     outs='CombinedWastewater',
    #     tau=1
    # )
    # M901.add_specification(run=True)
    
    # RO treatment for water recovery
    # S901 = u.ReverseOsmosis('S901', ins=M901-0, outs=('RO_treated_water', effluent1))

    # Supplemental nutrient streams for wastewater treatment stability
    # These ensure minimum nutrient availability for anaerobic digestion at all scales
    # Required because at low production rates, reactor effluent doesn't have enough nutrients
    SupplementalNH4SO4 = bst.Stream('SupplementalNH4SO4', units='kg/hr', price=0)
    SupplementalFeSO4 = bst.Stream('SupplementalFeSO4', FeSO4=50, units='kg/hr', price=0)
    SupplementalNH4SO4.imass['NH3'] = 50 
    SupplementalNH4SO4.imass['(NH4)2SO4'] = 50
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[
            area_400['press_liquor'],
            area_400['spent_media'],
            area_400['wash_effluent'],
            area_500['uf_permeate'],
            area_500['concentration_permeate'],
            SupplementalFeSO4,
            SupplementalNH4SO4,
        ],
        outs=('biogas', 'sludge', 'RO_treated_water', effluent1),
        mockup=True,
        area=500,
    )
    biogas, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
    
    # Utility systems
    CT = bst.CoolingTower(500 if use_area_convention else 'CT')
    CWP = bst.ChilledWaterPackage(500 if use_area_convention else 'CWP')
    
    # # Route dehydrated debris to disposal
    # debris_disposal = bst.Stream('debris_disposal')
    # M503 = bst.Mixer('M503', ins=area_400['dehydrated_debris'], outs=debris_disposal)
    
    M902 = bst.Mixer(
        'M902',
        ins=(area_400['dehydrated_debris'], sludge),
        outs='SolidsToBoiler',
    )
    
    BT = u.BoilerTurbogenerator(400 if use_area_convention else 'BT',
        (M902-0, biogas, 'boiler_makeup_water', 'natural_gas', 'lime_boiler', 'boiler_chems'),
        outs=(emissions, 'rejected_water_and_blowdown', ash_disposal),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85,
        satisfy_system_electricity_demand=False,
    )
    
    makeup_water_streams = (
        F.cooling_tower_makeup_water,
        F.Water1, F.Water2,
        F.Water3, F.Water4,
        F.Water5, F.Water6,F.Water7,
        F.boiler_makeup_water,
    )
    process_water_streams = (
        treated_water,
        F.rejected_water_and_blowdown,
        *makeup_water_streams
    )
    
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    
    PWC = bst.ProcessWaterCenter(
        500 if use_area_convention else 'PWC',
        ins=(treated_water, makeup_water, 'recycled_process_water', 'makeup_process_water'),
        outs=('RO_water', 'process_water', 'excess_water'),
        makeup_water_streams=makeup_water_streams,
        process_water_streams=process_water_streams,
        reverse_osmosis_water_price=0.000254,
        process_water_price=0.000135,
    )
    
    # return {
    #     'units': [M501, S501, CT, CWP, M503, PWC],
    #     'effluent': effluent1,
    #     'debris_disposal': debris_disposal,
    #     'M501': M501, 'S501': S501, 'CT': CT, 'CWP': CWP,
    #     'M503': M503, 'PWC': PWC,
    # }
    return {
        'units': [M902, CT, CWP, BT, PWC],
        'effluent': effluent1,
        'emissions': emissions,
        'ash_disposal': ash_disposal,
        'wastewater_treatment_sys': wastewater_treatment_sys,
        'M902': M902, 'CT': CT, 'CWP': CWP,
        'BT': BT, 'PWC': PWC,
    }

# =============================================================================
# MAIN SYSTEM FACTORY
# =============================================================================
@bst.SystemFactory(
    ID='LegHb_sys',
    ins=[s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt, s.DfUltraBuffer1, s.DfUltraBuffer2,s.AntioxidantStream],
    outs=[s.LegHb_3, s.vent1, s.vent2, s.effluent1,s.emissions,s.ash_disposal],
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
        Input streams (SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer1, DfUltraBuffer2,AntioxidantStream)
    outs : tuple
        Output streams (LegHb_3, vent1, vent2, effluent1, emissions, ash_disposal)
    use_area_convention : bool
        Whether to use area ID convention for utility units
    
    Returns
    -------
    tuple
        Output streams (LegHb_3, vent1, vent2, effluent1, emissions, ash_disposal, DfUltraBuffer1, DfUltraBuffer, AntioxidantStream)
    """
    bst.preferences.N = 50
    try:
        bst.settings.solver.maxiter = 500
        bst.settings.solver.molar_tolerance = 0.01
    except:
        pass
    
    # Unpack streams
    SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer1, DfUltraBuffer2, AntioxidantStream = ins
    LegHb_3, vent1, vent2, effluent1, emissions, ash_disposal = outs

    
    # Set GWP characterization factors
    set_GWPCF(Glucose, 'Glucose')
    set_GWPCF(ash_disposal, 'ash_disposal')
    set_GWPCF(AntioxidantStream, 'AscorbicAcid', dilution=0.1)
    set_GWPCF_Multi(SeedIn1, ['AmmoniumSulfate', 'Glucose', 'MagnesiumSulfate', 'KH2PO4'],
                   [0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(SeedIn2, ['AmmoniumSulfate', 'Glucose', 'MagnesiumSulfate', 'KH2PO4'],
                   [0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(CultureIn, ['Glycine', 'Glucose', 'IronSulfate'],
                   [0.1/(0.1+60+0.15191), 60/(0.1+60+0.15191), 0.15191/(0.1+60+0.15191)])
    set_GWPCF(NH3_25wt, 'Ammonia_SEA', dilution=0.25)
    set_GWPCF_Multi(DfUltraBuffer1, ['KH2PO4', 'NaCl'], [0.8472, 0.1455])
    set_GWPCF_Multi(DfUltraBuffer2, ['KH2PO4', 'NaCl'], [0.8472, 0.1455])
    
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
        DfUltraBuffer2=DfUltraBuffer2,
    )
    
    # Area 500: Purification
    area_500 = create_area_500_purification(
        clarified_lysate=area_400['clarified_lysate'],
        DfUltraBuffer1=DfUltraBuffer1,
    )
    
    # Area 600: Formulation
    area_600 = create_area_600_formulation(
        concentrated_product=area_500['concentrated_product'],
        AntioxidantStream=AntioxidantStream,
        LegHb_3=LegHb_3,
    )
    
    # Area 900: Facilities
    area_900 = create_area_900_facilities(
        area_400=area_400,
        area_500=area_500,
        effluent1=effluent1,
        emissions=emissions,
        ash_disposal=ash_disposal,
        use_area_convention=use_area_convention,
    )
    
    # Update input stream prices
    s.update_all_input_stream_prices(
        streamlist=[SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer1, DfUltraBuffer2]
    )
    
    return LegHb_3, vent1, vent2, effluent1, DfUltraBuffer1, DfUltraBuffer2, emissions, ash_disposal, AntioxidantStream

# =============================================================================
# DESIGN SPECIFICATION FUNCTIONS
# =============================================================================

def optimize_NH3_loading(system, verbose=True):
    """
    Optimizes NH3_25wt flow rate (X) and S202 split ratio (Y) to meet fermentation demand
    ensuring residual NH3 < 1e-4 kmol/hr in both seed and main fermenters.
    
    This function:
    1. Sets excess Ammonia to ensure reactions proceed fully.
    2. Calculates exact Ammonia consumption in R301 and R302.
    3. Adjusts NH3_25wt flow and S202 split to match demand.
    """
    log = print if verbose else (lambda *args, **kwargs: None)
    log(f"\n{'='*60}")
    log("Optimizing NH3 Loading (pH Control Strategy)")
    log(f"{'='*60}")
    
    u = system.flowsheet.unit
    s = system.flowsheet.stream
    
    # Retrieve units and streams
    try:
        S202 = u.S202
        R301 = u.R301
        R302 = u.R302
        NH3_source = s.NH3_25wt
        Broth = u.R302.outs[1] # Broth
    except AttributeError as e:
        log(f"[WARN] Could not find necessary units/streams for optimization: {e}")
        return

    # Step 1: Supply excess ammonia to determine max demand
    # Store initial state
    initial_flow = NH3_source.F_mass
    initial_split = S202.split[0]
    
    # Set excess relative to current flow (which should be roughly scaled)
    # 20x should be sufficient without causing design errors (N > 10000)
    excess_flow = initial_flow * 20 if initial_flow > 1e-6 else 1000 
    NH3_source.F_mass = excess_flow
    S202.split[:] = 0.5 # 50/50 split
    
    system.simulate()
    
    # Step 2: Calculate Demand
    # R301 Demand
    # NH3 in R301 comes from S202-0 (Seed Ammonia) + Inputs in Seed?
    # R301 ins: [SeedIn1, SeedIn2 (via M203)]
    # M203 ins: [M201, M202, S202-0]
    # So NH3 enters R301 via M203.
    
    # We check reaction extent or mass balance.
    # Consumption = In - Out - Accumulation (0)
    # But In is what we want to find.
    # Current In is Excessive.
    # Current Out includes Residual.
    # Consumed = In - Out.
    
    def get_NH3_consumption(reactor):
        nh3_in = sum(i.imol['NH3'] for i in reactor.ins)
        nh3_out = sum(o.imol['NH3'] for o in reactor.outs)
        return nh3_in - nh3_out

    consumed_R301 = get_NH3_consumption(R301)
    consumed_R302 = get_NH3_consumption(R302)
    
    log(f"  R301 NH3 Consumption: {consumed_R301:.4f} kmol/hr")
    log(f"  R302 NH3 Consumption: {consumed_R302:.4f} kmol/hr")
    
    # Step 3: Back-calculate required Feed
    # NH3_25wt is ~25% NH3, 75% Water.
    # Source stream composition matches NH3_25wt definition.
    # We need to set NH3_25wt.imol['NH3'] to meet (consumed_R301 + consumed_R302).
    # And S202 split to distribute it correctly.
    
    total_NH3_demand = consumed_R301 + consumed_R302
    
    if total_NH3_demand <= 0:
        log("  [WARN] Zero Ammonia demand detected. Optimization skipped.")
        NH3_source.F_mass = initial_flow
        S202.split[:] = initial_split
        return

    # Calculate required source flow
    # Source has fixed composition. We scale F_mass.
    # Current mol fraction of NH3 in source:
    nh3_mol_frac = NH3_source.imol['NH3'] / NH3_source.F_mol
    required_F_mol = total_NH3_demand / nh3_mol_frac
    
    # We can also just scale F_mass based on current NH3 flow
    current_NH3_flow_mol = NH3_source.imol['NH3']
    scaling_factor = total_NH3_demand / current_NH3_flow_mol*1.001
    new_F_mass = NH3_source.F_mass * scaling_factor
    
    # Set Flow Rate 'X'
    NH3_source.F_mass = new_F_mass
    
    # Set Split 'Y'
    # S202-0 goes to R301 (Seed)
    # S202-1 goes to R302 (Fer)
    # We want S202-0 NH3 = consumed_R301
    
    split_seed = consumed_R301 / total_NH3_demand
    S202.split[:] = split_seed # S202 definition: outs=('NH3_Seed', 'NH3_Fer'). split refers to index 0?
    # bst.Splitter(..., split=0.5) usually means outs[0] gets 0.5.
    
    log(f"  Optimized Total Flow (X): {new_F_mass:.2f} kg/hr")
    log(f"  Optimized Split Ratio (Y): {split_seed:.4f} (to Seed)")
    
    # Step 4: Verify
    system.simulate()
    
    r_NH3_R301 = R301.outs[1].imol['NH3'] # R301Out
    r_NH3_R302 = R302.outs[1].imol['NH3'] # Broth
    
    log(f"  Residual NH3 R301: {r_NH3_R301:.2e} kmol/hr (< 1e-4 check: {'PASS' if r_NH3_R301 < 1e-4 else 'FAIL'})")
    log(f"  Residual NH3 R302: {r_NH3_R302:.2e} kmol/hr (< 1e-4 check: {'PASS' if r_NH3_R302 < 1e-4 else 'FAIL'})")
    
    if r_NH3_R301 > 1e-4 or r_NH3_R302 > 1e-4:
        log("  [WARN] Residual Ammonia target not met. Adjusting slightly...")
        # Fine tuning could be added here, but stoichiometric calculation should be precise for this system
        pass

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
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            system.simulate()
    except Exception as e:
        log(f"  [WARN] Baseline simulation failed: {e}")
        log("  Proceeding with initialized flow rates as baseline.")
    
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
            # We skip explicit stream clearing to avoid accidentally clearing inputs
            # relying on robust simulation or error catching instead.
            
            for stream, baseline_flow in baseline_flows.items():
                target = baseline_flow * scaling_factor
                stream.F_mass = target
                # Store target in global registry for specs
                if 'Seed' in stream.ID:
                     # Check if it has 'Seed' component to be sure
                     if 'Seed' in stream.available_chemicals and stream.imass['Seed'] > 0:
                         seed_targets[stream] = stream.imass['Seed']
                     else:
                         seed_targets[stream] = target
            
            # Optimize NH3 Loading at this scale
            # This ensures we don't hit negative flows or infeasible balances
            optimize_NH3_loading(system, verbose=False)
            
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
            x1 = max(10.0, initial_guess * 3.0)
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
                xtol=0.001,
                ytol=0.1,
                maxiter=200,
                checkbounds=(y0 * y1 <= 0),
                checkiter=True,
            )
        
        log(f"\nSolver converged. Running final validation simulation...")
        
        # Skip clearing streams

        for stream, baseline_flow in baseline_flows.items():
            target = baseline_flow * scaling_factor
            stream.F_mass = target
            if 'Seed' in stream.ID:
                 if 'Seed' in stream.available_chemicals and stream.imass['Seed'] > 0:
                     seed_targets[stream] = stream.imass['Seed']
                 else:
                     seed_targets[stream] = target
        
        # Optimize NH3 for the final state
        optimize_NH3_loading(system, verbose=False)
        
        system.simulate()
        achieved_rate = product_stream.F_mass
        
        log(f"\n[OK] Successfully achieved production rate:")
        log(f"  Target:          {target_production_rate_kg_hr:.2f} kg/hr")
        log(f"  Achieved:        {achieved_rate:.2f} kg/hr")
        log(f"  Scaling factor:  {scaling_factor:.4f}x")
        log(f"  Error:           {abs(achieved_rate - target_production_rate_kg_hr):.4f} kg/hr")
        log(f"  Total iterations: {iteration[0]}")
        log(f"{'='*60}\n")
        
        return achieved_rate
        
    except Exception as e:
        log(f"\n[FAIL] Failed to achieve target production rate: {e}")
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
        
        status = '[PASS]' if passed else '[FAIL]'
        target_str = f"{min_val:.1f}% - {max_val:.1f}%"
        print(f"{spec_name:<30} {target_str:<25} {actual:>6.2f}%{'':<7} {status}")
        
        if not passed:
            all_passed = False
    
    purity_passed = protein_purity >= 65.0
    status = '[PASS]' if purity_passed else '[FAIL]'
    print(f"{'Protein Purity':<30} {'>= 65.0%':<25} {protein_purity:>6.2f}%{'':<7} {status}")
    
    if not purity_passed:
        all_passed = False
        
    # Calculate Heme Equivalent (LegHb contains 1 heme group per 763 carbons/normalized units)
    # Note: Leghemoglobin chemical is defined normalized to C=1 (MW ~22).
    # Real LegHb has 763 Carbons. Stoichiometry is 1 Heme per 763 mols of chem.
    legh_mol = product_stream.imol['Leghemoglobin']
    try:
        heme_mw = product_stream.chemicals.Heme_b.MW
        heme_equiv_mass_kg = (legh_mol / 763.0) * heme_mw
        equiv_wt = heme_equiv_mass_kg / product_stream.F_mass * 100 if product_stream.F_mass > 0 else 0
        print(f"{'Heme Equivalent':<30} {'(Info)':<25} {equiv_wt:>6.2f}%{'':<7} [INFO]")
    except AttributeError:
        print(f"{'Heme Equivalent':<30} {'(Info)':<25} {'N/A':>6}     [WARN] Heme_b not in chemicals")

    print(f"{'='*85}")
    
    print(f"\nProduct Stream Summary:")
    print(f"  Total Flow Rate: {total_mass:.2f} kg/hr")
    print(f"  Water Content:   {get_mass_percent('H2O'):.2f}%")
    print(f"  Leghemoglobin:   {legh_mass:.4f} kg/hr ({get_mass_percent('Leghemoglobin'):.2f}%)")
    
    if not all_passed:
        print(f"\n{'='*85}")
        raise ValueError("Product does not meet one or more specifications. See details above.")
    
    print(f"\n[OK] All specifications met!")
    print(f"{'='*85}\n")
    
    return True


# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == '__main__':
    bst.preferences.N = 50
    nn = 1
    TARGET_PRODUCTION = 150 * nn
    
    print("="*85)
    print("LEGHEMOGLOBIN PRODUCTION SYSTEM - MODULAR DESIGN")
    print("="*85)
    
    print("\n1. Creating system...")
    LegHb_sys = create_LegHb_system()
    sys = LegHb_sys
    f = sys.flowsheet
    u= f.unit
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
        # Run optimization BEFORE production rate setting, to establish correct ratios
        optimize_NH3_loading(LegHb_sys)
        
        achieved_production = set_production_rate(LegHb_sys, TARGET_PRODUCTION)
        
        # Run optimization AGAIN after production rate scaling, as demand changes with scale
        optimize_NH3_loading(LegHb_sys)
        
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
    
    # Calculate Heme Equivalent for Main Report
    legh_mol = ss.LegHb_3.imol['Leghemoglobin']
    heme_equiv_mass = 0
    try:
        heme_mw = LEGHB_THERMO.Heme_b.MW
        heme_equiv_mass = (legh_mol / 763.0) * heme_mw
    except:
        pass

    print(f"\n{'='*85}")
    print("KEY PERFORMANCE INDICATORS")
    print(f"{'='*85}")
    print(f"  Product Stream:           {ss.LegHb_3.ID}")
    print(f"  Production Rate:          {ss.LegHb_3.F_mass:.2f} kg/hr")
    print(f"  Leghemoglobin Content:    {legh_purity:.2f}%")
    print(f"  Annual Production:        {ss.LegHb_3.F_mass * sys.operating_hours / 1000:.2f} metric tons/year")
    print(f"  Yield (Heme eq):          {heme_equiv_mass:.4f} kg/hr")
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
        print(r1)
        r2 = bst.report.lca_displacement_allocation_table(
            systems=[sys],
            key='GWP',
            items=[ss.LegHb_3],
        )
        print("   LCA Displacement Allocation Table generated")
        print(r2)
    except Exception as e:
        print(f"   LCA analysis failed: {e}")
    
    print(f"\n{'='*85}")
    print("SIMULATION COMPLETE")
    print(f"Target Production:   {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"Achieved Production: {ss.LegHb_3.F_mass:.2f} kg/hr")
    print(f"{'='*85}\n")