# -*- coding: utf-8 -*-
"""
LegHb Production System - Config 1 NEW (Internalized Specifications)

Refactored from _config1.py to embed titer control (adjust_glucose_for_titer)
and NH3 optimization (optimize_NH3_loading) as unit add_specifications,
eliminating the need for external iterative convergence loops.

Key changes vs _config1.py:
    - adjust_glucose_for_titer logic → R302 add_specification
    - optimize_NH3_loading logic → S202 add_specification (analytical)
    - set_production_rate simplified (no internal optimize_NH3_loading calls)
    - No run_titer_convergence needed; system.simulate() handles everything

Process Areas (unchanged):
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

from biorefineries.prefers.v2._process_settings import set_GWPCF, GWP_CFs, set_GWPCF_Multi, load_process_settings
from biorefineries.prefers.v2.LegHb import _streams as s
import biosteam as bst
from thermosteam import Stream
from biosteam import F
import thermosteam as tmo
import numpy as np
from biorefineries.prefers.v2.LegHb import _chemicals as c
from biorefineries.prefers.v2 import _units as u
from biorefineries.prefers.v2 import _process_settings
from biorefineries.prefers.v2._process_settings import price

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
    # Initialize streams to reference values
    try:
        SeedIn1.imass['Seed'] = s.SeedSolution1['Seed']
        SeedIn2.imass['Seed'] = s.SeedSolution2['Seed']
        
        tmp = bst.Stream(ID='tmp_init', **s.SeedSolution2)
        ref_seed2_total = tmp.F_mass
        ratio_culture = (0.1 + 60 + 0.15191) / 1000
        CultureIn.imass['Culture'] = ref_seed2_total * ratio_culture
        tmp.empty()
    except:
        pass

    # Seed solution 1 preparation
    M201 = bst.MixTank('M201', ins=[SeedIn1, 'Water1'], outs='M201Out', tau=16)
    
    # Pre-calculate reference values for Seed1
    ref_stream1 = bst.Stream(**{**s.SeedSolution1, 'ID': 'RefSeed1'})
    ref_seed1_val = ref_stream1.imass['Seed']
    ref_water1_val = ref_stream1.imass['H2O']
    
    @M201.add_specification(run=True)
    def update_seed1_inputs():
        target = seed_targets.get(SeedIn1, ref_seed1_val)
        SeedIn1.imass['Seed'] = target
        if ref_seed1_val > 0:
            ratio = ref_water1_val / ref_seed1_val
            M201.ins[1].imass['H2O'] = target * ratio
        M201.ins[1].T = 25 + 273.15
    
    # Seed solution 2 + culture media preparation
    M202 = bst.MixTank('M202', ins=[SeedIn2, CultureIn, 'Water2'], outs='M202Out', tau=16)
    
    # Pre-calculate reference values for Seed2
    ref_stream2 = bst.Stream(**{**s.SeedSolution2, 'ID': 'RefSeed2'})
    ref_seed2_val = ref_stream2.imass['Seed']
    ref_water2_val = ref_stream2.imass['H2O']
    ref_total2_val = ref_stream2.F_mass
    
    @M202.add_specification(run=True)
    def update_culture_inputs():
        target = seed_targets.get(SeedIn2, ref_seed2_val)
        SeedIn2.imass['Seed'] = target
        if ref_seed2_val > 0:
            ratio_water = ref_water2_val / ref_seed2_val
            M202.ins[2].imass['H2O'] = target * ratio_water
            ratio_culture = (0.1 + 60 + 0.15191) / 1000
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
# AREA 300: CONVERSION (FERMENTATION) - WITH INTERNALIZED SPECS
# =============================================================================
def create_area_300_conversion(seed_in, glucose_in, ammonia_in, vent1, vent2, 
                                reactions, params=None, system_glucose_stream=None,
                                NH3_source_stream=None, S202_unit=None):
    """
    Create Area 300: Conversion (Fermentation)
    
    Purpose: Convert glucose to LegHb through microbial fermentation.
    
    NEW: Includes internalized specifications for:
        - Titer control (yield adjustment based on target titer) - replaces adjust_glucose_for_titer
        - NH3 optimization (analytical stoichiometric calculation) - replaces optimize_NH3_loading
    
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
    system_glucose_stream : Stream, optional
        Reference to the main system Glucose input stream for titer control
    NH3_source_stream : Stream, optional
        Reference to NH3_25wt input stream for NH3 optimization
    S202_unit : bst.Splitter, optional
        Reference to S202 splitter for NH3 split adjustment
    
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
        titer=params['titer_LegHb'],
        titer_IDs=('Leghemoglobin_In', 'Leghemoglobin'),
    )
    R302.target_titer = params['titer_LegHb']
    R302.target_productivity = params['productivity_LegHb']
    R302.target_yield = params['yield_LegHb']
    
    # Store reference streams for internalized specs
    R302._glucose_feed_stream = glucose_in  # T201Out (diluted)
    R302._system_glucose_stream = system_glucose_stream  # Main system input (Glucose)
    R302._nh3_source_stream = NH3_source_stream  # NH3_25wt
    R302._s202_unit = S202_unit  # S202 splitter
    
    # Store baseline values
    if system_glucose_stream is not None:
        R302._baseline_system_glucose = system_glucose_stream.imass['Glucose']
    R302._baseline_titer = params['titer_LegHb']
    R302._baseline_yield = params['yield_LegHb']
    
    # Elasticity factor for titer→yield correction (matches old adjust_glucose_for_titer)
    ELASTICITY = 1.3
    MAX_TITER_ITERATIONS = 15
    TITER_TOL = 0.001  # 0.1% relative tolerance
    
    # =========================================================================
    # INTERNALIZED SPECIFICATION: Titer Control + NH3 Optimization
    # (Replaces adjust_glucose_for_titer AND optimize_NH3_loading)
    # =========================================================================
    @R302.add_specification(run=False)
    def update_fermentation_and_nh3():
        """
        Combined specification that:
        1. Empirically measures NH3 demand (run R302 with excess, measure consumption)
           → sets correct NH3 in R302's ammonia inlet
           (replaces optimize_NH3_loading)
        2. Iteratively adjusts fermentation yield to hit target titer
           → uses elasticity-based correction with adaptive damping
           (replaces adjust_glucose_for_titer + run_titer_convergence)
        
        Each iteration: set yield → measure NH3 with excess → set exact NH3 →
        run R302 → check titer → adjust yield → repeat
        """
        target_titer = R302.titer if R302.titer is not None else R302.target_titer
        
        if target_titer is None or target_titer <= 0:
            R302._run()
            return
        
        baseline_titer = R302._baseline_titer
        baseline_yield = R302._baseline_yield
        
        if baseline_titer is None or baseline_titer <= 0 or baseline_yield is None or baseline_yield <= 0:
            R302._run()
            return
        
        # Use R302.target_yield as starting point if set (respects model parameters).
        # Fall back to _baseline_yield only if target_yield is not set.
        # This matches old adjust_glucose_for_titer which uses R302.target_yield.
        starting_yield = R302.target_yield
        if starting_yield is None or starting_yield <= 0:
            starting_yield = baseline_yield
        
        # Update tau based on titer/productivity relationship
        if R302.target_productivity and R302.target_productivity > 0:
            R302.tau = target_titer / R302.target_productivity
        
        # References
        nh3_source = R302._nh3_source_stream   # NH3_25wt (upstream of S202)
        s202 = R302._s202_unit                 # S202 splitter
        ammonia_in = R302.ins[2]               # NH3_Fer (R302's direct NH3 inlet)
        
        # Helper: get NH3 consumption of a reactor
        def get_nh3_consumption(reactor):
            nh3_in = sum(i.imol['NH3'] for i in reactor.ins)
            nh3_out = sum(o.imol['NH3'] for o in reactor.outs)
            return max(0, nh3_in - nh3_out)
        
        # Helper: set stream to supply target_nh3_mol kmol/hr NH3
        # with 25:75 NH3:H2O mass ratio (matching NH3_25wt composition)
        def set_nh3_stream(stream, target_nh3_mol):
            if target_nh3_mol <= 0:
                return
            nh3_mass = target_nh3_mol * 17.031  # kg/hr
            stream.imass['NH3'] = nh3_mass
            stream.imass['H2O'] = nh3_mass * 3.0  # 75/25 mass ratio
        
        # Helper: run R302 with excess NH3, return consumption in kmol/hr
        def measure_nh3_consumption(yield_val):
            """Run R302 with 20x excess NH3 to measure actual consumption."""
            R302.target_yield = yield_val
            fermentation_reaction[2].product_yield(
                'Leghemoglobin_In', basis='wt', product_yield=yield_val
            )
            # Save state and set excess
            saved_mol = ammonia_in.mol.copy()
            if ammonia_in.F_mol > 1e-10:
                ammonia_in.mol *= 20.0
            else:
                set_nh3_stream(ammonia_in, 50.0)
            R302._run()
            consumed = get_nh3_consumption(R302)
            # Restore (will be overwritten next)
            ammonia_in.mol = saved_mol
            return consumed
        
        # =====================================================================
        # Measure R301's NH3 consumption (R301 already ran upstream)
        # =====================================================================
        consumed_R301 = get_nh3_consumption(R301)
        
        # =====================================================================
        # Unified convergence loop: yield + NH3 interleaved
        # Each iteration: measure NH3 empirically → set exact → check titer → adjust yield
        #
        # Starting yield comes from R302.target_yield (which model parameters can
        # update), with elasticity correction applied if titer target differs from
        # the baseline titer. This matches old adjust_glucose_for_titer behavior
        # where current_yield = R302.target_yield and correction is incremental.
        # =====================================================================
        titer_ratio = target_titer / baseline_titer
        # Start from user-set yield with titer correction (incremental approach)
        if abs(titer_ratio - 1.0) > 1e-6:
            # Titer changed from baseline: apply elasticity correction from starting yield
            current_yield = starting_yield * titer_ratio ** (1.0 / ELASTICITY)
        else:
            # Titer at baseline: use starting yield directly
            current_yield = starting_yield
        current_yield = max(0.005, min(0.15, current_yield))
        consumed_R302 = 0.0
        
        for _iteration in range(MAX_TITER_ITERATIONS):
            # 1. Measure NH3 demand at current yield (run with excess)
            consumed_R302 = measure_nh3_consumption(current_yield)
            
            # 2. Set exact NH3 and run R302
            if consumed_R302 > 0:
                set_nh3_stream(ammonia_in, consumed_R302 * 1.001)
            R302.target_yield = current_yield
            fermentation_reaction[2].product_yield(
                'Leghemoglobin_In', basis='wt', product_yield=current_yield
            )
            R302._run()
            
            # 3. Check titer
            actual_titer = R302.actual_titer
            if not actual_titer or actual_titer <= 0:
                break
            
            error = abs(target_titer - actual_titer) / target_titer
            if error < TITER_TOL:
                break
            
            # 4. Adjust yield using elasticity correction
            correction_ratio = target_titer / actual_titer
            if error > 0.10:
                damping = 1.0   # Coarse mode
            else:
                damping = 0.5   # Fine mode
            step = (correction_ratio ** (1.0 / ELASTICITY)) ** damping
            current_yield *= step
            current_yield = max(0.005, min(0.15, current_yield))
        
        # =====================================================================
        # Update upstream NH3_25wt source and S202 split for consistency
        # (so future system.simulate() calls / set_production_rate use correct NH3)
        # =====================================================================
        if nh3_source is not None and s202 is not None:
            total_demand = consumed_R301 + consumed_R302
            
            if total_demand > 0:
                # Set NH3_25wt source to total demand * 1.001
                if nh3_source.F_mol > 1e-10:
                    nh3_mol_frac = nh3_source.imol['NH3'] / nh3_source.F_mol
                    required_total_nh3 = total_demand * 1.001
                    required_F_mol = required_total_nh3 / nh3_mol_frac
                    nh3_source.F_mass = nh3_source.F_mass * (required_F_mol / nh3_source.F_mol)
                else:
                    set_nh3_stream(nh3_source, total_demand * 1.001)
                
                # Set S202 split (outs[0] = NH3_Seed fraction)
                split_seed = consumed_R301 / total_demand
                s202.split[:] = split_seed
    
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
        split=c.chemical_groups['SolidsCentrifuge'],
        moisture_content=0.40,
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
        split=c.chemical_groups['SolidsCentrifuge'],
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
    C403 = u.Centrifuge(
        'C403',
        ins=H401-0,
        outs=('CellDebris', 'CrudeLysate'),
        split=c.chemical_groups['SolidsCentrifuge'],
        moisture_content=0.20,
    )
    
    # Depth filtration
    S403 = u.FiltrationAdv.from_preset(
        'MF',
        'S403',
        ins=C403-1,
        outs=('FilterCake', 'ClarifiedLysate'),
        solid_capture_efficiency=0.85,
        cake_moisture_content=0.30,
        solid_IDs=c.chemical_groups['FilteredSubstances'],
    )
    S403.add_specification(run=True)
    
    # Screw press for debris dewatering
    S404 = bst.ScrewPress(
        'S404',
        ins=(C403-0, S403-0),
        outs=('DehydratedDebris', 'PressLiquor'),
        split=0.99,
        moisture_content=0.01
    )
    
    return {
        'units': [C401, M401, H402, M402, C402, S401, H401, C403, S403, S404],
        'clarified_lysate': S403.outs[1],
        'spent_media': C401.outs[1],
        'wash_effluent': C402.outs[1],
        'dehydrated_debris': S404.outs[0],
        'press_liquor': S404.outs[1],
        'C401': C401, 'M401': M401, 'H402': H402,
        'M402': M402, 'C402': C402, 'S401': S401,
        'H401': H401, 'C403': C403, 'S403': S403, 'S404': S404,
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
    # DF buffer preparation
    M501 = bst.MixTank('M501', ins=(DfUltraBuffer1, 'Water5'), 
                       outs='DFBufferOut', tau=0.5)
    
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
    
    # UF Stage 1 - concentration + buffer exchange
    U501 = u.DiafiltrationAdv.from_preset(
        'UF',
        'U501',
        ins=(clarified_lysate, H501-0),
        outs=('UFConcentrate', 'UFPermeate'),
        TargetProduct_IDs=('Leghemoglobin',),
        Salt_IDs=tuple(c.chemical_groups['Salts']),
        OtherLargeMolecules_IDs=tuple(c.chemical_groups['OtherLargeMolecules']),
        TargetProduct_Retention=0.99,
        Salt_Retention=0.05,
    )
    U501.add_specification(run=True)
    
    # UF Stage 2 - Final concentration
    U502 = u.DiafiltrationAdv(
        'U502',
        ins=(U501-0, bst.Stream('Water6', H2O=0.001)),
        outs=('ConcentratedLegH', 'ConcentrationPermeate'),
        TargetProduct_IDs=('Leghemoglobin',),
        Salt_IDs=tuple(c.chemical_groups['Salts']),
        OtherLargeMolecules_IDs=tuple(c.chemical_groups['OtherLargeMolecules']),
        TMP_bar=3.0,
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
    AntioxidantStream : Stream
        Antioxidant input stream
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
    # Supplemental nutrient streams for wastewater treatment stability
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
        F.Water5, F.Water6, F.Water7,
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
    ins=[s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt, s.DfUltraBuffer1, s.DfUltraBuffer2, s.AntioxidantStream],
    outs=[s.LegHb_3, s.vent1, s.vent2, s.effluent1, s.emissions, s.ash_disposal],
    fthermo=lambda chemicals=None: LEGHB_THERMO,
)
def create_LegHb_system(ins, outs, use_area_convention=False):
    """
    Creates the LegHb (Leghemoglobin) production system.
    
    NEW: Titer control and NH3 optimization are embedded as unit specifications.
    No external adjust_glucose_for_titer or optimize_NH3_loading calls needed.
    
    Modular architecture using discrete Process Area functions:
        - Area 200: Media Preparation
        - Area 300: Conversion (Fermentation) [with internalized titer + NH3 specs]
        - Area 400: Recovery (Harvest & Cell Disruption)
        - Area 500: Purification (UF/DF)
        - Area 600: Formulation (Thermal & Product)
        - Area 900: Facilities & Utilities
    
    Parameters
    ----------
    ins : tuple
        Input streams (SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer1, DfUltraBuffer2, AntioxidantStream)
    outs : tuple
        Output streams (LegHb_3, vent1, vent2, effluent1, emissions, ash_disposal)
    use_area_convention : bool
        Whether to use area ID convention for utility units
    
    Returns
    -------
    tuple
        Output streams
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
    
    # Area 300: Conversion (Fermentation) - with internalized titer + NH3 specs
    area_300 = create_area_300_conversion(
        seed_in=area_200['seed_out'],
        glucose_in=area_200['glucose_out'],
        ammonia_in=area_200['ammonia_out'],
        vent1=vent1,
        vent2=vent2,
        reactions=reactions,
        params=params,
        system_glucose_stream=Glucose,
        NH3_source_stream=NH3_25wt,      # NEW: pass NH3 source for optimization
        S202_unit=area_200['S202'],       # NEW: pass S202 for split adjustment
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

def set_production_rate(system, target_production_rate_kg_hr, verbose=True):
    """
    Adjust system inputs to achieve target LegHb_3 production rate.
    
    Simplified compared to _config1.py:
    - Does NOT call optimize_NH3_loading (handled by R302 specification)
    - Uses flexsolve IQ_interpolation for scaling factor only
    
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
    
    log(f"\n{'='*60}")
    log(f"Setting Production Rate to {target_production_rate_kg_hr:.2f} kg/hr")
    log(f"{'='*60}")
    
    # Run initial simulation FIRST to establish proper system state.
    # R302 specification adjusts yield/NH3 during simulation, which changes
    # input stream flows. We must store baseline_flows AFTER simulation so
    # that IQ_interpolation scales from the correct post-spec state.
    log("Running initial simulation to establish baseline...")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            system.simulate()
    except Exception as e:
        log(f"  [WARN] Baseline simulation failed: {e}")
        log("  Proceeding with initialized flow rates as baseline.")
    
    # Store baseline flows AFTER simulation (post-spec state)
    baseline_flows = {}
    for stream in system.ins:
        if stream.F_mass > 0:
            baseline_flows[stream] = stream.F_mass
    
    if not baseline_flows:
        raise ValueError("No input streams with positive flow rates found")
    
    log(f"Baseline input streams stored: {len(baseline_flows)} streams")
    
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
                target = baseline_flow * scaling_factor
                stream.F_mass = target
                if 'Seed' in stream.ID:
                    if 'Seed' in stream.available_chemicals and stream.imass['Seed'] > 0:
                        seed_targets[stream] = stream.imass['Seed']
                    else:
                        seed_targets[stream] = target
            
            # NH3 optimization is handled automatically by R302 specification
            # No need for explicit optimize_NH3_loading call
            
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
            
            x0 = max(0.01, initial_guess * 0.3)
            x1 = max(initial_guess * 3.0, 5.0)
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
                maxiter=300,
                checkbounds=(y0 * y1 <= 0),
                checkiter=True,
            )
        
        log(f"\nSolver converged. Running final validation simulation...")
        
        for stream, baseline_flow in baseline_flows.items():
            target = baseline_flow * scaling_factor
            stream.F_mass = target
            if 'Seed' in stream.ID:
                if 'Seed' in stream.available_chemicals and stream.imass['Seed'] > 0:
                    seed_targets[stream] = stream.imass['Seed']
                else:
                    seed_targets[stream] = target
        
        # NH3 handled by specs - just simulate
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
        
    # Calculate Heme Equivalent
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
    print("LEGHEMOGLOBIN PRODUCTION SYSTEM - CONFIG1_NEW (Internalized Specs)")
    print("="*85)
    
    print("\n1. Creating system...")
    LegHb_sys = create_LegHb_system()
    sys = LegHb_sys
    f = sys.flowsheet
    u = f.unit
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
        target_titer = u.R302.titer if u.R302.titer is not None else u.R302.target_titer
        if target_titer is not None:
            u.R302.titer = float(target_titer)
            u.R302.target_titer = float(target_titer)
            
        print(f"   Target Titer: {target_titer} g/L")
        
        # Direct scaling - no run_titer_convergence needed
        # NH3 and titer control happen automatically via R302 specification
        set_production_rate(sys, TARGET_PRODUCTION, verbose=True)
        
        final_production = ss.LegHb_3.F_mass
        
        if abs(final_production - TARGET_PRODUCTION) > 1.0:
            print(f"\n   WARNING: Production rate drifted after final simulation!")
            print(f"   Target:  {TARGET_PRODUCTION:.2f} kg/hr")
            print(f"   Actual:  {final_production:.2f} kg/hr")
            
    except Exception as e:
        print(f"\n   Could not achieve target production: {e}")
        print("   Continuing with baseline production rate...")
        import traceback; traceback.print_exc()
    
    print(f"\n4. Verifying product specifications...")
    try:
        check_LegHb_specifications(ss.LegHb_3)
    except ValueError as e:
        print(f"\n   SPECIFICATION CHECK FAILED: {e}")
    
    print(f"\n5. System Summary")
    print("="*85)
    LegHb_sys.show()
    
    legh_purity = ss.LegHb_3.imass['Leghemoglobin'] / ss.LegHb_3.F_mass * 100
    
    # Calculate Heme Equivalent
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
    
    print(f"\n6. Broth Analysis")
    print(f"   target titer: {u.R302.titer:.2f} g/L")
    print(f"   Achieved titer: {u.R302.actual_titer:.2f} g/L")
    print(f"   Achieved titer2: {(ss.Broth.imass['Leghemoglobin']+ss.Broth.imass['Leghemoglobin_In']) / ss.Broth.F_vol:.2f} g/L")
    print(ss.Broth)
    
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
