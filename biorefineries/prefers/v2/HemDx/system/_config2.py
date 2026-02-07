# -*- coding: utf-8 -*-
"""
N-HemDx Full Production System - Config 2 NEW (Internalized Specifications)

Refactored from _config2.py to embed titer control (adjust_glucose_for_titer)
and NH3 optimization (optimize_NH3_loading) as unit add_specifications,
eliminating the need for external iterative convergence loops.

Config 2: Intracellular production (SF=0.1).
Supernatant filtration (S404) removed; supernatant goes directly to WWT via M504.

Key changes vs _config2.py:
    - adjust_glucose_for_titer logic → R302 add_specification
    - optimize_NH3_loading logic → R302 add_specification (empirical)
    - set_production_rate simplified (no internal optimize_NH3_loading calls)
    - No run_titer_convergence needed; system.simulate() handles everything

Process Areas:
    - Area 200: Media Preparation
    - Area 300: Conversion (Fermentation)
    - Area 400: Recovery (Clarification & Cell Disruption)
    - Area 500: Purification (Resin Adsorption)
    - Area 600: Concentration (NF/Diafiltration)
    - Area 700: Formulation (Complexation)
    - Area 800: Final Product
    - Area 900: Facilities & Utilities

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biorefineries.prefers.v2._process_settings import set_GWPCF, GWP_CFs, set_GWPCF_Multi, load_process_settings
from biorefineries.prefers.v2.HemDx import _streams as s
import biosteam as bst
from thermosteam import Stream
from biosteam import F
import thermosteam as tmo
import numpy as np
from biorefineries.prefers.v2.HemDx import _chemicals as c
from biorefineries.prefers.v2 import _units as u
from biorefineries.prefers.v2 import _process_settings
from biorefineries.prefers.v2.HemDx import _streams as s
import numpy as np

# Global registry for seed scaling targets
seed_targets = {}

# =============================================================================
# AREA 100: FEEDING
# =============================================================================
HEMDX_THERMO = c.create_chemicals_Hemodextrin()
bst.settings.set_thermo(HEMDX_THERMO, skip_checks=True)
bst.preferences.classic_mode()
try:
    bst.settings.chemicals.set_alias('Heme_b', 'Heme_b_In')
except ValueError:
    pass

__all__ = (
    'create_NHemDx_system',
    'create_area_200_media_prep',
    'create_area_300_conversion',
    'create_area_400_recovery',
    'create_area_500_purification',
    'create_area_600_concentration',
    'create_area_700_formulation',
    'create_area_900_facilities',
    'set_production_rate',
    'check_HemDx_specifications',
)

# =============================================================================
# PARAMETERS & REACTIONS
# =============================================================================
def get_fermentation_parameters():
    """Return fermentation parameters."""
    # Fermentation Parameters
    return {
        'theta_O2': 0.45,
        'agitation_power': 0.985,
        'design': 'Stirred tank',
        'method': "Riet",
        'T_operation': 273.15 + 30,
        'Q_O2_consumption': -110 * 4184,
        'dT_hx_loop': 8,
        'cooler_pressure_drop': 20684,
        'compressor_isentropic_efficiency': 0.85,
        'V_max': 500,
        # Kinetics
        'titer_Heme': 0.380,  # g/L
        'tau1': 72,  # hours growth
        'tau2': 90,  # hours production
        'P_Heme_1st': 1e-8,
        'P_pp_1st': 1e-8,
        'P_Heme_2nd': 4.2/1e3,
        'P_pp_2nd': 0.7/1e3,
        'SF': 0.1,  # Secretion fraction (Intracellular)
    }

def create_fermentation_reactions(params=None):
    if params is None: params = get_fermentation_parameters()
    
    tau1 = params['tau1']
    tau2 = params['tau2']
    tau = tau1 + tau2
    
    P_ComsuptionG1 = 500 * (0.5/3) / tau1
    P_cell_growth_rate1 = 150 * 0.35 / tau1
    Y_Heme1 = params['P_Heme_1st'] / P_ComsuptionG1
    Y_pp1 = params['P_pp_1st'] / P_ComsuptionG1
    Y_b1 = P_cell_growth_rate1 / P_ComsuptionG1

    P_ComsuptionG2 = 800 * (0.5/3) / tau2
    P_cell_growth_rate2 = 50 * 0.35 / tau2
    Y_Heme2 = params['P_Heme_2nd'] / P_ComsuptionG2
    Y_pp2 = params['P_pp_2nd'] / P_ComsuptionG2
    Y_b2 = P_cell_growth_rate2 / P_ComsuptionG2

    # Weighted averages
    Y_Heme = (params['P_Heme_1st']*tau1 + params['P_Heme_2nd']*tau2) / (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2)
    Y_pp = (params['P_pp_1st']*tau1 + params['P_pp_2nd']*tau2) / (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2)
    Y_b = (P_cell_growth_rate1*tau1 + P_cell_growth_rate2*tau2) / (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2)
    
    SF = params['SF']
    
    fermentation_reaction = bst.PRxn([
        bst.Rxn('1 Glucose + 0.70588 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b + 0.617647 O2 + 0.17647 H2SO4 + 4.05882 H2O',
                reactant='Glucose', X=Y_Heme*SF, check_atomic_balance=True),
        bst.Rxn('1 Glucose + 0.70588 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b_In + 0.617647 O2 + 0.17647 H2SO4 + 4.05882 H2O',
                reactant='Glucose', X=Y_Heme*(1-SF), check_atomic_balance=True),
        bst.Rxn('1 Glucose + 0.70588 NH3  -> 0.17647 ProtoporphyrinIX + 0.618 O2 + 4.06 H2O',
                reactant='Glucose', X=Y_pp*SF, correct_atomic_balance=True),
        bst.Rxn('1 Glucose + 0.70588 NH3 -> 0.17647 ProtoporphyrinIX_In + 0.618 O2 + 4.06 H2O',
                reactant='Glucose', X=Y_pp*(1-SF), correct_atomic_balance=True),
    ]) 
    fermentation_reaction[0].product_yield('Heme_b', basis='wt', product_yield=Y_Heme*SF)
    fermentation_reaction[1].product_yield('Heme_b_In', basis='wt', product_yield=Y_Heme*(1-SF))
    fermentation_reaction[2].product_yield('ProtoporphyrinIX', basis='wt', product_yield=Y_pp*SF)
    fermentation_reaction[3].product_yield('ProtoporphyrinIX_In', basis='wt', product_yield=Y_pp*(1-SF))

    neutralization_reaction = bst.Rxn(
        'H2SO4 + 2 NH3 -> (NH4)2SO4', reactant='H2SO4', X=0.99,
        check_atomic_balance=True
    )
    
    cell_growth_reactionCG1 = bst.Rxn(
        'Glucose + 1.32 NH3 + 1.32 H2O-> 6 Corynebacterium_glutamicum + 2.37 O2', 'Glucose', X=Y_b*2,
        correct_atomic_balance=True
    )
    cell_growth_reactionCG1.product_yield('Corynebacterium_glutamicum', basis='wt', product_yield=Y_b*2)
    
    respiration_reactionGC1 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose',X = 1 - Y_b*2,
        check_atomic_balance=True
    )
    
    cell_growth_reactionCG2 = bst.Rxn(
        'Glucose + 1.32 NH3 + 1.32 H2O-> 6 Corynebacterium_glutamicum + 2.37 O2', 'Glucose', X=Y_b*2,
        correct_atomic_balance=True
    )
    cell_growth_reactionCG2.product_yield('Corynebacterium_glutamicum', basis='wt', product_yield=Y_b*2)
    respiration_reactionGC2 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', (1 - Y_b - Y_pp - Y_Heme)*0.95,
        check_atomic_balance=True
    )
    
    respiration_reactionGC3 = bst.Rxn(
        'Glucose -> 2 EtOH + 2 CO2', 'Glucose', (1 - Y_b - Y_pp - Y_Heme)*0.05,
        check_atomic_balance=True
    )
    
    RXN = bst.ReactionSystem(
        fermentation_reaction,
        bst.PRxn([cell_growth_reactionCG2, respiration_reactionGC2, respiration_reactionGC3]),
    )
    
    return {
        'fermentation_reaction': fermentation_reaction,
        'cell_growth_reactionCG1': cell_growth_reactionCG1,
        'cell_growth_reactionCG2': cell_growth_reactionCG2,
        'respiration_reactionGC1': respiration_reactionGC1,
        'respiration_reactionGC2': respiration_reactionGC2,
        'respiration_reactionGC3': respiration_reactionGC3,
        'neutralization_reaction': neutralization_reaction,
        'RXN': RXN,
        'Y_Heme': Y_Heme,
        'Y_pp': Y_pp
    }

# =============================================================================
# AREA 200: MEDIA PREPARATION
# =============================================================================
def create_area_200_media_prep(SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt):
    # Initialize streams to reference values (run once during creation)
    try:
        SeedIn1.imass['Seed'] = s.SeedSolution1['Seed']
        SeedIn2.imass['Seed'] = s.SeedSolution2['Seed']
        
        tmp = bst.Stream(ID='tmp_init_h', **s.SeedSolution2)
        ref_seed2_total = tmp.F_mass
        ratio_culture = (0.1 + 60 + 0.15191) / 1000
        CultureIn.imass['Culture'] = ref_seed2_total * ratio_culture
        tmp.empty()
    except:
        pass

    M201 = bst.MixTank('M201', ins=[SeedIn1,'Water1'], outs='M201Out', tau=16)
    # Pre-calculate reference values for Seed1 to avoid re-registering stream in loop
    ref_stream1 = bst.Stream(**{**s.SeedSolution1, 'ID': 'RefSeed1'})
    ref_seed1_val = ref_stream1.imass['Seed']
    ref_water1_val = ref_stream1.imass['H2O']

    @M201.add_specification(run=True)
    def update_seed1_inputs():
        # Force SeedIn1 to target (from global dict or reference)
        target = seed_targets.get(SeedIn1, ref_seed1_val)
        
        SeedIn1.imass['Seed'] = target
        
        if ref_seed1_val > 0:
            ratio = ref_water1_val / ref_seed1_val
            M201.ins[1].imass['H2O'] = target * ratio
        
        M201.ins[1].T = 25+273.15
    
    M202 = bst.MixTank('M202', ins=[SeedIn2,CultureIn,'Water2'], outs='M202Out', tau=16)
    # Pre-calculate reference values for Seed2 to avoid re-registering stream in loop
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
            M202.ins[2].imass['H2O'] = target * ratio_water # Water2
            
            # CultureIn logic
            ratio_culture = (0.1 + 60 + 0.15191) / 1000
            
            # Predict current total solution mass
            current_total_mass_est = target * (ref_total2_val / ref_seed2_val)
            CultureIn.imass['Culture'] = current_total_mass_est * ratio_culture
        
        M202.ins[2].T = 25+273.15

    T202 = u.AmmoniaStorageTank('T202', ins=NH3_25wt, outs='T202Out')
    S202 = bst.Splitter('S202', ins=T202-0, outs=('NH3_Seed', 'NH3_Fer'), split=0.04)

    M203 = u.SeedHoldTank('M203', ins=[M201-0, M202-0, S202-0], outs='M203Out')
    
    M204 = bst.MixTank('M204', ins=[Glucose,'Water3'], outs='M204Out', tau=16)
    @M204.add_specification(run=True)
    def update_water_content():
        M204.ins[1].imass['H2O'] = Glucose.imass['Glucose']/2
        M204.ins[1].T = 25+273.15
    
    T201 = bst.StorageTank('T201', ins=M204-0, outs='T201Out', tau=16*4+72)
    
    return {
        'M203_out': M203-0,
        'T201_out': T201-0,
        'T202_out': S202-1,
        'S202': S202,
    }

# =============================================================================
# AREA 300: CONVERSION
# =============================================================================
def create_area_300_conversion(seed_in, glucose_in, ammonia_in, vent1, vent2, rxns, params, NH3_source_stream=None, S202_unit=None):
    R301 = u.SeedTrain(
        'R301',
        ins=[seed_in],
        outs=[vent1, 'R301Out'],
        reactions=bst.PRxn([rxns['cell_growth_reactionCG1'], rxns['respiration_reactionGC1'],]),
        saccharification=None,
        T=32+273.15,
    )
    R301.add_specification(run=True)
    
    R302 = u.AeratedFermentation(
        'R302',
        ins=[R301-1, glucose_in, ammonia_in, bst.Stream('FilteredAir', phase='g', P=2 * 101325)],
        outs=[vent2, 'Broth'],
        fermentation_reaction=rxns['fermentation_reaction'],
        cell_growth_reaction=rxns['cell_growth_reactionCG2'],
        respiration_reaction=bst.PRxn([rxns['respiration_reactionGC2'], rxns['respiration_reactionGC3']]),
        neutralization_reaction=rxns['neutralization_reaction'],
        design=params['design'], method=params['method'], theta_O2=params['theta_O2'],
        V_max=params['V_max'], Q_O2_consumption=params['Q_O2_consumption'],
        dT_hx_loop=params['dT_hx_loop'], T=params['T_operation'],
        batch=True, reactions=rxns['RXN'],
        kW_per_m3=params['agitation_power'],
        tau=params['tau1']+params['tau2'],
        cooler_pressure_drop=params['cooler_pressure_drop'],
        compressor_isentropic_efficiency=params['compressor_isentropic_efficiency'],
        P=1 * 101325,
        # Titer control parameters
        titer=None,  # Will be set below
        titer_IDs=('Heme_b', 'Heme_b_In'),  # HemDx products (secreted and intracellular)
    )
    # Expose parameters for Model analysis
    R302.reaction_params = params
    R302.fermentation_rxns_collection = rxns
    
    # Calculate target yields for R302 specification
    P_Heme = (params['P_Heme_1st']*params['tau1'] + params['P_Heme_2nd']*params['tau2']) / (params['tau1']+params['tau2'])
    Y_Heme = rxns['Y_Heme']
    Y_pp = rxns['Y_pp']
    SF = params['SF']
    
    # Set titer control attributes
    R302.target_titer = params['titer_Heme']
    R302.titer = R302.target_titer
    R302.target_productivity = P_Heme
    R302.target_yield = Y_Heme
    
    # Store reference streams for internalized specs
    R302._nh3_source_stream = NH3_source_stream  # NH3_25wt (upstream of S202)
    R302._s202_unit = S202_unit                  # S202 splitter
    
    # Store baseline values
    R302._baseline_titer = params['titer_Heme']
    R302._baseline_yield = Y_Heme
    
    # Elasticity factor for titer→yield correction
    ELASTICITY = 2.0    # Higher than LegHb (1.3) to dampen HemDx oscillation
    MAX_TITER_ITERATIONS = 15
    TITER_TOL = 0.005   # 0.5% relative tolerance
    
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
        2. Iteratively adjusts fermentation yield to hit target titer
           → uses elasticity-based correction with adaptive damping
        3. Updates all 4 fermentation reaction yields using SF split
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
        
        # Use R302.target_yield as starting point (respects model parameters)
        starting_yield = R302.target_yield
        if starting_yield is None or starting_yield <= 0:
            starting_yield = baseline_yield
        
        # Update tau based on titer/productivity relationship
        if R302.target_productivity and R302.target_productivity > 0:
            R302.tau = target_titer / R302.target_productivity
        
        # Get current SF from reaction_params (allows model to change SF dynamically)
        current_SF = R302.reaction_params.get('SF', SF)
        
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
        
        # Helper: update all 4 fermentation reaction yields using SF split
        def update_yields(yield_val):
            R302.target_yield = yield_val
            rxns['fermentation_reaction'][0].product_yield('Heme_b', basis='wt', product_yield=yield_val * current_SF)
            rxns['fermentation_reaction'][1].product_yield('Heme_b_In', basis='wt', product_yield=yield_val * (1 - current_SF))
            rxns['fermentation_reaction'][2].product_yield('ProtoporphyrinIX', basis='wt', product_yield=Y_pp * current_SF)
            rxns['fermentation_reaction'][3].product_yield('ProtoporphyrinIX_In', basis='wt', product_yield=Y_pp * (1 - current_SF))
        
        # Helper: run R302 with excess NH3, return consumption in kmol/hr
        def measure_nh3_consumption(yield_val):
            update_yields(yield_val)
            saved_mol = ammonia_in.mol.copy()
            if ammonia_in.F_mol > 1e-10:
                ammonia_in.mol *= 20.0
            else:
                set_nh3_stream(ammonia_in, 50.0)
            R302._run()
            consumed = get_nh3_consumption(R302)
            ammonia_in.mol = saved_mol
            return consumed
        
        # Measure R301's NH3 consumption (R301 already ran upstream)
        consumed_R301 = get_nh3_consumption(R301)
        
        # Starting yield with titer correction
        titer_ratio = target_titer / baseline_titer
        if abs(titer_ratio - 1.0) > 1e-6:
            current_yield = starting_yield * titer_ratio ** (1.0 / ELASTICITY)
        else:
            current_yield = starting_yield
        current_yield = max(0.0001, min(0.05, current_yield))
        consumed_R302 = 0.0
        
        for _iteration in range(MAX_TITER_ITERATIONS):
            # 1. Measure NH3 demand at current yield (run with excess)
            consumed_R302 = measure_nh3_consumption(current_yield)
            
            # 2. Set exact NH3 and run R302
            if consumed_R302 > 0:
                set_nh3_stream(ammonia_in, consumed_R302 * 1.001)
            update_yields(current_yield)
            R302._run()
            
            # 3. Check titer
            actual_titer = R302.actual_titer
            if not actual_titer or actual_titer <= 0:
                break
            
            error = abs(target_titer - actual_titer) / target_titer
            if error < TITER_TOL:
                break
            
            # 4. Adjust yield using elasticity correction with adaptive damping
            correction_ratio = target_titer / actual_titer
            if error > 0.10:
                damping = 1.0   # Coarse mode
            else:
                damping = 0.5   # Fine mode
            step = (correction_ratio ** (1.0 / ELASTICITY)) ** damping
            current_yield *= step
            current_yield = max(0.0001, min(0.05, current_yield))
        
        # Update upstream NH3_25wt source and S202 split
        if nh3_source is not None and s202 is not None:
            total_demand = consumed_R301 + consumed_R302
            
            if total_demand > 0:
                if nh3_source.F_mol > 1e-10:
                    nh3_mol_frac = nh3_source.imol['NH3'] / nh3_source.F_mol
                    required_total_nh3 = total_demand * 1.001
                    required_F_mol = required_total_nh3 / nh3_mol_frac
                    nh3_source.F_mass = nh3_source.F_mass * (required_F_mol / nh3_source.F_mol)
                else:
                    set_nh3_stream(nh3_source, total_demand * 1.001)
                
                split_seed = consumed_R301 / total_demand
                s202.split[:] = split_seed
    
    return R302-1

# =============================================================================
# AREA 400: RECOVERY (Config 2 - Intracellular)
# =============================================================================
def create_area_400_recovery(broth_in, DfUltraBuffer2):
    # C401: Primary Centrifuge (simplified using chemical_groups)
    C401 = u.Centrifuge(
        'C401',
        ins=broth_in,
        outs=('CellCream', 'Supernatant'),
        split=c.chemical_groups['SolidsCentrifuge'],
        moisture_content=0.40,
    )

    # Config 2: S404 removed - supernatant goes directly to M504 for WWT
    M504 = bst.Mixer('M504', ins=C401-1, outs='SupernatantToWWT')

    # Cell cream washing & disruption
    M401 = bst.MixTank('M401', ins=(DfUltraBuffer2, 'Water4'), outs='WashBufferOut', tau=0.5)

    @M401.add_specification(run=True)
    def update_wash_buffer():
        M401.ins[1].imass['H2O'] = C401.outs[0].F_mass * 1.5
        M401.ins[0].imol['DfUltraBuffer'] = M401.ins[1].imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000

    H402 = bst.HXutility('H402', ins=M401-0, outs='ColdWashBuffer', T=10 + 273.15, cool_only=True)
    M402 = bst.MixTank('M402', ins=(C401-0, H402-0), outs='WashedCellSlurry', tau=0.25)

    C402 = u.Centrifuge(
        'C402',
        ins=M402-0,
        outs=('WashedCellCream', 'WashEffluent'),
        split=c.chemical_groups['SolidsCentrifuge'],
        moisture_content=0.55,
    )

    S402 = u.CellDisruption(
        'S402', ins=C402-0, outs='CrudeHomogenate',
        Cell_ID='Corynebacterium_glutamicum',
        cell_disruption_efficiency=0.55, P_high=1000e5, P_low=101325,
        component_fractions={'Protein': 0.45, 'Cellulose': 0.22, 'Xylan': 0.15, 'OleicAcid': 0.08, 'RNA': 0.10}
    )

    H401 = bst.HXutility('H401', ins=S402-0, outs='CooledHomogenate', T=15 + 273.15, cool_only=True)

    C403 = u.Centrifuge(
        'C403', ins=H401-0, outs=('CellDebrisSolids', 'CrudeLysate'),
        split=c.chemical_groups['SolidsCentrifuge'],
        moisture_content=0.20,
    )

    S405 = u.FiltrationAdv.from_preset(
        'MF', 'S405', ins=C403-1, outs=('FilterCake', 'ClarifiedLysate'),
        solid_capture_efficiency=0.85, cake_moisture_content=0.30,
        solid_IDs=c.chemical_groups['FilteredSubstances'],
    )
    S405.add_specification(run=True)

    # Debris Dewatering (Config 2: no S404-0, only C403-0 and S405-0)
    M404 = bst.Mixer('M404', ins=(C403-0, S405-0), outs='CellDebrisRaw')
    S406 = bst.ScrewPress(
        'S406', ins=M404-0, outs=('DehydratedDebris', 'PressLiquor'),
        split=0.999, moisture_content=0.001
    )

    # Config 2: M405 only receives S405-1 (no S404-1)
    M405 = bst.Mixer('M405', ins=S405-1, outs='CombinedLysate')
    
    return {
        'CombinedLysate': M405-0,
        'DehydratedDebris': S406-0,
        'PressLiquor': S406-1,
        'WashEffluent': C402-1,
        'SupernatantToWWT': M504-0,  # Config 2: supernatant to WWT
    }

# =============================================================================
# AREA 500: PURIFICATION
# =============================================================================
def create_area_500_purification(lysate_in, NaCl_wash, NaOH_elute, Ethanol_regen):
    M501 = bst.MixTank('M501', ins=(NaCl_wash, 'Water5'), outs='WashBuffer', tau=1)
    M502 = bst.MixTank('M502', ins=(NaOH_elute, 'Water6'), outs='ElutionBuffer', tau=1)
    M503 = bst.MixTank('M503', ins=(Ethanol_regen, 'Water7'), outs='RegenBuffer', tau=1)

    H501 = bst.HXutility('H501', ins=M501-0, outs='H501Out', T=30+273.15, heat_only=True)
    H502 = bst.HXutility('H502', ins=M502-0, outs='H502Out', T=40+273.15, heat_only=True)
    H503 = bst.HXutility('H503', ins=M503-0, outs='H503Out', T=30+273.15, heat_only=True)
    
    U501 = u.ResinColumnAdv(
        'U501',
        ins=(lysate_in, H501-0, H502-0, H503-0),
        outs=('ResinFlowthrough', 'ResinEluate', 'ResinRegenWaste', 'ResinWash'),
        preset='Adsorption',
        TargetProduct_IDs=('Heme_b', 'Heme_b_In', 'ProtoporphyrinIX', 'ProtoporphyrinIX_In'),
        TargetProduct_Yield=0.99,
        wash_CV=3,
        elution_CV=0.05,
        regeneration_CV=0.05,
        NonTarget_Removal=0.99,
        Wash_Impurity_Carryover=0.02,
        Regen_Impurity_Carryover=0.01
    )

    @M501.add_specification(run=True)
    def update_wash_buffer_adsorption():
        M501.ins[1].imass['H2O'] = lysate_in.imass['H2O'] * U501.wash_CV
        M501.ins[0].imass['NaCl'] = M501.ins[1].imass['H2O'] * 0.002 / 0.998

    @M502.add_specification(run=True)
    def update_elution_buffer_adsorption():
        M502.ins[1].imass['H2O'] = lysate_in.imass['H2O'] * U501.elution_CV
        M502.ins[0].imol['NaOH'] = M502.ins[1].imass['H2O'] * 0.2 / 1000

    @M503.add_specification(run=True)
    def update_regen_buffer_adsorption():
        M503.ins[1].imass['H2O'] = lysate_in.imass['H2O'] * U501.regeneration_CV * 0.3
        M503.ins[0].imass['Ethanol'] = M503.ins[1].imass['H2O'] * U501.regeneration_CV * 0.7

    return {
        'ResinEluate': U501-1,
        'ResinFlowthrough': U501-0,
        'ResinRegenWaste': U501-2,
        'ResinWash': U501-3,
    }

# =============================================================================
# AREA 600: CONCENTRATION
# =============================================================================
def create_area_600_concentration(eluate_in, DfUltraBuffer1):
    M601 = bst.MixTank('M601', ins=(DfUltraBuffer1, 'Water8'), outs='DFBuffer', tau=0.5)
    
    U601 = u.DiafiltrationAdv.from_preset(
        'NF', 'U601',
        ins=(eluate_in, M601-0),
        outs=('ConcentratedHeme', 'NFPermeate'),
        TargetProduct_IDs=('Heme_b', 'Heme_b_In', 'ProtoporphyrinIX', 'ProtoporphyrinIX_In',),
        Salt_IDs=c.chemical_groups['Salts'],
        TargetProduct_Retention=0.98, Salt_Retention=0.10,
        diavolumes=5.0,
    )
    
    @U601.add_specification(run=True)
    def adjust_df_buffer():
        feed = U601.ins[0]
        buffer = M601.ins[1]
        buffer.imass['H2O'] = feed.imass['H2O'] * 2.0
        M601.ins[0].imol['DfUltraBuffer'] = buffer.imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000

    H601 = bst.HXutility('H601', ins=U601-0, outs='HemeConcentrate', T=30+273.15)
    
    return {
        'HemeConcentrate': H601-0,
        'NFPermeate': U601-1,
    }

# =============================================================================
# AREA 700: FORMULATION
# =============================================================================
def create_area_700_formulation(heme_concentrate, GammaCyclodextrinFeed, NicotinamideFeed):
    M701 = bst.MixTank('M701', ins=(GammaCyclodextrinFeed, 'Water9'), outs='GammaCDSolution', tau=0.5)
    M702 = bst.Mixer('M702', ins=(heme_concentrate, M701-0), outs='HemeWithCD')
    
    encapsulation_rxn = bst.Rxn(
        '0.0014711 Heme_b + 0.0197913 GammaCyclodextrin -> HemoDextrin',
        reactant='Heme_b', X=0.95, check_atomic_balance=False,
    )
    coordination_rxn = bst.Rxn(
        'HemoDextrin + 0.0029422 Nicotinamide -> N-HemoDextrin',
        reactant='HemoDextrin', X=0.95, check_atomic_balance=False,
    )
    hemdx_rxns = bst.ReactionSystem(encapsulation_rxn, coordination_rxn)
    
    M703 = bst.MixTank('M703', ins=(NicotinamideFeed, 'Water10'), outs='NicotinamideSolution', tau=0.5)
    M704 = bst.Mixer('M704', ins=(M702-0, M703-0), outs='HemeDxWithNic')

    R702 = u.HemDxCSTR(
        'R702',
        ins=M704-0,
        outs='NHemDx_Crude',
        reactions=hemdx_rxns,
        tau=3.0, T=40 + 273.15, P=101325,
    )
    
    @M701.add_specification(run=True)
    def dose_gamma_cd():
        heme_mol = heme_concentrate.imol['Heme_b'] + heme_concentrate.imol['Heme_b_In']
        GammaCyclodextrinFeed.imol['GammaCyclodextrin'] = heme_mol * (0.0197913 / 0.0014711)
        M701.ins[1].imass['H2O'] = GammaCyclodextrinFeed.imass['GammaCyclodextrin'] * 10
    
    @M703.add_specification(run=True)
    def dose_nicotinamide():
        heme_mol = heme_concentrate.imol['Heme_b'] + heme_concentrate.imol['Heme_b_In']
        nic_mol = heme_mol * 2.0 * 1.025
        NicotinamideFeed.imol['Nicotinamide'] = nic_mol
        M703.ins[1].imass['H2O'] = NicotinamideFeed.imass['Nicotinamide'] * 5

    return R702-0

# =============================================================================
# AREA 800: FINAL PRODUCT
# =============================================================================
def create_area_800_final_product(crude_product, product_out, AntioxidantStream):
    """
    Area 800: Finalization
    
    Processing:
        1. U801: Diafiltration to remove excess salts/nicotinamide and concentrate
        2. H802: Pasteurization (HTST)
        3. M802: Formulation with Antioxidant
        4. H803: Final Cooling
        5. T801: Final product storage
    """
    U801 = u.DiafiltrationAdv.from_preset(
        'UF', 'U801',
        ins=(crude_product, 'Water11'),
        outs=('ConcentratedFinal', 'FinalPermeate'),
        TargetProduct_IDs=('N-HemoDextrin', 'HemoDextrin', 'Nicotinamide', 'GammaCyclodextrin'),
        Salt_IDs=c.chemical_groups['Salts'],
        TargetProduct_Retention=0.99,
        Salt_Retention=0.10,
        diavolumes=5.0,
    )
    
    @U801.add_specification(run=True)
    def control_final_concentration():
        feed = U801.ins[0]
        U801.ins[1].imass['H2O'] = feed.F_mass * 3.0 
        target_wt = 0.075
        U801.FeedWater_Recovery_to_Permeate = 0.5 
        U801._run()
        
        n_hemdx_mass = U801.outs[0].imass['N-HemoDextrin']
        
        if n_hemdx_mass > 0:
            target_total_mass = n_hemdx_mass / target_wt
            current_solids = U801.outs[0].F_mass - U801.outs[0].imass['H2O']
            target_water_retentate = max(0, target_total_mass - current_solids)
            feed_water = feed.imass['H2O']
            
            if feed_water > 0:
                recovery = 1.0 - (target_water_retentate / feed_water)
            else:
                recovery = 0.0
            
            U801.FeedWater_Recovery_to_Permeate = recovery
            
        U801._run()

    H802 = bst.HXutility(
        'H802', 
        ins=U801-0, 
        outs='PasteurizedProduct', 
        T=74+273.15, 
        heat_only=True
    )
    
    M802 = bst.MixTank('M802', ins=(H802-0, AntioxidantStream, 'Water12'), outs='PreCooledProduct', tau=0.5)

    H803 = bst.HXutility(
        'H803',
        ins=M802-0,
        outs='CooledProduct',
        T=4+273.15,
        cool_only=True
    )

    T801 = bst.StorageTank(
        'T801',
        ins=H803-0,
        outs=product_out,
        tau=24*7,
    )
    
    return {
        'FinalPermeate': U801-1
    }

# =============================================================================
# AREA 900: FACILITIES
# =============================================================================
def create_area_900_facilities(waste_streams, dehydrated_debris, ProcessWaste, emissions, ash_disposal):
    
    SupplementalNH4SO4 = bst.Stream('SupplementalNH4SO4', units='kg/hr', price=0)
    SupplementalFeSO4 = bst.Stream('SupplementalFeSO4', FeSO4=50, units='kg/hr', price=0)
    SupplementalNH4SO4.imass['(NH4)2SO4'] = 50 
    SupplementalNH4SO4.imass['NH3'] = 0.5  
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[*waste_streams, SupplementalNH4SO4, SupplementalFeSO4],
        outs=('biogas', 'sludge', 'RO_treated_water', ProcessWaste),
        mockup=True,
        area=500,
    )
    biogas, sludge, ro_water, brine = wastewater_treatment_sys.outs
    
    CT = bst.CoolingTower('CT')
    CWP = bst.ChilledWaterPackage('CWP')
    
    M902 = bst.Mixer(
        'M902',
        ins=(dehydrated_debris, sludge),
        outs='SolidsToBoiler',
    )
    
    BT = u.BoilerTurbogenerator(
        'BT',
        (M902-0, biogas, 'boiler_makeup_water', 'natural_gas', 'lime_boiler', 'boiler_chems'),
        outs=(emissions, 'rejected_water_and_blowdown', ash_disposal),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85,
        satisfy_system_electricity_demand=False,
    )
    
    F = bst.main_flowsheet
    makeup_water_streams = (
        F.cooling_tower_makeup_water,
        F.Water1, F.Water2, F.Water3, F.Water4,
        F.Water5, F.Water6, F.Water7,
        F.Water8, F.Water9, F.Water10, F.Water11, F.Water12,
        F.boiler_makeup_water,
    )
    process_water_streams = (
        ro_water,
        F.rejected_water_and_blowdown,
        *makeup_water_streams
    )
    
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    
    PWC = bst.ProcessWaterCenter(
        'PWC',
        ins=(ro_water, makeup_water, 'recycled_process_water', 'makeup_process_water'),
        outs=('RO_water', 'process_water', 'excess_water'),
        makeup_water_streams=makeup_water_streams,
        process_water_streams=process_water_streams,
        reverse_osmosis_water_price=0.000254,
        process_water_price=0.000135,
    )

# =============================================================================
# QA & PRODUCTION CONTROL
# =============================================================================
def check_HemDx_specifications(product_stream):
    """
    Verify HemDx product quality specifications.
    
    Specs:
        - Salt < 2.0 wt%
        - Residual Cyclodextrin < 4.0 wt%
        - Residual Nicotinamide < 2.0 wt%
        - Intermediate HemDx < 2.0 wt%
        - N-HemDx ~ 7.5 wt%
    """
    def get_mass_percent(chem_ids):
        if isinstance(chem_ids, str): chem_ids = [chem_ids]
        mass = sum([product_stream.imass[i] for i in chem_ids if i in product_stream.chemicals])
        return mass / product_stream.F_mass * 100 if product_stream.F_mass > 0 else 0

    specs = {
        'Salt': {
            'chemicals': ['NaCl', 'NaOH', 'KH2PO4','SodiumAscorbate'],
            'limit': 2.0, 
            'type': 'max'
        },
        'Residual CD': {
            'chemicals': ['GammaCyclodextrin'],
            'limit': 4.0,
            'type': 'max'
        },
        'Residual Nicotinamide': {
            'chemicals': ['Nicotinamide'],
            'limit': 2.0,
            'type': 'max'
        },
        'Intermediate HemDx': {
            'chemicals': ['HemoDextrin'],
            'limit': 2.0,
            'type': 'max'
        },
        'N-HemoDextrin': {
            'chemicals': ['N-HemoDextrin'],
            'target': (6.5, 8.5),
            'type': 'range'
        }
    }

    all_passed = True
    print(f"\n{'Specification':<30} {'Limit/Target':<20} {'Actual':<15} {'Status'}")
    print("-" * 80)
    
    for name, data in specs.items():
        actual = get_mass_percent(data['chemicals'])
        
        if data['type'] == 'max':
            passed = actual < data['limit']
            limit_str = f"< {data['limit']}%"
        elif data['type'] == 'range':
            passed = data['target'][0] <= actual <= data['target'][1]
            limit_str = f"{data['target'][0]}-{data['target'][1]}%"
            
        status = "[PASS]" if passed else "[FAIL]"
        print(f"{name:<30} {limit_str:<20} {actual:>6.2f}%{'':<8} {status}")
        
        if not passed: all_passed = False

    n_hemdx_mol = product_stream.imol['N-HemoDextrin']
    hemdx_mol = product_stream.imol['HemoDextrin']
    heme_equiv_mass_kg = (n_hemdx_mol + hemdx_mol) * 0.0014711 * 616.487
    equiv_wt = heme_equiv_mass_kg / product_stream.F_mass * 100 if product_stream.F_mass > 0 else 0
    print(f"{'Heme Equivalent':<30} {'(Info)':<20} {equiv_wt:>6.2f}%{'':<8} [INFO]")
    
    print("-" * 80)
    if not all_passed:
        print("WARNING: Product Specification Check Failed.")
        return False
    return True


def set_production_rate(system, target_production_rate_kg_hr, verbose=True):
    """
    Adjust system inputs to achieve target production rate.
    
    Simplified compared to _config2.py:
    - Does NOT call optimize_NH3_loading (handled by R302 specification)
    - Uses flexsolve IQ_interpolation for scaling factor only
    - Stores baseline_flows AFTER initial simulation (post-spec state)
    
    Parameters
    ----------
    system : biosteam.System
        The NHemDx production system
    target_production_rate_kg_hr : float
        Target mass flow rate for NHemDx_Product stream [kg/hr]
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
    
    # Identify product stream
    try:
        product_stream = system.flowsheet.stream.NHemDx_Product
    except AttributeError:
        product_stream = next(s for s in system.outs if 'Product' in s.ID or 'HemDx' in s.ID)
    
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
            
            # Bound expansion logic
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

# =============================================================================
# SYSTEM FACTORY
# =============================================================================
@bst.SystemFactory(
    ID='NHemDx_sys',
    ins=[
        s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt,
        s.NaCl_wash, s.NaOH_elute, s.Ethanol_regen, s.DfUltraBuffer1, s.DfUltraBuffer2,
        s.GammaCyclodextrinFeed, s.NicotinamideFeed, s.AntioxidantStream
    ],
    outs=[
        s.vent1, s.vent2,
        s.NHemDx_Product,
        dict(ID='ProcessWaste'),
        dict(ID='emissions'),
        dict(ID='ash_disposal'),
    ],
    fthermo=lambda chemicals=None: HEMDX_THERMO,
)
def create_NHemDx_system(ins, outs, use_area_convention=False):
    bst.preferences.N = 50
    (SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt,
     NaCl_wash, NaOH_elute, Ethanol_regen, DfUltraBuffer1, DfUltraBuffer2,
     GammaCyclodextrinFeed, NicotinamideFeed, AntioxidantStream) = ins
    (vent1, vent2, NHemDx_Product, ProcessWaste, emissions, ash_disposal) = outs

    bst.settings.set_thermo(HEMDX_THERMO, skip_checks=True)
    load_process_settings()
    
    # Set GWP characterization factors
    set_GWPCF(Glucose, 'Glucose')
    set_GWPCF(NH3_25wt, 'Ammonia_SEA', dilution=0.25)
    set_GWPCF(NaCl_wash, 'NaCl')
    set_GWPCF(NaOH_elute, 'NaOH')
    set_GWPCF(ash_disposal, 'ash_disposal')
    set_GWPCF(Ethanol_regen, 'Ethanol')
    set_GWPCF(GammaCyclodextrinFeed, 'GammaCyclodextrin')
    set_GWPCF(NicotinamideFeed, 'Nicotinamide')
    set_GWPCF(AntioxidantStream, 'AscorbicAcid', dilution=0.1)
    set_GWPCF_Multi(SeedIn1, ['AmmoniumSulfate', 'Glucose', 'MagnesiumSulfate', 'KH2PO4'],
                   [0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(SeedIn2, ['AmmoniumSulfate', 'Glucose', 'MagnesiumSulfate', 'KH2PO4'],
                   [0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(CultureIn, ['Glycine', 'Glucose', 'IronSulfate'],
                   [0.1/(0.1+60+0.15191), 60/(0.1+60+0.15191), 0.15191/(0.1+60+0.15191)])
    set_GWPCF_Multi(DfUltraBuffer1, ['KH2PO4', 'NaCl'], [0.8472, 0.1455])
    set_GWPCF_Multi(DfUltraBuffer2, ['KH2PO4', 'NaCl'], [0.8472, 0.1455])
    
    # Create Reactions
    params = get_fermentation_parameters()
    rxns = create_fermentation_reactions(params)

    # Area 200
    A200 = create_area_200_media_prep(SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt)
    
    # Area 300
    broth = create_area_300_conversion(A200['M203_out'], A200['T201_out'], A200['T202_out'], vent1, vent2, rxns, params,
                                        NH3_source_stream=NH3_25wt, S202_unit=A200['S202'])
    
    # Area 400: Recovery
    area_400 = create_area_400_recovery(broth, DfUltraBuffer2)
    
    # Area 500: Purification
    area_500 = create_area_500_purification(
        lysate_in=area_400['CombinedLysate'],
        NaCl_wash=NaCl_wash,
        NaOH_elute=NaOH_elute,
        Ethanol_regen=Ethanol_regen,
    )
    
    # Area 600: Concentration
    area_600 = create_area_600_concentration(area_500['ResinEluate'], DfUltraBuffer1)
    
    # Area 700: Formulation
    area_700 = create_area_700_formulation(area_600['HemeConcentrate'], GammaCyclodextrinFeed, NicotinamideFeed)
    
    # Area 800: Final Product
    area_800 = create_area_800_final_product(area_700, NHemDx_Product, AntioxidantStream)
    
    # Area 900: Facilities
    create_area_900_facilities(
        waste_streams=(
            area_400['PressLiquor'], area_400['WashEffluent'],
            area_400['SupernatantToWWT'],  # Config 2: supernatant directly to WWT
            area_500['ResinFlowthrough'], area_500['ResinRegenWaste'], area_500['ResinWash'],
            area_600['NFPermeate'],
            area_800['FinalPermeate'],
        ),
        dehydrated_debris=area_400['DehydratedDebris'],
        ProcessWaste=ProcessWaste,
        emissions=emissions,
        ash_disposal=ash_disposal
    )
    
    # Update input stream prices
    s.update_all_input_stream_prices(
        streamlist=[SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer1, DfUltraBuffer2, GammaCyclodextrinFeed, NicotinamideFeed, AntioxidantStream]
    )

    return NHemDx_Product, vent1, vent2, ProcessWaste, emissions, ash_disposal, DfUltraBuffer1, DfUltraBuffer2, AntioxidantStream
    
if __name__ == '__main__':
    bst.preferences.N = 50
    TARGET_PRODUCTION = 150  # kg/hr
    
    print("="*85)
    print("N-HemDx FULL PRODUCTION SYSTEM - Config 2 NEW (Intracellular, Internalized Specs)")
    print("="*85)
    
    print("\n1. Creating system...")
    NHemDx_sys = create_NHemDx_system()
    sys = NHemDx_sys
    sys.operating_hours = 8000
    f = sys.flowsheet
    u = f.unit
    ss = f.stream

    print("\n2. Setting production rate (NEW internalized approach)...")
    set_production_rate(sys, TARGET_PRODUCTION, verbose=True)

    print(f"\n4. System Summary")
    print("="*85)
    sys.show()    
    
    print("\n5. Product Analysis:")
    product = ss.NHemDx_Product
    
    print(f"\n6. Generating system diagram...")
    sys.diagram(format='html', display=True)

    # Run QA Compliance Check
    check_HemDx_specifications(product)

    n_hemdx = product.imass['N-HemoDextrin']
    hemdx = product.imass['HemoDextrin']
    heme = product.imass['Heme_b'] + product.imass['Heme_b_In']
    
    hemdx_mol = product.imol['HemoDextrin']
    n_hemdx_mol = product.imol['N-HemoDextrin']
    heme_equiv_mol = (hemdx_mol + n_hemdx_mol) * 0.0014711
    heme_equiv_mass = heme_equiv_mol * HEMDX_THERMO['Heme_b'].MW

    print(f"\n{'='*85}")
    print("KEY PERFORMANCE INDICATORS")
    print(f"{'='*85}")
    print(f"  Product Stream:      {product.ID}")
    print(f"  Target Rate:         {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"  Actual Rate:         {product.F_mass:.2f} kg/hr")
    print(f"  Annual Production:   {product.F_mass * sys.operating_hours / 1000:.2f} metric tons/year")
    print(f"  Yield (Heme eq):     {heme_equiv_mass:.4f} kg/hr")
    print(f"{'='*85}\n")
    
    print(f"\n6. Performing LCA analysis...")
    try:
        r1 = bst.report.lca_inventory_table(
            systems=[sys],
            keys='GWP',
            items=[product],
        )
        print("   LCA Inventory Table generated")
        print(r1)
        r2 = bst.report.lca_displacement_allocation_table(
            systems=[sys],
            key='GWP',
            items=[product],
        )
        print("   LCA Displacement Allocation Table generated")
        print(r2)
    except Exception as e:
        print(f"   LCA analysis failed: {e}")

    print(f"\n{'='*85}")
    print("SIMULATION COMPLETE")
    print(f"Target Production:   {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"Achieved Production: {ss.NHemDx_Product.F_mass:.2f} kg/hr")
    print(f"{'='*85}\n")
    
    print(f"target titer: {u.R302.titer:.2f} g/L")
    print(f"Achieved titer: {u.R302.actual_titer:.2f} g/L")
    print(f"Achieved titer2: {(ss.Broth.imass['Heme_b']+ss.Broth.imass['Heme_b_In']) / ss.Broth.F_vol:.2f} g/L")
