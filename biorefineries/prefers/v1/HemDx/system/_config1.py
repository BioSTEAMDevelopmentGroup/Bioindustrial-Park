# -*- coding: utf-8 -*-
"""
N-HemDx Full Production System - Config 1 (Modular Refactored Version)

Modular refactoring of N-HemDx production system into standard BioSTEAM Process Areas.

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

from biorefineries.prefers.v1._process_settings import set_GWPCF, GWP_CFs, set_GWPCF_Multi, load_process_settings
from biorefineries.prefers.v1.HemDx import _streams as s
import biosteam as bst
from thermosteam import Stream
from biosteam import F
import thermosteam as tmo
import numpy as np
from biorefineries.prefers.v1.HemDx import _chemicals as c
from biorefineries.prefers.v1 import _units as u
from biorefineries.prefers.v1._process_settings import price

# =============================================================================
# MODULE INITIALIZATION
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
        'tau1': 72,  # hours growth
        'tau2': 90,  # hours production
        'P_Heme_1st': 1e-8,
        'P_pp_1st': 1e-8,
        'P_Heme_2nd': 4.2/1e3,
        'P_pp_2nd': 0.7/1e3,
        'SF': 0.45,  # Secretion fraction
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
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                reactant='Glucose', X=Y_Heme*SF, check_atomic_balance=True),
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b_In + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
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
        'H2SO4 + 2 NH3 -> (NH4)2SO4', reactant='H2SO4', X=1,
        check_atomic_balance=True
    )
    
    cell_growth_reactionCG = bst.Rxn(
        'Glucose + 1.32 NH3 + 1.32 H2O-> 6 Corynebacterium_glutamicum + 2.37 O2', 'Glucose', X=Y_b,
        correct_atomic_balance=True
    )
    cell_growth_reactionCG.product_yield('Corynebacterium_glutamicum', basis='wt', product_yield=Y_b)
    
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
        bst.PRxn([cell_growth_reactionCG, respiration_reactionGC2, respiration_reactionGC3]),
    )
    
    return {
        'fermentation_reaction': fermentation_reaction,
        'cell_growth_reactionCG': cell_growth_reactionCG,
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
    M201 = bst.MixTank('M201', ins=[SeedIn1,'Water1'], outs='M201Out', tau=16)
    @M201.add_specification(run=True)
    def update_seed1_inputs():
        target_stream = bst.Stream(**s.SeedSolution1)
        SeedIn1.imass['Seed'] = target_stream.imass['Seed']
        M201.ins[1].imass['H2O'] = target_stream.imass['H2O']
        M201.ins[1].T = 25+273.15
    
    M202 = bst.MixTank('M202', ins=[SeedIn2,CultureIn,'Water2'], outs='M202Out', tau=16)
    @M202.add_specification(run=True)
    def update_culture_inputs():
        target_stream = bst.Stream(**s.SeedSolution2)
        SeedIn2.imass['Seed'] = target_stream.imass['Seed']
        M202.ins[2].imass['H2O'] = target_stream.imass['H2O']
        M202.ins[2].T = 25+273.15
        CultureIn.imass['Culture'] = target_stream.imass['SeedSolution']*(0.1+60+0.15191)/1000

    M203 = u.SeedHoldTank('M203', ins=[M201-0, M202-0], outs='M203Out')
    
    M204 = bst.MixTank('M204', ins=[Glucose,'Water3'], outs='M204Out', tau=16)
    @M204.add_specification(run=True)
    def update_water_content():
        M204.ins[1].imass['H2O'] = Glucose.imass['Glucose']/2
        M204.ins[1].T = 25+273.15
    
    T201 = bst.StorageTank('T201', ins=M204-0, outs='T201Out', tau=16*4+72)
    T202 = u.AmmoniaStorageTank('T202', ins=NH3_25wt, outs='T202Out')
    
    return {
        'M203_out': M203-0,
        'T201_out': T201-0,
        'T202_out': T202-0
    }

# =============================================================================
# AREA 300: CONVERSION
# =============================================================================
def create_area_300_conversion(seed_in, glucose_in, ammonia_in, vent1, vent2, rxns, params):
    R301 = u.SeedTrain(
        'R301',
        ins=[seed_in],
        outs=[vent1, 'R301Out'],
        reactions=bst.PRxn([rxns['cell_growth_reactionCG'], rxns['respiration_reactionGC2'], rxns['respiration_reactionGC3']]),
        saccharification=None,
        T=32+273.15,
    )
    R301.add_specification(run=True)
    
    R302 = u.AeratedFermentation(
        'R302',
        ins=[R301-1, glucose_in, ammonia_in, bst.Stream('FilteredAir', phase='g', P=2 * 101325)],
        outs=[vent2, 'Broth'],
        fermentation_reaction=rxns['fermentation_reaction'],
        cell_growth_reaction=rxns['cell_growth_reactionCG'],
        respiration_reaction=rxns['respiration_reactionGC2'],
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
    )
    
    # Calculate target yields for R302 specification
    P_Heme = (params['P_Heme_1st']*params['tau1'] + params['P_Heme_2nd']*params['tau2']) / (params['tau1']+params['tau2'])
    Y_Heme = rxns['Y_Heme']
    Y_pp = rxns['Y_pp']
    SF = params['SF']
    
    R302.target_titer = (params['tau1']+params['tau2']) * P_Heme
    R302.target_productivity = P_Heme
    R302.target_yield = Y_Heme

    @R302.add_specification(run=True)
    def update_reaction_time_and_yield():
        R302.tau = R302.target_titer / R302.target_productivity
        rxns['fermentation_reaction'][0].product_yield('Heme_b', basis='wt', product_yield=R302.target_yield*SF)
        rxns['fermentation_reaction'][1].product_yield('Heme_b_In', basis='wt', product_yield=R302.target_yield*(1-SF))
        rxns['fermentation_reaction'][2].product_yield('ProtoporphyrinIX', basis='wt', product_yield=Y_pp*SF)
        rxns['fermentation_reaction'][3].product_yield('ProtoporphyrinIX_In', basis='wt', product_yield=Y_pp*(1-SF))
        
    return R302-1

# =============================================================================
# AREA 400: RECOVERY
# =============================================================================
def create_area_400_recovery(broth_in):
    # S401: Primary Centrifuge
    centrifuge_order = ('Corynebacterium_glutamicum', 'Heme_b_In', 'ProtoporphyrinIX_In')
    S401 = u.Centrifuge(
        'S401',
        ins=broth_in,
        outs=('CellCream', 'Supernatant'),
        moisture_content=0.4,
        split=(0.999, 0.999, 0.999),
        order=centrifuge_order,
    )
    @S401.add_specification(run=True)
    def update_centrifuge_splits():
        inlet = S401.ins[0]
        S401._run()
        outlets_solid = S401.outs[0]
        if inlet.imass['H2O'] > 0:
            waterresidual = outlets_solid.imass['H2O'] / inlet.imass['H2O']
        else:
            waterresidual = 0.
        order_chemicals = set(centrifuge_order)
        for chem in S401.chemicals:
            chem_id = chem.ID
            if chem_id in order_chemicals:
                S401.isplit[chem_id] = 0.999
            elif chem_id == 'H2O':
                pass
            else:
                S401.isplit[chem_id] = waterresidual
        S401._run()

    # Supernatant microfiltration
    S404 = u.Filtration.from_preset(
        'MF', 'S404', ins=S401-1, outs=('SupernatantCake', 'FilteredSupernatant'),
        solid_capture_efficiency=0.95, cake_moisture_content=0.30,
    )
    S404.add_specification(run=True)

    # Cell cream washing & disruption
    DfUltraBuffer = bst.Stream(**s.DfUltraBuffer)
    M401 = bst.MixTank('M401', ins=(DfUltraBuffer, 'WashWater'), outs='WashBufferOut', tau=0.5)

    @M401.add_specification(run=True)
    def update_wash_buffer():
        M401.ins[1].imass['H2O'] = S401.outs[0].F_mass * 1.5
        M401.ins[0].imol['DfUltraBuffer'] = M401.ins[1].imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000

    H402 = bst.HXutility('H402', ins=M401-0, outs='ColdWashBuffer', T=10 + 273.15, cool_only=True)
    M402 = bst.MixTank('M402', ins=(S401-0, H402-0), outs='WashedCellSlurry', tau=0.25)

    C402 = u.Centrifuge(
        'C402',
        ins=M402-0,
        outs=('WashedCellCream', 'WashEffluent'),
        split={
            'Corynebacterium_glutamicum': 0.999, 'Heme_b_In': 0.999, 'ProtoporphyrinIX_In': 0.999,
            'Heme_b': 0.85, 'ProtoporphyrinIX': 0.85, 'Protein': 0.99, 'Cellulose': 0.99,
            'Xylan': 0.99, 'OleicAcid': 0.99, 'RNA': 0.99,
        },
        moisture_content=0.55,
    )

    S402 = u.CellDisruption(
        'S402', ins=C402-0, outs='CrudeHomogenate',
        Cell_ID='Corynebacterium_glutamicum',
        cell_disruption_efficiency=0.55, P_high=1000e5, P_low=101325,
        component_fractions={'Protein': 0.45, 'Cellulose': 0.22, 'Xylan': 0.15, 'OleicAcid': 0.08, 'RNA': 0.10}
    )

    H401 = bst.HXutility('H401', ins=S402-0, outs='CooledHomogenate', T=15 + 273.15, cool_only=True)

    S403 = bst.SolidsCentrifuge(
        'S403', ins=H401-0, outs=('CellDebrisSolids', 'CrudeLysate'),
        split={'Corynebacterium_glutamicum': 0.98, 'Protein': 0.98, 'Cellulose': 0.98, 'Xylan': 0.98,
               'OleicAcid': 0.98, 'RNA': 0.85, 'Heme_b_In': 0.02, 'ProtoporphyrinIX_In': 0.02,
               'Heme_b': 0.02, 'ProtoporphyrinIX': 0.02},
        moisture_content=0.20,
    )

    S405 = u.Filtration.from_preset(
        'MF', 'S405', ins=S403-1, outs=('FilterCake', 'ClarifiedLysate'),
        solid_capture_efficiency=0.95, cake_moisture_content=0.30,
    )
    S405.add_specification(run=True)

    # Debris Dewatering
    M404 = bst.Mixer('M404', ins=(S403-0, S405-0, S404-0), outs='CellDebrisRaw')
    S406 = bst.ScrewPress(
        'S406', ins=M404-0, outs=('DehydratedDebris', 'PressLiquor'),
        split=0.999, moisture_content=0.001
    )

    M405 = bst.Mixer('M405', ins=(S405-1, S404-1), outs='CombinedLysate')
    
    return {
        'CombinedLysate': M405-0,
        'DehydratedDebris': S406-0,
        'PressLiquor': S406-1,
        'WashEffluent': C402-1,
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
    
    U501 = u.ResinColumn(
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
        M503.ins[0].imass['Ethanol'] = M503.ins[1].imass['H2O'] * U501.regeneration_CV * 0.7  # Fixed bug here (was M405 in orig)

    return {
        'ResinEluate': U501-1,
        'ResinFlowthrough': U501-0,
        'ResinRegenWaste': U501-2,
        'ResinWash': U501-3,
    }

# =============================================================================
# AREA 600: CONCENTRATION
# =============================================================================
def create_area_600_concentration(eluate_in, DFBufferSolute):
    M601 = bst.MixTank('M601', ins=(DFBufferSolute, 'DFBufferWater'), outs='DFBuffer', tau=0.5)
    
    U601 = u.Diafiltration.from_preset(
        'NF', 'U601',
        ins=(eluate_in, M601-0),
        outs=('ConcentratedHeme', 'NFPermeate'),
        TargetProduct_ID=('Heme_b', 'Heme_b_In', 'ProtoporphyrinIX', 'ProtoporphyrinIX_In'),
        Salt_ID=('NaCl', 'NaOH', 'Ethanol'),
        TargetProduct_Retention=0.98, Salt_Retention=0.10, concentration_factor=5.0,
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
    M701 = bst.MixTank('M701', ins=(GammaCyclodextrinFeed, 'Water8'), outs='GammaCDSolution', tau=0.5)
    M702 = bst.Mixer('M702', ins=(heme_concentrate, M701-0), outs='HemeWithCD')
    
    encapsulation_rxn = bst.Rxn(
        '0.0014711 Heme_b + 0.0197913 GammaCyclodextrin -> HemoDextrin',
        reactant='Heme_b', X=0.95, check_atomic_balance=False,
    )
    encapsulation_rxn_in = bst.Rxn(
        '0.0014711 Heme_b_In + 0.0197913 GammaCyclodextrin -> HemoDextrin',
        reactant='Heme_b_In', X=0.95, check_atomic_balance=False,
    )
    coordination_rxn = bst.Rxn(
        'HemoDextrin + 0.0029422 Nicotinamide -> N-HemoDextrin',
        reactant='HemoDextrin', X=0.95, check_atomic_balance=False,
    )
    hemdx_rxns = bst.ReactionSystem(encapsulation_rxn, encapsulation_rxn_in, coordination_rxn)
    
    M703 = bst.MixTank('M703', ins=(NicotinamideFeed, 'Water9'), outs='NicotinamideSolution', tau=0.5)
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
        # Access H601 output (heme_concentrate) indirectly or assume it's the input
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
def create_area_800_final_product(crude_product, product_out):
    """
    Area 800: Finalization
    
    Processing:
        1. U801: Diafiltration to remove excess salts/nicotinamide and concentrate
        2. M801: Formulation buffer addition (if needed for concentration adj)
        3. H802: Pasteurization (HTST)
        4. T801: Final product storage
    """
    # U801: Diafiltration (Concentration & Purification)
    # Remove excess Nicotinamide, salts, and Cyclodextrin
    M801 = bst.Mixer('M801', ins=(crude_product, 'FormulationBuffer'))
    
    U801 = u.Diafiltration.from_preset(
        'UF', 'U801',
        ins=(M801-0, 'FinalDFWater'),
        outs=('ConcentratedFinal', 'FinalPermeate'),
        TargetProduct_ID=('N-HemoDextrin', 'HemoDextrin'),
        Salt_ID=('Nicotinamide', 'GammaCyclodextrin', 'Sodium', 'Chloride', 'Ethanol'),
        TargetProduct_Retention=0.99,
        Salt_Retention=0.10,
        concentration_factor=2.0, # Initial guess, controlled by spec
    )
    
    # QA/Concentration Control
    # Target 7.5% N-HemDx concentration
    @U801.add_specification(run=True)
    def control_final_concentration():
        # 1. Run Baseline Diafiltration (Salt Removal)
        # Ensure sufficient DF water for washing (e.g., 3 diavolumes)
        feed = U801.ins[0]
        U801.ins[1].imass['H2O'] = feed.F_mass * 3.0 
        
        # 2. Target Concentration Calculation
        target_wt = 0.075 # 7.5%
        
        # To get the right solute mass, we assume retention maps are static or run once
        # Diafiltration retention params are fixed by init. 
        # But we need to run it to know how much solute is retained.
        # Run with a dummy water recovery first.
        U801.FeedWater_Recovery_to_Permeate = 0.5 
        U801._run()
        
        n_hemdx_mass = U801.outs[0].imass['N-HemoDextrin']
        
        if n_hemdx_mass > 0:
            target_total_mass = n_hemdx_mass / target_wt
            current_solids = U801.outs[0].F_mass - U801.outs[0].imass['H2O']
            # Target water in retentate
            target_water_retentate = max(0, target_total_mass - current_solids)
            
            feed_water = feed.imass['H2O']
            
            # Calculate Recovery Factor:
            # retentate_water = feed_water * (1 - R)
            # R = 1 - (retentate_water / feed_water)
            
            if feed_water > 0:
                recovery = 1.0 - (target_water_retentate / feed_water)
            else:
                recovery = 0.0 # Should not happen if feed is aqueous
            
            U801.FeedWater_Recovery_to_Permeate = recovery
            
        # 3. Final Run
        U801._run()

    # H802: Pasteurization (HTST) 74°C, 15s -> Cool to 4°C
    H802 = bst.HXutility(
        'H802', 
        ins=U801-0, 
        outs='PasteurizedProduct', 
        T=74+273.15, 
        heat_only=True
    )
    
    H803 = bst.HXutility(
        'H803',
        ins=H802-0,
        outs='CooledProduct',
        T=4+273.15,
        cool_only=True
    )

    # T801: Storage
    T801 = bst.StorageTank(
        'T801',
        ins=H803-0,
        outs=product_out,
        tau=24*7, # 1 week storage
    )
    
    return {
        'FinalPermeate': U801-1
    }

# =============================================================================
# AREA 900: FACILITIES
# =============================================================================
def create_area_900_facilities(waste_streams, dehydrated_debris, ProcessWaste, emissions, ash_disposal):
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=waste_streams,
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
        F.Water5, F.Water6, F.Water7,
        F.Water8, F.Water9,
        F.DFBufferWater,
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
        - Salt < 1.0 wt%
        - Residual Cyclodextrin < 1.0 wt%
        - Residual Nicotinamide < 1.0 wt%
        - Intermediate HemDx < 1.0 wt%
        - N-HemDx ~ 7.5 wt%
    """
    def get_mass_percent(chem_ids):
        if isinstance(chem_ids, str): chem_ids = [chem_ids]
        mass = sum([product_stream.imass[i] for i in chem_ids if i in product_stream.chemicals])
        return mass / product_stream.F_mass * 100 if product_stream.F_mass > 0 else 0

    specs = {
        'Salt': {
            'chemicals': ['NaCl', 'NaOH', 'Sodium', 'Chloride', 'Na+', 'Cl-'], # Covers ions if split
            'limit': 1.0, 
            'type': 'max'
        },
        'Residual CD': {
            'chemicals': ['GammaCyclodextrin'],
            'limit': 1.0,
            'type': 'max'
        },
        'Residual Nicotinamide': {
            'chemicals': ['Nicotinamide'],
            'limit': 1.0,
            'type': 'max'
        },
        'Intermediate HemDx': {
            'chemicals': ['HemoDextrin'],
            'limit': 1.0,
            'type': 'max'
        },
        'N-HemoDextrin': {
            'chemicals': ['N-HemoDextrin'],
            'target': (6.5, 8.5), # 7.5 +/- 1.0
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

    # Calculate Heme Equivalent
    n_hemdx_mol = product_stream.imol['N-HemoDextrin']
    hemdx_mol = product_stream.imol['HemoDextrin']
    heme_equiv_mass_kg = (n_hemdx_mol + hemdx_mol) * 0.0014711 * 616.487 # MW Heme_b
    equiv_wt = heme_equiv_mass_kg / product_stream.F_mass * 100 if product_stream.F_mass > 0 else 0
    print(f"{'Heme Equivalent':<30} {'(Info)':<20} {equiv_wt:>6.2f}%{'':<8} [INFO]")
    
    print("-" * 80)
    if not all_passed:
        # Don't raise error yet, just warn
        print("WARNING: Product Specification Check Failed.")
        return False
    return True

def set_production_rate(sys, target_rate):
    """
    Iteratively adjust system capacity to match target production rate.
    Adjusts: Glucose feed flow (which cascades to others via specs).
    """
    # Find Glucose stream
    glucose = next(i for i in sys.ins if i.ID == 'Glucose')
    product = next(i for i in sys.outs if i.ID == 'NHemDx_Product')
    
    print(f"\nAuto-adjusting production to {target_rate} kg/hr...")
    
    # Run optimization
    # Use simple feedback loop for speed/robustness vs heavy solver
    for i in range(5):
        sys.simulate()
        current_rate = product.F_mass
        error = target_rate - current_rate
        print(f"  Iter {i+1}: Glucose={glucose.F_mass:.2f}, Product={current_rate:.2f} (Target {target_rate})")
        
        if abs(error) < 1.0: 
            print("  Converged.")
            break
        
        # Correction factor (dampened)
        factor = target_rate / current_rate if current_rate > 0 else 1.1
        # Limit jump to 2x or 0.5x to avoid instability
        factor = max(0.5, min(2.0, factor))
        glucose.F_mass *= factor
        
    final_rate = product.F_mass
    print(f"Set Glucose flow to {glucose.F_mass:.2f} kg/hr -> Product: {final_rate:.2f} kg/hr")
    return final_rate

# =============================================================================
# SYSTEM FACTORY
# =============================================================================
@bst.SystemFactory(
    ID='NHemDx_sys',
    ins=[
        s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt,
        s.NaCl_wash, s.NaOH_elute, s.Ethanol_regen, s.DFBufferSolute,
        s.GammaCyclodextrinFeed, s.NicotinamideFeed,
    ],
    outs=[
        s.vent1, s.vent2,
        dict(ID='NHemDx_Product'),
        dict(ID='ProcessWaste'),
        dict(ID='emissions'),
        dict(ID='ash_disposal'),
    ],
    fthermo=lambda chemicals=None: HEMDX_THERMO,
)
def create_NHemDx_system(ins, outs, use_area_convention=False):
    bst.preferences.N = 50
    (SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt,
     NaCl_wash, NaOH_elute, Ethanol_regen, DFBufferSolute,
     GammaCyclodextrinFeed, NicotinamideFeed) = ins
    (vent1, vent2, NHemDx_Product, ProcessWaste, emissions, ash_disposal) = outs

    bst.settings.set_thermo(HEMDX_THERMO, skip_checks=True)
    load_process_settings()
    
    # Set GWP characterization factors
    # Knowns from LegHb/Process Settings
    set_GWPCF(Glucose, 'Glucose')
    set_GWPCF(NH3_25wt, 'Ammonia_US') # Assume US for now
    set_GWPCF(NaCl_wash, 'NaCl')
    set_GWPCF(NaOH_elute, 'NaOH')
    set_GWPCF(ash_disposal, 'ash_disposal')
    
    # Placeholder for unknowns (0.5 as requested)
    # Ethanol, DFBuffer, Cyclodextrin, Nicotinamide
    for s_unknown in [Ethanol_regen, DFBufferSolute, GammaCyclodextrinFeed, NicotinamideFeed]:
        s_unknown.characterization_factors['GWP'] = 0.5 # assumed temporarily
        
    # Complex feeds (Seed/Culture) - similar to LegHb but simplified/placeholder if composition varies
    # Or use set_GWPCF_Multi if we know composition. 
    # For now, default to placeholder or simple estimate
    SeedIn1.characterization_factors['GWP'] = 0.5 # assumed temporarily
    SeedIn2.characterization_factors['GWP'] = 0.5 # assumed temporarily
    CultureIn.characterization_factors['GWP'] = 0.5 # assumed temporarily

    # Create Reactions
    params = get_fermentation_parameters()
    rxns = create_fermentation_reactions(params)

    # Area 200
    A200 = create_area_200_media_prep(SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt)
    
    # Area 300
    broth = create_area_300_conversion(A200['M203_out'], A200['T201_out'], A200['T202_out'], vent1, vent2, rxns, params)
    
    # Area 400
    A400 = create_area_400_recovery(broth)
    
    # Area 500
    A500 = create_area_500_purification(A400['CombinedLysate'], NaCl_wash, NaOH_elute, Ethanol_regen)
    
    # Area 600
    A600 = create_area_600_concentration(A500['ResinEluate'], DFBufferSolute)
    
    # Area 700
    crude_product = create_area_700_formulation(A600['HemeConcentrate'], GammaCyclodextrinFeed, NicotinamideFeed)
    
    # Area 800
    # Area 800
    A800 = create_area_800_final_product(crude_product, NHemDx_Product)
    
    # Area 900
    waste_streams = [
        A500['ResinFlowthrough'],
        A500['ResinRegenWaste'],
        A500['ResinWash'],
        A600['NFPermeate'],
        A400['WashEffluent'],
        A400['PressLiquor'],
        A800['FinalPermeate'], # New waste stream
    ]
    
    create_area_900_facilities(
        waste_streams,
        A400['DehydratedDebris'],
        ProcessWaste,
        emissions,
        ash_disposal
    )

if __name__ == '__main__':
    bst.preferences.N = 50
    TARGET_PRODUCTION = 275  # kg/hr
    
    print("="*85)
    print("N-HemDx FULL PRODUCTION SYSTEM - Upstream + DSP Integration")
    print("="*85)
    
    print("\n1. Creating system...")
    NHemDx_sys = create_NHemDx_system()
    sys = NHemDx_sys
    sys.operating_hours = 8000
    
    print("\n2. Setting production rate...")
    try:
        set_production_rate(sys, TARGET_PRODUCTION)
    except Exception as e:
        print(f"   Production rate setting error: {e}")

    print(f"\n3. System Summary")
    print("="*85)
    sys.show()    
    
    print("\n" + "="*85)
    print("SIMULATION COMPLETE")
    print("="*85)
    
    sys.diagram(format='html')
    
    print("\n4. Product Analysis:")
    product = sys.flowsheet.stream.NHemDx_Product
    
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
    
    print(f"\n5. Performing LCA analysis...")
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
    print(f"{'='*85}\n")

