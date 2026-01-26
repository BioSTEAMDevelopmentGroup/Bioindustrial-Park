# -*- coding: utf-8 -*-
"""
N-HemDx Full Production System - Config 1 (Integrated Upstream + DSP)

Complete production system for N-HemDx (Nicotinamide-Stabilized Hemodextrin)
integrating upstream fermentation with downstream processing.

Process Areas:
    - Area 200: Media Preparation (Seed, Glucose, Ammonia)
    - Area 300: Conversion (Fermentation - Corynebacterium glutamicum)
    - Area 400: Clarification & Cell Disruption (Split-stream processing)
    - Area 500: Capture & Purification (Resin column for both streams)
    - Area 600: Concentration (NF/Diafiltration)
    - Area 700: Formulation (γ-CD + Nicotinamide complexation)
    - Area 800: Final Product
    - Area 900: Facilities & Utilities

Host: Corynebacterium glutamicum
Product: N-HemDx (Heme b + γ-Cyclodextrin + Nicotinamide complex)

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
)

# =============================================================================
# PROCESS PARAMETERS
# =============================================================================


class HemDxCSTR(bst.CSTR):
    _N_ins = 1
    _N_outs = 1

    def _init(self, reactions=None, **kwargs):
        super()._init(**kwargs)
        self.reactions = reactions

    def _run(self):
        out, = self.outs
        out.mix_from(self.ins, energy_balance=False)
        if self.reactions:
            self.reactions(out)
        out.T = self.T
        out.P = 101325 if self.P is None else self.P

# Fermentation Parameters (from upstream config)
theta_O2 = 0.45
agitation_power = 0.985
design = 'Stirred tank'
method = "Riet"
T_operation = 273.15 + 30
Q_O2_consumption = -110 * 4184
dT_hx_loop = 8
cooler_pressure_drop = 20684
compressor_isentropic_efficiency = 0.85
V_max = 500

# Fermentation kinetics
tau1 = 72  # hours growth
tau2 = 90  # hours production
tau = tau1 + tau2

P_Heme_1st = 1e-8
P_pp_1st = 1e-8
P_ComsuptionG1 = 500 * (0.5/3) / tau1
P_cell_growth_rate1 = 150 * 0.35 / tau1
Y_Heme1 = P_Heme_1st / P_ComsuptionG1
Y_pp1 = P_pp_1st / P_ComsuptionG1
Y_b1 = P_cell_growth_rate1 / P_ComsuptionG1

P_Heme_2nd = 4.2/1e3
P_pp_2nd = 0.7/1e3
P_ComsuptionG2 = 800 * (0.5/3) / tau2
P_cell_growth_rate2 = 50 * 0.35 / tau2
Y_Heme2 = P_Heme_2nd / P_ComsuptionG2
Y_pp2 = P_pp_2nd / P_ComsuptionG2
Y_b2 = P_cell_growth_rate2 / P_ComsuptionG2

P_Heme = (P_Heme_1st*tau1 + P_Heme_2nd*tau2) / tau
P_ComsuptionG = (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2) / tau
Y_Heme = (P_Heme_1st*tau1 + P_Heme_2nd*tau2) / (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2)
Y_pp = (P_pp_1st*tau1 + P_pp_2nd*tau2) / (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2)
Y_b = (P_cell_growth_rate1*tau1 + P_cell_growth_rate2*tau2) / (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2)
Y_r = 1 - Y_Heme - Y_pp - Y_b

SF = 0.45  # Secretion fraction of HemeB

# Formulation Parameters
NICOTINAMIDE_EXCESS = 5.0  # 5× molar excess
HPH_PRESSURE = 1000e5  # 1000 bar

# =============================================================================
# SYSTEM FACTORY
# =============================================================================

@bst.SystemFactory(
    ID='NHemDx_sys',
    ins=[
        s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt,
        # DSP inputs (solutes only; water is prepared internally)
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
def create_NHemDx_system(
        ins, outs,
        use_area_convention=False,
    ):
    """
    Creates the complete N-HemDx production system with upstream fermentation
    and downstream processing including formulation.
    """
    bst.preferences.N = 50
    
    # Unpack input streams
    (SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt,
     NaCl_wash, NaOH_elute, Ethanol_regen, DFBufferSolute,
     GammaCyclodextrinFeed, NicotinamideFeed) = ins

    # Unpack output streams
    (vent1, vent2, NHemDx_Product, ProcessWaste, emissions, ash_disposal) = outs

    
    # Ensure HemDx chemicals are active
    bst.settings.set_thermo(HEMDX_THERMO, skip_checks=True)
    load_process_settings()
    
    # =========================================================================
    # REACTIONS (from upstream config - DO NOT MODIFY)
    # =========================================================================
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
    
    # =========================================================================
    # AREA 200: MEDIA PREPARATION (DO NOT MODIFY)
    # =========================================================================
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

    R301 = u.SeedTrain(
        'R301',
        ins=[M203-0],
        outs=[vent1, 'R301Out'],
        reactions=bst.PRxn([cell_growth_reactionCG, respiration_reactionGC2, respiration_reactionGC3]),
        saccharification=None,
        T=32+273.15,
    )
    R301.add_specification(run=True)
    
    M204 = bst.MixTank('M204', ins=[Glucose,'Water3'], outs='M204Out', tau=16)
    
    @M204.add_specification(run=True)
    def update_water_content():
        M204.ins[1].imass['H2O'] = Glucose.imass['Glucose']/2
        M204.ins[1].T = 25+273.15
    
    T201 = bst.StorageTank('T201', ins=M204-0, outs='T201Out', tau=16*4+72)
    T202 = u.AmmoniaStorageTank('T202', ins=NH3_25wt, outs='T202Out')

    # =========================================================================
    # AREA 300: CONVERSION - FERMENTATION (DO NOT MODIFY)
    # =========================================================================
    R302 = u.AeratedFermentation(
        'R302',
        ins=[R301-1, T201-0, T202-0, bst.Stream('FilteredAir', phase='g', P=2 * 101325)],
        outs=[vent2, 'Broth'],
        fermentation_reaction=fermentation_reaction,
        cell_growth_reaction=cell_growth_reactionCG,
        respiration_reaction=respiration_reactionGC2,
        neutralization_reaction=neutralization_reaction,
        design='Stirred tank', method=method, theta_O2=theta_O2,
        V_max=V_max, Q_O2_consumption=Q_O2_consumption,
        dT_hx_loop=dT_hx_loop, T=T_operation,
        batch=True, reactions=RXN,
        kW_per_m3=agitation_power,
        tau=tau,
        cooler_pressure_drop=cooler_pressure_drop,
        compressor_isentropic_efficiency=compressor_isentropic_efficiency,
        P=1 * 101325,
    )
    R302.target_titer = tau*P_Heme
    R302.target_productivity = P_Heme
    R302.target_yield = Y_Heme

    @R302.add_specification(run=True)
    def update_reaction_time_and_yield():
        R302.tau = R302.target_titer / R302.target_productivity
        fermentation_reaction[0].product_yield('Heme_b', basis='wt', product_yield=R302.target_yield*SF)
        fermentation_reaction[1].product_yield('Heme_b_In', basis='wt', product_yield=R302.target_yield*(1-SF))
        fermentation_reaction[2].product_yield('ProtoporphyrinIX', basis='wt', product_yield=Y_pp*SF)
        fermentation_reaction[3].product_yield('ProtoporphyrinIX_In', basis='wt', product_yield=Y_pp*(1-SF))

    # =========================================================================
    # AREA 400: CLARIFICATION & CELL DISRUPTION (Split-Stream)
    # =========================================================================
    
    # S401: Primary Centrifuge - Split broth into supernatant and cell cream
    centrifuge_order = ('Corynebacterium_glutamicum', 'Heme_b_In', 'ProtoporphyrinIX_In')
    S401 = u.Centrifuge(
        'S401',
        ins=R302-1,
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

    # Supernatant microfiltration (deep removal of residual cells)
    S404 = u.Filtration.from_preset(
        'MF',
        'S404',
        ins=S401-1,
        outs=('SupernatantCake', 'FilteredSupernatant'),
        solid_capture_efficiency=0.95,
        cake_moisture_content=0.30,
    )
    S404.add_specification(run=True)

    # Cell cream washing & disruption sequence (LegHb-style)
    DfUltraBuffer = bst.Stream(**s.DfUltraBuffer)

    M401 = bst.MixTank('M401', ins=(DfUltraBuffer, 'WashWater'),
                       outs='WashBufferOut', tau=0.5)

    @M401.add_specification(run=True)
    def update_wash_buffer():
        M401.ins[1].imass['H2O'] = S401.outs[0].F_mass * 1.5
        M401.ins[0].imol['DfUltraBuffer'] = M401.ins[1].imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000

    H402 = bst.HXutility(
        'H402',
        ins=M401-0,
        outs='ColdWashBuffer',
        T=10 + 273.15,
        cool_only=True,
    )

    M402 = bst.MixTank('M402', ins=(S401-0, H402-0),
                       outs='WashedCellSlurry', tau=0.25)

    C402 = u.Centrifuge(
        'C402',
        ins=M402-0,
        outs=('WashedCellCream', 'WashEffluent'),
        split={
            'Corynebacterium_glutamicum': 0.999,
            'Heme_b_In': 0.999,
            'ProtoporphyrinIX_In': 0.999,
            'Heme_b': 0.85,
            'ProtoporphyrinIX': 0.85,
            'Protein': 0.99,
            'Cellulose': 0.99,
            'Xylan': 0.99,
            'OleicAcid': 0.99,
            'RNA': 0.99,
        },
        moisture_content=0.55,
    )

    # S402: Cell Disruption (HPH) - releases intracellular heme
    S402 = u.CellDisruption(
        'S402', ins=C402-0, outs='CrudeHomogenate',
        Cell_ID='Corynebacterium_glutamicum',
        cell_disruption_efficiency=0.55,
        P_high=HPH_PRESSURE,
        P_low=101325,
        component_fractions={
            'Protein': 0.45,
            'Cellulose': 0.22,
            'Xylan': 0.15,
            'OleicAcid': 0.08,
            'RNA': 0.10,
        }
    )

    H401 = bst.HXutility(
        'H401',
        ins=S402-0,
        outs='CooledHomogenate',
        T=15 + 273.15,
        cool_only=True,
    )

    S403 = bst.SolidsCentrifuge(
        'S403',
        ins=H401-0,
        outs=('CellDebrisSolids', 'CrudeLysate'),
        split={
            'Corynebacterium_glutamicum': 0.98,
            'Protein': 0.98,
            'Cellulose': 0.98,
            'Xylan': 0.98,
            'OleicAcid': 0.98,
            'RNA': 0.85,
            'Heme_b_In': 0.02,
            'ProtoporphyrinIX_In': 0.02,
            'Heme_b': 0.02,
            'ProtoporphyrinIX': 0.02,
        },
        moisture_content=0.20,
    )

    S405 = u.Filtration.from_preset(
        'MF',
        'S405',
        ins=S403-1,
        outs=('FilterCake', 'ClarifiedLysate'),
        solid_capture_efficiency=0.95,
        cake_moisture_content=0.30,
    )
    S405.add_specification(run=True)

    # Route debris to screw press for dewatering (like LegHb config1)
    M404 = bst.Mixer('M404', ins=(S403-0, S405-0, S404-0), outs='CellDebrisRaw')

    # Screw press for debris dewatering (replaces old S406 splitter)
    S406 = bst.ScrewPress(
        'S406',
        ins=M404-0,
        outs=('DehydratedDebris', 'PressLiquor'),
        split=0.999,
        moisture_content=0.001
    )

    # Combine clarified lysate + filtered supernatant
    M405 = bst.Mixer('M405', ins=(S405-1, S404-1), outs='CombinedLysate')

    # =========================================================================
    # AREA 500: CAPTURE & PURIFICATION (Single Resin Column)
    # =========================================================================
    
    M501 = bst.MixTank('M501', ins=(NaCl_wash, 'Water5'), outs='WashBuffer', tau=1)
    M502 = bst.MixTank('M502', ins=(NaOH_elute, 'Water6'), outs='ElutionBuffer', tau=1)
    M503 = bst.MixTank('M503', ins=(Ethanol_regen, 'Water7'), outs='RegenBuffer', tau=1)

    H501 = bst.HXutility('H501', ins=M501-0, outs='H501Out', T=30+273.15, heat_only=True)
    H502 = bst.HXutility('H502', ins=M502-0, outs='H502Out', T=40+273.15, heat_only=True)
    H503 = bst.HXutility('H503', ins=M503-0, outs='H503Out', T=30+273.15, heat_only=True)
    
    U501 = u.ResinColumn(
        'U501',
        ins=(M405-0, H501-0, H502-0, H503-0),
        outs=('ResinFlowthrough', 'ResinEluate', 'ResinRegenWaste', 'ResinWash'),
        preset='Adsorption',
        TargetProduct_IDs=('Heme_b', 'Heme_b_In', 'ProtoporphyrinIX', 'ProtoporphyrinIX_In'),
        TargetProduct_Yield=0.99,
        wash_CV=3,
        elution_CV=0.05,
        regeneration_CV=0.05,
    )

    @M501.add_specification(run=True)
    def update_wash_buffer_adsorption():
        M501.ins[1].imass['H2O'] = M405.outs[0].imass['H2O'] * U501.wash_CV
        M501.ins[0].imass['NaCl'] = M501.ins[1].imass['H2O'] * 0.002 / 0.998

    @M502.add_specification(run=True)
    def update_elution_buffer_adsorption():
        M502.ins[1].imass['H2O'] = M405.outs[0].imass['H2O'] * U501.elution_CV
        M502.ins[0].imol['NaOH'] = M502.ins[1].imass['H2O'] * 0.2 / 1000

    @M503.add_specification(run=True)
    def update_regen_buffer_adsorption():
        M503.ins[1].imass['H2O'] = M405.outs[0].imass['H2O'] * U501.regeneration_CV * 0.3
        M503.ins[0].imass['Ethanol'] = M405.outs[0].imass['H2O'] * U501.regeneration_CV * 0.7

    # =========================================================================
    # AREA 600: CONCENTRATION (NF/Diafiltration)
    # =========================================================================
    
    # DF buffer preparation (solute feed + water)
    M601 = bst.MixTank('M601', ins=(DFBufferSolute, 'DFBufferWater'), outs='DFBuffer', tau=0.5)
    
    U601 = u.Diafiltration.from_preset(
        'NF',
        ID='U601',
        ins=(U501-1, M601-0),
        outs=('ConcentratedHeme', 'NFPermeate'),
        TargetProduct_ID=('Heme_b', 'Heme_b_In', 'ProtoporphyrinIX', 'ProtoporphyrinIX_In'),
        Salt_ID=('NaCl', 'NaOH', 'Ethanol'),
        TargetProduct_Retention=0.98,
        Salt_Retention=0.10,
        concentration_factor=5.0,
    )
    
    @U601.add_specification(run=True)
    def adjust_df_buffer():
        feed = U601.ins[0]
        buffer = M601.ins[1]
        buffer.imass['H2O'] = feed.imass['H2O'] * 2.0
        M601.ins[0].imol['DfUltraBuffer'] = buffer.imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000


    H601 = bst.HXutility('H601', ins=U601-0, outs='HemeConcentrate', T=30+273.15)

    # =========================================================================
    # AREA 700: FORMULATION (γ-CD + Nicotinamide Complexation)
    # =========================================================================
    
    # M701: γ-Cyclodextrin solution preparation
    M701 = bst.MixTank('M701', ins=(GammaCyclodextrinFeed, 'Water8'), outs='GammaCDSolution', tau=0.5)

    # M702: γ-Cyclodextrin Addition
    M702 = bst.Mixer('M702', ins=(H601-0, M701-0), outs='HemeWithCD')
    
    # R702: Encapsulation + Coordination (single CSTR, series reactions)
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
    
    # M703: Nicotinamide solution preparation
    M703 = bst.MixTank('M703', ins=(NicotinamideFeed, 'Water9'), outs='NicotinamideSolution', tau=0.5)

    # M704: Nicotinamide Addition
    M704 = bst.Mixer('M704', ins=(M702-0, M703-0), outs='HemeDxWithNic')

    R702 = HemDxCSTR(
        'R702',
        ins=M704-0,
        outs='NHemDx_Crude',
        reactions=hemdx_rxns,
        tau=3.0,
        T=40 + 273.15,
        P=101325,
    )
    
    # Dynamic dosing specifications
    @M701.add_specification(run=True)
    def dose_gamma_cd():
        heme_mol = H601.outs[0].imol['Heme_b'] + H601.outs[0].imol['Heme_b_In']
        GammaCyclodextrinFeed.imol['GammaCyclodextrin'] = heme_mol * (0.0197913 / 0.0014711)
        M701.ins[1].imass['H2O'] = GammaCyclodextrinFeed.imass['GammaCyclodextrin'] * 10
    
    @M703.add_specification(run=True)
    def dose_nicotinamide():
        heme_mol = H601.outs[0].imol['Heme_b'] + H601.outs[0].imol['Heme_b_In']
        nic_mol = heme_mol * 2.0 * 1.025
        NicotinamideFeed.imol['Nicotinamide'] = nic_mol
        M703.ins[1].imass['H2O'] = NicotinamideFeed.imass['Nicotinamide'] * 5

    # =========================================================================
    # AREA 800: FINAL PRODUCT
    # =========================================================================
    
    H801 = bst.HXutility('H801', ins=R702-0, outs=NHemDx_Product, T=4+273.15, cool_only=True)

    # =========================================================================
    # AREA 900: WASTEWATER TREATMENT & FACILITIES (Following LegHb config1)
    # =========================================================================
    
    # Create wastewater treatment system following LegHb config1 pattern
    # Collects all waste streams and produces biogas, sludge, RO_treated_water
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[
            U501-0,    # ResinFlowthrough
            U501-2,    # ResinRegenWaste
            U501-3,    # ResinWash
            U601-1,    # NFPermeate
            C402-1,    # WashEffluent
            S406-1,    # PressLiquor (from ScrewPress)
        ],
        outs=('biogas', 'sludge', 'RO_treated_water', ProcessWaste),
        mockup=True,
        area=500,
    )
    biogas, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
    
    # Utility systems (following LegHb pattern)
    CT = bst.CoolingTower('CT')
    CWP = bst.ChilledWaterPackage('CWP')
    
    # Combine dehydrated debris with WWT sludge for boiler
    M902 = bst.Mixer(
        'M902',
        ins=(S406-0, sludge),  # DehydratedDebris + WWT sludge
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
    
    # Process Water Center (following LegHb pattern)
    F = bst.main_flowsheet
    makeup_water_streams = (
        F.cooling_tower_makeup_water,
        F.Water5, F.Water6, F.Water7,  # Buffer waters
        F.Water8, F.Water9,  # Formulation waters
        F.DFBufferWater,
        F.boiler_makeup_water,
    )
    process_water_streams = (
        treated_water,
        F.rejected_water_and_blowdown,
        *makeup_water_streams
    )
    
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    
    PWC = bst.ProcessWaterCenter(
        'PWC',
        ins=(treated_water, makeup_water, 'recycled_process_water', 'makeup_process_water'),
        outs=('RO_water', 'process_water', 'excess_water'),
        makeup_water_streams=makeup_water_streams,
        process_water_streams=process_water_streams,
        reverse_osmosis_water_price=0.000254,
        process_water_price=0.000135,
    )



# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == '__main__':
    bst.preferences.N = 50
    
    print("="*85)
    print("N-HemDx FULL PRODUCTION SYSTEM - Upstream + DSP Integration")
    print("="*85)
    
    print("\n1. Creating system...")
    NHemDx_sys = create_NHemDx_system()
    sys = NHemDx_sys
    f = sys.flowsheet
    u = f.unit
    ss = f.stream
    sys.operating_hours = 8000
    
    print("\n2. Running simulation...")
    try:
        sys.simulate()
        print("   Simulation completed successfully!")
    except Exception as e:
        print(f"   Simulation error: {e}")
        raise
    
    print(f"\n4. System Summary")
    print("="*85)
    sys.show()    
    
    print("\n5. Generating system diagram...")
    try:
        sys.diagram(format='html', display=True)
    except Exception as e:
        print(f"   Diagram error: {e}")
    
    print("\n" + "="*85)
    print("SIMULATION COMPLETE")
    print("="*85)

    print("\n3. Product Composition:")
    product = ss.NHemDx_Product
    n_hemdx = product.imass['N-HemoDextrin']
    hemdx = product.imass['HemoDextrin']
    heme = product.imass['Heme_b'] + product.imass['Heme_b_In']
    print(f"    N-HemoDextrin: {n_hemdx:.4f} kg/hr")
    print(f"    HemoDextrin:   {hemdx:.4f} kg/hr")
    print(f"    Free Heme:     {heme:.4f} kg/hr")
    print(f"    Total Product: {product.F_mass:.2f} kg/hr")

    hemdx_mol = product.imol['HemoDextrin']
    n_hemdx_mol = product.imol['N-HemoDextrin']
    heme_equiv_mol = (hemdx_mol + n_hemdx_mol) * 0.0014711
    heme_equiv_mass = heme_equiv_mol * HEMDX_THERMO['Heme_b'].MW
    print("\n4. Hemodextrin Equivalents:")
    print(f"    HemoDextrin:   {hemdx_mol:.6f} kmol/hr ({hemdx:.4f} kg/hr)")
    print(f"    N-HemoDextrin: {n_hemdx_mol:.6f} kmol/hr ({n_hemdx:.4f} kg/hr)")
    print(f"    Heme equiv.:   {heme_equiv_mol:.6f} kmol/hr ({heme_equiv_mass:.4f} kg/hr)")
    
    print(f"\n{'='*85}")
    print("KEY PERFORMANCE INDICATORS")
    print(f"{'='*85}")
    print(f"  Product Stream:      {product.ID}")
    print(f"  Production Rate:     {product.F_mass:.2f} kg/hr")
    print(f"  Annual Production:   {product.F_mass * sys.operating_hours / 1000:.2f} metric tons/year")
    print(f"{'='*85}\n")
    
    u.U501.show()