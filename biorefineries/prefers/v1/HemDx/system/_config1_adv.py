# -*- coding: utf-8 -*-
"""
N-HemDx Full Production System - Config 1 ADVANCED

This configuration overrides multiple areas with advanced units:
- Area 400 (Recovery): FiltrationAdv for S404, S405
- Area 500 (Purification): ResinColumnAdv for U501
- Area 600 (Concentration): DiafiltrationAdv for U601
"""

from biorefineries.prefers.v1.HemDx.system import _config1
from biorefineries.prefers.v1.HemDx.system._config1 import *
from biorefineries.prefers.v1.HemDx.system._config1 import (
    create_area_200_media_prep,
    create_area_300_conversion,
    create_area_400_recovery as _create_area_400_recovery_base,
    create_area_500_purification as _create_area_500_purification_base,
    create_area_600_concentration as _create_area_600_concentration_base,
    create_area_700_formulation,
    create_area_800_final_product,
    create_area_900_facilities,
    get_fermentation_parameters,
    create_fermentation_reactions,
    check_HemDx_specifications,
    optimize_NH3_loading,
    set_production_rate,
    load_process_settings,
    set_GWPCF, set_GWPCF_Multi,
    HEMDX_THERMO,
    s, c, u, bst, tmo, np
)

# Import new units
from biorefineries.prefers.v1._units_adv import ResinColumnAdv, DiafiltrationAdv, FiltrationAdv

# =============================================================================
# AREA 400: RECOVERY (OVERRIDE)
# =============================================================================
def create_area_400_recovery(broth_in, DfUltraBuffer2):
    """
    Create Area 400 (Recovery) using the advanced FiltrationAdv unit for S404/S405.
    
    Replaces base Filtration units with FiltrationAdv (backward mode).
    Maintains all other base units (Centrifuges, CellDisruption, etc.)
    """
    # S401: Primary Centrifuge (from base)
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

    # S404: Supernatant microfiltration (ADVANCED)
    S404 = FiltrationAdv.from_preset(
        'MF', 'S404',
        ins=S401-1,
        outs=('SupernatantCake', 'FilteredSupernatant'),
        solid_capture_efficiency=0.95,
        cake_moisture_content=0.30,
        parameter_mode='backward',
    )

    # Cell cream washing & disruption
    M401 = bst.MixTank('M401', ins=(DfUltraBuffer2, 'Water4'), outs='WashBufferOut', tau=0.5)

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

    # S405: Lysate clarification (ADVANCED)
    S405 = FiltrationAdv.from_preset(
        'MF', 'S405',
        ins=S403-1,
        outs=('FilterCake', 'ClarifiedLysate'),
        solid_capture_efficiency=0.95,
        cake_moisture_content=0.30,
        parameter_mode='backward',
    )

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
# AREA 500: PURIFICATION (OVERRIDE)
# =============================================================================
def create_area_500_purification(lysate_in, NaCl_wash, NaOH_elute, Ethanol_regen):
    """
    Create Area 500 (Purification) using the advanced ResinColumnAdv unit.
    
    Uses 'ratio' parameter mode: Performance ratios are set directly by the user.
    CVs are used for buffer sizing but do not affect performance ratios.
    """
    # Re-create streams/units locally as in original
    M501 = bst.MixTank('M501', ins=(NaCl_wash, 'Water5'), outs='WashBuffer', tau=1)
    M502 = bst.MixTank('M502', ins=(NaOH_elute, 'Water6'), outs='ElutionBuffer', tau=1)
    M503 = bst.MixTank('M503', ins=(Ethanol_regen, 'Water7'), outs='RegenBuffer', tau=1)

    H501 = bst.HXutility('H501', ins=M501-0, outs='H501Out', T=30+273.15, heat_only=True)
    H502 = bst.HXutility('H502', ins=M502-0, outs='H502Out', T=40+273.15, heat_only=True)
    H503 = bst.HXutility('H503', ins=M503-0, outs='H503Out', T=30+273.15, heat_only=True)
    
    # Use ResinColumnAdv with 'ratio' mode
    # Performance ratios are set directly; CVs are for buffer sizing only.
    U501 = ResinColumnAdv(
        'U501',
        ins=(lysate_in, H501-0, H502-0, H503-0),
        outs=('ResinFlowthrough', 'ResinEluate', 'ResinRegenWaste', 'ResinWash'),
        preset='Adsorption',
        TargetProduct_IDs=('Heme_b', 'Heme_b_In', 'ProtoporphyrinIX', 'ProtoporphyrinIX_In'),
        # --- Performance Ratios (Design Basis) ---
        TargetProduct_Yield=0.99,           # 99% yield
        Wash_Impurity_Carryover=0.02,       # 2% impurity to wash
        Regen_Impurity_Carryover=0.01,      # 1% impurity to regen

        # NonTarget_Removal=0.99, # Removed in standalone refactor             # 99% non-targets to flowthrough (legacy, used by base class)
        # --- CVs (for buffer sizing specifications) ---
        wash_CV=3.0,
        elution_CV=0.05, 
        regeneration_CV=0.05,
        # --- Advanced Parameters ---
        parameter_mode='ratio',             # Use ratios directly
        efficiency=1.0,                     # Efficiency already factored into ratios above
        operating_pressure_drop_bar=2.0
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
# AREA 600: CONCENTRATION (OVERRIDE)
# =============================================================================
def create_area_600_concentration(eluate_in, DfUltraBuffer1):
    """
    Create Area 600 (Concentration) using the advanced DiafiltrationAdv unit.
    
    Uses 'backward' parameter mode: Retention targets are set directly.
    Matches baseline U601 performance (NF preset, 98% target retention).
    """
    M601 = bst.MixTank('M601', ins=(DfUltraBuffer1, 'Water8'), outs='DFBuffer', tau=0.5)
    
    U601 = DiafiltrationAdv.from_preset(
        'NF', 'U601',
        ins=(eluate_in, M601-0),
        outs=('ConcentratedHeme', 'NFPermeate'),
        TargetProduct_IDs=('Heme_b', 'Heme_b_In', 'ProtoporphyrinIX', 'ProtoporphyrinIX_In'),
        Salt_IDs=('NaCl', 'NaOH', 'Ethanol', 'K2SO4', 'MgSO4', 'FeSO4', '(NH4)2SO4', 'KH2PO4'),
        # --- Backward Mode: Set retention targets directly ---
        parameter_mode='backward',
        TargetProduct_Retention=0.98,  # Match baseline
        Salt_Retention=0.10,           # Match baseline
        OtherLargeMolecules_Retention=0.98,
        DefaultSolutes_Retention=0.10,
        # --- Operating params ---
        diavolumes=5.0,
        FeedWater_Recovery_to_Permeate=0.80,
        efficiency=1.0,  # No scaling since retention is direct target
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
# SYSTEM FACTORY (OVERRIDE)
# =============================================================================
@bst.SystemFactory(
    ID='NHemDx_sys_adv',
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
    
    # Set GWPCF (Same as original)
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
    broth = create_area_300_conversion(A200['M203_out'], A200['T201_out'], A200['T202_out'], vent1, vent2, rxns, params)
    
    # Area 400: Recovery
    area_400 = create_area_400_recovery(broth, DfUltraBuffer2)
    
    # Area 500: Purification (USING UPGRADED VERSION)
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
    print("N-HemDx FULL PRODUCTION SYSTEM - ADVANCED (ResinColumn Upgrade)")
    print("="*85)
    
    print("\n1. Creating system...")
    NHemDx_sys = create_NHemDx_system()
    sys = NHemDx_sys
    sys.operating_hours = 8000
    
    print("\n2. Running baseline simulation...")
    try:
        sys.simulate()
        baseline_production = sys.flowsheet.stream.NHemDx_Product.F_mass
        print(f"   Baseline production rate: {baseline_production:.2f} kg/hr")
        
        # Check ResinColumnAdv status
        U501 = sys.flowsheet.unit.U501
        print(f"\n   ResinColumnAdv Status:")
        print(f"     Efficiency: {U501.efficiency}")
        print(f"     Operating dP: {U501.operating_pressure_drop_bar} bar")
        print(f"     Parameter Mode: {U501.parameter_mode}")
        print(f"     TargetProduct_Yield: {U501.TargetProduct_Yield:.4f}")
        print(f"     Wash_Impurity_Carryover: {U501.Wash_Impurity_Carryover:.4f}")
        print(f"     Regen_Impurity_Carryover: {U501.Regen_Impurity_Carryover:.4f}")
        
    except Exception as e:
        print(f"   Baseline simulation failed: {e}")
        import traceback
        traceback.print_exc()
        raise    
    
    print("\n3. Setting production rate...")
    try:
        set_production_rate(sys, TARGET_PRODUCTION)
    except Exception as e:
        print(f"   Production rate setting error: {e}")

    print("\n4. Product Analysis:")
    product = sys.flowsheet.stream.NHemDx_Product
    check_HemDx_specifications(product)
    
    print(f"\n{'='*85}")
    print("SIMULATION COMPLETE")
    print(f"{'='*85}\n")

    u = sys.flowsheet.unit
    
    print(f"\n--- S404 (FiltrationAdv) Final State ---")
    print(f"  Solid Capture Efficiency:   {u.S404.solid_capture_efficiency:.4f}")
    print(f"  Cake Moisture:              {u.S404.cake_moisture_content:.4f}")
    print(f"  Parameter Mode:             {u.S404.parameter_mode}")
    
    print(f"\n--- S405 (FiltrationAdv) Final State ---")
    print(f"  Solid Capture Efficiency:   {u.S405.solid_capture_efficiency:.4f}")
    print(f"  Cake Moisture:              {u.S405.cake_moisture_content:.4f}")
    print(f"  Parameter Mode:             {u.S405.parameter_mode}")
    
    print(f"\n--- U501 (ResinColumnAdv) Final State ---")
    print(f"  TargetProduct_Yield:        {u.U501.TargetProduct_Yield:.4f}")
    print(f"  Wash_Impurity_Carryover:    {u.U501.Wash_Impurity_Carryover:.4f}")
    print(f"  Regen_Impurity_Carryover:   {u.U501.Regen_Impurity_Carryover:.4f}")
    print(f"  Parameter Mode:             {u.U501.parameter_mode}")
    
    print(f"\n--- U601 (DiafiltrationAdv) Final State ---")
    print(f"  TargetProduct_Retention:    {u.U601.TargetProduct_Retention:.4f}")
    print(f"  Salt_Retention:             {u.U601.Salt_Retention:.4f}")
    print(f"  Diavolumes:                 {u.U601.diavolumes:.2f}")
    print(f"  Parameter Mode:             {u.U601.parameter_mode}")

