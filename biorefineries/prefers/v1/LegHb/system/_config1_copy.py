# -*- coding: utf-8 -*-
"""
Created on 2025-06-04 14:26:14

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

#from biorefineries.animal_bedding import _system
from ast import Yield
from biorefineries.prefers.v1._process_settings import set_GWPCF, GWP_CFs, set_GWPCF_Multi, load_process_settings
from biorefineries.prefers.v1.LegHb import _streams as s
import biosteam as bst
from pint import set_application_registry
from thermosteam import Stream
from biosteam import F
import thermosteam as tmo
import numpy as np
from biorefineries.prefers.v1.LegHb import _chemicals as c
from biorefineries.prefers.v1 import _units as u
from biorefineries.prefers.v1._process_settings import price

# %% Settings
LEGHB_THERMO = c.create_chemicals_LegHb()
bst.settings.set_thermo(LEGHB_THERMO, skip_checks=True)
bst.preferences.classic_mode()
SHOW_RXN = False

# %%
__all__ = (
    'create_LegHb_system',
    'set_production_rate',
    'check_LegHb_specifications',
)
# %%
@bst.SystemFactory(
    ID='LegHb_sys',
    ins=[s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt, s.DfUltraBuffer, s.DfUltraBuffer_wash],
    outs=[s.LegHb_3, s.vent1, s.vent2, s.effluent1],
    fthermo=lambda chemicals=None: LEGHB_THERMO,  # Force use of cached thermo instance
)
def create_LegHb_system(
        ins, outs,
        use_area_convention=False,
        # reactions_passed=None, # Placeholder if reactions were to be passed externally
        # One could add more configurable parameters here:
        # V_max_fermenter=500, target_titer=7.27, etc.
    ):
    """
    Creates the LegHb (Leghemoglobin) production system.
    This system is based on the process flow and parameters from test_LegHb.py.
    """
    bst.preferences.N=50
    # Update all input stream prices before system creation to avoid conflicts
    # This is now done at the module level in _streams.py to prevent ID conflicts
    
    # Uncomment the line below if you need to update prices dynamically:
    #s.update_all_input_stream_prices()
    # Unpack input streams
    SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, DfUltraBuffer_wash = ins

    # Unpack output streams
    (LegHb_3, vent1, vent2, effluent1) = outs
    
    set_GWPCF(Glucose, 'Glucose')
    set_GWPCF_Multi(SeedIn1, ['AmmoniumSulfate','Glucose','MagnesiumSulfate','KH2PO4'],[0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(SeedIn2, ['AmmoniumSulfate','Glucose','MagnesiumSulfate','KH2PO4'],[0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(CultureIn, ['Glycine','Glucose','IronSulfate'],[0.1/(0.1+60+0.15191), 60/(0.1+60+0.15191), 0.15191/(0.1+60+0.15191)])
    set_GWPCF(NH3_25wt, 'Ammonia_SEA',dilution=0.25)
    set_GWPCF_Multi(DfUltraBuffer, ['KH2PO4','NaCl'], [0.8472, 0.1455])
    set_GWPCF_Multi(DfUltraBuffer_wash, ['KH2PO4','NaCl'], [0.8472, 0.1455])

    load_process_settings()  # Load process settings to update prices and CFs
    """
    Fermentation Parameters
    -----------------------
    Validated against literature [Yao et al. 2025, World J Microbiol Biotechnol 41:404]
    
    Key findings from literature:
    - LegHb titer: Up to 7.27 g/L achieved in K. marxianus [Tian et al. 2024b]
    - Optimal temperature: 28-30°C for K. phaffii, 30°C for S. cerevisiae
    - Optimal pH: 5.0-5.5 for yeast hosts
    - Fed-batch with μ-STAT feeding strategy increases yields 67.6% vs DO-STAT
    - Heme supply is limiting factor; overexpression of Hem pathway genes beneficial
    - Cell growth yield Y_b: ~0.43 g/g glucose (typical for yeast)
    """
    theta_O2 = 0.5  # Dissolved oxygen concentration [% saturation]
    agitation_power = 0.985  # [kW/m³] - Ref: Typical for STR bioreactors
    design = 'Stirred tank'  # Reactor type
    method = "Riet"  # kLa correlation method
    T_operation = 273.15 + 32  # [K] - Ref: 28-32°C optimal [1]
    Q_O2_consumption = -110 * 4184  # [kJ/kmol] - Heat of aerobic metabolism
    dT_hx_loop = 8  # [°C] - Heat exchanger approach temperature
    cooler_pressure_drop = 20684  # [Pa]
    compressor_isentropic_efficiency = 0.85
    V_max = 500  # [m³] - Max vessel volume
    
    # LegHb production parameters - Ref: Yao et al. 2025 [1]
    # Titer 7.27 g/L achieved in K. marxianus (Tian et al. 2024b)
    titer_LegHb = 7.27  # [g/L] - Literature validated
    productivity_LegHb = 7.27 / 72  # [g/L/h] - 72h fermentation cycle
    Y_p = 7.27 * 4 / 1300  # [by wt] - ~2.2% product yield on glucose
    Y_b = 0.43  # [by wt] - Biomass yield, typical for yeast
    yield_LegHb = Y_p  # Product yield based on glucose for product formation
    
    # Heme B co-product (intracellular)
    titer_HemeB = 0.06632  # [g/L] - Heme accumulation
    productivity_HemeB = 0.06632 / 72  # [g/L/h]


    """
    Reactions
    """
    
    # fermentation_reaction = bst.PRxn([
    #     #           Reaction        Reactnat            Conversion           Check                  "
    #     bst.Rxn('8 Glucose + 6 NH3 + 1 FeSO4 + 10.5 O2 -> Heme_b + 1 (NH4)2SO4 + 37 H2O + 14 CO2',
    #                                 reactant = 'Glucose',X=LegHb_yield*0.05,check_atomic_balance=True),
    #     bst.Rxn('Glucose + (NH4)2SO4 + NH3 -> Globin + CO2 + H2O',
    #                                 reactant = 'Glucose', X= LegHb_yield*0.05,correct_atomic_balance=True),
    #     bst.Rxn('Glucose + FeSO4 + (NH4)2SO4 + NH3 -> Leghemoglobin + CO2  + H2O',
    #                                 reactant = 'Glucose', X=LegHb_yield,  correct_atomic_balance=True),
    #     ])
    fermentation_reaction = bst.PRxn([
        #           Reaction        Reactnat            Conversion           Check                  "
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                                    reactant = 'Glucose',X=yield_LegHb*0.01,check_atomic_balance=True),
        bst.Rxn('Glucose + 0.01646 (NH4)2SO4 + 1.61317 NH3 -> 6 Globin_In + 0.28807 O2 + 3.68724 H2O',
                                    reactant = 'Glucose', X= yield_LegHb*0.04,check_atomic_balance=True),
        bst.Rxn('Glucose + 0.00786 FeSO4 + 0.00786 (NH4)2SO4 + 1.58847 NH3 -> 6 Leghemoglobin_In + 0.30275 O2  + 3.70380 H2O',
                                    reactant = 'Glucose', X=yield_LegHb,  check_atomic_balance=True),
        ])
    fermentation_reaction[2].product_yield('Leghemoglobin_In', basis='wt', product_yield=yield_LegHb)

    neutralization_reaction = bst.Rxn(
        'H2SO4 + 2 NH3 -> (NH4)2SO4', reactant = 'NH3', X=1,
        check_atomic_balance=True
    )

    cell_growth_reaction = bst.Rxn(
        'Glucose + 0.8364 NH3 + 0.0108 (NH4)2SO4 -> 2.01 H2O + 0.106 O2 + 6 Pichia_pastoris', 'Glucose', X=(1-yield_LegHb*1.1)*Y_b,
        correct_atomic_balance=True
    )
    cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=(1-yield_LegHb*1.1)*Y_b)

    respiration_reaction1 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', 1 - Y_b,
        check_atomic_balance=True
    )

    respiration_reaction2 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', 1 - cell_growth_reaction.X - fermentation_reaction[2].X * 1.1,
        check_atomic_balance=True
    )

    bst.settings.chemicals.set_alias('Pichia_pastoris', 'cellmass')
    RXN = bst.ReactionSystem(
        fermentation_reaction,
        bst.PRxn([cell_growth_reaction, respiration_reaction2])
    )
    if SHOW_RXN:
        RXN.show()

    # =========================================================================
    # AREA 200: MEDIA PREPARATION
    # =========================================================================
    """
    Media Preparation Area
    
    Purpose: Prepare all feed solutions for fermentation
    
    Units:
    - M301: Seed solution 1 preparation (mixing seed with water)
    - M302: Seed solution 2 + culture media preparation
    - M303: Seed hold tank (combines both seed solutions)
    - M304: Glucose solution preparation (50% dilution with water)
    - T301: Glucose storage tank (96 hours residence time)
    - T302: Ammonia storage tank (25 wt% aqueous ammonia)
    
    Outputs:
    - M303Out: Combined seed solution → R301 (Seed Train)
    - T301Out: Diluted glucose solution → R302 (Main Fermentation)
    - T302Out: Ammonia solution → R302 (pH control)
    """
    
    M301 = bst.MixTank('M301', ins=[SeedIn1,'Water1'], outs='M301Out', tau=16)
    @M301.add_specification(run=True)
    def update_seed1_inputs():
        target_stream = bst.Stream(**{**s.SeedSolution1, 'ID': None})
        SeedIn1.imass['Seed'] = target_stream.imass['Seed']
        M301.ins[1].imass['H2O'] = target_stream.imass['H2O']
        M301.ins[1].T = 25+273.15
    
    M302 = bst.MixTank('M302', ins=[SeedIn2,CultureIn,'Water2'], outs='M302Out', tau=16)
    @M302.add_specification(run=True)
    def update_culture_inputs():
        target_stream = bst.Stream(**{**s.SeedSolution2, 'ID': None})
        SeedIn2.imass['Seed'] = target_stream.imass['Seed']
        M302.ins[2].imass['H2O'] = target_stream.imass['H2O']
        M302.ins[2].T = 25+273.15
        CultureIn.imass['Culture'] = target_stream.imass['SeedSolution']*(0.1+60+0.15191)/1000

    M303 = u.SeedHoldTank('M303', ins=[M301-0, M302-0], outs='M303Out')

    # =========================================================================
    # AREA 300: FERMENTATION / CONVERSION
    # =========================================================================
    """
    Fermentation/Conversion Area
    
    Purpose: Convert glucose to LegHb through microbial fermentation
    
    Process Flow:
    1. Seed Train (R301): Grow inoculum culture (28-32°C, pH 5.0-5.5)
    2. Main Fermentation (R302): Batch aerobic fermentation (72h cycle)
       - Target titer: 7.27 g/L LegHb [Yao et al. 2025]
       - Product yield: ~2.2% on glucose
       - Cell yield: 0.43 g/g glucose (typical yeast)
    
    Units:
    - R301: Seed train bioreactor
    - R302: Main aerobic fermentation vessel (500 m³ max)
    
    Outputs:
    - vent1: CO2 exhaust from seed train
    - vent2: CO2 exhaust from main fermentation
    - Broth: Cell suspension containing intracellular LegHb → Recovery
    """
    
    R301 = u.SeedTrain(
        'R301',
        ins=[M303-0],
        outs=[vent1, 'R301Out'],
        reactions=bst.PRxn([cell_growth_reaction, respiration_reaction1]),
        saccharification=None,
        T=32+273.15,
    )
    R301.add_specification(run=True)
    
    M304 = bst.MixTank('M304', ins=[Glucose,'Water3'], outs='M304Out', tau=16)
    
    @M304.add_specification(run=True)
    def update_water_content():
        M304.ins[1].imass['H2O'] = Glucose.imass['Glucose']/2
        M304.ins[1].T = 25+273.15
    
    T301 = bst.StorageTank('T301', ins=M304-0, outs='T301Out',tau=16*4+72)

    T302 = u.AmmoniaStorageTank('T302', ins=NH3_25wt, outs='T302Out')

    # checking more details...
    R302 = u.AeratedFermentation(
        'R302',
        ins=[R301-1, T301-0, T302-0, bst.Stream('FilteredAir', phase='g', P=2 * 101325)],
        outs=[vent2, 'Broth'],
        fermentation_reaction=fermentation_reaction,
        cell_growth_reaction=cell_growth_reaction,
        respiration_reaction=respiration_reaction2,
        neutralization_reaction=neutralization_reaction,
        design=design, method=method,theta_O2=theta_O2,
        V_max=V_max, Q_O2_consumption=Q_O2_consumption,
        dT_hx_loop=dT_hx_loop, T=T_operation,
        batch=True, reactions=RXN,
        kW_per_m3=agitation_power,
        tau=titer_LegHb/productivity_LegHb,
        cooler_pressure_drop=cooler_pressure_drop,
        compressor_isentropic_efficiency=compressor_isentropic_efficiency,
        P=1 * 101325,#optimize_power=True,
    )
    R302.target_titer = titer_LegHb # g / L
    R302.target_productivity = productivity_LegHb # g / L / h
    R302.target_yield = yield_LegHb  # wt %

    @R302.add_specification(run=True)
    def update_reaction_time_and_yield():
        R302.tau = R302.target_titer / R302.target_productivity
        fermentation_reaction[2].product_yield('Leghemoglobin_In', basis='wt', product_yield=R302.target_yield)


    # =========================================================================
    # AREA 400: RECOVERY (BIOMASS HARVEST & CELL DISRUPTION)
    # =========================================================================
    """
    Recovery Area - Intracellular LegHb Extraction
    
    Purpose: Release and clarify intracellular LegHb from yeast cells
    
    Process Flow:
    1. Biomass Harvest (C401): Centrifuge to concentrate cells (45-55% WCW)
    2. Cell Washing (M401, H401_wash, M402_wash, C402): Remove spent media
    3. Cell Disruption (S401): High-pressure homogenization to lyse cells
    4. Debris Removal (H401, S402): Cool and centrifuge cell debris
    5. Lysate Clarification (S403, S404): Depth filtration + screw press
    
    Units:
    - C401: Primary harvest centrifuge (98% cell capture)
    - M401, H401_wash, M402_wash: Cell wash system
    - C402: Wash centrifuge
    - S401: Cell disruption homogenizer (700-1200 bar)
    - H401: Post-disruption cooler (to 10°C)
    - S402: Debris centrifuge
    - S403: Depth filtration
    - S404: Screw press for dewatering debris
    
    Outputs:
    - ClarifiedLysate: Cell-free protein solution → Purification (Area 500)
    - SpentMedia, WashEffluent: To wastewater (Area 900)
    - DehydratedDebris: To disposal (Area 900)
    """
    
    # Biomass Harvest (Primary Separation)
    C401 = bst.SolidsCentrifuge(
        'C401',
        ins=R302-1,  # Fermentation broth
        outs=('CellCream', 'SpentMedia'),
        split={'cellmass': 0.98,  # High cell capture
               'Leghemoglobin_In': 0.98,  # Intracellular product stays with cells
               'Globin_In': 0.98,  # Intracellular protein stays with cells
               'Heme_b': 0.85,  # Heme mostly with cells
               'Glucan': 0.95,
               'Mannoprotein': 0.95,
               'Chitin': 0.95,
               'OleicAcid': 0.90,
               'RNA': 0.90,
               },  # Cells + intracellular product to cream
        moisture_content=0.55,  # ~45% dry solids = 55% moisture
    )
    
    # Wash buffer preparation for cell wash (1-2 diavolumes)
    M401 = bst.MixTank('M401', ins=(DfUltraBuffer_wash, 'WashWater'), 
                    outs='WashBufferOut', tau=0.5)
    
    @M401.add_specification(run=True)
    def update_wash_buffer():
        # 1.5 diavolumes of wash buffer
        M401.ins[1].imass['H2O'] = C401.outs[0].F_mass * 1.5
        M401.ins[0].imol['DfUltraBuffer'] = M401.ins[1].imass['H2O'] * (0.025+0.01+0.001)/1000
    
    # Cool wash buffer before use
    H401_wash = bst.HXutility(
        'H401_wash',
        ins=M401-0,
        outs='ColdWashBuffer',
        T=10+273.15,  # Cool to 10°C
        cool_only=True,
    )
    
    # Cell wash mixer - combines cell cream with wash buffer
    M402_wash = bst.MixTank('M402_wash', ins=(C401-0, H401_wash-0), 
                            outs='WashedCellSlurry', tau=0.25)
    
    # Second centrifuge for washed cell separation
    C402 = bst.SolidsCentrifuge(
        'C402', 
        ins=M402_wash-0,
        outs=('WashedCellCream', 'WashEffluent'),
        split={'cellmass': 0.98,
               'Leghemoglobin_In': 0.98,  # Intracellular product stays with cells
               'Globin_In': 0.98,  # Intracellular protein stays with cells
               'Heme_b': 0.85,
               'Glucan': 0.95,
               'Mannoprotein': 0.95,
               },
        moisture_content=0.55,
    )
    
    # =========================================================================
    # STEP 2: Cell Disruption (High-Pressure Homogenization)
    # =========================================================================
    # Multi-pass HPH at 800-1200 bar to achieve >90% cell lysis
    # Pressure generates ~2.5°C per 100 bar - requires cooling
    
    S401 = u.CellDisruption(
        'S401',
        ins=C402-0,  # Washed cell cream (not raw broth)
        outs='CrudeHomogenate',
        P_high=1000e5,  # 1000 bar = 100 MPa (spec: 800-1200 bar)
        cell_disruption_efficiency=0.90,  # >90% lysis per spec
    )
    
    # Immediate post-valve cooling to prevent LegHb denaturation
    H401 = bst.HXutility(
        'H401',
        ins=S401-0,
        outs='CooledHomogenate',
        T=15+273.15,  # Return to <15°C per spec
        cool_only=True,
    )
    
    # =========================================================================
    # STEP 3: Lysate Clarification (Debris Removal)
    # =========================================================================
    # High-speed centrifugation followed by depth filtration
    # Goal: Remove cell wall ghosts, nucleic acids, lipids
    
    # Primary clarification - high-speed centrifuge for bulk solids
    # Key: Remove as much Mannoprotein as possible to achieve ≥65% protein purity
    S402 = bst.SolidsCentrifuge(
        'S402',
        ins=H401-0,  # Cooled homogenate
        outs=('CellDebris', 'CrudeLysate'),
        split={'cellmass': 0.98,  # Undisrupted cells to debris (increased)
               'Glucan': 0.98,     # Cell wall components (increased)
               'Chitin': 0.98,     # Cell wall (increased)
               'OleicAcid': 0.98,  # Lipids mostly removed with debris (increased)
               'RNA': 0.85,        # Pelleted with debris (increased)
               'Leghemoglobin': 0.02,  # Minimal loss of product
               'Globin': 0.05,
               'Mannoprotein': 0.98,  # Cell wall protein - critical for purity!
               },
        moisture_content=0.20,  # Relatively dry sludge
    )
    
    # Secondary clarification - depth filtration (series: 30µm → 5µm → 0.5µm)
    # Modeled as MF to remove colloidal haze
    S403 = u.Filtration.from_preset(
        'MF',
        'S403',
        ins=S402-1,  # Crude lysate
        outs=('FilterCake', 'ClarifiedLysate'),
        solid_capture_efficiency=0.95,  # Capture remaining debris
        cake_moisture_content=0.30,
    )
    S403.add_specification(run=True)
    
    # Debris disposal via screw press for volume reduction
    S404 = bst.ScrewPress(
        'S404', 
        ins=(S402-0, S403-0),  # Combined debris streams
        outs=('DehydratedDebris', 'PressLiquor'),
        split=0.999, 
        moisture_content=0.20
    )
    
    # =========================================================================
    # AREA 500: PURIFICATION (UF/DF MEMBRANE SEPARATION)
    # =========================================================================
    """
    Purification Area - Membrane-based Protein Concentration & Buffer Exchange
    
    Purpose: Concentrate LegHb and remove impurities using tangential flow filtration
    
    Process Flow:
    1. DF Buffer Preparation (M403, H402): Prepare cold buffer for diafiltration
    2. UF/DF Stage 1 (U401): Concentrate + buffer exchange (5-7 diavolumes)
       - MWCO: 3-10 kDa (LegHb monomer = 16 kDa, high retention)
       - 99% LegHb retention, 95% salt removal
    3. UF/DF Stage 2 (U404): Final concentration to target protein content
       - 85% water recovery to permeate
       - Concentrated LegHb solution → Formulation (Area 600)
    
    Units:
    - M403: DF buffer mixing tank
    - H402: Buffer cooler (5°C for protein stability)
    - U401: Primary UF/DF unit (concentration + buffer exchange)
    - U404: Secondary UF unit (final concentration)
    
    Outputs:
    - ConcentratedLegH: Purified, concentrated LegHb solution → Formulation
    - UFPermeate, ConcentrationPermeate: To wastewater (Area 900)
    
    Config1 Pathway: Diafiltration only (no chromatography)
    Config2 Pathway: Adds ion exchange chromatography (ResinColumn U402)
    """
    
    # TFF with 3-10 kDa MWCO membrane
    # LegH monomer is ~16 kDa - high retention expected
    # VCF 5-10X, DF 5-7 diavolumes
    
    # DF buffer preparation with antioxidant
    M403 = bst.MixTank('M403', ins=(DfUltraBuffer, 'DFBufferWater'), 
                    outs='DFBufferOut', tau=0.5)
    
    @M403.add_specification(run=True)
    def update_df_buffer():
        # 6 diavolumes of DF buffer (middle of 5-7 range)
        feed_water = S403.outs[1].imass['H2O']
        M403.ins[1].imass['H2O'] = feed_water * 6
        M403.ins[0].imol['DfUltraBuffer'] = M403.ins[1].imass['H2O'] * (0.025+0.01+0.001)/1000
    
    # Cool DF buffer
    H402 = bst.HXutility(
        'H402',
        ins=M403-0,
        outs='ColdDFBuffer',
        T=5+273.15,  # Cold operation to preserve protein
        cool_only=True,
    )
    
    # UF/DF unit - combines concentration and buffer exchange
    U401 = u.Diafiltration.from_preset(
        'UF',
        'U401',
        ins=(S403-1, H402-0),  # Clarified lysate + DF buffer
        outs=('UFConcentrate', 'UFPermeate'),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
        TargetProduct_Retention=0.99,  # High LegH retention
        Salt_Retention=0.05,  # Salts wash through
    )
    U401.add_specification(run=True)
    
    # Final concentration to target protein content
    U404 = u.Diafiltration(
        'U404',
        ins=(U401-0, bst.Stream('U404_buffer', H2O=0.001)),  # Minimal buffer for DF interface
        outs=('ConcentratedLegH', 'ConcentrationPermeate'),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
        TMP_bar1=3,
        FeedWater_Recovery_to_Permeate=0.85,  # Moderate concentration
    )
    
    @U404.add_specification(run=True)
    def U404_adjust_water_recovery():
        """
        Concentrate the product via UF. The final dilution to meet total solids
        specification is handled in the M404 formulation unit.
        """
        U404.Salt_Retention = 0.95
        U404._run()
    
    # =========================================================================
    # AREA 600: FORMULATION (THERMAL TREATMENT & FINAL PRODUCT)
    # =========================================================================
    """
    Formulation Area - Thermal Stabilization & Product Finishing
    
    Purpose: Ensure microbiological safety and achieve final product specification
    
    Process Flow:
    1. HTST Pasteurization (H403, T401): Heat to 72°C for 30 seconds
       - Inactivates pathogens
       - Precipitates unstable host proteins
       - LegHb (stable heme-protein) remains soluble
    2. Formulation (M404): Add antioxidant + dilution water
       - Target: 20% total solids (below 24% limit)
       - Target: 7.5% LegHb (middle of 6-9% range)
       - Antioxidant: Sodium ascorbate for heme stability
    3. Final Cooling (H406): Rapid chill to 4°C for storage
    
    Units:
    - H403: Pasteurization heater (72°C)
    - T401: Hold time tank (30 seconds residence)
    - M404: Formulation mixer (antioxidant + dilution)
    - H406: Final cooler (4°C)
    
    Outputs:
    - LegHb_3: Final product stream (meets all specifications)
    """
    
    # Heat to pasteurization temperature
    H403 = bst.HXutility(
        'H403',
        ins=U404-0,
        outs='HeatedConcentrate',
        T=72+273.15,  # 72°C (middle of 70-75°C range)
        heat_only=True,
    )
    
    # Hold time tank (~30 seconds = 0.0083 hours)
    T401 = bst.MixTank(
        'T401', 
        ins=H403-0, 
        outs='PasteurizedConcentrate', 
        tau=30/3600  # 30 second hold time
    )
    
    # Antioxidant stream (SodiumAscorbate) 
    # Use existing SodiumAscorbate chemical from _chemicals.py
    AntioxidantStream = bst.Stream(
        'AntioxidantStream', 
        SodiumAscorbate=0.1,  # Initial value, adjusted by spec
        H2O=1.0,
        units='kg/hr',
        price=price.get('SodiumAscorbate', 5.0)  # Default price if not defined
    )
    
    # Dilution water stream for achieving final total solids target
    DilutionWater = bst.Stream(
        'DilutionWater',
        H2O=100.0,  # Initial value, adjusted by spec
        units='kg/hr',
    )
    
    # Formulation mixer - combines pasteurized concentrate, antioxidant, and dilution water
    M404 = bst.MixTank(
        'M404', 
        ins=(T401-0, AntioxidantStream, DilutionWater), 
        outs='FormulatedProduct',
        tau=0.1
    )
    
    # Target specifications for final product
    M404.target_total_solids_percent = 20.0  # Below 24% limit, target 20%
    M404.target_legh_percent = 7.5  # Target middle of 6-9% range
    
    @M404.add_specification(run=True)
    def update_formulation():
        """
        Adjust antioxidant and dilution water to achieve:
        - 0.1% w/w sodium ascorbate
        - Target total solids (default 20%, max 24%)
        - Target LegHb content (6-9% range)
        """
        import flexsolve as flx
        import warnings
        
        feed = T401.outs[0]
        feed_total_mass = feed.F_mass
        
        if feed_total_mass <= 0:
            DilutionWater.imass['H2O'] = 0
            AntioxidantStream.imass['SodiumAscorbate'] = 0
            AntioxidantStream.imass['H2O'] = 0
            return
        
        # Calculate feed composition
        feed_water_mass = feed.imass['H2O']
        feed_solids_mass = feed_total_mass - feed_water_mass
        feed_legh_mass = feed.imass['Leghemoglobin']
        
        # Target total solids and LegHb percentages
        target_solids_pct = getattr(M404, 'target_total_solids_percent', 20.0)
        target_legh_pct = getattr(M404, 'target_legh_percent', 7.5)
        
        # Calculate required dilution to achieve target total solids
        # total_solids_pct = solids / (solids + water + dilution) * 100
        # Rearranging: dilution = solids * (100/target_pct - 1) - water
        if target_solids_pct > 0:
            required_total_mass = feed_solids_mass * 100.0 / target_solids_pct
            required_dilution = required_total_mass - feed_total_mass
            
            # Also check LegHb target - use whichever requires more dilution
            if feed_legh_mass > 0 and target_legh_pct > 0:
                required_total_mass_legh = feed_legh_mass * 100.0 / target_legh_pct
                required_dilution_legh = required_total_mass_legh - feed_total_mass
                required_dilution = max(required_dilution, required_dilution_legh)
            
            # Ensure non-negative dilution
            required_dilution = max(0, required_dilution)
        else:
            required_dilution = 0
        
        # Set dilution water
        DilutionWater.imass['H2O'] = required_dilution
        
        # Set antioxidant (0.1% w/w of total final product)
        final_product_mass = feed_total_mass + required_dilution
        AntioxidantStream.imass['SodiumAscorbate'] = final_product_mass * 0.001
        AntioxidantStream.imass['H2O'] = AntioxidantStream.imass['SodiumAscorbate'] * 9  # 10% solution
    
    # Final rapid cooling to <4°C
    H406 = bst.HXutility(
        'H406',
        ins=M404-0,
        outs=LegHb_3,
        T=4+273.15,  # 4°C per spec (not 0°C as before)
        cool_only=True,
    )
    
    # =========================================================================
    # AREA 900: FACILITIES & UTILITIES
    # =========================================================================
    """
    Facilities Area - Utility Systems & Wastewater Treatment
    
    Purpose: Provide utilities (cooling, heating, power) and manage wastewater
    
    Subsystems:
    
    A. WASTEWATER TREATMENT & WATER RECOVERY:
       - M501: Collect all liquid waste streams
       - S501: RO treatment for water recovery → PWC recycling
       - M503_debris: Collect solid debris for disposal
    
    B. UTILITY SYSTEMS:
       - CT: Cooling tower (process cooling water)
       - CWP: Chilled water package (refrigeration for 4-10°C cooling)
       - PWC: Process water center (water distribution & recycling)
    
    Units:
    - M501: Wastewater collection mixer
    - S501: Reverse osmosis (water recovery)
    - M503_debris: Solid waste routing mixer
    - CT: Cooling tower (Area 500 ID if use_area_convention)
    - CWP: Chilled water package (Area 500 ID if use_area_convention)
    - PWC: Process water center (Area 500 ID if use_area_convention)
    
    Outputs:
    - effluent1: RO reject → disposal
    - debris_disposal: Solid waste to landfill
    
    Notes:
    - PWC integrates RO-treated water recycling into process water distribution
    """
    
    # Collect permeate streams for RO treatment and recycle
    
    M501 = bst.MixTank(
        'M501', 
        ins=(S404-1, C401-1, C402-1, U401-1, U404-1),  # All liquid waste streams
        outs='CombinedWastewater', 
        tau=1
    )
    M501.add_specification(run=True)
    
    # RO treatment for water recovery
    S501 = u.ReverseOsmosis('S501', ins=M501-0, outs=('RO_treated_water1', effluent1))
    
    # Utility systems
    CT = bst.CoolingTower(500 if use_area_convention else 'CT')
    CWP = bst.ChilledWaterPackage(500 if use_area_convention else 'CWP')
    
    # Route dehydrated debris (S404-0) to separate disposal stream
    # Note: Proteins lack complete gas-phase thermodynamic properties for boiler emissions
    debris_disposal = bst.Stream('debris_disposal')
    M503_debris = bst.Mixer('M503_debris', ins=S404-0, outs=debris_disposal)
    
    makeup_water_streams = (F.cooling_tower_makeup_water,
                            F.Water1, F.Water2,
                            F.Water3, F.WashWater,
                            F.DFBufferWater)
    process_water_streams = (F.S501.outs[1],
                            *makeup_water_streams)
    
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    
    PWC = bst.ProcessWaterCenter(500 if use_area_convention else 'PWC',
        ins=('recycled_RO_water', makeup_water, 'recycled_process_water', 'makeup_process_water'),
        outs=('RO_water', 'process_water', 'excess_water'),
        makeup_water_streams=makeup_water_streams,
        process_water_streams=process_water_streams,
        reverse_osmosis_water_price=0.000254,
        process_water_price=0.000135,
    )
    
    # Update input stream prices (only used streams)
    s.update_all_input_stream_prices(streamlist=[SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, DfUltraBuffer_wash])

    return LegHb_3, vent1, vent2, effluent1, DfUltraBuffer, DfUltraBuffer_wash


# %% Design Specification Functions


def set_production_rate(system, target_production_rate_kg_hr, verbose=True):
    """
    Adjust system inputs to achieve target LegHb_3 production rate using a global scaling factor.
    
    Parameters
    ----------
    system : biosteam.System
        The LegHb production system
    target_production_rate_kg_hr : float
        Target mass flow rate for LegHb_3 stream [kg/hr]
    
    Returns
    -------
    float
        Achieved production rate [kg/hr]
    """
    import flexsolve as flx
    import warnings
    
    log = print if verbose else (lambda *args, **kwargs: None)

    # Get product stream
    product_stream = system.flowsheet.stream.LegHb_3
    
    # Store baseline flow rates for all input streams
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
    
    # Run initial simulation to get baseline
    log("Running initial simulation to establish baseline...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        system.simulate()
    
    initial_production = product_stream.F_mass
    
    if initial_production <= 0:
        raise ValueError("Initial production rate is zero. Cannot scale system.")
    
    # Calculate initial scaling factor guess
    initial_guess = target_production_rate_kg_hr / initial_production
    
    log(f"  Initial production: {initial_production:.2f} kg/hr")
    log(f"  Target production:  {target_production_rate_kg_hr:.2f} kg/hr")
    log(f"  Initial scaling guess: {initial_guess:.4f}x")
    
    # Iteration tracking
    iteration = [0]
    
    # Define objective function that scales inputs, simulates, and returns error
    def objective_function(scaling_factor):
        """
        Scale all input streams, simulate system, return production error.
        
        Parameters
        ----------
        scaling_factor : float
            Global scaling factor for all inputs
            
        Returns
        -------
        float
            Error = (achieved_production - target_production) [kg/hr]
        """
        iteration[0] += 1
        
        try:
            # Scale all input streams
            for stream, baseline_flow in baseline_flows.items():
                stream.F_mass = baseline_flow * scaling_factor
            
            # Simulate system with scaled inputs (warnings suppressed by outer context)
            system.simulate()
            
            # Get achieved production rate
            achieved_rate = product_stream.F_mass
            error = achieved_rate - target_production_rate_kg_hr
            
            # Print progress (only every 5 iterations to reduce clutter)
            if iteration[0] % 5 == 0 or abs(error) < 1.0:
                log(f"    Iteration {iteration[0]}: scale={scaling_factor:.4f}x, "
                    f"production={achieved_rate:.2f} kg/hr, error={error:.4f} kg/hr")
            
            return error
            
        except Exception as e:
            log(f"    Iteration {iteration[0]} FAILED at scale={scaling_factor:.4f}x: {e}")
            # Return large error to guide solver away from infeasible region
            return 1e6 if scaling_factor > initial_guess else -1e6
    
    try:
        log(f"\nSolving for optimal scaling factor (intermediate warnings suppressed)...")
        
        # Suppress all warnings during the iterative solving process
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Build robust brackets around the initial guess
            x0 = max(0.01, initial_guess * 0.1)
            x1 = max(5.0, initial_guess * 2.0)
            y0 = objective_function(x0)
            y1 = objective_function(x1)
            expand_iter = 0
            while y0 * y1 > 0 and expand_iter < 8:
                expand_iter += 1
                # Expand toward the side with larger absolute error
                if abs(y1) <= abs(y0):
                    x0 = max(0.01, x0 / 2.0)
                    y0 = objective_function(x0)
                else:
                    x1 *= 2.0
                    y1 = objective_function(x1)
            log(f"  Search bounds: {x0:.4f}x to {x1:.4f}x (bracket {'OK' if y0*y1 <= 0 else 'NOT FOUND'})")

            # Use flexsolve directly (not as a system specification)
            scaling_factor = flx.IQ_interpolation(
                f=objective_function,
                x0=x0,               # Lower bound
                x1=x1,               # Upper bound
                x=initial_guess,     # Initial guess
                xtol=0.0001,         # Tolerance on scaling factor
                ytol=0.01,           # Tolerance on production rate [kg/hr]
                maxiter=100,
                checkbounds=(y0 * y1 <= 0),
                checkiter=True,
            )
        
        # After solver converges, run final simulation WITHOUT suppressing warnings
        # This allows any persistent issues to be displayed
        log(f"\nSolver converged. Running final validation simulation...")
        for stream, baseline_flow in baseline_flows.items():
            stream.F_mass = baseline_flow * scaling_factor
        
        # Final simulation with warnings enabled
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
        
        # Restore baseline flows
        log(f"\n  Restoring baseline input flows...")
        for stream, baseline_flow in baseline_flows.items():
            stream.F_mass = baseline_flow
        
        # Run one final simulation with baseline flows (warnings enabled)
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
    
    # Get total mass
    total_mass = product_stream.F_mass
    
    if total_mass <= 0:
        raise ValueError("Product stream has zero mass flow")
    
    # Helper function to safely get mass of chemicals
    def get_chemical_mass(chemical_id):
        """Safely get mass of a chemical, return 0 if not present"""
        try:
            return product_stream.imass[chemical_id]
        except:
            return 0
    
    # Calculate mass fractions (as percentages)
    def get_mass_percent(chemical_ids):
        """Get total mass percent for list of chemicals"""
        if isinstance(chemical_ids, str):
            chemical_ids = [chemical_ids]
        total = sum(get_chemical_mass(chem) for chem in chemical_ids)
        return total / total_mass * 100
    
    # Define specifications
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
    
    # Calculate actual values
    for spec_name, spec_data in specs.items():
        spec_data['actual'] = get_mass_percent(spec_data['chemicals'])
    
    # Calculate protein purity
    legh_mass = get_chemical_mass('Leghemoglobin')
    globin_mass = get_chemical_mass('Globin')
    mannoprotein_mass = get_chemical_mass('Mannoprotein')
    total_protein_mass = legh_mass + globin_mass + mannoprotein_mass
    
    protein_purity = (legh_mass / total_protein_mass * 100) if total_protein_mass > 0 else 0
    
    # Print results
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
    
    # Check protein purity
    purity_passed = protein_purity >= 65.0
    status = '✓ PASS' if purity_passed else '✗ FAIL'
    print(f"{'Protein Purity':<30} {'>= 65.0%':<25} {protein_purity:>6.2f}%{'':<7} {status}")
    
    if not purity_passed:
        all_passed = False
    
    print(f"{'='*85}")
    
    # Additional info
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


# %%
if __name__ == '__main__':
    # Set preferences
    bst.preferences.N = 50
    nn=1
    # Define target production rate
    TARGET_PRODUCTION = 275 * nn # kg/hr
    
    print("="*85)
    print("LEGHEMOGLOBIN PRODUCTION SYSTEM - DESIGN SPECIFICATION MODE")
    print("="*85)
    
    # Create the LegHb system
    print("\n1. Creating system...")
    LegHb_sys = create_LegHb_system()
    sys = LegHb_sys
    f = sys.flowsheet
    u = f.unit
    ss = f.stream
    sys.operating_hours = 8000
    
    # Run initial baseline simulation
    print("\n2. Running baseline simulation...")
    try:
        LegHb_sys.simulate()
        baseline_production = ss.LegHb_3.F_mass
        print(f"   Baseline production rate: {baseline_production:.2f} kg/hr")
    except Exception as e:
        print(f"   Baseline simulation failed: {e}")
        raise
    
    # Set production rate to target using design specification
    print(f"\n3. Applying design specification: TARGET_PRODUCTION = {TARGET_PRODUCTION} kg/hr")
    try:
        achieved_production = set_production_rate(LegHb_sys, TARGET_PRODUCTION)
        
        # Verify production rate is maintained
        LegHb_sys.simulate()
        final_production = ss.LegHb_3.F_mass
        
        if abs(final_production - TARGET_PRODUCTION) > 1.0:  # Allow 1 kg/hr tolerance
            print(f"\n   WARNING: Production rate drifted after final simulation!")
            print(f"   Target:  {TARGET_PRODUCTION:.2f} kg/hr")
            print(f"   Actual:  {final_production:.2f} kg/hr")
            
    except Exception as e:
        print(f"\n   Could not achieve target production: {e}")
        print("   Continuing with baseline production rate...")
        achieved_production = ss.LegHb_3.F_mass
    
    # Check product specifications
    print(f"\n4. Verifying product specifications...")
    try:
        check_LegHb_specifications(ss.LegHb_3)
    except ValueError as e:
        print(f"\n   SPECIFICATION CHECK FAILED: {e}")
        print("   System may require process parameter adjustments to meet specifications.")
    
    # Display system results
    print(f"\n5. System Summary")
    print("="*85)
    LegHb_sys.show()
    
    # Calculate key metrics
    legh_purity = ss.LegHb_3.imass['Leghemoglobin'] / ss.LegHb_3.F_mass * 100
    print(f"\n{'='*85}")
    print("KEY PERFORMANCE INDICATORS")
    print(f"{'='*85}")
    print(f"  Product Stream:           {ss.LegHb_3.ID}")
    print(f"  Production Rate:          {ss.LegHb_3.F_mass:.2f} kg/hr")
    print(f"  Leghemoglobin Content:    {legh_purity:.2f}%")
    print(f"  Annual Production:        {ss.LegHb_3.F_mass * sys.operating_hours / 1000:.2f} metric tons/year")
    print(f"{'='*85}\n")
    
    # Generate system diagram
    print(f"\n6. Generating system diagram...")
    LegHb_sys.diagram(format='html', display=True)
    
    # LCA analysis
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
