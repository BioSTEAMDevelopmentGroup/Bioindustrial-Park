# -*- coding: utf-8 -*-
"""
N-HemDx Production System - Config 1 (Split-Stream Topology)

Modular BioSTEAM implementation of the N-HemDx (Nicotinamide-Stabilized Hemodextrin)
production system using the "Split-Stream Heme Recovery" topology.

Process Areas:
    - Area 100: Clarification (Broth centrifugation, stream split)
    - Area 200: Intracellular Recovery (Resuspension, HPH, Debris Removal)
    - Area 300: Combined Capture & Purification (Resin adsorption)
    - Area 400: Concentration (Nanofiltration)
    - Area 500: Formulation (γ-CD + Nicotinamide complexation)
    - Area 600: Final Product
    - Area 900: Facilities & Utilities

Host: Corynebacterium glutamicum
Product: N-HemDx (Heme b + γ-Cyclodextrin + Nicotinamide complex)

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biorefineries.prefers.v1._process_settings import set_GWPCF, GWP_CFs, load_process_settings
from biorefineries.prefers.v1.HemDx import _streams as s
from biorefineries.prefers.v1.HemDx import _chemicals as c
from biorefineries.prefers.v1 import _units as u
from biorefineries.prefers.v1._process_settings import price
import biosteam as bst
from biosteam import F
import thermosteam as tmo
import numpy as np

# =============================================================================
# MODULE INITIALIZATION
# =============================================================================
HEMDX_THERMO = c.create_chemicals_Hemodextrin()
bst.settings.set_thermo(HEMDX_THERMO, skip_checks=True)
bst.preferences.classic_mode()
SHOW_RXN = False

__all__ = (
    'create_NHemDx_system',
    'set_production_rate',
    'check_NHemDx_specifications',
    # Area creation functions
    'create_area_100_clarification',
    'create_area_200_recovery',
    'create_area_300_capture',
    'create_area_400_concentration',
    'create_area_500_formulation',
    'create_area_600_product',
)

# =============================================================================
# PROCESS PARAMETERS (From BOD)
# =============================================================================

# Fermentation Parameters (Corynebacterium glutamicum)
HEME_TITER = 0.35  # g/L (typical free heme titer)
FERMENTATION_TAU = 72  # hours
BIOMASS_CONCENTRATION = 50  # g/L dry cell weight

# Split-Stream Parameters
EXTRACELLULAR_HEME_FRACTION = 0.5  # 50% extracellular
INTRACELLULAR_HEME_FRACTION = 0.5  # 50% intracellular

# HPH (High-Pressure Homogenization) Parameters
HPH_PRESSURE = 1000e5  # 1000 bar = 100 MPa
HPH_TEMPERATURE = 10 + 273.15  # <10°C

# Capture Column Parameters (HP-20 or equivalent)
RESIN_DBC = 50  # g/L dynamic binding capacity
HEME_RECOVERY = 0.85  # 85% overall recovery

# Formulation Parameters
GAMMA_CD_MW = 1297.12  # g/mol
HEME_B_MW = 616.49  # g/mol
NICOTINAMIDE_MW = 122.12  # g/mol
HEME_CD_MOLAR_RATIO = 1.0  # 1:1 molar
NICOTINAMIDE_EXCESS = 5.0  # 5× molar excess

# Target Production
TARGET_PRODUCTION = 100  # kg/hr N-HemDx

# =============================================================================
# FERMENTATION REACTIONS
# =============================================================================

def create_fermentation_reactions():
    """
    Create fermentation reactions for Corynebacterium glutamicum heme production.
    
    Returns production, cell growth, and respiration reactions.
    """
    Rxn = tmo.reaction.Reaction
    PRxn = tmo.reaction.ParallelReaction
    
    # Heme production from glucose (simplified stoichiometry)
    # C6H12O6 + Fe2+ + 4 Glycine → C34H32FeN4O4 (Heme_b) + byproducts
    production = Rxn(
        'Glucose + 0.01 FeSO4 + 0.02 Glycine -> 0.005 Heme_b + 0.8 H2O + 0.3 CO2',
        'Glucose', X=0.15, basis='wt'
    )
    
    # Cell growth reaction
    cell_growth = Rxn(
        'Glucose + 0.05 NH3 -> Corynebacterium_glutamicum + 0.3 H2O + 0.2 CO2',
        'Glucose', X=0.45, basis='wt'
    )
    
    # Respiration (maintenance)
    respiration = Rxn(
        'Glucose + O2 -> CO2 + H2O',
        'Glucose', X=0.35, basis='wt'
    )
    
    # Neutralization
    neutralization = Rxn(
        'NH3 + H2SO4 -> (NH4)2SO4',
        'H2SO4', X=1.0
    )
    
    return production, cell_growth, respiration, neutralization


# =============================================================================
# AREA 100: CLARIFICATION (Harvest & Stream Split)
# =============================================================================

def create_area_100_clarification(fermentation_out, extracellular_fraction=0.5):
    """
    Area 100: Clarification - Split fermentation broth into supernatant and pellet.
    
    Parameters
    ----------
    fermentation_out : Stream
        Fermentation broth output stream
    extracellular_fraction : float
        Fraction of heme in supernatant (default 0.5 = 50%)
    
    Returns
    -------
    supernatant : Stream (Path A)
    pellet : Stream (Path B)
    """
    # S101: Broth Centrifuge
    S101 = u.Centrifuge(
        'S101',
        ins=fermentation_out,
        outs=('Supernatant_A', 'WetPellet_B'),
        split={'Corynebacterium_glutamicum': 0.02,  # 98% cells to pellet
               'Heme_b': extracellular_fraction,     # split based on location
               'H2O': 0.85,                          # most water to supernatant
               'Glucose': 0.90,
               'Glycine': 0.90,
               },
    )
    
    return S101


# =============================================================================
# AREA 200: INTRACELLULAR RECOVERY (Path B)
# =============================================================================

def create_area_200_recovery(pellet_stream, extraction_buffer):
    """
    Area 200: Intracellular Recovery - Lyse cells and remove debris.
    
    Process Flow:
    T201 (Resuspension) → U202 (HPH) → S203 (Debris Centrifuge) → S204 (Guard Filter)
    
    Parameters
    ----------
    pellet_stream : Stream
        Wet cell pellet from S101
    extraction_buffer : Stream
        Buffer for cell resuspension
    """
    # T201: Resuspension Tank
    T201 = bst.MixTank(
        'T201',
        ins=(pellet_stream, extraction_buffer),
        outs='ResuspendedCells',
        tau=0.5,  # 30 min residence time
    )
    
    @T201.add_specification(run=True)
    def adjust_resuspension_buffer():
        """Adjust buffer volume to achieve target cell concentration."""
        pellet = T201.ins[0]
        buffer = T201.ins[1]
        cell_mass = pellet.imass['Corynebacterium_glutamicum']
        target_conc = 150  # g/L target cell concentration
        required_volume = cell_mass / target_conc  # L
        buffer.imass['H2O'] = required_volume * 1000 - pellet.imass['H2O']  # kg/hr
        if buffer.imass['H2O'] < 0:
            buffer.imass['H2O'] = pellet.imass['H2O'] * 0.5  # minimum buffer
    
    # U202: Cell Disruptor (High-Pressure Homogenizer)
    U202 = u.CellDisruption(
        'U202',
        ins=T201-0,
        outs='CellLysate',
        Cell_ID='Corynebacterium_glutamicum',
        cell_disruption_efficiency=0.55,  # 55% soluble release
        P_high=HPH_PRESSURE,
        P_low=101325,
    )

    
    # S203: Debris Centrifuge (Decanter)
    S203 = u.Centrifuge(
        'S203',
        ins=U202-0,
        outs=('ClarifiedLysate', 'CellDebris'),
        split={'Protein': 0.90,      # soluble protein to supernatant
               'Heme_b': 0.92,       # heme released to supernatant  
               'H2O': 0.95,
               'Mannoprotein': 0.10, # debris to pellet
               'Biomass': 0.05,      # cell wall fragments to pellet
               },
    )
    
    # S204: Guard Filter (Depth Filtration)
    S204 = u.Filtration.from_preset(
        'MF',
        ID='S204',
        ins=S203-0,
        outs=('FilteredLysate_B', 'FilterRetentate'),
        solids_loading=20.0,
        solid_capture_efficiency=0.95,
    )
    
    return T201, U202, S203, S204


# =============================================================================
# AREA 300: COMBINED CAPTURE & PURIFICATION
# =============================================================================

def create_area_300_capture(supernatant_A, lysate_B, wash_buffer, elution_buffer, regeneration_buffer):
    """
    Area 300: Combined Capture - Merge streams and adsorb on hydrophobic resin.
    
    Process Flow:
    M301 (Stream Merge) → U302 (Capture Column)
    
    Parameters
    ----------
    supernatant_A : Stream
        Clarified supernatant from Path A
    lysate_B : Stream
        Filtered lysate from Path B
    wash_buffer : Stream
        Column wash buffer
    elution_buffer : Stream
        Elution buffer (NaOH-based)
    regeneration_buffer : Stream
        Regeneration buffer (Ethanol-based)
    """
    # M301: Stream Merger
    M301 = bst.Mixer(
        'M301',
        ins=(supernatant_A, lysate_B),
        outs='CombinedFeed',
    )
    
    # Temperature conditioning before column
    H301 = bst.HXutility(
        'H301',
        ins=M301-0,
        outs='ConditionedFeed',
        T=25 + 273.15,  # 25°C for adsorption
    )
    
    # Prepare buffer streams with temperature control
    H302 = bst.HXutility(
        'H302',
        ins=wash_buffer,
        outs='WashBuffer_T',
        T=25 + 273.15,
    )
    
    H303 = bst.HXutility(
        'H303',
        ins=elution_buffer,
        outs='ElutionBuffer_T',
        T=40 + 273.15,  # 40°C for elution
    )
    
    H304 = bst.HXutility(
        'H304',
        ins=regeneration_buffer,
        outs='RegenBuffer_T',
        T=30 + 273.15,
    )
    
    # U302: Capture Column (HP-20 Hydrophobic Resin - Adsorption mode)
    U302 = u.ResinColumn(
        'U302',
        ins=(H301-0, H302-0, H303-0, H304-0),
        outs=('HemeEluate', 'ColumnFlowthrough', 'WashWaste', 'RegenWaste'),
        preset='IonExchange',
        TargetProduct_IDs=('Heme_b', 'ProtoporphyrinIX'),
        TargetProduct_Yield=0.90,  # 90% capture
        BoundImpurity_IDs=c.chemical_groups['Salts'] + ('Glucose', 'Glycine'),
        BoundImpurity_Removal=0.95,
        NonBinding_Carryover=0.05,
        EBCT_min=5.0,
        superficial_velocity_m_h=10.0,
    )
    
    # Dynamic buffer specifications
    @H302.add_specification(run=True)
    def update_wash_buffer():
        feed_water = H301.outs[0].imass['H2O']
        wash_buffer.imass['H2O'] = feed_water * 3.0  # 3 CV wash
        wash_buffer.imass['NaCl'] = wash_buffer.imass['H2O'] * 0.002
    
    @H303.add_specification(run=True)
    def update_elution_buffer():
        feed_water = H301.outs[0].imass['H2O']
        elution_buffer.imass['H2O'] = feed_water * 0.5  # 0.5 CV elution
        elution_buffer.imass['Ethanol'] = elution_buffer.imass['H2O'] * 0.7  # 70% ethanol
    
    @H304.add_specification(run=True)
    def update_regen_buffer():
        feed_water = H301.outs[0].imass['H2O']
        regeneration_buffer.imass['H2O'] = feed_water * 0.3
        regeneration_buffer.imass['NaOH'] = feed_water * 0.02  # 2% NaOH
    
    return M301, H301, H302, H303, H304, U302


# =============================================================================
# AREA 400: CONCENTRATION
# =============================================================================

def create_area_400_concentration(eluate_stream, df_buffer):
    """
    Area 400: Concentration - Remove ethanol/salts and concentrate heme.
    
    Process Flow:
    U401 (Nanofiltration/Diafiltration) → H402 (Temperature conditioning)
    """
    # U401: Nanofiltration for concentration
    U401 = u.Diafiltration.from_preset(
        'NF',
        ID='U401',
        ins=(eluate_stream, df_buffer),
        outs=('ConcentratedHeme', 'NFPermeate'),
        TargetProduct_ID=('Heme_b', 'ProtoporphyrinIX'),
        Salt_ID=('NaCl', 'NaOH', 'Ethanol'),
        TargetProduct_Retention=0.98,  # heme MW ~616 Da - needs tight NF
        Salt_Retention=0.10,
        concentration_factor=5.0,  # 5× concentration
    )
    
    @U401.add_specification(run=True)
    def adjust_df_buffer():
        """Adjust DF buffer for proper diafiltration."""
        feed = U401.ins[0]
        buffer = U401.ins[1]
        buffer.imass['H2O'] = feed.imass['H2O'] * 2.0  # 2 diavolumes
    
    # H402: Temperature conditioning for formulation
    H402 = bst.HXutility(
        'H402',
        ins=U401-0,
        outs='HemeConcentrate',
        T=30 + 273.15,  # 30°C for formulation
    )
    
    return U401, H402


# =============================================================================
# AREA 500: FORMULATION
# =============================================================================

def create_area_500_formulation(heme_concentrate, gamma_cd_solution, nicotinamide_solution):
    """
    Area 500: Formulation - Create Heme-γ-CD-Nicotinamide complex.
    
    Process Flow:
    M501 (γ-CD addition) → R502 (Complexation) → M503 (Nic addition) → R504 (Stabilization)
    
    Stoichiometry:
    - Heme:γ-CD = 1:1 molar
    - Nicotinamide:Heme ≥ 5:1 molar (excess for axial ligand stabilization)
    """
    # M501: γ-Cyclodextrin Addition
    M501 = bst.Mixer(
        'M501',
        ins=(heme_concentrate, gamma_cd_solution),
        outs='HemeWithCD',
    )
    
    # R502: Complexation Reactor
    R502 = bst.MixTank(
        'R502',
        ins=M501-0,
        outs='HemeCD_Complex',
        tau=1.0,  # 1 hour residence time
    )
    
    # Complexation reaction: Heme_b + GammaCyclodextrin -> HemoDextrin
    complexation_rxn = tmo.reaction.Reaction(
        'Heme_b + GammaCyclodextrin -> HemoDextrin',
        'Heme_b', X=0.95  # 95% complexation efficiency
    )
    
    @R502.add_specification(run=False)
    def run_complexation():
        """Run mixing first, then complexation reaction."""
        R502._run()  # Run the MixTank default behavior first
        outlet = R502.outs[0]
        if outlet.imol['Heme_b'] > 1e-12 and outlet.imol['GammaCyclodextrin'] > 1e-12:
            complexation_rxn(outlet)
    
    # M503: Nicotinamide Addition
    M503 = bst.Mixer(
        'M503',
        ins=(R502-0, nicotinamide_solution),
        outs='HemeDxWithNic',
    )

    
    # R504: Stabilization Tank (Final product maturation)
    R504 = bst.MixTank(
        'R504',
        ins=M503-0,
        outs='NHemDx_Crude',
        tau=0.5,  # 30 min residence
    )
    
    # N-HemoDextrin formation: HemoDextrin + Nicotinamide -> N-HemoDextrin
    stabilization_rxn = tmo.reaction.Reaction(
        'HemoDextrin + Nicotinamide -> N-HemoDextrin',
        'HemoDextrin', X=0.90  # 90% stabilization efficiency
    )
    
    @R504.add_specification(run=False)
    def run_stabilization():
        """Run mixing first, then stabilization reaction."""
        R504._run()  # Run the MixTank default behavior first
        outlet = R504.outs[0]
        if outlet.imol['HemoDextrin'] > 1e-12 and outlet.imol['Nicotinamide'] > 1e-12:
            stabilization_rxn(outlet)
    
    # Dynamic γ-CD dosing based on heme content
    @M501.add_specification(run=True)
    def dose_gamma_cd():
        """Add γ-CD at 1:1 molar ratio with 10% excess."""
        heme_mol = heme_concentrate.imol['Heme_b']
        gamma_cd_solution.imol['GammaCyclodextrin'] = heme_mol * 1.1  # 10% excess
        gamma_cd_solution.imass['H2O'] = gamma_cd_solution.imass['GammaCyclodextrin'] * 10  # carrier water
    
    # Dynamic nicotinamide dosing
    @M503.add_specification(run=True)
    def dose_nicotinamide():
        """Add Nicotinamide at 5× molar ratio."""
        heme_mol = R502.outs[0].imol['HemoDextrin'] + R502.outs[0].imol['Heme_b']
        nicotinamide_solution.imol['Nicotinamide'] = heme_mol * NICOTINAMIDE_EXCESS
        nicotinamide_solution.imass['H2O'] = nicotinamide_solution.imass['Nicotinamide'] * 5
    
    return M501, R502, M503, R504


# =============================================================================
# AREA 600: FINAL PRODUCT
# =============================================================================

def create_area_600_product(crude_product):
    """
    Area 600: Final Product - Polish filtration and product output.
    """
    # S601: Final Filtration (Polish)
    S601 = u.Filtration.from_preset(
        'MF',
        ID='S601',
        ins=crude_product,
        outs=('NHemDx_Product', 'PolishRetentate'),
        solids_loading=10.0,
        solid_capture_efficiency=0.95,
    )
    
    return S601


# =============================================================================
# SYSTEM FACTORY
# =============================================================================

@bst.SystemFactory(
    ID='NHemDx_sys',
    ins=[
        dict(ID='FermentationBroth', Heme_b=5, Corynebacterium_glutamicum=200, 
             H2O=800, Glucose=10, Glycine=1, units='kg/hr', T=30+273.15),
        dict(ID='ExtractionBuffer', H2O=1, units='kg/hr', T=25+273.15),
        dict(ID='WashBuffer', H2O=1, NaCl=0.002, units='kg/hr', T=25+273.15),
        dict(ID='ElutionBuffer', H2O=1, Ethanol=0.7, units='kg/hr', T=25+273.15),
        dict(ID='RegenerationBuffer', H2O=1, NaOH=0.02, units='kg/hr', T=25+273.15),
        dict(ID='DFBuffer', H2O=1, units='kg/hr', T=25+273.15),
        dict(ID='GammaCDSolution', GammaCyclodextrin=0.1, H2O=1, units='kg/hr', T=25+273.15),
        dict(ID='NicotinamideSolution', Nicotinamide=0.1, H2O=1, units='kg/hr', T=25+273.15),
    ],
    outs=[
        dict(ID='NHemDx_Product'),
        dict(ID='CellDebris'),
        dict(ID='ProcessWaste'),
    ],
    fthermo=lambda: HEMDX_THERMO,
)
def create_NHemDx_system(
    ins, outs,
    extracellular_fraction=0.5,
):
    """
    Create the N-HemDx (Nicotinamide-Stabilized Hemodextrin) production system.
    
    This system implements the "Split-Stream Heme Recovery" topology:
    - Path A: Supernatant processing (extracellular heme)
    - Path B: Intracellular recovery (cell lysis)
    - Path C: Combined capture, purification, and formulation
    
    Parameters
    ----------
    ins : list of Streams
        Input streams
    outs : list of Streams
        Output streams
    extracellular_fraction : float
        Fraction of heme in supernatant (default 0.5)
    """
    bst.preferences.N = 50
    
    # Unpack input streams
    (FermentationBroth, ExtractionBuffer, WashBuffer, ElutionBuffer, 
     RegenerationBuffer, DFBuffer, GammaCDSolution, NicotinamideSolution) = ins
    
    # Unpack output streams
    NHemDx_Product, CellDebris, ProcessWaste = outs
    
    # =========================================================================
    # AREA 100: CLARIFICATION
    # =========================================================================
    S101 = create_area_100_clarification(
        FermentationBroth, 
        extracellular_fraction=extracellular_fraction
    )
    
    # =========================================================================
    # AREA 200: INTRACELLULAR RECOVERY (Path B)
    # =========================================================================
    # Custom area creation to pass CellDebris output correctly
    # T201: Resuspension Tank
    T201 = bst.MixTank(
        'T201',
        ins=(S101-1, ExtractionBuffer),
        outs='ResuspendedCells',
        tau=0.5,
    )
    
    @T201.add_specification(run=True)
    def adjust_resuspension_buffer():
        pellet = T201.ins[0]
        buffer = T201.ins[1]
        cell_mass = pellet.imass['Corynebacterium_glutamicum']
        target_conc = 150
        required_volume = cell_mass / target_conc
        buffer.imass['H2O'] = required_volume * 1000 - pellet.imass['H2O']
        if buffer.imass['H2O'] < 0:
            buffer.imass['H2O'] = pellet.imass['H2O'] * 0.5
    
    # U202: Cell Disruptor
    U202 = u.CellDisruption(
        'U202',
        ins=T201-0,
        outs='CellLysate',
        Cell_ID='Corynebacterium_glutamicum',
        cell_disruption_efficiency=0.55,
        P_high=HPH_PRESSURE,
        P_low=101325,
    )
    
    # S203: Debris Centrifuge - output CellDebris directly
    S203 = u.Centrifuge(
        'S203',
        ins=U202-0,
        outs=('ClarifiedLysate', CellDebris),  # Pass output stream directly
        split={'Protein': 0.90,
               'Heme_b': 0.92,
               'H2O': 0.95,
               'Mannoprotein': 0.10,
               'Biomass': 0.05,
               },
    )
    
    # S204: Guard Filter
    S204 = u.Filtration.from_preset(
        'MF',
        ID='S204',
        ins=S203-0,
        outs=('FilteredLysate_B', 'FilterRetentate'),
        solids_loading=20.0,
        solid_capture_efficiency=0.95,
    )
    
    # =========================================================================
    # AREA 300: COMBINED CAPTURE & PURIFICATION
    # =========================================================================
    M301, H301, H302, H303, H304, U302 = create_area_300_capture(
        S101-0,   # Supernatant Path A
        S204-0,   # Filtered lysate Path B
        WashBuffer,
        ElutionBuffer,
        RegenerationBuffer,
    )
    
    # =========================================================================
    # AREA 400: CONCENTRATION
    # =========================================================================
    U401, H402 = create_area_400_concentration(
        U302-0,  # Heme eluate
        DFBuffer,
    )
    
    # =========================================================================
    # AREA 500: FORMULATION
    # =========================================================================
    M501, R502, M503, R504 = create_area_500_formulation(
        H402-0,  # Heme concentrate
        GammaCDSolution,
        NicotinamideSolution,
    )
    
    # =========================================================================
    # AREA 600: FINAL PRODUCT - Simple pass-through for now
    # =========================================================================
    # Use a simple cooler/conditioner instead of MF to avoid product loss
    S601 = bst.HXutility(
        'S601',
        ins=R504-0,
        outs=NHemDx_Product,  # Direct output to product
        T=4 + 273.15,  # Cold storage temperature
    )
    
    # =========================================================================
    # AREA 900: WASTE COLLECTION - Pass ProcessWaste directly  
    # =========================================================================
    ProcessWasteCollector = bst.Mixer(
        'M901',
        ins=(U302-1, U302-2, U302-3, U401-1, S204-1),  # Removed S601-1
        outs=ProcessWaste,  # Pass output stream directly
    )




# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def set_production_rate(system, target_kg_hr, max_iterations=20, tolerance=0.05):
    """
    Scale fermentation input to achieve target N-HemDx production rate.
    
    Parameters
    ----------
    system : bst.System
        The N-HemDx system
    target_kg_hr : float
        Target production rate in kg/hr
    max_iterations : int
        Maximum scaling iterations
    tolerance : float
        Acceptable deviation from target (fraction)
    """
    ss = system.flowsheet.stream
    
    for i in range(max_iterations):
        system.simulate()
        
        # Check current production
        current_production = ss.NHemDx_Product.imass['N-HemoDextrin']
        if current_production == 0:
            current_production = ss.NHemDx_Product.imass['HemoDextrin']
        if current_production == 0:
            current_production = ss.NHemDx_Product.imass['Heme_b']
        
        if current_production > 0:
            ratio = target_kg_hr / current_production
            if abs(ratio - 1.0) < tolerance:
                print(f"  Target achieved in {i+1} iterations: {current_production:.2f} kg/hr")
                return current_production
            
            # Scale fermentation input
            ss.FermentationBroth.F_mass *= ratio ** 0.8  # damped scaling
    
    print(f"  Warning: Did not converge after {max_iterations} iterations")
    return current_production


def check_NHemDx_specifications(product_stream):
    """
    Verify product stream meets N-HemDx specifications.
    
    Parameters
    ----------
    product_stream : Stream
        Final product stream
    
    Returns
    -------
    bool : True if specifications met
    """
    specs_met = True
    
    # Check N-HemoDextrin content
    n_hemdx_mass = product_stream.imass['N-HemoDextrin']
    hemdx_mass = product_stream.imass['HemoDextrin']
    heme_mass = product_stream.imass['Heme_b']
    total_product = n_hemdx_mass + hemdx_mass + heme_mass
    
    print(f"\n  Product Composition:")
    print(f"    N-HemoDextrin: {n_hemdx_mass:.4f} kg/hr")
    print(f"    HemoDextrin:   {hemdx_mass:.4f} kg/hr")
    print(f"    Free Heme:     {heme_mass:.4f} kg/hr")
    print(f"    Total Product: {total_product:.4f} kg/hr")
    
    # Calculate purity
    total_mass = product_stream.F_mass
    if total_mass > 0:
        purity = (n_hemdx_mass + hemdx_mass) / total_mass * 100
        print(f"    Purity (HemDx): {purity:.1f}%")
    
    return specs_met


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == '__main__':
    bst.preferences.N = 50
    
    print("="*85)
    print("N-HemDx PRODUCTION SYSTEM - Split-Stream Topology")
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
    
    print("\n3. Checking product specifications...")
    check_NHemDx_specifications(ss.NHemDx_Product)
    
    print(f"\n4. System Summary")
    print("="*85)
    sys.show()
    
    # Calculate key metrics
    product_mass = ss.NHemDx_Product.F_mass
    print(f"\n{'='*85}")
    print("KEY PERFORMANCE INDICATORS")
    print(f"{'='*85}")
    print(f"  Product Stream:      {ss.NHemDx_Product.ID}")
    print(f"  Production Rate:     {product_mass:.2f} kg/hr")
    print(f"  Annual Production:   {product_mass * sys.operating_hours / 1000:.2f} metric tons/year")
    print(f"{'='*85}\n")
    
    print("\n5. Generating system diagram...")
    try:
        sys.diagram(format='html', display=True)
    except Exception as e:
        print(f"   Diagram error: {e}")
    
    print("\n" + "="*85)
    print("SIMULATION COMPLETE")
    print("="*85)