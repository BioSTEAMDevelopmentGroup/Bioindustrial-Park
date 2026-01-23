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

# %%
__all__ = (
    'create_LegHb_system',
    'set_production_rate',
    'check_LegHb_specifications',
    # Modular area builders
    'get_fermentation_parameters',
    'create_fermentation_reactions',
    'create_area_200_media_prep',
    'create_area_300_conversion',
    'create_area_400_recovery',
    'create_area_500_purification',
    'create_area_600_formulation',
    'create_area_900_facilities',
)
# %%
def get_fermentation_parameters():
    """
    Return fermentation parameters for Config 2.
    """
    return {
        'theta_O2': 0.5,
        'agitation_power': 0.985,
        'design': 'Stirred tank',
        'method': 'Riet',
        'T_operation': 273.15 + 32,
        'Q_O2_consumption': -110 * 4184,
        'dT_hx_loop': 8,
        'cooler_pressure_drop': 20684,
        'compressor_isentropic_efficiency': 0.85,
        'V_max': 500,
        'titer_LegHb': 5,
        'productivity_LegHb': 5 / 72,
        'Y_p': 5 * 4 / 600,
        'Y_b': 0.53,
    }


def create_fermentation_reactions(params=None):
    """
    Create reaction systems for LegHb fermentation (Config 2).
    """
    if params is None:
        params = get_fermentation_parameters()

    yield_LegHb = params['Y_p']
    Y_b = params['Y_b']

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
        'H2SO4 + 2 NH3 -> (NH4)2SO4', reactant='NH3', X=1,
        check_atomic_balance=True
    )

    cell_growth_reaction = bst.Rxn(
        'Glucose + 0.8364 NH3 + 0.0108 (NH4)2SO4 -> 2.01 H2O + 0.106 O2 + 6 Pichia_pastoris', 'Glucose', X=Y_b,
        correct_atomic_balance=True
    )
    cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=1 - Y_b)

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

    return {
        'fermentation_reaction': fermentation_reaction,
        'neutralization_reaction': neutralization_reaction,
        'cell_growth_reaction': cell_growth_reaction,
        'respiration_reaction1': respiration_reaction1,
        'respiration_reaction2': respiration_reaction2,
        'RXN': RXN,
    }


def create_area_200_media_prep(SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt):
    """
    Area 200: Media Preparation (Config 2).
    """
    M301 = bst.MixTank('M301', ins=[SeedIn1, 'Water1'], outs='M301Out', tau=16)

    @M301.add_specification(run=True)
    def update_seed1_inputs():
        target_stream = bst.Stream(**{**s.SeedSolution1, 'ID': None})
        SeedIn1.imass['Seed'] = target_stream.imass['Seed']
        M301.ins[1].imass['H2O'] = target_stream.imass['H2O']
        M301.ins[1].T = 25 + 273.15

    M302 = bst.MixTank('M302', ins=[SeedIn2, CultureIn, 'Water2'], outs='M302Out', tau=16)

    @M302.add_specification(run=True)
    def update_culture_inputs():
        target_stream = bst.Stream(**{**s.SeedSolution2, 'ID': None})
        SeedIn2.imass['Seed'] = target_stream.imass['Seed']
        M302.ins[2].imass['H2O'] = target_stream.imass['H2O']
        M302.ins[2].T = 25 + 273.15
        CultureIn.imass['Culture'] = target_stream.imass['SeedSolution'] * (0.1 + 60 + 0.15191) / 1000

    M303 = u.SeedHoldTank('M303', ins=[M301-0, M302-0], outs='M303Out')

    M304 = bst.MixTank('M304', ins=[Glucose, 'Water3'], outs='M304Out', tau=16)

    @M304.add_specification(run=True)
    def update_water_content():
        M304.ins[1].imass['H2O'] = Glucose.imass['Glucose'] / 2
        M304.ins[1].T = 25 + 273.15

    T301 = bst.StorageTank('T301', ins=M304-0, outs='T301Out', tau=16*4+72)

    T302 = u.AmmoniaStorageTank('T302', ins=NH3_25wt, outs='T302Out')

    return {
        'units': [M301, M302, M303, M304, T301, T302],
        'seed_out': M303.outs[0],
        'glucose_out': T301.outs[0],
        'ammonia_out': T302.outs[0],
        'M301': M301, 'M302': M302, 'M303': M303,
        'M304': M304, 'T301': T301, 'T302': T302,
    }


def create_area_300_conversion(seed_in, glucose_in, ammonia_in, vent1, vent2, reactions, params=None):
    """
    Area 300: Fermentation / Conversion (Config 2).
    """
    if params is None:
        params = get_fermentation_parameters()

    fermentation_reaction = reactions['fermentation_reaction']
    cell_growth_reaction = reactions['cell_growth_reaction']
    respiration_reaction1 = reactions['respiration_reaction1']
    respiration_reaction2 = reactions['respiration_reaction2']
    neutralization_reaction = reactions['neutralization_reaction']
    RXN = reactions['RXN']

    R301 = u.SeedTrain(
        'R301',
        ins=[seed_in],
        outs=[vent1, 'R301Out'],
        reactions=bst.PRxn([cell_growth_reaction, respiration_reaction1]),
        saccharification=None,
        T=32 + 273.15,
    )
    R301.add_specification(run=True)

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
    R302.target_yield = params['Y_p']

    @R302.add_specification(run=True)
    def update_reaction_time_and_yield():
        R302.tau = R302.target_titer / R302.target_productivity
        fermentation_reaction[2].product_yield('Leghemoglobin_In', basis='wt', product_yield=R302.target_yield)

    return {
        'units': [R301, R302],
        'broth_out': R302.outs[1],
        'R301': R301, 'R302': R302,
    }


def create_area_400_recovery(broth_in, DfUltraBuffer_wash):
    """
    Area 400: Recovery (Config 2 - aligned with Config 1 recovery train).
    """
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

    M401 = bst.MixTank('M401', ins=(DfUltraBuffer_wash, 'WashWater'), outs='WashBufferOut', tau=0.5)

    @M401.add_specification(run=True)
    def update_wash_buffer():
        M401.ins[1].imass['H2O'] = C401.outs[0].F_mass * 1.5
        M401.ins[0].imol['DfUltraBuffer'] = M401.ins[1].imass['H2O'] * (0.025 + 0.01 + 0.001) / 1000

    H402 = bst.HXutility(
        'H402',
        ins=M401-0,
        outs='ColdWashBuffer',
        T=10 + 273.15,
        cool_only=True,
    )

    M402 = bst.MixTank('M402', ins=(C401-0, H402-0), outs='WashedCellSlurry', tau=0.25)

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

    S401 = u.CellDisruption(
        'S401',
        ins=C402-0,
        outs='CrudeHomogenate',
        P_high=1000e5,
        cell_disruption_efficiency=0.90,
    )

    H401 = bst.HXutility(
        'H401',
        ins=S401-0,
        outs='CooledHomogenate',
        T=15 + 273.15,
        cool_only=True,
    )

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

    S403 = u.Filtration.from_preset(
        'MF',
        'S403',
        ins=S402-1,
        outs=('FilterCake', 'ClarifiedLysate'),
        solid_capture_efficiency=0.95,
        cake_moisture_content=0.30,
    )
    S403.add_specification(run=True)

    S404 = bst.ScrewPress(
        'S404',
        ins=(S402-0, S403-0),
        outs=('DehydratedDebris', 'PressLiquor'),
        split=0.999,
        moisture_content=0.001,
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


def create_area_500_purification(clarified_lysate, DfUltraBuffer, IXEquilibriumBuffer, IXElutionBuffer, IXRegenerationSolution, DfNanoBuffer):
    """
    Area 500: UF/DF -> IEX -> DF (Config 2).
    """
    M401 = bst.MixTank('M401', ins=(DfUltraBuffer, 'Water4'), outs='M401Out', tau=1)

    @M401.add_specification(run=True)
    def update_DfUltraBuffer_initial():
        M401.ins[1].imass['H2O'] = clarified_lysate.imass['H2O'] * 4
        M401.ins[0].imol['DfUltraBuffer'] = (clarified_lysate.imass['H2O'] * 4) * (0.025 + 0.01 + 0.001) / 1000

    H401 = bst.HXutility('H401', ins=M401-0, outs='H401Out', T=5 + 273.15, cool_only=True)

    U401 = u.Diafiltration.from_preset(
        'UF',
        'U401',
        ins=(clarified_lysate, H401-0),
        outs=('U401Out', 'PermeateWasteUltra'),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
    )
    U401.add_specification(run=True)

    M402 = bst.MixTank('M402', ins=(IXEquilibriumBuffer, 'Water5'), outs='M402Out', tau=1)
    M403 = bst.MixTank('M403', ins=(IXElutionBuffer, 'Water6'), outs='M403Out', tau=1)
    M404 = bst.MixTank('M404', ins=(IXRegenerationSolution, 'Water7'), outs='M404Out', tau=1)

    H402 = bst.HXutility('H402', ins=M402-0, outs='H402Out', T=5 + 273.15, cool_only=True)
    H403 = bst.HXutility('H403', ins=M403-0, outs='H403Out', T=5 + 273.15, cool_only=True)
    H404 = bst.HXutility('H404', ins=M404-0, outs='H404Out', T=5 + 273.15, cool_only=True)

    U402 = u.ResinColumn(
        'U402',
        ins=(U401-0, H402-0, H403-0, H404-0),
        outs=('U402Out', 'FlowthroughWaste', 'WashWaste', 'RegenerationWaste'),
        preset='IonExchange',
        TargetProduct_IDs=c.chemical_groups['LegHbIngredients'],
        BoundImpurity_IDs=c.chemical_groups['BoundImpurities'],
    )

    @M402.add_specification(run=True)
    def update_IXEquilibriumBuffer_initial():
        M402.ins[1].imass['H2O'] = (U401-0).imass['H2O'] * U402.wash_CV
        M402.ins[0].imol['IXEquilibriumBuffer'] = ((U401-0).imass['H2O'] * U402.wash_CV) * (0.025 + 0.01 + 0.001) / 1000

    @M403.add_specification(run=True)
    def update_IXElutionBuffer_initial():
        M403.ins[1].imass['H2O'] = (U401-0).imass['H2O'] * U402.elution_CV
        M403.ins[0].imol['IXElutionBuffer'] = ((U401-0).imass['H2O'] * U402.elution_CV) * (0.025 + 1 + 0.1) / 1000

    @M404.add_specification(run=True)
    def update_IXRegenerationSolution_initial():
        M404.ins[1].imass['H2O'] = (U401-0).imass['H2O'] * U402.regeneration_CV
        M404.ins[0].imol['NaOH'] = ((U401-0).imass['H2O'] * U402.regeneration_CV) * (0.5) / 1000

    M405 = bst.MixTank('M405', ins=(DfNanoBuffer, 'Water8'), outs='M405Out', tau=1)

    @M405.add_specification(run=True)
    def update_DfNanoBuffer_initial():
        M405.ins[1].imass['H2O'] = (U402-0).imass['H2O'] * 2
        M405.ins[0].imol['DfNanoBuffer'] = ((U402-0).imass['H2O'] * 2) * (0.01 + 0.01) / 1000

    H405 = bst.HXutility('H405', ins=M405-0, outs='H405Out', T=5 + 273.15, cool_only=True)

    U403 = u.Diafiltration.from_preset(
        'NF',
        'U403',
        ins=(U402-0, H405-0),
        outs=('U403Out', 'PermeateWasteNano'),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
    )
    U403.add_specification(run=True)

    U404 = u.Diafiltration(
        'U404',
        ins=(U403-0, bst.Stream('U404_buffer', H2O=0.001)),
        outs=('U404Out', 'PermeateWater'),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
        TMP_bar1=3,
        FeedWater_Recovery_to_Permeate=0.75,
    )

    U404.target_total_solids_percent = 12.0
    U404.target_legh_percent = 7.5

    @U404.add_specification(run=True)
    def U404_adjust_water_recovery():
        import flexsolve as flx
        import warnings

        feed_stream = U404.ins[0]
        product_stream = U404.outs[0]

        if feed_stream.F_mass <= 0:
            U404._run()
            return

        feed_legh = feed_stream.imass['Leghemoglobin']
        if feed_legh <= 0:
            U404._run()
            return

        def set_water_recovery(recovery):
            recovery = max(0.05, min(0.99, recovery))
            U404.FeedWater_Recovery_to_Permeate = recovery
            U404.Salt_Retention = 0.05

        def calculate_legh_error(water_recovery):
            set_water_recovery(water_recovery)
            U404._run()
            product_total_mass = product_stream.F_mass
            if product_total_mass <= 0:
                return -100.0
            product_legh_mass = product_stream.imass['Leghemoglobin']
            actual_legh_percent = (product_legh_mass / product_total_mass) * 100
            return actual_legh_percent - U404.target_legh_percent

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                optimal_recovery = flx.IQ_interpolation(
                    f=calculate_legh_error,
                    x0=0.50,
                    x1=0.99,
                    x=0.90,
                    xtol=0.001,
                    ytol=0.1,
                    maxiter=50,
                )

            set_water_recovery(optimal_recovery)
            U404._run()

        except Exception:
            set_water_recovery(0.90)
            U404._run()

    T501 = u.SulfuricAcidStorageTank(
        'T501',
        ins=bst.Stream('SulfuricAcid', H2SO4=0.98, H2O=0.02, units='kg/hr', price=price['H2SO4']),
        outs='T501Out',
    )

    @T501.add_specification(run=True)
    def update_acid_flowrate():
        T501.ins[0].imol['H2SO4'] = U402.outs[3].imol['NaOH'] / 2 * 1.001
        T501.ins[0].T = 25 + 273.15

    M502 = bst.NeutralizationTank1('M502', ins=(U402-3, T501-0), outs='M502Out', T=20 + 273.15)

    return {
        'units': [M401, H401, U401, M402, M403, M404, H402, H403, H404, U402, M405, H405, U403, U404, T501, M502],
        'concentrated_product': U404.outs[0],
        'uf_permeate': U401.outs[1],
        'ix_flowthrough': U402.outs[1],
        'ix_wash': U402.outs[2],
        'ix_regen_waste': U402.outs[3],
        'nf_permeate': U403.outs[1],
        'final_permeate': U404.outs[1],
        'neutralized_regen_waste': M502.outs[0],
        'M401': M401, 'H401': H401, 'U401': U401,
        'M402': M402, 'M403': M403, 'M404': M404,
        'H402': H402, 'H403': H403, 'H404': H404,
        'U402': U402, 'M405': M405, 'H405': H405,
        'U403': U403, 'U404': U404, 'T501': T501, 'M502': M502,
    }


def create_area_600_formulation(concentrated_product, LegHb_3):
    """
    Area 600: Final formulation / cooling (Config 2).
    """
    DilutionWater = bst.Stream('DilutionWater', H2O=1.0, units='kg/hr')

    M406 = bst.MixTank('M406', ins=(concentrated_product, DilutionWater), outs='FormulatedProduct', tau=0.1)
    M406.target_legh_percent = 7.5
    M406.max_total_solids_percent = 24.0

    @M406.add_specification(run=True)
    def adjust_total_solids():
        feed = M406.ins[0]
        legh_mass = feed.imass['Leghemoglobin']
        if legh_mass <= 0:
            M406.ins[1].imass['H2O'] = 0
            return
        target_total = legh_mass / (M406.target_legh_percent / 100)
        required_water = max(target_total - feed.F_mass, 0)
        total_mass = feed.F_mass + required_water
        total_solids = (feed.F_mass - feed.imass['H2O']) / total_mass * 100
        if total_solids > M406.max_total_solids_percent:
            target_total = (feed.F_mass - feed.imass['H2O']) / (M406.max_total_solids_percent / 100)
            required_water = max(target_total - feed.F_mass, 0)
        M406.ins[1].imass['H2O'] = required_water

    H406 = bst.HXutility('H406', ins=M406-0, outs=LegHb_3, T=0 + 273.15, cool_only=True)

    return {
        'units': [M406, H406],
        'product': H406.outs[0],
        'M406': M406, 'H406': H406,
        'DilutionWater': DilutionWater,
    }


def create_area_900_facilities(area_400, area_500, effluent1, effluent2, effluent3, use_area_convention=False):
    """
    Area 900: Facilities & utilities (Config 2).
    """
    M905 = bst.MixTank(
        'M905',
        ins=(
            area_400['press_liquor'],
            area_400['spent_media'],
            area_400['wash_effluent'],
            area_500['uf_permeate'],
            area_500['ix_flowthrough'],
            area_500['ix_wash'],
            area_500['nf_permeate'],
            area_500['final_permeate'],
        ),
        outs='WWTCombined',
        tau=1,
    )

    @M905.add_specification(run=True)
    def update_wwt_neutralization():
        out = M905.outs[0]
        h2so4 = out.imol['H2SO4']
        if h2so4 > 0:
            out.imol['Na2SO4'] += h2so4
            out.imol['H2O'] += 2 * h2so4
            out.imol['H2SO4'] = 0

    S904 = bst.Splitter('S904', ins=M905-0, outs=('WWTFeed', 'AcidWaste'), split=1)

    @S904.add_specification(run=True)
    def route_acid_to_disposal():
        S904.isplit['H2SO4'] = 0.0

    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[S904-0],
        outs=('biogas', 'sludge', 'RO_treated_water', effluent1),
        mockup=True,
        area=500,
    )
    biogas, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs

    M904 = bst.Mixer('M904', ins=(area_500['neutralized_regen_waste'], S904-1), outs=effluent2)

    S902 = bst.Splitter('S902', ins=sludge, outs=('SludgeToBoiler', 'SludgePurge'), split=1)
    S903 = bst.Splitter('S903', ins=biogas, outs=('BiogasToBoiler', effluent3), split=1)

    CT = bst.CoolingTower(500 if use_area_convention else 'CT')
    CWP = bst.ChilledWaterPackage(500 if use_area_convention else 'CWP')

    M903 = bst.Mixer('M903', ins=(area_400['dehydrated_debris'], S902-0), outs='SolidsToBoiler')

    BT = u.BoilerTurbogenerator(400 if use_area_convention else 'BT',
        (M903-0, S903-0, 'boiler_makeup_water', 'natural_gas', 'lime_boiler', 'boiler_chems'),
        outs=('emissions', 'rejected_water_and_blowdown', 'ash_disposal'),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85,
        satisfy_system_electricity_demand=False,
    )

    makeup_water_streams = (
        F.cooling_tower_makeup_water,
        F.Water1, F.Water2,
        F.Water3, F.Water4,
        F.Water5, F.Water6,
        F.Water7, F.Water8,
        F.boiler_makeup_water,
    )
    process_water_streams = (
        treated_water,
        *makeup_water_streams,
    )

    makeup_water = bst.Stream('makeup_water', price=0.000254)

    PWC = bst.ProcessWaterCenter(500 if use_area_convention else 'PWC',
        ins=(treated_water, makeup_water, 'recycled_process_water', 'makeup_process_water'),
        outs=('RO_water', 'process_water', 'excess_water'),
        makeup_water_streams=makeup_water_streams,
        process_water_streams=process_water_streams,
        reverse_osmosis_water_price=0.000254,
        process_water_price=0.000135,
    )

    return {
        'units': [M905, S904, M904, S902, S903, M903, CT, CWP, BT, PWC],
        'effluent': effluent1,
        'wastewater_treatment_sys': wastewater_treatment_sys,
        'M905': M905, 'S904': S904, 'M904': M904,
        'S902': S902, 'S903': S903,
        'M903': M903, 'CT': CT, 'CWP': CWP, 'BT': BT, 'PWC': PWC,
    }

@bst.SystemFactory(
    ID='LegHb_sys',
    ins=[s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt,s.DfUltraBuffer, s.IXEquilibriumBuffer, s.IXElutionBuffer
         ,s.IXRegenerationSolution, s.DfNanoBuffer],
    outs=[s.LegHb_3, s.vent1, s.vent2, s.effluent1, s.effluent2, s.effluent3,
    ],
    fthermo=lambda chemicals=None: LEGHB_THERMO,
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
    bst.preferences.N = 50

    SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, IXEquilibriumBuffer, IXElutionBuffer, IXRegenerationSolution, DfNanoBuffer = ins
    LegHb_3, vent1, vent2, effluent1, effluent2, effluent3 = outs

    bst.settings.set_thermo(LEGHB_THERMO, skip_checks=True)

    set_GWPCF(Glucose, 'Glucose')
    set_GWPCF_Multi(SeedIn1, ['AmmoniumSulfate', 'Glucose', 'MagnesiumSulfate', 'KH2PO4'],
                   [0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(SeedIn2, ['AmmoniumSulfate', 'Glucose', 'MagnesiumSulfate', 'KH2PO4'],
                   [0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(CultureIn, ['Glycine', 'Glucose', 'IronSulfate'],
                   [0.1/(0.1+60+0.15191), 60/(0.1+60+0.15191), 0.15191/(0.1+60+0.15191)])
    set_GWPCF(NH3_25wt, 'Ammonia_SEA', dilution=0.25)
    set_GWPCF_Multi(DfUltraBuffer, ['KH2PO4', 'NaCl'], [0.8472, 0.1455])
    set_GWPCF_Multi(IXEquilibriumBuffer, ['KH2PO4', 'NaCl'], [0.8472, 0.1455])
    set_GWPCF_Multi(IXElutionBuffer, ['KH2PO4', 'NaCl', 'KCl'], [0.3616, 0.5911, 0.0792])
    set_GWPCF(IXRegenerationSolution, 'NaOH')
    set_GWPCF_Multi(DfNanoBuffer, ['Na2HPO4', 'NaH2PO4'], [0.5420, 0.4580])

    load_process_settings()

    params = get_fermentation_parameters()
    reactions = create_fermentation_reactions(params)

    area_200 = create_area_200_media_prep(SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt)

    area_300 = create_area_300_conversion(
        seed_in=area_200['seed_out'],
        glucose_in=area_200['glucose_out'],
        ammonia_in=area_200['ammonia_out'],
        vent1=vent1,
        vent2=vent2,
        reactions=reactions,
        params=params,
    )

    area_400 = create_area_400_recovery(
        broth_in=area_300['broth_out'],
        DfUltraBuffer_wash=DfUltraBuffer,
    )

    area_500 = create_area_500_purification(
        clarified_lysate=area_400['clarified_lysate'],
        DfUltraBuffer=DfUltraBuffer,
        IXEquilibriumBuffer=IXEquilibriumBuffer,
        IXElutionBuffer=IXElutionBuffer,
        IXRegenerationSolution=IXRegenerationSolution,
        DfNanoBuffer=DfNanoBuffer,
    )

    area_600 = create_area_600_formulation(
        concentrated_product=area_500['concentrated_product'],
        LegHb_3=LegHb_3,
    )

    area_900 = create_area_900_facilities(
        area_400=area_400,
        area_500=area_500,
        effluent1=effluent1,
        effluent2=effluent2,
        effluent3=effluent3,
        use_area_convention=use_area_convention,
    )

    s.update_all_input_stream_prices(
        streamlist=[SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, IXEquilibriumBuffer, IXElutionBuffer, IXRegenerationSolution, DfNanoBuffer]
    )

    return LegHb_3, vent1, vent2, effluent1, effluent2, effluent3, DfUltraBuffer, IXEquilibriumBuffer, IXElutionBuffer, IXRegenerationSolution, DfNanoBuffer

# %% Design Specification Functions


def set_production_rate(system, target_production_rate_kg_hr):
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
    
    # Get product stream
    product_stream = system.flowsheet.stream.LegHb_3
    
    # Store baseline flow rates for all input streams
    baseline_flows = {}
    for stream in system.ins:
        if stream.F_mass > 0:
            baseline_flows[stream] = stream.F_mass
    
    if not baseline_flows:
        raise ValueError("No input streams with positive flow rates found")
    
    print(f"\n{'='*60}")
    print(f"Setting Production Rate to {target_production_rate_kg_hr:.2f} kg/hr")
    print(f"{'='*60}")
    print(f"Baseline input streams stored: {len(baseline_flows)} streams")
    
    # Run initial simulation to get baseline
    print("Running initial simulation to establish baseline...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        system.simulate()
    
    initial_production = product_stream.F_mass
    
    if initial_production <= 0:
        raise ValueError("Initial production rate is zero. Cannot scale system.")
    
    # Calculate initial scaling factor guess
    initial_guess = target_production_rate_kg_hr / initial_production
    
    print(f"  Initial production: {initial_production:.2f} kg/hr")
    print(f"  Target production:  {target_production_rate_kg_hr:.2f} kg/hr")
    print(f"  Initial scaling guess: {initial_guess:.4f}x")
    print(f"  Search bounds: 0.1x to 5.0x baseline")
    
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
                print(f"    Iteration {iteration[0]}: scale={scaling_factor:.4f}x, "
                      f"production={achieved_rate:.2f} kg/hr, error={error:.4f} kg/hr")
            
            return error
            
        except Exception as e:
            print(f"    Iteration {iteration[0]} FAILED at scale={scaling_factor:.4f}x: {e}")
            # Return large error to guide solver away from infeasible region
            return 1e6 if scaling_factor > initial_guess else -1e6
    
    try:
        print(f"\nSolving for optimal scaling factor (intermediate warnings suppressed)...")
        
        # Suppress all warnings during the iterative solving process
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Use flexsolve directly (not as a system specification)
            scaling_factor = flx.IQ_interpolation(
                f=objective_function,
                x0=0.1,              # Lower bound
                x1=5.0,              # Upper bound
                x=initial_guess,     # Initial guess
                xtol=0.0001,         # Tolerance on scaling factor
                ytol=0.01,           # Tolerance on production rate [kg/hr]
                maxiter=100,
                checkbounds=True,
                checkiter=True,
            )
        
        # After solver converges, run final simulation WITHOUT suppressing warnings
        # This allows any persistent issues to be displayed
        print(f"\nSolver converged. Running final validation simulation...")
        for stream, baseline_flow in baseline_flows.items():
            stream.F_mass = baseline_flow * scaling_factor
        
        # Final simulation with warnings enabled
        system.simulate()
        
        achieved_rate = product_stream.F_mass
        
        print(f"\n✓ Successfully achieved production rate:")
        print(f"  Target:          {target_production_rate_kg_hr:.2f} kg/hr")
        print(f"  Achieved:        {achieved_rate:.2f} kg/hr")
        print(f"  Scaling factor:  {scaling_factor:.4f}x")
        print(f"  Error:           {abs(achieved_rate - target_production_rate_kg_hr):.4f} kg/hr")
        print(f"  Total iterations: {iteration[0]}")
        print(f"{'='*60}\n")
        
        return achieved_rate
        
    except Exception as e:
        print(f"\n✗ Failed to achieve target production rate: {e}")
        print(f"  Error type: {type(e).__name__}")
        print(f"  Total iterations attempted: {iteration[0]}")
        
        # Restore baseline flows
        print(f"\n  Restoring baseline input flows...")
        for stream, baseline_flow in baseline_flows.items():
            stream.F_mass = baseline_flow
        
        # Run one final simulation with baseline flows (warnings enabled)
        try:
            system.simulate()
            print(f"  System restored to baseline: {product_stream.F_mass:.2f} kg/hr")
        except Exception as restore_error:
            print(f"  Warning: Could not restore baseline simulation: {restore_error}")
        
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
        print(r1)
        r2 = bst.report.lca_displacement_allocation_table(
            systems=[sys],
            key='GWP',
            items=[ss.LegHb_3],
        )
        print(r2)
        print("   LCA Displacement Allocation Table generated")
    except Exception as e:
        print(f"   LCA analysis failed: {e}")
    
    print(f"\n{'='*85}")
    print("SIMULATION COMPLETE")
    print(f"Target Production:   {TARGET_PRODUCTION:.2f} kg/hr")
    print(f"Achieved Production: {ss.LegHb_3.F_mass:.2f} kg/hr")
    print(f"{'='*85}\n")
