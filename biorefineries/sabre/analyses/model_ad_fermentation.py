# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
biosteam.evaluation.Model construction for the standalone ad_fermentation
flowsheet (systems._ad_fermentation_system.create_ad_fermentation_system()).

Press/PressateConcentrator/Mill/AcidogenicAD/DigestateScrewPress/
VFAMicrofilter parameters are all reused directly from model_biostimulant.py
and model_ad_vfa.py -- ad_fermentation embeds create_ad_vfa_system()
internally (with add_product_splitter=False, see
systems/_ad_fermentation_system.py), so those units are the exact same
objects, not redefined here. add_fermentation_medium_tank_parameters(),
add_yarrowia_lipid_fermenter_parameters(), add_evaporator_parameters(),
add_oil_extraction_parameters(), add_oil_centrifuge_parameters(), and
add_biomass_recycle_splitter_parameters() are new here (the fermentation
train has no smaller sabre flowsheet to reuse them from) -- see
docs/superpowers/specs/2026-07-31-sabre-uncertainty-model-design.md for the
general design.
"""
import biosteam as bst
from biosteam.evaluation import Metric

from biorefineries.sabre.utils import load_assumptions
from biorefineries.sabre.systems._ad_fermentation_system import create_ad_fermentation_system
from biorefineries.sabre._tea import solve_product_msp
from biorefineries.sabre.analyses.model_utils import distribution_from_yaml, add_tea_parameters
from biorefineries.sabre.analyses.model_biostimulant import (
    add_press_parameters, add_pressate_concentrator_parameters,
)
from biorefineries.sabre.analyses.model_ad_vfa import (
    add_mill_parameters, add_acidogenic_ad_parameters,
    add_digestate_screw_press_parameters, add_vfa_microfilter_parameters,
)

__all__ = (
    'add_fermentation_medium_tank_parameters', 'add_yarrowia_lipid_fermenter_parameters',
    'add_evaporator_parameters', 'add_oil_extraction_parameters',
    'add_oil_centrifuge_parameters', 'add_biomass_recycle_splitter_parameters',
    'create_ad_fermentation_model',
)

_FERMENTATION_YAML = load_assumptions("fermentation.yaml")
_FERMENTATION_MEDIUM_TANK = _FERMENTATION_YAML["fermentation_medium_tank"]
_VFA_FERM = _FERMENTATION_YAML["vfa"]
_VFA_CASE = _VFA_FERM["cases"][_VFA_FERM["case"]]

_DOWNSTREAM_PROCESSING_YAML = load_assumptions("downstream_processing.yaml")
_OIL_EXTRACTION = _DOWNSTREAM_PROCESSING_YAML["oil_extraction"]

_TEA_PRICE = load_assumptions("tea.yaml")["price"]


def add_fermentation_medium_tank_parameters(model, T601):
    """
    Add FermentationMediumTank process parameters to `model`
    (data/fermentation.yaml `fermentation_medium_tank`): tau, kW_per_m3
    (the tank's own attribute name -- the yaml/constructor-kwarg name is
    `mixing_kW_per_m3`), magnesium_sulfate_dose_kg_per_m3.

    Excludes target_pH (recorded in design_results, not enforced --
    permanently inert, same category as PressateConcentrator's
    operating_pressure_bar in model_biostimulant.py), and
    ammonia_dose_kg_per_m3/phosphate_dose_kg_per_m3/base_dose_kg_per_m3
    (all baseline 0.0 in the current yaml -- a +/-50%-of-baseline rule on
    a baseline of exactly 0 gives a zero-width Uniform(0, 0), i.e. no
    perturbation at all, so these are excluded as degenerate rather than
    silently sampled as a no-op).
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_FERMENTATION_MEDIUM_TANK[attr])
        def setter(x, attr=attr):
            setattr(T601, attr, float(x))
        param(
            setter, element=T601, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    baseline = float(_FERMENTATION_MEDIUM_TANK["tau"])
    def set_tau(x):
        T601.tau = float(x)
    param(
        set_tau, element=T601, name='Fermentation medium tank residence time', units='hr',
        baseline=baseline, distribution=distribution_from_yaml(baseline),
    )

    baseline = float(_FERMENTATION_MEDIUM_TANK["mixing_kW_per_m3"])
    def set_kW_per_m3(x):
        T601.kW_per_m3 = float(x)
    param(
        set_kW_per_m3, element=T601, name='Fermentation medium tank mixing intensity', units='kW/m3',
        baseline=baseline, distribution=distribution_from_yaml(baseline),
    )

    _add('Fermentation medium tank MgSO4 dose', 'magnesium_sulfate_dose_kg_per_m3', units='kg/m3')

    return model


def add_yarrowia_lipid_fermenter_parameters(model, R601):
    """
    Add YarrowiaLipidFermenter process parameters to `model`
    (data/fermentation.yaml `vfa.cases.yarrowia_vfa_base`): tau,
    conversion, product_yield_kg_per_kg_vfa_consumed,
    biomass_yield_kg_per_kg_vfa_consumed, co2_yield_kg_per_kg_vfa_consumed,
    oxygen_kg_per_kg_vfa_consumed, Q_O2_consumption.

    Excludes target_pH (recorded in design_results, not enforced --
    permanently inert), and seed_water_kgph/seed_cellmass_kgph (both
    baseline 0.0 -- degenerate under the +/-50% rule, same reasoning as
    the fermentation medium tank's zero-baseline reagent doses).
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_VFA_CASE[attr])
        def setter(x, attr=attr):
            setattr(R601, attr, float(x))
        param(
            setter, element=R601, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Fermenter residence time', 'tau', units='hr')
    _add('Fermenter VFA conversion', 'conversion')
    _add('Fermenter product yield', 'product_yield_kg_per_kg_vfa_consumed', units='kg product/kg VFA')
    _add('Fermenter biomass yield', 'biomass_yield_kg_per_kg_vfa_consumed', units='kg biomass/kg VFA')
    _add('Fermenter CO2 yield', 'co2_yield_kg_per_kg_vfa_consumed', units='kg CO2/kg VFA')
    _add('Fermenter oxygen demand', 'oxygen_kg_per_kg_vfa_consumed', units='kg O2/kg VFA')
    _add('Fermenter heat of reaction', 'Q_O2_consumption', units='kJ/kmol O2')

    return model


def add_evaporator_parameters(model, Ev607):
    """
    Add the post-fermentation MultiEffectEvaporator's sabre-specific
    target concentration to `model` (data/downstream_processing.yaml
    `oil_extraction.target_oil_and_solids_content_g_per_L`) --
    `target_oil_and_solids_content` is set directly by
    systems/_ad_fermentation_system.py, not a MultiEffectEvaporator
    constructor argument, but it is a plain settable instance attribute.
    """
    param = model.parameter

    baseline = float(_OIL_EXTRACTION["target_oil_and_solids_content_g_per_L"])
    def setter(x):
        Ev607.target_oil_and_solids_content = float(x)
    param(
        setter, element=Ev607, name='Evaporator target oil and solids content', units='g/L',
        baseline=baseline, distribution=distribution_from_yaml(baseline),
    )

    return model


def add_oil_extraction_parameters(model, OE):
    """
    Add OilExtraction process parameters to `model`
    (data/downstream_processing.yaml `oil_extraction`):
    homogenization_kWh_per_kg_dry_biomass. Excludes
    oil_extraction_reagent_usd_per_kg_oil -- OE._design() reads it
    directly off the module-level yaml dict each call
    (`_OIL_EXTRACTION["oil_extraction_reagent_usd_per_kg_oil"]`), not off
    an instance attribute, so there is nothing on the unit itself to set;
    parameterizing it would require a code change to
    units/_downstream_processing.py, out of scope here (same category as
    VFAMicrofilter's membrane_cost_usd_per_m2 in model_ad_vfa.py, which is
    baked into a class-level `@cost` decorator instead).
    """
    param = model.parameter

    baseline = float(_OIL_EXTRACTION["homogenization_kWh_per_kg_dry_biomass"])
    def setter(x):
        OE.homogenization_kWh_per_kg_dry_biomass = float(x)
    param(
        setter, element=OE, name='Oil extraction electricity intensity', units='kWh/kg dry biomass',
        baseline=baseline, distribution=distribution_from_yaml(baseline),
    )

    return model


def add_oil_centrifuge_parameters(model, C603_2):
    """
    Add the post-extraction LiquidsSplitCentrifuge's split fractions to
    `model` (data/downstream_processing.yaml `oil_extraction`):
    oil_recovery (-> isplit['MicrobialOil']), oil_water_split (->
    isplit['Water']).
    """
    param = model.parameter

    baseline = float(_OIL_EXTRACTION["oil_recovery"])
    def set_oil_recovery(x):
        C603_2.isplit['MicrobialOil'] = float(x)
    param(
        set_oil_recovery, element=C603_2, name='Oil centrifuge oil recovery', units='',
        baseline=baseline, distribution=distribution_from_yaml(baseline),
    )

    baseline = float(_OIL_EXTRACTION["oil_water_split"])
    def set_oil_water_split(x):
        C603_2.isplit['Water'] = float(x)
    param(
        set_oil_water_split, element=C603_2, name='Oil centrifuge water split', units='',
        baseline=baseline, distribution=distribution_from_yaml(baseline),
    )

    return model


def add_biomass_recycle_splitter_parameters(model, S602):
    """
    Add the biomass recycle MockSplitter's process parameters to `model`
    (data/downstream_processing.yaml `oil_extraction`):
    recycle_total_fraction, recycle_cellmass_wt_frac.
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_OIL_EXTRACTION[attr])
        def setter(x, attr=attr):
            setattr(S602, attr, float(x))
        param(
            setter, element=S602, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Biomass recycle total fraction', 'recycle_total_fraction')
    _add('Biomass recycle cellmass wt fraction', 'recycle_cellmass_wt_frac')

    return model


def create_ad_fermentation_model(system=None):
    """
    Build the full biosteam.evaluation.Model for the standalone
    ad_fermentation flowsheet: reused Press/PressateConcentrator/Mill/
    AcidogenicAD/DigestateScrewPress/VFAMicrofilter parameters (from the
    embedded ad_vfa subsystem, built with add_product_splitter=False),
    the fermentation train's own process parameters (medium tank,
    fermenter, evaporator, oil extraction, oil centrifuge, biomass
    recycle splitter), the shared TEA parameters, and this system's own
    boundary-stream prices (sargassum feed, biostimulant product,
    permeate disposal, acidogenic solid digestate disposal, VFA cake
    disposal, ammonia/phosphate/NaOH/MgSO4 reagent prices, microbial oil
    product, wastewater disposal).

    Parameters
    ----------
    system : bst.System, optional
        Defaults to a fresh `create_ad_fermentation_system()`.
    """
    if system is None:
        system = create_ad_fermentation_system()
    system.simulate()

    flowsheet = system.flowsheet
    PR = flowsheet.unit.PR
    PC = flowsheet.unit.PC
    ML = flowsheet.unit.ML
    AD = flowsheet.unit.VFA_AD
    SP = flowsheet.unit.SP_VFA
    MF = flowsheet.unit.MF
    T601 = flowsheet.unit.T601
    R601 = flowsheet.unit.R601
    Ev607 = flowsheet.unit.Ev607
    OE = flowsheet.unit.OE
    C603_2 = flowsheet.unit.C603_2
    S602 = flowsheet.unit.S602

    feed = flowsheet.stream.sargassum_feed
    biostimulant_product = flowsheet.stream.biostimulant_product
    permeate = flowsheet.stream.permeate
    acidogenic_solid_digestate = flowsheet.stream.acidogenic_solid_digestate
    vfa_cake = flowsheet.stream.vfa_cake
    ammonia = flowsheet.stream.fermentation_ammonia
    phosphate = flowsheet.stream.fermentation_phosphate
    base = flowsheet.stream.fermentation_base
    mgso4 = flowsheet.stream.fermentation_mgso4
    microbial_oil = flowsheet.stream.microbial_oil
    wastewater = flowsheet.stream.wastewater

    tea = system.TEA

    def get_msp():
        return solve_product_msp(tea=tea, product_stream=microbial_oil)["usd_per_kg"]

    def get_fci():
        return tea.FCI

    def get_electricity():
        return sum(u.power_utility.rate for u in system.units)

    def get_production_rate():
        return microbial_oil.F_mass

    metrics = [
        Metric('MSP', get_msp, 'USD/kg'),
        Metric('Fixed capital investment', get_fci, 'USD'),
        Metric('Total electricity', get_electricity, 'kW'),
        Metric('Microbial oil production rate', get_production_rate, 'kg/hr'),
    ]

    model = bst.evaluation.Model(system, metrics)

    add_press_parameters(model, PR)
    add_pressate_concentrator_parameters(model, PC)
    add_mill_parameters(model, ML)
    add_acidogenic_ad_parameters(model, AD)
    add_digestate_screw_press_parameters(model, SP)
    add_vfa_microfilter_parameters(model, MF)
    add_fermentation_medium_tank_parameters(model, T601)
    add_yarrowia_lipid_fermenter_parameters(model, R601)
    add_evaporator_parameters(model, Ev607)
    add_oil_extraction_parameters(model, OE)
    add_oil_centrifuge_parameters(model, C603_2)
    add_biomass_recycle_splitter_parameters(model, S602)
    add_tea_parameters(model, tea)

    param = model.parameter

    def _add_price(name, stream, price_entry):
        baseline = float(price_entry["baseline"])
        def setter(x, stream=stream):
            stream.price = float(x)
        param(
            setter, element=stream, name=name, units='USD/kg',
            baseline=baseline, distribution=distribution_from_yaml(price_entry),
        )

    _add_price('Sargassum feed price', feed, _TEA_PRICE["sargassum"])
    _add_price('Biostimulant product price', biostimulant_product, _TEA_PRICE["biostimulant"])
    _add_price('Permeate disposal price', permeate, _TEA_PRICE["disposal_wastewater"])
    _add_price(
        'Acidogenic solid digestate disposal price', acidogenic_solid_digestate,
        _TEA_PRICE["disposal_solid"],
    )
    _add_price('VFA cake disposal price', vfa_cake, _TEA_PRICE["disposal_solid"])
    _add_price('Ammonia price', ammonia, _TEA_PRICE["ammonia"])
    _add_price('Phosphate price', phosphate, _TEA_PRICE["phosphate"])
    _add_price('NaOH price', base, _TEA_PRICE["naoh"])
    _add_price('MgSO4 price', mgso4, _TEA_PRICE["mgso4"])
    _add_price('Microbial oil product price', microbial_oil, _TEA_PRICE["microbial_oil"])
    _add_price('Wastewater disposal price', wastewater, _TEA_PRICE["disposal_wastewater"])

    return model
