# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
biosteam.evaluation.Model construction for the standalone ad_vfa flowsheet
(systems._ad_vfa_system.create_ad_vfa_system()).

add_mill_parameters(), add_acidogenic_ad_parameters(),
add_digestate_screw_press_parameters(), and add_vfa_microfilter_parameters()
add only process/cost parameters (no stream prices), so ad_biomethane's
(future) model file can reuse add_mill_parameters()/
add_digestate_screw_press_parameters() against its own embedded Mill/
DigestateScrewPress instances, and ad_fermentation's (future) model file can
reuse all four against the ad_vfa subsystem it embeds -- see
docs/superpowers/specs/2026-07-31-sabre-uncertainty-model-design.md for the
general design and biorefineries/sabre/systems/_ad_vfa_system.py for how
this flowsheet is built (it embeds create_biostimulant_system() internally,
so Press/PressateConcentrator parameters are reused directly from
model_biostimulant.py, not redefined here).
"""
import biosteam as bst
from biosteam.evaluation import Metric

from biorefineries.sabre.utils import load_assumptions
from biorefineries.sabre.systems._ad_vfa_system import create_ad_vfa_system
from biorefineries.sabre._tea import solve_product_msp
from biorefineries.sabre.analyses.model_utils import distribution_from_yaml, add_tea_parameters
from biorefineries.sabre.analyses.model_biostimulant import (
    add_press_parameters, add_pressate_concentrator_parameters,
)

__all__ = (
    'add_mill_parameters', 'add_acidogenic_ad_parameters',
    'add_digestate_screw_press_parameters', 'add_vfa_microfilter_parameters',
    'create_ad_vfa_model',
)

# GAL_TO_M3 matches units/_ad.py's own conversion constant -- AnaerobicDigester
# stores max_single_digester_volume_MG (yaml's million-US-gallon basis) as
# max_single_digester_volume_m3 (converted at __init__ time), so this
# module must apply the same conversion to compute the right baseline/bounds
# for that attribute.
GAL_TO_M3 = 0.003785411784

_MILL = load_assumptions("preprocessing.yaml")["mill"]

_AD_YAML = load_assumptions("ad.yaml")
_AD_SHARED = _AD_YAML["ad"]
_AD_ACIDOGENIC = {**_AD_SHARED, **_AD_SHARED.get("acidogenic", {})}
_AD_PERFORMANCE = _AD_YAML["ad_performance"]
_ADP_ACIDOGENIC = _AD_PERFORMANCE["acidogenic"]
_ADP_ACIDOGENIC_CASE = _ADP_ACIDOGENIC["cases"][_ADP_ACIDOGENIC["case"]]
_DIGESTATE_SCREW_PRESS = _AD_YAML["digestate_screw_press"]

_DOWNSTREAM_PROCESSING_YAML = load_assumptions("downstream_processing.yaml")
_VFA_MICROFILTER = _DOWNSTREAM_PROCESSING_YAML["vfa_microfilter"]
_VFA_PRODUCT_SPLITTER = _DOWNSTREAM_PROCESSING_YAML["vfa_product_splitter"]

_TEA_PRICE = load_assumptions("tea.yaml")["price"]


def add_mill_parameters(model, ML):
    """
    Add Mill (ML) process/cost parameters to `model` (data/preprocessing.yaml
    `mill`): loss_frac, power_kWh_per_dry_ton_dry. Excludes capex_model
    (categorical) and the ref_capacity_dry_ton_per_hr/purchase_cost_ref_usd/
    scale_exponent CAPEX anchor set, and F_BM. Reusable by any sabre
    flowsheet that embeds a Mill (ad_vfa, ad_biomethane, ad_fermentation).
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_MILL[attr])
        def setter(x, attr=attr):
            setattr(ML, attr, float(x))
        param(
            setter, element=ML, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Mill loss fraction', 'loss_frac')
    _add('Mill power intensity', 'power_kWh_per_dry_ton_dry', units='kWh/dry ton')

    return model


def add_acidogenic_ad_parameters(model, AD):
    """
    Add AcidogenicAD process parameters to `model` (data/ad.yaml `ad`
    shared/acidogenic sizing and `ad_performance.acidogenic`'s selected
    case): vs_destruction, vfa_kg_per_kg_vs, hrt_days, headspace_frac,
    max single digester volume (m3, converted from the yaml's MG basis),
    mixing_W_per_m3, target_temperature_K.

    Excludes vfa_split (a per-chemical composition that must sum to 1,
    like tea.yaml's construction_schedule), digestible_IDs (categorical),
    target_feed_moisture_frac (null in the current baseline),
    maintenance_usd_per_m3_yr (AcidogenicAD's own __init__ never forwards
    it to the base class, so it's not wired to yaml for this class), the
    heat-shock parameters (enable_heat_shock is False in the current
    baseline, so AcidogenicAD._design() never reads
    hs_target_temperature_K/hs_events_per_day/hs_heated_fraction_of_liquid/
    hs_duration_min -- dead code, same pattern as EV in model_biostimulant.py),
    and target_temperature_K -- fixed at baseline (not sampled) per user
    decision 2026-08-01: a generic +/-50%-of-absolute-Kelvin-value rule
    produced physically nonsensical ranges (e.g. 308.15 K baseline ->
    [154, 462] K) and caused ~20% of a trial N=20 sample to fail
    simulation ("no cooling agent that can cool under 157 K"). Every
    other Kelvin temperature parameter in this package
    (EnzymaticPretreatment.temperature_K, HeatingPretreatment.
    target_temperature_K, used by ad_biomethane) is excluded for the same
    reason.

    Reusable by ad_fermentation's (future) model file, which embeds ad_vfa's
    own system.
    """
    param = model.parameter

    def _add(name, attr, source, units=''):
        baseline = float(source[attr])
        def setter(x, attr=attr):
            setattr(AD, attr, float(x))
        param(
            setter, element=AD, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Acidogenic AD VS destruction', 'vs_destruction', _ADP_ACIDOGENIC_CASE)
    _add('Acidogenic AD VFA yield', 'vfa_kg_per_kg_vs', _ADP_ACIDOGENIC_CASE, units='kg VFA/kg VS')
    _add('Acidogenic AD HRT', 'hrt_days', _AD_ACIDOGENIC, units='days')
    _add('Acidogenic AD headspace fraction', 'gas_storage_frac_of_total_volume', _AD_ACIDOGENIC)
    _add('Acidogenic AD mixing intensity', 'mixing_W_per_m3', _AD_ACIDOGENIC, units='W/m3')

    baseline_MG = float(_AD_ACIDOGENIC["max_single_digester_volume_MG"])
    baseline_m3 = baseline_MG * 1e6 * GAL_TO_M3
    def set_max_single_digester_volume_m3(x):
        AD.max_single_digester_volume_m3 = float(x)
    param(
        set_max_single_digester_volume_m3, element=AD,
        name='Acidogenic AD max single digester volume', units='m3',
        baseline=baseline_m3, distribution=distribution_from_yaml(baseline_m3),
    )

    # target_temperature_K: fixed at baseline, not sampled -- see the
    # Excludes note above.

    return model


def add_digestate_screw_press_parameters(model, SP):
    """
    Add DigestateScrewPress process/cost parameters to `model`
    (data/ad.yaml `digestate_screw_press`): ts_capture_frac,
    cake_moisture_frac, capacity_tph_each, kWh_per_m3. Excludes
    include_polymer_dosing (bool; False in the current baseline, so
    polymer_dosing_cost_usd_each is dead code) and F_BM. Reusable by
    ad_biomethane's (future) model file (its own DigestateScrewPress on
    the methanogenic-digestate side) and by ad_fermentation's (future)
    model file (via ad_vfa reuse).
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_DIGESTATE_SCREW_PRESS[attr])
        def setter(x, attr=attr):
            setattr(SP, attr, float(x))
        param(
            setter, element=SP, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Screw press TS capture fraction', 'ts_capture_frac')
    _add('Screw press cake moisture fraction', 'cake_moisture_frac')
    _add('Screw press capacity per unit', 'capacity_tph_each', units='ton/h')
    _add('Screw press electricity intensity', 'kWh_per_m3', units='kWh/m3')

    return model


def add_vfa_microfilter_parameters(model, MF):
    """
    Add VFAMicrofilter process/cost parameters to `model`
    (data/downstream_processing.yaml `vfa_microfilter`):
    vfa_to_permeate_frac, water_to_permeate_frac, solids_to_permeate_frac,
    dissolved_other_to_permeate_frac, SEC_kWh_per_m3_feed,
    design_flux_L_m2_h. Excludes vfa_IDs (categorical) and
    membrane_cost_usd_per_m2/F_BM (baked into the class-level `@cost`
    decorator at class-definition time, not a settable instance attribute).
    Reusable by ad_fermentation's (future) model file, via ad_vfa reuse.
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_VFA_MICROFILTER[attr])
        def setter(x, attr=attr):
            setattr(MF, attr, float(x))
        param(
            setter, element=MF, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('MF VFA to permeate fraction', 'vfa_to_permeate_frac')
    _add('MF water to permeate fraction', 'water_to_permeate_frac')
    _add('MF solids to permeate fraction', 'solids_to_permeate_frac')
    _add('MF dissolved other to permeate fraction', 'dissolved_other_to_permeate_frac')
    _add('MF electricity intensity', 'SEC_kWh_per_m3_feed', units='kWh/m3')
    _add('MF design flux', 'design_flux_L_m2_h', units='L/m2-h')

    return model


def create_ad_vfa_model(system=None):
    """
    Build the full biosteam.evaluation.Model for the standalone ad_vfa
    flowsheet: reused Press/PressateConcentrator parameters (from the
    embedded biostimulant subsystem), Mill/AcidogenicAD/
    DigestateScrewPress/VFAMicrofilter process parameters, the product
    splitter's VFA recovery fraction, the shared TEA parameters, and this
    system's own boundary-stream prices (sargassum feed, biostimulant
    product, permeate disposal, acidogenic solid digestate disposal, VFA
    cake disposal, pure VFA product, VFA disposal).

    Parameters
    ----------
    system : bst.System, optional
        Defaults to a fresh `create_ad_vfa_system()` (standalone, with the
        product splitter enabled -- `pressed_cake`'s disposal price is
        NOT parameterized here: once embedded in this system it feeds the
        Mill instead of leaving as a leaf stream, so BioSTEAM's TEA no
        longer counts its price -- confirmed via `system.feeds`/
        `system.products`, which do not include `pressed_cake`).
    """
    if system is None:
        system = create_ad_vfa_system()
    system.simulate()

    flowsheet = system.flowsheet
    PR = flowsheet.unit.PR
    PC = flowsheet.unit.PC
    ML = flowsheet.unit.ML
    AD = flowsheet.unit.VFA_AD
    SP = flowsheet.unit.SP_VFA
    MF = flowsheet.unit.MF
    SP_PRODUCT = flowsheet.unit.SP_PRODUCT

    feed = flowsheet.stream.sargassum_feed
    biostimulant_product = flowsheet.stream.biostimulant_product
    permeate = flowsheet.stream.permeate
    acidogenic_solid_digestate = flowsheet.stream.acidogenic_solid_digestate
    vfa_cake = flowsheet.stream.vfa_cake
    pure_vfa = flowsheet.stream.pure_vfa
    vfa_disposal = flowsheet.stream.vfa_disposal

    tea = system.TEA

    def get_msp():
        return solve_product_msp(tea=tea, product_stream=pure_vfa)["usd_per_kg"]

    def get_fci():
        return tea.FCI

    def get_electricity():
        return sum(u.power_utility.rate for u in system.units)

    def get_production_rate():
        return pure_vfa.F_mass

    metrics = [
        Metric('MSP', get_msp, 'USD/kg'),
        Metric('Fixed capital investment', get_fci, 'USD'),
        Metric('Total electricity', get_electricity, 'kW'),
        Metric('VFA production rate', get_production_rate, 'kg/hr'),
    ]

    model = bst.evaluation.Model(system, metrics)

    add_press_parameters(model, PR)
    add_pressate_concentrator_parameters(model, PC)
    add_mill_parameters(model, ML)
    add_acidogenic_ad_parameters(model, AD)
    add_digestate_screw_press_parameters(model, SP)
    add_vfa_microfilter_parameters(model, MF)
    add_tea_parameters(model, tea)

    param = model.parameter

    # SP_PRODUCT's split dict is built at system-construction time as
    # {cid: vfa_recovery_frac for cid in MF.vfa_IDs} (see
    # systems/_ad_vfa_system.py) -- one uniform recovery fraction applied
    # to every VFA chemical ID, so a single Parameter sets all of them
    # together via SP_PRODUCT.isplit[cid], not one Parameter per chemical.
    vfa_ids = tuple(MF.vfa_IDs)
    baseline = float(_VFA_PRODUCT_SPLITTER["vfa_recovery_frac"])
    def set_vfa_recovery_frac(x, vfa_ids=vfa_ids):
        for cid in vfa_ids:
            SP_PRODUCT.isplit[cid] = float(x)
    param(
        set_vfa_recovery_frac, element=SP_PRODUCT,
        name='Product splitter VFA recovery fraction', units='',
        baseline=baseline, distribution=distribution_from_yaml(baseline),
    )

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
    _add_price('Pure VFA product price', pure_vfa, _TEA_PRICE["vfa"])
    _add_price('VFA disposal price', vfa_disposal, _TEA_PRICE["disposal_wastewater"])

    return model
