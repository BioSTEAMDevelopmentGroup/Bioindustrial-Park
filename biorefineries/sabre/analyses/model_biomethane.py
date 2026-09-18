# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
biosteam.evaluation.Model construction for the standalone ad_biomethane
flowsheet (systems._ad_biomethane_system.create_ad_biomethane_system()).

Press/PressateConcentrator parameters are reused directly from
model_biostimulant.py, and Mill/DigestateScrewPress parameters are reused
directly from model_ad_vfa.py (both embed create_biostimulant_system()/a
Mill/a DigestateScrewPress the same way ad_biomethane does -- see
docs/superpowers/specs/2026-07-31-sabre-uncertainty-model-design.md).
add_methanogenic_ad_parameters(), add_h2s_removal_parameters(),
add_biogas_upgrading_parameters(), add_peroxide_pretreatment_parameters(),
add_heating_pretreatment_parameters(), and
add_enzymatic_pretreatment_parameters() are new here and add only
process/cost parameters (no stream prices), so ad_fermentation's (future)
model file or any other future flowsheet reusing these same units can
import them too.

pretreatment_case selects which of PeroxidePretreatment (PX),
HeatingPretreatment (HT), and EnzymaticPretreatment (EZ) actually exist in
the built system (see systems/_ad_biomethane_system.py) --
create_ad_biomethane_model() detects which of these units are present via
`unit_id in flowsheet.unit` and only adds parameters for the ones that
exist, so it works unmodified for all five pretreatment_case options.
"""
import biosteam as bst
from biosteam.evaluation import Metric

from biorefineries.sabre.utils import load_assumptions
from biorefineries.sabre.systems._ad_biomethane_system import create_ad_biomethane_system
from biorefineries.sabre._tea import solve_product_msp, usd_per_mmbtu_to_usd_per_kg, CH4_MMBTU_PER_KG
from biorefineries.sabre.analyses.model_utils import distribution_from_yaml, add_tea_parameters
from biorefineries.sabre.analyses.model_biostimulant import (
    add_press_parameters, add_pressate_concentrator_parameters,
)
from biorefineries.sabre.analyses.model_ad_vfa import (
    add_mill_parameters, add_digestate_screw_press_parameters,
)

__all__ = (
    'add_methanogenic_ad_parameters', 'add_h2s_removal_parameters',
    'add_biogas_upgrading_parameters', 'add_peroxide_pretreatment_parameters',
    'add_heating_pretreatment_parameters', 'add_enzymatic_pretreatment_parameters',
    'create_ad_biomethane_model',
)

_AD_YAML = load_assumptions("ad.yaml")
_AD_SHARED = _AD_YAML["ad"]
_AD_METHANOGENIC = {**_AD_SHARED, **_AD_SHARED.get("methanogenic", {})}
_AD_PERFORMANCE = _AD_YAML["ad_performance"]
_ADP_METHANOGENIC = {**_AD_PERFORMANCE, **_AD_PERFORMANCE.get("methanogenic", {})}
_AD_COST = _AD_SHARED["cost"]
_H2S_REMOVAL = _AD_YAML["h2s_removal"]
_BIOGAS_UPGRADING = _AD_YAML["biogas_upgrading"]

_PRETREATMENT_AD = load_assumptions("pretreatment.yaml")["pretreatment_ad"]
_PEROXIDE = _PRETREATMENT_AD["peroxide"]["peroxide"]
_HEATING = _PRETREATMENT_AD["combined_PTE"]["heating"]
_ENZYMATIC = _PRETREATMENT_AD["enzymatic"]["enzymatic"]

_TEA_PRICE = load_assumptions("tea.yaml")["price"]


def add_methanogenic_ad_parameters(model, AD):
    """
    Add MethanogenicAD process parameters to `model` (data/ad.yaml `ad`
    shared/methanogenic sizing, `ad.cost`, and the selected pretreatment
    case's `ad_effects`): vs_destruction, ch4_kg_per_kg_vs_fed, hrt_days,
    headspace_frac, max single digester volume (m3, converted from the
    yaml's MG basis), mixing_W_per_m3, maintenance_usd_per_m3_yr.

    Excludes raw_biogas_molfrac (a per-chemical composition that must sum
    to 1, like tea.yaml's construction_schedule), digestible_IDs
    (categorical), biodegradability (a per-chemical dict with no yaml
    range/distribution of its own -- same treatment as vfa_split),
    target_feed_moisture_frac (null in the current baseline), and
    target_temperature_K (fixed at baseline -- see model_ad_vfa.py's
    add_acidogenic_ad_parameters docstring for why).

    Reusable by any future sabre flowsheet that embeds a MethanogenicAD.
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

    # vs_destruction/ch4_kg_per_kg_vs_fed come from the selected
    # pretreatment case's ad_effects (data/pretreatment.yaml
    # `pretreatment_ad.<case>.ad_effects`) -- read directly off the live
    # unit instance instead, since which case is active depends on how
    # the caller built the system, not a fixed module-level yaml lookup.
    def _add_from_ad(name, attr, units=''):
        baseline = float(getattr(AD, attr))
        def setter(x, attr=attr):
            setattr(AD, attr, float(x))
        param(
            setter, element=AD, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add_from_ad('Methanogenic AD VS destruction', 'vs_destruction')
    _add_from_ad('Methanogenic AD CH4 yield', 'ch4_kg_per_kg_vs_fed', units='kg CH4/kg VS')
    _add('Methanogenic AD HRT', 'hrt_days', _AD_METHANOGENIC, units='days')
    _add('Methanogenic AD headspace fraction', 'gas_storage_frac_of_total_volume', _AD_METHANOGENIC)
    _add('Methanogenic AD mixing intensity', 'mixing_W_per_m3', _AD_METHANOGENIC, units='W/m3')
    _add('Methanogenic AD maintenance', 'maintenance_usd_per_m3_yr', _AD_COST, units='USD/m3-yr')

    baseline_MG = float(_AD_METHANOGENIC["max_single_digester_volume_MG"])
    baseline_m3 = baseline_MG * 1e6 * 0.003785411784  # GAL_TO_M3, matches units/_ad.py
    def set_max_single_digester_volume_m3(x):
        AD.max_single_digester_volume_m3 = float(x)
    param(
        set_max_single_digester_volume_m3, element=AD,
        name='Methanogenic AD max single digester volume', units='m3',
        baseline=baseline_m3, distribution=distribution_from_yaml(baseline_m3),
    )

    return model


def add_h2s_removal_parameters(model, H2SR):
    """
    Add H2SRemoval process/cost parameters to `model` (data/ad.yaml
    `h2s_removal`): h2s_removal_efficiency, reagent_cost_usd_per_Nm3_raw.
    Excludes technology (categorical), the ref_flow_Nm3ph/
    ref_installed_cost_usd/scale_exponent CAPEX anchor set (baked into the
    class-level `@cost` decorator), and F_BM.
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_H2S_REMOVAL[attr])
        def setter(x, attr=attr):
            setattr(H2SR, attr, float(x))
        param(
            setter, element=H2SR, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('H2S removal efficiency', 'h2s_removal_efficiency')
    _add('H2S removal reagent cost', 'reagent_cost_usd_per_Nm3_raw', units='USD/Nm3')

    return model


def add_biogas_upgrading_parameters(model, UP):
    """
    Add BiogasUpgrading process/cost parameters to `model` (data/ad.yaml
    `biogas_upgrading`): methane_loss_frac, min_ch4_massfrac,
    electricity_kwh_per_Nm3_raw, capex_usd_per_Nm3ph_raw,
    maintenance_frac_of_capex_per_yr. Excludes technology (categorical)
    and F_BM.
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_BIOGAS_UPGRADING[attr])
        def setter(x, attr=attr):
            setattr(UP, attr, float(x))
        param(
            setter, element=UP, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Biogas upgrading methane loss fraction', 'methane_loss_frac')
    _add('Biogas upgrading min CH4 mass fraction', 'min_ch4_massfrac')
    _add('Biogas upgrading electricity intensity', 'electricity_kWh_per_Nm3_raw', units='kWh/Nm3')
    _add('Biogas upgrading capex per Nm3/h', 'capex_usd_per_Nm3ph_raw', units='USD/(Nm3/h)')
    _add('Biogas upgrading maintenance fraction of capex', 'maintenance_frac_of_capex_per_yr')

    return model


def add_peroxide_pretreatment_parameters(model, PX):
    """
    Add PeroxidePretreatment process/cost parameters to `model`
    (data/pretreatment.yaml `pretreatment_ad.peroxide.peroxide`):
    h2o2_wt_frac_on_dry_feed, residence_time_hr, capex_usd,
    h2o2_price_usd_per_kg, maintenance_frac_of_capex_per_yr. Excludes
    F_BM. Only present when `pretreatment_case` is 'peroxide',
    'combined_PE', or 'combined_PTE'.
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_PEROXIDE[attr])
        def setter(x, attr=attr):
            setattr(PX, attr, float(x))
        param(
            setter, element=PX, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Peroxide pretreatment H2O2 dose', 'h2o2_wt_frac_on_dry_feed')
    _add('Peroxide pretreatment residence time', 'residence_time_hr', units='hr')
    _add('Peroxide pretreatment capex', 'capex_usd', units='USD')
    _add('Peroxide pretreatment H2O2 price', 'h2o2_price_usd_per_kg', units='USD/kg')
    _add('Peroxide pretreatment maintenance fraction of capex', 'maintenance_frac_of_capex_per_yr')

    return model


def add_heating_pretreatment_parameters(model, HT):
    """
    Add HeatingPretreatment process/cost parameters to `model`
    (data/pretreatment.yaml `pretreatment_ad.combined_PTE.heating`):
    residence_time_hr, capex_usd, maintenance_frac_of_capex_per_yr.
    Excludes target_temperature_K (fixed at baseline -- see
    model_ad_vfa.py's add_acidogenic_ad_parameters docstring) and F_BM.
    Only present when `pretreatment_case` is 'combined_PTE'.
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_HEATING[attr])
        def setter(x, attr=attr):
            setattr(HT, attr, float(x))
        param(
            setter, element=HT, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Heating pretreatment residence time', 'residence_time_hr', units='hr')
    _add('Heating pretreatment capex', 'capex_usd', units='USD')
    _add('Heating pretreatment maintenance fraction of capex', 'maintenance_frac_of_capex_per_yr')

    return model


def add_enzymatic_pretreatment_parameters(model, EZ):
    """
    Add EnzymaticPretreatment process/cost parameters to `model`
    (data/pretreatment.yaml `pretreatment_ad.enzymatic.enzymatic`):
    residence_time_hr, enzyme_dose_kg_per_kg_dry_feed, treated_fraction,
    enzyme_recycle_factor, capex_usd, enzyme_price_usd_per_kg,
    maintenance_frac_of_capex_per_yr. Excludes temperature_K (fixed at
    baseline) and F_BM. Only present when `pretreatment_case` is
    'enzymatic', 'combined_PE', or 'combined_PTE'.
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_ENZYMATIC[attr])
        def setter(x, attr=attr):
            setattr(EZ, attr, float(x))
        param(
            setter, element=EZ, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Enzymatic pretreatment residence time', 'residence_time_hr', units='hr')
    _add('Enzymatic pretreatment enzyme dose', 'enzyme_dose_kg_per_kg_dry_feed', units='kg/kg dry feed')
    _add('Enzymatic pretreatment treated fraction', 'treated_fraction')
    _add('Enzymatic pretreatment enzyme recycle factor', 'enzyme_recycle_factor')
    _add('Enzymatic pretreatment capex', 'capex_usd', units='USD')
    _add('Enzymatic pretreatment enzyme price', 'enzyme_price_usd_per_kg', units='USD/kg')
    _add('Enzymatic pretreatment maintenance fraction of capex', 'maintenance_frac_of_capex_per_yr')

    return model


def create_ad_biomethane_model(system=None, pretreatment_case: str = 'press_mill_only'):
    """
    Build the full biosteam.evaluation.Model for the standalone
    ad_biomethane flowsheet: reused Press/PressateConcentrator parameters
    (from the embedded biostimulant subsystem), reused Mill/
    DigestateScrewPress parameters, MethanogenicAD/H2SRemoval/
    BiogasUpgrading process parameters, whichever pretreatment unit(s)
    `pretreatment_case` includes (PeroxidePretreatment/HeatingPretreatment/
    EnzymaticPretreatment -- detected from the built system, not assumed),
    the shared TEA parameters, and this system's own boundary-stream
    prices (sargassum feed, biostimulant product, permeate disposal,
    spent H2S media disposal, methanogenic solid digestate disposal,
    liquid digestate disposal, biomethane).

    Parameters
    ----------
    system : bst.System, optional
        Defaults to a fresh `create_ad_biomethane_system(pretreatment_case=
        pretreatment_case)`. If a system is passed in directly, this
        argument is ignored -- the pretreatment units actually present in
        `system` are auto-detected instead.
    pretreatment_case : str
        Forwarded to `create_ad_biomethane_system()` when `system` is not
        given. One of 'press_mill_only', 'enzymatic', 'peroxide',
        'combined_PE', 'combined_PTE' (data/pretreatment.yaml
        `pretreatment_ad`).
    """
    if system is None:
        system = create_ad_biomethane_system(pretreatment_case=pretreatment_case)
    system.simulate()

    flowsheet = system.flowsheet
    PR = flowsheet.unit.PR
    PC = flowsheet.unit.PC
    ML = flowsheet.unit.ML
    AD = flowsheet.unit.AD
    H2SR = flowsheet.unit.H2SR
    UP = flowsheet.unit.UP
    SP = flowsheet.unit.SP

    feed = flowsheet.stream.sargassum_feed
    biostimulant_product = flowsheet.stream.biostimulant_product
    permeate = flowsheet.stream.permeate
    spent_h2s_media = flowsheet.stream.spent_h2s_media
    methanogenic_solid_digestate = flowsheet.stream.methanogenic_solid_digestate
    liquid_digestate = flowsheet.stream.liquid_digestate
    biomethane = flowsheet.stream.biomethane

    tea = system.TEA

    def get_msp_usd_per_kg():
        return solve_product_msp(
            tea=tea, product_stream=biomethane, energy_content_mmbtu_per_kg=CH4_MMBTU_PER_KG,
        )["usd_per_kg"]

    def get_msp_usd_per_mmbtu():
        return solve_product_msp(
            tea=tea, product_stream=biomethane, energy_content_mmbtu_per_kg=CH4_MMBTU_PER_KG,
        )["usd_per_mmbtu"]

    def get_fci():
        return tea.FCI

    def get_electricity():
        return sum(u.power_utility.rate for u in system.units)

    def get_production_rate():
        return biomethane.F_mass

    metrics = [
        Metric('MSP', get_msp_usd_per_kg, 'USD/kg'),
        Metric('MSP (energy basis)', get_msp_usd_per_mmbtu, 'USD/mmbtu'),
        Metric('Fixed capital investment', get_fci, 'USD'),
        Metric('Total electricity', get_electricity, 'kW'),
        Metric('Biomethane production rate', get_production_rate, 'kg/hr'),
    ]

    model = bst.evaluation.Model(system, metrics)

    add_press_parameters(model, PR)
    add_pressate_concentrator_parameters(model, PC)
    add_mill_parameters(model, ML)
    add_methanogenic_ad_parameters(model, AD)
    add_h2s_removal_parameters(model, H2SR)
    add_biogas_upgrading_parameters(model, UP)
    add_digestate_screw_press_parameters(model, SP)
    add_tea_parameters(model, tea)

    if 'PX' in flowsheet.unit:
        add_peroxide_pretreatment_parameters(model, flowsheet.unit.PX)
    if 'HT' in flowsheet.unit:
        add_heating_pretreatment_parameters(model, flowsheet.unit.HT)
    if 'EZ' in flowsheet.unit:
        add_enzymatic_pretreatment_parameters(model, flowsheet.unit.EZ)

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
    _add_price('Spent H2S media disposal price', spent_h2s_media, _TEA_PRICE["disposal_solid"])
    _add_price(
        'Methanogenic solid digestate disposal price', methanogenic_solid_digestate,
        _TEA_PRICE["disposal_solid"],
    )
    _add_price('Liquid digestate disposal price', liquid_digestate, _TEA_PRICE["disposal_wastewater"])

    # Biomethane is priced on an energy basis (data/tea.yaml `price.
    # biomethane_mmbtu`, USD/mmbtu); the stream's own .price attribute is
    # USD/kg (see _tea.py's usd_per_mmbtu_to_usd_per_kg), so this
    # parameter samples in the yaml's native USD/mmbtu units and converts
    # in the setter, rather than sampling the already-converted USD/kg
    # value directly.
    biomethane_mmbtu_entry = _TEA_PRICE["biomethane_mmbtu"]
    mmbtu_per_kg = float(biomethane_mmbtu_entry["mmbtu_per_kg"])
    baseline = float(biomethane_mmbtu_entry["baseline"])
    def set_biomethane_price(x, mmbtu_per_kg=mmbtu_per_kg):
        biomethane.price = usd_per_mmbtu_to_usd_per_kg(float(x), mmbtu_per_kg)
    param(
        set_biomethane_price, element=biomethane, name='Biomethane price', units='USD/mmbtu',
        baseline=baseline, distribution=distribution_from_yaml(biomethane_mmbtu_entry),
    )

    return model
