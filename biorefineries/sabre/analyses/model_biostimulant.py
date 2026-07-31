# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
biosteam.evaluation.Model construction for the standalone biostimulant
flowsheet (systems._biostimulant_system.create_biostimulant_system()).

add_press_parameters() and add_pressate_concentrator_parameters() add only
process/cost parameters (no stream prices), so ad_biomethane/ad_vfa's own
(future) model files can import and reuse them against their own embedded
PR/PC instances -- see
docs/superpowers/specs/2026-07-31-sabre-uncertainty-model-design.md.
"""
import biosteam as bst
from biosteam.evaluation import Metric

from biorefineries.sabre.utils import load_assumptions
from biorefineries.sabre.systems._biostimulant_system import create_biostimulant_system
from biorefineries.sabre._tea import solve_product_msp
from biorefineries.sabre.analyses.model_utils import distribution_from_yaml, add_tea_parameters

__all__ = (
    'add_press_parameters', 'add_pressate_concentrator_parameters',
    'create_biostimulant_model',
)

_PRESS = load_assumptions("preprocessing.yaml")["press"]
_CONCENTRATOR = load_assumptions("biostimulant.yaml")["pressate_concentrator"]
_TEA_PRICE = load_assumptions("tea.yaml")["price"]


def add_press_parameters(model, PR):
    """
    Add Press (PR) process/cost parameters to `model` (data/preprocessing.yaml
    `press`): solids_capture_frac, cake_solids_wt_frac,
    solubles_to_pressate_frac, power_kWh_per_dry_ton_TS. Excludes
    capex_model/basis (categorical) and the ref_capacity/capex_installed_ref/
    scale_exponent CAPEX anchor set, and F_BM. Reusable by any sabre
    flowsheet that embeds a Press.
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_PRESS[attr])
        def setter(x, attr=attr):
            setattr(PR, attr, float(x))
        param(
            setter, element=PR, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('Press solids capture fraction', 'solids_capture_frac')
    _add('Press cake solids wt fraction', 'cake_solids_wt_frac')
    _add('Press solubles to pressate fraction', 'solubles_to_pressate_frac')
    _add('Press power intensity', 'power_kWh_per_dry_ton_TS', units='kWh/dry ton')

    return model


def add_pressate_concentrator_parameters(model, PC):
    """
    Add PressateConcentrator (PC) process/cost parameters to `model`
    (data/biostimulant.yaml `pressate_concentrator`): water_recovery_to_permeate,
    retained_solute_recovery_to_concentrate, nontarget_solute_recovery_to_permeate,
    design_flux_L_m2_h, electricity_kWh_per_m3_feed, capex_usd_per_m2,
    maintenance_frac_of_capex_per_yr. Excludes enabled/retained_solute_IDs
    (categorical), F_BM, and operating_pressure_bar (permanently inert --
    the unit's own docstring says it's design-basis only, never used in a
    cost/duty calculation). Reusable by any sabre flowsheet that embeds a
    PressateConcentrator.
    """
    param = model.parameter

    def _add(name, attr, units=''):
        baseline = float(_CONCENTRATOR[attr])
        def setter(x, attr=attr):
            setattr(PC, attr, float(x))
        param(
            setter, element=PC, name=name, units=units,
            baseline=baseline, distribution=distribution_from_yaml(baseline),
        )

    _add('PC water recovery to permeate', 'water_recovery_to_permeate')
    _add('PC retained solute recovery to concentrate', 'retained_solute_recovery_to_concentrate')
    _add('PC nontarget solute recovery to permeate', 'nontarget_solute_recovery_to_permeate')
    _add('PC design flux', 'design_flux_L_m2_h', units='L/m2-h')
    _add('PC electricity intensity', 'electricity_kWh_per_m3_feed', units='kWh/m3')
    _add('PC membrane capex per m2', 'capex_usd_per_m2', units='USD/m2')
    _add('PC maintenance fraction of capex', 'maintenance_frac_of_capex_per_yr')

    return model


def create_biostimulant_model(system=None):
    """
    Build the full biosteam.evaluation.Model for the standalone biostimulant
    flowsheet: Press + PressateConcentrator process/cost parameters, the
    shared TEA parameters, and this system's own boundary-stream prices
    (sargassum feed, pressed_cake disposal, permeate disposal).

    Parameters
    ----------
    system : bst.System, optional
        Defaults to a fresh `create_biostimulant_system()`.
    """
    if system is None:
        system = create_biostimulant_system()
    system.simulate()

    flowsheet = system.flowsheet
    PR = flowsheet.unit.PR
    PC = flowsheet.unit.PC
    feed = flowsheet.stream.sargassum_feed
    pressed_cake = flowsheet.stream.pressed_cake
    permeate = flowsheet.stream.permeate
    biostimulant_product = flowsheet.stream.biostimulant_product
    tea = system.TEA

    def get_msp():
        return solve_product_msp(tea=tea, product_stream=biostimulant_product)["usd_per_kg"]

    def get_fci():
        return tea.FCI

    def get_electricity():
        return sum(u.power_utility.rate for u in system.units)

    def get_production_rate():
        return biostimulant_product.F_mass

    metrics = [
        Metric('MSP', get_msp, 'USD/kg'),
        Metric('Fixed capital investment', get_fci, 'USD'),
        Metric('Total electricity', get_electricity, 'kW'),
        Metric('Biostimulant production rate', get_production_rate, 'kg/hr'),
    ]

    model = bst.evaluation.Model(system, metrics)

    add_press_parameters(model, PR)
    add_pressate_concentrator_parameters(model, PC)
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
    _add_price('Pressed cake disposal price', pressed_cake, _TEA_PRICE["disposal_solid"])
    _add_price('Permeate disposal price', permeate, _TEA_PRICE["disposal_wastewater"])

    return model
