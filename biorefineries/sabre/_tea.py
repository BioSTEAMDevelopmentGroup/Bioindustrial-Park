# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Techno-economic analysis for the SaBRe flowsheets.
"""
import biosteam as bst

from biorefineries.sabre.utils import OPERATING_DAYS_PER_YEAR, load_assumptions

__all__ = (
    'SaBReTEA', 'create_tea', 'solve_product_msp',
    'usd_per_mmbtu_to_usd_per_kg', 'usd_per_kg_to_usd_per_mmbtu',
    'disposal_avoided_credit_usd_per_yr', 'apply_revenue_credit',
)

_TEA_DEFAULTS = load_assumptions("tea.yaml")["tea"]
CH4_MMBTU_PER_KG = load_assumptions("tea.yaml")["price"]["biomethane_mmbtu"]["mmbtu_per_kg"]


def usd_per_mmbtu_to_usd_per_kg(usd_per_mmbtu: float, mmbtu_per_kg: float) -> float:
    """Convert an energy-basis price (USD/mmbtu) to a mass-basis price (USD/kg)."""
    return usd_per_mmbtu * mmbtu_per_kg


def usd_per_kg_to_usd_per_mmbtu(usd_per_kg: float, mmbtu_per_kg: float) -> float:
    """Convert a mass-basis price (USD/kg) to an energy-basis price (USD/mmbtu)."""
    return usd_per_kg / mmbtu_per_kg


class SaBReTEA(bst.TEA):
    """
    Simplified TEA

    Assumptions:
    - Fixed operating cost (FOC) is modeled as a fraction of FCI.
    - Variable operating cost (VOC) is taken from BioSTEAM system-level
      material and utility costs, plus any unit-level add_OPEX.
    """

    def __init__(
        self,
        *args,
        foc_frac_of_fci: float = 0.04,
        **kwargs,
    ):
        self.foc_frac_of_fci = float(foc_frac_of_fci)
        super().__init__(*args, **kwargs)

    def _FOC(self, FCI):
        return self.foc_frac_of_fci * FCI

    def _annual_unit_add_opex(self):
        """
        Sum all unit-level add_OPEX contributions and annualize them.
        These are stored as $/hr on the units.
        """
        hourly = 0.0
        for u in self.system.units:
            add_opex = getattr(u, "add_OPEX", None)
            if not add_opex:
                continue

            if isinstance(add_opex, dict):
                hourly += sum(float(v or 0.0) for v in add_opex.values())
            else:
                hourly += float(add_opex or 0.0)

        return hourly * self.operating_days * 24.0

    @property
    def VOC(self):
        """
        Override BioSTEAM's default VOC property so reagent add_OPEX is included
        in the cashflow / MSP calculation --> for reagent cost and other OPEX
        """
        return (
            self.material_cost
            + self.utility_cost
            + self._annual_unit_add_opex()
        )


def create_tea(
    sys,
    IRR: float = _TEA_DEFAULTS["IRR"],
    duration: tuple[int, int] = tuple(_TEA_DEFAULTS["duration"]),
    depreciation: str = _TEA_DEFAULTS["depreciation"],
    income_tax: float = _TEA_DEFAULTS["income_tax"],
    operating_days: int = OPERATING_DAYS_PER_YEAR,
    construction_schedule: tuple[float, ...] = tuple(_TEA_DEFAULTS["construction_schedule"]),
    WC_over_FCI: float = _TEA_DEFAULTS["WC_over_FCI"],
    finance_interest: float = _TEA_DEFAULTS["finance_interest"],
    finance_years: int = _TEA_DEFAULTS["finance_years"],
    finance_fraction: float = _TEA_DEFAULTS["finance_fraction"],
    startup_months: int = _TEA_DEFAULTS["startup_months"],
    startup_FOCfrac: float = _TEA_DEFAULTS["startup_FOCfrac"],
    startup_VOCfrac: float = _TEA_DEFAULTS["startup_VOCfrac"],
    startup_salesfrac: float = _TEA_DEFAULTS["startup_salesfrac"],
    foc_frac_of_fci: float = _TEA_DEFAULTS["foc_frac_of_fci"],
):
    """
    Create a TEA object for any SaBRe system.
    """
    tea = SaBReTEA(
        system=sys,

        # Economic targets
        IRR=IRR,
        duration=duration,
        depreciation=depreciation,
        income_tax=income_tax,
        operating_days=operating_days,

        # Capital deployment
        construction_schedule=construction_schedule,
        WC_over_FCI=WC_over_FCI,

        # Financing
        finance_interest=finance_interest,
        finance_years=finance_years,
        finance_fraction=finance_fraction,

        # Startup assumptions
        lang_factor=None,
        startup_months=startup_months,
        startup_FOCfrac=startup_FOCfrac,
        startup_VOCfrac=startup_VOCfrac,
        startup_salesfrac=startup_salesfrac,

        # Custom
        foc_frac_of_fci=foc_frac_of_fci,
    )
    return tea


def solve_product_msp(
    tea,
    product_stream,
    energy_content_mmbtu_per_kg: float | None = None,
):
    """
    Solve minimum selling price for a product stream, on a whole-stream
    basis. Used for biomethane by passing energy_content_mmbtu_per_kg
    =CH4_MMBTU_PER_KG (data/tea.yaml `price.biomethane_mmbtu.mmbtu_per_kg`).

    Returns a dict with:
    - usd_per_kg
    - annual_product_kg
    - usd_per_mmbtu        (if energy_content_mmbtu_per_kg is provided)
    - annual_product_mmbtu (if energy_content_mmbtu_per_kg is provided)
    """
    usd_per_kg = tea.solve_price(product_stream)
    annual_hours = tea.operating_days * 24
    annual_product_kg = float(product_stream.F_mass) * annual_hours

    result = {
        "usd_per_kg": usd_per_kg,
        "annual_product_kg": annual_product_kg,
        "usd_per_mmbtu": float("nan"),
        "annual_product_mmbtu": float("nan"),
    }

    if energy_content_mmbtu_per_kg is not None:
        result["usd_per_mmbtu"] = usd_per_kg_to_usd_per_mmbtu(usd_per_kg, energy_content_mmbtu_per_kg)
        result["annual_product_mmbtu"] = annual_product_kg * energy_content_mmbtu_per_kg

    return result


def disposal_avoided_credit_usd_per_yr(
    waste_stream, disposal_price_usd_per_kg: float, tea,
) -> float:
    """
    Annual value (USD/yr) of NOT paying to dispose of waste_stream at
    disposal_price_usd_per_kg (a baseline disposal cost -- see data/tea.yaml
    `price.disposal_solid`/`disposal_wastewater`, stored as a negative
    USD/kg) because it's diverted to a productive downstream use (e.g. as
    AD feedstock) instead of being landfilled/disposed. A tipping-fee
    credit for whichever system takes it on.
    """
    return -disposal_price_usd_per_kg * float(waste_stream.F_mass) * tea.operating_days * 24


def apply_revenue_credit(msp: dict, credit_usd_per_yr: float) -> dict:
    """
    Apply a fixed annual revenue credit (e.g. disposal_avoided_credit_usd_per_yr)
    to an already-solved solve_product_msp() result, returning a new dict.

    Valid because tea.solve_price() is linear in any additional fixed
    annual cashflow: adding credit_usd_per_yr of revenue lowers the solved
    product's own required price by exactly
    credit_usd_per_yr / annual_product_units, holding volumes fixed.
    """
    result = dict(msp)
    result["usd_per_kg"] -= credit_usd_per_yr / msp["annual_product_kg"]
    if msp["annual_product_mmbtu"] == msp["annual_product_mmbtu"]:  # not NaN
        result["usd_per_mmbtu"] -= credit_usd_per_yr / msp["annual_product_mmbtu"]
    return result
