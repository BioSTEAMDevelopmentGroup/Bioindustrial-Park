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
)

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
    IRR: float = 0.10,
    duration: tuple[int, int] = (2026, 2046),
    depreciation: str = "MACRS7",
    income_tax: float = 0.21,
    operating_days: int = OPERATING_DAYS_PER_YEAR,
    construction_schedule: tuple[float, ...] = (0.4, 0.6),
    WC_over_FCI: float = 0.05,
    finance_interest: float = 0.08,
    finance_years: int = 10,
    finance_fraction: float = 0.6,
    startup_months: int = 3,
    startup_FOCfrac: float = 1.0,
    startup_VOCfrac: float = 0.5,
    startup_salesfrac: float = 0.5,
    foc_frac_of_fci: float = 0.04,
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
    - usd_per_kg_stream
    - annual_stream_kg
    - usd_per_mmbtu        (if energy basis is provided)
    - annual_product_mmbtu (if energy basis is provided)
    """
    msp_usd_per_kg_stream = tea.solve_price(product_stream)
    annual_hours = tea.operating_days * 24
    annual_stream_kg = float(product_stream.F_mass) * annual_hours

    result = {
        "usd_per_kg_stream": msp_usd_per_kg_stream,
        "annual_stream_kg": annual_stream_kg,
        "usd_per_mmbtu": float("nan"),
        "annual_product_mmbtu": float("nan"),
    }

    if energy_content_mmbtu_per_kg is not None and energy_content_mmbtu_per_kg > 0:
        result["usd_per_mmbtu"] = usd_per_kg_to_usd_per_mmbtu(msp_usd_per_kg_stream, energy_content_mmbtu_per_kg)
        result["annual_product_mmbtu"] = annual_stream_kg * energy_content_mmbtu_per_kg

    return result
