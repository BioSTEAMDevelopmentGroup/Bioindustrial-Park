# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Cradle-to-gate life cycle assessment (GWP100) for the SaBRe flowsheets,
built on BioSTEAM's own impact accounting (stream and utility
characterization factors, `System.get_net_impact`).

Scope and conventions
---------------------
- Indicator: 100-year global warming potential, key 'GWP', kg CO2-eq.
- Functional unit: 1 kg of the product stream (for the 3-HP system, dried
  Ca(3HP)2), also reported per kg 3-HP.
- Boundary: feedstock, purchased chemicals and enzyme, natural gas, and
  grid electricity (consumed). Exported electricity displaces the grid
  (reported with and without this credit).
- Biogenic CO2 (from burning Sargassum solids and cell mass) is counted as
  zero; carbon stored in the product, and the product's end of life, are not
  counted. Fossil CO2 from natural gas combustion is counted.
- Not included: feedstock harvest and transport (the Sargassum burden is
  set in data/lca.yaml, default 0), capital goods and construction, and
  treatment and disposal of wastewater, mother-liquor purge and ash. The
  last of these has no treatment unit in the model, so those burdens are
  missing, not zero.
- Heating and cooling utilities carry no factor: steam comes from the
  boiler, whose fuel is already counted, so it would otherwise be double
  counted.

Characterization factors are placeholders reused from the lactic
biorefinery; see data/lca.yaml.
"""
import biosteam as bst

from biorefineries.sabre.utils import load_assumptions

__all__ = ('GWP_KEY', 'set_characterization_factors', 'get_gwp')

GWP_KEY = 'GWP'
_GWP = load_assumptions("lca.yaml")["gwp100"]

# Feed streams are matched by chemical content, since most feeds are unnamed.
# Chemical ID -> GWP (kg CO2e / kg chemical). Methane's combustion CO2 is
# added in `_chemical_cfs`.
_CHEMICAL_KEYS = {
    "Enzyme": "enzyme",
    "CalciumDihydroxide": "lime",
    "H2SO4": "sulfuric_acid",
    "Ammonia": "ammonia",
}
# Chemicals that carry no burden (process water, air) or whose burden is
# negligible and not tracked (boiler ash additive).
_IGNORED = {"Water", "Oxygen", "Nitrogen", "Ash"}


def _chemical_cfs(chemicals) -> dict:
    cfs = {ID: _GWP[key] for ID, key in _CHEMICAL_KEYS.items() if ID in chemicals.IDs}
    if "Methane" in chemicals.IDs:
        CO2_per_CH4 = chemicals.CO2.MW / chemicals.Methane.MW  # complete combustion
        cfs["Methane"] = _GWP["methane_upstream"] + CO2_per_CH4
    return cfs


def set_characterization_factors(
        sys, feedstock_ID: str = "sargassum_feed",
        feedstock_gwp_per_kg_wet: float | None = None) -> list:
    """
    Set GWP characterization factors on the system's feed streams and on
    electricity. Returns a list of (stream ID, composition) for feeds that
    were left without a factor.
    """
    bst.settings.define_impact_indicator(GWP_KEY, "kg*CO2e")
    kWh_cf = _GWP["electricity_kg_CO2e_per_kWh"]
    bst.PowerUtility.characterization_factors[GWP_KEY] = (kWh_cf, kWh_cf)
    if feedstock_gwp_per_kg_wet is None:
        feedstock_gwp_per_kg_wet = _GWP["sargassum_kg_CO2e_per_kg_wet"]

    unassigned = []
    for s in sys.feeds:
        if s.F_mass <= 0:
            continue
        if s.ID == feedstock_ID:
            s.characterization_factors[GWP_KEY] = feedstock_gwp_per_kg_wet
            continue
        cfs = _chemical_cfs(s.chemicals)
        total = float(s.F_mass)
        s.characterization_factors[GWP_KEY] = sum(
            float(s.imass[ID]) / total * cf for ID, cf in cfs.items()
        )
        leftover = {
            ID: float(s.imass[ID]) for ID in s.chemicals.IDs
            if ID not in cfs and ID not in _IGNORED and s.imass[ID] > 1e-6
        }
        if leftover:
            unassigned.append((s.ID or "<unnamed>", leftover))
    return unassigned


def get_gwp(sys, product, feedstock_ID: str = "sargassum_feed") -> dict:
    """
    GWP of `product` (a stream of `sys`), after
    `set_characterization_factors(sys)` and a simulation. Returns kg CO2e per
    kg of product, with the exported-electricity credit ("gwp_kg_per_kg") and
    without it ("gwp_no_credit_kg_per_kg"), per kg 3-HP (if the product is
    Ca(3HP)2), and the contributions per kg of product.
    """
    hours = sys.operating_hours
    annual_kg = float(product.F_mass) * hours
    key = GWP_KEY

    net = sys.get_net_impact(key)
    net_no_credit = sys.get_net_impact(key, displace=False)

    contrib = {"feedstock": 0.0, "enzyme": 0.0, "lime": 0.0, "natural_gas": 0.0, "other_feeds": 0.0}
    for s in sys.feeds:
        if s.F_mass <= 0 or key not in s.characterization_factors:
            continue
        impact = s.get_impact(key) * hours
        if s.ID == feedstock_ID:
            contrib["feedstock"] += impact
        elif "Enzyme" in s.chemicals.IDs and s.imass["Enzyme"] > 0.5 * s.F_mass:
            contrib["enzyme"] += impact
        elif "Methane" in s.chemicals.IDs and s.imass["Methane"] > 0.5 * s.F_mass:
            contrib["natural_gas"] += impact
        elif "CalciumDihydroxide" in s.chemicals.IDs and s.imass["CalciumDihydroxide"] > 0:
            contrib["lime"] += impact
        else:
            contrib["other_feeds"] += impact
    contrib["electricity_net"] = sys.get_net_electricity_impact(key)  # negative = credit

    result = {
        "gwp_kg_per_kg": net / annual_kg,
        "gwp_no_credit_kg_per_kg": net_no_credit / annual_kg,
        **{f"gwp_{k}_kg_per_kg": v / annual_kg for k, v in contrib.items()},
    }
    chem = product.chemicals
    if "Ca3HP2" in chem.IDs and "HP" in chem.IDs:
        to_3hp = chem.Ca3HP2.MW / (2 * chem.HP.MW)  # kg Ca(3HP)2 per kg 3-HP
        result["gwp_per_kg_3hp"] = result["gwp_kg_per_kg"] * to_3hp
        result["gwp_no_credit_per_kg_3hp"] = result["gwp_no_credit_kg_per_kg"] * to_3hp
    return result
