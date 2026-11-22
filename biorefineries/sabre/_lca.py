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
- Boundary: feedstock, purchased chemicals and enzyme, process and makeup
  water, facility chemicals (boiler, cooling tower), natural gas (upstream
  supply plus the fossil CO2 from burning it on site, in the dryer and, if
  used, the boiler), and grid electricity (consumed). Exported electricity
  displaces the grid (reported with and without this credit). Liquid waste
  leaving the plant and boiler ash carry a treatment/disposal burden.
- Biogenic CO2 (from burning Sargassum solids and cell mass) is counted as
  zero; carbon stored in the product, and the product's end of life, are not
  counted. Fossil CO2 from natural gas combustion is counted, on the natural
  gas inlet (upstream factor plus the stoichiometric CO2 of complete
  combustion), so the combustion gas outlets carry no factor of their own.
- Not included: feedstock harvest and transport (the Sargassum burden is
  set in data/lca.yaml, default 0) and capital goods and construction.
  Wastewater and ash burdens are per-kg placeholders (data/lca.yaml), not a
  treatment-plant model.
- Heating and cooling utilities carry no factor: steam comes from the
  boiler, whose fuel is already counted, so it would otherwise be double
  counted.
- Water is charged once, on the process water center (PWC) inlet; the boiler
  and cooling tower makeup streams, which the PWC supplies, carry no factor.

Characterization factors are placeholders reused from other biorefineries in
this repo; see data/lca.yaml.
"""
import biosteam as bst

from biorefineries.sabre.utils import load_assumptions

__all__ = ('GWP_KEY', 'set_characterization_factors', 'get_gwp', 'inventory_check')

GWP_KEY = 'GWP'
_GWP = load_assumptions("lca.yaml")["gwp100"]

# Feed streams are matched by chemical content, since most feeds are unnamed.
# Chemical ID -> GWP (kg CO2e / kg chemical). Methane's combustion CO2 is
# added in `_chemical_cfs`.
_CHEMICAL_KEYS = {
    "Enzyme": "enzyme",
    "CalciumDihydroxide": "lime",
    "Ammonia": "ammonia",
    "H2SO4": "sulfuric_acid",
    "NaCl": "sodium_chloride",
}
# Chemicals that carry no burden (process water, air) or whose burden is
# negligible and not tracked.
_IGNORED = {"Water", "Oxygen", "Nitrogen", "Ash"}


def _chemical_cfs(chemicals) -> dict:
    cfs = {ID: _GWP[key] for ID, key in _CHEMICAL_KEYS.items() if ID in chemicals.IDs}
    if "Methane" in chemicals.IDs:
        CO2_per_CH4 = chemicals.CO2.MW / chemicals.Methane.MW  # complete combustion
        cfs["Methane"] = _GWP["methane_upstream"] + CO2_per_CH4
    return cfs


def _units(sys) -> dict:
    return {u.ID: u for u in (*sys.units, *sys.facilities)}


def _mass(s, *IDs) -> float:
    """Mass flow (kg/h) of the listed chemicals, ignoring those not in the stream's set."""
    present = [i for i in IDs if i in s.chemicals.IDs]
    return float(s.imass[present].sum()) if present else 0.0


def _is_water_only(s) -> bool:
    return s.F_mass > 0 and _mass(s, "Water") > 0.999 * s.F_mass


def set_characterization_factors(
        sys, feedstock_ID: str = "sargassum_feed",
        feedstock_gwp_per_kg_wet: float | None = None) -> list:
    """
    Set GWP characterization factors on every feed stream of `sys` and on
    electricity. Streams are matched by role (boiler fuel, boiler chemicals,
    cooling tower chemicals, makeup water, dryer fuel) where the content is not
    enough to tell, and otherwise by chemical content. Zero-flow fuel inlets get
    their factor too, so a scenario that later burns natural gas in the boiler is
    counted. Returns a list of (stream ID, composition) for feeds left without a
    factor.
    """
    bst.settings.define_impact_indicator(GWP_KEY, "kg*CO2e")
    kWh_cf = _GWP["electricity_kg_CO2e_per_kWh"]
    bst.PowerUtility.characterization_factors[GWP_KEY] = (kWh_cf, kWh_cf)
    if feedstock_gwp_per_kg_wet is None:
        feedstock_gwp_per_kg_wet = _GWP["sargassum_kg_CO2e_per_kg_wet"]

    CO2_per_CH4 = sys.feeds[0].chemicals.CO2.MW / sys.feeds[0].chemicals.Methane.MW
    CH4_CF = _GWP["methane_upstream"] + CO2_per_CH4
    units = _units(sys)
    by_role = {}
    BT, CT = units.get("BT"), units.get("CT")
    if BT is not None:
        by_role[id(BT.fuel)] = CH4_CF
        by_role[id(BT.ins[5])] = _GWP["boiler_chemicals"]
        by_role[id(BT.makeup_water)] = 0.0  # supplied through the PWC inlet
    if CT is not None:
        by_role[id(CT.ins[0])] = 0.0  # supplied through the PWC inlet
        by_role[id(CT.ins[1])] = _GWP["cooling_tower_chemicals"]
    for u in sys.units:
        if hasattr(u, "natural_gas"):
            by_role[id(u.natural_gas)] = CH4_CF

    unassigned = []
    for s in sys.feeds:
        if id(s) in by_role:
            s.characterization_factors[GWP_KEY] = by_role[id(s)]
            continue
        if s.F_mass <= 0:
            continue
        if s.ID == feedstock_ID:
            s.characterization_factors[GWP_KEY] = feedstock_gwp_per_kg_wet
            continue
        if _is_water_only(s):
            s.characterization_factors[GWP_KEY] = _GWP["process_water"]
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


def _waste_streams(sys):
    """(liquid waste streams, boiler ash stream) leaving `sys`."""
    liquid = [
        s for s in sys.products
        if s.F_mass > 0 and s.phase == "l" and (s.price or 0.0) < 0.0
    ]
    BT = _units(sys).get("BT")
    ash = BT.outs[2] if BT is not None else None
    return liquid, ash


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

    liquid_waste, ash = _waste_streams(sys)
    wastewater = sum(float(s.F_mass) for s in liquid_waste) * hours * _GWP["wastewater_treatment"]
    ash_impact = (float(ash.F_mass) * hours * _GWP["ash_disposal"]) if ash is not None else 0.0

    net = sys.get_net_impact(key) + wastewater + ash_impact
    net_no_credit = sys.get_net_impact(key, displace=False) + wastewater + ash_impact

    contrib = {
        "feedstock": 0.0, "enzyme": 0.0, "lime": 0.0, "natural_gas": 0.0,
        "water": 0.0, "ammonia": 0.0, "sulfuric_acid": 0.0, "sodium_chloride": 0.0,
        "facility_chemicals": 0.0, "other_feeds": 0.0,
    }
    for s in sys.feeds:
        if key not in s.characterization_factors:
            continue
        impact = s.get_impact(key) * hours
        if impact == 0.0:
            continue
        if s.ID == feedstock_ID:
            contrib["feedstock"] += impact
        elif "Enzyme" in s.chemicals.IDs and s.imass["Enzyme"] > 0.5 * s.F_mass:
            contrib["enzyme"] += impact
        elif "Methane" in s.chemicals.IDs and s.imass["Methane"] > 0.5 * s.F_mass:
            contrib["natural_gas"] += impact
        elif "CalciumDihydroxide" in s.chemicals.IDs and s.imass["CalciumDihydroxide"] > 0:
            contrib["lime"] += impact
        elif "Ammonia" in s.chemicals.IDs and s.imass["Ammonia"] > 0.5 * s.F_mass:
            contrib["ammonia"] += impact
        elif "H2SO4" in s.chemicals.IDs and s.imass["H2SO4"] > 0:
            contrib["sulfuric_acid"] += impact
        elif "NaCl" in s.chemicals.IDs and s.imass["NaCl"] > 0.5 * s.F_mass:
            contrib["sodium_chloride"] += impact
        elif _is_water_only(s) and s.characterization_factors[key] == _GWP["process_water"]:
            contrib["water"] += impact
        else:
            contrib["facility_chemicals"] += impact
    contrib["wastewater_treatment"] = wastewater
    contrib["ash_disposal"] = ash_impact
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


def inventory_check(sys, product=None) -> list:
    """
    One row per feed and product stream with flow, price, GWP factor and how its
    cost and impact are accounted, so it is clear that nothing is left out.
    Run after `set_characterization_factors(sys)` and a simulation. Rows with
    status 'MISSING' are streams that carry flow but have neither a cost nor
    an impact and no stated reason.

    Costs that are not on the stream itself are booked by the TEA through
    utilities (natural gas, makeup water through the process water center,
    boiler ash disposal), and so are listed as such.
    """
    units = _units(sys)
    BT, CT, PWC = units.get("BT"), units.get("CT"), units.get("PWC")
    utility_cost_streams = set()
    if PWC is not None:
        utility_cost_streams.update(id(s) for s in PWC.ins)
    if BT is not None:
        utility_cost_streams.update((id(BT.makeup_water), id(BT.fuel)))
    if CT is not None:
        utility_cost_streams.add(id(CT.ins[0]))
    # Steam fed to a SteamMixer is booked as a heat utility (steam from the
    # boiler), so its cost and impact come through the boiler's fuel.
    steam_streams = set()
    for u in sys.units:
        if hasattr(u, "natural_gas"):
            utility_cost_streams.add(id(u.natural_gas))
        if isinstance(u, bst.SteamMixer):
            steam_streams.add(id(u.ins[1]))
    liquid_waste, ash = _waste_streams(sys)
    liquid_ids = {id(s) for s in liquid_waste}
    rows = []

    def add(kind, s, cost_basis, impact_basis, ok):
        origin = (s.source or s.sink)
        rows.append(dict(
            kind=kind, stream=s.ID or f"<{origin.ID}>",
            kg_per_h=float(s.F_mass), price_usd_per_kg=float(s.price or 0.0),
            gwp_kg_per_kg=float(s.characterization_factors.get(GWP_KEY, 0.0) or 0.0),
            cost_basis=cost_basis, impact_basis=impact_basis,
            status="ok" if ok else "MISSING",
        ))

    for s in sys.feeds:
        if s.F_mass <= 0:
            continue
        has_cf = GWP_KEY in s.characterization_factors
        air = _mass(s, "Oxygen", "Nitrogen") > 0.999 * s.F_mass
        if s.ID == "sargassum_feed":
            cost, imp = "assumed zero (data/tea.yaml)", "assumed zero (data/lca.yaml)"
        elif air:
            cost, imp = "free (air)", "zero (air)"
        elif id(s) in steam_streams:
            cost, imp = "heat utility (steam)", "boiler fuel (counted at BT)"
        else:
            cost = ("utility price" if id(s) in utility_cost_streams
                    else "stream price" if s.price else "")
            imp = "characterization factor" if has_cf else ""
        add("feed", s, cost, imp, bool(cost) and bool(imp))
    for s in sys.products:
        if s.F_mass <= 0:
            continue
        if product is not None and s is product:
            add("product", s, "sold (MSP solved)", "functional unit", True)
        elif s.source is not None and s.source.ID == "PS":
            add("product", s, "ships with product, unpriced", "ships with product", True)
        elif s.source is not None and s.source.ID == "PWC":
            add("internal", s, "water charged on PWC inlet", "counted on PWC inlet", True)
        elif s.phase == "g":
            add("gas", s, "none (gas outlet)", "none (fossil CO2 counted on the fuel inlet)", True)
        elif id(s) in liquid_ids:
            add("waste", s, "disposal price", "wastewater treatment factor", True)
        elif s is ash:
            add("waste", s, "utility price (ash disposal)", "ash disposal factor", True)
        else:
            add("product", s, "", "", False)
    return rows
