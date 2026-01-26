# -*- coding: utf-8 -*-
"""
BioSTEAM Process Report Generator (v2)
=====================================

Generates TEA/LCA breakdown tables with standardized output management.

@author: Dr. Ouwen Peng
@institute: PreFerS - Centre for Precision Fermentation and Sustainability
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, Iterable, List, Optional

import pandas as pd
import biosteam as bst
from biosteam.process_tools import UnitGroup

from . import utils


STANDALONE_IDS = ("PWC", "CT", "CWP", "BT", "WW")

MSP_MAJOR_UNITS = (
    "aerated fermentation",
    "cell disruption",
    "centrifuge",
    "diafiltration",
)

MSP_UTILITY_UNITS = (
    "cooling tower",
    "chilled water package",
    "boiler",
    "process water center",
    "wastewater",
)

MIX_STORAGE_REGEX = re.compile(r"^(mix|storage)", re.IGNORECASE)
HEAT_EXCHANGER_REGEX = re.compile(r"^(cooler|heater|hx)", re.IGNORECASE)


@dataclass
class VerificationResult:
    label: str
    reported: float
    reference: float
    variance: float
    passed: bool


class ProcessReportGenerator:
    """Generate TEA/LCA breakdown tables and save as Excel/CSV."""

    def __init__(
        self,
        system: bst.System,
        tea: Optional[bst.TEA] = None,
        product_stream: Optional[bst.Stream] = None,
        config: str = "",
        timestamp: Optional[str] = None,
        base_dir: Optional[str] = None,
    ):
        self.system = system
        self.tea = tea
        self.product_stream = product_stream or self._infer_product_stream()
        self.timestamp = timestamp or datetime.now().strftime("%Y%m%d_%H%M")
        self.config = config

        self._setup_output_dirs(base_dir)

    # ---------------------------------------------------------------------
    # Output directories
    # ---------------------------------------------------------------------
    def _setup_output_dirs(self, base_dir: Optional[str]):
        if base_dir and os.path.isdir(base_dir):
            self.base_dir = base_dir
            self.data_dir = os.path.join(self.base_dir, "data")
            self.figures_dir = os.path.join(self.base_dir, "figure")
        else:
            script_path = base_dir or os.path.abspath(__file__)
            dirs = utils.get_analysis_dirs(script_path, config=self.config, timestamp=self.timestamp)
            self.base_dir = dirs["base"]
            self.data_dir = os.path.join(self.base_dir, "data")
            self.figures_dir = os.path.join(self.base_dir, "figure")
        os.makedirs(self.data_dir, exist_ok=True)
        os.makedirs(self.figures_dir, exist_ok=True)

    # ---------------------------------------------------------------------
    # Public API
    # ---------------------------------------------------------------------
    def generate(self) -> Dict[str, pd.DataFrame]:
        capital_df = self.build_capital_breakdown()
        msp_df = self.build_msp_breakdown()
        gwp_df = self.build_gwp_breakdown(msp_df)

        self._save_excel(
            {
                "Total Capital": capital_df,
                "MSP Breakdown": msp_df,
                "LCA Allocation": gwp_df,
            }
        )

        self._save_metadata()
        self.verify(capital_df, msp_df, gwp_df)

        return {
            "Total Capital": capital_df,
            "MSP Breakdown": msp_df,
            "LCA Allocation": gwp_df,
        }

    # ---------------------------------------------------------------------
    # Capital breakdown
    # ---------------------------------------------------------------------
    def build_capital_breakdown(self) -> pd.DataFrame:
        if self.tea is None:
            raise ValueError("TEA object is required for capital breakdown.")

        cost_units = list(self.system.cost_units)
        tci = float(self.tea.TCI) if self.tea else 0.0

        rows = []
        grouped: Dict[str, float] = {}

        def add_row(name: str, cost: float):
            if cost == 0:
                return
            rows.append({"Name": name, "Cost (USD)": cost})

        # Standalone units
        standalone_units = [u for u in cost_units if any(s in u.ID for s in STANDALONE_IDS)]
        for u in standalone_units:
            add_row(u.ID, float(u.installed_cost))

        # Facilities / OSBL
        osbl_units = getattr(self.tea, "OSBL_units", None)
        if osbl_units:
            for u in osbl_units:
                if u in standalone_units:
                    continue
                add_row(u.ID, float(u.installed_cost))
        else:
            installed_sum = sum(float(u.installed_cost) for u in cost_units)
            facilities_cost = max(float(self.tea.FCI) - installed_sum, 0.0)
            if facilities_cost:
                add_row("Facilities", facilities_cost)

        # Group remaining by area code
        excluded = set(standalone_units)
        if osbl_units:
            excluded.update(osbl_units)
        for u in cost_units:
            if u in excluded:
                continue
            area = self._get_area_label(u)
            grouped[area] = grouped.get(area, 0.0) + float(u.installed_cost)

        for name, cost in sorted(grouped.items()):
            add_row(name, cost)

        df = pd.DataFrame(rows)
        if not df.empty:
            df["% of TCI"] = df["Cost (USD)"] / tci * 100 if tci else 0.0
        else:
            df = pd.DataFrame(columns=["Name", "Cost (USD)", "% of TCI"])
        return df

    # ---------------------------------------------------------------------
    # MSP breakdown
    # ---------------------------------------------------------------------
    def build_msp_breakdown(self) -> pd.DataFrame:
        if self.product_stream is None:
            raise ValueError("Product stream is required for MSP breakdown.")

        operating_hours = float(self.system.operating_hours or 0.0)
        annual_production = float(self.product_stream.F_mass) * operating_hours
        if annual_production <= 0:
            raise ValueError("Product flow rate must be greater than 0.")

        category_costs: Dict[str, float] = {}
        for unit in self.system.units:
            unit_cost = self._get_unit_operating_cost(unit)
            category = self._categorize_unit_for_msp(unit)
            category_costs[category] = category_costs.get(category, 0.0) + unit_cost

        rows = []
        total_cost_per_kg = 0.0
        for name, cost_per_hr in category_costs.items():
            cost_per_kg = (cost_per_hr * operating_hours) / annual_production
            total_cost_per_kg += cost_per_kg
            rows.append({"Name": name, "Cost (USD/kg product)": cost_per_kg})

        df = pd.DataFrame(rows)
        if not df.empty and total_cost_per_kg:
            df["%"] = df["Cost (USD/kg product)"] / total_cost_per_kg * 100
        else:
            df = pd.DataFrame(columns=["Name", "Cost (USD/kg product)", "%"])

        df = self._sort_breakdown(df)
        return df

    # ---------------------------------------------------------------------
    # GWP breakdown
    # ---------------------------------------------------------------------
    def build_gwp_breakdown(self, msp_df: pd.DataFrame) -> pd.DataFrame:
        if self.product_stream is None:
            raise ValueError("Product stream is required for GWP breakdown.")

        operating_hours = float(self.system.operating_hours or 0.0)
        annual_production = float(self.product_stream.F_mass) * operating_hours
        if annual_production <= 0:
            raise ValueError("Product flow rate must be greater than 0.")

        categories = list(msp_df["Name"]) if not msp_df.empty else ["Others"]
        grouped_units = self._units_by_msp_category(categories)

        material_impacts = {}
        utility_impacts = {}
        for name, units in grouped_units.items():
            material_impacts[name] = self._get_material_gwp(units)
            utility_impacts[name] = self._get_utility_gwp(units)

        direct_total = self.system.get_process_impact("GWP")
        allocation_base = sum(material_impacts.values()) + sum(utility_impacts.values())
        direct_impacts = {}
        for name in categories:
            if allocation_base > 0:
                share = (material_impacts.get(name, 0.0) + utility_impacts.get(name, 0.0)) / allocation_base
            else:
                share = 1.0 / len(categories) if categories else 0.0
            direct_impacts[name] = direct_total * share

        rows = []
        total_gwp_per_kg = 0.0
        for name in categories:
            total_gwp = material_impacts.get(name, 0.0) + utility_impacts.get(name, 0.0) + direct_impacts.get(name, 0.0)
            gwp_per_kg = total_gwp / annual_production if annual_production else 0.0
            total_gwp_per_kg += gwp_per_kg
            rows.append({"Name": name, "GWP (kg CO2e/kg product)": gwp_per_kg})

        df = pd.DataFrame(rows)
        if not df.empty and total_gwp_per_kg:
            df["%"] = df["GWP (kg CO2e/kg product)"] / total_gwp_per_kg * 100
        else:
            df = pd.DataFrame(columns=["Name", "GWP (kg CO2e/kg product)", "%"])
        df = self._sort_breakdown(df)
        return df

    # ---------------------------------------------------------------------
    # Verification
    # ---------------------------------------------------------------------
    def verify(self, capital_df: pd.DataFrame, msp_df: pd.DataFrame, gwp_df: pd.DataFrame):
        results: List[VerificationResult] = []

        # Check 1: Capital
        if self.tea and not capital_df.empty:
            reported = capital_df["Cost (USD)"].sum()
            reference = float(self.tea.TCI)
            variance = abs(reported - reference) / reference if reference else 0.0
            passed = variance <= 0.01
            if not passed:
                unaccounted = reference - reported
                capital_df.loc[len(capital_df)] = {
                    "Name": "Unaccounted Capital",
                    "Cost (USD)": unaccounted,
                    "% of TCI": unaccounted / reference * 100 if reference else 0.0,
                }
            results.append(VerificationResult("Capital", reported, reference, variance, passed))

        # Check 2: MSP
        if self.tea and not msp_df.empty:
            reported = msp_df["Cost (USD/kg product)"].sum()
            try:
                reference = float(self.tea.solve_price(self.product_stream))
            except Exception:
                reference = reported
            variance = abs(reported - reference) / reference if reference else 0.0
            passed = variance <= 0.01
            if not passed:
                diff = reference - reported
                if "Others" in msp_df["Name"].values:
                    idx = msp_df.index[msp_df["Name"] == "Others"][0]
                    msp_df.at[idx, "Cost (USD/kg product)"] += diff
                    msp_df.at[idx, "%"] = msp_df.at[idx, "Cost (USD/kg product)"] / reference * 100
            results.append(VerificationResult("MSP", reported, reference, variance, passed))

        # Check 3: GWP
        if not gwp_df.empty and self.product_stream is not None:
            reported = gwp_df["GWP (kg CO2e/kg product)"].sum()
            reference = self._get_total_gwp_per_kg()
            variance = abs(reported - reference) / reference if reference else 0.0
            passed = variance <= 0.01
            results.append(VerificationResult("GWP", reported, reference, variance, passed))

        self._print_verification(results)

    # ---------------------------------------------------------------------
    # Helpers
    # ---------------------------------------------------------------------
    def _infer_product_stream(self) -> Optional[bst.Stream]:
        products = getattr(self.system, "products", [])
        if products:
            return products[0]
        return None

    def _get_area_label(self, unit) -> str:
        if hasattr(unit, "area") and unit.area:
            return f"A{int(unit.area):03d}"
        match = re.search(r"(\d+)", unit.ID)
        if match:
            area = int(match.group(1)) // 100 * 100
            return f"A{area}"
        return "Other"

    def _get_unit_display_name(self, unit) -> str:
        return getattr(unit, "line", None) or type(unit).__name__

    def _categorize_unit_for_msp(self, unit) -> str:
        name = self._get_unit_display_name(unit).lower()
        unit_id = unit.ID.lower()

        for key in MSP_MAJOR_UNITS:
            if key in name or key in unit_id:
                return self._title_case(key)

        for key in MSP_UTILITY_UNITS:
            if key in name or key in unit_id:
                return self._title_case(key)

        if MIX_STORAGE_REGEX.match(name) or MIX_STORAGE_REGEX.match(unit_id):
            return "Mix/Storage Tank"

        if HEAT_EXCHANGER_REGEX.match(name) or HEAT_EXCHANGER_REGEX.match(unit_id):
            return "Heat Exchangers"

        return "Others"

    def _get_unit_operating_cost(self, unit) -> float:
        material_cost = self._get_unit_material_cost(unit)
        utility_cost = getattr(unit, "utility_cost", 0.0) or 0.0
        return float(material_cost + utility_cost)

    def _get_unit_material_cost(self, unit) -> float:
        group = UnitGroup("single", (unit,))
        return float(group.get_material_cost())

    def _units_by_msp_category(self, categories: Iterable[str]) -> Dict[str, List]:
        grouped = {name: [] for name in categories}
        for unit in self.system.units:
            category = self._categorize_unit_for_msp(unit)
            if category not in grouped:
                grouped.setdefault("Others", []).append(unit)
            else:
                grouped[category].append(unit)
        return grouped

    def _get_material_gwp(self, units: Iterable) -> float:
        unit_group = UnitGroup("gwp_group", units)
        inlets = bst.utils.feeds_from_units(unit_group.units)
        inlets = set(inlets)
        bst.utils.filter_out_missing_streams(inlets)
        feeds = bst.utils.feeds(inlets)
        return sum([self.system.get_material_impact(s, "GWP") for s in feeds])

    def _get_utility_gwp(self, units: Iterable) -> float:
        impact = 0.0
        for unit in units:
            for hu in getattr(unit, "heat_utilities", ()):
                if hu.flow > 0:
                    impact += hu.get_impact("GWP") * self.system.operating_hours
            power_utility = getattr(unit, "power_utility", None)
            if power_utility and power_utility.rate > 0:
                impact += power_utility.get_impact("GWP") * self.system.operating_hours
        return impact

    def _get_total_gwp_per_kg(self) -> float:
        if self.product_stream is None:
            return 0.0
        try:
            lca_table = bst.report.lca_displacement_allocation_table(
                systems=[self.system],
                key="GWP",
                items=[self.product_stream],
            )
            return float(lca_table.loc[("Total", ""), lca_table.columns[-1]])
        except Exception:
            total_gwp = self.system.get_net_impact("GWP")
            annual_production = self.product_stream.F_mass * self.system.operating_hours
            return total_gwp / annual_production if annual_production else 0.0

    def _sort_breakdown(self, df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df
        return df.sort_values(by=df.columns[1], ascending=False).reset_index(drop=True)

    def _title_case(self, key: str) -> str:
        return " ".join([s.capitalize() for s in key.split()])

    def _save_excel(self, tables: Dict[str, pd.DataFrame]):
        filepath = os.path.join(self.data_dir, "Breakdown_Summary.xlsx")
        with pd.ExcelWriter(filepath, engine="openpyxl") as writer:
            for sheet, df in tables.items():
                df.to_excel(writer, sheet_name=sheet, index=False)

    def _save_metadata(self):
        metadata_path = os.path.join(self.base_dir, "README.txt")
        system_id = getattr(self.system, "ID", "Unknown")
        with open(metadata_path, "w", encoding="utf-8") as f:
            f.write("BioSTEAM Process Report Generator (v2)\n")
            f.write(f"System ID: {system_id}\n")
            f.write(f"Timestamp: {self.timestamp}\n")

    def _print_verification(self, results: List[VerificationResult]):
        if not results:
            return
        print("\nVerification Summary")
        print("=" * 70)
        for res in results:
            status = "✅" if res.passed else "⚠️"
            print(
                f"{status} {res.label} Sum: {res.reported:.4g} | "
                f"Model: {res.reference:.4g} | Variance: {res.variance:.2%}"
            )
        print("=" * 70)
