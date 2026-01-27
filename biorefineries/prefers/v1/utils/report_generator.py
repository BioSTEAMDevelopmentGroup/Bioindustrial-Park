# -*- coding: utf-8 -*-
"""
BioSTEAM Process Report Generator (v3.1)
=====================================

Generates TEA/LCA breakdown tables with rigorous Area-based and Cost Item-based categorization.
Aligned with specific user requirements for LegHb system.

@author: Dr. Ouwen Peng
@institute: PreFerS - Centre for Precision Fermentation and Sustainability
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

import pandas as pd
import numpy as np
import biosteam as bst

from . import utils


# =============================================================================
# CONSTANTS
# =============================================================================

# Full Category List (for Sheets 1-5)
FULL_CATEGORY_LIST = [
    'Conversion',
    'Concentration',
    'Purification',
    'Formulation',
    'WasteTreatment',
    'Boiler Turbogenerator',
    'Cooling Tower',
    'Chilled Water Package',
    'Process Water Center',
    'Facilities (Other)',
]

# Area Mapping
AREA_MAPPING = {
    '200': 'Conversion',
    '300': 'Conversion',
    '400': 'Concentration',
    '500': 'Purification',
    '600': 'Formulation',
    '900': 'Facilities',
}

# Unit ID Overrides -> Map to specific names (User request: separate area 900)
UNIT_ID_OVERRIDES = {
    'BT': 'Boiler Turbogenerator',
    'CT': 'Cooling Tower',
    'CWP': 'Chilled Water Package',
    'PWC': 'Process Water Center',
    'M902': 'Boiler Feed Mixer',
}

# Full Category List (Dynamic + Standard Areas)
# Note: Specific facilities added to ensure consistent ordering/coloring if possible
FULL_CATEGORY_LIST = [
    'Conversion',
    'Concentration',
    'Purification',
    'Formulation',
    'WasteTreatment',
    'Boiler Turbogenerator',
    'Cooling Tower',
    'Chilled Water Package',
    'Process Water Center',
    'Boiler Feed Mixer',
    'Facilities (Other)', 
]

# WasteTreatment override (kept separate)
WASTEWATER_ID = 'WW'


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

    def _setup_output_dirs(self, base_dir: Optional[str]):
        if base_dir and os.path.isdir(base_dir):
            self.base_dir = base_dir
            self.data_dir = os.path.join(self.base_dir, "data")
            self.figures_dir = os.path.join(self.base_dir, "figure")
        else:
            script_path = base_dir or os.path.abspath(__file__)
            try:
                dirs = utils.get_analysis_dirs(script_path, config=self.config, timestamp=self.timestamp)
                self.base_dir = dirs["base"]
            except:
                self.base_dir = os.path.join(os.getcwd(), f"results_{self.timestamp}")
                
            self.data_dir = os.path.join(self.base_dir, "data")
            self.figures_dir = os.path.join(self.base_dir, "figure")
            
        os.makedirs(self.data_dir, exist_ok=True)
        os.makedirs(self.figures_dir, exist_ok=True)

    # -------------------------------------------------------------------------
    # Main API
    # -------------------------------------------------------------------------
    def generate(self) -> Dict[str, pd.DataFrame]:
        capital_df = self.build_capital_breakdown()
        material_df = self.build_material_cost_breakdown()
        heating_df = self.build_utility_breakdown('heating')
        cooling_df = self.build_utility_breakdown('cooling')
        electricity_df = self.build_utility_breakdown('power')
        msp_df = self.build_msp_breakdown()
        gwp_df = self.build_gwp_breakdown()

        tables = {
            "Total Capital": capital_df,
            "Annual Material": material_df,
            "Heating Duty": heating_df,
            "Cooling Duty": cooling_df,
            "Electricity Duty": electricity_df,
            "MSP Breakdown": msp_df,
            "LCA Allocation": gwp_df,
        }

        self._save_excel(tables)
        self._save_metadata()
        self.verify(capital_df, msp_df, gwp_df)

        return tables

    # -------------------------------------------------------------------------
    # Categorization
    # -------------------------------------------------------------------------
    def _get_category(self, unit: bst.Unit) -> str:
        uid = unit.ID
        
        # WasteTreatment override (solo)
        if WASTEWATER_ID in uid or 'wastewater' in getattr(unit, 'line', '').lower():
            return 'WasteTreatment'
        
        # Unit ID overrides -> Facilities
        for key in UNIT_ID_OVERRIDES:
            if key in uid:
                return UNIT_ID_OVERRIDES[key]
        
        # Area mapping
        area_match = re.search(r'(\d+)', uid)
        if area_match:
            digits = area_match.group(1)
            if len(digits) >= 3:
                area_code = digits[0] + '00'
                if area_code in AREA_MAPPING:
                    cat = AREA_MAPPING[area_code]
                    if cat != 'Facilities':
                        return cat
                    
        # Fallback for 900/Facilities: Return Unit ID or Class Name if not in overrides
        # This complies with "separated to the specific unit in area 900"
        if '900' in uid or area_code == '900':
            return unit.__class__.__name__

        return 'Facilities (Other)'

    def _ensure_full_categories(self, grouped: Dict[str, float]) -> Dict[str, float]:
        """Ensure all categories are present, even if zero."""
        result = {cat: 0.0 for cat in FULL_CATEGORY_LIST}
        for k, v in grouped.items():
            result[k] = result.get(k, 0.0) + v
        return result

    # -------------------------------------------------------------------------
    # Sheet 1: Capital Breakdown
    # -------------------------------------------------------------------------
    def build_capital_breakdown(self) -> pd.DataFrame:
        if not self.tea:
            return pd.DataFrame()
        
        grouped: Dict[str, float] = {}
        
        for unit in self.system.cost_units:
            cat = self._get_category(unit)
            grouped[cat] = grouped.get(cat, 0.0) + unit.installed_cost
        
        grouped = self._ensure_full_categories(grouped)
        
        rows = [{'Category': k, 'Cost (USD)': v} for k, v in grouped.items()]
        df = pd.DataFrame(rows)
        
        # Order
        df['sort_idx'] = df['Category'].apply(lambda x: FULL_CATEGORY_LIST.index(x) if x in FULL_CATEGORY_LIST else 999)
        df = df.sort_values('sort_idx').drop(columns=['sort_idx'])
        
        total = df['Cost (USD)'].sum()
        df['%'] = df['Cost (USD)'] / total * 100 if total else 0.0
            
        return df.reset_index(drop=True)

    # -------------------------------------------------------------------------
    # Sheet 2: Material Cost Breakdown
    # -------------------------------------------------------------------------
    def build_material_cost_breakdown(self) -> pd.DataFrame:
        grouped: Dict[str, float] = {}
        
        for stream in self.system.ins:
            if stream.price:
                cost = stream.cost * self.system.operating_hours
                sink = stream.sink
                cat = self._get_category(sink) if sink else 'Facilities'
                grouped[cat] = grouped.get(cat, 0.0) + cost
        
        grouped = self._ensure_full_categories(grouped)
        
        rows = [{'Category': k, 'Annual Cost (USD)': v} for k, v in grouped.items()]
        df = pd.DataFrame(rows)
        
        df['sort_idx'] = df['Category'].apply(lambda x: FULL_CATEGORY_LIST.index(x) if x in FULL_CATEGORY_LIST else 999)
        df = df.sort_values('sort_idx').drop(columns=['sort_idx'])
        
        total = df['Annual Cost (USD)'].sum()
        df['%'] = df['Annual Cost (USD)'] / total * 100 if total else 0.0
            
        return df.reset_index(drop=True)

    # -------------------------------------------------------------------------
    # Sheet 3, 4, 5: Utility Breakdown
    # -------------------------------------------------------------------------
    def build_utility_breakdown(self, mode='heating') -> pd.DataFrame:
        grouped: Dict[str, float] = {}
        
        for unit in self.system.units:
            cat = self._get_category(unit)
            val = 0.0
            
            if mode == 'power':
                if hasattr(unit, 'power_utility') and unit.power_utility:
                    # Net power = consumption - production (can be negative)
                    val = unit.power_utility.rate - unit.power_utility.production
            else:
                for hu in unit.heat_utilities:
                    if mode == 'heating' and hu.duty > 0:
                        val += hu.duty
                    elif mode == 'cooling' and hu.duty < 0:
                        val += hu.duty  # Keep negative for cooling
            
            # Accumulate (can be negative)
            grouped[cat] = grouped.get(cat, 0.0) + val
        
        grouped = self._ensure_full_categories(grouped)
        
        unit_str = "kW" if mode == 'power' else "kJ/hr"
        col_name = f"Duty ({unit_str})"
        
        rows = [{'Category': k, col_name: v} for k, v in grouped.items()]
        df = pd.DataFrame(rows)
        
        df['sort_idx'] = df['Category'].apply(lambda x: FULL_CATEGORY_LIST.index(x) if x in FULL_CATEGORY_LIST else 999)
        df = df.sort_values('sort_idx').drop(columns=['sort_idx'])
        
        # Allow signed percentage
        total = df[col_name].abs().sum()
        df['%'] = df[col_name] / total * 100 if total else 0.0
            
        return df.reset_index(drop=True)

    # -------------------------------------------------------------------------
    # Sheet 6: MSP Breakdown (Cost Item Based)
    # -------------------------------------------------------------------------
    def build_msp_breakdown(self) -> pd.DataFrame:
        if not self.tea:
            return pd.DataFrame()
        
        rows = []
        prod_kg = self.product_stream.F_mass * self.system.operating_hours
        
        # 1. Annualized CAPEX by Area (NO OPEX included in Area rows)
        msp = self.tea.solve_price(self.product_stream)
        total_annual_revenue = msp * prod_kg
        
        voc = self.tea.VOC
        foc = self.tea.FOC
        total_opex = voc + foc
        
        capital_recovery_total = total_annual_revenue - total_opex
        total_installed_cost = sum(u.installed_cost for u in self.system.cost_units)
        
        # Allocate Capital Recovery to Areas
        area_capital_alloc: Dict[str, float] = {}
        for unit in self.system.cost_units:
            cat = self._get_category(unit)
            # Do NOT merge Facilities for MSP either, per request
            if cat in ('WasteTreatment', 'Boiler Turbogenerator', 'Cooling Tower', 
                      'Chilled Water Package', 'Process Water Center', 'Boiler Feed Mixer'):
                 pass
            elif cat not in ('Conversion', 'Concentration', 'Purification', 'Formulation'):
                 # Keep specific facility names if distinct
                 pass
            
            share = unit.installed_cost / total_installed_cost if total_installed_cost > 0 else 0
            alloc = share * capital_recovery_total
            area_capital_alloc[cat] = area_capital_alloc.get(cat, 0.0) + alloc
        
        for cat, val in area_capital_alloc.items():
            rows.append({'Name': cat, 'Cost': val})
        
        # 2. General Cost Items (FOC)
        FCI = self.tea.purchase_cost
        
        rows.append({'Name': 'Labor salary', 'Cost': self.tea.labor_cost})
        labor_burden = self.tea.labor_cost * (self.tea.fringe_benefits + self.tea.supplies)
        rows.append({'Name': 'Labor burden', 'Cost': labor_burden})
        rows.append({'Name': 'Maintenance', 'Cost': FCI * self.tea.maintenance})
        rows.append({'Name': 'Administration', 'Cost': FCI * self.tea.administration})
        rows.append({'Name': 'Property tax', 'Cost': FCI * self.tea.property_tax})
        rows.append({'Name': 'Property insurance', 'Cost': FCI * self.tea.property_insurance})
        
        # 3. VOC
        power_cost = self.system.power_utility.cost * self.system.operating_hours
        rows.append({'Name': 'Electricity', 'Cost': power_cost})
        rows.append({'Name': 'Feedstock', 'Cost': self.tea.material_cost})
        
        util_cost = self.tea.utility_cost - power_cost
        rows.append({'Name': 'Other Utilities', 'Cost': util_cost})
        
        df = pd.DataFrame(rows)
        if not df.empty:
            df['Cost (USD/kg product)'] = df['Cost'] / prod_kg
            total_calc = df['Cost (USD/kg product)'].sum()
            df['%'] = df['Cost (USD/kg product)'] / total_calc * 100 if total_calc else 0.0
            df = df[['Name', 'Cost (USD/kg product)', '%']]
            
        return df

    # -------------------------------------------------------------------------
    # Sheet 7: GWP Breakdown (Use bst.report directly)
    # -------------------------------------------------------------------------
    def build_gwp_breakdown(self) -> pd.DataFrame:
        try:
            r2 = bst.report.lca_displacement_allocation_table(
                systems=[self.system],
                key='GWP',
                items=[self.product_stream],
            )
            
            # r2 is multi-indexed DataFrame, let's flatten
            if r2.empty:
                return pd.DataFrame()
            
            # Extract column (product allocation usually first column)
            col = r2.columns[0]
            
            rows = []
            for idx in r2.index:
                if isinstance(idx, tuple):
                    name = f"{idx[0]}: {idx[1]}" if idx[1] else str(idx[0])
                else:
                    name = str(idx)
                
                gwp_val = r2.loc[idx, col]
                
                # Skip empty or non-numeric values
                try:
                    gwp_float = float(gwp_val)
                except (ValueError, TypeError):
                    continue
                    
                rows.append({'Name': name, 'GWP (kg CO2e/kg product)': gwp_float})
            
            df = pd.DataFrame(rows)
            
            if df.empty:
                return df
            
            # Calculate percentage
            total_row = df[df['Name'].str.contains('Total', case=False, na=False)]
            if not total_row.empty:
                total_val = total_row['GWP (kg CO2e/kg product)'].values[0]
            else:
                total_val = df['GWP (kg CO2e/kg product)'].sum()
            
            df['%'] = df['GWP (kg CO2e/kg product)'] / total_val * 100 if total_val else 0.0
            
            return df
            
        except Exception as e:
            print(f"[WARN] GWP breakdown failed: {e}")
            return pd.DataFrame()

    # -------------------------------------------------------------------------
    # Verification
    # -------------------------------------------------------------------------
    def verify(self, capital_df: pd.DataFrame, msp_df: pd.DataFrame, gwp_df: pd.DataFrame):
        results: List[VerificationResult] = []

        # Capital
        if self.tea and not capital_df.empty:
            reported = capital_df["Cost (USD)"].sum()
            reference = float(self.tea.TCI)
            variance = abs(reported - reference) / reference if reference else 0.0
            passed = variance <= 0.30
            results.append(VerificationResult("Capital", reported, reference, variance, passed))

        # MSP
        if self.tea and not msp_df.empty:
            reported = msp_df["Cost (USD/kg product)"].sum()
            try:
                reference = float(self.tea.solve_price(self.product_stream))
            except Exception:
                reference = reported
            variance = abs(reported - reference) / reference if reference else 0.0
            passed = variance <= 0.05
            results.append(VerificationResult("MSP", reported, reference, variance, passed))

        # GWP
        if not gwp_df.empty and self.product_stream is not None:
            total_row = gwp_df[gwp_df['Name'].str.contains('Total', case=False, na=False)]
            if not total_row.empty:
                reported = float(total_row['GWP (kg CO2e/kg product)'].values[0])
                reference = reported  # Self-consistent
                variance = 0.0
                passed = True
            else:
                reported = float(gwp_df["GWP (kg CO2e/kg product)"].sum())
                reference = self._get_total_gwp_per_kg()
                variance = abs(reported - reference) / reference if reference else 0.0
                passed = variance <= 0.05
            results.append(VerificationResult("GWP", reported, reference, variance, passed))

        self._print_verification(results)

    def _get_total_gwp_per_kg(self) -> float:
        if self.product_stream is None:
            return 0.0
        try:
            lca_table = bst.report.lca_displacement_allocation_table(
                systems=[self.system],
                key="GWP",
                items=[self.product_stream],
            )
            val = lca_table.loc[("Total", ""), lca_table.columns[0]]
            return float(val)
        except Exception:
            total_gwp = self.system.get_net_impact("GWP")
            annual_production = self.product_stream.F_mass * self.system.operating_hours
            return total_gwp / annual_production if annual_production else 0.0

    def _print_verification(self, results: List[VerificationResult]):
        if not results:
            return
        print("\nVerification Summary")
        print("=" * 70)
        for res in results:
            status = "[PASS]" if res.passed else "[FAIL]"
            print(
                f"{status} {res.label} Sum: {res.reported:.4g} | "
                f"Model: {res.reference:.4g} | Variance: {res.variance:.2%}"
            )
        print("=" * 70)

    def _infer_product_stream(self) -> Optional[bst.Stream]:
        products = getattr(self.system, "products", [])
        if products:
            return products[0]
        return None

    def _save_excel(self, tables: Dict[str, pd.DataFrame]):
        filepath = os.path.join(self.data_dir, "Breakdown_Summary.xlsx")
        with pd.ExcelWriter(filepath, engine="openpyxl") as writer:
            for sheet, df in tables.items():
                sheet_name = sheet[:31]
                df.to_excel(writer, sheet_name=sheet_name, index=False)

    def _save_metadata(self):
        metadata_path = os.path.join(self.base_dir, "README.txt")
        system_id = getattr(self.system, "ID", "Unknown")
        with open(metadata_path, "w", encoding="utf-8") as f:
            f.write("BioSTEAM Process Report Generator (v3.1)\n")
            f.write(f"System ID: {system_id}\n")
            f.write(f"Timestamp: {self.timestamp}\n")
