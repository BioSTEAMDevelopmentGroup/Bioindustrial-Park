# -*- coding: utf-8 -*-
"""
Created on Sat July 20 17:50:00 2025

Microalgae biorefinery to produce medium chain fatty acids 
by anaerobic fermentation without external electron donor addition- TEA

References
----------
[1] BioSTEAM Documentation: 
    https://biosteam.readthedocs.io/en/latest/tutorial/Creating_a_System.html
[2] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
[3] 3-Hydroxypropionic acid biorefineries project:
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/HP
[4] Succinic projest
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/succinic

@author: Xingdong Shi
@version: 0.0.9
"""

import pandas as pd
import numpy as np
import biosteam as bst
from biorefineries.tea.cellulosic_ethanol_tea import CellulosicEthanolTEA as create_tea

def compute_labor_cost(dry_tpd: float,
                        base_tpd: float = 2205.0,
                        base_cost: float = 3212962.0,
                        exponent: float = 0.2,
                        floor_tpd: float = 100.0,
                        floor_cost: float = 0.5e6) -> float:
    if dry_tpd < floor_tpd:
        return floor_cost
    return base_cost * (dry_tpd / base_tpd) ** exponent

# Generic TEA builder for microalgae systems
def microalgae_tea(system):
    u = system.flowsheet.unit
    dry_tpd = u.U101.ins[0].F_mass * 24 / 1000
    tea = create_tea(system=system, IRR=0.10, duration=(2024, 2045),
        depreciation='MACRS7', income_tax=0.21, 
            operating_days=330,
        lang_factor= None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
            OSBL_units=(u.CT, u.CWP, u.ADP, u.PWC, u.BT601),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=max(0.5e6, compute_labor_cost(dry_tpd)),
            labor_burden=0.90, property_insurance=0.007, maintenance=0.03, boiler_turbogenerator=u.BT601,
        steam_power_depreciation='MACRS20')
    return tea

class BTElectricityRevenueConfig:
    def __init__(self, 
                 electricity_price=0.07,  # USD/kWh
                 include_in_material_cost=True,
                 include_in_operating_cost=False):
        self.electricity_price = electricity_price
        self.include_in_material_cost = include_in_material_cost
        self.include_in_operating_cost = include_in_operating_cost


bt_revenue_config = BTElectricityRevenueConfig()


def get_system_power_breakdown(system):
    power_breakdown = {
        'total_consumption': 0,
        'total_production': 0,
        'net_electricity': 0,
        'unit_details': {}
    }
    for unit in system.units:
        if hasattr(unit, 'power_utility') and unit.power_utility:
            pu = unit.power_utility
            consumption = pu.consumption
            production = pu.production
            if consumption != 0 or production != 0:
                power_breakdown['unit_details'][unit.ID] = {
                    'consumption': consumption,
                    'production': production,
                    'net': consumption - production
                }
                power_breakdown['total_consumption'] += consumption
                power_breakdown['total_production'] += production
    power_breakdown['net_electricity'] = (
        power_breakdown['total_consumption'] - power_breakdown['total_production']
    )
    return power_breakdown


def calculate_system_electricity_revenue(system, tea, config: BTElectricityRevenueConfig | None = None):   
    config = config or bt_revenue_config
    power_breakdown = get_system_power_breakdown(system)
    net_power_demand = power_breakdown['net_electricity']  # kW
    operating_hours = tea.operating_hours
    if net_power_demand < 0:
        surplus_electricity = abs(net_power_demand)  # kW
        annual_electricity_revenue = surplus_electricity * config.electricity_price * operating_hours
        return -annual_electricity_revenue 
    else:
        electricity_cost = net_power_demand * config.electricity_price * operating_hours
        return electricity_cost 


def create_tea_breakdown_data(system, tea, unit_groups, 
                              print_output: bool = False, 
                              fractions: bool = False,
                              config: BTElectricityRevenueConfig | None = None):
    config = config or bt_revenue_config
    system_electricity_revenue = calculate_system_electricity_revenue(system, tea, config)
    first_metrics_group = None
    for ug in unit_groups:
        if hasattr(ug, 'metrics') and ug.metrics:
            first_metrics_group = ug
            break
    metric_breakdowns = {metric.name: {} for metric in first_metrics_group.metrics}
    for ug in unit_groups:
        if hasattr(ug, 'metrics') and ug.metrics:
            for metric in ug.metrics:
                denominator = 1.0
                if fractions:
                    metric_name = metric.name
                    if 'cost' in metric_name.lower() or 'installed' in metric_name.lower():
                        denominator = tea.installed_equipment_cost / 1e6
                    elif 'cooling' in metric_name.lower():
                        denominator = 0
                        for unit in system.units:
                            if hasattr(unit, 'heat_utilities'):
                                for hu in unit.heat_utilities:
                                    if hasattr(hu, 'duty') and hu.duty > 0:
                                        denominator += hu.duty / 1e6
                    elif 'heating' in metric_name.lower():
                        denominator = 0
                        for unit in system.units:
                            if hasattr(unit, 'heat_utilities'):
                                for hu in unit.heat_utilities:
                                    if hasattr(hu, 'duty') and hu.duty < 0:
                                        denominator -= hu.duty / 1e6
                    elif 'electricity' in metric_name.lower() or 'power' in metric_name.lower():
                        if 'consumption' in metric_name.lower():
                            denominator = system.power_utility.consumption / 1e3
                        elif 'production' in metric_name.lower():
                            denominator = system.power_utility.production / 1e3
                    elif 'material' in metric_name.lower():
                        denominator = tea.material_cost / tea.operating_hours
                metric_value = metric() if callable(metric) else 0
                if ug.name == 'BT':
                    should_add_revenue = (
                        (config.include_in_material_cost and 'material' in metric.name.lower()) or
                        (config.include_in_operating_cost and ('operating' in metric.name.lower() or 'variable' in metric.name.lower()))
                    )
                    if should_add_revenue:
                        if fractions:
                            electricity_revenue_to_add = system_electricity_revenue
                        else:
                            electricity_revenue_to_add = system_electricity_revenue / tea.operating_hours
                        metric_value += electricity_revenue_to_add
                if ug.name == 'Storage':
                    continue
                elif ug.name == 'Other facilities':
                    storage_value = 0
                    for storage_ug in unit_groups:
                        if storage_ug.name == 'Storage':
                            if hasattr(storage_ug, 'metrics') and storage_ug.metrics:
                                storage_metric = storage_ug.metrics[ug.metrics.index(metric)] if ug.metrics.index(metric) < len(storage_ug.metrics) else None
                                storage_value = storage_metric() if callable(storage_metric) else 0
                            break
                    combined_name = 'Storage and ' + ug.name
                    metric_breakdowns[metric.name][combined_name] = (metric_value + storage_value) / denominator
                else:
                    metric_breakdowns[metric.name][ug.name] = metric_value / denominator
    if print_output and first_metrics_group:
        for metric in first_metrics_group.metrics:
            print(f"\n\n----- {metric.name} ({metric.units}) -----")
            metric_breakdown = metric_breakdowns.get(metric.name, {})
            for group_name, value in metric_breakdown.items():
                print(f"{group_name}: {value:.3f}")
    return metric_breakdowns


def create_tea_breakdown_dataframe(system, tea, unit_groups, fraction=True, 
                                   scale_fractions_to_positive_values=False,
                                   config: BTElectricityRevenueConfig | None = None):
    config = config or bt_revenue_config
    df = bst.UnitGroup.df_from_groups(
        unit_groups, 
        fraction=fraction,
        scale_fractions_to_positive_values=scale_fractions_to_positive_values
    )
    system_electricity_revenue = calculate_system_electricity_revenue(system, tea, config)
    if system_electricity_revenue != 0 and 'BT' in df.index:
        if 'Material cost' in df.columns:
            if fraction:
                total_annual_material_cost = tea.material_cost
                revenue_fraction = system_electricity_revenue / total_annual_material_cost
                df.loc['BT', 'Material cost'] += revenue_fraction
            else:
                revenue_hourly = system_electricity_revenue / tea.operating_hours
                df.loc['BT', 'Material cost'] += revenue_hourly
    return df


def get_unit_groups_for_system(system):
    try:
        return bst.UnitGroup.group_by_area(system)
    except Exception:
        try:
            from .model_utils import get_unit_groups  # type: ignore
            return get_unit_groups()
        except Exception:
            return []


def get_cost_breakdown_by_category(tea):
    breakdown = {}
    breakdown['Installed equipment cost'] = tea.installed_equipment_cost / 1e6
    breakdown['ISBL installed equipment cost'] = getattr(tea, 'ISBL_installed_equipment_cost', 0) / 1e6
    breakdown['OSBL installed equipment cost'] = getattr(tea, 'OSBL_installed_equipment_cost', 0) / 1e6
    breakdown['Fixed operating cost'] = tea.FOC / 1e6
    breakdown['Variable operating cost'] = tea.VOC / 1e6
    breakdown['Material cost'] = tea.material_cost / 1e6
    breakdown['Utility cost'] = tea.utility_cost / 1e6
    breakdown['Total capital investment'] = tea.TCI / 1e6
    breakdown['Fixed capital investment'] = tea.FCI / 1e6
    breakdown['Total depreciable capital'] = tea.TDC / 1e6
    breakdown['Direct permanent investment'] = tea.DPI / 1e6
    return breakdown


def analyze_microalgae_tea_breakdown():
    from .model_utils import get_unit_groups
    from .system import microalgae_mcca_sys
    microalgae_mcca_sys.simulate()
    unit_groups = get_unit_groups()
    tea_obj = microalgae_tea(microalgae_mcca_sys)
    _ = create_tea_breakdown_data(
        system=microalgae_mcca_sys,
        tea=tea_obj,
        unit_groups=unit_groups,
        print_output=True
    )
    df_breakdown = create_tea_breakdown_dataframe(
        system=microalgae_mcca_sys,
        tea=tea_obj,
        unit_groups=unit_groups,
        fraction=True,
        scale_fractions_to_positive_values=False
    )
    cost_breakdown = get_cost_breakdown_by_category(tea_obj)
    results = {
        'breakdown_dataframe': df_breakdown,
        'cost_breakdown': cost_breakdown,
        'unit_groups': unit_groups,
        'tea_system': microalgae_mcca_sys,
        'system_electricity_revenue': calculate_system_electricity_revenue(microalgae_mcca_sys, tea_obj, bt_revenue_config)
    }
    return results


if __name__ == "__main__":
    analyze_microalgae_tea_breakdown()




