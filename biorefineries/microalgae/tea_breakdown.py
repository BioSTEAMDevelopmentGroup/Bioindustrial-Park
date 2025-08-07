# -*- coding: utf-8 -*-
"""
Microalgae TEA Breakdown Analysis

@author: Xingdong Shi
"""

import pandas as pd
import numpy as np
import biosteam as bst
from .model_utils import get_unit_groups
from .system import microalgae_mcca_sys, microalgae_tea

class BTElectricityRevenueConfig:
    def __init__(self, 
                 electricity_price=0.07,  # USD/kWh
                 include_in_material_cost=True,
                 include_in_operating_cost=False):
        self.electricity_price = electricity_price
        self.include_in_material_cost = include_in_material_cost
        self.include_in_operating_cost = include_in_operating_cost

bt_revenue_config = BTElectricityRevenueConfig()

def create_tea_breakdown_data(unit_groups=None, tea=microalgae_tea, print_output=False, fractions=False):
    unit_groups = unit_groups or get_unit_groups()
    system_electricity_revenue = calculate_system_electricity_revenue()
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
                        for unit in microalgae_mcca_sys.units:
                            if hasattr(unit, 'heat_utilities'):
                                for hu in unit.heat_utilities:
                                    if hasattr(hu, 'duty') and hu.duty > 0:
                                        denominator += hu.duty / 1e6
                    elif 'heating' in metric_name.lower():
                        denominator = 0
                        for unit in microalgae_mcca_sys.units:
                            if hasattr(unit, 'heat_utilities'):
                                for hu in unit.heat_utilities:
                                    if hasattr(hu, 'duty') and hu.duty < 0:
                                        denominator -= hu.duty / 1e6
                    elif 'electricity' in metric_name.lower() or 'power' in metric_name.lower():
                        if 'consumption' in metric_name.lower():
                            denominator = microalgae_mcca_sys.power_utility.consumption / 1e3
                        elif 'production' in metric_name.lower():
                            denominator = microalgae_mcca_sys.power_utility.production / 1e3
                    elif 'material' in metric_name.lower():
                        denominator = tea.material_cost / tea.operating_hours
                metric_value = metric() if callable(metric) else 0
                if ug.name == 'BT':
                    should_add_revenue = (
                        (bt_revenue_config.include_in_material_cost and 'material' in metric.name.lower()) or
                        (bt_revenue_config.include_in_operating_cost and ('operating' in metric.name.lower() or 'variable' in metric.name.lower()))
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

def calculate_system_electricity_revenue():   
    power_breakdown = get_system_power_breakdown()
    net_power_demand = power_breakdown['net_electricity']  # kW
    operating_hours = microalgae_tea.operating_hours
    
    if net_power_demand < 0:
        surplus_electricity = abs(net_power_demand)  # kW
        annual_electricity_revenue = surplus_electricity * bt_revenue_config.electricity_price * operating_hours
        return -annual_electricity_revenue 
    else:
        electricity_cost = net_power_demand * bt_revenue_config.electricity_price * operating_hours
        return electricity_cost 

def create_tea_breakdown_dataframe(unit_groups=None, fraction=True, 
                                   scale_fractions_to_positive_values=False):
    unit_groups = unit_groups or get_unit_groups()
    df = bst.UnitGroup.df_from_groups(
        unit_groups, 
        fraction=fraction,
        scale_fractions_to_positive_values=scale_fractions_to_positive_values
    )
    
    system_electricity_revenue = calculate_system_electricity_revenue()
    if system_electricity_revenue != 0 and 'BT' in df.index:
        if 'Material cost' in df.columns:
            if fraction:
                total_annual_material_cost = microalgae_tea.material_cost
                revenue_fraction = system_electricity_revenue / total_annual_material_cost
                df.loc['BT', 'Material cost'] += revenue_fraction
            else:
                revenue_hourly = system_electricity_revenue / microalgae_tea.operating_hours
                df.loc['BT', 'Material cost'] += revenue_hourly
    return df

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

def get_system_power_breakdown():
    power_breakdown = {
        'total_consumption': 0,
        'total_production': 0,
        'net_electricity': 0,
        'unit_details': {}
    }
    
    for unit in microalgae_mcca_sys.units:
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
    power_breakdown['net_electricity'] = (power_breakdown['total_consumption'] - 
                                          power_breakdown['total_production'])
    return power_breakdown

def analyze_microalgae_tea_breakdown():
    microalgae_mcca_sys.simulate()
    system_electricity_revenue = calculate_system_electricity_revenue()
    power_breakdown = get_system_power_breakdown()
    unit_groups = get_unit_groups()
    breakdown_data = create_tea_breakdown_data(unit_groups, microalgae_tea, print_output=True)
    df_breakdown = create_tea_breakdown_dataframe(
        unit_groups,
        fraction=True,
        scale_fractions_to_positive_values=False
    )
    df_breakdown = create_tea_breakdown_dataframe(
        unit_groups,
        fraction=True,
        scale_fractions_to_positive_values=False
    )
    cost_breakdown = get_cost_breakdown_by_category(microalgae_tea)
    results = {
        'breakdown_data': breakdown_data,
        'breakdown_dataframe': df_breakdown,
        'cost_breakdown': cost_breakdown,
        'unit_groups': unit_groups,
        'tea_system': microalgae_mcca_sys,
        'power_breakdown': power_breakdown,
        'system_electricity_revenue': system_electricity_revenue
    }
    return results

if __name__ == "__main__":
    analyze_microalgae_tea_breakdown()


