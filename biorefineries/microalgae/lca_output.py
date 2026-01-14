# -*- coding: utf-8 -*-
"""
Microalgae LCA Output Module

基于succinic项目的lca_output.py模式，为microalgae项目提供LCA分析、报告和数据导出功能。

@author: Xingdong Shi
@version: 0.0.1
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# 添加项目路径以便导入
current_dir = Path(__file__).parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

try:
    from . import system
    from ._chemicals import chems
    from .lca import LCA, create_microalgae_lca
except ImportError:
    # 如果相对导入失败，尝试绝对导入
    try:
        from microalgae import system
        from microalgae._chemicals import chems
        from microalgae.lca import LCA, create_microalgae_lca
    except ImportError:
        print("警告：无法导入microalgae模块，请确保在正确的环境中运行")

def get_gwp_breakdown_data(lca):
    """
    获取GWP分解数据，按照microalgae项目的计算方式
    """
    material_GWP_breakdown = lca.get_material_impact_breakdown('GWP_100')
    
    gwp_breakdown = {
        'feedstock*': lca.feedstock_GWP,
        'H2SO4': material_GWP_breakdown.get('H2SO4', 0),
        'NaOH': material_GWP_breakdown.get('NaOH', 0),
        'NH4OH': material_GWP_breakdown.get('NH4OH', 0),
        'CalciumDihydroxide': material_GWP_breakdown.get('CalciumDihydroxide', 0),
        'Ethanol': material_GWP_breakdown.get('Ethanol', 0),
        'Octanol': material_GWP_breakdown.get('Octanol', 0),
        'GlucoAmylase': material_GWP_breakdown.get('GlucoAmylase', 0),
        'AlphaAmylase': material_GWP_breakdown.get('AlphaAmylase', 0),
        'CH4': material_GWP_breakdown.get('CH4', 0),
        'net electricity': lca.net_electricity_GWP,
        'direct non-biogenic emissions': lca.direct_non_biogenic_emissions_GWP,
    }
    
    # 计算正值的总和，然后归一化
    tot_positive_GWP = sum([v for v in gwp_breakdown.values() if v > 0])
    if tot_positive_GWP > 0:
        for k, v in gwp_breakdown.items():
            gwp_breakdown[k] = v / tot_positive_GWP
    else: # Fallback if no positive values
        tot_abs_GWP = sum([abs(v) for v in gwp_breakdown.values()])
        if tot_abs_GWP > 0:
            for k, v in gwp_breakdown.items():
                gwp_breakdown[k] = v / tot_abs_GWP
    
    return gwp_breakdown

def export_gwp_breakdown_data(lca, output_file=None):
    """
    导出GWP分解数据到CSV文件
    
    参数
    ----------
    lca : LCA
        LCA对象
    output_file : str, 可选
        输出文件路径，默认为'microalgae_gwp_breakdown.csv'
    """
    if output_file is None:
        output_file = 'microalgae_gwp_breakdown.csv'
    
    gwp_breakdown = get_gwp_breakdown_data(lca)
    
    # 创建DataFrame
    df = pd.DataFrame(list(gwp_breakdown.items()), columns=['Component', 'Percentage'])
    df['Percentage'] = df['Percentage'] * 100  # 转换为百分比
    
    # 保存到CSV
    df.to_csv(output_file, index=False)
    print(f"GWP分解数据已保存到: {output_file}")
    
    return df

def print_microalgae_lca_results(lca):
    """
    打印microalgae LCA结果
    
    参数
    ----------
    lca : LCA
        LCA对象
    """
    print("\n" + "="*60)
    print("MICROALGAE 生物炼制厂生命周期评估结果")
    print("="*60)
    
    print(f"\n功能单位: {lca.functional_unit}")
    print(f"系统碳平衡: {lca.system_carbon_balance:.4f}")
    
    # 主要影响类别结果
    print(f"\n----- 主要影响类别 -----")
    print(f"GWP (kg CO2e/{lca.functional_unit}): {lca.GWP:.4f}")
    print(f"FEC (MJ/{lca.functional_unit}): {lca.FEC:.4f}")
    print(f"WC (m³/{lca.functional_unit}): {lca.WC:.4f}")
    
    # GWP分解
    print(f"\n----- GWP分解 -----")
    print(f"物质影响: {lca.material_GWP:.4f} ({lca.material_GWP/lca.GWP*100 if lca.GWP else 0:.1f}%)")
    print(f"原料影响: {lca.feedstock_GWP:.4f} ({lca.feedstock_GWP/lca.GWP*100 if lca.GWP else 0:.1f}%)")
    print(f"电力影响: {lca.net_electricity_GWP:.4f} ({lca.net_electricity_GWP/lca.GWP*100 if lca.GWP else 0:.1f}%)")
    print(f"直接排放(非生物源): {lca.direct_non_biogenic_emissions_GWP:.4f} ({lca.direct_non_biogenic_emissions_GWP/lca.GWP*100 if lca.GWP else 0:.1f}%)")
    
    # 物质影响明细
    print(f"\n----- 物质影响明细 (kg CO2e/{lca.functional_unit}) -----")
    material_breakdown = lca.get_material_impact_breakdown('GWP_100')
    for component, impact in material_breakdown.items():
        if abs(impact) > 1e-6:
            print(f"{component:<20}: {impact:8.4f}")
    
    # 电力影响详情
    print(f"\n----- 电力影响详情 -----")
    elec_CF = lca.CFs['GWP_100']['Electricity']
    net_kWh = lca.net_electricity
    per_unit = net_kWh * elec_CF / lca.functional_quantity_per_h
    print(f"净电力消耗: {net_kWh:.3f} kWh/h")
    print(f"电力影响因子: {elec_CF} kg CO2e/kWh")
    print(f"电力影响: {per_unit:.4f} kg CO2e/{lca.functional_unit}")
    
    # 原料影响详情
    print(f"\n----- 原料影响详情 -----")
    for ID, (stream, mass_kind) in lca.complex_feeds.items():
        impact = lca.get_complex_feed_impact_by_ID('GWP_100', ID)
        print(f"{ID:<20}: {impact:8.4f} kg CO2e/{lca.functional_unit}")
    
    # 直接排放详情
    print(f"\n----- 直接排放详情 -----")
    print(f"总直接排放GWP: {lca.direct_emissions_GWP:.4f}")
    print(f"生物来源排放GWP: {lca.biogenic_emissions_GWP:.4f}")
    print(f"非生物来源直接排放GWP: {lca.direct_non_biogenic_emissions_GWP:.4f}")
    
    # 排放流详细
    print(f"\n排放流详细:")
    for stream in lca.emissions:
        C = stream.get_atomic_flow("C")
        if C > 0:
            stream_gwp = C * 44 / lca.functional_quantity_per_h
            if abs(stream_gwp) > 1e-6:
                print(f"  {stream.ID:<20}: {stream_gwp:8.4f} kg CO2e/{lca.functional_unit}")
    
    # GWP分解数据
    print(f"\n----- GWP分解数据 -----")
    gwp_breakdown = get_gwp_breakdown_data(lca)
    
    print("GWP分解数据:")
    print("="*50)
    total_percentage = 0
    for component, percentage in gwp_breakdown.items():
        print(f"{component:<30} {percentage*100:6.2f}%")
        if percentage > 0:  # 只计算正值
            total_percentage += percentage * 100
    print("-"*50)
    print(f"{'总计':<30} {total_percentage:6.2f}%")
    
    print("\n" + "="*60)

def create_microalgae_lca_scenario(system, main_product, main_product_chemical_IDs, boiler=None):
    """
    创建microalgae LCA场景
    
    参数
    ----------
    system : System
        微藻生产系统
    main_product : Stream
        主产品流
    main_product_chemical_IDs : list
        主产品化学品ID列表
    boiler : Unit, 可选
        锅炉设备单元
        
    返回
    ----------
    LCA
        生命周期评估对象
    """
    # 创建LCA对象
    lca = create_microalgae_lca(
        system=system,
        main_product=main_product,
        main_product_chemical_IDs=main_product_chemical_IDs,
        boiler=boiler,
        add_EOL_GWP=True
    )
    
    return lca

def export_lca_results_to_excel(lca, output_file=None):
    """
    导出LCA结果到Excel文件
    
    参数
    ----------
    lca : LCA
        LCA对象
    output_file : str, 可选
        输出文件路径，默认为'microalgae_lca_results.xlsx'
    """
    if output_file is None:
        output_file = 'microalgae_lca_results.xlsx'
    
    # 创建结果字典
    results = {
        '基本信息': {
            '功能单位': [lca.functional_unit],
            '系统碳平衡': [lca.system_carbon_balance],
            '主产品': [lca.main_product.ID],
        },
        '主要影响类别': {
            'GWP (kg CO2e/功能单位)': [lca.GWP],
            'FEC (MJ/功能单位)': [lca.FEC],
            'WC (m³/功能单位)': [lca.WC],
        },
        'GWP分解': {
            '物质影响 (kg CO2e/功能单位)': [lca.material_GWP],
            '原料影响 (kg CO2e/功能单位)': [lca.feedstock_GWP],
            '电力影响 (kg CO2e/功能单位)': [lca.net_electricity_GWP],
            '直接排放(非生物源) (kg CO2e/功能单位)': [lca.direct_non_biogenic_emissions_GWP],
        },
        '物质影响明细': lca.get_material_impact_breakdown('GWP_100'),
        'GWP分解百分比': get_gwp_breakdown_data(lca)
    }
    
    # 创建Excel文件
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # 基本信息
        pd.DataFrame(results['基本信息']).to_excel(writer, sheet_name='基本信息', index=False)
        
        # 主要影响类别
        pd.DataFrame(results['主要影响类别']).to_excel(writer, sheet_name='主要影响类别', index=False)
        
        # GWP分解
        pd.DataFrame(results['GWP分解']).to_excel(writer, sheet_name='GWP分解', index=False)
        
        # 物质影响明细
        material_df = pd.DataFrame(list(results['物质影响明细'].items()), 
                                 columns=['物质', '影响 (kg CO2e/功能单位)'])
        material_df.to_excel(writer, sheet_name='物质影响明细', index=False)
        
        # GWP分解百分比
        gwp_percent_df = pd.DataFrame(list(results['GWP分解百分比'].items()), 
                                     columns=['组分', '百分比'])
        gwp_percent_df['百分比'] = gwp_percent_df['百分比'] * 100
        gwp_percent_df.to_excel(writer, sheet_name='GWP分解百分比', index=False)
    
    print(f"LCA结果已导出到: {output_file}")
    return output_file

def compare_lca_scenarios(scenarios, output_file=None):
    """
    比较多个LCA场景
    
    参数
    ----------
    scenarios : dict
        场景字典，格式: {scenario_name: lca_object}
    output_file : str, 可选
        输出文件路径
    """
    if output_file is None:
        output_file = 'microalgae_lca_comparison.xlsx'
    
    # 创建比较数据
    comparison_data = {}
    
    for scenario_name, lca in scenarios.items():
        comparison_data[scenario_name] = {
            'GWP (kg CO2e/功能单位)': lca.GWP,
            'FEC (MJ/功能单位)': lca.FEC,
            'WC (m³/功能单位)': lca.WC,
            '物质影响 (kg CO2e/功能单位)': lca.material_GWP,
            '原料影响 (kg CO2e/功能单位)': lca.feedstock_GWP,
            '电力影响 (kg CO2e/功能单位)': lca.net_electricity_GWP,
            '直接排放(非生物源) (kg CO2e/功能单位)': lca.direct_non_biogenic_emissions_GWP,
        }
    
    # 创建DataFrame并保存
    df = pd.DataFrame(comparison_data).T
    df.to_excel(output_file)
    
    print(f"场景比较结果已保存到: {output_file}")
    return df

def generate_lca_report(lca, output_file=None):
    """
    生成详细的LCA报告
    
    参数
    ----------
    lca : LCA
        LCA对象
    output_file : str, 可选
        输出文件路径
    """
    if output_file is None:
        output_file = 'microalgae_lca_report.txt'
    
    report = lca.generate_report(output_file)
    
    if output_file:
        print(f"详细报告已保存到: {output_file}")
    
    return report

# 主程序
if __name__ == "__main__":
    import biosteam as bst
    
    # 初始化系统
    bst.settings.set_thermo(chems)
    microalgae_sys = system.create_microalgae_MCCA_production_sys()
    microalgae_sys.simulate()
    
    # 获取主产品流和主产品化学品ID
    main_product = microalgae_sys.flowsheet.stream.caproic_acid_product
    main_product_chemical_IDs = ['CaproicAcid']
    
    # 尝试找到锅炉涡轮发电机
    boiler = None
    for unit in microalgae_sys.units:
        # 优先查找锅炉涡轮发电机
        if 'BT' in unit.ID and hasattr(unit, 'turbogenerator_efficiency'):
            boiler = unit
            break
    
    # 如果没有找到锅炉涡轮发电机，查找普通锅炉
    if boiler is None:
        for unit in microalgae_sys.units:
            if hasattr(unit, 'heat_utilities') and unit.heat_utilities:
                boiler = unit
                break
    
    # 如果没有找到锅炉，使用第一个热交换器
    if boiler is None:
        for unit in microalgae_sys.units:
            if unit.__class__.__name__ == 'HXutility':
                boiler = unit
                break
    
    # 创建LCA对象
    lca = create_microalgae_lca_scenario(
        system=microalgae_sys,
        main_product=main_product,
        main_product_chemical_IDs=main_product_chemical_IDs,
        boiler=boiler,
        add_EOL_GWP=True
    )
    
    # 打印结果
    print_microalgae_lca_results(lca)
    
    # 导出GWP分解数据
    export_gwp_breakdown_data(lca)
    
    # 导出完整结果到Excel
    export_lca_results_to_excel(lca)
    
    # 生成详细报告
    generate_lca_report(lca) 