"""
Created on Thu Jun 20 15:30:00 2024

Visualization module for Microalgae biorefinery to produce medium chain fatty acids
This script creates technical, economic and environmental figures for the system

@author: Xingdong Shi
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import biosteam as bst
import thermosteam as tmo
from biosteam.evaluation import Model, Metric
from chaospy import distributions as shape

from .lca import create_microalgae_lca
from .system import create_microalgae_MCCA_production_sys, compute_labor_cost
from .tea import create_tea

# Set style for publication quality plots
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 7

def calculate_breakdown_data(tea_obj, system):
    """从TEA对象和系统中提取实际的成本分布数据"""
    
    try:
        # 总投资成本
        total_TCI = tea_obj.TCI
        total_VOC = tea_obj.VOC  
        total_FOC = tea_obj.FOC  
        
        # 按类别分组单元
        feedstock_units = []
        conversion_units = []
        separation_units = []
        wastewater_units = []
        facilities_units = []
        
        # 分类系统中的单元
        for unit in system.units:
            unit_name = unit.__class__.__name__
            unit_id = unit.ID
            
            if 'U1' in unit_id or 'feed' in unit_id.lower():
                feedstock_units.append(unit)
            elif 'R3' in unit_id or 'ferment' in unit_id.lower():
                conversion_units.append(unit)
            elif any(x in unit_id for x in ['S3', 'D3', 'F3', 'M3']):
                separation_units.append(unit)
            elif 'waste' in unit_id.lower() or 'W' in unit_id:
                wastewater_units.append(unit)
            else:
                facilities_units.append(unit)
        
        # 计算各类别的设备成本
        feedstock_cost = sum(getattr(unit, 'purchase_cost', 0) for unit in feedstock_units)
        conversion_cost = sum(getattr(unit, 'purchase_cost', 0) for unit in conversion_units)
        separation_cost = sum(getattr(unit, 'purchase_cost', 0) for unit in separation_units)
        wastewater_cost = sum(getattr(unit, 'purchase_cost', 0) for unit in wastewater_units)
        facilities_cost = sum(getattr(unit, 'purchase_cost', 0) for unit in facilities_units)
        
        # 归一化为百分比
        total_equipment_cost = feedstock_cost + conversion_cost + separation_cost + wastewater_cost + facilities_cost
        if total_equipment_cost > 0:
            equipment_costs = {
                'Feedstock': feedstock_cost / total_equipment_cost * 100,
                'Conversion': conversion_cost / total_equipment_cost * 100,
                'Separation': separation_cost / total_equipment_cost * 100,
                'Wastewater': wastewater_cost / total_equipment_cost * 100,
                'Cooling utility facilities': 0,
                'Boiler & turbogenerator': 0,
                'Other facilities': facilities_cost / total_equipment_cost * 100,
                'Heat exchanger network': 0,
                'Natural gas (steam)': 0,
                'Natural gas (drying)': 0,
                'Fixed operating costs': 0,
                'Credits': 0
            }
        else:
            # 回退到估算值
            equipment_costs = {k: np.nan for k in [
                'Feedstock','Conversion','Separation','Wastewater',
                'Cooling utility facilities','Boiler & turbogenerator','Other facilities',
                'Heat exchanger network','Natural gas (steam)','Natural gas (drying)',
                'Fixed operating costs','Credits']}
        
        # 计算公用设施消耗
        total_cooling = total_heating = total_electricity = 0
        
        for unit in system.units:
            # 冷却负荷
            if hasattr(unit, 'heat_utilities'):
                for hu in unit.heat_utilities:
                    if hu.duty < 0:  
                        total_cooling += abs(hu.duty)
                    else:  
                        total_heating += hu.duty
            
            # 电力消耗
            if hasattr(unit, 'power_utility') and unit.power_utility:
                total_electricity += abs(unit.power_utility.power)
        
        # 按单元类型分配公用设施负荷 (简化分配)
        cooling_duty = {k: np.nan for k in equipment_costs}
        
        heating_duty = {k: np.nan for k in equipment_costs}
        
        electricity = {k: np.nan for k in equipment_costs}
        
        # 运营成本分配 (基于FOC和VOC)
        operating_cost = {k: np.nan for k in equipment_costs}
        
        # 输出实际数值供参考
        print(f"实际TEA数据:")
        print(f"  总投资成本 (TCI): ${total_TCI/1e6:.1f} M")
        print(f"  年度变动成本 (VOC): ${total_VOC/1e6:.1f} M/y")
        print(f"  年度固定成本 (FOC): ${total_FOC/1e6:.1f} M/y")
        print(f"  总加热负荷: {total_heating/1e6:.1f} MW")
        print(f"  总冷却负荷: {total_cooling/1e6:.1f} MW")
        print(f"  总电力消耗: {total_electricity/1e3:.1f} MW")
        
    except Exception as e:
        print(f"获取实际成本数据时出错: {str(e)}")
    
    return {
        'Installed\nequipment\ncost': equipment_costs,
        'Cooling\nduty': cooling_duty,
        'Heating\nduty': heating_duty,
        'Electricity\nconsumption': electricity,
        'Operating\ncost': operating_cost
    }


def create_plots_from_monte_carlo(save_path=None):
    """Create and save publication-quality figures for microalgae MCCA system using Monte Carlo results."""

    
    try:
        bst.main_flowsheet.clear()
        bst.main_flowsheet.set_flowsheet('_')  
        
        from microalgae._chemicals import chems
        tmo.settings.set_thermo(chems)
    except:
        pass
    
    # 检查results目录是否存在，如果不存在则创建
    os.makedirs('results', exist_ok=True)
    
    # 定义结果文件路径
    mc_file_path = 'results/monte_carlo_results.csv'
    
    # 移除所有缓存检查，直接运行新模拟
    print("运行新的Monte Carlo模拟...")
    # 初始化系统
    microalgae_sys = create_microalgae_MCCA_production_sys()
    microalgae_sys.simulate()
    
    # 创建TEA
    u = microalgae_sys.flowsheet.unit
    s = microalgae_sys.flowsheet.stream
    dry_tpd = u.U101.ins[0].F_mass * 24 / 1000  # kg/h -> t/d
    microalgae_tea = create_tea(
        system=microalgae_sys, 
        IRR=0.10, 
        duration=(2024, 2045),
        depreciation='MACRS7', 
        income_tax=0.21, 
        operating_days=330,
        lang_factor=None, 
        construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, 
        startup_FOCfrac=1, 
        startup_salesfrac=0.5,
        startup_VOCfrac=0.75, 
        WC_over_FCI=0.05,
        finance_interest=0.08, 
        finance_years=10, 
        finance_fraction=0.4,
        OSBL_units=(u.CT, u.CWP, u.ADP, u.PWC),
        warehouse=0.04, 
        site_development=0.09, 
        additional_piping=0.045,
        proratable_costs=0.10, 
        field_expenses=0.10, 
        construction=0.20,
        contingency=0.10, 
        other_indirect_costs=0.10, 
        labor_cost=max(0.5e6, compute_labor_cost(dry_tpd)),
        labor_burden=0.90, 
        property_insurance=0.007, 
        maintenance=0.03, 
        boiler_turbogenerator=None,
        steam_power_depreciation='MACRS20'
    )
    
    # 设置主产品和相关变量，用于后续LCA计算
    main_product = s.caproic_acid_product
    main_product_chemical_IDs = ['CaproicAcid']
    
    # 尝试找到锅炉单元
    boiler = None
    for unit in microalgae_sys.units:
        if hasattr(unit, 'heat_utilities') and unit.heat_utilities:
            boiler = unit
            break
    
    # 如果没有找到锅炉，使用热交换器
    if boiler is None:
        for unit in microalgae_sys.units:
            if unit.__class__.__name__ == 'HXutility':
                boiler = unit
                break
    
    # 创建评估模型
    model = Model(microalgae_sys)
    
    # 定义setter函数
    def set_microalgae_price(price):
        u.U101.ins[0].price = price
        
    def set_caproic_acid_price(price):
        s.caproic_acid_product.price = price
        
    def set_caproic_acid_yield_factor(factor):
        u.R301.caproic_acid_yield_factor = factor
    
    def set_operating_days(days):
        microalgae_tea.operating_days = days
        
    def set_maintenance_factor(factor):
        microalgae_tea.maintenance = factor
    
    # 定义指标
    model.metrics = [
        Metric('ROI', lambda: microalgae_tea.ROI * 100, '%'),
        Metric('MSP', lambda: microalgae_tea.solve_price(s.caproic_acid_product), '$/kg'),
        Metric('TCI', lambda: microalgae_tea.TCI, '$'),
        Metric('VOC', lambda: microalgae_tea.VOC, '$'),
        Metric('FOC', lambda: microalgae_tea.FOC, '$'),
        Metric('GWP', lambda: create_microalgae_lca(system=microalgae_sys,
                                                   main_product=main_product,
                                                   main_product_chemical_IDs=main_product_chemical_IDs,
                                                   boiler=boiler).GWP, 'kg CO2 eq./kg')
    ]
    
    # 注册参数
    model.parameter(setter=set_microalgae_price, element=u.U101.ins[0], 
                  name='Microalgae Price', kind='isolated', 
                  units='$/kg', distribution=shape.Triangle(-0.2, -0.1, 0.1))
    
    model.parameter(setter=set_caproic_acid_price, element=s.caproic_acid_product,
                  name='Caproic Acid Price', kind='isolated',
                  units='$/kg', distribution=shape.Triangle(2.0, 2.89, 4.0))
    
    model.parameter(setter=set_caproic_acid_yield_factor, element=u.R301,
                  name='Caproic Acid Yield Factor', kind='isolated',
                  distribution=shape.Triangle(0.7, 1.0, 1.3))
    
    model.parameter(setter=set_operating_days, element=microalgae_tea,
                  name='Operating Days', kind='isolated',
                  distribution=shape.Triangle(300, 330, 350))
    
    model.parameter(setter=set_maintenance_factor, element=microalgae_tea,
                  name='Maintenance Factor', kind='isolated',
                  distribution=shape.Triangle(0.02, 0.03, 0.05))
    
    # 运行Monte Carlo模拟
    n_samples = 10
    
    # 获取参数分布
    param_distributions = {}
    for param in model.parameters:
        if hasattr(param, 'distribution'):
            param_distributions[param.index] = param.distribution
    
    # 获取指标名称映射
    metric_names = {}
    for i, metric in enumerate(model.metrics):
        metric_names[i] = metric.name
    
    # 获取基准样本作为模板
    baseline = model.get_baseline_sample()
    results = []
    
    # 生成样本并评估
    for i in range(n_samples):
        if i % 100 == 0 and i > 0:
            print(f"Completed {i}/{n_samples} samples")
            
        # 复制基准样本
        sample = baseline.copy()
        
        # 为每个参数生成随机值
        for param_name, distribution in param_distributions.items():
            # 简化的拉丁超立方抽样
            u_val = (i + 0.5) / n_samples
            if hasattr(distribution, 'ppf'):
                sample[param_name] = distribution.ppf(u_val)
            else:
                sample[param_name] = distribution.sample()
        
        # 评估样本
        try:
            result_series = model(sample)
            
            # 创建包含参数和指标值的字典
            result_dict = {}
            
            # 添加参数值
            for param_name in param_distributions.keys():
                result_dict[param_name] = sample[param_name]
            
            # 添加指标值
            for i, name in metric_names.items():
                result_dict[name] = result_series.iloc[i]

            # 保障 GWP 指标存在（某些情况下 Metric 计算可能失败）
            if 'GWP' not in result_dict or np.isnan(result_dict['GWP']):
                try:
                    lca_tmp = create_microalgae_lca(system=microalgae_sys,
                                                   main_product=main_product,
                                                   main_product_chemical_IDs=main_product_chemical_IDs,
                                                   boiler=boiler)
                    result_dict['GWP'] = lca_tmp.GWP
                except Exception as _:
                    result_dict['GWP'] = np.nan
            
            results.append(result_dict)
        except Exception as e:
            print(f"Sample {i} evaluation failed: {str(e)}")
    
    # 组合结果
    monte_carlo_results = pd.DataFrame(results)
    
    # 保存结果
    monte_carlo_results.to_csv(mc_file_path, index=False)
    
    # =========================================================================
    # FIGURE A: MFSP comparison across different scales (Monte Carlo based)
    # =========================================================================
    fig = plt.figure(figsize=(7, 8))
    
    # A: MFSP across different scales 
    ax1 = plt.subplot2grid((3, 2), (0, 0))

    # 使用 Monte Carlo 样本结果绘制 MFSP 分布 (均值 ± 标准差)
    if 'MSP' in monte_carlo_results.columns:
        mfsp_series = monte_carlo_results['MSP'].dropna()
        mean_mfsp = mfsp_series.mean()
        std_mfsp = mfsp_series.std()
        ax1.bar(['Monte Carlo'], [mean_mfsp], yerr=[std_mfsp], color='#A52A2A', alpha=0.7, width=0.5)
        # 价格范围背景
        ax1.axhspan(2.5, 3.5, alpha=0.3, color='gray', label='Bio-based MFSP range')
        ax1.axhspan(1.0, 2.5, alpha=0.2, color='gray', label='Market price range')
        ax1.set_ylim(0, max(5.0, mean_mfsp*1.5))
    else:
        ax1.text(0.5, 0.5, 'No MSP data', ha='center', va='center')

    ax1.set_ylabel('MFSP [$kg^{-1}$]')
    ax1.set_title('MFSP (Monte Carlo mean ± std)')
    ax1.set_xticks([])

    legend_elements = [
        plt.Rectangle((0,0), 1, 1, fc='gray', alpha=0.3, label='Bio-based MFSP range'),
        plt.Rectangle((0,0), 1, 1, fc='gray', alpha=0.2, label='Market price range')
    ]
    ax1.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=7)

    # Cost and utility breakdown chart - 使用计算数据
    ax2 = plt.subplot2grid((3, 2), (0, 1))

    # 通过计算获取实际的成本和公用设施数据
    # 获取分类数据
    breakdown_data = calculate_breakdown_data(microalgae_tea, microalgae_sys)

    # 定义类别
    categories = ['Installed\nequipment\ncost', 'Cooling\nduty', 'Heating\nduty', 'Electricity\nconsumption', 'Operating\ncost']

    # 创建数据字典
    data = {}
    for key in list(breakdown_data.values())[0].keys():  # 使用第一个字典的键
        data[key] = []
        for category in categories:
            data[key].append(breakdown_data[category][key])

    # 创建DataFrame用于绘图
    df = pd.DataFrame(data, index=categories)

    # 保存成本与公用设施百分比分布到 CSV
    try:
        df.to_csv(os.path.join('results', 'cost_utility_breakdown_percent.csv'))
    except Exception as _:
        pass

    # 绘制堆叠条形图
    bottom = np.zeros(len(categories))
    colors = ['#FADBD8', '#A52A2A', '#641E16', '#FFFFFF', '#EAECEE', '#D5D8DC', 
              '#85929E', '#2C3E50', '#1C2833', '#17202A', '#7D3C98', '#FFFFFF']
    hatches = ['', '', '', '', '', '', '', '', '', '', '', '////']

    for i, (col, color, hatch) in enumerate(zip(df.columns, colors, hatches)):
        ax2.bar(categories, df[col], bottom=bottom, label=col, color=color, hatch=hatch)
        bottom += df[col]

    # 计算每个类别的总和和单位 - 使用实际数据
    try:
        # 从TEA获取实际数值
        actual_TCI = microalgae_tea.TCI / 1e6  # M$
        actual_heating = sum(abs(hu.duty) for unit in microalgae_sys.units 
                           if hasattr(unit, 'heat_utilities') 
                           for hu in unit.heat_utilities if hu.duty > 0) / 1e6  # MW
        actual_cooling = sum(abs(hu.duty) for unit in microalgae_sys.units 
                           if hasattr(unit, 'heat_utilities') 
                           for hu in unit.heat_utilities if hu.duty < 0) / 1e6  # MW
        actual_electricity = sum(abs(unit.power_utility.power) for unit in microalgae_sys.units 
                               if hasattr(unit, 'power_utility') and unit.power_utility) / 1e3  # MW
        actual_operating = (microalgae_tea.VOC + microalgae_tea.FOC) / 1e6  # M$/y
        
        sum_values = [actual_TCI, actual_cooling, actual_heating, actual_electricity, actual_operating]
    except Exception as e:
        print(f"[WARNING] Failed to fetch aggregate totals: {e}")
        sum_values = [np.nan]*5
    
    units = ['[$10^6$ $]', '[MW]', '[MW]', '[MW]', '[$10^6$ $/y]']

    # 添加总和值和单位
    for i, category in enumerate(categories):
        if i < len(sum_values) and not np.isnan(sum_values[i]):
            plt.text(i, 105, f'sum: {sum_values[i]:.1f}', ha='center', fontsize=6)
            plt.text(i, 112, units[i], ha='center', fontsize=6)
    
    ax2.set_ylim(-40, 120)
    ax2.set_ylabel('Cost and Utility Breakdown [%]')
    # 将图例移到图右侧
    ax2.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), fontsize=5)
    
    # =========================================================================
    # FIGURE B: GWP comparison and breakdown (使用LCA计算)
    # =========================================================================
    
    # B: GWP across different scales - 使用LCA计算的实际数据
    ax3 = plt.subplot2grid((3, 2), (1, 0))

    # 使用与MFSP相同的条件计算GWP值
    # 这里复用之前设置的conditions
    if 'GWP' in monte_carlo_results.columns:
        gwp_series = monte_carlo_results['GWP'].dropna()
        mean_gwp = gwp_series.mean()
        std_gwp = gwp_series.std()
        ax3.bar(['Monte Carlo'], [mean_gwp], yerr=[std_gwp], color='#A52A2A', alpha=0.7, width=0.5)
        ax3.set_ylim(0, max(5.0, mean_gwp*1.3))
    else:
        ax3.text(0.5,0.5,'No GWP data',ha='center',va='center')
        ax3.set_ylim(0,5)

    # 添加参考线 - 行业基准值
    f1_value = 3.1  # 行业平均值示例
    f2_value = 1.7  # 行业最佳实践示例
    ax3.axhline(y=f1_value, linestyle='--', color='black', linewidth=0.5)  # f1
    ax3.axhline(y=f2_value, linestyle='--', color='black', linewidth=0.5)  # f2

    # 文字标注
    ax3.text(2.1, f1_value, '$f_1$', fontsize=8)
    ax3.text(2.1, f2_value, '$f_2$', fontsize=8)

    ax3.set_ylabel('GWP$_{100}$ [kg CO$_2$ eq. kg$^{-1}$]')

    # GWP breakdown chart - 使用LCA的实际数据
    ax4 = plt.subplot2grid((3, 2), (1, 1))

    # 从LCA获取GWP明细
    # 重新计算当前状态的LCA
    lca = create_microalgae_lca(
        system=microalgae_sys,
        main_product=main_product,
        main_product_chemical_IDs=main_product_chemical_IDs,
        boiler=boiler
    )

    # 将数据转换为百分比的字典
    gwp_total = abs(lca.GWP) if lca.GWP != 0 else 1.0  # 避免除零错误
    gwp_data = {
        'Feedstock': [max(0, lca.feedstock_GWP/gwp_total*100)],
        'Material Inputs': [max(0, lca.material_GWP/gwp_total*100)],
        'Electricity': [max(0, lca.net_electricity_GWP/gwp_total*100)],
        'Direct Emissions': [max(0, lca.direct_non_biogenic_emissions_GWP/gwp_total*100)],
    }

    # 如果有负值(碳捕获)，添加Credits类别
    if lca.feedstock_GWP < 0:
        gwp_data['Credits'] = [min(0, lca.feedstock_GWP/gwp_total*100)]  # 负值作为贷项
    else:
        gwp_data['Credits'] = [0]

    # 创建DataFrame以进行绘图
    gwp_df = pd.DataFrame(gwp_data)

    # 保存 GWP breakdown 百分比到 CSV
    try:
        gwp_df.to_csv(os.path.join('results', 'gwp_breakdown_percent.csv'), index=False)
    except Exception as _:
        pass

    # 绘制堆叠条形图 - GWP
    gwp_bottom = np.zeros(1)
    gwp_colors = ['#FADBD8', '#A52A2A', '#1C2833', '#85929E', '#FFFFFF']
    gwp_hatches = ['', '', '', '/////', '////']

    for i, (col, color, hatch) in enumerate(zip(gwp_df.columns, gwp_colors, gwp_hatches)):
        ax4.bar(['GWP$_{100}$ Breakdown'], gwp_df[col], bottom=gwp_bottom, label=col, color=color, hatch=hatch)
        gwp_bottom += gwp_df[col]
    
    ax4.set_ylim(-40, 100)
    ax4.set_ylabel('GWP$_{100}$ Breakdown [%]')
    # 将图例移到图右侧
    ax4.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), fontsize=5)
    
    # 添加脚注文本 - 基于实际LCA数据
    footnote = "*Feedstock (microalgae) farming, harvesting,\n"\
               "transportation, storage, and handling. Credit for\n"\
               f"fixed carbon: {lca.feedstock_GWP:.2f} kg CO2 eq./kg.\n"\
               "Sum of direct biogenic emissions: "\
               f"{lca.biogenic_emissions_GWP:.2f} kg CO2 eq./kg."
    plt.figtext(0.1, 0.02, footnote, fontsize=7, ha='left')
    
    # =========================================================================
    # FIGURES C-F: Contour plots for yield vs titer (使用实际计算数据)
    # =========================================================================

    # 定义产率和浓度的网格点 - 与组合图保持一致
    yields = np.linspace(30, 80, 10)  # % theoretical - 与组合图一致
    titers = np.linspace(30, 100, 10)  # g/L - 与组合图一致
    
    # 创建数据存储数组
    Z_neutral_mfsp = np.zeros((len(titers), len(yields)))
    Z_lowph_mfsp = np.zeros((len(titers), len(yields)))
    Z_neutral_gwp = np.zeros((len(titers), len(yields)))
    Z_lowph_gwp = np.zeros((len(titers), len(yields)))
    
    # 创建网格 - 与组合图一致
    Y, T = np.meshgrid(yields, titers)
    
    # 保存原始参数
    original_yield_factor = u.R301.caproic_acid_yield_factor if hasattr(u.R301, 'caproic_acid_yield_factor') else 1.0
    original_titer = u.R301.titer if hasattr(u.R301, 'titer') else 60.0

    print("开始计算等高线图数据...")
    # 计算每个网格点的MFSP和GWP值
    for i, titer in enumerate(titers):
        for j, yield_pct in enumerate(yields):
            try:
                # 改进产率因子转换公式 - 更合理地映射产率到产率因子
                yield_factor = 0.7 + (yield_pct - 30) / (80 - 30) * 0.6  # 映射到0.7-1.3范围
                
                # 设置参数
                if hasattr(u.R301, 'caproic_acid_yield_factor'):
                    u.R301.caproic_acid_yield_factor = yield_factor
                
                # 设置滴度参数
                if hasattr(u.R301, 'titer'):
                    u.R301.titer = titer
                    
                # 中性发酵条件
                microalgae_sys.simulate()
                
                # 创建LCA计算GWP
                try:
                    neutral_lca = create_microalgae_lca(
                        system=microalgae_sys,
                        main_product=main_product,
                        main_product_chemical_IDs=main_product_chemical_IDs,
                        boiler=boiler
                    )
                    neutral_gwp = neutral_lca.GWP
                except Exception as e:
                    # 如果LCA计算失败，使用简单估算模型
                    print(f"LCA计算失败：{str(e)}")
                    raise
                
                # 计算MFSP
                mfsp = microalgae_tea.solve_price(s.caproic_acid_product)
                Z_neutral_mfsp[i, j] = mfsp
                Z_neutral_gwp[i, j] = neutral_gwp
                
                # 低pH发酵条件 (假设pH对产率有5%的影响)
                if hasattr(u.R301, 'caproic_acid_yield_factor'):
                    u.R301.caproic_acid_yield_factor = yield_factor * 0.95
                
                microalgae_sys.simulate()
                
                # 创建LCA计算GWP
                try:
                    lowph_lca = create_microalgae_lca(
                        system=microalgae_sys,
                        main_product=main_product,
                        main_product_chemical_IDs=main_product_chemical_IDs,
                        boiler=boiler
                    )
                    lowph_gwp = lowph_lca.GWP
                except Exception as e:
                    # 如果LCA计算失败，使用简单估算模型
                    print(f"低pH LCA计算失败：{str(e)}")
                    raise
                    
                mfsp_lowph = microalgae_tea.solve_price(s.caproic_acid_product)
                Z_lowph_mfsp[i, j] = mfsp_lowph
                Z_lowph_gwp[i, j] = lowph_gwp
                
                # 仅显示部分计算结果以减少输出
                if (i == 0 or i == len(titers)-1) and (j == 0 or j == len(yields)-1):
                    print(f"计算完成: 浓度={titer:.1f}, 产率={yield_pct:.1f}%, 产率因子={yield_factor:.3f}, MFSP={mfsp:.2f}, GWP={neutral_gwp:.2f}")
                
            except Exception as e:
                print(f"[WARNING] Convergence failure (T={titer:.1f}, Y={yield_pct:.1f}%): {e}")
                Z_neutral_mfsp[i, j] = np.nan
                Z_lowph_mfsp[i, j] = np.nan
                Z_neutral_gwp[i, j] = np.nan
                Z_lowph_gwp[i, j] = np.nan

    # 恢复原始参数
    if hasattr(u.R301, 'caproic_acid_yield_factor'):
        u.R301.caproic_acid_yield_factor = original_yield_factor
    if hasattr(u.R301, 'titer'):
        u.R301.titer = original_titer
    microalgae_sys.simulate()

    # 添加平滑处理并在此处保存等高线数据矩阵
    try:
        from scipy.ndimage import gaussian_filter
        Z_neutral_mfsp = gaussian_filter(Z_neutral_mfsp, sigma=0.7)
        Z_lowph_mfsp = gaussian_filter(Z_lowph_mfsp, sigma=0.7)
        Z_neutral_gwp = gaussian_filter(Z_neutral_gwp, sigma=0.7)
        Z_lowph_gwp = gaussian_filter(Z_lowph_gwp, sigma=0.7)
        print("数据平滑处理已应用")
    except ImportError:
        print("无法导入scipy.ndimage，跳过数据平滑处理")

    # 保存等高线数据到 NPZ
    try:
        np.savez(os.path.join('results', 'contour_data.npz'),
                 yields=yields, titers=titers,
                 Z_neutral_mfsp=Z_neutral_mfsp, Z_lowph_mfsp=Z_lowph_mfsp,
                 Z_neutral_gwp=Z_neutral_gwp, Z_lowph_gwp=Z_lowph_gwp)
    except Exception as _:
        pass
    
    # Plot C: Neutral fermentation MFSP
    ax5 = plt.subplot2grid((3, 2), (2, 0))
    plt.title('Neutral Fermentation')
    cs = ax5.contourf(Y, T, Z_neutral_mfsp, 30, cmap='viridis_r')
    
    # 计算合适的轮廓线级别
    try:
        # 计算均匀分布的8个轮廓线级别，确保递增
        min_val = np.nanmin(Z_neutral_mfsp)
        max_val = np.nanmax(Z_neutral_mfsp)
        print(f"MFSP数据范围: {min_val:.3f} - {max_val:.3f}")
        
        if min_val >= max_val or np.isnan(min_val) or np.isnan(max_val):
            # 如果数据范围有问题，使用默认值
            levels = [1.0, 1.5, 2.0, 2.5, 3.0]
            print("使用默认级别")
        else:
            levels = np.linspace(min_val, max_val, 8)
            # 确保级别是递增的且唯一
            levels = np.unique(np.sort(levels))
            if len(levels) < 3:  # 如果唯一级别太少，使用默认值
                levels = np.linspace(min_val, max_val, 5)
            levels = np.round(levels, 3)
            
        print(f"MFSP等高线级别: {levels}")
        # 验证级别是递增的
        if len(levels) < 2 or not np.all(np.diff(levels) > 0):
            levels = [1.0, 1.5, 2.0, 2.5, 3.0]
            print("回退到默认级别")
    except Exception as e:
        print(f"级别计算错误: {e}")
        levels = [1.0, 1.5, 2.0, 2.5, 3.0]
    
    contours = ax5.contour(Y, T, Z_neutral_mfsp, levels=levels, colors='black', linewidths=0.5)
    ax5.clabel(contours, inline=True, fontsize=6, fmt='%.2f')
    ax5.set_xlabel('Yield [%theoretical]')
    ax5.set_ylabel('Titer [g L$^{-1}$]')
    
    # Add extra space on the right and place colorbar there
    fig.subplots_adjust(right=0.85)
    cbar = fig.colorbar(cs, ax=ax5, location='right', pad=0.02)
    
    # 标记基准案例点
    base_yield_pct = 50.0  # 基准产率百分比
    base_titer = original_titer  # 基准浓度
    ax5.plot(base_yield_pct, base_titer, 'D', color='white', markersize=5, markeredgecolor='black', markeredgewidth=0.5)
    
    # Adjust layout
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.4, wspace=0.3, left=0.1, right=0.9, top=0.95, bottom=0.1)
    
    # Save the figure if a path is provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    # Create additional GWP contour plots
    fig2 = plt.figure(figsize=(7, 3))
    
    # Plot E: Neutral fermentation GWP
    ax7 = plt.subplot(1, 2, 1)
    plt.title('Neutral Fermentation')
    cs = ax7.contourf(Y, T, Z_neutral_gwp, 30, cmap='viridis_r')
    
    # 计算合适的轮廓线级别
    try:
        # 计算均匀分布的8个轮廓线级别，确保递增
        min_val = np.nanmin(Z_neutral_gwp)
        max_val = np.nanmax(Z_neutral_gwp)
        print(f"GWP数据范围: {min_val:.3f} - {max_val:.3f}")
        
        if min_val >= max_val or np.isnan(min_val) or np.isnan(max_val):
            # 如果数据范围有问题，使用默认值
            levels = [1.5, 2.0, 2.5, 3.0, 3.5]
            print("使用默认级别")
        else:
            levels = np.linspace(min_val, max_val, 8)
            # 确保级别是递增的且唯一
            levels = np.unique(np.sort(levels))
            if len(levels) < 3:  # 如果唯一级别太少，使用默认值
                levels = np.linspace(min_val, max_val, 5)
            levels = np.round(levels, 3)
            
        print(f"GWP等高线级别: {levels}")
        # 验证级别是递增的
        if len(levels) < 2 or not np.all(np.diff(levels) > 0):
            levels = [1.5, 2.0, 2.5, 3.0, 3.5]
            print("回退到默认级别")
    except Exception as e:
        print(f"级别计算错误: {e}")
        levels = [1.5, 2.0, 2.5, 3.0, 3.5]
    
    contours = ax7.contour(Y, T, Z_neutral_gwp, levels=levels, colors='black', linewidths=0.5)
    ax7.clabel(contours, inline=True, fontsize=6, fmt='%.2f')
    ax7.set_xlabel('Yield [%theoretical]')
    ax7.set_ylabel('Titer [g L$^{-1}$]')
    
    # 添加基准点
    ax7.plot(base_yield_pct, base_titer, 'D', color='white', markersize=5, markeredgecolor='black', markeredgewidth=0.5)

    
    # Add colorbar outside of the two subplots
    fig2.subplots_adjust(right=0.85)
    cbar = fig2.colorbar(cs, ax=ax7, location='right', pad=0.02)
    
    # Adjust layout
    fig2.subplots_adjust(wspace=0.3, left=0.1, right=0.9)
    
    # Save the second figure if a path is provided
    if save_path:
        plt.savefig(save_path.replace('.png', '_gwp.png'), dpi=300, bbox_inches='tight')
    
        # ------------------------------------------------------------------
        # Additionally split every Axes object in fig and fig2 into its own
        # image file so that users can access each sub-plot independently.
        # ------------------------------------------------------------------
        base_dir = os.path.dirname(save_path) or '.'
        split_dir = os.path.join(base_dir, 'combined_subplots')
        os.makedirs(split_dir, exist_ok=True)

        # Helper to save each axis
        def _save_axes_individually(fig_obj, tag):
            fig_obj.canvas.draw()  # Ensure renderer exists
            renderer = fig_obj.canvas.get_renderer()
            for i, ax in enumerate(fig_obj.axes, start=1):
                # Use axis title as part of filename when available, else index
                title = ax.get_title() or f'{tag}_{i}'
                safe_title = ''.join(c if c.isalnum() else '_' for c in title)[:40]
                fname = f'{tag}_{safe_title or i}.png'
                path = os.path.join(split_dir, fname)
                bbox = ax.get_tightbbox(renderer).transformed(fig_obj.dpi_scale_trans.inverted())
                fig_obj.savefig(path, dpi=300, bbox_inches=bbox)

        _save_axes_individually(fig, 'fig1')
        _save_axes_individually(fig2, 'fig2')
    
    return fig, fig2

create_plots_from_monte_carlo('microalgae_mcca_summary.png')