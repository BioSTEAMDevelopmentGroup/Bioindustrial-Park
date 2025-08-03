# -*- coding: utf-8 -*-
"""
Created on Sat July 20 10:50:00 2025

Microalgae biorefinery to produce medium chain fatty acids 
by anaerobic fermentation without external electron donor addition- LCA analysis

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
@version: 0.0.1
"""

from thermosteam import Stream
import numpy as np
import pandas as pd
import biosteam as bst
from biosteam.units.decorators import cost
from ._chemicals import chems

__all__ = ['LCA']

class LCA:
    """
    微藻生物炼制厂生命周期环境影响评估类。
    
    参数
    ----------
    system : System
        微藻生产系统
    CFs : dict
        格式: {<impact_category_name>: {key: float}}
        其中key可以是流名称(对于complex_feeds)、化学品ID或'Electricity'
    main_product : Stream
        主产品流
    main_product_chemical_IDs : list
        主产品化学品ID列表
    boiler : Unit
        锅炉设备单元
    by_products : list, 可选
        副产品流列表
    complex_feeds : dict, 可选
        复合进料流，格式: {key: (Stream, mass_kind)}
        其中mass_kind必须是('wet', 'dry')之一
    cooling_tower : Unit, 可选
        冷却塔设备
    chilled_water_processing_units : list, 可选
        冷冻水处理单元列表
    has_turbogenerator : bool, 可选
        是否有涡轮发电机
    functional_unit : str, 可选
        功能单位，如 '1 kg', '1 metric ton', '1 MJ (by LHV)', '1 MJ (by HHV)'
    functional_quantity_per_h_fn : callable, 可选
        计算每小时功能数量的函数
    add_EOL_GWP : bool, 可选
        是否添加产品(主产品和副产品)的寿命终止温室气体排放
    input_biogenic_carbon_streams : list, 可选
        生物来源碳输入流列表
    """
    
    # 关键影响类别
    GWP_key = 'GWP_100'
    FEC_key = 'FEC'
    WC_key = 'WC'

    def __init__(self, 
                 system, 
                 CFs, 
                 main_product, 
                 main_product_chemical_IDs, 
                 boiler, 
                 by_products=[], 
                 complex_feeds={}, 
                 cooling_tower=None, 
                 chilled_water_processing_units=[],
                 has_turbogenerator=None, 
                 functional_unit=None, 
                 functional_quantity_per_h_fn=None, 
                 add_EOL_GWP=False, 
                 input_biogenic_carbon_streams=[]):
        
        self.system = system
        system._LCA = self
        self.chemicals = chemicals = main_product.chemicals
        self.units = self.system.units
        self.flowsheet = self.system.flowsheet
        self.streams = self.system.streams
        
        self.input_biogenic_carbon_streams = input_biogenic_carbon_streams
        
        self.add_EOL_GWP = add_EOL_GWP
        
        self.complex_feeds = complex_feeds
        
        kg_product_per_h_fn = lambda: self.main_product.imass[self.main_product_chemical_IDs].sum()
        MJ_LHV_product_per_h_fn = lambda: self.main_product.LHV
        MJ_HHV_product_per_h_fn = lambda: self.main_product.HHV
        
        if functional_unit and functional_quantity_per_h_fn:
            raise AttributeError('必须提供functional_unit或functional_quantity_per_h_fn，但不能同时提供两者。')
        
        elif not functional_unit: 
            self.functional_unit = functional_unit = '1 kg'
            self.functional_quantity_per_h_fn = kg_product_per_h_fn
            
        elif not functional_quantity_per_h_fn:
            self.functional_unit = functional_unit
            if functional_unit == '1 kg':
              self.functional_quantity_per_h_fn = kg_product_per_h_fn  
            elif functional_unit == '1 metric ton':
              self.functional_quantity_per_h_fn = lambda: kg_product_per_h_fn() * 1000.
            elif functional_unit == '1 MJ (by LHV)':
              self.functional_quantity_per_h_fn = MJ_LHV_product_per_h_fn  
            elif functional_unit == '1 MJ (by HHV)':
              self.functional_quantity_per_h_fn = MJ_HHV_product_per_h_fn  

        self.CFs = CFs
        
        self._chemical_IDs = [chem.ID for chem in chemicals]
        
        self._CF_streams = _CF_streams = {}
        for impact_category in CFs.keys():
            _CF_streams[impact_category] = ic_CF_stream = Stream(f'{system.ID}_{impact_category}_CF_stream')
            for k, v in CFs[impact_category].items():
                if not k in list(complex_feeds.keys()) + ['Electricity']:
                    try: 
                        ic_CF_stream.imass[k] = v
                    except:
                        pass # 假设其他complex_feed IDs存在于CFs.keys()中
                        
        self._LCA_stream = Stream(f'{system.ID}_LCA_stream')
        
        self.main_product_chemical_IDs = main_product_chemical_IDs
        self.main_product = main_product
        self.by_products = by_products
        
        self.chem_IDs = [i.ID for i in chemicals]
        
        if has_turbogenerator is None:
            has_turbogenerator = boiler.power_utility.production > 0.
        self.has_turbogenerator = has_turbogenerator
        
        self.BT = self.boiler = boiler
        try:
            self.natural_gas = self.BT.natural_gas
        except:
            # 假设biogas替代自然气
            self.natural_gas = next((stream for stream in self.system.streams 
                                    if 'biogas' in stream.ID.lower() or 
                                       (hasattr(stream, 'price') and stream.price and
                                        stream.imass['CH4'] > 0)), None)
            if not self.natural_gas:
                # 创建一个空的流作为自然气流的替代
                self.natural_gas = Stream('natural_gas', CH4=0)
        
        self.CT = self.cooling_tower = cooling_tower
        self.CWP_units = self.chilled_water_processing_units = chilled_water_processing_units
        
        self._CO2_MW = self.chemicals.CO2.MW
    
    @property
    def system_carbon_balance(self):
        """计算系统碳平衡"""
        total_C_in = sum([feed.get_atomic_flow('C') for feed in self.feeds])
        total_C_out = self.main_product.get_atomic_flow('C') +\
                      sum([i.get_atomic_flow('C') for i in self.by_products]) +\
                      sum([emission.get_atomic_flow('C') for emission in self.emissions])
        return total_C_out/total_C_in if total_C_in > 0 else 0
    
    @property
    def products(self):
        """获取所有产品流(主产品和副产品)"""
        return [self.main_product] + self.by_products
    
    @property
    def functional_quantity_per_h(self):
        """计算每小时功能数量"""
        return self.functional_quantity_per_h_fn()
    
    @property
    def emissions(self):
        """获取排放流列表"""
        emissions = list(self.system.products)
        for i in self.products: 
            if i in emissions:
                emissions.remove(i)
        return emissions
    
    @property
    def feeds(self): 
        """获取原料流列表"""
        return self.system.feeds
        
    @property
    def LCA_stream(self):
        """创建LCA评估用的综合流"""
        _LCA_stream = self._LCA_stream
        to_mix = list(self.feeds)
        for s, m_k in self.complex_feeds.values(): 
            if s in to_mix:
                to_mix.remove(s)
        _LCA_stream.mix_from(to_mix)
        return _LCA_stream
    
    def get_material_impact_array(self, impact_category):
        """获取特定影响类别的物质影响数组"""
        return self.LCA_stream.mass*self._CF_streams[impact_category].mass
    
    def get_material_impact(self, impact_category):
        """计算特定影响类别的物质影响"""
        return self.get_material_impact_array(impact_category).sum() / self.functional_quantity_per_h

    @property
    def net_electricity(self):
        """计算净电力消耗"""
        return self.system.power_utility.rate
        
    def get_net_electricity_impact(self, impact_category):
        """计算特定影响类别的净电力影响"""
        if 'Electricity' in self.CFs[impact_category]:
            return self.net_electricity * self.CFs[impact_category]['Electricity'] / self.functional_quantity_per_h
        return 0
    
    @property
    def EOL_GWP(self): 
        """计算寿命终止全球变暖潜能值"""
        return sum([i.get_atomic_flow('C') for i in [self.main_product] + self.by_products]) *\
            self._CO2_MW/self.functional_quantity_per_h
    
    @property
    def direct_emissions_GWP(self):
        """计算直接排放全球变暖潜能值"""
        return sum([stream.get_atomic_flow('C') for stream in self.emissions]) *\
            self._CO2_MW / self.functional_quantity_per_h
            
    @property
    def biogenic_emissions_GWP(self): 
        """计算生物来源直接排放全球变暖潜能值"""
        return sum([i.get_atomic_flow('C') for i in self.input_biogenic_carbon_streams]) *\
            self._CO2_MW/self.functional_quantity_per_h
    
    @property
    def direct_non_biogenic_emissions_GWP(self): 
        """计算非生物来源直接排放全球变暖潜能值"""
        return  self.direct_emissions_GWP -\
                self.biogenic_emissions_GWP +\
                int(self.add_EOL_GWP)*self.EOL_GWP
    
    def get_complex_feeds_impact(self, impact_category):
        """计算复合进料流在特定影响类别的影响"""
        tot_cfs_impact = 0.
        impact_CFs = self.CFs[impact_category]
        for k, (s, mass_kind) in self.complex_feeds.items():
            if k in impact_CFs:
                mass = s.F_mass
                if mass_kind=='dry': mass -= s.imass['Water']
                tot_cfs_impact += impact_CFs[k] * mass
        return tot_cfs_impact/self.functional_quantity_per_h
    
    def get_total_impact(self, impact_category):
        """计算特定影响类别的总影响"""
        tot_impact = self.get_complex_feeds_impact(impact_category) +\
                     self.get_material_impact(impact_category) +\
                     self.get_net_electricity_impact(impact_category)
        if impact_category in ('GWP', 'GWP_100', 'GWP100'):
            tot_impact += self.direct_non_biogenic_emissions_GWP
        return tot_impact

    
    def get_material_impact_breakdown(self, impact_category):
        """获取特定影响类别的物质影响明细"""
        # 初始化关键化学品的影响字典
        chemical_impact_dict = {'H2SO4':0, 'NaOH':0, 'NH4OH':0, 'CH4':0, 'CO2':0, 'Octanol':0}
        LCA_stream = self.LCA_stream
        ic_CF_stream = self._CF_streams[impact_category]
        functional_quantity_per_h = self.functional_quantity_per_h
        # 添加其他化学品的影响
        chemical_impact_dict_additional = {ID: LCA_stream.imass[ID] * ic_CF_stream.imass[ID] / functional_quantity_per_h
                                for ID in self.chem_IDs}
        # 移除值为0的项
        for k in list(chemical_impact_dict_additional.keys()):
            if chemical_impact_dict_additional[k] == 0.: del(chemical_impact_dict_additional[k])
        chemical_impact_dict.update(chemical_impact_dict_additional)
        return chemical_impact_dict
    
    def get_material_impact_breakdown_as_fraction_of_material_impact(self, impact_category):
        """获取物质影响明细，表示为物质影响的比例"""
        chemical_impact_dict = self.get_material_impact_breakdown(impact_category)
        tot_material_impact = self.get_material_impact(impact_category)
        if tot_material_impact > 0:
            for k,v in chemical_impact_dict.items():
                chemical_impact_dict[k] /= tot_material_impact
        return chemical_impact_dict
    
    def get_material_impact_breakdown_as_fraction_of_total_impact(self, impact_category):
        """获取物质影响明细，表示为总影响的比例"""
        chemical_impact_dict = self.get_material_impact_breakdown(impact_category)
        tot_impact = self.get_total_impact(impact_category)
        if tot_impact > 0:
            for k,v in chemical_impact_dict.items():
                chemical_impact_dict[k] /= tot_impact
        return chemical_impact_dict
    
    def get_natural_gas_impact(self, impact_category):
        """计算自然气在特定影响类别的影响"""
        if 'CH4' in self.CFs[impact_category]:
            return self.CFs[impact_category]['CH4']*self.natural_gas.F_mass/self.functional_quantity_per_h
        return 0
    
    @property
    def natural_gas_combustion_GWP(self):
        """计算自然气燃烧的全球变暖潜能值"""
        return (self.natural_gas.get_atomic_flow('C')) * self._CO2_MW / self.functional_quantity_per_h
    
    def get_complex_feed_impact_by_ID(self, impact_category, complex_feed_ID):
        """计算特定复合进料流在特定影响类别的影响"""
        if complex_feed_ID in self.complex_feeds and complex_feed_ID in self.CFs[impact_category]:
            s, mass_kind = self.complex_feeds[complex_feed_ID]
            mass = s.F_mass
            if mass_kind=='dry': mass -= s.imass['Water']
            return self.CFs[impact_category][complex_feed_ID] * mass
        return 0

    def get_material_impact_by_ID(self, impact_category, material_ID):
        """计算特定物质在特定影响类别的影响"""
        if material_ID in self.CFs[impact_category]:
            return self.CFs[impact_category][material_ID] * self.LCA_stream.imass[material_ID]
        return 0
    
    @property
    def electricity_demand(self): 
        """计算电力需求"""
        return sum([i.power_utility.consumption for i in self.system.units])

    @property
    def cooling_electricity_demand(self):
        """计算冷却系统电力需求"""
        return (self.CT.power_utility.rate if self.CT else 0) + \
               sum([i.power_utility.rate for i in self.CWP_units])
    
    def __repr__(self):
        """返回LCA对象的字符串表示"""
        return f'LCA对象，用于{self.system.ID}。影响类别: {list(self.CFs.keys())}。'
    
    def show(self):
        """打印LCA对象信息"""
        print(self.__repr__())
    _ipython_display_ = show

    def generate_report(self, output_file=None):
        """生成LCA报告"""
        report = []
        report.append(f"======= 微藻生物炼制厂生命周期评估报告 =======")
        report.append(f"系统ID: {self.system.ID}")
        report.append(f"功能单位: {self.functional_unit}")
        report.append(f"主产品: {self.main_product.ID}")
        
        report.append("\n----- 系统碳平衡 -----")
        report.append(f"碳平衡: {self.system_carbon_balance:.4f}")
        
        report.append("\n----- 影响类别结果 -----")
        for impact_category in self.CFs.keys():
            total_impact = self.get_total_impact(impact_category)
            report.append(f"{impact_category}: {total_impact:.4f} 单位/{self.functional_unit}")
            
            report.append(f"  - 物质影响: {self.get_material_impact(impact_category):.4f} ({self.get_material_impact(impact_category)/total_impact*100:.1f}%)")
            report.append(f"  - 电力影响: {self.get_net_electricity_impact(impact_category):.4f} ({self.get_net_electricity_impact(impact_category)/total_impact*100 if total_impact else 0:.1f}%)")
            report.append(f"  - 复合进料影响: {self.get_complex_feeds_impact(impact_category):.4f} ({self.get_complex_feeds_impact(impact_category)/total_impact*100 if total_impact else 0:.1f}%)")
            
            if impact_category in ('GWP', 'GWP_100', 'GWP100'):
                report.append(f"  - 直接排放GWP: {self.direct_emissions_GWP:.4f} ({self.direct_emissions_GWP/total_impact*100 if total_impact else 0:.1f}%)")
                report.append(f"  - 生物来源排放GWP: {self.biogenic_emissions_GWP:.4f}")
                report.append(f"  - 非生物来源直接排放GWP: {self.direct_non_biogenic_emissions_GWP:.4f}")
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write('\n'.join(report))
        
        return '\n'.join(report)
    
    def get_impact_dataframe(self):
        """返回包含所有影响结果的DataFrame"""
        data = {}
        
        for impact_category in self.CFs.keys():
            data[f'{impact_category}_总影响'] = [self.get_total_impact(impact_category)]
            data[f'{impact_category}_物质影响'] = [self.get_material_impact(impact_category)]
            data[f'{impact_category}_电力影响'] = [self.get_net_electricity_impact(impact_category)]
            data[f'{impact_category}_复合进料影响'] = [self.get_complex_feeds_impact(impact_category)]
            
            if impact_category in ('GWP', 'GWP_100', 'GWP100'):
                data[f'{impact_category}_直接排放'] = [self.direct_emissions_GWP]
                data[f'{impact_category}_生物来源排放'] = [self.biogenic_emissions_GWP]
                data[f'{impact_category}_非生物来源直接排放'] = [self.direct_non_biogenic_emissions_GWP]
        
        data['功能单位'] = [self.functional_unit]
        data['系统碳平衡'] = [self.system_carbon_balance]
        
        return pd.DataFrame(data)
    
    def export_results(self, file_path):
        """导出结果到Excel文件"""
        df = self.get_impact_dataframe()
        df.to_excel(file_path)
        return f"结果已导出到 {file_path}"
    
    def analyze_emissions_breakdown(self, threshold=1e-6):
        """
        分析emissions的详细组成
        
        参数
        ----------
        threshold : float
            显示阈值，只显示流量大于此值的组分
            
        返回
        ----------
        dict
            包含emissions分析结果的字典
        """
        # 获取所有emissions流
        emissions_streams = []
        for unit in self.system.units:
            if hasattr(unit, 'outs'):
                for out_stream in unit.outs:
                    if out_stream and out_stream not in self.products:
                        emissions_streams.append(out_stream)
        
        # 去重
        unique_emissions = []
        for stream in emissions_streams:
            if stream not in unique_emissions:
                unique_emissions.append(stream)
        
        results = {
            'total_streams': len(unique_emissions),
            'streams_analysis': []
        }
        
        print("=== Emissions详细分解分析 ===")
        print(f"总emissions流数量: {len(unique_emissions)}")
        
        for i, stream in enumerate(unique_emissions):
            print(f"\n--- Emissions流 {i+1}: {stream.ID} ---")
            print(f"来源单元: {stream.source.ID if stream.source else 'Unknown'}")
            print(f"总流量: {stream.F_mass:.3f} kg/h")
            print(f"总摩尔流量: {stream.F_mol:.3f} kmol/h")
            
            # 分析组分
            components = []
            total_mass = stream.F_mass
            total_mol = stream.F_mol
            
            print("\n组分详细:")
            print("组分名称".ljust(20) + "质量流量(kg/h)".ljust(15) + "摩尔流量(kmol/h)".ljust(15) + "质量分数(%)".ljust(12) + "摩尔分数(%)".ljust(12))
            print("-" * 80)
            
            for chem in stream.chemicals:
                mass_flow = stream.imass[chem.ID]
                mol_flow = stream.imol[chem.ID]
                
                if mass_flow > threshold:  # 只显示有意义的流量
                    mass_fraction = (mass_flow / total_mass * 100) if total_mass > 0 else 0
                    mol_fraction = (mol_flow / total_mol * 100) if total_mol > 0 else 0
                    
                    print(f"{chem.ID:<20} {mass_flow:>14.3f} {mol_flow:>14.3f} {mass_fraction:>11.2f} {mol_fraction:>11.2f}")
                    
                    components.append({
                        'chemical': chem.ID,
                        'mass_flow': mass_flow,
                        'mol_flow': mol_flow,
                        'mass_fraction': mass_fraction,
                        'mol_fraction': mol_fraction
                    })
            
            # 碳含量分析
            total_carbon = 0
            carbon_components = []
            
            for chem in stream.chemicals:
                if hasattr(chem, 'formula') and chem.formula:
                    carbon_atoms = 0
                    formula = chem.formula
                    if 'C' in formula:
                        c_index = formula.find('C')
                        if c_index + 1 < len(formula) and formula[c_index + 1].isdigit():
                            num_str = ""
                            i = c_index + 1
                            while i < len(formula) and formula[i].isdigit():
                                num_str += formula[i]
                                i += 1
                            carbon_atoms = int(num_str)
                        else:
                            carbon_atoms = 1
                    
                    if carbon_atoms > 0:
                        mol_flow = stream.imol[chem.ID]
                        carbon_mol = mol_flow * carbon_atoms
                        if carbon_mol > threshold:
                            total_carbon += carbon_mol
                            carbon_components.append({
                                'chemical': chem.ID,
                                'mol_flow': mol_flow,
                                'carbon_atoms': carbon_atoms,
                                'carbon_mol': carbon_mol
                            })
            
            if carbon_components:
                print(f"\n含碳组分分析:")
                print("组分名称".ljust(20) + "摩尔流量(kmol/h)".ljust(15) + "碳原子数".ljust(10) + "碳摩尔流量(kmol/h)".ljust(15))
                print("-" * 65)
                
                for comp in carbon_components:
                    print(f"{comp['chemical']:<20} {comp['mol_flow']:>14.3f} {comp['carbon_atoms']:>9} {comp['carbon_mol']:>14.3f}")
                
                print(f"\n总碳摩尔流量: {total_carbon:.3f} kmol/h")
                print(f"总碳质量流量: {total_carbon * 12.01:.3f} kg/h")
            
            # CO2分析
            if 'CO2' in stream.chemicals:
                co2_mol = stream.imol['CO2']
                co2_mass = stream.imass['CO2']
                if co2_mol > threshold:
                    co2_mass_fraction = (co2_mass / total_mass * 100) if total_mass > 0 else 0
                    co2_mol_fraction = (co2_mol / total_mol * 100) if total_mol > 0 else 0
                    
                    print(f"\nCO2分析:")
                    print(f"CO2摩尔流量: {co2_mol:.3f} kmol/h")
                    print(f"CO2质量流量: {co2_mass:.3f} kg/h")
                    print(f"CO2质量分数: {co2_mass_fraction:.2f}%")
                    print(f"CO2摩尔分数: {co2_mol_fraction:.2f}%")
                    print(f"CO2当量质量: {co2_mol * 44.01:.3f} kg/h")
            
            # 其他重要气体
            important_gases = ['H2O', 'N2', 'O2', 'SO2', 'NO', 'NO2', 'CH4', 'N2O']
            found_gases = []
            
            for gas in important_gases:
                if gas in stream.chemicals:
                    mol_flow = stream.imol[gas]
                    mass_flow = stream.imass[gas]
                    if mol_flow > threshold:
                        mass_fraction = (mass_flow / total_mass * 100) if total_mass > 0 else 0
                        mol_fraction = (mol_flow / total_mol * 100) if total_mol > 0 else 0
                        found_gases.append({
                            'gas': gas,
                            'mol_flow': mol_flow,
                            'mass_flow': mass_flow,
                            'mass_fraction': mass_fraction,
                            'mol_fraction': mol_fraction
                        })
            
            if found_gases:
                print(f"\n其他重要气体:")
                for gas_info in found_gases:
                    print(f"{gas_info['gas']}: {gas_info['mol_flow']:.3f} kmol/h, {gas_info['mass_flow']:.3f} kg/h, {gas_info['mass_fraction']:.2f}% mass, {gas_info['mol_fraction']:.2f}% mol")
            
            # 热力学性质
            print(f"\n热力学性质:")
            print(f"温度: {stream.T:.2f} K ({stream.T - 273.15:.2f} °C)")
            print(f"压力: {stream.P:.0f} Pa ({stream.P/101325:.3f} atm)")
            print(f"焓: {stream.H:.3f} kJ/h")
            
            # 保存分析结果
            stream_analysis = {
                'stream_id': stream.ID,
                'source_unit': stream.source.ID if stream.source else 'Unknown',
                'total_mass_flow': stream.F_mass,
                'total_mol_flow': stream.F_mol,
                'components': components,
                'carbon_analysis': {
                    'total_carbon_mol': total_carbon,
                    'total_carbon_mass': total_carbon * 12.01,
                    'carbon_components': carbon_components
                },
                'thermodynamic_properties': {
                    'temperature': stream.T,
                    'pressure': stream.P,
                    'enthalpy': stream.H
                }
            }
            
            results['streams_analysis'].append(stream_analysis)
        
        return results
    
    # 便捷属性 - GWP相关
    @property
    def material_GWP(self): 
        return self.get_material_impact(self.GWP_key)
 
    @property
    def feedstock_GWP(self): 
        return self.get_complex_feeds_impact(self.GWP_key)
     
    @property
    def net_electricity_GWP(self): 
        return self.get_net_electricity_impact(self.GWP_key)
     
    @property
    def GWP(self): 
        return self.get_total_impact(self.GWP_key)
     
    # 便捷属性 - FEC相关
    @property
    def material_FEC(self): 
        return self.get_material_impact(self.FEC_key)
     
    @property
    def feedstock_FEC(self):
        return self.get_complex_feeds_impact(self.FEC_key)
     
    @property
    def net_electricity_FEC(self): 
        return self.get_net_electricity_impact(self.FEC_key)
     
    @property
    def FEC(self): 
        return self.get_total_impact(self.FEC_key)
         
    # 便捷属性 - WC相关
    @property
    def material_WC(self): 
        return self.get_material_impact(self.WC_key)
     
    @property
    def feedstock_WC(self):
        return self.get_complex_feeds_impact(self.WC_key)
     
    @property
    def net_electricity_WC(self): 
        return self.get_net_electricity_impact(self.WC_key)
     
    @property
    def WC(self): 
        return self.get_total_impact(self.WC_key)


# 使用示例
def create_microalgae_lca(system, main_product, main_product_chemical_IDs, boiler):
    """
    创建微藻生物炼制厂的LCA对象
    
    参数
    ----------
    system : System
        微藻生产系统
    main_product : Stream
        主产品流
    main_product_chemical_IDs : list
        主产品化学品ID列表
    boiler : Unit
        锅炉设备单元
        
    返回
    ----------
    LCA
        生命周期评估对象
    """
    # 影响因子
    # 示例数据，实际应用中需要替换为准确的数据
    CFs = {
        'GWP_100':{
            'Electricity': 0.36, # [kg*CO2*eq / kWhr] From GREET; NG-Fired Simple-Cycle Gas Turbine CHP Plant
            # 0.66 is the GWP from producing diesel from GREET; Conventional diesel from crude oil for US Refineries.
            # Downstream fuel emissions are added in. Accounts for how biodiesel has less energy than diesel.
            'CH4': 0.33,  # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
            
            'H2SO4': 0.04447,  # kg CO2e/kg
            'NaOH': 2.01,  # GREET
            'NH4OH': 1.28304, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,   
            'CalciumDihydroxide': 1.29,  # GREET
            'Ethanol': 1.44, # from BDO project
            
            'GlucoAmylase': 6.16,  # from Succinic project
            'AlphaAmylase': 6.16,  # from Succinic project
            
            'Octanol': 1.9,  # kg CO2e/kg

            'Microalgae': 0.15,  # 
        },
        'FEC': {
            'Electricity': 10.0,  # MJ/kWh
            'H2SO4': 568.98/1e3,
            'NaOH': 29,  # MJ/kg
            'NH4OH': 42 * 0.4860, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,
            'CalciumDihydroxide': 4.0,  # MJ/kg
            'CH4': 55.0,  # MJ/kg
            'Ethanol': 16,
            'Octanol': 45.0,  # MJ/kg
            'GlucoAmylase': 90.0,  # MJ/kg
            'AlphaAmylase': 90.0,  # MJ/kg
            'Microalgae': 2.0,  # MJ/kg
        },
        'WC': {  # 水消耗
            'Electricity': 0.005, # m3 / kWh
            'H2SO4': 0.01,        # m3 / kg
            'NaOH': 0.02,         # m3 / kg
            'NH4OH': 0.015,       # m3 / kg
            'CalciumDihydroxide': 0.007, # m3 / kg
            'CH4': 0.0,           # m3 / kg
            'CO2': 0.0,           # m3 / kg
            'Octanol': 0.05,      # m3 / kg
            'GlucoAmylase': 0.1,  # m3 / kg
            'AlphaAmylase': 0.1,  # m3 / kg
            'Microalgae': 0.2,    # m3 / kg
        }
    }
    
    # 获取系统中的冷却塔和冷冻水处理单元
    cooling_tower = None
    chilled_water_processing_units = []
    
    for unit in system.units:
        if unit.__class__.__name__ == 'CoolingTower':
            cooling_tower = unit
        elif unit.__class__.__name__ == 'ChilledWaterPackage':
            chilled_water_processing_units.append(unit)
    
    # 尝试找到系统中的锅炉涡轮发电机
    if boiler is None:
        for unit in system.units:
            # 优先查找锅炉涡轮发电机
            if (isinstance(unit, bst.facilities.BoilerTurbogenerator) or
                ('BT' in unit.ID and hasattr(unit, 'turbogenerator_efficiency'))):
                boiler = unit
                break
        
        # 如果没有找到锅炉涡轮发电机，查找普通锅炉
        if boiler is None:
            for unit in system.units:
                if 'boiler' in unit.__class__.__name__.lower() or 'BT' in unit.ID:
                    boiler = unit
                    break
        
        # 如果没有找到锅炉，使用最大加热负荷的热交换器
        if boiler is None:
            heat_exchangers = []
            for unit in system.units:
                if unit.__class__.__name__ == 'HXutility' and hasattr(unit, 'duty') and unit.duty > 0:
                    heat_exchangers.append(unit)
            
            if heat_exchangers:
                # 选择具有最大加热负荷的热交换器
                boiler = max(heat_exchangers, key=lambda x: x.duty)
    
    # 创建复合进料字典
    complex_feeds = {
        'Microalgae': (system.feeds[0], 'dry'),  # 微藻进料
    }
    
    # 确定生物来源碳输入流
    # 找到微藻原料流
    microalgae_feed = None
    for feed in system.feeds:
        if 'microalgae' in feed.ID.lower():
            microalgae_feed = feed
            break
    
    # 如果没有找到microalgae流，使用第一个包含碳的原料流
    if microalgae_feed is None:
        for feed in system.feeds:
            if feed.get_atomic_flow('C') > 0:
                microalgae_feed = feed
                break
    
    input_biogenic_carbon_streams = [microalgae_feed] if microalgae_feed else []
    
    # 创建LCA对象
    microalgae_lca = LCA(
        system=system,
        CFs=CFs,
        main_product=main_product,
        main_product_chemical_IDs=main_product_chemical_IDs,
        boiler=boiler,
        by_products=[system.flowsheet.stream.butyric_acid_product,
                     system.flowsheet.stream.heptanoic_acid_product, 
                     system.flowsheet.stream.caprylic_acid_product], 
        complex_feeds=complex_feeds,
        cooling_tower=cooling_tower,
        chilled_water_processing_units=chilled_water_processing_units,
        has_turbogenerator=None,  # 自动确定    
        functional_unit='1 kg',  # 每千克产品
        add_EOL_GWP=True,
        input_biogenic_carbon_streams=input_biogenic_carbon_streams
    )
    
    return microalgae_lca

def print_non_biogenic_breakdown(lca, threshold=1e-4):
    print('\n-- 非生物源直接排放明细 (kg CO2-eq / 功能单位) --')
    tot = 0
    for s in lca.emissions:
        if s in lca.input_biogenic_carbon_streams:
            # 生物源流 => 扣除，不计入非生物
            continue
        C = s.get_atomic_flow('C')            # kmol/h
        if C == 0:
            continue
        impact = C * 44 / lca.functional_quantity_per_h
        if abs(impact) > threshold:
            print(f'{s.ID:<15} {impact:8.4f}')
        tot += impact
    if lca.add_EOL_GWP:
        print(f'End-of-Life      {lca.EOL_GWP:8.4f}')
        tot += lca.EOL_GWP
    print('-----------------------------------------')
    print(f'Total            {tot:8.4f}')


# 主程序运行示例
if __name__ == "__main__":
    import biosteam as bst
    from microalgae import system
    from microalgae._chemicals import chems
    
    # 初始化系统
    bst.settings.set_thermo(chems)
    microalgae_sys = system.create_microalgae_MCCA_production_sys()
    microalgae_sys.simulate()
    
    # 获取主产品流和主产品化学品ID
    main_product = microalgae_sys.flowsheet.stream.caproic_acid_product
    main_product_chemical_IDs = ['CaproicAcid']
    
    # 尝试找到锅炉单元
    boiler = None
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
    lca = create_microalgae_lca(
        system=microalgae_sys,
        main_product=main_product,
        main_product_chemical_IDs=main_product_chemical_IDs,
        boiler=boiler
    )
    
    # 输出LCA结果
    print("\n===== 微藻生物炼制厂LCA结果 =====")
    print(f"系统碳平衡: {lca.system_carbon_balance:.4f}")
    print(f"GWP (kg CO2e/{lca.functional_unit}): {lca.GWP:.4f}")
    print(f"FEC (MJ/{lca.functional_unit}): {lca.FEC:.4f}")
    print(f"WC (m³/{lca.functional_unit}): {lca.WC:.4f}")
    
    # 输出GWP明细
    print(f"----- GWP明细 -----")
    print(f"物质影响: {lca.material_GWP:.4f} ({lca.material_GWP/lca.GWP*100 if lca.GWP else 0:.1f}%)")
    print(f"电力影响: {lca.net_electricity_GWP:.4f} ({lca.net_electricity_GWP/lca.GWP*100 if lca.GWP else 0:.1f}%)")
    print(f"原料影响: {lca.feedstock_GWP:.4f} ({lca.feedstock_GWP/lca.GWP*100 if lca.GWP else 0:.1f}%)")
    print(f"直接排放(非生物源): {lca.direct_non_biogenic_emissions_GWP:.4f} ({lca.direct_non_biogenic_emissions_GWP/lca.GWP*100 if lca.GWP else 0:.1f}%)")
    
    # 分析emissions分解
    #print(f"\n{'='*60}")
    #lca.analyze_emissions_breakdown()
    
    # 生成并保存详细报告
    report = lca.generate_report()
    print("\n完整报告已生成。报表示例:")
    print(report[:500] + "...\n")
    
    # 导出结果到Excel
    # lca.export_results("microalgae_lca_results.xlsx")

    # 1) 物质影响 – 按单一化学品列出
    mat_break = lca.get_material_impact_breakdown('GWP_100')
    print('\n-- Material impact breakdown (kg CO2-eq / kg product) --')
    for k, v in mat_break.items():
        if abs(v) > 1e-6:         # 过滤 0
            print(f'{k:<15} {v:8.4f}')

    # 2) 电力影响 – 直接由净电力×CF 得到
    elec_CF   = lca.CFs['GWP_100']['Electricity']
    net_kWh   = lca.net_electricity              # kWh h⁻¹
    per_kg    = net_kWh * elec_CF / lca.functional_quantity_per_h
    print(f'\nElectricity ({net_kWh:.3f} kWh/h) × {elec_CF} = {per_kg:.4f} kg CO2-eq')

    # 3) 复合进料（feedstock）影响
    feed_imp  = lca.get_complex_feeds_impact('GWP_100')
    print(f'\nFeedstock impact = {feed_imp:.4f} kg CO2-eq')

    # 若想看各复合进料流的贡献:
    for ID, (s, _) in lca.complex_feeds.items():
        impact = lca.get_complex_feed_impact_by_ID('GWP_100', ID)
        print(f'  {ID:<15} {impact:8.4f}')

    # 4) 直接排放
    print(f'\nDirect non-biogenic emissions GWP = {lca.direct_non_biogenic_emissions_GWP:.4f}')
    print('\n  排放流详细:')
    for s in lca.emissions:
        C = s.get_atomic_flow("C")
        if C>0:
            stream_gwp = C * 44 / lca.functional_quantity_per_h   # 44 g CO2 / mol C
            if stream_gwp!=0:
                print(f'  {s.ID:<15}  {stream_gwp:8.4f} kg CO2-eq')

    # 若已把生物源排放列入 `input_biogenic_carbon_streams`,
    # 也可打印 biogenic 值:
    print(f'\nBiogenic CO2 credit = {lca.biogenic_emissions_GWP:.4f}')
    print_non_biogenic_breakdown(lca)

    # u = microalgae_sys.flowsheet.unit   
    # s = microalgae_sys.flowsheet.stream  

    # s.s49.show()    
    # print('From:', s.s49.source)   
    # print('To :', s.s49.sink)      

    # s.s51.show()
    # u.R501.show()
    # print('From :', s.s51.source)
    # print('To :', s.s51.sink)




