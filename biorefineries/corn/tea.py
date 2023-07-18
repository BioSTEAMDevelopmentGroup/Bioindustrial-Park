# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
#                     Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from biosteam.process_tools import UnitGroup
from biorefineries import corn
from biorefineries.tea import ConventionalEthanolTEA

__all__ = ('create_tea', 'default_tea_parameters', 'tea_summary')


default_tea_parameters = dict(
    IRR=0.15,
    duration=(2018, 2038),
    depreciation='MACRS7',
    income_tax=0.35,
    operating_days=330,
    lang_factor=4,
    construction_schedule=(0.4, 0.6),
    WC_over_FCI=0.05,
    labor_cost=2.3e6,
    fringe_benefits=0.4,
    property_tax=0.001,
    property_insurance=0.005,
    supplies=0.20,
    maintenance=0.01,
    administration=0.005,
)

def create_tea(system, cls=ConventionalEthanolTEA, **tea_parameters):
    tea_parameters = tea_parameters or default_tea_parameters
    return cls(system, **tea_parameters)

def tea_summary(*args, **kwargs):
    cn = corn.Biorefinery(*args, **kwargs)
    operating_days = cn.corn_tea.operating_days
    f = operating_days * 24
    ug = UnitGroup('biorefinery', cn.corn_tea.units)
    power_utilities = ug.power_utilities
    heat_utilities = ug.heat_utilities
    heating_utilities = [i for i in heat_utilities if i.duty * i.flow > 0]
    cooling_utilities = [i for i in heat_utilities if i.duty * i.flow < 0]
    FCI = cn.corn_tea.FCI
    FOC = cn.corn_tea.FOC
    maintenance = cn.corn_tea.maintenance * FCI
    
    dct = {
        'Material costs':
            {'Corn': f * cn.corn.cost ,
             'Denaturant': f * cn.denaturant.cost,
             'Enzymes': f * (cn.alpha_amylase.cost + cn.gluco_amylase.cost),
             'Yeast': f * cn.yeast.cost,
             'Other': f * sum([i.cost for i in (cn.ammonia, cn.sulfuric_acid, cn.lime)])},
         'Production':
             {'Ethanol [MMgal / yr]': f * cn.ethanol.F_mass / 2.987 / 1e6,
              'DDGS [MT / yr]': f * cn.DDGS.F_mass / 1e3,
              'Crude oil [MT / yr]': f * cn.crude_oil.F_mass / 1e3},
         'Sales':
             {'Ethanol': f * cn.ethanol.cost,
              'DDGS': f * cn.DDGS.cost,
              'Crude oil': f * cn.crude_oil.cost},
         'Utilities':
             {'Electricity': f * sum([i.cost for i in power_utilities]),
              'Steam': f * (sum([i.cost for i in heating_utilities]) + cn.steam.cost),
              'Natural gas': f * bst.stream_utility_prices['Natural gas']*cn.natural_gas.F_mass,
              'Cooling water': f * sum([i.cost for i in cooling_utilities])},
         'Labor and supplies': 
             {'Plant operations': cn.corn_tea.labor_cost * (1. + cn.corn_tea.fringe_benefits + cn.corn_tea.supplies),
              'Maintenance': maintenance},
         'Insurance, property tax, and administration':
             FCI*(cn.corn_tea.property_tax + cn.corn_tea.property_insurance + cn.corn_tea.administration),
         
         'Depreciation': cn.corn_tea.TDC / 10.,
         'Fixed operating cost': FOC,
         'Variable operating cost': cn.corn_tea.VOC,
         'Fixed capital cost': FCI,
         
    }
    sales = dct['Sales']
    dct['Co-product credit'] = credit = - sales['DDGS'] - sales['Crude oil']
    dct['Total production cost'] = x = sum([*dct['Material costs'].values(),
                                        *dct['Utilities'].values(),
                                        dct['Depreciation'],
                                        dct['Co-product credit'],
                                        dct['Fixed operating cost']]) 
    dct['Production cost [USD/gal]'] = x / dct['Production']['Ethanol [MMgal / yr]'] / 1e6
    dct['Operating cost [USD/gal]'] = (x - dct['Depreciation']) / dct['Production']['Ethanol [MMgal / yr]'] / 1e6
    return dct
