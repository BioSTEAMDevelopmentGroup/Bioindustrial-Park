# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biosteam.process_tools import UnitGroup
from biorefineries import corn as cn
from biorefineries.sugarcane import ConventionalEthanolTEA

__all__ = ('create_tea', 'tea_summary')

def create_tea(system, cls=ConventionalEthanolTEA):
    return cls(system, IRR=0.15,
               duration=(2018, 2038),
               depreciation='MACRS10', income_tax=0.35,
               operating_days=330, lang_factor=3,
               construction_schedule=(0.4, 0.6), WC_over_FCI=0.05,
               labor_cost=1e6, fringe_benefits=0.4,
               property_tax=0.001, property_insurance=0.005,
               supplies=0.20, maintenance=0.01, administration=0.005)

def tea_summary():
    operating_days = cn.corn_tea.operating_days
    f = operating_days * 24
    ug = UnitGroup('biorefinery', cn.corn_tea.units)
    power_utilities = ug.power_utilities
    heat_utilities = ug.heat_utilities
    heating_utilities = [i for i in heat_utilities if i.duty > 0]
    cooling_utilities = [i for i in heat_utilities if i.duty < 0]
    dct = {
        'Material costs':
            {'Corn': f * cn.corn.cost ,
             'Denaturant': f * cn.denaturant.cost,
             'Enzymes': f * (cn.alpha_amylase.cost + cn.gluco_amylase.cost),
             'Yeast': f * cn.yeast.cost,
             'Other': f * sum([i.cost for i in (cn.ammonia, cn.sulfuric_acid, cn.lime)])},
         'Utilities':
             {'Electricity':f * sum([i.cost for i in power_utilities]),
              'Steam': f * sum([i.cost for i in heating_utilities]),
              'Natural gas': f * cn.D610.natural_gas_cost},
              'Cooling water': f * sum([i.cost for i in cooling_utilities]),
         'Depreciation': cn.corn_tea.TDC / 10.,
         'Fixed operating costs': cn.corn_tea.VOC,
    }
    dct['Total production cost'] = sum([*dct['Material costs'].values(),
                                        *dct['Utilities'].values(),
                                        dct['Depreciation'],
                                        dct['Fixed operating costs']])
    return dct