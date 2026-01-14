#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model utilities for microalgae biorefinery uncertainty analysis

Based on succinic project's model_utils.py but adapted for microalgae system structure
"""
from pandas import DataFrame, read_excel
import chaospy as shape
import numpy as np
import biosteam as bst
from biosteam.evaluation import Model
from .system import microalgae_mcca_sys, microalgae_tea

class MicroalgaeModel(bst.Model):
    def __init__(self, system, metrics=None, specification=None, 
                 parameters=None, retry_evaluation=True, exception_hook='warn',
                 namespace_dict={}):
        Model.__init__(self, system=system, specification=specification, 
                     parameters=parameters, retry_evaluation=retry_evaluation, exception_hook=exception_hook)
        self.namespace_dict = namespace_dict
        # Set metrics after initialization
        if metrics is not None:
            self.metrics = metrics
    
    def load_parameter_distributions(self, distributions, namespace_dict=None):
        namespace_dict = namespace_dict or self.namespace_dict
            
        df = distributions
        if type(df) is not DataFrame:
            df = read_excel(distributions)
            
        create_function = self.create_function
        param = self.parameter
        
        for i, row in df.iterrows():
            name = row['Parameter name']
            element = row['Element']
            kind = row['Kind']
            units = row['Units']
            baseline = row['Baseline']
            shape_data = row['Shape']
            lower, midpoint, upper = row['Lower'], row['Midpoint'], row['Upper']
            load_statements = str(row['Load statement'])
            
            D = None
            if shape_data.lower() in ['triangular', 'triangle']:
                D = shape.Triangle(lower, midpoint, upper)
            elif shape_data.lower() in ['uniform']:
                D = shape.Uniform(lower, upper)
            
            if D is not None:
                param(name=name, 
                      setter=create_function(load_statements, namespace_dict), 
                      element=element, 
                      kind=kind, 
                      units=units,
                      baseline=baseline, 
                      distribution=D)
    
    def create_function(self, code, namespace_dict):
        def wrapper_fn(statement):
            def f(x):
                namespace_dict['x'] = x
                exec(statement, namespace_dict)
            return f
        return wrapper_fn(code) 

def create_unit_groups():
    units_dict = {unit.ID: unit for unit in microalgae_mcca_sys.units}
    unit_groups = []
    
    # Area 1: Microalgae cultivation and harvesting
    cultivation_units = [units_dict[uid] for uid in ['U101'] if uid in units_dict]
    if cultivation_units:
        unit_groups.append(bst.UnitGroup('Cultivation and harvesting', units=cultivation_units))
    
    # Area 2: Pretreatment and hydrolysis
    pretreatment_unit_ids = ['T201', 'P201', 'M201', 'P202', 'H201', 'R201', 
                            'T202', 'P203', 'R202', 'P204', 'H202', 'T203', 
                            'P205', 'T204', 'P206', 'T205', 'P207', 'S201', 
                            'M202', 'R203', 'H203', 'M203', 'R204', 'S202', 'P208']
    pretreatment_units = [units_dict[uid] for uid in pretreatment_unit_ids if uid in units_dict]
    if pretreatment_units:
        unit_groups.append(bst.UnitGroup('Pretreatment and hydrolysis', units=pretreatment_units))
    
    # Area 3: Conversion
    conversion_units = [units_dict[uid] for uid in ['H301', 'M301', 'T301', 'P301', 'R301', 'T302', 'S301'] if uid in units_dict]
    if conversion_units:
        unit_groups.append(bst.UnitGroup('Conversion', units=conversion_units))
    
    # Area 4: Separation
    separation_units = [units_dict[uid] for uid in ['M401', 'S402', 'D401', 'D402', 'D403', 'D404', 'D405'] if uid in units_dict]
    if separation_units:
        unit_groups.append(bst.UnitGroup('Separation', units=separation_units))
    
    # Area 5: Waste treatment - Anaerobic digestion
    waste_units = [units_dict[uid] for uid in ['M501', 'R501', 'M502', 'M503'] if uid in units_dict]
    if waste_units:
        unit_groups.append(bst.UnitGroup('Waste treatment and biogas', units=waste_units))
    
    # Area 6: Storage
    storage_unit_ids = ['T601', 'P601', 'T602', 'P602', 'T603', 'P603', 'T604', 'P604', 'T605', 'P605']
    storage_units = [units_dict[uid] for uid in storage_unit_ids if uid in units_dict]
    if storage_units:
        unit_groups.append(bst.UnitGroup('Storage', units=storage_units))

    # Wastewater    
    wastewater_units = [units_dict[uid] for uid in ['WastewaterT'] if uid in units_dict]
    remaining_units = [u for u in microalgae_mcca_sys.units 
                      if u not in [unit for group in unit_groups for unit in group.units]
                      and u.ID not in ['WastewaterT']]
    wastewater_units.extend(remaining_units)
    if wastewater_units:
        unit_groups.append(bst.UnitGroup('Wastewater treatment', units=wastewater_units))

    # Boiler & turbogenerator
    bt_units = [units_dict[uid] for uid in ['BT601'] if uid in units_dict]
    if bt_units:
        unit_groups.append(bst.UnitGroup('BT', units=bt_units))

    # Heat Exchanger Network
    hxn_units = [units_dict[uid] for uid in ['HXN601'] if uid in units_dict]
    if hxn_units:
        hxn_group = bst.UnitGroup('Heat exchange network', units=hxn_units)
        hxn_group.filter_savings = False
        unit_groups.append(hxn_group)

    # Other facilities
    facility_units = [units_dict[uid] for uid in ['CT', 'PWC', 'ADP', 'CWP'] if uid in units_dict]
    remaining_units = [u for u in microalgae_mcca_sys.units 
                      if u not in [unit for group in unit_groups for unit in group.units]
                      and u.ID not in ['CT', 'PWC', 'ADP', 'CWP']]
    facility_units.extend(remaining_units)
    if facility_units:
        unit_groups.append(bst.UnitGroup('Other facilities', units=facility_units))

    # Fixed Operating Costs
    foc_group = bst.UnitGroup('Fixed operating costs')
    unit_groups.append(foc_group)
    for ug in unit_groups:
        ug.autofill_metrics(shorthand=False, 
                           electricity_production=False, 
                           electricity_consumption=True,
                           material_cost=True)
    
    # Special metric configurations
    for ug in unit_groups:
        if ug.name in ['Storage', 'Other facilities'] and ug.metrics:
            for metric in ug.metrics:
                if metric.name.lower() == 'material cost':
                    metric.getter = lambda: 0.0
                    break
    
    # BT unit group
    for ug in unit_groups:
        if ug.name == 'BT' and ug.metrics:
            for i, metric in enumerate(ug.metrics):
                if metric.name.lower() == 'material cost':
                    ug.metrics[i] = bst.evaluation.Metric(
                        'Material cost', 
                        getter=lambda: microalgae_tea.utility_cost / microalgae_tea.operating_days / 24,
                        units='USD/hr',
                        element=None
                    )
                    break
    
    # Fixed Operating Costs
    for ug in unit_groups:
        if ug.name == 'Fixed operating costs' and ug.metrics:
            for i, metric in enumerate(ug.metrics):
                if metric.name.lower() == 'material cost':
                    ug.metrics[i] = bst.evaluation.Metric(
                        'Material cost', 
                        getter=lambda: microalgae_tea.FOC / microalgae_tea.operating_days / 24,
                        units='USD/hr',
                        element=None
                    )
                    break
    
    return unit_groups

def get_unit_groups():
    return create_unit_groups() 