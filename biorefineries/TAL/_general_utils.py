# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 15:36:06 2023

@author: sarangbhagwat
"""
from biosteam import UnitGroup
from biosteam.evaluation import Metric
from biosteam import HeatExchangerNetwork, DrumDryer
from flexsolve import IQ_interpolation

__all__ = {'call_all_specifications_or_run',
           'get_more_unit_groups',
           'add_metrics_to_unit_groups',
           'set_production_capacity',
           'TEA_breakdown',
           'update_facility_IDs',
           }

#%% For a given list of units, call all specifications of each unit or run each unit (in the presented order)
def call_all_specifications_or_run(units_to_run):
    if units_to_run.__class__ == list:
        for unit_to_run in units_to_run:
            if unit_to_run.specifications: 
                [i() for i in unit_to_run.specifications]
            else:
                unit_to_run._run()
    else:
        if units_to_run.specifications: 
            [i() for i in units_to_run.specifications]
        else:
            units_to_run._run()

#%% Get some more unit groups for a given system
def get_more_unit_groups(system,
                         groups_to_get=['wastewater',
                                        'storage',
                                        'boiler & turbogenerator',
                                        'cooling utility facilities',
                                        'other facilities',
                                        'heat exchanger network',]
                         ):
    unit_groups_temp = UnitGroup.group_by_area(system.units)
    u = system.flowsheet.unit
    
    def has_area_unit(units, area):
        for ui in units:
            if ui.ID[1] == str(area)[0]: return True
        return False
    
    unit_groups_ = []
    
    if 'wastewater' in groups_to_get:
        wastewater_group = [i for i in unit_groups_temp if has_area_unit(i.units, 500)][0]
        wastewater_group.name = 'wastewater'
        unit_groups_.append(wastewater_group)
        
    if 'storage' in groups_to_get:
        storage_group = [i for i in unit_groups_temp if has_area_unit(i.units, 600)][0]
        storage_group.name = 'storage'
        unit_groups_.append(storage_group)
    
    if 'boiler & turbogenerator' in groups_to_get:
        boiler_turbogenerator_group = UnitGroup('boiler & turbogenerator', 
                                                    units=(u.BT701,))
        unit_groups_.append(boiler_turbogenerator_group)
    
    if 'cooling utility facilities' in groups_to_get:
        cooling_utility_facilities_group = UnitGroup('cooling utility facilities', 
                                                         units=(u.CT801, u.CWP802,))
        unit_groups_.append(cooling_utility_facilities_group)

    if 'other facilities' in groups_to_get:
        other_facilities_group = [i for i in unit_groups_temp if has_area_unit(i.units, 900)][0]
        other_facilities_group.name = 'other facilities'
        unit_groups_.append(other_facilities_group)
    
    if 'heat exchanger network' in groups_to_get:
        heat_exchanger_network_group = UnitGroup('heat exchanger network', 
                                                         units=(u.HXN1001,))
        unit_groups_.append(heat_exchanger_network_group)
    
    unit_groups_.append(UnitGroup('natural gas (for steam generation)'))
    unit_groups_.append(UnitGroup('natural gas (for product drying)'))
    unit_groups_.append(UnitGroup('fixed operating costs'))
    
    return unit_groups_

#%% Add metrics to a given list of unit groups
def add_metrics_to_unit_groups(unit_groups, 
                               system,
                               TEA=None, 
                               BT=None,
                               hxn_class=HeatExchangerNetwork,
                               natural_gas_utilizing_non_BT_unit_classes = [DrumDryer,],
                               ):
    if not TEA:
        TEA = system.TEA
    if not BT:
        BT = system.flowsheet.unit.BT701
    
    natural_gas_utilizing_non_BT_system_units = [ui for ui in system.units
                                                 if ui.__class__ in natural_gas_utilizing_non_BT_unit_classes]
    for i in unit_groups: i.autofill_metrics(shorthand=False, 
                                             electricity_production=False, 
                                             electricity_consumption=True,
                                             material_cost=True)
    for i in unit_groups:
        if i.name == 'storage' or i.name=='other facilities' or i.name == 'cooling utility facilities':
            i.metrics[-1].getter = lambda: 0. # Material cost
        if i.name == 'cooling utility facilities':
            i.metrics[1].getter = lambda: 0. # Cooling duty
        if i.name == 'boiler & turbogenerator':
            i.filter_savings = False
            i.metrics[-1] = Metric('Material cost', 
                                                getter=lambda: TEA.utility_cost/TEA.operating_hours, 
                                                units='USD/hr',
                                                element=None) # Material cost
            # i.metrics[-2].getter = lambda: BT.power_utility.rate/1e3 # Electricity consumption [MW]
        if i.name == 'natural gas (for steam generation)':
            i.metrics[-1] = Metric('Material cost', 
                                                getter=lambda: BT.natural_gas_price * BT.natural_gas.F_mass, 
                                                units='USD/hr',
                                                element=None)
        if i.name == 'natural gas (for product drying)':
            pass
        if i.name == 'fixed operating costs':
            i.metrics[-1] = Metric('Material cost', 
                                                getter=lambda: TEA.FOC/TEA.operating_hours, 
                                                units='USD/hr',
                                                element=None)
        if i.name == 'natural gas (for product drying)':
            
            i.metrics[-1] = Metric('Material cost', 
                                                getter=lambda: sum([i.utility_cost-i.power_utility.cost\
                                                                    for i in natural_gas_utilizing_non_BT_system_units]), 
                                                units='USD/hr',
                                                element=None)
        if i.name == 'heat exchanger network':
            i.filter_savings = False
            assert isinstance(i.units[0], hxn_class)

#%% Set production capacity by adjusting feedstock capacity
def set_production_capacity(
                            desired_annual_production, # pure metric tonne /y
                            system,
                            product_stream, 
                            product_chemical_IDs,
                            feedstock_stream,
                            feedstock_F_mass_range, # wet-kg/h
                            method='analytical', # 'IQ_interpolation' or 'analytical'
                            TEA=None,
                            spec=None,
                            ):
    get_pure_product_mass = lambda: sum([product_stream.imass[i] for i in product_chemical_IDs])
    if not TEA:
        TEA = system.TEA
    
    if method=='analytical':
        system.simulate()
        feedstock_stream.F_mass *= desired_annual_production / (get_pure_product_mass() * system.TEA.operating_hours/1e3)
        if spec: spec.load_specifications(spec_1=spec.spec_1,
                                          spec_2=spec.spec_2,
                                          spec_3=spec.spec_3,)
        system.simulate()
        
    elif method=='IQ_interpolation':
        def obj_f_prod_cap(feedstock_F_mass):
            feedstock_stream.F_mass = feedstock_F_mass
            if spec: spec.load_specifications(spec_1=spec.spec_1,
                                              spec_2=spec.spec_2,
                                              spec_3=spec.spec_3,)
            system.simulate()
            return get_pure_product_mass() * system.TEA.operating_hours/1e3 - desired_annual_production
        IQ_interpolation(obj_f_prod_cap, feedstock_F_mass_range[0], feedstock_F_mass_range[1], ytol=5)
        
    else: raise RuntimeError("Invalid method for set_production_capacity: must be 'IQ_interpolation' or 'analytical'.")

#%% Get a breakdown of TEA results by unit group
def TEA_breakdown(unit_groups_dict,
                  print_output=False):
    unit_groups = list(unit_groups_dict.values())
    metric_breakdowns = {i.name: {} for i in unit_groups[0].metrics}
    for ug in unit_groups:
        for metric in ug.metrics:
            # storage_metric_val = None
            if not ug.name=='storage':
                if ug.name=='other facilities':
                    metric_breakdowns[metric.name]['storage and ' + ug.name] = metric() + unit_groups_dict['storage'].metrics[ug.metrics.index(metric)]()
                else:
                    metric_breakdowns[metric.name][ug.name] = metric()
                    
    # print metric_breakdowns
    if print_output:
        for i in unit_groups[0].metrics:
            print(f"\n\n----- {i.name} ({i.units}) -----")
            metric_breakdowns_i = metric_breakdowns[i.name]
            for j in metric_breakdowns_i.keys():
                print(f"{j}: {format(metric_breakdowns_i[j], '.3f')}")
    return metric_breakdowns

#%% Update a given system's facility IDs
def update_facility_IDs(system):
    u = system.flowsheet.unit
    u.BT701.ID = 'BT701'
    u.CT901.ID = 'CT801'
    u.CWP901.ID = 'CWP802'
    u.CIP901.ID = 'CIP901'
    u.ADP901.ID = 'ADP902'
    u.FWT901.ID = 'FWT903'
    u.PWC901.ID = 'PWC904'
