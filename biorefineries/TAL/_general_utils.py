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
                                        'storage & other facilities',
                                        'boiler & turbogenerator',
                                        'cooling utility facilities',
                                        'other facilities',
                                        'heat exchanger network',
                                        'natural gas (for steam generation)',
                                        'natural gas (for product drying)',
                                        'chilled brine',
                                        'fixed operating cost',
                                        'electricity consumption',
                                        'heating duty',
                                        'excess electricity',
                                        ]
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
        
    if 'storage & other facilities' in groups_to_get:
        storage_group = UnitGroup('storage & other facilities',
        units=[i for i in unit_groups_temp if has_area_unit(i.units, 600)][0].units+\
            [i for i in unit_groups_temp if has_area_unit(i.units, 900)][0].units,
            )
        # storage_group.name = 'storage & other facilities'
        unit_groups_.append(storage_group)
    
    if 'boiler & turbogenerator' in groups_to_get:
        boiler_turbogenerator_group = UnitGroup('boiler & turbogenerator', 
                                                    units=(u.BT701,))
        unit_groups_.append(boiler_turbogenerator_group)
    
    if 'cooling utility facilities' in groups_to_get:
        cooling_utility_facilities_group = UnitGroup('cooling utility facilities', 
                                                         units=(u.CT801, u.CWP802,))
        unit_groups_.append(cooling_utility_facilities_group)

    # if 'other facilities' in groups_to_get:
    #     other_facilities_group = [i for i in unit_groups_temp if has_area_unit(i.units, 900)][0]
    #     other_facilities_group.name = 'other facilities'
    #     unit_groups_.append(other_facilities_group)
    
    if 'heat exchanger network' in groups_to_get:
        heat_exchanger_network_group = UnitGroup('heat exchanger network', 
                                                         units=(u.HXN1001,))
        unit_groups_.append(heat_exchanger_network_group)
        
    if 'natural gas (for steam generation)' in groups_to_get:
        unit_groups_.append(UnitGroup('natural gas (for steam generation)'))
        
    if 'natural gas (for product drying)' in groups_to_get:
        unit_groups_.append(UnitGroup('natural gas (for product drying)'))
    
    if 'chilled brine' in groups_to_get:
        unit_groups_.append(UnitGroup('chilled brine'))
    
    if 'fixed operating cost' in groups_to_get:
        unit_groups_.append(UnitGroup('fixed operating cost'))
    # unit_groups_.append(UnitGroup('material cost'))
    
    if 'electricity consumption' in groups_to_get:
        unit_groups_.append(UnitGroup('electricity consumption'))
        
    if 'heating duty' in groups_to_get:
        unit_groups_.append(UnitGroup('heating duty'))
        
    # unit_groups_.append(UnitGroup('cooling duty'))
    
    if 'excess electricity' in groups_to_get:
        unit_groups_.append(UnitGroup('excess electricity'))
        
    return unit_groups_

# %% Update units of all metrics in a list of unit groups
def update_metric_units_of_unit_groups(unit_groups):
    for ui in unit_groups:
        for i in ui.metrics:
            name = i.name.lower()
            if name in ('heating duty', 'cooling duty', 'steam use'):
                # i.units = 'GJ/h'
                i.units = 'GJ\u00b7h\u207b\u00b9'
                # i.units = f'GJ\u00b7h{chr(0x2D)}\u00b9'
                # i.units = 'GJ\u00b7$\mathdefault{h^{-1}}$'
            elif name in ('electricity consumption', 'power consumption'):
                i.units = 'MW'
            elif name in ('operating cost', 'material cost'):
                # i.units = 'MM$/y'
                # i.units = 'MM$\u00b7y\u207b\u00b9'
                # i.units = f'GJ\u00b7h{chr(0x2D)}\u00b9'
                # i.units = 'USD\u00b7$\mathdefault{h^{-1}}$'
                
                # i.units = '$\u00b7h\u207b\u00b9'
                # i.units = r"$\mathrm{\$}\cdot\mathrm{h}^{-1}$"
                i.units = r"$\mathrm{\$}$" + '\u00b7h\u207b\u00b9'
                
            elif name in ('installed equipment cost'):
                i.units = 'MM$'
                
#%% Add metrics to a given list of unit groups
def add_metrics_to_unit_groups(unit_groups, 
                               system,
                               TEA=None, 
                               LCA=None,
                               BT=None,
                               hxn_class=HeatExchangerNetwork,
                               natural_gas_utilizing_non_BT_unit_classes = [DrumDryer,],
                               ):
    if not TEA:
        TEA = system.TEA
    if not LCA:
        LCA = system.LCA
    if not BT:
        BT = system.flowsheet.unit.BT701
    
    natural_gas_utilizing_non_BT_system_units = [ui for ui in system.units
                                                 if ui.__class__ in natural_gas_utilizing_non_BT_unit_classes]
    

    
        
    for i in unit_groups: i.autofill_metrics(shorthand=False, 
                                             electricity_production=False, 
                                             electricity_consumption=True,
                                             material_cost=True)
    
    for i in unit_groups:
        for j in i.metrics:
            if j.name.lower()=='material cost':
                j.name = 'Operating cost'
        i.metrics.append(Metric('Steam use', 
                                            getter=lambda: 0, 
                                            units='GJ/h',
                                            element=None))
        # i.metrics.append(Metric('Operating cost', 
        #                                     getter=lambda: 0, 
        #                                     units='MM$/y',
        #                                     element=None))
                         
    for i in unit_groups:
        # i.metrics[-3].getter = i.get_material_cost*TEA.annual_operating_hours
        
        if i.name == 'heat exchanger network':
            i.filter_savings = False
            assert isinstance(i.units[0], hxn_class)
        if i.name == 'cooling utility facilities':
            i.metrics[1].getter = lambda: 0. # Cooling duty
            i.metrics[3].getter = lambda: sum([ui.power_utility.rate for ui in i.units])/1e3 # Electricity consumption [MW]
        
        ############## Operating cost ##############
        if i.name == 'feedstock juicing':
            i.metrics[-2].getter = lambda: sum([sum([j.cost for j in ui.ins]) for ui in i.units])
        if i.name == 'storage & other facilities' or i.name == 'cooling utility facilities':
            i.metrics[-2].getter = lambda: 0.
        if i.name == 'natural gas (for steam generation)':
            i.metrics[-2].getter=lambda: BT.natural_gas_price * BT.natural_gas.F_mass
        if i.name == 'natural gas (for product drying)':
            i.metrics[-2].getter=lambda: sum([j.utility_cost-j.power_utility.cost
                                              for j in natural_gas_utilizing_non_BT_system_units])
        if i.name == 'chilled brine':
            i.metrics[-2].getter=lambda: sum([i.cost for i in system.heat_utilities 
                                              if i.agent.ID in ('chilled_brine', 'chilled brine')])
            
        ##############  ##############
        if i.name == 'fixed operating cost':
            i.metrics[-2].getter=lambda: TEA.FOC/TEA.operating_hours
            
        # if i.name == 'material cost':
        #     i.metrics[-2].getter=lambda: TEA.material_cost/TEA.operating_hours\
        #                                 +BT.natural_gas_price * BT.natural_gas.F_mass\
        #                                 +sum([j.utility_cost-j.power_utility.cost
        #                                      for j in natural_gas_utilizing_non_BT_system_units])\
        #                                 +sum([i.cost for i in system.heat_utilities 
        #                                       if i.agent.ID in ('chilled_brine', 'chilled brine')])
            
        if i.name == 'excess electricity':
            i.metrics[-2].getter=\
                lambda: TEA.utility_cost/TEA.operating_hours -\
                                        sum([i.cost for i in system.heat_utilities 
                                            if i.agent.ID in ('chilled_brine', 'chilled brine')])-\
                                        sum([j.utility_cost-j.power_utility.cost
                                              for j in natural_gas_utilizing_non_BT_system_units])
                
                # lambda: system.power_utility.cost
        ############## Electricity consumption ##############
        if i.name == 'cooling utility facilities':
            i.metrics[3].getter=lambda: sum([ui.power_utility.rate for ui in i.units])
            
        ############## Steam use ##############
        if i.name == 'heating duty':
            i.metrics[-1].getter=lambda: LCA.actual_steam_frac_heating*LCA.BT_steam_kJph_total/1e6
            
        # if i.name == 'cooling duty':
        #     i.metrics[-2].getter=lambda: LCA.actual_steam_frac_cooling*LCA.BT_steam_kJph_total/1e6
            
        if i.name == 'electricity consumption':
            i.metrics[-1].getter=lambda: (LCA.actual_steam_frac_electricity_non_cooling+\
                                          LCA.actual_steam_frac_cooling)*\
                                          LCA.BT_steam_kJph_total/1e6
            
        if i.name == 'excess electricity':
            i.metrics[-1].getter=lambda: LCA.actual_steam_frac_excess*LCA.BT_steam_kJph_total/1e6
        
        update_metric_units_of_unit_groups(unit_groups)
        
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
