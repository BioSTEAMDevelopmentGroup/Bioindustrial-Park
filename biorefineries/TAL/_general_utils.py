#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from biosteam import UnitGroup
from biosteam.evaluation import Metric
from biosteam import HeatExchangerNetwork, DrumDryer
from thermosteam import equilibrium
from flexsolve import IQ_interpolation
from numba import njit
from math import log
import numpy as np

__all__ = {'call_all_specifications_or_run',
           'get_more_unit_groups',
           'add_metrics_to_unit_groups',
           'set_production_capacity',
           'TEA_breakdown',
           'update_facility_IDs',
           'get_pH_polyprotic_acid_mixture',
           'get_major_units_df'
           }

#%% Estimate the pH of a mixture of polyprotic acids
def get_molarity(ID, stream):
    return stream.imol[ID]/stream.F_vol

@njit
def helper_acid_contribution_upto_triprotic(cH, acid_molarity, Ka1, Ka2, Ka3):
    return acid_molarity * ((Ka1*cH**2 + 2*cH*Ka1*Ka2 + 3*Ka1*Ka2*Ka3) / (cH**3 + Ka1*cH**2 + cH*Ka1*Ka2 + Ka1*Ka2*Ka3))

def get_acid_contribution_upto_triprotic(cH, acid_molarity, Kas):
    n_Kas = len(Kas)
    Ka1 = Kas[0]
    Ka2 = Kas[1] if n_Kas>1 else 0.
    Ka3 = Kas[2] if n_Kas>2 else 0.
    return helper_acid_contribution_upto_triprotic(cH, acid_molarity, Ka1, Ka2, Ka3)

def obj_f_cH_polyprotic_acid_mixture(cH, acid_molarities, Kas):
    return cH - (10**-14)/cH - sum([get_acid_contribution_upto_triprotic(cH, m, k) 
                                    for m,k in zip(acid_molarities, Kas)])

def get_cH_polyprotic_acid_mixture(stream, acid_IDs, Kas, activities):
    acid_molarities = [get_molarity(i, stream) for i in acid_IDs]
    gammas = np.ones(len(acid_molarities))
    
    if activities=='UNIFAC':
        stream_chems = stream.chemicals
        gamma_obj = equilibrium.UNIFACActivityCoefficients(stream_chems)
        indices = [stream_chems.index(i) for i in acid_IDs]
        gammas = gamma_obj(stream.mol[:], stream.T)[indices]
    elif activities=='Dortmund':
        stream_chems = stream.chemicals
        gamma_obj = equilibrium.DortmundActivityCoefficients(stream_chems)
        indices = [stream_chems.index(i) for i in acid_IDs]
        gammas = gamma_obj(stream.mol[:], stream.T)[indices]
        
    acid_molarities = np.multiply(acid_molarities, gammas)
    
    obj_f = lambda cH: 1000.* obj_f_cH_polyprotic_acid_mixture(cH, acid_molarities, Kas)
    return IQ_interpolation(obj_f, 10**(-14), 10**(-1), ytol=1e-6)

def get_pH_polyprotic_acid_mixture(stream, acid_IDs=[], Kas=[], activities='ideal'):
    implemented_activities = ('ideal', 'UNIFAC', 'Dortmund')
    if activities not in implemented_activities:
        raise ValueError(f"Parameter 'activities' must be one of {implemented_activities}, not {activities}.")
    return -log(get_cH_polyprotic_acid_mixture(stream, acid_IDs, Kas, activities), 10.)


#%% For a given list of units, call all specifications of each unit or run each unit (in the presented order)
def call_all_specifications_or_run(units_to_run):
    if units_to_run.__class__ in (list, tuple):
        for unit_to_run in units_to_run:
            if unit_to_run.specifications: 
                [i() for i in unit_to_run.specifications]
            else:
                unit_to_run.run()
    else:
        if units_to_run.specifications: 
            [i() for i in units_to_run.specifications]
        else:
            units_to_run.run()

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
                                        # 'chilled brine',
                                        'fixed operating cost',
                                        'electricity consumption',
                                        'heating duty',
                                        'excess electricity',
                                        ],
                         wastewater_areas=[500,],
                         storage_and_other_facilities_areas=[600,900],
                         has_brine_facility=True,
                         ):
    unit_groups_temp = UnitGroup.group_by_area(system.units)
    u = system.flowsheet.unit
    
    def has_area_unit(units, area):
        for ui in units:
            if ui.ID[1] == str(area)[0]: return True
        return False
    
    unit_groups_ = []
    
    if 'wastewater' in groups_to_get:
        wastewater_group = [i for i in unit_groups_temp if has_area_unit(i.units, wastewater_areas[0])][0]
        wastewater_group.name = 'wastewater'
        unit_groups_.append(wastewater_group)
        
    if 'storage & other facilities' in groups_to_get:
        storage_group = UnitGroup('storage & other facilities',
        units=[i for i in unit_groups_temp if has_area_unit(i.units, storage_and_other_facilities_areas[0])][0].units+\
            [i for i in unit_groups_temp if has_area_unit(i.units, storage_and_other_facilities_areas[1])][0].units,
            )
        # storage_group.name = 'storage & other facilities'
        unit_groups_.append(storage_group)
    
    if 'boiler & turbogenerator' in groups_to_get:
        boiler_turbogenerator_group = UnitGroup('boiler & turbogenerator', 
                                                    units=(u.BT701,))
        unit_groups_.append(boiler_turbogenerator_group)
    
    if 'cooling utility facilities' in groups_to_get:
        cuf_units = [u.CT801, u.CWP802]
        if has_brine_facility: cuf_units.append(u.CWP803)
        cooling_utility_facilities_group = UnitGroup('cooling utility facilities', 
                                                         units=cuf_units)
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
def update_metric_unit_labels_of_unit_groups(unit_groups):
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
                
                
                # i.units = r"$\mathrm{\$}$" + '\u00b7h\u207b\u00b9'
                i.units = r"$\mathrm{MM\$}$" + '\u00b7y\u207b\u00b9'
                
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
        if i.name == 'boiler & turbogenerator':
            i.metrics[-2].getter = lambda: sum([i.cost for i in BT.ins]) +\
                                        abs(BT.ash_disposal_price*BT.ash_disposal.F_mass)
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
                lambda: system.power_utility.cost
                # lambda: TEA.utility_cost/TEA.operating_hours -\
                #                         sum([i.cost for i in system.heat_utilities 
                #                             if i.agent.ID in ('chilled_brine', 'chilled brine')])-\
                #                         sum([j.utility_cost-j.power_utility.cost
                #                               for j in natural_gas_utilizing_non_BT_system_units])

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
        
        
        update_metric_unit_labels_of_unit_groups(unit_groups)
        
        
    # Rearrange metrics
    metrics_list_ordered = [
                            'Installed equipment cost', 
                            'Operating cost', 
                            'Steam use',
                            'Electricity consumption',
                            'Cooling duty', 
                            'Heating duty', 
                            ]
    for i in unit_groups:
        i.metrics.sort(key = lambda x: metrics_list_ordered.index(x.name))
        
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
    
    n_iterations_analytical = 1
    if method=='analytical':
        for i in range(n_iterations_analytical):
            if spec: spec.load_specifications(spec_1=spec.spec_1,
                                              spec_2=spec.spec_2,
                                              spec_3=spec.spec_3,)
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
                  print_output=False): # operating cost is in $/h, but unit name is MM$/y  by default for plotting purposes
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

#%%
import biosteam as bst
import pandas as pd
import re

def get_major_units_df(units, unit_groups, save_filename='major_units.xlsx', 
                       # non_biosteam_sources={},
                       remove_units_with_no_equipment=False):
    u = units
    IDs, lines, equipment, areas, sources = [], [], [], [], []
    for ui in u:
        IDs.append(ui.ID)
        lines.append(ui.line)
        equipment.append(list([i.lower() for i in ui.baseline_purchase_costs.keys()])
                         + list(ui.auxiliary_unit_names))
        
        found_area = False
        for j in unit_groups:
            if ui in j.units:
                areas.append(j.name)
                found_area = True
        if not found_area: areas.append('?')
        
        # found_source = False
        # for v in bst.units.__dict__.values(): 
        #     if type(v)==type:
        #         if v.__name__.lower() == ui.__class__.__name__.lower():
        #             sources.append('BioSTEAM')
        #             found_source = True
        #             break
        # if not found_source:
        #     for src_ID, src in non_biosteam_sources.items():
        #         for w in src.__dict__.values(): 
        #             if type(w)==type:
        #                 if w.__name__.lower() == ui.__class__.__name__.lower():
        #                     sources.append(src_ID)
        #                     found_source = True
        #                     break
        #         if found_source: break
        # if not found_source: sources.append('?')
        sources.append(ui.__class__.__module__)
    
    equipment_strings = [str(i) for i in equipment]
    equipment_strings = [re.sub('_', ' ', i) for i in equipment_strings]
    equipment_strings = [re.sub('\[', '', i) for i in equipment_strings]
    equipment_strings = [re.sub('\]', '', i) for i in equipment_strings]
    equipment_strings = [re.sub("'", '', i) for i in equipment_strings]
    
    exc_ind = [i for i in range(len(equipment)) if not equipment[i]] if remove_units_with_no_equipment else []
    
    areas, IDs, lines, equipment_strings, sources =\
        exclude_given_indices_from_list(areas, exc_ind),\
        exclude_given_indices_from_list(IDs, exc_ind),\
        exclude_given_indices_from_list(lines, exc_ind),\
        exclude_given_indices_from_list(equipment_strings, exc_ind),\
        exclude_given_indices_from_list(sources, exc_ind)
    
    df = pd.DataFrame(data={'Process': areas, 'ID': IDs, 'Unit': lines, 'Equipment': equipment_strings, 'Sources': sources})
    df.to_excel(save_filename)
    

    return df

def replace_first_instance_in_string(given_string, old_partial_string, new_partial_string):
    if old_partial_string in given_string:
        print(old_partial_string, given_string)
        partial_string_index = given_string.index(old_partial_string)
        len_partial_string = len(old_partial_string)
        return given_string[0:partial_string_index] +\
            new_partial_string +\
            given_string[partial_string_index+len_partial_string:]
    else:
        return given_string
    
def exclude_given_indices_from_list(given_list, given_indices):
    return [given_list[i] for i in range(len(given_list)) if not i in given_indices]

def identify_accumulating_and_depleting_streams(sys, threshold_F_mol=1):
    streams = list(sys.flowsheet.stream)
    stream_F_mol_1 = [i.F_mol for i in streams]
    sys.simulate()
    stream_F_mol_2 = [i.F_mol for i in streams]
    
    acc_indices = [i for i in range(len(streams)) if stream_F_mol_2[i] - stream_F_mol_1[i] > threshold_F_mol]
    
    dep_indices = [i for i in range(len(streams)) if stream_F_mol_1[i] - stream_F_mol_2[i] > threshold_F_mol]
    
    return [streams[i] for i in acc_indices], [streams[i] for i in dep_indices]
