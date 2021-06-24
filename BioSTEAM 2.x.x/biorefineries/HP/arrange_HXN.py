# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 13:57:14 2021

@author: sarangbhagwat
"""
#%% Imports
from biorefineries.HP.system_light_lle_vacuum_distillation import *
import itertools
from copy import copy 

#%% Setup
num_streams = len(HXN.streams)

original_heat_utils = HXN.original_heat_utils

cold_stream_indices = [index for index in range(num_streams) 
                       if original_heat_utils[index].duty>=0]

hot_stream_indices = [index for index in range(num_streams) 
                       if index not in cold_stream_indices]

is_cold = lambda i: i in cold_stream_indices
is_hot = lambda i: i in hot_stream_indices

new_HXs, new_HX_utils = HXN.new_HXs, HXN.new_HX_utils

stream_life_cycles = HXN.stream_life_cycles

individual_orders_dict = {i:[] for i in range(num_streams)}

for i in range(num_streams):
    slc = stream_life_cycles[i].life_cycle
    # position = 0
    for stage in slc:
        # individual_orders_dict[i][stage.unit.ID] = position
        # position += 1
        unit = stage.unit
        if not 'Util' in unit.ID and not 'util' in unit.ID:
            individual_orders_dict[i].append(unit.ID)

overall_order = [hx.ID for hx in new_HXs]

# all_possible_overall_orders = list(itertools.permutations(overall_order)) # inordinately large; n! permutations

def get_proposed_individual_order(i):
     proposed_individual_order = copy.copy(individual_orders_dict[i])
     proposed_individual_order.sort(key = lambda hx: overall_order.index(hx),
                                    reverse = is_hot(i))
     return proposed_individual_order