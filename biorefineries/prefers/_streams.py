# -*- coding: utf-8 -*-
"""
Created on 2025-05-07 18:26:22

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biosteam import stream_kwargs
from biorefineries.prefers._process_settings import price  # ADD THIS IMPORT
import biosteam as bst
import numpy as np  # ADD THIS IMPORT
LegH={}

titer = 7.27 # [g / L]
productivity = titer / 72 # [g / L / h]
LegH_yield = titer * 5 / 1300 # [by wt]

m=2 # scale up factor based on 1.5e5 kg/yr
# %% In
SeedIn = stream_kwargs('SeedIn', Seed=m*0.15*1e5/16, units='kg/hr', T=25+273.15)
CultureIn = stream_kwargs('CultureIn', Culture=m*1.5*1e5/16, units='kg/hr', T=25+273.15)

Glucose = stream_kwargs('Glucose', Glucose=m*1.3*1e5/72, units='kg/hr', T=25+273.15)
NH3_25wt = stream_kwargs('NH3_25wt', NH3_25wt=m*350, units='kg/hr', T=25+273.15)#, price=price['NH3_25wt'])

DfUltraBuffer = stream_kwargs('DfUltraBuffer', DfUltraBuffer=1, units='kg/hr', T=25+273.15)
IXEquilibriumBuffer = stream_kwargs('IXEquilibriumBuffer', IXEquilibriumBuffer=1, units='kg/hr', T=25+273.15)
IXElutionBuffer = stream_kwargs('IXElutionBuffer', IXElutionBuffer=1, units='kg/hr', T=25+273.15)
IXRegenerationSolution = stream_kwargs('IXRegenerationSolution', NaOH =1, units='kg/hr', T=25+273.15)
DfNanoBuffer = stream_kwargs('DfNanoBuffer', DfNanoBuffer=1, units='kg/hr', T=25+273.15)

# %% Out

LegH_1 = stream_kwargs('LegH_1')
LegH_2 = stream_kwargs('LegH_2')
LegH_3 = stream_kwargs('LegH_3',units='kg/hr', T=25+273.15, price = 25)
vent1 = stream_kwargs('vent1')
vent2 = stream_kwargs('vent2')
effluent1 = stream_kwargs('effluent1')
effluent2 = stream_kwargs('effluent2')
effluent3 = stream_kwargs('effluent3')
effluent4 = stream_kwargs('effluent4')
effluent5 = stream_kwargs('effluent5')
effluent6 = stream_kwargs('effluent6')

# %%
# Add missing pricing functions directly here since they're not in your _streams2.py
def set_stream_price_from_components(stream, chemical_prices=None):
    """Set stream price based on individual chemical component prices"""
    from biorefineries.prefers._process_settings import price
    if chemical_prices is None:
        chemical_prices = price
    
    try:
        # FIX: Check if stream is a proper Stream object
        if not hasattr(stream, 'imass'):
            print(f"Warning: {stream} is not a proper Stream object, skipping pricing")
            return 0
            
        mass = stream.imass
        total_cost = 0
        total_mass = stream.F_mass
        
        # Handle individual chemicals in groups
        for chemical_id in stream.chemicals.IDs:
            if mass[chemical_id] > 0:
                # Check if this is a grouped chemical with group composition
                chemical = stream.chemicals[chemical_id]
                if hasattr(chemical, '_group_wt_composition') and chemical._group_wt_composition:
                    # This is a group - calculate based on group composition
                    group_mass = mass[chemical_id]
                    for component_id, fraction in chemical._group_wt_composition.items():
                        if component_id in chemical_prices:
                            component_cost = group_mass * fraction * chemical_prices[component_id]
                            total_cost += component_cost
                else:
                    # Individual chemical
                    if chemical_id in chemical_prices:
                        total_cost += mass[chemical_id] * chemical_prices[chemical_id]
        
        stream.price = total_cost / total_mass if total_mass > 0 else 0
        return stream.price
        
    except Exception as e:
        print(f"Warning: Could not set price for {getattr(stream, 'ID', 'unknown')}: {e}")
        if hasattr(stream, 'price'):
            stream.price = 0
        return 0

def create_stream(stream_list):
    """Create a new stream with default properties"""
    streams = []
    for stream_kwargs in stream_list:
        s = bst.Stream(**stream_kwargs)
        streams.append(s)
    return streams

def update_all_input_stream_prices(streamlist):
    """Update prices for all input streams"""
    # FIX: Check if streams exist and are proper objects
    try:
        input_streams = streamlist
        # input_streams = create_stream([SeedIn, CultureIn, Glucose, NH3_25wt, 
        #                 BufferA, BufferB, BufferC])
        
        for stream in input_streams:
            # FIX: Verify each stream is a proper Stream object
            if hasattr(stream, 'imass') and hasattr(stream, 'ID'):
                old_price = getattr(stream, 'price', 0)
                new_price = set_stream_price_from_components(stream)
                print(f"{stream.ID}: ${old_price:.4f}/kg → ${new_price:.4f}/kg")
            else:
                print(f"Warning: Stream {stream} is not properly initialized")
    except Exception as e:
        print(f"Error updating stream prices: {e}")