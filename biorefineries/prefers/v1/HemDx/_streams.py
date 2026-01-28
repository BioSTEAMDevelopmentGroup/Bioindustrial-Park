# -*- coding: utf-8 -*-
"""
Created on 2025-05-07 18:26:22

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biosteam import stream_kwargs
from biorefineries.prefers.v1._process_settings import price  # ADD THIS IMPORT
import biosteam as bst
from biosteam.units import Fermentation
from httpx import stream
import numpy as np
from sympy import Trace  # ADD THIS IMPORT
LegHb={}


productivity = 4.2/1000 # [g / L / h]
tau = 90 # [h]
titer = productivity * tau # [g / L]
# OD600 180
# 1 OD600 = 0.35 g/L dry cell weight


m=2/5 # scale up factor based on 1.5e5 kg/yr
# %% In
# TraceMetalSolution = stream_kwargs('TraceMetalSolution', TraceMetalSolution=1, units='kg/hr', T=25+273.15)
# VitaminCSolution = stream_kwargs('VitaminCSolution', VitaminC=1, units='kg/hr', T=25+273.15)
# YPD = stream_kwargs('YPD', YPD_Medium=1, units='kg/hr', T=25+273.15)
# YPD_Medium = stream_kwargs('YPD_Medium', YPD_Medium=1, units='kg/hr', T=25+273.15)
# Supplemented = stream_kwargs('Supplemented', Supplemented=1, units='kg/hr', T=25+273.15)
# Fermentation_Medium = stream_kwargs('Supplemented_Medium', Fermentation_Medium=1, units='kg/hr', T=25+273.15)
# Feed1st = stream_kwargs('Feed1st', Fermentation_Medium=1, units='kg/hr', T=25+273.15)
# Feed2nd = stream_kwargs('Feed2nd', Fermentation_Medium=1, units='kg/hr', T=25+273.15)
# Feed1st_Solution = stream_kwargs('Feed1st_Solution', Fermentation_Medium=1, units='kg/hr', T=25+273.15)
# Feed2nd_Solution = stream_kwargs('Feed2nd_Solution', Fermentation_Medium=1, units='kg/hr', T=25+273.15)

SeedIn1 = stream_kwargs('SeedIn1', Seed=1, units='kg/hr', T=25+273.15)
SeedIn2 = stream_kwargs('SeedIn2', Seed=1, units='kg/hr', T=25+273.15)
CultureIn = stream_kwargs('CultureIn', Culture=1, units='kg/hr', T=25+273.15)
SeedSolution1 = stream_kwargs('SeedSolution1', SeedSolution=m*0.15*1e5/16, units='kg/hr', T=25+273.15)
SeedSolution2 = stream_kwargs('SeedSolution2', SeedSolution=m*1.5*1e5/16, units='kg/hr', T=25+273.15)
CultureSolution = stream_kwargs('CultureSolution', CultureSolution=m*1.5*1e5/16, units='kg/hr', T=25+273.15)

Glucose = stream_kwargs('Glucose', Glucose=m*1.3*1e5/(72+90), units='kg/hr', T=25+273.15)
NH3_25wt = stream_kwargs('NH3_25wt', NH3_25wt=m*230, units='kg/hr', T=25+273.15)#, price=price['NH3_25wt'])

DfUltraBuffer1 = stream_kwargs('DfUltraBuffer1', DfUltraBuffer=1, units='kg/hr', T=25+273.15)
DfUltraBuffer2 = stream_kwargs('DfUltraBuffer2', DfUltraBuffer=1, units='kg/hr', T=25+273.15)
NaCl_wash = stream_kwargs('NaCl_wash', NaCl=1, units='kg/hr', T=25+273.15)
NaOH_elute = stream_kwargs('NaOH_elute', NaOH=1, units='kg/hr', T=25+273.15)
Ethanol_regen = stream_kwargs('Ethanol_regen', Ethanol=1, units='kg/hr', T=25+273.15)

GammaCyclodextrinFeed = stream_kwargs('GammaCyclodextrinFeed', GammaCyclodextrin=1, units='kg/hr', T=25+273.15)
NicotinamideFeed = stream_kwargs('NicotinamideFeed', Nicotinamide=1, units='kg/hr', T=25+273.15)
AdsorbWash = stream_kwargs('AdsorbWash', H2O=1, units='kg/hr', T=25+273.15)
AdsorbElute = stream_kwargs('AdsorbElute', H2O=1, units='kg/hr', T=25+273.15)
AdsorbRegen = stream_kwargs('AdsorbRegen', H2O=1, units='kg/hr', T=25+273.15)
IXEquilibriumBuffer = stream_kwargs('IXEquilibriumBuffer', IXEquilibriumBuffer=1, units='kg/hr', T=25+273.15)
IXElutionBuffer = stream_kwargs('IXElutionBuffer', IXElutionBuffer=1, units='kg/hr', T=25+273.15)
IXRegenerationSolution = stream_kwargs('IXRegenerationSolution', NaOH=1, units='kg/hr', T=25+273.15)
DfNanoBuffer = stream_kwargs('DfNanoBuffer', DfNanoBuffer=1, units='kg/hr', T=25+273.15)
AntioxidantStream = stream_kwargs('AntioxidantStream', SodiumAscorbate=0.1, H2O=1.0, units='kg/hr', price=price.get('SodiumAscorbate', 5.0))

# %% Out

LegHb_1 = stream_kwargs('LegHb_1')
LegHb_2 = stream_kwargs('LegHb_2')
LegHb_3 = stream_kwargs('LegHb_3',units='kg/hr', T=25+273.15, price = 20)#11.006028727167918)
vent1 = stream_kwargs('vent1')
vent2 = stream_kwargs('vent2')
effluent1 = stream_kwargs('effluent1')
effluent2 = stream_kwargs('effluent2', units='kg/hr', price=-0.33)  # organic waste remove
effluent3 = stream_kwargs('effluent3')
NHemDx_Product = stream_kwargs('NHemDx_Product', units='kg/hr', T=25+273.15, price=20)#  20 for pharmaceutical grade



# %%
# Add missing pricing functions directly here since they're not in your _streams2.py
def set_stream_price_from_components(stream, chemical_prices=None):
    """Set stream price based on individual chemical component prices"""
    from biorefineries.prefers.v1._process_settings import price
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

def update_all_input_stream_prices(streamlist, verbose=True):
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
                if verbose:
                    print(f"{stream.ID}: ${old_price:.4f}/kg -> ${new_price:.4f}/kg")
            else:
                print(f"Warning: Stream {stream} is not properly initialized")
    except Exception as e:
        print(f"Error updating stream prices: {e}")