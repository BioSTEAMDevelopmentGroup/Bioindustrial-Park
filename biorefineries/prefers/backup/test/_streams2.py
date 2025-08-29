# -*- coding: utf-8 -*-
"""
Created on 2025-08-12 14:03:02

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biosteam import stream_kwargs
from biorefineries.prefers._process_settings import price  # ADD THIS IMPORT
import numpy as np  # ADD THIS IMPORT

# %% Manual price update function
def setup_stream_pricing():
    """Call this function to update all stream prices"""
    return update_all_input_stream_prices()

def set_stream_price_from_components(stream, chemical_prices=None):
    """
    Set stream price based on individual chemical component prices
    Automatically calculates price from grouped chemical compositions
    
    Parameters:
    -----------
    stream : tmo.Stream
        Stream to set price for
    chemical_prices : dict, optional
        Chemical prices, defaults to global price dict
    """
    if chemical_prices is None:
        chemical_prices = price
    
    try:
        mass = stream.imass
        total_cost = 0
        total_mass = stream.F_mass
        
        # Handle grouped chemicals - expand to individual components
        for chemical_id in stream.chemicals.IDs:
            if mass[chemical_id] > 0:
                # Check if this is a grouped chemical
                if hasattr(stream.chemicals, '_groups') and chemical_id in stream.chemicals._groups:
                    # This is a group - calculate based on group composition
                    group_mass = mass[chemical_id]
                    group_def = stream.chemicals._groups[chemical_id]
                    
                    # Calculate weighted cost for group components
                    group_total_fraction = sum(group_def.data.values())
                    for component_id, fraction in group_def.data.items():
                        if component_id in chemical_prices:
                            component_fraction = fraction / group_total_fraction
                            component_cost = group_mass * component_fraction * chemical_prices[component_id]
                            total_cost += component_cost
                else:
                    # Individual chemical
                    if chemical_id in chemical_prices:
                        total_cost += mass[chemical_id] * chemical_prices[chemical_id]
        
        stream.price = total_cost / total_mass if total_mass > 0 else 0
        return stream.price
        
    except Exception as e:
        print(f"Error setting price for stream {stream.ID}: {e}")
        stream.price = 0
        return 0

def create_auto_priced_stream(ID, auto_update=True, **kwargs):
    """
    Create a stream with automatic price calculation based on composition
    
    Parameters:
    -----------
    ID : str
        Stream ID
    auto_update : bool
        Whether to add specification for automatic price updates
    **kwargs : 
        Stream composition and properties
    """
    stream = stream_kwargs(ID, **kwargs)
    
    # Set initial price
    set_stream_price_from_components(stream)
    
    if auto_update:
        # Add specification for automatic price updates when stream changes
        @stream.add_specification
        def auto_price_update():
            set_stream_price_from_components(stream)
    
    return stream

def update_all_input_stream_prices():
    """
    Update prices for all input streams based on their composition
    """
    input_streams = [SeedIn, CultureIn, Glucose, NH3_18wt, 
                    WashingSolution, Elution, NFBuffer]
    
    price_info = {}
    for stream in input_streams:
        old_price = getattr(stream, 'price', 0)
        new_price = set_stream_price_from_components(stream)
        price_info[stream.ID] = {
            'old_price': f"${old_price:.4f}/kg",
            'new_price': f"${new_price:.4f}/kg"
        }
        print(f"{stream.ID}: ${old_price:.4f}/kg → ${new_price:.4f}/kg")
    
    return price_info


# ...existing code...

# Replace your current stream definitions with auto-priced versions:

# %% In - WITH AUTO-PRICING
SeedIn = create_auto_priced_stream('SeedIn', Seed=0.15*1e5/16, units='kg/hr', T=25+273.15)
CultureIn = create_auto_priced_stream('CultureIn', Culture=1.5*1e5/16, units='kg/hr', T=25+273.15)
Glucose = create_auto_priced_stream('Glucose', Glucose=1.3*1e5/72, units='kg/hr', T=25+273.15)
NH3_18wt = create_auto_priced_stream('NH3_18wt', NH3_18wt=1000, units='kg/hr', T=25+273.15)
WashingSolution = create_auto_priced_stream('WashingSolution', DiaBuffer=1, units='kg/hr', T=25+273.15)
Elution = create_auto_priced_stream('Elution', IEXBuffer=1, units='kg/hr', T=25+273.15)
NFBuffer = create_auto_priced_stream('NFBuffer', NanoBuffer=1, units='kg/hr', T=25+273.15)

# %% Manual price update function
def setup_stream_pricing():
    """Call this function to update all stream prices"""
    return update_all_input_stream_prices()