#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microalgae biorefinery analysis utilities

Helper functions for creating LCA objects and other analysis utilities.

@author: Xingdong Shi
@version: 0.0.1
"""

from . import lca

def create_microalgae_lca_simple(system, tea=None):
    """
    Simplified function to create microalgae LCA object with default parameters.
    
    Parameters
    ----------
    system : System
        Microalgae production system
    tea : TEA, optional
        TEA object (not used in current implementation)
        
    Returns
    -------
    LCA
        Life cycle assessment object
    """
    
    # Find main product
    main_product = None
    for stream in system.products:
        if stream.get_total_flow('kg/hr') > 0:
            main_product = stream
            break
    
    if main_product is None:
        # If no products found, try to find caproic acid stream
        for stream in system.streams:
            if 'caproic' in stream.ID.lower():
                main_product = stream
                break
    
    # Find boiler
    boiler = None
    for unit in system.units:
        if 'BT' in unit.ID or 'boiler' in unit.ID.lower():
            boiler = unit
            break
    
    # Default main product chemical IDs
    main_product_chemical_IDs = ['CaproicAcid']
    
    # Create LCA object
    microalgae_lca = lca.create_microalgae_lca(
        system=system,
        main_product=main_product,
        main_product_chemical_IDs=main_product_chemical_IDs,
        boiler=boiler
    )
    
    return microalgae_lca

def get_fermentation_unit(system):
    """
    Find the fermentation unit in the system.
    
    Parameters
    ----------
    system : System
        Microalgae production system
        
    Returns
    -------
    Unit or None
        Fermentation unit if found, None otherwise
    """
    for unit in system.units:
        if ('ferment' in unit.ID.lower() or 
            'MCCA' in unit.ID or 
            'R' in unit.ID and hasattr(unit, 'conversion')):
            return unit
    return None

def get_main_product_stream(system):
    """
    Find the main product stream in the system.
    
    Parameters
    ----------
    system : System
        Microalgae production system
        
    Returns
    -------
    Stream or None
        Main product stream if found, None otherwise
    """
    # First try to find streams with significant flow
    for stream in system.products:
        if stream.get_total_flow('kg/hr') > 0:
            return stream
    
    # If no products found, try to find caproic acid stream
    for stream in system.streams:
        if 'caproic' in stream.ID.lower() and stream.get_total_flow('kg/hr') > 0:
            return stream
    
    return None 