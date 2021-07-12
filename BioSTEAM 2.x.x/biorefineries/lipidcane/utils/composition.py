# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module defines the composition_balance function, which performs mass 
energy balance the given oil composition of lipid cane to determine its
composition.

"""
from thermosteam import Stream
from biosteam import main_flowsheet as f
from numpy import array

__all__ = ('set_lipid_fraction',
           'get_lipid_fraction')

def set_lipid_fraction(lipid_fraction, stream=None, data={}):
    """Adjust composition of lipid cane to achieve desired oil fraction (dry weight)."""
    if not stream: stream = f.stream.lipidcane
    carbs_IDs = ('Glucose', 'Sucrose')
    fiber_IDs = ('Lignin', 'Cellulose', 'Hemicellulose')
    lipid_IDs = ('PL', 'FFA', 'MAG', 'DAG', 'TAG')
    if not data:
        carbs   = Stream(None, thermo=stream.thermo)
        fiber   = Stream(None, thermo=stream.thermo)
        lipid   = Stream(None, thermo=stream.thermo)
        carbs.imol[carbs_IDs] = stream.imol[carbs_IDs]
        fiber.imol[fiber_IDs] = stream.imol[fiber_IDs]
        lipid.imol[lipid_IDs] = stream.imol[lipid_IDs]
        
        # Mass property arrays
        carbs_mass = stream.imass[carbs_IDs]
        fiber_mass = stream.imass[fiber_IDs]
        lipid_mass = stream.imass[lipid_IDs]
        
        # Net weight
        carbs_massnet = carbs_mass.sum()
        fiber_massnet = fiber_mass.sum()
        lipid_massnet = lipid_mass.sum()
        
        # Heats of combustion per kg
        LHV_carbs_kg = carbs.LHV/carbs_massnet
        LHV_lipid_kg = lipid.LHV/lipid_massnet
        data['LHV_lipid_over_carbs'] = LHV_lipid_kg / LHV_carbs_kg
        
        
        # Relative composition
        data['r_mass_carbs'] = array(carbs_mass/carbs_massnet)
        data['r_mass_fiber'] = array(fiber_mass/fiber_massnet)
        data['r_mass_lipid'] = array(lipid_mass/lipid_massnet)
        
    z_mass_lipid = lipid_fraction
    F_mass = stream.F_mass
    z_dry = 0.3
    z_mass_lipid = z_mass_lipid * z_dry 
    z_mass_carbs = 0.149 - z_mass_lipid * data['LHV_lipid_over_carbs']
    z_mass_fiber = z_dry - z_mass_carbs - z_mass_lipid - 0.006 - 0.015
    imass = stream.imass
    imass['Water'] = F_mass * 0.7
    imass['Ash'] = F_mass * 0.006
    imass['Solids'] = F_mass * 0.015
    imass[lipid_IDs] = data['r_mass_lipid'] * z_mass_lipid * F_mass
    imass[carbs_IDs] = data['r_mass_carbs'] * z_mass_carbs * F_mass
    imass[fiber_IDs] = data['r_mass_fiber'] * z_mass_fiber * F_mass
    if any(stream.mol < 0):
        raise ValueError(f'lipid cane oil composition of {z_mass_lipid/z_dry*100:.0f}% dry weight is infeasible')

def get_lipid_fraction(stream=None):
    if not stream: stream = f.stream.lipidcane
    imass = stream.imass
    F_dry_mass = stream.F_mass - imass['Water']
    return imass['Lipid'] / F_dry_mass