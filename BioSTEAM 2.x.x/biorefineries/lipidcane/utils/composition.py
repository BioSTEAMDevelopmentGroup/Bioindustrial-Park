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

def set_lipid_fraction(lipid_fraction, stream=None,
                       PL_fraction=0.,
                       FFA_fraction=0.,
                       z_mass_carbs_baseline=0.149,
                       z_mass_ash_baseline=0.006,
                       z_mass_solids_baseline=0.015):
    """Adjust composition of lipid cane to achieve desired oil fraction (dry weight)."""
    if not stream: stream = f.stream.lipidcane
    thermo = stream.thermo
    carbs_IDs = ('Glucose', 'Sucrose')
    fiber_IDs = ('Lignin', 'Cellulose', 'Hemicellulose')
    lipid_IDs = tuple([i for i in ('PL', 'FFA', 'MAG', 'DAG', 'TAG') if i in thermo.chemicals])
    carbs = Stream(None, thermo=thermo)
    fiber = Stream(None, thermo=thermo)
    lipid = Stream(None, TAG=1.0 - PL_fraction - FFA_fraction, units='kg/hr', thermo=thermo)
    if PL_fraction: lipid.imass['PL'] = PL_fraction
    if FFA_fraction: lipid.imass['FFA'] = FFA_fraction
    
    carbs.imol[carbs_IDs] = stream.imol[carbs_IDs]
    fiber.imol[fiber_IDs] = stream.imol[fiber_IDs]
    
    # Mass property arrays
    carbs_mass = stream.imass[carbs_IDs]
    fiber_mass = stream.imass[fiber_IDs]
    lipid_mass = lipid.imass[lipid_IDs]
    
    
    # Net weight
    carbs_massnet = carbs_mass.sum()
    fiber_massnet = fiber_mass.sum()
    lipid_massnet = lipid_mass.sum()
    
    # Heats of combustion per kg
    LHV_carbs_kg = carbs.LHV / carbs_massnet
    LHV_lipid_kg = lipid.LHV / lipid_massnet
    LHV_lipid_over_carbs = LHV_lipid_kg / LHV_carbs_kg
    
    # Relative composition
    r_mass_lipid = lipid_mass / lipid_massnet
    r_mass_carbs = carbs_mass / carbs_massnet
    r_mass_fiber = fiber_mass / fiber_massnet
    F_mass = stream.F_mass
    z_mass_lipid = lipid_fraction
    z_dry = 0.3
    z_mass_lipid = z_mass_lipid * z_dry 
    z_mass_carbs = z_mass_carbs_baseline - z_mass_lipid * LHV_lipid_over_carbs
    z_mass_fiber = z_dry - z_mass_carbs - z_mass_lipid - z_mass_ash_baseline - z_mass_solids_baseline
    imass = stream.imass
    imass['Water'] = F_mass * 0.7
    imass['Ash'] = F_mass * z_mass_ash_baseline
    imass['Solids'] = F_mass * z_mass_solids_baseline
    imass[lipid_IDs] = r_mass_lipid * z_mass_lipid * F_mass
    imass[carbs_IDs] = r_mass_carbs * z_mass_carbs * F_mass
    imass[fiber_IDs] = r_mass_fiber * z_mass_fiber * F_mass
    if any(stream.mol < 0):
        raise ValueError(f'lipid cane oil composition of {z_mass_lipid/z_dry*100:.0f}% dry weight is infeasible')

def get_lipid_fraction(stream=None):
    if not stream: stream = f.stream.lipidcane
    imass = stream.imass
    F_dry_mass = stream.F_mass - imass['Water']
    return imass['Lipid'] / F_dry_mass