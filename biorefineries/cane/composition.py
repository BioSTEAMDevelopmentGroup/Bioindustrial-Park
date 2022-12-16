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
import thermosteam as tmo
from thermosteam import Stream
from biosteam import main_flowsheet as f
import numpy as np

__all__ = (
    'set_lipid_fraction',
    'get_lipid_fraction',
    'set_sugarcane_composition',
    'convert_fiber_to_lignocelluosic_components',
    'set_composition',
    'get_composition'
)

def set_lipid_fraction(lipid_fraction, stream=None,
                       PL_fraction=0.,
                       FFA_fraction=0.,
                       z_mass_carbs_baseline=0.149,
                       z_mass_ash_baseline=0.006,
                       z_mass_solids_baseline=0.015,
                       z_mass_water_baseline=0.7):
    """Adjust composition of lipid cane to achieve desired oil fraction (dry weight)
    assuming that an increase in lipid content is accompanied by a decrease in sugar
    content using an energy balance and an increase in fiber using a mass balance 
    (lipid is more energy dense than sugar)."""
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
    z_dry = 1. - z_mass_water_baseline
    z_mass_lipid = z_mass_lipid * z_dry 
    z_mass_carbs = z_mass_carbs_baseline - z_mass_lipid * LHV_lipid_over_carbs
    z_mass_fiber = z_dry - z_mass_carbs - z_mass_lipid - z_mass_ash_baseline - z_mass_solids_baseline
    imass = stream.imass
    imass['Water'] = F_mass * z_mass_water_baseline
    imass['Ash'] = F_mass * z_mass_ash_baseline
    imass['Solids'] = F_mass * z_mass_solids_baseline
    imass[lipid_IDs] = r_mass_lipid * z_mass_lipid * F_mass
    imass[carbs_IDs] = r_mass_carbs * z_mass_carbs * F_mass
    imass[fiber_IDs] = r_mass_fiber * z_mass_fiber * F_mass
    if any(stream.mol < 0):
        raise ValueError(f'lipid cane oil composition of {z_mass_lipid/z_dry*100:.0f}% dry weight is infeasible')

def set_sugarcane_composition(stream, water, fiber, sugar):
    chemicals = stream.chemicals
    if 'Sugar' not in chemicals:
        IDs = ('Sucrose', 'Glucose')
        chemicals.define_group('Sugar', IDs, composition=stream.imass[IDs], wt=True)
    if 'Fiber' not in chemicals:
        IDs = ('Cellulose', 'Hemicellulose', 'Lignin')
        chemicals.define_group('Fiber', IDs, composition=stream.imass[IDs], wt=True)
    if 'Other' not in chemicals:
        IDs = ('Ash', 'Solids')
        chemicals.define_group('Other', IDs, composition=stream.imass[IDs], wt=True)
    F_mass = stream.F_mass
    other = 1 - water - fiber - sugar
    assert other > 0, (water, fiber, sugar, other)
    stream.imass['Water', 'Fiber', 'Sugar', 'Other'] = F_mass * np.array([water, fiber, sugar, other])

def convert_fiber_to_lignocelluosic_components(stream, ignore_acetate=False):
    chemicals = stream.chemicals
    if chemicals is convert_fiber_to_lignocelluosic_components.last_chemicals:
        prxn = convert_fiber_to_lignocelluosic_components.last_reaction
    else:
        cellulose_rxn = tmo.Reaction('Cellulose -> Glucan', 'Cellulose', 1.0,
                                     basis='wt', chemicals=chemicals)
        cellulose_rxn.basis = 'mol'
        # Bagasse composition https://www.sciencedirect.com/science/article/pii/S0144861710005072
        # South american; by HPLC
        # Glucan: 41.3%
        # Xylan: 24.9%
        # Galactan: 0.6%
        # Arabinan: 1.7%
        # Lignin: 23.2%
        # Acetyl: 3.0%
        if ignore_acetate:
            hemicellulose_rxn = tmo.Reaction(
                '27.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan', 'Hemicellulose', 
                1.0, basis='wt', chemicals=chemicals
            )
        else:
            hemicellulose_rxn = tmo.Reaction(
                '30.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan + 3 Acetate', 'Hemicellulose', 
                1.0, basis='wt', chemicals=chemicals
            )
        hemicellulose_rxn.basis = 'mol'
        convert_fiber_to_lignocelluosic_components.last_chemicals = chemicals
        convert_fiber_to_lignocelluosic_components.last_reaction = prxn = tmo.ParallelReaction(
            [cellulose_rxn, hemicellulose_rxn]
        )
    prxn(stream)

convert_fiber_to_lignocelluosic_components.last_chemicals = None

def set_composition(
        stream,
        # By weight (e.g., g water / g cane)
        moisture,
        # By dry weight (e.g., g lipid / g dry cane)
        lipid, 
        fiber, 
        solids,
        ash,
        # By lipid weight (e.g., g TAG / g Lipid)
        TAG,
        FFA,
        PL,
    ):
    """Set the composition of a sugarcane stream by dry weight 
    (lipid, fiber, ash, and solids), total moisture content (moisture),
    and lipid weight (TAG, FFA, PL)."""
    imass = stream.imass
    dry_content = 1 - moisture
    F_mass = stream.F_mass
    total_dry = F_mass * dry_content
    total_lipid = total_dry * lipid
    imass['TAG'] = total_lipid * TAG
    imass['FFA'] = total_lipid * FFA
    imass['PL'] = total_lipid * PL
    imass['Fiber'] = total_dry * fiber
    imass['Water'] = F_mass * moisture
    imass['Ash'] = total_dry * ash
    imass['Solids'] = total_dry * solids
    imass['Sugar'] = 0.
    sugar = F_mass - stream.F_mass
    if sugar < 0.: raise ValueError('composition is infeasibile')
    imass['Sugar'] = sugar

def get_composition(stream):
    total = stream.F_mass
    imass = stream.imass
    water = imass['Water']
    total_dry = total - water
    total_lipid = imass['Lipid']
    return dict(
        moisture=water / total,
        lipid=total_lipid / total_dry,
        sugar=imass['Sugar'] / total_dry,
        fiber=imass['Fiber'] / total_dry,
        solids=imass['Solids'] / total_dry,
        ash=imass['Ash'] / total_dry,
        TAG=imass['TAG'] / total_lipid if total_lipid else 0.,
        FFA=imass['FFA'] / total_lipid if total_lipid else 0.,
        PL=imass['PL'] / total_lipid if total_lipid else 0.,
    )

def get_lipid_fraction(stream=None):
    if not stream: stream = f.stream.lipidcane
    imass = stream.imass
    F_dry_mass = stream.F_mass - imass['Water']
    return imass['Lipid'] / F_dry_mass