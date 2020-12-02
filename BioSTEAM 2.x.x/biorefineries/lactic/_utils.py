#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu>,
# Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Tue May 19 14:17:17 2020

Modified from the biorefineries constructed in [1] and [2] for the production of
lactic acid from lignocellulosic feedstocks

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040
    
[2] Li et al., Tailored Pretreatment Processes for the Sustainable Design of
    Lignocellulosic Biorefineries across the Feedstock Landscape. Submitted,
    2020.

@author: yalinli_cabbi
"""


# %% Setup

import numpy as np
import pandas as pd
import thermosteam as tmo
from biorefineries.lactic._chemicals import chems
_kg_per_ton = 907.18474

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
CEPCI = {1997: 386.5,
         1998: 389.5,
         2007: 525.4,
         2009: 521.9,
         2010: 550.8,
         2011: 585.7,
         2012: 584.6,
         2013: 567.3,
         2014: 576.1,
         2016: 541.7}


# %% 

# =============================================================================
# Function to get feedstock flow by giving dry weight composition and moisture content
# =============================================================================

def get_feedstock_flow(dry_composition, moisture_content, dry_flow):
    dry_array = chems.kwarray(dry_composition)
    wet_flow = dry_flow / (1-moisture_content)
    moisture_array = chems.kwarray(dict(Water=moisture_content))
    feedstock_flow = wet_flow * (dry_array*(1-moisture_content)+moisture_array)
    return feedstock_flow

dry_composition = dict(
    Glucan=0.3505, Xylan=0.1953, Lignin=0.1576, Ash=0.0493, Acetate=0.0181,
    Protein=0.0310, Arabinan=0.0238, Galactan=0.0143, Mannan=0.0060, 
    Sucrose=0.0077, Extractives=0.1465, SuccinicAcid=0)

moisture_content = 0.2
dry_feedstock_flow = 2205 * _kg_per_ton / 24     
baseline_feedflow = get_feedstock_flow(dry_composition, moisture_content, 
                                       dry_feedstock_flow)


# %%

# =============================================================================
# Function to calculate lactic acid titer in g/L (kg/m3) of a stream
# =============================================================================

def compute_lactic_titer(stream, V=None):
    lactic_mass = 2*stream.imol['CalciumLactate']*chems.LacticAcid.MW \
        + stream.imass['LacticAcid']
    if V:
        vol = V
    else: vol = stream.F_vol
    if vol != 0:
        titer = lactic_mass / vol
    else: titer = 0
    return titer

def set_yield(lactic_yield, R301, R302):
    if lactic_yield > 1:
        raise ValueError(f'Lactic acid yield of {lactic_yield:.2f} is infeasible')
    R301_X = R301.cofermentation_rxns.X
    R301_X[0] = R301_X[3] = lactic_yield
    R301_X[2] = R301_X[5] = R301._X[2]
    R301_X[1] = R301_X[4] = min(R301._X[1], 1-1e-6-R301_X[0]-R301_X[2])
    if R301_X[1] < 0:
        R301_X[2] = R301_X[5] = min(R301._X[2], -(R301_X[1]+1e-6), 1-1e-6-R301_X[0])
        R301_X[1] = R301_X[4] = 0
    R302_X = R302.cofermentation_rxns.X
    R302_X[0] = R302_X[3] = R301_X[0] * R302.ferm_ratio
    R302_X[1] = R302_X[4] = min(R301_X[1]*R302.ferm_ratio, 1-1e-6-R302_X[0]-R302_X[2])


# %% 

# =============================================================================
# Functions to compute chemical loading and adjust recycle flows to maintain
# a certain ratio for Esterification and Hydrolysis reactor
# =============================================================================

def compute_extra_chemical(feed, recycle, reactants_ID, chemical_ID, ratios):
    reactants_in_feed = feed.imol[reactants_ID]
    reactants_in_recycle = recycle.imol[reactants_ID]
    chemical_needed = (ratios*(reactants_in_feed+reactants_in_recycle)).sum()
    chemical_extra = (feed.imol[chemical_ID]+recycle.imol[chemical_ID]) - chemical_needed
    return chemical_extra

def adjust_recycle(feed, recycle, reactants_ID, chemical_ID, ratios):
    feed_chemical_needed = (feed.imol[reactants_ID]*ratios).sum() \
        - feed.imol[chemical_ID]
    
    recycle_chemical_extra = recycle.imol[chemical_ID] \
        - (recycle.imol[reactants_ID]*ratios).sum()
    
    split = feed_chemical_needed / recycle_chemical_extra
    effluent = feed.copy()
    recycle_recycled = recycle.copy()
    recycle_recycled.mol *= split
    recycle_discarded = recycle.copy()
    recycle_discarded.mol *= (1 - split)
    effluent.mix_from([feed, recycle_recycled])
    
    return effluent, recycle_discarded


# %% 

# =============================================================================
# Function to compute chemical oxygen demand (COD) of a given stream based on:
#     C_nH_aO_bN_c + (n+a/4-b/2-3/4c)O2 -> nCO2 + (a/2-3/2c)H2O +cNH3
# =============================================================================

def compute_COD(IDs, stream):
    unit_COD = []
    if not iter(IDs):
        raise TypeError(f'{IDs.__class__} is not iterable')
    if isinstance(IDs, str):
        IDs = (IDs,)
    for i in IDs:
        if i not in stream.chemicals.IDs: continue
        chemical = getattr(stream.chemicals, i)
        atoms = {}
        for j in ('C', 'H', 'O', 'N'): atoms[j] = 0
        atoms.update(chemical.atoms)
        COD_ratio = atoms['C'] + atoms['H']/4 - atoms['O']/2 - 3/4*atoms['N']
        unit_COD.append(COD_ratio)
    unit_COD = np.asarray(unit_COD)
    return unit_COD.sum()*32


# %% 

# =============================================================================
# Function to find the split ratios for Splitters, assume 0 for chemicals not
# specified in splits
# =============================================================================

def find_split(IDs, flow0, flow1, chemical_groups):
    # Add 1e-6 to avoid flow0 and flow1 both being 0
    flow0 = np.asarray(flow0) + 1e-6
    flow1 = np.asarray(flow1) + 1e-6    
    splits = flow0/(flow0 + flow1)
    thermo = tmo.settings.get_thermo()
    chemicals = thermo.chemicals
    array = np.zeros(chemicals.size)  
    for ID, split in zip(IDs, splits):
        if ID in chemical_groups:
            array[chemicals.get_index(chemical_groups[ID])] = split
        else:
            array[chemicals.index(ID)] = split
    # WWTsludge is removed from the cell mass group 
    array[chemicals.index('WWTsludge')] = array[chemicals.index('FermMicrobe')]
    return array

IDs = ('Ethanol', 'H2O', 'Glucose', 'Xylose', 'OtherSugars',
    'SugarOligomers', 'OrganicSolubleSolids', 'InorganicSolubleSolids', 'Ammonia', 'AceticAcid', 
    'SulfuricAcid', 'Furfurals', 'OtherOrganics', 'CO2', 'CH4',
    'O2', 'N2', 'COSOxNOxH2S', 'Glucan', 'Xylan', 
    'OtherStructuralCarbohydrates', 'Acetate', 'Lignin', 'Protein', 'CellMass',
    'OtherInsolubleSolids')

streams = {}

streams['stream_535'] = (177, 329030, 502, 1022, 2094,
                         1552, 15808, 2513, 0, 0,
                         0, 513, 1348, 0, 0,
                         0, 0, 0, 25, 8,
                         2, 0, 250, 69, 19,
                         92)

streams['stream_571'] = (6, 12797, 19, 49, 81,
                         60, 612, 97, 0, 0,
                         0, 19, 52, 0, 0,
                         1, 1, 0, 1230, 415,
                         94, 0, 12226, 3376, 925,
                         4489)

streams['stream_611'] = (15, 356069, 42, 85, 175,
                         130, 2387, 110, 633, 5, 
                         0, 70, 113, 181, 3, 
                         1, 0, 300, 6, 2, 
                         0, 0, 64, 18, 280,
                         23)
streams['stream_612'] = (1, 27158, 3, 7, 13,
                         10, 182, 8, 48, 0,
                         0, 5, 9, 14, 0,
                         0, 0, 23, 19, 6,
                         1, 0, 186, 51, 813,
                         68)

streams['stream_616'] = (1, 109098, 3, 6, 13,
                         9, 187, 1068, 46, 0,
                         0, 5, 8, 14, 0, 
                         1, 1, 31, 1, 0,
                         0, 0, 13, 3, 80, 
                         5)

streams['stream_623'] = (0, 7708, 0, 0, 1,
                         1, 13, 75, 3, 0,
                         0, 0, 1, 1, 0,
                         0, 0, 2, 25, 8,
                         2, 0, 250, 52, 1523,
                         92)

streams['stream_624'] = (0, 381300, 0, 1, 1,
                         1, 79, 4828, 3, 0,
                         0, 0, 1, 6, 0,
                         3, 5, 44, 0, 0,
                         0, 0, 0, 0, 0,
                         0)

streams['stream_625'] = (1, 2241169, 2, 3, 7,
                         6, 466, 28378, 16, 0,
                         0, 3, 7, 38, 0,
                         17, 32, 259, 194, 65,
                         15, 0, 1925, 90, 19778,
                         707)

splits_df = pd.DataFrame.from_dict(streams)
splits_df.index = IDs




