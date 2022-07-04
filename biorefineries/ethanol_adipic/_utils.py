#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
References
----------
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf

[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
    https://doi.org/10.2172/1483234

[3] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design,
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty.
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
    https://doi.org/10.1021/acssuschemeng.9b07040
'''


# %% Setup

import numpy as np
import thermosteam as tmo
from biorefineries.ethanol_adipic._chemicals import chems, chemical_groups
from biorefineries.lactic._utils import (
    auom, CEPCI, get_feedstock_flow, dry_composition, baseline_feedflow,
    compute_COD, splits_df
    )

__all__ = ('auom', 'CEPCI', 'get_feedstock_flow', 'dry_composition',
           'baseline_feedflow', 'compute_muconic_titer',
           #'set_yield',
           'compute_COD')

_kg_per_ton = auom('ton').conversion_factor('kg')
ethanol_V = chems.Ethanol.V('l', 298.15, 101325) # molar volume in m3/mol
ethanol_MW = chems.Ethanol.MW
_ethanol_kg_2_gal = auom('gal').conversion_factor('liter')/ethanol_V*ethanol_MW/1e6


# %%

# =============================================================================
# Function to convert ethanol fraction from weight basis to mol basis in an
# ethanol-water mixture
# =============================================================================

def convert_ethanol_wt_2_mol(wt):
    Ethanol_MW = chems.Ethanol.MW
    Water_MW = chems.Water.MW
    return (wt/Ethanol_MW) / (wt/Ethanol_MW + (1-wt)/Water_MW)


# %%

# =============================================================================
# Function to compute titer of muconic acid fermentation in g/L (kg/m3)
# =============================================================================

def compute_muconic_titer(stream, V=None):
    muconic_mass = stream.imol['MonoSodiumMuconate']*chems.MuconicAcid.MW \
        + stream.imass['MuconicAcid']
    if V: vol = V
    else: vol = stream.F_vol
    if vol != 0:
        titer = muconic_mass / vol
    else: titer = 0
    return titer

#!!! Maybe need a function to set yield for muconic acid?


# %%

# =============================================================================
# Function to find the split ratios for Splitters, assume 0 for chemicals not specified in splits,
# stream numbers refer to those in ref [1]
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
    array[chemicals.index('WWTsludge')] = array[chemicals.index('Z_mobilis')]
    return array

# 1 is water, changed by moisture content rather than using data from ref [1]
cell_mass_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
cell_mass_solid = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
cell_mass_filtrate = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
cell_mass_split = find_split(cell_mass_index, cell_mass_solid, cell_mass_filtrate,
                             chemical_groups)

# Anaerobic digestion
AD_split = find_split(splits_df.index,
                      splits_df['stream_611'],
                      splits_df['stream_612'],
                      chemical_groups)

# Membrane bioreactor
MB_split = find_split(splits_df.index,
                      splits_df['stream_624'],
                      splits_df['stream_625'],
                      chemical_groups)