#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Fri Jun 26 11:48:00 2020

References:
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

@author: yalinli_cabbi
"""


# %% Setup

import numpy as np
import thermosteam as tmo
from biorefineries.ethanol_adipic._chemicals import chems
from biorefineries.lactic import (
    CEPCI, # Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
    get_feedstock_flow, dry_composition, baseline_feedflow,
    compute_COD,
    splits_df
    )


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


