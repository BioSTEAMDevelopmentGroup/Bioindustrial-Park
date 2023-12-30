#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 09:43:38 2023

@author: wenjun
"""


# %% Setup

import numpy as np
import thermosteam as tmo
import biorefineries.SAF.chemicals as chems


# %%

# =============================================================================
# Function to convert ethanol fraction from weight basis to mol basis in an
# ethanol-water mixture
# =============================================================================

def convert_ethanol_wt_2_mol(wt):
    Ethanol_MW = chems.Ethanol.MW
    Water_MW = chems.Water.MW
    return (wt/Ethanol_MW) / (wt/Ethanol_MW + (1-wt)/Water_MW)