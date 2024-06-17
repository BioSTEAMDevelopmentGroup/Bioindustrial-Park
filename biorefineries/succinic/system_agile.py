#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat
"""

import biosteam as bst
from warnings import filterwarnings; filterwarnings('ignore')
from biorefineries.succinic import system_sc

s = system_sc.s
u = system_sc.u

# Create the agile system object
agile_sys = bst.AgileSystem()

# Create sugarcane system
sc_succinic_sys = system_sc.succinic_sys

# Create new soghum system; which is exactly the same as sugarcane but with a different feedstock.
# Systems remember the configuration they had when they were created, which makes it
# possible to change configurations without affecting old systems.
sugarcane = u.U201.ins[0]
u.U201.ins[0] = None # Remove stream without warnings
u.U201.ins[0] = sorghum = bst.Stream(
    'sorghum',
    Water=0.7,
    Glucose=0.0111,
    Sucrose=0.126,
    Ash=0.006,
    Cellulose=0.0668,
    Hemicellulose=0.0394,
    Lignin=0.0358,
    Solids=0.015,
    total_flow=333334,
    price=s.sugarcane.price,
    units='kg/hr',
)
sorghum_sys = bst.System.from_units('sorghum_sys', sc_succinic_sys.units)

# # Define operation parameters.
# @agile_sys.operation_parameter
# def set_fermentation_efficiency(X):
#     u.R302.X = X

# If mode dependent, the mode is also passed to the setter when the agile system is simulated.
@agile_sys.operation_parameter(mode_dependent=True)
def set_feedstock_flow_rate(flow_rate, mode):
    # We can define a "feedstock" attribute later.
    mode.feedstock.F_mass = flow_rate


# Note how the name of the parameter defaults to the name used in the function (e.g., X, flow_rate).
sugarcane_mode = agile_sys.operation_mode(
    system=sc_succinic_sys, operating_hours=24*200, 
    # X=0.90, 
    feedstock=sugarcane,
    flow_rate=96000.,
)

# Assume fermentation of sorghum is less efficient for demonstration purposes.
# Assume a similar flow rate is available for sorghum.
sorghum_mode = agile_sys.operation_mode(
    system=sorghum_sys, operating_hours=24*60, 
    # X=0.85,
    feedstock=sorghum,
    flow_rate=96000.,
)