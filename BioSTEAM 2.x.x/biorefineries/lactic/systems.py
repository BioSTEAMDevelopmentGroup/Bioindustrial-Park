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

'''
References
----------
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of 
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic 
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764; 
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf
[2] Kuo et al., Production of Optically Pure L-Lactic Acid from Lignocellulosic
    Hydrolysate by Using a Newly Isolated and d-Lactate Dehydrogenase
    Gene-Deficient Lactobacillus Paracasei Strain.
    Bioresource Technology 2015, 198, 651–657.
    https://doi.org/10.1016/j.biortech.2015.09.071.
[3] Aden et al., Process Design Report for Stover Feedstock: Lignocellulosic
    Biomass to Ethanol Process Design and Economics Utilizing Co-Current Dilute
    Acid Prehydrolysis and Enzymatic Hydrolysis for Corn Stover; NREL/TP-510-32438;
    National Renewable Energy Lab (NREL), 2002.
    https://doi.org/10.2172/1218326.

Naming conventions:
    D = Distillation column
    E = Evaporator
    F = Flash tank
    H = Heat exchange
    M = Mixer
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units
    PS = Process specificiation, not physical units, but for adjusting streams

Processes:
    100: Feedstock preprocessing
    200: Pretreatment
    300: Conversion
    400: Separation
    500: Wastewater treatment
    600: Facilities

'''


# %%

import biosteam as bst
from biorefineries.lactic._chemicals import chems
from biorefineries.lactic._settings import CFs
from biorefineries.lactic._utils import find_split, splits_df
from biorefineries.lactic._chemicals import chemical_groups

flowsheet = bst.Flowsheet('')
feedstock = flowsheet.stream.feedstock
U101 = flowsheet.unit.feedstock
# Feedstock flow rate in dry U.S. ton per day
get_flow_tpd = lambda: \
    (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185*(1-U101.diversion_to_CHP)

# 1 is water, changed by moisture content rather than using data from ref [1]
cell_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
cell_solid = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
cell_filtrate = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
cell_split = find_split(cell_index, cell_solid, cell_filtrate, chemical_groups)

# Moisture content (20%) and gypsum removal (99.5%) on Page 24 of ref [3]
gypsum_index = cell_index + ['Gypsum']
gypsum_solid = cell_solid + [0.995]
gypsum_filtrate = cell_filtrate + [0.005]
gypsum_split = find_split(gypsum_index, gypsum_solid, gypsum_filtrate, chemical_groups)
















