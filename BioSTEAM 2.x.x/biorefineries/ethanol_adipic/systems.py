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
    https://doi.org/10.2172/1218326
[4] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234
'''

import biosteam as bst
from biorefineries.ethanol_adipic._chemicals import chems
from biorefineries.ethanol_adipic._utils import auom
from biorefineries.ethanol_adipic._settings import set_feedstock_price, \
    price, CFs
from biorefineries.ethanol_adipic._processes import (
    create_preprocessing_process,
    )

bst.settings.set_thermo(chems)
bst.CE = 541.7 # year 2016
_kg_per_ton = auom('ton').conversion_factor('kg')


# %%

# =============================================================================
# Different depot systems
# =============================================================================

CPP_flowsheet, CPP_cost = create_preprocessing_process(kind='CPP', with_AFEX=False)
CPP_AFEX_flowsheet, CPP_AFEX_cost = create_preprocessing_process(kind='CPP', with_AFEX=True)
HMPP_flowsheet, HMPP_cost = create_preprocessing_process(kind='HMPP', with_AFEX=False)
HMPP_AFEX_flowsheet, HMPP_AFEX_cost = create_preprocessing_process(kind='HMPP', with_AFEX=True)

CPP_feedstock = CPP_flowsheet.stream.feedstock
CPP_AFEX_feedstock = CPP_AFEX_flowsheet.stream.feedstock
HMPP_feedstock = HMPP_flowsheet.stream.feedstock
HMPP_AFEX_feedstock = HMPP_AFEX_flowsheet.stream.feedstock

# # If want to use the default preprocessing price ($24.35/Mg)
# set_feedstock_price(feedstock)
# # If want to use the price in ref [2], note that the price here is $/dry U.S. ton
# CPP_feedstock.price = price['Feedstock']


# %%

def create_acid_sys(feedstock):
    flowsheet = bst.Flowsheet('Acid')
    s = flowsheet.stream
    u = flowsheet.unit
    
    get_flow_tpd = \
        lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/_kg_per_ton









# Acid:
M601 = bst.units.Mixer('M601', ins=(H201-0, D402_P-0, S401-1, ''))
# (last one is blowdown)






acid_feedstock = CPP_feedstock.copy('acid_feedstock')
acid_flowsheet = create_acid_sys(acid_feedstock)














