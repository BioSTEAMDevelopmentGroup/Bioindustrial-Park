#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the oilcane biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/oilcane
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biosteam import Flowsheet, main_flowsheet
from biorefineries.oilcane import create_chemicals, load_process_settings, price
from biorefineries.wwt import add_wwt_chemicals, create_wastewater_system




# %%

# =============================================================================
# Function to make the system
# =============================================================================

new_oc_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()
bst.settings.set_thermo(new_oc_chems)