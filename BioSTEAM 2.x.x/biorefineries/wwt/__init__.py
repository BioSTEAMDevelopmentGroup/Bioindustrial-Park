#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

# =============================================================================
# Initializing here to avoid circular importing
# =============================================================================

# Units of measure
from thermosteam.units_of_measure import AbsoluteUnitsOfMeasure as auom
_lb_per_kg = 0.4536 # auom('lb').conversion_factor('kg')
_GDP_2007to2016 = 1.160

# Harmonized prices
# References
# ----------
# [1] Hossain et al. Techno-Economic Evaluation of Heat Integrated
# Second Generation Bioethanol and Furfural Coproduction.
# Biochemical Engineering Journal 2019, 144, 89–103.
# https://doi.org/10.1016/j.bej.2019.01.017.

# [2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
# Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
# NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
# https://doi.org/10.2172/1483234.

# [3] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
# Valorization of Dilute Organic Carbon Waste Streams.
# Energy Environ. Sci. 2016, 9 (3), 1102–1112.
# https://doi.org/10.1039/C5EE03715H.

# .. note::
#
#     The cornstover biorefinery uses 2011 USD whereas sugarcane and lipidcane
#     biorefineries use 2013 USD, thus ideally prices in the the `new_price` dct
#     should be adjusted accordingly. However, since the calculated MESPs are only
#     used for comparison between biorefineries with the new wastewater treatment
#     process vs. without/with the original process, this will not affect the
#     conclusions.

#!!! Want to adjust to $2016
new_price = { # $/kg unless otherwise noted
    'Wastewater': -0.03, # ref [1], negative value for cost from product,
    'NaOCl': 0.14, # $/L
    'CitricAcid': 0.22, # $/L
    'Bisulfite': 0.08, # $/L
    'Caustics': 0.2627, # la._settings.price['NaOH]/2 as the caustic is 50% NaOH/water
    'Polymer': 2.6282 / _lb_per_kg / _GDP_2007to2016, # ref [2]
    }

ethanol_density_kggal = 2.9867 # cs.ethanol_density_kggal
import biosteam as bst
def print_MESP(ethanol, tea, tea_name=''):
    bst.settings.set_thermo(ethanol.chemicals)
    tea.system.simulate()
    ethanol.price = tea.solve_price(ethanol)
    ethanol_price_gal = ethanol.price * ethanol_density_kggal
    pre_tea = ' of ' if tea_name else ''
    print(f'MESP{pre_tea}{tea_name} is ${ethanol_price_gal:.2f}/gal.')
    return ethanol_price_gal


# Path
import os
wwt_path = os.path.dirname(__file__)
results_path = os.path.join(wwt_path, 'results')
del os


# %%

# =============================================================================
# Importing from this modules
# =============================================================================

# Importing order matters
from . import (
    _chemicals,
    utils,
    _internal_circulation_rx,
    _wwt_pump,
    _polishing_filter,
    _membrane_bioreactor,
    _sludge_handling,
    _wwt_sys,
    _lca,
    )

# These orders don't matter
from .utils import *
from ._chemicals import *
from ._internal_circulation_rx import *
from ._wwt_pump import *
from ._polishing_filter import *
from ._membrane_bioreactor import *
from ._sludge_handling import *
from ._wwt_sys import *
from ._lca import *


__all__ = (
    # Unit conversion
    'auom', 'ethanol_density_kggal', 'new_price',
    # Path
    'wwt_path', 'results_path',
    *_chemicals.__all__,
    *_internal_circulation_rx.__all__,
    *_wwt_pump.__all__,
    *_polishing_filter.__all__,
    *_membrane_bioreactor.__all__,
    *_sludge_handling.__all__,
    *_wwt_sys.__all__,
    *_lca.__all__,
    *utils.__all__,
    )