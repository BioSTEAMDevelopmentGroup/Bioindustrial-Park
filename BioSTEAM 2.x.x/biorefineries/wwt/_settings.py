#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
Settings for cornstover (cs), sugarcane (sc), and lipidcane (lc) biorefineries.

References
----------
[1] Hossain et al. Techno-Economic Evaluation of Heat Integrated
Second Generation Bioethanol and Furfural Coproduction.
Biochemical Engineering Journal 2019, 144, 89–103.
https://doi.org/10.1016/j.bej.2019.01.017.

[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
https://doi.org/10.2172/1483234.

[3] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
Valorization of Dilute Organic Carbon Waste Streams.
Energy Environ. Sci. 2016, 9 (3), 1102–1112.
https://doi.org/10.1039/C5EE03715H.

.. note::
    The cornstover biorefinery uses 2011 USD whereas sugarcane and lipidcane
    biorefineries use 2013 USD, thus ideallly prices in the the `new_price` dct
    should be adjusted accordingly. However, since the calculated MESPs are only
    used for comparison between biorefineries with the new wastewater treatment
    process vs. without/with the original process, this will not affect the
    conclusions.
'''

__all__ = (
    'cs_price',
    'lc_price',
    'new_price',
    'load_cs_settings',
    'load_sc_settings',
    'load_lc_settings'
    )

from biorefineries import cornstover as cs
from biorefineries.cornstover._process_settings import (
    price as cs_price,
    load_process_settings as load_cs_settings,
    )

from biorefineries.sugarcane._process_settings import \
    load_process_settings as load_sc_settings

from biorefineries.lipidcane._process_settings import (
    price as lc_price,
    load_process_settings as load_lc_settings
    )

# from biorefineries.utils import auom
from utils import auom
_lb_per_kg = auom('lb').conversion_factor('kg')
_GDP_2007to2016 = 1.160

new_price = { # $/kg unless otherwise noted
    'Wastewater': -0.03, # ref [1], negative value for cost from product,
    'NaOCl': 0.14, # $/L
    'CitricAcid': 0.22, # $/L
    'Bisulfite': 0.08, # $/L
    'Caustics': cs.caustic.price,
    'Polymer': 2.6282 / _lb_per_kg / _GDP_2007to2016, # ref [2]
    }