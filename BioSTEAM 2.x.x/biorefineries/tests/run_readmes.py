# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

import os
from doctest import testfile
from biorefineries import (
    sugarcane as sc,
    lipidcane as lc,
    cornstover as cs,
    LAOs as laos,
    ethanol_adipic,
    lactic,
    )

get_readme = lambda module: os.path.join(os.path.dirname(module.__file__), 'README.rst')
testfile(get_readme(sc))
testfile(get_readme(lc))
testfile(get_readme(cs))
testfile(get_readme(laos))
testfile(get_readme(ethanol_adipic))
testfile(get_readme(lactic))