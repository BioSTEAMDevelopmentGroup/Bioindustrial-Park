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
    lactic,
    ethanol_adipic,
    )

__all__ = (
    'test_sugarcane_readme',
    'test_lipidcane_readme',
    'test_cornstover_readme',
    'test_laos_readme',
    'test_lactic_readme',
    'test_ethanol_adipic_readme',
)

get_readme = lambda module: os.path.join(os.path.dirname(module.__file__), 'README.rst')
test_sugarcane_readme = lambda: testfile(get_readme(sc))
test_lipidcane_readme = lambda: testfile(get_readme(lc))
test_cornstover_readme = lambda: testfile(get_readme(cs))
test_laos_readme = lambda: testfile(get_readme(laos))
test_lactic_readme = lambda: testfile(get_readme(lactic))
test_ethanol_adipic_readme = lambda: testfile(get_readme(ethanol_adipic))

if __name__ == '__main__':
    test_sugarcane_readme()
    test_lipidcane_readme()
    test_cornstover_readme()
    test_laos_readme()
    test_lactic_readme()
    test_ethanol_adipic_readme()
