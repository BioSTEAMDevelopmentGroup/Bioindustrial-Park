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


__all__ = (
    'run_sugarcane_readme',
    'run_lipidcane_readme',
    'run_cornstover_readme',
    'run_laos_readme',
    'run_lactic_readme',
    'run_ethanol_adipic_readme',
)

get_readme = lambda module: os.path.join(os.path.dirname(module.__file__), 'README.rst')

def run_sugarcane_readme():
    from biorefineries import sugarcane as sc
    testfile(get_readme(sc))
    
def run_lipidcane_readme():
    from biorefineries import lipidcane as lc
    testfile(get_readme(lc))

def run_cornstover_readme():
    from biorefineries import cornstover as cs
    testfile(get_readme(cs))

def run_laos_readme():
    from biorefineries import LAOs as laos
    testfile(get_readme(laos))
    
def run_lactic_readme():
    from biorefineries import lactic
    testfile(get_readme(lactic))
    
def run_ethanol_adipic_readme():
    from biorefineries import ethanol_adipic
    testfile(get_readme(ethanol_adipic))

if __name__ == '__main__':
    run_sugarcane_readme()
    run_lipidcane_readme()
    run_cornstover_readme()
    run_laos_readme()
    run_lactic_readme()
    run_ethanol_adipic_readme()
