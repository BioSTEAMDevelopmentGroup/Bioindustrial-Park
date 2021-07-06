# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import sys, biosteam

version_components = sys.version.split('.')
PY_MAJOR, PY_MINOR = int(version_components[0]), int(version_components[1])
PY37 = (PY_MAJOR, PY_MINOR) >= (3, 7)
version_components = biosteam.__version__.split('.')
BST_MAJOR, BST_MINOR = int(version_components[0]), int(version_components[1])
BST222 = (BST_MAJOR, BST_MINOR) >= (2, 22)
del sys, version_components, PY_MAJOR, PY_MINOR, BST_MAJOR, BST_MINOR

__version__ = '2.22.0'


