# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biorefineries import lipidcane as lc

unit_groups = {
    'Oil separation': (lc.U101, lc.U102, lc.U103, lc.T201, 
                       lc.U201, lc.U202, lc.M201, lc.S201, 
                       lc.T202, lc.H201, lc.T203, lc.P201, 
                       lc.T204, lc.P202, lc.T205, lc.M202, 
                       lc.H202, lc.T206, lc.C201, lc.C202, 
                       lc.P203, lc.T207, lc.T207_2, lc.H203, 
                       lc.S202, lc.T208, lc.C203, lc.F201,
                       lc.U202),
    'Biodiesel production': (lc.T403, lc.P403, lc.R401, lc.C401, lc.R402,
                             lc.PS5, lc.C402, lc.T406, lc.P409, lc.C404, 
                             lc.T407, lc.P410, lc.H402, lc.D401, lc.D402, 
                             lc.H404, lc.P412, lc.T405, lc.P406, lc.C403, 
                             lc.F401, lc.H401, lc.P408, lc.P407, lc.T409, 
                             lc.P405, lc.B401, lc.H403, lc.P411, lc.T401, 
                             lc.P401, lc.T402, lc.P402, lc.T404, lc.P404, 
                             lc.S401),
    'Fermentation': (lc.S301, lc.F301, lc.P306, lc.M301, lc.H301, lc.T305,
                     lc.R301, lc.T301, lc.C301, lc.M302),
    'Ethanol separation': (lc.P301, lc.H302, lc.D302, lc.P302, lc.M303, 
                           lc.D303, lc.H303, lc.U301, lc.H304, lc.T302, 
                           lc.P304, lc.PS4, lc.T303, lc.P305, lc.T304, 
                           lc.D301, lc.T408, lc.P303, lc.M305),
    'Facilities': (lc.BT, lc.CWP, lc.CT, lc.PWC)
}

