#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

from biorefineries.cornstover import CellulosicEthanolTEA

__all__ = ('EthanolAdipicTEA',)

class EthanolAdipicTEA(CellulosicEthanolTEA):

    # For uncertainty analysis
    _TCI_ratio_cached = 1