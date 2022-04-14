# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import flexsolve as flx
from scipy.optimize import minimize_scalar
import numpy as np
from thermosteam.exceptions import InfeasibleRegion
from biorefineries.lipidcane import set_lipid_fraction as set_oil_fraction
import flexsolve as flx
import thermosteam as tmo
import biosteam as bst

__all__ = ('OilExtractionSpecification', 'MockExtractionSpecification')

class MockExtractionSpecification:
    
    def load_bagasse_efficiency(self, bagasse_efficiency):
        pass
      
    def load_efficiency(self, efficiency):
        pass
    
    def load_oil_content(self, oil_content):
        pass
    
    def load_specifications(self, 
            efficiency=None,
            bagasse_efficiency=None,
            oil_content=None,
        ):
        pass

class OilExtractionSpecification:
    """
    Create a OilExtractionSpecification object for setting process 
    specifications related to oil extraction.
    
    Parameters
    ----------
    system : System
        System associated to feedstocks.
    feedstocks : Stream
        Oilcane feedstocks.
    isplit_efficiency : ChemicalIndexer
        Defines extraction efficiency as a material split.
    oil_content : float 
        Oil content of feedstocks [dry wt. %].
        
    Notes
    -----
    Use the `load_specifications` method to load specifications.
    
    """
    
    __slots__ = (
        'system',
        'feedstocks',
        'isplit_efficiency',
        'isplit_bagasse_efficiency',
        'efficiency',
        'bagasse_efficiency',
        'oil_content',
        'locked_oil_content',
        'isplit_efficiency_is_reversed',
        'FFA_content',
        'PL_content',
    )
    
    def __init__(self, system, feedstocks, isplit_efficiency, isplit_bagasse_efficiency=None, 
                 oil_content=0.1, FFA_content=0.1, PL_content=0.1):
        self.system = system #: [System] System associated to feedstocks
        self.feedstocks = feedstocks #: Stream Oil feedstocks
        self.isplit_efficiency = isplit_efficiency #: [SplitIndexer] Defines extraction efficiency as a material split.
        self.isplit_bagasse_efficiency = isplit_bagasse_efficiency #: [SplitIndexer] Defines extraction efficiency from bagasse as a material split.
        self.oil_content = oil_content #: [float] Oil content of feedstocks [dry wt. %].
        self.locked_oil_content = False
        self.FFA_content = FFA_content
        self.PL_content = PL_content
        
    def MFPP(self):
        return self.system.TEA.solve_price(self.feedstocks)
    
    def load_oil_content(self, oil_content):
        if self.locked_oil_content: return
        for i in self.feedstocks: 
            set_oil_fraction(
                oil_content, i, PL_fraction=self.PL_content, 
                FFA_fraction=self.FFA_content
            )
        self.oil_content = oil_content
      
    def load_efficiency(self, efficiency):
        if self.isplit_efficiency is None: return
        self.efficiency = efficiency
        self.isplit_efficiency['Oil'] = 1. - efficiency
    
    def load_bagasse_efficiency(self, bagasse_efficiency):
        if self.isplit_bagasse_efficiency is None: return
        self.bagasse_efficiency = self.isplit_bagasse_efficiency['Oil'] = bagasse_efficiency
    
    def load_specifications(self, 
            efficiency=None,
            bagasse_efficiency=None,
            oil_content=None,
        ):
        """
        Load oil extraction specifications.

        Parameters
        ----------
        efficiency : float, optional
            Oil extraction efficiency.
        bagasse_efficiency : float, optional
            Oil extraction efficiency from bagasse.
        oil_content : float, optional
            Oil content of feedstocks [dry wt. %].

        """
        if efficiency is None: efficiency = self.efficiency
        if oil_content is None: oil_content = self.oil_content
        self.load_efficiency(efficiency)
        self.load_oil_content(oil_content)
        self.load_bagasse_efficiency(bagasse_efficiency)
    