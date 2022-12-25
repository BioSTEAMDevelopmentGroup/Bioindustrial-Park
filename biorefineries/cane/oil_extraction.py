# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

__all__ = (
    'OilExtractionSpecification', 
    'MockExtractionSpecification',
)

class MockExtractionSpecification:
    
    def load_saccharification_oil_recovery(self, bagasse_recovery):
        pass
      
    def load_crushing_mill_oil_recovery(self, recovery):
        pass
    
    def load_specifications(self, 
            recovery=None,
            bagasse_recovery=None,\
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
    isplit_crushing_mill : ChemicalIndexer
        Defines extraction recovery as a material split.
        
    Notes
    -----
    Use the `load_specifications` method to load specifications.
    
    """
    
    __slots__ = (
        'feedstocks',
        'isplit_crushing_mill',
        'isplit_saccharification',
        'crushing_mill_oil_recovery',
        'saccharification_oil_recovery',
        'isplit_recovery_is_reversed',
    )
    
    def __init__(self, system, isplit_crushing_mill, isplit_saccharification=None):
        self.isplit_crushing_mill = isplit_crushing_mill #: [SplitIndexer] Defines extraction recovery as a material split.
        self.isplit_saccharification = isplit_saccharification #: [SplitIndexer] Defines extraction recovery from bagasse as a material split.
        if isplit_crushing_mill is not None: 
            self.crushing_mill_oil_recovery = isplit_crushing_mill['Oil'].mean()
        else:
            self.crushing_mill_oil_recovery = None
        if isplit_saccharification is not None: 
            self.saccharification_oil_recovery = 1 - isplit_saccharification['Oil'].mean()
        else:
            self.saccharification_oil_recovery = None
      
    def load_crushing_mill_oil_recovery(self, recovery):
        if self.isplit_crushing_mill is None: return
        self.crushing_mill_oil_recovery = recovery
        self.isplit_crushing_mill['Oil'] = 1. - recovery
    
    def load_saccharification_oil_recovery(self, recovery):
        if self.isplit_saccharification is None: return
        self.saccharification_oil_recovery = recovery
        self.isplit_saccharification['Oil'] = 1. - recovery
    
    def load_specifications(self, 
            crushing_mill_oil_recovery=None,
            saccharification_oil_recovery=None,
        ):
        """
        Load oil extraction specifications.

        Parameters
        ----------
        recovery : float, optional
            Oil extraction recovery.
        bagasse_recovery : float, optional
            Oil extraction recovery from bagasse.
        oil_content : float, optional
            Oil content of feedstocks [dry wt. %].

        """
        if crushing_mill_oil_recovery is None: 
            crushing_mill_oil_recovery = self.crushing_mill_oil_recovery
        if saccharification_oil_recovery is None:
            saccharification_oil_recovery = self.saccharification_oil_recovery
        self.load_crushing_mill_oil_recovery(crushing_mill_oil_recovery)
        self.load_saccharification_oil_recovery(saccharification_oil_recovery)
    
