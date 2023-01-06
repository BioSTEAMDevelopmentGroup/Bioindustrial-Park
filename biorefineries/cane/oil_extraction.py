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
    
    def load_bagasse_oil_recovery(self, bagasse_recovery):
        pass
      
    def load_crushing_mill_oil_recovery(self, recovery):
        pass
    
    def load_microbial_oil_recovery(self, recovery):
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
    crushing_mill : ChemicalIndexer
        Defines extraction recovery as a material split.
        
    Notes
    -----
    Use the `load_specifications` method to load specifications.
    
    """
    
    __slots__ = (
        'feedstocks',
        'crushing_mill',
        'pressure_filter',
        'cellmass_centrifuge',
        'crushing_mill_oil_recovery',
        'bagasse_oil_recovery',
        'microbial_oil_recovery',
    )
    
    def __init__(self, system, crushing_mill, pressure_filter=None, cellmass_centrifuge=None, microbial_oil_recovery=None):
        self.crushing_mill = crushing_mill #: [Unit] Defines extraction recovery as a material split.
        self.pressure_filter = pressure_filter #: [Unit] Defines extraction recovery from bagasse as a material split.
        self.cellmass_centrifuge = cellmass_centrifuge
        if crushing_mill is not None: 
            self.crushing_mill_oil_recovery = crushing_mill.isplit['Oil'].mean()
        else:
            self.crushing_mill_oil_recovery = None
        if pressure_filter is not None: 
            self.bagasse_oil_recovery = 1 - pressure_filter.isplit['Oil'].mean()
            if pressure_filter in cellmass_centrifuge.get_downstream_units(facilities=False):
                # Microbial oil is removed further downstream
                @pressure_filter.add_specification(run=True)
                def adjust_oil_recovery():
                    cellmass = self.cellmass_centrifuge.outs[0]
                    pressure_filter = self.pressure_filter
                    bagasse_oil_recovery = self.bagasse_oil_recovery
                    microbial_oil_recovery = self.microbial_oil_recovery
                    f = cellmass.imass['Oil'] / pressure_filter.ins[0].imass['Oil']
                    pressure_filter.isplit['Oil'] = (1.
                        - f * microbial_oil_recovery
                        - (1. - f) * bagasse_oil_recovery
                    )
            else:
                @pressure_filter.add_specification(run=True)
                def adjust_oil_recovery():
                    self.pressure_filter.isplit['Oil'] = 1. - self.bagasse_oil_recovery
        else:
            self.bagasse_oil_recovery = None
        self.microbial_oil_recovery = microbial_oil_recovery
      
    def load_crushing_mill_oil_recovery(self, recovery):
        if self.crushing_mill is None: return
        self.crushing_mill_oil_recovery = recovery
        self.crushing_mill.isplit['Oil'] = 1. - recovery
    
    def load_bagasse_oil_recovery(self, recovery):
        if self.pressure_filter is None: return
        self.bagasse_oil_recovery = recovery 
    
    def load_microbial_oil_recovery(self, recovery):
        if self.pressure_filter is None: return
        self.microbial_oil_recovery = recovery 
    
    def load_specifications(self, 
            crushing_mill_oil_recovery=None,
            bagasse_oil_recovery=None,
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
        if bagasse_oil_recovery is None:
            bagasse_oil_recovery = self.bagasse_oil_recovery
        self.load_crushing_mill_oil_recovery(crushing_mill_oil_recovery)
        self.load_bagasse_oil_recovery(bagasse_oil_recovery)
    
