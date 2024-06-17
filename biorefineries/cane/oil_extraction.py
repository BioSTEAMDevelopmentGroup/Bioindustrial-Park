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
      
    def load_juicing_oil_recovery(self, recovery):
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
        'juice',
        'screw_press',
        'crushing_mill_oil_recovery',
        'bagasse_oil_recovery',
        'microbial_oil_recovery',
        'juice_oil_recovery',
        'microbial_oil_centrifuge',
        'bioreactor',
    )
    
    def __init__(self, system, crushing_mill, pressure_filter=None, cellmass_centrifuge=None, 
                 microbial_oil_centrifuge=None, juice=None, screw_press=None, 
                 bioreactor=None, microbial_oil_recovery=None, juice_oil_recovery=None):
        self.crushing_mill = crushing_mill #: [Unit] Defines extraction recovery as a material split.
        self.pressure_filter = pressure_filter #: [Unit] Defines extraction recovery from bagasse as a material split.
        self.screw_press = screw_press
        self.cellmass_centrifuge = cellmass_centrifuge
        self.microbial_oil_centrifuge = microbial_oil_centrifuge
        self.bioreactor = bioreactor
        self.juice = juice
        if crushing_mill is None: 
            self.crushing_mill_oil_recovery = None
        else:
            self.crushing_mill_oil_recovery = crushing_mill.isplit['Oil'].mean()
        if juice is None:
            self.juice_oil_recovery = juice_oil_recovery
        else:
            self.juice_oil_recovery = 1. if juice_oil_recovery is None else juice_oil_recovery
        if pressure_filter is None:
            self.bagasse_oil_recovery = None
        else:
            self.bagasse_oil_recovery = 1. - pressure_filter.isplit['Oil'].mean()
        if microbial_oil_centrifuge:
            if juice is None: raise ValueError('must pass both juice stream and pressure filter unit')
            self.microbial_oil_recovery = 0.7 if microbial_oil_recovery is None else microbial_oil_recovery
            # Microbial oil is removed further downstream
            @pressure_filter.add_specification(run=True)
            def adjust_bagasse_oil_recovery():
                pressure_filter = self.pressure_filter
                feedstock = self.crushing_mill.ins[0]
                juice = self.juice
                bagasse_oil = feedstock.imass['Oil'] - juice.imass['Oil']
                total_oil = pressure_filter.ins[0].imass['Oil']
                bagasse_oil_fraction = bagasse_oil / total_oil
                pressure_filter.isplit['Oil'] = (1.
                    - bagasse_oil_fraction * self.bagasse_oil_recovery
                )
                
            @cellmass_centrifuge.add_specification(run=True)
            def adjust_microbial_oil_split():
                bioreactor = self.bioreactor
                cellmass_oil = bioreactor.outs[1].imass['Oil'] - sum([i.imass['Oil'] for i in bioreactor.ins])
                cellmass_centrifuge = self.cellmass_centrifuge
                total_oil = cellmass_centrifuge.ins[0].imass['Oil']
                cellmass_centrifuge.aqueous_isplit['Oil'] = cellmass_oil / total_oil
                cellmass_centrifuge.solids_isplit['Oil', 'Cellmass'] = 1.
                
            @microbial_oil_centrifuge.add_specification(run=True)
            def adjust_microbial_oil_recovery():
                cellmass_centrifuge = self.cellmass_centrifuge
                cellmass_oil = cellmass_centrifuge.outs[2].imass['Oil']
                total_oil = microbial_oil_centrifuge.ins[0].imass['Oil']
                cellmass_oil_fraction = cellmass_oil / total_oil
                microbial_oil_centrifuge.isplit['Oil'] = (1.
                    - cellmass_oil_fraction * self.microbial_oil_recovery
                )
        elif screw_press:
            # Microbial oil is removed separately from bagasse oil
            self.microbial_oil_recovery = 1. - screw_press.isplit['Oil'].mean()
            
            @cellmass_centrifuge.add_specification(run=True)
            def adjust_oil_recovery():
                cellmass_centrifuge = self.cellmass_centrifuge
                total_oil = cellmass_centrifuge.ins[0].imass['Oil']
                cellmass_oil = bioreactor.outs[1].imass['Oil'] - sum([i.imass['Oil'] for i in bioreactor.ins])
                cellmass_centrifuge.aqueous_split[:] = 1
                cellmass_centrifuge.aqueous_isplit['Oil', 'Water'] = [cellmass_oil / total_oil, 0.9999]
                cellmass_centrifuge.solids_isplit['Oil', 'Cellmass'] = 1.
        else:
            # No microbial oil
            self.microbial_oil_recovery = None
      
    def load_juicing_oil_recovery(self, recovery):
        if self.crushing_mill is None: return
        self.crushing_mill_oil_recovery = recovery
        self.crushing_mill.isplit['Oil'] = 1. - recovery
    
    def load_bagasse_oil_recovery(self, recovery):
        if self.pressure_filter is None: return
        self.bagasse_oil_recovery = recovery 
        self.pressure_filter.isplit['Oil'] = 1. - self.bagasse_oil_recovery
    
    def load_microbial_oil_recovery(self, recovery):
        if self.cellmass_centrifuge is None: return
        self.microbial_oil_recovery = recovery 
        if self.screw_press: self.screw_press.isplit['Oil'] = 1. - recovery
    
    def load_specifications(self, 
            crushing_mill_oil_recovery=None,
            bagasse_oil_recovery=None,
            microbial_oil_recovery=None,
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
        if microbial_oil_recovery is None:
            microbial_oil_recovery = self.microbial_oil_recovery
        self.load_juicing_oil_recovery(crushing_mill_oil_recovery)
        self.load_bagasse_oil_recovery(bagasse_oil_recovery)
        self.load_microbial_oil_recovery(microbial_oil_recovery)
    
