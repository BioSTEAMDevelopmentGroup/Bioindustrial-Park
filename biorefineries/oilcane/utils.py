# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from biorefineries.cane.composition import (
    set_lipid_fraction as set_oil_fraction,
    get_composition,
)
from biorefineries.oilcane._load_data import load_composition as load_feedstock_composition
from biorefineries import oilcane as oc
import chaospy as cp

__all__ = (
    'OilExtractionSpecification', 
    'MockExtractionSpecification',
    'CaneCompositionSpecification',
    'set_composition_parameters',
    'set_line_composition_parameters',
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
    
    
class CaneCompositionSpecification:
    """
    Create a CaneCompositionSpecification object for setting the 
    composition of cane feedstocks.
    
    Parameters
    ----------
    feedstock : Stream
        Cane feedstock.
    oil : float 
        Oil content of feedstocks [dry wt. %].
    moisture : float 
        Moisture content of feedstocks [wt. %].
    fiber : float 
        Fiber content of feedstocks [dry wt. %].
    FFA : float 
        FFA content of feedstocks [oil wt. %].
    PL : float 
        PL content of feedstocks [oil wt. %].
    
    """
    __slots__ = (
        'feedstock',
        'oil',
        'FFA',
        'PL',
        'moisture',
        'fiber',
    )
    
    def __init__(self, feedstock, oil=None, moisture=None, fiber=None, FFA=None, PL=None):
        composition = get_composition(feedstock) # Get default values
        self.feedstock = feedstock #: [Stream] Cane feedstock.
        self.oil = composition['lipid'] if oil is None else oil #: [float] Oil content of feedstocks [dry wt. %].
        self.moisture = composition['moisture'] if moisture is None else moisture
        self.fiber = composition['fiber'] if fiber is None else fiber
        self.FFA = composition['FFA'] if FFA is None else FFA
        self.PL = composition['PL'] if PL is None else PL
        
    def load_oil_content(self, oil=None, **kwargs):
        """
        Set the oil content [dry wt. %] of cane feedstocks assuming that an 
        increase in oil content is accompanied by a decrease in sugar content 
        using an energy balance and an increase in fiber using a mass balance 
        (oil is more energy dense than sugar).
        
        """
        if oil is None: oil = self.oil
        set_oil_fraction(
            oil, self.feedstock, PL_fraction=self.PL, 
            FFA_fraction=self.FFA, **kwargs,
        )
        self.oil = oil
        
    def load_composition(self, oil=None, moisture=None, fiber=None, FFA=None, PL=None):
        if oil is None: oil = self.oil
        if moisture is None: moisture = self.moisture
        if fiber is None: fiber = self.fiber
        if FFA is None: FFA = self.FFA
        if PL is None: PL = self.PL
        load_feedstock_composition(self.feedstock, oil, moisture, fiber, FFA, PL)
        try:
            assert (self.feedstock.mol >= 0.).all()
        except:
            breakpoint()
    
    
def bounded_gaussian_distribution_from_mean_and_std(mean, std):
    if mean == 0: mean = 1e-6
    if std == 0: std = mean / 1e3
    normal = cp.Normal(mean, std)
    trunc = cp.Trunc(normal, lower=mean - 2 * std, upper=mean + 2 * std)
    return trunc

def set_line_composition_parameters(line):
    df = oc.get_composition_data()
    line = df.loc[line]
    oil = line['Stem Oil (dw)']
    water = line['Water (wt)']
    fiber = line['Fiber (dw)']
    biomass = line['Dry biomass yield (MT/hc)']
    set_composition_parameters(
        oil['Mean'],
        oil['STD'],
        fiber['Mean'],
        fiber['STD'],
        water['Mean'],
        water['STD'],
        biomass['Mean'],
        biomass['STD'],
    )
    
def set_composition_parameters(
        mean_oil_content,
        std_oil_content,
        mean_fiber_content,
        std_fiber_content,
        mean_moisture_content,
        std_moisture_content,
        mean_biomass_yield,
        std_biomass_yield,
    ):
    try:
        _reset_composition_parameters_for_cane_line(
            mean_oil_content,
            std_oil_content,
            mean_fiber_content,
            std_fiber_content,
            mean_moisture_content,
            std_moisture_content,
            mean_biomass_yield,
            std_biomass_yield,
        )
    except:
        _add_model_composition_parameters_for_cane_line(
            mean_oil_content,
            std_oil_content,
            mean_fiber_content,
            std_fiber_content,
            mean_moisture_content,
            std_moisture_content,
            mean_biomass_yield,
            std_biomass_yield,
        )

def _add_model_composition_parameters_for_cane_line(
        mean_oil_content,
        std_oil_content,
        mean_fiber_content,
        std_fiber_content,
        mean_moisture_content,
        std_moisture_content,
        mean_biomass_yield,
        std_biomass_yield,
        model=None,
    ):
    if model is None: model = oc.model
    removed = {oc.set_cane_oil_content, oc.set_relative_sorghum_oil_content}
    model.parameters = tuple([i for i in model.parameters if i not in removed])
    parameter = model.parameter
    @parameter(
        element=oc.feedstock, 
        baseline=mean_oil_content,
        distribution=bounded_gaussian_distribution_from_mean_and_std(
            mean_oil_content, std_oil_content
        ),
        units='dw %',
    )
    def set_cane_oil_content(oil_content):
        oc.composition_specification.oil = oil_content
    
    @parameter(
        element=oc.feedstock, 
        baseline=mean_fiber_content,
        distribution=bounded_gaussian_distribution_from_mean_and_std(
            mean_fiber_content, std_fiber_content
        ),
        units='dw %',
    )
    def set_cane_fiber_content(fiber_content):
        oc.composition_specification.fiber = fiber_content
    
    @parameter(
        element=oc.feedstock, 
        baseline=mean_moisture_content,
        distribution=bounded_gaussian_distribution_from_mean_and_std(
            mean_moisture_content, std_moisture_content
        ),
        units='wt %',
    )
    def set_cane_moisture_content(moisture_content):
        oc.composition_specification.moisture = moisture_content
        oc.composition_specification.load_composition()
        
    @parameter(
        element=oc.feedstock, 
        baseline=mean_biomass_yield,
        distribution=bounded_gaussian_distribution_from_mean_and_std(
            mean_biomass_yield, std_biomass_yield
        ),
        units='dry MT/hc',
    )
    def set_cane_biomass_yield(biomass_yield):
        oc.dry_biomass_yield = biomass_yield
        
    oc.set_cane_oil_content = set_cane_oil_content
    oc.set_cane_fiber_content = set_cane_fiber_content
    oc.set_cane_moisture_content = set_cane_moisture_content
    oc.set_cane_biomass_yield = set_cane_biomass_yield
    
def _reset_composition_parameters_for_cane_line(
        mean_oil_content,
        std_oil_content,
        mean_fiber_content,
        std_fiber_content,
        mean_moisture_content,
        std_moisture_content,
        mean_biomass_yield,
        std_biomass_yield,
        parameters=None, 
    ):
    if parameters is None: 
        parameters = (
            oc.set_cane_oil_content, 
            oc.set_cane_fiber_content, 
            oc.set_cane_moisture_content,
            oc.set_cane_biomass_yield,
        )
    oil, fiber, moisture, biomass = parameters
    oil.baseline = mean_oil_content
    oil.distribution = bounded_gaussian_distribution_from_mean_and_std(
        mean_oil_content, std_oil_content
    )
    fiber.baseline = mean_fiber_content
    fiber.distribution = bounded_gaussian_distribution_from_mean_and_std(
        mean_fiber_content, std_fiber_content
    )
    moisture.baseline = mean_moisture_content
    moisture.distribution = bounded_gaussian_distribution_from_mean_and_std(
        mean_moisture_content, std_moisture_content
    )
    biomass.baseline = mean_biomass_yield
    biomass.distribution = bounded_gaussian_distribution_from_mean_and_std(
        mean_biomass_yield, std_biomass_yield
    )
    