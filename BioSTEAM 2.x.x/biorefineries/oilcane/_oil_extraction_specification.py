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

def evaluate_across_oil_content(spec, efficiency, oil_retention, 
                                  metrics, oil_content):
    spec.load_specifications(efficiency=efficiency, oil_retention=oil_retention)
    return spec._evaluate_across_oil_content(metrics, oil_content)
    
evaluate_across_oil_content = np.vectorize(
    evaluate_across_oil_content, 
    excluded=['spec', 'metrics', 'oil_content'],
    signature='(),(),(),(m),(p)->(m,p)'
)

def evaluate_across_oil_retention(spec, efficiency, oil_content, 
                                    metrics, oil_retention):
    spec.load_specifications(efficiency=efficiency, oil_content=oil_content)
    return spec._evaluate_across_oil_retention(metrics, oil_retention)
    
evaluate_across_oil_retention = np.vectorize(
    evaluate_across_oil_retention, 
    excluded=['spec', 'metrics', 'oil_retention'],
    signature='(),(),(),(m),(p)->(m,p)'
)       

class MockExtractionSpecification:
    
    def load_oil_content(self, oil_content):
        pass
      
    def load_efficiency(self, efficiency):
        pass
    
    def load_oil_retention(self, oil_retention):
        pass
    
    def load_specifications(self, 
            efficiency=None,
            oil_retention=None,
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
    isplit_oil_retention : ChemicalIndexer
        Defines bagasse oil retention as a material split.
    oil_content : float 
        Oil content of feedstocks [dry wt. %].
        
    Notes
    -----
    Use the `load_specifications` method to load specifications.
    Setting attributes (i.e. yield_, titer, productivity, production) sets 
    defaults for this method, but does not actually load any specifications.
    
    """
    
    __slots__ = (
        'system',
        'feedstocks',
        'isplit_efficiency',
        'isplit_oil_retention',
        'efficiency',
        'oil_retention',
        'oil_content',
        'locked_oil_content',
        'isplit_efficiency_is_reversed',
        'FFA_content',
        'PL_content',
    )
    
    def __init__(self, system, feedstocks, isplit_efficiency, isplit_oil_retention, 
                 isplit_efficiency_is_reversed=False, efficiency=0.9, 
                 oil_retention=0.9, oil_content=0.1, FFA_content=0.1,
                 PL_content=0.1):
        self.system = system #: [System] System associated to feedstocks
        self.feedstocks = feedstocks #: Stream Oil feedstocks
        self.isplit_efficiency = isplit_efficiency #: [ChemicalIndexer] Defines extraction efficiency as a material split.
        self.isplit_oil_retention = isplit_oil_retention #: [ChemicalIndexer] Defines bagasse oil retention as a material split.
        self.isplit_efficiency_is_reversed = isplit_efficiency_is_reversed #: [bool] Whether oil extraction efficiency is 1 - split.
        self.efficiency = efficiency #: [float] Oil extraction efficiency b
        self.oil_retention = oil_retention #: [float] Oil extraction oil retention
        self.oil_content = oil_content #: [float] Oil content of feedstocks [dry wt. %].
        self.locked_oil_content = False
        self.FFA_content = FFA_content
        self.PL_content = PL_content
        
    def MFPP(self):
        return self.system.TEA.solve_price(self.feedstocks)
        
    def dMFPP_over_doil_content_at_efficiency(self, efficiency, doil_content=0.01):
        oil_content = self.oil_content
        system = self.system
        self.load_efficiency(efficiency)   
        self.load_oil_content(oil_content + doil_content)
        system.simulate()
        MFPP1 = self.MFPP()
        self.load_oil_content(oil_content)
        system.simulate()
        MFPP0 = self.MFPP()
        return (MFPP1 - MFPP0) / doil_content
        
    def solve_MFPP_inflection(self, oil_retention=None, oil_content=None):
        if oil_content is not None: self.load_oil_content(oil_content)
        if oil_retention is not None: self.load_oil_retention(oil_retention)
        f = self.dMFPP_over_doil_content_at_efficiency
        x0 = 0.1
        x1 = 1.
        y0 = f(x0)
        y1 = f(x1)
        if not y0 < 0. < y1: return np.nan
        return flx.IQ_interpolation(f, 0.1, 1., y0, y1, xtol=1e-4, ytol=1e-4, checkiter=False,
                                    maxiter=10)
    
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
        if self.isplit_efficiency_is_reversed:
            self.isplit_efficiency['Oil'] = 1. - efficiency
        else:
            self.isplit_efficiency['Oil'] = efficiency
    
    def load_oil_retention(self, oil_retention):
        if self.isplit_efficiency is None: return
        self.oil_retention = oil_retention
        self.isplit_oil_retention['Oil'] = oil_retention    
    
    def load_specifications(self, 
            efficiency=None,
            oil_retention=None,
            oil_content=None,
        ):
        """
        Load oil extraction specifications.

        Parameters
        ----------
        efficiency : float, optional
            Oil extraction efficiency a.
        oil_retention : float, optional
            Oil extraction efficiency b.
        oil_content : float, optional
            Oil content of feedstocks [dry wt. %].

        """
        if efficiency is None: efficiency = self.efficiency
        if oil_retention is None: oil_retention = self.oil_retention
        if oil_content is None: oil_content = self.oil_content
        self.load_oil_retention(oil_retention)
        self.load_efficiency(efficiency)
        self.load_oil_content(oil_content)

    def _evaluate_across_oil_content(self, metrics, oil_content):
        """
        Evaluate metrics across oil content and return an array with the all
        metric results.
        
        Parameters
        ----------
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        oil_content : array_like[P elements]
            Oil contents to evaluate.
        
        Returns
        -------
        results : array[M x P]
            All metric results.
        
        """
        M = len(metrics)
        P = len(oil_content)
        data = np.zeros([M, P])
        for i in range(P):
            self.load_oil_content(oil_content[i])
            self.system.simulate()
            data[:, i] = [j() for j in metrics]
        return data
    
    def _evaluate_across_oil_retention(self, metrics, oil_retention):
        """
        Evaluate metrics across oil retention and return an array with the all
        metric results.
        
        Parameters
        ----------
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        oil_retention: array_like[P elements]
            Oil retentions to evaluate.
        
        Returns
        -------
        results : array[M x P]
            All metric results.
        
        """
        M = len(metrics)
        P = len(oil_retention)
        data = np.zeros([M, P])
        for i in range(P):
            self.load_oil_retention(oil_retention[i])
            self.system.simulate()
            data[:, i] = [j() for j in metrics]
        return data

    def evaluate_across_oil_content(self, 
            efficiency, oil_retention, 
            metrics, oil_content
        ):
        """
        Evaluate metrics at given efficiency and oil retention across a set of 
        oil contents. Return an array with the all metric results.
            
        Parameters
        ----------
        efficiency : float
            Oil extraction efficiency.
        oil_retention : float
            Fraction of oil retained in bagasse.
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        oil_content : array_like[P elements]
            Oil content of feedstocks [dry wt. %].
        
        Returns
        -------
        results : array[shape x M x P]
            All metric results at given specifications.
        
        Notes
        -----
        This method is vectorized along titer and yield. If, for example,
        the parameters had the following dimensions:
            
        efficiency [Y x T], oil_retention [Y x T], metrics [M], oil_content [P]
        
        This method would return an array with the following dimensions:
        
        results [Y x T x M x P]
        
        """
        return evaluate_across_oil_content(self, 
                                             efficiency, oil_retention, 
                                             metrics, oil_content)
    def evaluate_across_oil_retention(self, 
            efficiency, oil_content, 
            metrics, oil_retention
        ):
        """
        Evaluate metrics at given efficiency and oil retention across a set of 
        oil contents. Return an array with the all metric results.
            
        Parameters
        ----------
        efficiency : float
            Oil extraction efficiency.
        oil_content : float
            Oil content of feedstocks [dry wt. %].
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        oil_retention : array_like[P elements]
            Fraction of oil retained in bagasse.
        
        Returns
        -------
        results : array[shape x M x P]
            All metric results at given specifications.
        
        Notes
        -----
        This method is vectorized along efficiency and oil_content. If, for example,
        the parameters had the following dimensions:
            
        efficiency [Y x T], oil_content [Y x T], metrics [M], oil_retention [P]
        
        This method would return an array with the following dimensions:
        
        results [Y x T x M x P]
        
        """
        return evaluate_across_oil_retention(self, 
                                               efficiency, oil_content, 
                                               metrics, oil_retention)
