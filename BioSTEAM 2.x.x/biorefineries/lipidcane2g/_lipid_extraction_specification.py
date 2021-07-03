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
from biorefineries.lipidcane import set_lipid_fraction
import flexsolve as flx
import thermosteam as tmo

__all__ = ('LipidExtractionSpecification', 'MockExtractionSpecification')

def evaluate_across_lipid_content(spec, efficiency, lipid_retention, 
                                  metrics, lipid_content):
    spec.load_specifications(efficiency=efficiency, lipid_retention=lipid_retention)
    return spec._evaluate_across_lipid_content(metrics, lipid_content)
    
evaluate_across_lipid_content = np.vectorize(
    evaluate_across_lipid_content, 
    excluded=['spec', 'metrics', 'lipid_content'],
    signature='(),(),(),(m),(p)->(m,p)'
)

def evaluate_across_lipid_retention(spec, efficiency, lipid_content, 
                                    metrics, lipid_retention):
    spec.load_specifications(efficiency=efficiency, lipid_content=lipid_content)
    return spec._evaluate_across_lipid_retention(metrics, lipid_retention)
    
evaluate_across_lipid_retention = np.vectorize(
    evaluate_across_lipid_retention, 
    excluded=['spec', 'metrics', 'lipid_retention'],
    signature='(),(),(),(m),(p)->(m,p)'
)       

class MockExtractionSpecification:
    
    def load_lipid_content(self, lipid_content):
        pass
      
    def load_efficiency(self, efficiency):
        pass
    
    def load_lipid_retention(self, lipid_retention):
        pass
    
    def load_specifications(self, 
            efficiency=None,
            lipid_retention=None,
            lipid_content=None,
        ):
        pass

class LipidExtractionSpecification:
    """
    Create a LipidExtractionSpecification object for setting process 
    specifications related to lipid extraction.
    
    Parameters
    ----------
    system : System
        System associated to feedstock.
    feedstock : Stream
        Lipidcane feedstock.
    isplit_efficiency : ChemicalIndexer
        Defines extraction efficiency as a material split.
    isplit_lipid_retention : ChemicalIndexer
        Defines bagasse lipid retention as a material split.
    lipid_content : float 
        Lipid content of feedstock [dry wt. %].
        
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
        'isplit_lipid_retention',
        'efficiency',
        'lipid_retention',
        'lipid_content',
        'locked_lipid_content',
        'isplit_efficiency_is_reversed',
    )
    
    def __init__(self, system, feedstocks, isplit_efficiency, isplit_lipid_retention, 
                 isplit_efficiency_is_reversed=False, efficiency=0.9, 
                 lipid_retention=0.9, lipid_content=0.1):
        self.system = system #: [System] System associated to feedstock
        self.feedstocks = feedstocks #: tuple[Stream] Lipid feedstocks
        self.isplit_efficiency = isplit_efficiency #: [ChemicalIndexer] Defines extraction efficiency as a material split.
        self.isplit_lipid_retention = isplit_lipid_retention #: [ChemicalIndexer] Defines bagasse lipid retention as a material split.
        self.isplit_efficiency_is_reversed = isplit_efficiency_is_reversed #: [bool] Whether lipid extraction efficiency is 1 - split.
        self.efficiency = efficiency #: [float] Lipid extraction efficiency b
        self.lipid_retention = lipid_retention #: [float] Lipid extraction lipid retention
        self.lipid_content = lipid_content #: [float] Lipid content of feedstock [dry wt. %].
        self.locked_lipid_content = False
        
    def MFPP(self):
        return self.system.TEA.solve_price(self.feedstocks)
        
    def dMFPP_over_dlipid_content_at_efficiency(self, efficiency, dlipid_content=0.01):
        lipid_content = self.lipid_content
        system = self.system
        self.load_efficiency(efficiency)   
        self.load_lipid_content(lipid_content + dlipid_content)
        system.simulate()
        MFPP1 = self.MFPP()
        self.load_lipid_content(lipid_content)
        system.simulate()
        MFPP0 = self.MFPP()
        return (MFPP1 - MFPP0) / dlipid_content
        
    def solve_MFPP_inflection(self, lipid_retention=None, lipid_content=None):
        if lipid_content is not None: self.load_lipid_content(lipid_content)
        if lipid_retention is not None: self.load_lipid_retention(lipid_retention)
        f = self.dMFPP_over_dlipid_content_at_efficiency
        x0 = 0.1
        x1 = 1.
        y0 = f(x0)
        y1 = f(x1)
        if not y0 < 0. < y1: return np.nan
        return flx.IQ_interpolation(f, 0.1, 1., y0, y1, xtol=1e-4, ytol=1e-4, checkiter=False,
                                    maxiter=10)
    
    def load_lipid_content(self, lipid_content):
        if self.locked_lipid_content: return
        for i in self.feedstocks: set_lipid_fraction(lipid_content, i)
        self.lipid_content = lipid_content
      
    def load_efficiency(self, efficiency):
        if self.isplit_efficiency is None: return
        self.efficiency = efficiency
        if self.isplit_efficiency_is_reversed:
            self.isplit_efficiency['Lipid'] = 1. - efficiency
        else:
            self.isplit_efficiency['Lipid'] = efficiency
    
    def load_lipid_retention(self, lipid_retention):
        if self.isplit_efficiency is None: return
        self.lipid_retention = lipid_retention
        self.isplit_lipid_retention['Lipid'] = lipid_retention    
    
    def load_specifications(self, 
            efficiency=None,
            lipid_retention=None,
            lipid_content=None,
        ):
        """
        Load lipid extraction specifications.

        Parameters
        ----------
        efficiency : float, optional
            Lipid extraction efficiency a.
        lipid_retention : float, optional
            Lipid extraction efficiency b.
        lipid_content : float, optional
            Lipid content of feedstock [dry wt. %].

        """
        if efficiency is None: efficiency = self.efficiency
        if lipid_retention is None: lipid_retention = self.lipid_retention
        if lipid_content is None: lipid_content = self.lipid_content
        self.load_lipid_retention(lipid_retention)
        self.load_efficiency(efficiency)
        self.load_lipid_content(lipid_content)

    def _evaluate_across_lipid_content(self, metrics, lipid_content):
        """
        Evaluate metrics across lipid content and return an array with the all
        metric results.
        
        Parameters
        ----------
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        lipid_content : array_like[P elements]
            Lipid contents to evaluate.
        
        Returns
        -------
        results : array[M x P]
            All metric results.
        
        """
        M = len(metrics)
        P = len(lipid_content)
        data = np.zeros([M, P])
        for i in range(P):
            self.load_lipid_content(lipid_content[i])
            self.system.simulate()
            data[:, i] = [j() for j in metrics]
        return data
    
    def _evaluate_across_lipid_retention(self, metrics, lipid_retention):
        """
        Evaluate metrics across lipid retention and return an array with the all
        metric results.
        
        Parameters
        ----------
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        lipid_retention: array_like[P elements]
            Lipid retentions to evaluate.
        
        Returns
        -------
        results : array[M x P]
            All metric results.
        
        """
        M = len(metrics)
        P = len(lipid_retention)
        data = np.zeros([M, P])
        for i in range(P):
            self.load_lipid_retention(lipid_retention[i])
            self.system.simulate()
            data[:, i] = [j() for j in metrics]
        return data

    def evaluate_across_lipid_content(self, 
            efficiency, lipid_retention, 
            metrics, lipid_content
        ):
        """
        Evaluate metrics at given efficiency and lipid retention across a set of 
        lipid contents. Return an array with the all metric results.
            
        Parameters
        ----------
        efficiency : float
            Lipid extraction efficiency.
        lipid_retention : float
            Fraction of lipid retained in bagasse.
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        lipid_content : array_like[P elements]
            Lipid content of feedstock [dry wt. %].
        
        Returns
        -------
        results : array[shape x M x P]
            All metric results at given specifications.
        
        Notes
        -----
        This method is vectorized along titer and yield. If, for example,
        the parameters had the following dimensions:
            
        efficiency [Y x T], lipid_retention [Y x T], metrics [M], lipid_content [P]
        
        This method would return an array with the following dimensions:
        
        results [Y x T x M x P]
        
        """
        return evaluate_across_lipid_content(self, 
                                             efficiency, lipid_retention, 
                                             metrics, lipid_content)
    def evaluate_across_lipid_retention(self, 
            efficiency, lipid_content, 
            metrics, lipid_retention
        ):
        """
        Evaluate metrics at given efficiency and lipid retention across a set of 
        lipid contents. Return an array with the all metric results.
            
        Parameters
        ----------
        efficiency : float
            Lipid extraction efficiency.
        lipid_content : float
            Lipid content of feedstock [dry wt. %].
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        lipid_retention : array_like[P elements]
            Fraction of lipid retained in bagasse.
        
        Returns
        -------
        results : array[shape x M x P]
            All metric results at given specifications.
        
        Notes
        -----
        This method is vectorized along efficiency and lipid_content. If, for example,
        the parameters had the following dimensions:
            
        efficiency [Y x T], lipid_content [Y x T], metrics [M], lipid_retention [P]
        
        This method would return an array with the following dimensions:
        
        results [Y x T x M x P]
        
        """
        return evaluate_across_lipid_retention(self, 
                                               efficiency, lipid_content, 
                                               metrics, lipid_retention)
