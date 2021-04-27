# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import flexsolve as flx
import numpy as np
from thermosteam.exceptions import InfeasibleRegion
from biorefineries.lipidcane import set_lipid_fraction
import thermosteam as tmo

__all__ = ('LipidExtractionSpecification',)

def evaluate_across_specifications(spec, efficiency_a, efficiency_b, 
                                   metrics, lipid_content):
    spec.load_specifications(efficiency_a, efficiency_b)
    return spec.evaluate_across_lipid_content(metrics, lipid_content)
    
evaluate_across_specifications = np.vectorize(
    evaluate_across_specifications, 
    excluded=['spec', 'metrics', 'lipid_content'],
    signature='(),(),(),(m),(p)->(m,p)'
)

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
    isplit_a : ChemicalIndexer
        Defines extraction efficiency as a material split.
    isplit_b : ChemicalIndexer
        Defines extraction efficiency as a material split.
    
    Notes
    -----
    Use the `load_specifications` method to load specifications.
    Setting attributes (i.e. yield_, titer, productivity, production) sets 
    defaults for this method, but does not actually load any specifications.
    
    """
    
    __slots__ = (
        'system',
        'feedstock',
        'isplit_a',
        'isplit_b',
        'efficiency_a',
        'efficiency_b',
        'lipid_content',
    )
    
    def __init__(self, system, feedstock, isplit_a, isplit_b, 
                 efficiency_a=0.9, efficiency_b=0.9, lipid_content=0.1):
        self.system = system #: [System] System associated to feedstock
        self.feedstock = feedstock #: [Stream] Lipidcane feedstock
        self.isplit_a = isplit_a #: [ChemicalIndexer] Defines extraction efficiency as a material split.
        self.isplit_b = isplit_b #: [ChemicalIndexer] Defines extraction efficiency as a material split.
        self.efficiency_a = efficiency_a #: [float] Lipid extraction efficiency b
        self.efficiency_b = efficiency_b #: [float] Lipid extraction efficiency b
        self.lipid_content = lipid_content #: [float] Lipid content of feedstock [dry wt. %].
        
    def load_lipid_content(self, lipid_content):
        set_lipid_fraction(lipid_content, self.feedstock)
        self.lipid_content = lipid_content
      
    def load_efficiency_a(self, efficiency_a):
        self.efficiency_a = efficiency_a
        self.isplit_a['Lipid'] = efficiency_a
    
    def load_efficiency_b(self, efficiency_b):
        self.efficiency_b = efficiency_b
        self.isplit_b['Lipid'] = efficiency_b    
    
    def load_specifications(self, 
            efficiency_a=None,
            efficiency_b=None,
            lipid_content=None,
        ):
        """
        Load lipid extraction specifications.

        Parameters
        ----------
        efficiency_a : float, optional
            Lipid extraction efficiency a.
        efficiency_b : float, optional
            Lipid extraction efficiency b.
        lipid_content : float, optional
            Lipid content of feedstock [dry wt. %].

        """
        if efficiency_a is None: efficiency_a = self.efficiency_a
        if efficiency_b is None: efficiency_b = self.efficiency_b
        if lipid_content is None: lipid_content = self.lipid_content
        self.load_efficiency_a(efficiency_a)
        self.load_efficiency_b(efficiency_b)
        self.load_lipid_content(lipid_content)

    def evaluate_across_lipid_content(self, metrics, lipid_content):
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

    def evaluate_across_specifications(self, 
            efficiency_a, efficiency_b, 
            metrics, lipid_content
        ):
        """
        Evaluate metrics at given titer and yield across a set of 
        productivities. Return an array with the all metric results.
            
        Parameters
        ----------
        efficiency_a : float, optional
            Lipid extraction efficiency a.
        efficiency_b : float, optional
            Lipid extraction efficiency b.
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
            
        efficiency_a [Y x T], efficiency_b [Y x T], metrics [M], lipid_content [P]
        
        This method would return an array with the following dimensions:
        
        results [Y x T x M x P]
        
        """
        return evaluate_across_specifications(self, 
                                              efficiency_a, efficiency_b, 
                                              metrics, lipid_content)

