# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

__all__ = ('dextrose_equivalent',)

def dextrose_equivalent(n):
    """
    Return the dextrose equivalent of starch given the degree of polymerization.

    Parameters
    ----------
    n : int
        Degree of polymerization.

    Returns
    -------
    DE : float
        Dextrose equivalent.

    Notes
    -----
    The dextrose equivalent (DE) is a measure of the amount of reducing sugars 
    present in a sugar product, expressed as a percentage on a dry basis 
    relative to dextrose. For polymerized glucose (i.e. dissacharides, 
    oligosaccharides, and polysaccharides), the dextrose equivalent is given 
    by the following formula:
        
    .. math::
        DE = 100% \frac{180}{180n - 18(n - 1)}
        
    This formula (and this function) is not valid for sugars with linkages other 
    than 1-4 and 1-6 glycosidic linkages.
    
    """
    if n < 1: 
        raise ValueError('degree of polymerization, n, must be greater or equal to 1')
    n = float(n)
    MW_glucose = 180.
    MW_water = 18.
    return 100. * MW_glucose / (MW_glucose * n - MW_water * (n - 1))