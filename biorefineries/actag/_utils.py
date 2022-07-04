# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 17:16:05 2021

@author: Latin Fire
"""

__all__ = ('product_viscosity', )

from math import log, exp

def product_viscosity(acTAG, TAG, nu_acTAG=20.3, nu_TAG=30.6, rule='log'):
    total = acTAG + TAG
    acTAG /= total
    TAG /= total
    if rule == 'log':
        return exp(acTAG*log(nu_acTAG) + TAG*log(nu_TAG))
    elif rule == 'linear':
        return acTAG*nu_acTAG + TAG*nu_TAG
    else:
        raise ValueError("rule must be either log or linear")