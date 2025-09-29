# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 13:14:15 2025

@author: saran
"""

from . import (
    # system, 
    units,
    )

__all__ = [
            # *system.__all__,
           *units.__all__,
]

from .units import *
# from .system import *

# def load(*args, **kwargs):
#     br = Biorefinery(*args, **kwargs)
#     globals().update(br.__dict__)
#     globals().update({
#         'biorefinery': br,
#         'system': br.system,
#         'tea': br.TEA,
#     })
