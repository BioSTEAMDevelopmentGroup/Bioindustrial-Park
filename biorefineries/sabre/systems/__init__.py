# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autofunction:: biorefineries.sabre.systems.create_biomethane_system
.. autofunction:: biorefineries.sabre.systems.create_vfa_system
.. autofunction:: biorefineries.sabre.systems.create_microbial_oil_system
.. autofunction:: biorefineries.sabre.systems.create_biostimulant_system
.. autofunction:: biorefineries.sabre.systems.create_ad_integrated_system
.. autofunction:: biorefineries.sabre.systems.create_3hp_system

Note: create_3hp_system() rebuilds the process-global chemical set with
sabre._chemicals.create_chemicals(include_hp3=True) every time it's called
-- see _3hp_system.py's own module docstring. Don't call it in the same
Python process as any of the other systems above; not wired into
sabre.__init__.load()'s dispatcher for this reason.
"""
from . import _biomethane_system
from . import _vfa_system
from . import _microbial_oil_system
from . import _biostimulant_system
from . import _ad_integrated_system
from . import _3hp_system

from ._biomethane_system import *
from ._vfa_system import *
from ._microbial_oil_system import *
from ._biostimulant_system import *
from ._ad_integrated_system import *
from ._3hp_system import *

__all__ = (
    *_biomethane_system.__all__,
    *_vfa_system.__all__,
    *_microbial_oil_system.__all__,
    *_biostimulant_system.__all__,
    *_ad_integrated_system.__all__,
    *_3hp_system.__all__,
)
