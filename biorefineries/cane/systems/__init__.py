# -*- coding: utf-8 -*-
"""

.. contents:: :local:

Oilcane
-------
.. automodule:: biorefineries.cane.systems.oilcane

Sugarcane
---------
.. automodule:: biorefineries.cane.systems.sugarcane

Biodiesel
---------
.. autofunction:: biorefineries.cane.systems.create_oilcane_to_biodiesel_1g
.. autofunction:: biorefineries.cane.systems.create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation
.. autofunction:: biorefineries.cane.systems.create_acTAG_separation_system
.. autofunction:: biorefineries.cane.systems.create_oilcane_to_biodiesel_and_actag_combined_1_and_2g_post_fermentation_oil_separation

Fermentation
------------
.. automodule:: biorefineries.cane.systems.fermentation.add_urea_MgSO4_nutrients
.. automodule:: biorefineries.cane.systems.fermentation.add_urea_nutrient
.. automodule:: biorefineries.cane.systems.fermentation.create_sucrose_fermentation_system

Juicing
-------
.. autofunction:: biorefineries.cane.systems.juicing.create_feedstock_handling_system
.. autofunction:: biorefineries.cane.systems.juicing.create_juicing_system_without_treatment
.. autofunction:: biorefineries.cane.systems.juicing.create_juicing_system_up_to_clarification
.. autofunction:: biorefineries.cane.systems.juicing.create_juicing_system

Bagasse
-------
.. autofunction:: biorefineries.cane.systems.bagasse.create_bagasse_pelleting_system
.. autofunction:: biorefineries.cane.systems.bagasse.create_bagasse_drying_system

Lipid extraction
----------------
.. autofunction:: biorefineries.cane.systems.lipid_extraction.create_lipid_exctraction_system
.. autofunction:: biorefineries.cane.systems.lipid_extraction.create_post_fermentation_oil_separation_system

Pretreatment
------------
.. autofunction:: biorefineries.cane.systems.pretreatment.create_cane_combined_1_and_2g_pretreatment

Sugar
-----
.. autofunction:: biorefineries.cane.systems.sugar.create_sugar_crystallization_system

"""
from . import bagasse
from . import fermentation
from . import juicing 
from . import sugar
from . import oilcane 
from . import sugarcane 
from . import lipidcane
from . import biodiesel
from . import cellulosic_biodiesel
from . import biodiesel_actag
from . import cellulosic_biodiesel_actag

__all__ = (
    *bagasse.__all__,
    *fermentation.__all__,
    *juicing.__all__,
    *sugar.__all__,
    *oilcane.__all__,
    *sugarcane.__all__, 
    *lipidcane.__all__,
    *biodiesel.__all__,
    *cellulosic_biodiesel.__all__,
    *biodiesel_actag.__all__,
    *cellulosic_biodiesel_actag.__all__,
)

from .bagasse import *
from .fermentation import *
from .juicing import *
from .sugar import *
from .oilcane import *
from .sugarcane import *
from .lipidcane import *
from .biodiesel import *
from .cellulosic_biodiesel import *
from .biodiesel_actag import *
from .cellulosic_biodiesel_actag import *