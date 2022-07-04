=====================================================================
fattyalcohols: Production of fatty alcohols by fermentation of sugars
=====================================================================

.. figure:: ./images/LAOs_areas.png

This library models the production of fatty alcohols produced by fermentation.
The library is not yet complete as it only includes the fermentation section
without any separation or purification of products. A tentative date for 
completion is November 2020.

Getting Started
---------------

For now, only the fermentation section can be modeled:

.. code-block:: python

    >>> import biosteam as bst
    >>> from biorefineries import fattyalcohols as fa
    >>> chemicals = fa.create_chemicals()
    >>> bst.settings.set_thermo(chemicals)
    >>> fattyalcohol_production_sys = fa.create_fattyalcohol_production_sys()
    >>> fattyalcohol_production_sys.show()
    System: fattyalcohol_production_sys
     path: (M101, T102, P102, T101, T104, T103,
            T106, M102, H101, R101, T105,
            P104, C101, P107)
    


