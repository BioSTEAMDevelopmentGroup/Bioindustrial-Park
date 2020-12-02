=================================================================
lipidcane: Co-production of Ethanol and Biodiesel from Lipid-cane
=================================================================

.. figure:: ./images/lipidcane_areas.png

The lipidcane biorefinery design for the co-production of ethanol and biodiesel
follows all assumptions from a study by Huang et. al. executed with SuperPro 
Designer [1]_. As discussed in the original BioSTEAM manuscript [2]_, the 
lipid-cane biorefinery can be divided into six main areas: feedstock handling, 
oil/sugar separation, ethanol production, biodiesel production, boiler and 
turbo-generator system, and utilities. During oil and sugar 
separation (Area 200), the cane juice is extracted and the oil is separated 
from the juice by a settler. The bagasse is burned to produce steam 
and electricity for the plant (Area 500), with excess electricity sold to the 
grid. The sugar solution is fermented and distilled (Area 300). The oil is 
trans-esterified with methanol and sodium methoxide catalyst to produce 
biodiesel and glycerol (Area 400). The glycerol is distilled to 80 wt % and 
sold as crude glycerol. The utility area includes on-site recirculation of 
cooling water and chilled water (Area 600).

Getting Started
---------------

To load the biorefinery, simply import it. All data and variables
are lazy loaded by the module:

.. code-block:: python

    >>> from biorefineries import lipidcane as lc
    >>> # This is optional; it forces the biorefinery to load
    >>> # Otherwise, first time accessing will take a bit to load.
    >>> lc.load()
    >>> lc.chemicals # All chemicals used by the biorefinery.
    CompiledChemicals([Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, CH4, Ash, Cellulose, Hemicellulose, Flocculant, Lignin, Solids, Yeast, CaO, Biodiesel, Methanol, Glycerol, HCl, NaOH, NaOCH3, Lipid])
    >>> lc.lipidcane_sys.show() # The complete biorefinery system
    System: lipidcane_sys
     path: (U101, U102, U103, T201,
            juice_extraction_sys, T202, H201,
            T203, P201, T204, T205, P202,
            juice_separation_sys, T207,
            T207_2, H203, T208, C203, F201,
            T403, P403, R401, C401, R402,
            C402, methanol_recycle_sys, F401,
            P407, T409, H401, P408, P405,
            B401, H403, P411, T401, P401,
            T402, P402, T404, P404, S401,
            S202, ethanol_production_sys,
            T408, U202)
     facilities: (CWP, BT, CT, PWC)
    >>> lc.lipidcane_tea.show() # The TEA object
    ConventionalEthanolTEA: lipidcane_sys
     NPV: 2 USD at 17.8% IRR
    >>> lc.flowsheet # The complete flowsheet
    <Flowsheet: lipidcane>
    >>> lc.R301.show() # Any unit operations and streams can be accessed through the module
    Fermentation: R301
    ins...
    [0] d235  from  HXutility-H301
        phase: 'l', T: 295.15 K, P: 101325 Pa
        flow (kmol/hr): Water    4.11e+03
                        Glucose  10.5
                        Sucrose  62.7
                        H3PO4    0.85
    [1] d260  from  MixTank-T305
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  1.37e+03
                        Yeast  1.03e+04
    outs...
    [0] CO2  to  VentScrubber-D301
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    1.83
                        Ethanol  0.456
                        CO2      245
    [1] d236  to  StorageTank-T301
        phase: 'l', T: 296.06 K, P: 101325 Pa
        flow (kmol/hr): Water    5.42e+03
                        Ethanol  244
                        Glucose  13.6
                        H3PO4    0.85
                        Yeast    1.03e+04


References
----------
.. [1] Huang, H.; Long, S.; Singh, V. Techno-Economic Analysis of Biodiesel and 
    Ethanol Co-Production from Lipid-Producing Sugarcane: Biodiesel and Ethanol
    Co-Production from Lipid-Producing Sugarcane. Biofuels, Bioprod. Bioref. 
    2016, 10 (3), 299–315. https://doi.org/10.1002/bbb.1640.

.. [2] Cortes-Peña, Y.; Kumar, D.; Singh, V.; Guest, J. S.
    BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and 
    Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020. https://doi.org/10.1021/acssuschemeng.9b07040.


