=================================================================
sugarcane: Ethanol production from Sugarcane
=================================================================

The sugarcane biorefinery design for the production of ethanol follows all 
assumptions from a study by Huang et. al. executed with SuperPro 
Designer [1]_. As discussed in the original BioSTEAM manuscript [2]_, the 
sugarcane biorefinery can be divided into five main areas: feedstock handling, 
juice separation, ethanol production, boiler and turbo-generator system, and 
utilities. 

Getting Started
---------------

To load the biorefinery, simply import it. All data and variables
are lazy loaded by the module:

.. code-block:: python

    >>> from biorefineries import sugarcane as sc
    >>> # This is optional; it forces the biorefinery to load
    >>> # Otherwise, first time accessing will take a bit to load.
    >>> sc.load()
    >>> sc.chemicals # All chemicals used by the biorefinery.
    CompiledChemicals([Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, CH4, Ash, Cellulose, Hemicellulose, Flocculant, Lignin, Solids, Yeast, CaO])
    >>> sc.sugarcane_sys.show() # The complete biorefinery system
    System: sugarcane_sys
     path: (U101, U102, U103, T201,
            juice_extraction_sys, T202, H201,
            T203, P201, T204, T205, P202,
            juice_separation_sys, S202,
            ethanol_production_sys, U202)
     facilities: (CWP, BT, CT, PWC)
    >>> sc.sugarcane_tea.show() # The TEA object
    ConventionalEthanolTEA: sugarcane_sys
     NPV: 0 USD at 10.3% IRR
    >>> sc.flowsheet # The complete flowsheet
    <Flowsheet: sugarcane>
    >>> sc.R301.show() # Any unit operations and streams can be accessed through the module
    Fermentation: R301
    ins...
    [0] s30  from  HXutility-H301
        phase: 'l', T: 295.15 K, P: 101325 Pa
        flow (kmol/hr): Water    4.34e+03
                        Glucose  21.1
                        Sucrose  126
                        H3PO4    0.85
    [1] s55  from  MixTank-T305
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  1.37e+03
                        Yeast  1.03e+04
    outs...
    [0] CO2  to  VentScrubber-D301
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    4.1
                        Ethanol  1.73
                        CO2      491
    [1] s31  to  StorageTank-T301
        phase: 'l', T: 295.98 K, P: 101325 Pa
        flow (kmol/hr): Water    5.58e+03
                        Ethanol  489
                        Glucose  27.3
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


