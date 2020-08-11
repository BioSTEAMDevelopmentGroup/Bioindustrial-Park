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
    >>> # First time accessing takes a bit to load the chemicals and system.
    >>> sc.chemicals # All chemicals used by the biorefinery.
    CompiledChemicals([Water, Methanol, Ethanol, Glycerol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, Biodiesel, CH4, Ash, Cellulose, Hemicellulose, Flocculant, Lignin, Solids, DryYeast, CaO, HCl, NaOH, NaOCH3, Lipid])
    >>> sc.sugarcane_sys.show() # The complete biorefinery system
    System: lipidcane_sys
     path: (U101, U102, U103, T201, SYS1, T202,
            H201, T203, P201, T204, T205,
            P202, SYS2, T207, T207_2, H203,
            T208, C203, F201, T403, P403,
            R401, C401, R402, C402, SYS3,
            F401, H401, P408, P407, T409,
            P405, B401, H403, P411, T401,
            P401, T402, P402, T404, P404,
            S401, S202, S301, F301, P306,
            M301, H301, T305, R301, T301,
            C301, M302, P301, SYS4, SYS5,
            H304, T302, P304, T303, P305,
            T304, D301, T408, P303, M305, U202)
     facilities: (BT, CWP, CT, PWC)
    >>> sc.sugarcane_tea.show() # The TEA object
    ConventionalEthanolTEA: lipidcane_sys
     NPV: 1.84e+07 USD at 17.8% IRR
     ROI: 0.176 1/yr
     PBP: 5.4 yr
    >>> sc.flowsheet # The complete flowsheet
    <Flowsheet: lipidcane>
    >>> sc.R301.show() # Any unit operations and streams can be accessed through the module
    Fermentation: R301
    ins...
    [0] d235  from  HXutility-H301
        phase: 'l', T: 295.15 K, P: 101325 Pa
        flow (kmol/hr): Water    4.11e+03
                        Glucose  10.5
                        Sucrose  62.7
                        H3PO4    0.85
    [1] d259  from  MixTank-T305
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     1.37e+03
                        DryYeast  1.03e+04
    outs...
    [0] CO2  to  VentScrubber-D301
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    1.83
                        Ethanol  0.456
                        CO2      245
    [1] d236  to  StorageTank-T301
        phase: 'l', T: 296.06 K, P: 101325 Pa
        flow (kmol/hr): Water     5.42e+03
                        Ethanol   244
                        Glucose   13.6
                        H3PO4     0.85
                        DryYeast  1.03e+04


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


