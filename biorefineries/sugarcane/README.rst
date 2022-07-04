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
    >>> sc.load() # This is optional; it forces the biorefinery to load
    >>> sc.chemicals # All chemicals used by the biorefinery.
    CompiledChemicals([Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, N2, CH4, Ash, Cellulose, Hemicellulose, Flocculant, Lignin, Solids, Yeast, CaO])
    >>> sc.sugarcane_sys.show(data=False) # The complete biorefinery system
    System: sugarcane_sys
    ins...
    [0] sugarcane
    [1] enzyme
    [2] H3PO4
    [3] lime
    [4] polymer
    [5] denaturant
    outs...
    [0] ethanol
    [1] vinasse
    [2] wastewater
    [3] emissions
    [4] ash_disposal
    >>> sc.sugarcane_tea.show() # The TEA object
    ConventionalEthanolTEA: sugarcane_sys
     NPV: -3 USD at 11.9% IRR
    >>> sc.flowsheet # The complete flowsheet
    <Flowsheet: sugarcane>
    >>> sc.R301.show() # Any unit operations and streams can be accessed through the module
    Fermentation: R301
    ins...
    [0] s35  from  HXutility-H301
        phase: 'l', T: 295.15 K, P: 101325 Pa
        flow (kmol/hr): Water    8.71e+03
                        Glucose  21.1
                        Sucrose  126
                        H3PO4    0.85
    [1] s41  from  MockSplitter-S302
        phase: 'l', T: 294.19 K, P: 101325 Pa
        flow (kmol/hr): Water  25.3
                        Yeast  27.5
    outs...
    [0] CO2  to  VentScrubber-D301
        phase: 'g', T: 294.19 K, P: 101325 Pa
        flow (kmol/hr): Water    10.7
                        Ethanol  5.52
                        CO2      491
    [1] s36  to  StorageTank-T301
        phase: 'l', T: 294.19 K, P: 101325 Pa
        flow (kmol/hr): Water    8.6e+03
                        Ethanol  486
                        Glucose  8.19
                        H3PO4    0.85
                        Yeast    166


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


