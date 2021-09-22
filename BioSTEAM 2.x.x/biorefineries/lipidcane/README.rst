=================================================================
lipidcane: Co-production of Ethanol and Biodiesel from Lipid-cane
=================================================================

.. figure:: ./images/lipidcane_areas.png

The lipid-cane biorefinery design for the co-production of ethanol and biodiesel
follows all assumptions from a study by Huang et al. executed with SuperPro 
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

This biorefinery module makes preliminary assumptions on the lipid mass balance
throughout the process and assumes that the lipid is 100% tri-acylglyceride. 
For a more accurate assessment of lipid-cane biorefinery configuration, 
checkout the biorefineries.oilcane module, which relieves assumptions around 
the lipid mass balance and compostion. Note that the name "oilcane" is prefered 
over lipid-cane, as it resonates better with non-scientific audience and is more
consistent with how we talk about vegetable oils.

Getting Started
---------------

To load the biorefinery, simply import it and run the load function.

.. code-block:: python

    >>> from biorefineries import lipidcane as lc
    >>> # This is optional; it forces the biorefinery to load
    >>> # Otherwise, first time accessing will take a bit to load.
    >>> lc.load()
    >>> lc.chemicals # All chemicals used by the biorefinery.
    CompiledChemicals(
        [Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, N2,
         CH4, Ash, Cellulose, Hemicellulose, Flocculant, Lignin, Solids, Yeast,
         CaO, Biodiesel, Methanol, Glycerol, HCl, NaOH, NaOCH3, Phosphatidylinositol,
         OleicAcid, MonoOlein, DiOlein, TriOlein, Acetone]
    )
    >>> lc.lipidcane_sys.show() # The complete biorefinery system
    System: lipidcane_sys
    ins...
    [0] lipidcane
    [1] enzyme
    [2] H3PO4
    [3] lime
    [4] polymer
    [5] denaturant
    outs...
    [0] ethanol
    [1] biodiesel
    [2] crude_glycerol
    [3] wastewater
    [4] emissions
    [5] ash_disposal
    >>> lc.lipidcane_tea.show() # The TEA object
    ConventionalEthanolTEA: lipidcane_sys
     NPV: -0 USD at 20.7% IRR
    >>> lc.flowsheet # The complete flowsheet
    <Flowsheet: lipidcane>
    >>> lc.R301.show() # Unit operations and streams can be accessed through the module
    Fermentation: R301
    ins...
    [0] s164  from  HXutility-H301
        phase: 'l', T: 295.15 K, P: 101325 Pa
        flow (kmol/hr): Water    4.34e+03
                        Glucose  10.5
                        Sucrose  62.7
                        H3PO4    0.85
    [1] s170  from  MockSplitter-S302
        phase: 'l', T: 294.19 K, P: 101325 Pa
        flow (kmol/hr): Water  12.6
                        Yeast  13.7
    outs...
    [0] CO2  to  VentScrubber-D301
        phase: 'g', T: 294.19 K, P: 101325 Pa
        flow (kmol/hr): Water    5.35
                        Ethanol  2.75
                        CO2      245
    [1] s165  to  StorageTank-T301
        phase: 'l', T: 294.19 K, P: 101325 Pa
        flow (kmol/hr): Water    4.28e+03
                        Ethanol  242
                        Glucose  4.08
                        H3PO4    0.85
                        Yeast    82.8


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


