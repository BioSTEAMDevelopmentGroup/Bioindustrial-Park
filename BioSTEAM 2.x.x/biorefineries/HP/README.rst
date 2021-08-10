=======================================================================================
HP: Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass
=======================================================================================

Process flow scheme
------------------------------
.. figure:: ./images/system_diagram.png


The biorefinery is developed for Bhagwat et al. [1]_ for the production of 
acrylic acid via 3-hydroxypropionic acid (3-HP) from lignocellulosic biomass
in a biorefinery comprising sulfuric acid pretreatment, fermentation, 
and separation processes. Part of the script is adopted from [2]_ and [3]_.


Getting started
---------------
Two configurations are now included in this biorefinery, one with corn stover
('lignocellulosic') and one with sugarcane ('sugarcane') as the feedstock.
You can choose which system to load.

.. code-block:: python

    >>> from biorefineries import HP
    >>> HP.load_system('lignocellulosic')
    >>> # Once loaded, you can have a quick glance of the results
    >>> # Note that the first simulation will take a longer time
    >>> HP.simulate_and_print()
    
    ---------- Simulation Results ----------
    MPSP is $1.831/kg
    GWP is 3.898 kg CO2-eq/kg AA
    FEC is 51.43 MJ/kg AA or 1.64 MJ/MJ AA
    ----------------------------------------
    >>> HP.chems.show()
    >>> # All chemicals used in the bioreinfery
    CompiledChemicals([H2O, O2, N2, CH4, CarbonMonoxide, CO2, NH3, NitricOxide, NO2, H2S, SO2, HCl, H2SO4, HNO3, NaOH, AmmoniumHydroxide, CalciumDihydroxide, AmmoniumSulfate, NaNO3, Na2SO4, CaSO4, Ethanol, AmmoniumAcetate, CalciumLactate, CalciumAcetate, CalciumSuccinate, AceticAcid, Acetate, AcrylicAcid, Glucose, Decanol, TOA, AQ336, Octanol, Hexanol, 3-methyl-1-butanol, Octanediol, Butyl acetate, Toluene, Isobutyraldehyde, DPHP, GlucoseOligomer, Extract, Xylose, XyloseOligomer, Sucrose, Cellobiose, Mannose, MannoseOligomer, Galactose, GalactoseOligomer, Arabinose, ArabinoseOligomer, SolubleLignin, Protein, Enzyme, FermMicrobe, WWTsludge, Furfural, Acetoin, HMF, Xylitol, Glycerol, LacticAcid, HP, MethylHP, SuccinicAcid, MethylAcetate, EthylLactate, MethylSuccinate, Glucan, Mannan, Galactan, Xylan, Arabinan, Lignin, P4O10, Ash, Tar, TiO2, CSL, BoilerChems, BaghouseBag, CoolingTowerChems, DAP, Methanol, Denaturant, DenaturedEnzyme, FermMicrobeXyl, H3PO4, Hemicellulose, Flocculant, Solids, Yeast, CaO])
    >>> # To load an alternative feedstock biorefinery, you can just use the load_system function
    >>> HP.load_system('sugarcane')
    >>> HP.simulate_and_print()
    
    ---------- Simulation Results ----------
    MPSP is $1.912/kg
    GWP is 3.129 kg CO2-eq/kg AA
    FEC is 40.61 MJ/kg AA or 1.29 MJ/MJ AA
    ----------------------------------------

Systems
-------
The `Flowsheet`, `System`, and `TEA`, objects are `flowsheet`, `system`, and `HP_tea`, respectively.

Naming conventions:
    D = Distillation column

    E = Evaporator
    
    F = Flash tank

    H = Heat exchange

    M = Mixer

    P = Pump (including conveying belt)

    R = Reactor

    S = Splitter (including solid/liquid and liquid/liquid separator)

    T = Tank or bin for storage

    U = Other units

    PS = Process specificiation, not physical units, but for adjusting streams

Processes:
    100: Preprocessing

    200: Pretreatment

    300: Conversion

    400: Separation

    500: Wastewater

    600: Facilities

.. code-block:: python

    >>> # You can directly access the Flowsheet, System, and TEA objects
    >>> HP.load_system('lignocellulosic')
    >>> HP.simulate_and_print()
    >>> HP.flowsheet
    <Flowsheet: HP>
    >>> HP.system
    System: HP_sys
    >>> HP.HP_tea.show()
    CellulosicEthanolTEA: HP_sys
     NPV: -0 USD at 10.0% IRR
    >>> # You can use the flowsheet to access streams, unit operations, and subsystems
    >>> u = HP.flowsheet.unit
    >>> u.R302.show()
    CoFermentation: R302
    ins...
    [0] to_cofermentation  from  Splitter-S302
        phase: 'l', T: 303.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                1.78e+04
                        AmmoniumHydroxide  3.36
                        AmmoniumSulfate    16.8
                        AceticAcid         20.9
                        Glucose            148
                        GlucoseOligomer    6.26
                        Extract            60.7
                        ...
    [1] s267  from  SeedHoldTank-T301
        phase: 'l', T: 303.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                1.39e+03
                        AmmoniumHydroxide  0.253
                        AmmoniumSulfate    1.27
                        AceticAcid         3.44
                        Glucose            4.91e-09
                        GlucoseOligomer    0.471
                        Extract            4.57
                        ...
    [2] CSL  from  CSLstorageTank-T603
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CSL  101
    [3] fermentation_lime  from  LimeStorageBin-T604
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CalciumDihydroxide  149
    outs...
    [0] fermentation_effluent  to  SolidsCentrifuge-S401
        phase: 'l', T: 303.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                 2e+04
                        AmmoniumHydroxide   3.62
                        CalciumDihydroxide  2.6
                        AmmoniumSulfate     18.1
                        CalciumLactate      120
                        CalciumAcetate      26
                        AceticAcid          5.19e-08
                        ...
    [1] CO2_fermentation
        phase: 'g', T: 303.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  524

    

Analyses
--------
Multiple analysis modules were used to evaluate the biorefinery's
different aspects for [1]_, including: full Monte Carlo simulation (1_full_evaluation.py in ./analyses),
titer-yield-productivity analysis for the fermentation performance space (TRY_analysis.py), targeted improvements (system_targeted_improvements.py),
and evaluating feedstocks of varying carbohydrate contents, sugar contents, and prices (TRY_analysis.py).

Note that results used in the manuscript [1]_ were generated using biosteam (HP_sys branch),
thermosteam (HP_sys branch), and dependencies (`commit fcc4b50 <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/commit/fcc4b5080c097f2fa74804344b11c02ed51399b8>`_).

To reproduce the results, directly run the script of interest, and results will
be saved as Excel files in the same directory path as the module.


References
----------

.. [1] Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. 
    ACS Sustainable Chem. Eng. 2021.
    Submitted August 10, 2021.

.. [2] Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass.
    ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. 
    `<https://doi.org/10.1021/acssuschemeng.0c08055>`_
     
.. [3] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design,
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    `<https://doi.org/10.1021/acssuschemeng.9b07040>`_







