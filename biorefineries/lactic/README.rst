==============================================================
lactic: Production of Lactic Acid from Lignocellulosic Biomass
==============================================================

The biorefinery is developed for Li et al. [1]_ for the production of lactic acid
via dilute sulfuric acid pretreatment, fermentation, and separation processes
from lignocellulosic biomass. Part of the script is adopted from [2]_.


Getting Started
---------------
Two configurations are now included in this biorefinery, one with simultaneous
saccharification and co-fermentation (SSCF), one with separate hydrolysis and
co-fermentation (SHF). You can choose which system to load.

.. code-block:: python

    >>> # Loading the biorefinery is easy and you can have a quick glance of the results
    >>> # Note that the first simulation will take a longer time
    >>> from biorefineries import lactic as la
    >>> la.load(kind='SSCF', print_results=True) # `kind` can be "SSCF" or "SHF"

    ---------- Simulation Results ----------
    MPSP is $1.414/kg
    GWP is 4.567 kg CO2-eq/kg lactic acid
    FEC is 58.67 MJ/kg lactic acid
    ------------------------------------------
    >>> la.chemicals.show()
    >>> # All chemicals used in the biorefinery
    CompiledChemicals([H2O, O2, N2, CH4, CO, CO2, NH3, NO, NO2, H2S, SO2, H2SO4, HNO3, NaOH, NH4OH, CalciumDihydroxide, AmmoniumSulfate, NaNO3, Na2SO4, CaSO4, Ethanol, AceticAcid, Glucose, GlucoseOligomer, Extractives, Xylose, XyloseOligomer, Sucrose, Cellobiose, Mannose, MannoseOligomer, Galactose, GalactoseOligomer, Arabinose, ArabinoseOligomer, SolubleLignin, Protein, Enzyme, FermMicrobe, WWTsludge, Furfural, HMF, Xylitol, LacticAcid, SuccinicAcid, EthylAcetate, EthylLactate, EthylSuccinate, Acetate, AmmoniumAcetate, CalciumLactate, CalciumAcetate, CalciumSuccinate, Glucan, Mannan, Galactan, Xylan, Arabinan, Lignin, P4O10, Ash, Tar, CSL, BoilerChems, Polymer, BaghouseBag, CoolingTowerChems])
    >>> # To load the alternative configuration, you can just use the load function
    >>> la.load('SHF')

    ---------- Simulation Results ----------
    MPSP is $1.454/kg
    GWP is 4.489 kg CO2-eq/kg lactic acid
    FEC is 57.52 MJ/kg lactic acid
    ------------------------------------------


Systems
-------
The `Flowsheet`, `System`, and `TEA`, objects are `flowsheet`, `lactic_sys`, and `lactic_tea`, respectively.

Naming conventions:
    D = Distillation column

    E = Evaporator

    F = Flash tank

    H = Heat exchange

    M = Mixer

    P = Pump (including conveying belt)

    R = Reactor

    S = Splitter (including solid/liquid separator)

    T = Tank or bin for storage

    U = Other units

    PS = Process specification, not physical units, but for adjusting streams

Processes:
    100: Preprocessing

    200: Pretreatment

    300: Conversion

    400: Separation

    500: Wastewater

    600: Facilities

.. code-block:: python

    >>> # You can directly access the Flowsheet, System, and TEA objects
    >>> la.flowsheet
    <Flowsheet: SHF>
    >>> la.lactic_sys.show()
    System: lactic_sys
    ins...
    [0] feedstock
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O          1.16e+03
                        Extractives  62.2
                        Sucrose      1.87
                        Protein      113
                        Acetate      25.1
                        Glucan       180
                        Mannan       3.08
                        ...
    [1] water_M201
        phase: 'l', T: 387.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  2.05e+03
    [2] water_M202
        phase: 'l', T: 368.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  7.83e+03
    [3] steam_M203
        phase: 'g', T: 506.15 K, P: 1.041e+06 Pa
        flow (kmol/hr): H2O  1.75e+03
    [4] water_M205
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  723
    [5] enzyme_M301
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Enzyme  24
    [6] water_M301
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  7.69e+03
    [7] water_R301
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [8] s98
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [9] water_R403
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [10] s99
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [11] caustic_R502
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NaOH  56.3
    [12] polymer_R502
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Polymer  2.14
    [13] air_R502
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  1.83e+03
                        N2  6.9e+03
    [14] sulfuric_acid
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O    104
                        H2SO4  255
    [15] ammonia
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NH4OH  103
    [16] CSL
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CSL  92.3
    [17] lime
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CalciumDihydroxide  243
    [18] ethanol
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  13.6
    [19] firewater_in
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  445
    [20] plant_air_in
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  9.01e+03
                        N2  3.87e+04
    [21] lime_CHP
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CalciumDihydroxide  1.68
    [22] boiler_chems
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): BoilerChems  0.00433
    [23] baghouse_bag
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): BaghouseBag  2.32
    [24] natural_gas
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CH4  1.25e+03
    [25] boiler_makeup_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  1.15e+03
    [26] CIP_chems_in
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  8.05
    [27] return_cooling_water
        phase: 'l', T: 310.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  1.79e+06
    [28] cooling_tower_chems
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CoolingTowerChems  5.94
    [29] CT_makeup_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  2.77e+03
    [30] system_makeup_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  4.91e+03
    outs...
    [0] gypsum
        phase: 'l', T: 356.64 K, P: 101325 Pa
        flow (kmol/hr): H2O                471
                        H2SO4              10.7
                        AmmoniumSulfate    0.712
                        CaSO4              212
                        AceticAcid         30.1
                        Glucose            0.788
                        GlucoseOligomer    0.241
                        ...
    [1] vent_R502
        phase: 'g', T: 306.67 K, P: 101325 Pa
        flow (kmol/hr): H2O  227
                        O2   1.77e+03
                        N2   6.9e+03
                        CO2  63.3
    [2] brine
        phase: 'l', T: 306.67 K, P: 101325 Pa
        flow (kmol/hr): H2O                251
                        NaOH               55.1
                        CaSO4              7.15e-09
                        AceticAcid         0.0632
                        Glucose            2.98e-08
                        GlucoseOligomer    0.00254
                        Extractives        0.0229
                        ...
    [3] lactic_acid
        phase: 'l', T: 345 K, P: 101325 Pa
        flow (kmol/hr): H2O           86.1
                        Ethanol       0.538
                        AceticAcid    0.155
                        Furfural      0.185
                        HMF           0.00504
                        LacticAcid    256
                        EthylLactate  13.1
    [4] firewater_out
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  445
    [5] plant_air_out
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  9.01e+03
                        N2  3.87e+04
    [6] vent_CHP
        phase: 'g', T: 539.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  3.98e+03
                        N2   44.7
                        CO2  3.14e+03
                        NH3  64.4
                        SO2  1.89
    [7] ash
        phase: 's', T: 539.15 K, P: 101325 Pa
        flow (kmol/hr): NaOH                1.17
                        CalciumDihydroxide  0.279
                        AmmoniumSulfate     0.11
                        CaSO4               23.8
                        CalciumLactate      7.04
                        CalciumAcetate      1.16
                        Ash                 73.2
                        ...
    [8] boiler_blowdown
        phase: 'l', T: 373.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  1.14e+03
    [9] CIP_chems_out
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  8.05
    [10] process_cooling_water
        phase: 'l', T: 301.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  1.79e+06
    [11] cooling_tower_blowdown
        phase: 'l', T: 301.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  2.77e+03
    [12] process_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  2.39e+04
    [13] discharged_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    >>> la.lactic_tea.show()
    LacticTEA: lactic_sys
     NPV: -850 USD at 10.0% IRR
    >>> # You can directly access streams, unit operations, and subsystems
    >>> # or use the flowsheet
    >>> # Flowsheet is recommended when you have switched configurations
    >>> # to avoid retrieving objects from the preceding configuration
    >>> la.R301 is la.flowsheet.unit.R301
    True
    >>> la.R301.show()
    CoFermentation: R301
    ins...
    [0] to_fermenter  from  Splitter-S302
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                1.87e+04
                        NH4OH              3.62
                        AmmoniumSulfate    18.1
                        AceticAcid         10.3
                        Glucose            159
                        GlucoseOligomer    6.73
                        Extractives        59.9
                        ...
    [1] s68  from  SeedHoldTank-T301
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                1.41e+03
                        NH4OH              0.272
                        AmmoniumSulfate    1.36
                        AceticAcid         4.53
                        Glucose            1.93
                        GlucoseOligomer    0.507
                        Extractives        4.51
                        ...
    [2] CSL_R301  from  CSLstorage-T604
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CSL  92.3
    [3] lime_R301  from  LimeStorage-T605
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CalciumDihydroxide  243
    [4] water_R301
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [5] s67  from  Pump-E301_P
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow: 0
    outs...
    [0] fermentation_effluent  to  Pump-R301_P1
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                 1.92e+04
                        NH4OH               3.62
                        CalciumDihydroxide  22.1
                        AmmoniumSulfate     18.1
                        Glucose             22.4
                        GlucoseOligomer     6.73
                        Extractives         59.9
                        ...
    [1] sidedraw  to  Pump-R301_P2
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                1.41e+03
                        NH4OH              0.272
                        AmmoniumSulfate    1.36
                        AceticAcid         1.04
                        Glucose            11.3
                        GlucoseOligomer    0.507
                        Extractives        4.51
                        ...


Analyses
--------
Multiple analysis modules (in `./analyses <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/f56692d3bc06527b57dc77ed7cb929a40b59bc4d/BioSTEAM%202.x.x/biorefineries/lactic/analyses>`_) were used to evaluate the biorefinery from different aspects for [1]_, including: full Monte Carlo simulation, titer-yield-productivity analysis for the fermentation performance space, and evaluate feedstocks of varying carbohydrate contents and prices.

Analyses (in ./analyses) performed for the paper [1]_ and archived results (see `commit f56692d <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/commit/f56692d3bc06527b57dc77ed7cb929a40b59bc4d>`_) were generated using biosteam v2.20.21, thermosteam v0.20.26, and dependencies.

To reproduce the results, set up the environment to be using the biosteam, thermosteam, and this biorefinery at `commit f56692d <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/commit/f56692d3bc06527b57dc77ed7cb929a40b59bc4d>`_, then directly run the script of interest, results will be saved as Excel files in the same directory path as the module.


References
----------
.. [1] Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass.
    ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351.
     `<https://doi.org/10.1021/acssuschemeng.0c08055>`_

.. [2] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design,
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty.
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
    `<https://doi.org/10.1021/acssuschemeng.9b07040>`_