==================================================================================
ethanol_adipic: Production of Ethanol and Adipic Acid from Lignocellulosic Biomass
==================================================================================

Simplified process flow schemes
-------------------------------
.. figure:: ./images/acid_pretreatment_biorefinery.png

    *Acid-pretreatment biorefinery*


.. figure:: ./images/AFEX_pretreatment_biorefinery.png

    *AFEX-pretreatment biorefinery*


.. figure:: ./images/base_pretreatment_biorefinery.png

    *Base-pretreatment biorefinery*


This module contains several options for preprocessing and pretreatment of
lignocellulosic biomass. Preprocessing options at biomass depot can be either
conventional or high-moisture pelleting process (CPP/HMPP). Pretreatment of the
feedstocks can be ammonia fiber expansion (AFEX, at the biomass depot), acid
(dilute sulfuric acid, at the biorefinery), or base (deacetylation and mechanical
refining, DMR, at the biorefinery).

Carbohydrates in these biorefineries are converted to ethanol; unconverted
carbohydrates and lignin in all biorefineries can be burned in the boiler for
energy production. When AFEX- and base-pretreatment are used, it is also
feasible to include a process to convert those into adipic acid (default).


Getting Started
---------------
First choose which system configuration to load, the default is acid-pretreatment
biorefinery preprocessed by HMPP.

.. code-block:: python

    >>> from biorefineries import ethanol_adipic as ea
    >>> ea.load_system('acid', 'HMPP')
    >>> # Once loaded, you can have a quick glance of the results
    >>> # Note that the first simulation will take a longer time
    >>> ea.simulate_and_print()
    
    ---------- Acid Biorefinery ----------
    MESP: $3.40/gal
    GWP: 1.016 kg CO2-eq/gal ethanol without feedstock
    --------------------------------------
    >>> ea.chems.show()
    CompiledChemicals([H2O, O2, N2, CH4, CO, CO2, NH3, NO, NO2, H2S, SO2, H2, H2SO4, HNO3, NH4OH, CalciumDihydroxide, AmmoniumSulfate, NaNO3, CaSO4, NaOH, Na2SO4, DAP, AceticAcid, Glucose, GlucoseOligomer, Extractives, Xylose, XyloseOligomer, Sucrose, Cellobiose, Mannose, MannoseOligomer, Galactose, GalactoseOligomer, Arabinose, ArabinoseOligomer, SolubleLignin, Protein, Enzyme, FermMicrobe, WWTsludge, Furfural, HMF, Xylitol, LacticAcid, SuccinicAcid, Ethanol, Glycerol, P_putida, P_putidaGrow, Denaturant, AdipicAcid, MuconicAcid, MonoSodiumMuconate, Acetate, AmmoniumAcetate, Glucan, Mannan, Galactan, Xylan, Arabinan, Lignin, P4O10, Ash, Tar, CSL, BoilerChems, Polymer, BaghouseBag, CoolingTowerChems])
    >>> # To load the alternative configuration, you can just use the load_system function
    >>> ea.load_system('AFEX', 'CPP_AFEX')
    >>> # Note the default GWP doesn't include impacts from feedstocks,
    >>> # you can include feedstock impacts by indicating the which depot
    >>> # configuration is used for feedstock preprocessing
    >>> ea.simulate_and_print('CPP_AFEX')
    
    ---------- AFEX Biorefinery ----------
    MESP: $2.83/gal
    GWP: -3.583 kg CO2-eq/gal ethanol without feedstock
    GWP: -1.418 kg CO2-eq/gal ethanol with feedstock
    --------------------------------------


Systems
-------
The `Flowsheet`, `System`, and `TEA`, objects are `flowsheet`, `tea`, and `biorefinery`, respectively.

Naming conventions:
    D = Distillation column
    
    F = Flash tank
    
    H = Heat exchange
    
    M = Mixer
    
    P = Pump
    
    R = Reactor
    
    S = Splitter (including solid/liquid separator)
    
    T = Tank or bin for storage
    
    U = Other units
    
    PS = Process specification, not physical units, but for adjusting streams

Processes:
    100: Feedstock preprocessing
    
    200: Pretreatment
    
    300: Carbohydrate conversion
    
    400: Carbohydrate product separation
    
    500: Wastewater treatment
    
    600: Facilities
    
    700: Lignin conversion and separation


.. code-block:: python

    >>> # You can directly access the Flowsheet, System, and TEA objects
    >>> ea.load_system('acid', 'HMPP')
    >>> ea.simulate_and_print()
    >>> ea.flowsheet
    <Flowsheet: acid>
    >>> ea.biorefinery.show()
    System: biorefinery
    ins...
    [0] water_M201
        phase: 'l', T: 387.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  2.04e+03
    [1] feedstock
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O          1.16e+03
                        Extractives  62.2
                        Sucrose      1.87
                        Protein      113
                        Acetate      25.1
                        Glucan       180
                        Mannan       3.08
                        ...
    [2] water_M202
        phase: 'l', T: 368.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  7.82e+03
    [3] steam_M203
        phase: 'g', T: 506.15 K, P: 1.041e+06 Pa
        flow (kmol/hr): H2O  1.75e+03
    [4] water_M205
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  723
    [5] enzyme_M301
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Enzyme  23.9
    [6] water_M301
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  7.69e+03
    [7] water_U401
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  1.43e+03
    [8] caustic_R502  from  ReversedSplitter-T605_S
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NaOH  56.3
    [9] polymer_R502
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Polymer  2.42
    [10] air_R502
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  2.07e+03
                        N2  7.79e+03
    [11] denaturant
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Denaturant  4.57
    [12] sulfuric_acid
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O    7.7
                        H2SO4  18.8
    [13] ammonia
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NH4OH  103
    [14] caustic
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [15] CSL
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CSL  36.2
    [16] DAP
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): DAP  1.19
    [17] firewater_in
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  445
    [18] plant_air_in
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  9.01e+03
                        N2  3.87e+04
    [19] lime_CHP
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CalciumDihydroxide  22.2
    [20] boiler_chems
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): BoilerChems  0.00135
    [21] baghouse_bag
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): BaghouseBag  0.722
    [22] natural_gas
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [23] boiler_makeup_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  440
    [24] CIP_chems_in
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  8.05
    [25] return_cooling_water
        phase: 'l', T: 310.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  4.11e+05
    [26] cooling_tower_chems
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CoolingTowerChems  1.37
    [27] CT_makeup_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  637
    [28] system_makeup_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  769
    outs...
    [0] U401_vent
        phase: 'g', T: 304.51 K, P: 101325 Pa
        flow (kmol/hr): O2   5.22
                        CO2  456
    [1] vent_R502
        phase: 'g', T: 308.05 K, P: 101325 Pa
        flow (kmol/hr): H2O  257
                        O2   2.01e+03
                        N2   7.79e+03
                        CO2  59.6
    [2] S506_vent
        phase: 'g', T: 373.12 K, P: 101325 Pa
        flow (kmol/hr): H2O  274
                        SO2  0.00755
    [3] residuals_to_CHP
        phase: 'l', T: 373.12 K, P: 101325 Pa
        flow (kmol/hr): NaNO3              3.52
                        NaOH               51.6
                        DAP                1.94
                        AceticAcid         0.0271
                        Glucose            5.22e-09
                        GlucoseOligomer    0.00358
                        Extractives        0.0322
                        ...
    [4] ethanol
        phase: 'l', T: 349.14 K, P: 101325 Pa
        flow (kmol/hr): H2O         6.85
                        Ethanol     464
                        Denaturant  4.57
    [5] s292  to  ReversedSplitter-T605_S
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [6] s296
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [7] s298
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [8] firewater_out
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  445
    [9] sodium_sulfate
        phase: 's', T: 373.12 K, P: 101325 Pa
        flow: 0
    [10] plant_air_out
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  9.01e+03
                        N2  3.87e+04
    [11] vent_CHP
        phase: 'g', T: 539.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  1.36e+03
                        N2   49.1
                        CO2  1.76e+03
                        NH3  63.2
                        SO2  1.55
    [12] ash
        phase: 's', T: 539.15 K, P: 101325 Pa
        flow (kmol/hr): CalciumDihydroxide  3.69
                        AmmoniumSulfate     0.0558
                        NaNO3               0.0721
                        CaSO4               18.5
                        NaOH                1.06
                        DAP                 0.129
                        Ash                 73.2
                        ...
    [13] CIP_chems_out
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  8.05
    [14] process_cooling_water
        phase: 'l', T: 301.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  4.11e+05
    [15] process_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  2.25e+04
    [16] discharged_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    >>> ea.tea.show()
    EthanolAdipicTEA: biorefinery
     NPV: -579 USD at 10.0% IRR
    >>> # You can use the flowsheet to access streams, unit operations, and subsystems
    >>> u = ea.flowsheet.unit
    >>> u.R301.show()
    SaccharificationAndCoFermentation: R301
    ins...
    [0] s268  from  EnzymeHydrolysateMixer-M301
        phase: 'l', T: 346.95 K, P: 101325 Pa
        flow (kmol/hr): H2O                1.96e+04
                        NH4OH              3.76
                        AmmoniumSulfate    18.8
                        DAP                0.981
                        AceticAcid         20.6
                        Glucose            19.7
                        GlucoseOligomer    0.54
                        ...
    [1] s269  from  SeedHoldTank-T301
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                2.16e+03
                        NH4OH              0.417
                        AmmoniumSulfate    2.09
                        DAP                0.43
                        AceticAcid         2.28
                        Glucose            0.832
                        GlucoseOligomer    0.78
                        ...
    [2] CSL_R301  from  ReversedSplitter-T606_S
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CSL  29.6
    [3] DAP_R301  from  ReversedSplitter-T607_S
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): DAP  0.981
    outs...
    [0] R301_g  to  Mixer-M401
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): H2O           19.7
                        O2            4.68
                        CO2           411
                        AceticAcid    0.0405
                        Furfural      0.0223
                        HMF           0.0095
                        LacticAcid    1.14e-05
                        ...
    [1] effluent  to  Mixer-M402
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                1.94e+04
                        NH4OH              3.76
                        AmmoniumSulfate    18.8
                        DAP                2.07
                        AceticAcid         20.6
                        Glucose            2.91
                        GlucoseOligomer    7.02
                        ...
    [2] side_draw  to  SeedTrain-R302
        phase: 'l', T: 321.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                2.16e+03
                        NH4OH              0.417
                        AmmoniumSulfate    2.09
                        DAP                0.448
                        AceticAcid         2.28
                        Glucose            16.6
                        GlucoseOligomer    0.78
                        ...
