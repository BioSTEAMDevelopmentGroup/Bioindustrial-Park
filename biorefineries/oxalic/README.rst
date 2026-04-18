=======================================================================================
Oxalic: Production of Oxalic Acid from Sugarcane
=======================================================================================

This biorefinery is developed for Lu, et.al. [1]_ for the production of 
oxalic acid from sugarcane. It includes juicing, fermentation, 
and separation processes. Part of the script is adopted from [2]_, [3]_, and [4]_.

Installation
------------
In an environment with Python v3.9.21, do the following:
    (1) Clone this repository (this may require a few minutes' time; git clone https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park)
    (2) Add /Bioindustrial-Park/biorefineries to your Python paths
    (3) Install the following packages in sequence (this may require a few minutes' time):
	    (i) biosteam==2.46.0
	    (ii) thermosteam==0.44.0
	    (iii) contourplots==0.3.4
	    (iv) numpy==1.26.4

Getting started
---------------
Simply import the biorefinery.

.. code-block:: python

    >>> from biorefineries import oxalic
    >>> oxalic.load()
    
    >>> # Once loaded, you can have a quick glance of the results of first simulation
    >>> # (this is not baseline referenced in the main text)
    >>> # Note that the first simulation will take a longer time
    >>> oxalic.simulate_and_print()
    
    ---------- Simulation Results ----------
    MPSP is $1.77/kg AA
    GWP-100a is -1.822 kg CO2-eq/kg AA
    FEC is -49.221 MJ/kg AA
    GWP-100a without electricity credit is 4.131 kg CO2-eq/kg AA
    FEC without electricity credit is 24.270 MJ/kg AA
    ----------------------------------------
    >>> oxalic.system.feeds[0].chemicals.show()
    >>> # All chemicals used in the biorefinery
    CompiledChemicals([H2O, O2, N2, CH4, CarbonMonoxide, CO2, NH3, NitricOxide, NO2, H2S, SO2, HCl, H2SO4, HNO3, NaOH, AmmoniumHydroxide, CalciumDihydroxide, NaNO3, Na2SO4, CaSO4, MagnesiumChloride, ZincSulfate, Ethanol, CalciumLactate, CalciumOxalate, CalciumAcetate, SodiumLactate, CalciumSuccinate, AceticAcid, AcrylicAcid, Glucose, Decanol, Dodecanol, TOA, AQ336, Octanol, Hexanol, Octanediol, Toluene, Isobutyraldehyde, DPHP, GlucoseOligomer, Extract, Xylose, XyloseOligomer, Sucrose, Cellobiose, Mannose, MannoseOligomer, Galactose, GalactoseOligomer, Arabinose, ArabinoseOligomer, SolubleLignin, Protein, Enzyme, FermMicrobe, WWTsludge, Furfural, Acetoin, HMF, Xylitol, Glycerol, LacticAcid, HP, MethylHP, SuccinicAcid, OxalicAcid, MethylAcetate, EthylLactate, MethylSuccinate, Glucan, Mannan, Galactan, MEA, Xylan, Arabinan, P4O10, Tar, TiO2, CSL, BoilerChems, BaghouseBag, CoolingTowerChems, DAP, Methanol, Denaturant, DenaturedEnzyme, FermMicrobeXyl, H3PO4, Cellulose, Hemicellulose, CaO, Solids, Flocculant, Lignin, Acetate, AmmoniumSulfate, AmmoniumAcetate, Cellulase, Ash, Starch, Fiber, SolubleProtein, InsolubleProtein, TriOlein, Yeast, Octane])

Systems
-------
The `Flowsheet`, `System`, and `TEA`, objects are `flowsheet`, `system`, and `succinic_tea`, respectively.

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

    200: Juicing

    300: Fermentation

    400: Separation

    500: Wastewater

    600-900: Facilities

.. code-block:: python

    >>> # You can directly access the Flowsheet, System, and TEA objects
    >>> oxalic.load()
    >>> oxalic.simulate_and_print()
    >>> oxalic.flowsheet
    <Flowsheet: oxalic>
    >>> oxalic.system
    System: oxalic_sys
    Highest convergence error among components in recycle
    stream P203-0 after 3 loops:
    - flow rate   7.78e+00 kmol/hr (0.028%)
    - temperature 9.68e-05 K (2.6e-05%)
    

Analyses
--------
Multiple analysis modules were used to evaluate biorefineries for [1]_, including
full Monte Carlo simulation (oxalic.run_uncertainty_analysis)
and titer-yield-productivity analysis for the fermentation performance space (oxalic.run_TRY_analysis).

Note that results used in the manuscript [1]_ were generated using biosteam==2.46.0,
thermosteam==0.44.0, contourplots==0.3.4, numpy==1.26.4, and dependencies (`commit 3a8acb2 <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/commit/3a8acb269f1bf4bf30787d684d89eed23cf076ae>`_).

To reproduce the results, directly run the script of interest, and results will
be saved as Excel files and figures in analyses/results.


References
----------
.. [1] Lu et al., Metabolic Engineering of Issatchenkia orientalis to Cost-effectively Produce Bio-oxalic Acid for Rare Earth Elements Recovery. 
    Nat. Commun. 2025. Submitted May, 2025.

.. [2] Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. 
    ACS Sustainable Chem. Eng. 2021, 9 (49) 16659–16669.
    `<https://doi.org/10.1021/acssuschemeng.1c05441>`_

.. [3] Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass.
    ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. 
    `<https://doi.org/10.1021/acssuschemeng.0c08055>`_
     
.. [4] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design,
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    `<https://doi.org/10.1021/acssuschemeng.9b07040>`_
