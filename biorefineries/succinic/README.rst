=======================================================================================
Succinic: Production of Succinic Acid from Sugarcane
=======================================================================================

Process flowsheet
------------------------------
.. figure:: ./images/system_diagram.png


These biorefineries were developed to support studies by Tran, Mishra, Bhagwat et al. [1]_ and Tran et al. [3]_ for the production of 
succinic acid from sugarcane. They include juicing, fermentation, 
and separation processes. Part of the script is adopted from [2]_.

Installation
------------
In an environment with Python v3.9.15, do the following:
    (1) Clone this repository (this may require a few minutes' time; git clone https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park)
    (2) Add /Bioindustrial-Park/biorefineries to your Python paths
    (3) Install the following packages in sequence (this may require a few minutes' time):
	    (i) biosteam==2.36.4
	    (ii) thermosteam==0.33.0
	    (iii) contourplots==0.2.1
	    (iv) numpy==1.23.5

Getting started
---------------
Simply import the biorefinery and choose from three scenarios to load ('lab_batch', 'lab_fed-batch', 'pilot_batch').

.. code-block:: python

    >>> from biorefineries import succinic
    >>> succinic.load('pilot_batch')
    >>> # Once loaded, you can have a quick glance of the results
    >>> # Note that the first simulation will take a longer time
    >>> succinic.simulate_and_print()
    
    ---------- Simulation Results ----------
    MPSP is $1.367/kg
    GWP100a is 1.673 kg CO2-eq./kg
    FEC is -0.210 MJ/kg
    ----------------------------------------
    >>> succinic.system.chemicals.show()
    >>> # All chemicals used in the bioreinfery
    CompiledChemicals([H2O, H2, CH4, CarbonMonoxide, CO2, NH3, NitricOxide, NO2, H2S, SO2, H2SO4, HCl, HNO3, NaOH, KOH, KCl, AmmoniumHydroxide, CalciumDihydroxide, DiammoniumSulfate, MagnesiumSulfate, NaNO3, Na2SO4, CaSO4, Ethanol, Acetate, AmmoniumAcetate, DiammoniumSuccinate, CalciumLactate, CalciumAcetate, CalciumSuccinate, AceticAcid, Glucose, Fructose, GlucoseOligomer, Extract, Xylose, XyloseOligomer, Sucrose, Cellobiose, Mannose, MannoseOligomer, Galactose, GalactoseOligomer, Arabinose, ArabinoseOligomer, SolubleLignin, Protein, Enzyme, FermMicrobe, WWTsludge, Furfural, Hexanol, HMF, PotassiumSorbate, TAL, Pyrone, DHL, PyruvicAcid, SuccinicAcid, LacticAcid, SorbicAcid, PSA, PolyPSA, ButylSorbate, HMTHP, HMDHP, PD, VitaminA, VitaminD2, Xylitol, Glucan, Mannan, Galactan, MEA, Xylan, Arabinan, Lignin, P4O10, Ash, Tar, CSL, BoilerChems, BaghouseBag, CoolingTowerChems, DAP, Methanol, Denaturant, DenaturedEnzyme, MethylLactate, FermMicrobeXyl, H3PO4, Cellulose, Hemicellulose, CaO, Solids, Flocculant, CO2_compressible, O2, N2])
    >>> # To load an alternative biorefinery, simply use the load function
    >>> succinic.load('lab_fed-_batch')
    >>> succinic.simulate_and_print()
    
    ---------- Simulation Results ----------
    MPSP is $1.063/kg
    GWP100a is 0.932 kg CO2-eq./kg
    FEC is -5.360 MJ/kg
    ----------------------------------------

Systems
-------
The `Flowsheet`, `System`, and `TEA`, objects are `flowsheet`, `system`, and `succinic_tea`, respectively.

Naming conventions:
    D = Distillation column

    E = Evaporator
    
    F = Flash tank

    C = Crystallizer

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
    >>> succinic.load('pilot_batch')
    >>> succinic.simulate_and_print()
    >>> succinic.flowsheet
    <Flowsheet: succinic>
    >>> succinic.system
    System: succinic_sys
    Highest convergence error among components in recycle
    stream M505-0 after 1 loops:
    - flow rate   2.02e-03 kmol/hr (0.0018%)
    - temperature 2.21e-07 K (7.2e-08%)
    >>> succinic.succinic_tea.show()
    CellulosicEthanolTEA: succinic_sys
    NPV: 49 USD at 10.0% IRR
    >>> # You can use the flowsheet to access streams, unit operations, and subsystems
    >>> u = succinic.flowsheet.unit
    >>> u.C401.show(N=100)
    SuccinicAcidCrystallizer: C401
    ins...
    [0] s28  from  Pump-F401_P
        phase: 'l', T: 333.21 K, P: 101325 Pa
        flow (kmol/hr): H2O           1.34e+03
                        Ethanol       4.81e-08
                        AceticAcid    4.39e-07
                        Sucrose       0.0354
                        FermMicrobe   0.438
                        PyruvicAcid   1.06
                        SuccinicAcid  65.1
                        H3PO4         1.05
                        Flocculant    8.3
    outs...
    [0] C401_0  to  PressureFilter-S402
        phases: ('l', 's'), T: 273.4 K, P: 101325 Pa
        flow (kmol/hr): (l) H2O           1.34e+03
                            Ethanol       4.81e-08
                            AceticAcid    4.39e-07
                            Sucrose       0.0354
                            FermMicrobe   0.438
                            PyruvicAcid   1.06
                            SuccinicAcid  28.8
                            H3PO4         1.05
                            Flocculant    8.3
                        (s) SuccinicAcid  36.3

    
Analyses
--------
Multiple analysis modules were used to evaluate biorefineries for [1]_ and [3]_, including
full Monte Carlo simulation (succinic.run_uncertainty_analysis)
and titer-yield-productivity analysis for the fermentation performance space (succinic.run_TRY_analysis).
The crystallization simulation model was fit to experimental data (crystallization_curvefit.py).

Note that results used in the manuscript [1]_ were generated using biosteam==2.36.4,
thermosteam==0.33.0, contourplots==0.2.1, numpy==1.23.5, and dependencies (`commit 616b812 <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/616b8122863e966aa2f97e2ef7b0e85544731392/biorefineries/succinic>`_).

Note that results used in the manuscript [3]_ were generated using biosteam==2.46.1,
thermosteam==0.45.0, contourplots==0.2.1, numpy==1.26.4, and dependencies (`this branch <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/succinic_2025_no_FGI/biorefineries/succinic>`_).

To reproduce the results, directly run the script of interest, and results will
be saved as Excel files and figures in analyses/results.


References
----------

.. [1] Tran, V. G.\ :sup:`*`; Mishra, S.\ :sup:`*`; Bhagwat, S. S.\ :sup:`*`; Shafaei, S.; Shen, Y.; Allen, J. L.; Crosly, B. A.; Tan, S.; Fatma, Z.; Rabinowitz, J. D.; Guest, J. S.; Singh, V. S.; Zhao, H. An end-to-end pipeline for succinic acid production at an industrially relevant scale using Issatchenkia orientalis. Nature Communications. 2023, 14 (1), 6152. \ :sup:`*` These authors contributed equally to this work.
    `<https://doi.org/10.1038/s41467-023-41616-9>`_

.. [2] Bhagwat, S. S.; Li, Y.; Cortés-Peña, Y. R.; Brace, E. C.; Martin, T. A.; Zhao, H.; Guest, J. S. Sustainable production of acrylic acid via 3-hydroxypropionic acid from lignocellulosic biomass. ACS Sustainable Chemistry & Engineering. 2021, 9 (49), 16659–16669.
    `<https://doi.org/10.1021/acssuschemeng.1c05441>`_
     
.. [3] Tran, V. G.; Tan, S.; Xu, H; Weilandt, D.R.; Li, X.; Bhagwat, S. S.; Zhu, Z.; Guest, J. S.; Rabinowitz, J.D.; Zhao, H. Decompartmentalization of the yeast mitochondrial metabolism to improve chemical production in Issatchenkia orientalis. Nature Communications. Accepted June 25, 2025.
    `<https://doi.org/10.1021/acssuschemeng.9b07040>`_

