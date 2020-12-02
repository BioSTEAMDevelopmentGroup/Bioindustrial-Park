==============================================================
cornstover: 2nd Generation Ethanol Production from Corn Stover
==============================================================

.. figure:: ./images/cornstover_areas.png

The corn stover biorefinery design for the production of cellulosic ethanol 
follows all assumptions from a study by Humbird et. al. executed with Aspen 
plus [1]_. As discussed in the original BioSTEAM manuscript [2]_, the 
corn stover biorefinery can be divided into eight areas: feedstock handling, 
pretreatment and conditioning, enzyme production, recovery, wastewater 
treatment, storage, boiler and turbo-generator system, and utilities. In 
contrast to the original report by Humbird et al., the on-site enzyme 
production model is not included. Instead, the cellulase mixture is bought at 
the same concentration of 50 g·L–1 and at a price of 0.424 $·kg–1, which is 
the back-calculated cost of on-site enzyme production as reported in the 
benchmark study. During pretreatment, the corn stover is presteamed and 
loaded into pretreatment reactors with dilute sulfuric acid under high 
pressures and temperatures to break down the biomass and facilitate enzymatic 
hydrolysis (Area 200). The hydrolysate is then saccharified and fermented in 
the fermentation area (Area 300). In the recovery area, the beer is distilled, 
and the lignin-rich solid fraction is filter-pressed (Area 400). The stillage 
is sent to wastewater treatment to produce biogas (Area 500), while the lignin 
is sent to the boiler-turbo-generator section (Area 700) to produce steam and 
electricity.

Getting Started
---------------

To load the biorefinery, simply import it. All data and variables
are lazy loaded by the module:

.. code-block:: python

    >>> from biorefineries import cornstover as cs
    >>> # This is optional; it forces the biorefinery to load
    >>> # Otherwise, first time accessing will take a bit to load.
    >>> cs.load()
    >>> cs.chemicals # All chemicals used by the biorefinery.
    CompiledChemicals([Water, Ethanol, AceticAcid, Furfural, Glycerol, H2SO4, LacticAcid, SuccinicAcid, P4O10, HNO3, Denaturant, DAP, AmmoniumAcetate, AmmoniumSulfate, NaNO3, Oil, HMF, N2, NH3, O2, CH4, H2S, SO2, CO2, NO2, NO, CO, Glucose, Xylose, Sucrose, CaSO4, Mannose, Galactose, Arabinose, CellulaseNutrients, Extract, Acetate, Tar, CaO, Ash, NaOH, Lignin, SolubleLignin, GlucoseOligomer, GalactoseOligomer, MannoseOligomer, XyloseOligomer, ArabinoseOligomer, Z_mobilis, T_reesei, Biomass, Cellulose, Protein, Enzyme, Glucan, Xylan, Xylitol, Cellobiose, CSL, DenaturedEnzyme, Arabinan, Mannan, Galactan, WWTsludge, Cellulase])
    >>> cs.cornstover_sys.show() # The complete biorefinery system
    System: cornstover_sys
     path: (pretreatment_sys, fermentation_sys,
            purification, S401, M601, WWTC,
            R601, aerobic_digestion_sys, S604)
     facilities: (M501, CWP, BT, CT, PWC, ADP, CIP_package,
                  S301, S302, DAP_storage, CSL_storage,
                  FT, blowdown_mixer)
    >>> cs.cornstover_tea.show() # The TEA object
    CombinedTEA: cornstover_sys, Area700
     NPV: -193,156 USD at 10.0% IRR
    >>> cs.flowsheet # The complete flowsheet
    <Flowsheet: cornstover>
    >>> cs.R301.show() # Any unit operations and streams can be accessed through the module
    SaccharificationAndCoFermentation: R301
    ins...
    [0] d329  from  Mixer-M302
        phase: 'l', T: 319.61 K, P: 101325 Pa
        flow (kmol/hr): Water              1.92e+04
                        Ethanol            49.9
                        AceticAcid         21.5
                        Furfural           0.327
                        Glycerol           0.21
                        H2SO4              20.9
                        SuccinicAcid       0.409
                        ...
    [1] CSL2  from  ReversedSplitter-S302
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CSL  22.6
    [2] DAP2  from  ReversedSplitter-S301
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): DAP  0.779
    outs...
    [0] d330  to  Mixer-M304
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water         17.5
                        Ethanol       4.53
                        AceticAcid    0.0394
                        Furfural      0.00191
                        Glycerol      1.34e-07
                        LacticAcid    1.14e-05
                        SuccinicAcid  0.000737
                        ...
    [1] d331  to  Mixer-M401
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water              1.72e+04
                        Ethanol            460
                        AceticAcid         19.3
                        Furfural           0.294
                        Glycerol           1.84
                        H2SO4              18.8
                        LacticAcid         17.5
                        ...
    [2] d332  to  Mixer-M303
        phase: 'l', T: 321.15 K, P: 101325 Pa
        flow (kmol/hr): Water              1.91e+03
                        Ethanol            4.99
                        AceticAcid         2.15
                        Furfural           0.0327
                        Glycerol           0.021
                        H2SO4              2.09
                        SuccinicAcid       0.0409
                        ...


References
----------
.. [1] Humbird, D.; Davis, R.; Tao, L.; Kinchin, C.; Hsu, D.; Aden, A.; Schoen, 
    P.; Lukas, J.; Olthof, B.; Worley, M.; Sexton, D.; Dudgeon, D. Process 
    Design and Economics for Biochemical Conversion of Lignocellulosic Biomass 
    to Ethanol: Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn 
    Stover, Technical Report NREL/TP-5100-47764; DOE: NREL, 2011.

.. [2] Cortes-Peña, Y.; Kumar, D.; Singh, V.; Guest, J. S.
    BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and 
    Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020. https://doi.org/10.1021/acssuschemeng.9b07040.


