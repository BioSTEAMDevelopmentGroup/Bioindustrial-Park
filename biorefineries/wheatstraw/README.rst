==============================================================
wheatstraw: 2nd Generation Ethanol Production from Wheat Straw
==============================================================

The wheatstraw biorefinery design for the production of ethanol
follows all assumptions from a study by Humbird et. al. executed with Aspen 
plus [1]_. As discussed in our paper [2]_, the wheatstraw biorefinery can be 
divided into eight areas: pretreatment, fermentation, 
ethanol purification, biogas production, wastewater treatment, storage, 
boiler and turbo-generator system, and utilities. During pretreatment, the 
wheatstraw is soaked with sulfuric acid and water, filtered and presteamed 
under high pressures and temperatures to break down the biomass and facilitate 
enzymatic hydrolysis. The liquid fraction is sent to the biogas reactor, and the 
solid fraction is sent to the fermentor, where both saccharification and fermentation
are carried out. In the ethanol purification, the liquid outlet stream from the
fermentor is distilled in several stages at differnt pressures so that the heat needed
in the bolier of a column can be provided by the heat produced in the condenser of the
other column. The bottom fraction is filter-pressed. The solid fraction is sent to the 
boiler turbogenerator to produce steam and electricity, and the liquid fraction is sent
to the biogas production area. The liquid from the biogas reactor is sent to the wastewater
treatment plant and the treated water in recirculated into the process.
The results shown in [2]_ correspond to an implementation of the biorefinery in the Biosteam 
version 2.1.9. The update of the biorefinery to the latest version (v2.20.19) led to a slight
change in the results due to the fact that the boiler turbogenerator is modelled differently 
in the latest version.


References
----------
.. [1] Humbird, D.; Davis, R.; Tao, L.; Kinchin, C.; Hsu, D.; Aden, A.; Schoen, 
    P.; Lukas, J.; Olthof, B.; Worley, M.; Sexton, D.; Dudgeon, D. Process 
    Design and Economics for Biochemical Conversion of Lignocellulosic Biomass 
    to Ethanol: Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn 
    Stover, Technical Report NREL/TP-5100-47764; DOE: NREL, 2011.

.. [2] Sanchis-Sebastiá, M.; Gomis-Fons, J.; Galbe, M.; Wallberg, O. Techno-Economic
       Evaluation of Biorefineries Based on Low-Value Feedstocks Using the BioSTEAM 
       Software: A Case Study for Animal Bedding. Processes 2020, 8, 904. 
       https://doi.org/10.3390/pr8080904