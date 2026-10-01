=======================================================================================
Plasma-based nitric acid production biorefinery
=======================================================================================

This biorefinery is developed for Nguyen, et.al. [1]_ for the production of 
plasma_based nitric acid.

Installation
------------
In an environment with Python v3.13.7, do the following:
    (1) Clone this repository (this may require a few minutes' time; git clone https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park)
    (2) Add /Bioindustrial-Park/biorefineries to your Python paths
    (3) Install the following packages in sequence (this may require a few minutes' time):
	    (i) biosteam==2.52.14
	    (ii) thermosteam==0.52.10


Systems
-------
The system includes such units: 


    C101: IsothermalCompressor

    R101: PlasmaReactor
    
    U101 = PowerUnit


Analyses
--------
For full Monte Carlo simulation for TEA/LCA, directly run _model, and results will
be saved as Excel files in results modules.

Note that results used in the manuscript [1]_ were generated using biosteam==2.52.14,
thermosteam==0.52.10, and dependencies (`commit fb23c23 <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/commit/fb23c23328c287188e05c2e699c32cd4fece3a11>`_).


References
----------
.. [1] Nguyen et al., Integrating Plasma-based Nitrogen Fixation with Microbial Fermentation for Bio-based Manufacturing. 
    Nat. Commun. 2026. Submitted 2025.


