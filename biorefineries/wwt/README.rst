=========================================================================
wwt: High-Rate Wastewater Treatment Process for Energy and Water Recovery
=========================================================================

Description
-----------

This repository contains the modules for the design of a high-rate wastewater treatment (WWT) and energy/water recovery process as described in Li et al. [1]_ This module has been added as a whole to BioSTEAM's existing wastewater treatment configurations (the "high-rate") one.

For the paper, this process was implemented for seven different biorefineries:
	
- CN (corn-to-ethanol) as in the `corn <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/corn>`_ module, implementation and comparison codes in `corn.py <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/wwt/corn.py>`_
- SC-1G (sugarcane juice-to-ethanol) as in the `oilcane <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/oilcane>`_ module (S1 configuration), implementation and comparison codes in `sugarcane1g.py <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/wwt/sugarcane1g.py>`_
- OC-1G (oilcane juice-to-ethanol and oil-to-biodiesel) as in the `oilcane <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/oilcane>`_ module (O1 configuration), implementation and comparison codes in `oilcane1g.py <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/wwt/oilcane1g.py>`_
- CS (corn stover-to-ethanol) as in the `cornstover <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/cornstover>`_ module, implementation and comparison codes in `cornstover.py <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/wwt/cornstover.py>`_
- SC-2G (sugarcane juice and bagasse-to-ethanol) as in the `oilcane <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/oilcane>`_ module (S2 configuration), implementation and comparison codes in `sugarcane2g.py <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/wwt/sugarcane2g.py>`_
- OC-2G (oilcane juice and bagasse-to-ethanol and oil-to-biodiesel) as in the `oilcane <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/oilcane>`_ module (O2 configuration), implementation and comparison codes in `oilcane2g.py <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/wwt/oilcane2g.py>`_
- LA (corn stover-to-lactic acid) as in the `lactic <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/lactic>`_ module, implementation and comparison codes in `lactic.py <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/wwt/lactic.py>`_

Codes for a quick comparison of key biorefinery metrics (baseline) can be found in `comparison.py <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/wwt/comparison.py>`_ (together with the environment used to generate the results; codes for full comparison in the respective modules as listed above). Archived results generated from this module can be found in the `biorefinery-WWT-archived-results repository <https://github.com/yalinli2/biorefinery-WWT-archived-results>`_


References
----------
.. [1] Li et al., Design of a High-Rate Wastewater Treatment Process for Energy and Water Recovery at Biorefineries. ACS Sustainable Chem. Eng. 2023, 11 (9), 3861–3872. https://doi.org/10.1021/acssuschemeng.2c07139.