# -*- coding: utf-8 -*-
"""
Created on Sat June 28 10:00:00 2025

Microalgae biorefinery to produce medium chain fatty acids 
by anaerobic fermentation without external electron donor addition- feed streams

References
----------
[1] BioSTEAM Documentation: 
    https://biosteam.readthedocs.io/en/latest/tutorial/Creating_a_System.html
[2] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
[3] 3-Hydroxypropionic acid biorefineries project:
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/HP
[4] Succinic projest
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/succinic

@author: Xingdong Shi
@version: 0.0.1
"""

from biosteam import stream_kwargs

microalgae_feed = stream_kwargs('microalgae', 
                                Microalgae=5000, 
                                units='kg/hr',
                                price = 0.15  # this is a baseline value for TEA (cost for cultivation, harvesting and trasporting)
)
