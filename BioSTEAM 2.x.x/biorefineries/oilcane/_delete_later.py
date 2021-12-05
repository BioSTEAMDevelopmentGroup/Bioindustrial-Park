# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 22:30:49 2021

@author: yrc2
"""

from biorefineries import oilcane as oc

cfdiff = 3.3111 * (
    oc.GWP_characterization_factors['biodiesel']
    - oc.GWP_characterization_factors['biodiesel displacement'] 
)
configurations = ['O1', 'O2', 'O1*', 'O2*']
for name in configurations:
    df = oc.get_monte_carlo(name)
    file_name = oc.monte_carlo_file(name)
    df[oc.GWP_ethanol_displacement.index] += (
        df[oc.biodiesel_production.index] * cfdiff
        / df[oc.ethanol_production.index]
    )
    df.to_excel(file_name)
    