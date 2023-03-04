# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2023-, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from biorefineries.cornstover import ethanol_density_kggal
from biorefineries.miscanthus import create_ethanol_system

sys = create_ethanol_system()
sys.simulate()
tea = sys.TEA

f = sys.flowsheet
ethanol = f.stream.ethanol

get_MESP = lambda: tea.solve_price(ethanol)*ethanol_density_kggal # from $/kg to $/gallon
get_GWP = lambda: sys.get_net_impact('GWP')/sys.operating_hours/ethanol.F_mass*ethanol_density_kggal

print(f'price: {get_MESP()}')
print(f'GWP: {get_GWP()}')