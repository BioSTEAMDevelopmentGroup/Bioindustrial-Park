# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2023-, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from biorefineries.miscanthus import create_lactic_system

sys = create_lactic_system()
tea = sys.TEA

f = sys.flowsheet
lactic_acid = f.stream.lactic_acid

sys.simulate()
get_MPSP = lambda: tea.solve_price(lactic_acid)
get_GWP = lambda: sys.get_net_impact('GWP')/sys.operating_hours/lactic_acid.F_mass

print(f'price: {get_MPSP()}')
print(f'GWP: {get_GWP()}')