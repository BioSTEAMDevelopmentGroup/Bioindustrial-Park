# -*- coding: utf-8 -*-
"""
Created on 2025-05-07 18:26:22

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biosteam import stream_kwargs

LegH={}

# %% In
SeedIn = stream_kwargs('SeedIn', Seed=0.15*1e5/16, units='kg/hr', T=25+273.15,price=0.6*1.3*1e5/72)
CultureIn = stream_kwargs('CultureIn', Culture=1.5*1e5/16, units='kg/hr', T=25+273.15,price=0.6*1.3*1e5/72)

Glucose = stream_kwargs('Glucose', Glucose=1.3*1e5/72, units='kg/hr', T=25+273.15, price=0.6*1.3*1e5/72)
NH3_18wt = stream_kwargs('NH3_18wt', NH3_18wt=1000, units='kg/hr', T=25+273.15,price=0.6*1.3*1e5/72)

WashingSolution = stream_kwargs('WashingSolution', DiaBuffer=1, units='kg/hr', T=25+273.15, price=0.6*1.3*1e5/72)
Elution = stream_kwargs('Elution', IEXBuffer=1 , units='kg/hr', T=25+273.15,price=0.6*1.3*1e5/72)
NFBuffer = stream_kwargs('NFBuffer', NanoBuffer=1, units='kg/hr', T=25+273.15,price=0.6*1.3*1e5/72)

# %% Out
LegH_Ingredients = stream_kwargs('LegH_Ingredients',units='kg/hr', T=25+273.15, price=100)
vent1 = stream_kwargs('vent1')
vent2 = stream_kwargs('vent2')
effluent1 = stream_kwargs('effluent1')
effluent2 = stream_kwargs('effluent2')
effluent3 = stream_kwargs('effluent3')
effluent4 = stream_kwargs('effluent4')
effluent5 = stream_kwargs('effluent5')
effluent6 = stream_kwargs('effluent6')


# SeedIn = stream_kwargs('SeedIn')
# CultureIn = stream_kwargs('CultureIn')
# Glucose = stream_kwargs('Glucose')
# FilteredAir = stream_kwargs('FilteredAir', phase='g')
# NaOH = stream_kwargs('NaOH', price=0.14952852932689323)
# ammonia = stream_kwargs('ammonia', phase='l', P=12 * 101325, price=0.4485966110937889)

# sulfuric_acid = stream_kwargs(
#     'sulfuric_acid',
#     P=5.4*101325,
#     T=294.15,
#     Water=130,
#     H2SO4=1800,
#     units='kg/hr',
#     price=0.0897171175961359
# )

# denaturant = stream_kwargs(
#     'denaturant',
#     Octane=1,
#     price=0.756,
# )

# DAP = stream_kwargs('DAP', price=0.9869213628968231)
# CSL = stream_kwargs('CSL', price=0.0568241480781522)
# vent = stream_kwargs('vent')
# cellulase = stream_kwargs('cellulase', price=0.212)