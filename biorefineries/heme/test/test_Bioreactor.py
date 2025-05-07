# -*- coding: utf-8 -*-
"""
Created on 2025-04-30 11:03:04

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""
# %%
import biosteam as bst
from biorefineries.heme._chemicals import create_chemical_LegH
# from biorefineries.heme._units import Fermentation

# %%
bst.nbtutorial()
bst.settings.set_thermo(create_chemical_LegH(),skip_checks=True)
bst.preferences.N=50

# %%

# Culture_Media = bst.Stream('Culture_Media', Culture_Media=10, units='kg/hr', T=32+273.15)
# Trace_Metal_Solution = bst.Stream('Trace_Metal_Solution', Trace_Metal_Solution=10, units='kg/hr', T=32+273.15)
# Antibiotics = bst.Stream('Antibiotics', Antibiotics=1, units='kg/hr', T=32+273.15)
# Feed1 = bst.Stream('Feed1', Feed1=1000, units='kg/hr', T=32+273.15)
# Feed2 = bst.Stream('Feed2', Feed2=1000, units='kg/hr', T=32+273.15)

Seed = bst.Stream('Seed', Seed=1500, units='kg/hr', T=32+273.15)
Culture = bst.Stream('Culture', Culture=1500, units='kg/hr', T=32+273.15)
FilteredAir = bst.Stream('FilteredAir', air=100,units='kg/hr', T=32+273.15)
Glucose = bst.Stream('Glucose', Glucose=1300, units='kg/hr', T=32+273.15)
NH3 = bst.Stream('NH3', NH3=500, units='kg/hr', T=32+273.15)

# Seed = bst.Stream('Seed', Seed, units='kg/hr', T=32+273.15)
# Culture = bst.Stream('Culture', Culture, units='kg/hr', T=32+273.15)
# FilteredAir = bst.Stream('FilteredAir', air ,units='kg/hr', T=32+273.15)
# Glucose = bst.Stream('Glucose', Glucose, units='kg/hr', T=32+273.15)
# NH3 = bst.Stream('NH3', NH3=4, units='kg/hr', T=32+273.15)


# Set all stream properties as water but retain composition
# for stream in [Culture_Media, Trace_Metal_Solution, Antibiotics, Feed1, Feed2]:
#     stream.copy_thermal_condition(bst.Stream('temp', Water=1, T=stream.T, P=stream.P))

# %%
theta_O2 = 0.5 # Dissolved oxygen concentration [% saturation]
agitation_power = 0.985 # [kW / m3]
design = 'Stirred tank' # Reactor type
method = "Riet" # Name of method

T_operation = 273.15 + 32 # [K]
Q_O2_consumption = -110 * 4184 # [kJ/kmol]
dT_hx_loop = 8 # [degC]

cooler_pressure_drop = 20684 # [Pa]
compressor_isentropic_efficiency = 0.85

V_max = 500 # [m3]
titer = 7.27 # [g / L]
productivity = 0.1 # [g / L / h]
LegH_yield = 0.18 # [by wt]
Y_b = 0.43 # [by wt]



# %%
Production = bst.SRxn([
        bst.Rxn('8 Glucose + 4 NH3 + 1 FeSO4 + 10.5 O2 -> Heme_b + 1 H2SO4 + 37 H2O + 14 CO2',reactant = 'Glucose',X=0.1,check_atomic_balance=True),
        bst.Rxn('Glucose + NH3 + H2SO4+ O2 -> Protein + H2O',reactant = 'Glucose', X= 0.1,correct_atomic_balance=True),
        bst.Rxn('Heme_b + Protein -> Leghemoglobin + H2O', reactant = 'Heme_b', X=0.95, correct_atomic_balance=True),
        ])

bst.settings.chemicals.set_alias('Yeast', 'cellmass')

Production[2].product_yield('Leghemoglobin', basis='wt', product_yield=LegH_yield)
growth = bst.Reaction(
    'Glucose -> H2O + CO2 + Yeast', 'Glucose', 1,
    correct_atomic_balance=True
)
growth.product_yield('Yeast', basis='wt', product_yield=Y_b)
respiration = bst.Reaction(
    'Glucose + O2 -> CO2 + H2O', 'Glucose', 1. - growth.X,
    correct_atomic_balance=True
)
RXN = bst.ReactionSystem(
    Production,
    bst.ParallelReaction([growth, respiration])
)
RXN.show()

# %%
effluent = bst.Stream('effluent')
AB1 = bst.AeratedBioreactor('AB1',
    ins=[Culture, Glucose, NH3, bst.Stream('air', phase='g')],
    outs=('vent','product'),
    design='Stirred tank', method=method,
    V_max=V_max, Q_O2_consumption=Q_O2_consumption,
    dT_hx_loop=dT_hx_loop, T=T_operation,
    batch=True, reactions=RXN,
    kW_per_m3=agitation_power,
    tau=titer/productivity,
    cooler_pressure_drop=cooler_pressure_drop,
    compressor_isentropic_efficiency=compressor_isentropic_efficiency,
    optimize_power=True,
)
AB1.target_titer = titer # g / L
AB1.target_productivity = productivity # g / L / h
AB1.target_yield = LegH_yield  # wt %

@AB1.add_specification(run=True)
def update_reaction_time_and_yield():
    AB1.tau = AB1.target_titer / AB1.target_productivity
    Production[2].product_yield('Leghemoglobin', basis='wt', product_yield=AB1.target_yield)


# %%
AB1.simulate()
AB1.show()
AB1.diagram(format='html')
# %%
AB1.design_results
# %%
AB1.purchase_costs
# %%
# AB1.heat_utilities

# %%
# %%

# R1 = bst.AeratedBioreactor(
#     'R1', ins=[Culture, Glucose, NH3, FilteredAir], outs=('vent', 'product'), tau=12, V_max=500,
#     reactions=Production
# )
# # RuntimeError: Liquid molar volume method 'NEGLECT_P' is not valid at T=305.15 K and P=101325.0 Pa for component with CASRN '7733-02-0'

# # %%
# R1.simulate()
# R1.show()
# # %%
# R1.design_results

# R1.purchase_costs

# R1.heat_utilities

# # %%
# R1.diagram(format='html')
# # %%
