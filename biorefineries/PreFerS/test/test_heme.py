# -*- coding: utf-8 -*-
"""
Created on 2025-05-06 18:26:54

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

# %%
import biosteam as bst
import thermosteam as tmo
import numpy as np
from biorefineries.PreFerS import _chemicals, _units,_streams

# %% Settings
bst.nbtutorial()
_chemicals.__all__
_units.__all__
_streams
bst.settings.set_thermo(_chemicals.create_chemicals_LegH(), skip_checks=True)
bst.preferences.N=50

# %% Streams

theta_O2 = 0.5 # Dissolved oxygen concentration [% saturation]
agitation_power = 0.985 # [kW / m3]
design = 'Stirred tank' # Reactor type
method = "Riet" # Name of method

T_operation = 273.15 + 32 # [K]
Q_O2_consumption = -110 * 4184 # [kJ/kmol]
dT_hx_loop = 8 # [degC]

cooler_pressure_drop = 20684 # [Pa]
compressor_isentropic_efficiency = 0.85

V_max = 500 # [m3] 5 L * 1e5
titer = 7.27 # [g / L]
productivity = 7.27 / 72 # [g / L / h]
LegH_yield = 17.27 * 5 * 1e5/1000/1300 # [by wt]
Y_b = 0.43 # [by wt]


# %% Reactions
conversion = productivity*V_max/(1.3*1e5/72)

fermentation_reaction = bst.PRxn([
        bst.Rxn('8 Glucose + 4 NH3 + 1 FeSO4 + 10.5 O2 -> Heme_b + 1 H2SO4 + 37 H2O + 14 CO2',reactant = 'Glucose',X=conversion*0.05,check_atomic_balance=True),
        bst.Rxn('Glucose + NH3 + H2SO4 + O2 -> Globin + H2O',reactant = 'Glucose', X= conversion*0.05,correct_atomic_balance=True),
        bst.Rxn('Glucose + FeSO4 + NH3 + H2SO4 + O2 -> Leghemoglobin + H2O', reactant = 'Glucose', X=conversion,  correct_atomic_balance=True),
        ])
fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=LegH_yield)

cell_growth_reaction = bst.Rxn(
    'Glucose -> H2O + CO2 + Yeast', 'Glucose', 1,
    correct_atomic_balance=True
)
cell_growth_reaction.product_yield('Yeast', basis='wt', product_yield=Y_b)

respiration_reaction = bst.Rxn(
    'Glucose + O2 -> CO2 + H2O', 'Glucose', 1. - cell_growth_reaction.X,
    correct_atomic_balance=True
)

bst.settings.chemicals.set_alias('Yeast', 'cellmass')

RXN = bst.ReactionSystem(
    fermentation_reaction,
    bst.PRxn([cell_growth_reaction, respiration_reaction])
)
RXN.show()

# %% Process
Mix1Out = bst.Stream('Mix1Out')
SeedIn = bst.Stream('SeedIn', Seed=0.15*1e5/16, units='kg/hr', T=25+273.15)
CultureIn = bst.Stream('CultureIn', Culture=1.5*1e5/16, units='kg/hr', T=25+273.15)

MX1 = bst.Mixer(
    ins=[SeedIn, CultureIn],
    outs=[Mix1Out],
)
# %%

SeedOut = bst.Stream('SeedOut')
ST1 = _units.SeedTrain(
    ins=[Mix1Out],
    outs=['vent', SeedOut],
    reactions=None,
    saccharification=None,
    T=32+273.15,
)
# %%
Glucose = bst.Stream('Glucose', Glucose=1.3*1e5/72, units='kg/hr', T=25+273.15)
_18wtNH3 = bst.Stream('_18wtNH3', _18wtNH3=2000, units='kg/hr', T=25+273.15)
effluent = bst.Stream('effluent')

AB1 = _units.AeratedFermentation(
    ins=[SeedOut, Glucose, _18wtNH3, bst.Stream('FilteredAir', phase='g', P = 2 * 101325)],
    outs=['vent', 'product'],
    fermentation_reaction=fermentation_reaction,
    cell_growth_reaction=cell_growth_reaction,
    respiration_reaction=respiration_reaction,
    design='Stirred tank', method=method,
    V_max=V_max, Q_O2_consumption=Q_O2_consumption,
    dT_hx_loop=dT_hx_loop, T=T_operation,
    batch=True, reactions=RXN,
    kW_per_m3=agitation_power,
    tau=titer/productivity,
    cooler_pressure_drop=cooler_pressure_drop,
    compressor_isentropic_efficiency=compressor_isentropic_efficiency,
    P=1 * 101325, #optimize_power=True,
)

AB1.target_titer = titer # g / L
AB1.target_productivity = productivity # g / L / h
AB1.target_yield = LegH_yield  # wt %

@AB1.add_specification(run=True)
def update_reaction_time_and_yield():
    AB1.tau = AB1.target_titer / AB1.target_productivity
    fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=AB1.target_yield)

# %%
LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
LegH_sys.simulate()
LegH_sys.diagram(format='html', kind='cluster', number=True)
LegH_sys.show()

# %%
