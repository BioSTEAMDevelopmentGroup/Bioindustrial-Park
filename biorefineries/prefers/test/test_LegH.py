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
from biorefineries.prefers import _chemicals, _units

# %% Settings
# bst.nbtutorial()
# _chemicals.__all__
# _units.__all__
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
productivity = titer / 72 # [g / L / h]
LegH_yield = titer * 5 / 1300 # [by wt]
Y_b = 0.43 # [by wt]


# %% Reactions

fermentation_reaction = bst.PRxn([
        #           Reaction        Reactnat            Conversion           Check                  
        bst.Rxn('8 Glucose + 4 NH3 + 1 FeSO4 + 10.5 O2 -> Heme_b + 1 H2SO4 + 37 H2O + 14 CO2',
                                    reactant = 'Glucose',X=LegH_yield*0.05,check_atomic_balance=True),
        bst.Rxn('Glucose + (NH4)2SO4 + O2 -> Globin + NH3 + H2O',
                                    reactant = 'Glucose', X= LegH_yield*0.05,correct_atomic_balance=True),
        bst.Rxn('Glucose + FeSO4 + (NH4)2SO4 + O2 -> Leghemoglobin + NH3 + H2O', 
                                    reactant = 'Glucose', X=LegH_yield,  correct_atomic_balance=True),
        ])
fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=LegH_yield)

neutralization_reaction = bst.Rxn(
    'H2SO4 + NH3 -> (NH4)2SO4', reactant = 'H2SO4', X=1,
    correct_atomic_balance=True
)

cell_growth_reaction = bst.Rxn(
    'Glucose -> H2O + CO2 + Pichia_pastoris', 'Glucose', X=Y_b,
    correct_atomic_balance=True
)
cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=Y_b)

respiration_reaction = bst.Rxn(
    'Glucose + O2 -> CO2 + H2O', 'Glucose', 1. - cell_growth_reaction.X - fermentation_reaction[2].X,
    correct_atomic_balance=True
)

bst.settings.chemicals.set_alias('Pichia_pastoris', 'cellmass')

RXN = bst.ReactionSystem(
    fermentation_reaction,
    bst.PRxn([cell_growth_reaction, respiration_reaction])
)
RXN.show()

# %% Upstream Process
Mix1Out = bst.Stream('Mix1Out')
SeedIn = bst.Stream('SeedIn', Seed=0.15*1e5/16, units='kg/hr', T=25+273.15)
CultureIn = bst.Stream('CultureIn', Culture=1.5*1e5/16, units='kg/hr', T=25+273.15)

MX1 = bst.Mixer(
    ins=[SeedIn, CultureIn],
    outs=[Mix1Out],
)

SeedOut = bst.Stream('SeedOut')
vent1 = bst.Stream('vent1')
ST1 = _units.SeedTrain('ST1',
    ins=[Mix1Out],
    outs=[vent1, SeedOut],
    reactions=bst.PRxn([cell_growth_reaction, respiration_reaction]),
    saccharification=None,
    T=32+273.15,
)

Y_b = 0.25 # [by wt]
RXN.show()

Glucose = bst.Stream('Glucose', Glucose=1.3*1e5/72, units='kg/hr', T=25+273.15)
NH3_18wt = bst.Stream('NH3_18wt', NH3_18wt=1000, units='kg/hr', T=25+273.15)
vent2 = bst.Stream('vent2')
AB1Out = bst.Stream('AB1Out')

AB1 = _units.AeratedFermentation('AB1',
    ins=[SeedOut, Glucose, NH3_18wt, bst.Stream('FilteredAir', phase='g', P = 2 * 101325)],
    outs=[vent2, AB1Out],
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
# LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
# LegH_sys.simulate()
# %% Downstream process
CD1Out = bst.Stream('CD1Out')
CD1 = _units.CellDisruption('CD1',
    ins=AB1Out,
    outs=CD1Out,
)
# LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
# LegH_sys.simulate()
# %%
PC1Out = bst.Stream('PC1Out')
effluent1 = bst.Stream('effluent1')
PC1 = _units.ProteinCentrifuge('PC1',
    ins = CD1Out,
    outs = (effluent1, PC1Out),
    moisture_content = 0.5,
    split = (1, 0.99, 1, 1),
    order = ('Glucose','cellmass', 'LeghemoglobinIntre','GlobinIntre'),
    solids = ('Glucose','cellmass', 'LeghemoglobinIntre','GlobinIntre'),
    centrifuge_type = 'reciprocating_pusher'
)

# ???
EV1Out = bst.Stream('EV1Out')
effluent2 = bst.Stream('effluent2')
EV1 = _units.Evaporator('EV1',
    ins = PC1Out,
    outs = (EV1Out, effluent2),
    P = (101325, 73581, 50892, 32777, 20000),
    V = 0.1,
    V_definition = 'First-effect',
) # ???
# LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
# LegH_sys.simulate()
# This refers to how much diafiltration buffer you use relative to the volume of your protein solution (the retentate) during constant volume diafiltration.

# Recommended: 3 to 5 diavolumes.
# 1 diavolume (DV) = Volume of diafiltration buffer added is equal to the volume of the protein solution being processed (the retentate volume, which is kept constant).
# 3 DV will result in approximately 95% exchange of the original buffer/small solutes.
# 5 DV will result in approximately 99.3% exchange.
# For very thorough exchange (>99.9%): You might consider 6 to 7 diavolumes.
# 10 DV is not recommended as it will result in a very low concentration of the protein solution and a very high volume of diafiltration buffer.

# %%
DF1Out = bst.Stream('DF1Out')
WashingSolution1 = bst.Stream('WashingSolution1', units='kg/hr', T=25+273.15)
# WashingSolution1 = bst.Stream('WashingSolution1', DiaBuffer=EV1Out.imass['H2O']*4, units='kg/hr', T=25+273.15)
effluent3 = bst.Stream('effluent3')
DF1 = _units.Diafiltration( 'DF1',
    ins = (EV1Out, WashingSolution1),
    outs = (DF1Out, effluent3),
    TargetProduct_ID = 'Leghemoglobin',
    Salt_ID = _chemicals.chemical_groups['Salts'],
    OtherLargeMolecules_ID = _chemicals.chemical_groups['OtherLargeMolecules'],
    DefaultSolutes_ID = _chemicals.chemical_groups['DefaultSolutes'],
)

@DF1.add_specification(run=True)
def update_WashingSolution1():
    WashingSolution1 = DF1.ins[1]
    EV1Out = EV1.outs[0]
    WashingSolution1.imass['DiaBuffer'] = EV1Out.imass['H2O']*4

LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
LegH_sys.simulate()

# %%
IEX1Out = bst.Stream('IEX1Out')
Elution = bst.Stream('Elution', units='kg/hr', T=25+273.15)
# Elution = bst.Stream('Elution', IEXBuffer=DF1Out.imass['H2O']/2 ,units='kg/hr', T=25+273.15)
effluent4 = bst.Stream('effluent4')

IEX1 = _units.IonExchange( 'IEX1',
    ins = (DF1Out, Elution),
    outs = (IEX1Out, effluent4),
    TargetProduct_ID = 'Leghemoglobin',
    BoundImpurity_ID=_chemicals.chemical_groups['BoundImpurities'],
    ElutionBuffer_Defining_Component_ID =_chemicals.chemical_groups['ElutionBuffer'],
)

@IEX1.add_specification(run=True)
def update_Elution():
    Elution = IEX1.ins[1]
    DF1Out = DF1.outs[0]
    Elution.imass['IEXBuffer'] = DF1Out.imass['H2O']/2

# LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
# LegH_sys.simulate()
# %%
# NF1Out = bst.Stream('NF1Out')
# NFBuffer = bst.Stream('NFBuffer', NanoBuffer=IEX1Out.imass['H2O']*5, units='kg/hr', T=25+273.15)
# effluent5 = bst.Stream('effluent5')
# NF1 = _units.NanofiltrationDF(
#     ins = (IEX1Out, NFBuffer),
#     outs = (NF1Out, effluent5),
#     TargetProduct_ID = 'Leghemoglobin',
#     Salt_ID = _chemicals.chemical_groups['Salts'],
#     OtherSmallSolutes_ID = _chemicals.chemical_groups['DefaultSolutes'],
# )
# LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
# LegH_sys.simulate()
# %%
NF1Out = bst.Stream('NF1Out')
NFBuffer = bst.Stream('NFBuffer', units='kg/hr', T=25+273.15)
# NFBuffer = bst.Stream('NFBuffer', 
#                     NanoBuffer=1.1*0.05*1000*
#                     IEX1Out.imass['Leghemoglobin']/(0.25*bst.Chemical('TrehaloseDH',search_ID='6138-23-4', phase='l', default=True).MW),
#                     units='kg/hr', T=25+273.15)
effluent5 = bst.Stream('effluent5')
NF1 = _units.Diafiltration('NF1',
    ins = (IEX1Out, NFBuffer),
    outs = (NF1Out, effluent5),
    TargetProduct_ID = 'Leghemoglobin',
    membrane_cost_USD_per_m2=10000, # Nanomembrane cost
    Salt_ID = _chemicals.chemical_groups['Salts'],
    OtherLargeMolecules_ID = _chemicals.chemical_groups['OtherLargeMolecules'],
    DefaultSolutes_ID = _chemicals.chemical_groups['DefaultSolutes'],
    TargetProduct_Retention=0.995, Salt_Retention=0.1,
    OtherLargeMolecules_Retention=0.99, DefaultSolutes_Retention=0.15,
    FeedWater_Recovery_to_Permeate=0.2,
    TMP_bar= 5
)

@NF1.add_specification(run=True)
def update_NFBuffer():
    NFBuffer = NF1.ins[1]
    IEX1Out = IEX1.outs[0]
    NFBuffer.imass['NanoBuffer']=1.1*0.05*1000*IEX1Out.imass['Leghemoglobin']/(0.25*bst.Chemical('TrehaloseDH',search_ID='6138-23-4', phase='l', default=True).MW),

# LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
# LegH_sys.simulate()
# %%
SD1Out = bst.Stream('SD1Out')
effluent6 = bst.Stream('effluent6')
SD1 = bst.SprayDryer(
    ins=NF1Out,
    outs=(effluent6,SD1Out),
    moisture_content=0.05,  # 5% moisture content in the final product
)

# %%
LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
LegH_sys.simulate()
LegH_sys.diagram(kind=0,format='html',label=False)
LegH_sys.show()

# %%
