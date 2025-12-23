# -*- coding: utf-8 -*-
"""
Created on 2025-12-04 13:25:17

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

#from biorefineries.animal_bedding import _system
from ast import Yield
from math import tau
from biorefineries.prefers._process_settings import set_GWPCF, GWP_CFs,set_GWPCF_Multi,load_process_settings
from biorefineries.prefers.systems import LegH
from biorefineries.prefers.systems.LegH import _streams as s
import biosteam as bst
from pint import set_application_registry
from regex import P
from thermosteam import Stream
from biosteam import F
import thermosteam as tmo
import numpy as np
from biorefineries.prefers import _chemicals as c, _units as u
from biorefineries.prefers._process_settings import price  # ADD THIS IMPORT

# %% Settings
bst.settings.set_thermo(c.create_chemicals_Hemodextrin(), skip_checks=True)
bst.preferences.classic_mode()

# %%
__all__ = (
    'create_Heme_system',
)
# %%
@bst.SystemFactory(
    ID='Heme_sys',
    ins=[s.TraceMetalSolution, s.VitaminCSolution, s.YPD, s.Supplemented, s.Feed1st, s.Feed2nd, s.NH3_25wt,s.DfUltraBuffer, s.IXEquilibriumBuffer, s.IXElutionBuffer
        ,s.IXRegenerationSolution, s.DfNanoBuffer],
    outs=[s.Hemodextrin, s.vent1, s.vent2, s.effluent1, s.effluent2, s.effluent3
    ],
    fthermo=c.create_chemicals_Hemodextrin, # Pass the function itself
)
def create_Heme_system(
        ins, outs,
        use_area_convention=False,
        # reactions_passed=None, # Placeholder if reactions were to be passed externally
        # One could add more configurable parameters here:
        # V_max_fermenter=500, target_titer=7.27, etc.
    ):
    """
    Creates the Hemedextrin (Protein Free) production system.
    This system is based on the process flow and parameters from test_LegH.py.
    """
    bst.preferences.N=50
    
    # Unpack input streams
    SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, IXEquilibriumBuffer, IXElutionBuffer, IXRegenerationSolution, DfNanoBuffer = ins

    # Unpack output streams
    (Heme_3, vent1, vent2, effluent1, effluent2, effluent3) = outs    
    
    ## LCA Settings ##
    
    load_process_settings()  # Load process settings to update prices and CFs
    """
    Fermentation Parameter
    """
    theta_O2 = 0.45 # Dissolved oxygen concentration [% saturation] >= 40%
    agitation_power = 0.985 # [kW / m3]
    design = 'Stirred tank' # Reactor type
    method = "Riet" # Name of method
    T_operation = 273.15 + 30 # [K]
    Q_O2_consumption = -110 * 4184 # [kJ/kmol]
    dT_hx_loop = 8 # [degC]
    cooler_pressure_drop = 20684 # [Pa]
    compressor_isentropic_efficiency = 0.85
    V_max = 500 # [m3] #here cause pressure vessel design problem

    # OD600 180
    # 1 OD600 = 0.35 g/L dry cell weight
    # first 10 hour, 10 g/L/h
    Glucose_Utilization = 0.5614 # [by wt] based on glucose fed
    # 72 hours growth + 90 hours fermentation
    # 1st    
    tau1= 72 # [h] growth time
    P_Heme_1st = 1e-6 # [g / L / h]
    P_pp_1st = 1e-6 # [g / L / h]
    P_ComsuptionG1 = 500 * (0.5/3) / tau1 # [g / L / h] glucose consumption rate
    P_cell_growth_rate1 = 150 * 0.35 / tau1 # [g / L / h] cell growth rate
    Y_Heme1 = P_Heme_1st / P_ComsuptionG1 # product yield on glucose
    Y_pp1 = P_pp_1st / P_ComsuptionG1 # proph
    Y_b1 = P_cell_growth_rate1 / P_ComsuptionG1 # biomass yield on glucose
    Y_r1 = 1 - Y_Heme1 - Y_pp1 - Y_b1 # respiration yield on glucose
    
    # 2nd
    tau2= 90 # [h] fermentation time
    P_Heme_2nd = 4.2 # [g / L / h]
    P_pp_2nd = 0.7 # [g / L / h]
    P_ComsuptionG2 = 800 * (0.5/3) / tau2 # [g / L / h] glucose consumption rate
    P_cell_growth_rate2 = 50 * 0.35 / tau2 # [g / L / h] cell growth rate
    Y_Heme2 = P_Heme_2nd / P_ComsuptionG2 # product yield on glucose
    Y_pp2 = P_pp_2nd / P_ComsuptionG2 # proph
    Y_b2 = P_cell_growth_rate2 / P_ComsuptionG2 # biomass yield on glucose
    Y_r2 = 1 - Y_Heme2 - Y_pp2 - Y_b2 # respiration yield on glucose
    
    # Overall
    tau = tau1 + tau2 # [h] total time
    P_Heme = (P_Heme_1st*tau1 + P_Heme_2nd*tau2) / tau # [g / L / h]
    P_pp = (P_pp_1st*tau1 + P_pp_2nd*tau2) / tau # [g / L / h]
    P_ComsuptionG = (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2) / tau # [g / L / h]
    P_cell_growth_rate = (P_cell_growth_rate1*tau1 + P_cell_growth_rate2*tau2) / tau # [g / L / h]
    Y_Heme = (P_Heme_1st*tau1 + P_Heme_2nd*tau2) / (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2)
    Y_pp = (P_pp_1st*tau1 + P_pp_2nd*tau2) / (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2)
    Y_b = (P_cell_growth_rate1*tau1 + P_cell_growth_rate2*tau2) / (P_ComsuptionG1*tau1 + P_ComsuptionG2*tau2)
    Y_r = 1 - Y_Heme - Y_pp - Y_b # respiration yield on glucose    


    """
    Reactions
    """
    fermentation_reaction = bst.PRxn([
        #           Reaction        Reactnat            Conversion           Check                  "
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                                    reactant = 'Glucose',X=Y_Heme*0.09,check_atomic_balance=True),
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b_Intra + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                                    reactant = 'Glucose', X= Y_Heme*0.09,check_atomic_balance=True),
        ])

    
    
    
    fermentation_reaction = bst.PRxn([
        #           Reaction        Reactnat            Conversion           Check                  "
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                                    reactant = 'Glucose',X=yield_HemeB*0.09,check_atomic_balance=True),
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b_Intra + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                                    reactant = 'Glucose', X= yield_HemeB*0.09,check_atomic_balance=True),
        ])
    fermentation_reaction[0].product_yield('Heme_b', basis='wt', product_yield=yield_HemeB*0.09)
    
    cell_growth_reaction = bst.Rxn(
        'Glucose + 0.8364 NH3 + 0.0108 (NH4)2SO4 -> 2.01 H2O + 0.106 O2 + 6 Pichia_pastoris', 'Glucose', X=(1-yield_HemeB*1.1)*Y_b,
        correct_atomic_balance=True
    )
    cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=(1-yield_HemeB*1.1)*Y_b)
    
    respiration_reaction1 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', 1 - Y_b,
        check_atomic_balance=True
    )

    respiration_reaction2 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', 1 - cell_growth_reaction.X - fermentation_reaction[2].X * 1.1,
        check_atomic_balance=True
    )

    bst.settings.chemicals.set_alias('Pichia_pastoris', 'cellmass')
    RXN = bst.ReactionSystem(
        fermentation_reaction,
        bst.PRxn([cell_growth_reaction, respiration_reaction2])
    )
    RXN.show()
    
    
    """
    Upstream Process
    """
    M301 = bst.MixTank('M301', ins=[SeedIn1,'Water1'], outs='M301Out', tau=16)
    @M301.add_specification(run=True)
    def update_seed1_inputs():
        target_stream = bst.Stream(**s.SeedSolution1)
        SeedIn1.imass['Seed'] = target_stream.imass['Seed']
        M301.ins[1].imass['H2O'] = target_stream.imass['H2O']
        M301.ins[1].T = 25+273.15
    
    M302 = bst.MixTank('M302', ins=[SeedIn2,CultureIn,'Water2'], outs='M302Out', tau=16)
    @M302.add_specification(run=True)
    def update_culture_inputs():
        target_stream = bst.Stream(**s.SeedSolution2)
        SeedIn2.imass['Seed'] = target_stream.imass['Seed']
        M302.ins[2].imass['H2O'] = target_stream.imass['H2O']
        M302.ins[2].T = 25+273.15
        CultureIn.imass['Culture'] = target_stream.imass['SeedSolution']*(0.1+60+0.15191)/1000

    M303 = u.SeedHoldTank('M303', ins=[M301-0, M302-0], outs='M303Out')

    R301 = u.SeedTrain(
        'R301',
        ins=[M303-0],
        outs=[vent1, 'R301Out'],
        reactions=bst.PRxn([cell_growth_reaction, respiration_reaction1]),
        saccharification=None,
        T=32+273.15,
    )
    R301.add_specification(run=True)
    
    M304 = bst.MixTank('M304', ins=[Glucose,'Water3'], outs='M304Out', tau=16)
    
    @M304.add_specification(run=True)
    def update_water_content():
        M304.ins[1].imass['H2O'] = Glucose.imass['Glucose']/2
        M304.ins[1].T = 25+273.15
    
    T301 = bst.StorageTank('T301', ins=M304-0, outs='T301Out',tau=16*4+72)

    T302 = u.AmmoniaStorageTank('T302', ins=NH3_25wt, outs='T302Out')

    # checking more details...
    R302 = u.AeratedFermentation(
        'R302',
        ins=[R301-1, T301-0, T302-0, bst.Stream('FilteredAir', phase='g', P=2 * 101325)],
        outs=[vent2, 'Broth'],
        fermentation_reaction=fermentation_reaction,
        cell_growth_reaction=cell_growth_reaction,
        respiration_reaction=respiration_reaction2,
        #neutralization_reaction=neutralization_reaction,
        design='Stirred tank', method=method,theta_O2=theta_O2,
        V_max=V_max, Q_O2_consumption=Q_O2_consumption,
        dT_hx_loop=dT_hx_loop, T=T_operation,
        batch=True, reactions=RXN,
        kW_per_m3=agitation_power,
        tau=titer_HemeB/productivity_HemeB,
        cooler_pressure_drop=cooler_pressure_drop,
        compressor_isentropic_efficiency=compressor_isentropic_efficiency,
        P=1 * 101325,#optimize_power=True,
    )
    R302.target_titer = titer_HemeB # g / L
    R302.target_productivity = productivity_HemeB # g / L / h
    R302.target_yield = yield_HemeB  # wt %

    @R302.add_specification(run=True)
    def update_reaction_time_and_yield():
        R302.tau = R302.target_titer / R302.target_productivity
        fermentation_reaction[2].product_yield('Heme_b', basis='wt', product_yield=R302.target_yield)

    """
    Downstream process
    """
    
    
    
    
    