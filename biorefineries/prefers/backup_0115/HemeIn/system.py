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
import re
from altair import Y
from biorefineries.prefers.process_settings import set_GWPCF, GWP_CFs, set_GWPCF_Multi, load_process_settings
from biorefineries.prefers.HemeIn import streams as s
import biosteam as bst
from docutils import Component
from pint import set_application_registry
from regex import B, P
from thermosteam import Stream
from biosteam import F
import thermosteam as tmo
import numpy as np
from biorefineries.prefers.HemeIn import chemicals as c
from biorefineries.prefers import units as u
from biorefineries.prefers.process_settings import price

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
    ins=[s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose,s.NH3_25wt,],
    outs=[s.vent1, s.vent2,
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
    SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt = ins

    # Unpack output streams
    (vent1, vent2) = outs    
    
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
    P_Heme_1st = 1e-8 # [g / L / h]
    P_pp_1st = 1e-8 # [g / L / h]
    P_ComsuptionG1 = 500 * (0.5/3) / tau1 # [g / L / h] glucose consumption rate
    P_cell_growth_rate1 = 150 * 0.35 / tau1 # [g / L / h] cell growth rate
    Y_Heme1 = P_Heme_1st / P_ComsuptionG1 # product yield on glucose
    Y_pp1 = P_pp_1st / P_ComsuptionG1 # proph
    Y_b1 = P_cell_growth_rate1 / P_ComsuptionG1 # biomass yield on glucose
    Y_r1 = 1 - Y_Heme1 - Y_pp1 - Y_b1 # respiration yield on glucose
    
    # 2nd
    tau2= 90 # [h] fermentation time
    P_Heme_2nd = 4.2/1e3 # [g / L / h]
    P_pp_2nd = 0.7/1e3 # [g / L / h]
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
    
    SF = 0.45 # Secretion fraction of HemeB

    """
    
    Reactions
    """
    fermentation_reaction = bst.PRxn([
        #           Reaction        Reactnat            Conversion           Check                  "
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                                    reactant = 'Glucose',X=Y_Heme*SF,check_atomic_balance=True),
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b_In + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                                    reactant = 'Glucose', X= Y_Heme*(1-SF),check_atomic_balance=True),
        bst.Rxn('1 Glucose + 0.70588 NH3  -> 0.17647 ProtoporphyrinIX + 0.618 O2 + 4.06 H2O',
                                    reactant = 'Glucose', X=Y_pp*SF,correct_atomic_balance=True),
        bst.Rxn('1 Glucose + 0.70588 NH3 -> 0.17647 ProtoporphyrinIX_In + 0.618 O2 + 4.06 H2O',
                                    reactant = 'Glucose', X=Y_pp*(1-SF),correct_atomic_balance=True),
        ])
    fermentation_reaction[0].product_yield('Heme_b', basis='wt', product_yield=Y_Heme*SF)
    fermentation_reaction[1].product_yield('Heme_b_In', basis='wt', product_yield=Y_Heme*(1-SF))
    fermentation_reaction[2].product_yield('ProtoporphyrinIX', basis='wt', product_yield=Y_pp*SF)
    fermentation_reaction[3].product_yield('ProtoporphyrinIX_In', basis='wt', product_yield=Y_pp*(1-SF))
    
    neutralization_reaction = bst.Rxn(
        'H2SO4 + 2 NH3 -> (NH4)2SO4', reactant = 'H2SO4', X=1,
        check_atomic_balance=True
    )
    
    cell_growth_reactionCG = bst.Rxn(
        'Glucose + 1.32 NH3 + 1.32 H2O-> 6 Corynebacterium_glutamicum + 2.37 O2', 'Glucose', X=Y_b,
        correct_atomic_balance=True
    )
    cell_growth_reactionCG.product_yield('Corynebacterium_glutamicum', basis='wt', product_yield=Y_b)
    
    respiration_reactionGC1 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', 1 - Y_b,
        check_atomic_balance=True
    )

    respiration_reactionGC2 = bst.Rxn(
        'Glucose + 6 O2 -> 6 CO2 + 6 H2O', 'Glucose', (1 - Y_b - Y_pp - Y_Heme)*0.95,
        check_atomic_balance=True
    )
    
    respiration_reactionGC3 = bst.Rxn(
        'Glucose -> 2 EtOH + 2 CO2', 'Glucose', (1 - Y_b - Y_pp - Y_Heme)*0.05,
        check_atomic_balance=True
    )
    
    RXN = bst.ReactionSystem(
        fermentation_reaction,
        bst.PRxn([cell_growth_reactionCG, respiration_reactionGC2, respiration_reactionGC3]),
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
        reactions=bst.PRxn([cell_growth_reactionCG, respiration_reactionGC2, respiration_reactionGC3]),
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
        cell_growth_reaction=cell_growth_reactionCG,
        respiration_reaction=respiration_reactionGC2,
        neutralization_reaction=neutralization_reaction,
        design='Stirred tank', method=method,theta_O2=theta_O2,
        V_max=V_max, Q_O2_consumption=Q_O2_consumption,
        dT_hx_loop=dT_hx_loop, T=T_operation,
        batch=True, reactions=RXN,
        kW_per_m3=agitation_power,
        tau=tau,
        cooler_pressure_drop=cooler_pressure_drop,
        compressor_isentropic_efficiency=compressor_isentropic_efficiency,
        P=1 * 101325,
    )
    R302.target_titer = tau*P_Heme # g / L
    R302.target_productivity = P_Heme # g / L / h
    R302.target_yield = Y_Heme  # wt %

    @R302.add_specification(run=True)
    def update_reaction_time_and_yield():
        R302.tau = R302.target_titer / R302.target_productivity
        fermentation_reaction[0].product_yield('Heme_b', basis='wt', product_yield=R302.target_yield*SF)
        fermentation_reaction[1].product_yield('Heme_b_In', basis='wt', product_yield=R302.target_yield*(1-SF))
        fermentation_reaction[2].product_yield('ProtoporphyrinIX', basis='wt', product_yield=Y_pp*SF)
        fermentation_reaction[3].product_yield('ProtoporphyrinIX_In', basis='wt', product_yield=Y_pp*(1-SF))

    """
    Downstream process
    """
    
    # Define order chemicals outside the unit for easy access in specification
    centrifuge_order = ('Corynebacterium_glutamicum', 'Heme_b_In', 'ProtoporphyrinIX_In')
    S401 = u.Centrifuge(
        'S401',
        ins=R302-1,
        outs=('Deposite', 'Supernatant'),
        moisture_content = 0.4,
        split = (0.999, 0.999, 0.999),
        order = centrifuge_order,
    )
    @S401.add_specification(run=True)
    def update_centrifuge_splits():
        inlet = S401.ins[0]
        
        # First run to get initial outlet composition
        S401._run()
        
        outlets_solid = S401.outs[0]
        
        # Calculate water residual ratio (fraction of water going to solids)
        if inlet.imass['H2O'] > 0:
            waterresidual = outlets_solid.imass['H2O'] / inlet.imass['H2O']
        else:
            waterresidual = 0.
        
        # Get all chemical IDs
        order_chemicals = set(centrifuge_order)
        
        # Update splits: target chemicals go to solids (0.999), others follow water
        for chem in S401.chemicals:
            chem_id = chem.ID
            if chem_id in order_chemicals:
                S401.isplit[chem_id] = 0.999  # These go to solids stream
            elif chem_id == 'H2O':
                pass  # Water is handled by moisture_content
            else:
                S401.isplit[chem_id] = waterresidual  # Follow water distribution
        
        # Re-run with updated splits
        S401._run()

    
    S402 = u.CellDisruption(
        'S402', ins=S401-0, outs='DisruptedDeposite', 
        Cell_ID ='Corynebacterium_glutamicum',
        component_fractions={
                            'Protein': 0.45,
                            'Cellulose': 0.22,
                            'Xylan': 0.15,
                            'OleicAcid': 0.08,
                            'RNA':0.10,
        }
    )

    S403 = bst.Centrifuge(
        'S403', ins=S402-0, outs=('Debris','IntraProduct'),
        moisture_content=0.20,
        split=(1, 1, 1, 1, 1, 1, 1),
        order=('Heme_b_In', 'ProtoporphyrinIX_In','Corynebacterium_glutamicum','Protein','Cellulose', 'Xylan','OleicAcid','RNA'),
    )
    
    M401 = bst.MixTank('M401', ins=('','Water5'), 
                    outs='RAEquilibrium', tau=1)
    M402 = bst.MixTank('M402', ins=('','Water6'),
                    outs='RAElution', tau=1)
    M403 = bst.MixTank('M403', ins=('','Water6'),
                    outs='RARegeneration', tau=1)
    
    H401 = bst.HXutility(
        'H401',
        ins=M401-0,
        outs='H401Out',
        T=30+273.15,  # Cool to 30°C
        heat_only=True,
    )
    H402 = bst.HXutility(
        'H402',
        ins=M402-0,
        outs='H402Out',
        T=40+273.15,  # Cool to 40°C
        heat_only=True,
    )
    H403 = bst.HXutility(
        'H403',
        ins=M403-0,
        outs='H403Out',
        T=30+273.15,  # Cool to 30°C
        heat_only=True,
    )
    
    U401 = u.ResinAdsorption(
        'U401',
        ins = (S401-1, H401-0, H402-0, H403-0),
        outs = ('ExtraProduct','FlowthroughWaste','WashWaste','RegenerationWaste'),
        TargetProduct_IDs = ('Heme_b','ProtoporphyrinIX'),
        TargetProduct_Yield=0.98,
        BoundImpurity_IDs=(c.chemical_groups['Salts'],'Glucose','Glycine','Ethanol'),
        BoundImpurity_Removal=0.999,
        NonBinding_Carryover=0.05,
        wash_CV=3,
        elution_CV=0.05,
        regeneration_CV=0.05,
    ) 
    
    @M401.add_specification(run=True)
    def update_RAEquilibrium_initial():
        M401.ins[1].imass['H2O'] = (S401-1).imass['H2O']*U401.wash_CV
        M401.ins[0].imass['NaCl'] = (S401-1).imass['H2O']*U401.wash_CV*0.002/0.998
    @M402.add_specification(run=True)
    def update_RAElution_initial():
        M402.ins[1].imass['H2O'] = (S401-1).imass['H2O']*U401.elution_CV
        M402.ins[0].imol['NaOH'] = (S401-1).imass['H2O']*U401.elution_CV*0.2/1000
    @M403.add_specification(run=True)
    def update_RARegeneration_initial():
        # M403.ins[1].imass['H2O'] = (S401-1).imass['H2O']*U401.regeneration_CV
        # M403.ins[0].imol['NaOH'] = (S401-1).imass['H2O']*U401.regeneration_CV*1/1000
        M403.ins[1].imass['H2O'] = (S401-1).imass['H2O']*U401.regeneration_CV*0.3
        M403.ins[0].imass['Ethanol'] = (S401-1).imass['H2O']*U401.regeneration_CV*0.7
    

if __name__ == '__main__':
    bst.preferences.N = 50
    Hm_sys = create_Heme_system()
    sys = Hm_sys
    f = sys.flowsheet
    u = f.unit
    ss = f.stream
    sys.operating_hours = 8000
    sys.simulate()
    sys.diagram(format='html',display=True)