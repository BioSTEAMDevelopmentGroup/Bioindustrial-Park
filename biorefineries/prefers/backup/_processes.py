# -*- coding: utf-8 -*-
"""
Created on 2025-08-15 18:39:59

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

# -*- coding: utf-8 -*-
# In biorefineries/legh/_processes.py

import biosteam as bst
from thermosteam import Stream
# Assuming you have these custom modules set up
from biorefineries.prefers import _chemicals as c
from biorefineries.prefers import _units as u

__all__ = (
    'create_LegH_conversion_process',
    'create_LegH_concentration_process',
    'create_LegH_purification_process',
)

# %% 1. Conversion Process

def create_LegH_conversion_process(ins, outs):
    """
    Creates the Leghemoglobin (LegH) conversion process area.
    This includes the seed train and the main fermentation.
    """
    # Unpack input and output streams
    SeedIn, CultureIn, Glucose, NH3_25wt = ins
    vent1, vent2, fermentation_broth = outs

    # --- Fermentation Parameters (centralized from original code) ---
    theta_O2 = 0.5
    agitation_power = 0.985
    T_operation = 273.15 + 32
    Q_O2_consumption = -110 * 4184
    dT_hx_loop = 8
    cooler_pressure_drop = 20684
    compressor_isentropic_efficiency = 0.85
    V_max = 100
    titer = 7.27
    productivity = titer / 72
    LegH_yield = titer * 5 / 1300
    Y_b = 0.43

    # --- Reactions ---
    fermentation_reaction = bst.PRxn([
        bst.Rxn('8 Glucose + 6 NH3 + 1 FeSO4 + 10.5 O2 -> Heme_b + 1 (NH4)2SO4 + 37 H2O + 14 CO2',
                reactant='Glucose', X=LegH_yield * 0.05, check_atomic_balance=True),
        bst.Rxn('Glucose + (NH4)2SO4 + NH3 -> Globin + CO2 + H2O',
                reactant='Glucose', X=LegH_yield * 0.05, correct_atomic_balance=True),
        bst.Rxn('Glucose + FeSO4 + (NH4)2SO4 + NH3 -> Leghemoglobin + CO2  + H2O',
                reactant='Glucose', X=LegH_yield, correct_atomic_balance=True),
    ])
    fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=LegH_yield)

    neutralization_reaction = bst.Rxn(
        'H2SO4 + NH3 -> (NH4)2SO4', reactant='H2SO4', X=1, correct_atomic_balance=True
    )

    cell_growth_reaction = bst.Rxn(
        'Glucose -> H2O + CO2 + Pichia_pastoris', 'Glucose', X=Y_b, correct_atomic_balance=True
    )
    cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=Y_b)

    respiration_reaction1 = bst.Rxn(
        'Glucose + O2 -> CO2 + H2O', 'Glucose', 1. - cell_growth_reaction.X, correct_atomic_balance=True
    )

    respiration_reaction2 = bst.Rxn(
        'Glucose + O2 -> CO2 + H2O', 'Glucose', 1. - cell_growth_reaction.X - fermentation_reaction[2].X,
        correct_atomic_balance=True
    )

    bst.settings.chemicals.set_alias('Pichia_pastoris', 'cellmass')
    RXN = bst.ReactionSystem(
        fermentation_reaction,
        bst.PRxn([cell_growth_reaction, respiration_reaction2])
    )

    # --- Unit Operations ---
    M301 = bst.Mixer('M301', ins=[SeedIn, CultureIn])

    R301 = u.SeedTrain(
        'R301',
        ins=[M301-0],
        outs=[vent1, 'SeedOut'],
        reactions=bst.PRxn([cell_growth_reaction, respiration_reaction1]),
        saccharification=None,
        T=T_operation,
    )

    R302 = u.AeratedFermentation(
        'R302',
        ins=[R301-1, Glucose, NH3_25wt, bst.Stream('FilteredAir', phase='g', P=2 * 101325)],
        outs=[vent2, fermentation_broth],
        reactions=RXN, # Pass the complete reaction system here
        fermentation_reaction=fermentation_reaction,
        cell_growth_reaction=cell_growth_reaction,
        respiration_reaction=respiration_reaction2,
        neutralization_reaction=neutralization_reaction,
        design='Stirred tank', method="Riet", theta_O2=theta_O2,
        V_max=V_max, Q_O2_consumption=Q_O2_consumption,
        dT_hx_loop=dT_hx_loop, T=T_operation,
        batch=True,
        kW_per_m3=agitation_power,
        tau=titer/productivity,
        cooler_pressure_drop=cooler_pressure_drop,
        compressor_isentropic_efficiency=compressor_isentropic_efficiency,
        P=1 * 101325,
    )
    R302.target_titer = titer
    R302.target_productivity = productivity
    R302.target_yield = LegH_yield

    @R302.add_specification(run=True)
    def update_reaction_time_and_yield():
        R302.tau = R302.target_titer / R302.target_productivity
        fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=R302.target_yield)

# %% 2. Concentration Process

def create_LegH_concentration_process(ins, outs):
    """
    Creates the LegH concentration process area.
    This includes cell disruption, centrifugation, and evaporation.
    """
    # Unpack input and output streams
    fermentation_broth, = ins
    effluent1, effluent2, concentrated_legh = outs

    # --- Unit Operations ---
    S401 = u.CellDisruption(
        'S401',
        ins=fermentation_broth,
        outs='CD1Out',
    )

    S402 = u.ProteinCentrifuge(
        'S402',
        ins=S401-0,
        outs=(effluent1, 'S402Out'),
        moisture_content=0.3,
        split=(1, 0.99, 1, 1),
        order=('Glucose', 'cellmass', 'LeghemoglobinIntre', 'GlobinIntre'),
    )

    E401 = u.Evaporator(
        'E401',
        ins=S402-1,
        outs=(concentrated_legh, effluent2),
        P=(101325, 73581, 50892),
        V=0.3,
        V_definition='Overall',
    )

# %% 3. Purification Process

def create_LegH_purification_process(ins, outs):
    """
    Creates the LegH purification process area.
    This includes diafiltration, ion exchange, and spray drying.
    """
    # Unpack input and output streams
    concentrated_legh, BufferA, BufferB, BufferC = ins
    LegH_3, effluent3, effluent4, effluent5, effluent6 = outs

    # --- Unit Operations ---
    U401 = u.Diafiltration(
        'U401',
        ins=(concentrated_legh, BufferA),
        outs=('U401Out', effluent3),
        TargetProduct_ID='Leghemoglobin',
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
        DefaultSolutes_ID=c.chemical_groups['DefaultSolutes'],
    )

    U402 = u.IonExchange(
        'U402',
        ins=(U401-0, BufferB),
        outs=('U402Out', effluent4),
        TargetProduct_ID='Leghemoglobin',
        BoundImpurity_ID=c.chemical_groups['BoundImpurities'],
        ElutionBuffer_Defining_Component_ID=c.chemical_groups['ElutionBuffer'],
    )

    U403 = u.Diafiltration(
        'U403',
        ins=(U402-0, BufferC),
        outs=('U403Out', effluent5),
        TargetProduct_ID='Leghemoglobin',
        membrane_cost_USD_per_m2=10000,
        Salt_ID=c.chemical_groups['Salts'],
        OtherLargeMolecules_ID=c.chemical_groups['OtherLargeMolecules'],
        DefaultSolutes_ID=c.chemical_groups['DefaultSolutes'],
        TargetProduct_Retention=0.995, Salt_Retention=0.1,
        OtherLargeMolecules_Retention=0.99, DefaultSolutes_Retention=0.15,
        FeedWater_Recovery_to_Permeate=0.2,
        TMP_bar=5
    )

    S403 = bst.SprayDryer(
        'S403',
        ins=U403-0,
        outs=(effluent6, LegH_3),
        moisture_content=0.05,
    )