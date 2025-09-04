# -*- coding: utf-8 -*-
"""
Created on 2025-06-04 14:26:14

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

#from biorefineries.animal_bedding import _system
import biosteam as bst
from thermosteam import Stream

import thermosteam as tmo
import numpy as np
from biorefineries.prefers import _chemicals as c, _units as u ,_streams as s
from biorefineries.prefers._process_settings import price  # ADD THIS IMPORT

# %% Settings
bst.settings.set_thermo(c.create_chemicals_LegH(), skip_checks=True)
bst.preferences.classic_mode()

# %%
__all__ = (
    'create_LegH_system',
)
# %%
@bst.SystemFactory(
    ID='LegH_sys',
    ins=[s.SeedIn, s.CultureIn, s.Glucose, s.NH3_25wt,s.BufferA, s.BufferB, s.BufferC
        ],
    outs=[s.LegH_3, s.vent1, s.vent2, s.effluent1, #s.effluent2,
        s.effluent3, s.effluent4, s.effluent5, s.effluent6],
    fthermo=c.create_chemicals_LegH, # Pass the function itself
)
def create_LegH_system(
        ins, outs,
        # reactions_passed=None, # Placeholder if reactions were to be passed externally
        # One could add more configurable parameters here:
        # V_max_fermenter=500, target_titer=7.27, etc.
    ):
    """
    Creates the LegH (Leghemoglobin) production system.
    This system is based on the process flow and parameters from test_LegH.py.
    """
    bst.preferences.N=50
    # Update all input stream prices before system creation to avoid conflicts
    # This is now done at the module level in _streams.py to prevent ID conflicts
    
    # Uncomment the line below if you need to update prices dynamically:
    #s.update_all_input_stream_prices()
    # Unpack input streams
    SeedIn, CultureIn, Glucose, NH3_25wt, BufferA, BufferB, BufferC = ins

    # Unpack output streams
    (LegH_3, vent1, vent2, effluent1, #effluent2, 
    effluent3, effluent4, effluent5, effluent6) = outs
    """
    Fermentation Parameter
    """
    theta_O2 = 0.5 # Dissolved oxygen concentration [% saturation]
    agitation_power = 0.985 # [kW / m3]
    design = 'Stirred tank' # Reactor type
    method = "Riet" # Name of method
    T_operation = 273.15 + 32 # [K]
    Q_O2_consumption = -110 * 4184 # [kJ/kmol]
    dT_hx_loop = 8 # [degC]
    cooler_pressure_drop = 20684 # [Pa]
    compressor_isentropic_efficiency = 0.85
    V_max = 150 # [m3] #here cause pressure vessel design problem
    titer = 7.27 # [g / L]
    productivity = titer / 72 # [g / L / h]
    Y_p = titer * 5 / 1300 # [by wt] 3 wt%
    Y_b = 0.43 # [by wt] 
    LegH_yield = Y_p/0.517
    """
    Reactions
    """
    
    # fermentation_reaction = bst.PRxn([
    #     #           Reaction        Reactnat            Conversion           Check                  "
    #     bst.Rxn('8 Glucose + 6 NH3 + 1 FeSO4 + 10.5 O2 -> Heme_b + 1 (NH4)2SO4 + 37 H2O + 14 CO2',
    #                                 reactant = 'Glucose',X=LegH_yield*0.05,check_atomic_balance=True),
    #     bst.Rxn('Glucose + (NH4)2SO4 + NH3 -> Globin + CO2 + H2O',
    #                                 reactant = 'Glucose', X= LegH_yield*0.05,correct_atomic_balance=True),
    #     bst.Rxn('Glucose + FeSO4 + (NH4)2SO4 + NH3 -> Leghemoglobin + CO2  + H2O',
    #                                 reactant = 'Glucose', X=LegH_yield,  correct_atomic_balance=True),
    #     ])
    fermentation_reaction = bst.PRxn([
        #           Reaction        Reactnat            Conversion           Check                  "
        bst.Rxn('1 Glucose + 1.05882 NH3 + 0.17647 FeSO4  -> 0.17647 Heme_b + 0.617647 O2 + 0.17647 (NH4)2SO4 + 4.05882 H2O',
                                    reactant = 'Glucose',X=LegH_yield*0.05,check_atomic_balance=True),
        bst.Rxn('Glucose + 0.01646 (NH4)2SO4 + 1.61317 NH3 -> 6 Globin + 0.28807 O2 + 3.68724 H2O',
                                    reactant = 'Glucose', X= LegH_yield*0.05,check_atomic_balance=True),
        bst.Rxn('Glucose + 0.00786 FeSO4 + 0.00786 (NH4)2SO4 + 1.58847 NH3 -> 6 Leghemoglobin + 0.30275 O2  + 3.70380 H2O',
                                    reactant = 'Glucose', X=LegH_yield,  check_atomic_balance=True),
        ])
    fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=LegH_yield)

    neutralization_reaction = bst.Rxn(
        'H2SO4 + 2 NH3 -> (NH4)2SO4', reactant = 'H2SO4', X=1,
        check_atomic_balance=True
    )

    cell_growth_reaction = bst.Rxn(
        'Glucose + 0.8364 NH3 + 0.0108 (NH4)2SO4 -> 2.01 H2O + 0.106 O2 + 6 Pichia_pastoris', 'Glucose', X=(1-LegH_yield*1.1)*Y_b,
        correct_atomic_balance=True
    )
    cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=(1-LegH_yield*1.1)*Y_b)

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
    M301 = bst.Mixer('M301', ins=[SeedIn,CultureIn], outs='M301Out')
    R301 = u.SeedTrain(
        'R301',
        ins=[M301-0],
        outs=[vent1, 'R301Out'],
        reactions=bst.PRxn([cell_growth_reaction, respiration_reaction1]),
        saccharification=None,
        T=32+273.15,
    )
    R301.add_specification(run=True)

    M302 = bst.Mixer('M302', ins=[R301-1, Glucose, NH3_25wt], outs='M302Out')

    R302 = u.AeratedFermentation(
        'R302',
        ins=[M302-0, bst.Stream('FilteredAir', phase='g', P = 2 * 101325)],
        outs=[vent2, 'R302Out'],
        fermentation_reaction=fermentation_reaction,
        cell_growth_reaction=cell_growth_reaction,
        respiration_reaction=respiration_reaction2,
        neutralization_reaction=neutralization_reaction,
        design='Stirred tank', method=method,theta_O2=theta_O2,
        V_max=V_max, Q_O2_consumption=Q_O2_consumption,
        dT_hx_loop=dT_hx_loop, T=T_operation,
        batch=True, reactions=RXN,
        kW_per_m3=agitation_power,
        tau=titer/productivity,
        cooler_pressure_drop=cooler_pressure_drop,
        compressor_isentropic_efficiency=compressor_isentropic_efficiency,
        P=1 * 101325,#optimize_power=True,
    )
    R302.target_titer = titer # g / L
    R302.target_productivity = productivity # g / L / h
    R302.target_yield = LegH_yield  # wt %

    @R302.add_specification(run=True)
    def update_reaction_time_and_yield():
        R302.tau = R302.target_titer / R302.target_productivity
        fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=R302.target_yield)
    """
    Downstream process
    """
    S401 = u.CellDisruption(
        'S401',
        ins=R302-1,
        outs='S401Out',
    )

    S402 = u.ProteinCentrifuge(
        'S402',
        ins = S401-0,
        outs = (effluent1, 'S402Out'),
        moisture_content = 0.20,  # 20% moisture content in the final product
        split = (0, 0.2, 1, 1),
        order = ('Glucose','cellmass', 'LeghemoglobinIntre','GlobinIntre'),
    )
    S402.add_specification(run=True)
    

    # E401 = u.Evaporator(
    #     'E401',
    #     ins = S402-1,
    #     outs = ('E401Out',effluent2),
    #     P = (101325, 73581, 50892, 32777, 20000),  # Reduced to 3 effects to increase vessel size
    #     V = 0.1,  # Increased vapor fraction to create larger vessels
    #     V_definition = 'First-effect',  # Changed to overall for better load distribution
    # )
    # E401.add_specification(run=True)

    U401 = u.Diafiltration(
        'U401',
        ins = (S402-1, BufferA),
               #bst.Stream('BufferA', BufferA=(E401-0).imass['H2O']*4, units='kg/hr', T=25+273.15)),
        outs = ('U401Out',effluent3),
        TargetProduct_ID = 'Leghemoglobin',
        Salt_ID = c.chemical_groups['Salts'],
        OtherLargeMolecules_ID = c.chemical_groups['OtherLargeMolecules'],
        TMP_bar =5
    )
    @U401.add_specification(run=True)
    def update_BufferA():
        BufferA.imass['H2O'] = (S402-1).imass['H2O']*8#4
        BufferA.T = 25+273.15


    U402 = u.IonExchange(
        'U402',
        ins = (U401-0,BufferB,) ,
            #bst.Stream('BufferB', BufferB=(U401-0).imass['H2O']/2, units='kg/hr', T=25+273.15)),
        outs = ('U402Out',effluent4),
        TargetProduct_IDs = c.chemical_groups['LegHIngredients'],
        BoundImpurity_IDs=c.chemical_groups['BoundImpurities'],
        #ElutionBuffer_Defining_Component_ID =c.chemical_groups['ElutionBuffer'],
    )
    @U402.add_specification(run=True)
    def update_BufferB():
        BufferB.imass['H2O'] = (U401-0).imass['H2O']*5#/2
        BufferB.T = 25+273.15


    U403 = u.Diafiltration(
        'U403',
        ins = (U402-0, BufferC),
            #    bst.Stream('BufferC', 
            #         BufferC=1.1*0.05*1000*
            #         (U402-0).imass['Leghemoglobin']/(0.25*bst.Chemical('TrehaloseDH',search_ID='6138-23-4', phase='l', default=True).MW),
            #         units='kg/hr', T=25+273.15)),
        outs = ('U403Out',effluent5),
        TargetProduct_ID = 'Leghemoglobin',
        membrane_cost_USD_per_m2=1000, # Nanomembrane cost
        Salt_ID = c.chemical_groups['Salts'],
        OtherLargeMolecules_ID = c.chemical_groups['OtherLargeMolecules'],
        TargetProduct_Retention=0.995, Salt_Retention=0.1,
        OtherLargeMolecules_Retention=0.995, DefaultSolutes_Retention=0.15,
        FeedWater_Recovery_to_Permeate=0.2,
        TMP_bar= 20
    )
    @U403.add_specification(run=True)
    def update_BufferC():
        BufferC.imass['H2O'] = 4*1.1*0.20*1000*(U402-0).imass['Leghemoglobin']/(0.25*bst.Chemical('TrehaloseDH',search_ID='6138-23-4', phase='l', default=True).MW)
        BufferC.T = 25+273.15


    S403 = bst.SprayDryer(
        'S403', 
        ins=U403-0,
        outs=(effluent6, LegH_3),
        moisture_content=0.9,  # 90% moisture content in the final product
    )
    
    s.update_all_input_stream_prices(streamlist=[SeedIn, CultureIn, Glucose, NH3_25wt, BufferA, BufferB, BufferC])

    return LegH_3, vent1, vent2, effluent1, effluent3, effluent4, effluent5, effluent6, BufferA, BufferB, BufferC
#effluent2
# %%
if __name__ == '__main__':
    # # Create the LegH system
    bst.preferences.N=50
    LegH_sys = create_LegH_system()
    sys = LegH_sys
    f = sys.flowsheet
    u = f.unit

    LegH_sys.simulate()
    LegH_sys.show()
    LegH_sys.diagram(format='html',display=True,)


# %%
