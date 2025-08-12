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
    # 'create_LegH_Pretreatment_system',
    # 'create_LegH_Conversion_system',
    # 'create_LegH_CellRemoval_Concentration_system',
    # 'create_LegH_Purification_system',
    'create_LegH_system',
)
# %%

# @bst.SystemFactory(
#     ID='LegH_Conversion_sys',
#     ins=[s.SeedIn, s.CultureIn, s.Glucose, s.NH3_18wt],
#     outs=[s.LegH_1, s.vent1, s.vent2,],
#     fthermo=c.create_chemicals_LegH,
# )
# def create_LegH_Conversion_system(
#         ins, outs,
#         # reactions_passed=None, # Placeholder if reactions were to be passed externally
#         # One could add more configurable parameters here:
#         # V_max_fermenter=500, target_titer=7.27, etc.
#     ):
#     """
#     Creates the LegH (Leghemoglobin) conversion system.
#     This system is based on the process flow and parameters from test_LegH.py.
#     """
#     # Unpack input streams
#     SeedIn, CultureIn, Glucose, NH3_18wt = ins

#     # Unpack output streams
#     LegH_1, vent1, vent2 = outs

#     """
#     Fermentation Parameter
#     """
#     theta_O2 = 0.5 # Dissolved oxygen concentration [% saturation]
#     agitation_power = 0.985 # [kW / m3]
#     design = 'Stirred tank' # Reactor type
#     method = "Riet" # Name of method
#     T_operation = 273.15 + 32 # [K]
#     Q_O2_consumption = -110 * 4184 # [kJ/kmol]
#     dT_hx_loop = 8 # [degC]
#     cooler_pressure_drop = 20684 # [Pa]
#     compressor_isentropic_efficiency = 0.85
#     V_max = 500 # [m3]
#     titer = 7.27 # [g / L]
#     productivity = titer / 72 # [g / L / h]
#     LegH_yield = titer * 5 / 1300 # [by wt]
#     Y_b = 0.43 # [by wt]

#     """
#     Reactions
#     """
    
#     fermentation_reaction = bst.PRxn([
#         #           Reaction        Reactnat            Conversion           Check                  "
#         bst.Rxn('8 Glucose + 4 NH3 + 1 FeSO4 + 10.5 O2 -> Heme_b + 1 H2SO4 + 37 H2O + 14 CO2',
#                                     reactant = 'Glucose',X=LegH_yield*0.05,check_atomic_balance=True),
#         bst.Rxn('Glucose + (NH4)2SO4 + O2 -> Globin + NH3 + H2O',
#                                     reactant = 'Glucose', X= LegH_yield*0.05,correct_atomic_balance=True),
#         bst.Rxn('Glucose + FeSO4 + (NH4)2SO4 + O2 -> Leghemoglobin + NH3 + H2O', 
#                                     reactant = 'Glucose', X=LegH_yield,  correct_atomic_balance=True),
#         ])
#     fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=LegH_yield)

#     neutralization_reaction = bst.Rxn(
#         'H2SO4 + NH3 -> (NH4)2SO4', reactant = 'H2SO4', X=1,
#         correct_atomic_balance=True
#     )

#     cell_growth_reaction = bst.Rxn(
#         'Glucose -> H2O + CO2 + Pichia_pastoris', 'Glucose', X=Y_b,
#         correct_atomic_balance=True
#     )
#     cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=Y_b)

#     respiration_reaction = bst.Rxn(
#         'Glucose + O2 -> CO2 + H2O', 'Glucose', 1. - cell_growth_reaction.X - fermentation_reaction[2].X,
#         correct_atomic_balance=True
#     )

#     bst.settings.chemicals.set_alias('Pichia_pastoris', 'cellmass')
#     RXN = bst.ReactionSystem(
#         fermentation_reaction,
#         bst.PRxn([cell_growth_reaction, respiration_reaction])
#     )
#     RXN.show()
#     """
#     Seed Train
#     """
#     M101= bst.Mixer(ins=[SeedIn,CultureIn])
#     R101 = u.SeedTrain(
#         'R101',
#         ins=[M101-0],
#         outs=[vent1, 'SeedOut'],
#         reactions=bst.PRxn([cell_growth_reaction, respiration_reaction]),
#         saccharification=None,
#         T=32+273.15,
#     )
    
#     R102 = u.AeratedFermentation(
#         'R102',
#         ins=[R101-1, Glucose, NH3_18wt, bst.Stream('FilteredAir', phase='g', P = 2 * 101325)],
#         outs=[vent2, LegH_1],
#         fermentation_reaction=fermentation_reaction,
#         cell_growth_reaction=cell_growth_reaction,
#         respiration_reaction=respiration_reaction,
#         neutralization_reaction=neutralization_reaction,
#         design='Stirred tank', method=method,theta_O2=theta_O2,
#         V_max=V_max, Q_O2_consumption=Q_O2_consumption,
#         dT_hx_loop=dT_hx_loop, T=T_operation,
#         batch=True, reactions=RXN,
#         kW_per_m3=agitation_power,
#         tau=titer/productivity,
#         cooler_pressure_drop=cooler_pressure_drop,
#         compressor_isentropic_efficiency=compressor_isentropic_efficiency,
#         P=1 * 101325,#optimize_power=True,
#     )
#     R102.target_titer = titer # g / L
#     R102.target_productivity = productivity # g / L / h
#     R102.target_yield = LegH_yield  # wt %

#     @R102.add_specification(run=True)
#     def update_reaction_time_and_yield():
#         R102.tau = R102.target_titer / R102.target_productivity
#         fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=R102.target_yield)

#     return LegH_1, vent1, vent2

# @bst.SystemFactory(
#     ID='LegH_CellRemoval_Concentration_sys',
#     ins=[s.LegH_1],
#     outs=[s.LegH_2, s.effluent1, s.effluent2,],
#     fthermo=c.create_chemicals_LegH, 
# )
# def create_LegH_CellRemoval_Concentration_system(
#         ins, outs,
#         # reactions_passed=None, # Placeholder if reactions were to be passed externally
#         # One could add more configurable parameters here:
#         # V_max_fermenter=500, target_titer=7.27, etc.
#     ):
#     """
#     Creates the LegH (Leghemoglobin) cell removal system.
#     This system is based on the process flow and parameters from test_LegH.py.
#     """
#     # Unpack input streams
#     LegH_1 = ins

#     # Unpack output streams
#     LegH_2, effluent1, effluent2 = outs

    
#     U201 = u.CellDisruption(
#         'U201',
#         ins=LegH_1,
#         outs='U201Out',
#     )

#     C201 = u.ProteinCentrifuge(
#         'C201',
#         ins = U201-0,
#         outs = (effluent1, 'C201Out'),
#         moisture_content = 0.5,
#         split = (1, 0.99, 1, 1),
#         order = ('Glucose','cellmass', 'LeghemoglobinIntre','GlobinIntre'),
#     )

#     E201 = u.Evaporator(
#         'E201',
#         ins = C201-1,
#         outs = ('E201Out',effluent2),
#         P = (101325, 73581, 50892, 32777, 20000),
#         V = 0.1,
#         V_definition = 'First-effect',
#     )

#     return LegH_2, effluent1, effluent2

# @bst.SystemFactory(
#     ID='LegH_Purification_sys',
#     ins=[s.SeedIn, s.WashingSolution, s.Elution, s.NFBuffer
#         ],
#     outs=[s.LegH_3, s.effluent3, s.effluent4, s.effluent5, s.effluent6], 
#     fthermo=c.create_chemicals_LegH,
# )
# def create_LegH_Purification_system(
#     ins, outs,
#     # reactions_passed=None, # Placeholder if reactions were to be passed externally
#     # One could add more configurable parameters here:
#     # V_max_fermenter=500, target_titer=7.27, etc.
# ):
#     '''
#     creates the LegH (Leghemoglobin) purification system.
#     '''
#     # Unpack input streams
#     SeedIn, WashingSolution, Elution, NFBuffer = ins
#     # Unpack output streams
#     LegH_3, effluent3, effluent4, effluent5, effluent6 = outs
#     #

#     U301 = u.Diafiltration(
#         'U301',
#         ins = (LegH_3, bst.Stream('WashingSolution', DiaBuffer=(LegH_3).imass['H2O']*4, units='kg/hr', T=25+273.15)),
#         outs = ('DF1Out',effluent3),
#         TargetProduct_ID = 'Leghemoglobin',
#         Salt_ID = c.chemical_groups['Salts'],
#         OtherLargeMolecules_ID = c.chemical_groups['OtherLargeMolecules'],
#         DefaultSolutes_ID = c.chemical_groups['DefaultSolutes'],
#     )
#     LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
#     LegH_sys.simulate()
#     U302 = u.IonExchange(
#         'U302',
#         ins = (U301-0, 
#             bst.Stream('Elution', IEXBuffer=(U301-0).imass['H2O']/2, units='kg/hr', T=25+273.15)),
#         outs = ('IEX1Out',effluent4),
#         TargetProduct_ID = 'Leghemoglobin',
#         BoundImpurity_ID=c.chemical_groups['BoundImpurities'],
#         ElutionBuffer_Defining_Component_ID =c.chemical_groups['ElutionBuffer'],
#     )
#     LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
#     LegH_sys.simulate()
#     U303 = u.Diafiltration(
#         'U303',
#         ins = (U302-0, bst.Stream('NFBuffer', 
#                     NanoBuffer=1.1*0.05*1000*
#                     (U302-0).imass['Leghemoglobin']/(0.25*bst.Chemical('TrehaloseDH',search_ID='6138-23-4', phase='l', default=True).MW),
#                     units='kg/hr', T=25+273.15)),
#         outs = ('NF1Out',effluent5),
#         TargetProduct_ID = 'Leghemoglobin',
#         membrane_cost_USD_per_m2=10000, # Nanomembrane cost
#         Salt_ID = c.chemical_groups['Salts'],
#         OtherLargeMolecules_ID = c.chemical_groups['OtherLargeMolecules'],
#         DefaultSolutes_ID = c.chemical_groups['DefaultSolutes'],
#         TargetProduct_Retention=0.995, Salt_Retention=0.1,
#         OtherLargeMolecules_Retention=0.99, DefaultSolutes_Retention=0.15,
#         FeedWater_Recovery_to_Permeate=0.2,
#         TMP_bar= 5
#     )
#     LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
#     LegH_sys.simulate()
#     E304 = bst.SprayDryer(
#         'E304',
#         ins=U303-0,
#         outs=(effluent6, LegH_3),
#         moisture_content=0.05,  # 5% moisture content in the final product
#     )
#     LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
#     LegH_sys.simulate()

#     WashingSolution = (1-U301)
#     Elution = (1-U302)
#     NFBuffer = (1-U303)

#     return LegH_3, effluent3, effluent4, effluent5, effluent6, WashingSolution, Elution, NFBuffer




@bst.SystemFactory(
    ID='LegH_sys',
    ins=[s.SeedIn, s.CultureIn, s.Glucose, s.NH3_25wt,
        ],
    outs=[s.LegH_3, s.vent1, s.vent2, s.effluent1, s.effluent2,
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
    
    # Update all input stream prices before simulation
    # NOW update stream prices (after they're properly created)

    s.update_all_input_stream_prices()

    # Unpack input streams
    SeedIn, CultureIn, Glucose, NH3_25wt = ins

    # Unpack output streams
    (LegH_3, vent1, vent2, effluent1, effluent2, 
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
    V_max = 500 # [m3]
    titer = 7.27 # [g / L]
    productivity = titer / 72 # [g / L / h]
    LegH_yield = titer * 5 / 1300 # [by wt]
    Y_b = 0.43 # [by wt]

    """
    Reactions
    """
    
    fermentation_reaction = bst.PRxn([
        #           Reaction        Reactnat            Conversion           Check                  "
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

    """
    Upstream Process
    """
    M1= bst.Mixer(ins=[SeedIn,CultureIn])
    ST1 = u.SeedTrain(
        'ST1',
        ins=[M1-0],
        outs=[vent1, 'SeedOut'],
        reactions=bst.PRxn([cell_growth_reaction, respiration_reaction]),
        saccharification=None,
        T=32+273.15,
    )
    
    AB1 = u.AeratedFermentation(
        'AB1',
        ins=[ST1-1, Glucose, NH3_25wt, bst.Stream('FilteredAir', phase='g', P = 2 * 101325)],
        outs=[vent2, 'AB1Out'],
        fermentation_reaction=fermentation_reaction,
        cell_growth_reaction=cell_growth_reaction,
        respiration_reaction=respiration_reaction,
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
    AB1.target_titer = titer # g / L
    AB1.target_productivity = productivity # g / L / h
    AB1.target_yield = LegH_yield  # wt %

    @AB1.add_specification(run=True)
    def update_reaction_time_and_yield():
        AB1.tau = AB1.target_titer / AB1.target_productivity
        fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=AB1.target_yield)
    """
    Downstream process
    """
    CD1 = u.CellDisruption(
        'CD1',
        ins=AB1-1,
        outs='CD1Out',
    )

    PC1 = u.ProteinCentrifuge(
        'PC1',
        ins = CD1-0,
        outs = (effluent1, 'PC1Out'),
        moisture_content = 0.5,
        split = (1, 0.99, 1, 1),
        order = ('Glucose','cellmass', 'LeghemoglobinIntre','GlobinIntre'),
    )

    EV1 = u.Evaporator(
        'EV1',
        ins = PC1-1,
        outs = ('EV1Out',effluent2),
        P = (101325, 73581, 50892, 32777, 20000),
        V = 0.1,
        V_definition = 'First-effect',
    )

    LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
    LegH_sys.simulate()

    DF1 = u.Diafiltration(
        'DF1',
        ins = (EV1-0, bst.Stream('BufferA', BufferA=(EV1-0).imass['H2O']*4, units='kg/hr', T=25+273.15)),
        outs = ('DF1Out',effluent3),
        TargetProduct_ID = 'Leghemoglobin',
        Salt_ID = c.chemical_groups['Salts'],
        OtherLargeMolecules_ID = c.chemical_groups['OtherLargeMolecules'],
        DefaultSolutes_ID = c.chemical_groups['DefaultSolutes'],
    )
    LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
    LegH_sys.simulate()
    IEX1 = u.IonExchange(
        'IEX1',
        ins = (DF1-0, 
            bst.Stream('BufferB', BufferB=(DF1-0).imass['H2O']/2, units='kg/hr', T=25+273.15)),
        outs = ('IEX1Out',effluent4),
        TargetProduct_ID = 'Leghemoglobin',
        BoundImpurity_ID=c.chemical_groups['BoundImpurities'],
        ElutionBuffer_Defining_Component_ID =c.chemical_groups['ElutionBuffer'],
    )
    LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
    LegH_sys.simulate()
    NF1 = u.Diafiltration(
        'NF1',
        ins = (IEX1-0, bst.Stream('BufferC', 
                    BufferC=1.1*0.05*1000*
                    (IEX1-0).imass['Leghemoglobin']/(0.25*bst.Chemical('TrehaloseDH',search_ID='6138-23-4', phase='l', default=True).MW),
                    units='kg/hr', T=25+273.15)),
        outs = ('NF1Out',effluent5),
        TargetProduct_ID = 'Leghemoglobin',
        membrane_cost_USD_per_m2=10000, # Nanomembrane cost
        Salt_ID = c.chemical_groups['Salts'],
        OtherLargeMolecules_ID = c.chemical_groups['OtherLargeMolecules'],
        DefaultSolutes_ID = c.chemical_groups['DefaultSolutes'],
        TargetProduct_Retention=0.995, Salt_Retention=0.1,
        OtherLargeMolecules_Retention=0.99, DefaultSolutes_Retention=0.15,
        FeedWater_Recovery_to_Permeate=0.2,
        TMP_bar= 5
    )
    
    LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
    LegH_sys.simulate()

    SD1 = bst.SprayDryer(
        'SD1',
        ins=NF1-0,
        outs=(effluent6, LegH_3),
        moisture_content=0.05,  # 5% moisture content in the final product
    )

    LegH_sys = bst.main_flowsheet.create_system('LegH_sys')
    LegH_sys.simulate()
    
    BufferA = (1-DF1)
    BufferB = (1-IEX1)
    BufferC = (1-NF1)

    return LegH_3, vent1, vent2, effluent1, effluent2, effluent3, effluent4, effluent5, effluent6, BufferA, BufferB, BufferC

if __name__ == '__main__':
    # # Create the LegH system
    LegH_sys = create_LegH_system()
    sys = LegH_sys
    f = sys.flowsheet
    u = f.unit
    SD1 = u.SD1

    LegH_sys.simulate()
    LegH_sys.diagram(format='html')
    LegH_sys.show()
    # Check stream prices
    print(f"\nStream Prices:")
    print(f"SeedIn price: ${f.SeedIn.price:.4f}/kg")
    print(f"CultureIn price: ${f.CultureIn.price:.4f}/kg")
    print(f"Glucose price: ${f.Glucose.price:.4f}/kg")
    print(f"NH3_25wt price: ${f.NH3_25wt.price:.4f}/kg")

# %%
