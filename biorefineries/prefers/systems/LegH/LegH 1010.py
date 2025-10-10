# -*- coding: utf-8 -*-
"""
Created on 2025-06-04 14:26:14

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

#from biorefineries.animal_bedding import _system
from biorefineries.prefers._process_settings import set_GWPCF, GWP_CFs,set_GWPCF_Multi,load_process_settings
import biosteam as bst
from pint import set_application_registry
from thermosteam import Stream
from biosteam import F
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
    ins=[s.SeedIn1, s.SeedIn2, s.CultureIn, s.Glucose, s.NH3_25wt,s.DfUltraBuffer, s.IXEquilibriumBuffer, s.IXElutionBuffer
         ,s.IXRegenerationSolution, s.DfNanoBuffer],
    outs=[s.LegH_3, s.vent1, s.vent2, s.effluent1, s.effluent2, s.effluent3,
    ],
    fthermo=c.create_chemicals_LegH, # Pass the function itself
)
def create_LegH_system(
        ins, outs,
        use_area_convention=False,
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
    SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, IXEquilibriumBuffer, IXElutionBuffer, IXRegenerationSolution, DfNanoBuffer = ins

    # Unpack output streams
    (LegH_3, vent1, vent2, effluent1, effluent2, effluent3) = outs
    
    set_GWPCF(Glucose, 'Glucose')
    set_GWPCF_Multi(SeedIn1, ['AmmoniumSulfate','Glucose','MagnesiumSulfate','KH2PO4'],[0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(SeedIn2, ['AmmoniumSulfate','Glucose','MagnesiumSulfate','KH2PO4'],[0.5/(0.5+1+0.05+0.3), 1/(0.5+1+0.05+0.3), 0.05/(0.5+1+0.05+0.3), 0.3/(0.5+1+0.05+0.3)])
    set_GWPCF_Multi(CultureIn, ['Glycine','Glucose','IronSulfate'],[0.1/(0.1+60+0.15191), 60/(0.1+60+0.15191), 0.15191/(0.1+60+0.15191)])
    set_GWPCF(NH3_25wt, 'Ammonia_SEA',dilution=0.25)
    set_GWPCF_Multi(DfUltraBuffer, ['KH2PO4','NaCl'], [0.8472, 0.1455])
    set_GWPCF_Multi(IXEquilibriumBuffer, ['KH2PO4','NaCl'], [0.8472, 0.1455])
    set_GWPCF_Multi(IXElutionBuffer, ['KH2PO4','NaCl','KCl'], [0.3616,0.5911, 0.0792])
    set_GWPCF(IXRegenerationSolution, 'NaOH')
    set_GWPCF_Multi(DfNanoBuffer, ['Na2HPO4','NaH2PO4'], [0.5420, 0.4580])

    load_process_settings()  # Load process settings to update prices and CFs
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
    V_max = 500 # [m3] #here cause pressure vessel design problem
    titer = 7.27 # [g / L]
    productivity = titer / 72 # [g / L / h]
    Y_p = titer * 5 / 1300 # [by wt] 3 wt%
    Y_b = 0.43 # [by wt] 
    LegH_yield = Y_p#/0.517 # yield based on glucose utilized for product formation
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
                                    reactant = 'Glucose',X=LegH_yield*0.01,check_atomic_balance=True),
        bst.Rxn('Glucose + 0.01646 (NH4)2SO4 + 1.61317 NH3 -> 6 Globin + 0.28807 O2 + 3.68724 H2O',
                                    reactant = 'Glucose', X= LegH_yield*0.04,check_atomic_balance=True),
        bst.Rxn('Glucose + 0.00786 FeSO4 + 0.00786 (NH4)2SO4 + 1.58847 NH3 -> 6 Leghemoglobin + 0.30275 O2  + 3.70380 H2O',
                                    reactant = 'Glucose', X=LegH_yield,  check_atomic_balance=True),
        ])
    fermentation_reaction[2].product_yield('Leghemoglobin', basis='wt', product_yield=LegH_yield)

    neutralization_reaction = bst.Rxn(
        'H2SO4 + 2 NH3 -> (NH4)2SO4', reactant = 'NH3', X=1,
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

    R302 = u.AeratedFermentation(
        'R302',
        ins=[R301-1, T301-0, T302-0, bst.Stream('FilteredAir', phase='g', P=2 * 101325)],
        outs=[vent2, 'Broth'],
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
        outs='DisruptedBroth',
    )

    S402 = u.ProteinCentrifuge(
        'S402',
        ins = S401-0,
        outs = ('Deposit', 'Supernatant'),
        moisture_content = 0.20,  # 20% moisture content in the final product
        split = (0, 0.995, 1, 1),
        order = ('Glucose','cellmass', 'LeghemoglobinIntre','GlobinIntre'),
    )
    S402.add_specification(run=True)

    #S403 = bst.SolidsSeparator('S403', ins=S402-0, outs=(effluent1, 'S403Out'), split=(1) , moisture_content=0.01)
    S403 = bst.ScrewPress('S403', ins=S402-0, outs=('CellMassWaste', 'S403Out'), split=(0.999) , moisture_content=0.01)
    
    S404 = bst.Splitter('S404', ins=S402-1, outs=('ResidualCellMass', 'S404Out'), split=(0.999) , order=('cellmass',))

    # E401 = u.Evaporator(
    #     'E401',
    #     ins = S402-1,
    #     outs = ('E401Out',effluent2),
    #     P = (101325, 73581, 50892, 32777, 20000),  # Reduced to 3 effects to increase vessel size
    #     V = 0.1,  # Increased vapor fraction to create larger vessels
    #     V_definition = 'First-effect',  # Changed to overall for better load distribution
    # )
    # E401.add_specification(run=True)

    M401 = bst.MixTank('M401', ins=(DfUltraBuffer, 'Water4'), 
                    outs='M401Out', tau=1) 
    @M401.add_specification(run=True)
    def update_DfUltraBuffer_initial():
        M401.ins[1].imass['H2O'] = (S402-1).imass['H2O']*4
        M401.ins[0].imol['DfUltraBuffer'] = ((S402-1).imass['H2O']*4)*(0.025+0.01+0.001)/1000

    H401 = bst.HXutility(
        'H401',
        ins=M401-0,
        outs='H401Out',
        T=5+273.15,  # Cool to 5°C
        cool_only=True,
    )

    U401 = u.Diafiltration(
        'U401',
        ins = (S404-1, H401-0),
        outs = ('U401Out','PermeateWasteUltra'),
        TargetProduct_ID = 'Leghemoglobin',
        Salt_ID = c.chemical_groups['Salts'],
        OtherLargeMolecules_ID = c.chemical_groups['OtherLargeMolecules'],
        TMP_bar1 = 3 ,#2~4
        TMP_bar2 = 2 ,#1.5~3
    )
    U401.add_specification(run=True)

    #S404 = bst.SolidsSeparator('S404', ins=U401-1, outs=(effluent2, 'S404Out'), split=(1) , moisture_content=0.01)
    #S405 = u.ReverseOsmosis('S405', ins=U401-1, outs=('S405Out',effluent2))

    M402 = bst.MixTank('M402', ins=(IXEquilibriumBuffer,'Water5'), 
                    outs='M402Out', tau=1)
    M403 = bst.MixTank('M403', ins=(IXElutionBuffer,'Water6'),
                    outs='M403Out', tau=1)
    M404 = bst.MixTank('M404', ins=(IXRegenerationSolution,'Water7'), 
                    outs='M404Out', tau=1)

    H402 = bst.HXutility(
        'H402',
        ins=M402-0,
        outs='H402Out',
        T=5+273.15,  # Cool to 5°C
        cool_only=True,
    )
    H403 = bst.HXutility(
        'H403',
        ins=M403-0,
        outs='H403Out',
        T=5+273.15,  # Cool to 5°C
        cool_only=True,
    )
    H404 = bst.HXutility(
        'H404', 
        ins=M404-0,
        outs='H404Out',
        T=5+273.15,  # Cool to 5°C
        cool_only=True,
    )

    U402 = u.IonExchangeCycle(
        'U402',
        ins = (U401-0, H402-0, H403-0, H404-0),
        outs = ('U402Out','FlowthroughWaste','WashWaste','RegenerationWaste'),
        TargetProduct_IDs = c.chemical_groups['LegHIngredients'],
        BoundImpurity_IDs=c.chemical_groups['BoundImpurities'],
    )    

    @M402.add_specification(run=True)
    def update_IXEquilibriumBuffer_initial():
        M402.ins[1].imass['H2O'] = (U401-0).imass['H2O']*U402.wash_CV
        M402.ins[0].imol['IXEquilibriumBuffer'] = ((U401-0).imass['H2O']*U402.wash_CV)*(0.025+0.01+0.001)/1000
    @M403.add_specification(run=True)
    def update_IXElutionBuffer_initial():
        M403.ins[1].imass['H2O'] = (U401-0).imass['H2O']*U402.elution_CV
        M403.ins[0].imol['IXElutionBuffer'] = ((U401-0).imass['H2O']*U402.elution_CV)*(0.025+1+0.1)/1000
    @M404.add_specification(run=True)
    def update_IXRegenerationSolution_initial():
        M404.ins[1].imass['H2O'] = (U401-0).imass['H2O']*U402.regeneration_CV
        M404.ins[0].imol['NaOH'] = ((U401-0).imass['H2O']*U402.regeneration_CV)*(0.5)/1000

    #S405 = bst.SolidsSeparator('S405', ins=U402-1, outs=(effluent3, 'S405Out'), split=(1) , moisture_content=0.01)
    #S406 = u.ReverseOsmosis('S406', ins=U402-1, outs=('S406Out',effluent3))

    M405 = bst.MixTank('M405', ins=(DfNanoBuffer,'Water8'), 
                    outs='M405Out', tau=1)
    @M405.add_specification(run=True)
    def update_DfNanoBuffer_initial():
        M405.ins[1].imass['H2O'] = (U402-0).imass['H2O']*2
        M405.ins[0].imol['DfNanoBuffer'] = ((U402-0).imass['H2O']*2)*(0.01+0.01)/1000

    H405 = bst.HXutility(
        'H405',
        ins=M405-0,
        outs='H405Out',
        T=5+273.15,  # Cool to 5°C
        cool_only=True,
    )

    U403 = u.Diafiltration(
        'U403',
        ins = (U402-0, H405-0),
        outs = ('U403Out','PermeateWasteNano'),
        TargetProduct_ID = 'Leghemoglobin',
        membrane_cost_USD_per_m2=250, # Nanomembrane cost
        Salt_ID = c.chemical_groups['Salts'],
        OtherLargeMolecules_ID = c.chemical_groups['OtherLargeMolecules'],
        TargetProduct_Retention=0.995, Salt_Retention=0.001,
        OtherLargeMolecules_Retention=0.995, DefaultSolutes_Retention=0.015,
        FeedWater_Recovery_to_Permeate=0.2,
        TMP_bar1= 15 ,# 10-25
        TMP_bar2= 4  ,# 3-6
        membrane_flux_LMH=25, # 10-40
    )
    U403.add_specification(run=True)

    #S406 = bst.SolidsSeparator('S406', ins=U403-1, outs=(effluent4, 'S406Out'), split=(1) , moisture_content=0.01)
    #S407 = u.ReverseOsmosis('S407', ins=U403-1, outs=('S407Out',effluent4))

    S408 = bst.SprayDryer(
        'S408', 
        ins=U403-0,
        outs=('EvaporatedWater', 'S408Out'),
        moisture_content=0.9,  # 90% moisture content in the final product
    )
    S408.add_specification(run=True)

    M501 = bst.MixTank('M501', ins=(S403-1,U401-1,U402-1,U402-2,U403-1), 
                    outs='M501Out', tau=1)
    M501.add_specification(run=True)

    T501 = u.SulfuricAcidStorageTank('T501', 
                        ins=bst.Stream('SulfuricAcid', H2SO4=0.98,H2O=0.02, 
                        units='kg/hr', price=price['H2SO4']), 
                        outs='T501Out')
    @T501.add_specification(run=True)
    def update_acid_flowrate():
        T501.ins[0].imol['H2SO4'] = U402.outs[3].imol['NaOH']/2*1.001
        T501.ins[0].T = 25+273.15

    M502 = bst.NeutralizationTank1('M502', ins=(U402-3,T501-0), outs='M502Out', T=20+273.15)

    S501 = u.ReverseOsmosis('S501', ins=M501-0, outs=('RO_treated_water1', effluent1))
    # effluent2 to neutralization and then to biological treatment

    S502 = bst.Splitter('S502', ins=M502-0, outs=('SaltWaste', effluent2), split=[0.99,0.99,0.99], order=['H2O','Na2SO4','NaHSO4'])

    S503 = u.ReverseOsmosis('S503', ins=S502-0, outs=('RO_treated_water2', effluent3))

    H406 = bst.HXutility(
        'H406',
        ins=S408-1,
        outs=LegH_3,
        T=0+273.15,  # Cool to 0°C
        cool_only=True,
    )

    # # ### Facilities ###
    CT = bst.CoolingTower(500 if use_area_convention else 'CT')
    CWP = bst.ChilledWaterPackage(500 if use_area_convention else 'CWP')

    #ADP = bst.AirDistributionPackage(500 if use_area_convention else 'ADP')

    BT = bst.BoilerTurbogenerator(400 if use_area_convention else 'BT',
        (S403-0, 'gas_to_boiler', 'boiler_makeup_water', 'natural_gas', 'lime_boiler', 'boiler_chems'),
        outs=('emissions', 'rejected_water_and_blowdown', 'ash_disposal'),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85,
        satisfy_system_electricity_demand=False,
    )

    makeup_water_streams = (F.cooling_tower_makeup_water,
                            F.Water1, F.Water2,
                            F.Water3, F.Water4,
                            F.Water5, F.Water6,
                            F.Water7, F.Water8,
                            F.boiler_makeup_water)
    process_water_streams = (F.S501.outs[1],
                            F.EvaporatedWater,
                            *makeup_water_streams)

    makeup_water = bst.Stream('makeup_water', price=0.000254)

    PWC = bst.ProcessWaterCenter(500 if use_area_convention else 'PWC',
        ins=('recycled_RO_water', makeup_water, 'recycled_process_water', 'makeup_process_water'),
        outs=('RO_water', 'process_water', 'excess_water'),
        makeup_water_streams=makeup_water_streams,
        process_water_streams=process_water_streams,
        reverse_osmosis_water_price=0.000254,  # USD/kg
        process_water_price=0.000135,  # USD/kg
    )
    # HXN = bst.HeatExchangerNetwork(600 if use_area_convention else 'HXN')
    # load_process_settings()  # Load process settings to update prices and CFs
    s.update_all_input_stream_prices(streamlist=[SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt, DfUltraBuffer, IXEquilibriumBuffer, IXElutionBuffer, IXRegenerationSolution, DfNanoBuffer])

    return LegH_3, vent1, vent2, effluent1, effluent2,effluent3, DfUltraBuffer, IXEquilibriumBuffer, IXElutionBuffer, IXRegenerationSolution, DfNanoBuffer

# %%
if __name__ == '__main__':
    # # Create the LegH system
    bst.preferences.N=50
    LegH_sys = create_LegH_system()
    sys = LegH_sys
    f = sys.flowsheet
    u = f.unit
    ss = f.stream
    sys.operating_hours = 8000
    LegH_sys.simulate()
    LegH_sys.show()
    LegH_sys.diagram(format='html',display=True,)
    # %%
    ss.LegH_3.imass['Leghemoglobin']/ss.LegH_3.F_mass*100
    # %%
    r1 = bst.report.lca_inventory_table(
        systems=[sys],
        key='GWP',
        items=[ss.LegH_3], # For including products without characterization factors
    )
    #r1.to_excel('LegH_LCA_inventory.xlsx',index=True)
    # %%
    r2 = bst.report.lca_displacement_allocation_table(
        systems=[sys],
        key='GWP',
        items=[ss.LegH_3], # For dividing yearly impact by ethanol production
    )
    
    #r2.to_excel('LegH_LCA_displacement_allocation.xlsx',index=True)

# %%
