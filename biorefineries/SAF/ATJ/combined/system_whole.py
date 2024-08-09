#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 10:51:39 2024

@author: wenjun
"""

import numpy as np
import thermosteam as tmo
import biosteam as bst
from biosteam import Flowsheet as F 
from biosteam import main_flowsheet
from biosteam import units, Stream, SystemFactory
from biosteam.process_tools import UnitGroup
from biorefineries.SAF._chemicals import chems
from biorefineries.SAF import _units
from biorefineries.cellulosic import units
from biorefineries.SAF.utils import convert_ethanol_wt_2_mol
from biorefineries.SAF._process_settings import add_utility_agent, price
from biorefineries.SAF.ATJ.system_ethanol import F, create_cellulosic_ethanol_system, create_cellulosic_ethanol_chemicals, swg
# from biorefineries.cellulosic.chemicals import create_cellulosic_ethanol_chemicals
# from biorefineries.cellulosic.systems import create_cellulosic_ethanol_system
 
add_utility_agent()

bst.main_flowsheet.set_flowsheet(F)



__all__ = ('sys')

#%%

@SystemFactory(
    ID='upgrading_sys',
    ins=[dict(ID='NaOH'),
         dict(ID='Syndol_catalyst'),
         dict(ID='first_catalyst'),
         dict(ID='second_catalyst'),
         dict(ID='ethanol'),
         dict(ID='hydrogen'),
         dict(ID='Como_catalyst'),],
    outs=[dict(ID='F401_to_WWT'),
          dict(ID='D401_heavy_impurities'),
          dict(ID='D402_top_product'),
          dict(ID='CH4_C2H6'),
          dict(ID='gasoline'),
          dict(ID='jet_fuel'),
          dict(ID='diesel'),
          dict(ID='spent_catalyst_R401'),
          dict(ID='spent_catalyst_R402'),
          dict(ID='spent_catalyst_R403'),
          dict(ID='spent_catalyst_R404'),])
          
def upgrading_sys(ins,outs):
    NaOH,Syndol_catalyst,first_catalyst,second_catalyst,ethanol,hydrogen,Como_catalyst = ins
    F401_1_to_WWT,D401_1_heavy_impurities,D402_1_top_product,CH4_C2H6,gasoline,jet_fuel,diesel,spent_catalyst_R401_1,spent_catalyst_R402_1,spent_catalyst_R403_1,spent_catalyst_R404_1 = outs
    
    P401_1 = bst.Pump('P401_1',ins=ethanol,P=4.5*101325)
    M400_1 = bst.Mixer('M400_1', (P401_1-0,''))
    H400_1 = bst.HXutility('H400',ins=M400_1-0,V=1,rigorous=True)

    R401_1 = _units.AdiabaticFixedbedDehydrationReactor('R401_1', ins=(H400_1-0,Syndol_catalyst),outs=('',spent_catalyst_R401_1))

    # Depressurize to 1 bar before quenching
    V401_1 = bst.IsenthalpicValve('V401_1', ins=R401_1-0,P=101325)

    # Simulate as a quench tower
    H402_1 = bst.HXutility('H402_1', ins=V401_1-0,T=40+273.15,outs='crude_ethylene',rigorous=True)
    S401_1 = bst.PhaseSplitter('S401_1', ins=H402_1-0,outs=('vapor','liquid'))

    # Recover ethylene from water in S401-1 and S402-1
    M401_1 = bst.Mixer('M401_1',ins=(S401_1-1,''),rigorous=False)
    
    # Split flash setting based on eq.plot_vle_binary_phase_envelope(['Ethylene', 'Water'], P=101325) NO VACCUM (T cannot meet)
    F401_1 = bst.SplitFlash('F401_1', ins=M401_1-0, outs=('recovered_ethylene',F401_1_to_WWT), 
                          T=230, P=101325, split=dict(Ethylene=0.999,
                                                      Water=0.0001))
    # Only consider ethylene and water in recovered_ethylene, ignore impurities to simplify
    M402_1 = bst.Mixer('M402_1',ins=(S401_1-0,F401_1-0),outs='')

    # Reduce temperature to 15 C
    # H403 = bst.HXutility('H403', ins=M402-0,outs='',T=15+273.15)

    # Pressurize to 22 bar before purification;
    # Based on Bioethylene Production from Ethanol: A Review and Techno-economical Evaluation
    C401_1 = bst.MultistageCompressor('C401_1', ins=M402_1-0,n_stages=3,pr=2.8,vle=True)

    # Condense water
    S402_1 = bst.PhaseSplitter('S402_1', ins=C401_1-0,outs=('vapor','liquid'))


    # Ethylene recyle from water loop 
    S402_1-1-1-M401_1

    # Remove CO2
    U402_1 = _units.CausticTower('U402_1', ins=(S402_1-0, NaOH),P=22*101325)

    # Remove water and ethanol
    U403_1 = bst.MolecularSieve('U403_1', ins=U402_1-0,outs=('water_and_ethanol','dried_ethylene'),
                              split=dict(Water=1,
                                         Ethanol=1))
    U403_1.outs[0].vle(T=U403_1.outs[0].T,P=U403_1.outs[0].P)
    U403_1-0-1-M400_1
   
    # Reduce temperature before cryogenic distillation low temperature + high pressure
    # Temperature is based on Bioethylene Production from Ethanol: A Review and Techno-economical Evaluation
    H404_1 = bst.HXutility('H404',ins=U403_1-1,T=-25+273.15)

    # Ethylene column by cryogenic distillation, remove heavy keys (Propylene,
    # Butadiene, Ethane, Diethyl ether, and Acetaldehyde) 
    D401_1 = bst.BinaryDistillation('D401_1', ins=H404_1-0, 
                                  outs=('D401_1_top_product', D401_1_heavy_impurities),
                                  LHK=('Ethylene', 'Ethane'),
                                  P=22*101325, 
                                  Lr=0.9999, Hr=0.9999, 
                                  k=1.2)
    # Stripper column by cryogenic distillation, remove light keys (H2, CarbonMonoxide, CH4) 
    # M403 = bst.Mixer('M403',ins=(D401-0,''))
    D402_1 = bst.BinaryDistillation('D402_1', ins=D401_1-0, 
                                  outs=(D402_1_top_product, 'D402_1_pure_ethylene'),
                                  LHK=('CH4', 'Ethylene'), 
                                  P=22*101325,
                                  Lr=0.99, Hr=0.999,
                                  k=1.2)
    
    H405_1 = bst.HXutility('H405_1', ins=D402_1-1, T=10+273.15)

    # Ethylene purity is 100%, ignore impurity for simplification but considered in LCA
    S403_1 = bst.Splitter('S403_1', ins=H405_1-0, outs=('pure_ethylene',CH4_C2H6),
                        split=dict(Ethylene=1))

    # First oligomerization
    R402_1 = _units.Oligomerization1_Reactor('R402_1', ins=(S403_1-0,first_catalyst),outs=('',spent_catalyst_R402_1))

    # Second oligomerization
    R403_1 = _units.Oligomerization2_Reactor('R403_1',ins=(R402_1-0,'',second_catalyst),outs=('',spent_catalyst_R403_1))

    # Recycle light olefins to R403
    D403_1 = bst.BinaryDistillation('D403_1', ins=R403_1-0,outs=('light_olefins','heavy_olefins'),
                                  LHK=('C7H14','C8H16'),
                                  Lr=0.9999,
                                  Hr=0.9999,
                                  k=1.2)
    D403_1-0-1-R403_1

    # Hydrogenation
    P402_1 = bst.Pump('P402_1', ins=D403_1-1, outs='')
    R404_1 = _units.HydrogenationReactor('R404_1', ins=(P402_1-0,hydrogen,Como_catalyst),outs=(spent_catalyst_R404_1,''))
    
    @R404_1.add_specification(run=True)
    def correct_hydrogen_flow():
        feed=R404_1.ins[0]
        h2=R404_1.ins[1]
        h2.imol['H2']=feed.F_mol

    # Fractionation

    D404_1 = bst.BinaryDistillation('D404_1', ins=R404_1-1,outs=('C6_C8_distillate','more_than_C9_bottoms'),
                                  LHK=('C8H18','C9H20'),
                                  Lr=0.9999,
                                  Hr=0.9999,
                                  k=2,
                                  is_divided=False)

    D405_1 = bst.BinaryDistillation('D405_1', ins=D404_1-1,outs=('C9_C16_distillate','more_than_C16_bottoms'),
                                  LHK=('C16H34','C18H38'),
                                  Lr=0.999,
                                  Hr=0.999,
                                  k=1.2,
                                  is_divided=False)

    # Standard temperature for storing gasoline is around 15 C
    H406_1 = bst.HXutility('H406_1', ins=D404_1-0,outs='gasoline1',T=15+273.15,rigorous=True)

    # Standard temperature for storing jet fuel is around 15 C
    H407_1 = bst.HXutility('H407_1', ins=D405_1-0,outs='jet_fuel1',T=15+273.15,rigorous=True)

    # Standard temperature for storing diesel is around 20 C
    H408_1 = bst.HXutility('H408_1', ins=D405_1-1,outs='diesel1',T=20+273.15,rigorous=True)
    
    gasoline_storage = bst.StorageTank('gasoline_storage', ins=H406_1-0, outs=gasoline, tau=7*24, vessel_type='Floating roof', vessel_material='Carbon steel')
    
    jet_storage = bst.StorageTank('jet_storage', ins=H407_1-0, outs=jet_fuel, tau=7*24, vessel_type='Floating roof', vessel_material='Carbon steel')
    
    diesel_storage = bst.StorageTank('diesel_storage', ins=H408_1-0, outs=diesel, tau=7*24, vessel_type='Floating roof', vessel_material='Carbon steel')
    

#%%

# Complete energycane to SAF system
@SystemFactory(ID='SAF_sys')
def SAF_sys(ins,outs,material_storage,product_storage,WWTC,BoilerTurbo,hydrogenation_distillation,h2_purchase,opportunity_cost):
    #%%
    chem = create_cellulosic_ethanol_chemicals()
    chem.set_synonym('Extract','Extractives')
    bst.settings.set_thermo(chem, cache= True)
    
    sys_ethanol = create_cellulosic_ethanol_system('sys_switchgrass',ins = swg)

    F.M701.denaturant_fraction = 0. # 0 for SAF
    sys_ethanol.simulate()
    
    tmo.settings.set_thermo(chems)
    eth_2 = tmo.Stream('eth_2')
    J1 = bst.units.Junction('J1', sys_ethanol-0, tmo.Stream(''))
    J1.prioritize = False
    J1.simulate()
    
    #%%
    SAF_sys = upgrading_sys(ins = (Stream(ID='NaOH',
                                          NaOH=0.5,
                                          Water=0.5,
                                          units='kg/hr',
                                          price=price['NaOH']),
                                   Stream(ID='Syndol_catalyst',
                                          Ash=1,
                                          phase='s',
                                          units='kg/hr',
                                          price=price['Syndol catalyst']),      
                                   Stream(ID='first_catalyst',
                                          Ash=1,
                                          phase='s',
                                          units='kg/hr',
                                          price=price['Ni-loaded aluminosilicate catalyst']),     
                                   Stream(ID='second_catalyst',
                                          Ash=1,
                                          phase='s',
                                          units='kg/hr',
                                          price=price['Aluminosilicate catalyst']),    
                                   J1-0,
                                   Stream(ID='hydrogen',
                                          phase='g',
                                          units='kg/hr',
                                          price=price['h2']),
                                   Stream('Como_catalyst',
                                          Ash=1,
                                          phase='s',
                                          units='kg/hr',
                                          price=price['Como catalyst'])),
                            outs = (Stream(ID='F401_to_WWT'),
                                    Stream(ID='D401_heavy_impurities'),
                                    Stream(ID='D402_top_product'),
                                    Stream(ID='CH4_C2H6'),
                                    Stream(ID='gasoline',
                                           price=price['gasoline']),
                                    Stream(ID='jet_fuel',
                                           price=price['jet fuel']),
                                    Stream(ID='diesel',
                                           price=price['diesel']),
                                    Stream(ID='spent_catalyst_R401',
                                           Ash=1),
                                    Stream(ID='spent_catalyst_R402',
                                           Ash=1),
                                    Stream(ID='spent_catalyst_R403',
                                           Ash=1),
                                    Stream(ID='spent_catalyst_R404',
                                           Ash=1)
                                    ))
                          
                                  
    #==============================================================================
    #                      Area 500 WWT
    #==============================================================================
    M501_1 = bst.Mixer(ID='M501',
                     ins=(
                          # F.pretreatment_to_WWT,
                          # F.U302_to_WWT,
                          # F.D302_to_WWT,
                          F.F401_to_WWT))
    WWT801 = bst.create_conventional_wastewater_treatment_system(ID='WWT801',
                                                              ins=M501_1.outs[0],
                                                              outs=('biogas_from_wastewater_treatment',
                                                                    'sludge_from_wastewater_treatment',
                                                                    'RO_treated_water_from_wastewater_treatment',
                                                                    'brine_from_wastewater_treatment'), # brine goes to waste
                                                              autopopulate=False,
                                                              NaOH_price=price['NaOH'])
    #==============================================================================
    #                      Area 600 Facilities
    #==============================================================================
    # Scale ADP, CIP, FWT by flow 

    # Streams
    natural_gas_1 = Stream('natural_gas_1',price=price['natural gas'])
    lime_boiler_1 = Stream('lime_boiler_1', price=price['lime'])
    boiler_chems_1 = Stream('boiler_chems_1', price=price['boiler chems'])
    ash_1 = Stream('ash_1', price=price['ash disposal'])
    cooling_tower_chems_1 = Stream('cooling_tower_chems_1', price=price['cooling tower chems'])
    system_makeup_water_1 = Stream('system_makeup_water_1', price=price['water'])

    plant_air_in_1 = bst.Stream('plant_air_in', phase='g', N2=0.79, O2=0.21, units='kg/hr')
    CIP_chems_in_1 = Stream('CIP_chems_in', units='kg/hr')
    fire_water_in_1 = Stream('fire_water_in', units='kg/hr')
    
    # Stream mixer to boiler turbogenerator for burning
    M901_1 = bst.Mixer('M901_1',
                     ins=(
                          # F.bagasse_to_CHP,
                          # F.U302_cell_mass,
                          F.D401_heavy_impurities,
                          F.sludge_from_wastewater_treatment), # WWT.outs[1]
                     outs='total_effluent_to_be_burned')  
    
    M902_1 = bst.Mixer('M902_1',
                    ins=(F.D402_top_product,
                         F.biogas_from_wastewater_treatment), # WWT.outs[0]
                    outs='total_gas_to_be_burned')

    BT_1 = bst.BoilerTurbogenerator('BT_1',
                                  ins=(F.total_effluent_to_be_burned, # M901 outs[0]
                                       F.total_gas_to_be_burned, # M902 outs[0]
                                       'boiler_makeup_water_1',
                                       natural_gas_1,
                                       lime_boiler_1,
                                       boiler_chems_1),
                                  outs=('gas_emission',
                                        'boiler_blowdown_water',
                                         ash_1),
                                  satisfy_system_electricity_demand=True,
                                  natural_gas_price=0, # price separately set by the stream
                                  ash_disposal_price=0, # price separately set by the stream
                                  ) # price separately set by the stream
                                  
    BT_1.register_alias('CHP_1')

    # Chilled water package for cooling requirements
    # CWP_1 = bst.ChilledWaterPackage('CWP_1')
                                     
    CT_1 = bst.CoolingTower('CT_1')
    CT_1.ins[-1].price = price['cooling tower chems']
                                                         
    # All water used in the system, here only consider water usage
    system_makeup_water_streams = (F.boiler_makeup_water_1,
                                   CT_1.ins[1]) # Second ins of CT (cooling_tower_makeup_water)
                             
    system_process_water_streams = (
                                    # F.warm_process_water,
                                    # F.water_M301,
                                    # F.water_U301
                                    )
                             
    PWC_1 = bst.ProcessWaterCenter('PWC_1',
                                  ins=(F.RO_treated_water_from_wastewater_treatment, # WWT.outs[2],
                                      '',
                                      system_makeup_water_1),
                                  outs=('','process_water', 'discharged_water'),
                                  makeup_water_streams=system_makeup_water_streams,
                                  process_water_streams=system_process_water_streams)

    ADP_1 = bst.AirDistributionPackage('ADP', ins=plant_air_in_1, outs='plant_air_out_1')
    ADP_1.plant_air_over_feedstock = 0.8
    @ADP_1.add_specification(run=True)
    def adjust_plant_air():
        plant_air_in_1.imass['N2'] = F.ethanol.F_mass * ADP_1.plant_air_over_feedstock
        

    CIP_1 = bst.CIPpackage('CIP', ins=CIP_chems_in_1, outs='CIP_chems_out_1')
    CIP_1.CIP_over_feedstock = 0.00121
    @CIP_1.add_specification(run=True)
    def adjust_CIP():
        CIP_chems_in_1.imass['H2O'] = F.ethanol.F_mass * CIP_1.CIP_over_feedstock


    FWT_1 = _units.FireWaterTank('FWT', ins=fire_water_in_1, outs='fire_water_out_1')
    FWT_1.fire_water_over_feedstock = 0.08
    @FWT_1.add_specification(run=True)
    def adjust_fire_water():
        fire_water_in_1.imass['Water'] = F.ethanol.F_mass * FWT_1.fire_water_over_feedstock
        
    
    HXN_1 = bst.HeatExchangerNetwork('HXN_1',cache_network=True)
    
    
    
    
    if material_storage == False:
        @F.DAP_storage.add_specification(run=True)
        def DAP_storage_no_cost():
            F.DAP_storage._design = lambda:0
            F.DAP_storage._cost = lambda:0
        
        @F.CSL_storage.add_specification(run=True)
        def CSL_storage_no_cost():
            F.CSL_storage._design = lambda:0
            F.CSL_storage._cost = lambda:0
    else: 
        pass
    
    if product_storage == False:
        @F.ethanol_storage.add_specification(run=True)
        def ethanol_storage_no_cost():
            F.ethanol_storage._design = lambda:0
            F.ethanol_storage._cost = lambda:0
        @F.gasoline_storage.add_specification(run=True)
        def gasoline_storage_no_cost():
            F.gasoline_storage._design = lambda:0
            F.gasoline_storage._cost = lambda:0
        @F.jet_storage.add_specification(run=True)
        def jet_storage_no_cost():
            F.jet_storage._design = lambda:0
            F.jet_storage._cost = lambda:0
        @F.diesel_storage.add_specification(run=True)
        def diesel_storage_no_cost():
            F.diesel_storage._design = lambda:0
            F.diesel_storage._cost = lambda:0
    else:
        pass
    
    
    
    if hydrogenation_distillation == False:
        # R404
        R404_cost =F.R404._cost
        R404_design = F.R404._design
        def R404_no_cost():
            R404_cost()
            bpc = F.R404.baseline_purchase_costs
            pc = F.R404.purchase_costs
            ic = F.R404.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0. 
        def R404_no_cost_design():
            R404_design()
            bpc = F.R404.baseline_purchase_costs
            pc = F.R404.purchase_costs
            ic = F.R404.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0.
        F.R404._cost = R404_no_cost
        F.R404._design = R404_no_cost_design 
        def R404_summary(design_kwargs=None, cost_kwargs=None, lca_kwargs=None):
            F.R404._check_run()
            F.R404._design()
            F.R404._cost()
            for i in F.R404.auxiliary_units:
                F.R404.heat_utilities += i.heat_utilities
                F.R404.power_utility.rate += i.power_utility.rate
            F.R404._lca()
            F.R404._load_operation_costs()
        F.R404._summary = R404_summary
        
        # D404
        D404_cost =F.D404._cost
        D404_design = F.D404._design
        def D404_no_cost():
            D404_cost()
            bpc = F.D404.baseline_purchase_costs
            pc = F.D404.purchase_costs
            ic = F.D404.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0. 
        def D404_no_cost_design():
            D404_design()
            bpc = F.D404.baseline_purchase_costs
            pc = F.D404.purchase_costs
            ic = F.D404.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0.
        F.D404._cost = D404_no_cost
        F.D404._design = D404_no_cost_design 
        def D404_summary(design_kwargs=None, cost_kwargs=None, lca_kwargs=None):
            F.D404._check_run()
            F.D404._design()
            F.D404._cost()
            for i in F.D404.auxiliary_units:
                F.D404.heat_utilities += i.heat_utilities
                F.D404.power_utility.rate += i.power_utility.rate
            F.D404._lca()
            F.D404._load_operation_costs()
        F.D404._summary = D404_summary
        
        # D405
        D405_cost =F.D405._cost
        D405_design = F.D405._design
        def D405_no_cost():
            D405_cost()
            bpc = F.D405.baseline_purchase_costs
            pc = F.D405.purchase_costs
            ic = F.D405.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0. 
        def D405_no_cost_design():
            D405_design()
            bpc = F.D405.baseline_purchase_costs
            pc = F.D405.purchase_costs
            ic = F.D405.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0.
        F.D405._cost = D405_no_cost
        F.D405._design = D405_no_cost_design 
        def D405_summary(design_kwargs=None, cost_kwargs=None, lca_kwargs=None):
            F.D405._check_run()
            F.D405._design()
            F.D405._cost()
            for i in F.D405.auxiliary_units:
                F.D405.heat_utilities += i.heat_utilities
                F.D405.power_utility.rate += i.power_utility.rate
            F.D405._lca()
            F.D405._load_operation_costs()
        F.D405._summary = D405_summary
        
        # H406
        H406_cost =F.H406._cost
        H406_design = F.H406._design
        def H406_no_cost():
            H406_cost()
            bpc = F.H406.baseline_purchase_costs
            pc = F.H406.purchase_costs
            ic = F.H406.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0. 
        def H406_no_cost_design():
            H406_design()
            bpc = F.H406.baseline_purchase_costs
            pc = F.H406.purchase_costs
            ic = F.H406.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0.
        F.H406._cost = H406_no_cost
        F.H406._design = H406_no_cost_design 
        def H406_summary(design_kwargs=None, cost_kwargs=None, lca_kwargs=None):
            F.H406._check_run()
            F.H406._design()
            F.H406._cost()
            for i in F.H406.auxiliary_units:
                F.H406.heat_utilities += i.heat_utilities
                F.H406.power_utility.rate += i.power_utility.rate
            F.H406._lca()
            F.H406._load_operation_costs()
        F.H406._summary = H406_summary
        
        # H407
        H407_cost =F.H407._cost
        H407_design = F.H407._design
        def H407_no_cost():
            H407_cost()
            bpc = F.H407.baseline_purchase_costs
            pc = F.H407.purchase_costs
            ic = F.H407.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0. 
        def H407_no_cost_design():
            H407_design()
            bpc = F.H407.baseline_purchase_costs
            pc = F.H407.purchase_costs
            ic = F.H407.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0.
        F.H407._cost = H407_no_cost
        F.H407._design = H407_no_cost_design 
        def H407_summary(design_kwargs=None, cost_kwargs=None, lca_kwargs=None):
            F.H407._check_run()
            F.H407._design()
            F.H407._cost()
            for i in F.H407.auxiliary_units:
                F.H407.heat_utilities += i.heat_utilities
                F.H407.power_utility.rate += i.power_utility.rate
            F.H407._lca()
            F.H407._load_operation_costs()
        F.H407._summary = H407_summary
        
        # H408
        H408_cost =F.H408._cost
        H408_design = F.H408._design
        def H408_no_cost():
            H408_cost()
            bpc = F.H408.baseline_purchase_costs
            pc = F.H408.purchase_costs
            ic = F.H408.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0. 
        def H408_no_cost_design():
            H408_design()
            bpc = F.H408.baseline_purchase_costs
            pc = F.H408.purchase_costs
            ic = F.H408.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0.
        F.H408._cost = H408_no_cost
        F.H408._design = H408_no_cost_design 
        def H408_summary(design_kwargs=None, cost_kwargs=None, lca_kwargs=None):
            F.H408._check_run()
            F.H408._design()
            F.H408._cost()
            for i in F.H408.auxiliary_units:
                F.H408.heat_utilities += i.heat_utilities
                F.H408.power_utility.rate += i.power_utility.rate
            F.H408._lca()
            F.H408._load_operation_costs()
        F.H408._summary = H408_summary
    
    else:
        pass


    # Operating cost (heat and steam demand shouldn't be accounted since we have considered the opportunity cost)
    # if hydrogenation_distillation == False:
    #     @F.R404.add_specification(run=True)
    #     def R404_no_cost():
    #         F.R404._design = lambda:0
    #         F.R404._cost = lambda:0
    #     @F.D404.add_specification(run=True)
    #     def D404_no_cost():
    #         F.D404._design = lambda:0
    #         F.D404._cost = lambda:0
    #     @F.D405.add_specification(run=True)
    #     def D405_no_cost():
    #         F.D405._design = lambda:0
    #         F.D405._cost = lambda:0
    #     @F.H406.add_specification(run=True)
    #     def H406_no_cost():
    #         F.H406._design = lambda:0
    #         F.H406._cost = lambda:0
    #     @F.H407.add_specification(run=True)
    #     def H407_no_cost():
    #         F.H407._design = lambda:0
    #         F.H407._cost = lambda:0
    #     @F.H408.add_specification(run=True)
    #     def H408_no_cost():
    #         F.H408._design = lambda:0
    #         F.H408._cost = lambda:0
    # else:
    #     pass
    
    
    
    if WWTC == False:
        WWTC_cost =F.WWTC._cost
        WWTC_design = F.WWTC._design
        def WWTC_no_cost():
            WWTC_cost()
            bpc = F.WWTC.baseline_purchase_costs
            pc = F.WWTC.purchase_costs
            ic = F.WWTC.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0. 
        def WWTC_no_cost_design():
            WWTC_design()
            bpc = F.WWTC.baseline_purchase_costs
            pc = F.WWTC.purchase_costs
            ic = F.WWTC.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0.
        F.WWTC._cost = WWTC_no_cost
        F.WWTC._design = WWTC_no_cost_design 
        def WWTC_summary(design_kwargs=None, cost_kwargs=None, lca_kwargs=None):
            F.WWTC._check_run()
            F.WWTC._design()
            F.WWTC._cost()
            for i in F.WWTC.auxiliary_units:
                F.WWTC.heat_utilities += i.heat_utilities
                F.WWTC.power_utility.rate += i.power_utility.rate
            F.WWTC._lca()
            F.WWTC._load_operation_costs()
        F.WWTC._summary = WWTC_summary
    else:
        pass
    
    
    
    
    
    if h2_purchase == False:
        @F.U405.add_specification(run=True)
        def ng_flow():
            F.R404.run()
            F.U405.ins[0].imol['CH4'] = F.R404.ins[1].imol['H2']/3 # assume CH4+H2O=CO+3H2
    else:
        F.U405.ins[0].imol['CH4'] = 0
            
    
    
    if BoilerTurbo == False:
        BT_cost = F.BT._cost
        BT_design = F.BT._design
        def BT_no_cost():
            BT_cost()
            bpc = F.BT.baseline_purchase_costs
            pc = F.BT.purchase_costs
            ic = F.BT.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0. 
        def BT_no_cost_design():
            BT_design()
            bpc = F.BT.baseline_purchase_costs
            pc = F.BT.purchase_costs
            ic = F.BT.installed_costs
            for i in bpc.keys(): bpc[i] = 0. 
            for i in pc.keys(): pc[i] = 0. 
            for i in ic.keys(): ic[i] = 0.
        F.BT._design = BT_no_cost_design
        F.BT._cost = BT_no_cost
        def BT_summary(design_kwargs=None, cost_kwargs=None, lca_kwargs=None):
            F.BT._check_run()
            F.BT._design()
            F.BT._cost()
            for i in F.BT.auxiliary_units:
                F.BT.heat_utilities += i.heat_utilities
                F.BT.power_utility.rate += i.power_utility.rate
            F.BT._lca()
            F.BT._load_operation_costs()
        F.BT._summary = BT_summary
    else:
        pass
    
    
    if opportunity_cost == True:
        @F.U404.add_specification(run=True)
        def oc_flow():
            oc = F.U404.ins[0]
            oc.reset_flow(Water=1, phase='l', units='kg/hr', total_flow=1.07*F.R404.ins[0].F_vol*845) # 1.07 is volume gain from https://www.eia.gov/energyexplained/oil-and-petroleum-products/refining-crude-oil-inputs-and-outputs.php; assume 845 is density but will offset when calculating toal price
    else:
        pass
    