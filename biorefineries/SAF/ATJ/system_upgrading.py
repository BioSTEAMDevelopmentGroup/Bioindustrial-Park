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
def SAF_sys(ins,outs):
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
    
#%% System setup and process groups
sys = SAF_sys()
