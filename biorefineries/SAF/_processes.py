#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 13:13:15 2023

@author: wenjun
"""
#%%
'''
Naming conventions:
    D = Distillation column
    F = Flash tank
    H = Heat exchanger
    M = Mixer
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including PhaseSplitter)
    T = Tank or bin for storage
    U = Other units
    PS = Process specifications, not physical units, but for adjusting streams

Processes:
    100: Preprocessing
    200: Pretreatment
    300: Fermentation
    400: Upgrading
    500: Wastewater
    600: Facilities
'''


#%% Set up

import numpy as np
import biosteam as bst
import thermosteam as tmo
from thermosteam import Stream, reaction as rxn
from biosteam import System, SystemFactory
from biorefineries.SAF import _chemicals as chems
from biorefineries.SAF.utils import convert_ethanol_wt_2_mol
from biorefineries.SAF._process_settings import price
from biorefineries.SAF import _units 
from biorefineries.cellulosic import units

F=bst.Flowsheet('energycane_SAF')
bst.main_flowsheet.set_flowsheet(F)
u = F.unit 
s = F.stream
#%%
#==============================================================================
#                       Area 100 Preprocessing process
#==============================================================================

# Streams
energycane = Stream('energycane',
                    Water=0.6,
                    Sucrose=0.077,
                    Glucose=0.007,
                    Fructose=0.006,
                    Ash=0.029,
                    Glucan=0.129,
                    Xylan=0.07,
                    Arabinan=0.008,
                    Lignin=0.071,
                    Extract=0.004,
                    total_flow=333333.33,
                    units='kg/hr',
                    price=price['Feedstock'])
# Assume processing capacity of wet basis=1,600,000 metric tons (MT)/year, with 200 operating days,8000 MT/day,operating hours of 24 hr/day.
# Content is based on 'Techno-economic feasibility analysis of engineered energycane-based biorefinery co-producing biodiesel and ethanol'. 
# Bagasse content is based on 'Enzyme hydrolysis and ethanol fermentation of dilute ammonia pretreated energy cane'.

water_for_imbibition = Stream('water_for_imbibition',Water=90000,units='kg/hr',T=25+273.15,price=price['Water']) # From 'Techno-economic feasibility analysis of engineered energycane-based biorefinery co-producing biodiesel and ethanol'.

imbibition_water_H = bst.HXutility('water_for_imbibition', 
                               ins=water_for_imbibition,
                               outs='imbibition_water', 
                               T=60+273.15)

imbibition_water = imbibition_water_H.outs[0]

water_for_rvf = Stream('water_for_rvf',Water=1000,units='kg/hr',T=25+273.15,price=price['Water'])

rvf_wash_water_H = bst.HXutility('rvf_wash_water_H',
                                 ins=water_for_rvf,
                                 outs='rvf_wash_water',
                                 T=90+273.15)  

rvf_wash_water = rvf_wash_water_H.outs[0]  
               
H3PO4 = Stream('H3PO4',
               H3PO4=74.23,
               Water=13.1,
               units='kg/hr',
               price=price['H3PO4'])
                   
lime = Stream('lime',
              CaO=333,
              Water=2200,
              units='kg/hr', 
              price=price['Lime'])
                  
polymer = Stream('polymer',
                     Flocculant=0.83, 
                     units='kg/hr',
                     price=price['Flocculant'])

U101 = bst.ConveyingBelt('U101', ins=energycane)
U102 = bst.MagneticSeparator('U102', ins=U101-0, outs=('A','B'))
U103 = bst.Shredder('U103', ins=U102-0,outs='')
U104 = bst.CrushingMill('U104', ins=(U103-0,''),
                        split=dict(Ash=0.92,
                                   Glucan=0.92,
                                   Xylan=0.92,
                                   Arabinan=0.92,
                                   Lignin=0.92,
                                   Sucrose=0.04,
                                   Glucose=0.04,
                                   Fructose=0.04,
                                   Extract=0.04,
                                   Solids=1), 
                        moisture_content=0.5)
# Correct imbibition water
@U101.add_specification(run=True)
def correct_imbibition_water():
    feed = U101.ins[0]
    # moisture = feed.imass['Water']
    # dry = feed.F_mass-moisture
    water_for_imbibition.imass['Water'] = 0.25 * feed.F_mass

# Mix in water
M101 = bst.Mixer('M101',ins=('',imbibition_water),outs=1-U104)

# Screen out fibers
S101 = bst.VibratingScreen('S101',ins=U104-1,outs=('', 0-M101),
                           split=dict(Ash=0.35,
                                      Glucan=0.35,
                                      Xylan=0.35,
                                      Arabinan=0.35,
                                      Lignin=0.35,
                                      Solids=0,
                                      Sucrose=0.88,
                                      Glucose=0.88,
                                      Fructose=0.88,
                                      Extract=0.88,
                                      Water=0.88))

# Store juice
T101 = bst.StorageTank('T101',ins=S101-0,outs='untreated_juice',
                       tau=4,vessel_material='Carbon steel')

U105 = bst.ConveyingBelt('U105',ins=U104-0,outs='bagasse')

# Split bagasse
S102 = bst.Splitter('S102',ins=U105-0,outs=['bagasse_to_ethanol','bagasse_to_CHP'],split=0.8)

# Heat up juice #
H101 = bst.HXutility('H101',ins=T101-0,T=343.15)
        
# Mix with acid #
T102 = bst.MixTank('T102',ins=[H101-0,H3PO4])
        
# Pump solution #
P101 = bst.Pump('P101',ins=T102-0)
        
# Mix lime solution #
T104 = bst.MixTank('T104',ins=(P101-0,lime),tau=0.1)
P102 = bst.Pump('P102',ins=T104-0)
        
# Mix recycle with rvf #
M103 = bst.Mixer('M103',ins=(P102-0, ''))
        
# Heat #
H102 = bst.HXutility('H102',M103-0,T=105+273.15)
        
# Mix in flocculant polymer #
T105 = bst.MixTank('T105',ins=(H102-0,polymer),tau=0.1)
        
# Clarify #
C101 = bst.Clarifier('C101',ins=T105-0,outs=('clarified_juice',''),
                     split=dict(Ash=0,
                                CaO=0,
                                Glucan=0,
                                Xylan=0,
                                Arabinan=0,
                                Lignin=0,
                                Flocculant=0.522,
                                H3PO4=0.522,
                                Glucose=0.522,
                                Sucrose=0.522,
                                Fructose=0.522,
                                Extract=0.522,
                                Water=0.522))
# Remove solids as filter cake #
C102 = bst.RVF('C102',ins=(C101-1,rvf_wash_water),outs=('filter_cake',''),
                       moisture_content=0.80,
                       split=dict(Ash=0.85,
                                  CaO=0.85,
                                  Glucan=0.85,
                                  Xylan=0.85,
                                  Arabinan=0.85,
                                  Lignin=0.85,
                                  Glucose=0.01,
                                  Sucrose=0.01,
                                  Fructose=0.01,
                                  Extract=0.01))

P103 = bst.Pump('P103',ins=C102-1,outs=1-M103)
        
# Specifications dependent on energy cane flow rate
@U104.add_specification(run=True) 
def correct_imbibition_water():
    feed = U104.ins[0]
    lime.imass['CaO', 'Water'] = 0.001 * feed.F_mass * np.array([0.046, 0.954])
    H3PO4.imass['H3PO4', 'Water'] = 0.00025 * feed.F_mass
    
# Specifications within P102
@P102.add_specification(run=True)
def correct_wash_water():
    P102._run()
    solids = P102.outs[0].imol['Ash', 'CaO', 'Glucan',
                                'Xylan','Arabinan', 'Lignin'].sum()
    rvf_wash_water.imol['Water'] = 0.0574 * solids
           
# Juice screening process
S103 = bst.VibratingScreen('S103',ins=C101-0,outs=('screened_juice','fiber_fines'),
                           split=dict(Ash=0.998,
                                      CaO=0.998,
                                      H3PO4=1.0,
                                      Flocculant=0.0,
                                      Glucan=0.0,
                                      Xylan=0.0,
                                      Arabinan=0.0,
                                      Lignin=0.0,
                                      Glucose=0.998,
                                      Fructose=0.998,
                                      Sucrose=0.998,
                                      Extract=0.998,
                                      Water=0.998))
S103.mesh_opening=2

H103 = bst.HXutility('H103',ins=S103-0,T=25+273.15)   
   
# Juice concentrating process 
S104 = bst.Splitter('S104',ins=H103-0,split=0.9)
E101 = bst.MultiEffectEvaporator('E101', ins=S104-0,outs=('solids','condensate'),
                                 V=0.1, V_definition='First-effect',
                                 P=(101325, 73581, 50892, 32777))
M104 = bst.Mixer('M104', ins=(E101-0,S104-1),outs='')
H104 = bst.HXutility('H104', ins=M104-0, outs='concentrated_juice', T=32+273.15)
    
E101.target_sugar_concentration=0.23
@E101.add_bounded_numerical_specification(x0=0,x1=1,xtol=1e-5,ytol=1e-2)
def sugar_concentration_at_fraction_evaporation(V):
    E101.V=V
    E101.run_until(M104,inclusive=True)
    Glucose_concentration = M104.outs[0].get_mass_fraction('Glucose')
    Sucrose_concentration = M104.outs[0].get_mass_fraction('Sucrose')
    Fructose_concentration = M104.outs[0].get_mass_fraction('Fructose')
    Extract_concentration = M104.outs[0].get_mass_fraction('Extract')
    sugar_concentration = Glucose_concentration + Sucrose_concentration + Fructose_concentration + Extract_concentration
    return E101.target_sugar_concentration - sugar_concentration

# # create system object
# sys_preprocessing = bst.main_flowsheet.create_system('sys_preprocessing')

# # simulate system
# sys_preprocessing.simulate()

# sys_preprocessing.diagram(format='png')

# for i in sys_preprocessing.units: i.show('cwt100')

   

#%%
#==============================================================================
#                        Area 200 Pretreatment process
#==============================================================================

solids_loading = 0.305
nonsolids = ['Water']
T_pretreatment_reactor = 130.+273.15

# Heat process water
water_for_pretreatment = Stream('water_for_pretreatment', 
                            T=25+273.15,
                            Water=1,
                            price=price['Water'])
warm_process_water_H = bst.HXutility('warm_process_water_H',
                                     ins=water_for_pretreatment,
                                     outs='',
                                     T=368.15,)
warm_process_water_P = bst.Pump('warm_process_water_P',
                                ins=warm_process_water_H-0,
                                outs='warm_process_water',
                                P=4.7*101325)

# Streams
warm_process_water = warm_process_water_P.outs[0]

pretreatment_steam = Stream('pretreatment_steam',
                            phase='g',
                            T=268+273.15,
                            P=13*101325,
                            Water=24534+3490,
                            units='kg/hr')

P = pretreatment_steam.chemicals['H2O'].Psat(T_pretreatment_reactor + 25)
M201 = bst.SteamMixer('M201',ins=(S102-0,pretreatment_steam,warm_process_water),
                      P=P,solids_loading=solids_loading)
R201 = _units.PretreatmentReactor('R201',ins=M201-0,T=T_pretreatment_reactor)
P201 = units.BlowdownDischargePump('P201',ins=R201-1)
F201 = units.PretreatmentFlash('F201', ins=P201-0,P=101325,Q=0)
M202 = bst.Mixer('M202',ins=(R201-0,F201-0))
units.WasteVaporCondenser('H201',ins=M202-0,outs='pretreatment_to_WWT',V=0)

P202 = units.HydrolyzatePump('P202',ins=F201-1)

U201 = bst.HammerMill('U201', ins=P202-0,outs='pretreated_bagasse')

# sys_pretreatment=bst.main_flowsheet.create_system('sys_pretreatment')
# sys_pretreatment.simulate()

#%%
#==============================================================================
#                        Area 300 Fermentation process
#==============================================================================

# Stream
enzyme_M301 = Stream('enzyme_M301',units='kg/hr',price=price['Enzyme'])
water_M301 = Stream('water_M301',units='kg/hr',price=price['Water'])

CSL_R301 = Stream('CSL_R301',units='kg/hr',price=price['CSL'])
CSL_R302 = Stream('CSL_R302',units='kg/hr',price=price['CSL'])

DAP_R301 = Stream('DAP_R301',units='kg/hr',price=price['DAP'])
DAP_R302 = Stream('DAP_R302',units='kg/hr',price=price['DAP'])

water_U301 = Stream('water_U301',units='kg/hr',price=price['Water'])

# SSCF
nonsolids=chems.default_nonsolids
insoluble_solids=chems.default_insoluble_solids
ignored=chems.default_ignored

insoluble_solids_loading=10.3

M301 = _units.EnzymeHydrolysateMixer('M301',ins=(U201-0,enzyme_M301,water_M301),
                                     T=48+273.15,enzyme_loading=20,solids_loading=0.2)
R301 = _units.SaccharificationAndCoFermentation('R301',ins=(M301-0,'',CSL_R301,DAP_R301,H104-0), 
                                                outs=('R301_g','effluent','side_draw'))
R302 = _units.SeedTrain('R302', ins=(R301-2,CSL_R302,DAP_R302),
                       outs=('R302_g','seed'))
T301 = _units.SeedHoldTank('T301', ins=R302-1,outs=1-R301)

M302 = bst.Mixer('M302',ins=(R301-0,R302-0),outs='fermentation_vapor')

U301 = bst.VentScrubber('U301', ins=(water_U301, M302-0), 
                        outs=('U301_vent', 'U301_recycled'),
                        gas=('CO2', 'NH3', 'O2'))

@U301.add_specification(run=True)
def update_U301_water():
    water_U301 = U301.ins[0]
    M302.run_until(U301,inclusive=True)
    water_U301.imass['Water']=26836/21759 * M302.F_mass_in

M303 = bst.Mixer('M303',ins=(R301-1,U301-1),outs='')
T302 = _units.BeerTank('T302', ins=M303-0,outs='')

# Heat up crude beer by exchanging heat with stillage
H302 = bst.HXprocess('H302', ins=(T302-0,''),outs='',
                     phase0='l',phase1='l',U=1.28)

# Remove solids from fermentation broth, based on the pressure filter in ref [1]
# Moisture content is 35% in ref [1] but 25% in ref [2], used 35% to be conservative
U302 = _units.CellMassFilter('U302', ins=H302-1,outs=('U302_cell_mass','U302_to_WWT'),
                             moisture_content=0.35,split=0.99)

# Beer column
xbot = convert_ethanol_wt_2_mol(0.00001) 
ytop = convert_ethanol_wt_2_mol(0.5)

D301 = bst.BinaryDistillation('D301',ins=H302-0,outs=('',''),k=1.25,Rmin=0.6,
                              P=101325,y_top=ytop,x_bot=xbot,
                              LHK=('Ethanol','Water'),
                              tray_material='Stainless steel 304',
                              vessel_material='Stainless steel 304')
D301.reboiler.U = 1.85
D301_P = bst.Pump('D301_P', ins=D301-1, outs=1-H302)
D301_P.BM = 3.1

# Mix recycled ethanol
M304 = bst.Mixer('M304', ins=(D301-0, ''))

ytop = convert_ethanol_wt_2_mol(0.915)
D302 = bst.BinaryDistillation('D302', ins=M304-0,outs=('',''), k=1.25, Rmin=0.6,
                              P=101325, y_top=ytop, x_bot=xbot,
                              LHK=('Ethanol', 'Water'),
                              tray_material='Stainless steel 304',
                              vessel_material='Stainless steel 304',
                              is_divided=True)
                                    
D302.reboiler.U = 1.85
D302_P = bst.Pump('D302_P', ins=D302-1, outs='D302_to_WWT')
D302_P.BM = 3.1
D302_H = bst.HXutility('D302_H', ins=D302-0, outs='', T=115+283.15, V=1)

# Molecular sieve, split based on streams 515 and 511 in ref [1]
split_ethanol = 1 - 21673/27022
split_water = 1 - 108/2164
S301 = bst.MolecularSieve('S301', ins=D302_H-0, outs=(1-M304, ''),
                          split=(split_ethanol, split_water),
                          order=('Ethanol', 'Water'))
                                
# Condense ethanol product
S301_H = bst.HXutility('S301_H', ins=S301-1, outs='ethanol_to_storage',
                       V=0, T=20+273.15)
T303 = bst.StorageTank('T303',S301_H-0,tau=12,
                       vessel_type='Floating roof',
                       vessel_material='Carbon steel')

# Ethanol need not be denatured 

# sys_SSCF=bst.main_flowsheet.create_system('sys_SSCF')
# sys_SSCF.simulate()

#%%
#==============================================================================
#                        Area 400 Upgrading process
#==============================================================================

# Add a heating agent
import qsdsan as qs
DPO_chem = qs.Chemical('DPO_chem', search_ID='101-84-8')
BIP_chem = qs.Chemical('BIP_chem', search_ID='92-52-4')
DPO = qs.Component.from_chemical('DPO', chemical=DPO_chem, particle_size='Soluble',
                                 degradability='Slowly', organic=True)
BIP = qs.Component.from_chemical('BIP', chemical=BIP_chem, particle_size='Soluble',
                                 degradability='Slowly', organic=True)
HTF_thermo = bst.Thermo((DPO, BIP,))
HTF = bst.UtilityAgent('HTF', DPO=0.735, BIP=0.265, T=673.15, P=10.6*101325, phase='g',
                       # 400 C (673.15 K) and 138 psig (951477 pa) are max temp and pressure for HTF
                       thermo=HTF_thermo,
                       # T_limit = 495 F (530.372 K) is the highest temp that vapor can exist
                       regeneration_price=1) # Lang
                       # use default heat transfer efficiency (1)
# Temperature and pressure: https://www.dow.com/content/dam/dcc/documents/\
# en-us/app-tech-guide/176/176-01334-01-dowtherm-heat-transfer-fluids-\
# engineering-manual.pdf?iframe=true (accessed on 11-16-2022)
bst.HeatUtility.heating_agents.append(HTF)
bst.CE = qs.CEPCI_by_year[2020] # use 2020$ to match up with latest PNNL report

# Add a cooling agent
Decamethyltetrasiloxane_chem = qs.Chemical('Decamethyltetrasiloxane_chem',search_ID='141-62-8')
Octamethyltrisiloxane_chem = qs.Chemical('Octamethyltrisiloxane',search_ID='107-51-7')
DEC = qs.Component.from_chemical('DEC', chemical=Decamethyltetrasiloxane_chem, particle_size='Soluble',
                                 degradability='Slowly', organic=True)
OCT = qs.Component.from_chemical('OCT', chemical=Octamethyltrisiloxane_chem, particle_size='Soluble',
                                 degradability='Slowly', organic=True)
LTF_thermo = bst.Thermo((DEC,OCT,))
LTF = bst.UtilityAgent('LTF', DEC=0.4, OCT=0.6, T=173.15, P=5.2*101325, phase='l',
                       thermo = LTF_thermo,
                       regeneration_price=1)
bst.HeatUtility.cooling_agents.append(LTF)


# Dehydration
# ethanol_to_storage = Stream('ethanol_to_storage',
#                             phase='l', T=298.15, P=101325,
#                             H2O=10.3,
#                             Ethanol=694,
#                             units='kmol/hr') 

# Stream
NaOH = Stream('NaOH',NaOH=0.5,Water=0.5,units='kg/hr',price=price['NaOH'])
Syndol_catalyst = Stream('Syndol_catalyst',units='kg/hr',price=price['Syndol catalyst'])

P401 = bst.Pump('P401',ins=T303-0,P=4.5*101325)
M400 = bst.Mixer('M401', (P401-0,''))
H400 = bst.HXutility('H400',ins=P401-0,V=1,rigorous=True)

R401 = _units.AdiabaticFixedbedDehydrationReactor('R401', ins=(H400-0,Syndol_catalyst),outs=('','spent_catalyst'))

# Depressurize to 1 bar before quenching
V401 = bst.IsenthalpicValve('V401', ins=R401-0,P=101325)

# Simulate as a quench tower
H402 = bst.HXutility('H402', ins=V401-0,T=50+273.15,outs='crude_ethylene',rigorous=True)
S401 = bst.PhaseSplitter('S401', ins=H402-0,outs=('vapor','liquid'))

S401.target_ethylene_cwt=0.98
@H402.add_bounded_numerical_specification(x0=10+273.15,x1=80+273.15,xtol=1e-5,ytol=1e-2)
def ethylene_cwt_at_T(T):
    H402.T=T
    H402.run_until(S401,inclusive=True)
    ethylene_cwt=S401.outs[0].get_mass_fraction('Ethylene')
    return S401.target_ethylene_cwt - ethylene_cwt

# Recover ethylene from water in S401-1 and S402-1
M401 = bst.Mixer('M401',ins=(S401-1,''),rigorous=False)
F401 = bst.Flash('F401',ins=M401-0,outs=('recovered_ethylene','F401_to_WWT'),T=353,P=101325)

F401.target_ethylene_recovery_ratio=0.99
@F401.add_bounded_numerical_specification(x0=280,x1=373,xtol=1e-5,ytol=1e-2)
def ethylene_recovery_at_T(T):
    F401.T=T
    F401.run()
    ethylene_recovery=F401.outs[0].imass['Ethylene']
    ethylene_recovery_ratio=ethylene_recovery/F401.ins[0].imass['Ethylene']
    return F401.target_ethylene_recovery_ratio-ethylene_recovery_ratio

# Only consider ethylene and water in recovered_ethylene, ignore impurities to simplify
M402 = bst.Mixer('M402',ins=(S401-0,F401-0),outs='',rigorous=False)

# Reduce temperature to 15 C
# H403 = bst.HXutility('H403', ins=M402-0,outs='',T=15+273.15)

# Pressurize to 22 bar before purification;
# Based on Bioethylene Production from Ethanol: A Review and Techno-economical Evaluation
C401 = bst.MultistageCompressor('C401', ins=M402-0,n_stages=3,pr=2.8,eta=0.72,vle=True)

# Condense water
S402 = bst.PhaseSplitter('S402', ins=C401-0,outs=('vapor','liquid'))


# Ethylene recyle from water loop 
S402-1-1-M401

# Remove CO2
U402 = _units.CausticTower('U402', ins=(S402-0, NaOH),P=22*101325)

# Remove water and ethanol
U403 = bst.MolecularSieve('U403', ins=U402-0,outs=('water_and_ethanol','dried_ethylene'),
                          split=dict(Water=1,
                                     Ethanol=1))
U403.outs[0].vle(T=U403.outs[0].T,P=U403.outs[0].P)
U403-0-1-M400
# Reduce temperature before cryogenic distillation low temperature + high pressure
# Temperature is based on Bioethylene Production from Ethanol: A Review and Techno-economical Evaluation
H404 = bst.HXutility('H404',ins=U403-1,T=-25+273.15)

# Ethylene column by cryogenic distillation, remove heavy keys (Propylene,
# Butadiene, Ethane, Diethyl ether, and Acetaldehyde) 
D401 = bst.BinaryDistillation('D401', ins=H404-0, 
                              outs=('D401_top_product', 'D401_heavy_impurities'),
                              LHK=('Ethylene', 'Ethane'),
                              P=22*101325, 
                              Lr=0.9999, Hr=0.9999, 
                              k=1.2)
# Stripper column by cryogenic distillation, remove light keys (H2, CarbonMonoxide, CH4) 
# M403 = bst.Mixer('M403',ins=(D401-0,''))
D402 = bst.BinaryDistillation('D402', ins=D401-0, 
                              outs=('D402_top_product', 'D402_pure_ethylene'),
                              LHK=('CH4', 'Ethylene'), 
                              P=22*101325,
                              Lr=0.99, Hr=0.999,
                              k=1.2)
# S403 = bst.Splitter('S403',ins=D402-0,outs=('recycled_D402_top_product', 'wasted_D402_top_product'),split=0.5)
# S403-0-1-M403
H405 = bst.HXutility('H405', ins=D402-1, T=10+273.15)

# Ethylene purity is 100%, ignore impurity for simplification
S403 = bst.Splitter('S403', ins=H405-0, outs=('pure_ethylene','CH4_C2H6'),
                    split=dict(Ethylene=1))

# S403.outs[0].as_stream()

# Oligomerization
first_catalyst = Stream('first_catalyst',phase='s',units='kg/hr',price=price['Ni-loaded aluminosilicate catalyst'])
second_catalyst = Stream('second_catalyst',phase='s',units='kg/hr',price=price['Aluminosilicate catalyst'])

# First oligomerization

# Ethylene = Stream('Ethylene',Ethylene=682,units='kmol/hr',phase='g')
R402 = _units.EthyleneDimerizationReactor('R402', ins=(S403-0,first_catalyst),outs=('','spent_catalyst'))

# Second oligomerization
R403 = _units.OligomerizationReactor('R403',ins=(R402-0,'',second_catalyst),outs=('','spent_catalyst'))

# Recycle light olefins to R403
D403 = bst.BinaryDistillation('D403', ins=R403-0,outs=('light_olefins','heavy_olefins'),
                              LHK=('C6H12','C7H14'),
                              Lr=0.9999,
                              Hr=0.9999,
                              k=1.2)
D403-0-1-R403

# Hydrogenation

# Stream
Pd_Al2O3_catalyst = Stream('Pd_Al2O3_catalyst',phase='s',units='kg/hr')
como_catalyst = Stream('como_catalyst',phase='s',units='kg/hr',price=price['Como catalyst'])
hydrogen = Stream('hydrogen',units='kg/hr',price=price['H2'])

R404 = _units.HydrogenationReactor('R404', ins=(D403-1,hydrogen,Pd_Al2O3_catalyst))
@R404.add_specification(run=True)
def correct_hydrogen_flow():
    feed=R404.ins[0]
    h2=R404.ins[1]
    h2.imol['H2']=feed.F_mol


# Fractionation

D404 = bst.BinaryDistillation('D404', ins=R404-0,outs=('C6_C8_distillate','more_than_C9_bottoms'),
                              LHK=('C8H18','C9H20'),
                              Lr=0.9999,
                              Hr=0.9999,
                              k=1.2,
                              is_divided=True)

D405 = bst.BinaryDistillation('D405', ins=D404-1,outs=('C9_C16_distillate','more_than_C16_bottoms'),
                              LHK=('C16H34','C18H38'),
                              Lr=0.999,
                              Hr=0.999,
                              k=1.2,
                              is_divided=True)

# Standard temperature for storing gasoline is around 15 C
H406 = bst.HXutility('H406', ins=D404-0,outs='gasoline',T=15+273.15,rigorous=True)

# Standard temperature for storing jet fuel is around 15 C
H407 = bst.HXutility('H407', ins=D405-0,outs='jet_fuel',T=15+273.15,rigorous=True)

# Standard temperature for storing diesel is around 20 C
H408 = bst.HXutility('H408', ins=D405-1,outs='diesel',T=20+273.15,rigorous=True)

# sys_upgrading = bst.main_flowsheet.create_system('sys_upgrading')

# sys_upgrading.simulate()

#%%
#==============================================================================
#                      Area 500 Wastewater treatment (WWT)
#==============================================================================

M501 = bst.Mixer(ID='M501',
                 ins=(s.pretreatment_to_WWT,
                      s.U302_to_WWT,
                      s.D302_to_WWT,
                      s.F401_to_WWT))
WWT = bst.create_conventional_wastewater_treatment_system(ID='WWT',
                                                          ins=M501.outs[0],
                                                          outs=('biogas_from_wastewater_treatment',
                                                                'sludge_from_wastewater_treatment',
                                                                'RO_treated_water_from_wastewater_treatment',
                                                                'brine_from_wastewater_treatment'), # brine goes to waste
                                                          autopopulate=False)
                                                             
#%%
#==============================================================================
#                      Area 600 Facilities
#==============================================================================

# Scale ADP, CIP, FWT by flow 
get_flow_tpd = lambda: (energycane.F_mass-energycane.imass['H2O'])*24/907.185

# Streams
natural_gas = Stream('natural_gas',price=price['Natural gas'])
lime_boiler = Stream('lime_boiler', price=price['Lime'])
boiler_chems = Stream('boiler_chems', price=price['Boiler chems'])
ash = Stream('ash', price=price['Ash disposal'])
cooling_tower_chems = Stream('cooling_tower_chems', price=price['Cooling tower chems'])
system_makeup_water = Stream('system_makeup_water', price=price['Water'])

plant_air_in = bst.Stream('plant_air_in',phase='g',units='kg/hr',
                          N2=0.79*1372608*get_flow_tpd()/2205,
                          O2=0.21*1372608*get_flow_tpd()/2205)
CIP_chems_in = Stream('CIP_chems_in',Water=145*get_flow_tpd()/2205, units='kg/hr')
fire_water_in = Stream('fire_water_in', Water=8021*get_flow_tpd()/2205, units='kg/hr')

# Stream mixer to boiler turbogenerator for burning
M601 = bst.Mixer('M601',
                 ins=(s.bagasse_to_CHP,
                      s.filter_cake,
                      s.fiber_fines,
                      s.U302_cell_mass,
                      s.sludge_from_wastewater_treatment), # WWT.outs[1]
                 outs='total_effluent_to_be_burned')    
M602 = bst.Mixer('M602',
                ins=(F.D402_top_product,
                     F.biogas_from_wastewater_treatment), # WWT.outs[0]
                outs='total_gas_to_be_burned')
BT = bst.BoilerTurbogenerator('BT',
                              ins=(s.total_effluent_to_be_burned, # M601 outs[0]
                                   s.total_gas_to_be_burned, # WWT.outs[0],
                                   'boiler_makeup_water',
                                    natural_gas,
                                    lime_boiler,
                                    boiler_chems),
                              outs=('gas_emission',
                                    'boiler_blowdown_water',
                                     ash),
                              satisfy_system_electricity_demand=True,
                              natural_gas_price=0, # price separately set by the stream
                              ash_disposal_price=0) # price separately set by the stream
                              
BT.register_alias('CHP')

# Chilled water package for cooling requirements
CWP = bst.ChilledWaterPackage('CWP')
                                 
CT = bst.CoolingTower('CT')
CT.ins[-1].price = price['Cooling tower chems']
                                                     
# All water used in the system, here only consider water usage
system_makeup_water_streams = (s.boiler_makeup_water,
                               CT.ins[1]) # Second ins of CT (cooling_tower_makeup_water)
                         
system_process_water_streams = (s.water_for_imbibition,
                                s.water_for_rvf,
                                s.water_for_pretreatment,
                                s.water_M301,
                                s.water_U301)
                         
                         
PWC = bst.ProcessWaterCenter('PWC',
                              ins=(s.RO_treated_water_from_wastewater_treatment, # WWT.outs[2],
                                  '',
                                  s.condensate, # E101.outs[1],
                                  system_makeup_water),
                              outs=('','process_water', 'discharged_water'),
                              makeup_water_streams=system_makeup_water_streams,
                              process_water_streams=system_process_water_streams)

ADP = bst.AirDistributionPackage('ADP', ins=plant_air_in, outs='plant_air_out')

CIP = bst.CIPpackage('CIP', ins=CIP_chems_in, outs='CIP_chems_out')

FWT = _units.FireWaterTank('FWT', ins=fire_water_in, outs='fire_water_out')


sys = bst.main_flowsheet.create_system('sys')

sys.simulate()

# HXN = bst.HeatExchangerNetwork('HXN',T_min_app = 5.)
 # sys.save_report(file='report.xlsx',dpi='900')

