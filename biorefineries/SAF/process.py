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
    E = Evaporator
    F = Flash tank
    H = Heat exchanger
    M = Mixer
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units
    PS = Process specifications, not physical units, but for adjusting streams

Processes:
    100: Preprocessing
    200: Pretreatment
    300: Conversion
    400: Separation
    500: Wastewater
    600: Facilities
'''


#%% Set up

import numpy as np
import biosteam as bst
import thermosteam as tmo
from biosteam import Stream, System
from biorefineries.SAF import chemicals as chems
from biorefineries.SAF.utils import convert_ethanol_wt_2_mol
# from biorefineries.SAF.process_settings import _load_process_settings
from biosteam import Flowsheet as F
from biosteam import main_flowsheet,units,SystemFactory
from thermosteam import Stream, reaction as rxn

from biorefineries.SAF import units as _units
from biorefineries.cellulosic import units

F=bst.Flowsheet('energycane_SAF')
bst.main_flowsheet.set_flowsheet(F)

#%%

### A100 Preprocessing process

# Streams
energycane=Stream('energycane',
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
                  price=0.035)
# Assume processing capacity of wet basis=1,600,000 metric tons (MT)/year, with 200 operating days,8000 MT/day,operating hours of 24 hr/day.
# Content based on ''Techno-economic feasibility analysis of engineered energycane-based biorefinery co-producing biodiesel and ethanol''. 
# Bagasse content is based on Enzyme hydrolysis and ethanol fermentation of dilute ammonia pretreated energy cane''.
imbibition_water=Stream('imbibition_water',
                        Water=90000, 
                        units='kg/hr',
                        T = 60+273.15)
H3PO4=Stream('H3PO4',units='kg/hr')
lime=Stream('lime', CaO=333.0,
            Water=2200.0, 
            units='kg/hr', 
            price=0.077)
polymer=Stream('polymer',Flocculant=0.83, 
               units='kg/hr')
rvf_wash_water=bst.Stream('rvf_wash_water',
                          Water=16670,units='kg/hr',
                          T=363.15)

U101=bst.ConveyingBelt('U101', ins=energycane)
U102=bst.MagneticSeparator('U102', U101-0, outs=('A','B'))
U103=bst.Shredder('U103', U102-0)
U104=bst.CrushingMill('U104', ins=(U103-0,''),
                      split=dict(Ash=0.92,
                                 Glucan=0.92,
                                 Glucose=0.04,
                                 Xylan=0.92,
                                 Arabinan=0.92,
                                 Lignin=0.92,
                                 Sucrose=0.04,
                                 Extract=0,
                                 Solids=1), 
                      moisture_content=0.5)

                             
@U104.add_specification(run=True) 
def correct_imbibition_water():
     feed = U104.ins[0]
     moisture=feed.imass['Water']
     dry=feed.F_mass-moisture
     imbibition_water.imass['Water'] = max(2.68333333 * dry - moisture, 0)
     
# Mix in water
M101 = bst.Mixer('M101', ('',imbibition_water),1-U104)

# Screen out fibers
S101 = bst.VibratingScreen('S101', U104-1, ('', 0-M101),
                           split=dict(Ash=0.35,
                                      Glucan=0.35,
                                      Glucose=0.88,
                                      Xylan=0.35,
                                      Arabinan=0.88,
                                      Lignin=0.35,
                                      Solids=0,
                                      Sucrose=0.88,
                                      Extract=0,
                                      Water=0.88))



# Store juice
T101=bst.StorageTank('T101', S101-0, outs='untreated_juice',
                     tau=4,vessel_material='Carbon steel')
U105=bst.ConveyingBelt('U105', U104-0,outs='bagasse')

### Create_juicing_system_up_to_clarification_process

# Heat up juice #
H101=bst.HXutility('H101', ins=T101-0, T=343.15)
        
# Mix with acid #
T102=bst.MixTank('T102', ins=[H101-0,H3PO4])
        
#Pump solution #
P101=bst.Pump('P101', T102-0)
        
# Mix lime solution #
T104=bst.MixTank('T104', ins=(P101-0,lime),tau=0.1)
P102=bst.Pump('P102',ins=T104-0)
        
# Mix recycle with rvf #
M103=bst.Mixer('M103',ins=(P102-0, ''))
        
# Heat #
H102=bst.HXutility('H102', M103-0,T=105+273.15)
        
# Mix in flocculant polymer #
T105=bst.MixTank('T105', ins=(H102-0,polymer),tau=0.1)
        
# Clarify #
C101=bst.Clarifier('C101', T105-0,outs=('clarified_juice',''),
                   split=dict(Ash=0,
                              CaO=0,
                              Glucan=0,
                              Flocculant=0.522,
                              Glucose=0.522,
                              Xylan=0,
                              Arabinan=0,
                              Lignin=0,
                              H3PO4=0.522,
                              Sucrose=0.522,
                              Extract=0.522,
                              Water=0.522))
# Remove solids as filter cake #
C102=bst.RVF('C102',(C101-1,rvf_wash_water),('filter_cake',''),
                     moisture_content=0.80,
                     split=dict(Ash=0.85,
                                CaO=0.85,
                                Glucan=0.85,
                                Glucose=0.01,
                                Xylan=0.85,
                                Lignin=0.85,
                                Arabinan=0.85,
                                Sucrose=0.01,
                                Extract=0.01))
P103=bst.Pump('P103',C102-1,1-M103)
        
# Specifications dependent on energy cane flow rate
@U104.add_specification(run=True)
def correct_flows():
    feed = U104.ins[0]
    F_mass = feed.F_mass
    # correct lime, phosphoric acid
    lime.imass['CaO', 'Water'] = 0.001 * F_mass * np.array([0.046, 0.954])
    H3PO4.imass['H3PO4', 'Water'] = 0.00025 * F_mass

            
# Specifications within a system
@P102.add_specification
def correct_wash_water():
    P102._run()
    solids = P102.outs[0].imol['Ash', 'CaO', 'Glucan',
                                'Xylan','Arabinan', 'Lignin'].sum()
    rvf_wash_water.imol['Water'] = 0.0574 * solids
           
        
### Create_clarified_juice_screening_process
S102=bst.VibratingScreen('S102', C101-0,outs=('screened_juice','fiber_fines'),
                                 split=dict(Ash=0.998,
                                            CaO=0.998,
                                            Glucan=0.0,
                                            Flocculant=0.0,
                                            Glucose=0.998,
                                            Xylan=0.0,
                                            Arabinan=0.0,
                                            Lignin=0.0,
                                            H3PO4=1.0,
                                            Sucrose=0.998,
                                            Water=0.998,
                                            Extract=0.998))
    
S102.mesh_opening=2

H103=bst.HXutility('H103', ins=S102-0,T=25+273.15)      
### Create_screened_juice_concentrating_process 
S103=bst.Splitter('S103', H103-0,
                  split=0.57)
E101=bst.MultiEffectEvaporator('E101', ins=S103-0,outs=('solids','condensate'),
                               V=0.1, V_definition='First-effect',
                               P=(101325, 73581, 50892, 32777))
M104=bst.Mixer('M104', (E101-0,S103-1))
H104=bst.HXutility('H104', M104-0, outs='concentrated_juice', T=32+273.15)
    
E101.target_sugar_concentration=0.2
@E101.add_bounded_numerical_specification(x0=0,x1=1,xtol=1e-5,ytol=1e-2)
def sugar_concentration_at_fraction_evaporation(V):
    E101.V=V
    E101.run_until(M104,inclusive=True)
    Glucose_mass = M104.outs[0].get_mass_fraction('Glucose')
    Sucrose_mass = M104.outs[0].get_mass_fraction('Sucrose')
    sugar_concentration = Glucose_mass + Sucrose_mass
    return E101.target_sugar_concentration - sugar_concentration

# sys_preprocessing = bst.main_flowsheet.create_system('sys_preprocessing')
# sys_preprocessing.simulate()

#%%

### A200 Pretreatment process

solids_loading=0.305
nonsolids=['Water']
T_pretreatment_reactor = 130.+273.15

# Streams

warm_process_water = Stream('warm_process_water', 
                            T=368.15,
                            P=4.7*101325,
                            Water=1)
pretreatment_steam = Stream('pretreatment_steam',
                            phase='g',
                            T=268+273.15,
                            P=13*101325,
                            Water=24534+3490,
                            units='kg/hr')

P=pretreatment_steam.chemicals['H2O'].Psat(T_pretreatment_reactor + 25)
M201=bst.SteamMixer('M201',ins=(U105-0,pretreatment_steam,warm_process_water),
                    P=P,solids_loading=solids_loading,
                    )
R201=units.PretreatmentReactorSystem('R201', M201-0,T=T_pretreatment_reactor)
P201=units.BlowdownDischargePump('P201',R201-1)
F201=units.PretreatmentFlash('F201', P201-0,P=101325,Q=0)
M202=bst.Mixer('M202',(R201-0,F201-0))
units.WasteVaporCondenser('H201', M202-0,'pretreatment_wastewater',V=0)

P202=units.HydrolyzatePump('P202', F201-1)

U201=bst.HammerMill('U201', ins=P202-0,outs='pretreated_bagasse')

    

sys_pretreatment=bst.main_flowsheet.create_system('sys_pretreatment')
sys_pretreatment.simulate()

#%%

# A300 Fermentation process

# Stream

enzyme_M301 = Stream('enzyme_M301',units='kg/hr')
water_M301=Stream('water_M301',units='kg/hr')

CSL_R301 = Stream('CSL_R301',units='kg/hr')
CSL_R302 = Stream('CSL_R302',units='kg/hr')

DAP_R301 = Stream('DAP_R301',units='kg/hr')
DAP_R302 = Stream('DAP_R302',units='kg/hr')

water_U301 = Stream('water_U301',units='kg/hr')

# SSCF

nonsolids=chems.default_nonsolids
insoluble_solids=chems.default_insoluble_solids
ignored=chems.default_ignored

insoluble_solids_loading=10.3


M301 = _units.EnzymeHydrolysateMixer('M301', ins=(U201-0,enzyme_M301,water_M301),
                                   enzyme_loading=20,solids_loading=0.2)

R301 = _units.SimultaneousSaccharificationAndCoFermentation('R301', ins=(M301-0,'',CSL_R301,DAP_R301,H104-0),
                                              outs=('R301_g','effluent','side_draw'))
R302 = _units.SeedTrain('R302', ins=(R301-2,CSL_R302,DAP_R302),
                       outs=('R302_g','seed'))
T301 = _units.SeedHoldTank('T301', ins=R302-1,outs=1-R301)

M302 = bst.Mixer('M302',ins=(R301-0,R302-0),outs='fermentation_vapor')

def update_U301_water():
    M302._run()
    water_U301.imass['Water']=26836/21759 * M302.F_mass_in
M302.specifications=update_U301_water()

U301 = bst.VentScrubber('U301', ins=(water_U301, M302-0), 
                        outs=('U301_vent', 'U301_recycled'),
                        gas=('CO2', 'NH3', 'O2'))
M303=bst.Mixer('M303',ins=(R301-1,U301-1))
T302=_units.BeerTank('T302', ins=M303-0)

# Heat up crude beer by exchanging heat with stillage
H302=bst.HXprocess('H302', ins=(T302-0,''), 
                   phase0='l',phase1='l',U=1.28)

# Remove solids from fermentation broth, based on the pressure filter in ref [1]
# Moisture content is 35% in ref [1] but 25% in ref [2], used 35% to be conservative
U302=_units.CellMassFilter('U302', ins=H302-1,outs=('U302_cell_mass','U302_to_WWT'),
                           moisture_content=0.35,split=0.99)

# Beer column
xbot = convert_ethanol_wt_2_mol(0.00001)
ytop = convert_ethanol_wt_2_mol(0.5)

D301 = bst.BinaryDistillation('D301',ins=H302-0,k=1.25,Rmin=0.6,
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
D302 = bst.BinaryDistillation('D302', ins=M304-0, k=1.25, Rmin=0.6,
                                    P=101325, y_top=ytop, x_bot=xbot,
                                    LHK=('Ethanol', 'Water'),
                                    tray_material='Stainless steel 304',
                                    vessel_material='Stainless steel 304',
                                    is_divided=True)
D302.reboiler.U = 1.85
D302_P = bst.Pump('D302_P', ins=D302-1, outs='D302_to_WWT')
D302_P.BM = 3.1

D302_H = bst.HXutility('D302_H', ins=D302-0, T=115+283.15, V=1)

# Molecular sieve, split based on streams 515 and 511 in ref [1]
split_ethanol = 1 - 21673/27022
split_water = 1 - 108/2164
S301 = bst.MolecularSieve('S301', ins=D302_H-0, outs=(1-M304, ''),
                                split=(split_ethanol, split_water),
                                order=('Ethanol', 'Water'))
# Condense ethanol product
S301_H = bst.HXutility('S301_H', ins=S301-1, outs='ethanol_to_storage',
                             V=0, T=80+273.15)



sys_SSCF=bst.main_flowsheet.create_system('sys_SSCF')

sys_SSCF.simulate()

#%%

### A400 Upgrading

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

ethanol_to_storage = Stream('ethanol_to_storage',
                            phase='l', T=298.15, P=101325,
                            H2O=10.3,
                            Ethanol=694,
                            units='kmol/hr') 
# Stream

NaOH=Stream('NaOH',NaOH=0.5,Water=0.5,units='kg/hr')
Syndol_catalyst=Stream('Syndol_catalyst',units='kg/hr')

P401=bst.Pump('P401',ins=ethanol_to_storage,P=4.5*101325)
H400=bst.HXutility('H400', ins=P401-0,V=1,rigorous=True)

R401=_units.AdiabaticFixedbedDehydrationReactor('R401', ins=(H400-0,Syndol_catalyst),outs=('','spent_catalyst'))

# Depressurize to 1 bar before quenching
V401=bst.IsenthalpicValve('V401', ins=R401-0,P=101325)
H402=bst.HXutility('H402', ins=V401-0,T=50+273.15,outs='crude_ethylene',rigorous=True)
S401=bst.PhaseSplitter('S401', ins=H402-0,outs=('vapor','liquid'))

S401.target_ethylene_cwt=0.921
@H402.add_bounded_numerical_specification(x0=20+273.15,x1=80+273.15,xtol=1e-5,ytol=1e-2)
def ethylene_cwt_at_T(T):
    H402.T=T
    H402.run_until(S401,inclusive=True)
    ethylene_cwt=S401.outs[0].get_mass_fraction('Ethylene')
    return S401.target_ethylene_cwt - ethylene_cwt

# Pressurize to 27 bar before purification;
# Based on Bioethylene Production from Ethanol: A Review and Techno-economical Evaluation
C401=bst.MultistageCompressor('C401', ins=S401-0,n_stages=3,pr=3,eta=0.72,vle=True)
# Condense water
S402=bst.PhaseSplitter('S402', ins=C401-0,outs=('vapor','liquid'))
# Remove CO2
U402=_units.CausticTower('U402', ins=(S402-0, NaOH),P=27*101325)
# Remove water and ethanol
U403=bst.MolecularSieve('U403', ins=U402-0,outs=('water_and_ethanol','dried_ethylene'),
                        split=dict(Water=1,
                                   Ethanol=1))

# Reduce temperature before cryogenic distillation low temperature + high pressure
# Temperature is based on Bioethylene Production from Ethanol: A Review and Techno-economical Evaluation
H404=bst.HXutility('H404',ins=U403-1,T=-25+273.15,rigorous=True)

# D401 cannot work!!
D401=bst.BinaryDistillation('D401',ins=H404-0,outs=('distillate','bottoms'),
                            LHK=('Ethylene','Ethane'),
                            Lr=681.5/682,
                            Hr=1.7/1.86,
                            product_specification_format='Recovery',
                            k=1.2,
                            is_divided=True,
                            partial_condenser=True,
                            )
# D402=_units.Stripper`Column('D402', ins=D401-0,outs=('distillate','bottoms'),
#                           LHK=(llight_keys,'ethylene'),
#                           is_divided=False,
#                           y_top=1,
#                           x_bottom=0,
#                           k=2,
#                           )


# Oligomerization

first_catalyst=Stream('first_catalyst',units='kg/hr')
second_catalyst=Stream('second_catalyst',units='kg/hr')

ethylene=Stream('ethylene',Ethylene=682,units='kmol/hr',phase='g',T=25+273.15)

# First oligomerization
C402=bst.IsentropicCompressor('C402',ins=ethylene,P=21*101325,eta=1, compressor_type='Centrifugal')
R402=_units.EthyleneDimerizationReactor('R402', ins=(C402-0,first_catalyst),outs=('','spent_catalyst'))

# Second oligomerization
R403=_units.OligomerizationReactor('R403',ins=(R402-0,'',second_catalyst))

D403=bst.BinaryDistillation('D403', ins=R403-0,outs=('distillate','bottoms'),
                            LHK=('1-hexene','1-heptene'),
                            y_top=0.999,
                            x_bot=0.001,
                            k=3)
H405=bst.HXutility('H405', ins=D403-0,outs=1-R403,V=0,rigorous=True)
# Light olefins recycled to the 2nd oligomerization reactor


# Hydrogenation

# Stream

Pd_Al2O3_catalyst=Stream('Pd_Al2O3_catalyst',units='kg/hr')
como_catalyst=Stream('como_catalyst',units='kg/hr')
h2=Stream('h2',H2=1,P=34.5*101325,units='kg/hr')

R403=_units.HydrogenationReactor('R403', ins=(D403-1,h2,Pd_Al2O3_catalyst))




# Fractionation

D404=bst.BinaryDistillation('D404', ins=R403-0,outs=('distillate','bottoms'),
                            LHK=('C8H18','C9H20'),
                            Lr=0.999,
                            Hr=0.999,
                            k=3,
                            is_divided=True)

D405=bst.BinaryDistillation('D405', ins=D404-1,outs=('distillate','bottoms'),
                            LHK=('C16H34','C18H38'),
                            Lr=0.999,
                            Hr=0.999,
                            k=3,
                            is_divided=True)
  
sys=bst.main_flowsheet.create_system('sys')

sys.simulate()
