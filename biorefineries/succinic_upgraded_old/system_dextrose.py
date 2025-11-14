#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 19:22:32 2025

@author: princyk2
"""



import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from numba import njit
from biorefineries.succinic._process_specification import ProcessSpecification
from biosteam import System
from thermosteam import Stream
from biosteam.process_tools import UnitGroup
from biosteam import SystemFactory
from biorefineries.succinic import units
from biorefineries.succinic_upgraded import units_adsorption
from biorefineries.succinic import facilities
from biorefineries.succinic_upgraded.process_settings import price, CFs, chem_index, _GDP_2007_to_2010
from biorefineries.succinic.utils import find_split, splits_df, baseline_feedflow


from biorefineries.succinic.chemicals_data import chems, chemical_groups, \
                                                                soluble_organics, combustibles
# TODO  changes chems to succininc_upgraded                                                              # 
from biorefineries.succinic.tea import TemplateTEA as SuccinicTEA
from biorefineries.succinic.lca import SuccinicLCA
# from biosteam.units import SingleComponentAdsorptionColumn
from biorefineries.succinic.crystallization_curvefit import Ct_given_C0
# from biorefineries.make_a_biorefinery.auto_waste_management import AutoWasteManagement
from biorefineries.cellulosic import create_facilities
from contourplots import stacked_bar_plot
import math
import matplotlib.pyplot as plt
from biorefineries import corn
# from biorefineries.succinic import system_glucose
# from biorefineries.succinic.system_glucose import chems

from biorefineries.succinic.system_glucose import (s,
    U101, M201, U302, S301, H302, P303, F301, H301,
    M304, M304_H, S302, CO2_fermentation,makeup_MEA_A301,sulfuric_acid_R401, base_neutralization,A301_makeup_water,natural_gas_drying,air_lagoon,aerobic_caustic, base_fermentation,\
    base_neutralization, S303, CSL, R303, M305,sulfuric_acid_acidulation,sulfuric_acid_R401,
    T301, A301, K301, H401, M401, R401, R401_P, S405, M404,
    F401, F401_P, C401, S402, F402, F402_P, C402, S403, F403, F403_P,C403,
    S404, M402, M403, S406, M501, M510,T601_P,
     get_succinic_acid_solubility_gpL, HeatExchangerNetwork,\
        T601,T602,create_succinic_sys
)

  
        
# U101, M201, U302, S308,H302,P303,F301,H301,H304,M304_H,S302,\
 # S301,R302,R303, M305,T301, A301, K301, K302
 
IQ_interpolation = flx.IQ_interpolation
HeatExchangerNetwork = bst.facilities.HeatExchangerNetwork

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# bst.speed_up()
flowsheet = bst.Flowsheet('succinic')
bst.main_flowsheet.set_flowsheet(flowsheet)
bst.preferences.update(flow = 'kg/hr',composition = True)
bst.preferences.N = 300
# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0

# Baseline cost year is 2016

bst.CE = bst.units.design_tools.CEPCI_by_year[2019]

# Set default thermo object for the system
tmo.settings.set_thermo(chems)

# System settings
System.converge_method = 'wegstein'
# princy changed the function from (default_converge_method) to converge_method as it not defined in the systen ffile
# System.default_converge_method = 'fixed-point'
# System.default_converge_method = 'aitken'

System.default_maxiter = 100
System.default_molar_tolerance = 0.1
System.default_relative_molar_tolerance = 0.0001 # supersedes absolute tolerance
System.strict_convergence = True # True => throw exception if system does not converge; False => continue with unconverged system
exp = math.exp
feedstock_ID = 'Glucose'

l_phenylalanine=bst.Stream('l-phenylalanine')



# %% System
@SystemFactory(ID = 'succinic_sys')
def create_succinic_sys(ins, outs):
    def fix_split(isplit, ID):
        isplit['SuccinicAcid', 'LacticAcid', 'Ethanol'] = isplit[ID]
    process_groups = []
    
    f = bst.main_flowsheet
    u = f.unit
    s = f.stream
 
     # %% Feedstock
     
    feedstock = bst.Stream('glucose_feedstock', Glucose=1., Water=1., units='kmol/h')
    feedstock.price = price['Glucose']*0.909 # dextrose monohydrate stream is 90.9 wt% glucose
    feedstock.F_mass = 200_000 # initial value; updated by spec.set_production_capacity
    # feedstock_sys=bst.System('feedstock_sys', path=[U101, M201, U302])
  
    # U101-0-M201-0-U302
    
    # %%
    # =============================================================================
    # Conversion units
    # =============================================================================
    M304.water_to_sugar_mol_ratio = 200  
    dilution_water = Stream('dilution_water', price=price['Makeup water'])                                                          
    U302-0-0-S301
   
    S301-0-0-H302-0-0-P303
    S301-1-0-F301
    F301-0-0-H301
    (H301-0,dilution_water)-M304
    M304-0-M304_H
    M304_H-0-0-S302
    (CO2_fermentation)-S303
    
    
    M304-0-M304_H
    M304_H-0-S302
    
    
      # Corn steep liquor as nitrogen nutrient for microbes, flow updated in R301
 
    magnesiumchloride_hexahydrate_fermentation =Stream('magnesiumchloride_hexahydrate_fermentation', units='kg/hr',
                                             price=price['MagnesiumChlorideHexahydrate'],) 
                                       
                                          
    zincsulfate_heptahydrate_fermentation =Stream('zincsulfate_heptahydrate_fermentation', units='kg/hr', 
                                      price=price['Zinc Sulfate Heptahydrate'],)
    diammonium_sulfate_fermentation = Stream('diammonium_sulfate_fermentation', units='kg/hr',
                              price=price['Diammonium sulfate'],)
    
    magnesium_sulfate_fermentation = Stream('magnesium_sulfate_fermentation', units='kg/hr',
                              price=price['Magnesium sulfate'],)
     # adsorption_chemicals
     # nutruients
    
    SuccinicAcid=Stream('SuccinicAcid', units='kg/hr')
    
    
    R302 = units_adsorption.CoFermentation('R302', 
                                    ins=(S302-1, 'seed', CSL, S303-1, base_neutralization, 
                                         # M300_R-0, 
                                         '', # CO2 recycle (0)
                                         diammonium_sulfate_fermentation, magnesium_sulfate_fermentation,
                                        # magnesiumchloride_hexahydrate_fermentation, 
                                        # zincsulfate_heptahydrate_fermentation, 
                                         '',
                                         ), # air (K302-0)
                                    outs=('fermentation_broth', 'fermentation_vent'),
                                    mode='batch',
                                    V_wf=0.8,
                                    neutralization=False,
                                    pH_control=True,
                                    neutralizing_agent='Lime',
                                    base_neutralizes_product=False, # if True, acidulation is performed downstream
                                    )
    
    #TODO didnt change the pH control used Lime instead of NH4OH
    #inlcuded teh model of tactical CO2 dosing
   
    @R302.add_specification(run=False)
    def R302_spec(): # note: effluent always has 0 CSL
         # R302.show(N=100)
         R302._run()
         if R302.ins[2].F_mol:
             R302.ins[2].F_mass*=1./(1-S302.split[0])
         R302._run()
         S303.simulate()
         K302.specifications[0]()   
    (S302-0,S303-0)-R303
  
    (R302-1,R303-1)-M305
    R303-0-0-T301
    T301-0-1-R302
    (M305-0, makeup_MEA_A301,A301_makeup_water)-A301
    
    A301-1-0-K301
  
    K301-0-5-R302
   
    

  
  
    K302 = units.IsothermalCompressor('K302', ins='atmospheric_air', outs=('pressurized_air'), 
                                     P=3e7, # set to R302.air_pressure when simulated
                                     # vle=True,
                                     eta=0.6,
                                     driver='Electric motor',
                                     )
    
    K302-0-8-R302
    @K302.add_specification(run=False)
    def K302_spec():
        K302.P = R302.air_pressure
        K302.ins[0].phase = 'g'
        K302.ins[0].mol[:] = K302.outs[0].mol[:]
        K302._run()
     

    # %% seperation units
    
    R302-0-0-H401
    
    H401.Ts_to_try = np.linspace(R302.T, 90.+273.15, 50)
    @H401.add_specification(run=False)
    def H401_spec():
        Ts = H401.Ts_to_try
        outstream = H401.outs[0]
        for T in Ts:
            H401.T = T
            H401._run()
            if outstream.imass['SuccinicAcid']/outstream.F_vol - get_succinic_acid_solubility_gpL(H401.T) < 0.:
                break
    H401-0-0-M401
  
    S401 = bst.SolidsCentrifuge('S401', ins=M401-0, 
                            outs=('S401_solid_waste', 'S401_1'),
                            solids=['FermMicrobe',
                                   
                                    'Extract',
                                    'Lignin',
                                    'Protein', 'Arabinose', 'Ash', 
                                    'AmmoniumSulfate',
                                    'Glucan', 
                                    'Xylan', 
                                    'SolubleLignin',
                                     'XyloseOligomer', 'GlucoseOligomer','Cellobiose'], 
                            split={'FermMicrobe':0.995,
                                   'Extract': 0.995, #!!! this assumption is purely to keep the same separation process as for the sugarcane->succinic biorefinery, 
                                                     #    update if using this model in a study focused on cornstover -> succinic biorefineries
                                   'Lignin': 0.995,  #!!! same as Extract
                                   'Protein': 0.995, #!!! same as Extract
                                   'Arabinose': 0.995, #!!! same as Extract
                                   'Ash': 0.995,     #!!! same as Extract
                                   'AmmoniumSulfate': 0.995, #!!! same as Extract
                                    'Glucan': 0.995, #!!! same as Extract
                                    'Xylan': 0.995, #!!! same as Extract
                                   'SolubleLignin': 0.995, #!!! same as Extract
                                     #!!! same as Extract
                                    'XyloseOligomer': 0.995, 
                                    'GlucoseOligomer': 0.995,
                                  
                                         #!!! same as Extract
                                    'Cellobiose': 0.995, #!!! same as Extract
                                    })
    
        # FermMicrobe=bst.Stream('FermMicrobe', )
        # N_flow = FermMicrobe.get_atomic_flow('N')
        
        
    @S401.add_specification(run=False)
    def S401_succinic_acid_precipitation_spec():
          
        S401._run()
        instream = S401.ins[0]
        tot_SA = instream.imass['SuccinicAcid']
        dissolved_SA = get_succinic_acid_solubility_gpL(instream.T) * instream.F_vol
        S401.outs[1].imass['SuccinicAcid'] = min(tot_SA, dissolved_SA)
        S401.outs[0].imass['SuccinicAcid'] = max(0, tot_SA - dissolved_SA)
        liquid = S401.outs[1]
    
    # Isolate FermMicrobe
        microbe_stream = liquid.copy()
        microbe_stream.empty()
        microbe_stream.imass['FermMicrobe'] = liquid.imass['FermMicrobe']
        # Calculate atomic nitrogen flow
        N_flow = microbe_stream.get_atomic_flow('N')
        # Assume 1 N atom per mol of amino acid
        l_phenylalanine = N_flow / 1
        # Remove FermMicrobe and add amino acid
        liquid.imass['FermMicrobe'] = 0
        liquid.imol['l-phenylalanine'] += l_phenylalanine
        # Store for downstream use
        S401.l_phenylalanine = l_phenylalanine
        S401._run()
            
    #         Assumptions
    # 1 N per amino acid: This is a reasonable average for most amino acids.
    
    # No carbon constraint: Allows flexibility in modeling mixed carbon sources.
                
            
            # aminoacid defined as 'l-phenylalanine'] in the chemical data because the adsoprtion model we used based on 'l-phenylalanine'] as adsorbates
    H402 = units_adsorption.HXutility('H402', ins=S401-1, outs = ('broth_to_adsorbtion',), T=20. + 273.15)
    make_up_water = Stream('make_up_water', water=1,units='kg/hr', price=price['Makeup water'],)
    # the temperature is updated to 40 ˚C but #TODO temp is changed from operational temp 40 to 20 ˚C becuase data for 20 ˚ C is only available)
    recycled_water_A401 =bst.Stream(ID='recycled_water_A401',units='kg/hr', price=price['Makeup water'],)
    recycled_water_F400_P=bst.Stream(ID='recycled_water_F400_P',units='kg/hr', price=price['Makeup water'],)
   
    NaOH=bst.Stream('NaOH', units='kg/hr',  price=price['NaOH'])
    T400 = units_adsorption.CausticStorage('T400', ins = NaOH, outs = 'fresh_NaOH')
    @T400.add_specification(run=False)
    def T400_spec():
        T400.ins[0].imass['NaOH']=A401.outs[1].imass['NaOH']
        
    # TODO price of NaoH
    # M400=bst.Mixer('M400',ins= (make_up_water,'' ), outs= 'mixed_water_for_1.0M_NaOH_solution')
    M400=bst.Mixer('M400',ins= (make_up_water,T400-0,''  ), outs= '1.0M_NaOH_solution')

         # Kmol/hr × MW of phenylalanine (Kg/Kmol) to Kg/hr
    # solubility of l phenylalanine in 1 N NaOH solution is 0.0465 Kg/L~46.5Kg/m3
       # Calculate required NaOH solution volume (mL/hr)
       
       # 1N NaoH solution NaOH mass fraction=0.018
       
    
    # M400_N = bst.Mixer('M400_N', ins=(M400-0,T400-0), outs=('1.0M_NaOH_solution'))
    
        # M400 spec added according to the amount of NaOH solution required to desorb 'aminoacid' [3]
    #refer source 4: file:///Users/princyk2/Downloads/RobinM_Smith_Chemicalprocessdesignandint.pdf
    #page No 191 - from Robin Smith - the isotherm model- Freundlich isotherm is mainly selected for the adsoprtion from liquids
    #certainily if the identity of solute is unknown and when measuring the colur intensity as a proxy for concentration
    #his is particularly relevant in applications like removing color from liquids through adsorption. In such cases, the "concentration of solute" can be measured by a colorimeter. 
    #The Freundlich equation can then be expressed in terms of arbitrary units of color intensity, provided that this color scale varies linearly with the actual concentration of the solute responsible for the colo 
    make_up_AC=bst.Stream('Make_Up_adsorbent',units='kg/hr', price=price['ActivatedCarbon'])
    @M400.add_specification(run=False)
    def M400_spec():
        
     
        # MW_phenylalanine = 165.19  # kg/kmol
        # adsorbate_flow= S401.l_phenylalanine*165.19
        # volume_per_Kg = 1 / 46.5    #m3/Kg
        # required_NaOH_solution_volume = volume_per_Kg * adsorbate_flow*1040
        # density of 1 N NaOH=1040kg/m3  
 
        A401._run()
        # F400_P._run()
        M400.ins[0].phase = 'l'
        M400.ins[2].phase = 'l'
        M400.ins[1].imass['NaOH']=A401.outs[1].imass['NaOH']
        T400.specifications[0]()
    
        Required_Water_mass =   A401.outs[1].imass['Water']
        M400.outs[0].imass['Water'] = Required_Water_mass
        recycled_water = M400.ins[2].imass['Water']
        total_water_req = Required_Water_mass
        amt_of_fresh_water = max(0, total_water_req - recycled_water)
        M400.ins[0].imass['Water'] = amt_of_fresh_water
    
        M400._run()
       
        
       

        
       
  
    
    A401=units_adsorption.SingleComponentAdsorptionColumn('A401',ins=(H402-0,M400-0,make_up_AC), outs=('broth_post_adsorption','spent_fluid','spent_adsorbent' ),
                                        adsorbate = ['l-phenylalanine'],
                                                          # [3]contains **{'l-phenylalanine'} adsorption study so for this study i used **{'l-phenylalanine} as aminoacid present in the dextrose 
                                                          # I didnt find any sources that evident to prove vitamins A and D present in th dextrose.
                                        cycle_time= 3, superficial_velocity=14,
                                        # assume superficial velocity as in the range between ( velocity fro liquids typically between 0.001-0.004 m/s)
                                        
                                        
                                        # regeneration_superficial_velocity = None, #TODO: this does not affect the design and cost
                                        isotherm_model='Langmuir', isotherm_args=(5.4,64.4),
                                    # (5.4,64.4)
                                        #  adsorption model-[9]
                                        regeneration_fluid={'T':20+273.15,'NaOH': 0.018 , 'Water': 1-0.018}
                                                    
                                                            #4.08e+03 }
                                                            # 'NaOH': 1.8,'water': 98.2}
                                                            ,#TODO: change water
                                        rho_adsorbent= 2000,#TODO: this was based on 
                                        particle_diameter=0.004,#refre[11]GAC is composed of larger granules, usually ranging from 0.4 mm to 4 mm in size.avg 2.2 mm
                                        # density of DARCO-https://www.usbio.net/biochemicals/136542/activated-charcoal-darco-g-60
                                                      # K [L/mg], q_max [mg/g]) and cycle time , the data is obtained from https://www.sciencedirect.com/science/article/pii/S1665642316300773
                                        # superficial_velocity is assumed to be 4 from the book(page no:119) rule of thumb(file:///Users/princyk2/Downloads/Rules%20of%20Thumb%20in%20Engineering%20Practice%20-%202007%20-%20Woods%20(1).pdf)
                                        # cycle time is total 8 hours (adsoprtion+regenration)
                                        # regeneration temp is 20 ˚C
                                        k=1.2,
                                        void_fraction = 0.5, 
                                        k_regeneration=360,
                                        # _estimate_regeneration_time=, #TODO: does not affect the cost and design, and can be assumed anything else if needed
                                        # Mass transfer coefficient of the regeration solvent [1/h]
            regeneration_isotherm_args= 'Langmuir',
            regeneration_isotherm_model= (1,1000) ,adsorbent='ActivatedCarbon',) # [1 / hr]k=0.3,)
   
        
       
        
        
   
    @A401.add_specification(run=False)
    def A401_spec(): 
     
        A401._run() 
        # A401.outs[1].imass['NaOH']=A401.ins[1].imass['NaOH']#actual Naoh and water flow rate required by the adsorptioncolumn based on chemistry
        # A401.outs[1].imass['Water'] = A401.ins[1].imass['Water']
#adsorption unit module calculte regeneration velocity not based on chemistry the calculation is based on inventories, length, area adn regneration vleocity, cycle time =regenration time 
        loss_effectiveness=A401.loss_effectiveness #refer[]
        no_of_cycles=A401.no_cycles
        
        total_loss_fraction=loss_effectiveness/(A401.cycle_time*no_of_cycles)
        Total_mass_AC = float(A401.area) * float(A401.column_length) * A401.rho_adsorbent
     
       
        make_up_AC_flow= Total_mass_AC * total_loss_fraction
        A401.ins[2].imass['ActivatedCarbon']=make_up_AC_flow
        A401._run() 
        # M400.specifications[0]()
     
    S400_W = bst.Splitter('S400_W', ins= A401-1, outs=('NaOH_distillation', 'water_wwt'),split=0.25)
    
    F400 = bst.MultiEffectEvaporator('F400', ins=S400_W-0, outs=('F400_l', 'F400_g'),
                                              P = (101325, 73581, 50892, 32777, 20000), V_definition='Overall', V=0.1) 
    
    F400_P =bst.Pump('F400_P', ins=F400-1, outs= recycled_water_F400_P, P=101325)
    
    F400_P-0-2-M400          
          
            
    (A401-0,sulfuric_acid_R401)-R401  
    R401-0-0-R401_P
   
        
  
    R401_P-0-0-S405
   
    S405-1-0-M401
    M401-0-0-F401
    F401.target_concentration = 0.24
    F401.add_specification()
    F401-0-0-F401_P
    C401-0-0-F401_P
    C401-0-0-S402
    S402.recovery = 0.85
    S401-1-0-F402
    F402-0-0-F403_P
    F403_P-0-0-C403
    C403-0-0-S404
    S404.recovery = 0.85
    (S402-0, S403-0, S404-0)-M402
    F404 = bst.DrumDryer('F404', ins=(M402-0, 'dryer_air_in', natural_gas_drying,),
                         outs=('dry_solids', 'hot_air', 'dryer_emissions'),
                         split={'SuccinicAcid': 1e-4,
                                    'FermMicrobe': 0.}
                         )
    F404.ins[2].price = 0.2527
    (F401-1, F402-1, F403-1)-M403
    S404-1-0-S406
    S406-0-1-M404
          # %% 
      
      # =============================================================================
      # Wastewater treatment units
      # =============================================================================
      # Mix waste liquids for treatment 
    # dilution_water = Stream('dilution_water', price=price['Makeup water'])
    (P303-0, F301-1, S406-1,S400_W-1, M403-0)-M501
   

    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=M501-0,
        mockup=True,
        area=500,
    )
    (S401-0,A401-2,F400-0)-M510
     
    MX = bst.Mixer(400, ['', ''])
    imbibition_water = Stream('imbibition_water', price=price['Makeup water'])
    rvf_wash_water = Stream('rvf_wash_water', price=price['Makeup water'])
    dilution_water = Stream('dilution_water', price=price['Makeup water'])
    s = flowsheet.stream
    create_facilities(
        solids_to_boiler=M510-0,
        gas_to_boiler=wastewater_treatment_sys-0,
        process_water_streams=[
           s.imbibition_water,
         s.rvf_wash_water,
         s.dilution_water,
         # s.makeup_water,
         # s.fire_water,
         # s.boiler_makeup_water,
         # s.CIP,
         # s.recirculated_chilled_water,
         # s.s.3,
         # s.cooling_tower_makeup_water,
         # s.cooling_tower_chemicals,
         ],
        feedstock=s.glucose_feedstock,
        RO_water=wastewater_treatment_sys-2,
        recycle_process_water=MX-0,
        BT_area=700,
        area=900,
    )
    
    #%% 
    
     # =============================================================================
     # Facilities units
     # =============================================================================
    
     # # 7-day storage time, similar to ethanol's in Humbird et al.   
     # Product storage
    F404-0-0-T601
    T601.line = 'SuccinicAcidStorageTank'
    T601_P-0-0-T601
    
    # # Storage tanks for other process inputs (besides water) limestorage didnt change to maintain the consistency in the unit
    T602 = units.LimeStorage('T604', ins=base_fermentation, 
                             outs=base_neutralization,
                             # tau=7.*24., V_wf=0.9,
                                          # vessel_type='Floating roof',
                                          # vessel_material='Stainless steel',
                                          )
    
    T603 = units.SulfuricAcidStorage('T602', ins=sulfuric_acid_acidulation,
                                     outs=sulfuric_acid_R401)
    
   

    
    
    
    
    # M510-0 is the liquid/solid mixture, R501-0 is the biogas, blowdown is discharged
    
    HXN = HeatExchangerNetwork('HXN1001',
                                              # ignored=[H401, H402],
                                              cache_network=True,
                                              )
    def HXN_no_run_cost():
        HXN.heat_utilities = []
        HXN._installed_cost = 0.
       
#%%
# SuccinicAcid=bst.Stream('SuccinicAcid',units='kg/hr')

succinic_sys = create_succinic_sys()
u = flowsheet.unit
s = flowsheet.stream
feedstock = s.glucose_feedstock
product_stream = s.SuccinicAcid

feeds = succinic_sys.feeds

products = [product_stream] # Don't include gypsum since we want to include carbon impurities in GWP calculation

u.BT701.ID = 'BT701'
u.CT901.ID = 'CT801'
u.CWP901.ID = 'CWP802'
u.CIP901.ID = 'CIP901'
u.ADP901.ID = 'ADP902'
u.FWT901.ID = 'FWT903'
u.PWC901.ID = 'PWC904'


BT = flowsheet('BT')

BT.natural_gas_price = 0.2527
BT.ins[4].price = price['Lime']

globals().update(flowsheet.to_dict())
unit_groups = bst.UnitGroup.group_by_area(succinic_sys.units)

# # %% keeep this commented till other results are fixed
# # =============================================================================
# # TEA
# # =============================================================================

# # Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legislation)

# get_flow_dry_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
# succinic_tea = SuccinicTEA(system=succinic_sys, IRR=0.10, duration=(2016, 2046),
#         depreciation='MACRS7', income_tax=0.21, 
#         # operating_days=0.9*365,
#         operating_days = 200,
#         lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
#         startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
#         startup_VOCfrac=0.75, WC_over_FCI=0.05,
#         finance_interest=0.08, finance_years=10, finance_fraction=0.4,
#         # biosteam Splitters and Mixers have no cost, 
#         # cost of all wastewater treatment units are included in WWT_cost,
#         # BT is not included in this TEA
#         OSBL_units=(
#                     u.U501,
#                     # u.T601, u.T602, 
#                     # u.T601, u.T602, u.T603, u.T604,
#                     # u.T606, u.T606_P,
#                     u.BT701, u.CT801, u.CWP802, u.CIP901, u.ADP902, u.FWT903, u.PWC904,
#                     ),
#         warehouse=0.04, site_development=0.09, additional_piping=0.045,
#         proratable_costs=0.10, field_expenses=0.10, construction=0.20,
#         contingency=0.10, other_indirect_costs=0.10, 
#         labor_cost=3212962*get_flow_dry_tpd()/2205,
#         labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
#         steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT701)

# seed_train_system = bst.System('seed_train_system', path=(u.S302, u.R303, u.T301))

# theoretical_max_g_succinic_acid_per_g_glucose = (2*chems.SuccinicAcid.MW)/(chems.Glucose.MW) # ~1.311 g-succinic-acid/g-glucose

# spec = ProcessSpecification(
#     evaporator = u.F301,
#     pump = None,
#     mixer = u.M304,
#     heat_exchanger = u.M304_H,
#     seed_train_system = [],
#     seed_train = u.R303,
#     reactor= u.R302,
#     reaction_name='fermentation_reaction',
#     substrates=('Xylose', 'Glucose'),
#     products=('SuccinicAcid',),
    
#     # spec_1=0.19,
#     spec_1=0.7,
#     spec_2=100.,
#     spec_3=1.0,
#     xylose_utilization_fraction = 0.80,
#     feedstock = feedstock,
#     dehydration_reactor = None,
#     byproduct_streams = [],
#     HXN = u.HXN1001,
#     maximum_inhibitor_concentration = 1.,
#     pre_conversion_units = succinic_sys.split(u.M304.ins[0])[0],
    
#     # set baseline fermentation performance here
#     # baseline_yield = 0.4478, # 0.587 g/g-glucose
#     # baseline_yield = 0.3738, # 0.49 g/g-glucose
    
#     baseline_yield = 0.65/theoretical_max_g_succinic_acid_per_g_glucose, # 49.6 % of theoretical maximum
#     # 65 % actual yield
#     baseline_titer = 78.37,
#     baseline_productivity = 0.95,
#     neutralization=False,
#     baseline_pyruvic_acid_yield = 0.006,
#     baseline_cell_mass_yield = 0.186,
#     #TODO check if baseline pyruvic and cell mass is needed

#     feedstock_mass = feedstock.F_mass,
#     pretreatment_reactor = None)


# def M304_titer_obj_fn(water_to_sugar_mol_ratio):
#     # M304, R302 = u.M304, u.R302
#     u.M304.water_to_sugar_mol_ratio = water_to_sugar_mol_ratio
#     u.M304.specifications[0]()
#     u.M304_H._run()
#     u.S302._run()
#     u.R303.specifications[0]()
#     u.T301._run()
#     # R302.simulate()
#     u.R302.specifications[0]()
#     return u.R302.effluent_titer - u.R302.titer_to_load

# def F301_titer_obj_fn(V):
#     u.F301.V = V
#     u.F301._run()
#     u.H301._run()
#     u.M304.specifications[0]()
#     u.M304_H._run()
#     u.S302._run()
#     u.R303.specifications[0]()
#     u.T301._run()
#     # R302.simulate()
#     u.R302.specifications[0]()
#     return u.R302.effluent_titer - u.R302.titer_to_load


# def load_titer_with_glucose(titer_to_load, set_F301_V=None):
#     # clear_units([V301, K301])
#     F301_ub = 0.8
#     F301_lb = 0. if set_F301_V is None else set_F301_V
#     M304_lb, M304_ub = 0., 100_000.  # for low-titer high-yield combinations, if infeasible, use a higher upper bound
    
#     spec.spec_2 = titer_to_load
#     u.R302.titer_to_load = titer_to_load
#     F301_titer_obj_fn(F301_lb)
    
#     if M304_titer_obj_fn(M304_lb) < 0.: # if there is too low a conc even with no dilution
#         # IQ_interpolation(F301_titer_obj_fn, F301_lb, F301_ub, ytol=1e-3)
#         spec.titer_inhibitor_specification.check_sugar_concentration() 
#     # elif F301_titer_obj_fn(1e-4)>0: # if the slightest evaporation results in too high a conc
#     elif M304_titer_obj_fn(M304_ub) > 0.:
#         IQ_interpolation(M304_titer_obj_fn, 
#                          M304_lb,
#                          M304_ub, 
#                          ytol=1e-3)
#     else:
#         F301_titer_obj_fn(F301_lb)
#         IQ_interpolation(M304_titer_obj_fn, 
#                          M304_lb, 
#                          M304_ub, 
#                          ytol=1e-3)
        
#     if not feedstock_ID=='Glucose' and (set_F301_V is None): spec.titer_inhibitor_specification.check_sugar_concentration()
    
# spec.load_spec_1 = spec.load_yield
# spec.load_spec_3 = spec.load_productivity
# spec.load_spec_2 = load_titer_with_glucose





# # %% Diagram
# import biosteam as bst
# bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = True
# succinic_sys.diagram('cluster')

# #%% Define unit groups

# area_names = [
#     'feedstock',
#     'conversion',
#     'separation',
#     # 'recovery',
#     'wastewater system',
#     'storage',
#     'boiler & turbogenerator',
#     'cooling utility facilities',
#     'other facilities',
#     'heat exchanger network',
# ]

# # area_names = [
# #     'feedstock',
# #     'conversion',
# #     'separation',
# #     'wastewater',
# #     'storage',
# #     'boiler & turbogenerator',
# #     'cooling utility facilities',
# #     'other facilities',
# #     'heat exchanger network',
# # ]
# # I changed area_name 'wastewater in the system glucose --> recovery 
# unit_groups = bst.UnitGroup.group_by_area(succinic_sys.units)


# unit_groups.append(bst.UnitGroup('natural gas (for steam generation)'))
# unit_groups.append(bst.UnitGroup('natural gas (for product drying)'))

# unit_groups.append(bst.UnitGroup('fixed operating costs'))


# for i, j in zip(unit_groups, area_names): i.name = j
# for i in unit_groups: i.autofill_metrics(shorthand=False, 
#                                          electricity_production=False, 
#                                          electricity_consumption=True,
#                                          material_cost=True)
# for i in unit_groups:
#     if i.name == 'storage' or i.name=='other facilities' or i.name == 'cooling utility facilities':
#         i.metrics[-1].getter = lambda: 0. # Material cost
#     if i.name == 'cooling utility facilities':
#         i.metrics[1].getter = lambda: 0. # Cooling duty
#     if i.name == 'boiler & turbogenerator':
#         i.metrics[-1] = bst.evaluation.Metric('Material cost', 
#                                             getter=lambda: succinic_tea.utility_cost/succinic_tea.operating_hours, 
#                                             units='USD/hr',
#                                             element=None) # Material cost
#         # i.metrics[-2].getter = lambda: BT.power_utility.rate/1e3 # Electricity consumption [MW]
        
# for HXN_group in unit_groups:
#     if HXN_group.name == 'heat exchanger network':
#         HXN_group.filter_savings = False
#         # assert isinstance(HXN_group.units[0], HeatExchangerNetwork) #!!! Note groups aren't working for non-sugarcane succinic acid biorefineries, fix if needed


# unit_groups[-3].metrics[-1] = bst.evaluation.Metric('Material cost', 
#                                                     getter=lambda: BT.natural_gas_price * BT.natural_gas.F_mass, 
#                                                     units='USD/hr',
#                                                     element=None)

# unit_groups[-2].metrics[-1] = bst.evaluation.Metric('Material cost', 
#                                                     getter=lambda: u.F404.ins[2].cost, 
#                                                     units='USD/hr',
#                                                     element=None)

# unit_groups[-1].metrics[-1] = bst.evaluation.Metric('Material cost', 
#                                                     getter=lambda: succinic_tea.FOC/succinic_tea.operating_hours, 
#                                                     units='USD/hr',
#                                                     element=None)

# for i in unit_groups:
#     i.metric(i.get_net_electricity_production,
#                 'Net electricity production',
#                 'kW')




# unit_groups_dict = {}
# for i in unit_groups:
#     unit_groups_dict[i.name] = i

# BT = u.BT701
# CT = u.CT801
# CWP = u.CWP802
# HXN = u.HXN1001

# # HXN1001.force_ideal_thermo = True


# # %%
# # =============================================================================
# # Simulate system and get results
# # =============================================================================

# num_sims = 3
# num_solve_tea = 3
# def get_product_stream_MPSP():
#     for i in range(num_sims):
#         succinic_sys.simulate()
#     for i in range(num_solve_tea):
#         product_stream.price = succinic_tea.solve_price(product_stream)
#         # product_stream.price = 2.7 
#         # to check the NPV value
        
        
#     return product_stream.price * product_stream.F_mass / product_stream.imass['SuccinicAcid']

# def simulate_and_print():
#     MPSP = get_product_stream_MPSP()
#     GWP = succinic_LCA.GWP
#     FEC = succinic_LCA.FEC
#     print('\n---------- Simulation Results ----------')
#     print(f'MPSP is ${MPSP:.3f}/kg')
#     print(f'GWP100a is {GWP:.3f} kg CO2-eq./kg')
#     print(f'FEC is {FEC:.3f} MJ/kg')
#     print('----------------------------------------\n')

# get_product_stream_MPSP()
# succinic_tea.labor_cost=3212962*get_flow_dry_tpd()/2205
# get_product_stream_MPSP()




# # # #%% LCA

# succinic_LCA = SuccinicLCA(system=succinic_sys, 
#                  CFs=CFs, 
#                  feedstock=s.glucose_feedstock, 
#                  feedstock_ID=feedstock_ID,
#                  input_biogenic_carbon_streams=[feedstock, s.CSL],
#                  main_product=product_stream, 
#                  main_product_chemical_IDs=['SuccinicAcid',], 
#                  by_products=[], 
#                  cooling_tower=u.CT801, 
#                  chilled_water_processing_units=[u.CWP802,], 
#                  boiler=u.BT701, has_turbogenerator=True,
#                  add_EOL_GWP=True,
#                  )

# #%% Simulate and print
# spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)
# get_product_stream_MPSP()

# #%% TEA breakdown (old)

# def TEA_breakdown(print_output=False, fractions=False):
#     metric_breakdowns = {i.name: {} for i in unit_groups[0].metrics}
#     for ug in unit_groups:
#         for metric in ug.metrics:
#             denominator = 1.
#             if fractions:
#                 metric_name = metric.name
#                 if metric_name=='Inst. eq. cost':
#                     denominator = succinic_tea.installed_equipment_cost/1e6
#                 elif metric_name=='Cooling':
#                     denominator=0
#                     for i in list(CT.heat_utilities) + list(u.CWP802.heat_utilities):
#                         if i.flow<0 and i.duty>0:
#                             denominator += i.duty/1e6
#                 elif metric_name=='Heating':
#                     for i in list(BT.heat_utilities):
#                         if i.flow<0 and i.duty<0:
#                             denominator -= i.duty/1e6
#                 elif metric_name=='Elec. cons.':
#                     denominator = succinic_sys.power_utility.consumption/1e3
#                 elif metric_name=='Elec. prod.':
#                     denominator = succinic_sys.power_utility.production/1e3
#                 elif metric_name=='Mat. cost':
#                     denominator = (succinic_tea.material_cost/succinic_tea.operating_hours) + BT.natural_gas.F_mass*BT.natural_gas_price
#             # storage_metric_val = None
#             if not ug.name=='storage':
#                 if ug.name=='other facilities':
#                     metric_breakdowns[metric.name]['storage and ' + ug.name] = (metric() + unit_groups_dict['storage'].metrics[ug.metrics.index(metric)]())/denominator
#                 else:
#                     metric_breakdowns[metric.name][ug.name] = metric()/denominator
                    
                    
#             if ug.name=='natural gas (for steam generation)':
#                 if metric.name=='Mat. cost' or metric.name=='Material cost':
#                     metric_breakdowns[metric.name][ug.name] = BT.natural_gas.F_mass*BT.natural_gas_price/denominator
            
            
#             # else:
#             #     storage_metric_val = metric()
                
#     # print and return metric_breakdowns
#     if print_output:
#         for i in unit_groups[0].metrics:
#             print(f"\n\n----- {i.name} ({i.units}) -----")
#             metric_breakdowns_i = metric_breakdowns[i.name]
#             for j in metric_breakdowns_i.keys():
#                 print(f"{j}: {format(metric_breakdowns_i[j], '.3f')}")
#     return metric_breakdowns

# TEA_breakdown(print_output=True, fractions=True)# -*- coding: utf-8 -*-


# # %% Carbon balance

# def get_carbon_flow_by_area():
#     c_flow_units = []
#     for ug in unit_groups:
#         sys_ID = ''.join(e for e in ug.name if e.isalnum())
#         ug_sys = ug.to_system(sys_ID)
#         ug_ins = [stream for stream in ug_sys.streams if stream.source not in ug_sys.units]
#         ug_outs = [stream for stream in ug_sys.streams if stream.sink not in ug_sys.units]
#         if sum([stream.get_atomic_flow('C') for stream in ug_ins])>0 or sum([stream.get_atomic_flow('C') for stream in ug_outs])>0:
#             c_flow_units.append(ug_sys)
#     # print(c_flow_units)
#     source=[]
#     target=[]
#     value=[]
#     carbon_feeds = [stream for stream in succinic_sys.feeds if stream.get_atomic_flow('C')>0]
#     label=[unit.ID  + ' process' for unit in c_flow_units] + [stream.ID for stream in carbon_feeds] + ['SuccinicAcid'] + ['Emissions']
#     for unit in c_flow_units:
#         ug_outs = [stream for stream in unit.streams if stream.sink not in unit.units]
#         # print(ug_outs)
#         for stream in ug_outs:
#             if stream.get_atomic_flow('C'):
#                 if stream.sink:
#                     for ug in c_flow_units:
#                         if stream.sink in ug.units:
#                             target.append(c_flow_units.index(ug))
#                     source.append(c_flow_units.index(unit))
#                     value.append(stream.get_atomic_flow('C'))
#                 else:
#                     if not stream.ID==product_stream.ID:
#                         source.append(c_flow_units.index(unit))
#                         target.append(len(label)-1)
#                         value.append(stream.get_atomic_flow('C'))
#                     else:
#                         source.append(c_flow_units.index(unit))
#                         target.append(len(label)-2)
#                         value.append(stream.get_atomic_flow('C'))
#     for stream in carbon_feeds:
#         source.append(label.index(stream.ID))
#         for ug in c_flow_units:
#             if stream.sink in ug.units:
#                 target.append(c_flow_units.index(ug))
#         value.append(stream.get_atomic_flow('C'))
#     return source, target, value, label


# def get_carbon_flow():
#     c_flow_units = []
#     for unit in u:
#         if sum([stream.get_atomic_flow('C') for stream in unit.ins])>0 or sum([stream.get_atomic_flow('C') for stream in unit.outs])>0:
#             c_flow_units.append(unit)
#     source=[]
#     target=[]
#     value=[]
#     carbon_feeds = [stream for stream in succinic_sys.feeds if stream.get_atomic_flow('C')>0]
#     label=[unit.ID  + ': ' + unit.__class__.__name__ for unit in c_flow_units] + [stream.ID for stream in carbon_feeds] + ['SuccinicAcid'] + ['Emissions']
#     for unit in c_flow_units:
#         for stream in unit.outs:
#             if stream.get_atomic_flow('C'):
#                 if stream.sink:
#                     source.append(c_flow_units.index(unit))
#                     target.append(c_flow_units.index(stream.sink))
#                     value.append(stream.get_atomic_flow('C'))
#                 else:
#                     if not stream.ID==product_stream.ID:
#                         source.append(c_flow_units.index(unit))
#                         target.append(len(label)-1)
#                         value.append(stream.get_atomic_flow('C'))
#                     else:
#                         source.append(c_flow_units.index(unit))
#                         target.append(len(label)-2)
#                         value.append(stream.get_atomic_flow('C'))
#     for stream in carbon_feeds:
#         source.append(label.index(stream.ID))
#         target.append(c_flow_units.index(stream.sink))
#         value.append(stream.get_atomic_flow('C'))
#     return source, target, value, label

# def plot_carbon_flow():
#     import plotly.graph_objs as go#create links
#     from plotly.offline import init_notebook_mode,  plot
#     init_notebook_mode()
#     source, target, value, label = get_carbon_flow_by_area()
#     link = dict(source=source, target=target, value=value, 
#     color=['rgba(0,0,0, 0.5)'] * len(source))#create nodes
#     node = dict(label=label, pad=15, thickness=8, 
#                 color = ['blue' if ':' in label[i] else 'gray' for i in range(len(label))],
#                )#create a sankey object
#     chart = go.Sankey(link=link, node=node, arrangement="freeform")#build a figure
#     fig = go.Figure(chart)
#     plot(fig)
#     fig.show()
#     fig.write_image("fig1.png")



# # %% Plot figures
# plot = True
# if plot:
    
#     ## TEA breakdown figure
#     df_TEA_breakdown = bst.UnitGroup.df_from_groups(
#         unit_groups, fraction=True,
#         scale_fractions_to_positive_values=False
#     )
    
    
#     # df_TEA_breakdown['Net electricity production']*=-1
#     # df_TEA_breakdown = df_TEA_breakdown.rename(columns={'Net electricity production': 'Net electricity demand'})
    
    
#     stacked_bar_plot(dataframe=df_TEA_breakdown, 
#                      # y_ticks=[-200, -175, -150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175], 
#                      y_ticks=[-125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175, 200, 225], 
#                      y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
#                      y_units = "%", 
#                      yticks_fontsize=10,xticks_fontsize=10,
                     
#                      colors=['#7BBD84', '#F7C652', '#63C6CE', '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#F8858A', 'red', 'magenta'],
#                      filename='TEA_breakdown_stacked_bar_plot')
    