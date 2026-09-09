#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:17:56 2026

@author: princyk2
"""

#TODO after drying might need heat exchanger when cools to 70 ° C
#TODO find out the how water split during extraction
import biosteam as bst
from biosteam import preferences
from biorefineries.cellulosic import create_facilities
from biosteam import SystemFactory
from biorefineries.covercress import units_oil 
from biorefineries.covercress.chemicals_data import create_covercress_chemicals
from biorefineries.covercress.process_settings import price, CFs, chem_index, _GDP_2007_to_2010,_lps
from biosteam.units import LLEUnit
from thermosteam import indexer, equilibrium, settings
__all__=('create_covercress_oil_sys')

chemicals = create_covercress_chemicals()
bst.settings.set_thermo(chemicals)
T_PREHEAT_K = 35.0 + 273.15 #-refer[4]
T_DRYER_FLAKES_K = 85.0 + 273.15  #-refer[4] design target; align with DrumDryer API
debris_frac = 0.015 #-refer[4]
hexane_mass_ratio = 1.5*0.66 #TODO (update)check exact hexane ratio  from the refer[ppt] assume 3:1 v/w for 1 Kg 3 L of hexane required and at 25 °C density of hexane is 0.66 Kg/L thus the solvent ratio from the formula--> solvent required in Kg= feed mass*3L/Kg*density of hexane
# Extraction specs
solvent_ID = 'Hexane'
# Mass fraction of *lipids* in miscella (miscella = lipids + solvent in this model).
# 0.25 => 25 % "oil" in miscella, 75 % hexane by mass.
miscella_lipid_mass_frac = 0.25
# Cake-side extraction parameter (SolventExtraction.residual_oil_mass_fraction).
residual_oil_in_cake = 0.01
# Wet-cake hexane holdup (SolventExtraction.solvent_mass_fraction_in_cake); set None to disable.
solvent_mass_frac_in_cake = 0.30
#TODO temp set is 300 ° k for inlet covercress but later update from 110 C to 300 K by setting heat utility
total_feed = 100_00.0          # kg/hr total at gate

feed = bst.Stream(
    'flake_feed',
    covercress=total_feed * (1 - debris_frac),   # 98,500 kg/hr
    Debris=total_feed * debris_frac,
    units='kg/hr',
    T=300.0,
    P=101325.0,
)
solids_retentate={'Glucose','Protein','Ash','Lignin','Sucrose','Cellulose','Hemicellulose'}
lipids = (
     'Lipid', 'Oil',
     'PL', 'Phospholipid',
     'FFA', 'OleicAcid',
     'MAG', 'MonoOlein', 'DAG', 'DiOlein', 'TAG', 'TriOlein','Tristearin',
     'LinoleicAcid', 'LinolenicAcid','Stericacid','OOO', 'LLL', 'LnLnLn', 'SSS'
 )

def screw_press_permeate(chemicals):
    # Build a split dictionary which return the rententate split fraction
    lipid_permeate = 0.7
    water_permeate = 0.2 #TODO later can be changed the fraction of water goes into the oil
    
   
    s={k:1 for k in solids_retentate if k in chemicals}
    print("dictionary of all chemicl")
    for k in solids_retentate:
        if k in chemicals:
            s[k]=0.0   
    for ID in lipids:
        if ID in chemicals:
            s[ID] = lipid_permeate
        if 'Water' in chemicals:
            s['Water'] = water_permeate
        elif 'H2O' in chemicals:
            s['H2O'] = water_permeate
    
    return s
@SystemFactory(ID='oil_extraction_sys', ins=['raw_seed'], outs=[
        dict(ID='crude_oil', units='kg/hr', price=price['crude_oil']),
        dict(ID='protein_powder', units='kg/hr', price=price['protein_isolate']),
    ],)
# @SystemFactory(ID='oil_extraction_sys', ins=['raw_seed'], )
def create_covercress_oil_sys(ins,outs,moisture_content=0.03):# maxmium oil recovery obtained when the cooked seeds was at 3 %
    raw_seed, = ins
    crude_oil, protein_isolate = outs 
    # creating a function which returns the output of screw press ie:cake and liquid
    #TODO check for %wrt to seed rather than direct number.
    # seed_mass = 985.0   # kg/hr
    # foreign_mass = 15.0 # kg/hr  (1.5% of 1000)
    # T_DRYER_FLAKES_K= 273+85
    # raw_seed = bst.Stream(
    #     'raw_seed',
    #     covercress=seed_mass,
    #     Debris=foreign_mass,
    #     units='kg/hr',
        
    # )
    # raw_seed,= ins
    # miscella_combined=outs
    chemicals=feed.chemicals
    split=screw_press_permeate(chemicals)
    
    
    # crude_oil = bst.Stream('crude_oil', units='kg/hr', price=price['crude_oil'])
    # protein_isolate   = bst.Stream('protein_isolate', units='kg/hr', price=price['protein_isolate'])  
    C101 = units_oil.ConveyingBelt('C101', ins=raw_seed, outs='to_screen')
    U102 = units_oil.VibratingScreen(
        'U102', ins=C101-0, outs=('screened_seed', 'overs'),
        split= 1,  # set per screening assumption assumed here =1
    )
    U103 = units_oil.SeedCleaningSystem(
        'U103',
        ins=U102-0,
        outs=('clean_seed', 'cleaning_reject'),
        split={'covercress': 1.0, 'Debris': 0.0},  # refre[4]
    )
    
    U104 = units_oil.Cooking(
        'U104',
        ins=(U103-0, 'conditioner_air', 'conditioner_natural_gas'),
        outs=('conditioned_seed', 'conditioner_outlet_air', 'conditioner_emissions'),
        moisture_content=0.09,
        split=0.0,
        T=273+35,
        utility_agent='Steam',
    )#refer 4
    #refer [6] and refer [2] in the canola processing then operation wired from conditioninf->flaking>cooking, 
    # conditionng -->soften the outer shell and moisture levels for optimal preparation performance.in refer[2] the MC durinf storage is 9 % so I assumed once after seed cleaning we adjust the moisture to 9 % and then after cooking final MC reduced to 3-4 % refer[2]
    U105=units_oil.SeedFlaking("U105", ins=U104-0, outs= "seed_flakes")
  
    U106= units_oil.SeedConditioning(
        'U106',
        ins=(U105-0, 'dryer_air', 'dryer_natural_gas'),
        outs=('flakes_dried', 'dryer_outlet_air', 'dryer_emissions'),
        moisture_content=0.03,
        split=0.0,
        T=273+85,
        utility_agent='Steam',
    ) #refer4
    #TODO check later where to direct emsiions from dryer 
    U107 = bst.ScrewPress('U107', U106-0,outs=('pressed_oil','press_cake'), 
                         split=split)
    H101 = bst.HXutility(
        'H101',
        ins=U107-1,
        outs='press_cake_cooled',
        T=60+273,
        
    )
    C102 = units_oil.ConveyingBelt(
        'C102',
        ins=H101-0,
        outs='cake_to_extraction',
    )
    
    # solvent_recycle = bst.Stream('Hexane_recycle', phase='l', T=273+60, P=101325)
    fresh_solvent =bst.Stream('Hexane_1', phase='l',units='kg/hr',  price=price['Hexane'])
    
    #TODO temp refer[3] for extraction is 63 °C update it accordingly
    # Split recovered hexane between the two extractors (tune 0.5 / 0.5 or by demand).
    # S401 = bst.Splitter(
    #     'S401',
    #     ins=solvent_recycle,
    #     outs=('hexane_recycle_1', 'hexane_recycle_2'),
    #     split=0.5,
    # )
    T101 = units_oil.HexaneStorageTank('T101', ins=fresh_solvent, outs='hexane_from_storage')
    # @T101.add_specification(run=False)
    # def T101_hexane_makeup():
        
    #     cake = C102.outs[0]  # extractor cake feed
    #     solvent_to_cake = 1  # kg hexane / kg cake refer[14] and also refer [20]
    #     H_need = solvent_to_cake * cake.F_mass
    #     recycle_solvent = M401.ins[0].imass['Hexane']
    #     cake_hex = cake.imass['Hexane']
    #     fresh = max(0.0, H_need - recycle_solvent - cake_hex)
    #     T101.ins[0].imass['Hexane'] = fresh
    #     T101._run()
    #     H102._run()

 

    # M401 = bst.Mixer('M401', ins=('solvent_recycle', H102-0), outs=('to_M102',))

    # @M401.add_specification(run=False)
    # def M401_mix_hexane():
    #     """Mix recycle + hot makeup; hexane mass is set upstream on T101.ins[0]."""
     
            
    #     cake = C102.outs[0]  # extractor cake feed
    #     solvent_to_cake = 1  # kg hexane / kg cake refer[14] and also refer [20]
    #     H_need = solvent_to_cake * cake.F_mass
    #     recycle=M401.ins[0].imass['Hexane']
    #     cake_hex = cake.imass['Hexane']
    #     fresh = max(0.0, H_need - recycle  - cake_hex)
    #     T101.ins[0].imass['Hexane'] = fresh
    #     # T101._run()
    #     # H102._run()
        # M401._run()
    M102=bst.Mixer('M102', ins=(C102-0,T101-0,'' ),outs= 'U201_feed')
    @M102.add_specification(run=False)
    def adjust_solvent():
        feed, fresh_solvent,solvent_recycle = M102.ins
        required_solvent = feed.F_mass * 1
        recycled_solvent = M102.ins[2].imass['Hexane']
        if recycled_solvent > required_solvent:
            M102.ins[2].imass['Hexane'] = required_solvent
            solvent_recycle.F_mass = required_solvent
            recycled_solvent = solvent_recycle.imass['Hexane']
        fresh_solvent.imass['Hexane'] = max(0, required_solvent - recycled_solvent)
        M102._run()
#     # M102=bst.Mixer('M102', ins=(fresh_solvent,'' ),outs= 'U201_hexane_feed')
#     # @M102.add_specification
#     # def adjust_solvent():
#     #     fresh_solvent,solvent_recycle = M102.ins
#     #     required_solvent = feed.F_mass * 1
#     #     recycled_solvent = M102.ins[1].imass['Hexane']
#     #     if recycled_solvent > required_solvent:
#     #         solvent_recycle.F_mass = required_solvent
#     #         recycled_solvent = solvent_recycle.imass['Hexane']
#     #     fresh_solvent.imass['Hexane'] = max(0, required_solvent - recycled_solvent)
#     #     M102._run()
#     # M102_C=bst.Mixer('M102_C', ins=(C102-0, M102-0))
    H102 = bst.HXutility(
           'H102',
           ins=M102-0,
           outs='U201_hot_feed',
           T=60+273,
       )
    # =============================================================================
    #OIL extraction
    # =============================================================================
    U201 = units_oil.CovercressExtractor(
        'U201',
        ins=H102-0,
        outs=('miscella_2', 'cake_2'),
        top_chemical='Hexane',
    )
    # E401=bst.MultiEffectEvaporator('E401',ins=U201-0, outs= ('oil', 'Hexane'),V=0.1,V_definition='First-effect',
    #                              P=(101325, 73581, 50892, 32777, 20000)) 
        
        

#     # makeup_water = bst.Stream('makeup_water', phase='l', T=95 + 273.15)
#     # M103=bst.Mixer('M103',ins=(U201-1, makeup_water),outs='F101_feed')
#     # @M103.add_specification(run=False)
#     # def M103_adjust_moisture():
#     #     MC_present=U201.outs[1].imass['Water']
#     #     # M103.ins[1].imass['Water'] = max(0,0.43 * (U201.outs[1].F_mass - 0.99*U201.outs[1].imass['Hexane']) - 1.42 * MC_present)
#     #     M103.ins[1].imass['Water'] = max(0,0.28 * (U201.outs[1].F_mass - 0.99*U201.outs[1].imass['Hexane']) - 1.28 * MC_present)
#     #     M103._run()
          
    F101 = units_oil.Desolventizer(
           'F101',
           ins=(U201-1, 'dryer_gas', 'Steam'),
           outs=('toasted_meal', 'desolventizer_vapor', 'emissions'),
           split={solvent_ID: 0.99},
           moisture_content=0.01,
           T=110.0 + 273.15,
           utility_agent='Steam',
       )
    # F102=units_oil.Drying('F102', ins=(F101-0,'dryer_gas', 'Steam'), outs=('dried_meal','hot_air', 'emissions'),split={solvent_ID: 0.99},
    # moisture_content=0.065,
    # T=110.0 + 273.15,
    # utility_agent='Steam')

    H103 = units_oil.MealCooler(
        'H103',
        ins=F101-0,                   
        outs='meal_cooled',
        T=30.0 + 273.15,cool_only=True,
    ) 
    
    E101=bst.MultiEffectEvaporator('E101', ins=U201-0, outs=('E101_l','E101_g'), P=(101325, 73581, 50892, 32777, 20000), V=0.5, V_definition='First-effect',chemical='Hexane',flash=False)
    E101.target_conc=0.98 #refer[20]
    @E101.add_bounded_numerical_specification(x0=1e-3, x1=1.0 - 1e-3, xtol=1e-4, ytol=1e-4, x=0.5)
    def E101_V(V):
        E101.V=V
        E101._run()
        out_mass = E101.outs[0]
        mass_oil =E101.outs[0].imass['OOO']+E101.outs[0].imass['LLL']+ E101.outs[0].imass['LnLnLn']+E101.outs[0].imass['SSS']
        return mass_oil/out_mass.F_mass-E101.target_conc
    # E101.run_after_specifications = True
    E101_P0=bst.Pump('E101_P0',ins=E101-0,P=101325,)
    E101_P1=bst.Pump('E101_P1',ins=E101-1,P=101325,)
    steam = bst.Stream(
    'steam',
     Water=10,            
    phase='g',
    T=110+273,             # 110 °C, or up to 413.15 for 140 °C used low pressure steam
                 # match E102.P (150 mbar)
    units='kg/hr',
    price=_lps.regeneration_price,   # from process setting
 )
    import numpy as np
    K_strip = np.array([
    # 2, #refer https://pubs.acs.org/doi/10.1021/acsomega.1c07044
    27, 
    #27,   # Hexane  strongly vapor refer [22] for activity coefficent and then calculated the partition coefficient
    11.56, #NBwhen i fix the K no need to bound specs for steam flowrate so now stripper does the job
    #141,   # Water   strongly vapor #TODO refernce required waterhas theoreticaly high activity coefficnet for imscible mixture  and following fatty acids has very 
    1e-5,   # OOO  strongly liquid
    1e-5,   # LLL
    1e-5,   # LnLnLn # Tb of Linolenic is less than all TAG which is greater than >600 , reason why Linolenic is identified in the overhead
    1e-5,   # SSS
])
    E102 = bst.Stripper(
    'E102',
     N_stages=7,# refer[23]        
    ins=[E101_P0-0,steam], # refer 18
    outs=['strip_overhead', 'stripped_oil'],#solute =water because already majority of hexane removed uisng MEE so other majority component is water that must removed so we use solute =water to calculate the internal calcs like stage efficiency etc also also its boiling boint is >hexane thus it would be partitoned in significant amount is both phases  
    P=6000,T=110+273, #110 as per literature and activity coefficent for 110 °C is 0.572 
    solute="Water",  partition_data={
          'IDs': ['Hexane', 'Water', 'OOO', 'LLL', 'LnLnLn', 'SSS'],'K':K_strip,
        'vapor_chemicals': ['Hexane', 'Water'],
          'liquid_chemicals': ['OOO', 'LLL', 'LnLnLn', 'SSS'],
      },
)    
    E102.x_hex_target=0.003 #refer[excel] 0.3 wt % residual hexane
    # @E102.add_specification(run=True)
    # def adjust_steam():
    #     feed, strip_gas = E102.ins
    #     strip_gas.phase='g'
    #     feed.phase='l'
    #     steam.imass['Water'] = 1*E102.ins[0].imass['Hexane']   # kg/hrTODO refer [21] but later update why we assume 1:1 because there is only traces of hexane present which leads to the required MC -0.1-0.3 %
    
    # @E102.add_bounded_numerical_specification( # why we use this specs because when we increase the feedflowrate incoming steam is not enough to remove traces of hexane to maintain 0.3 wt % in the stripped oil but when i fix the K now this bounded specs is not required
    #       x0=0.0,        # minimum steam kg/hr
    #       x1=100,      # maximum steam kg/hr (increase if needed)
    #       x=10,         # initial guess
    #       xtol=1e-3,
    #       ytol=1e-4,
    #   )
    # def E102_steam_to_remove_hexane(Wsteam):
    #     feed, steam = E102.ins
       
    #     steam.imass['Water'] = Wsteam
    #     steam.phase = 'g'
    #     steam.P = E102.P
    #     E102._run()
    #     oil = E102.outs[1]
    #     x_hex = oil.imass['Hexane']/oil.F_mass
    #     return x_hex - E102.x_hex_target
    @E102.add_specification(run=False)
    def E102_oil_to_product():
        E102._run()
        vapor, oil = E102.outs
        lipids = ('OOO', 'LLL', 'LnLnLn', 'SSS')
        vapor.phase = 'g'
        oil.phase = 'l'
        for ID in lipids:
            oil.imass[ID] += vapor.imass[ID]
            vapor.imass[ID] = 0. # why we force all oil in the vpor is 0 because it found linolenic acid is present in the overhead because we treat TAG as freefatty acids so it has low boiling point compare to TAG (to match the real world vacuum stripping of miscelle)
        # P = vapor.P
        # vapor.vle(P=P, V=1.0)
        # oil.vle(P=P, V=0.0)  # calculates vle at given V and P refer tutorial VLE
    E102_P0=bst.Pump('E102_P0',ins=E102-0,P=101325)
    E102_P1=bst.Pump('E102_P1',ins=E102-1,P=101325)
    # =============================================================================
    #Solvent Recovery
    # =============================================================================
        
#     # E102 = bst.Flash('E102', ins=E101-0, outs=('E102_vap', 'stripped_oil'), P=15_000, T=383.15,vacuum_system_preference='Steam-jet ejector')
    M104=bst.Mixer('M104', ins=(F101-1,E101_P1-0,E102_P0-0),outs='hexane_vapour_with_water_and_noncondensates')
    
    # H104=bst.HXutility('H104',ins=M104-0,outs='cooled_hexane',T=40+273,)
    F104 = bst.Flash('F104', ins=M104-0, outs=('vent_gas', 'condensate'), P=101325,T=273+50) #TODOshould we recycle this vent gas ? and check for appropriate units because in industry they use minerla oil to abosorb the gases
     # flash drum hexane requires more vapour pressure than water to reach equilibrium that’s reason hexane in liquid. The Psat for water is <P sat for hexane so to attain equilibrium more hexane vapours isrequired thus water attain equilibrium faster with less Psat compared to hexane (at equili, Psat =vapour press)						
    # H105=bst.HXutility('H105',ins=F104-1,outs='cooled_hexane',T=60+273)
    F104_P=bst.Pump('F104_P',ins=F104-1,outs='solvent_recycle')
    F104_P-0-2-M102
    # the water hexane mixture is very immiscible so it has very large partiion coeffiecnt thus it leaves all water molecules in the vapour state
    # D104 = bst.LLEUnit('D104', ins=F104-1,
    #                    outs=('water_rich', solvent_recycle),  
    #                    top_chemical='Hexane')
    
    
   
   
    # =============================================================================
    #Protein extraction feeed-dried cake
    # =============================================================================
    #protein extraction unit -->2-stage single-solvent extraction to manufacture the actual protein product refer[27]
    
    G201=units_oil.HammerMill('G201', ins=H103-0,outs='grinded_meal')
    fresh_nacl  = bst.Stream('fresh_nacl',  phase='s', units='kg/hr',price=price['NaCl'])   # 
    fresh_NaOH=bst.Stream('fresh_NaOH',units='kg/hr',price=price['NaOH'])
    makeup_water_1 = bst.Stream('makeup_water_1', phase='l',  units='kg/hr', price=price['Makeup water'])
    M105 = bst.Mixer('M105', ins=(makeup_water_1, fresh_nacl), outs='0.1M_NaCl_extraction')
    @M105.add_specification(run=False)
    def M105_make_saline():
        water, nacl = M105.ins
        meal=G201.outs[0]
        solid = meal.F_mass 
        vol_L = 0.01 * solid           # 1 g : 10 mL  ->  0.006 g naCLfor 1 L solution refer[27] and also 10 L/kg = 0.01 m3/kg  -> m3/hr
        nacl.imass['NaCl']  = 0.006 * vol_L     # 0.1 M
        water.imass['Water'] = 1.0 * vol_L - nacl.F_mass   # ~1 kg/L solution density
        M105._run()
    
    
    
    H104=bst.HXutility('H104',ins=M105-0,outs='heated_NaCl',T=50+273)
   
    U202=bst.MixerSettler('U202', ins=(G201-0,H104-0), outs=('protein_extract', 'spent_solids'), mixer_data={'tau': 2.0,               # 2 h stir
                'agitator_kW_per_m3': 1.0,},model='split',settler_data={
                    'split': {
                        'Protein': 0.47,              # 47.0% of feed protein -> supernatant_1, protein recovery into supernatant (set from TT8/Y1126) from refr [27]recovered proterin from the TT8 variety is 50.1 % and also from this 15 fold of extracted protein entered into first stage
                        'Water': 0.90, 'NaCl': 0.90,  #TODO check water splits into extract
                        'Glucose': 0.90, 'Sucrose': 0.90,#TODO checksplits into extract refer[26] - some amount of carboyhydraes are present and also the extract contains water, saline ixture ieNaCl and some glucose and sucrose since it is water and saline soluble
                        'Cellulose': 0.0, 'Hemicellulose': 0.0,
                        'Lignin': 0.0, 'Ash': 0.0,    # insolubles -> spent solids
                        'OOO': 0.01, 'LLL': 0.01,
                        'LnLnLn': 0.01, 'SSS': 0.01,
                    },
                },
            )#refer 24
#refer tutorial mixer settler with model=split does not use lle its just replicates the centrifuge
    # M502=bst.Mixer('M502', ins=(U202-0, U202-1),outs='slurry')
    M110=bst.Mixer('M110',ins=(U202-0,U202-1),outs='mixed_slurry_centrifuge_1')
    U203 = bst.SolidsCentrifuge(
       'U203',
       ins=(M110-0),                                 # after the extraction the output is in th eform of slurry which then goes into centrifugation and separate into solids
       outs=('spent_solids', 'protein_supernatant'), # [0] solids-rich, [1] liquid-rich
       split={
           # insolubles -> spent solids
           'Cellulose': 1.0, 'Hemicellulose': 1.0, 'Lignin': 1.0, 'Ash': 1.0,
           # protein: 31% retained in spent solids => 69% to supernatant refer[27]
           'Protein': 0.53,          
           'NaCl': 0, 'NaOH': 0, 'Glucose': 0.0, 'Sucrose': 0.0, 'Hexane': 0.0,
           'OOO': 0.0, 'LLL': 0.0, 'LnLnLn': 0.0, 'SSS': 0.0,
          
       },
       moisture_content=0.01,   #TODO check the moisture contetnt                    # wet-cake moisture 
       solids=('Cellulose', 'Hemicellulose', 'Lignin', 'Ash', 'Protein'),
       centrifuge_type='scroll_solid_bowl',
   )
    
    
    
    extraction_water = bst.Stream('extraction_water', phase='l', units='kg/hr', price=price['Makeup water'])# TODO later update the price of fresh solvent
    # M107=bst.Mixer('M107',ins=(extraction_water)) #TODO if recycle needed
    # @M107.add_specification(run=False)
    # def M107_make_wash():
    #     wextraction_water = M107.ins
    #     pellet=U203.outs[0].F_mass
    #     water_required = 0.01 * pellet       # 1:10 ratio (use pellet.F_mass - imass['Water'] for dry basis) 10 L/kg = 0.01 m3/kg  -> m3/hr

    #      # PURE water (no NaCl in the wash)
    #                   # pre-heated to 50 C
        
    #     # extraction_water= max(0,water_required)#TODO if we need to recylce
    #     M107._run()
                   
    H105=bst.HXutility('H105',ins= extraction_water,outs='heated_water',T=50+273)
    @H105.add_specification(run=False)
    def H105_make_wash():
        
        pellet=U203.outs[0].F_mass
        H105.ins[0].imass['Water'] = 0.01 * pellet      # 1:10 ratio (use pellet.F_mass - imass['Water'] for dry basis) 10 L/kg = 0.01 m3/kg  -> m3/hr
        H105._run()
        # U204.ins[1].imass['Water'] = 1.0 * vol_L 
    U204 = bst.MixerSettler(
        'U204',
        ins=(U203-0, H105-0),                     # [0] pellet from stage 1, [1] fresh water
        outs=('supernatant_2', 'spent_solids_final'), # [0] extract (collected), [1] raffinate (solids)
        model='split',
        mixer_data={'tau': 1.0, 'agitator_kW_per_m3': 1.0},   # 1 h stir; 250 rpm -> specific power
        settler_data={'split': {
            'Protein': 0.059,                            # # 5.9% of the entering (pellet) protein -> supernatant_2 also refer word file-protein extraction
            'Water':0.90,'NaCl':0.90,'NaOH':0.90,'Glucose':0.90,'Sucrose':0.90,'Hexane':0.90,
            'Cellulose':0.0,'Hemicellulose':0.0,'Lignin':0.0,'Ash':0.0,
            'OOO':0.0,'LLL':0.0,'LnLnLn':0.0,'SSS':0.0,
        }},
    )
    # @U204.add_specification(run=False)
    # def U204_make_wash():
    #    solids,extraction_water = U204.ins
    #    pellet=U204.ins[0].F_mass
    #    extraction_water= 0.01 * pellet       # 1:10 ratio (use pellet.F_mass - imass['Water'] for dry basis) 10 L/kg = 0.01 m3/kg  -> m3/hr

    #     # PURE water (no NaCl in the wash)
    #                  # pre-heated to 50 C
       
    #    # extraction_water= max(0,water_required)#TODO if we need to recylce
   
    # @U204.add_specification(run=True)
    # def U204_make_wash():
    #     pellet, water = U204.ins
    #     vol_L = 0.01 * pellet.F_mass        # 1:10 ratio (use pellet.F_mass - imass['Water'] for dry basis) 10 L/kg = 0.01 m3/kg  -> m3/hr

    #     U204.ins[1].imass['Water'] = 1.0 * vol_L  # PURE water (no NaCl in the wash)
     
                      # pre-heated to 50 C
    M111=bst.Mixer('M111',ins=(U204-0,U204-1),outs='mixed_slurry_centrifuge_2')
    U205 = bst.SolidsCentrifuge(
      'U205',
      ins=(M111-0),                                 # after the extraction the output is in th eform of slurry which then goes into centrifugation and separate into solids
      outs=('spent_solids', 'protein_supernatant'), # [0] solids-rich, [1] liquid-rich
      split={
          # insolubles -> spent solids
          'Cellulose': 1.0, 'Hemicellulose': 1.0, 'Lignin': 1.0, 'Ash': 1.0,
          # protein: 31% retained in spent solids --> 69% to supernatant refer[27]
          'Protein': 0.94,## 94.1% stays in final solids --> 5.9% to supernatant_2          
          'NaCl': 0, 'NaOH': 0, 'Glucose': 0.0, 'Sucrose': 0.0, 'Hexane': 0.0,
          'OOO': 0.0, 'LLL': 0.0, 'LnLnLn': 0.0, 'SSS': 0.0,
         
      },
      moisture_content=0.01,   #TODO check the moisture contetnt                    # wet-cake moisture 
      solids=('Cellulose', 'Hemicellulose', 'Lignin', 'Ash', 'Protein'),
      centrifuge_type='scroll_solid_bowl',
  )
    # pool the two supernatants
    M108 = bst.Mixer('M108', ins=(U203-1, U205-1), outs='pooled_supernatant')
    makeup_water_2 = bst.Stream('makeup_water_2', phase='l',  units='kg/hr', price=price['Makeup water'])
    M106=bst.Mixer('M106',ins=(makeup_water_2,fresh_NaOH), outs='10%_NaOH')
    @M106.add_specification(run=False)
    def M106_make_NaOH():
        water, NaOH = M106.ins
        slurry_vol_L = M108.outs[0].F_vol * 1000 #convert m3/h to L/h
        mass_NaOH_required=2.7e-5*slurry_vol_L    # refer calc 2.77731E-05 mass of naoH required to neutralize 1 L of solution
        total_solution_required = mass_NaOH_required/0.1  # 10 wt % of NaOH is solution so 0.1= mass of NaOH/mass of solution
        M106.ins[1].imass['NaOH']=mass_NaOH_required
        M106.ins[0].imass['Water']=0.90*total_solution_required
        M106._run()
   
    T201=bst.MixTank('T201', ins=(M108-0,M106-0), outs='neutralized_slurry')# to maintian pH the extraction nad centrifugation referred from refer [27] but this neutralization referred from 
    #Pilot plant scale-up of protein extraction by salt solubilization coupled with UF so we refer muconic membrane unit from ethanol_adipic because its uses ultrafilter as cost line
    U206 = units_oil.ProteinUltrafiltration('U206', ins=T201-0,
                                         outs=('UF_permeate', 'UF_retentate'))
    U207 = units_oil.Lyophilizer(
    'U207',
    ins=U206-1,                              # UF retentate (or U207-1 after diafiltration)
    outs=(protein_isolate, 'fd_vapor'),
    moisture_content=0.001,                   # very low  freeze-dried powder
    T=233.15,                                # -40 °C
)
 # =============================================================================
 # mix extracted and stripped oil
  # =============================================================================s
    M109= bst.Mixer('M109',ins=(U107-0,E102_P1-0),outs=crude_oil)
  
  
    # =============================================================================
    # Wastewater treatment & facilities 
    # =============================================================================

    M501 = bst.Mixer('M501', ins=(
        U206-0,          # UF permeate (saline protein wastewater)
        U207-1,          # lyophilizer vapor (mostly water)
    ), outs='process_wastewater')

    @M501.add_specification(run=False)
    def M501_spec():
        for i in M501.ins:
            i.phase = 'l'
        M501._run()

    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=M501-0,
        mockup=True,
        area=500,
    )

    sludge = wastewater_treatment_sys.get_outlet('sludge')
    biogas = wastewater_treatment_sys.get_outlet('biogas')
    RO_water = wastewater_treatment_sys.get_outlet('RO_treated_water')

    M510 = bst.Mixer('M510', ins=(
        U103-1,          # cleaning reject (Debris)
        U102-1,          # screen overs
        U205-0,          # final protein-extraction spent solids (fiber/meal)
        sludge,          # WWT sludge → boiler fuel (oilcane pattern)
    ), outs='wastes_to_boiler')

    @M510.add_specification(run=False)
    def M510_spec():
        for i in M510.ins:
            i.phase = 'l'   # same as succinic
        M510._run()

    create_facilities(
        solids_to_boiler=M510-0,
        gas_to_boiler=biogas,
        process_water_streams=[
            makeup_water_1,    # saline makeup (M105)
            makeup_water_2,    # NaOH neutralization makeup (M106)
            extraction_water,  # stage-2 wash water (H105)
        ],
        feedstock=raw_seed,    # scales CIP, plant air, fire water
        RO_water=RO_water,
        recycle_process_water='',   # no process-water recycle yet (succinic uses empty MX)
        blowdown_to_wastewater=M501-0,  # BT + CT blowdown back to WWT
        BT_area=700,
        area=900,
    )
    
    HXN = bst.HeatExchangerNetwork('HXN1001',
                                              ignored=[],
                                              cache_network=True,
                  )
    
   



oil_extraction_sys = create_covercress_oil_sys(ins=[feed], mockup=False)
oil_extraction_sys.simulate()

u = oil_extraction_sys.flowsheet.unit
s=oil_extraction_sys.flowsheet.stream
