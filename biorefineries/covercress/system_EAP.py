#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 19:06:37 2026

@author: princyk2
"""


# protein extraction with enzyme

#TODO after drying might need heat exchanger when cools to 70 ° C
#TODO find out the how water split during extraction
import biosteam as bst
from biosteam import preferences
import numpy as np
from biorefineries.cellulosic import create_facilities
from biosteam import SystemFactory
from biorefineries.covercress import units_oil 
from biorefineries.covercress.chemicals_cc import create_covercress_chemicals
from biorefineries.covercress.process_settings import price, CFs, chem_index, _GDP_2007_to_2010,_lps
from biosteam.units import LLEUnit
from thermosteam import indexer, equilibrium, settings
from biorefineries.corn.units import SlurryMixTank
__all__=('create_covercress_oil_sys')
bst.preferences.update(flow = 'kg/hr',composition = True)
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
total_feed = 13_310        # kg/hr total at gate assumed pennycress refer[4]

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
     'LinoleicAcid', 'LinolenicAcid','Stericacid','OOO', 'LLL', 'LnLnLn', 'SSS','NHP',
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
@SystemFactory(
    ID='oil_extraction_sys',
    ins=['raw_seed'],
    outs=[
        dict(ID='SAF', units='kg/hr', price=price['SAF']),
        dict(ID='protein_isolate', units='kg/hr', price=price['protein_isolate']),
        #  coproducts:
        dict(ID='naphtha_product', units='kg/hr', price=price['Naptha']),
        dict(ID='green_diesel', units='kg/hr', price=price['Green_diesel']),
        dict(ID='propane', units='kg/hr', price=price['Propane']),
    ],
)
# @SystemFactory(ID='oil_extraction_sys', ins=['raw_seed'], 
#                outs=[dict(ID='SAF', units='kg/hr', price=price['SAF']),
#                dict(ID='protein_powder', units='kg/hr', price=price['protein_isolate']),
#                ],
#                      )
        
 
# @SystemFactory(ID='oil_extraction_sys', ins=['raw_seed'], )
def create_covercress_oil_sys(ins,outs,moisture_content=0.03):# maxmium oil recovery obtained when the cooked seeds was at 3 %
    raw_seed, = ins
    # protein_isolate=outs
    SAF, protein_isolate, naphtha_product, green_diesel, propane = outs
    chemicals=feed.chemicals
    split=screw_press_permeate(chemicals)
    
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
    T00=units_oil.SeedStorage('T00', ins=U103-0,outs='seed_for_processing')
    U104 = units_oil.Cooking(
        'U104',
        ins=(T00-0, 'conditioner_air', 'conditioner_natural_gas'),
        outs=('conditioned_seed', 'conditioner_outlet_air', 'conditioner_emissions'),
        moisture_content=0.09,
        split=0.0,
        T=273+35,
        utility_agent='Steam',
    )#refer 4
    #refer [6] and refer [2] in the canola processing then operation wired from conditioninf->flaking>cooking, 
    # conditionng -->soften the outer shell and moisture levels for optimal preparation performance.in refer[2] the MC durinf storage is 9 % so I assumed once after seed cleaning we adjust the moisture to 9 % and then after cooking final MC reduced to 3-4 % refer[2]
    U105=units_oil.SeedFlaking("U105", ins=U104-0, outs= "seed_flakes", )
  
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
    U107 =bst.ScrewPress('U107', U106-0,outs=('pressed_oil','press_cake'),  #for cost refer[seider pg No482]
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
    
  
    fresh_solvent =bst.Stream('Hexane_1', phase='l',units='kg/hr',  price=price['Hexane'])
 
    T101 = units_oil.HexaneStorageTank('T101', ins=fresh_solvent, outs='hexane_from_storage')
  
    M102 = bst.Mixer('M102', ins=(T101-0,'solvent_recycle'), outs='U201_feed')

    @M102.add_specification(run=False)   # not run=True
    def adjust_solvent():
        fresh_out, solvent_recycle = M102.ins
        required_solvent = C102.outs[0].F_mass * 1.0
        if required_solvent <= 0:
            required_solvent = raw_seed.F_mass * hexane_mass_ratio
    
        recycled = solvent_recycle.imass['Hexane']
        if recycled > required_solvent:
            solvent_recycle.imass['Hexane'] = required_solvent
            recycled = required_solvent
    
        makeup = max(0.0, required_solvent - recycled)
        T101.ins[0].imass['Hexane'] = makeup   # Hexane_1
        T101._run()                            # then tank outlet
        M102._run()                            # mix once — NOT M102.run()
       


    H102 = bst.HXutility(
           'H102',
           ins=M102-0,
           outs='U201_hot_feed',
           T=60+273,
       )
    # =============================================================================
    #OIL extraction from cake 
    # =============================================================================
    U201 = bst.MultiStageMixerSettlers(
        'U201',
        N_stages=6,   #Refer [25]for Hexane  assumes overall efficiency is 50 %                   # real extractors ~6-8 ideal stages when stage =1, maximum Hexane in miscelle is 7 %
        ins=[C102-0, H102-0],
        # feed_stages=[0, -1],              # solids in top, fresh solvent bottom (countercurrent) 
        outs=['miscella', 'spent_cake'],
        # phases=('s', 'l'),
        partition_data={
            'IDs': ('OOO', 'LLL', 'LnLnLn', 'SSS','Hexane','Water','NHP'),
            'K': np.array([6.41E+05,6.41E+05,6.41E+05,6.41E+05,1.8,0.25,0]),  #as per calc hexane partila coefficent is 1.8 for maintian 52 % hexane in misceele  # ideal leaching: oil conc. equal in miscella & retained liquid
            # 'extract_chemicals': ['Hexane'],       # solvent reports to miscella
            'raffinate_chemicals': ['Cellulose','Hemicellulose','Lignin',
                                    'Ash','Protein','Glucose','Sucrose'],
        },
    )





    S100 = bst.Splitter('S100', ins=U201-0, outs=('hexane_loss_m','miscella_to_recovery'), split={'Hexane': 0.0})
    S101 = bst.Splitter('S101', ins=U201-1, outs=('hexane_loss_c','cake_to_DT'), split={'Hexane': 0.0})
    M100 = bst.Mixer('M100', ins=(S100-0, S101-0), outs='hexane_loss_total')
    alpha = 0.5
    Hex = chemicals.index('Hexane')
    @S101.add_specification(run=False)
    def set_hexane_loss():
        target_loss = 0.5 * feed.F_mass / 1000.0
        loss_m = alpha * target_loss
        loss_c = (1 - alpha) * target_loss
        hex_m = S100.ins[0].imass['Hexane']
        hex_c = S101.ins[0].imass['Hexane']
        fm = 0.0 if hex_m <= loss_m else min(1.0 - 1e-6, loss_m / hex_m)
        fc = 0.0 if hex_c <= loss_c else min(1.0 - 1e-6, loss_c / hex_c)
        S100.split[Hex] = fm
        S101.split[Hex] = fc
        S100._run()
        S101._run()
        M100._run()   
    F101 = units_oil.Desolventizer(
           'F101',
           ins=(S101-1, 'dryer_gas', 'Steam'),
           outs=('toasted_meal', 'desolventizer_vapor', 'emissions'),
           split={solvent_ID: 0.99},
           moisture_content=0.01,
           T=110.0 + 273.15,
           utility_agent='Steam',
       )


    H103 = units_oil.MealCooler(
        'H103',
        ins=F101-0,                   
        outs='meal_cooled',
        T=30.0 + 273.15,cool_only=True,
    ) 
    
    E101=bst.MultiEffectEvaporator('E101', ins=S100-1, outs=('E101_l','E101_g'), P=(101325, 73581, 50892, 32777, 20000), V=0.5, V_definition='First-effect',chemical='Hexane',flash=False)
    E101.target_conc=0.98 #refer[20]
    @E101.add_bounded_numerical_specification(x0=1e-3, x1=1.0 - 1e-3, xtol=1e-4, ytol=1e-4, x=0.5)
    def E101_V(V):
        if E101.ins[0].imol['Hexane'] <= 0:
            return 0.0
        E101.V=V
        E101._run()
        out_mass = E101.outs[0]
        mass_oil =E101.outs[0].imass['OOO']+E101.outs[0].imass['LLL']+ E101.outs[0].imass['LnLnLn']+E101.outs[0].imass['SSS']+E101.outs[0].imass['NHP']
        return mass_oil/out_mass.F_mass-E101.target_conc
    # E101.run_after_specifications = True
    E101_P0=bst.Pump('E101_P0',ins=E101-0,P=101325,)
    E101_P1=bst.IsothermalCompressor('E101_P1',ins=E101-1,P=101325,)
    steam = bst.Stream(
    'steam',
     Water=10,            
    phase='g',
    T=110+273,             # 110 °C, or up to 413.15 for 140 °C used low pressure steam
                 # match E102.P (150 mbar)
    units='kg/hr',
    price=_lps.regeneration_price,   # from process setting
 )
   
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
    1e-5,   # TODO update NHP K
])
    E102 = bst.Stripper(
    'E102',
     N_stages=7,# refer[23]        
    ins=[E101_P0-0,steam], # refer 18
    outs=['strip_overhead', 'stripped_oil'],#solute =water because already majority of hexane removed uisng MEE so other majority component is water that must removed so we use solute =water to calculate the internal calcs like stage efficiency etc also also its boiling boint is >hexane thus it would be partitoned in significant amount is both phases  
    P=6000,T=110+273, #110 as per literature and activity coefficent for 110 °C is 0.572 
    solute="Water",  partition_data={
          'IDs': ['Hexane', 'Water', 'OOO', 'LLL', 'LnLnLn', 'SSS','NHP'],'K':K_strip,
        'vapor_chemicals': ['Hexane', 'Water'],
          'liquid_chemicals': ['OOO', 'LLL', 'LnLnLn', 'SSS','NHP'],
      },
)    
    E102.x_hex_target=0.003 #refer[excel] 0.3 wt % residual hexane
 
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
    E102_P0=bst.IsothermalCompressor('E102_P0',ins=E102-0,P=101325)
    H1=bst.HXutility('H1',ins= E102-1, T=273+30)
    E102_P1=bst.Pump('E102_P1',ins=H1-0,P=101325)
    # =============================================================================
    #Solvent Recovery
    # =============================================================================
        

    M104=bst.Mixer('M104', ins=(F101-1,E101_P1-0,E102_P0-0),outs='hexane_vapour_with_water_and_noncondensates')
    
   #cooling of vapours to separate hexane 
    F104 = bst.Flash('F104', ins=M104-0, outs=('vent_gas', 'condensate'), P=101325,T=273+50) #TODOshould we recycle this vent gas ? and check for appropriate units because in industry they use minerla oil to abosorb the gases
    
    F104_P=bst.Pump('F104_P',ins=F104-1,outs='solvent_recycle')
    F104_P-0-1-M102

    
    
   
   
    # =============================================================================
    #Protein extraction feeed-dried cake
    # =============================================================================
    #protein extraction unit -->2-stage single-solvent extraction to manufacture the actual protein product refer[27]
    
    G201=units_oil.HammerMill('G201', ins=H103-0,outs='grinded_meal')
    fresh_nacl  = bst.Stream('fresh_nacl',  phase='s', units='kg/hr',price=price['NaCl'])   # 
    fresh_NaOH=bst.Stream('fresh_NaOH',units='kg/hr',price=price['NaOH'])
    makeup_water_1 = bst.Stream('makeup_water_1', phase='l',  units='kg/hr', price=price['Makeup water'])
    fresh_enzyme = bst.Stream('fresh_enzyme', units='kg/hr', price=price['enzyme'])

   #refer [42] for the enzyme asssited extraction
    # Heat water to enzyme optimum (use 50 °C; set 60 °C if needed)
   
    M105 = bst.Mixer('M105', ins=(makeup_water_1, ''), outs='extraction_water')
    @M105.add_specification(run=False)
    def M105_make_extraction_water():
        """
        Total water to extraction = 25 kg / kg meal protein (lab: 2 g protein / 50 mL).
        Recycle from clarified UF permeate (S206-0); makeup covers the deficit.
        """
        makeup, recycle = M105.ins
        meal = G201.outs[0]
        P = meal.imass['Protein']                 # kg/hr
        water_needed = 25.0 * P                   # kg/hr
        recycled_water = recycle.imass['Water']    # kg/hr (0 on first pass)
        M105.ins[0].imass['Water'] = max(0.0, water_needed - recycled_water)
        # leave recycle as-is (set by S206 / U208 upstream)
        M105._run()
    H104 = bst.HXutility('H104', ins=M105-0, outs='heated_water', T=50 + 273.15)
   
 # Slurry MixTank: meal + water 
    T201 = bst.MixTank(
        'T201',
        ins=(G201-0, H104-0),
        outs='meal_water_slurry',
        tau=0.5,                 # ~30 min dissolve (lab)
        
    )
    T201.line = 'Meal–water mixing'
    
        
    # ---------------------------------------------------------------------
    # 2) SlurryMixTank: pH adjust (NaOH)
    # ---------------------------------------------------------------------
    U202 = SlurryMixTank(
        'U202',
        ins=(T201-0, fresh_NaOH),
        outs='pH_adjusted_slurry',
        tau=0.25,
        
    )
    U202.line = 'pH adjustment'
    
    @U202.add_specification(run=False)
    def U202_pH_adjust():
        slurry, naoh = U202.ins
        P = G201.outs[0].imass['Protein']
        M105._run()
        H104._run()
        T201._run()
        naoh.imass['NaOH'] = 0.01 * P              #
        U202._run()
    
    # ---------------------------------------------------------------------
    # 3) SlurryMixTank: protease addition + incubation
    # ---------------------------------------------------------------------
  
    U202E = SlurryMixTank(
        'U202E',
        ins=(U202-0, fresh_enzyme),
        outs='enzyme_extraction_slurry',
        tau=3.0,                 # 3 h incubation
                   # keep 50–60 °C here (not 90)
    )
    U202E.line = 'Enzyme-assisted extraction'
    
    @U202E.add_specification(run=False)
    def U202E_enzyme_dose():
        slurry, enzyme = U202E.ins
        P = G201.outs[0].imass['Protein']
        enzyme.imass['Enzyme'] = 0.05 * P          # 5% on protein
        U202E._run()
   
    # Enzyme kill: heat slurry to 90 °C (10 min at lab;)
    H106 = bst.HXutility('H106', ins=U202E-0, outs='slurry', T=90 + 273.15)
    # Optional cool before UF (membranes cant hold 90 °C)
    H107 = bst.HXutility('H107', ins=H106-0, outs='cooled_supernatant', T=40 + 273.15)
    # Centrifuge: 90% protein -> supernatant  ->split to solids = 0.10
    U203 = bst.SolidsCentrifuge(
        'U203',
        ins=H107-0,
        outs=('spent_solids', 'protein_supernatant'),  # [0] solids, [1] liquid
        split={
            'Cellulose': 1.0, 'Hemicellulose': 1.0, 'Lignin': 1.0, 'Ash': 1.0,
            'Protein': 0.10,   # 10% left in solids -> 90% extracted yield calcuted by u.U207.outs[0]/u.G201.outs[0] is ~45 %
            'Enzyme': 0.0,    # denatured enzyme stays with liquid (to UF permeate)
            'NaOH': 0.0, 'HCl': 0.0, 'Glucose': 0.0, 'Sucrose': 0.0, 'Hexane': 0.0,
            'OOO': 0.0, 'LLL': 0.0, 'LnLnLn': 0.0, 'SSS': 0.0, 'NHP': 0.0,
        },
        moisture_content=None,   # wet spent solids; adjust if needed
        solids=('Cellulose', 'Hemicellulose', 'Lignin', 'Ash', 'Protein'),
        centrifuge_type='scroll_solid_bowl',
        strict_moisture_content=False,
    )
   
   
   
    # Same downstream as before
    U206 = units_oil.ProteinUltrafiltration(
        'U206',
        ins=U203-1,
        outs=('UF_permeate', 'UF_retentate'),
    )
    
    U208 = bst.SolidsCentrifuge(
    'U208',
    ins=U206-0,
    outs=('permeate_solids', 'clarified_water'),
    split={
        # to solids [0] → boiler
        'Enzyme': 1.0,
        'Protein': 1.0,
        'Cellulose': 1.0, 'Hemicellulose': 1.0, 'Lignin': 1.0, 'Ash': 1.0,
        # stay in water [1] → recycle
        'NaOH': 0.0, 'Glucose': 0.0, 'Sucrose': 0.0, 'Hexane': 0.0, 'Water': 0.0,
    },
    moisture_content=None,
    solids=('Enzyme', 'Protein', 'Cellulose', 'Hemicellulose', 'Lignin', 'Ash'),
    centrifuge_type='scroll_solid_bowl',
    strict_moisture_content=False,
)
    S206 = bst.Splitter(
    'S206',
    ins=U208-1,
    outs=('protein_water_recycle', 'protein_water_purge'),
    split=0.90,   # 90% recycle / 10% purge (NaOH + sugar bleed)
)
    # H108=bst.HXutility('H108', ins=S206-0, outs='cooled_swater', T=50 + 273.15)
    S206-0-1-M105
    # U207 = units_oil.Lyophilizer(
    #     'U207',
    #     ins=U206-1,
    #     outs=(protein_isolate, 'fd_vapor'),
    #     moisture_content=0.001,
    #     T=233.15,
    # )
#     U207 = bst.SolidsCentrifuge(
#     'U207',
#     ins=U206-1,                                    # UF retentate
#     outs=(protein_isolate, 'protein_centrate'),    # [0] product, [1] water
#     split={
#         'Protein': 0.995,      # protein isolate (solids)
#         'Enzyme': 1.0,
#         'Cellulose': 1.0, 'Hemicellulose': 1.0, 'Lignin': 1.0, 'Ash': 1.0,
#         'NaOH': 0.0, 'Glucose': 0.0, 'Sucrose': 0.0, 'Hexane': 0.0, 'Water': 0.0, #glucose, sucose are water soluble
#     },
#     moisture_content=0.10,    # wet protein cake; or None if it crashes
#     solids=('Protein', 'Enzyme', 'Cellulose', 'Hemicellulose', 'Lignin', 'Ash'),
#     centrifuge_type='scroll_solid_bowl',
#     strict_moisture_content=False,
# )
#     U207.line = 'Protein dewatering centrifuge'
    U207 = units_oil.Drying(
    'U207',
    ins=(U206-1, 'dryer_air', 'dryer_natural_gas'),
    outs=(protein_isolate, 'dryer_outlet_air', 'dryer_emissions'),
    moisture_content=0.05,
    split=0.0,
    T=80 + 273.15,
    utility_agent='Steam',
)
 # =============================================================================
 # mix extracted and stripped oil
  # =============================================================================s
    M109= bst.Mixer('M109',ins=(U107-0,E102_P1-0),outs='mixed_crude_oil')
    
    #Degumming parameters 
    # CITRIC_G_PER_100G_OIL = 0.63      # 0.05 – 1.2 refer[29] avg-0.63
    # WATER_WT_ON_OIL = 0.02           # 0.01 – 0.03
    # CITRIC_SOLUTION_WT = 0.30        # 30 wt% citric in solution refer[4]
                   
    # DEGUM_TAU_HR = 20. / 60. #refer [4]
    # NHP_WT_ON_OIL = 0.01             # impurity in crude oil oleochemicals
    citric_acid_feed = bst.Stream('citric_acid_feed', units='kg/hr',price=price['citric_acid'])
    degum_water_feed = bst.Stream('degum_water_feed', units='kg/hr', price=price['Makeup water'])
    T100=bst.StorageTank('T100',ins=citric_acid_feed)
    T100.line = 'Citric acid solution storage'
    H201=bst.HXutility('H201', M109-0, T=273.15 + 80)  # K 80 °C refer [4]
    M201 = bst.Mixer('M201', ins=(T100-0, degum_water_feed), outs='degumming_liquor')

    @M201.add_specification(run=False)
    def dose_degumming_chemicals():
        oil = H201.outs[0]
        nhp_mol = oil.imol['NHP']
        oil_mass = sum(oil.imass[i] for i in ('OOO', 'LLL', 'LnLnLn', 'SSS'))
        mw_citric = oil.chemicals['Citric_acid'].MW

        # Eq. (1): 1 citric per 6 NHP — stoichiometric dose 
        citric_mol = nhp_mol / 6
        citric_mass = citric_mol * mw_citric

        # 30 wt% citric solution
        M201.ins[0].imass['Citric_acid'] = citric_mass # we calc citric acid mass per stiochiometry
       
        M201.ins[0].imass['Water'] = citric_mass * (1 - 0.30) / 0.30
        

        # Hydration wash water (2 wt% on oil)
        M201.ins[1].imass['Water'] = 0.02 * oil_mass #TODO update the water composition
        M201._run()

    #  Mix feeds 
    M202 = bst.Mixer('M202', ins=(H201-0, M201-0), outs='degumming_feed')

    #  Reactor 
    degumming_rxn = bst.Reaction(
    'Citric_acid + 6 NHP + 0.2 Water -> 6.2 PL',  # refre [4] but water moles #TODO update 0.2 moles because non hydrous PL to PL 6>6.2
    reactant='NHP',
    X=1.0,
    check_atomic_balance=False,
)
    
    R201 = bst.SinglePhaseReactor(
        'R201',
        ins=M202-0,
        outs='degumming_slurry',
        reaction=degumming_rxn,
        T=273.15 + 80,
        batch=False,
        tau=20 / 60,   # 20 min refer[4]
    )

    #  Centrifuge
    C201 = bst.LiquidsSplitCentrifuge(
        'C201',
        ins=R201-0,
        outs=('crude_oil_purified', 'degumming_gums'),
        split={
        # outs[0] crude_oil (light / purified oil)
        'OOO': 0.99, 'LLL': 0.99, 'LnLnLn': 0.99, 'SSS': 0.99,
        'Hexane': 0.99,
        # → outs[1] degumming_gums (heavy)
        'NHP': 0.001,   # reaction product PL;
        'PL': 0.00,
       
        'Citric_acid': 0.05,
        'Water': 0.02,
    },
    )
   
    # =============================================================================
    # HEFA
    # =============================================================================
    # =============================================================================
    #catlyst recovery 
    # =============================================================================
    cat_fresh = bst.Stream('cat_fresh', units='kg/hr',price=price['Pt/Al2O3'],phase='s')
    cat_recycle = bst.Stream('cat_recycle', phase='s', units='kg/hr',)
    catalyst_lifetime_hr = 8400  # refer [34] catalyst life is 8400h
    
    M300=bst.Mixer('M300', ins=(cat_fresh,cat_recycle), outs='catalyst_HDO & HI', )
    @M300.add_specification(run=False)
    def adjust_catalyst_makeup():
        feed_R301 = HX301-1
        feed_R302 = HX304-1

        m_R301 = R301._required_catalyst(feed_R301)
        m_R302 = R302._required_catalyst(feed_R302)
        m_total = m_R301 + m_R302
        holdup_hr = R301.tau + R302.tau 

        M_inventory = m_total * holdup_hr 
        makeup = M_inventory / catalyst_lifetime_hr

        recycled = S301.outs[0].imass['Pt/Al2O3']
        cat_fresh.imass['Pt/Al2O3'] = max(makeup,m_total-recycled)   
        cat_fresh.phase = 's'
        cat_recycle.phase = 's'

        M300._run()

    
    hydrogen_fresh = bst.Stream(
        'hydrogen_fresh',
        Hydrogen=1.0,          # initial guess [kmol/hr]; 
        phase='g',
        units='kmol/hr',
        P=3e6,                   # 30 bar — matches R301 & R401
        price=price['Hydrogen'],
    )
    T102 = units_oil.HydrogenGasStorageTank(
        'T102',
        ins=hydrogen_fresh,
        outs=('H2'),
        tau=7 * 24,              # 7 d inventory — adjust from ref [4]#TODO
        V_wf=0.9,
        vessel_material='Carbon steel',
    )
    
    T102.line = 'H2 storage tank'
  
    M301=bst.Mixer('M301',ins=(T102-0,'',''),outs='H2_supply',)
    @M301.add_specification(run=False)
    def _adjust_H2_flow():
        fresh_h2, h2_rec_hdo, h2_rec_hck = M301.ins
        n_target_R301 = (1 + R301.excess_H2) * R301._required_H2_mol(R301.ins[0])
        n_target_R302 = (1 + R302.excess_H2) * R302._required_H2_mol(R302.ins[0]) 
        n_target=n_target_R301+n_target_R302
        # n_recycle = h2_rec_hdo.imol['Hydrogen'] + h2_rec_hck.imol['Hydrogen']
        n_recycle = S304.outs[1].imol['Hydrogen'] + F306.outs[1].imol['Hydrogen']
        hydrogen_fresh.imol['Hydrogen'] = max(0., n_target - n_recycle)
        hydrogen_fresh.phase = 'g'
        M301._run()
    
    S300=bst.Splitter('S300',ins=M301-0, outs=('H2_HDO','H2_HI'), split=0.5)
    @S300.add_specification(run=False)
    def _adjust_H2_split():
        n_R301 = (1 + R301.excess_H2) * R301._required_H2_mol(R301.ins[0])
        n_R302 = (1 + R302.excess_H2) * R302._required_H2_mol(R302.ins[0])
        n_tot = n_R301 + n_R302
        S300.split = n_R301 / n_tot if n_tot > 0 else 0.5
        S300._run()
        for outs in S300.outs:
            outs.phase = 'g'
        
    K301=bst.IsothermalCompressor('K301',ins=S300-0,P=4e6)
    K302=bst.IsothermalCompressor('K302',ins=S300-1,P=3e6)
    
  
    
    HX_stream_1=bst.Stream('HX_stream_1') #referred Glycerolsis setup with HX from oilcane
    P102=bst.Pump('P102',ins=C201-0, P=4e6)
    HX301=bst.HXprocess('HX301',ins=(HX_stream_1,(P102)-0,),outs=('cooled_R301_effluent', 'R301_preheated_feed'), phase0='l', phase1='l') ##why HX as in the refer [4] process heat recovery is opted to preheat the feed to HDO reactor by utilizing the heat from the outelt
    #phase 0 and 1 =l because we assume no rigours VLE refer tutorial , this is just heat transfer between liquids no flash integration
    # =============================================================================
    # HDO 
    # =============================================================================
    S300_C = bst.Splitter(
    'S300_C',
    ins=M300-0,
    outs=('cat_to_R301', 'cat_to_R302'),
    split=0.5,   # 
)
    @S300_C.add_specification(run=False)
    def _adjust_catalsyt_split():
        feed_R301 = HX301-1
        feed_R302 = HX304-1
        m_R301 = R301._required_catalyst(feed_R301)
        m_R302 = R302._required_catalyst(feed_R302)
        m_total = m_R301 + m_R302
        S300_C.split = m_R301 / m_total if m_total > 0 else 0.5
        S300_C._run()
        for outs in S300_C.outs:
            outs.phase = 's'
           
    R301 = units_oil.Hydrodeoxygenation(
        'R301',
        ins=(HX301-1, K301-0, S300_C-0),
        outs='light_hydrocarbons and vapours',
        T=583.15, P=4e6, tau=2.0,   # refer[4]
        catalyst_wt_frac=0.05,
        excess_H2=0.3,
      #TODO check for exact range of excess H2 check for other lit and set uncertainity range
)
    V301=bst.IsenthalpicValve('V301',ins=R301-0, outs=HX_stream_1, P=1e6, vle=True)
    
    # =============================================================================
    # GasSeparator
    # =============================================================================
    F301 = bst.SplitFlash(
    'F301',
    ins=HX301-0,
    outs=('hdo_gases', 'hdo_praffin'),
    split={
        # mostly to vapor (lights)
        'Hydrogen': 0.995,
        'Propane': 0.95,
        'CO': 0.995,
        'CO2': 0.995, #Recovery in 3-phase sepa'rator 99% refer [30]
        # at 25 bar, water stays liquid at high pressure : The boiling point increases to 223.8°C
        'Water': 0.995,
        # mostly to liquid (paraffins + cat)
        'Octadecane': 0.0,
        'Heptadecane': 0.0,
        'Stericacid': 0.0,
        'OOO': 0.0, 'LLL': 0.0, 'LnLnLn': 0.0, 'SSS': 0.0,
        'Pt/Al2O3': 0.0,
        'NHP': 0.01, 'PL': 0.0,
        'Citric_acid': 0.0,
        'Hexane': 0.995, #why split=1 because in R302 for hexane traces R302 failed to compute molar volume of hexane as it in crtical temp no longer can be in liquid
    },
    P=1e6,              # 
    T=273.15 + 100,       # TODO check later at what temp and pressure flash removes gases from paraffins
    Q=0,
)
    F301.line = 'HDO 3-phase separator'
     
    # =============================================================================
    # Dewatering refer [34]
    # =============================================================================
    F302 = bst.Flash(  #to maintian VLE
    'F302',
    ins=F301-0,
    outs=('gas_to_recovery', 'hdo_wastewater'),
    T=273.15 + 28,  # TODO check later at what temp and pressure flash removes gases from paraffins
    P=1e6,
)
     
    # =============================================================================
    # Propane recovery-1 refer[34] #TODOcheck later do we have another option in biosteam to separate gases 
    # =============================================================================
    # V302_g=bst.IsenthalpicValve('V302',ins=F302-0, P=101325, vle=False)
    # V302_l=bst.IsenthalpicValve('V302',ins=F302-1, P=101325, vle=False)
    # S303 =units_oil.PropaneDistillation(
    # 'S303',
    # ins=V302_g-0,
    # outs=('gas_to_recovery', 'propane_'),split={
    #     'Propane': 0.01,     # 90% C3 product (outs[1]) — 
    #     'CO2': 0.9995,         
    #     'CO': 0.995,
    #     'Hydrogen': 0.995,    
    #     'Water': 0.995,
    #     'Hexane': 1, 
    # },
    
# )

#refer gas fermnetaion isentropic compressor is used as gases increase their temp while compreses so cant assume the isothermal isentroipic is adiabatic and then we use HX to cool this matches the reality 
    K1=bst.IsentropicCompressor('K1', ins=F302-0, P=2.5e6, vle=True) #why K1 To condense LPG you need refrigeration (you use T=273-40) and often high pressure (you use K1 to 25 bar before flash).
    T1 = bst.HXutility('K100', ins=K1-0, outs = ['s2'],T=273-40,)
    S303 = bst.Flash(
    'S303',
    ins=T1-0,
    outs=('fuel_gas_overhead', 'lpg_condensate'),
    P=2.5e6,          # match compressor outlet
  
  Q=0
            # adiabatic flash — VLE at fixed T, P
)
#     #why not shortcut because H₂, CO, CO₂ do not form a liquid phase at any practical column pressure/temperature. They are non-condensables in the model.
# Propane is only significantly in the liquid phase if you compress and refrigerate (your K1 → T1 at −40 °C, 25 bar).
# ShortcutColumn still tries to build stage-by-stage VLE with relative volatility. If a species is only vapor (H₂) or wrongly treated as liquid (your comment on D301 crash), thermodynamics break → infeasible region / enthalpy errors.
# the feed is non-condensable gas + one very light condensable. BioSTEAM’s ShortcutColumn is built for liquid hydrocarbon distillation; 
    K1_l=bst.IsenthalpicValve('K1_l',ins=S303-1,P=1e5)
    K1_g=bst.IsenthalpicValve('K1_g',ins=S303-0,P=1e5)
    S304=bst.Splitter('S304', ins=K1_g-0, outs=('acid_gas','H2_recycle'), split={'Water':0.995, 'CO':0.999,'CO2':0.999,'Hydrogen':0.10})
   
    
    S304-1-1-M301
   
    
   
    # V303=bst.IsenthalpicValve('V303',ins=F301-1,  P=101325, vle=False)
    
    
   
    


#     F303 = bst.SplitFlash(
#     'F303',
#     ins=V303-0,
#     outs=('hdo_water_vapor', 'paraffins_to_isom'),
#     split={
#         'Water': 0.999,
#         'Hydrogen': 0.99,
#         'Propane': 0.99,
#         'CO': 0.99,
#         'CO2': 0.99,
#         'Hexane': 1,
#         'Octadecane': 0.0,
#         'Heptadecane': 0.0,
#         'Stericacid': 0.0,
#         'Pt/Al2O3': 0.00,
#         'OOO': 0.0, 'LLL': 0.0, 'LnLnLn': 0.0, 'SSS': 0.0,
#         'NHP': 0.0, 'PL': 0.0, 'Citric_acid': 0.0,
#     },
#     P=101325,              # ~0.1 bar — biodiesel F401; refer[biodiesel -->systems
#     T=273.15 + 120,        # ~120 °C vacuum dryer — trefer [4]/[30]
#     Q=0,
# )
    # F303.line = 'Paraffin dewatering'
    # F303_1 = bst.SplitFlash(
    # 'F303_1',
    # ins=F303-1,
    # outs=('lights_to_boiler', 'liquid_paraffins_to_pump'),
    # split={
    #     'CO2': 0.999, 'CO': 0.999, 'Hydrogen': 0.01,
    #     'Propane': 0.0, 'Water': 0.0,
    #     # heavy paraffins →0 (stay liquid)
        
    # },
    # P=101325., T=273.15 + 120,
# )
#     P301=bst.IsothermalCompressor('P301',ins=F303-0, P=101325) #since the input is gas so we use isothermal compressor
    # HX303 = bst.HXutility('HX303', ins=F303-1, T=273.15+25, rigorous=False)
    
    S305 = bst.Flash( #why splitter because for R302 the liquid molar enthalpy model failed to  find molar enthalpy for hexane at 523
    'S305',
    ins=F301-1,
    outs=('hexane_lights_to_boiler', 'paraffins_to_pump'),
    # split={
    #     'Hexane': 1,      # # why 1 because for R302 the liquid molar enthalpy model failed to  find molar enthalpy for hexane at 523
    #     'Water':0.995,
    #     'CO2': 0.999,
    #     'CO': 0.999,
    #     'Hydrogen': 0.99,
    #     'Propane': 0.99,
    #     # paraffins default 0 stay on outs[1]
    # },
    P=1e6,              # ~0.1 bar — biodiesel F401; refer[biodiesel -->systems
       # T=273.15 + 120,        # ~120 °C vacuum dryer — trefer [4]/[30]
       Q=0,
   
)
    cat_hck = bst.Stream('cat_hck', units='kg/hr', price=price['Pt/Al2O3'])
    # H2=bst.HXutility('H2',ins=S305-1, T=273+40)
  
    P201 = bst.Pump('P201', ins=S305-1, P=3e6) #TODO or valve check later
    
#     M302 = bst.Mixer('M302', ins=(F303_1-1, F302-1), outs='paraffin_feed_hydrocracking')
    
    # HX304 = bst.HXutility('HX304', ins=P201-0, T=273.15 + 250, rigorous=False)
    HX_stream_2=bst.Stream('HX_stream_2') #referred Glycerolsis setup with HX from oilcane
    HX304=bst.HXprocess('HX304',ins=(HX_stream_2,P201-0,),outs=('cooled_R302_effluent', 'R302_preheated_feed'), phase0='l', phase1='l') 
#     # =============================================================================
#     # Hydrocracking/Isomerization
#     # =============================================================================
    R302 = units_oil.Hydrocracking(
    'R302',
    ins=(HX304-1, K302-0, S300_C-1),     # h2 from T602
    outs='hck_effluent',
    T=273.15 + 250, # as refernce [4] suggests 355 C but decane and heptane forms gas at crtical temp which is for hexane-507 K and for decane 618 so for safe choice 250-260 C
    P=3e6, #refer [4]
    tau=1.5,#TODO update 
    catalyst_wt_frac=0.02,
    excess_H2=0.3,
)
    R302.line = 'hydrocracking reactor'
    # @R302.add_specification(run=False)
    # def R302_run():
    #     R302._run()
    #     s = R302.outs[0]
    #     n = s.imol['Hexane']
    #     if n:
    #         s.phases = ('l', 'g')
    #         s.imol['l', 'Hexane'] = 0
    #         s.imol['g', 'Hexane'] = n
    
    V304=bst.IsenthalpicValve('V304',ins=R302-0,outs=HX_stream_2, P=101325)
    V304.line='Depressurizer_R302_out'
    # HX305 = bst.HXutility('HX305', ins=HX304-0, T=273.15 + 98, rigorous=False)
    
    F304 = bst.SplitFlash(
        'F304',
        ins=HX304-0,
        outs=('hck_flash_vapor', 'hck_flash_liquid'),
        split={
            # vapor @ 98 °C, ~1 atm non condensables 
            'Hydrogen': 1,#0.999,why split=1 because when these gases entered into D301 column crash occurs due to presence of these gases in may be liquid state as we unlock the phases of these gases
            'CO2': 1,#0.999,
            'CO': 1,#0.999,
            'Water': 0.95,
            'Propane': 0.95,
            'Heptane': 0.15,     # small flash strip only (Tb ~ 98 °C)#Species	Normal Tb (°C)but according to lit Remove non-condensables (H₂, CO₂, CO), water, propane (and similar lights)
            # Octane # ~126 # Mostly liquid, some vapor
# Propane Tb-42# Essentially all vapor
# Heptane # ~98 # Near boiling → large vapor fraction 
# Nonane # ~151 # Mostly liquid
# Decane# ~174# Almost all liquid
            'Octane': 0.00,
            'Nonane': 0.0,
            'Decane': 0.0,
            # liquid / solids
            'Pt/Al2O3': 0.0,
            'Undecane': 0.00, 'Dodecane': 0.00, 'Tridecane': 0.00,
            'Tetradecane': 0.0, 'Pentadecane': 0.00,
            'Hexadecane': 0.0, 'Heptadecane': 0.0, 'Octadecane': 0.0,
            'IsoOctane': 0.0, 'IsoDecane': 0.00, 'IsoHeptadecane': 0.0,
           
        },
        P=101325,
        T=273.15 + 98,
        Q=0,
    )
    F304.line = 'Post-HydroCrack_gas-liquid separator'
    # HX305 = bst.HXutility(
    # 'HX305',
    # ins=F304-0,
    # T=273.15 + 40,
    # V=0,
    # rigorous=False,
# )
#     # =============================================================================
#     # Propane recovery-2
#     # =============================================================================
    # F305 = units_oil.PropaneDistillation(
    #     'F305',
    #     ins=F304-0,          # or F304-0 if no HX401
    #     outs=('propane', 'h2_co2_co_gas'),
    #     split={
    #         # outs[0] lpg_product (mostly C3)
    #         'Propane': 0.999,      # 92% propane to product
    #         'Hydrogen': 0.0,     # little H2 in LPG
    #         'CO2': 0.0,
    #         'CO': 0.0,
    #         'Water': 0.00,
    #         # lights stay with H2 stream (outs[1])
    #         'Heptane': 0.00,
    #         'Octane': 0.0,
    #     },
    # )
    # F305.line = 'Propane recovery (LPG product)'
    K2=bst.IsentropicCompressor('K2', ins=F304-0, P=2.5e6, vle=True)
    T2 = bst.HXutility('K100', ins=K2-0, outs = ['s2'],T=273-40,)
    F305 = bst.Flash(
    'F305',
    ins=T2-0,
    outs=('fuel_gas_overhead', 'lpg_condensate'),
    P=2.5e6,          # match compressor outlet
  
 T=273.15 - 40,#make it isothermal if its adiabatic the propane separatin was not adequate in adiabatic condition temp of the oulet increases makes the propane in vapour state.
            # Why flash , because flash establishes single stage equilibrium for even gases unlike shortcut column ie, the compenents should be in the both phases during VLE for flash its not essentila
)
    K3_l=bst.IsenthalpicValve('K3_l',ins=F305-1,P=1e5)
    K3_g=bst.IsenthalpicValve('K3_g',ins=F305-0,P=1e5)
#     # =============================================================================
#     # Hydrogen recovery
#     # =============================================================================
    F306 = bst.Splitter(
    'F306',
    ins=K3_g-0,
    outs=('acid_gas_purge', 'h2_recycle_raw'),
    split={
        # outs[0] acid gas / flare
        'CO2': 0.999,          # 98% CO2 rejected
        'CO': 0.999,
        'Water': 0.90,
        'Hydrogen': 0.03,     # ~3% H2 purge (avoids CO buildup in loop)
        'Propane': 0.995,      # remaining propane split 
        'Heptane': 0.99, },
        # outs[1] h2_recycle_raw: mostly H2
        # P=101325,
        # T=273.15 + 98,
        # Q=0,
   
)
    F306.line = 'CO2/CO removal & H2 recovery'
    F306-1-2-M301
    product_storage_T=273+30
    _PRODUCT_TAU_HR = 7 * 24 
    M303=bst.Mixer('M303', ins=(K3_l-0,K1_l-0), outs='propane_mix')
    HX302 = bst.HXutility('HX302', ins=M303-0,outs='propane_cooled', T=273.15 -40, rigorous=False) # biosteam uses the ethylen as cooling agent so CWP cost  is 0
    T600 = bst.StorageTank(
        'T600', ins=HX302-0, outs=propane,
        tau=_PRODUCT_TAU_HR, V_wf=0.9,
    )
    S301 = bst.SolidsCentrifuge(
    'S301',
    ins=F304-1,
    outs=('spent_catalyst', 'decat_liquid'),
    split={'Pt/Al2O3': 1},moisture_content=None,#assume all calatayst is separated if any traces is present in the downstream biosteam model crashes because the model fails to find gas enthalpy of catalsyt
) 
    # HX306=bst.HXutility('HX306', ins=S301-0, outs='cooled_catalsyt', T=273+25)
    # HX306.outs[0].phase = 's'
    S301-0-1-M300
    D301 = bst.ShortcutColumn(
    'D301',
    ins=S301-1,  # or S401-0 after catalyst removal
    outs=('light_ends', 'd401_bottoms'),
    LHK=('Propane', 'Water'), #refer from fig 1 [4] most of propane is distillaed out in first colum
    Lr=0.99,   # most heptane to bottoms (into naphtha/SAF pool), not light_ends even with Lr and Hr =0.99 no change in results
    
    Hr=0.95,   # most octane to bottoms too — split C8 later in D302
    P=101325,
    k=1.5,
)
    D301.check_LHK = False #to ignore the validation check by biosteam whether there any interediate Tb between light nad heavy
    D302 = bst.ShortcutColumn(
    'D302',
    ins=D301-1,                    # d401_bottoms
    outs=('D302_distillate', 'd302_bottoms'),
    LHK=('Heptane', 'Octane'), #refer table 2  [4] naptha C5-C7 LK = last carbon number in that product C7 (Heptane).HK = first carbon number of the next product  C8 (Octane) for SAF.
    Lr=0.99,    # 51.62% IsoNonane distillate --> this is from refer[converion pathway], TODO 
    Hr=0.9456,    # 94.56% Decane bottoms (5.44% in distillate) 
    P=101325,
    k=1.5,
)
    D302.check_LHK = False
#     D302 = bst.ShortcutColumn(
#     'D302',
#     ins=D301-1,                    # d401_bottoms
#     outs=(naphtha_product, 'd302_bottoms'),
#     LHK=('IsoNonane', 'Decane'),
#     Lr=0.5162,    # 51.62% IsoNonane distillate
#     Hr=0.9456,    # 94.56% Decane bottoms (5.44% in distillate) 
#     P=101325,
#     k=1.5,
# )
#     D302.check_LHK = False
    D303 = bst.ShortcutColumn(
    'D303',
    ins=D302-1,                    # d402_bottoms
    outs=('D303_distilalte', 'D303_bottoms'),
    LHK=('IsoPentadecane','Heptadecane'),  # C16 / C17 cut when i check the entire components in the light and heavy keys in LHK L have to be higher component that needs to be in the overhead
    Lr=0.99,                         # C16 and lighter  distillate (SAF) C8-C16
    Hr=0.92,                         # C17 bottoms (diesel)C17-C22 
    P=101325,                          # mild vacuum (optional)
    k=1.5,
)
    D303.check_LHK = False
    # product_storage_T=273+30
    # _PRODUCT_TAU_HR = 7 * 24           
    HX305=bst.HXutility(
        'HX305', ins=D302-0, outs='naphtha_product_cooled',
        T=product_storage_T, V=0, rigorous=False,
    )
    HX306 = bst.HXutility(
        'HX306', ins=D303-0, outs='SAF_cooled',
        T=product_storage_T, V=0, rigorous=False,
    )
    HX307 = bst.HXutility(
        'HX307', ins=D303-1, outs='green_diesel_cooled',
        T=product_storage_T, V=0, rigorous=False,
    )
    T601 = bst.StorageTank(
        'T601', ins=HX305-0, outs=naphtha_product,
        tau=_PRODUCT_TAU_HR, V_wf=0.9,
    )
    T601.line = 'Naphtha product storage'
    T602 = bst.StorageTank(
        'T602', ins=HX306-0, outs=SAF,
        tau=_PRODUCT_TAU_HR, V_wf=0.9,
    )
    T602.line = 'SAF product storage'
    T603 = bst.StorageTank(
        'T603', ins=HX307-0, outs=green_diesel,
        tau=_PRODUCT_TAU_HR, V_wf=0.9,
    )
    T603.line = 'Green diesel product storage'

    # =============================================================================
    # Wastewater treatment & facilities 
    # =============================================================================

    M501 = bst.Mixer('M501', ins=('',S206-1,
        # U206-0,          # UF permeate (saline protein wastewater)
        U207-1,          # lyophilizer vapor (mostly water)
        C201-1,          #degumming _gums from centrifuge
     
        F104-0,
         F302-1, # waster water from the splitter which splits water from the vaoprs
       
    ), outs='process_wastewater')

    @M501.add_specification(run=False)
    def M501_spec():
        for i in M501.ins:
            i.phase = 'l'
        M501._run()

    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=M501-0, # same like SA blowdown mixer is wired to M501  M501-0 waster water from units
        mockup=True,
        area=500,
    )
    for ui in wastewater_treatment_sys.units:
        if type(ui).__name__ == 'SludgeCentrifuge': # in the module conventional in WWS_sys there is default sludgecentrifuge moisture ocntent i s0.79 but for this sys BIOSTEAM crash dude to infeabile moisture content (could be becuase of mor eoil and sucrose) so we set the MC to 0.55
           ui.strict_moisture_content = False
           ui.moisture_content = 0.55   # 
    sludge = wastewater_treatment_sys.get_outlet('sludge')
    biogas = wastewater_treatment_sys.get_outlet('biogas')
    RO_water = wastewater_treatment_sys.get_outlet('RO_treated_water')
    MX = bst.Mixer('MX401', ['', ''])
    M510 = bst.Mixer('M510', ins=(
        U103-1,          # cleaning reject (Debris)
        U102-1,          # screen overs
        U203-0,          # final protein-extraction spent solids (fiber/meal)
        U208-0,
        sludge,          # WWT sludge → boiler fuel (oilcane pattern)
       # hexane from S305
        
    ), outs='wastes_to_boiler')

    @M510.add_specification(run=False)
    def M510_spec():
        for i in M510.ins:
            i.phase = 'l'   # same as succinic
        M510._run()
      
    # =============================================================================
    # i acid gas mainly H2, propane are combustible which  combusted in flare before venting making  refer [33]
    # =============================================================================
    # M511 = bst.Mixer(
    #     'M511',
    #     ins=(biogas),#F306-0,S304-0,S305-0,),
    #     outs='fuel_gas_to_boiler',
    # )
    # @M511.add_specification(run=False)
    # def M511_spec():
    #     M511.ins[0].phase = 'g'
    #     M511.ins[1].phase = 'g'
    #     M511._run()
    #     M511.outs[0].phase = 'g'
    create_facilities(
        solids_to_boiler=M510-0,
        gas_to_boiler=biogas, #M511-0,
        process_water_streams=[
            makeup_water_1,    # saline makeup (M105)
            # makeup_water_2,    # NaOH neutralization makeup (M106)
            # extraction_water,  # stage-2 wash water (H105)
            degum_water_feed,
        ],
        feedstock=raw_seed,    # scales CIP, plant air, fire water
        RO_water=RO_water,
        recycle_process_water=MX-0,   # no process-water recycle yet (succinic uses empty MX) #TODO update
       blowdown_to_wastewater=M501.ins[0],  # Waste water from BT + CT blowdown back to WWT
        BT_area=700,
        area=900,
    )
    
    HXN = bst.HeatExchangerNetwork('HXN1001',
                                            ignored=[HX301,HX304],  # ignored=[HX301,],#HX304],
                                              cache_network=True,
                  )
    
   



oil_extraction_sys = create_covercress_oil_sys(ins=[feed], mockup=False)
oil_extraction_sys.simulate()

u = oil_extraction_sys.flowsheet.unit
s=oil_extraction_sys.flowsheet.stream

# =============================================================================
# TEA — SAF MPSP & feedstock MFPP
# =============================================================================
from biorefineries.tea.cellulosic_ethanol_tea import create_cellulosic_ethanol_tea

num_sims = 5
num_solve_tea = 3

feedstock = feed
SAF = s.SAF
protein_isolate = s.protein_isolate
naphtha_product = s.naphtha_product
green_diesel = s.green_diesel
propane = s.propane

BT = next((ui for ui in oil_extraction_sys.units
           if isinstance(ui, bst.BoilerTurbogenerator)), None)

get_flow_dry_tpd = lambda: feedstock.F_mass * 24 / 907.185

def set_market_prices():
    """Reset only feed + products (other feeds keep price= from Stream creation)."""
    feedstock.price = price['Feedstock']
    SAF.price = price['SAF']
    protein_isolate.price = price['protein_isolate']
    naphtha_product.price = price['Naptha']
    green_diesel.price = price['Green_diesel']
    propane.price = price['Propane']

    if BT is not None:
        BT.natural_gas_price = price['Natural gas']
        if len(BT.ins) > 4:
            BT.ins[4].price = price['Lime']


OSBL_units = bst.get_OSBL(oil_extraction_sys.cost_units)

saf_tea = create_cellulosic_ethanol_tea(
    oil_extraction_sys,
    OSBL_units=OSBL_units,
)
saf_tea.IRR = 0.10
saf_tea.operating_days =100 # availbility of covercress is in may-june so max 60 days
saf_tea.income_tax = 0.21
saf_tea.duration = (2016, 2046)
saf_tea.labor_cost = 3212962 * get_flow_dry_tpd() / 2205
if BT is not None:
    saf_tea.boiler_turbogenerator = BT


def solve_scenario(stream_to_solve, label=''):
   
    for _ in range(num_sims):
        oil_extraction_sys.simulate()
    for _ in range(num_solve_tea):
        stream_to_solve.price = saf_tea.solve_price(stream_to_solve)
    print(f'\n{label} ')
    print(f'  Stream     : {stream_to_solve.ID}')
    print(f'  Break-even : ${stream_to_solve.price:.4f}/kg')
    print(f'  NPV        : ${saf_tea.NPV:,.0f}')
    return stream_to_solve.price


# def get_SAF_MPSP():
#     return solve_scenario(SAF, label='Scenario A: SAF MPSP')
SAF_liquid_density= 0.78 
def saf_price_per_L(price_per_kg, rho_kg_per_L=SAF_liquid_density):
    """Convert SAF price from $/kg to $/L at liquid reference density."""
    return price_per_kg * rho_kg_per_L
def get_SAF_MPSP():
    mpsp_kg = solve_scenario(SAF, label='Scenario A: SAF MPSP')
    mpsp_L = saf_price_per_L(mpsp_kg)
    print(f'  Break-even : ${mpsp_L:.4f}/L  (${mpsp_kg:.4f}/kg)')
    return mpsp_kg, mpsp_L   # or return mpsp_L only if you prefer
def get_feedstock_MFPP():
    saved_saf_price = SAF.price
   
    for _ in range(num_sims):
        oil_extraction_sys.simulate()
    for _ in range(num_solve_tea):
        feedstock.price = saf_tea.solve_price(feedstock)
    mfpp = feedstock.price
    feedstock.price = price['Feedstock']
    SAF.price = saved_saf_price
    print('\nScenario B: feedstock MFPP (separate from VOC above)')
    print(f'  Stream              : {feedstock.ID}')
    print(f'  Max feedstock price : ${mfpp:.4f}/kg')
    print(f'  Reference feedstock : ${price["Feedstock"]:.4f}/kg')
    print(f'  NPV @ MFPP solve    : ${saf_tea.NPV:,.0f}')
    return mfpp



def print_tea_summary(MSP_SAF):
    MSP_SAF_L = saf_price_per_L(MSP_SAF)
    market_SAF_L = saf_price_per_L(price['SAF'])
    print('\n TEA SUMMARY ')
    print(f'FCI                 : ${saf_tea.FCI/1e6:.2f} MM$')
    print(f'Material cost       : ${saf_tea.material_cost/1e6:.2f} MM$/yr')
    print(f'Utility cost        : ${saf_tea.utility_cost/1e6:.2f} MM$/yr')
    print(f'AOC                 : ${saf_tea.AOC/1e6:.2f} MM$/yr')
    print(f'\nSAF MPSP            : ${MSP_SAF_L:.4f}/L')
    print(f'Market SAF          : ${market_SAF_L:.4f}/L')
    # print(f'  Profitable?       {MSP_SAF <= price["SAF"]}')
    # print(f'\nMax feedstock price : ${MFPP:.4f}/kg')
    print(f'Reference feedstock : ${price["Feedstock"]:.4f}/kg')
    # print(f'  Affordable?       {MFPP >= price["Feedstock"]}')
 # def print_tea_summary(MSP_SAF, MFPP):
 #     MSP_SAF_L = saf_price_per_L(MSP_SAF)
 #     market_SAF_L = saf_price_per_L(price['SAF'])
 #     print('\n TEA SUMMARY ')
 #     print(f'FCI                 : ${saf_tea.FCI/1e6:.2f} MM$')
 #     print(f'Material cost       : ${saf_tea.material_cost/1e6:.2f} MM$/yr')
 #     print(f'Utility cost        : ${saf_tea.utility_cost/1e6:.2f} MM$/yr')
 #     print(f'AOC                 : ${saf_tea.AOC/1e6:.2f} MM$/yr')
 #     print(f'\nSAF MPSP            : ${MSP_SAF_L:.4f}/L')
 #     print(f'Market SAF          : ${market_SAF_L:.4f}/L')
 #     # print(f'  Profitable?       {MSP_SAF <= price["SAF"]}')
 #     # print(f'\nMax feedstock price : ${MFPP:.4f}/kg')
 #     print(f'Reference feedstock : ${price["Feedstock"]:.4f}/kg')
 #     # print(f'  Affordable?       {MFPP >= price["Feedstock"]}')
      


def run_tea():
    set_market_prices() #for the basleine prices as inital guess
    oil_extraction_sys.simulate()
    print(f'FCI: ${saf_tea.FCI/1e6:.2f} MM$')

    MSP_SAF_kg, MSP_SAF_L = get_SAF_MPSP()
    MFPP = get_feedstock_MFPP()
    # print_tea_summary(MSP_SAF_kg, MFPP)
    # return MSP_SAF_L, MFPP   # primary result in $/L
    print_tea_summary(MSP_SAF_kg)
    op_hr = saf_tea.operating_hours
    print(f'\nTotal material: ${saf_tea.material_cost/1e6:.3f} MM$/yr  '
          f'(${saf_tea.material_cost/op_hr:.1f} USD/hr)')
    print('\n--- Material cost by feed stream ---')
    for s in sorted(oil_extraction_sys.feeds, key=lambda x: -x.cost):
        if s.price and s.F_mass > 0 and abs(s.cost) > 1e-6:
            print(f'  {s.ID:22s}  {s.F_mass:10.1f} kg/hr  '
                  f'${s.price:8.4f}/kg  ${s.cost:8.2f}/hr  '
                  f'${s.cost*op_hr/1e6:6.3f} MM$/yr')
    MFPP = get_feedstock_MFPP()
    return MSP_SAF_L, MFPP
    


MSP_SAF, MFPP = run_tea()
# MSP_SAF = run_tea()
from biorefineries.covercress import tea as cc_tea
cc_tea.setup(oil_extraction_sys, saf_tea, BT)
cc_tea.TEA_breakdown(print_output=True, fractions=False)
cc_tea.export_excel(fraction=False)
cc_tea.plot_breakdown(fraction=True, show=False)
# %% LCA — SAF GWP (succinic pattern)
from biorefineries.covercress.lca import CovercressLCA

# Functional unit: kg SAF/h (all hydrocarbons in SAF stream)
SAF_product_IDs = [ID for ID in SAF.chemicals.IDs if SAF.imass[ID] > 1e-9]

CT = next((ui for ui in oil_extraction_sys.units if isinstance(ui, bst.CoolingTower)), None)
CWP = next((ui for ui in oil_extraction_sys.units if isinstance(ui, bst.ChilledWaterPackage)), None)

covercress_LCA = CovercressLCA(
    system=oil_extraction_sys,
    CFs=CFs,
    feedstock=feedstock,
    feedstock_ID='Feedstock',
    input_biogenic_carbon_streams=[feedstock],   # biogenic C from seed
    main_product=SAF,
    main_product_chemical_IDs=SAF_product_IDs,
    by_products=[protein_isolate, naphtha_product, green_diesel, propane],
    boiler=BT,
    cooling_tower=CT,
    chilled_water_processing_units=[CWP] if CWP else [],
    has_turbogenerator=BT.power_utility.production > 0 if BT else False,
    add_EOL_GWP=True,   # coproduct carbon credited at end-of-life (succinic default)
)


def print_carbon_and_gwp(lca=None):
    """Collect C from feeds, products, emissions; print GWP breakdown."""
    if lca is None:
        lca = covercress_LCA
    CO2_MW = lca.chemicals.CO2.MW
    fq = lca.functional_quantity_per_h  # kg SAF/h

    def co2eq_rate(stream):
        return stream.get_atomic_flow('C') * CO2_MW

    C_in = sum(s.get_atomic_flow('C') for s in lca.feeds)
    C_products = sum(s.get_atomic_flow('C') for s in lca.products)
    C_emissions = sum(s.get_atomic_flow('C') for s in lca.emissions)

    print('\n--- Carbon balance ---')
    print(f'C in (feeds)      : {C_in:.1f} mol-C/h')
    print(f'C in products     : {C_products:.1f} mol-C/h')
    print(f'C in emissions    : {C_emissions:.1f} mol-C/h')
    print(f'Closure (out/in)  : {(C_products + C_emissions) / C_in:.4f}')
    print(f'SAF flow          : {fq:.2f} kg/h')

    print('\n--- Emission streams (CO2-eq from C) ---')
    for stream in lca.emissions:
        if stream.get_atomic_flow('C') > 1e-6:
            print(f'  {stream.ID:30s}  {co2eq_rate(stream):10.1f} kg CO2-eq/h')

    print('\n--- GWP per kg SAF ---')
    print(f'  Total GWP              : {lca.GWP:.3f} kg CO2-eq/kg')
    print(f'  Feedstock              : {lca.feedstock_GWP:.3f}')
    print(f'  Materials              : {lca.material_GWP:.3f}')
    print(f'  Net electricity        : {lca.net_electricity_GWP:.3f}')
    print(f'  Direct emissions (gross): {lca.direct_emissions_GWP:.3f}')
    print(f'  Non-biogenic (+ EOL)   : {lca.direct_non_biogenic_emissions_GWP:.3f}')
    print(f'  System C balance (LCA) : {lca.system_carbon_balance:.4f}')


oil_extraction_sys.simulate()
print_carbon_and_gwp()
print(f'\nSAF GWP = {covercress_LCA.GWP:.3f} kg CO2-eq/kg SAF')