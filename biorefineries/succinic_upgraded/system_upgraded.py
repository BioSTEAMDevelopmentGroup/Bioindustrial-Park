#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 17:39:07 2025

@author: princyk2
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 19:22:32 2025

@author: princyk2
"""


import math
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
IQ_interpolation = flx.IQ_interpolation
from biorefineries.succinic.system_glucose import get_succinic_acid_solubility_gpL
from biorefineries.succinic_upgraded import units_adsorption_upgraded
from biorefineries.succinic import facilities
from biorefineries.succinic_upgraded.process_settings import price, CFs, chem_index, _GDP_2007_to_2010
from biorefineries.succinic.utils import find_split, splits_df, baseline_feedflow

from biorefineries.succinic.tea import TemplateTEA as SuccinicTEA
from biorefineries.succinic.lca import SuccinicLCA
from biorefineries.cellulosic import create_facilities
from biorefineries.succinic.chemicals_data import chems, chemical_groups, \
                                                              soluble_organics, combustibles
# TODO  changes chems to succininc_upgraded                                                              # 

from biorefineries import succinic
bst.preferences.update(flow = 'kg/hr',composition = True)
HeatExchangerNetwork = bst.facilities.HeatExchangerNetwork 
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
bst.preferences.N = 300
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0
bst.CE = bst.units.design_tools.CEPCI_by_year[2019]
tmo.settings.set_thermo(chems)

feedstock_ID = 'Glucose'

# System settings
System.converge_method = 'wegstein'
System.default_maxiter = 100
System.default_molar_tolerance = 0.1
System.default_relative_molar_tolerance = 0.0001 # supersedes absolute tolerance
System.strict_convergence = True 

succinic.load(mode='pilot_batch_upgraded_glucose',machine='mac')
base_sys = succinic.system

System.default_relative_molar_tolerance = 0.0001 # supersedes absolute tolerance
System.strict_convergence = True # True => throw exception if system does not converge; False => continue with unconverged system
flowsheet = bst.Flowsheet('succinic_upgraded')
bst.main_flowsheet.set_flowsheet(flowsheet)

@SystemFactory(ID = 'succinic_sys_new')
def create_new_adsorption_units_only_sys(ins, outs):
    def fix_split(isplit, ID):
        isplit['SuccinicAcid', 'LacticAcid', 'Ethanol'] = isplit[ID]
    process_groups = []
    f1=bst.main_flowsheet
    u_1 = f1.unit
    s_1 = f1.stream
    # Access units directly from base_sys (loaded system)
    feedstock = bst.Stream('glucose_feedstock', Glucose=1., Water=1., units='kmol/h')
    product_stream = bst.Stream('SuccinicAcid')
    feedstock.price = price['Glucose']*0.909 # dextrose monohydrate stream is 90.9 wt% glucose
    feedstock.F_mass = 200_000 # initial value; updated by spec.set_production_capacity
     
    # Simplify base_sys.flowsheet.unit access
    u = base_sys.flowsheet.unit
    s = base_sys.flowsheet.stream
    
    # Connect upstream units using pipe notation
    feedstock-0-u.U101
    u.U101-0-0-u.M201
    u.M201-0-0-u.U302
    u.U302-0-0-u.S301
    u.S301-0-0-u.H302
    u.H302-0-0-u.P303
    u.S301-1-0-u.F301
    u.F301-0-0-u.H301
    (u.H301-0, s.dilution_water)-u.M304
    u.M304-0-u.M304_H
    u.M304_H-0-0-u.S302
    (s.CO2_fermentation)-u.S303

    seed = bst.Stream('seed', units='kg/hr')
    A301_makeup_water = bst.Stream('A301_makeup_water', units='kg/hr')
    CSL = Stream('CSL', units='kg/hr', price=price['CSL'])
    # Register CSL in flowsheet for material cost tracking
    s_1.data['CSL'] = CSL
    base_neutralization = Stream('base_neutralization', units='kg/hr',)
    magnesiumchloride_hexahydrate_fermentation =Stream('magnesiumchloride_hexahydrate_fermentation', units='kg/hr',
                                     price=price['MagnesiumChlorideHexahydrate'],) 
    # Register nutrients in flowsheet for material cost tracking
    s_1.data['magnesiumchloride_hexahydrate_fermentation'] = magnesiumchloride_hexahydrate_fermentation
    zincsulfate_heptahydrate_fermentation =Stream('zincsulfate_heptahydrate_fermentation', units='kg/hr', 
                              price=price['Zinc Sulfate Heptahydrate'],)
    s_1.data['zincsulfate_heptahydrate_fermentation'] = zincsulfate_heptahydrate_fermentation
    ActivatedCarbon= Stream('ActivatedCarbon',units='kg/hr', price=price['ActivatedCarbon'])
    s_1.data['ActivatedCarbon'] = ActivatedCarbon
   
    R302_1 = units_adsorption_upgraded.CoFermentation('R302_1', 
                                    ins=(
                                        u.S302-1, 'seed', CSL, u.S303-1, base_neutralization, 
                                         
                                         # M300_R-0, 
                                         u.K301-0, # CO2 recycle (K301-0)
                                         magnesiumchloride_hexahydrate_fermentation,
                                         zincsulfate_heptahydrate_fermentation, 
                                        # diammonium_sulfate_fermentation,, 
                                       # magnesium_sulfate_fermentation,# 
                                         u.K302-0,
                                         ), # air (K302-0)
                                    outs=('fermentation_broth', 'fermentation_vent'),
                                    mode='batch',
                                    V_wf=0.33,
                                    neutralization=False,
                                    pH_control=True,
                                    neutralizing_agent='Lime',
                                    base_neutralizes_product=False, # if True, acidulation is performed downstream
                                    )
     #TODO didnt change the pH control used Lime instead of NH4OH
     #inlcuded teh model of tactical CO2 dosing
  
    @R302_1 .add_specification(run=False)
    def R302_1_spec(): # note: effluent always has 0 CSL
     # R302.show(N=100)
        R302_1 ._run()
        if R302_1 .ins[2].F_mol:
            R302_1 .ins[2].F_mass*=1./(1-u.S302.split[0])
            R302_1 ._run()
            u.S303.simulate()
            u.K302.specifications[0]()
    (u.S302-0, u.S303-0)-u.R303
    
    (R302_1 -1, u.R303-1)-u.M305
    u.R303-0-0-u.T301
    u.T301-0-1-R302_1 
    makeup_MEA_A301 = Stream('makeup_MEA_A301', units='kg/hr', price=price['Monoethanolamine'])
    # Register MEA in flowsheet for material cost tracking
    s_1.data['makeup_MEA_A301'] = makeup_MEA_A301
    # Note: u.makeup_MEA_A301 might be used, but we need to use the new one
    # Check if base_sys has makeup_MEA_A301 and update reference
    A301_1 = bst.AmineAbsorption('A301_1', ins=(u.M305-0, makeup_MEA_A301, 'A301_makeup_water'), outs=('absorption_vent', 'captured_CO2'),
                               CO2_recovery=0.52)
   
    
    def A301_1_obj_f(CO2_recovery):
        A301_1.CO2_recovery = CO2_recovery
        A301_1._run()
        u.K301.specifications[0]()
        R302_1.specifications[0]()
        return R302_1.fresh_CO2_required
    
    @A301_1.add_specification(run=False)
    def A301_1_spec():
        # A301.outs[1].phase='g'
        A301_1._run()
        if A301_1_obj_f(1-1e-3)>0.:
            pass
        else:
            IQ_interpolation(A301_1_obj_f, 1e-3, 1-1e-3, x=0.5, ytol=1e-4)
        A301_1.outs[1].phase='g'
           
    A301_1-1-0-u.K301
       
    R302_1 -0-0-u.H401
    u.H401-0-0-u.M401
    S401_1= bst.SolidsCentrifuge('S401_1', ins=u.M401-0, 
                        outs=('S401_1_solid_waste', 'S401_1_1'),
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
     
    @S401_1.add_specification(run=False)
    def S401_1_succinic_acid_precipitation_spec():
        S401_1._run()
        
        instream = S401_1.ins[0]
        tot_SA = instream.imass['SuccinicAcid']
        dissolved_SA = get_succinic_acid_solubility_gpL(instream.T) * instream.F_vol
        S401_1.outs[1].imass['SuccinicAcid'] = min(tot_SA, dissolved_SA)
        S401_1.outs[0].imass['SuccinicAcid'] = max(0, tot_SA - dissolved_SA)
        liquid = S401_1.outs[1]
        microbe_stream = liquid.copy()

        microbe_stream.imass['FermMicrobe'] = liquid.imass['FermMicrobe']

        N_flow = microbe_stream.get_atomic_flow('N')

        l_phenylalanine = N_flow / 1

        liquid.imass['FermMicrobe'] = 0
        liquid.imol['l-phenylalanine'] += l_phenylalanine

        S401_1.l_phenylalanine= l_phenylalanine
    
    H402 = units_adsorption_upgraded.HXutility('H402', ins=S401_1-1, outs = ('broth_to_adsorbtion',), T=20. + 273.15)
    
    #TODO temp is changed from operational temp 40 to 20 ˚C becuase data for 20 ˚ C is only available)
    recycled_water_A401 =bst.Stream(ID='recycled_water_A401',units='kg/hr', price=price['Makeup water'],)
    recycled_water_F400_P=bst.Stream(ID='recycled_water_F400_P',units='kg/hr', price=price['Makeup water'],)
    make_up_water = Stream('make_up_water', water=1,units='kg/hr', price=price['Makeup water'],)   
    NaOH_1=bst.Stream('NaOH_1', units='kg/hr',  price=price['NaOH'])
    T400 = units_adsorption_upgraded.CausticStorage('T400', ins = NaOH_1, outs = 'fresh_NaOH')
    @T400.add_specification(run=True)
    def T400_spec():
        T400.ins[0].imass['NaOH']=A401.outs[1].imass['NaOH']
    
    M400=bst.MixTank('M400',ins= (make_up_water,T400-0, ''), outs= '1.0M_NaOH_solution')
    
    
    
    make_up_AC=bst.Stream('Make_Up_adsorbent',units='kg/hr', price=price['ActivatedCarbon'])

    
    A401=units_adsorption_upgraded.SingleComponentAdsorptionColumn('A401',ins=(H402-0,M400-0,make_up_AC), outs=('broth_post_adsorption','spent_fluid','spent_adsorbent' ),
                                adsorbate = ['l-phenylalanine'],
                                                  # [3]contains **{'l-phenylalanine'} adsorption study so for this study i used **{'l-phenylalanine} as aminoacid present in the dextrose 
                                                  # I didnt find any sources that evident to prove vitamins A and D present in th dextrose.
                                cycle_time= 3, superficial_velocity=14,
                                # assume superficial velocity as in the range between ( velocity fro liquids typically between 0.001-0.004 m/s)
                                
                                
                                # regeneration_superficial_velocity = None, #TODO: this does not affect the design and cost
                                isotherm_model='Langmuir', isotherm_args=(5.4,64.4),
                            # (5.4,64.4)
                                #  adsorption model-[9]
                                regeneration_fluid={'T':20+273.15,'NaOH': 0.018 , 'Water': 1-0.018},
                                            
                                                   
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
    regeneration_isotherm_args=(1,1000),
    regeneration_isotherm_model='Langmuir', adsorbent='ActivatedCarbon',) # [1 / hr]k=0.3,)
       
    
       
    @M400.add_specification(run=True)
    def M400_spec():
        T400.run()
        M400.ins[0].phase = 'l'
        M400.ins[2].phase = 'l'
        M400.ins[1].imass['NaOH']=A401.outs[1].imass['NaOH']
        
        
        Required_Water_mass =   A401.outs[1].imass['Water']
        M400.outs[0].imass['Water'] = Required_Water_mass
        recycled_water = M400.ins[2].imass['Water']
        total_water_req = Required_Water_mass
        amt_of_fresh_water = max(0, total_water_req - recycled_water)
        M400.ins[0].imass['Water'] = amt_of_fresh_water
        
        M400.run()    
    
       
    @A401.add_specification(run=True)
    def A401_spec(): 
        A401._run()
        #adsorption unit module calculte regeneration velocity not based on chemistry the calculation is based on inventories, length, area adn regneration vleocity, cycle time =regenration time 
        loss_effectiveness=A401.loss_effectiveness #refer[]
        no_of_cycles=A401.no_cycles
        
        total_loss_fraction=loss_effectiveness/(A401.cycle_time*no_of_cycles)
        Total_mass_AC = float(A401.area) * float(A401.column_length) * A401.rho_adsorbent
         
           
        make_up_AC_flow= Total_mass_AC * total_loss_fraction
        A401.ins[2].imass['ActivatedCarbon']=make_up_AC_flow
        A401.run() 
    
    
    S400_W = bst.Splitter('S400_W', ins= A401-1, outs=('NaOH_distillation', 'water_wwt'),split=0.25)
    
    F400 = bst.MultiEffectEvaporator('F400', ins=S400_W-0, outs=('F400_l', 'F400_g'),
                                      P = (101325, 73581, 50892, 32777, 20000), V_definition='Overall', V=0.1) 
    
    F400_P =bst.Pump('F400_P', ins=F400-1, outs= recycled_water_F400_P, P=101325)
    
    F400_P-0-2-M400 
   
    # sys_adsorption=succinic_sys.from_units(ID='sys_adsorption',units=(u_1.R302_1, u_1.A301_1, u_1.S401_1, u_1.H402, u_1.T400, u_1.M400, u_1.A401, u_1.S400_W, u_1.F400, u_1.F400_P))
                             
   
    succinic_downstream = bst.System.from_units(ID = 'downstream_sys', units= [u.R401, u.R401_P, u.S405, u.M404, u.F401, u.F401_P, u.C401,\
                                                                        u.S402, u.F402, u.F402_P, u.C402, u.S403, u.F403, u.F403_P,\
                                                                        u.C403, u.S404, u.M402, u.F404, u.M403,\
                                                                        u.S406, u.M501, u.M502, u.U501, u.R501,\
                                                                        u.M503, u.R502, u.R503,\
                                                                        u.S501, u.M504, u.C501, u.M505, u.U502,\
                                                                        u.M510, u.M405, u.T601, u.T601_P,\
                                                                        u.T602, u.T603,\
                                                                        u.BT701,\
                                                                        u.CT801, u.CWP802, u.CIP901,\
                                                                         u.ADP902, u.FWT903,  u.PWC904])
      
                                    
                                                                  
          
  
     
    HXN = HeatExchangerNetwork('HXN1001',
                                              # ignored=[H401, H402],
                                              cache_network=True,
                                              )
    def HXN_no_run_cost():
        HXN.heat_utilities = []
        HXN._installed_cost = 0.
 
    (A401-0, s.sulfuric_acid_R401)-u.R401
    
 
  
    # Create adsorption subsystem with new units
    sys_adsorption = bst.System.from_units(ID='sys_adsorption', units=[R302_1, A301_1, S401_1, H402, T400, M400, A401, S400_W, F400, F400_P])
    
    # Create upstream system using pipe notation path
    succinic_upstream = bst.System(ID='upstream_sys', path=(u.U101, u.M201, u.U302, u.S301, u.H302, u.P303,
                                                             u.F301, u.H301, u.M304, u.M304_H, u.S302, u.S303, u.R303,
                                                             u.M305, u.T301, u.K301, u.K302, u.H401, u.M401))
    
    # Collect all units from subsystems and new units
    combined_units = list(succinic_upstream.units)
    combined_units.extend(sys_adsorption.units)
    combined_units.extend(succinic_downstream.units)
    combined_units.append(HXN)
    
    # Remove duplicates while preserving order
    seen = set()
    unique_units = []
    for unit in combined_units:
        if unit.ID not in seen:
            seen.add(unit.ID)#record hte ID of unit
            unique_units.append(unit)# keep the unit
    
    # Build the final system from all units
    # BioSTEAM will automatically detect subsystems
    succinic_sys_new = bst.System.from_units(
        ID='succinic_sys_new',
        units=unique_units
    )
    
    # Don't overwrite flowsheet.system as it's the system registry
    # The system is automatically registered via SystemFactory
    return succinic_sys_new
 

succinic_sys_new = create_new_adsorption_units_only_sys()

succinic_sys_new.simulate()

f1 = bst.main_flowsheet

def get_all_units_from_system(system):  
    units = list(system.units)
    # for subsystem in system.subsystems:
    #     units.extend(get_all_units_from_system(subsystem))
    return units

# Register all units from the new system (succinic_sys_new)
for unit in get_all_units_from_system(succinic_sys_new):
    f1.unit.data[unit.ID] = unit

# Register all units from base_sys that aren't already registered
# This makes them accessible via u_1 (f1.unit) for consistency
for unit in get_all_units_from_system(base_sys):
    if unit.ID not in f1.unit.data:
        f1.unit.data[unit.ID] = unit

u_1 = f1.unit
s_1 = f1.stream
feedstock = s_1.glucose_feedstock
# Product stream should be the output of T601_P (pump after storage tank)

product_stream = u_1.T601_P.outs[0]
# 
# Register the product stream with ID 'SuccinicAcid' in the flowsheet for easy access
s_1.data['SuccinicAcid'] = product_stream
product_stream.register_alias('SuccinicAcid')
feeds = succinic_sys_new.feeds
products = [product_stream]
 # Don't include gypsum since we want to include carbon impurities in GWP calculation


u_1.BT701.ID = 'BT701'
# u_1.CT901.ID = 'CT801'
u_1.CT801.ID = 'CT801'
u_1.CWP802.ID = 'CWP802'
u_1.CIP901.ID = 'CIP901'
u_1.ADP902.ID = 'ADP902'
u_1.FWT903.ID = 'FWT903'
u_1.PWC904.ID = 'PWC904'

# Access BT from unit registry instead of using flowsheet('BT')
# since f1.system is now a System object, not a registry
BT = u_1.BT701  # BT701 is the boiler unit



BT.natural_gas_price = 0.2527
BT.ins[4].price = price['Lime']

globals().update(f1.to_dict())


# %% keeep this commented till other results are fixed
# =============================================================================
# TEA
# =============================================================================

# Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legislation)

get_flow_dry_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
succinic_tea = SuccinicTEA(system=succinic_sys_new , IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, 
        # operating_days=0.9*365,
        operating_days = 200,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(
                    u_1.U501,
                    # f_old.T601, f_old.T602, 
                    # f_old.T601, f_old.T602, f_old.T603, f_old.T604,
                    # f_old.T606, f_old.T606_P,
                    u_1.BT701, u_1.CT801, u_1.CWP802, u_1.CIP901, u_1.ADP902, u_1.FWT903, u_1.PWC904,
                    ),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_dry_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u_1.BT701)

seed_train_system = bst.System('seed_train_system', path=(u_1.S302, u_1.R303, u_1.T301))

theoretical_max_g_succinic_acid_per_g_glucose = (2*chems.SuccinicAcid.MW)/(chems.Glucose.MW) # ~1.311 g-succinic-acid/g-glucose

spec = ProcessSpecification(
    evaporator = u_1.F301,
    pump = None,
    mixer = u_1.M304,
    heat_exchanger = u_1.M304_H,
    seed_train_system = [],
    seed_train = u_1.R303,
    reactor= u_1.R302_1,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('SuccinicAcid',),
    
    # spec_1=0.19,
    spec_1=0.7,
    spec_2=100.,
    spec_3=1.0,
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = None,
    byproduct_streams = [],
    HXN = u_1.HXN1001,
    maximum_inhibitor_concentration = 1.,
    pre_conversion_units = succinic_sys_new.split(u_1.H301.outs[0])[0],
    
    # set baseline fermentation performance here
    # baseline_yield = 0.4478, # 0.587 g/g-glucose
    # baseline_yield = 0.3738, # 0.49 g/g-glucose
    # baseline_yield = 0.4733/theoretical_max_g_succinic_acid_per_g_glucose, # 36.1% of theoretical maximum
    # baseline_titer = 63.1,
    # baseline_productivity = 0.6573,# for glucose
    
    baseline_yield = 0.65/theoretical_max_g_succinic_acid_per_g_glucose, # 49.6 % of theoretical maximum
    # 65 % actual yield
    baseline_titer = 78.37,
    baseline_productivity = 0.96,#fro dextrose
    neutralization=False,
    baseline_pyruvic_acid_yield = 0.006,
    baseline_cell_mass_yield = 0.186,
    #TODO check if baseline pyruvic and cell mass is needed

    feedstock_mass = feedstock.F_mass,
    pretreatment_reactor = None)


def M304_titer_obj_fn(water_to_sugar_mol_ratio):
    # M304, R302 = u_1.M304, u_1.R302_1
    u_1.M304.water_to_sugar_mol_ratio = water_to_sugar_mol_ratio
    u_1.M304.specifications[0]()
    u_1.M304_H._run()
    u_1.S302._run()
    u_1.R303.specifications[0]()
    u_1.T301._run()
    # R302_1.simulate()
    u_1.R302_1.specifications[0]()
    return u_1.R302_1.effluent_titer - u_1.R302_1.titer_to_load

def F301_titer_obj_fn(V):
    u_1.F301.V = V
    u_1.F301._run()
    u_1.H301._run()
    u_1.M304.specifications[0]()
    u_1.M304_H._run()
    u_1.S302._run()
    u_1.R303.specifications[0]()
    u_1.T301._run()
    # R302_1.simulate()
    u_1.R302_1.specifications[0]()
    return u_1.R302_1.effluent_titer - u_1.R302_1.titer_to_load


def load_titer_with_glucose(titer_to_load, set_F301_V=None):
    # clear_units([V301, K301])
    F301_ub = 0.8
    F301_lb = 0.5 if set_F301_V is None else set_F301_V
    M304_lb, M304_ub = 0., 100_000.  # for low-titer high-yield combinations, if infeasible, use a higher upper bound
    
    spec.spec_2 = titer_to_load
    u_1.R302_1.titer_to_load = titer_to_load
    F301_titer_obj_fn(F301_lb)
    
    if M304_titer_obj_fn(M304_lb) < 0.: # if there is too low a conc even with no dilution
        # IQ_interpolation(F301_titer_obj_fn, F301_lb, F301_ub, ytol=1e-3)
        spec.titer_inhibitor_specification.check_sugar_concentration() 
    # elif F301_titer_obj_fn(1e-4)>0: # if the slightest evaporation results in too high a conc
    elif M304_titer_obj_fn(M304_ub) > 0.:
        IQ_interpolation(M304_titer_obj_fn, 
                         M304_lb,
                         M304_ub, 
                         ytol=1e-3)
    else:
        F301_titer_obj_fn(F301_lb)
        IQ_interpolation(M304_titer_obj_fn, 
                         M304_lb, 
                         M304_ub, 
                         ytol=1e-3)
        
    if not feedstock_ID=='Glucose' and (set_F301_V is None): spec.titer_inhibitor_specification.check_sugar_concentration()
    
spec.load_spec_1 = spec.load_yield
spec.load_spec_3 = spec.load_productivity
spec.load_spec_2 = load_titer_with_glucose





# %% Diagram
import biosteam as bst
bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = True
succinic_sys_new .diagram('cluster')

#%% Define unit groups

area_names = [
    'feedstock_handling',
    'feedstock_preprocessing',
    #'feedstock',
    'conversion',
    'separation',
    # 'recovery',
    'wastewater system',
    'storage',
    'boiler & turbogenerator',
    'cooling utility facilities',
    'other facilities',
    'heat exchanger network',
]

# area_names = [
#     'feedstock',
#     'conversion',
#     'separation',
#     'wastewater',
#     'storage',
#     'boiler & turbogenerator',
#     'cooling utility facilities',
#     'other facilities',
#     'heat exchanger network',
# ]
# I changed area_name 'wastewater in the system glucose --> recovery 
unit_groups = bst.UnitGroup.group_by_area(succinic_sys_new .units)


unit_groups.append(bst.UnitGroup('natural gas (for steam generation)'))
unit_groups.append(bst.UnitGroup('natural gas (for product drying)'))

unit_groups.append(bst.UnitGroup('fixed operating costs'))


for i, j in zip(unit_groups, area_names): i.name = j
for i in unit_groups: i.autofill_metrics(shorthand=False, 
                                         electricity_production=False, 
                                         electricity_consumption=True,
                                         material_cost=True)
for i in unit_groups:
    if i.name == 'storage' or i.name=='other facilities' or i.name == 'cooling utility facilities':
        i.metrics[-1].getter = lambda: 0. # Material cost
    if i.name == 'cooling utility facilities':
        i.metrics[1].getter = lambda: 0. # Cooling duty
    if i.name == 'boiler & turbogenerator':
        i.metrics[-1] = bst.evaluation.Metric('Material cost', 
                                            getter=lambda: succinic_tea.utility_cost/succinic_tea.operating_hours, 
                                            units='USD/hr',
                                            element=None) # Material cost
        # i.metrics[-2].getter = lambda: BT.power_utility.rate/1e3 # Electricity consumption [MW]
        
for HXN_group in unit_groups:
    if HXN_group.name == 'heat exchanger network':
        HXN_group.filter_savings = False
        # assert isinstance(HXN_group.units[0], HeatExchangerNetwork) #!!! Note groups aren't working for non-sugarcane succinic acid biorefineries, fix if needed


unit_groups[-3].metrics[-1] = bst.evaluation.Metric('Material cost', 
                                                    getter=lambda: BT.natural_gas_price * BT.natural_gas.F_mass, 
                                                    units='USD/hr',
                                                    element=None)

unit_groups[-2].metrics[-1] = bst.evaluation.Metric('Material cost', 
                                                    getter=lambda: u_1.F404.ins[2].cost, 
                                                    units='USD/hr',
                                                    element=None)

unit_groups[-1].metrics[-1] = bst.evaluation.Metric('Material cost', 
                                                    getter=lambda: succinic_tea.FOC/succinic_tea.operating_hours, 
                                                    units='USD/hr',
                                                    element=None)

# for i in unit_groups:
#     i.metric(i.get_net_electricity_production,
#                 'Net electricity production',
#                 'kW')




unit_groups_dict = {}
for i in unit_groups:
    unit_groups_dict[i.name] = i

# Override material cost calculation for conversion group to include external feeds
# (CSL, nutrients, MEA) that don't have source units
# conversion_group = unit_groups_dict.get('conversion')
# if conversion_group:
#     def conversion_material_cost():
#         """Calculate material cost from all input streams in conversion units."""
#         total_cost = 0
        
#         # Sum material costs from all input streams
#         for unit in conversion_group.units:
#             for stream in unit.ins:
#                 if stream and hasattr(stream, 'price') and hasattr(stream, 'F_mass'):
#                     if stream.price and stream.F_mass:
#                         total_cost += stream.price * stream.F_mass
        
#         return total_cost
    
#     conversion_group.metrics[-1].getter = conversion_material_cost

BT = u_1.BT701
CT = u_1.CT801
CWP = u_1.CWP802
HXN = u_1.HXN1001

# HXN1001.force_ideal_thermo = True
# %%
# =============================================================================
# Simulate system and get results
# =============================================================================

num_sims = 3
num_solve_tea = 3
def get_product_stream_MPSP():
    for i in range(num_sims):
        succinic_sys_new .simulate()
    for i in range(num_solve_tea):
        # Use product_stream (T601_P.outs[0]) which is the actual system product
        # flowsheet.stream.SuccinicAcid should now point to the same stream
        product_stream.price = succinic_tea.solve_price(product_stream)
        # product_stream.price = 2.7 
        # to check the NPV value
        
        
    return product_stream.price * product_stream.F_mass / product_stream.imass['SuccinicAcid']

def simulate_and_print():
    MPSP = get_product_stream_MPSP()
    GWP = succinic_LCA.GWP
    FEC = succinic_LCA.FEC
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${MPSP:.3f}/kg')
    print(f'GWP100a is {GWP:.3f} kg CO2-eq./kg')
    print(f'FEC is {FEC:.3f} MJ/kg')
    print('----------------------------------------\n')

get_product_stream_MPSP()
succinic_tea.labor_cost=3212962*get_flow_dry_tpd()/2205
get_product_stream_MPSP()




# # #%% LCA

# Get the actual feedstock from system feeds to ensure it's in the system
system_feedstock = feedstock if feedstock in succinic_sys_new.feeds else (succinic_sys_new.feeds[0] if succinic_sys_new.feeds else s_1.glucose_feedstock)
succinic_LCA = SuccinicLCA(system=succinic_sys_new , 
                 CFs=CFs, 
                 feedstock=system_feedstock, 
                 feedstock_ID=feedstock_ID,
                 input_biogenic_carbon_streams=[system_feedstock, s_1.CSL] if hasattr(s_1, 'CSL') else [system_feedstock],
                 main_product=product_stream, 
                 main_product_chemical_IDs=['SuccinicAcid',], 
                 by_products=[], 
                 cooling_tower=u_1.CT801, 
                 chilled_water_processing_units=[u_1.CWP802,], 
                 boiler=u_1.BT701, has_turbogenerator=True,
                 add_EOL_GWP=True,
                 )

#%% Simulate and print
spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)
get_product_stream_MPSP()

# Verify conversion material cost after simulation( why we used this code snippet here is the material cost calclation was not reflecting the entire materials price.)
# if conversion_group:
#     # Re-apply the override after simulation to ensure it uses updated flows
#     def conversion_material_cost_updated():
#         total_cost = 0
#         for unit in conversion_group.units:
#             for i, stream in enumerate(unit.ins):
#                 if stream and hasattr(stream, 'price') and stream.price and hasattr(stream, 'F_mass'):
#                     stream_cost = stream.price * stream.F_mass
#                     if stream_cost > 0.01:
#                         total_cost += stream_cost
        
#         # Add Lime cost for pH control in R302_1 (Lime is consumed during unit operation)
#         try:
#             if hasattr(u_1, 'R302_1'):
#                 r302_1 = u_1.R302_1
#                 # Check if pH control is used with Lime
#                 if (hasattr(r302_1, 'pH_control') and r302_1.pH_control and 
#                     hasattr(r302_1, 'neutralizing_agent') and r302_1.neutralizing_agent in ['Lime', 'CaO']):
#                     # Calculate Lime consumption from succinic acid production
#                     # Formula from unit code: base.imol['Lime'] = mol_lime_per_acid_pH_control * (SuccinicAcid/ X[1]) * safety_factor
#                     if r302_1.outs and len(r302_1.outs) > 0:
#                         effluent = r302_1.outs[0]
#                         if hasattr(effluent, 'imol') and 'SuccinicAcid' in effluent.chemicals:
#                             try:
#                                 succinic_mol = effluent.imol['SuccinicAcid']
#                                 if succinic_mol > 0:
#                                     # Get reaction parameters
#                                     if hasattr(r302_1, 'lime_pH_control_rxns') and hasattr(r302_1.lime_pH_control_rxns, 'X'):
#                                         X_lime = r302_1.lime_pH_control_rxns.X[1] if len(r302_1.lime_pH_control_rxns.X) > 1 else 1.0
#                                     else:
#                                         X_lime = 1.0
                                    
#                                     mol_lime_per_acid = getattr(r302_1, 'mol_lime_per_acid_pH_control', 1.0)
#                                     safety_factor = getattr(r302_1, 'neutralization_safety_factor', 1.0)
                                    
#                                     # Calculate Lime moles consumed
#                                     lime_mol = mol_lime_per_acid * (succinic_mol / X_lime) * safety_factor
                                    
#                                     if lime_mol > 0 and 'Lime' in chems and 'Lime' in price:
#                                         lime_mass = lime_mol * chems.Lime.MW
#                                         lime_cost = price['Lime'] * lime_mass
#                                         if lime_cost > 0.01:
#                                             total_cost += lime_cost
#                             except (KeyError, AttributeError, TypeError):
#                                 pass
#         except Exception:
#             pass  # If R302_1 doesn't exist or can't calculate Lime, continue
        
#         return total_cost
    # conversion_group.metrics[-1].getter = conversion_material_cost_updated

#%% TEA breakdown (old)

def TEA_breakdown(print_output=False, fractions=False):
    metric_breakdowns = {i.name: {} for i in unit_groups[0].metrics}
    for ug in unit_groups:
        for metric in ug.metrics:
            denominator = 1.
            if fractions:
                metric_name = metric.name
                if metric_name=='Inst. eq. cost':
                    denominator = succinic_tea.installed_equipment_cost/1e6
                elif metric_name=='Cooling':
                    denominator=0
                    for i in list(CT.heat_utilities) + list(u_1.CWP802.heat_utilities):
                        if i.flow<0 and i.duty>0:
                            denominator += i.duty/1e6
                elif metric_name=='Heating':
                    for i in list(BT.heat_utilities):
                        if i.flow<0 and i.duty<0:
                            denominator -= i.duty/1e6
                elif metric_name=='Elec. cons.':
                    denominator = succinic_sys_new .power_utility.consumption/1e3
                elif metric_name=='Elec. prod.':
                    denominator = succinic_sys_new .power_utility.production/1e3
                elif metric_name=='Mat. cost':
                    denominator = (succinic_tea.material_cost/succinic_tea.operating_hours) + BT.natural_gas.F_mass*BT.natural_gas_price
            # storage_metric_val = None
            if not ug.name=='storage':
                if ug.name=='other facilities':
                    metric_breakdowns[metric.name]['storage and ' + ug.name] = (metric() + unit_groups_dict['storage'].metrics[ug.metrics.index(metric)]())/denominator
                else:
                    metric_breakdowns[metric.name][ug.name] = metric()/denominator
                    
                    
            if ug.name=='natural gas (for steam generation)':
                if metric.name=='Mat. cost' or metric.name=='Material cost':
                    metric_breakdowns[metric.name][ug.name] = BT.natural_gas.F_mass*BT.natural_gas_price/denominator
            
            
            # else:
            #     storage_metric_val = metric()
                
    # print and return metric_breakdowns
    if print_output:
        for i in unit_groups[0].metrics:
            print(f"\n\n----- {i.name} ({i.units}) -----")
            metric_breakdowns_i = metric_breakdowns[i.name]
            for j in metric_breakdowns_i.keys():
                print(f"{j}: {format(metric_breakdowns_i[j], '.3f')}")
    return metric_breakdowns

TEA_breakdown(print_output=True, fractions=True)# -*- coding: utf-8 -*-


# %% Carbon balance

def get_carbon_flow_by_area():
    c_flow_units = []
    for ug in unit_groups:
        sys_ID = ''.join(e for e in ug.name if e.isalnum())
        ug_sys = ug.to_system(sys_ID)
        ug_ins = [stream for stream in ug_sys.streams if stream.source not in ug_sys.units]
        ug_outs = [stream for stream in ug_sys.streams if stream.sink not in ug_sys.units]
        if sum([stream.get_atomic_flow('C') for stream in ug_ins])>0 or sum([stream.get_atomic_flow('C') for stream in ug_outs])>0:
            c_flow_units.append(ug_sys)
    # print(c_flow_units)
    source=[]
    target=[]
    value=[]
    carbon_feeds = [stream for stream in succinic_sys_new .feeds if stream.get_atomic_flow('C')>0]
    label=[unit.ID  + ' process' for unit in c_flow_units] + [stream.ID for stream in carbon_feeds] + ['SuccinicAcid'] + ['Emissions']
    for unit in c_flow_units:
        ug_outs = [stream for stream in unit.streams if stream.sink not in unit.units]
        # print(ug_outs)
        for stream in ug_outs:
            if stream.get_atomic_flow('C'):
                if stream.sink:
                    for ug in c_flow_units:
                        if stream.sink in ug.units:
                            target.append(c_flow_units.index(ug))
                    source.append(c_flow_units.index(unit))
                    value.append(stream.get_atomic_flow('C'))
                else:
                    if not stream.ID==product_stream.ID:
                        source.append(c_flow_units.index(unit))
                        target.append(len(label)-1)
                        value.append(stream.get_atomic_flow('C'))
                    else:
                        source.append(c_flow_units.index(unit))
                        target.append(len(label)-2)
                        value.append(stream.get_atomic_flow('C'))
    for stream in carbon_feeds:
        source.append(label.index(stream.ID))
        for ug in c_flow_units:
            if stream.sink in ug.units:
                target.append(c_flow_units.index(ug))
        value.append(stream.get_atomic_flow('C'))
    return source, target, value, label


def get_carbon_flow():
    c_flow_units = []
    for unit in u_1:
        if sum([stream.get_atomic_flow('C') for stream in unit.ins])>0 or sum([stream.get_atomic_flow('C') for stream in unit.outs])>0:
            c_flow_units.append(unit)
    source=[]
    target=[]
    value=[]
    carbon_feeds = [stream for stream in succinic_sys_new .feeds if stream.get_atomic_flow('C')>0]
    label=[unit.ID  + ': ' + unit.__class__.__name__ for unit in c_flow_units] + [stream.ID for stream in carbon_feeds] + ['SuccinicAcid'] + ['Emissions']
    for unit in c_flow_units:
        for stream in unit.outs:
            if stream.get_atomic_flow('C'):
                if stream.sink:
                    source.append(c_flow_units.index(unit))
                    target.append(c_flow_units.index(stream.sink))
                    value.append(stream.get_atomic_flow('C'))
                else:
                    if not stream.ID==product_stream.ID:
                        source.append(c_flow_units.index(unit))
                        target.append(len(label)-1)
                        value.append(stream.get_atomic_flow('C'))
                    else:
                        source.append(c_flow_units.index(unit))
                        target.append(len(label)-2)
                        value.append(stream.get_atomic_flow('C'))
    for stream in carbon_feeds:
        source.append(label.index(stream.ID))
        target.append(c_flow_units.index(stream.sink))
        value.append(stream.get_atomic_flow('C'))
    return source, target, value, label

def plot_carbon_flow():
    import plotly.graph_objs as go#create links
    from plotly.offline import init_notebook_mode,  plot
    init_notebook_mode()
    source, target, value, label = get_carbon_flow_by_area()
    link = dict(source=source, target=target, value=value, 
    color=['rgba(0,0,0, 0.5)'] * len(source))#create nodes
    node = dict(label=label, pad=15, thickness=8, 
                color = ['blue' if ':' in label[i] else 'gray' for i in range(len(label))],
               )#create a sankey object
    chart = go.Sankey(link=link, node=node, arrangement="freeform")#build a figure
    fig = go.Figure(chart)
    plot(fig)
    fig.show()
    fig.write_image("fig1.png")



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
    
#     #
# # After creating the plot and before saving or showing it
# #                    Adjust y-axis tick font size
    

