# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:09:15 2022

@author: sarangbhagwat
"""
import numpy as np
import biosteam as bst
import thermosteam as tmo
from biosteam import SystemFactory
# from biorefineries.TAL.system_TAL_adsorption_glucose import create_TAL_sys as create_sys_TAL_adsorption_glucose
from biorefineries.TAL.system_TAL_adsorption_glucose import TAL_sys as sys_TAL
from biorefineries.TAL import units
from biorefineries.TAL.chemicals_data import TAL_chemicals
from biorefineries.TAL.process_settings import price, CFs
from biorefineries.cornstover import CellulosicEthanolTEA as TALTEA

from warnings import filterwarnings
filterwarnings('ignore')

Stream = tmo.Stream
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# from lactic.hx_network import HX_Network

# # Do this to be able to show more streams in a diagram
# bst.units.Mixer._graphics.edge_in *= 2
bst.speed_up()
flowsheet = bst.Flowsheet('TAL')
bst.main_flowsheet.set_flowsheet(flowsheet)

# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0

# Baseline cost year is 2016
bst.CE = 541.7
# _labor_2007to2016 = 22.71 / 19.55

# Set default thermo object for the system
tmo.settings.set_thermo(TAL_chemicals)


@SystemFactory(ID = 'sorbic_acid_sys')
def create_sys_sorbic_acid_adsorption_glucose(ins, outs):
    # sys_TAL = create_sys_TAL_adsorption_glucose()
    units_TAL = sys_TAL.flowsheet.unit
    streams_TAL = sys_TAL.flowsheet.stream
    
    Hydrogen = Stream('Hydrogen', units = 'kg/hr')
    KOH = Stream('KOH', units = 'kg/hr')
    HCl = Stream('HCl', units = 'kg/hr')
    
    F401, M401, M501 = units_TAL['F401'], units_TAL['M401'], units_TAL['M501']
    excluded_units = [units_TAL['F402'], units_TAL['H403'], units_TAL['T620'], units_TAL['T620_P']]
    
    R401 = units.HydrogenationReactor('R401', ins = (F401-1, '', Hydrogen), outs = 'HMTHP',
                                      vessel_material='Stainless steel 316',)
    
    
    R402 = units.DehydrationRingOpeningReactor('R402', ins = (R401-0, ''), outs = 'SA',
                                               vessel_material='Stainless steel 316',
                                               tau = 12)
    
    R403 = units.HydrolysisReactor('R403', ins = (R402-0, '', KOH), outs = 'KSA',
                                               vessel_material='Stainless steel 316')
    
    splits_S406 = np.zeros(len(TAL_chemicals))
    
    splits_S406[TAL_chemicals.index('KSA')] = 0.98
    splits_S406[TAL_chemicals.index('Water')] = 0.2
    
    S406 = bst.units.SolidsCentrifuge('S406', ins=R403-0, outs=('K_sorbate', ''),
                                # moisture_content=0.50,
                                split=splits_S406,
                                solids = ['KSA'])
    # !!! splits for moisture content (negligible in feed), hexanol content 
    def S406_spec():
        try:
            S406._run()
        except:
            moisture_content = S406.moisture_content
            # S406.ins[0].imol['Water'] = 0.
            S406.moisture_content/= 5.
            S406._run()
            S406.moisture_content = moisture_content
    
    S406.specification = S406_spec
    
    F401 = bst.Flash('F401', ins=S406-1, outs=('F401_hexanol_to_recycle', 'F401_organics_to_boiler'),
                      V = 0.9999, P=101325.) # need to eventually change this to a rotary evaporator
 
    Ethanol_Tb = TAL_chemicals['Ethanol'].Tb
    def F401_spec():
        F401_V = 0.
        F401_ins_0 = F401.ins[0]
        F401_ins_0_F_mol = F401_ins_0.F_mol
        for chem_i in F401_ins_0.vle_chemicals:
            if chem_i.Tb<=Ethanol_Tb:
                F401_V += F401_ins_0.imol[chem_i.ID]/F401_ins_0_F_mol
        F401.V = min(1.2*F401_V, 1-1e-6)
        F401._run()
        F401.outs[1].imol['TAL'] += F401.outs[0].imol['TAL']
        F401.outs[0].imol['TAL'] = 0.
        
    F401.specification = F401_spec
    
    F401-1-1-R401
    
    F401_H = bst.HXutility('F401_H', ins=F401-0, T = 313.15, rigorous=True,)
    

    F401_P = bst.Pump('F401_P', ins=F401_H-0, P=101325.)
    X401 = bst.Splitter('X401', ins=F401_P-0, outs=('recycled_F401_0', 'wasted_F401_0'), split = 0.9999)
    
    # X401-0-1-R402 # recycle HMTHP to R402
    X401-0-2-M401 # recycle Hexanol to M401

    S407 = units.Crystallization('S407', ins = (S406-0, HCl, '', ''), outs = ('wet_SorbicAcid_crystals', 'KCl'))
    
    def S407_spec():
        S407._run()
        S408.run()
    S407.specification = S407_spec
    
    R404 = units.HClKOHRecovery('R404', ins = (S407-1, 'water'),
                                outs = ('HCl_recycle', 'KOH_recycle'))
    
    
    R404-0-2-S407
    R404-1-1-R403
    
    S408 = units.Decantation('S408', ins=S407-0,
                             outs=('dissolved_SorbicAcid_to_evaporative_crystallization_or_WWT', 'SorbicAcid_crystals'),
                             forced_recovery=None)
    
    #%% Wastewater treatment and CHP
    M501.ins = (
        units_TAL.AC401-0,
        units_TAL.S402-1,
                X401-1, 
                S408-0,
                )
    # M505.ins = (
    #             S503-1,    
    #             S401-0, 
    #             )
    #%% Storage
    hydrogen_fresh = Stream('hydrogen_fresh', price=price['Hydrogen'])
    KOH_fresh = Stream('KOH_fresh', price=price['KOH'])
    HCl_fresh = Stream('HCl_fresh', price=price['HCl'])
    SA = Stream('SA', units='kg/hr', price=price['SA'])
    
    
    T607 = units.DPHPStorageTank('T607', ins=hydrogen_fresh, outs = Hydrogen)
    T607.line = 'Hydrogen storage tank'
    
    T608 = units.DPHPStorageTank('T608', ins=HCl_fresh, outs = HCl,
                                 vessel_material = 'Stainless steel')
    T608.line = 'HCl storage tank'
    
    T609 = units.DPHPStorageTank('T609', ins=KOH_fresh, outs = KOH,
                                 vessel_material = 'Stainless steel')
    T609.line = 'KOH storage tank'
    
    
    T621 = units.TALStorageTank('T621', ins=S408-1, tau=7*24, V_wf=0.9,
                                          vessel_type='Floating roof',
                                          vessel_material='Stainless steel')
    
    T621.line = 'SAStorageTank'
    T621_P = units.TALPump('T621_P', ins=T621-0, outs=SA)


#%% Create system
SA_sys = create_sys_sorbic_acid_adsorption_glucose()
u, s = SA_sys.flowsheet.unit, SA_sys.flowsheet.stream
feedstock = s.feedstock

get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
SA = s['SA']

TAL_tea = TALTEA(system=SA_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(u.U101, u.WWTcost501,
                    # u.T601, u.T602, 
                    u.T603, u.T604, u.T620,
                    # u.T606, u.T606_P,
                    u.CWP802, u.CT801, u.PWC904, u.CIP901, u.ADP902, u.FWT903, u.BT701),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT701)


def get_SA_MPSP():
    for i in range(3):
        SA_sys.simulate()
    for i in range(3):
        SA.price = TAL_tea.solve_price(SA)
    return SA.price*SA.F_mass/SA.imass['TAL']

def simulate_and_print():
    get_SA_MPSP()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${get_SA_MPSP():.3f}/kg')


simulate_and_print()