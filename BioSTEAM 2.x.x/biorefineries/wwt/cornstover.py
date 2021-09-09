#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# Part of this module is based on the cornstover biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/cornstover
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


import biosteam as bst
from biosteam import main_flowsheet as F
from biorefineries import (
    cornstover as cs,
    sugarcane as sc
    )

# from biorefineries.wwt import (
#     create_cs_chemicals,
#     cs_price, load_cs_settings,
#     create_wastewater_treatment_system,
#     get_cs_GWP,
#     )
# from biorefineries.utils import get_MESP
from _chemicals import create_cs_chemicals
from _settings import cs_price, load_cs_settings
from _wwt_sys import create_wastewater_treatment_system
from _lca import get_cs_GWP
from utils import get_MESP


# %%

# =============================================================================
# Function to make the system
# =============================================================================

load_cs_settings()
chems = create_cs_chemicals()
bst.settings.set_thermo(chems)

@bst.SystemFactory(
    ID='cornstover_sys',
    ins=[*cs.create_dilute_acid_pretreatment_system.ins,
          dict(ID='denaturant',
               Octane=1,
               price=cs_price['Denaturant'])],
    outs=[dict(ID='ethanol', price=cs_price['Ethanol'])],
)
def create_cs_system(ins, outs, include_blowdown_recycle=True, **wwt_kwargs):
    feedstock, denaturant = ins
    ethanol, = outs
    f = bst.main_flowsheet
    s = f.stream
    u = f.unit
    U101 = cs.FeedStockHandling('U101', feedstock)
    U101.cost_items['System'].cost = 0.

    pretreatment_sys = cs.create_dilute_acid_pretreatment_system(
        ins=U101-0,
        mockup=True
    )

    fermentation_sys = cs.create_cellulosic_fermentation_system(
        ins=pretreatment_sys-0,
        mockup=True,
    )

    ethanol_purification_sys = sc.create_ethanol_purification_system(
        ins=[fermentation_sys-1, denaturant],
        outs=[ethanol],
        IDs={'Beer pump': 'P401',
             'Beer column heat exchange': 'H401',
             'Beer column': 'D402',
             'Beer column bottoms product pump': 'P402',
             'Distillation': 'D403',
             'Distillation bottoms product pump': 'P403',
             'Ethanol-denaturant mixer': 'M701',
             'Recycle mixer': 'M402',
             'Heat exchanger to superheat vapor to molecular sieves': 'H402',
             'Molecular sieves': 'U401',
             'Ethanol condenser': 'H403',
             'Ethanol day tank': 'T701',
             'Ethanol day tank pump': 'P701',
             'Denaturant storage': 'T702',
             'Denaturant pump': 'P702',
             'Product tank': 'T703'},
        mockup=True,
    )
    ethanol, stillage, stripper_bottoms_product = ethanol_purification_sys.outs
    recycled_water = bst.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    S401 = bst.PressureFilter('S401', (stillage, recycled_water))
    if include_blowdown_recycle:
        blowdown_to_wastewater = bst.Stream('blowdown_to_wastewater')
    else:
        blowdown_to_wastewater = None

    create_wastewater_treatment_system(
        ins=[S401-1, pretreatment_sys-1, blowdown_to_wastewater],
        outs=['biogas', 'sludge_S603', 'recycled_water', 'brine'],
        mockup=True,
        R601_kwargs={'method': 'lumped'},
        R602_kwargs={'HRT': 35},
        **wwt_kwargs,
    )

    M501 = bst.Mixer('M501', (u.S603-1, S401-0))

    cs.create_facilities(
        solids_to_boiler=M501-0,
        gas_to_boiler=s.biogas,
        process_water_streams=(s.stripping_water,
                               s.warm_process_water_1,
                               s.warm_process_water_2,
                               s.pretreatment_steam,
                               s.saccharification_water),
        feedstock=feedstock,
        RO_water=u.S604-0,
        recycle_process_water=stripper_bottoms_product,
        blowdown_to_wastewater=blowdown_to_wastewater,
    )


# %%

# =============================================================================
# Create the system
# =============================================================================

flowsheet = bst.Flowsheet('wwt_cornstover')
F.set_flowsheet(flowsheet)
cornstover_sys = create_cs_system(include_blowdown_recycle=True)

u = F.unit
wwt_units = [i for i in u if i.ID[1:3]=='60']
OSBL_units = (*wwt_units, u.CWP, u.CT, u.PWC, u.ADP,
              u.T701, u.T702, u.P701, u.P702, u.M701, u.FT,
              u.CSL_storage, u.DAP_storage, u.BT)
cornstover_tea = cs.create_tea(cornstover_sys, OSBL_units)
ethanol = F.stream.ethanol

# Compare MESP
assert(cornstover_tea.IRR==cs.cornstover_tea.IRR==0.1)
print(f'\n\nIRR = {cornstover_tea.IRR:.0%}')
MESP_old = get_MESP(cs.ethanol, cs.cornstover_tea, 'old cs sys')
MESP_new = get_MESP(ethanol, cornstover_tea, 'new cs sys')


# %%

# =============================================================================
# Util functions for result comparison
# =============================================================================

# Old system
cs_wwt_units = [i for i in cs.cornstover_sys.units if i.ID[1:3]=='60']
old_capex = cs.WWTC.installed_cost / 1e6
old_capex_ratio = old_capex/(cs.cornstover_tea.installed_equipment_cost/1e6)

old_power_wwt = cs.WWTC.power_utility.rate
old_power_tot = sum(i.power_utility.consumption for i in cs.cornstover_sys.units)
old_power_ratio = old_power_wwt / old_power_tot

old_power_net = sum(i.power_utility.rate for i in cs.cornstover_sys.units)

# New system
new_capexes = {i.ID: i.installed_cost/1e6 for i in wwt_units}
new_capex = sum(i for i in new_capexes.values())
new_capex_ratio = new_capex/(cornstover_tea.installed_equipment_cost/1e6)

new_powers_wwt = {i.ID: i.power_utility.rate/1e3 for i in wwt_units}
new_power_wwt = sum(i for i in new_powers_wwt.values())
new_power_tot = sum(i.power_utility.consumption for i in cornstover_sys.units)/1e3
new_power_ratio = new_power_wwt / new_power_tot

new_power_net = sum(i.power_utility.rate for i in cornstover_sys.units)/1e3

s = F.stream
# Hf in kJ/hr
net_e = s.ethanol.Hf/3600/1e3 + u.BT.power_utility.rate/1e3
net_e_ratio = net_e/(s.cornstover.Hf/3600/1e3)

# Water
# water_usage = (u.PWC.F_mass_in+s.cornstover.imass['Water'])/s.ethanol.F_mass # 18 kg/kg
water_usage = u.PWC.F_mass_in/s.ethanol.F_mass # 17 kg/kg
# water_consumption = (u.PWC.ins[1].F_mass-u.PWC.outs[1].F_mass+s.cornstover.imass['Water']) / \
#     s.ethanol.F_mass # about 1.7 kg/kg
water_consumption = (u.PWC.ins[1].F_mass-u.PWC.outs[1].F_mass) / \
    s.ethanol.F_mass # about 0.8 kg/kg

# Wastewater
from utils import compute_stream_COD
wastewater = u.M601.F_mass_in/s.ethanol.F_mass # about 18 kg/kg
COD = compute_stream_COD(u.M601.outs[0]) # about 54 g COD/L

# OLR
from utils import auom
_ft_to_m3 = auom('ft3').conversion_factor('m3')
old_OLR_R601 = compute_stream_COD(cs.R601.ins[0])*cs.R601.ins[0].F_vol / \
    (4*31*1e6*auom('gal').conversion_factor('m3')) * 24 # 1.4 g COD/L/d
old_OLR_R602 = compute_stream_COD(cs.R602.ins[0])*cs.R602.ins[0].F_vol / \
    (3*25*115*344*_ft_to_m3) * 24 # 0.7 g COD/L/d

new_OLR_R601 = compute_stream_COD(u.R601.ins[0])*u.R601.ins[0].F_vol / \
    (u.R601.Vliq) * 24 # 30 g COD/L/d
V_R602 = u.R602.D_tank*u.R602.W_tank*u.R602.L_CSTR*u.R602.N_train * \
    _ft_to_m3
new_OLR_R602 = compute_stream_COD(u.R602.ins[0])*u.R602.ins[0].F_vol / \
    V_R602 * 24 # 10.3 g COD/L/d
# new_OLR_R603 = compute_stream_COD(u.R603._inf)*u.R603._inf.F_vol / \
#     (u.R603.design_results['Volume [ft3]']*_ft_to_m3) * 24 # 2.25 g COD/L/d
new_OLR_R603 = u.R603.OLR # 2.25 g COD/L/d


# Adjust methane production
def adjust_methane():
    from utils import get_digestable_chemicals
    flowsheet_ch4 = bst.Flowsheet('wwt_cornstover_ch4')
    F.set_flowsheet(flowsheet_ch4)
    cornstover_sys_ch4 = create_cs_system(include_blowdown_recycle=True,
                                          skip_R603=True)
    u_ch4 = F.unit
    wwt_units = [i for i in u if i.ID[1:3]=='60']
    OSBL_units = (*wwt_units, u_ch4.CWP, u_ch4.CT, u_ch4.PWC, u_ch4.ADP,
                  u_ch4.T701, u_ch4.T702, u_ch4.P701, u_ch4.P702, u_ch4.M701, u_ch4.FT,
                  u_ch4.CSL_storage, u_ch4.DAP_storage, u_ch4.BT)

    cornstover_tea_ch4 = cs.create_tea(cornstover_sys_ch4, OSBL_units)

    BD_dct = {k.ID: 1. for k in get_digestable_chemicals(chems)}
    u_ch4.R601.biodegradability = BD_dct
    u_ch4.R601._refresh_rxns(X_biogas=0.86, X_growth=0.05)
    # About 5481, vs. 5681 from Humbird, about 3-4% less,
    # might be due to the different CH4/CO2 ratios in the biogas production reaction,
    # Humbird assumed 51%:49% (1.04) CH4:CO2 on a molar basis, here is about 1
    # print(u.R601.outs[0].imass['CH4'])
    return flowsheet_ch4, cornstover_sys_ch4, cornstover_tea_ch4

flowsheet_ch4, cornstover_sys_ch4, cornstover_tea_ch4 = adjust_methane()

get_MESP(F.stream.ethanol, cornstover_tea_ch4, 'new cs sys with higher BD')


# %%

get_ratio = lambda cornstover, ethanol: \
    ethanol.F_mass/cs.ethanol_density_kggal/(cornstover.F_mass-cornstover.imass['Water'])

lca_streams_cs = [
    cs.denaturant,
    cs.cellulase,
    cs.sulfuric_acid,
    cs.DAP,
    cs.CSL,
    cs.ammonia,
    cs.FGD_lime,
    cs.caustic,
    cs.emissions,
    ]

# One stream that representing all chemicals
lca_stream_cs = bst.Stream('lca_stream_cs')
lca_stream_cs.mix_from(lca_streams_cs)

GWP_cs = get_cs_GWP(lca_stream_cs, cs.flowsheet,
                    get_ratio(cs.cornstover, cs.ethanol))

# Get the gross results
print(f'Total GWP for original cs sys is {GWP_cs:.2f} kg CO2-eq/gal ethanol.')

s_ch4 = flowsheet_ch4.stream
lca_streams_ch4 = [
    s_ch4.denaturant,
    s_ch4.cellulase,
    s_ch4.sulfuric_acid,
    s_ch4.DAP,
    s_ch4.CSL,
    s_ch4.ammonia,
    s_ch4.FGD_lime,
    s_ch4.emissions,
    ]

# One stream that representing all chemicals
lca_stream_ch4 = bst.Stream('lca_stream_ch4')
lca_stream_ch4.mix_from(lca_streams_ch4)

GWP_ch4 = get_cs_GWP(lca_stream_ch4, flowsheet_ch4,
                    get_ratio(s_ch4.cornstover, s_ch4.ethanol))

# Get the gross results
print(f'Total GWP for new cs sys with higher BD is {GWP_ch4:.3f} kg CO2-eq/gal ethanol.')




# %%

# # C/N ratio
# from _utils import get_CN_ratio
# R601_CN = get_CN_ratio(u.R601._inf)