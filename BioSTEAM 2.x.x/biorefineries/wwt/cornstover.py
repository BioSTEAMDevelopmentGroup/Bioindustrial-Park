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

# from biorefineries.wwt import (
#     sc, cs,
#     create_cs_chemicals,
#     cs_price, load_cs_settings,
#     create_wastewater_treatment_system,
#     get_cs_GWP,
#     )
# from biorefineries.utils import get_MESP, get_digestable_chemicals
from __init__ import sc, cs, ethanol_density_kggal
from _chemicals import create_cs_chemicals
from _settings import cs_price, load_cs_settings
from _wwt_sys import create_wastewater_treatment_system
from _lca import get_cs_GWP
from utils import print_MESP, get_digestable_chemicals


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
def create_cs_system(ins, outs, include_blowdown_recycle=True,
                     default_BD=True, wwt_kwargs={}):
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

    if default_BD:
        skip_R603 = wwt_kwargs.get('skip_R603')
        wwt_kwargs['skip_R603'] = True if skip_R603 is None else skip_R603

    create_wastewater_treatment_system(
        ins=[S401-1, pretreatment_sys-1, blowdown_to_wastewater],
        outs=['biogas', 'sludge_S603', 'recycled_water', 'brine'],
        mockup=True,
        R601_kwargs={'method': 'lumped'},
        R602_kwargs={'HRT': 35},
        **wwt_kwargs,
    )

    # Using the default organic biodegradability as in Humbird et al.
    # Methane production is  5481, vs. 5681 from Humbird (u.R601.outs[0].imass['CH4']),
    # which is about 3-4% less,
    # This might be due to the different CH4/CO2 ratios in the biogas production reaction,
    # Humbird assumed 51%:49% (1.04) CH4:CO2 on a molar basis, here is about 1
    if default_BD:
        BD_dct = {k.ID: 1. for k in get_digestable_chemicals(chems)}
        u.R601.biodegradability = BD_dct
        u.R601._refresh_rxns(X_biogas=0.86, X_growth=0.05)

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
# Create the system, do TEA and LCA
# =============================================================================

flowsheet = bst.Flowsheet('cs_wwt')
F.set_flowsheet(flowsheet)
sys_wwt = create_cs_system(include_blowdown_recycle=True, default_BD=False)

u = F.unit
wwt_units = [i for i in u if i.ID[1:3]=='60']
OSBL_units = (*wwt_units, u.CWP, u.CT, u.PWC, u.ADP,
              u.T701, u.T702, u.P701, u.P702, u.M701, u.FT,
              u.CSL_storage, u.DAP_storage, u.BT)
tea_wwt = cs.create_tea(sys_wwt, OSBL_units)

s = F.stream
ethanol = s.ethanol

# Get the ratio of ethanol-to-dry-cornstover
get_ratio = lambda cornstover, ethanol: \
    ethanol.F_mass/ethanol_density_kggal/(cornstover.F_mass-cornstover.imass['Water'])

# Input chemicals/generated emissions with impacts
lca_streams_cs = [
    cs.CSL,
    cs.DAP,
    cs.FGD_lime,
    cs.ammonia,
    cs.caustic,
    cs.cellulase,
    cs.denaturant,
    cs.sulfuric_acid,
    cs.emissions,
    ]

lca_streams_wwt = [
    s.CSL,
    s.DAP,
    s.FGD_lime,
    s.ammonia,
    s.bisulfite_R602,
    s.cellulase,
    s.citric_R602,
    s.denaturant,
    s.naocl_R602,
    s.sulfuric_acid,
    s.emissions,
    ]

# Use one stream to represent all chemicals/emissions
lca_stream_cs = bst.Stream('lca_stream_cs')
lca_stream_wwt = bst.Stream('lca_stream_wwt')


# Util functions
def get_GWP(system):
    f = system.flowsheet
    u = f.unit
    s = f.stream
    if hasattr(u, 'WWTC'):
        streams = lca_streams_cs
        stream = lca_stream_cs
    else:
        streams = lca_streams_wwt
        stream = lca_stream_wwt
    stream.mix_from(streams)
    return get_cs_GWP(stream, f, get_ratio(s.cornstover, s.ethanol))


def get_WWT_CAPEX(system, ratio=False):
    u = system.flowsheet.unit
    if hasattr(u, 'WWTC'): # Humbird et al.
        capex = u.WWTC.installed_cost
    else:
        wwt_units = [i for i in u if i.ID[1:3]=='60']
        capex = sum(i.installed_cost for i in wwt_units)
    if not ratio:
        return capex/1e6 # in MM$
    return capex/system.TEA.installed_equipment_cost

def get_WWT_power(system, ratio=False):
    u = system.flowsheet.unit
    if hasattr(u, 'WWTC'): # Humbird et al.
        power = cs.WWTC.power_utility.rate
    else:
        wwt_units = [i for i in u if i.ID[1:3]=='60']
        power = sum(i.power_utility.rate for i in wwt_units)
    if not ratio:
        return power # in kWh
    tot = sum(i.power_utility.rate for i in u)
    return power/tot


# %%

if __name__ == '__main__':
    # Print MESP and GWP
    print(f'\n\nIRR = {tea_wwt.IRR:.0%}')
    print_MESP(cs.ethanol, cs.cornstover_tea, 'old cs sys')
    print_MESP(ethanol, tea_wwt, 'new cs sys')
    print(f'Total GWP for the original cs sys is {get_GWP(cs.cornstover_sys):.2f} kg CO2-eq/gal ethanol.')
    print(f'Total GWP for new cs sys is {get_GWP(sys_wwt):.2f} kg CO2-eq/gal ethanol.')


# %%

# =============================================================================
# Legacy codes
# =============================================================================

# from utils import compute_stream_COD

# # C/N ratio
# from _utils import get_CN_ratio
# R601_CN = get_CN_ratio(u.R601._inf)


# # Water usage
# # water_usage = (u.PWC.F_mass_in+s.cornstover.imass['Water'])/s.ethanol.F_mass # 18 kg/kg
# water_usage = u.PWC.F_mass_in/s.ethanol.F_mass # 17 kg/kg
# # water_consumption = (u.PWC.ins[1].F_mass-u.PWC.outs[1].F_mass+s.cornstover.imass['Water']) / \
# #     s.ethanol.F_mass # about 1.7 kg/kg
# water_consumption = (u.PWC.ins[1].F_mass-u.PWC.outs[1].F_mass) / \
#     s.ethanol.F_mass # about 0.8 kg/kg


# # Wastewater generated
# wastewater = u.M601.F_mass_in/s.ethanol.F_mass # about 18 kg/kg
# COD = compute_stream_COD(u.M601.outs[0]) # about 54 g COD/L


# # OLR
# from utils import auom
# _ft3_to_m3 = auom('ft3').conversion_factor('m3')
# old_OLR_R601 = compute_stream_COD(cs.R601.ins[0])*cs.R601.ins[0].F_vol / \
#     (4*31*1e6*auom('gal').conversion_factor('m3')) * 24 # 1.4 g COD/L/d
# old_OLR_R602 = compute_stream_COD(cs.R602.ins[0])*cs.R602.ins[0].F_vol / \
#     (3*25*115*344*_ft3_to_m3) * 24 # 0.7 g COD/L/d

# new_OLR_R601 = compute_stream_COD(u.R601.ins[0])*u.R601.ins[0].F_vol / \
#     (u.R601.Vliq) * 24 # 30 g COD/L/d
# V_R602 = u.R602.D_tank*u.R602.W_tank*u.R602.L_CSTR*u.R602.N_train * \
#     _ft3_to_m3
# new_OLR_R602 = compute_stream_COD(u.R602.ins[0])*u.R602.ins[0].F_vol / \
#     V_R602 * 24 # 10.3 g COD/L/d
# # new_OLR_R603 = compute_stream_COD(u.R603._inf)*u.R603._inf.F_vol / \
# #     (u.R603.design_results['Volume [ft3]']*_ft3_to_m3) * 24 # 2.25 g COD/L/d
# new_OLR_R603 = u.R603.OLR # 2.25 g COD/L/d


# # Energy recovery
# # Hf in kJ/hr
# net_e = s.ethanol.Hf/3600/1e3 + u.BT.power_utility.rate/1e3
# net_e_ratio = net_e/(s.cornstover.Hf/3600/1e3)