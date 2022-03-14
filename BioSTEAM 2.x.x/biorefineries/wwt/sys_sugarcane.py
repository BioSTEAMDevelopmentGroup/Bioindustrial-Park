#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the sugarcane biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/sugarcane
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


import biosteam  as bst
from biosteam import main_flowsheet

# from biorefineries.wwt import (
#     sc
#     create_sc_chemicals,
#     new_price, load_cs_settings,
#     create_wastewater_treatment_system
#     )
# from biorefineries.utils import get_MESP
from __init__ import sc
from _chemicals import create_sc_chemicals
from _settings import new_price, load_sc_settings
from _wwt_sys import create_wastewater_system
from utils import print_MESP


# %%

# =============================================================================
# Function to make the system
# =============================================================================

load_sc_settings()
chems = create_sc_chemicals()
bst.settings.set_thermo(chems)

@bst.SystemFactory(
    ID='sugarcane_sys',
    ins=[*sc.create_juicing_system_with_fiber_screener.ins,
         sc.create_ethanol_purification_system.ins[1]], # denaturant
    outs=[sc.create_ethanol_purification_system.outs[0], # ethanol
          dict(ID='emissions'),
          dict(ID='ash_disposal')]

)
def create_sc_system(ins, outs, **wwt_kwargs):
    s = main_flowsheet.stream
    u = main_flowsheet.unit

    sugarcane, enzyme, H3PO4, lime, polymer, denaturant = ins
    ethanol, emissions, ash_disposal = outs

    feedstock_handling_sys = sc.create_feedstock_handling_system(
        ins=[sugarcane],
        outs=[''],
        mockup=True,
    )

    juicing_sys = sc.create_juicing_system_with_fiber_screener(
        ins=[feedstock_handling_sys-0, enzyme, H3PO4, lime, polymer],
        mockup=True
    )

    ethanol_production_sys = sc.create_sucrose_to_ethanol_system(
        ins=(juicing_sys-0, denaturant), outs=(ethanol, 'vinasse'),
        mockup=True
    )

    M305 = bst.units.Mixer(
        'M305',
        ins=(juicing_sys-2, *ethanol_production_sys-[2, 3]),
        outs='wastewater'
    )

    ### Wastewater treatment ###
    create_wastewater_system(
        ins=[ethanol_production_sys-1, M305-0],
        outs=['biogas', 'sludge_S603', 'recycled_water', 'brine'],
        mockup=True,
        R601_kwargs={'method': 'lumped'},
        **wwt_kwargs,
    )

    ### Facilities ###
    bst.units.BoilerTurbogenerator('BT',
        (juicing_sys-1, u.M602-0, 'boiler_makeup_water', 'natural_gas', '', ''),
        outs=(emissions, 'rejected_water_and_blowdown', ash_disposal),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85)
    bst.units.CoolingTower('CT')
    bst.units.ChilledWaterPackage('CWP')

    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    process_water_streams = (s.imbibition_water,
                             s.rvf_wash_water,
                             s.stripping_water,
                             *makeup_water_streams)
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    bst.units.ProcessWaterCenter('PWC',
                                 (u.S604-0, # recycled wastewater from reverse osmosis
                                 makeup_water),
                                 (),
                                 None,
                                 makeup_water_streams,
                                 process_water_streams)

    plant_air = bst.Stream('plant_air', N2=83333, units='kg/hr')
    ADP = bst.facilities.AirDistributionPackage('ADP', plant_air)
    @ADP.add_specification(run=True)
    def adjust_plant_air():
        plant_air.imass['N2'] = 0.8 * feedstock_handling_sys.ins[0].F_mass

    F301 = u.F301
    D303 = u.D303
    bst.HeatExchangerNetwork('HXN', units=[F301, D303])


# %%

# =============================================================================
# Create the system
# =============================================================================

flowsheet = bst.Flowsheet('wwt_sugarcane')
main_flowsheet.set_flowsheet(flowsheet)
sugarcane_sys = create_sc_system(skip_IC=True, skip_AnMBR=True)
u = main_flowsheet.unit
s = main_flowsheet.stream

sugarcane_sys.simulate()
sugarcane_tea = sc.create_tea(sugarcane_sys)

ethanol = s.ethanol

# Compare MESP
original_IRR = 0.1267
sugarcane_tea.IRR = sc.sugarcane_tea.IRR = original_IRR
print(f'\n\nIRR = {original_IRR:.0%}')
sc.wastewater.price = 0.
MESP_old = print_MESP(sc.ethanol, sc.sugarcane_tea, 'old sc sys w/o ww cost')
sc.wastewater.price = new_price['Wastewater']
MESP_old = print_MESP(sc.ethanol, sc.sugarcane_tea, 'old sc sys w ww cost')
MESP_new = print_MESP(ethanol, sugarcane_tea, 'new sc sys')

sugarcane_tea.IRR = sc.sugarcane_tea.IRR = 0.1
assert(sugarcane_tea.IRR==sc.sugarcane_tea.IRR)
print(f'\n\nIRR = {sugarcane_tea.IRR:.0%}')
sc.wastewater.price = 0.
MESP_old = print_MESP(sc.ethanol, sc.sugarcane_tea, 'old sc sys w/o ww cost')
sc.wastewater.price = new_price['Wastewater']
MESP_old = print_MESP(sc.ethanol, sc.sugarcane_tea, 'old sc sys w ww cost')
MESP_new = print_MESP(ethanol, sugarcane_tea, 'new sc sys')


# %%

# =============================================================================
# Sum up results
# =============================================================================

wwt_units = main_flowsheet.wastewater_treatment_system.units
new_capexes = {i.ID: i.installed_cost/1e6 for i in wwt_units}
new_capex = sum(i for i in new_capexes.values())

new_powers_wwt = {i.ID: i.power_utility.rate/1e3 for i in wwt_units}
new_power_wwt = sum(i for i in new_powers_wwt.values())/1e3
new_power_tot = sum(i.power_utility.consumption for i in sugarcane_sys.units)
new_power_ratio = new_power_wwt / new_power_tot

new_power_net = sum(i.power_utility.rate for i in sugarcane_sys.units)/1e3


# Hf in kJ/hr
net_e = s.ethanol.Hf/3600/1e3 + u.BT.power_utility.rate/1e3
net_e_ratio = net_e/(s.sugarcane.Hf/3600/1e3)

# # A lot (similar to ethanol) goes to filter_cake
# s.filter_cake.Hf/3600/1e3

# Water
water_usage = u.PWC.F_mass_in/s.ethanol.F_mass # 11 kg/kg
water_consumption = (u.PWC.ins[1].F_mass-u.PWC.outs[1].F_mass) / \
    s.ethanol.F_mass # about -5 kg/kg since sugarcane has 70% of water

# Wastewater
from utils import compute_stream_COD
wastewater = u.M601.F_mass_in/s.ethanol.F_mass # about 12 kg/kg
COD = compute_stream_COD(u.M601.outs[0]) # about 5.8 g COD/L

# Ratios for IC design
# from _utils import get_CN_ratio
# R601_CN = get_CN_ratio(u.R601._inf)