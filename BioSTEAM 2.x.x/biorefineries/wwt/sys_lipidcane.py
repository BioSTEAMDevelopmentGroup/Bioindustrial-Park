#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the lipidcane biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/lipidcane
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam  as bst
from biosteam import main_flowsheet
from biorefineries import lipidcane as lc
from biorefineries.sugarcane import create_sucrose_to_ethanol_system
from biorefineries.lipidcane import create_chemicals, load_process_settings, price
from biorefineries.wwt import (
    add_wwt_chemicals, create_wastewater_system,
    print_MESP, new_price,
    )


# %%

# =============================================================================
# Function to make the system
# =============================================================================

new_lc_chems = add_wwt_chemicals(create_chemicals())
bst.settings.set_thermo(new_lc_chems)
load_process_settings()

@bst.SystemFactory(
    ID='lipidcane_sys',
    ins=[*lc.create_juicing_and_lipid_extraction_system.ins,
         create_sucrose_to_ethanol_system.ins[1]],
    outs=[dict(ID='ethanol', price=price['Ethanol']),
          dict(ID='biodiesel', price=price['Biodiesel']),
          dict(ID='crude_glycerol', price=price['Crude glycerol']),
          dict(ID='emissions'),
          dict(ID='ash_disposal')]
)
def create_lc_system(ins, outs, **wwt_kwargs):
    s = main_flowsheet.stream
    u = main_flowsheet.unit

    lipidcane, enzyme, H3PO4, lime, polymer, denaturant = ins
    ethanol, biodiesel, crude_glycerol, emissions, ash_disposal = outs

    feedstock_handling_sys = lc.create_feedstock_handling_system(
        ins=lipidcane,
        mockup=True,
    )

    ### Oil and juice separation ###
    juicing_and_lipid_extraction_sys = lc.create_juicing_and_lipid_extraction_system(
        ins=[feedstock_handling_sys-0, enzyme, H3PO4, lime, polymer],
        mockup=True,
    )

    ### Ethanol section ###
    ethanol_production_sys = create_sucrose_to_ethanol_system(
        ins=[juicing_and_lipid_extraction_sys-0, denaturant],
        outs=[ethanol, 'vinasse'],
        mockup=True,
    )

    ### Biodiesel section ###
    # Fresh degummed oil
    oil = juicing_and_lipid_extraction_sys-1
    lc.create_transesterification_and_biodiesel_separation_system(
        ins=oil,
        outs=[biodiesel, crude_glycerol],
        mockup=True,
    )

    ### Wastewater treatment ###
    M305 = bst.Mixer(
        'M305',
        ins=(juicing_and_lipid_extraction_sys-4,
             juicing_and_lipid_extraction_sys-3,
             *ethanol_production_sys-[2, 3]),
        outs='wastewater',
    )

    create_wastewater_system(
        ins=[ethanol_production_sys-1, M305-0],
        outs=['biogas', 'sludge_S603', 'recycled_water', 'brine'],
        mockup=True,
        R601_kwargs={'method': 'lumped'},
        **wwt_kwargs,
    )

    ### Facilities ###
    bst.units.BoilerTurbogenerator('BT',
                                   (juicing_and_lipid_extraction_sys-2, u.M602-0,
                                    'boiler_makeup_water', 'natural_gas', '', ''),
                                   (emissions, 'rejected_water_and_blowdown', ash_disposal),
                                   boiler_efficiency=0.80,
                                   turbogenerator_efficiency=0.85)
    bst.units.CoolingTower('CT')
    bst.units.ChilledWaterPackage('CWP')

    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    process_water_streams = (s.imbibition_water,
                             s.biodiesel_wash_water,
                             s.oil_wash_water,
                             s.rvf_wash_water,
                             s.stripping_water,
                             s.lipid_wash_water,
                             s.dilution_water,
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

flowsheet = bst.Flowsheet('wwt_lipidcane')
main_flowsheet.set_flowsheet(flowsheet)
lipidcane_sys = create_lc_system(skip_IC=True, skip_AnMBR=True)
u = main_flowsheet.unit
s = main_flowsheet.stream

lipidcane_sys.simulate()
lipidcane_tea = lc.create_tea(lipidcane_sys)
ethanol = s.ethanol

# Compare MESP
original_IRR = 0.2108
lipidcane_tea.IRR = lc.lipidcane_tea.IRR = original_IRR
print(f'\n\nIRR = {original_IRR:.0%}')
lc.wastewater.price = 0.
MESP_old = print_MESP(lc.ethanol, lc.lipidcane_tea, 'old lc sys w/o ww cost')
lc.wastewater.price = new_price['Wastewater']
MESP_old = print_MESP(lc.ethanol, lc.lipidcane_tea, 'old lc sys w ww cost')
MESP_new = print_MESP(ethanol, lipidcane_tea, 'new lc sys')

lipidcane_tea.IRR = lc.lipidcane_tea.IRR = 0.1
ethanol.price = lc.ethanol.price = price['Ethanol']
assert(lipidcane_tea.IRR==lc.lipidcane_tea.IRR)
print(f'\n\nIRR = {lipidcane_tea.IRR:.0%}')
lc.wastewater.price = 0.
MESP_old = print_MESP(lc.ethanol, lc.lipidcane_tea, 'old lc sys w/o ww cost')
lc.wastewater.price = new_price['Wastewater']
MESP_old = print_MESP(lc.ethanol, lc.lipidcane_tea, 'old lc sys w ww cost')
MESP_new = print_MESP(ethanol, lipidcane_tea, 'new lc sys')


# %%

# =============================================================================
# Sum up results
# =============================================================================

wwt_units = main_flowsheet.wastewater_treatment_system.units
new_capexes = {i.ID: i.installed_cost/1e6 for i in wwt_units}
new_capex = sum(i for i in new_capexes.values())

new_powers_wwt = {i.ID: i.power_utility.rate/1e3 for i in wwt_units}
new_power_wwt = sum(i for i in new_powers_wwt.values())
new_power_tot = sum(i.power_utility.consumption for i in lipidcane_sys.units)/1e3
new_power_ratio = new_power_wwt / new_power_tot

new_power_net = sum(i.power_utility.rate for i in lipidcane_sys.units)/1e3

# Hf in kJ/hr
net_e = (s.biodiesel.Hf+s.ethanol.Hf)/3600/1e3 + u.BT.power_utility.rate/1e3
net_e_ratio = net_e/(s.lipidcane.Hf/3600/1e3)

# # A lot (twice as ethanol) goes to filter_cake
# s.filter_cake.Hf/3600/1e3

# For energy-based allocation between ethanol and biodiesel
ratio = s.ethanol.Hf/(s.ethanol.Hf+s.biodiesel.Hf)
# Water
water_usage = u.PWC.F_mass_in*ratio/s.ethanol.F_mass # 19 kg/kg
water_consumption = (u.PWC.ins[1].F_mass-u.PWC.outs[1].F_mass) * \
    ratio/s.ethanol.F_mass # about -8.8 kg/kg since lipidcane has 70% of water

# Wastewater
from utils import compute_stream_COD
wastewater = u.M601.F_mass_in/s.ethanol.F_mass # about 26 kg/kg
COD = compute_stream_COD(u.M601.outs[0]) # only about 3.6 g COD/L

# Ratios for IC design
# from _utils import get_CN_ratio
# R601_CN = get_CN_ratio(u.R601._inf)