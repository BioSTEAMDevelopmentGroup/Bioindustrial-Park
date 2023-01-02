#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from biosteam import SystemFactory, F
from ... import streams as s
from biorefineries.ethanol import (
    create_ethanol_purification_system,
)
from ..fermentation import create_sucrose_fermentation_system
from ..juicing import create_feedstock_handling_system, create_juicing_system

__all__ = (
    'create_sucrose_to_ethanol_system',
    'create_sugarcane_to_ethanol_system',
)

# %% Sugarcane to ethanol

@SystemFactory(
    ID='sucrose_to_ethanol_sys',
    ins=[s.screened_juice, s.denaturant],
    outs=[s.ethanol, s.stillage, s.recycle_process_water, s.evaporator_condensate]
)
def create_sucrose_to_ethanol_system(ins, outs, add_urea=False):
    screened_juice, denaturant = ins
    ethanol, stillage, recycle_process_water, evaporator_condensate = outs
    
    beer = bst.Stream()
    
    create_sucrose_fermentation_system(
        ins=screened_juice,
        outs=[beer, evaporator_condensate],
        mockup=True,
        add_urea=add_urea,
    )
    create_ethanol_purification_system(
        ins=[beer, denaturant], 
        outs=[ethanol, stillage, recycle_process_water],
        mockup=True
    )

@SystemFactory(
    ID='sugarcane_sys', 
    ins=[s.sugarcane, s.H3PO4, s.lime, s.polymer, s.denaturant], 
    outs=[s.ethanol, s.vinasse, s.wastewater, s.emissions, s.ash_disposal]
)
def create_sugarcane_to_ethanol_system(ins, outs, 
                                       use_area_convention=False,
                                       pellet_bagasse=None,
                                       dry_bagasse=None):
    sugarcane, H3PO4, lime, polymer, denaturant = ins
    ethanol, vinasse, wastewater, emissions, ash_disposal = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        area=100 if use_area_convention else None,
        ins=[sugarcane],
        outs=[''],
        mockup=True,
    )
    juicing_sys = create_juicing_system(
        area=200 if use_area_convention else None,
        ins=[feedstock_handling_sys-0, H3PO4, lime, polymer],
        pellet_bagasse=pellet_bagasse,
        dry_bagasse=dry_bagasse,
        mockup=True
    )
    ethanol_production_sys, edct = create_sucrose_to_ethanol_system(
        area=300 if use_area_convention else None,
        udct=True,
        ins=(juicing_sys-0, denaturant), outs=(ethanol, vinasse),
        mockup=True
    )
    M305 = bst.Mixer(400 if use_area_convention else 'M305', 
        ins=(juicing_sys-2, *ethanol_production_sys-[2, 3]),
        outs=wastewater
    )
    
    ### Facilities ###    
    
    BT = bst.BoilerTurbogenerator(400 if use_area_convention else 'BT',
        (juicing_sys-1, '', 'boiler_makeup_water', 'natural_gas', '', ''),
        outs=(emissions, 'rejected_water_and_blowdown', ash_disposal),
        boiler_efficiency=0.80,
        turbogenerator_efficiency=0.85
    )
    CT = bst.CoolingTower(500 if use_area_convention else 'CT')
    makeup_water_streams = (F.cooling_tower_makeup_water,
                            F.boiler_makeup_water)
    process_water_streams = (F.imbibition_water,
                             F.rvf_wash_water,
                             F.stripping_water,
                             *makeup_water_streams)
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    CWP = bst.ChilledWaterPackage(500 if use_area_convention else 'CWP')
    PWC = bst.ProcessWaterCenter(500 if use_area_convention else 'PWC',
                                   (bst.Stream(), makeup_water),
                                   (),
                                   None,
                                   makeup_water_streams,
                                   process_water_streams)
    
    F301 = edct['F301']
    D303 = edct['D303']
    HXN = bst.HeatExchangerNetwork(600 if use_area_convention else 'HXN',
                                   units=[F301, D303.condenser])
    
    # if vinasse_to_wastewater:
    #     plant_air = bst.Stream('plant_air', N2=83333, units='kg/hr')
    #     ADP = bst.facilities.AirDistributionPackage('ADP', plant_air)
    #     @ADP.add_specification(run=True)
    #     def adjust_plant_air():
    #         plant_air.imass['N2'] = 0.8 * feedstock_handling_sys.ins[0].F_mass
            
    #     wastewater_treatment_sys = bst.create_wastewater_treatment_system(
    #         ins=[vinasse],
    #         mockup=True,
    #     )