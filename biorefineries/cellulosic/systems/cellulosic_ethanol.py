# -*- coding: utf-8 -*-
"""
.. autofunction:: biorefineries.cellulosic.systems.cellulosic_ethanol.create_cellulosic_ethanol_system

"""

import biosteam as bst
import thermosteam as tmo
from thermosteam import Stream
from .pretreatment import create_dilute_acid_pretreatment_system
from .fermentation import create_cellulosic_fermentation_system
from biorefineries.ethanol import create_ethanol_purification_system
from biorefineries.cellulosic import streams as s
from biorefineries.cellulosic import units


__all__ = (
    'create_cellulosic_ethanol_system',
    'create_facilities',
)

def create_facilities(
        solids_to_boiler, 
        gas_to_boiler,
        process_water_streams,
        feedstock,
        RO_water='',
        recycle_process_water='',
        blowdown_to_wastewater=None,
        BT_area=None,
        area=None,
        BTkw=None
    ):
    
    BT = bst.facilities.BoilerTurbogenerator(BT_area or area or 'BT',
                                             **(BTkw or {}),
                                             ins=(solids_to_boiler,
                                                  gas_to_boiler, 
                                                  'boiler_makeup_water',
                                                  'natural_gas',
                                                  'FGD_lime',
                                                  'boilerchems'))
    
    bst.facilities.ChilledWaterPackage(area or 'CWP')
    CT = bst.facilities.CoolingTower(area or 'CT')
    
    process_water_streams = (*process_water_streams,
                              BT-1, CT-1)
            
    makeup_water = Stream('makeup_water', price=0.00025353160151260924)
    
    bst.facilities.ProcessWaterCenter(area or 'PWC',
        (RO_water, makeup_water, recycle_process_water),
        (),
        None,
        (BT-1, CT-1),
        process_water_streams
    )
    CIP = Stream('CIP', Water=126, units='kg/hr')
    CIP_package = bst.facilities.CIPpackage(area or 'CIP_package', CIP)
    CIP_package.CIP_over_feedstock = 0.00121
    @CIP_package.add_specification(run=True)
    def adjust_CIP(): CIP.imass['Water'] = feedstock.F_mass * CIP_package.CIP_over_feedstock
    
    plant_air = Stream('plant_air', N2=83333, units='kg/hr')
    ADP = bst.facilities.AirDistributionPackage(area or 'ADP', plant_air)
    ADP.plant_air_over_feedstock = 0.8
    @ADP.add_specification(run=True)
    def adjust_plant_air(): plant_air.imass['N2'] = feedstock.F_mass * ADP.plant_air_over_feedstock
        
    fire_water = Stream('fire_water', Water=8343, units='kg/hr')
    FT = bst.FireWaterTank(area or 'FT', fire_water)
    FT.fire_water_over_feedstock = 0.08
    @FT.add_specification(run=True)
    def adjust_fire_water(): fire_water.imass['Water'] = feedstock.F_mass * FT.fire_water_over_feedstock
    
    ### Complete system
    if blowdown_to_wastewater:
        bst.BlowdownMixer(area or 'blowdown_mixer', (BT-1, CT-1), blowdown_to_wastewater)

@bst.SystemFactory(
    ID='cornstover_sys',
    ins=[*create_dilute_acid_pretreatment_system.ins,
          s.denaturant],
    outs=[s.ethanol],
)
def create_cellulosic_ethanol_system(
        ins, outs,
        include_blowdown_recycle=None,
        WWT_kwargs=None,
    ):
    feedstock, sulfuric_acid, ammonia, denaturant = ins
    ethanol, = outs
    U101 = units.FeedStockHandling('U101', feedstock)
    U101.cost_items['System'].cost = 0.
    pretreatment_sys = create_dilute_acid_pretreatment_system(
        ins=[U101-0, sulfuric_acid, ammonia],
        mockup=True
    )
    fermentation_sys = create_cellulosic_fermentation_system(
        ins=pretreatment_sys-0,
        mockup=True,
    )
    ethanol_purification_sys, udct = create_ethanol_purification_system(
        ins=[fermentation_sys-1, denaturant],
        outs=[ethanol],
        udct=True,
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
    udct['H401'].dT = 10
    udct['D402'].k = 1.4
    udct['D403'].k = 1.4
    ethanol, stillage, recycle_process_water = ethanol_purification_sys.outs
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    S401 = bst.PressureFilter('S401', (stillage, recycled_water))

    bst.create_all_facilities(
        feedstock, 
        blowdown_recycle=include_blowdown_recycle,
        WWT_kwargs=WWT_kwargs,
        HXN=False,
        recycle_process_water_streams=[recycle_process_water],
    )
