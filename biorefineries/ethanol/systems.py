# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from biosteam import SystemFactory
from . import streams as s

__all__ = (
    'create_beer_distillation_system',
    'create_ethanol_purification_system_after_beer_column',
    'create_ethanol_purification_system',
)

@SystemFactory(
    ID='beer_distillation_sys',
    ins=[s.beer],
    outs=[s.distilled_beer, s.stillage] # TODO: Find good selling price for stillage/vinasse and possibly yeast
)
def create_beer_distillation_system(ins, outs,
                                    beer_column_heat_integration=True,
                                    IDs={}):
    beer, = ins
    distilled_beer, stillage = outs
    
    P301 = bst.Pump(IDs.get('Beer pump', 'P301'), P=2. * 101325)
    
    # Beer column
    x_bot = 3.91e-06
    y_top = 0.28
    D302 = bst.BinaryDistillation(IDs.get('Beer column', 'D302'), P=2. * 101325,
                                    outs=(distilled_beer, ''),
                                y_top=y_top, x_bot=x_bot, k=1.1, Rmin=0.001,
                                LHK=('Ethanol', 'Water'))
    D302.tray_material = 'Stainless steel 304'
    D302.vessel_material = 'Stainless steel 304'
    D302.boiler.U = 1.85
    P302 = bst.Pump(IDs.get('Beer column bottoms product pump', 'P302'), P=101325.)
    
    # Heat up before beer column
    # Exchange heat with stillage    
    if beer_column_heat_integration:
        H302 = bst.HXprocess(IDs.get('Beer column heat exchange', 'H302'), 
                               outs=('', stillage),
                               phase0='l', phase1='l', U=1.28)
        (beer-P301-0, P302-0)-H302-0-D302-1-P302
    else:
        beer-P301-0-D302-1-P302
        P302.outs[0] = stillage
        
@SystemFactory(
    ID='ethanol_purification_from_distilled_beer_sys',
    ins=[s.distilled_beer, s.denaturant],
    outs=[s.ethanol, s.recycle_process_water]
)
def create_ethanol_purification_system_after_beer_column(ins, outs, IDs={}):
    distilled_beer, denaturant = ins
    ethanol, recycle_process_water = outs
  
    
    # Mix ethanol Recycle (Set-up)
    M303 = bst.Mixer(IDs.get('Recycle mixer', 'M303'))
    
    D303 = bst.BinaryDistillation(IDs.get('Distillation', 'D303'), 
                                x_bot=3.9106e-06, y_top=0.80805, k=1.25, Rmin=0.01,
                                LHK=('Ethanol', 'Water'),
                                tray_material='Stainless steel 304',
                                vessel_material='Stainless steel 304',
                                P=1013250.,
                                is_divided=True)
    D303.boiler.U = 1.85
    P303 = bst.Pump(IDs.get('Distillation bottoms product pump', 'P303'), 
                      outs=recycle_process_water)
    
    # Superheat vapor for mol sieve
    H303 = bst.HXutility(IDs.get('Heat exchanger to superheat vapor to molecular sieves', 'H303'),
                           T=115+273.15, V=1, heat_only=True)
    
    # Molecular sieve
    U301 = bst.MolecularSieve(IDs.get('Molecular sieves', 'U301'),
                                split=(2165.14/13356.04, 1280.06/1383.85),
                                order=('Ethanol', 'Water'))
    
    # Condense ethanol product
    H304 = bst.HXutility(IDs.get('Ethanol condenser', 'H304'), 'S149', V=0, T=340.)
    T302 = bst.StorageTank(IDs.get('Ethanol day tank', 'T302'), tau=12,
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel')
    P304 = bst.Pump(IDs.get('Ethanol day tank pump', 'P304'))
    
    # Storage for gasoline
    T303 = bst.StorageTank(IDs.get('Denaturant storage', 'T303'), tau=7*24,
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel')
    P305 = bst.Pump(IDs.get('Denaturant pump', 'P305'))
    
    # Denatured ethanol product
    M304 = bst.Mixer(IDs.get('Ethanol-denaturant mixer', 'M304'))
    T304 = bst.StorageTank(IDs.get('Product tank', 'T304'),
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel',
                             tau=6.5*24, outs=ethanol)
    
    ### Ethanol system set-up ###
    (distilled_beer, U301-0)-M303-0-D303-0-H303-U301
    D303-1-P303
    M304.denaturant_fraction = 0.022
    
    @M304.add_specification(run=True, impacted_units=[T303])
    def adjust_denaturant():
        pure_ethanol = M304.ins[1]
        denaturant.imol['Octane'] = M304.denaturant_fraction * pure_ethanol.F_mass / 114.232
    
    (denaturant-T303-P305-0, U301-1-H304-0-T302-0-P304-0)-M304-T304

        
@SystemFactory(
    ID='ethanol_purification_sys',
    ins=[s.beer, s.denaturant],
    outs=[s.ethanol, s.stillage, s.recycle_process_water]
)
def create_ethanol_purification_system(ins, outs,
                                       beer_column_heat_integration=True,
                                       IDs={}):
    beer, denaturant = ins
    ethanol, stillage, recycle_process_water = outs
    distilled_beer = bst.Stream('')
    create_beer_distillation_system(
        ins=beer,
        outs=[distilled_beer, stillage],
        beer_column_heat_integration=beer_column_heat_integration,
        IDs=IDs,
        mockup=True,
    )
    create_ethanol_purification_system_after_beer_column(
        ins=[distilled_beer, denaturant],
        outs=[ethanol, recycle_process_water],
        IDs=IDs,
        mockup=True,
    )
    