#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 11:05:24 2017

@author: Yoel
"""
import biosteam as bst
from biosteam import units
from ._process_settings import price

__all__ = ('create_ethanol_production_system',
           'mass2molar_ethanol_fraction')

# %% Helpful functions

def mass2molar_ethanol_fraction(x):
    """Return ethanol mol fraction in a ethanol water mixture"""
    return x/46.06844 / (x/46.06844 + (1-x)/18.01528)

# %% Pretreatment section

def create_ethanol_production_system(ID='ethanol_production_sys',
                                     sugar_solution=None):
    ### Streams ###
    
    # Fresh water
    stripping_water = bst.Stream('stripping_water', Water=5000, units='kg/hr')
    
    # Gasoline
    denaturant = bst.Stream('denaturant', Octane=230.69,
                            units='kg/hr', price=price['Gasoline'])
    
    if not sugar_solution: # Feedstock
        sugar_solution = bst.Stream('sugar_solution',
            Glucose = 3802,
            Sucrose = 4.309e+04,
            Water   = 2.59e+05,
            H3PO4   = 83.33,
            units = 'kg/hr',
            T = 372,
        )
    
    # Yeast
    yeast = bst.Stream('yeast', Water=24700, DryYeast=10300, units='kg/hr')
    
    # Ethanol product
    ethanol = bst.Stream('ethanol', price=price['Ethanol'])
    
    ### Units ###
    
    # Split sugar solution
    S301 = units.Splitter('S301',
                        split=0.265)
    
    # Concentrate sugars
    F301 = units.MultiEffectEvaporator('F301',
                                       P=(101325, 73581, 50892, 32777),
                                       V=0.95) # fraction evaporated
    F301.components['condenser'].U = 1.85
    # Note: value of steam ~ 6.86 for the following 
    # (101325, 73580.467, 50891.17, 32777.406, 19999.925, 11331.5),
    
    # Mix sugar solutions
    M301 = units.Mixer('M301')
    
    # Cool for fermentation
    H301 = units.HXutility('H301', T=295.15)
    
    # Ethanol Production
    R301 = units.Fermentation('R301', outs=('CO2', ''), tau=9, efficiency=0.90, N=4) 
    T301 = units.StorageTank('T301', tau=4, vessel_material='Carbon steel')
    T301.line = 'Beer tank'
    
    D301 = units.VentScrubber('D301', ins=(stripping_water, R301-0), gas=('CO2',))
    
    # Separate 99% of yeast
    C301 = units.SolidsCentrifuge('C301', outs=('', 'recycle_yeast'),
                                split=(1, 0.99999, 1, 0.96, 0.01),
                                order=('Ethanol', 'Glucose', 'H3PO4', 
                                       'Water', 'DryYeast'),
                                solids=('DryYeast',))
    
    # Mix in Water
    M302 = units.Mixer('M302')
    P301 = units.Pump('P301')
    
    # Heat up before beer column
    # Exchange heat with stillage
    H302 = units.HXprocess('H302', outs=('', 'stillage'),
                          phase0='l', phase1='l', U=1.28)
    
    # Beer column
    xbot = mass2molar_ethanol_fraction(0.00001)
    ytop = mass2molar_ethanol_fraction(0.574)
    D302 = units.BinaryDistillation('D302', P=101325,
                                y_top=ytop, x_bot=xbot, k=1.25,
                                LHK=('Ethanol', 'Water'))
    D302.tray_material = 'Stainless steel 304'
    D302.vessel_material = 'Stainless steel 304'
    D302.boiler.U = 1.85
    P302 = units.Pump('P302')
    
    # Mix ethanol Recycle (Set-up)
    M303 = units.Mixer('M303')
    
    ytop = mass2molar_ethanol_fraction(0.9061726)
    D303 = units.BinaryDistillation('D303', P=101325,
                                y_top=ytop, x_bot=xbot, k=1.25,
                                LHK=('Ethanol', 'Water'),
                                tray_material='Stainless steel 304',
                                vessel_material='Stainless steel 304',
                                is_divided=True)
    D303.boiler.U = 1.85
    P303 = units.Pump('P303')
    
    # Superheat vapor for mol sieve
    H303 = units.HXutility('H303', T=115+273.15, V=1)
    
    # Molecular sieve
    U301 = units.MolecularSieve('U301',
                                split=(2165.14/13356.04, 1280.06/1383.85),
                                order=('Ethanol', 'Water'))
    
    # Condense ethanol product
    H304 = units.HXutility('H304', 'S149', V=0, T=340.)
    T302 = units.StorageTank('T302', tau=7*24,
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel')
    P304 = units.Pump('P304')
    
    # Storage for gasoline
    T303 = units.StorageTank('T303', tau=7*24,
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel')
    P305 = units.Pump('P305')
    
    # Denatured ethanol product
    T304 = units.MixTank('T304', outs=ethanol)
    T304.tau = 0.10
    
    # Recycle water to Process Condensate Tank
    M305 = units.Mixer('M305', outs='recycle_water')
    
    # Yeast mixing
    T305 = units.MixTank('T305')
    T305.tau = 0.1
    yeast-T305
    
    # Multi-effect evaporator pumps
    P306 = units.Pump('P306')
    
    
    ### Ethanol system set-up ###
    
    sugar_solution-S301-1-F301-0-P306
    (S301-0, P306-0)-M301-H301
    (H301-0, yeast-T305-0)-R301-1-T301-0-C301
    (C301-0, D301-1)-M302-P301
    (P301-0, P302-0)-H302-0-D302-1-P302
    (D302-0, U301-0)-M303-0-D303-0-H303-U301
    D303-1-P303
    
    pure_ethanol = P304.outs[0]
    def adjust_denaturant():
        denaturant.imol['Octane'] = 0.021*pure_ethanol.F_mass/114.232
        
    PS4 = bst.ProcessSpecification('PS4', specification=adjust_denaturant)
        
    U301-1-H304-0-T302-0-P304-0-PS4
    denaturant-T303-P305
    (P305-0, PS4-0)-T304
    (P303-0, F301-1)-M305
    
    ### System ###
    
    return bst.System(ID, 
                [S301, 
                 F301, 
                 P306, 
                 M301, 
                 H301, 
                 T305,
                 R301,
                 T301, 
                 C301, 
                 M302, 
                 P301,
                 bst.System('beer_column_heat_integration',
                     [H302,
                      D302,
                      P302],
                     recycle=P302-0),
                 bst.System('ethanol_recycle_from_molecular_sieves',
                     [M303,
                      D303,
                      H303,
                      U301],
                     recycle=U301-0),
                 H304,
                 T302, 
                 P304,
                 PS4, 
                 T303, 
                 P305, 
                 T304, 
                 D301,
                 P303, 
                 M305])
