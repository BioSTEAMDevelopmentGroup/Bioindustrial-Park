#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 11:05:24 2017

@author: Yoel
"""
import numpy as np
import biosteam as bst
from biosteam import units
from biorefineries.lipidcane.tea import LipidcaneTEA
from biorefineries.sugarcane.chemicals import sugarcane_chemicals
from biorefineries.sugarcane.process_settings import price

__all__ = ('sugarcane_sys', 'sugarcane_tea', 'sugarcane', 'sugar_cane')

# %% Pretreatment section

bst.main_flowsheet.set_flowsheet('sugarcane')
bst.settings.set_thermo(sugarcane_chemicals)

### Streams ###

F_mass_sugarcane = 333334
z_mass_sugarcane = sugarcane_chemicals.kwarray(
    dict(Glucose=0.0120811,
         Lignin=0.0327653,
         Solids=0.015,
         Sucrose=0.136919,
         Ash=0.006,
         Cellulose=0.0611531,
         Hemicellulose=0.036082,
         Water=0.7)
)

sugarcane = sugar_cane = bst.Stream('sugar_cane',
                                    flow=F_mass_sugarcane * z_mass_sugarcane,
                                    units='kg/hr',
                                    price=price['Sugar cane'])

enzyme = bst.Stream('enzyme',
                    Cellulose=100, Water=900, units='kg/hr',
                    price=price['Protease'])

imbibition_water = bst.Stream('imbibition_water',
                              Water=87023.35, units='kg/hr',
                              T = 338.15)

H3PO4 = bst.Stream('H3PO4',
                   H3PO4=74.23, Water=13.10, units='kg/hr',
                   price=price['H3PO4'])  # to T203

lime = bst.Stream('lime',
                  CaO=333.00, Water=2200.00, units='kg/hr',
                  price=price['Lime'])  # to P5

polymer = bst.Stream('polymer',
                     Flocculant=0.83, units='kg/hr',
                     price=price['Polymer'])  # to T205

rvf_wash_water = bst.Stream('rvf_wash_water',
                            Water=16770, units='kg/hr',
                            T=363.15)  # to C202

oil_wash_water = bst.Stream('oil_wash_water',
                            Water=1350, units='kg/hr',
                            T=358.15)  # to T207

### Unit operations ###

bst.Stream.ticket_name = 'd'
bst.Stream.ticket_number = 100

# Feed the shredder
U101 = units.ConveyingBelt('U101', ins=sugar_cane)

# Separate metals
U102 = units.MagneticSeparator('U102', ins=U101-0)

# Shredded cane
U103 = units.Shredder('U103', ins=U102-0)

bst.Stream.ticket_number = 200

# Hydrolyze starch
T201 = units.EnzymeTreatment('T201', T=323.15)  # T=50

# Finely crush lipid cane
U201 = units.CrushingMill('U201',
                          split=dict(Ash=0.92,
                                     Cellulose=0.92,
                                     Glucose=0.04,
                                     Hemicellulose=0.92,
                                     Lignin=0.92,
                                     Sucrose=0.04,
                                     Solids=1),
                          moisture_content=0.5)

# Convey out bagasse
U202 = units.ConveyingBelt('U202', ins=U201.outs[0], outs='Bagasse')

# Mix in water
M201 = units.Mixer('M201')

# Screen out fibers
S201 = units.VibratingScreen('S201',
                             split=dict(Ash=0.35,
                                        Cellulose=0.35,
                                        Glucose=0.88,
                                        Hemicellulose=0.35,
                                        Lignin=0.35,
                                        Solids=0,
                                        Sucrose=0.88,
                                        Water=0.88))

# Store juice before treatment
T202 = units.StorageTank('T202', tau=4, vessel_material='Carbon steel')

# Heat up before adding acid
H201 = units.HXutility('H201', T=343.15)

# Mix in acid
T203 = units.MixTank('T203')

# Pump acid solution
P201 = units.Pump('P201')

# Mix lime solution
T204 = units.MixTank('T204', tau=0.25)
P202 = units.Pump('P202')

# Blend acid lipid solution with lime
T205 = units.MixTank('T205')

# Mix recycle
M202 = units.Mixer('M202')

# Heat before adding flocculant
H202 = units.HXutility('H202', T=372.15)

# Mix in flocculant
T206 = units.MixTank('T206')
T206.tau = 0.10

# Separate residual solids
C201 = units.Clarifier('C201',
                       split=dict(Ash=0,
                                  CaO=0,
                                  Cellulose=0,
                                  Flocculant=0.522,
                                  Glucose=0.522,
                                  Hemicellulose=0,
                                  Lignin=0,
                                  H3PO4=0.522,
                                  Sucrose=0.522,
                                  Water=0.522))

# Remove solids as filter cake
C202 = units.RVF('C202', 
                 outs=('filte_cake', ''),
                 moisture_content=0.80,
                 split=dict(Ash=0.85,
                            CaO=0.85,
                            Cellulose=0.85,
                            Glucose=0.01,
                            Hemicellulose=0.85,
                            Lignin=0.85,
                            Sucrose=0.01))
P203 = units.Pump('P203')


# Screen out small fibers from sugar stream
S202 = units.VibratingScreen('S202', outs=('', 'fiber_fines'),
                             split=dict(Ash=1.0,
                                        CaO=1.0,
                                        Cellulose=1.0,
                                        Flocculant=0.0,
                                        Glucose=0.998,
                                        Hemicellulose=1.0,
                                        Lignin=1.0,
                                        H3PO4=1.0,
                                        Sucrose=0.998,
                                        Water=0.998))
sugar = S202-0
S202.mesh_opening = 2

### Process specifications ###

# Specifications dependent on lipid cane flow rate
F_mass_last_sugarcane = int(sugar_cane.F_mass)
def correct_flows():
    global F_mass_last_sugarcane
    F_mass = sugar_cane.F_mass
    if int(F_mass) != F_mass_last_sugarcane:
        # correct enzyme, lime, phosphoric acid, and imbibition water
        enzyme.imass['Cellulose', 'Water'] = 0.003 * F_mass * np.array([0.1, 0.9])
        lime.imass['CaO', 'Water'] = 0.001 * F_mass * np.array([0.046, 0.954])
        H3PO4.imass['H3PO4', 'Water'] = 0.00025 * F_mass
        imbibition_water.imass['Water'] = 0.25* F_mass
        F_mass_last_sugarcane = int(F_mass)

PS2 = bst.ProcessSpecification('PS2', specification=correct_flows)

# Specifications within a system
def correct_wash_water():
    solids = P202.outs[0].imol['Ash', 'CaO', 'Cellulose',
                               'Hemicellulose', 'Lignin'].sum()
    rvf_wash_water.imol['Water'] = 0.0574 * solids

PS3 = bst.ProcessSpecification('PS3', specification=correct_wash_water)


### System set-up ###

U103-0-PS2
(PS2-0, enzyme)-T201
(T201-0, M201-0)-U201-1-S201-0-T202
(S201-1, imbibition_water)-M201
crushing_mill_recycle_sys = bst.System('crushing_mill_recycle_sys',
                               path=(U201, S201, M201),
                               recycle=M201-0)

T202-0-H201
(H201-0, H3PO4)-T203-P201
(P201-0, lime-T204-0)-T205-P202
P202-0-PS3
(PS3-0, P203-0)-M202-H202
(H202-0, polymer)-T206-C201
(C201-1, rvf_wash_water)-C202-1-P203
C201-0-S202



# %% Ethanol section

### Utilities ###

MW_etoh = sugarcane_chemicals.Ethanol.MW
MW_water = sugarcane_chemicals.Water.MW

def mass2molar_ethanol_fraction(ethanol_mass_fraction):
    """Return ethanol mol fraction in a ethanol water mixture"""
    x = ethanol_mass_fraction
    return x/MW_etoh / (x/MW_etoh + (1-x)/MW_water)

### Input streams ###

# Fresh water
stripping_water = bst.Stream('stripping_water', Water=5000, units='kg/hr')

# Gasoline
denaturant = bst.Stream('denaturant', Octane=230.69,
                        units='kg/hr', price=price['Gasoline'])

# Yeast
yeast = bst.Stream('yeast', Water=24700, DryYeast=10300, units='kg/hr')

# From Pretreatment section
sugar_solution = S202-0

# Ethanol product
ethanol = bst.Stream('ethanol', price=price['Ethanol'])

### Units ###

# Split sugar solution
S301 = units.Splitter('S301',
                    split=0.265)

# Concentrate sugars
F301 = units.MultiEffectEvaporator('F301',
                                   P=(101325, 73581, 50892, 32777, 20000),
                                   V=0.95248) # fraction evaporated
F301.components['condenser'].U = 1.85
# Note: value of steam ~ 6.86 for the following 
# (101325, 73580.467, 50891.17, 32777.406, 19999.925, 11331.5),

# Mix sugar solutions
M301 = units.Mixer('M301')

# Cool for fermentation
H301 = units.HXutility('H301', T=295.15)

# Ethanol Production
R301 = units.Fermentation('R301', outs=('CO2', ''), tau=10, efficiency=0.90, N=6) 
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
                     fluid_type='ss', U=1.28)

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

# denaturantd ethanol product
T304 = units.MixTank('T304', outs=ethanol)
T304.tau = 0.10

# Add denaturant
M304 = units.Mixer('M304')

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


# %% Facilities

emission = bst.Stream('emission')
stream = bst.main_flowsheet.stream

# Stream.default_ID_number = 500

BT = units.BoilerTurbogenerator('BT',
                                ins=U202-0, # Bagasse from conveyor belt
                                outs=emission,
                                boiler_efficiency=0.80,
                                turbogenerator_efficiency=0.85)

bst.Stream.ticket_number = 600

CT = units.CoolingTower('CT')
makeup_water_streams = (stream.cooling_tower_makeup_water,
                        stream.boiler_makeup_water)
process_water_streams = (stream.imbibition_water,
                         *makeup_water_streams)
makeup_water = bst.Stream('makeup_water', price=0.000254)

CWP = units.ChilledWaterPackage('CWP')
PWC = units.ProcessWaterCenter('PWC',
                               (M305-0, makeup_water),
                               (),
                               makeup_water_streams,
                               process_water_streams)

# %% Perform TEA

sugarcane_sys = bst.main_flowsheet.create_system('sugarcane_sys')
sugarcane_tea = LipidcaneTEA(system=sugarcane_sys, IRR=0.15,
                             duration=(2018, 2038),
                             depreciation='MACRS7', income_tax=0.35,
                             operating_days=200, lang_factor=3,
                             construction_schedule=(0.4, 0.6), WC_over_FCI=0.05,
                             labor_cost=2.5e6, fringe_benefits=0.4,
                             property_tax=0.001, property_insurance=0.005,
                             supplies=0.20, maintenance=0.01, administration=0.005)
sugarcane_sys.simulate()
sugarcane_tea.IRR = sugarcane_tea.solve_IRR()
