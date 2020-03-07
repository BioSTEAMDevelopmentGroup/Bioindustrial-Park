# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 23:03:19 2019

This module defines the composition_balance function, which performs mass energy balance the given oil composition of lipid cane to determine its composition.

@author: yoelr
"""
import biosteam as bst
import thermosteam as tmo
from thermosteam import Stream
from biorefineries.lipidcane.chemicals import pretreatment_chemicals
from biorefineries.lipidcane.system import lipid_cane
from array_collections import tuple_array

__all__ = ('set_lipid_fraction',)

tmo.settings.set_thermo(pretreatment_chemicals)
getattr_ = getattr
carbs_IDs = ('Glucose', 'Sucrose')
fiber_IDs = ('Lignin', 'Cellulose', 'Hemicellulose')
lipid_IDs = ('Lipid',)
water_IDs = ('Water',)

lc = lipid_cane
carbs   = Stream('Carbs')
fiber   = Stream('Fiber')
lipid   = Stream('Lipid')
carbs.imol[carbs_IDs] = lc.imol[carbs_IDs]
fiber.imol[fiber_IDs] = lc.imol[fiber_IDs]
lipid.imol[lipid_IDs] = lc.imol[lipid_IDs]
streams = (carbs, fiber, lipid)

# Mass property arrays
carbs_mass = lc.imass[carbs_IDs]
fiber_mass = lc.imass[fiber_IDs]
lipid_mass = lc.imass[lipid_IDs]

# Net weight
carbs_massnet = carbs_mass.sum()
fiber_massnet = fiber_mass.sum()
lipid_massnet = lipid_mass.sum()

# Heats of combustion per kg
LHV_carbs_kg = carbs.LHV/carbs_massnet
LHV_lipid_kg = lipid.LHV/lipid_massnet

# Composition
carbs_massfrac = tuple_array(carbs_mass/carbs_massnet)
fiber_massfrac = tuple_array(fiber_mass/fiber_massnet)

water_mass = lc.imass['Water']
solids_mass = lc.imass['Solids']
ash_mass = lc.imass['Ash']

def set_lipid_fraction(lipid_fraction):
    """Adjust composition of lipid cane to achieve desired oil fraction (dry weight)."""
    z_mass_lipid = lipid_fraction
    F_mass = lc.F_mass
    z_dry = 0.3
    z_mass_new_lipid = z_mass_lipid * z_dry 
    z_mass_new_carbs = 0.149 - z_mass_new_lipid * (LHV_lipid_kg / LHV_carbs_kg)
    new_fiber_massfrac = z_dry - z_mass_new_carbs - z_mass_new_lipid - 0.006 - 0.015
    
    water_mass.value = F_mass*0.7
    ash_mass.value = F_mass*0.006
    solids_mass.value = F_mass*0.015
    lipid_mass[:] = z_mass_new_lipid * F_mass
    carbs_mass[:] = carbs_massfrac * z_mass_new_carbs * F_mass
    fiber_mass[:] = fiber_massfrac * new_fiber_massfrac * F_mass
    if any(lc.mol < 0):
        raise ValueError(f'lipid cane oil composition of {z_mass_lipid*100:.0f}% dry weight is infeasible')


### Old attempt to include energy balance. It does not work because fiber is more energy dense than sugar (to my surprise) ###

# LHV_0 = lc.LHV
# def LHV_error(carbs_massnet, dryweight_no_oil):
#     """Return the difference between the original heat of combustion and the heat of combustion at the new carbohydrate composition"""
#     carbs_mass[:] = carbs_massfrac * carbs_massnet
#     fiber_mass[:] = fiber_massfrac * (dryweight_no_oil - carbs_massnet)
#     return LHV_0 - lc.LHV

# def composition_balance(oilfrac):
#     """Adjust composition of lipid cane to achieve desired oil fraction (dry weight)."""
#     global carbs_massnet
#     arg = oilfrac
#     oil_mass[:] = mass_oil = oilfrac * dryweight
#     dryweight_no_oil = dryweight - mass_oil
#     carbs_massnet = newton(LHV_error, carbs_massnet, args=(dryweight_no_oil,))
#     if any(lc.mol < 0):
#         raise ValueError(f'Lipid cane oil composition of {arg*100:.0f}% dry weight is infeasible')
#     return carbs_massnet
     
# composition_balance(0.1)
# check = lambda: oil_mass.sum()/dryweight
    