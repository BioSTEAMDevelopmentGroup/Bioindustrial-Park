# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""

import biosteam as bst
from thermosteam import Stream
from biorefineries.cellulosic import chemicals as c
from biorefineries.cellulosic import streams as s
from biorefineries.cellulosic import units
import numpy as np

__all__ = (
    'create_saccharification_system',
)

@bst.SystemFactory(
    ID='saccharification_sys',
    ins=[s.pretreated_biomass,
         s.cellulase,
         s.saccharification_water],
    outs=[s.slurry],
)
def create_saccharification_system(
        ins, outs, solids_loading=None, 
        nonsolids=None,
        insoluble_solids=None,
        insoluble_solids_loading=None,
        ignored=None,
        Saccharification=None,
        saccharification_reactions=None,
    ):
    pretreated_biomass, cellulase, saccharification_water = ins
    slurry, = outs
    if nonsolids is None: nonsolids = c.default_nonsolids
    if insoluble_solids is None: insoluble_solids = c.default_insoluble_solids
    if insoluble_solids_loading is None: insoluble_solids_loading = 0.103
    if solids_loading is None: solids_loading = 0.2
    if ignored is None: ignored = c.default_ignored
    M301 = units.EnzymeHydrolysateMixer('M301', (pretreated_biomass, cellulase, saccharification_water))
    H301 = units.HydrolysateHeatExchanger('H301', M301-0, T=48+273.15)
    M301.solids_loading = solids_loading
    M301.insoluble_solids_loading = insoluble_solids_loading
    M301.enzyme_loading = 0.02 # (20 g enzyme / 1000 g cellulose) 
    M301.enzyme_concentration = 0.05 # (50 g cellulase / 1000g cellulose mixture)
    M301.loading_basis = lambda: 1.2 * (pretreated_biomass.imass['Glucan'])
    
    @M301.add_specification
    def update_cellulase_loading():
        pretreated_biomass, cellulase, water = M301.ins
        # Note: An additional 10% is produced to produce sophorose
        # Humbird (2011) pg. 37 
        enzyme_concentration = M301.enzyme_concentration
        enzyme_over_cellulose = M301.enzyme_loading / enzyme_concentration
        z_mass_cellulase_mixture = np.array([1 - enzyme_concentration , enzyme_concentration], float)
        cellulase.imass['Water', 'Cellulase'] = (
            enzyme_over_cellulose
            * z_mass_cellulase_mixture
            * M301.loading_basis()
        )
    
    @M301.add_specification(run=True)
    def update_moisture_content():
        pretreated_biomass, cellulase, saccharification_water, *other = M301.ins
        chemicals = M301.chemicals
        s_mix = Stream.sum([pretreated_biomass, cellulase], None, M301.outs[0].thermo, energy_balance=False)
        mass = s_mix.mol * chemicals.MW
        solids_loading = M301.solids_loading
        insoluble_solids_loading = M301.insoluble_solids_loading
        indices = chemicals.available_indices(nonsolids)
        mass_moisture = mass[indices].sum()
        indices = chemicals.available_indices(insoluble_solids)
        mass_insoluble = mass[indices].sum()
        indices = chemicals.available_indices(ignored)
        if indices:
            mass_ignored = mass[indices].sum()
        else:
            mass_ignored = 0.
        total_mass = mass.sum() - mass_ignored
        water_over_solids = (1 - solids_loading) / solids_loading
        water_over_insoluble_solids = (1 - insoluble_solids_loading) / insoluble_solids_loading
        M301.required_saccharification_water = max(
            water_over_solids * (total_mass - mass_moisture) - mass_moisture,
            water_over_insoluble_solids * mass_insoluble - (total_mass - mass_insoluble)
        )
        saccharification_water.imass['Water'] = max(M301.required_saccharification_water, 0.)
    
    if Saccharification is None:
        Saccharification = units.ContinuousPresaccharification
    Saccharification('R301', H301-0, slurry, reactions=saccharification_reactions)
