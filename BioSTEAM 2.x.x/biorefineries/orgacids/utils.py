#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 14:17:17 2020

@author: yalinli_cabbi
"""


# %% Setup

import numpy as np
import pandas as pd
import thermosteam as tmo
from orgacids.chemicals import orgacids_chemicals
_kg_per_ton = 907.18474


# %% Function to find the split ratios for Splitters, assume 0 for chemicals not specified in splits

def find_split(IDs, flow0, flow1, chemical_groups):
    flow0 = np.asarray(flow0)
    splits = flow0/(flow0 + np.asarray(flow1))
    thermo = tmo.settings.get_thermo()
    chemicals = thermo.chemicals
    array = np.zeros(chemicals.size)  
    for ID, split in zip(IDs, splits):
        if ID in chemical_groups:
            array[chemicals.get_index(chemical_groups[ID])] = split
        else:
            array[chemicals.index(ID)] = split
    return array


# %% Function to get feedstock flow by giving dry weight composition and moisture content

def get_feedstock_flow(dry_composition, moisture_content, dry_flow):
    dry_array = orgacids_chemicals.kwarray(dry_composition)
    wet_flow = dry_flow / (1-moisture_content)
    moisture_array = orgacids_chemicals.kwarray(dict(Water=moisture_content))
    feedstock_flow = wet_flow * (dry_array*(1-moisture_content)+moisture_array)
    return feedstock_flow

dry_composition = dict(
    Glucan=0.3505, Xylan=0.1953, Lignin=0.1576, Ash=0.0493, Acetate=0.0181,
    Protein=0.0310, Arabinan=0.0238, Galactan=0.0143, Mannan=0.0060, 
    Sucrose=0.0077, Extract=0.1465, SuccinicAcid=0)

moisture_content = 0.2
dry_feedstock_flow = 2205 * _kg_per_ton / 24     
baseline_feedflow = get_feedstock_flow(dry_composition, moisture_content, 
                                       dry_feedstock_flow)


# %% Function to output chemical properties

def get_chemical_properties(chemicals, T, P, output=False):
    formulas = [chemical.formula for chemical in chemicals]
    MWs = [chemical.MW for chemical in chemicals]
    Hfs = [chemical.Hf for chemical in chemicals]
    HHVs = [chemical.HHV for chemical in chemicals]
    LHVs = [chemical.LHV for chemical in chemicals]
    phases = []
    Tbs = []
    Psats = []
    Vs = []
    Cns = []
    mus = []
    kappas = []
    
    for chemical in chemicals:
        if chemical.locked_state:
            phases.append(chemical.phase_ref)
            Tbs.append('NA')
            try: Psats.append(chemical.Psat(T=T, P=P))
            except: Psats.append('')
            try: Vs.append(chemical.V(T=T, P=P))
            except: Vs.append('')
            try: Cns.append(chemical.Cn(T=T))
            except: Cns.append('')
            try: mus.append(chemical.mu(T=T, P=P))
            except: mus.append('')
            try: kappas.append(chemical.kappa(T=T, P=P))
            except: kappas.append('')
        else:
            ref_phase = chemical.get_phase(T=T, P=P)
            phases.append(f'variable, ref={ref_phase}')
            Tbs.append(chemical.Tb)
            try: Psats.append(chemical.Psat(T=T, P=P))
            except: Psats.append('')
            try: Vs.append(chemical.V(ref_phase, T=T, P=P))
            except: Vs.append('')
            try: Cns.append(chemical.Cn(ref_phase, T=T))
            except: Cns.append('')
            try: mus.append(chemical.mu(ref_phase, T=T, P=P))
            except: mus.append('')
            try: kappas.append(chemical.kappa(ref_phase, T=T, P=P))
            except: kappas.append('')
    
    properties = pd.DataFrame(
        {'ID': chemicals.IDs,
          'formula': formulas,
          'MW': MWs,
          'HHV': HHVs,
          'LHV': LHVs,
          'Hf': Hfs,
          'phase': phases,
          'boiling point': Tbs,
          'Psat': Psats,
          'V': Vs,
          'Cn': Cns,
          'mu': mus,
          'kappa': kappas}
        )
    
    if output:
        properties.to_excel('chemical_properties.xlsx', sheet_name='properties')


# %% Function to get the sugar content of a stream

def get_sugar_conc(stream, sugars=()):
    fermentable_sugar = 0
    for sugar in sugars:
        fermentable_sugar += stream.imass[sugar]
    fermentable_sugar_conc = fermentable_sugar/stream.F_vol    
    return fermentable_sugar_conc





