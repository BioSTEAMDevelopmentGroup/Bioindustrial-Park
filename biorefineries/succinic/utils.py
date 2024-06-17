#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat
"""


# %% Setup

import numpy as np
import pandas as pd
import thermosteam as tmo
from biorefineries.succinic.chemicals_data import chems
_kg_per_ton = 907.18474

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
CEPCI = {1997: 386.5,
         1998: 389.5,
         2007: 525.4,
         2009: 521.9,
         2010: 550.8,
         2011: 585.7,
         2012: 584.6,
         2013: 567.3,
         2014: 576.1,
         2016: 541.7}

# %% 

# =============================================================================
# Function to find the split ratios for Splitters, assume 0 for chemicals not specified in splits
# =============================================================================

def find_split(IDs, flow0, flow1, chemical_groups):
    # Add 1e-6 to avoid flow0 and flow1 both being 0
    flow0 = np.asarray(flow0) + 1e-6
    flow1 = np.asarray(flow1) + 1e-6    
    splits = flow0/(flow0 + flow1)
    thermo = tmo.settings.get_thermo()
    chemicals = thermo.chemicals
    array = np.zeros(chemicals.size)  
    for ID, split in zip(IDs, splits):
        if ID in chemical_groups:
            array[chemicals.get_index(chemical_groups[ID])] = split
        else:
            array[chemicals.index(ID)] = split
    # WWTsludge is removed from the cell mass group 
    array[chemicals.index('WWTsludge')] = array[chemicals.index('FermMicrobe')]
    return array

IDs = ('Ethanol', 'H2O', 'Glucose', 'Xylose', 'OtherSugars',
    'SugarOligomers', 'OrganicSolubleSolids', 'InorganicSolubleSolids', 'Ammonia', 'AceticAcid', 
    'SulfuricAcid', 'Furfurals', 'OtherOrganics', 'CO2', 'CH4',
    'O2', 'N2', 'COSOxNOxH2S', 'Glucan', 'Xylan', 
    'OtherStructuralCarbohydrates', 'Acetate', 'Lignin', 'Protein', 'CellMass',
    'OtherInsolubleSolids')

streams = {}

streams['stream_535'] = (177, 329030, 502, 1022, 2094,
                         1552, 15808, 2513, 0, 0,
                         0, 513, 1348, 0, 0,
                         0, 0, 0, 25, 8,
                         2, 0, 250, 69, 19,
                         92)

streams['stream_571'] = (6, 12797, 19, 49, 81,
                         60, 612, 97, 0, 0,
                         0, 19, 52, 0, 0,
                         1, 1, 0, 1230, 415,
                         94, 0, 12226, 3376, 925,
                         4489)

streams['stream_611'] = (15, 356069, 42, 85, 175,
                         130, 2387, 110, 633, 5, 
                         0, 70, 113, 181, 3, 
                         1, 0, 300, 6, 2, 
                         0, 0, 64, 18, 280,
                         23)
streams['stream_612'] = (1, 27158, 3, 7, 13,
                         10, 182, 8, 48, 0,
                         0, 5, 9, 14, 0,
                         0, 0, 23, 19, 6,
                         1, 0, 186, 51, 813,
                         68)

streams['stream_616'] = (1, 109098, 3, 6, 13,
                         9, 187, 1068, 46, 0,
                         0, 5, 8, 14, 0, 
                         1, 1, 31, 1, 0,
                         0, 0, 13, 3, 80, 
                         5)

streams['stream_623'] = (0, 7708, 0, 0, 1,
                         1, 13, 75, 3, 0,
                         0, 0, 1, 1, 0,
                         0, 0, 2, 25, 8,
                         2, 0, 250, 52, 1523,
                         92)

streams['stream_624'] = (0, 381300, 0, 1, 1,
                         1, 79, 4828, 3, 0,
                         0, 0, 1, 6, 0,
                         3, 5, 44, 0, 0,
                         0, 0, 0, 0, 0,
                         0)

streams['stream_625'] = (1, 2241169, 2, 3, 7,
                         6, 466, 28378, 16, 0,
                         0, 3, 7, 38, 0,
                         17, 32, 259, 194, 65,
                         15, 0, 1925, 90, 19778,
                         707)

stream_626 = (0,) + (376324,) + (0,) * (len(IDs)-2)
streams['stream_626'] = stream_626

streams['stream_627'] = (0, 4967, 0, 1, 1,
                         1, 79, 2828, 3, 0,
                         0, 0, 1, 3, 0,
                         0, 0, 44, 0, 0,
                         0, 0, 0, 0, 0,
                         0)

splits_df = pd.DataFrame.from_dict(streams)
splits_df.index = IDs


# %% 

# =============================================================================
# Function to get feedstock flow by giving dry weight composition and moisture content
# =============================================================================

def get_feedstock_flow(dry_composition, moisture_content, dry_flow):
    dry_array = chems.kwarray(dry_composition)
    wet_flow = dry_flow / (1-moisture_content)
    moisture_array = chems.kwarray(dict(Water=moisture_content))
    feedstock_flow = wet_flow * (dry_array*(1-moisture_content)+moisture_array)
    return feedstock_flow

# !!! This is the dry composition of corn stover; update to the composition you need
dry_composition = dict(
    Glucan=0.3505, Xylan=0.1953, Lignin=0.1576, Ash=0.0493, Acetate=0.0181,
    Protein=0.0310, Arabinan=0.0238, Galactan=0.0143, Mannan=0.0060, 
    Sucrose=0.0077, Extract=0.1465, SuccinicAcid=0)

moisture_content = 0.2 #!!! This is the moisture content of preprocessed corn stover; update to the moisture content you need
dry_feedstock_flow = 2205. * _kg_per_ton / 24. # !!! Changing this will alter the total flow without altering the composition
 
baseline_feedflow = get_feedstock_flow(dry_composition, moisture_content, 
                                       dry_feedstock_flow)

# %% 

# =============================================================================
# Functions to compute chemical loading and adjust recycle flows to maintain
# a certain ratio for Esterification and Hydrolysis reactor
# =============================================================================

def compute_extra_chemical(feed, recycle, reactants_ID, chemical_ID, ratios):
    reactants_in_feed = feed.imol[reactants_ID]
    reactants_in_recycle = recycle.imol[reactants_ID]
    chemical_needed = (ratios*(reactants_in_feed+reactants_in_recycle)).sum()
    chemical_extra = (feed.imol[chemical_ID]+recycle.imol[chemical_ID]) - chemical_needed
    return chemical_extra

def adjust_recycle(feed, recycle, reactants_ID, chemical_ID, ratios):
    feed_chemical_needed = (feed.imol[reactants_ID]*ratios).sum() \
        - feed.imol[chemical_ID]
    
    recycle_chemical_extra = recycle.imol[chemical_ID] \
        - (recycle.imol[reactants_ID]*ratios).sum()
    
    split = feed_chemical_needed / recycle_chemical_extra
    effluent = feed.copy()
    recycle_recycled = recycle.copy()
    recycle_recycled.mol *= split
    recycle_discarded = recycle.copy()
    recycle_discarded.mol *= (1 - split)
    effluent.mix_from([feed, recycle_recycled])
    
    return effluent, recycle_discarded


# %% 

# =============================================================================
# Function to output chemical properties
# =============================================================================

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


# %% 

# =============================================================================
# Function for quick result checking
# =============================================================================

def get_sugar_conc(stream, sugars=()):
    fermentable_sugar = 0
    for sugar in sugars:
        fermentable_sugar += stream.imass[sugar]
    fermentable_sugar_conc = fermentable_sugar/stream.F_vol    
    return fermentable_sugar_conc


#%% Estimate MPSP at a given IRR from results for IRR at a set product selling price (assuming uniform cash flows across all years)

def get_MPSP_at_IRR(new_IRR, # fraction; e.g., 10% is 0.10
                      old_IRR, # fraction; e.g., 10% is 0.10
                      old_selling_price, # USD/kg
                      yearly_sales_mass_of_main_product, #kg/y
                      t, # y
                      initial_investment, # USD
                      old_NPV=0., # USD
                      ):
    
    t = int(round(t))
    
    yearly_cash_flow_without_product_sale = \
        + yearly_sales_mass_of_main_product * old_selling_price \
        - ((old_NPV + initial_investment) / (sum([1/((1+old_IRR)**(i+1)) for i in range(t)])))
    
    print(yearly_cash_flow_without_product_sale)
    MPSP_new = \
        (1/yearly_sales_mass_of_main_product) *\
        (yearly_cash_flow_without_product_sale \
        +  (initial_investment / (sum([1/((1+new_IRR)**(i+1)) for i in range(t)])))) \
        
    return round(MPSP_new, 3) # USD/kg


#%% Estimate MPSP assuming all installed costs are 0

def get_MPSP_without_capital_cost(product_stream, TEA, flowsheet, sys, product_chem_ID=None, N_sims=1, N_tea_solves=20):
    for i in range(N_sims):
        sys.simulate()
    for u in flowsheet.unit:
        for k in u.installed_costs.keys():
            u.installed_costs[k] = 0.
    MPSP = None
    for i in range(N_tea_solves):
        product_stream.price = TEA.solve_price(product_stream)
        if product_chem_ID:
            MPSP = product_stream.price*product_stream.F_mass/product_stream.imass[product_chem_ID]
        else:
            MPSP = product_stream.price
    return MPSP


#%% Model utils

