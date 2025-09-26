# -*- coding: utf-8 -*-
"""
Created on 2025-07-07 14:55:35

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import biosteam as bst

__all__ = ('load_process_settings','price','set_GWPCF','set_FECCF','GWP_CFs','FEC_CFs','set_GWPCF_Multi')

# %% Prices

# =============================================================================
# Prices for techno-economic analysis (TEA), all in $/kg (electricity in $/kWh)
# and from ref [1] if not noted
# =============================================================================

price = {
    'H2O': 0.2/1e3,  # $/kg 0.15 to 0.5 /1e3
    'H2SO4': 0.12,  # $/kg 0.1 to 0.3
    'NH3': 0.46,  # $/kg 0.45 to 0.65
    'NH4OH': 0.4, # $/kg 0.3 to 0.5
    'NH3_25wt': 0.46,  # $/kg 0.3 to 0.5
    'NH4SO4': 0.16,  # $/kg 0.15 to 0.45
    'Glucose': 0.42,  # $/kg 0.5 to 0.9
    'MgSO4': 0.175,  # $/kg 0.3 to 0.6
    'NaOH': 0.28,  # $/kg 0.13 to 0.4
    'KH2PO4': 1.31,  # $/kg 1.2 to 2
    'NaH2PO4': 0.73,  # $/kg 0.5 to 0.82
    'Na2HPO4': 0.55,  # $/kg 0.35 to 0.65
    'K2HPO4': 1.7,  # $/kg 1.3 to 2.2
    'FeSO4': 0.055,  # $/kg 0.4 to 0.8
    'NH42HPO4': 0.569,  # $/kg 0.4 to 0.7
    'NaCl': 0.078,  # $/kg 0.1 to 0.25
    'KCl': 0.23,  # $/kg 0.1 to 0.25
    'EDTA': 2.43,  # $/kg 2 to 4
    'K2HPO4': 1.7,  # $/kg 1.3 to 2.2
    'Glycine': 1.75,  # $/kg 1.5 to 2.5
    'TrehaloseDH': 1.5,  # $/kg 15 to 40
    'SodiumAscorbate': 4.35,  # $/kg 5 to 10
    'Glycerol': 0.45,  # $/kg 0.3 to 0.7
    'sugar': 0.45,  # $/kg
    'Electricity': 0.065,  # $/kWh 0.035 to 0.25
    'Low pressure steam': 0.30626,  # $/kg
    'Cooling water': 0,  # $/kg
}
# %%
# =============================================================================
# Characterization factors (CFs) for life cycle analysis (LCA), all from ref [5] if not noted otherwise
# =============================================================================

# Individual characterization factor dictionaries for each impact category
##### 100-year global warming potential (GWP) in kg CO2-eq/kg unless noted otherwise #####
GWP_CFs = {
    'ElectricitySG': (0.55, 0.55), # from ecoinvent 3.11 cutoff, electricity 0.543 to 0.553
    'Electricity': (0.55, 0.55), # from ecoinvent 3.11 cutoff, electricity, GLO
    'Glucose': 1.61, # from ecoinvent 3.11 cutoff, market for glucose, GLO
    'Sugar_Cane': 0.835, # from ecoinvent 3.11 cutoff, market for sugar, GLO
    'Sugar_Beet': 0.564, # from ecoinvent 3.11 cutoff, market for sugar, GLO
    'IronSulfate': 0.278, # from ecoinvent 3.11 cutoff, market for iron sulfate , RoW
    'AmmoniumSulfate': 0.179, # from ecoinvent 3.11 cutoff, market for ammonium sulfate, RoW
    'MagnesiumSulfate': 0.884, # from ecoinvent 3.11 cutoff, market for magnesium sulfate, GLO
    'Ammonia_SEA': 2.84, # from ecoinvent 3.11 cutoff, market for ammonia, SEA
    'Ammonia_CN': 5.07, # from ecoinvent 3.11 cutoff, market for ammonia, CN
    'Ammonia_US': 2.70, # from ecoinvent 3.11 cutoff, market for ammonia, US
    'AmmoniumSulfate': 0.857, # from ecoinvent 3.11 cutoff, market for ammonium sulfate, RoW
    'NaCl': 0.268, # from ecoinvent 3.11 cutoff, market for sodium chloride powder, GLO
    'Na2HPO4': 2.29, # from ecoinvent 3.11 cutoff, market for disodium phosphate, GLO
    'NaH2PO4': 2.85, # from ecoinvent 3.11 cutoff, market for monosodium phosphate, RoW
    'K2HPO4': 2.29, # Assumed same as Na2HPO4
    'KH2PO4': 2.85, # Assumed same as NaH2PO4
    'KCl': 0.494, # from ecoinvent 3.11 cutoff, market for potassium chloride, RoW
    'NaOH': 1.41, # from ecoinvent 3.11 cutoff, market for sodium hydroxide, Row
    'CO2': 1.0,  # Direct CO2 emissions
}

# ##### Fossil energy consumption (FEC), in MJ/kg of material unless noted otherwise #####
# FEC_CFs = {
#     'Electricity': (5.926, 5.926), # assume production==consumption, both in MJ/kWh
#     'H2SO4': 568.98/1e3,
#     'NaOH': 29,
#     'NH4OH': 42 * 0.4860, # chemicals.NH3.MW/chemicals.NH4OH.MW,
#     'NH3_25wt': 42 * 0.25,
#     'CH4': 55.5,  # Natural gas FEC (placeholder - update with actual value)
#     # Add more chemicals as needed
# }

# # Consolidated CFs dictionary in the format expected by BioSTEAM LCA module
# # Format: {<impact_category>: {<chemical_ID_or_stream_name>: <CF_value>}}
# CFs = {
#     'GWP_100': {
#         # Extract single values from tuples for electricity, use first value
#         'Electricity': GWP_CFs['Electricity'][0] if isinstance(GWP_CFs['Electricity'], tuple) else GWP_CFs['Electricity'],
#         # Add all other GWP characterization factors
#         **{k: v for k, v in GWP_CFs.items() if k != 'Electricity'},
#         # Placeholder for complex feeds (feedstock streams) - to be updated per biorefinery
#         # 'Glucose': 0.0,  # kg CO2-eq per kg wet feedstock (uncomment and update as needed)
#         # 'Corn': 0.0,     # kg CO2-eq per kg wet feedstock (uncomment and update as needed)
#     },
#     'FEC': {
#         # Extract single values from tuples for electricity
#         'Electricity': FEC_CFs['Electricity'][0] if isinstance(FEC_CFs['Electricity'], tuple) else FEC_CFs['Electricity'],
#         # Add all other FEC characterization factors
#         **{k: v for k, v in FEC_CFs.items() if k != 'Electricity'},
#         # Placeholder for complex feeds (feedstock streams) - to be updated per biorefinery
#         # 'Glucose': 0.0,  # MJ per kg wet feedstock (uncomment and update as needed)
#         # 'Corn': 0.0,     # MJ per kg wet feedstock (uncomment and update as needed)
#     },
# }


# %% Process settings
def set_GWPCF(obj, name='', dilution=None):
    if not dilution: obj.characterization_factors['GWP'] = GWP_CFs[name]
    else: obj.characterization_factors['GWP'] = GWP_CFs[name] * dilution

def set_FECCF(obj, name='', dilution=None):
    if not dilution: obj.characterization_factors['FEC'] = FEC_CFs[name]
    else: obj.characterization_factors['FEC'] = FEC_CFs[name] * dilution

def set_GWPCF_Multi(obj, names=[], dilutions=None):
    """
    Set GWP characterization factors from multiple sources.
    
    Args:
        obj: Object to set characterization factors for
        names: List of names corresponding to GWP_CFs keys
        dilutions: List of dilution factors (same length as names) or single value for all
    """
    if not names:
        return
    
    if dilutions is None:
        # No dilution, just sum the CFs
        total_cf = sum(GWP_CFs[name] for name in names)
    elif isinstance(dilutions, (int, float)):
        # Single dilution value for all
        total_cf = sum(GWP_CFs[name] * dilutions for name in names)
    else:
        # List of dilutions - this handles array input
        if len(dilutions) != len(names):
            raise ValueError("Length of dilutions must match length of names")
        total_cf = sum(GWP_CFs[name] * dil for name, dil in zip(names, dilutions))
    
    obj.characterization_factors['GWP'] = total_cf

def load_process_settings():
    settings = bst.settings
    bst.process_tools.default_utilities()
    settings.CEPCI = 798.8 # 2024
    settings.electricity_price = 0.065 
    steam_utility = settings.get_agent('low_pressure_steam')
    steam_utility.heat_transfer_efficiency = 0.9
    steam_utility.regeneration_price = 0.30626
    steam_utility.T = 529.2
    steam_utility.P = 44e5
    settings.get_agent('cooling_water').regeneration_price = 0
    settings.get_agent('chilled_water').heat_transfer_price = 0
    bst.PowerUtility.price = price['Electricity']
    set_GWPCF(bst.PowerUtility, 'Electricity')
    #set_FECCF(bst.PowerUtility, 'Electricity')
    bst.settings.define_impact_indicator(key='GWP', units='kg*CO2e')