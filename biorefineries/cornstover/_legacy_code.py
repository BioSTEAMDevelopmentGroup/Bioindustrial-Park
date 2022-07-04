# -*- coding: utf-8 -*-
"""
This module contains legacy code for the cellulosic ethanol biorefinery
based on the 2011 NREL report [1]_.

References
----------
.. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
    Dudgeon, D. (2011). Process Design and Economics for Biochemical 
    Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
    Pretreatment and Enzymatic Hydrolysis of Corn Stover
    (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269

"""

# %% Legacy code

### DO NOT DELELTE:
### This code is used to generate splits for unit operations in wastewater 
### treatment of a cellulosic ethanol biorefinery.
import thermosteam as tmo
import numpy as np
from biorefineries.cornstover import chemicals
tmo.settings.set_thermo(chemicals)
chemical_groups = dict(
        OtherSugars = ('Arabinose',
                        'Mannose',
                        'Galactose',
                        'Cellobiose',
                        'Sucrose'),
        SugarOligomers = ('GlucoseOligomer',
                          'XyloseOligomer',
                          'GalactoseOligomer',
                          'ArabinoseOligomer',
                          'MannoseOligomer'),
        OrganicSolubleSolids = ('AmmoniumAcetate',
                                'SolubleLignin',
                                'Extract', 
                                'LacticAcid', 
                                'Cellulase'),
        InorganicSolubleSolids = ('AmmoniumSulfate',
                                  'DAP',
                                  'NaOH',
                                  'HNO3',
                                  'NaNO3'),
        Furfurals = ('Furfural',
                      'HMF'),
        OtherOrganics = ('Glycerol',
                          'Denaturant',
                          'Oil',
                          'SuccinicAcid',
                          'Xylitol'),
        COxSOxNOxH2S = ('NO',
                        'NO2',
                        'SO2',
                        'CO',
                        'H2S'),
        Protein = ('Protein',
                    'Enzyme',
                    'DenaturedEnzyme'),
        CellMass = ('WWTsludge',
                    'Z_mobilis',
                    'T_reesei'),
        OtherInsolubleSolids = ('Tar',
                                'Ash',
                                'Lime'),
        OtherStructuralCarbohydrates = ('Arabinan', 
                                        'Mannan', 
                                        'Galactan')
)

def find_split(IDs, flow0, flow1):
    flow0 = np.asarray(flow0)
    splits = flow0/(flow0 + np.asarray(flow1))
    chemicals = tmo.settings.get_chemicals()
    array = np.zeros(chemicals.size)
    for ID, split in zip(IDs, splits):
        if ID in chemical_groups:
            array[chemicals.get_index(chemical_groups[ID])] = split
        else:
            array[chemicals.index(ID)] = split
    return array

splits = [
    ('Ethanol', 1, 15),
    ('Water', 27158, 356069),
    ('Glucose', 3, 42),
    ('Xylose', 7, 85),
    ('OtherSugars', 13, 175),
    ('SugarOligomers', 10, 130),
    ('OrganicSolubleSolids', 182, 2387),
    ('InorganicSolubleSolids', 8, 110),
    ('Ammonia', 48, 633),
    ('AceticAcid', 0, 5),
    ('Furfurals', 5, 70),
    ('OtherOrganics', 9, 113),
    ('Cellulose', 19, 6),
    ('Xylan', 6, 2),
    ('OtherStructuralCarbohydrates', 1, 0),
    ('Lignin', 186, 64),
    ('Protein', 51, 18),
    ('CellMass', 813, 280),
    ('OtherInsolubleSolids', 68, 23)
]

anaerobic_bioreactor_sludge_splits = find_split(*zip(*splits))

splits = [
    ('Glucose', 19, 502),
    ('Xylose', 40, 1022),
    ('OtherSugars', 81, 2175),
    ('SugarOligomers', 60, 1552),
    ('OrganicSolubleSolids', 612, 15808),
    ('InorganicSolubleSolids', 97, 2513),
    ('Furfurals', 19, 513),
    ('OtherOrganics', 52, 1348),
    ('Glucan', 1230, 25),
    ('Xylan', 415, 8),
    ('OtherStructuralCarbohydrates', 94, 2),
    ('Lignin', 12226, 250),
    ('Protein', 3376, 69),
    ('CellMass', 925, 19),
    ('OtherInsolubleSolids', 4489, 92)
]

pressure_filter_splits = find_split(*zip(*splits))

splits = [
    ('Ethanol', 0, 1),
    ('Water', 381300, 2241169),
    ('Glucose', 0, 2),
    ('Xylose', 1, 3),
    ('OtherSugars', 1, 7),
    ('SugarOligomers', 1, 6),
    ('OrganicSolubleSolids', 79, 466),
    ('InorganicSolubleSolids', 4828, 28378),
    ('Ammonia', 3, 16),
    ('Furfurals', 0, 3),
    ('OtherOrganics', 1, 7),
    ('CarbonDioxide', 6, 38),
    ('O2', 3, 17),
    ('N2', 5, 32),
    ('Cellulose', 0, 194),
    ('Xylan', 0, 65),
    ('OtherStructuralCarbohydrates', 0, 15),
    ('Lignin', 0, 1925),
    ('Protein', 0, 90),
    ('CellMass', 0, 19778),
    ('OtherInsolubleSolids', 0, 707)
]

membrane_bioreactor_splits = find_split(*zip(*splits))

centrifuge_species = ('Water', 'Glucose', 'Xylose', 'OtherSugars',
                      'SugarOligomers', 'OrganicSolubleSolids',
                      'InorganicSolubleSolids', 'Ammonia', 'Furfurals', 
                      'OtherOrganics', 'CO2', 'COxSOxNOxH2S', 'Cellulose',
                      'Xylan', 'OtherStructuralCarbohydrates', 'Lignin',
                      'Protein', 'CellMass', 'OtherInsolubleSolids')
S623_flow = np.array([7708, 0, 0, 1, 1, 13, 75, 3, 0, 1, 1, 2, 25, 8, 2, 250, 52, 1523, 92])
S616_flow = np.array([109098, 3, 6, 13, 9, 187, 1068, 46, 5, 8, 14, 31, 1, 0, 0, 13, 3, 80, 5])

sludge_centrifuge_splits = find_split(centrifuge_species, S616_flow, S623_flow)