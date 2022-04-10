"""
Created on Sat Jan 16 19:09:17 2021

@author: sarangbhagwat
"""

#%% Imports
from biorefineries.TAL.system_TAL_adsorption_glucose import R302
from biorefineries.make_a_biorefinery.analyses.solvents_barrage import run_solvents_barrage

#%% Run
run_solvents_barrage(stream=R302.outs[0], # Stream from which you wish to extract the solute
                     solute_ID='TAL', # solute chemical ID
                     impurity_IDs=['VitaminA', 'VitaminD2'], # List of IDs of impurities in "stream" that you want get partitioning results for, other than water; note that all chemicals in the stream will affect LLE interaction effects, regardless of which chemicals are present in impurity_IDs
                     T=50+273.15, # Temperature (K) at which you wish to run the solvents barrage; temperature (K) of "stream" by default
                     stream_modifiers='baseline_stream', # String: 'baseline_stream' to analyze the "stream" passed in arguments; 'impurity_free_stream' to remove the impurities listed in impurity_IDs before performing analyses; 'solute_in_pure_water' to analyze simply for the solute in pure water
                     solvent_IDs = [
                                     # 'AceticAcid',
                                     'Pentadecanol',
                                     'Tridecanol',
                                     'Ethanol',
                                     'Methanol',
                                     'Propyl acetate',
                                     'Butyl acetate',
                                     'Hexanol',
                                     'Hexane',
                                     'Cyclohexanol',
                                     'Cyclohexanone',
                                     'Heptanol',
                                     'Octanol',
                                     '1,8-Octanediol',
                                     '2-Ethyl hexanol',
                                     'Nonanol',
                                     'Decanol',
                                     'Dodecanol',
                                     '117-81-7', # CAS number for Dioctyl (Diethylhexyl) phthalate
                                     'Diethyl sebacate', # No Psat, Hvap
                                     # 'Glycerol',
                                     'Toluene',
                                     'Trioctylamine',
                                     'Isoamyl alcohol',
                                     '5137-55-3', # CAS number for Aliquat 336
                                     # 'Water',
                                     'Benzene',
                                     '143-28-2', # CAS number for Oleyl alcohol
                                     'Tetrahydrofuran'
                                     ]
                     )