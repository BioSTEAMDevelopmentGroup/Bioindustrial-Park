# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 19:09:17 2021

@author: sarangbhagwat
"""
# %% Imports and chemicals initialization

import numpy as np
import thermosteam as tmo
import biosteam as bst
import copy 

from biorefineries.HP.chemicals_data import HP_chemicals
from biorefineries.HP.system_light_lle_vacuum_distillation import S404
# tmo.settings.set_thermo(['Water', 'octanol', 'hexanol', 'butyl acetate', HP_chemicals['Xylose'], HP_chemicals['Glucose'], HP_chemicals['Triacetic acid lactone'], 'isoamyl alcohol'])

# tmo.settings.set_thermo(HP_chemicals)

Water = HP_chemicals['Water']
# Glucose = HP_chemicals['Glucose']
# HP = HP_chemicals['Triacetic acid lactone']
Glycerol = HP_chemicals['Glycerol']
AQ336 = HP_chemicals['AQ336']
Octanol = HP_chemicals['Octanol']
Hexanol = tmo.Chemical('Hexanol')
Heptanol = tmo.Chemical('Heptanol')
Butyl_acetate = tmo.Chemical('Butyl acetate')
Propyl_acetate = tmo.Chemical('Propyl acetate')
Isoamyl_alcohol = tmo.Chemical('Isoamyl alcohol')
TOA = tmo.Chemical('Trioctylamine')
Decanol = tmo.Chemical('Decanol')
Dodecanol = tmo.Chemical('Dodecanol')
Nonanol = tmo.Chemical('Nonanol')
te_hexanol = tmo.Chemical('2-Ethyl hexanol')
Cyclohexanol = tmo.Chemical('Cyclohexanol')
Cyclohexanone = tmo.Chemical('Cyclohexanone')
Dioctyl_phthalate = tmo.Chemical('117-81-7')
Diethyl_sebacate = tmo.Chemical('Diethyl sebacate')
Diethyl_sebacate.copy_models_from(Water, ['Psat', 'Hvap'])
Octanediol = HP_chemicals['Octanediol']
Toluene = HP_chemicals['Toluene']
# AceticAcid = HP_chemicals['AceticAcid']
# Glycerol = tmo.Chemical('Glycerol')
# H2SO4 = tmo.Chemical('H2SO4')
# TAL = tmo.Chemical('Triacetic acid lactone')
# Furfural = tmo.Chemical('Furfural')

# TAL.copy_models_from(Furfural, ['Psat', 'Hvap', 'V'])
# TAL.Hfus = 30883.6698 # !!! from solubility modeling method 4(ii)
# TAL.Tm = 185 + 273.15
# TAL.Tb = 239.1 + 273.15

# Glycerol.copy_models_from(Furfural, ['V',])

#%% Set solute and solvents
solute = HP_chemicals['HP']
solvents = [Propyl_acetate, Butyl_acetate, Hexanol, Cyclohexanol, Cyclohexanone, Heptanol, Octanol, Octanediol, te_hexanol, Nonanol, Decanol, Dodecanol, Isoamyl_alcohol, Dioctyl_phthalate, Diethyl_sebacate, Glycerol, Toluene]

# test_solvents_chemicals = copy.deepcopy(HP_chemicals)
# test_solvents_chemicals.extend(solvents)

#%% Set thermo
# tmo.settings.set_thermo(solvents + ['Water', 'H2SO4', Glycerol, solute, HP_chemicals['Xylose'], HP_chemicals['Glucose'], HP_chemicals['AceticAcid']])
tmo.settings.set_thermo(list(HP_chemicals.tuple) + solvents)
# %% Streams initialization

T = 80. + 273.15
# process_stream = tmo.Stream('process_stream',
#                             Water = 4.88e4,
#                             units = 'kg/hr',
#                             T = T)
# process_stream.imass[solute.ID] = 2.12e+04
# process_stream.imol['AceticAcid'] = 5.
# process_stream.imol['Glycerol'] = 5.

process_stream = S404.ins[0].copy()
process_stream.T = T
solvent_stream = tmo.Stream('solvent_stream',
                            T = T)

mixed_stream = tmo.Stream('mixed_stream')
extract_phase = 'l'
raffinate_phase = 'L'

# %%% Functions

def get_K(chem_ID, stream, phase_1, phase_2):
    return (stream[phase_1].imol[chem_ID]/stream[phase_1].F_mol)/(stream[phase_2].imol[chem_ID]/stream[phase_2].F_mol)

def set_solvent(solvent_chemical, solvent_mol=1314, solvent_stream=solvent_stream):
    solvent_stream.empty()
    if type(solvent_chemical) is str:
        solvent_chemical = tmo.Chemical(solvent_chemical)
    solvent_ID = solvent_chemical.ID
    solvent_stream.imol[solvent_ID] = solvent_mol
    
def run_single_test(solvent_chemical, solvent_mol=1314, process_stream=process_stream, solvent_stream=solvent_stream):
    mixed_stream.empty()
    set_solvent(solvent_chemical, solvent_mol, solvent_stream)
    solvent_ID = solvent_chemical.ID
    mixed_stream.mix_from([process_stream, solvent_stream])
    # print(mixed_stream.lle_chemicals)
    mixed_stream.lle(T=T, top_chemical = solvent_ID)
    # mixed_stream.lle(T=T)
    # mixed_stream.show(N=100, composition=True)
    K_HP_in_solvent = get_K(solute.ID, mixed_stream, extract_phase, raffinate_phase)
    # K_Glucose_in_solvent = get_K('Glucose', mixed_stream, extract_phase, raffinate_phase)
    print(K_HP_in_solvent)
    K_Water_in_solvent = get_K('Water', mixed_stream, extract_phase, raffinate_phase)
    print(K_Water_in_solvent)
    K_solvent_in_Water = get_K(solvent_ID, mixed_stream, raffinate_phase, extract_phase)
    print(K_solvent_in_Water)
    # K_Water_in_Water = get_K('Water', mixed_stream, raffinate_phase, extract_phase)
    T_diff =  solute.Tb - solvent_chemical.Tb
    print(T_diff)
    print('\n')
    K_acetic_acid_in_solvent = 0.000001
    K_glycerol_in_solvent = 0.000001
    try:
        K_acetic_acid_in_solvent = get_K('AceticAcid', mixed_stream, raffinate_phase, extract_phase)
    except:
        pass
    try:
        K_glycerol_in_solvent = get_K('Glycerol', mixed_stream, raffinate_phase, extract_phase)
    except:
        pass
    return (mixed_stream.copy(), K_HP_in_solvent, K_Water_in_solvent, K_solvent_in_Water, T_diff, K_acetic_acid_in_solvent, K_glycerol_in_solvent)

def get_results_meeting_constraints(results_list, constraints):
    filtered_results = []
    for result in results_list:
        values = result[1]
        meets_constraints = True
        for i in range(1,len(constraints)):
            if not (values[i]>constraints[i-1][0] and values[i]<constraints[i-1][1]):
                meets_constraints = False
                # print(i)
                # print(values[i],constraints[i-1][0],values[i],constraints[i-1][1])
        if meets_constraints:
            filtered_results.append(result)
    return filtered_results
# %% Run test barrage
results_dict = {}
for solvent in solvents:
    results_dict[solvent.common_name] = run_single_test(solvent)

results_list = list(results_dict.items())
# results_list is a list of tuples each of the following format:
# (solvent_common_name, (mixed_stream, K_HP_in_solvent, K_Water_in_solvent, K_solvent_in_Water, T_diff))

# %% Sort results by criterion
key = None

criterion = 'Partition of solute into extract'
# criterion = 'Partition of water into extract'
# criterion = 'Partition of solvent into raffinate'
# criterion = 'Difference between solute and solvent boiling temperatures'
# criterion = 'Absolute difference between solute and solvent boiling temperatures'

if criterion == 'Partition of solute into extract':
    key = lambda i: i[1][1]
elif criterion == 'Partition of water into extract':
    key = lambda i: i[1][2]
elif criterion == 'Partition of solvent into raffinate':
    key = lambda i: i[1][3]
elif criterion == 'Difference between solvent and solute boiling temperatures':
    key = lambda i: i[1][4]
elif criterion == 'Absolute difference between solvent and solute boiling temperatures':
    key = lambda i: abs(i[1][4])

reverse = True

results_list.sort(key=key, reverse=reverse)

print(f"\n\nResult with highest {criterion.lower()}: \n")
print(results_list[0])
results_list[0][1][0].show()

# %% Filter results through constraints

inf = float('inf')

range_K_HP_in_solvent = (20., inf)
range_K_Water_in_solvent = (0., 0.2)
range_K_solvent_in_Water = (0., 0.02)
range_T_diff = (25, inf)
# range_abs_T_diff = (0., inf)

constraints = (range_K_HP_in_solvent, range_K_Water_in_solvent, range_K_solvent_in_Water,
               range_T_diff)

filtered_results = get_results_meeting_constraints(results_list, constraints)

print(f"\n\nResults meeting specified constraints ({len(filtered_results)} results out of {len(results_list)} total): \n")
for result in filtered_results:
    print(result)
    result[1][0].show()
    print('\n')
    
# %% Unit initialization and tests
solvent_to_run = '1-octanol'
set_solvent(solvent_to_run)
partition_data = dict(IDs=(solute.ID, 'Water', solvent_to_run,
                            'AceticAcid', 'Glycerol',), 
            K=np.array([1./results_dict[solvent_to_run][1],
                        1/results_dict[solvent_to_run][2],
                        results_dict[solvent_to_run][3],
                       results_dict[solvent_to_run][5],
                       results_dict[solvent_to_run][6]]),
            phi = 0.5)

MS = bst.units.MultiStageMixerSettlers('MS', ins = (process_stream, solvent_stream),
                                      outs = ('raffinate', 'extract'),
                                      N_stages = 5, partition_data = partition_data,) 

MS.simulate()
MS.show(N=100, composition=True, flow='kg/hr')


# %% For 3-HP only

import pandas as pd
from datetime import datetime

compiled_results_dict = {}
for result in results_list:
    compiled_results_dict[result[0]] = rn = {}
    rn['K_solute_in_extract'] = result[1][1]
    rn['K_water_in_extract'] = result[1][2]
    rn['K_solvent_in_water'] = result[1][3]
    rn['K_acetic_acid_in_extract'] = result[1][5]
    rn['K_glycerol_in_extract'] = result[1][6]
    rn['Tb_solute - Tb_solvent'] = result[1][4]

results_df = pd.DataFrame(compiled_results_dict)

dateTimeObj = datetime.now()
file_to_save = 'Solvents_barrage_3-HP_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)

map_dict = {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H', 8: 'I', 
9: 'J', 10: 'K', 11: 'L', 12: 'M', 13: 'N', 14: 'O', 15: 'P', 16: 'Q', 
17: 'R', 18: 'S', 19: 'T', 20: 'U', 21: 'V', 22: 'W', 23: 'X', 24: 'Y', 25: 'Z'}

with pd.ExcelWriter(file_to_save+'.xlsx') as writer:
    final_results_df = results_df.transpose()
    final_results_df.to_excel(writer, sheet_name='All solvent candidates')
    workbook  = writer.book
    worksheet = writer.sheets['All solvent candidates']
    # wrap_format = workbook.add_format({'text_wrap': True})
    # worksheet.set_column('A:A', None, wrap_format)
    # worksheet.set_row('A:A', None, wrap_format)
    # writer.save()
    # Add a header format.
    header_format = workbook.add_format({
        'bold': True,
        'text_wrap': True,
        'valign': 'top',
        # 'fg_color': '#D7E4BC',
        'border': 1})
    decimal_2_format = workbook.add_format({'num_format': '#,##0.00'})
    decimal_3_format = workbook.add_format({'num_format': '#,##0.000'})
    worksheet.set_column('A:A', 18, decimal_2_format)
    
    # Write all cells in the specified number format.
    for i in range(len(final_results_df.columns.values)):
        worksheet.set_column(map_dict[i+1]+':'+map_dict[i+1], 18, decimal_2_format)
    # worksheet.set_row(0, 28, decimal_2_format)
    
    # Write the column headers with the defined format.
    for col_num, value in enumerate(final_results_df.columns.values):
        worksheet.write(0, col_num + 1, value, header_format)
    # Write the row headers with the defined format.
    for row_num, value in enumerate(final_results_df.index.values):
        worksheet.write(row_num + 1, 0, value, header_format)
    