# -*- coding: utf-8 -*-
# Copyright (C) 2022-2023, Sarang Bhagwat <sarangb2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Created on Sat Jan 16 19:09:17 2021

@author: sarangbhagwat
"""
# %% Imports and chemicals initialization

import numpy as np
import thermosteam as tmo
import biosteam as bst
import copy 
from matplotlib import pyplot as plt

#%% Default solvent IDs list
solvent_IDs = [
                # 'AceticAcid',
                'Pentadecanol',
                'Tridecanol',
                'Ethanol',
                'Methanol',
                'Propyl acetate',
                'Butyl acetate',
                'Hexanol',
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
    #%% Initialize chemicals
    
def run_solvents_barrage(stream, # Stream from which you wish to extract the solute
                         solute_ID, # solute chemical ID
                         impurity_IDs, # List of IDs of impurities in "stream" that you want get partitioning results for, other than water; note that all chemicals in the stream will affect LLE interaction effects, regardless of which chemicals are present in impurity_IDs
                         T=None, # Temperature (K) at which you wish to run solvents barrage; temperature (K) of "stream" by default
                         solvent_IDs=solvent_IDs, # List of solvents to run the barrage for; defaults to a list of 25 common organic solvents
                         stream_modifiers='baseline_stream', # 'baseline_stream' to analyze the "stream" passed in arguments; 'impurity_free_stream' to remove the impurities listed in impurity_IDs before performing analyses; 'solute_in_pure_water' to analyze simply for the solute in pure water
                         show_all_mixer_settlers=False,
                         plot_Ks=True,
                         save_excel=True,
                         save_K_plots=True,
                         criterion = 'Partition of solute into extract',
                         print_result_with_optimal_criterion=False): 
    test_env_chems = tmo.Chemicals([])
    borrowed_chemicals = stream.chemicals
    if not T:
        T=stream.T
    def chemical_database(ID, phase=None, **kwargs):
        chemical = tmo.Chemical(ID, **kwargs)
        if not chemical.CAS in [c.CAS for c in test_env_chems]:
            if phase:
                chemical.at_state(phase)
                chemical.phase_ref = phase
            test_env_chems.append(chemical)
            if not chemical.Psat:
                chemical.copy_models_from(tmo.Chemical('H2O'), ('Psat',))
            if not chemical.Hvap:
                chemical.copy_models_from(tmo.Chemical('H2O'), ('Hvap',))
            if chemical.CAS == '143-28-2':
                chemical.copy_models_from(tmo.Chemical('Octyldodecanol'),('Psat', 'mu', 'sigma'))
                chemical.Tb = tmo.Chemical('Octyldodecanol').Tb
        # database_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
        # test_env_chems.set_synonym(ID, chemical.common_name)
        return chemical
    
    
    #%% Select solute, impurities, and solvents
    solute = borrowed_chemicals[solute_ID]
    # solute_ID = 'TAL'
    # solute_ID = solute.ID
    # impurity_IDs = ['VitaminA', 'VitaminD2']
    
    
    #%% Load solvents
    def load_solvents_and_get_dict(solvent_IDs):
        solvents_dict = {}
        for solv_ID in solvent_IDs:
            solvents_dict[solv_ID] = chemical_database(solv_ID)
        return solvents_dict
    
    solvents_dict = load_solvents_and_get_dict(solvent_IDs)
    
    solvents = list(solvents_dict.values())
    
    #%% Load borrowed chemicals
    borrowed_chemicals_list = list(borrowed_chemicals.tuple)
    def load_chemicals(borrowed_chemicals_list):
        for chem in borrowed_chemicals_list:
            # try:
            #     chemical_database(chem.CAS)
            # except:
            if not chem.CAS in [c.CAS for c in test_env_chems]:
                test_env_chems.append(borrowed_chemicals[chem.ID])
            
    load_chemicals(borrowed_chemicals_list)
    test_env_chems.compile()
    
    def formatted_name(name):
        capitalized_name_chars = list(name) 
        for i in range(len(name)):
            if i==0:
                capitalized_name_chars[i]=name[i].upper()
            elif name[i]==' ' or name[i]=='-' and not i==len(name)-1:
                capitalized_name_chars[i+1]=name[i+1].upper()
        formatted_name_chars = [n for n in capitalized_name_chars \
                                if not n==' ' \
                                    ]
            # and not n=='-' and not n.isnumeric()]
        formatted_name = ''
        for char in formatted_name_chars:
            formatted_name+=char
        return formatted_name
    
    # %% Add chemical synonyms
    def set_formatted_common_and_iupac_names_as_IDs(test_env_chems):
        
        for chem in test_env_chems:
            formatted_common_name = ''
            if chem.common_name==None:
                chem.common_name = chem.CAS
            formatted_common_name = formatted_name(chem.common_name)
            
            try:
                test_env_chems[formatted_common_name]
                # print(f'\n{test_env_chems[formatted_common_name]} is already set for a chemical with CAS number {test_env_chems[formatted_common_name].CAS}.')
            except:
                test_env_chems.set_synonym(chem.ID, formatted_common_name)
            
            # test_env_chems.set_synonym(chem.ID, formatted_name(chem.common_name))
            
            # if not chem.iupac_name=='()':
            #     test_env_chems.set_synonym(chem.ID, chem.iupac_name)
    
    set_formatted_common_and_iupac_names_as_IDs(test_env_chems)
    
    
    # %% Streams initialization
  
    stream_modifiers = 'baseline_stream'
    process_stream_orig = tmo.Stream('process_stream_orig')
    process_stream_orig.copy_like(stream)
    
    # process_stream_orig.mix_from([u.S402.ins[0], u.S402.ins[1]])
    
    if stream_modifiers == 'baseline_stream':
        pass
    
    elif stream_modifiers == 'impurity_free_stream':
        for i_ID in impurity_IDs:
            process_stream_orig.imol[i_ID] = 0.
    
    elif stream_modifiers == 'solute_in_pure_water':
        process_stream_orig.empty()
        process_stream_orig.imol[solute_ID] = stream.imol[solute_ID]
        process_stream_orig.imol['Water'] = stream.outs[0].imol['Water']
    
    tmo.settings.set_thermo(test_env_chems)
    
    process_stream = tmo.Stream('process_stream',
                                T = T)
    
    for c in process_stream_orig.chemicals:
        try: 
            process_stream.imol[c.iupac_name] = process_stream_orig.imol[c.ID]
        except:
            try:
                process_stream.imol[c.common_name] = process_stream_orig.imol[c.ID]
            except:
                process_stream.imol[c.ID] = process_stream_orig.imol[c.ID]
    
    
    solvent_stream = tmo.Stream('solvent_stream',
                                T = T)
    
    mixed_stream = tmo.Stream('mixed_stream')
    extract_phase = 'l'
    raffinate_phase = 'L'
    
    # %%% Functions
    
    def get_K(chem_ID, stream, phase_1, phase_2):
        # try:
        return (stream[phase_1].imol[chem_ID]/stream[phase_1].F_mol)/max(1e-6, (stream[phase_2].imol[chem_ID]/stream[phase_2].F_mol))
        # except:
        #     import pdb
        #     pdb.set_trace()
    def set_solvent(solvent_chemical, solvent_mol=1314, solvent_stream=solvent_stream):
        solvent_stream.empty()
        if type(solvent_chemical) is str:
            solvent_chemical = test_env_chems[solvent_chemical]
        solvent_ID = formatted_name(solvent_chemical.common_name)
        solvent_stream.imol[solvent_ID] = solvent_mol
        
    def run_single_test(solvent_chemical, solvent_mol=1314, process_stream=process_stream, solvent_stream=solvent_stream):
        mixed_stream.empty()
        set_solvent(solvent_chemical, solvent_mol, solvent_stream)
        solvent_ID = formatted_name(solvent_chemical.common_name)
        mixed_stream.mix_from([process_stream, solvent_stream])
        # print(mixed_stream.lle_chemicals)
        mixed_stream.lle(T=T, top_chemical = solvent_chemical.ID)
        # mixed_stream.lle(T=T)
        # mixed_stream.show(N=100, composition=True)
        K_solute_in_solvent = get_K(solute.ID, mixed_stream, extract_phase, raffinate_phase)
        # K_Glucose_in_solvent = get_K('Glucose', mixed_stream, extract_phase, raffinate_phase)
        # print(K_solute_in_solvent)
        K_Water_in_solvent = get_K('Water', mixed_stream, extract_phase, raffinate_phase)
        # print(K_Water_in_solvent)
        K_solvent_in_Water = get_K(solvent_ID, mixed_stream, raffinate_phase, extract_phase)
        # print(K_solvent_in_Water)
        # K_Water_in_Water = get_K('Water', mixed_stream, raffinate_phase, extract_phase)
        
        
        # T_diff =  solute.Tb - solvent_chemical.Tb
        
        T_diff =   solvent_chemical.Tb - 100. # solvent - water
        
        
        # print(T_diff)
        # print('\n')
        K_impurity_1_in_solvent = 0.000001
        K_impurity_2_in_solvent = 0.000001
        try:
            K_impurity_1_in_solvent = get_K(impurity_IDs[0], mixed_stream, extract_phase, raffinate_phase)
        except:
            pass
        try:
            K_impurity_2_in_solvent = get_K(impurity_IDs[1], mixed_stream, extract_phase, raffinate_phase)
        except:
            pass
        return (mixed_stream.copy(), K_solute_in_solvent, K_Water_in_solvent, K_solvent_in_Water, T_diff, K_impurity_1_in_solvent, K_impurity_2_in_solvent,
                solvent_chemical.Tb-solute.Tb)
    
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
    
    def get_LHK(chem_IDs):
        chem1 = tmo.Chemical(chem_IDs[0])
        chem2 = tmo.Chemical(chem_IDs[1])
        if chem1.Tb < chem2.Tb:
            return (chem_IDs[0], chem_IDs[1])
        if chem1.Tb > chem2.Tb:
            return (chem_IDs[1], chem_IDs[0])
        else:
            raise RuntimeError('get_LHK: chemicals have the same Tb')
            
    def get_solute_solvent_distillation_steam(stream, LHK, Lr, Hr, k=1.2):
        DSS = bst.units.BinaryDistillation('DSS', ins=stream, outs=('DSS_g', 'DSS_l'),
                                            LHK=LHK,
                                            product_specification_format='Recovery',
                                            Lr=Lr, Hr=Hr, k=1.2, P = 101325.,
                                            vessel_material = 'Stainless steel 316')
        steam_flow = 0.
        try:
            DSS.simulate()
            steam_flow = sum ([i.flow for i in DSS.heat_utilities if i.duty>0])
        except: steam_flow = np.inf
        
        return steam_flow
    
    def forms_azeotrope_with_water(chem_ID):
        d_az_in = tmo.Stream('d_az_in')
        d_az_in.imol['Water'] = 100
        d_az_in.imol[chem_ID] = 100
        LHK = (chem_ID, 'Water') if test_env_chems[chem_ID].Tb<test_env_chems['Water'].Tb else ('Water', chem_ID)
        d_az = bst.BinaryDistillation('d_az', ins=d_az_in, LHK=LHK, k=1.2, Lr=0.99, Hr=0.99)
        try:
            d_az.simulate()
            return False
        except:
            return True
    # %% Run test barrage
    results_dict = {}
    solvent_mol = process_stream.ivol['Water']
    for solvent in solvents:
        results_dict[formatted_name(solvent.common_name)] = run_single_test(solvent_chemical=solvent, solvent_mol=solvent_mol)
    
    results_list = list(results_dict.items())
    # results_list is a list of tuples each of the following format:
    # (solvent_common_name, (mixed_stream, K_solute_in_solvent, K_Water_in_solvent, K_solvent_in_Water, T_diff))
    
    # %% Sort results by criterion
    key = None
    reverse = True
    
    # criterion = 'Partition of water into extract'
    # criterion = 'Partition of solvent into raffinate'
    # criterion = 'Difference between solute and solvent boiling temperatures'
    # criterion = 'Absolute difference between solute and solvent boiling temperatures'
    
    if criterion == 'Partition of solute into extract':
        key = lambda i: i[1][1]
        reverse = True
    elif criterion == 'Partition of water into extract':
        key = lambda i: i[1][2]
        reverse = True
    elif criterion == 'Partition of solvent into raffinate':
        key = lambda i: i[1][3]
        reverse = False
    elif criterion == 'Difference between solvent and solute boiling temperatures':
        key = lambda i: i[1][4]
        reverse = True
    elif criterion == 'Absolute difference between solvent and water boiling temperatures':
        key = lambda i: i[1][4]
        reverse = True
        
    elif criterion == 'Partition of impurity 1 into extract':
        key = lambda i: i[1][5]
    
    elif criterion == 'Partition of impurity 2 into extract':
        key = lambda i: i[1][6]
        reverse = False
        
    elif criterion == 'Ease of solute-solvent distillative separation':
        key = lambda i: get_solute_solvent_distillation_steam(stream=i[1][0],
                                                              LHK=get_LHK((solute_ID, i[0])), 
                                                              Lr=0.995, Hr=0.995, k=1.2, P=101325./10.)
        reverse = False
    # reverse = True
    
    results_list.sort(key=key, reverse=reverse)
    
    if print_result_with_optimal_criterion:
        print(f"\n\nResult with optimal {criterion.lower()}: \n")
        print(results_list[0])
        results_list[0][1][0].show()
    
    # %% Filter results through constraints
    
    inf = float('inf')
    
    range_K_solute_in_solvent = (20., inf)
    range_K_Water_in_solvent = (0., 0.2)
    range_K_solvent_in_Water = (0., 0.02)
    range_T_diff = (25, inf)
    # range_abs_T_diff = (0., inf)
    
    constraints = (range_K_solute_in_solvent, range_K_Water_in_solvent, range_K_solvent_in_Water,
                   range_T_diff)
    
    filtered_results = get_results_meeting_constraints(results_list, constraints)
    
    # print(f"\n\nResults meeting specified constraints ({len(filtered_results)} results out of {len(results_list)} total): \n")
    # for result in filtered_results:
    #     print(result)
    #     result[1][0].show()
    #     print('\n')
        
    # %% Unit initialization and tests
    solvent_to_run = results_list[0][0]
    set_solvent(solvent_to_run)
    partition_data = dict(IDs=(solute_ID, 'Water', solvent_to_run,
                                'GlucoseOligomer', 'SolubleLignin',), 
                K=np.array([1./results_dict[solvent_to_run][1],
                            1/results_dict[solvent_to_run][2],
                            results_dict[solvent_to_run][3],
                           results_dict[solvent_to_run][5],
                           results_dict[solvent_to_run][6]]),
                phi = 0.5)
    
    MS = bst.units.MultiStageMixerSettlers('MS', ins = (process_stream, solvent_stream),
                                          outs = ('raffinate', 'extract'),
                                          N_stages = 8, partition_data = partition_data,
                                          thermo=tmo.settings.get_thermo())
    
    # S404 = bst.units.MultiStageMixerSettlers('S404', ins = (F401_H-0, M401_H-0),
    #                                          outs = ('raffinate', 'extract'),
    #                                          N_stages = 15, partition_data = Kds,) 
    
    if show_all_mixer_settlers:
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
        rn['K_solvent_in_raffinate'] = result[1][3]
        rn[f'K_{impurity_IDs[0]}_in_extract'] = result[1][5]
        rn[f'K_{impurity_IDs[1]}_in_extract'] = result[1][6]
        rn['Tb_solvent - Tb_water'] = result[1][4]
        rn['Tb_solvent - Tb_solute'] = result[1][7]
        rn['Forms azeotrope with water'] = forms_azeotrope_with_water(result[0])
        
    results_df = pd.DataFrame(compiled_results_dict)
    
    dateTimeObj = datetime.now()
    file_to_save = 'solvents_barrage-' + str(T-273.15) + '_Celsius-' + stream_modifiers + '-%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)
    
    map_dict = {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H', 8: 'I', 
    9: 'J', 10: 'K', 11: 'L', 12: 'M', 13: 'N', 14: 'O', 15: 'P', 16: 'Q', 
    17: 'R', 18: 'S', 19: 'T', 20: 'U', 21: 'V', 22: 'W', 23: 'X', 24: 'Y', 25: 'Z'}
    
    final_results_df = None
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
    
    #%% Plot results
    # for i in range(len(final_results_df.columns)):
    #     final_results_df.sort_values(by=final_results_df.columns[i], inplace=True)
    #     ax = final_results_df.plot.barh( y=i, figsize=(40,30), fontsize=52., ylabel=final_results_df.columns[i])
    #     # ax.xlabel('xlabel', fontsize=18)
    #     plt.ylabel(final_results_df.columns[i], fontsize=52)
    #     fig = ax.get_figure()
    #     fig.savefig(file_to_save+'___'+str(final_results_df.columns[i])+'.png', bbox_inches='tight')
    
    
    final_results_df.sort_values(by=final_results_df.columns[0], inplace=True)
    ax1 = final_results_df.plot.barh( y=final_results_df.columns[0:3], figsize=(40,30), fontsize=52., ylabel=final_results_df.columns[i])
    # ax.xlabel('xlabel', fontsize=18)
    plt.xlabel('Partition coefficients at ' + str(T-273.15) + ' deg Celsius [(mol/mol)/(mol/mol)]', fontsize=52)
    plt.legend(loc='lower right', prop={'size': 52})
    fig1 = ax1.get_figure()
    fig1.savefig(file_to_save+'___'+ 'Ks-solute-solvent-water'+'.png', bbox_inches='tight')
    
    
    ax2 = final_results_df.plot.barh( y=final_results_df.columns[3:5], figsize=(40,30), fontsize=52., ylabel=final_results_df.columns[i])
    # ax.xlabel('xlabel', fontsize=18)
    plt.xlabel('Partition coefficients at ' + str(T-273.15) + ' deg Celsius [(mol/mol)/(mol/mol)]', fontsize=52)
    plt.legend(loc='lower right', prop={'size': 52})
    fig2 = ax2.get_figure()
    fig2.savefig(file_to_save+'___'+ 'Ks-impurities-solvent'+'.png', bbox_inches='tight')
    
    # %% Get rough solubility in organic solvent vs T