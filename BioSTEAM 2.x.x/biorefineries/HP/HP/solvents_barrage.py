# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 19:09:17 2021

@author: sarangbhagwat
"""
# %% Imports and chemicals initialization

import numpy as np
# import thermosteam as tmo
# import biosteam as bst

from biorefineries.HP.chemicals_data import HP_chemicals

# tmo.settings.set_thermo(['Water', 'octanol', 'hexanol', 'butyl acetate', HP_chemicals['Xylose'], HP_chemicals['Glucose'], HP_chemicals['Triacetic acid lactone'], 'isoamyl alcohol'])

class Metric():
    
    def __init__(self, metric_name, metric_getter=lambda:None):
        self.name = metric_name
        self.getter = metric_getter
    
    def get():
        return self.getter()

# class Metrics():
    
#     def __init__(self, name='Metrics', metrics_dict={}):
#         self.name = name,
#         self.metrics = metrics_dict
    
#     def add_metric(metric_name, metric_getter=lambda:None):
#         self.metrics[metric_name] = Metric(metric_name, metric_getter)
#         print (f'Added metric "{metric_name}" to {self.name}.')

class Criterion():
    
    def __init__(self, criterion_name, criterion_metric=None, criterion_range=(-np.inf, np.inf)):
        self.name = criterion_name
        self.metric = criterion_metric
        self.range = criterion_range
        
    def get_metric():
        return self.metric.get()
    
    def test():
        gotten_metric = self.get_metric()
        lb, ub = self.range
        return gotten_metric>lb and gotten_metric<ub
    
class Criteria():
    
    def __init__(self, name='Criteria', criteria_list=[]):
        self.name = name,
        self.criteria = []
    
    def add_criterion(criterion_name, criterion_metric=None, criterion_range=(-np.inf, np.inf)):
        self.criteria.append(Criterion(criterion_name, criterion_metric, criterion_range))
        print (f'Added criterion to {self.name} for "{criterion_metric.name}" to be within range {criterion_range}.')
    
    def test_all():
        criteria = self.criteria
        return [criterion.test() for criterion in criteria]
    
class SolventsBarrage():
    

    def __init__(self, solute, water, solvents_list, T_range, process_stream,
                 extract_phase, raffinate_phase,
                 additional_criteria=[],
                 get_concentration=lambda:None,
                 get_K=lambda:None,
                 get_name_from_chemical=lambda:None,
                 get_chemical_from_name=lambda:None,
                 get_all_chemicals_in_stream=lambda:None):
        self.solute = solute
        self.solute_name = solute_name = get_name_from_chemical(solute)
        self.water = water
        self.water_name = water_name = get_name_from_chemical(water)
        
        self.solvents_list = solvents_list
        self.T_range = T_range
        self.extract_phase = extract_phase
        self.raffinate_phase = raffinate_phase
        
        self.partition_criteria = partition_criteria = Criteria('Partition criteria')
        
        get_K_solute_in_extract = lambda: get_K(self.solute_name, self.mixed_stream, 
                                                self.extract_phase, self.raffinate_phase)
        
        partition_criteria.add_criterion('K_solute_in_extract', get_K(solute, ))
        self.process_stream = None
        self.solvent_streams_dict = {}
        self.mixed_streams_dict = {}
        self.current_mixed_stream = None
        self.results = {}
        
        self.get_concentration = get_concentration
        self.get_K = get_K
        self.get_name_from_chemical = get_name_from_chemical
        self.get_chemical_from_name = get_chemical_from_name
        self.get_all_chemicals_in_stream = get_all_chemicals_in_stream


        solute_name = get_name_from_chemical(self.solute)
        water_name = get_name_from_chemical(self.water)

    def get_impurities(stream):
        impurities_dict = {}
        chemicals_in_stream = get_all_chemicals_in_stream(stream)
        solute_name = self.solute_name
        water_name = self.water_name
        
        for chemical in chemicals_in_stream:
            chemical_name = get_name_from_chemical(chemical)
            if not chemical_name==water_name and not chemical_name==solute_name:
                chemical_concentration = get_concentration(chemical, stream)
                if chemical_concentration>0.:
                    impurities_dict[chemical_name] = chemical_concentration
        
        return impurities_dict
    
        # return {get_name_from_chemical(chemical): get_concentration(chemical) \
        #         for chemical in stream.chemicals \
        #         if not (get_name_from_chemical(chemical)==solute \
        #                 and not get_name_from_chemical(chemical)==get_name_from_chemical(water)) \
        #             and get_concentration(get_name_from_chemical(chemical))>0.}
            
    def set_process_stream(process_stream):
        self.process_stream = process_stream
        self.process_stream_impurities = get_impurities(process_stream)
        
    def set_solvent_stream(solvent_name, solvent_mol, water_mol, solvent_stream=None):
        if solvent_stream:
            self.solvent_stream = solvent_stream
        else:
            
    def get_results(mixed_stream, )
    def run_single_test(solvent_chemical, solvent_mol, T):
        process_stream_impurities = self.process_stream_impurities
        solvent_stream = get_solvent_stream(solvent_chemical, solvent_mol)
        mixed_stream = get_mixed_stream(self.process_stream, solvent_stream, T)
        extract_phase = self.extract_phase
        raffinate_phase = self.raffinate_phase
        
        # Ks = [get_K(solute, mixed_stream, extract_phase, raffinate_phase)] \
        #     + [get_K(chemical, stream, phase_1, phase_2) for chemical in process_stream_impurities.keys()]
#%% Set solute and solvents
solute = TAL
solvents = [Propyl_acetate, Butyl_acetate, Hexanol, Cyclohexanol, Cyclohexanone, Heptanol, Octanol, Octanediol, te_hexanol, Nonanol, Decanol, Dodecanol, Isoamyl_alcohol, Dioctyl_phthalate, Diethyl_sebacate, Glycerol]

#%% Set thermo
tmo.settings.set_thermo(solvents + ['Water', 'H2SO4', Arabitol, solute, HP_chemicals['Xylose'], HP_chemicals['Glucose'], HP_chemicals['AceticAcid']])

# %% Streams initialization

T = 303
process_stream = tmo.Stream('process_stream',
                            Water = 224000.,
                            units = 'kmol/hr',
                            T = T)
process_stream.imol[solute.ID] = 250.*12.2/7.7
process_stream.imol['AceticAcid'] = 5.
process_stream.imol['Arabitol'] = 5.

solvent_stream = tmo.Stream('solvent_stream',
                            T = T)

mixed_stream = tmo.Stream('mixed_stream')
extract_phase = 'l'
raffinate_phase = 'L'

# %%% Functions

def get_K(chem_ID, stream, phase_1, phase_2):
    return (stream[phase_1].imol[chem_ID]/stream[phase_1].F_mol)/(stream[phase_2].imol[chem_ID]/stream[phase_2].F_mol)

def set_solvent(solvent_chemical, solvent_mol=1350, solvent_stream=solvent_stream):
    solvent_stream.empty()
    if type(solvent_chemical) is str:
        solvent_chemical = tmo.Chemical(solvent_chemical)
    solvent_ID = solvent_chemical.ID
    solvent_stream.imol[solvent_ID] = solvent_mol
    
def run_single_test(solvent_chemical, solvent_mol=1350, process_stream=process_stream, solvent_stream=solvent_stream):
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
    K_arabitol_in_solvent = 0.000001
    try:
        K_acetic_acid_in_solvent = get_K('AceticAcid', mixed_stream, raffinate_phase, extract_phase)
    except:
        pass
    try:
        K_arabitol_in_solvent = get_K('Arabitol', mixed_stream, raffinate_phase, extract_phase)
    except:
        pass
    return (mixed_stream.copy(), K_HP_in_solvent, K_Water_in_solvent, K_solvent_in_Water, T_diff, K_acetic_acid_in_solvent, K_arabitol_in_solvent)

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
partition_data = dict(IDs=('Triacetic acid lactone', 'Water', solvent_to_run,
                            'AceticAcid', 'Arabitol',), 
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