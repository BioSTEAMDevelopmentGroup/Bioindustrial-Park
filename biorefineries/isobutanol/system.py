#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import nskinetics as nsk
import biosteam as bst
import thermosteam as tmo
import numpy as np
from matplotlib import pyplot as plt
from biorefineries import corn
from biorefineries.isobutanol import units
from nskinetics.examples.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r, reset_kinetic_reaction_system
from scipy.optimize import differential_evolution

from warnings import filterwarnings
filterwarnings('ignore')

MultiEffectEvaporator = bst.MultiEffectEvaporator

__all__ = ('corn_EtOH_IBO_sys',)

corn_chems_compiled = corn.chemicals.create_chemicals()
chems = [c for c in corn_chems_compiled]
chems.append(tmo.Chemical('Isobutanol'))

solvent_chem = 'Isopentyl acetate' # 1
# solvent_chem = 'Valeraldehyde' # 2
# solvent_chem = '2-ethyl hexanol' # 3

chems.append(tmo.Chemical(solvent_chem))
chems.append(tmo.Chemical('AceticAcid'))
chems.append(tmo.Chemical('Acetaldehyde'))
tmo.settings.set_thermo(chems)

#%%
# corn.load(chemicals=chems)

# corn.system.simulate()
# corn.system.diagram('cluster')

settings = corn.process_settings.BiorefinerySettings()

#%%

corn_EtOH_sys = corn.systems.create_system(biorefinery_settings=settings)
corn_EtOH_sys.simulate()
f = corn_EtOH_sys.flowsheet

parameters = settings.process_parameters

#%% Add splitter and feed & spike evaporators, mixers, and heat exchangers

## Splitter
S301 = bst.Splitter('S301', ins=f.E402-0,
                    outs = ('fermentation_initial_feed', 'fermentation_spike'),
                    split = 0.8, # initial value, updated in FeedStrategySpecification object
                    )

## Initial feed evaporator, pumps, mixer, and hx
F301 = bst.MultiEffectEvaporator('F301', ins=S301-0, outs=('F301_l', 'F301_g'),
                                        P = (101325, 73581, 50892, 32777, 20000), V = 0.1,
                                        flash=False)
F301.V = 0.8 # initial value, updated in FeedStrategySpecification object
F301_design = F301._design
F301_cost = F301._cost

@F301.add_specification(run=False)
def F301_spec():
    feed = F301.ins[0]
    if feed.F_mol:
        # and feed.imass['Water']/feed.F_mass > 0.2:
        F301._run()
        F301._design = F301_design
        F301._cost = F301_cost
    else:
        F301.outs[1].empty()
        F301.outs[0].copy_like(feed)
        F301._design = lambda:0
        F301._cost = lambda:0

F301_P0 = bst.units.Pump('F301_P0', ins=F301-0, outs='', P=101325.)
F301_P1 = bst.units.Pump('F301_P1', ins=F301-1, outs='', P=101325.)

M301 = bst.units.Mixer('M301', ins=(F301_P0-0, 'dilution_water'))
M301.water_to_sugar_mol_ratio = 5. # initial value, updated in FeedStrategySpecification object

@M301.add_specification(run=False)
def adjust_M301_water():
    M301_ins_1 = M301.ins[1]
    M301_ins_1.imol['Water'] = M301.water_to_sugar_mol_ratio * M301.ins[0].imol[V406.sugar_IDs].sum()
    M301._run()
    
H301 = bst.units.HXutility('H301', ins=M301-0, outs=('glucose_initial_feed',), T=32+273.15, rigorous=True)

## Spike evaporator, pumps, mixer, and hx
F302 = bst.MultiEffectEvaporator('F302', ins=S301-1, outs=('F302_l', 'F302_g'),
                                        P = (101325, 73581, 50892, 32777, 20000), V = 0.1,
                                        flash=False)
F302.V = 0.8 # initial value, updated in FeedStrategySpecification object
F302_design = F302._design
F302_cost = F302._cost

@F302.add_specification(run=False)
def F302_spec():
    feed = F302.ins[0]
    if feed.F_mol:
    # and feed.imass['Water']/feed.F_mass > 0.2:
        F302._run()
        F302._design = F302_design
        F302._cost = F302_cost
    else:
        F302.outs[1].empty()
        F302.outs[0].copy_like(feed)
        F302._design = lambda:0
        F302._cost = lambda:0

def F302_design_extended():
    try:
        F302_design()
    except:
        F302._V_overall(F302._V_first_effect)
        F302_design()
        
F302_P0 = bst.units.Pump('F302_P0', ins=F302-0, outs='', P=101325.)
F302_P1 = bst.units.Pump('F302_P1', ins=F302-1, outs='', P=101325.)

M302 = bst.units.Mixer('M302', ins=(F302_P0-0, 'dilution_water'))
M302.water_to_sugar_mol_ratio = 5. # initial value

@M302.add_specification(run=False)
def adjust_M302_water():
    M302_ins_1 = M302.ins[1]
    M302_ins_1.imol['Water'] = M302.water_to_sugar_mol_ratio * M302.ins[0].imol[V406.sugar_IDs].sum()
    M302._run()
    
H302 = bst.units.HXutility('H302', ins=M302-0, outs=('glucose_spike_feed',), T=32+273.15, rigorous=True)

#%%
V405_old = f.V405

# V406 = units.SSFEtOHIBO(ID='V405', ins=(f.E402-0, f.P404-0), outs=('CO2', ''), V=1.9e3,
#                             kinetic_reaction_system=te_r,
#                             n_simulation_steps=1000,
#                             f_reset_kinetic_reaction_system=reset_kinetic_reaction_system,
#                             map_chemicals_nsk_to_bst = {'s_glu': 'Glucose',
#                                                         'x': 'Yeast',
#                                                         's_EtOH': 'Ethanol',
#                                                         's_IBO': 'Isobutanol'}
#                             )

V406 = nsk.units.NSKFermentation('V406', 
                                 ins=(H301-0, f.P404-0, H302-0), 
                                 kinetic_reaction_system=te_r,
                                 n_simulation_steps=2000,
                                 map_chemicals_nsk_to_bst = {'[s_glu]': 'Glucose',
                                                             '[x]': 'Yeast',
                                                             '[s_EtOH]': 'Ethanol',
                                                             '[s_IBO]': 'Isobutanol',
                                                             '[s_acetate]': 'AceticAcid',
                                                             # '[s_acetald]': 'Acetaldehyde',
                                                             },
                                 track_vars = ['y_EtOH_glu_added', 
                                               'y_EtOH_glu_consumed',
                                               'curr_n_glu_spikes',
                                               # 'tot_mass_glu', 
                                               'prod_EtOH',
                                               'curr_tot_vol_glu_feed_added',
                                               'curr_env',],
                                 f_reset_kinetic_reaction_system=reset_kinetic_reaction_system,
                                 tau=3*24,
                                 tau_max=3*24,
                                 sugar_IDs=('Glucose',),
                                 # tau_update_policy=None,
                                 tau_update_policy=('max', '[s_EtOH]'),
                                 perform_hydrolysis=False)

V406-0-1-f.V409
V406-1-0-f.P406

yeast = f.yeast
gluco_amylase = f.gluco_amylase
@V406.add_specification(run=True)
def correct_saccharification_feed_flows():
    mash = V406.ins[0]
    mash_flow = mash.F_mass
    mash_dry_flow = mash_flow - mash.imass['Water']
    yeast.F_mass = parameters['yeast_loading'] * mash_flow
    gluco_amylase.F_mass = parameters['saccharification_gluco_amylase_loading'] * mash_dry_flow

# V406.simulate()

#%%

f.S1.outs[0].disconnect_sink()

#%%
corn_EtOH_IBO_sys_no_IBO_recovery = bst.System.from_units('corn_EtOH_IBO_sys_no_IBO_recovery', 
                                          units = [i for i in corn_EtOH_sys.units 
                                                   if not (i.ID=='V405')]
                                                  + [S301,
                                                     F301, F301_P0, F301_P1, M301, H301,
                                                     F302, F302_P0, F302_P1, M302, H302,
                                                     V406])
corn_EtOH_IBO_sys_no_IBO_recovery.simulate()

#%% Add isobutanol recovery system - stage 1/2
stillage = f.stillage

S404 = bst.Splitter('S404', ins=stillage, outs=('to_IBO_recovery', 'direct_to_DDGS_recovery'), 
                    split=0.999)

makeup_isopentyl_acetate = tmo.Stream('makeup_isopentyl_acetate')

M401 = bst.Mixer('M401', ins=('', makeup_isopentyl_acetate), outs=('isopentyl_acetate_solvent'))

M401_design = M401._design
M401_cost = M401._cost

bypass_IBO_separation_conditions = [lambda: True] # if any return True, don't try to recover Isobutanol

@M401.add_specification(run=False)
def M401_adjust_makeup_solvent():
    if not np.any([i() for i in bypass_IBO_separation_conditions]):
        M401._design = M401_design
        M401._cost = M401_cost
        req = S401.mol_solvent_per_mol_carrier*S401.ins[0].imol['Water']
        recycle = M401.ins[1]
        makeup = M401.ins[0]
        makeup.imol[solvent_chem] = max(0, req-recycle.imol[solvent_chem])
        M401._run()
    else:
        M401._design = lambda: 0
        M401._cost = lambda: 0
        M401.ins[1].empty()
        M401._run()

# solvent_extraction_thermo = tmo.Thermo(chemicals=[i for i in chems if i.ID in ('Water', 'Isobutanol', solvent_chem)])

S401 = bst.MultiStageMixerSettlers('S401', 
                                    ins=(S404-0, M401.outs[0]), 
                                    outs=('S401_extract', 'S401_raffinate'), N_stages=5,
                                    top_chemical=solvent_chem,
                                    )
S401.mol_solvent_per_mol_carrier = 0.2

S401_design = S401._design
S401_cost = S401._cost

@S401.add_specification(run=False)
def S401_partial_chems():
    if not np.any([i() for i in bypass_IBO_separation_conditions]):
        S401._design = S401_design
        S401._cost = S401_cost
        M401.simulate()
        feed = S401.ins[0]
        raffinate = S401.outs[1]
        chems_included_lle = ('Water', 'Isobutanol', solvent_chem)
        chems_excluded_lle = {}
        for i in feed.chemicals:
            if i.ID not in chems_included_lle:
                chems_excluded_lle[i.ID] = feed.imol[i.ID]
                feed.imol[i.ID] = 0.0
        
        S401._run()
        
        for k, v in chems_excluded_lle.items():
            for stream in (feed, raffinate):
                stream.imol[k] = v
        
    else:
        S401._design = lambda: 0
        S401._cost = lambda: 0
        S401.ins[1].empty()
        S401.outs[0].empty() # extract
        S401.outs[1].copy_like(S401.ins[0]) # raffinate
    

# M401.simulate()
# S401.simulate()
# S401.show(N=100)

M402 = bst.Mixer('M402', ins=(S401.extract, ''),)

D401 = bst.BinaryDistillation('D401', ins=M402-0, outs=('D401_t', 'D401_b'), LHK=('Isobutanol', solvent_chem), 
                              Lr=0.9999, Hr=0.9999, 
                              k=1.2, P=101325.0,
                              partial_condenser=False)
D401-1-0-M401 # recycle

D401_design = D401._design
D401_cost = D401._cost
@D401.add_specification(run=False)
def D401_bypass_spec():
    if D401.ins[0].F_mol:
        D401._design = D401_design
        D401._cost = D401_cost
        D401._run()
    else:
        D401._design = lambda: 0
        D401._cost = lambda: 0
        
        S401.outs[0] = S401.ins[0].copy()
        
# M402.simulate()
# D401.simulate()
# D401.show(N=100)

D401_0_P = bst.Pump('D401_0_P', ins=D401-0, P=101325.)

# D401_0_P.simulate()

#%% Continue adding isobutanol recovery system - stage 2/2
stage_2_feed = D401_0_P-0

S402 = bst.units.MolecularSieve('S402', ins=stage_2_feed, 
                                # split=(1280.06/1383.85, 2165.14/13356.04),
                                split=(0.99, 2165.14/13356.04), # !!! water split assumed
                                order=('Water', 'Isobutanol'))

S403 = bst.units.Splitter('S403', ins=S402-0, outs=('S403_recycle', 'S403_purge'), split=0.0)

S403-0-1-M402

H401 = bst.HXutility('H401', ins=S402-1, T=273.15+25, rigorous=True)

#%% add HX to cool and reconnect to DDGS units

H402 = bst.HXutility('H402', ins=S401-1, T=360.15, rigorous=True)

M403 = bst.Mixer('M403', ins=(S404-1, H402-0), outs='mixed_stream_to_DDGS_recovery')

M403-0-0-f.V601

#%% Add storage for isobutanol product
V514 = bst.StorageTank('V514', ins=H401-0, outs=('isobutanol'), tau=7*24)

V514.isobutanol_price = 1.43

@V514.add_specification(run=False)
def V514_update_IBO_price():
    V514._run()
    ibo = V514.outs[0]
    if ibo.F_mol: ibo.price = V514.isobutanol_price * ibo.imass['Isobutanol']/ibo.F_mass

#%% Remove all existing HXprocess units

def reconnect_without_HXprocess_unit(HXprocess_unit):
    for i in [0,1]:
        instream = HXprocess_unit.ins[i]
        outstream = HXprocess_unit.outs[i]
        
        insource = instream.source
        insource_index = insource.outs.index(instream)
        
        outsink = outstream.sink
        outsink_index = outsink.ins.index(outstream)
        
        insource-insource_index-outsink_index-outsink
     
HXprocess = bst.units.HXprocess
HXprocess_units = []
for i in corn_EtOH_sys.units:
    if isinstance(i, HXprocess):
        reconnect_without_HXprocess_unit(i)
        HXprocess_units.append(i)
        
#%% Create corn to ethanol + isobutanol system
keep_non_rigorous = [f.HX101, f.HX500, f.HX501]
for i in corn_EtOH_IBO_sys_no_IBO_recovery.units + []:
    if isinstance(i, bst.HXutility) and not i in keep_non_rigorous: 
        i.rigorous = True
    
recovery_units = [S404,
                  M401, S401, M402, D401, D401_0_P, S402, 
                  S403, H401, 
                  M403, H402, 
                  V514]

HXN = bst.HeatExchangerNetwork('HXN1001', ignored=keep_non_rigorous)


corn_EtOH_IBO_sys = bst.System.from_units('corn_EtOH_IBO_sys', 
                                          units = [i for i in corn_EtOH_IBO_sys_no_IBO_recovery.units + recovery_units + [HXN]
                                                   if not i in HXprocess_units]
                                          )
corn_EtOH_IBO_sys.simulate()

#%% Set prices
f.isobutanol.price = 1.43
f.makeup_isopentyl_acetate.price = 3.2 # https://www.alibaba.com/product-detail/High-Quality-Colorless-Liquid-99-min_1600206242747.html?spm=a2700.galleryofferlist.normal_offer.d_price.2ed613a0wyq5n8

#%% Create TEA object

corn_EtOH_IBO_sys._TEA = corn_EtOH_IBO_sys_tea = corn.tea.create_tea(corn_EtOH_IBO_sys)

#%% Set baseline specifications

baseline_spec = {'target_conc_sugars': 100,
                 'threshold_conc_sugars': 10,
                 'conc_sugars_feed_spike': 600,
                 'tau_max': 10*24,
                 'max_n_glu_spikes': 10,}

#% Create fed-batch strategy specification object
fbs_spec = nsk.units.FedBatchStrategySpecification(
    target_conc_sugars=100,
    threshold_conc_sugars=10,
    conc_sugars_feed_spike=600,
    tau_max=72,
    max_n_glu_spikes=10,
    fermentation_reactor=V406,
    splitter=S301,
    feed_evaporator=F301,
    feed_mixer=M301,
    feed_units_sequential=[F301, F301_P0, F301_P1, M301, H301],
    spike_units_sequential=[F302, F302_P0, F302_P1, M302, H302],
    spike_evaporator=F302,
    spike_mixer=M302,
    sugar_IDs=['Glucose',],
    baseline_specifications=baseline_spec,
    )

#%%

def get_purity_adj_price(stream, chem_IDs):
    return stream.price * stream.F_mass/sum([stream.imass[ID] for ID in chem_IDs])

def load_simulate_get_EtOH_MPSP(target_conc_sugars=None,
    threshold_conc_sugars=None,
    conc_sugars_feed_spike=None,
    tau_max=None,
    max_n_glu_spikes=None,
    n_sims=3,
    n_tea_solves=3,
    plot=False,
    ):
    
    if target_conc_sugars is None:
        target_conc_sugars = fbs_spec.target_conc_sugars
    
    if threshold_conc_sugars is None:
        threshold_conc_sugars = fbs_spec.threshold_conc_sugars
    
    if conc_sugars_feed_spike is None:
        conc_sugars_feed_spike = fbs_spec.conc_sugars_feed_spike
    
    if tau_max is None:
        tau_max = fbs_spec.tau_max
    
    if max_n_glu_spikes is None:
        max_n_glu_spikes = fbs_spec.max_n_glu_spikes
        
    ethanol = f.ethanol
    
    for i in range(n_sims):
        fbs_spec.load_specifications(target_conc_sugars=target_conc_sugars,
        threshold_conc_sugars=threshold_conc_sugars,
        conc_sugars_feed_spike=conc_sugars_feed_spike,
        tau_max=tau_max,
        max_n_glu_spikes=max_n_glu_spikes)
        
        corn_EtOH_IBO_sys.simulate()
        
    for i in range(n_tea_solves):
        ethanol.price = corn_EtOH_IBO_sys_tea.solve_price(ethanol)

    if plot:
        plot_results()
    
    return get_purity_adj_price(ethanol, ['Ethanol'])

def plot_results():
    # for i in V406.map_chemicals_nsk_to_bst.keys():
    plt.plot(V406.results_dict['time'], V406.results_dict['[x]'], label='cell mass')
    plt.plot(V406.results_dict['time'], V406.results_dict['[s_glu]'], label='glucose')
    plt.plot(V406.results_dict['time'], V406.results_dict['[s_EtOH]'], label='ethanol')
    plt.plot(V406.results_dict['time'], V406.results_dict['[s_acetate]'], label='acetate')
    plt.plot(V406.results_dict['time'], V406.results_dict['[s_IBO]'], label='isobutanol')
    plt.legend()
    plt.xlabel('Time [h]')
    plt.ylabel('Concentration [g/L]')
    plt.show()

def reset_and_reload(**curr_spec):
    print('Resetting cache and emptying recycles ...')
    corn_EtOH_IBO_sys.reset_cache()
    corn_EtOH_IBO_sys.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    # curr_spec = {i: fbs_spec.__getattribute__(i) for i in baseline_spec.keys()}
    corn_EtOH_IBO_sys.simulate()
    load_simulate_get_EtOH_MPSP(**fbs_spec.baseline_specifications)
    print('Loading and simulating with required specifications ...')
    # load_simulate_get_EtOH_MPSP(**curr_spec)
    load_simulate_get_EtOH_MPSP(**curr_spec)
    
def reset_and_switch_solver(solver_ID, **curr_spec):
    corn_EtOH_IBO_sys.reset_cache()
    corn_EtOH_IBO_sys.empty_recycles()
    corn_EtOH_IBO_sys.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    corn_EtOH_IBO_sys.simulate()
    load_simulate_get_EtOH_MPSP(**curr_spec)

# F403 = u.F403
def run_bugfix_barrage(**curr_spec):
    try:
        reset_and_reload(**curr_spec)
    except Exception as e:
        print(str(e))
        if 'length' in str(e).lower():
            try:
                corn_EtOH_IBO_sys.reset_cache()
                corn_EtOH_IBO_sys.empty_recycles()
                corn_EtOH_IBO_sys.simulate()
                load_simulate_get_EtOH_MPSP()
            except:
                print(str(e))
                raise e
        
        # elif 'invalid value encountered' in str(e).lower():
        #     print('\n\n\n\n\n\n\n\nAAAAAAAAAAAAAAAA\n\n\n\n\n\n\n')
        #     try:
        #         for i in corn_EtOH_IBO_sys.units:
        #             if isinstance(i, MultiEffectEvaporator):
        #                 for j in i.evaporators:
        #                     try:
        #                         j.outs[0].T = j.T
        #                         j.outs[1].T = j.T
        #                     except:
        #                         pass
        #         load_simulate_get_EtOH_MPSP()
                
        #     except:
        #         print(str(e))
        #         raise e
                
        else:
            try:
                reset_and_switch_solver('fixedpoint', **curr_spec)
            except Exception as e:
                print(str(e))
                try:
                    reset_and_switch_solver('aitken', **curr_spec)
                except Exception as e:
                    print(str(e))
                    print("Bugfix barrage failed.\n")
                    # breakpoint()
                    raise e


#%%
def model_specification(**kwargs):
    curr_spec = {k: v for k,v in fbs_spec.current_specifications.items()}
    curr_spec.update(kwargs)
    try:
        load_simulate_get_EtOH_MPSP(**curr_spec)
    except Exception as e:
        str_e = str(e).lower()
        print('Error in model spec: %s'%str_e)
        # raise e
        if 'specifications do not meet required condition' in str_e:
            # flowsheet('AcrylicAcid').F_mass /= 1000.
            raise e
        else:
            run_bugfix_barrage(**curr_spec)


def optimize_tau_for_MPSP(threshold_s_EtOH=5, **kwargs):
    original_run_type = V406.run_type
    V406.run_type = 'index saved results by tau'
    results = V406.results_dict
    where_greq_threshold = np.where(V406.results_dict['[s_EtOH]']>=5)[0]
    taus = results['time']
    bounds_tau = (taus[where_greq_threshold[0]], taus[where_greq_threshold[-1]])
    def f(x):
        V406.tau = x[0]
        try:
            # corn_EtOH_IBO_sys.simulate()
            model_specification(**kwargs)
            return get_purity_adj_price(ethanol, ['Ethanol'])
        except:
            return np.inf
    res = differential_evolution(f, bounds=(bounds_tau,), atol=1e-2)
    V406.run_type = original_run_type
    return res.x[0]
    

def optimize_max_n_glu_spikes_for_MPSP(bounds=(0, 10), 
                                       optimize_tau=True, 
                                       show_progress=False, 
                                       threshold_s_EtOH=5):
    curr_spec = {k: v for k,v in fbs_spec.baseline_specifications.items()}
    original_run_type = V406.run_type
    V406.run_type = 'simulate kinetics'
    r = V406.kinetic_reaction_system._te
    results = []
    for max_n_glu_spikes in range(bounds[0], bounds[1]+1):
        curr_curr_spec = {k:v for k,v in curr_spec.items()}
        curr_curr_spec.update({'max_n_glu_spikes':max_n_glu_spikes,})
        model_specification(**curr_curr_spec)
        if optimize_tau: 
            tau = optimize_tau_for_MPSP(threshold_s_EtOH=threshold_s_EtOH)
        opt_MPSP = get_purity_adj_price(ethanol, ['Ethanol'])
        results.append((fbs_spec.current_specifications, 
                        V406.tau,
                        V406.results_specific_tau_dict['curr_n_glu_spikes'], # num spikes at selected tau
                        V406.results_dict['curr_n_glu_spikes'][-1], # num spikes at tau_max
                        opt_MPSP))
        if show_progress: print(results[-1])
        if len(results)>1 and results[-1][3]==results[-2][3]: # if this iteration resulted in the same number of glucose spikes at tau_max as the last one (i.e., increasing the max number won't make a difference), break
            break
        
    results.sort(key=lambda x: x[4], reverse=False) # sort by MPSP
    opt_curr_spec, opt_tau, opt_curr_n_glu_spikes, opt_curr_n_glu_spikes_at_tau_max, opt_MPSP = results[0]
    V406.run_type = original_run_type
    model_specification(**opt_curr_spec)
    if show_progress: 
        print('optimum :', results[0])
        plt.scatter([i[0]['max_n_glu_spikes'] for i in results], [i[4] for i in results])
    return opt_curr_spec['max_n_glu_spikes']


#%% Initialize 
r = V406.kinetic_reaction_system._te
corn_EtOH_IBO_sys.simulate()

#%% Baseline -- simulate and solve TEA

ethanol = f.ethanol

simulate_baseline = True
if simulate_baseline:
    model_specification(**fbs_spec.baseline_specifications,
        n_sims=3,
        n_tea_solves=3,
        plot=True,
        )
    # print(get_purity_adj_price(ethanol, ['Ethanol']))
    
# optimize_max_n_glu_spikes_for_MPSP(baseline_spec)