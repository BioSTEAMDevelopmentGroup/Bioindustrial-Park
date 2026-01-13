# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 18:10:51 2026

@author: saran
"""

import numpy as np
import nskinetics as nsk
import biosteam as bst
from biorefineries import isobutanol
from matplotlib import pyplot as plt

system = corn_EtOH_IBO_sys = isobutanol.system.corn_EtOH_IBO_sys
load_simulate_get_EtOH_MPSP = isobutanol.system.load_simulate_get_EtOH_MPSP
model_specification = isobutanol.system.model_specification
baseline_spec = isobutanol.system.baseline_spec
tea = corn_EtOH_IBO_sys_tea = isobutanol.system.corn_EtOH_IBO_sys_tea
fbs_spec = isobutanol.system.fbs_spec
optimize_tau_for_MPSP = isobutanol.system.optimize_tau_for_MPSP
optimize_max_n_glu_spikes_for_MPSP = isobutanol.system.optimize_max_n_glu_spikes_for_MPSP

f = system.flowsheet
V406 = f.V406

# %%
# conc_sugars_feed_spikes = np.linspace(150, 600, 20)
# MPSPs = []
# ethanol = f.ethanol
# for c in conc_sugars_feed_spikes:
#     curr_spec = {k:v for k,v in baseline_spec.items()}
#     curr_spec.update({'conc_sugars_feed_spike':c,})
#     model_specification(
#     **curr_spec,
#     n_sims=3,
#     n_tea_solves=3,
#     plot=True,
#     )
#     MPSPs.append(ethanol.price * ethanol.F_mass/ethanol.imass['Ethanol'])

# plt.plot(conc_sugars_feed_spikes, MPSPs)
# plt.xlabel('Glucose spike feed concentration [g/L]')
# plt.ylabel('MPSP [$/kg]')
# plt.show()

# %%
target_conc_sugarses = np.linspace(100, 500, 10)
MPSPs = []
taus = []
max_yields = []
final_n_glu_spikes = []
titers = []
ethanol = f.ethanol
for c in target_conc_sugarses:
    curr_spec = {k:v for k,v in baseline_spec.items()}
    curr_spec.update({'target_conc_sugars':c,})
    model_specification(
    **curr_spec,
    n_sims=3,
    n_tea_solves=3,
    plot=True,
    )
    # optimize_tau_for_MPSP(**curr_spec)
    optimize_max_n_glu_spikes_for_MPSP(curr_spec=curr_spec, optimize_tau=True)
    MPSPs.append(ethanol.price * ethanol.F_mass/ethanol.imass['Ethanol'])
    taus.append(V406.tau)
    max_yields.append(V406.results_specific_tau_dict['y_EtOH_glu_added'])
    final_n_glu_spikes.append(V406.results_specific_tau_dict['curr_n_glu_spikes'])
    titers.append(V406.results_specific_tau_dict['[s_EtOH]'])
    # plt.plot(V406.results_dict['time'], V406.results_dict['y_EtOH_glu_added'])
    # plt.xlabel(f'Time [h]; target conc={c}')
    # plt.ylabel('Yield [g Ethanol produced / g Glucose added]')
    # plt.show()
    print(c, MPSPs[-1], taus[-1], max_yields[-1])
    # try:
    #     if MPSPs[-1] > MPSPs[-2]:
    #         breakpoint()
    # except:
    #     pass

plt.plot(target_conc_sugarses, MPSPs)
plt.xlabel('Target glucose concentration [g/L]')
plt.ylabel('MPSP [$/kg]')
plt.show()

plt.plot(target_conc_sugarses, taus)
plt.xlabel('Target glucose concentration [g/L]')
plt.ylabel('tau [h]')
plt.show()

plt.plot(target_conc_sugarses, max_yields)
plt.xlabel('Target glucose concentration [g/L]')
plt.ylabel('max yield [g Ethanol produced / g Glucose added]')
plt.show()


plt.plot(target_conc_sugarses, final_n_glu_spikes)
plt.xlabel('Target glucose concentration [g/L]')
plt.ylabel('final n glu spikes [-]')
plt.show()

plt.plot(target_conc_sugarses, titers)
plt.xlabel('Target glucose concentration [g/L]')
plt.ylabel('Titer [g/L]')
plt.show()

# %%
threshold_conc_sugarses = np.linspace(0, 80, 20)
MPSPs = []
taus = []
max_yields = []
final_n_glu_spikes = []
ethanol = f.ethanol
for c in threshold_conc_sugarses:
    curr_spec = {k:v for k,v in baseline_spec.items()}
    curr_spec.update({'threshold_conc_sugars':c,})
    model_specification(
    **curr_spec,
    n_sims=3,
    n_tea_solves=3,
    plot=False,
    )
    MPSPs.append(ethanol.price * ethanol.F_mass/ethanol.imass['Ethanol'])
    taus.append(V406.tau)
    max_yields.append(V406.results_specific_tau_dict['y_EtOH_glu_added'])
    final_n_glu_spikes.append(V406.results_specific_tau_dict['curr_n_glu_spikes'])
    plt.plot(V406.results_dict['time'], V406.results_dict['y_EtOH_glu_added'])
    plt.xlabel(f'Time [h]; target conc={c}')
    plt.ylabel('Yield [g Ethanol produced / g Glucose added]')
    plt.show()
    print(c, MPSPs[-1], taus[-1], max_yields[-1])
    # try:
    #     if MPSPs[-1] > MPSPs[-2]:
    #         breakpoint()
    # except:
    #     pass

plt.plot(threshold_conc_sugarses, MPSPs)
plt.xlabel('Threshold glucose concentration [g/L]')
plt.ylabel('MPSP [$/kg]')
plt.show()

plt.plot(threshold_conc_sugarses, taus)
plt.xlabel('Threshold glucose concentration [g/L]')
plt.ylabel('tau [h]')
plt.show()

plt.plot(threshold_conc_sugarses, max_yields)
plt.xlabel('Threshold glucose concentration [g/L]')
plt.ylabel('max yield [g Ethanol produced / g Glucose added]')
plt.show()


plt.plot(threshold_conc_sugarses, final_n_glu_spikes)
plt.xlabel('Threshold glucose concentration [g/L]')
plt.ylabel('final n glu spikes [-]')
plt.show()


# %%
taus =  np.linspace(5, 200, 10)
MPSPs = []
max_yields = []
final_n_glu_spikes = []
ethanol = f.ethanol
V406.run_type = 'index saved results by tau'
for t in taus:
    curr_spec = {k:v for k,v in baseline_spec.items()}
    # curr_spec.update({'threshold_conc_sugars':t,})
    V406.tau = t
    model_specification(
    **curr_spec,
    n_sims=3,
    n_tea_solves=3,
    plot=False,
    )
    MPSPs.append(ethanol.price * ethanol.F_mass/ethanol.imass['Ethanol'])
    max_yields.append(V406.results_specific_tau_dict['y_EtOH_glu_added'])
    final_n_glu_spikes.append(V406.results_specific_tau_dict['curr_n_glu_spikes'])
    # plt.plot(V406.results_dict['time'], V406.results_dict['y_EtOH_glu_added'])
    # plt.xlabel(f'Time [h]; tau={t}')
    # plt.ylabel('Yield [g Ethanol produced / g Glucose added]')
    # plt.show()
    print(t, MPSPs[-1], taus[-1], max_yields[-1], V406.outs[1].imass['Ethanol']/(sum([i.imass['Glucose'] for i in V406.ins])))
    # try:
    #     if MPSPs[-1] > MPSPs[-2]:
    #         breakpoint()
    # except:
    #     pass

V406.run_type = 'simulate kinetics'

plt.plot(taus, MPSPs)
plt.xlabel('Tau [h]')
plt.ylabel('MPSP [$/kg]')
plt.show()


plt.plot(taus, max_yields)
plt.xlabel('Tau [h]')
plt.ylabel('max yield [g Ethanol produced / g Glucose added]')
plt.show()


plt.plot(taus, final_n_glu_spikes)
plt.xlabel('Tau [h]')
plt.ylabel('final n glu spikes [-]')
plt.show()
