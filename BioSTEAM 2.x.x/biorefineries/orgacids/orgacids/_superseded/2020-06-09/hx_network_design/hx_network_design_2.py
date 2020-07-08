# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 13:47:23 2020

@author: saran
"""

# -*- coding: utf-8 -*-
"""
Created on Sat May  2 16:44:24 2020

@author: sarangbhagwat
"""
import numpy as np
# import copy
import biosteam as bst
from stream_life_cycle import HotStreamNetwork, ColdStreamNetwork

# def get_TinTout()
# def get_streams(hx_utils):
    
#     hxs = [hx_util.heat_exchanger for hx_util in hx_util]
#     streams = [[hx.ins[0], hx.outs[0]] for hx in hxs]
#     # T_vector = [(hx_uti)]
    
def temperature_interval_pinch_analysis(hus, T_min_app = 10, find = None):
    
    
    # hx_utils = [hu for hu in hus if abs(hu.T_in - hu.T_out)>0.01]
    hx_utils = [hu for hu in hus if abs(hu.heat_exchanger.ins[0].T - hu.heat_exchanger.outs[0].T)>0.01]
    hus_heating = [hu for hu in hx_utils if hu.duty > 0]
    hus_cooling = [hu for hu in hx_utils if hu.duty < 0]
    hxs_heating = [hu.heat_exchanger for hu in hx_utils if hu.duty > 0]
    hxs_cooling = [hu.heat_exchanger for hu in hx_utils if hu.duty < 0]
    
    streams_heating = [hx.ins[0].copy() for hx in hxs_heating]
    streams_cooling = [hx.ins[0].copy() for hx in hxs_cooling]
    
    
    init_heat_util, init_cool_util = 0,0
    init_Q = 0
    
    #!!! change all hx.ins[0] refs to streams refs
    hx_utils_rearranged = hus_heating + hus_cooling
    hxs = hxs_heating + hxs_cooling
    streams = streams_heating + streams_cooling
    
    for stream in streams:
        stream.vle(H = stream.H, P = stream.P)
        
    len_hxs_heating = len(hxs_heating)
    is_cold_stream_index = lambda x: x<len_hxs_heating
    T_in_arr = np.array([stream.T for stream in streams])
    T_out_arr = np.array([hx.outs[0].T for hx in hxs])    
    
    # hot_stream_indices, cold_stream_indices = [], []

            
    T_hot_side_arr = np.array([stream.T for stream in streams_heating] + [hx.outs[0].T for hx in hxs_cooling])
    T_cold_side_arr = np.array([hx.outs[0].T for hx in hxs_heating] + [stream.T for stream in streams_cooling])
    
    
    adj_T_in_arr = T_in_arr.copy()
    adj_T_in_arr[:len(hxs_heating)] -= T_min_app
    
    adj_T_out_arr = T_out_arr.copy()
    adj_T_out_arr[:len(hxs_heating)] -= T_min_app
    
    T_changes_tuples = list(zip(adj_T_in_arr, adj_T_out_arr))
    
    all_Ts_descending = [*adj_T_in_arr, *adj_T_out_arr]
    all_Ts_descending.sort(reverse=True)
    
    stream_indices_for_T_intervals = {
        (all_Ts_descending[i], all_Ts_descending[i+1]):[] for i in range(len(all_Ts_descending)-1)
    }
    

    H_for_T_intervals = dict.fromkeys(stream_indices_for_T_intervals, 0)
    
    cold_indices = list(range(len_hxs_heating))
    hot_indices = list(range(len_hxs_heating, len(hxs)))

        
    for i in range(len(all_Ts_descending)-1):
        T_start = all_Ts_descending[i]
        T_end = all_Ts_descending[i+1]
        dT = T_start - T_end
        
        for stream_index in range(len(T_changes_tuples)):
            # print('Stream %s'%stream_index)
            if (T_changes_tuples[stream_index][0]>= T_start and T_changes_tuples[stream_index][1]<= T_end)\
                or (T_changes_tuples[stream_index][1]>= T_start and T_changes_tuples[stream_index][0]<= T_end):
        
        
                stream_indices_for_T_intervals[(T_start, T_end)].append(stream_index)
                multiplier = 1
                if is_cold_stream_index(stream_index):
                    multiplier = -1
                    
                    
                stream = streams[stream_index].copy()
                
                stream.vle(T = T_start, P = stream.P)
                H1 = stream.H

                stream.vle(T = T_end, P = stream.P)
                H2 = stream.H
                
                H = multiplier*(H1 - H2)
                # print(H)
                H_for_T_intervals[(T_start, T_end)] += H
                    
                
    res_H_vector = []
    
    prev_res_H = 0
    
    for interval, H in H_for_T_intervals.items():
        res_H_vector.append(prev_res_H + H)
        prev_res_H = res_H_vector[len(res_H_vector)-1]
        
        
    # assert not res_H_vector == []
    hot_util_load = - min(res_H_vector)
    
    assert hot_util_load>= 0
    # print(res_H_vector)
    # print(all_Ts_descending)
    pinch_cold_stream_T = all_Ts_descending[res_H_vector.index(-hot_util_load)+1] # the lower temperature of the temperature interval for which the res_H is minimum
    print(pinch_cold_stream_T)
    pinch_hot_stream_T = pinch_cold_stream_T + T_min_app
    cold_util_load = res_H_vector[len(res_H_vector)-1] + hot_util_load
    
    # assert cold_util_load>=0
    
    pinch_T_arr = []
    for i in range(len(T_in_arr)):
        if not is_cold_stream_index(i):
            if T_in_arr[i]<pinch_hot_stream_T:
                pinch_T_arr.append(T_in_arr[i])
            elif T_out_arr[i]>pinch_hot_stream_T:
                pinch_T_arr.append(T_out_arr[i])
            else:
                
                pinch_T_arr.append(pinch_hot_stream_T)
                
        else:
            if T_in_arr[i]>pinch_cold_stream_T:
                pinch_T_arr.append(T_in_arr[i])
            elif T_out_arr[i]<pinch_cold_stream_T:
                pinch_T_arr.append(T_out_arr[i])
            else:
                pinch_T_arr.append(pinch_cold_stream_T)
                
    pinch_T_arr = np.array(pinch_T_arr) 
    # print('HEEEEEEEEEEEERE')
    # print(pinch_T_arr)
    # print(hot_util_load, cold_util_load)
    return(pinch_T_arr, hot_util_load, cold_util_load, T_in_arr, T_out_arr, T_hot_side_arr, T_cold_side_arr, \
           hus_heating, hus_cooling, hxs_heating, hxs_cooling, hxs, hot_indices, cold_indices, streams, hx_utils_rearranged)
    
        
def load_duties(streams, pinch_T_arr, T_out_arr, indices, is_cold, Q_hot_side, Q_cold_side):
    
    # T_transient = list(pinch_T_arr)
    # T_transient[index] = T_in_arr[index]
    for index in indices:
        # need to account for streams that don't need to be at pinch
        # also need to use SLC's T_transient
        stream = streams[index].copy()
        # stream.vle(H = stream.H, P = stream.P)
        stream_in = stream.copy()
        H_in = stream.H
        T_in = stream.T
        stream.T = pinch_T_arr[index]
        stream.vle(T = pinch_T_arr[index], P = stream.P)
        stream_pinch = stream.copy()
        H_pinch = stream.H
        stream.T = T_out_arr[index]
        stream.vle(T = T_out_arr[index], P = stream.P)
        H_out = stream.H
        
        if not is_cold(index):
            dH1 = abs(H_pinch - H_in)
            dH2 = abs(H_out - H_pinch)
            
            if abs(dH1)<0.01: dH1 = 0 
            if abs(dH2)<0.01: dH2 = 0
            
            Q_hot_side[index] = ['cool', dH1]
            Q_cold_side[index] = ['cool', dH2]

        else:
            dH1 = H_out - H_pinch
            dH2 = H_pinch - H_in
            
            if abs(dH1)<0.01: dH1 = 0 
            if abs(dH2)<0.01: dH2 = 0
            
            Q_hot_side[index] = ['heat', dH1]
            Q_cold_side[index] = ['heat', dH2]
            
            print('\n-----')
            print(T_in, pinch_T_arr[index], T_out_arr[index])
            print(dH1, dH2)
            
    print(Q_hot_side)
    print(Q_cold_side)
    
def get_T_transient(pinch_T_arr, indices, T_in_arr):
    T_transient = pinch_T_arr.copy()
    T_transient[indices] = T_in_arr[indices]
    return T_transient

def synthesize_network(hus, T_min_app = 10, find = None):
    
    ID_original =  bst.main_flowsheet.ID
    bst.main_flowsheet.set_flowsheet('%s-HXN'%ID_original)
        
    pinch_T_arr, hot_util_load, cold_util_load, T_in_arr, T_out_arr, T_hot_side_arr, T_cold_side_arr, \
           hus_heating, hus_cooling, hxs_heating, hxs_cooling, hxs, hot_indices, cold_indices, streams, hx_utils_rearranged = \
           temperature_interval_pinch_analysis(hus, T_min_app = 10, find = None)
           
    # hsn = HotStreamNetwork()
    # csn = ColdStreamNetwork()
    # clcs = [csn.stream_life_cycle_from_hx(hx) for hx in hxs_heating]
    # hlcs = [hsn.stream_life_cycle_from_hx(hx) for hx in hxs_cooling]
    duties = np.array([abs(hx.Q) for hx in hxs])
    
    print(T_out_arr - T_in_arr)
    
    C_flow_vector = duties/np.abs(T_in_arr - T_out_arr)
    
    # C_flow_vector = [a/(abs(b-c)) for a,b,c in zip(duties, T_in_arr, T_out_arr)]

    Q_hot_side = {}
    Q_cold_side = {}

    
    T_transient_hot_side = get_T_transient(pinch_T_arr, hot_indices, T_in_arr)
    T_transient_cold_side = get_T_transient(pinch_T_arr, cold_indices, T_in_arr)
    
    indices = hot_indices + cold_indices
        
    stream_HXs_dict = {i:[] for i in indices}
    
    is_cold = lambda x: x in cold_indices
    load_duties(streams, pinch_T_arr, T_out_arr, indices, is_cold, Q_hot_side, Q_cold_side)
    
    # Aim is to have a network in which we only heat on the hs and only cool on the cs
    matches_hs = {i: [] for i in cold_indices}
    matches_cs = {i: [] for i in hot_indices}
    
    candidate_hot_streams = list(hot_indices)
    candidate_cold_streams = list(cold_indices)
    HXs_hot_side = []
    HXs_cold_side = []
    
    # ------------- Cold side design ------------- #     
    # print('\n\n--- Cold side design ---')
    # print(T_transient_cold_side)
    # print(pinch_T_arr)

    unavailables = set([i for i in hot_indices if T_out_arr[i] >= pinch_T_arr[i]])
    unavailables.update([i for i in cold_indices if T_in_arr[i] >= pinch_T_arr[i]])
    
    # for i in range(len(pinch_T_arr)):
    #     if i in hot_indices:
    #         if T_out_arr[i] >= pinch_T_arr[i]:
    #             unavailables.append(i)
    #     elif i in cold_indices:
    #         if T_in_arr[i] >= pinch_T_arr[i]:
    #             unavailables.append(i)
    # print(unavailables)      
    
    
    
    for hot in hot_indices:
        stream_quenched = False
        original_hot_stream = streams[hot]
        
        for cold in cold_indices:
            original_cold_stream = streams[cold]
            if C_flow_vector[hot]>= C_flow_vector[cold] \
                and T_transient_cold_side[hot] > T_transient_cold_side[cold] + T_min_app \
                and (hot not in unavailables) and (cold not in unavailables) \
                and (cold not in matches_cs[hot]) and (cold in candidate_cold_streams):

                Q_hstr = Q_cold_side[hot][1]
                Q_cstr = Q_cold_side[cold][1]
                
                Q_res = Q_cstr - Q_hstr
                
                hot_stream = original_hot_stream.copy()
                hot_stream.vle(T = T_transient_cold_side[hot], P = hot_stream.P)
                
                cold_stream = original_cold_stream.copy()
                cold_stream.vle(T = T_transient_cold_side[cold], P = cold_stream.P)
                
                if abs(T_transient_cold_side[cold] - pinch_T_arr[cold])<= 0.01:
                    continue
                
                ID = 'HX_%s_%s'%(hot, cold)
                new_HX = bst.units.HXprocess(ID = ID, ins = (hot_stream, cold_stream), outs = (),
                         T_lim0 = T_out_arr[hot], T_lim1 = pinch_T_arr[cold])
                new_HX.simulate()
                HXs_cold_side.append(new_HX)
                
                stream_HXs_dict[hot].append(new_HX)
                stream_HXs_dict[cold].append(new_HX)
                
                Q_cold_side[hot][1] -= new_HX.Q
                Q_cold_side[cold][1] -= new_HX.Q
                
                T_transient_cold_side[hot] = new_HX.outs[0].T
                T_transient_cold_side[cold] = new_HX.outs[1].T
                
                stream_quenched = T_transient_cold_side[hot] <= T_out_arr[hot]
                
                matches_cs[hot].append(cold)
                # !!! TODO: remove matches_cs, matches_hs and dependencies
                # since we're already making HXs_cold_side and HXs_hot_side
                
                if stream_quenched:
                    break

        
    # Hot side design
    # print('\n\n--- Hot side design ---')
    # print(T_transient_hot_side)
    # print(pinch_T_arr)
    # unavailables = []
    unavailables = set([i for i in hot_indices if T_in_arr[i] <= pinch_T_arr[i]])
    unavailables.update([i for i in cold_indices if T_out_arr[i] <= pinch_T_arr[i]])
    
    # for i in range(len(pinch_T_arr)):
    #     if i in hot_indices:
    #         if T_in_arr[i] <= pinch_T_arr[i]:
    #             unavailables.append(i)
    #     elif i in cold_indices:
    #         if T_out_arr[i] <= pinch_T_arr[i]:
    #             unavailables.append(i)
    # print(unavailables)            
    for cold in cold_indices:
        stream_quenched = False
        original_cold_stream = streams[cold]
        
        for hot in candidate_hot_streams:
            original_hot_stream = streams[hot]
            # print(hot,cold)
            # print(C_flow_vector)
            if C_flow_vector[cold]>= C_flow_vector[hot] \
                and T_transient_hot_side[hot] > T_transient_hot_side[cold] + T_min_app \
                and (hot not in unavailables) and (cold not in unavailables) \
                and (hot not in matches_hs[cold]) and (hot in candidate_hot_streams):
            # if T_transient_hot_side[hot] - T_transient_hot_side[cold] > T_min_app:
                
                # Q_hs = C_flow_vector[hot] * (T_transient_hot_side[hot] - pinch_T_arr[hot])
                # Q_cs = C_flow_vector[cold] * (T_out_arr[cold] - T_transient_hot_side[cold])
                
                Q_hstr = Q_hot_side[hot][1]
                Q_cstr = Q_hot_side[cold][1]
                    
                Q_res = Q_cstr - Q_hstr
                
                cold_stream = original_cold_stream.copy()
                cold_stream.vle(T = T_transient_hot_side[cold], P = cold_stream.P)
                
                hot_stream = original_hot_stream.copy()
                hot_stream.vle(T = T_transient_hot_side[hot], P = hot_stream.P)
                
                if abs(T_transient_hot_side[hot] - pinch_T_arr[hot])<= 0.01:
                    continue
                ID = 'HX_%s_%s'%(cold, hot)
                new_HX = bst.units.HXprocess(ID = ID, ins = (cold_stream, hot_stream), outs = (),
                         T_lim0 = T_out_arr[cold], T_lim1 = pinch_T_arr[hot])
                new_HX.simulate()
                HXs_hot_side.append(new_HX)
                
                stream_HXs_dict[hot].append(new_HX)
                stream_HXs_dict[cold].append(new_HX)
                
                Q_hot_side[hot][1] -= new_HX.Q
                Q_hot_side[cold][1] -= new_HX.Q
                
                T_transient_hot_side[cold] = new_HX.outs[0].T
                T_transient_hot_side[hot] = new_HX.outs[1].T
                
                
                stream_quenched = T_transient_hot_side[cold] >= T_out_arr[cold]

                matches_hs[cold].append(hot)

                if stream_quenched:
                    break


    
    
    # Offset heating requirement on cold side
    
    for cold in cold_indices:
        original_cold_stream = streams[cold]
        if Q_cold_side[cold][0]=='heat' and Q_cold_side[cold][1]>0:
            for hot in hot_indices:
                original_hot_stream = streams[hot]
                if Q_cold_side[hot][0]=='cool' and Q_cold_side[hot][1]>0\
                    and T_transient_cold_side[hot] - T_transient_cold_side[cold] >= T_min_app:
                    # print('\n\n\n\n\nARGH\n\n\n\n')
                    # if Q_cold_side[hot][1]>Q_cold_side[cold][1]:
                    #     Q_cold_side[hot][1] -= Q_cold_side[cold][1]
                    #     Q_cold_side[cold][1] = 0
                    # else:
                    #     # print('\n')
                    #     # print(cold, Q_cold_side[cold][1], Q_cold_side[hot][1])
                    #     Q_cold_side[cold][1] -= Q_cold_side[hot][1]
                    #     Q_cold_side[hot][1] = 0
                    #     # print('\n')
                    #     # print(cold, Q_cold_side[cold][1], Q_cold_side[hot][1])
                    hot_stream = original_hot_stream.copy()
                    hot_stream.vle(T = T_transient_cold_side[hot], P = hot_stream.P)
                    
                    cold_stream = original_cold_stream.copy()
                    cold_stream.vle(T = T_transient_cold_side[cold], P = cold_stream.P)
                    
                    
                    if abs(T_transient_cold_side[cold] - pinch_T_arr[cold])<= 0.01:
                        continue
                    
                    ID = 'HX_%s_%s'%(hot, cold)
                    new_HX = bst.units.HXprocess(ID = ID, ins = (hot_stream, cold_stream), outs = (),
                             T_lim0 = T_out_arr[hot], T_lim1 = pinch_T_arr[cold])
                    new_HX.simulate()
                    HXs_cold_side.append(new_HX)
                    
                    stream_HXs_dict[hot].append(new_HX)
                    stream_HXs_dict[cold].append(new_HX)
                
                    Q_cold_side[hot][1] -= new_HX.Q
                    Q_cold_side[cold][1] -= new_HX.Q
                    
                    T_transient_cold_side[hot] = new_HX.outs[0].T
                    T_transient_cold_side[cold] = new_HX.outs[1].T
                    
                    matches_cs[hot].append(cold)
                    
                    
    # Offset cooling requirement on hot side
    
    for hot in hot_indices:
        original_hot_stream = streams[hot]
        if Q_hot_side[hot][0]=='cool' and Q_hot_side[hot][1]>0:
            for cold in cold_indices:
                original_cold_stream = streams[cold]
                if Q_hot_side[cold][0]=='heat' and Q_hot_side[cold][1]>0\
                    and T_transient_hot_side[hot] - T_transient_hot_side[cold] >= T_min_app:
                    # if Q_hot_side[hot][1]>Q_hot_side[cold][1]:
                    #     Q_hot_side[hot][1] -= Q_hot_side[cold][1]
                    #     Q_hot_side[cold][1] = 0
                    # else:
                    #     Q_hot_side[cold][1] -= Q_hot_side[hot][1]
                    #     Q_hot_side[hot][1] = 0
                    #     # print('\n')
                    #     # print(cold, Q_hot_side[cold][1], Q_hot_side[hot][1])
                    cold_stream = original_cold_stream.copy()
                    cold_stream.vle(T = T_transient_hot_side[cold], P = cold_stream.P)
                    
                    hot_stream = original_hot_stream.copy()
                    hot_stream.vle(T = T_transient_hot_side[hot], P = hot_stream.P)
                    
                    if abs(T_transient_hot_side[hot] - pinch_T_arr[hot])<= 0.01:
                        continue
                    
                    ID = 'HX_%s_%s'%(cold, hot)
                    new_HX = bst.units.HXprocess(ID = ID, ins = (cold_stream, hot_stream), outs = (),
                             T_lim0 = T_out_arr[cold], T_lim1 = pinch_T_arr[hot])
                    new_HX.simulate()
                    HXs_hot_side.append(new_HX)
                    
                    stream_HXs_dict[hot].append(new_HX)
                    stream_HXs_dict[cold].append(new_HX)
                    
                    Q_hot_side[hot][1] -= new_HX.Q
                    Q_hot_side[cold][1] -= new_HX.Q
                    
                    T_transient_hot_side[cold] = new_HX.outs[0].T
                    T_transient_hot_side[hot] = new_HX.outs[1].T
                    
                    
                    stream_quenched = T_transient_hot_side[cold] >= T_out_arr[cold]                    
                    matches_hs[cold].append(hot)
                    
     
    new_HX_utils = []
    
    
    for hot in hot_indices:
        new_HX_util = None
        if T_transient_cold_side[hot] > T_out_arr[hot]:
            hot_stream = streams[hot].copy()
            hot_stream.vle(T = T_transient_cold_side[hot], P = hot_stream.P)
            ID = 'Util_%s_cold_side'%(hot)
            new_HX_util = bst.units.HXutility(ID = ID, ins = hot_stream, T = T_out_arr[hot], rigorous = True)
            new_HX_util.simulate()
            new_HX_utils.append(new_HX_util)
            stream_HXs_dict[hot].append(new_HX_util)
                
            
        if T_transient_hot_side[hot] > pinch_T_arr[hot]:
            hot_stream = streams[hot].copy()
            hot_stream.vle(T = T_transient_hot_side[hot], P = hot_stream.P)
            ID = 'Util_%s_hot_side'%(hot)
            new_HX_util = bst.units.HXutility(ID = ID, ins = hot_stream, T = pinch_T_arr[hot], rigorous = True)   
            new_HX_util.simulate()
            new_HX_utils.append(new_HX_util)
            stream_HXs_dict[hot].append(new_HX_util)
            
    for cold in cold_indices:
        if T_transient_hot_side[cold] < T_out_arr[cold]:
            cold_stream = streams[cold].copy()
            cold_stream.vle(T = T_transient_hot_side[cold], P = cold_stream.P)
            ID = 'Util_%s_hot_side'%(cold)
            new_HX_util = bst.units.HXutility(ID = ID, ins = cold_stream, T = T_out_arr[cold], rigorous = True)
            new_HX_util.simulate()
            new_HX_utils.append(new_HX_util)
            stream_HXs_dict[cold].append(new_HX_util)
            
        if T_transient_cold_side[cold] < pinch_T_arr[cold]:
            cold_stream = streams[cold].copy()
            cold_stream.vle(T = T_transient_cold_side[cold], P = cold_stream.P)
            ID = 'Util_%s_cold_side'%(cold)
            new_HX_util = bst.units.HXutility(ID = ID, ins = cold_stream, T = pinch_T_arr[cold], rigorous = True)
            new_HX_util.simulate()
            new_HX_utils.append(new_HX_util)
            stream_HXs_dict[cold].append(new_HX_util)
        
    # act_cold_util_load = sum([j for i, j in Q_cold_side.values() if i == 'cool']) + \
    #     sum([j for i, j in Q_hot_side.values() if i == 'cool'])
        
    # act_hot_util_load = sum([j for i, j in Q_hot_side.values() if i == 'heat']) + \
    #     sum([j for i, j in Q_cold_side.values() if i == 'heat'])
    
    
    new_hus = bst.process_tools.heat_exchanger_utilities_from_units(new_HX_utils)
    
    act_cold_util_load = sum([abs(hu.duty) for hu in new_hus if hu.duty<0])
    
    act_hot_util_load = sum([hu.duty for hu in new_hus if hu.duty>0])
    
    act_cold_util_cost = sum([hu.cost for hu in new_hus if hu.duty<0])
    
    act_hot_util_cost = sum([hu.cost for hu in new_hus if hu.duty>0])
    
    try: amh = act_hot_util_load/hot_util_load
    except: amh = None
    
    try: amc = act_cold_util_load/cold_util_load
    except: amc = None
    print('\n')
    print('Min Heat Req. = %s, Act. Heat Req. = %s, Act/Min = %s'%\
          (hot_util_load, act_hot_util_load,amh))
    print('Min Cool Req. = %s, Act. Cool Req. = %s, Act/Min = %s'%\
          (cold_util_load,act_cold_util_load,amc))
    
    orig_heat_util = sum([hu.duty for hu in hus_heating])
    orig_cool_util = sum([abs(hu.duty) for hu in hus_cooling])
    
    print('\nAct/Original (QCool) = %s, Act/Original (QHeat) = %s'%\
          (act_cold_util_load/orig_cool_util, act_hot_util_load/orig_heat_util))
    
    orig_heat_util_cost = sum([hu.cost for hu in hus_heating])
    orig_cool_util_cost = sum([hu.cost for hu in hus_cooling])
        
    print('\nAct/Original (CostUtilCool) = %s, Act/Original (CostUtilHeat) = %s'%\
          (act_cold_util_cost/orig_cool_util_cost, act_hot_util_cost/orig_heat_util_cost))
        
    Q_prev_heating = sum([hx.Q for hx in hxs_heating])
    Q_prev_cooling = sum([abs(hx.Q) for hx in hxs_cooling])
    
    Q_prev = Q_prev_heating + Q_prev_cooling
    
    Q_HXp = sum([hx.Q for hx in HXs_hot_side]) + sum([hx.Q for hx in HXs_cold_side])
    
    Q_new_heating = sum([hx_util.Q for hx_util in new_HX_utils if hx_util.ins[0].T < hx_util.outs[0].T])
    Q_new_cooling = sum([abs(hx_util.Q) for hx_util in new_HX_utils if hx_util.ins[0].T > hx_util.outs[0].T])
    Q_new_utils = Q_new_heating + Q_new_cooling
        
    Q_new = 2*Q_HXp + Q_new_utils
    print('\n2*Q_HXp/Q_prev = %s'%(2*Q_HXp/(Q_prev)))
    print('\nQ balance: Q_new/Q_prev = %s'%(Q_new/Q_prev))
    
    print('Q_new_heating/Q_prev_heating = %s, Q_new_cooling/Q_prev_cooling = %s'%(\
            Q_new_heating/Q_prev_heating,Q_new_cooling/Q_prev_cooling ))
    assert act_hot_util_load>=0.9*hot_util_load
    assert act_cold_util_load>=0.9*cold_util_load

    
        
    return matches_hs, matches_cs, Q_hot_side, Q_cold_side, unavailables, act_hot_util_load,\
        act_cold_util_load, HXs_hot_side, HXs_cold_side, new_HX_utils, hxs, T_in_arr,\
            T_out_arr, pinch_T_arr, C_flow_vector, hx_utils_rearranged, streams, stream_HXs_dict,\
                hot_indices, cold_indices
