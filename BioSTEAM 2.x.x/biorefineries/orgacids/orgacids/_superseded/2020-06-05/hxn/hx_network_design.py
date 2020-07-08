# -*- coding: utf-8 -*-
"""
Created on Sat May  2 16:44:24 2020

@author: sarangbhagwat
"""


# import numpy as np
import copy

__all__ = ('Temperature_Interval_Method', 'Design_Network')


def Temperature_Interval_Method(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T = 10, duties = None):
    hot_indices, cold_indices = [], []
    
    for i in range(len(T_in_vector)):
        dT = T_in_vector[i] - T_out_vector[i]
        assert not dT == 0
        
        if dT > 0:
            hot_indices.append(i)
        else:
            cold_indices.append(i)
    
    hot_side_temperatures = []
    cold_side_temperatures = []
    
    for i in range(len(T_in_vector)):
        if i in hot_indices:
            hot_side_temperatures.append(T_in_vector[i])
            cold_side_temperatures.append(T_out_vector[i])
        else:
            hot_side_temperatures.append(T_out_vector[i])
            cold_side_temperatures.append(T_in_vector[i])
            
    # hot_indices = [0,1]
    # cold_indices = [2,3]
    
    
    dVF_vector = [a-b for a,b in zip(VF_in_vector,VF_out_vector)]
    
    if duties == None:
        C_flow_vector = [a*b for a,b in zip(Cp_vector,Flow_vector)]
    else:
        C_flow_vector = [a/(abs(b-c)) for a,b,c in zip(duties, T_in_vector, T_out_vector)]
    
    LH_flow_vector = [a*b*c for a,b,c in zip(LH_vector, dVF_vector, Flow_vector)]
    
    # hots_adjusted = False
    
    def give_adjusted_hots(T_in_vector, T_out_vector, hot_indices, min_app_T):
        adj_T_in_vector, adj_T_out_vector = [], []
        for i in range(len(T_in_vector)):
            
    
            if i in hot_indices:
                adj_T_in_vector.append(T_in_vector[i] - min_app_T)
                adj_T_out_vector.append(T_out_vector[i] - min_app_T)
            else:
                adj_T_in_vector.append(T_in_vector[i])
                adj_T_out_vector.append(T_out_vector[i])
        # hots_adjusted = True
        return adj_T_in_vector, adj_T_out_vector
    
    adj_T_in_vector, adj_T_out_vector = \
        give_adjusted_hots(T_in_vector, T_out_vector, hot_indices, min_app_T)
    
    all_Ts_descending = [*adj_T_in_vector, *adj_T_out_vector]
    all_Ts_descending.sort(reverse=True)
    
    T_changes_tuples = []
    
    for i in range(len(adj_T_in_vector)):
        T_changes_tuples.append((adj_T_in_vector[i], adj_T_out_vector[i]))
        
    
    stream_indices_for_T_intervals = {}
    H_for_T_intervals = {}
    
    
    for i in range(len(all_Ts_descending)-1):
        T_start = all_Ts_descending[i]
        T_end = all_Ts_descending[i+1]
        stream_indices_for_T_intervals[(T_start, T_end)] = []   
        H_for_T_intervals[(T_start, T_end)] = 0 
        
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
                if stream_index in cold_indices:
                    multiplier = -1
                # print('Stream %s'%stream_index)
                H = multiplier * (C_flow_vector[stream_index] * dT + LH_flow_vector[stream_index])
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
    pinch_cold_side_T = all_Ts_descending[res_H_vector.index(-hot_util_load)+1] # the lower temperature of the temperature interval for which the res_H is minimum
    
    cold_util_load = res_H_vector[len(res_H_vector)-1] + hot_util_load
    
    # assert cold_util_load>=0
    
    pinch_temperatures = []
    for i in range(len(T_in_vector)):
        if i in hot_indices:
            if T_in_vector[i]<pinch_cold_side_T + min_app_T:
                pinch_temperatures.append(T_in_vector[i])
            else:
                
                pinch_temperatures.append(pinch_cold_side_T + min_app_T)
        else:
            assert i in cold_indices
            if T_in_vector[i]>pinch_cold_side_T + min_app_T:
                pinch_temperatures.append(T_in_vector[i])
            else:
                pinch_temperatures.append(pinch_cold_side_T)
            
    # print('HEEEEEEEEEEEERE')
    # print(pinch_temperatures)
    # print(hot_util_load, cold_util_load)
    return(pinch_temperatures, hot_util_load, cold_util_load, hot_side_temperatures, cold_side_temperatures, hot_indices, cold_indices)
    


def Design_Network(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T = 10, duties = None, no_hot_cooling_reqs = True):
    
    pinch_temperatures, hot_util_load, cold_util_load, hot_side_temperatures, cold_side_temperatures, hot_indices, cold_indices = Temperature_Interval_Method(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T, duties)
    # if no_hot_cooling_reqs:
        
    #     for hot_index in hot_indices:
    #         T_out_vector[hot_index] = pinch_temperatures[hot_index]
    #     pinch_temperatures, hot_util_load, cold_util_load, hot_side_temperatures, cold_side_temperatures, hot_indices, cold_indices = Temperature_Interval_Method(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T, duties)    
    
    dVF_vector = [a-b for a,b in zip(VF_in_vector,VF_out_vector)]
    
    if duties == None:
        C_flow_vector = [a*b for a,b in zip(Cp_vector,Flow_vector)]
    else:
        C_flow_vector = [a/(abs(b-c)) for a,b,c in zip(duties, T_in_vector, T_out_vector)]
        
    LH_flow_vector = [a*b*c for a,b,c in zip(LH_vector, dVF_vector, Flow_vector)]
    
    matches_hs = {}
    
    matches_cs = {}
    
    
    Q_hot_side = {}
    Q_cold_side = {}
    
    
        
    candidate_hot_streams = copy.deepcopy(hot_indices)
    candidate_cold_streams = copy.deepcopy(cold_indices)
    
    
    T_transient_hot_side = copy.deepcopy(pinch_temperatures)
    
    for hot in hot_indices:
        T_transient_hot_side[hot] = T_in_vector[hot]
        Q_hot_side[hot] = ['cool', C_flow_vector[hot] * (T_in_vector[hot] - pinch_temperatures[hot])]
        Q_cold_side[hot] = ['cool', C_flow_vector[hot] * (pinch_temperatures[hot] - T_out_vector[hot])]
        # print('\n')
        # print(hot, Q_hot_side[hot])
        # print(hot, Q_cold_side[hot])
        matches_cs[hot] = []
    # print('-----------------------------------------------------')
    T_transient_cold_side = copy.deepcopy(pinch_temperatures)
    for cold in cold_indices:
        T_transient_cold_side[cold] = T_in_vector[cold]
        Q_hot_side[cold] = ['heat', C_flow_vector[cold] * (T_out_vector[cold] - pinch_temperatures[cold])]
        Q_cold_side[cold] = ['heat', C_flow_vector[cold] * (pinch_temperatures[cold] - T_in_vector[cold])]
        matches_hs[cold] = []
        # print('\n')
        # print(cold, Q_hot_side[cold])
        # print(cold, Q_cold_side[cold])
    # print(pinch_temperatures)
    # print(hot_util_load, cold_util_load)
    # print(Q_hot_side)
    # print(Q_cold_side)
    # print('\n\n')
    # print(T_transient_hot_side)
    # print(T_transient_cold_side)
    
    
    # Cold side design     
    # print('\n\n--- Cold side design ---')
    # print(T_transient_cold_side)
    # print(pinch_temperatures)
    unavailables = []
    for i in range(len(pinch_temperatures)):
        if i in hot_indices:
            if T_out_vector[i] >= pinch_temperatures[i]:
                unavailables.append(i)
        elif i in cold_indices:
            if T_in_vector[i] >= pinch_temperatures[i]:
                unavailables.append(i)
    # print(unavailables)      
    for hot in hot_indices:
        stream_quenched = False
        
        for cold in candidate_cold_streams:
            # print(cold)
            if C_flow_vector[hot]>= C_flow_vector[cold] \
                and T_transient_cold_side[hot] > T_transient_cold_side[cold] \
                and (hot not in unavailables) and (cold not in unavailables) \
                and (cold not in matches_cs[hot]):
            # if T_transient_cold_side[hot] - T_transient_cold_side[cold] > min_app_T:
                
                # Q_hot = C_flow_vector[hot] * (T_transient_cold_side[hot] - T_out_vector[hot])
                # Q_cold = C_flow_vector[cold] * (pinch_temperatures[cold] - T_transient_cold_side[cold])
                # Q_res = Q_cold - Q_hot
                
                Q_hstr = Q_cold_side[hot][1]
                Q_cstr = Q_cold_side[cold][1]
                    
                Q_res = Q_cstr - Q_hstr
                # print(hot, cold, Q_hstr, Q_cstr)
                if Q_res>0:
                    
                    # hot stream cooled to output temp
                    Q_cold_side[cold][1] -= Q_hstr # reduces some heating duty reqd by cold stream
                    Q_cold_side[hot][1] = 0
                    T_transient_cold_side[cold] = (Q_hstr/C_flow_vector[cold] + T_transient_cold_side[cold])
                    T_transient_cold_side[hot] = T_out_vector[hot]
                    # candidate_hot_streams.remove(hot)
                    
                elif Q_res<0:
                    
                    # cold stream heated to pinchtemp
                    candidate_cold_streams.remove(cold)
                    Q_cold_side[hot][1] -= Q_cstr # reduces some cooling duty reqd by hot stream
                    Q_cold_side[cold][1] = 0
                    T_transient_cold_side[cold] = pinch_temperatures[cold]
                    T_transient_cold_side[hot] = -(Q_cstr/C_flow_vector[hot] - T_transient_cold_side[hot])
                    stream_quenched = True
                    
                else:
                    candidate_cold_streams.remove(cold)
                    T_transient_cold_side[cold] = T_out_vector[cold]
                    T_transient_cold_side[hot] = pinch_temperatures[hot]
                    Q_cold_side[hot][1] = 0
                    Q_cold_side[cold][1] = 0
                    stream_quenched = True
                    
                #     print('Perfect match! %s <-> %s'%(cold, hot))
                
                matches_cs[hot].append(cold)
                # match_found = True
                # print('Match: %s --> %s (Q_res = %s)'%(hot, cold, Q_res))
                # print('Match: %s (brought to %s K) --> %s (brought to %s K) [Q_res = %s]'\
                      # %(hot, T_transient_cold_side[hot], cold, T_transient_cold_side[cold], Q_res))
                
                if stream_quenched:
                    break
            # if match_found:
            #     break
        # if not match_found:
        #     # Q_cold_side[hot] = ['cool', C_flow_vector[hot] * (T_transient_cold_side[hot] - T_out_vector[hot])]
        #     Q_cold_side[hot] = ['cool', C_flow_vector[hot] * (T_in_vector[hot] - pinch_temperatures[hot])]
        
    # Hot side design
    # print('\n\n--- Hot side design ---')
    # print(T_transient_hot_side)
    # print(pinch_temperatures)
    unavailables = []
    for i in range(len(pinch_temperatures)):
        if i in hot_indices:
            if T_in_vector[i] <= pinch_temperatures[i]:
                unavailables.append(i)
        elif i in cold_indices:
            if T_out_vector[i] <= pinch_temperatures[i]:
                unavailables.append(i)
    # print(unavailables)            
    for cold in cold_indices:
        stream_quenched = False
        
        for hot in candidate_hot_streams:
            # print(hot,cold)
            # print(C_flow_vector)
            if C_flow_vector[cold]>= C_flow_vector[hot] \
                and T_transient_hot_side[hot] > T_transient_hot_side[cold] \
                and (hot not in unavailables) and (cold not in unavailables) \
                and (hot not in matches_hs[cold]):
            # if T_transient_hot_side[hot] - T_transient_hot_side[cold] > min_app_T:
                
                # Q_hs = C_flow_vector[hot] * (T_transient_hot_side[hot] - pinch_temperatures[hot])
                # Q_cs = C_flow_vector[cold] * (T_out_vector[cold] - T_transient_hot_side[cold])
                
                Q_hstr = Q_hot_side[hot][1]
                Q_cstr = Q_hot_side[cold][1]
                    
                Q_res = Q_cstr - Q_hstr
                
                if Q_res>0:
                    # hot stream cooled to pinch temp
                    
                    candidate_hot_streams.remove(hot)
                    Q_hot_side[cold][1] -= Q_hstr # reduces some heating duty reqd by cold stream
                    Q_hot_side[hot][1] = 0
                    T_transient_hot_side[cold] = (Q_hstr/C_flow_vector[cold] + T_transient_hot_side[cold])
                    T_transient_hot_side[hot] = pinch_temperatures[hot]
                    
                    
                elif Q_res<0:
                    # cold stream heated to output temp

                    Q_hot_side[hot][1] -= Q_cstr # reduces some cooling duty reqd by hot stream
                    Q_hot_side[cold][1] = 0
                    T_transient_hot_side[cold] = T_out_vector[cold]
                    T_transient_hot_side[hot] = -(Q_cstr/C_flow_vector[hot] - T_transient_hot_side[hot])
                    stream_quenched = True
                    
                else:
                    candidate_hot_streams.remove(hot)
                    T_transient_hot_side[cold] = T_out_vector[cold]
                    T_transient_hot_side[hot] = pinch_temperatures[hot]
                    Q_hot_side[cold][1] = 0
                    Q_hot_side[hot][1] = 0
                    stream_quenched = True
                #     print('Perfect match! %s <-> %s'%(cold, hot))
                
                matches_hs[cold].append(hot)
                # match_found = True
                # print('Match: %s (brought to %s K) <-- %s (brought to %s K) [Q_res = %s]'\
                #       %(cold, T_transient_hot_side[cold], hot, T_transient_hot_side[hot], Q_res))
                
                if stream_quenched:
                    break
            # if stream_quenched:
            #     break
        # if not match_found:
            # Q_hot_side[cold] = ['heat', C_flow_vector[cold] * (T_out_vector[cold] - T_transient_hot_side[cold])]
            # Q_hot_side[cold] = ['heat', C_flow_vector[cold] * (pinch_temperatures[cold] - T_transient_hot_side[cold])]
    

    
    
    # Offset heating requirement on cold side
    
    for cold in cold_indices:
        if Q_cold_side[cold][0]=='heat' and Q_cold_side[cold][1]>0:
            for hot in hot_indices:
                if Q_cold_side[hot][0]=='cool' and Q_cold_side[hot][1]>0\
                    and T_transient_cold_side[hot] - T_transient_cold_side[cold] >= min_app_T:
                    # print('\n\n\n\n\nARGH\n\n\n\n')
                    if Q_cold_side[hot][1]>Q_cold_side[cold][1]:
                        Q_cold_side[hot][1] -= Q_cold_side[cold][1]
                        Q_cold_side[cold][1] = 0
                    else:
                        # print('\n')
                        # print(cold, Q_cold_side[cold][1], Q_cold_side[hot][1])
                        Q_cold_side[cold][1] -= Q_cold_side[hot][1]
                        Q_cold_side[hot][1] = 0
                        # print('\n')
                        # print(cold, Q_cold_side[cold][1], Q_cold_side[hot][1])
                    matches_cs[hot].append(cold)
                    
                    
    # Offset cooling requirement on hot side
    
    for hot in hot_indices:
        if Q_hot_side[hot][0]=='cool' and Q_hot_side[hot][1]>0:
            for cold in cold_indices:
                if Q_hot_side[cold][0]=='heat' and Q_hot_side[cold][1]>0\
                    and T_transient_hot_side[hot] - T_transient_hot_side[cold] >= min_app_T:
                    if Q_hot_side[hot][1]>Q_hot_side[cold][1]:
                        Q_hot_side[hot][1] -= Q_hot_side[cold][1]
                        Q_hot_side[cold][1] = 0
                    else:
                        Q_hot_side[cold][1] -= Q_hot_side[hot][1]
                        Q_hot_side[hot][1] = 0
                        # print('\n')
                        # print(cold, Q_hot_side[cold][1], Q_hot_side[hot][1])
                    
                    matches_hs[cold].append(hot)
                    
       
    act_hot_util_load, act_cold_util_load = 0, 0
    for key, item in Q_hot_side.items():
        if item[0]=='heat':
            act_hot_util_load += item[1]
        elif item[0] == 'cool':
            act_cold_util_load += item[1]
            
    for key, item in Q_cold_side.items():
        if item[0]=='heat':
            act_hot_util_load += item[1]
        elif item[0] == 'cool':
            act_cold_util_load += item[1]
    
    print('\n')
    print(act_hot_util_load,hot_util_load)
    print(act_cold_util_load,cold_util_load)
    
    # assert act_hot_util_load>=hot_util_load
    # assert act_cold_util_load>=cold_util_load
    return matches_hs, matches_cs, Q_hot_side, Q_cold_side, unavailables, act_hot_util_load, act_cold_util_load
