# -*- coding: utf-8 -*-
"""
Created on Thu May 14 02:50:49 2020

@author: saran
"""

def Design_Network_old(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T = 10, duties = None, no_hot_cooling_reqs = True):
    
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
        matches_cs[hot] = []
        
    T_transient_cold_side = copy.deepcopy(pinch_temperatures)
    for cold in cold_indices:
        T_transient_cold_side[cold] = T_in_vector[cold]
        matches_hs[cold] = []
        
    # print('\n\n')
    # print(T_transient_hot_side)
    # print(T_transient_cold_side)
            
    # Hot side design
    print('\n\n--- Hot side design ---')
    unavailables = []
    for i in range(len(pinch_temperatures)):
        if i in hot_indices:
            if T_in_vector[i] == pinch_temperatures[i]:
                unavailables.append(i)
        elif i in cold_indices:
            if T_out_vector[i] == pinch_temperatures[i]:
                unavailables.append(i)
                
    for cold in cold_indices:
        match_found = False
        
        for hot in candidate_hot_streams:
            # print(hot,cold)
            # print(C_flow_vector)
            if C_flow_vector[cold]>= C_flow_vector[hot] \
                and T_transient_hot_side[hot] > T_transient_hot_side[cold] \
                and (hot not in unavailables) and (cold not in unavailables):
            # if T_transient_hot_side[hot] - T_transient_hot_side[cold] > min_app_T:
                
                Q_hot = C_flow_vector[hot] * (T_transient_hot_side[hot] - pinch_temperatures[hot])
                Q_cold = C_flow_vector[cold] * (T_out_vector[cold] - T_transient_hot_side[cold])
                Q_res = Q_cold - Q_hot
                
                if Q_res>0:
                    Q_hot_side[cold] = ['heat', Q_res]
                    T_transient_hot_side[cold] = (Q_hot/C_flow_vector[cold] + T_transient_hot_side[cold])
                    T_transient_hot_side[hot] = pinch_temperatures[hot]
                    candidate_hot_streams.remove(hot)
                    
                elif Q_res<0:
                    Q_hot_side[hot] = ['cool', -Q_res]
                    T_transient_hot_side[cold] = T_out_vector[cold]
                    T_transient_hot_side[hot] = -(Q_cold/C_flow_vector[hot] - pinch_temperatures[hot])
                    
                else:
                    candidate_hot_streams.remove(hot)
                    T_transient_hot_side[cold] = T_out_vector[cold]
                    T_transient_hot_side[hot] = pinch_temperatures[hot]
                #     print('Perfect match! %s <-> %s'%(cold, hot))
                
                matches_hs[cold].append(hot)
                match_found = True
                print('Match: %s <-- %s (Q_res = %s)'%(cold, hot, Q_res))
    
            if match_found:
                break
        if not match_found:
            # Q_hot_side[cold] = ['heat', C_flow_vector[cold] * (T_out_vector[cold] - T_transient_hot_side[cold])]
            Q_hot_side[cold] = ['heat', C_flow_vector[cold] * (pinch_temperatures[cold] - T_in_vector[cold])]
    
    # Cold side design     
    print('\n\n--- Cold side design ---')
    unavailables = []
    for i in range(len(pinch_temperatures)):
        if i in hot_indices:
            if T_out_vector[i] == pinch_temperatures[i]:
                unavailables.append(i)
        elif i in cold_indices:
            if T_in_vector[i] == pinch_temperatures[i]:
                unavailables.append(i)
                
    for hot in hot_indices:
        match_found = False
        
        for cold in candidate_cold_streams:
            
            if C_flow_vector[hot]>= C_flow_vector[cold] \
                and T_transient_cold_side[hot] > T_transient_cold_side[cold] \
                and (hot not in unavailables) and (cold not in unavailables):
            # if T_transient_cold_side[hot] - T_transient_cold_side[cold] > min_app_T:
                
                Q_hot = C_flow_vector[hot] * (T_transient_cold_side[hot] - T_out_vector[hot])
                Q_cold = C_flow_vector[cold] * (pinch_temperatures[cold] - T_transient_cold_side[cold])
                Q_res = Q_cold - Q_hot
                
                if Q_res>0:
                    Q_cold_side[cold] = ['heat', Q_res]
                    T_transient_cold_side[cold] = -(Q_hot/C_flow_vector[cold] - pinch_temperatures[cold])
                    T_transient_cold_side[hot] = T_out_vector[hot]
                    
                elif Q_res<0:
                    Q_cold_side[hot] = ['cool', -Q_res]
                    T_transient_cold_side[cold] = pinch_temperatures[cold]
                    T_transient_cold_side[hot] = -(Q_cold/C_flow_vector[hot] - T_transient_cold_side[hot])
                    candidate_cold_streams.remove(cold)
                else:
                    candidate_cold_streams.remove(cold)
                #     print('Perfect match! %s <-> %s'%(cold, hot))
                
                matches_cs[hot].append(cold)
                match_found = True
                print('Match: %s --> %s (Q_res = %s)'%(hot, cold, Q_res))
            
            if match_found:
                break
        if not match_found:
            # Q_cold_side[hot] = ['cool', C_flow_vector[hot] * (T_transient_cold_side[hot] - T_out_vector[hot])]
            Q_cold_side[hot] = ['cool', C_flow_vector[hot] * (T_in_vector[hot] - pinch_temperatures[hot])]
    
    # Offset heating requirement on cold side
    q = copy.deepcopy(Q_cold_side)
    to_delete = []
    # items = copy.deepcopy(Q_cold_side.items())
    for stream1, util1 in q.items():
        if not stream1 in to_delete:
            
            if Q_cold_side[stream1][0]=='heat':
                for stream2, util2 in q.items():
                    if (not stream2 in to_delete) and (not stream1 in to_delete):
                        if Q_cold_side[stream2][0] == 'cool':
                            if T_transient_cold_side[stream2] - T_transient_cold_side[stream1] >= min_app_T:
                                matches_cs[stream2].append(stream1)
                                if Q_cold_side[stream2][1]>Q_cold_side[stream1][1]:
                                    Q_cold_side[stream2][1] -= Q_cold_side[stream1][1]
                                    if stream1 in Q_cold_side.keys():
                                        del(Q_cold_side[stream1])
                                        to_delete.append(stream1)
                                elif Q_cold_side[stream1][1]>Q_cold_side[stream2][1]:
                                    # print(q[stream1])
                                    # print(util2[1])
                                    # print(Q_cold_side[stream1][1])
                                    
                                    # print(Q_cold_side[stream2][1])
                                    Q_cold_side[stream1][1] -= Q_cold_side[stream2][1]
                                    if stream2 in Q_cold_side.keys():
                                        del(Q_cold_side[stream2])
                                        to_delete.append(stream2)
                                        # del(q[stream2])
                                else:
                                    if stream1 in Q_cold_side.keys() and stream2 in Q_cold_side.keys():
                                        del(Q_cold_side[stream1])
                                        # del(q[stream1])
                                        del(Q_cold_side[stream2])
                                        
                                        to_delete.append(stream1)
                                        to_delete.append(stream2)
                                        # del(q[stream2])
                            
    # Offset cooling requirement on hot side
    q = copy.deepcopy(Q_hot_side)
    to_delete = []
    # items = copy.deepcopy(Q_hot_side.items())
    for stream1, util1 in q.items():
        if not stream1 in to_delete:
            if Q_hot_side[stream1][0]=='cool':
                for stream2, util2 in q.items():
                    if (not stream2 in to_delete) and (not stream1 in to_delete):
                        if Q_hot_side[stream2][0] == 'heat':
                            if T_transient_hot_side[stream1] - T_transient_hot_side[stream2] >= min_app_T:
                                matches_hs[stream2].append(stream1)
                                if Q_hot_side[stream2][1]>Q_hot_side[stream1][1]:
                                    Q_hot_side[stream2][1] -= Q_hot_side[stream1][1]
                                    del(Q_hot_side[stream1])
                                    to_delete.append(stream1)
                                elif Q_hot_side[stream1][1]>Q_hot_side[stream2][1]:
                                    Q_hot_side[stream1][1] -= Q_hot_side[stream2][1]
                                    del(Q_hot_side[stream2])
                                    to_delete.append(stream2)
                                else:
                                    if stream1 in Q_hot_side.keys() and stream2 in Q_hot_side.keys():
                                        del(Q_hot_side[stream1])
                                        del(Q_hot_side[stream2])
                                        to_delete.append(stream1)
                                        to_delete.append(stream2)
    # print('\n\n')
    # print(T_transient_hot_side)
    # print(T_transient_cold_side)
    # Matchless cold side
    # for hot in hot_indices:
    #     if hot not in match.keys():
    #         for cold in candidate_cold_streams:
    #             if T_transient_cold_side[hot] - T_transient_cold_side[cold] > min_app_T: # must be >, not >=
                    
    #                 Q_hot = C_flow_vector[hot] * (T_transient_cold_side[hot] - T_out_vector[hot])
    #                 Q_cold = C_flow_vector[cold] * (pinch_temperatures[cold] - T_transient_cold_side[cold])
    #                 Q_res = Q_cold - Q_hot
                    
    #                 if Q_res>0:
    #                     Q_cold_side[cold] = ('heat', Q_res)
    #                     T_transient_cold_side[cold] = -(Q_hot/C_flow_vector[cold]) - T_transient_cold_side[cold] 
    #                     T_transient_cold_side[hot] = pinch_temperatures[hot]
    #                     candidate_cold_streams.remove(cold)
    #                 elif Q_res<0:
    #                     T_transient_cold_side[cold] = pinch_temperatures[cold]
    #                     Q_cold_side[hot] = ('cool', -Q_res)
    #                     T_transient_cold_side[hot] = -(Q_cold/C_flow_vector[hot]) - T_transient_cold_side[hot]
                    
    #                 else:
    #                     candidate_cold_streams.remove(cold)
    #                 #     print('Perfect match! %s <-> %s'%(cold, hot))
                    
    #                 matches_cs[hot] = cold
    #                 match_found = True
                    
    #                 print('Match: %s <-> %s (Q_res = %s)'%(cold, hot, Q_res))
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
            
    print(act_hot_util_load,hot_util_load)
    print(act_cold_util_load,cold_util_load)
    
    # assert act_hot_util_load>=hot_util_load
    # assert act_cold_util_load>=cold_util_load
    return matches_hs, matches_cs, Q_hot_side, Q_cold_side, unavailables, act_hot_util_load, act_cold_util_load

