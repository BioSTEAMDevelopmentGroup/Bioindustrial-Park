#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 08:24:22 2020

@author: sarangbhagwat
"""

import copy
import biosteam as bst
from biosteam import HeatUtility, Facility
from biosteam.units import HXprocess
from orgacids.hxn.hx_network_design import Temperature_Interval_Method, Design_Network

__all__ = ('HX_Network',)


class HX_Network(Facility):
    network_priority = -1
    # _N_ins = 6
    # _N_outs = 3
    _N_heat_utilities = 1
    _units= {'Flow rate': 'kg/hr',
             'Work': 'kW'} 
    
    #!!! Yalin doesn't think this is needed
    # from biosteam import HeatUtility, settings, Stream
    
    def __init__(self, ID='', ins=None, outs=(), get_HXN=None, get_streams = None, min_app_T = 10, sys = []):
        
        Facility.__init__(self, ID, ins, outs)
        self.get_HXN = get_HXN
        self.min_app_T = min_app_T
        self.sys = sys
    def _run(self):
        pass
    def _cost(self):
        # def _design(self):
                
        min_app_T = self.min_app_T
        sys = self.system
        Flow_vector = []
        
        T_in_vector = [] 
        VF_in_vector = []
        
        
        T_out_vector = []
        VF_out_vector = []
        
        Cp_vector = []
        LH_vector = []
        Tb_vector = []
        Td_vector = []
        duties = []
        
        hx_utils = bst.process_tools.heat_exchanger_utilities_from_units(sys.units)
        
        hx_utils.sort(key = lambda x: x.duty)
        stream_dict = {}
        i = 0
        
        init_hot_util = 0
        init_cold_util = 0
        
        # def increment_P_till_target_V(stream, V):
        #     # og_P = stream.P
        #     # stream.vle(P = og_P, T = stream.T)
        #     # while stream.vapor_fraction > 0:
        #     #     stream.vle(P = stream.P + og_P/10, T = stream.T)
        #     T_new = stream.T
        #     T = 0
        #     H = stream.H
        #     stream.phase = 'l'
        #     # print(stream.P)
        #     # print(T_new, T)
        #     while abs(T_new - T) > 0.01:
                
        #         try: 
        #             bp = stream.bubble_point_at_T()
        #         except: 
        #             stream.T = 500
        #             bp = stream.bubble_point_at_T()
                    
        #         # print('these')
        #         # print(bp.P)
        #         stream.P = bp.P
        #         stream.H = H
                
        #         T = T_new
        #         T_new = stream.T
                
                
        P_changes_dict = {}
        T_out_changes_dict = {}
        
        h0i = []
        c0i = []
        dT0 = []
        hx_utils_applicable = []
        hx_utils_applicable2 = []
        vaporizes = []
        
        original_purchase_costs = []
        original_installation_costs = []
        original_utility_costs = []
        original_heat_utilities = []
        for hx_util in hx_utils:
            applicable = True
            hx = hx_util.heat_exchanger
            instr = hx.ins[0].copy()
            outstr = hx.outs[0].copy()
            T_in = instr.T
            T_out = outstr.T
            
            # if outstr.vapor_fraction > 0:
            #     # outstr.vle(V = 0, H = outstr.H)
            #     # T_out = outstr.T
            #     increment_P_till_target_V(outstr, 0)
                
            #     T_out_changes_dict[i] = [hx.outs[0].T, outstr.T]
            #     T_out = outstr.T
            P_changes_dict[i] = [instr.P, outstr.P]
            
            # if T_in < T_out:
                # P_changes_dict[i] = [instr.P, 2*outstr.P]
                
                
                    
                
            # if T_in > T_out and hx.ins[0].vapor_fraction > 0:
            #     applicable = False
            
            
            if abs(T_in - T_out)>min_app_T and applicable:
                # and hx.ins[0].F_mass*hx.ins[0].Cp*1000*(abs(T_in-T_out)) <= abs(hx_util.duty):
                if T_in > T_out:
                    h0i.append(i)
                elif T_in < T_out:
                    c0i.append(i)
                    if hx.outs[0].vapor_fraction > 0:
                        vaporizes.append(i)
                else: 
                    dT0.append(i)
                    
                stream_dict[hx] = i
                
                
                duty = hx_util.duty
                duties.append(abs(duty))
                if duty>0:
                    # original_heat_utilities.append([duty,0])
                    init_hot_util+= duty
                else:
                    init_cold_util -= duty
                    # original_heat_utilities.append([0,-duty])
                assert len(hx.ins) == 1
                T_in_vector.append(T_in)
                assert len(hx.outs) == 1
                T_out_vector.append(T_out)
                Flow_vector.append(hx.ins[0].F_mass)
                VF_in_vector.append(0)
                VF_out_vector.append(0)
                Cp_vector.append(hx.ins[0].Cp) # *1000 for g-> kg, /1000 for J-> kJ
                LH_vector.append(0)
                
                hx_utils_applicable.append(hx_util)
                # hx_utils_applicable2.append(hx_util)
                original_purchase_costs.append(hx.purchase_cost)
                original_installation_costs.append(hx.installation_cost)
                original_utility_costs.append(hx.utility_cost)
                # Tb_vector.append(hx.ins[0].bubble_point_at_P().T)
                # Td_vector.append(hx.ins[0].dew_point_at_P().T)
                i+=1
        
        
        
        # duties = None
        
        # t_pinches, h, c, hst, cst, hi, ci= Temperature_Interval_Method(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T)
        min_T_in = min(T_in_vector)
        
        min_T_out = min(T_out_vector)
        
        sum_by_agent1 = HeatUtility.sum_by_agent(hx_utils_applicable)
        
        # waste_stream_indices = []
        # for waste_stream in M501.ins:
        #     if waste_stream.F_mass>0 and waste_stream.T - min_T_in >= min_app_T:
                
        #         Cp_vector.append(waste_stream.Cp)
                
        #         stream_dict[waste_stream.ID] = i
        #         T_in_vector.append(waste_stream.T)
        #         waste_stream_indices.append(len(T_in_vector) -1)
        #         T_out_vector.append(min_T_out)
        #         Flow_vector.append(waste_stream.F_mass)
        #         VF_in_vector.append(0)
        #         VF_out_vector.append(0)
                
        #         LH_vector.append(0)
        #         # Tb_vector.append(waste_stream.bubble_point_at_P().T)
        #         # Td_vector.append(waste_stream.dew_point_at_P().T)
        #         i+=1
                
        # for unit in orgacids_sys.units:
        t_pinches, h, c, hst, cst, hi, ci= Temperature_Interval_Method(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T)
        
        pinch = min(t_pinches)
        
        # for i in range(len(T_in_vector)):
            # if (T_in_vector[i] < T_out_vector[i] and T_out_vector[i] >= pinch):
            #     # T_in_vector[i] = pinch
            #     T_out_vector[i] = pinch + 10
            # elif (T_in_vector[i] > T_out_vector[i] and T_out_vector[i] <= pinch):
            #     # T_in_vector[i] = pinch + 10
            #     T_out_vector[i] = pinch 
            
            # if i in waste_stream_indices:
            #     T_out_vector[i] = pinch + min_app_T
                    
        # print(init_hot_util, init_cold_util)
        
        t_pinches2, h, c, hst, cst, hi, ci= Temperature_Interval_Method(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T)
        
        matches_hs, matches_cs, Q_hot_side, Q_cold_side, unavailables, hot_util_load, cold_util_load = Design_Network(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T=min_app_T, duties=duties)
        
        #############################################

        
        rev_stream_dict = {}
        new_HXs = []
        new_HXs_ins = []
        new_streams_out_dict = {}
        
        for exchanger, stream_index in stream_dict.items():
            rev_stream_dict[stream_index] = exchanger
        # for hx in rev_stream_dict.values():
        #     # hx.show()
        
        for stream, matches in matches_hs.items():
            HX1 = rev_stream_dict[stream]
            i = 0
            
            for match in matches:
                HX2 = rev_stream_dict[match]
                ID = 'HX_%s_%s'%(stream, match)
                # for i in range(len(ID)):
                #     if ID[i] =='.':
                #         ID[i] = '_'
                # print(ID)
                # import pdb
                # pdb.set_trace()
                in1 = HX2.ins[0].copy()
                # in1.vle
                in1.P = P_changes_dict[match][1]
                fluid_type = 'ss'
                
                if not in1.vapor_fraction == 0:
                    fluid_type = 'ls'
                
                if i == 0:
                    in0 = HX1.ins[0].copy()
                    in0.P = P_changes_dict[stream][1]
                    # in0.vle
                    new_HXs.append(HXprocess(ID=ID, ins = (in0, in1), outs = ('out0', 'out1')))
                    # print(new_HXs[len(new_HXs) - 1].ID)
                    # new_HXs[len(new_HXs) - 1]._run()
                    # new_HXs[len(new_HXs) - 1].show()
                    
                    new_HXs[len(new_HXs) - 1].simulate()
                    
                else:
                    in0 = new_HXs[len(new_HXs) - 1].outs[0].copy()
                    in0.P = P_changes_dict[stream][1]
                    # in0.vle
                    new_HXs.append(HXprocess(ID=ID, ins = (in0, in1), outs = ('out0', 'out1')))
                    # print(new_HXs[len(new_HXs) - 1].ID)
                    # # new_HXs[len(new_HXs) - 1]._run()
                    # # new_HXs[len(new_HXs) - 1].show()
                    new_HXs[len(new_HXs) - 1].simulate()
                
                new_HXs_ins.append([stream, match])
                
                new_streams_out_dict[match] = new_HXs[len(new_HXs) - 1].outs[1]
                i += 1
                
        
        ### !!! Future implementation: inherit attributes from old HXs
        
        for stream, matches in matches_cs.items():
            HX1 = rev_stream_dict[stream]
            i = 0
            
            for match in matches:
                HX2 = rev_stream_dict[match]
                ID = 'HX_%s_%s'%(stream, match)
                # for i in range(len(ID)):
                #     if ID[i] =='.':
                #         ID[i] = '_'
                # print(ID)
                # import pdb
                # pdb.set_trace()
                in1 = HX2.ins[0].copy()
                # in1.vle
                in1.P = P_changes_dict[match][1]
                
                if i == 0:
                    in0 = HX1.ins[0].copy()
                    in0.P = P_changes_dict[stream][1]
                    # in0.vle
                    fluid_type = 'ss'
                    if not in1.vapor_fraction == 0:
                        fluid_type = 'ls'
                    new_HXs.append(HXprocess(ID=ID, ins = (in0, in1), outs = ('out0', 'out1')))
                    # print(new_HXs[len(new_HXs) - 1].ID)
                    # new_HXs[len(new_HXs) - 1]._run()
                    # new_HXs[len(new_HXs) - 1].show()
                    new_HXs[len(new_HXs) - 1].simulate()
                    
                else:
                    in0 = new_HXs[len(new_HXs) - 1].outs[0].copy()
                    in0.P = P_changes_dict[stream][1]
                    # in0.vle
                    fluid_type = 'ss'
                    if not in1.vapor_fraction == 0:
                        fluid_type = 'ls'
                    new_HXs.append(HXprocess(ID=ID, ins = (in0, in1), outs = ('out0', 'out1')))
                    # print(new_HXs[len(new_HXs) - 1].ID)
                    # new_HXs[len(new_HXs) - 1]._run()
                    # new_HXs[len(new_HXs) - 1].show()
                    new_HXs[len(new_HXs) - 1].simulate()
                
                new_HXs_ins.append([stream, match])
                new_streams_out_dict[match] = new_HXs[len(new_HXs) - 1].outs[1]
                i += 1
    
        
        
        new_hx_utils = {}
        
        for i in range(len(new_HXs)):
            new_HX = new_HXs[0]
            streams = new_HXs_ins[i]
            
            for j in range(len(streams)):
                stream = streams[j]
                
                if stream in ci:
                    if T_out_vector[stream]>new_HX.outs[j].T:
                        ID = 'Util_%s'%(stream)
                        new_hx_utils[stream] = bst.units.HXutility(ID = ID, ins = new_HX.outs[j].copy(), T = T_out_vector[stream], rigorous = True)
                        # new_hx_utils[stream]._run()
                        
                        # new_hx_utils[stream].show()
                        new_hx_utils[stream].simulate()
                if stream in hi:
                    if T_out_vector[stream]<new_HX.outs[j].T:
                        ID = 'Util_%s'%(stream)
                        new_hx_utils[stream] = bst.units.HXutility(ID = ID, ins = new_HX.outs[j].copy(), T = T_out_vector[stream], rigorous = True)
                        new_hx_utils[stream].simulate()
                    
        new_purchase_costs_HXp = []
        new_purchase_costs_HXu = copy.deepcopy(original_purchase_costs)
        new_installation_costs_HXp = []
        new_installation_costs_HXu = copy.deepcopy(original_installation_costs)
        
        new_utility_costs = copy.deepcopy(original_utility_costs)
        # new_heat_utilities = copy.deepcopy(original_heat_utilities)
        # for i in range(len(new_HXu)):
        heat_utilities = []
        # new_heat_utilities2 = copy.deepcopy(hx_utils_applicable)
        # for hx_util in hx_utils_applicable:
        #     heat_utilities
            
        for s1 in new_hx_utils.keys():
            # s1, s2 = new_HXs_ins[i]
            # new_HX = new_HXs[i]
            new_installation_costs_HXu[s1] = new_hx_utils[s1].installation_cost
            new_purchase_costs_HXu[s1] = new_hx_utils[s1].purchase_cost
            new_utility_costs[s1] = new_hx_utils[s1].utility_cost
            # s1_util = bst.process_tools.heat_exchanger_utilities_from_units([new_hx_utils[s1]])[0].duty
            # if s1_util>0:
            #     new_heat_utilities[s1] = [s1_util, 0]
            # else:
            #     new_heat_utilities[s1] = [0, -s1_util]
            
            # heat_utilities[s1] = new_hx_utils[s1].heat_utilities[0]
            
            # hx_utils_applicable[s1].scale((abs(s1_util)-abs(hx_utils_applicable[s1].duty))/abs(hx_utils_applicable[s1].duty))
            # heat_utilities[len(heat_utilities)-1].duty -= 
            # new_purchase_costs_HXu[s2] = new_hx_utils[s2].purchase_cost
            # new_utility_costs[s2] = new_hx_utils[s2].utility_cost
            # s2_util = bst.process_tools.heat_exchanger_utilities_from_units([new_hx_utils[s2]])[0].duty
            # if s2_util>0:
            #     new_heat_utilities[s2] = [s2_util, 0]
            # else:
            #     new_heat_utilities[s2] = [0, -s2_util]
        
        
        for new_HX in new_HXs:
            new_purchase_costs_HXp.append(new_HX.purchase_cost)
            new_installation_costs_HXp.append(new_HX.installation_cost)
        init_heating_sum, init_cooling_sum = init_hot_util, init_cold_util
        
        # new_heating_sum, new_cooling_sum = 0, 0
        # for nhu in new_heat_utilities:
        #     new_heating_sum += nhu[0]
        #     new_cooling_sum += nhu[1]
        # for i in range(len(hx_utils_applicable)):
        self.purchase_costs['Heat exchangers'] = (sum(new_purchase_costs_HXp) + sum(new_purchase_costs_HXu)) \
            - (sum(original_purchase_costs))
        
        # for i in range(len(hx_utils_applicable)):
        #     og_duty = new_heat_utilities[i].duty
        #     new_duty = hx_utils_applicable[i].duty - hx_utils_applicable[s1].duty
        #     hx_utils_applicable[i].scale(new_duty/og_duty)
        
        #####
        
        

        hu_sums1 = HeatUtility.sum_by_agent(hx_utils_applicable)
        
        hx_utils_applicable2 = [HeatUtility() for i in range(len(hx_utils_applicable))]
        # QTs1 = []
        QTs2 = []
        # for hx in hx_utils_applicable:
        #     QTs1.append(hx.Q, hx.T)                                                  
        # QTs = [(100, 400), (-100, 400), (100, 290), (-100, 290)]
        # for hu, (Q, T) in zip(hx_utils_applicable2, QTs1): hu(Q, T)
        
        
        for stream, hx_util in new_hx_utils.items():
            # print(hx_util)
            QTs2.append((hx_util.Q, hx_util.T))
        for hu, (Q, T) in zip(hx_utils_applicable2, QTs2): hu(Q, T)
        hu_sums2 = HeatUtility.sum_by_agent(hx_utils_applicable2)
        
        
        # to change sign on duty without switching heat/cool (i.e. negative costs):
        for hu in hu_sums1: hu.reverse()
        #
        
        # new_heat_utilties = [HeatUtility() for i in range(4)]
        # QTs = [(100, 400), (-100, 400), (200, 290), (-150, 400)]
        # for hu, (Q, T) in zip(new_heat_utilties, QTs): hu(Q, T)
        # new_hu_sums = HeatUtility.sum_by_agent(new_heat_utilties)
        
        hus_final = tuple(HeatUtility.sum_by_agent(hu_sums1 + hu_sums2))

        #####
        sum_by_agent2 =HeatUtility.sum_by_agent(hx_utils_applicable)
        
        self._installation_cost = (sum(new_installation_costs_HXp) + sum(new_installation_costs_HXu)) \
            - (sum(original_installation_costs))
        
        
        self.heat_utilities = hus_final
        self.new_HXs = new_HXs
        self.new_HX_utils = new_hx_utils
        self.orig_heat_utils = hx_utils_applicable
        self.original_purchase_costs = original_purchase_costs
        self.original_utility_costs = hu_sums1
        self.new_purchase_costs_HXp = new_purchase_costs_HXp
        self.new_purchase_costs_HXu = new_purchase_costs_HXu
        self.new_utility_costs = hu_sums2
        self.sum_by_agent1 = sum_by_agent1
        self.sum_by_agent2 = sum_by_agent2
        # self.heat_utilities = (new_heating_sum - init_heating_sum, new_cooling_sum - init_cooling_sum)
        # self.utility_costs['Utilities'] = sum(new_utility_costs) - sum(original_utility_costs)
    
    @property
    def installation_cost(self):
        return self._installation_cost
    
    def _design(self): pass
    