#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 09:53:03 2020

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

@author: yalinli_cabbi
"""

'''
TODO:
    Might want a better way to deal with insufficient BT steam generation
'''



# %% Setup

from biosteam import HeatUtility, Facility
from biosteam.units.decorators import cost
from biosteam.units import HXprocess
from thermosteam import Stream
import biosteam as bst
import copy

from hx_network_design import hx_network_design

Temperature_Interval_Method = hx_network_design.Temperature_Interval_Method
Design_Network = hx_network_design.Design_Network
# from orgacids.system import orgacids_sys

__all__ = ('OrganicAcidsCIP', 'OrganicAcidsADP', 'OrganicAcidsCT', 'OrganicAcidsPWC', 
           'OrganicAcidsBT')

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
# Year  1997    1998    2009    2010    2016
# CE    386.5   389.5   521.9   550.8   541.7


# %% Clean-in-place package

@cost(basis='Flow rate', ID='Clean-in-place package', units='kg/hr',
      cost=421000, S=63, CE=521.9, BM=1.8, n=0.6)
class OrganicAcidsCIP(Facility):
    network_priority = 2
    line = 'Clean-in-place package'


# %% Air distribution package

# Size based on stream 950 in Humbird et al.
@cost(basis='Flow rate', ID='Plant air compressor', units='kg/hr',
      kW=111.855, cost=28000, S=1372608, CE=550.8, n=0.6, BM=1.6)
@cost(basis='Flow rate', ID='Plant air reciever', units='kg/hr',
      cost=16000, S=1372608, CE=521.9, n=0.6, BM=3.1)
@cost(basis='Flow rate', ID='Instrument air dryer', units='kg/hr',
      cost=15000, S=1372608, CE=521.9, n=0.6, BM=1.8)
class OrganicAcidsADP(Facility): 
    line = 'Air distribution package'
    network_priority = 2
    

# %% Cooling tower

@cost('Flow rate', 'Cooling tower', units= 'kg/hr',
      kW=559.275, cost=1375000, S=10037820, CE=550.8, n=0.6, BM=1.5)
@cost('Flow rate', 'Cooling water pump', units='kg/hr',
      kW=1118.55, cost=283671, S=10982556,  CE=550.8, n=0.8, BM=3.1)
class OrganicAcidsCT(Facility):
    """
    Create a cooling tower with capital cost and power based on the flow rate 
    of cooling water as in [1]_.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    network_priority = 1
    _N_ins = 2
    _N_outs = 2    
    _N_heat_utilities = 1

    # Page 55 of Humbird et al.
    blowdown = 0.00005+0.0015
    
    def __init__(self, ID='', ins=None, outs=()):
        Facility.__init__(self, ID, ins, outs)
        # self.system_cooling_utilities = system_cooling_utilities
        
        
    def _run(self):
        return_cw, makeup_water = self.ins
        process_cw, blowdown = self.outs
        system_cooling_utilities = self.system_cooling_utilities = {}

        # Based on stream 945 in Humbird et al.
        return_cw.T = 37 + 273.15
        # Based on streams 940/944 in Humbird et al.
        process_cw.T = blowdown.T = 28 + 273.15
        
        for u in self.system.units:
            if u is self: continue
            if hasattr(u, 'heat_utilities'):
                for cw in u.heat_utilities:
                    # hu.flow>0 is to avoid having BT's steam counted here
                    # if cw.duty <0 and cw.flow >0:
                    if cw.duty * cw.flow <0:
                        # Note that technically, CT can only produce cooling water,
                        # but here includes all cool_utilities to check if
                        # cooling water is the only cooling utility in the system
                        system_cooling_utilities[f'{u.ID} - {cw.ID}'] = cw
        
        hu_cooling = HeatUtility().sum_by_agent(system_cooling_utilities.values())[0]
        hu_cooling.reverse()
        self.hu_cooling = hu_cooling
        self.heat_utilities = (hu_cooling,)
        self.system_cooling_demand = -hu_cooling.duty
        
        # Total amount of cooling water needed in the whole system (kmol/hr)
        total_cooling_water = self.total_cooling_water = -hu_cooling.flow
        return_cw.imol['H2O'] = process_cw.imol['H2O'] = total_cooling_water
        makeup_water.imol['H2O'] = total_cooling_water * self.blowdown
        blowdown.imol['H2O'] = makeup_water.imol['H2O']

        self.design_results['Flow rate'] = total_cooling_water


# %% Process water center

@cost(basis='Total water flow rate', ID='Tank', units='kg/hr',
      # Size basis changed from 451555 as the original design is not sufficient
      CE=521.9, cost=250000, S=106453, n=0.7, BM=1.7)
@cost(basis='Total water flow rate', ID='PWC circulating pump', units='kg/hr',
      CE=550.8, kW=55.9275, cost=15292, S=518924, n=0.8, BM=3.1)
@cost(basis='Makeup/discharged water flow rate', ID='Makeup/discharged water pump', units='kg/hr',
      CE=550.8, kW=14.914, cost=6864, S=155564, n=0.8, BM=3.1)
class OrganicAcidsPWC(Facility):
    """
    Create a ProcessWaterCenter object that takes care of balancing the amount
    of water required for the process. The capital cost and power are based on 
    the flow rate of process and makeup water as in [1]_.
    
    Parameters
    ----------
    ins :
        [0] Recycled water.
        
        [1] Makeup water (>0 when recycled water < process water).
    outs :
        [0] Process water.
        
        [1] Discharged water (>0 when recycled water > process water).
    process_water_streams : streams, optional
        All process water streams.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    network_priority = 1
    _N_ins = 2
    _N_outs = 2
    _units= {'Total water flow rate': 'kg/hr',
             'Makeup/discharged water flow rate': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), process_water_streams=None):
        Facility.__init__(self, ID, ins, outs)
        self.process_water_streams = process_water_streams

    def _run(self):
        recycled, makeup = self.ins
        process, discharged = self.outs
        process_water_streams = self.process_water_streams        
        
        process.imol['Water'] = sum(i.imol['Water'] for i in process_water_streams)           
        discharged.imol['Water'] = (recycled.imol['Water']+makeup.imol['Water']) - \
            process.imol['Water']
        if discharged.imol['Water'] < 0:
            makeup.imol['Water'] -= discharged.imol['Water']
            discharged.imol['Water'] = 0
        else:
            makeup.imol['Water'] = 0
        total_stream = process.copy()
        total_stream.mix_from([process, makeup])
            
        Design = self.design_results
        Design['Total water flow rate'] = total_stream.F_mass
        Design['Makeup/discharged water flow rate'] = max(makeup.imass['Water'],
                                                          discharged.imass['Water'])
        

# %% Boiler and turbogenerator

@cost(basis='Flow rate', ID='Boiler', units='kg/hr',
      kW=2752, cost=28550000, S=238686, CE=550.8, n=0.6, BM=1.8)
@cost(basis='Work', ID='Turbogenerator', units='kW',
      cost=9500000, S=42200, CE=550.8, n=0.6, BM=1.8)
@cost(basis='Flow rate', ID='Hot process water softener system', units='kg/hr',
      cost=78000, S=235803, CE=550.8, n=0.6, BM=1.8)
@cost(basis='Flow rate', ID='Deaerator', units='kg/hr',
      cost=305000, S=235803, CE=550.8, n=0.6, BM=3)
@cost(basis='Flow rate', ID='Amine addition pkg', units='kg/hr',
      cost=40000, S=235803, CE=550.8, n=0, BM=1.8)
class OrganicAcidsBT(Facility):
    """
    Create a BoilerTurbogenerator object that will calculate electricity
    generation from burning the feed. It also takes into account how much
    steam is being produced, and the required cooling utility of the turbo
    generator. Combustion reactions are populated based on molecular formula
    of the ins. Purchase natural gas if BT cannot meet system heating/steam demand.
    
    Parameters
    ----------
    ins : 
        [0] Liquid/solid wastes to be burned.
        
        [1] Gas waste to be burned.
        
        [2] FGD lime.
        
        [3] Boiler chemicals.
        
        [4] Baghouse bag.
                
        [5] Supplementary natural gas in the case that waste energy is not enough
            to meet system steam demand

        [6] Makeup water.
        
    outs :
        [0] Gas emission.
        
        [1] Ash residues.
        
        [2] Blowdown water.
   
    boiler_efficiency : float
        Fraction of heat transfered to steam.
        
    turbo_generator_efficiency : float
        Fraction of steam heat converted to electricity.
        
    agent : UtilityAgent
        Steam produced.
        
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
        
    .. [2] Cortes-Peña, Y.; Kumar, D.; Singh, V.; Guest, J. S. 
        BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, 
        and Techno-Economic Analysis of Biorefineries under Uncertainty. 
        ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
        https://doi.org/10.1021/acssuschemeng.9b07040.
        
    """
    network_priority = 0
    _N_ins = 7
    _N_outs = 3
    _units= {'Flow rate': 'kg/hr',
             'Work': 'kW'} 
    
    blowdown = 0.03
    
    def __init__(self, ID='', ins=None, outs=(), *, B_eff=0.8,
                 TG_eff=0.85, combustibles, ratio, side_streams_to_heat=()):
        Facility.__init__(self, ID, ins, outs)
        self.B_eff = B_eff
        self.TG_eff = TG_eff
        self.combustibles = combustibles
        self.ratio = ratio
        self.side_streams_to_heat = side_streams_to_heat
        self.side_streams_lps = None

    def _run(self): pass


    def _design(self):
        feed_solids, feed_gases, lime, boiler_chems, bag, natural_gas, makeup_water \
            = self.ins
        emission, ash, blowdown_water = self.outs
        ratio = self.ratio
        side_streams_to_heat = self.side_streams_to_heat
        side_streams_lps = self.side_streams_lps
        system_heating_utilities = self.system_heating_utilities = {}
        lps = HeatUtility.get_heating_agent('low_pressure_steam')
        hps = HeatUtility.get_heating_agent('high_pressure_steam')
        
        # Use combustion reactions to create outs
        combustion_rxns = self.chemicals.get_combustion_reactions()
        combustable_feeds = Stream(None)
        emission.mol = feed_solids.mol + feed_gases.mol
        combustable_feeds.copy_flow(emission, tuple(self.combustibles), remove=True)
        combustion_rxns.force_reaction(combustable_feeds.mol)
        emission.mol += combustable_feeds.mol

        # FGD lime scaled based on SO2 generated,
        # 20% stoichiometetric excess based on P52 of Humbird et al.
        # 895 kg/hr based on Table 30 in Humbird et al.
        lime.imol['CalciumDihydroxide'] = emission.imol['SO2'] * 1.2
        # boiler_chems and bag are scaled based on plant size,
        # 1.23 is $2007/hour and 2.2661 is $/lb from Table 30 in Humbird et al.
        # 2.20462 is kg to lb
        boiler_chems.imass['BoilerChemicals'] = 1.23 / 2.2661 / 2.20462 * ratio
        bag.imass['BaghouseBag'] = ratio
        
        # 92% of SO2 removed by lime
        CaSO4_mol = emission.imol['SO2'] * 0.92
        emission.imol['SO2'] *= (1 - 0.92)

        # Air/O2 usage not rigorously modeled
        emission.imol['O2'] = 0
        
        ash.empty()
        for chemical in emission.chemicals:
            if chemical.locked_state == 'l' or chemical.locked_state == 's':
                ash.imol[chemical.ID] = emission.imol[chemical.ID]
                emission.imol[chemical.ID] = 0
                
        ash.mol += boiler_chems.mol
        ash.imol['CaSO4'] = CaSO4_mol
        ash.imol['Lime'] += lime.F_mol - CaSO4_mol
        
        emission.phase = 'g'
        ash.phase = 's'
        # Assume T of emission and ash are the same as hps, which
        # has highest T amont all heating agents
        emission.T = ash.T = hps.T

        # Total heat generated by the boiler (kJ/hr)
        H_in = feed_solids.H + feed_gases.H
        H_out = emission.H + ash.H
        # Water evaporation energy is already accounted for in emission enthalpy
        heat_from_combustion = -(feed_solids.HHV+feed_gases.HHV)
        heat_generated = self.heat_generated = \
            (H_in+heat_from_combustion)*self.B_eff - H_out
        
        for u in self.system.units:
            if u is self: continue
            if hasattr(u, 'heat_utilities'):
                for hu in u.heat_utilities:
                    # hu.flow>0 is to avoid having CT's cooling water counted here
                    # if hu.
                    # if hu.duty >0 and hu.flow >0:
                    if hu.duty * hu.flow >0:
                        system_heating_utilities[f'{u.ID} - {hu.ID}'] = hu
        
        # Use lps to account for the energy needed for the side steam
        if side_streams_to_heat:
            if not side_streams_lps:
                side_streams_lps = self.side_streams_lps = HeatUtility()
                side_streams_lps.load_agent(lps)
            side_streams_lps(duty=sum([i.H for i in side_streams_to_heat]), 
                             T_in=298.15)
            system_heating_utilities['BT - side_streams_lps'] = side_streams_lps

        system_heating_demand = self.system_heating_demand = \
            sum([i.duty for i in system_heating_utilities.values()])
            
        BT_heat_surplus = self.BT_heat_surplus = heat_generated - system_heating_demand
        self.system_steam_demand = sum([i.flow for i in system_heating_utilities.values()])
        
        # hu_BT = set()
        hu_cooling = HeatUtility()
        # BT can meet system heating/steam demand
        if BT_heat_surplus >0:
            # 3600 is conversion of kJ/hr to kW (kJ/s)
            electricity_generated = self.electricity_generated = \
                BT_heat_surplus * self.TG_eff / 3600
            
            # Take the opposite for cooling duty (i.e., cooling duty should be negative)
            # this is to condense the unused steam
            cooling_need = self.cooling_need = -(BT_heat_surplus-electricity_generated)
            # hu_cooling = HeatUtility()
            hu_cooling(duty=cooling_need, T_in=lps.T)
            natural_gas.empty()
            # hu_BT.add(hu_cooling)
            
        # BT cannot meet system heating/steam demand, supplement with natural gas
        else:
            CH4_LHV = natural_gas.chemicals.CH4.LHV
            natural_gas.imol['CH4'] = BT_heat_surplus / (CH4_LHV*self.B_eff)
            emission.imol['CO2'] += natural_gas.imol['CH4']
            emission.imol['H2O'] += 2 * natural_gas.imol['CH4']
            electricity_generated = self.electricity_generated = 0
            

        # for i in system_heating_utilities.values():
        #     j = HeatUtility()
        #     j.copy_like(i)
        #     j.reverse()
        #     hu_BT.add(j)
        heating_utilities = HeatUtility.sum_by_agent(system_heating_utilities.values())
        for i in heating_utilities:
            i.reverse()
        
        
        # heating_utilities.reverse()
        self.heat_utilities = [hu_cooling, *heating_utilities]
        
        self.heat_utilities = tuple(self.heat_utilities)
        
        total_steam = sum([i.flow for i in system_heating_utilities.values()])
        blowdown_water.imol['H2O'] = total_steam * self.blowdown
        blowdown_water.T = 373.15

        # Additional need for making lime slurry
        makeup_water.imol['H2O'] = blowdown_water.imol['H2O'] + lime.F_mol/0.2*0.8

        # BT_utilities = self.BT_utilities = HeatUtility.sum_by_agent(hu_BT)
        # self.heat_utilities = tuple(BT_utilities)
        Design = self.design_results        
        Design['Flow rate'] = total_steam
        Design['Work'] = electricity_generated

    def _end_decorated_cost_(self):
        self.power_utility(self.power_utility.rate - self.electricity_generated)
        

class HX_Network(Facility):
    # network_priority = 3
    network_priority = -1
    # _N_ins = 6
    # _N_outs = 3
    _N_heat_utilities = 1
    _units= {'Flow rate': 'kg/hr',
             'Work': 'kW'} 
    
    from biosteam import HeatUtility, settings, Stream
    
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
    def _design(self):pass
    