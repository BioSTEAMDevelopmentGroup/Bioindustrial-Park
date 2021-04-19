# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 03:19:53 2021

@author: yrc2
"""
# from biosteam import Unit
# from biosteam.units.decorators import cost
# from biosteam.units.design_tools import size_batch
# import thermosteam as tmo
# import biosteam as bst

# Rxn = tmo.reaction.Reaction
# ParallelRxn = tmo.reaction.ParallelReaction

# # %% Constants

# _gal2m3 = 0.003785
# _gpm2m3hr = 0.227124
# # _m3hr2gpm = 4.40287
# _hp2kW = 0.7457
# _Gcal2kJ = 4184e3

# # %% Consolidated bioprocess

# @cost('Flow rate', 'Recirculation pumps', kW=30, S=340*_gpm2m3hr,
#       cost=47200, n=0.8, BM=2.3, CE=522, N='N_recirculation_pumps')
# @cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900,
#       S=-5*_Gcal2kJ, n=0.7, BM=2.2, N='N_reactors') # Based on a similar heat exchanger
# @cost('Batch duty', 'Fermentor batch cooler', CE=522, cost=86928,
#       S=-5*_Gcal2kJ, n=0.7, BM=1.8) # Based on a similar heat exchanger
# @cost('Reactor volume', 'Agitators', CE=522, cost=52500,
#       S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_reactors')
# @cost('Reactor volume', 'Reactors', CE=522, cost=844000,
#       S=1e6*_gal2m3, n=0.5, BM=1.5, N='N_reactors')
# class SaccharificationAndCoFermentation(Unit):
#     _N_ins = 2
#     _N_outs = 2
#     _ins_size_is_fixed = False
#     _N_heat_utilities = 1
        
#     #: Saccharification temperature (K)
#     T_saccharification = 48+273.15
    
#     #: Fermentation temperature (K)
#     T_fermentation = 32+273.15
    
#     #: Saccharification time (hr)
#     tau_saccharification = 60
    
#     #: Co-Fermentation time (hr)
#     tau_cofermentation = 36
    
#     #: Unload and clean up time (hr)
#     tau_0 = 4
    
#     #: Working volume fraction (filled tank to total tank volume)
#     V_wf = 0.9
    
#     #: Number of reactors
#     N_reactors = 12
    
#     #: Number of recirculation pumps
#     N_recirculation_pumps = 5
    
#     _units = {'Flow rate': 'm3/hr',
#               'Reactor volume': 'm3',
#               'Batch duty': 'kJ/hr',
#               'Reactor duty': 'kJ/hr'}
    
#     def __init__(self, ID='', ins=None, outs=(), P=101325):
#         Unit.__init__(self, ID, ins, outs)
#         self.P = P
        
#         self.sucrose_hydrolysis_reaction = Rxn(
#             'Sucrose + Water -> 2Glucose', 'Sucrose', 1.00
#         )
        
#         self.loss = ParallelRxn([
#     #   Reaction definition               Reactant    Conversion
#     Rxn('Glucose -> 2 LacticAcid',       'Glucose',   0.0300),
#     Rxn('3 Xylose -> 5 LacticAcid',      'Xylose',    0.0300),
#     Rxn('3 Arabinose -> 5 LacticAcid',   'Arabinose', 0.0300),
#     Rxn('Galactose -> 2 LacticAcid',     'Galactose', 0.0300),
#     Rxn('Mannose -> 2 LacticAcid',       'Mannose',   0.0300),])
    
#         self.cofermentation = ParallelRxn([
#     #   Reaction definition                                          Reactant    Conversion
#     Rxn('Glucose -> 2 Ethanol + 2 CO2',                             'Glucose',   0.9500),
#     Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O', 'Glucose',   0.0200),
#     Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',                       'Glucose',   0.0040),
#     Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',                   'Glucose',   0.0060),
#     Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                            'Xylose',    0.8500),
#     Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
#                                                                     'Xylose',    0.0190),
#     Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',                  'Xylose',    0.0030),
#     Rxn('Xylose + H2O -> Xylitol + 0.5 O2',                         'Xylose',    0.0460),
#     Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',              'Xylose',    0.0090),
#     ])
    
#         self.CSL_to_constituents = Rxn(
#             'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL', 1.0000, basis='wt',
#         )
#         self.CSL_to_constituents.basis = 'mol'
        
#     def _run(self):
#         vent, effluent = self.outs
#         vent.P = effluent.P = self.P
#         vent.T = effluent.T = self.T_fermentation
#         vent.phase = 'g'
#         effluent.mix_from(self.ins, energy_balance=False)
#         self.loss(effluent)
#         self.sucrose_hydrolysis_reaction(effluent)
#         self.cofermentation(effluent)
#         self.CSL_to_constituents(effluent)
#         vent.receive_vent(effluent)

#     def _design(self):
#         effluent = self.outs[1]
#         v_0 = effluent.F_vol
#         Design = self.design_results
#         Design['Flow rate'] = v_0 / self.N_recirculation_pumps
#         tau = self.tau_saccharification + self.tau_cofermentation
#         Design.update(size_batch(v_0, tau, self.tau_0, self.N_reactors, self.V_wf))
#         hu_fermentation, = self.heat_utilities
#         ei = effluent.chemicals.index('Ethanol')
#         ethanol = (sum([i.mol[ei] for i in self.outs])
#                    - sum([i.mol[ei] for i in self.ins]))
#         Design['Batch duty'] = batch_duty = self.H_out - self.H_in
#         Design['Reactor duty'] = reactor_duty = self.Hf_out - self.Hf_in
#         hu_fermentation(reactor_duty + batch_duty, effluent.T)
   