#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 07:16:33 2020

@author: yalinli_cabbi
"""

'''
Note: moved from units.py to streamline the code and reduce confusion
'''


# %%

# =============================================================================
# Old flash unit modified from bst.units.Flash, add agitator,
# not used as the cost is too low compared to Humbird
# =============================================================================

# @cost(basis='Flow rate', ID='Agitator', units='kg/hr',	
#       kW=170, cost=90000, S=252891, CE=521.9, n=0.5, BM=1.5)	
# class OrganicAcidsFlash(Flash):	
    	
#     def _design(self):
#         super()._design()
#         self._decorated_design()
    	
#     def _cost(self):	
#         super()._cost()
#         rate = self.power_utility.rate
#         self._decorated_cost()
#         if self.P and self.P < 101325:	
#             self.power_utility += rate


# %%

# =============================================================================
# Old codes in current Esterification unit
# =============================================================================

    # def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
    #               T=351.15, P=101325, tau=None, tau_max=5, 
    #               V_wf=0.8, length_to_diameter=2, kW_per_m3=0.985,
    #               cat_load=0.039,  X1=None, X2=None, ethanol2LA=1.2, 
    #               assumeX2equalsX1=True, allow_higher_T=False,
    #               vessel_material='Stainless steel 316',
    #               vessel_type='Vertical'):
        
    #     Unit.__init__(self, ID, ins, outs)
        
    #     # self._methanol_composition = chemicals.kwarray(
    #     #         dict(Methanol=1-catalyst_molfrac,
    #     #              NaOCH3=catalyst_molfrac))
    #     self.T = T #: Operating temperature (K)
    #     self.P = P
    #     self.tau_max = tau_max
    #     self.V_wf = V_wf
    #     self.length_to_diameter = length_to_diameter
    #     self.kW_per_m3 = kW_per_m3
    #     self.cat_load = cat_load # 2.4%
    #     self.X1, self.X2, self._tau = X1, X2, tau
    #     self.ethanol2LA = ethanol2LA # 1.134 #!!! isn't the default 1.2?
    #     self.assumeX2equalsX1 = assumeX2equalsX1
    #     self.allow_higher_T = allow_higher_T
    #     self.vessel_material = vessel_material
    #     self.vessel_type = vessel_type
        
    #     self.heat_exchanger = HXutility(None, None, None, T=T)
    #     # self.heat_utilities = hx.heat_utilities

    # def compute_coefficients(self, T):
    #     K = self.K = exp(2.9625 - 515.13/T)
    #     kc = self.kc = 2.70 * (1e8) * exp(-6011.55/T)
    #     KW = self.KW = 15.19 * exp(12.01/T)
    #     KEt = self.KEt = 1.22 * exp(359.63/T)
    #     return K, kc, KW, KEt

    # def compute_X1_and_tau(self, time_step=1): # time_step in min
    #     T = self.T
    #     lle_chemicals = self.chemicals.lle_chemicals
    #     reactives = self.reactives
    #     # LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
    #     #     lle_chemicals.get_index(reactives())
    #     K, kc, KW, KEt = self.compute_coefficients(T)

    #     cat_load = self.cat_load
    #     feed, ethanol1, recycled_LA, supplement_ethanol, ethanol2 = self.ins
        
    #     time_max = self.tau_max * 60
    #     # ethanol = ethanol1.copy()
    #     # ethanol.mix_from([ethanol,ethanol2])

    #     # lle_chemicals = self.ins[0].lle_chemicals
    #     IDs = tuple([i.ID for i in lle_chemicals])
        
    #     f_gamma = tmo.equilibrium.DortmundActivityCoefficients(lle_chemicals)
        
    #     temp_flow = feed.copy() 
    #     # recycle_LA = recycled_LA.copy()
    #     # temp_flow.mix_from([temp_flow, ])
        
    #     # temp_flow.imol['Ethanol'] = ethanol_reqd
    #     # temp_flow.mol += self.ins[1].mol.copy()
    #     ethanol = self.ethanol
    #     temp_flow.mix_from([temp_flow, recycled_LA, ethanol, supplement_ethanol])
    #     # temp_flow.imol['Ethanol'] = ethanol_reqd
    #     #mcat = cat_load * (self.ins[0].F_mass + ethanol_reqd * self.ins[0].chemicals.LacticAcid.MW)
    #     mcat = cat_load * temp_flow.F_mass
    #     self.mcat = mcat
    #     tau = time_step
    #     LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
    #         [IDs.index(ID) for ID in reactives]
            
    #     gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
    #     gammas_reactives = np.array([gammas[LA_index], 
    #                                   gammas[EtOH_index], 
    #                                   gammas[H2O_index], 
    #                                   gammas[EtLA_index]])
            
    #     curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 
    #                                                 'Water', 'EthylLactate'))
    #     normalized_mol = temp_flow.get_normalized_mol(IDs)
    #     curr_conc = np.array([normalized_mol[LA_index], 
    #                           normalized_mol[EtOH_index], 
    #                           normalized_mol[H2O_index], 
    #                           normalized_mol[EtLA_index]])
    #     activities = gammas_reactives * curr_conc
    #     # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
    #     r = kc * (activities[1]*activities[0]-(activities[3]*activities[2]/K)) \
    #         /(1+KEt*activities[3]+KW*activities[2])**2
    #     dX = r * time_step * mcat / 1000 # r is in mol g-1 min-1
        
    #     new_flows = [1, 1, 1, 1]
    #     LA_initial = temp_flow.imol['LacticAcid']
        
    #     # print(LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index)
    #     # print(lle_chemicals)
    #     # # print(IDs)
    #     # print(temp_flow.get_normalized_mol(IDs))
    #     # print(temp_flow.mol)
    #     # print(IDs,gammas)
    #     # print(activities)
    #     # print(r)
    #     # print(dX)
    #     while dX/LA_initial>1e-4:

    #         # print('tau = %s' % tau)
    #         # # print(curr_flow)
    #         # # print(curr_conc)
    #         # print(r)
    #         # print(curr_flow)
    #         # print(dX)
    #         # print(activities)
    #         # print(new_flows)
            
    #         if curr_flow[0] < dX or curr_flow[1] < dX:
    #             dX = min(curr_flow[0], curr_flow[1])
                
    #         new_flows = [curr_flow[0]-dX, 
    #                       curr_flow[1]-dX, 
    #                       curr_flow[2]+dX, 
    #                       curr_flow[3]+dX]
            
    #         temp_flow.set_flow(new_flows, 'kmol/hr', ('LacticAcid', 'Ethanol', 
    #                                                   'Water', 'EthylLactate'))
            
    #         # Zhao et al. 2008 reported 96% conversion of NH4LA -> BuLA in 6h
    #         if new_flows[0]<=0 or new_flows[1]<=0 or tau>time_max-time_step: 
    #             # dX = min(new_flows)
    #             break
            
    #         # print(new_flows)
            
    #         curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 
    #                                                     'Water', 'EthylLactate'))
    #         normalized_mol = temp_flow.get_normalized_mol(IDs)
    #         curr_conc = np.array([normalized_mol[LA_index], 
    #                               normalized_mol[EtOH_index], 
    #                               normalized_mol[H2O_index], 
    #                               normalized_mol[EtLA_index]])
    #         gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
    #         gammas_reactives = np.array([gammas[LA_index], 
    #                                       gammas[EtOH_index], 
    #                                       gammas[H2O_index], 
    #                                       gammas[EtLA_index]])
    #         activities = gammas_reactives * curr_conc
            
    #         # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
    #         r = kc * (activities[1]*activities[0]-(activities[3]*activities[2]/K))/(1+KEt*activities[3]+KW*activities[2])**2
    #         # print(LA_initial)
    #         dX = r * time_step * mcat / 1000  # r is in mol g-1 min-1
            
    #         tau += time_step
        
    #     LA_in_feeds = feed.imol['LacticAcid']+recycled_LA.imol['LacticAcid']
    #     X1 = (LA_in_feeds-temp_flow.imol['LacticAcid']) / LA_in_feeds
    #     self._tau = tau / 60 # convert min to hr
    #     self.X1 = X1
    #     return X1, tau


# %%

# =============================================================================
# Unused Hydrolysis class
# =============================================================================

# # Cost copied from Transesterification
# @cost('Volume', 'Reactor', 
#       CE=525.4, cost=15000, n=0.55, kW=1.5, BM=4.3, ub = 4000, N = 'N_reactors')
# # @cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
# #       cost=85000, S=-8, CE=550.8, n=0.7, BM=2.2)
# class Hydrolysis(Unit):
#     """
#     Create an esterification reactor that converts 'LacticAcid', 'AceticAcid', and 'Ethanol'
#     to 'EthylLactate', 'EthylAcetate', and 'Water'. Finds the amount of catalyst 'Amberlyst-15'
#     required and consumes it.
    
#     Parameters
#     ----------
#     ins : stream sequence
#         * [0] Organic acids feed
#         * [1] Ethanol feed (includes catalyst)
#     outs : stream
#         Reactor effluent.
#     efficiency : float
#         Efficiency of conversion (on a 'LacticAcid' basis).
#     ethanol2LA : float
#         Ethanol feed to LacticAcid molar ratio.
#     T : float
#         Operating temperature [K].
    
#     """
#     _bounds = {'Volume': (0.1, 20)}
#     _units = {'Volume': 'm^3'}
#     _tau = 1
#     _N_ins = 2
#     _N_outs = 1
#     _N_heat_utilities = 1

#     def _more_design_specs(self):
#         return (('Residence time', self._tau, 'hr'),
#                 ('Conversion efficiency', self.X1, ''),
#                 ('Working volume fraction', 0.8, ''))

#     def __init__(self, ID='', ins=None, outs=(), thermo = None, *, X1 = None, tau = None, X2 = None, cat_load=0.024, T=353.40, assumeX2equalsX1 = True):
#         Unit.__init__(self, ID, ins, outs, thermo)
#         #: [:class:`~thermosteam.ParallelReaction`] Transesterification and catalyst consumption reaction
        

#         # self._methanol_composition = chemicals.kwarray(
#         #         dict(Methanol=1-catalyst_molfrac,
#         #              NaOCH3=catalyst_molfrac))
#         self.reactives = \
#                 ('LacticAcid', 'Ethanol', 'EthylLactate', 'H2O', 'AceticAcid', 'EthylAcetate')
#         self.assumeX2equalsX1 = assumeX2equalsX1
#         self.X1, self.X2, self._tau = X1, X2, tau
        
#         self.T = T #: Operating temperature (K). # 353.40 K
        
#         # self.ethanol2LA = ethanol2LA # 1.134
#         # self.ethanol_reqd = 
        
#         self.cat_load = cat_load # 2.4%
#          # kg
        
#         self.K = exp(2.9625 - 515.13/T)
#         self.kc = 2.70 * (1e8) * exp(-6011.55/T)
#         self.KW = 15.19 * exp(12.01/T)
#         self.KEt = 1.22 * exp(359.63/T)
        
#         self.heat_exchanger = hx = HXutility(None, None, None, T=T)
#         self.heat_utilities = hx.heat_utilities
        
        
#     def compute_X1_and_tau(self, time_step=450):
#         lle_chemicals = self.chemicals.lle_chemicals
#         reactives = self.reactives
#         # LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
#         #     lle_chemicals.get_index(reactives())
#         K = self.K
#         KW = self.KW
#         KEt = self.KEt
#         kc = self.kc
#         T = self.T
#         cat_load = self.cat_load
#         # ethanol2LA = self.ethanol2LA
#         # self.ethanol_reqd = ethanol2LA * self.ins[0].imol['LacticAcid']
#         # ethanol_reqd = self.ethanol_reqd
#         mcat = cat_load * self.ins[0].F_mass
        
        
#         # lle_chemicals = self.ins[0].lle_chemicals
#         IDs = tuple([i.ID for i in lle_chemicals])
        
#         f_gamma = tmo.equilibrium.DortmundActivityCoefficients(lle_chemicals)
        
        
#         temp_flow = self.ins[0].copy()
#         # temp_flow.imol['Ethanol'] = ethanol_reqd
#         tau = time_step
#         LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
#             [IDs.index(ID) for ID in reactives]
            

#         gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
#         gammas_reactives = np.array([gammas[LA_index], gammas[EtOH_index], gammas[H2O_index], gammas[EtLA_index]])
            
#         curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
#         normalized_mol = temp_flow.get_normalized_mol(IDs)
#         curr_conc = np.array([normalized_mol[LA_index], normalized_mol[EtOH_index], normalized_mol[H2O_index], normalized_mol[EtLA_index]])
#         activities = gammas_reactives*curr_conc
#         # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
#         r = kc* (activities[1]*activities[0] - (activities[3]*activities[2]/K))/(1 + KEt*activities[3] + KW*activities[2])**2
#         dX = r*time_step*mcat/1000
        
#         new_flows = [1, 1, 1, 1]
#         EtLA_initial = temp_flow.imol['EthylLactate']
        
#         print(LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index)
#         print(lle_chemicals)
#         # print(IDs)
#         print(temp_flow.get_normalized_mol(IDs))
#         print(temp_flow.mol)
#         print(IDs,gammas)
#         print(dX)
#         while dX/EtLA_initial<-1e-5:
            
            
#             print('tau = %s' % tau)
#             print(curr_flow)
#             print(curr_conc)
#             print(r)
#             print(dX)
#             print(new_flows)
#             new_flows = [curr_flow[0] - dX, curr_flow[1] - dX, curr_flow[2] + dX, curr_flow[3] + dX]
#             if new_flows[0]<=0 or new_flows[1]<=0 or tau>7200-time_step: # Zhao et al. 2008 reported 96% conversion of NH4LA -> BuLA in 6h
#                 break
#             print(new_flows)
#             temp_flow.set_flow(new_flows, 'kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
            
            
#             activities = gammas_reactives*curr_conc
#             # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
#             r = kc* (activities[1]*activities[0] - (activities[3]*activities[2]/K))/(1 + KEt*activities[3] + KW*activities[2])**2
#             # print(LA_initial)
#             dX = r*time_step*mcat/1000
#             curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
#             normalized_mol = temp_flow.get_normalized_mol(IDs)
#             curr_conc = np.array([normalized_mol[LA_index], normalized_mol[EtOH_index], normalized_mol[H2O_index], normalized_mol[EtLA_index]])
#             gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
#             gammas_reactives = np.array([gammas[LA_index], gammas[EtOH_index], gammas[H2O_index], gammas[EtLA_index]])
#             tau += time_step
            
#         X1 = (self.ins[0].imol['EthylLactate'] - temp_flow.imol['EthylLactate']) / self.ins[0].imol['EthylLactate']
#         self._tau = tau
#         self.X1 = X1
#         return X1, tau

#     @property
#     def tau(self):
#         """Residence time (hr)."""
#         return self._tau
#     # @tau.setter
#     # def tau(self, tau):
#     #     self._tau = tau
    
#     # @property
#     def efficiency(self):
#         """Esterification conversion efficiency."""
#         return self.reaction.X[0]
#     # @efficiency.setter
#     # def efficiency(self, efficiency):
#     #     self.reaction.X[0] = efficiency
    
#     # @property
#     # def methanol2lipid(self):
#     #     """Methanol feed to lipid molar ratio."""
#     #     return self._methanol2lipid
#     # @methanol2lipid.setter
#     # def methanol2lipid(self, ratio):
#     #     self._methanol2lipid = ratio

#     # @property
#     # def catalyst_molfrac(self):
#     #     """Catalyst molar fraction in methanol feed."""
#     #     return self._methanol_composition[self._catalyst_index]
#     # @catalyst_molfrac.setter
#     # def catalyst_molfrac(self, molfrac):
#     #     meoh = self._methanol_composition
#     #     meoh[self._catalyst_index] = molfrac
#     #     meoh[self._methanol_index] = 1-molfrac

#     def _run(self):
#         if self.X1==None:
#             X1, tau = self.compute_X1_and_tau()
#         else:
#             X1, tau = self.X1, self._tau
#         if self.assumeX2equalsX1:
#             X2 = X1
#             self.X2 = X2
#         else:
#             X2 = self.X2
        
#         self.reaction = ParallelRxn([
#             #   Reaction definition                               Reactant  Conversion
#             Rxn('EthylLactate + H2O -> LacticAcid + Ethanol', 'EthylLactate', X1),
#             Rxn('EthylAcetate + H2O -> AceticAcid + Ethanol', 'EthylAcetate', X2)
#             ])
#         feed, ethanol = self.ins
#         product, = self.outs
#         # ethanol.imol['Ethanol'] = self.ethanol_reqd
        
#         product.mol[:] = feed.mol + ethanol.mol
#         self.reaction(product.mol)
#         product.T = self.T
        
#     # def _design(self):
#     #     effluent = self._outs[0]
#     #     self.design_results['Volume'] = self._tau * effluent.F_vol / 0.8
#     #     self.heat_utilities[0](self.Hnet, effluent.T)
#     def _design(self):
        
#         product = self.outs[0]
#         super()._design()
        
#         self.design_results['Volume'] = (self._tau/60) * product.F_vol / 0.8
#         # self.heat_utilities[0](self.Hnet, product.T)
        
#     def _end_decorated_cost_(self):
#         # super()._cost()
#         N = self.design_results['N_reactors']
#         # self.heat_exchanger.simulate_as_auxiliary_exchanger(self.Hnet, self.outs[0])
#         # self.purchase_costs['Heat exchanger'] = self.heat_exchanger.purchase_costs['Heat exchanger'] 
        
#         hx_effluent = self.outs[0].copy()
#         hx_effluent.phase = 'l'
#         hx_effluent.mol[:] /= N
#         heat_exchanger = self.heat_exchanger
#         heat_exchanger.simulate_as_auxiliary_exchanger(self.Hnet/N, hx_effluent)
#         hu_bioreactor, = self.heat_utilities
#         hu_hx, = heat_exchanger.heat_utilities
#         hu_bioreactor.copy_like(hu_hx)
#         hu_bioreactor.scale(N)
#         self.purchase_costs['Heat exchangers'] = self.heat_exchanger.purchase_costs['Heat exchanger'] * N  
        

# %% 

# =============================================================================
# Old esterification codes
# =============================================================================

# # Cost copied from PretreatmentFlash
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=511000, S=264116, CE=521.9, n=0.7, BM=2)
# @cost(basis='Flow rate', ID='Agitator', units='kg/hr',
#       kW=170, cost=90000, S=252891, CE=521.9, n=0.5, BM=1.5)
# @cost(basis='Flow rate', ID='Pump', units='kg/hr',
#       kW=55.9275, cost=30000, S=204390, CE=521.9, n=0.8, BM=2.3)
# class EsterificationReactor(Unit):
#     _N_ins = 3
#     _N_outs = 1

#     def __init__(self, ID='', ins=None, outs=()):
#         Unit.__init__(self, ID, ins, outs)
#         # self.chems= chemicals
#         self.esterification_rxns = ParallelRxn([
#             #   Reaction definition                               Reactant  Conversion
#             Rxn('LacticAcid + Methanol -> MethylLactate + H2O', 'LacticAcid', 0.95),
#             Rxn('AceticAcid + Methanol -> MethylAcetate + H2O', 'AceticAcid', 0.95)
#             ])

#     def _run(self):
#         feed, recycled_methanol, supplement_methanol = self.ins
#         effluent = self.outs[0]
        
        
#         # chemicals = self.chems
#         # chemicals = effluent.Chemicals
#         # f_gamma = tmo.equilibrium.activity_coefficients.UNIFACActivityCoefficiencts(chemicals)
        
#         # gammas = np.ones(len(chemicals))
#         # err = 1
        
#         # # LA_index = chemicals.index('LacticAcid')
#         # # EtOH_index = chemicals.index('Ethanol')
#         # # EtLA_index = chemicals.index('EthylLactate')
#         # # H2O_index = chemicals.index('Water')
        
#         # LA_index = 1
#         # EtOH_index = 2
#         # EtLA_index = 3
#         # H2O_index = 0
#         # T = effluent.T
#         # init = 0.5*effluent.imol['LacticAcid']/effluent.F_mol
#         # K_eq = 2.9625 - 515.13/T # Thermodynamic Equilibrium and Reaction Kinetics for the Esterification of
#         #                                     # Lactic Acid with Ethanol Catalyzed by Acid Ion-Exchange Resin
#         #                                     # (2008) Carla S. M. Pereira,† Sima˜o P. Pinho,‡ Viviana M. T. M. Silva,§ and Alı´rio E. Rodrigues*,†
#         # while (err>1e-4):
#         #     # EtLA * (Water + EtLA) = K_eq/K_gamma
#         #     K_gamma = ((gammas[EtLA_index]*gammas[H2O_index])/(gammas[LA_index]*gammas[EtOH_index]))
#         #     water_0 = effluent.imol['Water']/effluent.F_mol
#         #     balance = lambda EtLA: K_eq/K_gamma - EtLA*(water_0 + EtLA)
            
#         #     EtLA = fsolve(balance, [init, init, init])
#         #     print(EtLA)
#         #     effluent.imol['LacticAcid'] -= EtLA*effluent.F_mol
#         #     effluent.imol['Ethanol'] -= EtLA*effluent.F_mol
#         #     effluent.imol['EthylLactate'] += EtLA*effluent.F_mol
#         #     effluent.imol['Water'] += EtLA*effluent.F_mol
            
#         #     gammas = f_gamma(effluent.mol[:])
            
            
#         # Add 5% extra
#         methanol_needed = (feed.imol['LacticAcid']/self.esterification_rxns.X[0] \
#                             +feed.imol['AceticAcid']/self.esterification_rxns.X[1]) \
#                           * 1.05 - (feed.imol['Methanol']+recycled_methanol.imol['Methanol'])
#         supplement_methanol.imol['Methanol'] = max(0, methanol_needed)
#         effluent.mix_from([feed, recycled_methanol, supplement_methanol])
#         self.esterification_rxns(effluent.mol)




# # Cost copied from Transesterification
# @cost('Volume', 'Reactor',
#       CE=525.4, cost=15000, n=0.55, kW=1.5, BM=4.3,)

# class Esterification(Unit):
#     """
#     Create an esterification reactor that converts 'LacticAcid', 'AceticAcid', and 'Ethanol'
#     to 'EthylLactate', 'EthylAcetate', and 'Water'. Finds the amount of catalyst 'Amberlyst-15'
#     required and consumes it.
    
#     Parameters
#     ----------
#     ins : stream sequence
#         * [0] Organic acids feed
#         * [1] Ethanol feed (includes catalyst)
#     outs : stream
#         Reactor effluent.
#     efficiency : float
#         Efficiency of conversion (on a 'LacticAcid' basis).
#     ethanol2LA : float
#         Ethanol feed to LacticAcid molar ratio.
#     T : float
#         Operating temperature [K].
    
#     """
#     _bounds = {'Volume': (0.1, 20)}
#     _units = {'Volume': 'm^3'}
#     _tau = 1
#     _N_ins = 3
#     _N_outs = 1
#     _N_heat_utilities = 1

#     def _more_design_specs(self):
#         return (('Residence time', self._tau, 'hr'),
#                 ('Conversion efficiency', self.X1, ''),
#                 ('Working volume fraction', 0.8, ''))

#     def __init__(self, ID='', ins=None, outs=(), thermo = None, *, X1 = None, tau = None, X2 = None, ethanol2LA=1.8, cat_load=0.024, T=353.40, assumeX2equalsX1 = True):
#         Unit.__init__(self, ID, ins, outs, thermo)
#         #: [:class:`~thermosteam.ParallelReaction`] Transesterification and catalyst consumption reaction
        

#         # self._methanol_composition = chemicals.kwarray(
#         #         dict(Methanol=1-catalyst_molfrac,
#         #              NaOCH3=catalyst_molfrac))
#         self.reactives = \
#                 ('LacticAcid', 'Ethanol', 'EthylLactate', 'H2O', 'AceticAcid', 'EthylAcetate')
#         self.assumeX2equalsX1 = assumeX2equalsX1
#         self.X1, self.X2, self._tau = X1, X2, tau
        
#         self.T = T #: Operating temperature (K). # 353.40 K
        
#         self.ethanol2LA = ethanol2LA # 1.134
#         # self.ethanol_reqd = 
        
#         self.cat_load = cat_load # 2.4%
#          # kg
        
#         self.K = exp(2.9625 - 515.13/T)
#         self.kc = 2.70 * (1e8) * exp(-6011.55/T)
#         self.KW = 15.19 * exp(12.01/T)
#         self.KEt = 1.22 * exp(359.63/T)
    
#     def compute_X1_and_tau(self, time_step=1): # time_step in min
#         lle_chemicals = self.chemicals.lle_chemicals
#         reactives = self.reactives
#         # LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
#         #     lle_chemicals.get_index(reactives())
#         K = self.K
#         KW = self.KW
#         KEt = self.KEt
#         kc = self.kc
#         T = self.T
#         cat_load = self.cat_load
#         ethanol2LA = self.ethanol2LA
#         feed, ethanol, supplement_ethanol = self.ins
#         self.ethanol_reqd = ethanol2LA * (feed.imol['LacticAcid'])
#         ethanol_reqd = self.ethanol_reqd
        
        
        
        
#         # lle_chemicals = self.ins[0].lle_chemicals
#         IDs = tuple([i.ID for i in lle_chemicals])
        
#         f_gamma = tmo.equilibrium.DortmundActivityCoefficients(lle_chemicals)
        
        
#         temp_flow = feed.copy() 
#         # temp_flow.mol += + recycled_LA.mol
#         # temp_flow.imol['Ethanol'] = ethanol_reqd
#         # temp_flow.mol += self.ins[1].mol.copy()
#         temp_flow.imol['Ethanol'] = ethanol_reqd
#         #mcat = cat_load * (self.ins[0].F_mass + ethanol_reqd * self.ins[0].chemicals.LacticAcid.MW)
#         mcat = cat_load * temp_flow.F_mass
        
#         tau = time_step
#         LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
#             [IDs.index(ID) for ID in reactives]
            
                
#         gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
#         gammas_reactives = np.array([gammas[LA_index], gammas[EtOH_index], gammas[H2O_index], gammas[EtLA_index]])
            
#         curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
#         normalized_mol = temp_flow.get_normalized_mol(IDs)
#         curr_conc = np.array([normalized_mol[LA_index], normalized_mol[EtOH_index], normalized_mol[H2O_index], normalized_mol[EtLA_index]])
#         activities = gammas_reactives*curr_conc
#         # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
#         r = kc* (activities[1]*activities[0] - (activities[3]*activities[2]/K))/(1 + KEt*activities[3] + KW*activities[2])**2
#         dX = r*time_step*mcat/(1000) # r is in mol g-1 min-1
        
#         new_flows = [1, 1, 1, 1]
#         LA_initial = temp_flow.imol['LacticAcid']
        
#         print(LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index)
#         print(lle_chemicals)
#         # print(IDs)
#         print(temp_flow.get_normalized_mol(IDs))
#         print(temp_flow.mol)
#         print(IDs,gammas)
#         print(activities)
#         print(r)
#         print(dX)
#         while dX/LA_initial>1e-4:
            
            
#             print('tau = %s' % tau)
#             # print(curr_flow)
#             # print(curr_conc)
#             print(r)
#             print(curr_flow)
#             print(dX)
#             print(activities)
#             # print(new_flows)
            
#             if curr_flow[0]<dX or curr_flow[1]<dX:
#                 dX = min(curr_flow[0], curr_flow[1])
                
#             new_flows = [curr_flow[0] - dX, curr_flow[1] - dX, curr_flow[2] + dX, curr_flow[3] + dX]
#             temp_flow.set_flow(new_flows, 'kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
            
#             if new_flows[0]<=0 or new_flows[1]<=0 or tau>360-time_step: # Zhao et al. 2008 reported 96% conversion of NH4LA -> BuLA in 6h
#                 # dX = min(new_flows)
#                 break
            
#             # print(new_flows)
            
            
#             curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
#             normalized_mol = temp_flow.get_normalized_mol(IDs)
#             curr_conc = np.array([normalized_mol[LA_index], normalized_mol[EtOH_index], normalized_mol[H2O_index], normalized_mol[EtLA_index]])
#             gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
#             gammas_reactives = np.array([gammas[LA_index], gammas[EtOH_index], gammas[H2O_index], gammas[EtLA_index]])
#             activities = gammas_reactives*curr_conc
            
#             # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
#             r = kc* (activities[1]*activities[0] - (activities[3]*activities[2]/K))/(1 + KEt*activities[3] + KW*activities[2])**2
#             # print(LA_initial)
#             dX = r*time_step*mcat/(1000) # r is in mol g-1 min-1
            
#             tau += time_step
            
#         X1 = (feed.imol['LacticAcid'] - temp_flow.imol['LacticAcid']) / feed.imol['LacticAcid']
#         self._tau = tau
#         self.X1 = X1
#         return X1, tau

#     @property
#     def tau(self):
#         """Residence time (hr)."""
#         return self._tau
#     # @tau.setter
#     # def tau(self, tau):
#     #     self._tau = tau
    
#     # @property
#     def efficiency(self):
#         """Esterification conversion efficiency."""
#         return self.reaction.X[0]
#     # @efficiency.setter
#     # def efficiency(self, efficiency):
#     #     self.reaction.X[0] = efficiency
    
#     # @property
#     # def methanol2lipid(self):
#     #     """Methanol feed to lipid molar ratio."""
#     #     return self._methanol2lipid
#     # @methanol2lipid.setter
#     # def methanol2lipid(self, ratio):
#     #     self._methanol2lipid = ratio

#     # @property
#     # def catalyst_molfrac(self):
#     #     """Catalyst molar fraction in methanol feed."""
#     #     return self._methanol_composition[self._catalyst_index]
#     # @catalyst_molfrac.setter
#     # def catalyst_molfrac(self, molfrac):
#     #     meoh = self._methanol_composition
#     #     meoh[self._catalyst_index] = molfrac
#     #     meoh[self._methanol_index] = 1-molfrac

#     def _run(self):
#         if self.X1==None:
#             X1, tau = self.compute_X1_and_tau()
#         else:
#             X1, tau = self.X1, self._tau
#         if self.assumeX2equalsX1:
#             X2 = X1
#             self.X2 = X2
#         else:
#             X2 = self.X2
            
#         self.reaction = ParallelRxn([
#             #   Reaction definition                               Reactant  Conversion
#             Rxn('LacticAcid + Ethanol -> EthylLactate + H2O', 'LacticAcid', X1),
#             Rxn('AceticAcid + Ethanol -> EthylAcetate + H2O', 'AceticAcid', X2)
#             ])
#         feed, ethanol, supplement_ethanol = self.ins # broth, recycled ethanol, recycled LA, supplementary ethanol
#         product, = self.outs
        
#         self.added_ethanol = max(0,self.ethanol_reqd - ethanol.imol['Ethanol'])
        
#         supplement_ethanol.imol['Ethanol'] = max(0,  self.added_ethanol)
        
#         product.imol['Ethanol'] += self.added_ethanol
#         product.mol[:] = feed.mol + ethanol.mol + supplement_ethanol.mol
#         self.reaction(product.mol)
#         product.T = self.T
        
#     def _design(self):
#         effluent = self._outs[0]
#         self.design_results['Volume'] = self._tau * effluent.F_vol / 0.8
#         self.heat_utilities[0](self.Hnet, effluent.T)

# %%

# =============================================================================
# Old LiquidLiquidExtraction
# =============================================================================

# class LiquidLiquidExtraction(Unit):
#     _N_ins = 2
#     _N_outs = 2
    
#     def __init__(self, ID='', ins=None, outs=(), *, solute, feed_solvent, extraction_solvent, solute_target_conc):
#         Unit.__init__(self, ID, ins, outs)
#         self.solute = solute
#         self.feed_solvent = feed_solvent
#         self.extraction_solvent = extraction_solvent
#         self.Ye = solute_target_conc
        
#     def calculate_number_of_stages():
        
#     def _run(self):
#         water, vent_entry = self.ins
#         vent_exit, bottoms = self.outs
        
#         vent_exit.copy_like(vent_entry)
#         bottoms.copy_flow(vent_exit, self.gas, remove=True, exclude=True)
#         bottoms.mol += water.mol
        
#     def _design(self):
#         Design = self.design_results
#         Design['Vent flow rate'] = self.ins[1].F_mass
#         Design['Bottom flow rate'] = self.outs[1].F_mass