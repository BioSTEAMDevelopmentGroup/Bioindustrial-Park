#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import numpy as np
import biosteam as bst
import thermosteam as tmo
import nskinetics as nsk

TelluriumReactionSystem = nsk.TelluriumReactionSystem
ReactionSystem = nsk.ReactionSystem
cost = bst.decorators.cost
Splitter = bst.Splitter

__all__ = ('SimultaneousSaccharificationFermentationEthanolIsobutanol',
           'SSFEtOHIBO', 'MolecularSieve')

class SimultaneousSaccharificationFermentationEthanolIsobutanol(bst.BatchBioreactor):
    """
    Create a SimultaneousSaccharificationFermentation unit operation that 
    models the simultaneous saccharification and fermentation in the conventional
    dry-grind ethanol process.
    
    Parameters
    ----------
    ins : streams
        Inlet fluids.
    outs : stream
        Outlet fluid.
    yield_: float
        Yield of glucose to ethanol as a fraction of the theoretical yield.
    
    Notes
    -----
    This unit operation doesn't actually model the saccharification process.
    The reactor is modeled by the stoichiometric conversion of glucose to
    ethanol by mol:
        
    .. math:: 
        Glucose -> 2Ethanol + 2CO_2
    
    Yeast is assumed to be produced from any remaining glucose:
        Glucose -> Yeast
    
    A compound with name 'Yeast' must be present. Note that only glucose is 
    taken into account for conversion. Cleaning and unloading time,
    `tau_0`, fraction of working volume, `V_wf`, and number of reactors,
    `N_reactors`, are attributes that can be changed. Cost of a reactor
    is based on the NREL batch fermentation tank cost assuming volumetric
    scaling with a 6/10th exponent [1]_. 
    
    References
    ----------
    .. [1] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden
        National. Renewable Energy Laboratory Golden, Colorado. P. Schoen,
        J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group
        Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics
        for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid
        Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical
        Report NREL/TP-5100-47764
    
    
    """
    
    def __init__(self, ID, 
                 ins=None, outs=(), thermo=None, 
                 tau=60.0, 
                 kinetic_reaction_system=None, 
                 map_chemicals_nsk_to_bst={}, 
                 n_simulation_steps=1000,
                 f_reset_kinetic_reaction_system=None,
                 N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36,
                 # glucose_depletion=0.99,
                 yield_ethanol=0.80, yield_isobutanol=0.15, V_wf=0.83):
        
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
            tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax
        )
        self.fermentation_reactions = tmo.ParallelReaction([tmo.Rxn('Glucose -> 2Ethanol + 2CO2',  'Glucose', yield_ethanol),
                                             tmo.Rxn('Glucose -> Isobutanol + 2CO2',  'Glucose', yield_isobutanol),],
                                             )
        self.growth = tmo.Rxn('Glucose -> Yeast',  'Glucose', 1.0)
        self.V_wf = V_wf
        
        # self.glucose_depletion = glucose_depletion
        
        self.kinetic_reaction_system = kinetic_reaction_system
        if isinstance(kinetic_reaction_system, TelluriumReactionSystem):
            self.simulate_kinetics = self._nsk_te_simulate_kinetics
        elif isinstance(kinetic_reaction_system, ReactionSystem):
            self.simulate_kinetics = self._nsk_simulate
        
        self.map_chemicals_nsk_to_bst = map_chemicals_nsk_to_bst
        self.n_simulation_steps = n_simulation_steps
        self.f_reset_kinetic_reaction_system = f_reset_kinetic_reaction_system if f_reset_kinetic_reaction_system is not None else lambda model: model.reset()
        
    def _nsk_te_simulate_kinetics(self, feed, tau, feed_spike_condition=None, plot=False):
        self.tau = tau
        map_chemicals_nsk_to_bst = self.map_chemicals_nsk_to_bst
        kinetic_reaction_system = self.kinetic_reaction_system
        self.f_reset_kinetic_reaction_system(kinetic_reaction_system)
        te_r = kinetic_reaction_system._te
        
        # get unit conversion factors and unit-based material indexers
        
        time_units = kinetic_reaction_system._units['time']
        if time_units.lower() in ('min', 'm'):
            time_conv_factor = 60.0
        elif time_units.lower() in ('sec', 's'):
            time_conv_factor = 3600.0
        elif time_units.lower() in ('hr', 'h'):
            time_conv_factor = 1.0
            
        conc_units = kinetic_reaction_system._units['conc']
        if conc_units in ('M', 'mol/L', 'kg/m3', 'kg/m^3'):
            material_indexer = 'imol'
            volume_indexer = 'F_vol'
        elif conc_units in ('g/L', 'kg/m3', 'kg/m^3'):
            material_indexer = 'imass'
            volume_indexer = 'F_vol'
            
        self._nsk_initial_concentration = initial_concentrations = {}
        for c_nsk, c_bst in map_chemicals_nsk_to_bst.items():
            exec(f'te_r.{c_nsk} = feed.{material_indexer}[c_bst]/feed.{volume_indexer}')
            exec(f'initial_concentrations[c_nsk] = te_r.{c_nsk}')
        
        te_r.simulate(0, tau*time_conv_factor, self.n_simulation_steps)
        
        if plot: te_r.plot()
        
        effluent = feed.copy()
        initially_zero = []
        initially_nonzero = []
        # breakpoint()
        for c_nsk, c_bst in map_chemicals_nsk_to_bst.items():
            if initial_concentrations[c_nsk] > 0.0:
                exec(f'effluent.{material_indexer}[c_bst] *= te_r.{c_nsk}/initial_concentrations[c_nsk]')
                initially_nonzero.append((c_nsk, c_bst))
            else:
                initially_zero.append((c_nsk, c_bst))
        
        for c_nsk, c_bst in initially_zero:
            exec(f'material_factor = feed.{material_indexer}[initially_nonzero[0][1]]/initial_concentrations[initially_nonzero[0][0]]; effluent.{material_indexer}[c_bst] = material_factor * te_r.{c_nsk}')
                
        return effluent
    
    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins)
        self.mixed_feed = mixed_feed = effluent.copy()
        kinetic_reaction_system = self.kinetic_reaction_system
        if kinetic_reaction_system is not None:
            effluent.copy_like(self.simulate_kinetics(feed=effluent, tau=self._tau))
            effluent.imol['CO2'] = effluent.imol['Ethanol']
            # breakpoint()
            # !!! make ammonia yeast-dependent
        else:
            self.fermentation_reactions(effluent)
            self.growth(effluent)
            effluent.imass['Yeast'] += effluent.imass['NH3']
            effluent.imol['NH3'] = 0.
        effluent.empty_negative_flows()
        vent.empty()
        # vent.receive_vent(effluent)
        vent.receive_vent(effluent, energy_balance=False)
        
    @property
    def Hnet(self):
        # X = self.fermentation_reactions.X
        fermentation_reactions = self.fermentation_reactions
        glucose = self.ins[0].imol['Glucose']
        # return self.reaction.dH * glucose + self.growth.dH * (1 - X) + self.H_out - self.H_in
        # return np.sum([i.dH for i in fermentation_reactions]) * glucose + np.sum([self.growth.dH * (1 - i.X) + self.H_out - self.H_in for i in fermentation_reactions])
        return np.sum(fermentation_reactions[0].dH) * glucose + np.sum([self.growth.dH * (1 - i.X) + self.H_out - self.H_in for i in fermentation_reactions]) #!!! reactions[0]

SSFEtOHIBO = SimultaneousSaccharificationFermentationEthanolIsobutanol

#%%
@cost('Flow rate', 'Column', kW=151, BM=1.8,
      cost=2601000, CE=521.9, S=22687, n=0.6)
class MolecularSieve(Splitter):
    """
    Create an isobutanol/water molecular sieve for bioisobutanol plants.
    The molecular sieve is modeled as a component wise separator. Costing
    is based on scaling by the 6/10ths rule from an NREL TEA report [1]_.
    
    Parameters
    ----------
    ins : 
        * [0] Feed (gas)
    outs : 
        * [0] Split stream (gas)
        * [1] Remainder stream (gas)
    split : array_like
            Componentwise split to the 0th output stream

    References
    ----------
    .. [1] Process Design and Economics for Biochemical Conversion of
        Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and
        Enzymatic Hydrolysis of Corn Stover. D. Humbird, R. Davis, L.
        Tao, C. Kinchin, D. Hsu, and A. Aden (National Renewable Energy
        Laboratory Golden, Colorado). P. Schoen, J. Lukas, B. Olthof,
        M. Worley, D. Sexton, and D. Dudgeon (Harris Group Inc. Seattle,
        Washington and Atlanta, Georgia)
    
    """
    _units = {'Flow rate': 'kg/hr'}
    def _init(self, split, order=None, P=None, approx_duty=True):
        Splitter._init(self, order=order, split=split)
        self.P = None
        self.approx_duty = approx_duty
        
    def _run(self):
        Splitter._run(self)
        P = self.P
        if P is None: P = self.ins[0].P
        for i in self.outs: i.P = P

    def _design(self):
        self.design_results['Flow rate'] = flow = self._outs[1].F_mass
        if self.approx_duty:
            T = self.outs[0].T
            self.add_heat_utility(1429.65 * flow, T)
            self.add_heat_utility(-55.51 * flow, T)