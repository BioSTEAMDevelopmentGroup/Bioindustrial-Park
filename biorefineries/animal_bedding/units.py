# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 00:55:03 2019

Equipment from the Humbird 2011 Report.
    
Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A., Dudgeon, D. (2011). Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269

@author: yoelr
"""
import os
import sys
import flexsolve as flx
from thermosteam import MultiStream
from biosteam import Unit
from biosteam.units.decorators import cost, design
from biosteam.units.design_tools import size_batch
from biosteam.units.design_tools.specification_factors import  material_densities_lb_per_in3
from biosteam.units.design_tools import column_design
import thermosteam as tmo
import biosteam as bst

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# %% Add excel unit operations

from biosteam.units.factories import xl2mod
path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '_humbird2011.xlsx')
xl2mod(path, sys.modules[__name__])
del sys, xl2mod, os, path

# %% Constants

_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
# _m3hr2gpm = 4.40287
_hp2kW = 0.7457
_Gcal2kJ = 4184e3




# %% Pre-washing

@cost('Flow rate', 'Drum filter',
      cost=30000, CE=525.4, S=10, n=1, kW=1.63, BM=1, N='N_filters')
class DrumFilter(bst.Splitter):
    _units = {'Flow rate': 'm3/hr', 'Cost': 'USD'}   
    _N_ins = 1
    _N_outs= 2
    #: Number of parallel filters
    N_filters = 190 #total
    N_spare = 4
    def __init__(self, ID='', ins=None, outs=(), *, order=None, flux=1220.6, split, moisture_content):
        bst.Splitter.__init__(self, ID, ins, outs, order=order, split=split)
        #: Moisture content of retentate
        self.moisture_content = moisture_content
        self.flux = flux #kg/hr/m2
        assert self.isplit['Water'] == 0, 'cannot define water split, only moisture content'
    
    def run_split_with_solids(self,mc):
        """Splitter mass and energy balance function with mixing all input streams."""
        top, bot = self.outs
        feed = self.ins[0]
        top.copy_like(feed)
        bot.copy_like(top)
        top_mass = top.mass
        F_mass_solids = sum(top_mass*self.split)
        TS=1-mc
        F_mass_tot = F_mass_solids/TS
        F_mass_wat = F_mass_tot - F_mass_solids
        top_mass[:] *= self.split
        top.imass['Water']=F_mass_wat
        bot.mass[:] -= top_mass
        
    def _run(self):
        
        self.run_split_with_solids(self.moisture_content)
        retentate, permeate = self.outs
        if permeate.imass['Water'] < 0:
            import warnings
            warnings.warn(f'not enough water for {repr(self)}')
            
    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_vol/(self.N_filters-self.N_spare)


@cost('Tank volume', 'Tanks', cost=3840e3/8, S=250e3*_gal2m3, 
      CE=522, n=0.7, BM=2.0, N='N_tanks')
@cost('Tank volume', 'Agitators', CE=522, cost=52500,
      S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_tanks')
@cost('Flow rate', 'Transfer pumps', kW=58, S=352*_gpm2m3hr,
      cost=47200/5, CE=522, n=0.8, BM=2.3, N='N_tanks')
class WashingTank(bst.Unit):
    """Washing tank system where part of the manure organics and inorganics 
        dissolve   
    **Parameters**   
        **reactions:** [ReactionSet] Washing reactions.
    **ins**    
        [0] Feed        
        [1] Process water        
    **outs**   
        [0] Washed feed

    """
#    purchase_cost = installation_cost = 0
    _N_ins = 2
    _N_outs = 1
    N_tanks = 1
    tau_tank = 5/60/N_tanks #Residence time in each tank
    V_wf = 0.9
    _units = {'Flow rate': 'm3/hr',
              'Tank volume': 'm3'}
    def __init__(self, ID='', ins=None, outs=(), *, reactions):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions  
        
    def _run(self):
        feed, process_water = self.ins
        washed, = self.outs   
        washed.mix_from([feed,process_water])
        self.reactions(washed.mass)

    def _design(self):
        effluent = self.outs[0]
        v_0 = effluent.F_vol
        Design = self.design_results
        Design['Tank volume'] = v_0*self.tau_tank/self.V_wf
        Design['Flow rate'] = v_0


# %% Pretreatment
class SteamMixer(Unit):
    """
    **ins**
    
        [0] Feed
        
        [1] Steam
    
    **outs**
    
        [0] Mixed
    
    """
    _N_outs = 1
    _N_ins = 2
    _N_heat_utilities = 1
    def __init__(self, ID='', ins=None, outs=(), *, P):
        super().__init__(ID, ins, outs)
        self.P = P
    
    @staticmethod
    def _P_at_flow(mol_water, P, steam, mixed, feed):
        steam.imol['7732-18-5'] = mol_water
        mixed.mol[:] = steam.mol + feed.mol
        mixed.H = feed.H + mol_water * 40798
        P_new = mixed.chemicals.Water.Psat(mixed.T)
        return P - P_new
    
    def _run(self):
        feed, steam = self._ins
        steam_mol = steam.F_mol
        mixed = self.outs[0]
        steam_mol = flx.aitken_secant(self._P_at_flow,
                                  steam_mol, steam_mol+0.1, 
                                  1e-4, 1e-4,
                                  args=(self.P, steam, mixed, feed))
        mixed.P = self.P      
        hu = self.heat_utilities[0]
        hu(steam.Hvap, mixed.T)

    @property
    def installation_cost(self): return 0
    @property
    def purchase_cost(self): return 0
    def _design(self): pass
    def _cost(self): pass
        

@cost('Dry flow rate', units='kg/hr', S=83333, CE=522,
      cost=19812400, n=0.6, kW=4578, BM=1.5)
class PretreatmentReactorSystem(Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = bst.Flash._graphics
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self._multistream = MultiStream(None)
        self.reactions = ParallelRxn([
    #            Reaction definition                 Reactant    Conversion
    Rxn('Glucan + H2O -> Glucose',                   'Glucan',   0.0510),#
    Rxn('Glucan + H2O -> GlucoseOligomer',           'Glucan',   0.0750),#
    Rxn('Glucan -> HMF + 2 H2O',                     'Glucan',   0.0043),#
    Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.5190),#
    Rxn('Xylan + H2O -> XyloseOligomer',             'Xylan',    0.2610),#
    Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.0860),#
    Rxn('Arabinan + H2O -> Arabinose',               'Arabinan', 1.0000),#
    Rxn('Arabinan + H2O -> ArabinoseOligomer',       'Arabinan', 0.0000),#
    Rxn('Arabinan -> Furfural + 2 H2O',              'Arabinan', 0.0000),#
    Rxn('Acetate -> AceticAcid',                     'Acetate',  1.0000),
    Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.0470),
    Rxn('Extract -> ExtractVol',                     'Extract',  0.7000),
    Rxn('Extract -> ExtractNonVol',                  'Extract',  0.3000)])
        vapor, liquid = self.outs
        vapor.phase = 'g'
    
    def _run(self):
        ms = self._multistream
        feed = self.ins[0]
        vapor, liquid = self.outs
        liquid.copy_like(feed)
        self.reactions(liquid) #self.reactions.adiabatic_reaction(liquid) 
        ms.copy_like(liquid)
        H = ms.H + liquid.Hf - feed.Hf
        ms.vle(T=190+273.15, H=H)
        vapor.mol[:] = ms.imol['g']
        liquid.mol[:] = ms.imol['l']
        vapor.T = liquid.T = ms.T
        vapor.P = liquid.P = ms.P



@cost('Flow rate', 'Sieve filter',
      cost=14800, CE=551, S=0.2273, n=0.64, BM=1)
class SieveFilter(bst.Splitter):
    _units = {'Flow rate': 'm3/hr'}   
    _N_ins = 1
    _N_outs= 2
    
    def __init__(self, ID='', ins=None, outs=(), *, order=None, WIS=False, split, moisture_content):
        bst.Splitter.__init__(self, ID, ins, outs, order=order, split=split)
        #: Moisture content of retentate
        #If WIS=True, the moisture content is 1-WIS content. Otherwise is 1- totals solid content
        self.WIS=WIS
        self.moisture_content = moisture_content
        assert self.isplit['Water'] == 0, 'cannot define water split, only moisture content'
    
    def run_split_with_solids(self,mc):
        """Splitter mass and energy balance function with mixing all input streams."""
        top, bot = self.outs
        feed = self.ins[0]
        top.copy_like(feed)
        bot.copy_like(top)
        top_mass = top.mass
        F_mass_ins = sum(top_mass*self.split)
        F_mass_sol = top.F_mass - F_mass_ins
        F_mass_wat = top.imass['Water']
        x_sol = mc*F_mass_ins/(F_mass_wat-mc*F_mass_sol)
        self.split[self.split==0] = x_sol
        top_mass[:] *= self.split
        bot.mass[:] -= top_mass
    
    def run_split_with_solidsWIS(self,mc):
        """Splitter mass and energy balance function with mixing all input streams."""
        top, bot = self.outs
        feed = self.ins[0]
        top.copy_like(feed)
        bot.copy_like(top)
        top_mass = top.mass
        WIS=1-mc
        F_mass_ins = sum(top_mass*self.split)
        F_mass_tot_out = F_mass_ins/WIS
        F_mass_sol_out = F_mass_tot_out - F_mass_ins
        F_mass_sol_in = feed.F_mass - F_mass_ins
    
        x_sol = F_mass_sol_out/F_mass_sol_in
        self.split[self.split==0] = x_sol
        top_mass[:] *= self.split
        bot.mass[:] -= top_mass
        
    def _run(self):
        if self.WIS:
            self.run_split_with_solidsWIS(self.moisture_content)
        else:
            self.run_split_with_solids(self.moisture_content)
        retentate, permeate = self.outs
        if permeate.imass['Water'] < 0:
            import warnings
            warnings.warn(f'not enough water for {repr(self)}')
            
    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_vol

        
   
@cost('Area', 'Pressure filter',
      cost=1220, CE=551, S=0.092903, n=0.71, BM=1.7)
class PressureFilter(bst.Splitter):
    _units = {'Flow rate': 'kg/hr','Area': 'm2','Power': 'kW'}   
    _N_ins = 1
    _N_outs= 2
    _pressure = 13e5
    _efficiency = 0.95
    
    def __init__(self, ID='', ins=None, outs=(), *, order=None, WIS=False, flux=1220.6, split, moisture_content):
        bst.Splitter.__init__(self, ID, ins, outs, order=order, split=split)
        #: Moisture content of retentate
        #If WIS=True, the moisture content is 1-WIS content. Otherwise is 1- totals solid content
        self.WIS=WIS
        self.moisture_content = moisture_content
        self.flux = flux #kg/hr/m2
        assert self.isplit['Water'] == 0, 'cannot define water split, only moisture content'
    
    def run_split_with_solids(self,mc):
        """Splitter mass and energy balance function with mixing all input streams."""
        top, bot = self.outs
        feed = self.ins[0]
        top.copy_like(feed)
        bot.copy_like(top)
        top_mass = top.mass
        F_mass_ins = sum(top_mass*self.split)
        F_mass_sol = top.F_mass - F_mass_ins
        F_mass_wat = top.imass['Water']
        x_sol = mc*F_mass_ins/(F_mass_wat-mc*F_mass_sol)
        self.split[self.split==0] = x_sol
        top_mass[:] *= self.split
        bot.mass[:] -= top_mass
    
    def run_split_with_solidsWIS(self,mc):
        """Splitter mass and energy balance function with mixing all input streams."""
        top, bot = self.outs
        feed = self.ins[0]
        top.copy_like(feed)
        bot.copy_like(top)
        top_mass = top.mass
        WIS=1-mc
        F_mass_ins = sum(top_mass*self.split)
        F_mass_tot_out = F_mass_ins/WIS
        F_mass_sol_out = F_mass_tot_out - F_mass_ins
        F_mass_sol_in = feed.F_mass - F_mass_ins
    
        x_sol = F_mass_sol_out/F_mass_sol_in
        self.split[self.split==0] = x_sol
        top_mass[:] *= self.split
        bot.mass[:] -= top_mass
        
    def _run(self):
        if self.WIS:
            self.run_split_with_solidsWIS(self.moisture_content)
        else:
            self.run_split_with_solids(self.moisture_content)
        retentate, permeate = self.outs
        if permeate.imass['Water'] < 0:
            import warnings
            warnings.warn(f'not enough water for {repr(self)}')
            
    def _design(self):
        self.design_results['Flow rate'] = abs(self.outs[1].F_mass)
        self.design_results['Area'] = abs(self.outs[1].F_mass)/self.flux #m2
        feed = self.ins[0]        
        light_ind = feed.chemicals._light_indices  
        l = [a for a in feed.vol[light_ind] if not a==0]
        v_0 = feed.F_vol - sum(l)
        self.design_results['Power'] = power_rate =  (v_0/3600)*self._pressure/self._efficiency/1000 #kW
        pu = self.power_utility
        pu(rate=power_rate)


# %% Saccharification and fermentation
        
@cost('Flow rate', 'Pumps',
      S=43149, CE=522, cost=24800, n=0.8, kW=40, BM=2.3)
@cost('Stage #1 reactor volume', 'Stage #1 reactors',
      cost=37700, S=20*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #2 reactor volume', 'Stage #2 reactors',
      cost=58300, S=200*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #3 reactor volume', 'Stage #3 reactors',
      cost=78800, S=2e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 reactors',
      cost=176e3, S=20e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 agitators',
      cost=26e3/2, S=20e3*_gal2m3, kW=7.5, CE=522, n=0.5, BM=1.5)
@cost('Stage #5 reactor volume', 'Stage #5 reactors',
      cost=590e3, S=200e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #5 reactor volume', 'Stage #5 agitators',
      cost=43e3/2, S=200e3*_gal2m3, kW=10, CE=522, n=0.5, BM=1.5)
class SeedTrain(Unit):
    _N_ins = 1
    _N_outs= 2
    _N_heat_utilities = 1
    
    _units= {'Flow rate': 'kg/hr',
             'Stage #1 reactor volume': 'm3',
             'Stage #2 reactor volume': 'm3',
             'Stage #3 reactor volume': 'm3',
             'Stage #4 reactor volume': 'm3',
             'Stage #5 reactor volume': 'm3'}
    
    @property
    def N_stages(self): 
        """Number of stages in series."""
        return 5
    
    #: Number of parallel seed trains
    N_trains = 2
    
    #: Cycle time for each batch (hr)
    tau_batch = 24
    
    @property
    def tau_turnover(self):
        """Turnover time (hr) calculated by batch time divided by number of trains."""
        return self.tau_batch/self.N_trains
    
    #: Operating temperature (K)
    T = 32+273.15
    
    # #: wt % media (e.g. corn steep liquor) in each stage 
    # media_loading = 0.50
    
    # #: Diammonium phosphate loading in g/L of fermentation broth
    # DAP = 0.67 
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = ParallelRxn([
    #   Reaction definition                                                         Reactant    Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',                                             'Glucose',   0.0400),##0.04
    Rxn('Glucose + 0.62 NH3 + 2.17 O2 -> 3.65 S_cerevisiae + 3.71 H2O + 2.34 CO2',  'Glucose',   0.2400),##0.90
    Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',                                       'Glucose',   0.0040),
    Rxn('2 Glucose + 1.5 O2 -> 3 SuccinicAcid + 3 H2O',                             'Glucose',   0.0060)
    ])
    #20.65 g/mol S_cerivesiae
    def _run(self):
        feed, = self.ins
        vent, effluent= self.outs
        effluent.copy_flow(feed)
        self.reactions.force_reaction(effluent.mol) # TODO: Ignore negative O2; probably bug in _system.py
        effluent.T = self.T
        vent.phase = 'g'
        vent.copy_flow(effluent, ('CO2', 'NH3', 'O2', 'N2'), remove=True)

    def _design(self): 
        maxvol = self.outs[1].F_vol*self.tau_turnover
        vol = maxvol*10**-self.N_stages
        Design = self.design_results
        for i in range(1, self.N_stages+1):
            Design[f'Stage #{i} reactor volume'] = vol
            vol *= 10 
        Design['Flow rate'] = sum([i.F_mass for i in self.outs])
        self.heat_utilities[0](self.Hnet, self.T)

    def _cost(self):
        N = self.N_trains
        D = self.design_results
        C = self.purchase_costs
        kW = 0
        for i, x in self.cost_items.items():
            S = D[x._basis]
            q = S/x.S
            C[i] = N*bst.CE/x.CE*x.cost*q**x.n
            kW += N*x.kW*q
        self.power_utility(kW)
        

 

@cost('Flow rate', 'Recirculation pumps', kW=30, S=340*_gpm2m3hr,
      cost=47200, n=0.8, BM=2.3, CE=522, N='N_recirculation_pumps')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900,
      S=5*_Gcal2kJ, n=0.7, BM=2.2, N='N_reactors') # Based on a similar heat exchanger
@cost('Reactor volume', 'Agitators', CE=522, cost=52500,
      S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_reactors')
@cost('Reactor volume', 'Reactors', CE=522, cost=844000,
      S=1e6*_gal2m3, n=0.5, BM=1.5, N='N_reactors')
@cost('Flow rate', 'Transfer pumps', kW=58, S=352*_gpm2m3hr,
      cost=47200/5, CE=522, n=0.8, BM=2.3, N='N_transfer_pumps')
@cost('Tank volume', 'Tanks', cost=3840e3/8, S=250e3*_gal2m3, 
      CE=522, n=0.7, BM=2.0, N='N_tanks')
class SaccharificationAndCoFermentation(Unit):
    _N_ins = 3
    _N_outs = 3
    _N_heat_utilities = 2
    
    #: Saccharification temperature (K)
    T_saccharification = 48+273.15
    
    #: Fermentation temperature (K)
    T_fermentation = 32+273.15
    
    #: Residence time of countinuous saccharification tanks (hr)
    tau_tank = 24
    
    #: Saccharification time (hr)
    tau_saccharification = 60
    
    #: Co-Fermentation time (hr)
    tau_cofermentation = 36
    
    #: Unload and clean up time (hr)
    tau_0 = 4
    
    #: Working volume fraction (filled tank to total tank volume)
    V_wf = 0.9
    
    #: Number of reactors
    N_reactors = 12
    
    #: Number of continuous saccharification tanks
    N_tanks = 8
    
    #: Number of transfer pumps
    N_transfer_pumps = 5
    
    #: Number of recirculation pumps
    N_recirculation_pumps = 5
    
    _units = {'Flow rate': 'm3/hr',
              'Tank volume': 'm3',
              'Reactor volume': 'm3',
              'Reactor duty': 'kJ/hr'}
    

    
    def __init__(self, ID='', ins=None, outs=(), saccharified_slurry_split = 0.1, P=101325):
        Unit.__init__(self, ID, ins, outs)
        self.P = P
            # Split to outs[2]
        self.saccharified_slurry_split = saccharified_slurry_split
        #: [ParallelReaction] Enzymatic hydrolysis reactions including from downstream batch tank in co-fermentation.
        self.saccharification = ParallelRxn([
    #   Reaction definition                   Reactant    Conversion
    Rxn('Glucan -> GlucoseOligomer',          'Glucan',             0.0400),
    Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',             0.0120),
    Rxn('Glucan + H2O -> Glucose',            'Glucan',             0.7800),#changed
    Rxn('GlucoseOligomer + H2O -> Glucose',   'GlucoseOligomer',    1.0000),
    Rxn('Cellobiose + H2O -> Glucose',        'Cellobiose',         1.0000),
    Rxn('Xylan -> XyloseOligomer',            'Xylan',              0.0400),
    Rxn('Xylan + H2O -> Xylose',              'Xylan',              0.9000)
    ])
    
        self.loss = ParallelRxn([
    #   Reaction definition               Reactant    Conversion
    Rxn('Glucose -> 2 LacticAcid',       'Glucose',   0.0300),
    Rxn('3 Xylose -> 5 LacticAcid',      'Xylose',    0.0300),
    Rxn('3 Arabinose -> 5 LacticAcid',   'Arabinose', 0.0300),
    Rxn('Galactose -> 2 LacticAcid',     'Galactose', 0.0300),
    Rxn('Mannose -> 2 LacticAcid',       'Mannose',   0.0300),])
    
        self.cofermentation = ParallelRxn([
    #   Reaction definition                                                         Reactant    Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',                                             'Glucose',   0.8400),#changed
    Rxn('Glucose + 0.62 NH3 + 2.17 O2 -> 3.65 S_cerevisiae + 3.71 H2O + 2.34 CO2',  'Glucose',   0.0130),
    Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',                                       'Glucose',   0.0040),
    Rxn('2 Glucose + 1.5 O2 -> 3 SuccinicAcid + 3 H2O',                             'Glucose',   0.0060),    
    ])
        
        self.saccharified_stream = tmo.Stream(None)
    
    def _run(self):
        feed, DAP, air = self.ins
        vent, effluent, sidedraw = self.outs
        vent.P = effluent.P = sidedraw.P = self.P
        ss = self.saccharified_stream
        ss.T = sidedraw.T = self.T_saccharification
        vent.T = effluent.T = self.T_fermentation
        vent.phase = 'g'
        ss.copy_flow(feed)
        self.saccharification(ss.mol)
        sidedraw.mol[:] = ss.mol * self.saccharified_slurry_split
        effluent.mol[:] = ss.mol - sidedraw.mol + DAP.mol + air.mol
        self.loss(effluent.mol)
        self.cofermentation.force_reaction(effluent.mol)
        vent.receive_vent(effluent)
    
    def _design(self):
        effluent = self.outs[1]
        v_0 = effluent.F_vol
        Design = self.design_results
        Design['Tank volume'] = v_0*self.tau_tank/self.V_wf/self.N_tanks
        Design['Flow rate'] = v_0/self.N_transfer_pumps
        tau = self.tau_saccharification + self.tau_cofermentation
        Design.update(size_batch(v_0, tau, self.tau_0, self.N_reactors, self.V_wf))
        hu_cooling, hu_fermentation = self.heat_utilities
        mixture = self.thermo.mixture
        ss = self.saccharified_stream
        mol = ss.mol
        hu_cooling(mixture.H('l', mol, self.T_fermentation, 101325.)
                   - mixture.H('l', mol, self.T_saccharification, 101325.),
                   self.T_fermentation)
        ei = effluent.chemicals.index('Ethanol')
        ethanol = (sum([i.mol[ei] for i in self.outs])
                   - sum([i.mol[ei] for i in self.ins]))
        duty = ethanol*-5568
        hu_fermentation(duty, effluent.T)
        Design['Reactor duty'] = -duty



# %% Ethanol purification
class DistillationColumn(bst.BinaryDistillation):

    def __init__(self, ID='', ins=None, outs=(),*,P=101325, energy_integration=False, LHK, k, y_top,x_bot):
        
        bst.BinaryDistillation.__init__(self, ID=ID, ins=ins,
                                        P=P, y_top=y_top, x_bot=x_bot,
                                        k=k, LHK=LHK)

        self.energy_integration=energy_integration
        
    def _run(self):
        bst.BinaryDistillation._run(self)
            
    def _design(self):
        bst.BinaryDistillation._design(self)       
        H=self.get_design_result('Height','ft')
        Di=self.get_design_result('Diameter','ft')
        Po = self.P * 0.000145078 #to psi
        Pgauge = Po - 14.69
        if Pgauge<0.0: Po=-Pgauge+14.69
        Design = self.design_results
        Design['Wall thickness'] = tv = column_design.compute_tower_wall_thickness(Po, Di, H)
        rho_M = material_densities_lb_per_in3[self.vessel_material]
        Design['Weight'] = column_design.compute_tower_weight(Di, H, tv, rho_M)
        W = Design['Weight'] # in lb
        L = Design['Height']*3.28 # in ft
        Cost = self.purchase_costs
        F_VM = self._F_VM
        Cost['Tower'] = column_design.compute_purchase_cost_of_tower(Di, L, W, F_VM)
        
        if self.energy_integration:  
            self.boiler.heat_utilities[0].flow=0
            self.boiler.heat_utilities[0].cost=0
        self._simulate_components()    
# %% Biogas production

@cost('Reactor cooling', 'Heat exchangers', CE=522, cost=23900,
      S=5*_Gcal2kJ, n=0.7, BM=2.2, N='N_reactors') # Based on a similar heat exchanger
@cost('Reactor volume', 'Agitators', CE=522, cost=52500,
      S=1e6*_gal2m3, n=0.5, kW=30, BM=1.5, N='N_reactors')
@cost('Reactor volume', 'Reactors', CE=522, cost=12.8e6/4,
      S=5450, n=0.5, BM=1, N='N_reactors')
@cost('Flow rate', 'Transfer pumps', kW=4.4, S=352*_gpm2m3hr,
      cost=47200/5, CE=522, n=0.8, BM=2.3, N='N_transfer_pumps')
@cost('Raw biogas', 'Biogas purification', kW=0.25, S=0.0446,
      cost=1430, CE=522, n=1, BM=1)
class AnaerobicDigestion(bst.Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **reactions:** [ReactionSet] Anaerobic digestion reactions.
        
        **sludge_split:** [Array] Split between waste water and sludge
        
    **ins**
    
        [0] Waste water
        
        [1] Cool well water
        
    **outs**
    
        [0] Biogas
        
        [1] Waste water
        
        [2] Sludge
        
        [3] Hot well water
    
    """
    _N_ins = 1
    _N_outs = 3
    _N_heat_utilities = 1
    #: Residence time of countinuous anaerobic digesters (hr)
    tau = 30*24
      
    #: Working volume fraction (filled tank to total tank volume)
    V_wf = 0.8
    
    #: Number of reactors
    N_reactors = 4
    
    #: Number of transfer pumps
    N_transfer_pumps = 4
    
    #: Number of recirculation pumps
    N_recirculation_pumps = 4
    
    _units = {'Flow rate': 'm3/hr',
              'Reactor volume': 'm3',
              'Reactor cooling': 'kJ/hr',
              'Raw biogas':  'kmol/hr',
              'Pure biogas':  'kmol/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), *, reactions, sludge_split,T):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions
        self.sludge_split = sludge_split
        self.multi_stream = tmo.MultiStream()
        self.T=T       
        
    def _run(self):
        feed, = self.ins
        biogas, waste, sludge = self.outs
        biogas.phase = 'g'        
        biogas.T = waste.T = sludge.T = self.T  
        sludge.copy_flow(feed)
        self.reactions(sludge.mol)
        self.multi_stream.copy_flow(sludge)
        self.multi_stream.vle(P=101325, H=self.multi_stream.H)
        biogas.mol[:] = self.multi_stream.imol['g']
        liquid_mol = self.multi_stream.imol['l']
        sludge.mol[:] = liquid_mol * self.sludge_split
        waste.mol[:] = liquid_mol - sludge.mol
        biogas.receive_vent(waste, accumulate=True)

    def _design(self):
        feed = self.ins[0]
        biogas = self.outs[0]
        light_ind = feed.chemicals._light_indices  
        l = [a for a in feed.vol[light_ind] if not a==0]
        v_0 = feed.F_vol - sum(l)
        Design = self.design_results
        Design['Reactor volume'] = v_0*self.tau/self.V_wf/self.N_reactors
        Design['Flow rate'] = v_0/self.N_transfer_pumps
        Design['Raw biogas'] = biogas.F_mol - biogas.imol['Water']
        Design['Pure biogas'] = biogas.imol['CH4']
        hu_cooling, = self.heat_utilities
        H_at_35C = feed.thermo.mixture.H(mol=feed.mol, phase='l', T=self.T, P=101325)
        duty = H_at_35C - feed.H
        hu_cooling(duty,self.T)
        Design['Reactor cooling'] = abs(duty)
# %% Waste water treatment

@cost('Flow rate', 'Waste water system', units='kg/hr', CE=551,
      cost=50280080., n=0.6, BM=1, kW=7139/1.05, S=393100)
class WasteWaterSystemCost(bst.Unit): pass


class AnaerobicDigestionWWT(bst.Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **reactions:** [ReactionSet] Anaerobic digestion reactions.
        
        **sludge_split:** [Array] Split between waste water and sludge
        
    **ins**
    
        [0] Waste water
        
        [1] Cool well water
        
    **outs**
    
        [0] Biogas
        
        [1] Waste water
        
        [2] Sludge
        
        [3] Hot well water
    
    """
    purchase_cost = installation_cost = 0
    _N_ins = 2
    _N_outs = 4
    def __init__(self, ID='', ins=None, outs=(), *, reactions, sludge_split,T):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions
        self.sludge_split = sludge_split
        self.multi_stream = tmo.MultiStream()
        self.T=T       
        
    def _run(self):
        feed, cool_water = self.ins
        biogas, waste, sludge, hot_water = self.outs
        hot_water.link_with(cool_water, TP=False)
        hot_water.T = feed.T - 5
        H_at_35C = feed.thermo.mixture.H(mol=feed.mol, phase='l', T=self.T, P=101325)
        cool_water.mol[:] *= (feed.H - H_at_35C)/(hot_water.H - cool_water.H)
        biogas.phase = 'g'        
        biogas.T = waste.T = sludge.T = self.T  
        sludge.copy_flow(feed)
        self.reactions(sludge.mol)
        self.multi_stream.copy_flow(sludge)
        self.multi_stream.vle(P=101325, H=self.multi_stream.H)
        biogas.mol[:] = self.multi_stream.imol['g']
        liquid_mol = self.multi_stream.imol['l']
        sludge.mol[:] = liquid_mol * self.sludge_split
        waste.mol[:] = liquid_mol - sludge.mol
        biogas.receive_vent(waste, accumulate=True)
        
    
class AerobicDigestionWWT(bst.Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **reactions:** [ReactionSet] Anaerobic digestion reactions.
        
        **sludge_split:** [Array] Split between waste water and sludge
        
    **ins**
    
        [0] Waste water
        
        [1] Air
        
        [2] Caustic
        
    **outs**
    
        [0] Vent
        
        [1] Treated waste water

    """    
    _N_ins = 3
    _N_outs = 2
    purchase_cost = installation_cost = 0
    evaporation = 4/355
    
    def __init__(self, ID='', ins=None, outs=(), *, reactions):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions
    
    def _run(self):
        waste, air, caustic = self._ins
        vent, water = self.outs
        vent.phase = 'g'
        water.copy_like(waste)
        water.mol[:] += air.mol
        water.mol[:] += caustic.mol
        self.reactions(water.mol)
        vent.copy_flow(water, ('CO2', 'O2', 'N2'))
        vent.imol['7732-18-5'] = water.imol['7732-18-5'] * self.evaporation
        water.mol[:] -= vent.mol

@cost('Flow rate', units='kg/hr',
      S=63, cost=421e3, CE=522, BM=1.8, n=0.6)
class CIPpackage(bst.Facility):
    line = 'CIP Package'
    network_priority = 0
    _N_ins = 1
    _N_outs = 1
        
