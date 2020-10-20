# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx

__all__ = (
    'HammerMill', 
    'CornStorage',
    'CleaningSystem',
    'Liquefaction',
    'JetCooker',
    'SimultaneousSaccharificationFermentation', 'SSF',
    'DDGS_Dryer',
)


CAS_water = '7732-18-5'

class HammerMill(bst.Unit):
    pass


class CornStorage(bst.Unit):
    pass


class CleaningSystem(bst.Unit):
    pass


class JetCooker(bst.Unit):
    """
    ins : stream sequence
    
        [0] Feed
        
        [1] Steam
    
    outs : stream
        Mixed product.
    
    """
    _N_outs = 1
    _N_ins = 2
    _N_heat_utilities = 1
    
    def __init__(self, ID="", ins=None, outs=(), thermo=None, T=483.15):
        super().__init__(ID, ins, outs, thermo)
        self.T = T
    
    @staticmethod
    def _T_objective_function(steam_mol, T, steam, effluent, feed):
        steam.imol[CAS_water] = steam_mol
        effluent.mol[:] = steam.mol + feed.mol
        effluent.H = feed.H + steam.H
        return effluent.T - T
    
    def _run(self):
        feed, steam = self._ins
        steam_mol = steam.F_mol
        effluent, = self.outs
        steam_mol = flx.aitken_secant(self._T_objective_function,
                                      steam_mol, steam_mol+1., 
                                      1e-4, 1e-4,
                                      args=(self.T, steam, effluent, feed),
                                      checkroot=False)
        effluent.P = steam.P / 2.
        hu, = self.heat_utilities
        hu(steam.H, effluent.T)


class DDGS_Dryer(bst.Unit):
    
    def __init__(self, ID="", ins=None, outs=(), thermo=None,
                 moisture_content=0.5):
        super().__init__(ID, ins, outs, thermo)
        self.moisture_content = moisture_content
        
    def _run(self):
        wet_solids, air = self.ins
        vent, dry_solids = self.outs
        vent.copy_like(air)
        dry_solids.copy_like(wet_solids)
        tmo.separations.adjust_moisture_content(dry_solids, vent, self.moisture_content)


class Liquefaction(bst.Unit):
    """
    Create a Liquefaction unit operation that models the converion
    for Starch to Glucose oligomers.
    
    Parameters
    ----------
    ins : stream
        Inlet fluid.
    outs : stream
        Outlet fluid.
    yield_: float
        Yield of starch to glucose as a fraction of the theoretical yield.
    
    Notes
    -----
    The conversion of Starch to Glucose oligomers is modeled according to the
    following stoichiometry:
        
    Starch + H2O -> Glucose
    
    Where starch is a chemical with formula C6H10O5 that represents linked 
    glucose monomers (dehydrated from linkage).
    
    The dextrose equivalent, and for that manner the degree of polymerization,
    is not taken into account in this unit. However, the conversion is equivalent
    to the conversion of starch to fermentable saccharides, which is what matters
    downstream.
    
    References
    ----------
    TODO
    
    """
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, yield_=0.95):
        super().__init__(ID, ins, outs, thermo)
        self.reaction = tmo.reaction.Reaction('Starch + H2O -> Glucose', 'Starch', yield_)
        
    @property
    def yield_(self):
        return self.reaction.X
    @yield_.setter
    def yield_(self, X):
        self.reaction.X = X

    def _run(self):
        effluent, = self.outs
        effluent.mix_from(self.ins)
        self.reaction(effluent)
        

class SimultaneousSaccharificationFermentation(bst.BatchBioreactor):
    """
    Create a SimultaneousSaccharificationFermentation unit operation that 
    models the simultaneous saccharification and fermentation in the conventional
    dry-grind enthanol process.
    
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
    
    Yeast is assumed to be produced from the any remaining glucose:
        Glucoes -> Yeast
    
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
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 tau,  N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36,
                 yield_=0.9):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
            tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax
        )
        self.reaction = tmo.reaction.Reaction('Glucose -> 2Ethanol + 2CO2',  'Glucose', yield_)
    
    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins)
        self.reaction(effluent)
        mass_glucose = effluent.imass['Glucose']
        effluent.imass['Yeast'] += mass_glucose
        mass_glucose.value = 0.
        vent.receive_vent(effluent)
    
SSF = SimultaneousSaccharificationFermentation