# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""TankPurchaseCostAlgorithm
"""
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
from biosteam.units.decorators import cost, copy_algorithm
from biosteam.units.design_tools import CEPCI_by_year, cylinder_diameter_from_volume, cylinder_area
from biosteam import tank_factory
import numpy as np

__all__ = (
    'GrainHandling',
    'HammerMill', 
    'CornStorage',
    'CleaningSystem',
    'MilledCornSurgeTank',
    'MilledCornHopper',
    'MilledCornWeighTank',
    'LimeHopper',
    'AmmoniaTank',
    'AlphaAmylaseTank',
    'SlurryMixTank',
    'Liquefaction',
    'JetCooker',
    'CookedSlurrySurgeTank',
    'GlucoAmylaseTank',
    'SulfuricAcidTank',
    'Saccharification',
    'YeastTank',
    'WetDDGSConveyor',
    'DDGSDryer',
    'Liquefaction',
    'SimultaneousSaccharificationFermentation', 'SSF',
    'ThermalOxidizer',
    'DDGSHandling',
    'DDGSCentrifuge',
    'PlantAir_CIP_WasteWater_Facilities',
)

CE2007 = CEPCI_by_year[2007]

CAS_water = '7732-18-5'
# Material factors. Applicable to some units only.
MF38 = 0.38 
MF90 = 0.90
MF52 = 0.52

@cost('Flow rate', units='kg/hr', CE=CE2007, cost=MF38 * 120800, S=46350.72, kW=97, n=0.6)
class GrainHandling(bst.Unit): pass

CornStorage = tank_factory('CornStorage', 
    CE=CE2007, cost=MF38 * 979300., S=185400, tau=259.2, n=1.0, V_wf=0.9, V_max=3e5, 
    V_units='m3'
)

@cost('Flow rate', units='kg/hr', CE=CE2007, cost=60300., S=45350., n=0.6, ub=7.2e5)
class CleaningSystem(bst.Splitter): pass
    

@cost('Flow rate', units='kg/hr', CE=CE2007, cost=MF38 * 98200., S=46211.33, n=0.6, ub=720000., kW=314.237)
class HammerMill(bst.Unit): pass


MilledCornSurgeTank = tank_factory('MilledCornSurgeTank', 
    CE=CE2007, cost=MF38 * 32500., S=76.90, tau=2.0, n=0.6, V_wf=0.9, V_max=200., V_units='m3'
)

MilledCornHopper = tank_factory('MilledCornHopper', 
    CE=CE2007, cost=50700., S=100.93, tau=2.0, n=0.6, V_wf=0.9, V_max=1e3, V_units='m3'
)

MilledCornWeighTank = tank_factory('MilledCornWeighTank',
    CE=CE2007, cost=MF38 * 43600., S=76.90, tau=2.0, n=0.6, V_wf=0.9, V_max=20e3, V_units='m3'
)

LimeHopper = tank_factory('LimeHopper', 
    CE=CE2007, cost=9100., S=4.02, tau=46.3, n=0.6, V_wf=0.75, V_max=100., V_units='m3'
)

AmmoniaTank = tank_factory('AmmoniaTank', 
    CE=CE2007, cost=MF38 * 28400., S=8.77, tau=100., n=0.6, V_wf=0.90, V_max=100., V_units='m3'
)

AlphaAmylaseTank = tank_factory('AlphaAmylaseTank', 
    CE=CE2007, cost=49700., S=12.77, tau=336., n=0.6, V_wf=0.90, V_max=100., V_units='m3'
)

SlurryMixTank = tank_factory('SlurryMixTank', mixing=True,
    CE=CE2007, cost=133300., S=45.29, tau=0.25, n=0.6, V_wf=0.65, V_max=80., kW_per_m3=1.1607597, V_units='m3',
)


@cost('Flow rate', units='kg/hr', CE=CE2007, cost=14000, n=0.6, S=150347)
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
        steam.imol[CAS_water] = abs(steam_mol)
        effluent.mol[:] = steam.mol + feed.mol
        effluent.H = feed.H + steam.H
        return effluent.T - T
    
    def _run(self):
        feed, steam = self._ins
        steam_mol = feed.F_mol / 100.
        effluent, = self.outs
        effluent.T = self.T
        steam_mol = flx.aitken_secant(self._T_objective_function,
                                      steam_mol, 1/8 * steam_mol + 1., 
                                      1e-4, 1e-4,
                                      args=(self.T, steam, effluent, feed),
                                      checkroot=False)
        effluent.P = steam.P / 2.
        hu, = self.heat_utilities
        hu(steam.H, effluent.T)

CookedSlurrySurgeTank = tank_factory('CookedSlurrySurgeTank',
    CE=CE2007, cost=MF90 * 173700., S=14.16, tau=0.25, n=0.6, V_wf=0.90, V_max=100., V_units='m3',
)

GlucoAmylaseTank = tank_factory('AlphaAmylaseTank', 
    CE=CE2007, cost=84200., S=17.57, tau=336., n=0.6, V_wf=0.90, V_max=100., V_units='m3'
)

SulfuricAcidTank = tank_factory('AlphaAmylaseTank', 
    CE=CE2007, cost=MF38 * 19300., S=18.87, tau=336., n=0.6, V_wf=0.90, V_max=283.17, V_units='m3'
)

Saccharification = tank_factory('AlphaAmylaseTank', kW_per_m3=0.036, mixing=True,
    CE=CE2007, cost=MF90 * 102700., S=52.10, tau=1./3., n=0.6, V_wf=0.90, V_min=20., V_max=610., V_units='m3'
)

YeastTank = tank_factory('YeastTank', kW_per_m3=0.5,
    CE=CE2007, cost=114700., S=2.97, tau=40., n=0.6, V_wf=0.90, V_max=80., V_units='m3'
)

@cost('Flow rate', S=24158., n=0.6, units='kg/hr', CE=CE2007, kW=13.1, cost=55700.)
class WetDDGSConveyor(bst.Unit): pass

@cost('Peripheral drum area', CE=CE2007, ub=7854.0,
      S=1235.35, units='m2', n=0.6, cost=MF52 * 2268000., kW=938.866)
class DDGSDryer(bst.Unit):
    """
    Create a DDGSDryer to dry dried distillers grains with solubles (DDGS) by passing
    hot air (heated by burning natural gas).
    
    Parameters
    ----------
    ins : stream sequence
        [0] Wet DDGS.
        [1] Air.
        [2] Natural gas.
    outs : stream sequence
        [0] Dried DDGS
        [1] Hot air
        [2] Emissions
    split : dict[str, float]
        Component splits to hot air (stream [1]).
    R : float, optional
        Flow of hot air over evaporation. Defaults to 1.4 wt gas / wt evap.
    H : float, optional
        Specific evaporation rate [kg/hr/m3]. Defaults to 20. 
    length_to_diameter : float, optional
        Note that the drum is horizontal. Defaults to 25.
    T : float, optional
        Operating temperature [K]. Defaults to 343.15.
    natural_gas_price : float
        Price of natural gas [USD/kg]. Defaults to 0.218.
    moisture_content : float
        Moisutre content of DDGS [wt / wt]. Defaults to 0.10.
        
    Notes
    -----
    The flow rate for air in the inlet is varied to meet the `R` specification
    (i.e. flow of hot air over flow rate evaporated). The flow rate of inlet natural
    gas is also altered to meet the heat demand.
    
    
    """
    _units = {'Evaporation': 'kg/hr',
              'Peripheral drum area': 'm2',
              'Diameter': 'm'}
    _N_ins = 3
    _N_outs = 3
    
    @property
    def isplit(self):
        """[ChemicalIndexer] Componentwise split of feed to 0th outlet stream."""
        return self._isplit
    @property
    def split(self):
        """[Array] Componentwise split of feed to 0th outlet stream."""
        return self._isplit._data
    
    @property
    def natural_gas(self):
        """[Stream] Natural gas to satisfy steam and electricity requirements."""
        return self.ins[2]
    
    @property
    def utility_cost(self):
        return super().utility_cost + self.natural_gas_cost
    
    @property
    def natural_gas_cost(self):
        return self.natural_gas_price * self.natural_gas.F_mass
    
    def __init__(self, ID="", ins=None, outs=(), thermo=None, *,
                 split, R=1.4, H=20., length_to_diameter=25, T=343.15,
                 natural_gas_price=0.289, moisture_content=0.10):
        super().__init__(ID, ins, outs, thermo)
        self._isplit = self.chemicals.isplit(split)
        self.T = T
        self.R = R
        self.H = H
        self.length_to_diameter = length_to_diameter
        self.natural_gas_price = natural_gas_price
        self.moisture_content = moisture_content
        
    def _run(self):
        wet_solids, air, natural_gas = self.ins
        dry_solids, hot_air, emissions = self.outs
        tmo.separations.split(wet_solids, hot_air, dry_solids, self.split)
        tmo.separations.adjust_moisture_content(dry_solids, hot_air, self.moisture_content)
        design_results = self.design_results
        design_results['Evaporation'] = evaporation = hot_air.F_mass
        air.imass['N2', 'O2'] = np.array([0.78, 0.32]) * self.R * evaporation
        hot_air.mol += air.mol
        dry_solids.T = hot_air.T = self.T
        duty = (dry_solids.H + hot_air.H) - (wet_solids.H + air.H)
        natural_gas.empty()
        CO2 = CH4 = - duty / self.chemicals.CH4.LHV
        H2O = 2. * CH4
        natural_gas.imol['CH4'] = CH4
        emissions.imol['CO2', 'H2O'] = [CO2, H2O]
        emissions.T = self.T + 30.
        emissions.phase = air.phase = natural_gas.phase = hot_air.phase = 'g'
        
    def _design(self):
        length_to_diameter = self.length_to_diameter
        design_results = self.design_results
        design_results['Volume'] = volume = design_results['Evaporation'] / self.H 
        design_results['Diameter'] = diameter = cylinder_diameter_from_volume(volume, length_to_diameter)
        design_results['Length'] = length = diameter * length_to_diameter
        design_results['Peripheral drum area'] = cylinder_area(diameter, length)

LiquefactionTank = tank_factory('LiquefactionTank', 
    CE=CE2007, cost=160900., S=141.3, tau=0.9, n=0.6, V_wf=0.90, V_max=500., kW_per_m3=0.6,
    mixing=True,
)

class Liquefaction(LiquefactionTank):
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
    def __init__(self, *args, yield_=1.0, **kwargs):
        super().__init__(*args, **kwargs)
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
                 tau=60.,  N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36,
                 yield_=0.9, V_wf=0.83):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
            tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax
        )
        self.reaction = tmo.reaction.Reaction('Glucose -> 2Ethanol + 2CO2',  'Glucose', yield_)
        self.V_wf = V_wf
    
    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins)
        self.reaction(effluent)
        effluent.imass['Yeast'] += effluent.imass['Glucose', 'NH3'].sum()
        effluent.imass['Glucose', 'NH3'] = 0.
        vent.receive_vent(effluent)
    
SSF = SimultaneousSaccharificationFermentation

@copy_algorithm(bst.SolidLiquidsSplitCentrifuge, run=False)
class DDGSCentrifuge(bst.Splitter): pass
    
class ThermalOxidizer(bst.Unit):
    """
    Create a ThermalOxidizer that burns any remaining combustibles.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Feed gas
        [1] Air
    outs : stream
        Emissions.
    
    Notes
    -----
    Adiabatic operation is assummed. Simulation and cost is based on [1]_.
    
    References
    ----------
    .. [1] Kwiatkowski, J. R.; McAloon, A. J.; Taylor, F.; Johnston, D. B. 
        Modeling the Process and Costs of Fuel Ethanol Production by the Corn 
        Dry-Grind Process. Industrial Crops and Products 2006, 23 (3), 288–296.
        https://doi.org/10.1016/j.indcrop.2005.08.004.

    """
    _N_ins = 2
    _N_outs = 1
    max_volume = 20. # m3
    
    def __init__(self, *args, tau=0.00014, kW_per_m3=18.47, V_wf=0.95, **kwargs):
        bst.Unit.__init__(self, *args, **kwargs)
        self.tau = tau
        self.kW_per_m3 = kW_per_m3
        self.V_wf = V_wf

    def _run(self):
        feed, air = self.ins
        emissions, = self.outs
        emissions.copy_like(feed)
        combustion_rxns = self.chemicals.get_combustion_reactions()
        combustion_rxns.force_reaction(emissions)
        O2 = max(-emissions.imol['O2'], 0.)
        emissions.copy_like(feed)
        air.imol['N2', 'O2'] = [0.78/0.32 * O2, O2]
        emissions.mol += air.mol
        combustion_rxns.adiabatic_reaction(emissions)
        
    def _design(self):
        design_results = self.design_results
        volume = self.tau * self.outs[0].F_vol / self.V_wf
        V_max = self.max_volume
        design_results['Number of vessels'] = N = np.ceil(volume / V_max)
        design_results['Vessel volume'] = volume / N
        design_results['Total volume'] = volume
        
    def _cost(self):
        design_results = self.design_results
        total_volume = design_results['Total volume']
        N = design_results['Number of vessels']
        vessel_volume = design_results['Vessel volume']
        self.power_utility.consumption = total_volume * self.kW_per_m3
        purchase_costs = self.purchase_costs
        purchase_costs['Vessels'] = N * 918300. * (vessel_volume / 13.18)**0.6
        

@cost('Flow rate', units='kg/hr', CE=CE2007, cost=122800, S=15303.5346, kW=37.3, n=0.6)
class DDGSHandling(bst.Unit): pass


class PlantAir_CIP_WasteWater_Facilities(bst.Facility):
    network_priority = 0
    
    def __init__(self, ID, corn):
        self.corn = corn
        super().__init__(ID)
        
    def _run(self):
        pass
        
    def _cost(self):
        C = self.purchase_costs
        C['Facilities'] = 6e5 * (self.corn.F_mass / 46211.6723)**0.6