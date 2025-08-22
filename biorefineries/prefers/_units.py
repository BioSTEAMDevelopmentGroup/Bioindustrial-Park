# -*- coding: utf-8 -*-
"""
Created on 2025-04-18 15:20:45

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""
# %%
import biosteam as bst
import thermosteam as tmo

from biosteam.units.decorators import cost, copy_algorithm
from biosteam.units.design_tools import CEPCI_by_year, cylinder_diameter_from_volume, cylinder_area
from biosteam import tank_factory
from thermosteam import MultiStream
from biosteam import Unit
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch
import thermosteam as tmo
import biosteam as bst
import numpy as np

import numpy as np
import flexsolve as flx

__all__ = (
    ##### Upstream #####    
    # Feedstock Preparation
    #'PretreatmentReactorSystem',

    ## glucose as beginning

    # Bioreactor Fermentaion
    'SeedTrain',
    'AeratedFermentation',

    ##### Downstream #####
    'CellDisruption',
    'ProteinCentrifuge',
    'Evaporator',
    'DiaFiltration',
    'IonExchange',
    'NanofiltrationDF',
    'SprayDrying',
    #'ScrewPress'
)
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# %% Constants
_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
# _m3hr2gpm = 4.40287
_hp2kW = 0.7457
_Gcal2kJ = 4184e3

# %% 
####################
##### UpStream #####
####################

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
    _N_outs = 2
    _ins_size_is_fixed = False
    
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
    tau_batch = 16
    
    @property
    def tau_turnover(self):
        """Turnover time (hr) calculated by batch time divided by number of trains."""
        return self.tau_batch/self.N_trains
    
    #: Operating temperature (K)
    T = 32+273.15
    
    # #: Diammonium phosphate loading in g/L of fermentation broth
    # DAP = 0.67 
    
    def _init(self, reactions=None, saccharification=None, T=None):
        chemicals = self.chemicals

        if reactions is None:
            self.reactions = ParallelRxn([
        #   Reaction definition                             Reactant    Conversion
        Rxn('Glucose -> 2 Ethanol + 2 CO2',                 'Glucose',   0.9000, chemicals),
        Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O',
                                                            'Glucose',   0.0400, chemicals),
        Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',           'Glucose',   0.0040, chemicals),
        Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',       'Glucose',   0.0060, chemicals),
        # Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                'Xylose',    0.8000, chemicals),
        # Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
        #                                                     'Xylose',    0.0400, chemicals),
        # Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',      'Xylose',    0.0030, chemicals),
        # Rxn('Xylose + H2O -> Xylitol + 0.5 O2',             'Xylose',    0.0460, chemicals),
        # Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',  'Xylose',    0.0090, chemicals)
            ])
            self.glucose_to_ethanol = self.reactions[0]
            # self.xylose_to_ethanol = self.reactions[4]
            self.glucose_to_byproducts = self.reactions[1:4]
            # self.xylose_to_byproducts = self.reactions[5:]
        else:
            self.reactions = reactions
        
        if callable(saccharification):
            self.saccharification = saccharification
        elif saccharification:
            self.saccharification = ParallelRxn([
                Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400, chemicals),
                Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120, chemicals),
                Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000, chemicals),
                Rxn('Cellobiose + H2O -> 2Glucose',       'Cellobiose',  1.0000, chemicals)]
            )
        else:
            self.saccharification = None
        if T: self.T = T
    
    _setup = Unit._setup
    
    def _run(self):
        vent, effluent= self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        if self.saccharification:
            self.saccharification(effluent)
        self.reactions.force_reaction(effluent)
        effluent.empty_negative_flows()
        effluent.T = self.T
        vent.phase = 'g'
        vent.copy_flow(effluent, ('CO2', 'O2'), remove=True)

    def _design(self): 
        maxvol = self.outs[1].F_vol*self.tau_turnover
        vol = maxvol*10**-self.N_stages
        Design = self.design_results
        for i in range(1, self.N_stages+1):
            Design[f'Stage #{i} reactor volume'] = vol
            vol *= 10 
        Design['Flow rate'] = sum([i.F_mass for i in self.outs])
        self.add_heat_utility(self.Hnet, self.T)

    def _cost(self):
        pass
        N = self.N_trains
        D = self.design_results
        C = self.baseline_purchase_costs
        kW = 0
        for i, x in self.cost_items.items():
            S = D[x._basis]
            q = S/x.S
            C[i] = N*bst.CE/x.CE*x.cost*q**x.n
            kW += N*x.kW*q
        self.power_utility(kW)

class AeratedFermentation(bst.AeratedBioreactor):
    # V_max_default = 500
    def _init(
            self, 
            fermentation_reaction, 
            cell_growth_reaction, 
            respiration_reaction,
            neutralization_reaction,
            dT_hx_loop=8,
            Q_O2_consumption=-460240, # [kJ/kmol] equivalent to 110 kcal / mol as in https://www.academia.edu/19636928/Bioreactor_Design_for_Chemical_Engineers
            batch=True,
            **kwargs,
        ):
        bst.AeratedBioreactor._init(self, batch=batch, dT_hx_loop=dT_hx_loop, 
                                    Q_O2_consumption=Q_O2_consumption,
                                    optimize_power=True, **kwargs)
        chemicals = self.chemicals
        # self.hydrolysis_reaction = Rxn('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00, chemicals)
        self.V_max_default = 500
        self.fermentation_reaction = fermentation_reaction
        self.cell_growth_reaction = cell_growth_reaction
        self.respiration_reaction = respiration_reaction
        self.neutralization_reaction = neutralization_reaction

    def _run_vent(self, vent, effluent):
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)
        assert not effluent.imol['CO2', 'O2', 'N2'].any()
    
    def _run_reactions(self, effluent):
        # self.hydrolysis_reaction.force_reaction(effluent)
        # self.lipid_reaction.force_reaction(effluent)
        if effluent.imol['H2O'] < 0.: effluent.imol['H2O'] = 0.
        self.fermentation_reaction.force_reaction(effluent)
        self.cell_growth_reaction.force_reaction(effluent)
        self.respiration_reaction.force_reaction(effluent)
        self.neutralization_reaction.force_reaction(effluent)

class PSA(bst.Flash):
    _units= {'Liquid flow': 'kg/hr'}
    
    def _run(self):
        influent = self.ins[0]
        vapor, liquid = self.outs
        
        ms = tmo.MultiStream('ms')
        ms.copy_like(influent)
        ms.vle(P=101325, H=ms.H)
        
        vapor.mol = ms.imol['g']
        vapor.phase = 'g'
        liquid.mol = ms.imol['l']
        vapor.T = liquid.T = ms.T
        vapor.P = liquid.P = ms.P
        
    def _design(self):
        self.design_results['Liquid flow'] = self.outs[1].F_mass

# %%
######################
##### Downstream #####
######################

class CellDisruption(Unit):
    """
    Cell disruption unit for breaking open cells to release intracellular components.
    Providing various methods such as 'HighPressure', 'Sonication', #'Enzymatic'

    Parameters
    ----------
    ins : Sequence[Stream]
        [0] Feed stream containing cells to be disrupted.
    outs : Sequence[Stream]
        [0] Disrupted cell stream containing released intracellular components.
    cell_disruption_method : str, optional
        Method used for cell disruption (e.g., 'HighPressure', 'Sonication', 'Enzymatic').
    cell_disruption_efficiency : float, optional
        Efficiency of the cell disruption process, defined as the fraction of cells disrupted.
    T : float, optional
        Operating temperature in Kelvin. Defaults to 310.15 K (37 °C).
    """

    _N_ins = 1
    _N_outs = 1

    _graphics = bst.Pump._graphics
    
    # Bare-module factors for equipment cost estimation
    _F_BM_default = {
        'High-Pressure Pump': 1.89,      # Typical bare-module factor for high-pressure pumps
        'High-Pressure Valve': 1.35,     # Typical bare-module factor for high-pressure valves
    }

    def __init__(self, ID='', 
                ins=None, outs=(),
                intracellular_compounds=None,
                cell_disruption_method='HighPressure',
                cell_disruption_efficiency=0.95,
                T=310.15,
                Parameters=None, **kwargs):
        bst.Unit.__init__(self, ID, ins, outs, thermo=None, **kwargs)

        self.intracellular_compounds = intracellular_compounds if intracellular_compounds is not None else ['Globin', 'Leghemoglobin']
        self.cell_disruption_method = cell_disruption_method
        self.cell_disruption_efficiency = cell_disruption_efficiency
        self.T = T
        self.Parameters = Parameters

        # Set default parameters based on the cell disruption method
        if self.cell_disruption_method == 'HighPressure':
            self.Parameters = Parameters if Parameters is not None else {
                'P High': 150e5,  # High pressure in Pascals
                'P Low': 101325,  # Low pressure in Pascals (atmospheric)
            }

        elif self.cell_disruption_method == 'BallMilling':
            self.Parameters = Parameters if Parameters is not None else {
                'Flow rate': 100e-3,  # Flow rate in m^3/s
                'Media': 'Water',
            }
        elif self.cell_disruption_method == 'Sonication':
            self.Parameters = Parameters if Parameters is not None else {
                'Frequency': 20e3,  # Sonication frequency in Hz
                'Amplitude': 100e-6,  # Sonication amplitude in meters
            }
        elif self.cell_disruption_method == 'Enzymatic':
            self.Parameters = Parameters if Parameters is not None else {
                'Enzyme': 'Cellulase',
                'Concentration': 1e-3,  # Enzyme concentration in kg/m^3
            }
        else:
            raise ValueError(f"Unsupported cell disruption method: {self.cell_disruption_method}")
        
    def _run(self):
        """ Run the cell disruption process. """
        self.outs[0].copy_like(self.ins[0])

        # Apply the cell disruption efficiency
        for compound in self.intracellular_compounds:
            self.outs[0].imol[compound+'Intre'] = (1-self.cell_disruption_efficiency) * self.ins[0].imol[compound]
            self.outs[0].imol[compound] = self.cell_disruption_efficiency * self.ins[0].imol[compound]

    def _design(self):
        """ Design the cell disruption unit based on the method and parameters. """
        Design = self.design_results
        Design['Type'] = self.cell_disruption_method
        if self.cell_disruption_method == 'HighPressure':
            # Create a temporary stream for internal calculations
            internal_stream = self.ins[0].copy()

            # Initialize and simulate the internal pump and valve
            self.pump = bst.Pump(None, None, P=self.Parameters['P High'])
            self.valve = bst.Flash(None, None, T=self.T, P=self.Parameters['P Low'])

            self.pump.ins[0] = internal_stream
            self.pump.simulate()
            self.valve.ins[0] = self.pump.outs[0]
            self.valve.simulate()

            # Design the internal units to get their costs
            self.pump._design()
            self.valve._design()

        # # incomplete yet below
        # elif self.cell_disruption_method is 'Sonication':
        #     # Design based on sonication parameters
        #     Design['Disruption efficiency'] = self.cell_disruption_efficiency
        #     Design['Frequency'] = self.Parameters['Frequency']
        #     Design['Amplitude'] = self.Parameters['Amplitude']
        # elif self.cell_disruption_method is 'Enzymatic':
        #     # Design based on enzymatic parameters
        #     Design['Disruption efficiency'] = self.cell_disruption_efficiency
        #     Design['Enzyme'] = self.Parameters['Enzyme']
        #     Design['Concentration'] = self.Parameters['Concentration']

    def _cost(self):
        """ Cost the cell disruption unit based on the method and parameters. """
        # --- Aggregate Purchase Costs ---
        # To avoid confusion, we prefix the keys.
        self.baseline_purchase_costs['High-Pressure Pump'] = self.pump.purchase_cost
        self.baseline_purchase_costs['High-Pressure Valve'] = self.valve.purchase_cost

        # --- Aggregate Utility Costs ---
        # The main cost is the electricity for the pump.
        self.power_utility.rate = self.pump.power_utility.rate

        # Add any heat utilities from the valve (e.g., Joule-Thomson cooling)
        for hu in self.valve.heat_utilities:
            self.heat_utilities.append(hu)
    
    
class ProteinCentrifuge(bst.SolidsCentrifuge): pass


class Evaporator(bst.MultiEffectEvaporator): pass


# class Diafiltration(bst.Unit):
#     """
#     Diafiltration unit for separation of solutes based on size, typically
#     retaining larger molecules (like proteins) while allowing smaller ones
#     (like salts and water) to pass through the permeate. Includes continuous
#     addition of wash solution (diafiltration buffer).

#     Parameters
#     ----------
#     ins : Sequence[Stream]
#         [0] Feed stream to be processed.
#         [1] Wash solution (diafiltration buffer).
#     outs : Sequence[Stream]
#         [0] Permeate stream (water, salts, smaller molecules).
#         [1] Retentate stream (concentrated target product).
#     TargetProduct_ID : str, optional
#         Chemical ID of the target product to be retained.
#         Defaults to 'TargetProduct'.
#     Salt_ID : str or list[str], optional
#         Chemical ID(s) for salts that are mostly permeated.
#         Defaults to 'Salt'.
#     OtherLargeMolecules_ID : str or list[str], optional
#         Chemical ID(s) for other large molecules that are mostly retained.
#         Defaults to 'OtherLargeMolecules'.
#     TargetProduct_Retention : float, optional
#         Fraction of the target product retained in the retentate.
#         Defaults to 0.98.
#     Salt_Retention : float, optional
#         Fraction of salt(s) retained in the retentate.
#         Defaults to 0.05.
#     OtherLargeMolecules_Retention : float, optional
#         Fraction of other large molecules retained in the retentate.
#         Defaults to 0.99.
#     DefaultSolutes_Retention : float, optional
#         Fraction of any other solutes (not specified above) retained.
#         Defaults to 0.05.
#     FeedWater_Recovery_to_Permeate : float, optional
#         Fraction of water from the initial feed stream that is recovered in the permeate.
#         The rest of the water (including all wash water) is balanced.
#         Defaults to 0.75.
#     membrane_flux_LMH : float, optional
#         Average design membrane flux in Liters per square meter per hour (LMH).
#         Defaults to 50.0.
#     TMP_bar : float, optional
#         Transmembrane pressure in bar. Used for pump power calculation.
#         Defaults to 2.0.
#     pump_efficiency : float, optional
#         Efficiency of the pump(s) for the diafiltration system.
#         Defaults to 0.75.
#     membrane_cost_USD_per_m2 : float, optional
#         Replacement cost of membranes in USD per square meter.
#         Defaults to 200.0.
#     membrane_lifetime_years : float, optional
#         Expected lifetime of the membranes in years.
#         Defaults to 2.0.
#     module_cost_factor : float, optional
#         Cost factor for membrane system purchase cost calculation.
#         Assumes Cost = factor * (Area_m2 ** exponent).
#         Defaults to 2500.0 (USD for a base CEPCI).
#     module_cost_exponent : float, optional
#         Exponent for membrane system purchase cost calculation based on area.
#         Defaults to 0.7.
#     base_CEPCI : float, optional
#         The Chemical Engineering Plant Cost Index (CEPCI) for which the
#         `module_cost_factor` is valid.
#         Defaults to 500.
#     """
    
#     _N_ins = 2
#     _N_outs = 2
    
#     # Bare-module factors for equipment cost estimation
#     _F_BM_default = {
#         'Membrane System': 1.65,          # Typical bare-module factor for membrane modules
#         'Membrane replacement': 1.0,      # Operating cost item (no installation factor)
#         'Pump': 1.89,                     # Typical bare-module factor for centrifugal pumps
#     }

#     water_ID = 'H2O'
#     _default_TargetProduct_ID = 'TargetProduct'
#     _default_Salt_ID = 'Salt'
#     _default_OtherLargeMolecules_ID = 'OtherLargeMolecules'

#     _default_TargetProduct_Retention = 0.98
#     _default_Salt_Retention = 0.05
#     _default_OtherLargeMolecules_Retention = 0.99
#     _default_DefaultSolutes_Retention = 0.08
#     _default_FeedWater_Recovery_to_Permeate = 0.75

#     _default_membrane_flux_LMH = 50.0
#     _default_TMP_bar = 2.0
#     _default_membrane_cost_USD_per_m2 = 200.0
#     _default_membrane_lifetime_years = 2.0
#     _default_module_cost_factor = 2500.0 # e.g., USD for CEPCI=500
#     _default_module_cost_exponent = 0.7
#     _default_base_CEPCI = 500.0

#     _units = {
#         'Membrane Area': 'm2',
#         'TargetProduct_Retention': '%',
#         'Salt_Retention': '%',
#         'OtherLargeMolecules_Retention': '%',
#         'DefaultSolutes_Retention': '%',
#         'FeedWater_Recovery_to_Permeate': '%',
#         'membrane_flux_LMH': 'LMH',
#         'TMP_bar': 'bar',
#         'pump_efficiency': '%',
#         'membrane_cost_USD_per_m2': '$/m2',
#         'membrane_lifetime_years': 'years',
#         'module_cost_factor': '$/m2^exponent',
#         'module_cost_exponent': '',
#         'base_CEPCI': '',
#     }
    
#     def __init__(self, ID='', ins=None, outs=None, thermo=None,
#                 TargetProduct_ID=None, Salt_ID=None, OtherLargeMolecules_ID=None,
#                 TargetProduct_Retention=None, Salt_Retention=None,
#                 OtherLargeMolecules_Retention=None, DefaultSolutes_Retention=None,
#                 FeedWater_Recovery_to_Permeate=None,
#                 membrane_flux_LMH=None, TMP_bar=None,
#                 membrane_cost_USD_per_m2=None, membrane_lifetime_years=None,
#                 module_cost_factor=None, module_cost_exponent=None, base_CEPCI=None,
#                  **kwargs):
#         super().__init__(ID, ins, outs, thermo)

#         self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
#         self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
#         self.OtherLargeMolecules_ID = OtherLargeMolecules_ID if OtherLargeMolecules_ID is not None else self._default_OtherLargeMolecules_ID

#         self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
#         self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
#         self.OtherLargeMolecules_Retention = OtherLargeMolecules_Retention if OtherLargeMolecules_Retention is not None else self._default_OtherLargeMolecules_Retention
#         self.DefaultSolutes_Retention = DefaultSolutes_Retention if DefaultSolutes_Retention is not None else self._default_DefaultSolutes_Retention
#         self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate

#         self.membrane_flux_LMH = membrane_flux_LMH if membrane_flux_LMH is not None else self._default_membrane_flux_LMH
#         self.TMP_bar = TMP_bar if TMP_bar is not None else self._default_TMP_bar
#         self.membrane_cost_USD_per_m2 = membrane_cost_USD_per_m2 if membrane_cost_USD_per_m2 is not None else self._default_membrane_cost_USD_per_m2
#         self.membrane_lifetime_years = membrane_lifetime_years if membrane_lifetime_years is not None else self._default_membrane_lifetime_years
#         self.module_cost_factor = module_cost_factor if module_cost_factor is not None else self._default_module_cost_factor
#         self.module_cost_exponent = module_cost_exponent if module_cost_exponent is not None else self._default_module_cost_exponent
#         self.base_CEPCI = base_CEPCI if base_CEPCI is not None else self._default_base_CEPCI
        
#         self.power_utility = bst.PowerUtility()

#     def _run(self):
#         feed, wash_solution = self.ins
#         retentate, permeate = self.outs

#         # Assume isothermal operation and outlet pressures equal to feed pressure
#         # Actual pressure drop for pumping is handled by TMP_bar in _design
#         retentate.T = permeate.T = feed.T
#         retentate.P = permeate.P = feed.P # This is outlet P, TMP is an internal driving force

#         permeate.empty()
#         retentate.empty()

#         # --- Water Balance ---
#         feed_water_mass = feed.imass[self.water_ID]
#         wash_water_mass = wash_solution.imass[self.water_ID]
#         total_incoming_water = feed_water_mass + wash_water_mass

#         retentate_water_from_feed = feed_water_mass * (1.0 - self.FeedWater_Recovery_to_Permeate)
#         retentate.imass[self.water_ID] = max(0.0, retentate_water_from_feed)
        
#         permeate_water_total = total_incoming_water - retentate.imass[self.water_ID]
#         permeate.imass[self.water_ID] = max(0.0, permeate_water_total)

#         # Ensure overall water balance due to max(0,...) or floating point precision
#         current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
#         if abs(current_total_water_out - total_incoming_water) > 1e-9: # Tolerance
#             if total_incoming_water >= retentate.imass[self.water_ID]:
#                 permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
#             else: # This case should ideally not be reached with the logic above
#                 retentate.imass[self.water_ID] = total_incoming_water
#                 permeate.imass[self.water_ID] = 0.0

#         # --- Solute Balance ---
#         for chem in self.chemicals:
#             ID = chem.ID
#             if ID == self.water_ID:
#                 continue

#             mass_in_feed = feed.imass[ID]
#             mass_in_wash = wash_solution.imass[ID] # Solutes might be in wash solution
#             total_mass_in = mass_in_feed + mass_in_wash

#             if total_mass_in <= 1e-12: # Effectively zero
#                 retentate.imass[ID] = 0.0
#                 permeate.imass[ID] = 0.0
#                 continue

#             current_retention = self.DefaultSolutes_Retention
#             if self.TargetProduct_ID:
#                 if isinstance(self.TargetProduct_ID, str) and ID == self.TargetProduct_ID:
#                     current_retention = self.TargetProduct_Retention
#                 if isinstance(self.TargetProduct_ID, list) and ID in self.TargetProduct_ID:
#                     current_retention = self.TargetProduct_Retention
#             elif self.Salt_ID:
#                 if isinstance(self.Salt_ID, str) and ID == self.Salt_ID:
#                     current_retention = self.Salt_Retention
#                 elif isinstance(self.Salt_ID, list) and ID in self.Salt_ID: # Handles list of salt IDs
#                     current_retention = self.Salt_Retention
#             elif self.OtherLargeMolecules_ID:
#                 if isinstance(self.OtherLargeMolecules_ID, str) and ID == self.OtherLargeMolecules_ID:
#                     current_retention = self.OtherLargeMolecules_Retention
#                 elif isinstance(self.OtherLargeMolecules_ID, list) and ID in self.OtherLargeMolecules_ID: # Handles list
#                     current_retention = self.OtherLargeMolecules_Retention
            
#             retentate_mass_solute = total_mass_in * current_retention
#             retentate.imass[ID] = max(0.0, retentate_mass_solute)
            
#             permeate_mass_solute = total_mass_in - retentate.imass[ID] # Permeate by difference for mass balance
#             permeate.imass[ID] = max(0.0, permeate_mass_solute)
            
#             # Final check for solute mass balance due to max(0,...) or floating point nuances
#             current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
#             mass_balance_error = current_total_solute_out - total_mass_in
#             # Check relative and absolute error
#             if abs(mass_balance_error) > (1e-9 * abs(total_mass_in) + 1e-12):
#                 permeate.imass[ID] -= mass_balance_error # Adjust permeate
#                 if permeate.imass[ID] < 0:
#                     retentate.imass[ID] += permeate.imass[ID] # Add deficit to retentate
#                     permeate.imass[ID] = 0.0
#                     if retentate.imass[ID] < 0: # Should not happen if total_mass_in >=0
#                         retentate.imass[ID] = 0.0
#                         # Consider logging a warning here if mass is lost.
#                         # print(f"Warning: Mass balance issue for {ID} in {self.ID}. Input: {total_mass_in:.2e}, Output: {current_total_solute_out - mass_balance_error:.2e}")

#     def _design(self):
#         Design = self.design_results
#         # --- Membrane Area Calculation ---
#         permeate_stream = self.outs[0]
#         # Calculate permeate volumetric flow rate (L/hr)
#         # F_mass (kg/hr) / rho (kg/m^3) -> m^3/hr. Then * 1000 for L/hr.
#         if permeate_stream.isempty() or permeate_stream.rho == 0:
#             permeate_vol_L_per_hr = 0.0
#         else:
#             permeate_vol_L_per_hr = (permeate_stream.F_mass / permeate_stream.rho) * 1000.0
            
#         if self.membrane_flux_LMH > 0 and permeate_vol_L_per_hr > 0:
#             membrane_area_m2 = permeate_vol_L_per_hr / self.membrane_flux_LMH
#         else:
#             membrane_area_m2 = 0.0
#         Design['Membrane Area'] = membrane_area_m2
#         # Design['TargetProduct_Retention'] = self.TargetProduct_Retention * 100.0
#         # Design['Salt_Retention'] = self.Salt_Retention * 100.0
#         # Design['OtherLargeMolecules_Retention'] = self.OtherLargeMolecules_Retention * 100.0
#         # Design['DefaultSolutes_Retention'] = self.DefaultSolutes_Retention * 100.0
#         # Design['FeedWater_Recovery_to_Permeate'] = self.FeedWater_Recovery_to_Permeate * 100.0
#         Design['membrane_flux_LMH'] = self.membrane_flux_LMH
#         Design['TMP_bar'] = self.TMP_bar
#         Design['membrane_cost_USD_per_m2'] = self.membrane_cost_USD_per_m2
#         Design['membrane_lifetime_years'] = self.membrane_lifetime_years
#         # Design['module_cost_factor'] = self.module_cost_factor
#         # Design['module_cost_exponent'] = self.module_cost_exponent
#         # Design['base_CEPCI'] = self.base_CEPCI

#         # --- Pump Power Calculation ---
#         # Total volumetric flow to be pumped (feed + wash solution) in m^3/hr
#         internal_stream = self.ins[0].copy()
#         self.pump = bst.Pump(None, None, P=self.TMP_bar * 1e5)
#         self.pump.ins[0] = internal_stream
#         self.pump.simulate()
#         self.pump._design()  # Design the pump to get its cost
#         self.power_utility = self.pump.power_utility  # Use the pump's power utility
#         Design['pump_efficiency'] = self.pump.design_results['Efficiency'] * 100.0 

#     def _cost(self):
#         super()._cost()  # Call parent cost method to initialize purchase_costs and add_OPEX
#         # --- Capital Cost (Purchase Cost) ---
#         area_m2 = self.design_results.get('Membrane Area', 0.0)

#         if area_m2 > 0 and self.module_cost_factor > 0 and self.base_CEPCI > 0:
#             # Calculate base purchase cost using the power law
#             base_purchase_cost = self.module_cost_factor * (area_m2 ** self.module_cost_exponent)
#             # Adjust cost from base_CEPCI to current BioSTEAM CEPCI (bst.CE)
#             current_purchase_cost = base_purchase_cost * (bst.CE / self.base_CEPCI)
#             self.baseline_purchase_costs['Membrane System'] = current_purchase_cost
#         else:
#             self.baseline_purchase_costs['Membrane System'] = 0.0

#         # --- Annual Operating Cost (OPEX) for Membrane Replacement ---
#         # This is added to `add_OPEX` for the TEA.
#         if (self.membrane_lifetime_years > 0 and
#             self.membrane_cost_USD_per_m2 > 0 and
#             area_m2 > 0):
#             annual_replacement_cost = (area_m2 * self.membrane_cost_USD_per_m2) / self.membrane_lifetime_years
#             self.baseline_purchase_costs['Membrane replacement'] = annual_replacement_cost
#         else:
#             self.baseline_purchase_costs['Membrane replacement'] = 0.0
#         #self.power_utility.cost = self.power_utility.rate * bst.annual_hours * bst.electricity_cost_per_kWh
#         self.baseline_purchase_costs['Pump'] = self.pump.purchase_cost
#         # Tank costs are not included here.



class Diafiltration(bst.Unit):
    """
    Diafiltration unit for separation of solutes based on size, typically
    retaining larger molecules (like proteins) while allowing smaller ones
    (like salts and water) to pass through the permeate. Includes continuous
    addition of wash solution (diafiltration buffer).

    Parameters
    ----------
    ins : Sequence[Stream]
        [0] Feed stream to be processed.
        [1] Wash solution (diafiltration buffer).
    outs : Sequence[Stream]
        [0] Retentate stream (concentrated target product).
        [1] Permeate stream (water, salts, smaller molecules).
    TargetProduct_ID : str, optional
        Chemical ID of the target product to be retained.
        Defaults to 'TargetProduct'.
    Salt_ID : str or list[str], optional
        Chemical ID(s) for salts that are mostly permeated.
        Defaults to 'Salt'.
    OtherLargeMolecules_ID : str or list[str], optional
        Chemical ID(s) for other large molecules that are mostly retained.
        Defaults to 'OtherLargeMolecules'.
    TargetProduct_Retention : float, optional
        Fraction of the target product retained in the retentate.
        Defaults to 0.98.
    Salt_Retention : float, optional
        Fraction of salt(s) retained in the retentate.
        Defaults to 0.05.
    OtherLargeMolecules_Retention : float, optional
        Fraction of other large molecules retained in the retentate.
        Defaults to 0.99.
    DefaultSolutes_Retention : float, optional
        Fraction of any other solutes (not specified above) retained.
        Defaults to 0.05.
    FeedWater_Recovery_to_Permeate : float, optional
        Fraction of water from the initial feed stream that is recovered in the permeate.
        The rest of the water (including all wash water) is balanced.
        Defaults to 0.75.
    membrane_flux_LMH : float, optional
        Average design membrane flux in Liters per square meter per hour (LMH).
        Defaults to 50.0.
    TMP_bar : float, optional
        Transmembrane pressure in bar. Used for pump power calculation.
        Defaults to 2.0.
    pump_efficiency : float, optional
        Efficiency of the pump(s) for the diafiltration system.
        Defaults to 0.75.
    membrane_cost_USD_per_m2 : float, optional
        Replacement cost of membranes in USD per square meter.
        Defaults to 200.0.
    membrane_lifetime_years : float, optional
        Expected lifetime of the membranes in years.
        Defaults to 2.0.
    module_cost_factor : float, optional
        Cost factor for membrane system purchase cost calculation.
        Assumes Cost = factor * (Area_m2 ** exponent).
        Defaults to 2500.0 (USD for a base CEPCI).
    module_cost_exponent : float, optional
        Exponent for membrane system purchase cost calculation based on area.
        Defaults to 0.7.
    base_CEPCI : float, optional
        The Chemical Engineering Plant Cost Index (CEPCI) for which the
        `module_cost_factor` is valid.
        Defaults to 500.
    """
    
    _N_ins = 2
    _N_outs = 2
    
    # Bare-module factors for equipment cost estimation
    _F_BM_default = {
        'Membrane System': 1.65,          # Typical bare-module factor for membrane modules
        'Membrane replacement': 1.0,      # Operating cost item (no installation factor)
        'Pump': 1.89,                     # Typical bare-module factor for centrifugal pumps
    }

    water_ID = 'H2O'
    _default_TargetProduct_ID = 'TargetProduct'
    _default_Salt_ID = 'Salt'
    _default_OtherLargeMolecules_ID = 'OtherLargeMolecules'

    _default_TargetProduct_Retention = 0.98
    _default_Salt_Retention = 0.05
    _default_OtherLargeMolecules_Retention = 0.99
    _default_DefaultSolutes_Retention = 0.08
    _default_FeedWater_Recovery_to_Permeate = 0.75

    _default_membrane_flux_LMH = 50.0
    _default_TMP_bar = 2.0
    _default_membrane_cost_USD_per_m2 = 200.0
    _default_membrane_lifetime_years = 1.0
    _default_module_cost_factor = 2500.0 # e.g., USD for CEPCI=500
    _default_module_cost_exponent = 0.7
    _default_base_CEPCI = 500.0

    _units = {
        'Membrane Area': 'm2',
        'TargetProduct_Retention': '%',
        'Salt_Retention': '%',
        'OtherLargeMolecules_Retention': '%',
        'DefaultSolutes_Retention': '%',
        'FeedWater_Recovery_to_Permeate': '%',
        'membrane_flux_LMH': 'LMH',
        'TMP_bar': 'bar',
        'pump_efficiency': '%',
        'membrane_cost_USD_per_m2': '$/m2',
        'membrane_lifetime_years': 'years',
        'module_cost_factor': '$/m2^exponent',
        'module_cost_exponent': '',
        'base_CEPCI': '',
    }
    
    def __init__(self, ID='', ins=None, outs=None, thermo=None,
                TargetProduct_ID=None, Salt_ID=None, OtherLargeMolecules_ID=None,
                TargetProduct_Retention=None, Salt_Retention=None,
                OtherLargeMolecules_Retention=None, DefaultSolutes_Retention=None,
                FeedWater_Recovery_to_Permeate=None,
                membrane_flux_LMH=None, TMP_bar=None,
                membrane_cost_USD_per_m2=None, membrane_lifetime_years=None,
                module_cost_factor=None, module_cost_exponent=None, base_CEPCI=None,
                 **kwargs):
        super().__init__(ID, ins, outs, thermo)

        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
        self.OtherLargeMolecules_ID = OtherLargeMolecules_ID if OtherLargeMolecules_ID is not None else self._default_OtherLargeMolecules_ID

        self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
        self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
        self.OtherLargeMolecules_Retention = OtherLargeMolecules_Retention if OtherLargeMolecules_Retention is not None else self._default_OtherLargeMolecules_Retention
        self.DefaultSolutes_Retention = DefaultSolutes_Retention if DefaultSolutes_Retention is not None else self._default_DefaultSolutes_Retention
        self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate

        self.membrane_flux_LMH = membrane_flux_LMH if membrane_flux_LMH is not None else self._default_membrane_flux_LMH
        self.TMP_bar = TMP_bar if TMP_bar is not None else self._default_TMP_bar
        self.membrane_cost_USD_per_m2 = membrane_cost_USD_per_m2 if membrane_cost_USD_per_m2 is not None else self._default_membrane_cost_USD_per_m2
        self.membrane_lifetime_years = membrane_lifetime_years if membrane_lifetime_years is not None else self._default_membrane_lifetime_years
        self.module_cost_factor = module_cost_factor if module_cost_factor is not None else self._default_module_cost_factor
        self.module_cost_exponent = module_cost_exponent if module_cost_exponent is not None else self._default_module_cost_exponent
        self.base_CEPCI = base_CEPCI if base_CEPCI is not None else self._default_base_CEPCI
        
        self.power_utility = bst.PowerUtility()

    def _run(self):
        feed, wash_solution = self.ins
        retentate, permeate = self.outs

        # Assume isothermal operation
        retentate.T = permeate.T = feed.T
        retentate.P = permeate.P = feed.P

        permeate.empty()
        retentate.copy_like(feed) # Start with feed composition in retentate

        # --- Water Balance ---
        # This part of your logic is kept, as it defines the final retentate volume.
        feed_water_mass = feed.imass[self.water_ID]
        wash_water_mass = wash_solution.imass[self.water_ID]
        total_incoming_water = feed_water_mass + wash_water_mass

        retentate_water_from_feed = feed_water_mass * (1.0 - self.FeedWater_Recovery_to_Permeate)
        retentate.imass[self.water_ID] = max(0.0, retentate_water_from_feed)
        
        permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
        if permeate.imass[self.water_ID] < 0: # Safety check
            permeate.imass[self.water_ID] = 0.0
            retentate.imass[self.water_ID] = total_incoming_water

        # --- Solute Balance ---
        # Create dictionaries for easier lookup of retention values
        retention_map = {}
        # Ensure IDs are in lists for consistent processing
        target_ids = self.TargetProduct_ID if isinstance(self.TargetProduct_ID, list) else [self.TargetProduct_ID]
        large_mol_ids = self.OtherLargeMolecules_ID if isinstance(self.OtherLargeMolecules_ID, list) else [self.OtherLargeMolecules_ID]
        salt_ids = self.Salt_ID if isinstance(self.Salt_ID, list) else [self.Salt_ID]

        for chem_id in target_ids: retention_map[chem_id] = self.TargetProduct_Retention
        for chem_id in large_mol_ids: retention_map[chem_id] = self.OtherLargeMolecules_Retention
        for chem_id in salt_ids: retention_map[chem_id] = self.Salt_Retention

        for chem in self.chemicals:
            ID = chem.ID
            if ID == self.water_ID:
                continue

            total_mass_in = feed.imass[ID] + wash_solution.imass[ID]
            if total_mass_in < 1e-12: continue

            # Get the retention for the current chemical, or use the default
            current_retention = retention_map.get(ID, self.DefaultSolutes_Retention)

            # --- NEW LOGIC: Differentiate between retained and permeable species ---
            # We use a threshold (e.g., 0.5) to decide the model chemistry.
            # High retention means the species is mostly retained (e.g., protein, cells).
            if current_retention > 0.5:
                # MODEL 1: RETAINED SPECIES
                # The final amount is a fraction of the total input.
                retentate.imass[ID] = total_mass_in * current_retention
                permeate.imass[ID] = total_mass_in - retentate.imass[ID]
            else:
                # MODEL 2: PERMEABLE SPECIES (BUFFER EXCHANGE)
                # The final concentration in the retentate matches the wash buffer.
                wash_water = wash_solution.imass[self.water_ID]
                if wash_water > 1e-9:
                    # Concentration in wash buffer (per kg of water)
                    conc_in_wash = wash_solution.imass[ID] / wash_water
                    # Set retentate mass based on this concentration and the retentate's water content
                    retentate.imass[ID] = conc_in_wash * retentate.imass[self.water_ID]
                else:
                    # If wash solution has no water, this solute can't be in the retentate liquid phase
                    retentate.imass[ID] = 0.0

                # Permeate is the remainder to ensure mass balance
                permeate.imass[ID] = total_mass_in - retentate.imass[ID]

            # Final safety check to prevent negative mass flows
            if permeate.imass[ID] < 0:
                retentate.imass[ID] += permeate.imass[ID] # Adjust retentate
                permeate.imass[ID] = 0.0
            if retentate.imass[ID] < 0:
                permeate.imass[ID] += retentate.imass[ID] # Adjust permeate
                retentate.imass[ID] = 0.0

    def _design(self):
        Design = self.design_results
        # --- Membrane Area Calculation ---
        permeate_stream = self.outs[1]
        # Calculate permeate volumetric flow rate (L/hr)
        # F_mass (kg/hr) / rho (kg/m^3) -> m^3/hr. Then * 1000 for L/hr.
        if permeate_stream.isempty() or permeate_stream.rho == 0:
            permeate_vol_L_per_hr = 0.0
        else:
            permeate_vol_L_per_hr = (permeate_stream.F_mass / permeate_stream.rho) * 1000.0
            
        if self.membrane_flux_LMH > 0 and permeate_vol_L_per_hr > 0:
            membrane_area_m2 = permeate_vol_L_per_hr / self.membrane_flux_LMH
        else:
            membrane_area_m2 = 0.0
        Design['Membrane Area'] = membrane_area_m2
        # Design['TargetProduct_Retention'] = self.TargetProduct_Retention * 100.0
        # Design['Salt_Retention'] = self.Salt_Retention * 100.0
        # Design['OtherLargeMolecules_Retention'] = self.OtherLargeMolecules_Retention * 100.0
        # Design['DefaultSolutes_Retention'] = self.DefaultSolutes_Retention * 100.0
        # Design['FeedWater_Recovery_to_Permeate'] = self.FeedWater_Recovery_to_Permeate * 100.0
        Design['membrane_flux_LMH'] = self.membrane_flux_LMH
        Design['TMP_bar'] = self.TMP_bar
        Design['membrane_cost_USD_per_m2'] = self.membrane_cost_USD_per_m2
        Design['membrane_lifetime_years'] = self.membrane_lifetime_years
        # Design['module_cost_factor'] = self.module_cost_factor
        # Design['module_cost_exponent'] = self.module_cost_exponent
        # Design['base_CEPCI'] = self.base_CEPCI

        # --- Pump Power Calculation ---
        # Total volumetric flow to be pumped (feed + wash solution) in m^3/hr
        internal_stream = self.ins[0].copy() + self.ins[1].copy()
        self.pump = bst.Pump(None, None, P=self.TMP_bar * 1e5)
        self.pump.ins[0] = internal_stream
        self.pump.simulate()
        self.pump._design()  # Design the pump to get its cost
        self.power_utility = self.pump.power_utility + self.Membrane_replacement # Use the pump's power utility
        Design['pump_efficiency'] = self.pump.design_results['Efficiency'] * 100.0

    def _cost(self):
        super()._cost()  # Call parent cost method to initialize purchase_costs and add_OPEX
        # --- Capital Cost (Purchase Cost) ---
        area_m2 = self.design_results.get('Membrane Area', 0.0)

        if area_m2 > 0 and self.module_cost_factor > 0 and self.base_CEPCI > 0:
            # Calculate base purchase cost using the power law
            base_purchase_cost = self.module_cost_factor * (area_m2 ** self.module_cost_exponent)
            # Adjust cost from base_CEPCI to current BioSTEAM CEPCI (bst.CE)
            current_purchase_cost = base_purchase_cost * (bst.CE / self.base_CEPCI)
            self.baseline_purchase_costs['Membrane System'] = current_purchase_cost
        else:
            self.baseline_purchase_costs['Membrane System'] = 0.0

        # --- Annual Operating Cost (OPEX) for Membrane Replacement ---
        # This is added to `add_OPEX` for the TEA.
        if (self.membrane_lifetime_years > 0 and
            self.membrane_cost_USD_per_m2 > 0 and
            area_m2 > 0):
            annual_replacement_cost = (area_m2 * self.membrane_cost_USD_per_m2) / self.membrane_lifetime_years
            self.baseline_purchase_costs['Membrane replacement'] = annual_replacement_cost
        else:
            self.baseline_purchase_costs['Membrane replacement'] = 0.0
        #self.power_utility.cost = self.power_utility.rate * bst.annual_hours * bst.electricity_cost_per_kWh
        self.baseline_purchase_costs['Pump'] = self.pump.purchase_cost
        # Tank costs are not included here.



class IonExchange(bst.Unit): 

    _N_ins = 2  # Feed (conditioned) and Elution Buffer Profile
    _N_outs = 2 # Product (in elution buffer) and Waste Stream
    
    # Bare-module factors for equipment cost estimation
    _F_BM_default = {
        'IEX Resin': 1.0,                 # Resin is typically purchased without installation factor
        'IEX Column Hardware': 2.6,       # Typical bare-module factor for pressure vessels/columns
    }
    
    _units = {'resin_DBC_g_L': 'g/L', 
            'load_safety_factor': 'unitless',
            'operating_pressure_bar': 'bar', 
            'pump_efficiency': 'unitless',
            'resin_cost_USD_per_L': 'USD/L', 
            'column_hardware_cost_factor': 'unitless',
            'column_hardware_cost_exponent': 'unitless', 
            'resin_lifetime_years': 'years'}
    # --- Default Values for Parameters ---

    _default_water_recovery = 0.15 # 15% recovery of feed water into product stream

    _default_water_ID = 'H2O' # Using your preferred ID for water
    _default_TargetProduct_ID = 'Leghemoglobin' # As per your context

    _default_TargetProduct_Yield = 0.95      # e.g., 95% recovery of loaded target in product pool

    # Impurities specifically targeted for removal by binding differently than the product
    _default_BoundImpurity_ID = ('HostCellProtein', 'DNA', 'Endotoxin') # Example impurity IDs
    # Fraction of these 'Bound Impurities' (from feed) that are successfully removed from the product path (i.e., go to waste)
    _default_BoundImpurity_Removal_Efficiency = 0.99 # e.g., 99% removal (2 LRV)

    # For other solutes in the feed (not TargetProduct, not BoundImpurity, not ElutionBufferSalt)
    # This fraction of these other solutes (from feed) ends up in the product stream. The rest goes to waste.
    _default_NonBinding_Solutes_Carryover_to_Product = 0.05 
    
    # Key components defining the elution buffer matrix (e.g., the elution salt from ins[1])
    # If these components are also present in the feed (ins[0]), the feed's portion is assumed to go to waste.
    _default_ElutionBuffer_Defining_Component_IDs_tuple = ('NaCl', 'KCl') # Example elution salts

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                water_recovery=None,
                water_ID=None,
                TargetProduct_ID=None,
                TargetProduct_Yield=None,
                BoundImpurity_ID=None,
                BoundImpurity_Removal_Efficiency=None,
                NonBinding_Solutes_Carryover_to_Product=None,
                ElutionBuffer_Defining_Component_ID=None,
                **kwargs):
        super().__init__(ID, ins, outs, thermo, **kwargs)

        # --- Initialize the new water recovery attribute ---
        self.water_recovery = water_recovery if water_recovery is not None else self._default_water_recovery

        # Set other operational parameters
        self.water_ID = water_ID if water_ID is not None else self._default_water_ID
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.TargetProduct_Yield = TargetProduct_Yield if TargetProduct_Yield is not None else self._default_TargetProduct_Yield
        self.BoundImpurity_ID = BoundImpurity_ID if BoundImpurity_ID is not None else self._default_BoundImpurity_ID
        self.BoundImpurity_Removal_Efficiency = BoundImpurity_Removal_Efficiency if BoundImpurity_Removal_Efficiency is not None else self._default_BoundImpurity_Removal_Efficiency
        self.NonBinding_Solutes_Carryover_to_Product = NonBinding_Solutes_Carryover_to_Product if NonBinding_Solutes_Carryover_to_Product is not None else self._default_NonBinding_Solutes_Carryover_to_Product
        self.ElutionBuffer_Defining_Component_ID = ElutionBuffer_Defining_Component_ID if ElutionBuffer_Defining_Component_ID is not None else self._default_ElutionBuffer_Defining_Component_ID

    def _run(self):
        feed = self.ins[0]
        elution_buffer_profile = self.ins[1]
        product_stream = self.outs[0]
        waste_stream = self.outs[1]

        # Initialize output streams
        product_stream.T = elution_buffer_profile.T
        waste_stream.T = feed.T
        product_stream.mass = elution_buffer_profile.mass.copy()
        waste_stream.empty()

        # --- Updated Water Balance Logic ---
        # Partition water from the feed based on the recovery ratio
        feed_water_mass = feed.imass[self.water_ID]
        water_to_product = feed_water_mass * self.water_recovery
        water_to_waste = feed_water_mass - water_to_product

        product_stream.imass[self.water_ID] += water_to_product
        waste_stream.imass[self.water_ID] += water_to_waste

        # Process all other components from the feed stream
        for chem in self.chemicals:
            ID = chem.ID
            # Skip water as it's already handled
            if ID == self.water_ID:
                continue

            feed_mass_component = feed.imass[ID]

            if feed_mass_component <= 1e-12:
                continue

            # Route components based on their type
            if ID == self.TargetProduct_ID:
                mass_to_product = feed_mass_component * self.TargetProduct_Yield
                product_stream.imass[ID] += mass_to_product
                waste_stream.imass[ID] = feed_mass_component - mass_to_product

            elif self.ElutionBuffer_Defining_Component_ID and ID in self.ElutionBuffer_Defining_Component_ID:
                # Elution components from the feed go to waste
                waste_stream.imass[ID] = feed_mass_component

            elif self.BoundImpurity_ID and ID in self.BoundImpurity_ID:
                mass_removed_to_waste = feed_mass_component * self.BoundImpurity_Removal_Efficiency
                waste_stream.imass[ID] = mass_removed_to_waste
                product_stream.imass[ID] += (feed_mass_component - mass_removed_to_waste)

            else:
                # Partition all other non-binding solutes
                mass_to_product = feed_mass_component * self.NonBinding_Solutes_Carryover_to_Product
                product_stream.imass[ID] += mass_to_product
                waste_stream.imass[ID] = feed_mass_component - mass_to_product


    def _design(self):
        """
        Placeholder for Ion Exchange Column design.
        Key design parameters would include resin volume, column dimensions.
        """
        # Example: Calculate resin volume based on Dynamic Binding Capacity (DBC)
        # Parameters needed for this (to be added to __init__ if implementing):
        self.resin_DBC_g_L = 50 # g of TargetProduct per L of resin (e.g.)
        self.DBC_safety_factor = 0.8 # Operate at 80% of DBC
        self.load_flow_rate_CV_hr = 5 # Column Volumes per hour for loading
        self.num_cycles_per_year = 300 # For equipment sizing based on annual throughput

        target_product_in_feed_kg_hr = self.ins[0].imass[self.TargetProduct_ID] # If continuous average
        if hasattr(self, 'resin_DBC_g_L') and self.resin_DBC_g_L > 0:
            effective_DBC_kg_L = (self.resin_DBC_g_L / 1000.0) * self.DBC_safety_factor
            if effective_DBC_kg_L > 0:
                # This calculation depends on whether flow rate is per hour or per batch
                # For a batch process, target_product_in_feed_kg_hr would be kg/batch
                resin_volume_L = (target_product_in_feed_kg_hr / self.load_flow_rate_CV_hr) / effective_DBC_kg_L
                self.design_results['Resin Volume (L)'] = resin_volume_L
        
        # --- Pump Power Calculation ---
        # Total volumetric flow to be pumped (feed + wash solution) in m^3/hr
        internal_stream = self.ins[0].copy()
        self.pump = bst.Pump(None, None, P=2 * 1e5)
        self.pump.ins[0] = internal_stream
        self.pump.simulate()
        self.pump._design()  # Design the pump to get its cost
        self.power_utility = self.pump.power_utility  # Use the pump's power utility
        self.design_results['pump_efficiency'] = self.pump.design_results['Efficiency'] * 100.0


    def _cost(self):
        """
        Placeholder for Ion Exchange Column and Resin cost.
        Costs would be based on resin volume, column hardware.
        """
        # super()._cost()  # Ensure parent class cost logic runs and attributes are initialized
        if 'Resin Volume (L)' in self.design_results:
            resin_volume_L = self.design_results['Resin Volume (L)']
            # Example cost factors (these would ideally be class attributes or from a config)
            cost_per_L_resin = 1500 # $/L (highly variable based on resin type)
            column_hardware_factor = 0.6 # Hardware cost as a fraction of resin cost
            
            self.baseline_purchase_costs['IEX Resin'] = resin_volume_L * cost_per_L_resin
            self.baseline_purchase_costs['IEX Column Hardware'] = self.baseline_purchase_costs['IEX Resin'] * column_hardware_factor



class SprayDryer(bst.SprayDryer):pass
    # def _cost(self):
    #     super()._cost()  # Ensure parent class cost logic runs and attributes are initialized
    #     if 'Evaporation rate' in self.design_results:
    #         evaporation_rate = self.design_results['Evaporation rate']
    #         # Example cost factors (these would ideally be class attributes or from a config)
    #         # self.purchase_cost = 1000 * evaporation_rate ** 0.5  # Placeholder formula
    #         # self.power_utility(50 * evaporation_rate)  # Placeholder power calculation

