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
    V_max_default = 500
    def _init(
            self, 
            fermentation_reaction, 
            cell_growth_reaction, 
            respiration_reaction,
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
        self.fermentation_reaction = fermentation_reaction
        self.cell_growth_reaction = cell_growth_reaction
        self.respiration_reaction = respiration_reaction
    
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
        if self.cell_disruption_method is 'HighPressure':
            self.Parameters = Parameters if Parameters is not None else {
                'P High': 150e5,  # High pressure in Pascals
                'P Low': 101325,  # Low pressure in Pascals (atmospheric)
            }

        elif self.cell_disruption_method is 'BallMilling':
            self.Parameters = Parameters if Parameters is not None else {
                'Flow rate': 100e-3,  # Flow rate in m^3/s
                'Media': 'Water',
            }
        elif self.cell_disruption_method is 'Sonication':
            self.Parameters = Parameters if Parameters is not None else {
                'Frequency': 20e3,  # Sonication frequency in Hz
                'Amplitude': 100e-6,  # Sonication amplitude in meters
            }
        elif self.cell_disruption_method is 'Enzymatic':
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
        self.purchase_costs['High-Pressure Pump'] = self.pump.purchase_cost
        self.purchase_costs['High-Pressure Valve'] = self.valve.purchase_cost

        # --- Aggregate Utility Costs ---
        # The main cost is the electricity for the pump.
        self.power_utility.rate = self.pump.power_utility.rate

        # Add any heat utilities from the valve (e.g., Joule-Thomson cooling)
        for hu in self.valve.heat_utilities:
            self.heat_utilities.append(hu)
    
    
class ProteinCentrifuge(bst.SolidsCentrifuge): pass


class Evaporator(bst.MultiEffectEvaporator): pass


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
        [0] Permeate stream (water, salts, smaller molecules).
        [1] Retentate stream (concentrated target product).
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
    _N_ins = 2  # Feed and Wash Solution
    _N_outs = 2 # Permeate and Retentate

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
    _default_membrane_lifetime_years = 2.0
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

        # Assume isothermal operation and outlet pressures equal to feed pressure
        # Actual pressure drop for pumping is handled by TMP_bar in _design
        retentate.T = permeate.T = feed.T
        retentate.P = permeate.P = feed.P # This is outlet P, TMP is an internal driving force

        permeate.empty()
        retentate.empty()

        # --- Water Balance ---
        feed_water_mass = feed.imass[self.water_ID]
        wash_water_mass = wash_solution.imass[self.water_ID]
        total_incoming_water = feed_water_mass + wash_water_mass

        retentate_water_from_feed = feed_water_mass * (1.0 - self.FeedWater_Recovery_to_Permeate)
        retentate.imass[self.water_ID] = max(0.0, retentate_water_from_feed)
        
        permeate_water_total = total_incoming_water - retentate.imass[self.water_ID]
        permeate.imass[self.water_ID] = max(0.0, permeate_water_total)

        # Ensure overall water balance due to max(0,...) or floating point precision
        current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
        if abs(current_total_water_out - total_incoming_water) > 1e-9: # Tolerance
            if total_incoming_water >= retentate.imass[self.water_ID]:
                permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
            else: # This case should ideally not be reached with the logic above
                retentate.imass[self.water_ID] = total_incoming_water
                permeate.imass[self.water_ID] = 0.0

        # --- Solute Balance ---
        for chem in self.chemicals:
            ID = chem.ID
            if ID == self.water_ID:
                continue

            mass_in_feed = feed.imass[ID]
            mass_in_wash = wash_solution.imass[ID] # Solutes might be in wash solution
            total_mass_in = mass_in_feed + mass_in_wash

            if total_mass_in <= 1e-12: # Effectively zero
                retentate.imass[ID] = 0.0
                permeate.imass[ID] = 0.0
                continue

            current_retention = self.DefaultSolutes_Retention
            if self.TargetProduct_ID:
                if isinstance(self.TargetProduct_ID, str) and ID == self.TargetProduct_ID:
                    current_retention = self.TargetProduct_Retention
                if isinstance(self.TargetProduct_ID, list) and ID in self.TargetProduct_ID:
                    current_retention = self.TargetProduct_Retention
            elif self.Salt_ID:
                if isinstance(self.Salt_ID, str) and ID == self.Salt_ID:
                    current_retention = self.Salt_Retention
                elif isinstance(self.Salt_ID, list) and ID in self.Salt_ID: # Handles list of salt IDs
                    current_retention = self.Salt_Retention
            elif self.OtherLargeMolecules_ID:
                if isinstance(self.OtherLargeMolecules_ID, str) and ID == self.OtherLargeMolecules_ID:
                    current_retention = self.OtherLargeMolecules_Retention
                elif isinstance(self.OtherLargeMolecules_ID, list) and ID in self.OtherLargeMolecules_ID: # Handles list
                    current_retention = self.OtherLargeMolecules_Retention
            
            retentate_mass_solute = total_mass_in * current_retention
            retentate.imass[ID] = max(0.0, retentate_mass_solute)
            
            permeate_mass_solute = total_mass_in - retentate.imass[ID] # Permeate by difference for mass balance
            permeate.imass[ID] = max(0.0, permeate_mass_solute)
            
            # Final check for solute mass balance due to max(0,...) or floating point nuances
            current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
            mass_balance_error = current_total_solute_out - total_mass_in
            # Check relative and absolute error
            if abs(mass_balance_error) > (1e-9 * abs(total_mass_in) + 1e-12):
                permeate.imass[ID] -= mass_balance_error # Adjust permeate
                if permeate.imass[ID] < 0:
                    retentate.imass[ID] += permeate.imass[ID] # Add deficit to retentate
                    permeate.imass[ID] = 0.0
                    if retentate.imass[ID] < 0: # Should not happen if total_mass_in >=0
                        retentate.imass[ID] = 0.0
                        # Consider logging a warning here if mass is lost.
                        # print(f"Warning: Mass balance issue for {ID} in {self.ID}. Input: {total_mass_in:.2e}, Output: {current_total_solute_out - mass_balance_error:.2e}")

    def _design(self):
        Design = self.design_results
        # --- Membrane Area Calculation ---
        permeate_stream = self.outs[0]
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
        # Design['membrane_lifetime_years'] = self.membrane_lifetime_years
        # Design['module_cost_factor'] = self.module_cost_factor
        # Design['module_cost_exponent'] = self.module_cost_exponent
        # Design['base_CEPCI'] = self.base_CEPCI

        # --- Pump Power Calculation ---
        # Total volumetric flow to be pumped (feed + wash solution) in m^3/hr
        internal_stream = self.ins[0].copy()
        self.pump = bst.Pump(None, None, P=self.TMP_bar * 1e5)
        self.pump.ins[0] = internal_stream
        self.pump.simulate()
        self.pump._design()  # Design the pump to get its cost
        self.power_utility = self.pump.power_utility  # Use the pump's power utility
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
            self.purchase_costs['Membrane System'] = current_purchase_cost
        else:
            self.purchase_costs['Membrane System'] = 0.0

        # --- Annual Operating Cost (OPEX) for Membrane Replacement ---
        # This is added to `add_OPEX` for the TEA.
        if (self.membrane_lifetime_years > 0 and
            self.membrane_cost_USD_per_m2 > 0 and
            area_m2 > 0):
            annual_replacement_cost = (area_m2 * self.membrane_cost_USD_per_m2) / self.membrane_lifetime_years
            self.purchase_costs['Membrane replacement'] = annual_replacement_cost
        else:
            self.purchase_costs['Membrane replacement'] = 0.0
        #self.power_utility.cost = self.power_utility.rate * bst.annual_hours * bst.electricity_cost_per_kWh
        self.purchase_costs['Pump'] = self.pump.purchase_cost
        # Tank costs are not included here.

# class Diafiltration(bst.Unit):
#     _N_ins = 2  # Feed and Wash Solution
#     _N_outs = 2 # Permeate and Retentate

#     water_ID='H2O'
#     _default_TargetProduct_ID = 'Leghemoglobin'
#     _default_Salt_ID = 'Salts'
#     _default_OtherLargeMolecules_ID = 'OtherLargeMolecules'
#     _default_DefaultSolutes_ID = 'DefaultSolutes'
#     _default_TargetProduct_Retention = 0.98
#     _default_Salt_Retention = 0.05
#     _default_OtherLargeMolecules_Retention = 0.99
#     _default_DefaultSolutes_Retention = 0.05
#     _default_FeedWater_Recovery_to_Permeate = 0.75

#     def __init__(self, ID='', ins=None, outs=None, thermo=None,
#             # Membrane properties
#             TargetProduct_ID=None,
#             Salt_ID=None,
#             OtherLargeMolecules_ID=None,
#             DefaultSolutes_ID=None,
#             TargetProduct_Retention=None,
#             Salt_Retention=None,
#             OtherLargeMolecules_Retention=None,
#             DefaultSolutes_Retention=None,
#             FeedWater_Recovery_to_Permeate=None, **kwargs):
#         super().__init__(ID, ins,outs,thermo)

#         # set membrane properties
#         self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
#         self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
#         self.OtherLargeMolecules_ID = OtherLargeMolecules_ID if OtherLargeMolecules_ID is not None else self._default_OtherLargeMolecules_ID
#         self.DefaultSolutes_ID = DefaultSolutes_ID if DefaultSolutes_ID is not None else self._default_DefaultSolutes_ID
#         self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
#         self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
#         self.OtherLargeMolecules_Retention = OtherLargeMolecules_Retention if OtherLargeMolecules_Retention is not None else self._default_OtherLargeMolecules_Retention
#         self.DefaultSolutes_Retention = DefaultSolutes_Retention if DefaultSolutes_Retention is not None else self._default_DefaultSolutes_Retention
#         self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate

#     def _run(self):
#         feed, wash_solution = self.ins
#         retentate, permeate = self.outs

#         retentate.T = permeate.T = feed.T
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

#         current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
#         if abs(current_total_water_out - total_incoming_water) > 1e-9: # Tolerance for floating point
#             if total_incoming_water >= retentate.imass[self.water_ID]:
#                 permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
#             else:
#                 retentate.imass[self.water_ID] = total_incoming_water
#                 permeate.imass[self.water_ID] = 0.0

#         # --- Solute Balance ---
#         for chem in self.chemicals:
#             ID = chem.ID
#             if ID == self.water_ID:
#                 continue

#             mass_in_feed = feed.imass[ID]
#             mass_in_wash = wash_solution.imass[ID]
#             total_mass_in = mass_in_feed + mass_in_wash

#             if total_mass_in <= 1e-12: # Effectively zero mass
#                 retentate.imass[ID] = 0.0
#                 permeate.imass[ID] = 0.0
#                 continue

#             current_retention = self.DefaultSolutes_Retention

#             if ID == self.TargetProduct_ID and self.TargetProduct_ID:
#                 current_retention = self.TargetProduct_Retention
#             elif self.Salt_ID and ID in self.Salt_ID:
#                 current_retention = self.Salt_Retention
#             elif self.OtherLargeMolecules_ID and ID in self.OtherLargeMolecules_ID:
#                 current_retention = self.OtherLargeMolecules_Retention
            
#             retentate_mass_solute = total_mass_in * current_retention
#             permeate_mass_solute = total_mass_in * (1.0 - current_retention)
            
#             retentate.imass[ID] = max(0.0, retentate_mass_solute)
#             # Permeate by difference for better mass balance, after retentate is set
#             permeate.imass[ID] = max(0.0, total_mass_in - retentate.imass[ID])


#             # Final check to ensure solute mass balance due to max(0,...) or floating point nuances
#             current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
#             mass_balance_error = current_total_solute_out - total_mass_in
#             if abs(mass_balance_error) > 1e-9 * total_mass_in and abs(mass_balance_error) > 1e-12 :
#                 # If there's a significant discrepancy, adjust one stream (e.g. permeate)
#                 # This logic can be refined, but aims to conserve mass.
#                 permeate.imass[ID] -= mass_balance_error 
#                 if permeate.imass[ID] < 0: # If adjustment makes it negative
#                     retentate.imass[ID] += permeate.imass[ID] # Add the negative part to retentate
#                     permeate.imass[ID] = 0.0
#                     if retentate.imass[ID] < 0: # Should not happen if total_mass_in >=0
#                         retentate.imass[ID] = 0.0
#                         # At this point, mass is lost if total_mass_in was positive. This indicates an issue.
#                         # For robustness, one might assign all to one stream if total_mass_in > 0.
#                         # print(f"Warning: Mass balance issue for {ID} in {self.ID}")


#     def _design(self):
#         """Placeholder for design calculations (e.g., membrane area)."""
#         # Example: Calculate membrane area based on permeate flux
#         # self.membrane_flux_rate = 50 # L/(m^2*hr), would be a parameter
#         # permeate_vol_flow_m3hr = self.outs[0].F_vol / 1000 # Assuming F_vol is in L/hr
        
#         # if permeate_vol_flow_m3hr > 0 and hasattr(self, 'membrane_flux_rate') and self.membrane_flux_rate > 0:
#         #     # Flux often given in L/m2/h or GFD. Ensure units are consistent.
#         #     # Example: flux in L/m2/h, F_vol in L/hr
#         #     required_area_m2 = self.outs[0].F_vol / self.membrane_flux_rate
#         #     self.design_results['Membrane Area (m^2)'] = required_area_m2
#         pass

#     def _cost(self):
#         """Placeholder for cost calculations."""
#         # if 'Membrane Area (m^2)' in self.design_results:
#         #     area = self.design_results['Membrane Area (m^2)']
#         #     # Example: CEPCI = bst.CE # Chemical Engineering Plant Cost Index
#         #     # self.purchase_costs['Membrane Module'] = cost_correlation_function(area) * CEPCI / 500
#         pass


# class IonExchange(bst.Unit):
#     """
#     Ion Exchange Chromatography unit for purification of a target product.
#     The unit separates the target product from impurities based on differential
#     binding to an ion exchange resin, followed by elution.

#     Parameters
#     ----------
#     ins : Sequence[Stream]
#         [0] Feed stream (conditioned) containing the target product and impurities.
#         [1] Elution Buffer Profile stream, defining the matrix of the product stream
#             (e.g., buffer composition and volume for elution).
#     outs : Sequence[Stream]
#         [0] Product stream, containing the eluted target product in elution buffer.
#         [1] Waste stream, containing unbound components, impurities, and feed buffer.
#     water_ID : str, optional
#         Chemical ID for water. Defaults to 'H2O'.
#     TargetProduct_ID : str, optional
#         Chemical ID of the target product. Defaults to 'Leghemoglobin'.
#     TargetProduct_Yield : float, optional
#         Fraction of the target product from the feed that is recovered in the product stream.
#         Defaults to 0.95.
#     BoundImpurity_ID : tuple[str], optional
#         Tuple of chemical IDs for impurities that bind to the resin but are
#         separated from the target product (e.g., eluted in waste).
#         Defaults to ('HostCellProtein', 'DNA', 'Endotoxin').
#     BoundImpurity_Removal_Efficiency : float, optional
#         Fraction of specified bound impurities (from feed) that are removed to the waste stream.
#         Defaults to 0.99.
#     NonBinding_Solutes_Carryover_to_Product : float, optional
#         Fraction of other solutes in the feed (not target, not specified bound impurities,
#         not elution buffer defining components) that carry over to the product stream.
#         Defaults to 0.05.
#     ElutionBuffer_Defining_Component_ID : tuple[str], optional
#         Tuple of chemical IDs for key components defining the elution buffer matrix
#         (e.g., elution salts from `ins[1]`). If these are also in `ins[0]`, the feed's
#         portion is assumed to go to waste. Defaults to ('NaCl', 'KCl').
#     resin_DBC_g_L : float, optional
#         Dynamic Binding Capacity of the resin in grams of TargetProduct per Liter of resin.
#         Defaults to 50.0.
#     load_safety_factor : float, optional
#         Safety factor applied to DBC for resin volume calculation (e.g., 0.8 for 80% utilization).
#         Defaults to 0.8.
#     operating_pressure_bar : float, optional
#         Assumed operating pressure for the column system for pump power calculation (bar).
#         Defaults to 3.0.
#     pump_efficiency : float, optional
#         Overall efficiency of the pumps for the IEX system.
#         Defaults to 0.75.
#     resin_cost_USD_per_L : float, optional
#         Cost of the ion exchange resin in USD per Liter.
#         Defaults to 1500.0.
#     resin_lifetime_years : float, optional
#         Expected lifetime of the resin in years, for replacement cost calculation.
#         Defaults to 1.0.
#     column_hardware_cost_factor : float, optional
#         Factor for estimating column hardware purchase cost based on resin volume.
#         Cost_hardware = factor * (ResinVolume_L ** exponent).
#         Defaults to 20000.0.
#     column_hardware_cost_exponent : float, optional
#         Exponent for column hardware cost correlation. Defaults to 0.6.
#     base_CEPCI : float, optional
#         The Chemical Engineering Plant Cost Index (CEPCI) for which the
#         cost factors are valid. Defaults to 500.0.

#     """
#     _N_ins = 2  # Feed (conditioned) and Elution Buffer Profile
#     _N_outs = 2 # Product (in elution buffer) and Waste Stream

#     # --- Default Values for Parameters ---
#     _default_water_ID = 'H2O'
#     _default_TargetProduct_ID = 'Leghemoglobin'
#     _default_TargetProduct_Yield = 0.95
#     _default_BoundImpurity_ID = ('HostCellProtein', 'DNA', 'Endotoxin')
#     _default_BoundImpurity_Removal_Efficiency = 0.99
#     _default_NonBinding_Solutes_Carryover_to_Product = 0.05
#     _default_ElutionBuffer_Defining_Component_IDs_tuple = ('NaCl', 'KCl')

#     # Design and Costing Defaults
#     _default_resin_DBC_g_L = 50.0  # g target / L resin
#     _default_load_safety_factor = 0.8
#     _default_operating_pressure_bar = 3.0 # For pumping energy
#     _default_pump_efficiency = 0.75
#     _default_resin_cost_USD_per_L = 1500.0
#     _default_resin_lifetime_years = 1.0 # Can be highly variable
#     _default_column_hardware_cost_factor = 20000.0 # For Cost = A * V_R^B
#     _default_column_hardware_cost_exponent = 0.6
#     _default_base_CEPCI = 500.0


#     def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
#                 water_ID=None, TargetProduct_ID=None, TargetProduct_Yield=None,
#                 BoundImpurity_ID=None, BoundImpurity_Removal_Efficiency=None,
#                 NonBinding_Solutes_Carryover_to_Product=None,
#                 ElutionBuffer_Defining_Component_ID=None,
#                 # Design and Costing Parameters
#                 resin_DBC_g_L=None, load_safety_factor=None,
#                 operating_pressure_bar=None, pump_efficiency=None,
#                 resin_cost_USD_per_L=None, resin_lifetime_years=None,
#                 column_hardware_cost_factor=None, column_hardware_cost_exponent=None,
#                 base_CEPCI=None,
#                  **kwargs):
#         super().__init__(ID, ins, outs, thermo, **kwargs)

#         # Set operational parameters
#         self.water_ID = water_ID if water_ID is not None else self._default_water_ID
#         self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
#         self.TargetProduct_Yield = TargetProduct_Yield if TargetProduct_Yield is not None else self._default_TargetProduct_Yield
#         self.BoundImpurity_ID = BoundImpurity_ID if BoundImpurity_ID is not None else self._default_BoundImpurity_ID
#         self.BoundImpurity_Removal_Efficiency = BoundImpurity_Removal_Efficiency if BoundImpurity_Removal_Efficiency is not None else self._default_BoundImpurity_Removal_Efficiency
#         self.NonBinding_Solutes_Carryover_to_Product = NonBinding_Solutes_Carryover_to_Product if NonBinding_Solutes_Carryover_to_Product is not None else self._default_NonBinding_Solutes_Carryover_to_Product
#         self.ElutionBuffer_Defining_Component_ID = ElutionBuffer_Defining_Component_ID if ElutionBuffer_Defining_Component_ID is not None else self._default_ElutionBuffer_Defining_Component_IDs_tuple

#         # Set design and costing parameters
#         self.resin_DBC_g_L = resin_DBC_g_L if resin_DBC_g_L is not None else self._default_resin_DBC_g_L
#         self.load_safety_factor = load_safety_factor if load_safety_factor is not None else self._default_load_safety_factor
#         self.operating_pressure_bar = operating_pressure_bar if operating_pressure_bar is not None else self._default_operating_pressure_bar
#         self.pump_efficiency = pump_efficiency if pump_efficiency is not None else self._default_pump_efficiency
#         self.resin_cost_USD_per_L = resin_cost_USD_per_L if resin_cost_USD_per_L is not None else self._default_resin_cost_USD_per_L
#         self.resin_lifetime_years = resin_lifetime_years if resin_lifetime_years is not None else self._default_resin_lifetime_years
#         self.column_hardware_cost_factor = column_hardware_cost_factor if column_hardware_cost_factor is not None else self._default_column_hardware_cost_factor
#         self.column_hardware_cost_exponent = column_hardware_cost_exponent if column_hardware_cost_exponent is not None else self._default_column_hardware_cost_exponent
#         self.base_CEPCI = base_CEPCI if base_CEPCI is not None else self._default_base_CEPCI

#         self.power_utility = bst.PowerUtility()

#     def _run(self):
#         feed = self.ins[0]
#         elution_buffer_profile = self.ins[1]
#         product_stream = self.outs[0]
#         waste_stream = self.outs[1]

#         product_stream.T = elution_buffer_profile.T
#         waste_stream.T = feed.T

#         product_stream.mass = elution_buffer_profile.mass.copy()
#         waste_stream.empty()

#         for chem in self.chemicals:
#             ID = chem.ID
#             feed_mass_component = feed.imass[ID]

#             if feed_mass_component <= 1e-12: # Effectively zero in feed
#                 # Ensure waste stream has this component initialized if it's not already
#                 if ID not in waste_stream.imol: waste_stream.imass[ID] = 0.0
#                 continue

#             if ID == self.TargetProduct_ID:
#                 mass_to_product = feed_mass_component * self.TargetProduct_Yield
#                 product_stream.imass[ID] += mass_to_product
#                 waste_stream.imass[ID] = feed_mass_component - mass_to_product
            
#             elif ID == self.water_ID:
#                 waste_stream.imass[ID] = feed_mass_component # Feed water to waste
            
#             elif self.ElutionBuffer_Defining_Component_ID and ID in self.ElutionBuffer_Defining_Component_ID:
#                 waste_stream.imass[ID] = feed_mass_component # Feed buffer components to waste
            
#             elif self.BoundImpurity_ID and ID in self.BoundImpurity_ID:
#                 mass_removed_to_waste = feed_mass_component * self.BoundImpurity_Removal_Efficiency
#                 waste_stream.imass[ID] = mass_removed_to_waste
#                 product_stream.imass[ID] += (feed_mass_component - mass_removed_to_waste)
            
#             else: # Other solutes
#                 mass_to_product = feed_mass_component * self.NonBinding_Solutes_Carryover_to_Product
#                 product_stream.imass[ID] += mass_to_product
#                 waste_stream.imass[ID] = feed_mass_component - mass_to_product
        
#         # Ensure all chemical amounts are non-negative
#         for chem_obj in self.chemicals:
#             idx = chem_obj.ID
#             if product_stream.imass[idx] < 0: product_stream.imass[idx] = 0.0
#             if waste_stream.imass[idx] < 0: waste_stream.imass[idx] = 0.0

#     # def _design(self):
#     #     # --- Resin Volume Calculation ---
#     #     target_product_load_kg_hr = self.ins[0].imass[self.TargetProduct_ID]
        
#     #     if (self.resin_DBC_g_L > 0 and
#     #         self.load_safety_factor > 0 and
#     #         target_product_load_kg_hr > 0):
#     #         # Effective DBC in kg/L
#     #         effective_DBC_kg_L = (self.resin_DBC_g_L / 1000.0) * self.load_safety_factor
#     #         # Resin volume (L) needed to process the continuous flow of target product
#     #         resin_volume_L = target_product_load_kg_hr / effective_DBC_kg_L
#     #     else:
#     #         resin_volume_L = 0.0
#     #     self.design_results['Resin Volume (L)'] = resin_volume_L

#     #     # --- Pump Power Calculation ---
#     #     # Total volumetric flow through pumps (feed + elution buffer)
#     #     # Assuming densities are available and non-zero
#     #     vol_feed_m3_hr = self.ins[0].F_vol if not self.ins[0].isempty() else 0.0
#     #     vol_elution_m3_hr = self.ins[1].F_vol if not self.ins[1].isempty() else 0.0
#     #     total_vol_m3_hr = vol_feed_m3_hr + vol_elution_m3_hr

#     #     if total_vol_m3_hr > 0 and self.operating_pressure_bar > 0 and self.pump_efficiency > 0:
#     #         op_pressure_Pa = self.operating_pressure_bar * 1e5 # Convert bar to Pascals
#     #         total_vol_m3_s = total_vol_m3_hr / 3600.0
#     #         # Power (kW) = Volumetric flow (m^3/s) * Pressure (Pa) / efficiency / 1000 (W to kW)
#     #         power_kW = (total_vol_m3_s * op_pressure_Pa) / self.pump_efficiency / 1000.0
#     #     else:
#     #         power_kW = 0.0
#     #     self.power_utility.rate = power_kW

#     # def _cost(self):
#     #     resin_volume_L = self.design_results.get('Resin Volume (L)', 0.0)

#     #     # --- Resin Purchase Cost ---
#     #     if resin_volume_L > 0 and self.resin_cost_USD_per_L > 0:
#     #         base_resin_cost = resin_volume_L * self.resin_cost_USD_per_L
#     #         # Adjust cost from base_CEPCI to current BioSTEAM CEPCI
#     #         current_resin_cost = base_resin_cost * (bst.CE / self.base_CEPCI)
#     #         self.purchase_costs['IEX Resin'] = current_resin_cost
#     #     else:
#     #         self.purchase_costs['IEX Resin'] = 0.0

#     #     # --- Column Hardware Purchase Cost ---
#     #     if resin_volume_L > 0 and self.column_hardware_cost_factor > 0:
#     #         base_hardware_cost = self.column_hardware_cost_factor * (resin_volume_L ** self.column_hardware_cost_exponent)
#     #         current_hardware_cost = base_hardware_cost * (bst.CE / self.base_CEPCI)
#     #         self.purchase_costs['IEX Column Hardware'] = current_hardware_cost
#     #     else:
#     #         self.purchase_costs['IEX Column Hardware'] = 0.0

#     #     # --- Annual Operating Cost (OPEX) for Resin Replacement ---
#     #     if (resin_volume_L > 0 and
#     #         self.resin_cost_USD_per_L > 0 and
#     #         self.resin_lifetime_years > 0):
#     #         # Cost of one fill of resin (using current cost, assuming replacement cost similar to initial)
#     #         resin_fill_cost = self.purchase_costs.get('IEX Resin', 0.0) # Use already CEPCI-adjusted cost
#     #         annual_resin_replacement_cost = resin_fill_cost / self.resin_lifetime_years
#     #         self.add_OPEX['IEX Resin Replacement'] = annual_resin_replacement_cost
#     #     else:
#     #         self.add_OPEX['IEX Resin Replacement'] = 0.0

#     #     # Electricity cost is handled by self.power_utility through the TEA.
#     #     # Costs of buffers (ins[1]) are handled by their upstream sourcing in the flowsheet.



class IonExchange(bst.Unit): 

    
    _N_ins = 2  # Feed (conditioned) and Elution Buffer Profile
    _N_outs = 2 # Product (in elution buffer) and Waste Stream
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
        if 'Resin Volume (L)' in self.design_results:
            resin_volume_L = self.design_results['Resin Volume (L)']
            # Example cost factors (these would ideally be class attributes or from a config)
            cost_per_L_resin = 1500 # $/L (highly variable based on resin type)
            column_hardware_factor = 0.6 # Hardware cost as a fraction of resin cost
            
            self.purchase_costs['IEX Resin'] = resin_volume_L * cost_per_L_resin
            self.purchase_costs['IEX Column Hardware'] = self.purchase_costs['IEX Resin'] * column_hardware_factor


# class NanofiltrationDF(bst.Unit):
#     _N_ins = 2  # ins[0] is Feed, ins[1] is Diafiltration Water
#     _N_outs = 2 # outs[0] is Retentate (Product), outs[1] is Permeate (Waste)

#     # --- Default Values ---
#     _default_water_ID = 'H2O'
#     _default_TargetProduct_ID = 'Leghemoglobin', 
#     _default_TargetProduct_MembraneRetention = 0.995

#     # NEW: Defaults for Additives (mirrors TargetProduct's retention default)
#     _default_Additive_ID = ('Pichia_pastoris','TrehaloseDH','SodiumAscorbate') # Default to empty tuple, user can specify additives
#     _default_Additive_MembraneRetention = 0.2 # Same high retention as target product

#     _default_Salt_ID = ('NaCl', 'KCl') 
#     _default_Salt_MembraneRetention_NF = 0.10

#     _default_OtherSmallSolutes_ID = ('ResidualSugars', 'SmallPeptides')
#     _default_OtherSmallSolutes_MembraneRetention_NF = 0.15

#     _default_DefaultUnspecifiedSolute_MembraneRetention_NF = 0.10
#     _default_FinalRetentateWater_MassRatio_to_FeedWater = 0.20
    
#     _min_retentate_volume_L_for_N_calc = 1e-6

#     def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
#                 water_ID=None,
#                 TargetProduct_ID=None,
#                 TargetProduct_MembraneRetention=None,
#                 # NEW: Additive parameters
#                 Additive_ID=None, # Expects tuple of chemical IDs
#                 Additive_MembraneRetention=None,
#                 Salt_ID=None, 
#                 Salt_MembraneRetention_NF=None,
#                 OtherSmallSolutes_ID=None, 
#                 OtherSmallSolutes_MembraneRetention_NF=None,
#                 DefaultUnspecifiedSolute_MembraneRetention_NF=None,
#                 FinalRetentateWater_MassRatio_to_FeedWater=None,
#                  **kwargs):
#         super().__init__(ID, ins, outs, thermo, **kwargs)

#         # Assign parameters
#         self.water_ID = water_ID if water_ID is not None else self._default_water_ID
#         self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
#         self.TargetProduct_MembraneRetention = TargetProduct_MembraneRetention if TargetProduct_MembraneRetention is not None else self._default_TargetProduct_MembraneRetention

#         # NEW: Assign Additive parameters
#         self.Additive_ID = Additive_ID if Additive_ID is not None else self._default_Additive_ID
#         self.Additive_MembraneRetention = Additive_MembraneRetention if Additive_MembraneRetention is not None else self._default_Additive_MembraneRetention

#         self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
#         self.Salt_MembraneRetention_NF = Salt_MembraneRetention_NF if Salt_MembraneRetention_NF is not None else self._default_Salt_MembraneRetention_NF

#         self.OtherSmallSolutes_ID = OtherSmallSolutes_ID if OtherSmallSolutes_ID is not None else self._default_OtherSmallSolutes_ID
#         self.OtherSmallSolutes_MembraneRetention_NF = OtherSmallSolutes_MembraneRetention_NF if OtherSmallSolutes_MembraneRetention_NF is not None else self._default_OtherSmallSolutes_MembraneRetention_NF

#         self.DefaultUnspecifiedSolute_MembraneRetention_NF = DefaultUnspecifiedSolute_MembraneRetention_NF if DefaultUnspecifiedSolute_MembraneRetention_NF is not None else self._default_DefaultUnspecifiedSolute_MembraneRetention_NF
#         self.FinalRetentateWater_MassRatio_to_FeedWater = FinalRetentateWater_MassRatio_to_FeedWater if FinalRetentateWater_MassRatio_to_FeedWater is not None else self._default_FinalRetentateWater_MassRatio_to_FeedWater

#     def _run(self):
#         feed = self.ins[0]
#         diafiltration_water_input = self.ins[1] # Buffer stream, may contain solutes
#         retentate = self.outs[0]
#         permeate = self.outs[1]

#         retentate.T = permeate.T = feed.T
#         permeate.empty()
#         retentate.empty()

#         # --- Water Balance --- (Your existing water balance logic is fine)
#         feed_water_mass = feed.imass[self.water_ID]
#         df_water_mass = diafiltration_water_input.imass[self.water_ID]
#         total_incoming_water = feed_water_mass + df_water_mass

#         if feed_water_mass < 1e-9:
#             retentate.imass[self.water_ID] = 0.0
#             permeate.imass[self.water_ID] = df_water_mass
#         else:
#             retentate.imass[self.water_ID] = feed_water_mass * self.FinalRetentateWater_MassRatio_to_FeedWater
#             retentate.imass[self.water_ID] = max(0.0, retentate.imass[self.water_ID])
            
#             permeate_water_from_feed = feed_water_mass - retentate.imass[self.water_ID]
#             permeate.imass[self.water_ID] = permeate_water_from_feed + df_water_mass
#             permeate.imass[self.water_ID] = max(0.0, permeate.imass[self.water_ID])

#         current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
#         if abs(current_total_water_out - total_incoming_water) > 1e-9:
#             if total_incoming_water >= retentate.imass[self.water_ID]:
#                 permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
#             else:
#                 retentate.imass[self.water_ID] = total_incoming_water
#                 permeate.imass[self.water_ID] = 0.0
        
#         # --- Diafiltration Efficiency Calculation --- (Your existing logic is fine)
#         rho_water_kg_L_approx = 1.0 
#         V_wash_L = df_water_mass / rho_water_kg_L_approx
#         V_retentate_during_DF_L = retentate.imass[self.water_ID] / rho_water_kg_L_approx
        
#         N_diavolumes = 0.0 
#         # Corrected access to _min_retentate_volume_L_for_N_calc
#         if V_retentate_during_DF_L > self._min_retentate_volume_L_for_N_calc and V_wash_L > 0:
#             N_diavolumes = V_wash_L / V_retentate_during_DF_L
        
#         # --- Solute Balance ---
#         for chem in self.chemicals:
#             ID = chem.ID
            
#             total_initial_solute_mass = feed.imass[ID] + diafiltration_water_input.imass[ID]

#             if total_initial_solute_mass <= 1e-12:
#                 retentate.imass[ID] = 0.0
#                 permeate.imass[ID] = 0.0
#                 continue

#             if ID == self.water_ID:
#                 continue

#             current_membrane_retention = self.DefaultUnspecifiedSolute_MembraneRetention_NF
#             # Flag for components that are highly retained and not subject to diafiltration wash-out formula
#             is_highly_retained_component = False 

#             if self.TargetProduct_ID and ID in self.TargetProduct_ID:
#                 current_membrane_retention = self.TargetProduct_MembraneRetention
#                 is_highly_retained_component = True
#             # NEW: Check for Additives
#             elif self.Additive_ID and ID in self.Additive_ID:
#                 current_membrane_retention = self.Additive_MembraneRetention
#                 is_highly_retained_component = True
#             elif self.Salt_ID and ID in self.Salt_ID:
#                 current_membrane_retention = self.Salt_MembraneRetention_NF
#             elif self.OtherSmallSolutes_ID and ID in self.OtherSmallSolutes_ID:
#                 current_membrane_retention = self.OtherSmallSolutes_MembraneRetention_NF
            
#             retentate_mass = 0.0
#             if is_highly_retained_component: # For TargetProduct and Additives
#                 retentate_mass = total_initial_solute_mass * current_membrane_retention
#             elif N_diavolumes > 0: # Apply diafiltration effect for salts and other small (washable) solutes
#                 passage_coefficient = 1.0 - current_membrane_retention
#                 fraction_remaining = np.exp(-N_diavolumes * passage_coefficient)
#                 retentate_mass = total_initial_solute_mass * fraction_remaining
#             else: # No diafiltration (N_diavolumes = 0), simple NF concentration for washable solutes
#                 # Solute retention applies to the total amount of that solute present.
#                 retentate_mass = total_initial_solute_mass * current_membrane_retention
            
#             retentate.imass[ID] = max(0.0, retentate_mass)
#             permeate.imass[ID] = max(0.0, total_initial_solute_mass - retentate.imass[ID])
            
#             # Final check for solute mass balance
#             current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
#             mass_balance_error = current_total_solute_out - total_initial_solute_mass
            
#             relative_error_check = False
#             if total_initial_solute_mass > 1e-12:
#                 relative_error_check = abs(mass_balance_error) > (1e-9 * total_initial_solute_mass)

#             if relative_error_check or abs(mass_balance_error) > 1e-12 :
#                 permeate.imass[ID] -= mass_balance_error
#                 if permeate.imass[ID] < 0.0:
#                     retentate.imass[ID] += permeate.imass[ID] 
#                     permeate.imass[ID] = 0.0
#                     if retentate.imass[ID] < 0.0: retentate.imass[ID] = 0.0

#     # def _design(self):
#     #     """Placeholder for Nanofiltration system design."""
#     #     # Calculation of required membrane area based on total permeate flow rate
#     #     # and typical NF membrane flux (e.g., L/(m^2*h) or GFD).
#     #     # permeate_flow_L_hr = self.outs[1].F_vol # F_vol is in L/hr
#     #     # membrane_flux_Lm2h = 20 # Example, would be a parameter
#     #     # if membrane_flux_Lm2h > 0 and permeate_flow_L_hr > 0:
#     #     #     self.design_results['Membrane Area (m^2)'] = permeate_flow_L_hr / membrane_flux_Lm2h
#     #     pass

#     # def _cost(self):
#     #     """Placeholder for Nanofiltration system capital cost."""
#     #     # Based on membrane area, pumps, instrumentation.
#     #     # if 'Membrane Area (m^2)' in self.design_results:
#     #     #     area = self.design_results['Membrane Area (m^2)']
#     #     #     # Cost = A * area^B (example correlation)
#     #     #     # self.purchase_costs['NF System'] = ...
#     #     pass



class SprayDryer(bst.SprayDryer): pass



# class Nanofiltration(bst.Unit):
#     """
#     Nanofiltration unit operation.

#     Separates components based on specified retention coefficients and
#     achieves a defined water recovery to the permeate. This class can be
#     configured for various nanofiltration applications by adjusting the
#     retention coefficients for different chemical species or groups.

#     Parameters
#     ----------
#     ins : streams
#         [0] Feed solution
#         [1] Wash solution (optional, can be an empty stream or water if
#             diafiltration-like operation is desired, otherwise flow can be set to zero)
#     outs : streams
#         [0] Retentate
#         [1] Permeate
#     TargetProduct_ID : str, optional
#         ID of the primary target product to be retained.
#         Defaults to 'Leghemoglobin'.
#     Salt_ID : str or tuple[str], optional
#         ID(s) of salts. For nanofiltration, this might represent
#         divalent salts (highly retained) or monovalent salts (partially retained).
#         Defaults to 'Salts'.
#     OtherLargeMolecules_ID : str or tuple[str], optional
#         ID(s) of other large molecules to be retained.
#         Defaults to 'OtherLargeMolecules'.
#     DefaultSolutes_ID : str or tuple[str], optional
#         ID(s) for solutes not otherwise categorized.
#         This is not used in the current solute logic but is kept for consistency.
#         The retention for any solute not matching other categories will fall
#         under `DefaultSolutes_Retention`.
#     TargetProduct_Retention : float, optional
#         Retention coefficient for the target product (fraction retained).
#         Defaults to 0.98 (98%).
#     Salt_Retention : float, optional
#         Retention coefficient for salts.
#         For NF, typical divalent salt retention: 0.8-0.99.
#         Typical monovalent salt retention: 0.2-0.7.
#         Defaults to 0.05 (5%), which might be low for typical NF salts; adjust as needed.
#     OtherLargeMolecules_Retention : float, optional
#         Retention coefficient for other large molecules.
#         Defaults to 0.99 (99%).
#     DefaultSolutes_Retention : float, optional
#         Retention coefficient for any other solutes not specified.
#         Defaults to 0.05 (5%).
#     FeedWater_Recovery_to_Permeate : float, optional
#         Fraction of water in the feed stream that is recovered in the permeate.
#         Note: This is based on feed water only. Wash water partitioning is
#         derived from total water balance.
#         Defaults to 0.75 (75%).
#     membrane_flux : float, optional
#         Average membrane flux in L/m^2/hr (LMH). Used for sizing.
#         Defaults to 30 LMH.
#     membrane_cost : float, optional
#         Cost of membrane per square meter (e.g., USD/m^2).
#         Defaults to 200 USD/m^2.
#     operating_pressure_bar : float, optional
#         Operating transmembrane pressure in bar. Used for energy calculation.
#         Defaults to 15 bar.
#     pump_efficiency : float, optional
#         Efficiency of the pump supplying the pressure.
#         Defaults to 0.75.

#     """
#     _N_ins = 2  # Feed and optional Wash Solution
#     _N_outs = 2 # Permeate and Retentate

#     # Default chemical identifiers - user should customize these for their system
#     water_ID = 'H2O'
#     _default_TargetProduct_ID = 'Protein' # Generic placeholder
#     _default_Addictive_ID = ('TrehaloseDH','SodiumAscorbate') # Placeholder for specific addictive
#     _default_Salt_ID = ('NaCl', 'KCl') # Example monovalent salts
#     _default_DivalentSalt_ID = ('CaCl2', 'MgSO4') # Example divalent salts
#     _default_OrganicMolecule_ID = 'Glucose' # Example organic molecule
#     # Note: The original _default_DefaultSolutes_ID is not directly used in selection logic below,
#     # it's implicitly handled by DefaultSolutes_Retention.

#     # Default retention values - user MUST adjust these for their specific NF membrane and application
#     _default_TargetProduct_Retention = 0.995       # High retention for a target macromolecule
#     _default_Addictive_Retention = 0.90         # Retention for a specific addictive
#     _default_Salt_Retention = 0.40                # Moderate retention for generic/monovalent salts
#     _default_DivalentSalt_Retention = 0.95        # Higher retention for divalent salts
#     _default_OrganicMolecule_Retention = 0.90     # Retention for a specific organic molecule
#     _default_DefaultSolutes_Retention = 0.10      # Low retention for other small, uncharged solutes

#     _default_FeedWater_Recovery_to_Permeate = 0.80 # Typical for NF concentration

#     # Design and cost parameters
#     _default_membrane_flux = 30.0 # L/m^2/hr (LMH) - typical NF range 10-70 LMH
#     _default_membrane_cost = 200.0 # USD/m^2
#     _default_operating_pressure_bar = 15.0 # bar - typical NF range 5-30 bar
#     _default_pump_efficiency = 0.75

#     def __init__(self, ID='', ins=None, outs=None, thermo=None,
#                 # --- Membrane separation properties ---
#                 TargetProduct_ID=None,
#                 Addictive_ID=None, # Specific category for NF
#                 Salt_ID=None, # Can be used for monovalent salts
#                 DivalentSalt_ID=None, # Specific category for NF
#                 OrganicMolecule_ID=None, # Specific category for NF
#                 TargetProduct_Retention=None,
#                 Addictive_Retention=None,
#                 Salt_Retention=None,
#                 DivalentSalt_Retention=None,
#                 OrganicMolecule_Retention=None,
#                 DefaultSolutes_Retention=None,
#                 FeedWater_Recovery_to_Permeate=None,
#                 # --- Design and cost parameters ---
#                 membrane_flux=None,
#                 membrane_cost=None,
#                 operating_pressure_bar=None,
#                 pump_efficiency=None,
#                  **kwargs):
#         super().__init__(ID, ins, outs, thermo)

#         # --- Set solute/membrane interaction properties ---
#         self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
#         self.Addictive_ID = Addictive_ID if Addictive_ID is not None else self._default_Addictive_ID
#         self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
#         self.DivalentSalt_ID = DivalentSalt_ID if DivalentSalt_ID is not None else self._default_DivalentSalt_ID
#         self.OrganicMolecule_ID = OrganicMolecule_ID if OrganicMolecule_ID is not None else self._default_OrganicMolecule_ID
#         # Ensure Salt_ID, DivalentSalt_ID, OrganicMolecule_ID are tuples for consistent 'in' checking
#         if isinstance(self.Salt_ID, str): self.Salt_ID = (self.Salt_ID,)
#         if isinstance(self.DivalentSalt_ID, str): self.DivalentSalt_ID = (self.DivalentSalt_ID,)
#         if isinstance(self.OrganicMolecule_ID, str): self.OrganicMolecule_ID = (self.OrganicMolecule_ID,)


#         self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
#         self.Addictive_Retention = Addictive_Retention if Addictive_Retention is not None else self._default_Addictive_Retention
#         self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
#         self.DivalentSalt_Retention = DivalentSalt_Retention if DivalentSalt_Retention is not None else self._default_DivalentSalt_Retention
#         self.OrganicMolecule_Retention = OrganicMolecule_Retention if OrganicMolecule_Retention is not None else self._default_OrganicMolecule_Retention
#         self.DefaultSolutes_Retention = DefaultSolutes_Retention if DefaultSolutes_Retention is not None else self._default_DefaultSolutes_Retention
#         self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate

#         # --- Set design and cost parameters ---
#         self.membrane_flux = membrane_flux if membrane_flux is not None else self._default_membrane_flux
#         self.membrane_cost = membrane_cost if membrane_cost is not None else self._default_membrane_cost
#         self.operating_pressure_bar = operating_pressure_bar if operating_pressure_bar is not None else self._default_operating_pressure_bar
#         self.pump_efficiency = pump_efficiency if pump_efficiency is not None else self._default_pump_efficiency


#     def _run(self):
#         feed, wash_solution = self.ins
#         retentate, permeate = self.outs

#         retentate.T = permeate.T = feed.T # Assume isothermal operation
#         permeate.P = retentate.P = feed.P # Pressure drop across membrane handled by pump energy
#         permeate.phase = retentate.phase = feed.phase

#         permeate.empty()
#         retentate.empty()

#         # --- Water Balance ---
#         feed_water_mass = feed.imass[self.water_ID]
#         wash_water_mass = wash_solution.imass[self.water_ID] if wash_solution else 0.0
#         total_incoming_water = feed_water_mass + wash_water_mass

#         # Water recovery is based on FEED water going to permeate
#         permeate_water_from_feed = feed_water_mass * self.FeedWater_Recovery_to_Permeate
        
#         # All wash water is assumed to be available for partitioning based on overall hydraulics
#         # A simple assumption: wash water splits in the same ratio as feed water, OR
#         # more realistically, wash water contributes to the flux.
#         # For diafiltration mode, wash water primarily exits via permeate.
#         # Here, we'll ensure overall water balance.
#         # Retentate water from feed is what's not recovered
#         retentate_water_from_feed = feed_water_mass * (1.0 - self.FeedWater_Recovery_to_Permeate)

#         # Total water in retentate = unrecovered feed water
#         # If there's wash water, it needs to be distributed.
#         # A common model for continuous diafiltration might assume wash water helps maintain volume
#         # while salts are washed out. For a general NF, if wash is used, it also passes through.
#         # Let's assume wash water primarily goes to permeate if system is not water limited for flux
#         # and retentate primarily holds back unrecovered feed water.

#         retentate.imass[self.water_ID] = max(0.0, retentate_water_from_feed)
        
#         # Permeate water is the sum of recovered feed water and all wash water, adjusted for retentate.
#         permeate_water_total = total_incoming_water - retentate.imass[self.water_ID]
#         permeate.imass[self.water_ID] = max(0.0, permeate_water_total)

#         # Ensure water balance rigorously
#         current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
#         if abs(current_total_water_out - total_incoming_water) > 1e-9 * total_incoming_water: # Tolerance
#             # Adjust permeate to close balance, assuming retentate water is correctly set by recovery
#             if total_incoming_water >= retentate.imass[self.water_ID]:
#                 permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
#             else: # Should not happen if recovery < 1
#                 retentate.imass[self.water_ID] = total_incoming_water
#                 permeate.imass[self.water_ID] = 0.0


#         # --- Solute Balance ---
#         for chem in self.chemicals:
#             ID = chem.ID
#             if ID == self.water_ID:
#                 continue

#             mass_in_feed = feed.imass[ID] if feed else 0.0
#             mass_in_wash = wash_solution.imass[ID] if wash_solution else 0.0
#             total_mass_in = mass_in_feed + mass_in_wash

#             if total_mass_in <= 1e-12: # Effectively zero mass
#                 retentate.imass[ID] = 0.0
#                 permeate.imass[ID] = 0.0
#                 continue

#             # Determine retention based on solute type
#             current_retention = self.DefaultSolutes_Retention # Default assumption

#             if self.TargetProduct_ID and ID == self.TargetProduct_ID:
#                 current_retention = self.TargetProduct_Retention
#             elif self.Addictive_ID and ID in self.Addictive_ID:
#                 current_retention = self.Addictive_Retention
#             elif self.DivalentSalt_ID and ID in self.DivalentSalt_ID:
#                 current_retention = self.DivalentSalt_Retention
#             elif self.Salt_ID and ID in self.Salt_ID: # Monovalent or general salts
#                 current_retention = self.Salt_Retention
#             elif self.OrganicMolecule_ID and ID in self.OrganicMolecule_ID:
#                 current_retention = self.OrganicMolecule_Retention
            
#             retentate_mass_solute = total_mass_in * current_retention
#             # Permeate mass is by difference to ensure mass balance for the solute
#             permeate_mass_solute = total_mass_in - retentate_mass_solute
            
#             retentate.imass[ID] = max(0.0, retentate_mass_solute)
#             permeate.imass[ID] = max(0.0, permeate_mass_solute)
            
#             # Final check to ensure solute mass balance due to max(0,...) or floating point nuances
#             current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
#             mass_balance_error = current_total_solute_out - total_mass_in
            
#             if abs(mass_balance_error) > 1e-9 * total_mass_in and abs(mass_balance_error) > 1e-12 :
#                 # If there's a significant discrepancy, adjust permeate preferably
#                 permeate.imass[ID] -= mass_balance_error 
#                 if permeate.imass[ID] < 0:
#                     retentate.imass[ID] += permeate.imass[ID] # Add the negative part (error) to retentate
#                     permeate.imass[ID] = 0.0
#                     if retentate.imass[ID] < 0: # Should not happen if total_mass_in >=0
#                         # This implies an issue, potentially mass was lost.
#                         # Fallback: put all solute in retentate if permeate becomes negative
#                         # or if total_mass_in was positive and now retentate is negative.
#                         if total_mass_in > 0:
#                             retentate.imass[ID] = total_mass_in
#                             permeate.imass[ID] = 0.0
#                         else:
#                             retentate.imass[ID] = 0.0
#                         # print(f"Warning: Corrected mass balance issue for {ID} in {self.ID}. Check inputs/retentions.")


#     def _design(self):
#         # """
#         # Design calculations for the nanofiltration unit.
#         # Calculates membrane area and pumping power requirement.
#         # """
#         # self.design_results = {}
        
#         # # --- Membrane Area Calculation ---
#         # # Permeate flow rate in L/hr (F_vol is m3/hr by default in BioSTEAM)
#         # permeate_L_per_hr = self.outs[1].F_vol * 1000 

#         # if permeate_L_per_hr > 1e-6 and self.membrane_flux > 1e-6:
#         #     required_area_m2 = permeate_L_per_hr / self.membrane_flux
#         #     self.design_results['Membrane Area (m^2)'] = required_area_m2
#         # else:
#         #     self.design_results['Membrane Area (m^2)'] = 0.0

#         # # --- Pumping Power Calculation ---
#         # # Power = (Volumetric Flow Rate * Pressure) / Efficiency
#         # # Feed volumetric flow rate (m^3/s)
#         # # Note: Using total feed flow (feed + wash) to the membrane system
#         # total_feed_vol_m3_hr = self.ins[0].F_vol
#         # if self.ins[1] and self.ins[1].F_mass > 0: # if wash stream exists and has flow
#         #     total_feed_vol_m3_hr += self.ins[1].F_vol

#         # feed_vol_m3_s = total_feed_vol_m3_hr / 3600.0
        
#         # # Pressure in Pascals (1 bar = 100,000 Pa)
#         # pressure_Pa = self.operating_pressure_bar * 1e5

#         # if feed_vol_m3_s > 1e-9 and pressure_Pa > 0 and self.pump_efficiency > 1e-3:
#         #     power_kW = (feed_vol_m3_s * pressure_Pa) / self.pump_efficiency / 1000.0 # Power in kW
#         #     self.power_utility.rate = power_kW
#         # else:
#         #     self.power_utility.rate = 0.0
#         pass


#     def _cost(self):
#         """
#         # Cost calculations for the nanofiltration unit.
#         # Estimates capital cost based on membrane area.
#         # Operating costs include electricity for pumping (handled by power_utility).
#         # """
#         # if 'Membrane Area (m^2)' in self.design_results:
#         #     area = self.design_results['Membrane Area (m^2)']
#         #     if area > 0:
#         #         # Capital cost for membrane modules
#         #         # This is a simple linear relationship. More complex correlations can be used.
#         #         # CEPCI for scaling can be applied here if base year of cost data is known.
#         #         # Example: self.purchase_costs['Membrane Modules'] = self.membrane_cost * area * (bst.CE / 500)
#         #         self.purchase_costs['Membrane Modules'] = self.membrane_cost * area
#         #     else:
#         #         self.purchase_costs['Membrane Modules'] = 0.0
#         # else:
#         #     self.purchase_costs['Membrane Modules'] = 0.0
        
#         # Add other costs like housing, pumps, instrumentation if needed
#         # For example, a factor of the membrane module cost or a separate correlation.
#         # self.purchase_costs['System Auxiliaries'] = 0.5 * self.purchase_costs['Membrane Modules']

#         # Operating costs (electricity) are handled by power_utility defined in _design.
#         # Membrane replacement costs could be added as an operating expense based on lifetime.
#         pass

# class ScrewPress(bst.ScrewPress): pass
