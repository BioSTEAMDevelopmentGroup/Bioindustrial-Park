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
import numpy as np

from biosteam.units.decorators import cost, copy_algorithm
from biosteam.units.design_tools import CEPCI_by_year, cylinder_diameter_from_volume, cylinder_area
from biosteam import tank_factory
from thermosteam import MultiStream
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch
import biosteam as bst
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
class SeedTrain(bst.Unit):
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
    
    _setup = bst.Unit._setup
    
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
        
        # Ensure water doesn't go negative
        if effluent.imol['H2O'] < 0.: 
            effluent.imol['H2O'] = 0.
        
        # Check if we have sufficient reactants before running reactions
        glucose_available = effluent.imol['Glucose']
        if glucose_available < 1e-6:  # Very small amount of glucose
            return
            
        # Store initial state to prevent excessive consumption
        initial_glucose = effluent.imol['Glucose']
        initial_O2 = effluent.imol['O2']
        
        self.fermentation_reaction.force_reaction(effluent)
        self.cell_growth_reaction.force_reaction(effluent)
        self.respiration_reaction.force_reaction(effluent)
        self.neutralization_reaction.force_reaction(effluent)
        
        # # Ensure no negative flows after reactions
        # for chemical in effluent.chemicals:
        #     if effluent.imol[chemical.ID] < 0:
        #         effluent.imol[chemical.ID] = 0
                
        # # Ensure minimum oxygen consumption for aeration calculation
        # # If oxygen is being produced instead of consumed, set a minimum consumption
        # O2_change = effluent.imol['O2'] - initial_O2
        # if O2_change >= 0:  # Oxygen is being produced or unchanged
        #     # Set a minimum oxygen consumption to ensure proper aeration
        #     min_O2_consumption = initial_glucose * 0.1  # 10% of glucose as minimum O2 consumption
        #     effluent.imol['O2'] = max(0, initial_O2 - min_O2_consumption)

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


# https://www.alibaba.com/product-detail/1000L-2000L-SUS-High-Shear-Maquina_1601244807309.html?spm=a2700.galleryofferlist.topad_classic.d_title.358b13a0lSlmVy&priceId=fa9c440338e5471698929aa241dc50fa
# 1000L 7000 USD single price
class CellDisruption(bst.Unit):
    """
    Cell disruption unit using high-pressure homogenization.
    Converts biomass into its constituent components using a manual mass balance.
    Models both power consumption and thermal effects of pressure drop.
    """
    _N_ins = 1
    _N_outs = 1
    _graphics = bst.Pump._graphics
    
    _F_BM_default = {
        'High-Pressure Homogenizer': 3.5, 
    }

    def __init__(self, ID='', ins=None, outs=(), 
                 cell_disruption_efficiency=0.55, # 50~60% typical as soluble ,others are debris
                 component_fractions=None,
                 P_high=150e5, P_low=101325):
        super().__init__(ID, ins, outs)
        
        self.cell_disruption_efficiency = cell_disruption_efficiency
        self.P_high = P_high
        self.P_low = P_low
        
        if component_fractions is None:
            self.component_fractions = {
                'Mannoprotein': 0.40, 'Glucan': 0.50, 'OleicAcid': 0.06,
                'Chitin': 0.03, 'RNA': 0.01,
                # 'Mannoprotein': 0.40, 'Glucan': 0.50/6, 'OleicAcid': 0.06/18,
                # 'Chitin': 0.03/8, 'RNA': 0.01/4,
            }
        else:
            self.component_fractions = component_fractions
        
        total_fraction = sum(self.component_fractions.values())
        if not np.isclose(total_fraction, 1.0):
            raise ValueError(f"Component fractions must sum to 1.0, but they sum to {total_fraction}.")
            
    def _run(self):
        """Simulate the mass balance of the disruption."""
        feed = self.ins[0]
        outlet = self.outs[0]
        outlet.copy_like(feed)
        
        disrupted_mass = feed.imass['Pichia_pastoris'] * self.cell_disruption_efficiency
        
        outlet.imass['Pichia_pastoris'] -= disrupted_mass
        
        for component, fraction in self.component_fractions.items():
            outlet.imass[component] += fraction * disrupted_mass

        # The final temperature will be determined by the valve simulation in _design
        # For now, we set pressure and leave temperature as is.
        pumpout = bst.Stream('')
        pumpout.copy_like(outlet)
        pumpout.P = self.P_high
        temp_valve = bst.IsenthalpicValve(
            ID='valve',
            ins=pumpout,
            P=self.P_low, 
            vle=False,
        )
        temp_valve.simulate()
        outlet.copy_like(temp_valve.outs[0])

        # outlet.P= self.P_low

    
    def _design(self):
        """Design the homogenizer, calculating power and thermal effects."""
        feed = self.ins[0]
        
        # --- Create temporary, local unit operations for calculation ---
        # This is the key: they are not attached to `self`.
        
        # 1. Multi-stage pump system to achieve high pressure
        # Calculate number of stages needed (typical compression ratio ~3-5 per stage)
        compression_ratio_per_stage = 4.0
        n_stages = max(1, int(np.ceil(np.log(self.P_high / feed.P) / np.log(compression_ratio_per_stage))))
        
        # Create pressure stages
        pressure_stages = [feed.P * (compression_ratio_per_stage ** i) for i in range(1, n_stages + 1)]
        pressure_stages[-1] = self.P_high  # Ensure final pressure is exact
        
        # Create and simulate multi-stage pumps
        total_power = 0
        current_stream = feed.copy()
        
        for target_pressure in pressure_stages:
            temp_pump = bst.Pump(
            ins=current_stream.copy(),
            P=target_pressure
            )
            temp_pump.simulate()
            total_power += temp_pump.power_utility.rate
            current_stream = temp_pump.outs[0].copy()
        
        # Store total power consumption
        self.power_utility.rate = total_power
        
        # 2. Simulate the valve to find the thermal effect (e.g., cooling)
        temp_valve = bst.IsenthalpicValve(
            ID=f'_{self.ID}_valve',
            ins=temp_pump.outs[0],
            P=self.P_low, 
            vle=False,
        )
        temp_valve.simulate()

        # --- Capture the results for the main CellDisruption unit ---
        #self.power_utility.rate = temp_pump.power_utility.rate
        self.heat_utilities = temp_valve.heat_utilities
        
        # # Update the outlet stream's temperature to the final calculated temperature
        # self.outs[0].T = temp_valve.outs[0].T
        
        # # --- Costing based on the calculated power ---
        # power_kW = self.power_utility.rate
        # Inguva, P., Grasselli, S. & Heng, P. W. S. High pressure homogenization – 
        # An update on its usage and understanding. Chemical Engineering Research and Design 202, 284–302 (2024).
        dP = (self.P_high - self.P_low)/1e5
        power_kW = (feed.F_vol * 1000) * dP / (36000 * 0.85)  # Assuming 85% efficiency
        #self.power_utility.rate = power_kW

        # purchase_cost = 20000 * (power_kW / 10)**0.64 if power_kW > 0 else 0
        purchase_cost = 90000 * (feed.F_vol*1000)**0.5 * (dP/1000)**1.5

        self.design_results['Power (kW)'] = power_kW
        self.design_results['Purchase cost (USD)'] = purchase_cost
        
    def _cost(self):
        """Cost the homogenizer."""
        self.baseline_purchase_costs['High-Pressure Homogenizer'] = self.design_results['Purchase cost (USD)']


class ProteinCentrifuge(bst.SolidsCentrifuge): pass


class Diafiltration(bst.Unit):
    """
    Diafiltration unit for separation of solutes based on size, typically
    retaining larger molecules (like proteins) while allowing smaller ones
    (like salts and water) to pass through the permeate. Includes continuous
    addition of wash solution (diafiltration buffer).

    (Docstrings for parameters remain the same as your original code)
    """
    _N_ins = 2
    _N_outs = 2
    
    # --- All default values, _F_BM_default, and _units remain the same ---
    _F_BM_default = {
        'Membrane System': 1.65,
        'Membrane replacement': 1.65,#1.0,
        'Pump': 1.89,
    }
    water_ID = 'H2O'
    _default_TargetProduct_ID = 'TargetProduct'
    _default_Salt_ID = 'Salt'
    _default_OtherLargeMolecules_ID = 'OtherLargeMolecules'
    _default_TargetProduct_Retention = 0.98
    _default_Salt_Retention = 0.05
    _default_OtherLargeMolecules_Retention = 0.98
    _default_DefaultSolutes_Retention = 0.08
    _default_FeedWater_Recovery_to_Permeate = 0.75
    _default_membrane_flux_LMH = 40.0 # ultra: 20~80, nano: 10~40
    _default_TMP_bar1 = 2.0 # ultra: 2~4, nano: 10 ~25
    _default_TMP_bar2 = 2.0 # ultra: 1.5~3, nano: 3~6
    _default_membrane_cost_USD_per_m2 = 200.0 
    # food ultra: 50-250, nano: 70~300
    # for pharm: ultra: 1500~3000, nano: 2500~4500
    _default_membrane_lifetime_years = 1.0
    _default_module_cost_factor = 2500.0
    _default_module_cost_exponent = 0.7
    _default_base_CEPCI = 500.0
    _default_reciculation_ratio = 10.0  # Recirculation ratio to reduce fouling 5~20
    _default_equipment_lifetime_years = 20.0  # Typical equipment lifetime for CAPEX calculation
    _units = {
        'Membrane Area': 'm2',
        'TargetProduct_Retention': '%',
        'Salt_Retention': '%',
        'OtherLargeMolecules_Retention': '%',
        'DefaultSolutes_Retention': '%',
        'FeedWater_Recovery_to_Permeate': '%',
        'membrane_flux_LMH': 'LMH',
        'TMP_bar1': 'bar',
        'TMP_bar2': 'bar',
        'pump_efficiency': '%',
        'membrane_cost_USD_per_m2': '$/m2',
        'membrane_lifetime_years': 'years',
        'equipment_lifetime_years': 'years',
        'module_cost_factor': '$/m2^exponent',
        'module_cost_exponent': '0.6',
        'base_CEPCI': '500',
    }
    
    def __init__(self, ID='', ins=None, outs=None, thermo=None,
                 TargetProduct_ID=None, Salt_ID=None, OtherLargeMolecules_ID=None,
                 TargetProduct_Retention=None, Salt_Retention=None,
                 OtherLargeMolecules_Retention=None, DefaultSolutes_Retention=None,
                 FeedWater_Recovery_to_Permeate=None,
                 membrane_flux_LMH=None, TMP_bar1=None, TMP_bar2=None,
                 membrane_cost_USD_per_m2=None, membrane_lifetime_years=None,
                 equipment_lifetime_years=None,
                 module_cost_factor=None, module_cost_exponent=None, base_CEPCI=None,
                 reciculation_ratio=None,
                 **kwargs):
        super().__init__(ID, ins, outs, thermo)
        # --- All __init__ logic remains the same ---
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
        self.OtherLargeMolecules_ID = OtherLargeMolecules_ID if OtherLargeMolecules_ID is not None else self._default_OtherLargeMolecules_ID
        self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
        self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
        self.OtherLargeMolecules_Retention = OtherLargeMolecules_Retention if OtherLargeMolecules_Retention is not None else self._default_OtherLargeMolecules_Retention
        self.DefaultSolutes_Retention = DefaultSolutes_Retention if DefaultSolutes_Retention is not None else self._default_DefaultSolutes_Retention
        self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate
        self.membrane_flux_LMH = membrane_flux_LMH if membrane_flux_LMH is not None else self._default_membrane_flux_LMH
        self.TMP_bar1 = TMP_bar1 if TMP_bar1 is not None else self._default_TMP_bar1
        self.TMP_bar2 = TMP_bar2 if TMP_bar2 is not None else self._default_TMP_bar2
        self.membrane_cost_USD_per_m2 = membrane_cost_USD_per_m2 if membrane_cost_USD_per_m2 is not None else self._default_membrane_cost_USD_per_m2
        self.membrane_lifetime_years = membrane_lifetime_years if membrane_lifetime_years is not None else self._default_membrane_lifetime_years
        self.equipment_lifetime_years = equipment_lifetime_years if equipment_lifetime_years is not None else self._default_equipment_lifetime_years
        self.module_cost_factor = module_cost_factor if module_cost_factor is not None else self._default_module_cost_factor
        self.module_cost_exponent = module_cost_exponent if module_cost_exponent is not None else self._default_module_cost_exponent
        self.base_CEPCI = base_CEPCI if base_CEPCI is not None else self._default_base_CEPCI
        self.reciculation_ratio = reciculation_ratio if reciculation_ratio is not None else self._default_reciculation_ratio
        self.power_utility = bst.PowerUtility()

    def _run(self):
        feed, wash_solution = self.ins
        retentate, permeate = self.outs

        # Temperature and pressure logic remains the same
        retentate.T = permeate.T = feed.T
        retentate.P = permeate.P = feed.P

        permeate.empty()
        retentate.empty() # Start with empty streams for clarity

        # --- Water Balance (Unchanged) ---
        feed_water_mass = feed.imass[self.water_ID]
        wash_water_mass = wash_solution.imass[self.water_ID]
        total_incoming_water = feed_water_mass + wash_water_mass

        retentate_water_from_feed = feed_water_mass * (1.0 - self.FeedWater_Recovery_to_Permeate)
        retentate.imass[self.water_ID] = max(0.0, retentate_water_from_feed)
        
        permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
        if permeate.imass[self.water_ID] < 0:
            permeate.imass[self.water_ID] = 0.0
            retentate.imass[self.water_ID] = total_incoming_water

        # --- Solute Balance ---
        # Build a retention map from configured IDs. Accept str, list, tuple, or set.
        # Normalize inputs and validate against available chemical IDs.
        retention_map = {}

        # Helper to normalize IDs to a flat list
        def _to_id_list(x):
            if x is None:
                return []
            if isinstance(x, (list, tuple, set)):
                return [i for i in x if i is not None]
            return [x]

        target_ids = _to_id_list(self.TargetProduct_ID)
        large_mol_ids = _to_id_list(self.OtherLargeMolecules_ID)
        salt_ids = _to_id_list(self.Salt_ID)

        # Use a set of available chemical IDs for robust membership tests
        available_ids = {chem.ID for chem in self.chemicals}

        # Only add to retention_map if the chemical ID exists in thermo chemicals
        # If an ID is not found, skip and optionally warn once per call
        import warnings as _warnings

        def _map_ids(ids, retention, label):
            missing = []
            for chem_id in ids:
                if chem_id in available_ids:
                    retention_map[chem_id] = retention
                else:
                    missing.append(chem_id)
            if missing:
                _warnings.warn(
                    f"Diafiltration: {label} IDs not found in thermo chemicals and will use default retention: {missing}",
                    RuntimeWarning,
                )

        _map_ids(target_ids, self.TargetProduct_Retention, "TargetProduct")
        _map_ids(large_mol_ids, self.OtherLargeMolecules_Retention, "OtherLargeMolecules")
        _map_ids(salt_ids, self.Salt_Retention, "Salt")

        for chem in self.chemicals:
            ID = chem.ID
            if ID == self.water_ID:
                continue

            total_mass_in = feed.imass[ID] + wash_solution.imass[ID]
            if total_mass_in < 1e-12: continue

            current_retention = retention_map.get(ID, self.DefaultSolutes_Retention)

            # ####################################################################
            # POINT OF CHANGE: UNIFIED SOLUTE BALANCE LOGIC
            # ####################################################################
            # The previous if/else logic is removed. All solutes are now handled
            # by this single, more realistic model.
            
            # The mass of any solute in the retentate is its total input mass
            # multiplied by its specific retention factor.
            retentate.imass[ID] = total_mass_in * current_retention
            
            # The mass in the permeate is simply the remainder, ensuring mass balance.
            permeate.imass[ID] = total_mass_in - retentate.imass[ID]
            # ####################################################################
            
            # Final safety check remains the same
            if permeate.imass[ID] < 0:
                retentate.imass[ID] += permeate.imass[ID]
                permeate.imass[ID] = 0.0
            if retentate.imass[ID] < 0:
                permeate.imass[ID] += retentate.imass[ID]
                retentate.imass[ID] = 0.0

    def _design(self):
        # --- No changes to _design method ---
        Design = self.design_results
        permeate_stream = self.outs[1]
        if permeate_stream.isempty() or permeate_stream.rho == 0:
            permeate_vol_L_per_hr = 0.0
        else:
            permeate_vol_L_per_hr = permeate_stream.F_vol *1000 #(permeate_stream.F_mass / permeate_stream.rho) 
        if self.membrane_flux_LMH > 0 and permeate_vol_L_per_hr > 0:
            membrane_area_m2 = permeate_vol_L_per_hr / self.membrane_flux_LMH
        else:
            membrane_area_m2 = 0.0
        Design['Membrane Area'] = membrane_area_m2
        Design['membrane_flux_LMH'] = self.membrane_flux_LMH
        Design['TMP_bar1'] = self.TMP_bar1
        Design['TMP_bar2'] = self.TMP_bar2
        Design['membrane_cost_USD_per_m2'] = self.membrane_cost_USD_per_m2
        Design['membrane_lifetime_years'] = self.membrane_lifetime_years
        Design['equipment_lifetime_years'] = self.equipment_lifetime_years
        internal_stream = self.ins[0].copy() + self.ins[1].copy()
        self.pump1 = bst.Pump(None, None, P=self.TMP_bar1 * 1e5)
        self.pump1.ins[0] = internal_stream
        self.pump1.simulate()
        self.pump1._design()
        self.pump2 = bst.Pump(None, None, P=self.TMP_bar2 * 1e5)
        self.pump2.ins[0] = internal_stream * self.reciculation_ratio 
        self.pump2.simulate()
        self.pump2._design()
        Design['pump1_efficiency'] = 0.85 * 100.0
        Design['pump2_efficiency'] = 0.85 * 100.0
        self.power_utility.rate = self.pump1.power_utility.rate / (Design['pump1_efficiency']/100) + self.pump2.power_utility.rate / (Design['pump2_efficiency']/100)

    def _cost(self):
        # --- MODIFIED: Properly calculate total membrane replacement cost over equipment lifetime ---
        area_m2 = self.design_results.get('Membrane Area', 0.0)
        if area_m2 > 0 and self.module_cost_factor > 0 and self.base_CEPCI > 0:
            base_purchase_cost = self.module_cost_factor * (area_m2 ** self.module_cost_exponent)
            current_purchase_cost = base_purchase_cost * (bst.CE / self.base_CEPCI)
            self.baseline_purchase_costs['Membrane System'] = current_purchase_cost
        else:
            self.baseline_purchase_costs['Membrane System'] = 0.0
            
        # Calculate total membrane replacement cost over equipment lifetime
        if (self.membrane_lifetime_years > 0 and
            self.membrane_cost_USD_per_m2 > 0 and
            area_m2 > 0 and
            self.equipment_lifetime_years > 0):
            
            # Calculate number of membrane replacements over equipment lifetime
            num_replacements = self.equipment_lifetime_years / self.membrane_lifetime_years
            
            # Total membrane replacement cost over equipment lifetime
            total_replacement_cost = num_replacements * area_m2 * self.membrane_cost_USD_per_m2
            
            self.baseline_purchase_costs['Membrane replacement'] = total_replacement_cost
        else:
            self.baseline_purchase_costs['Membrane replacement'] = 0.0
            
        self.baseline_purchase_costs['Pump'] = self.pump1.purchase_cost + self.pump2.purchase_cost # Assume 2 pumps for operation feed & recirculation


class IonExchange(bst.Unit):
    """
    Simulates an Ion Exchange (IEX) chromatography column for protein purification.
    The unit models the binding of one or more target products and impurities 
    from a feed stream onto a resin, followed by elution into a separate elution buffer.
    
    (Docstrings for other parameters remain the same)
    """
    _N_ins = 2
    _N_outs = 2
    
    _F_BM_default = {
        'IEX Column': 2.5,
        'IEX Resin replacement': 2.5, #1.0,
        'Pump': 1.89,
    }

    _units = {
        'Resin Volume': 'L',
        'TargetProduct_Yield': '%',
        'BoundImpurity_Removal': '%',
        'NonBinding_Carryover': '%',
        'resin_DBC_g_L': 'g/L',
        'resin_cost_USD_per_L': 'USD/L',
        'resin_lifetime_years': 'years',
        'equipment_lifetime_years': 'years',
    }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 # MODIFIED: Now accepts a string, list, or tuple of IDs
                 TargetProduct_IDs=('Leghemoglobin'),
                 TargetProduct_Yield=0.95,
                 BoundImpurity_IDs=('Heme_b','Chitin','Globin','Mannoprotein','Glucan'),
                 BoundImpurity_Removal=0.97,
                 NonBinding_Carryover=0.04,
                 resin_DBC_g_L=50.0, # 30~120
                 load_safety_factor=0.8,
                 resin_cost_USD_per_L=7,#1500.0, 
                    # strong acid cation resin: 2~5
                    # weak acid cation resin: 3 ~ 6
                    # strong base anion resin: 7~15
                    # weak base anion resin: 5~10
                 # 4000~10000
                 resin_lifetime_years=1.0,
                 equipment_lifetime_years=20.0,
                 column_hardware_cost_factor=3000.0,
                 column_hardware_cost_exponent=0.6,
                 base_CEPCI=500.0,
                 **kwargs):
        super().__init__(ID, ins, outs, thermo, **kwargs)

        # MODIFIED: Ensure TargetProduct_IDs is always a tuple for consistency
        if isinstance(TargetProduct_IDs, str):
            self.TargetProduct_IDs = (TargetProduct_IDs,)
        else:
            self.TargetProduct_IDs = tuple(TargetProduct_IDs)
            
        self.TargetProduct_Yield = TargetProduct_Yield
        self.BoundImpurity_IDs = BoundImpurity_IDs
        self.BoundImpurity_Removal = BoundImpurity_Removal
        self.NonBinding_Carryover = NonBinding_Carryover
        self.resin_DBC_g_L = resin_DBC_g_L
        self.load_safety_factor = load_safety_factor
        self.resin_cost_USD_per_L = resin_cost_USD_per_L
        self.resin_lifetime_years = resin_lifetime_years
        self.equipment_lifetime_years = equipment_lifetime_years
        self.column_hardware_cost_factor = column_hardware_cost_factor
        self.column_hardware_cost_exponent = column_hardware_cost_exponent
        self.base_CEPCI = base_CEPCI
        self.power_utility = bst.PowerUtility()

    def _run(self):
        feed, elution_buffer = self.ins
        product, waste = self.outs

        product.copy_like(elution_buffer)
        waste.copy_like(feed)
        
        # MODIFIED: Streamlined logic into a single loop
        for chem_ID in feed.chemicals.IDs:
            if feed.imass[chem_ID] < 1e-12:
                continue
            
            solute_in_feed = feed.imass[chem_ID]

            if chem_ID in self.TargetProduct_IDs:
                # Handle target products
                solute_to_product = solute_in_feed * self.TargetProduct_Yield
                product.imass[chem_ID] += solute_to_product
                waste.imass[chem_ID] -= solute_to_product
            elif chem_ID in self.BoundImpurity_IDs:
                # Handle bound impurities
                solute_to_product = solute_in_feed * (1.0 - self.BoundImpurity_Removal)
                product.imass[chem_ID] += solute_to_product
                waste.imass[chem_ID] -= solute_to_product
            else:
                # Handle all other non-binding solutes
                solute_to_product = solute_in_feed * self.NonBinding_Carryover
                product.imass[chem_ID] += solute_to_product
                waste.imass[chem_ID] -= solute_to_product
        
        product.T = elution_buffer.T
        waste.T = feed.T

    def _design(self):
        Design = self.design_results
        
        # MODIFIED: Sum the mass of ALL target products for sizing
        total_target_mass_kg_hr = sum(self.ins[0].imass[ID] for ID in self.TargetProduct_IDs)
        target_mass_g_hr = total_target_mass_kg_hr * 1000.0
        
        if self.resin_DBC_g_L > 0 and self.load_safety_factor > 0:
            effective_DBC = self.resin_DBC_g_L * self.load_safety_factor
            resin_volume_L = target_mass_g_hr / effective_DBC if effective_DBC > 0 else 0.0
        else:
            resin_volume_L = 0.0
        
        Design['Resin Volume'] = resin_volume_L
        Design['TargetProduct_Yield'] = self.TargetProduct_Yield * 100
        Design['BoundImpurity_Removal'] = self.BoundImpurity_Removal * 100
        Design['NonBinding_Carryover'] = self.NonBinding_Carryover * 100
        Design['resin_lifetime_years'] = self.resin_lifetime_years
        Design['equipment_lifetime_years'] = self.equipment_lifetime_years

        # Pump logic remains the same
        internal_stream = self.ins[0].copy() + self.ins[1].copy() # * 15 # consider 15 CV times the elution buffer for cleaning
        self.pump = bst.Pump(None, P= 4.0 * 1e5) # assume 4 bar pressure drop
        self.pump.ins[0] = internal_stream
        self.pump.simulate()
        self.pump._design()
        self.power_utility.rate = self.pump.power_utility.rate

    def _cost(self):
        # MODIFIED: Calculate total resin replacement cost over equipment lifetime
        Costs = self.baseline_purchase_costs
        Design = self.design_results
        resin_volume_L = Design.get('Resin Volume', 0.0)
        
        if resin_volume_L > 0:
            base_cost = self.column_hardware_cost_factor * (resin_volume_L ** self.column_hardware_cost_exponent)
            Costs['IEX Column'] = base_cost * (bst.CE / self.base_CEPCI)
        else:
            Costs['IEX Column'] = 0.0

        # Calculate total resin replacement cost over equipment lifetime
        if (self.resin_lifetime_years > 0 and 
            self.resin_cost_USD_per_L > 0 and 
            resin_volume_L > 0 and 
            self.equipment_lifetime_years > 0):
            
            # Calculate number of resin replacements over equipment lifetime
            num_replacements = self.equipment_lifetime_years / self.resin_lifetime_years
            
            # Total resin replacement cost over equipment lifetime
            total_replacement_cost = num_replacements * resin_volume_L * self.resin_cost_USD_per_L
            
            Costs['IEX Resin replacement'] = total_replacement_cost
        else:
            Costs['IEX Resin replacement'] = 0.0
            
        Costs['Pump'] = self.pump.purchase_cost


class SprayDryer(bst.SprayDryer):pass

