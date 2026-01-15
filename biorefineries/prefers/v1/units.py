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
from matplotlib.pyplot import cool
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
                 Cell_ID='Pichia_pastoris', 
                 cell_disruption_efficiency=0.55, # 50~60% typical as soluble ,others are debris
                 component_fractions=None,
                 P_high=150e5, P_low=101325):
        super().__init__(ID, ins, outs)
        
        self.cell_disruption_efficiency = cell_disruption_efficiency
        self.P_high = P_high
        self.P_low = P_low
        self.Cell_ID= Cell_ID
        
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
        
        disrupted_mass = feed.imass[self.Cell_ID] * self.cell_disruption_efficiency
        
        outlet.imass[self.Cell_ID] -= disrupted_mass
        
        for component, fraction in self.component_fractions.items():
            outlet.imass[component] += fraction * disrupted_mass


        # Convert intracellular molecules to extracellular form based on disruption efficiency
        for chem in self.chemicals:
            chem_id = chem.ID
            if chem_id.endswith('_In'):
                # Get the extracellular molecule name by removing '_In' suffix
                extracellular_id = chem_id[:-3]
                if extracellular_id in self.chemicals:
                    # Calculate mass released from disrupted cells
                    intracellular_mass = feed.imass[chem_id]
                    released_mass = intracellular_mass * self.cell_disruption_efficiency
                    
                    # Transfer released mass to extracellular form
                    outlet.imass[extracellular_id] += released_mass
                    # Remaining intracellular mass stays in undisrupted cells
                    outlet.imass[chem_id] = intracellular_mass - released_mass
                    
        
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


class Centrifuge(bst.SolidsCentrifuge):pass

class ReverseOsmosis(bst.Unit):
    """
    Create a reverse osmosis unit operation for recovering water from brine.
    The model is based on a fraction of water recovered.
    
    Parameters
    ----------
    ins : 
        Inlet fluid to be split.
    outs : 
        * [0] Filtered water
        * [1] Brine
    water_recovery : float, optional
        Water recovered to 0th stream. Defaults to 0.987
    
    """
    _N_ins = 1
    _N_outs = 2
    
    _F_BM_default = {
        'Pump': 2.3,
        'Membrane replacement': 1.65,
        'Hardware': 2.0,
    }

    _units = {
        'Flow rate': 'm3/hr',
        'Power': 'kW',
        'Pump Pressure': 'bar',
        'Water recovery': '%',
        'Membrane Area': 'm2',
        'Membrane flux': 'L/m2/hr',
        'Membrane lifetime': 'years',
        'Plant lifetime': 'years',
        'Membrane cost per m2': 'USD/m2',
    }

    def _init(self, water_recovery=0.987, membrane_flux=40, membrane_cost_per_m2=300, 
                membrane_lifetime_years=3, plant_lifetime_years=20, operating_pressure_bar=30):
        self.water_recovery = water_recovery
        self.membrane_flux = membrane_flux  # L/m2/hr
        self.membrane_cost_per_m2 = membrane_cost_per_m2  # USD/m2
        self.membrane_lifetime_years = membrane_lifetime_years
        self.plant_lifetime_years = plant_lifetime_years
        self.operating_pressure_bar = operating_pressure_bar
    
    @property
    def RO_treated_water(self):
        return self.outs[0]
    
    def _run(self):
        feed, = self.ins
        water, brine = self.outs
        water.copy_thermal_condition(feed)
        brine.copy_like(feed)
        water_index = self.chemicals.index('H2O')
        water_flow = brine.mol[water_index]
        water_recovered = self.water_recovery * water_flow
        water.mol[water_index] = water_recovered
        brine.mol[water_index] = water_flow - water_recovered

    def _design(self):
        Design = self.design_results
        feed = bst.Stream()
        feed.copy_like(self.ins[0])
        
        # Record basic operating parameters
        Design['Flow rate'] = feed.F_vol  # m3/hr
        Design['Water recovery'] = self.water_recovery * 100  # %
        Design['Pump Pressure'] = self.operating_pressure_bar  # bar
        Design['Membrane flux'] = self.membrane_flux  # L/m2/hr
        Design['Membrane lifetime'] = self.membrane_lifetime_years  # years
        Design['Plant lifetime'] = self.plant_lifetime_years  # years
        Design['Membrane cost per m2'] = self.membrane_cost_per_m2  # USD/m2

        # Calculate membrane area based on flux and flow rate
        flow_rate_L_hr = feed.F_vol * 1000  # convert m3/hr to L/hr
        #flow_rate_L_hr = feed.ivol['H2O'] * 1000  # convert m3/hr to L/hr
        if self.membrane_flux > 0:
            membrane_area_m2 = flow_rate_L_hr / self.membrane_flux
        else:
            membrane_area_m2 = 0.0
        Design['Membrane Area'] = membrane_area_m2  # m2

        # Create temporary pump for power calculation
        self.pump = bst.Pump(
            ins=feed,
            outs=(),
            P=self.operating_pressure_bar * 1e5  # Convert bar to Pa
        )
        self.pump.simulate()
        self.pump._design()
        
        # Record power consumption
        power_kW = self.pump.power_utility.rate
        Design['Power'] = power_kW  # kW
        self.power_utility.rate = power_kW

    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        
        # Cost the pump
        if hasattr(self, 'pump') and self.pump:
            self.pump._cost()
            C['Pump'] = self.pump.purchase_cost
        else:
            C['Pump'] = 0.0
        
        # Cost membrane replacement based on membrane area
        membrane_area = D.get('Membrane Area', 0.0)
        if membrane_area > 0 and self.membrane_lifetime_years > 0:
            # Calculate total membrane replacement cost over plant lifetime
            num_replacements = self.plant_lifetime_years / self.membrane_lifetime_years
            total_membrane_cost = membrane_area * self.membrane_cost_per_m2 * num_replacements
            C['Membrane replacement'] = total_membrane_cost
        else:
            C['Membrane replacement'] = 0.0
        
        # Cost other hardware based on flow rate
        feed_flow_rate = D.get('Flow rate', 0.0)  # m3/hr


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
    _default_membrane_cost_USD_per_m2 = 150.0
    # Industry	Membrane Type	Approximate Value Range (USD/m²)
    # Food & Beverage	Ultrafiltration	$100 - $400
    #                   Nanofiltration	$150 - $500
    # Pharmaceutical	Ultrafiltration	$1,200 - $3,500+
    #                   Nanofiltration	$2,000 - $5,000+
    # https://www.biopharminternational.com/view/economic-analysis-single-use-tangential-flow-filtration-biopharmaceutical-applications


    # food ultra: 50-250, nano: 70~300
    # for pharm: ultra: 1500~3000, nano: 2500~4500
    # Membrane Type	Industrial / Water	Food & Beverage	Biopharmaceutical (TFF Cassettes)
    # Ultrafiltration (UF)	$80 - $350	$250 - $900	$1,500 - $5,000+
    # Nanofiltration (NF)	$100 - $400	$300 - $1,200	$2,000 - $7,000+
    _default_membrane_lifetime_years = 2.0 #1~3
    _default_module_cost_factor = 25000.0
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
        
        _mixed_stream = bst.Stream()
        _mixed_stream.mol = feed.mol + wash_solution.mol
        _mixed_stream.H = feed.H + wash_solution.H
        
        # Assign the correct temperature and pressure
        retentate.T = permeate.T = _mixed_stream.T
        retentate.P = permeate.P = _mixed_stream.P
        # retentate.T = permeate.T = feed.T
        # retentate.P = permeate.P = feed.P

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

class Ultrafiltration(bst.Unit):
    """
    Single-pass ultrafiltration unit with one feed stream.
    Retains the Diafiltration retention logic but removes diafiltration buffer
    and circulation pump handling.
    """
    _N_ins = 1
    _N_outs = 2

    _F_BM_default = {
        'Membrane System': 1.65,
        'Membrane replacement': 1.65,
        'Pump': 1.89,
    }
    water_ID = Diafiltration.water_ID
    _default_TargetProduct_ID = Diafiltration._default_TargetProduct_ID
    _default_Salt_ID = Diafiltration._default_Salt_ID
    _default_OtherLargeMolecules_ID = Diafiltration._default_OtherLargeMolecules_ID
    _default_TargetProduct_Retention = Diafiltration._default_TargetProduct_Retention
    _default_Salt_Retention = Diafiltration._default_Salt_Retention
    _default_OtherLargeMolecules_Retention = Diafiltration._default_OtherLargeMolecules_Retention
    _default_DefaultSolutes_Retention = Diafiltration._default_DefaultSolutes_Retention
    _default_FeedWater_Recovery_to_Permeate = Diafiltration._default_FeedWater_Recovery_to_Permeate
    _default_membrane_flux_LMH = Diafiltration._default_membrane_flux_LMH
    _default_TMP_bar = Diafiltration._default_TMP_bar1
    _default_membrane_cost_USD_per_m2 = Diafiltration._default_membrane_cost_USD_per_m2
    _default_membrane_lifetime_years = Diafiltration._default_membrane_lifetime_years
    _default_module_cost_factor = Diafiltration._default_module_cost_factor
    _default_module_cost_exponent = Diafiltration._default_module_cost_exponent
    _default_base_CEPCI = Diafiltration._default_base_CEPCI
    _default_equipment_lifetime_years = Diafiltration._default_equipment_lifetime_years

    _units = Diafiltration._units

    def __init__(self, ID='', ins=None, outs=None, thermo=None,
                 TargetProduct_ID=None, Salt_ID=None, OtherLargeMolecules_ID=None,
                 TargetProduct_Retention=None, Salt_Retention=None,
                 OtherLargeMolecules_Retention=None, DefaultSolutes_Retention=None,
                 FeedWater_Recovery_to_Permeate=None,
                 membrane_flux_LMH=None, TMP_bar=None,
                 membrane_cost_USD_per_m2=None, membrane_lifetime_years=None,
                 equipment_lifetime_years=None,
                 module_cost_factor=None, module_cost_exponent=None, base_CEPCI=None,
                 **kwargs):
        super().__init__(ID, ins, outs, thermo, **kwargs)
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
        self.OtherLargeMolecules_ID = OtherLargeMolecules_ID if OtherLargeMolecules_ID is not None else self._default_OtherLargeMolecules_ID
        self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
        #self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
        self.OtherLargeMolecules_Retention = OtherLargeMolecules_Retention if OtherLargeMolecules_Retention is not None else self._default_OtherLargeMolecules_Retention
        self.DefaultSolutes_Retention = DefaultSolutes_Retention if DefaultSolutes_Retention is not None else self._default_DefaultSolutes_Retention
        self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate
        self.Salt_Retention = 1 - self.FeedWater_Recovery_to_Permeate
        self.membrane_flux_LMH = membrane_flux_LMH if membrane_flux_LMH is not None else self._default_membrane_flux_LMH
        self.TMP_bar = TMP_bar if TMP_bar is not None else self._default_TMP_bar
        self.membrane_cost_USD_per_m2 = membrane_cost_USD_per_m2 if membrane_cost_USD_per_m2 is not None else self._default_membrane_cost_USD_per_m2
        self.membrane_lifetime_years = membrane_lifetime_years if membrane_lifetime_years is not None else self._default_membrane_lifetime_years
        self.equipment_lifetime_years = equipment_lifetime_years if equipment_lifetime_years is not None else self._default_equipment_lifetime_years
        self.module_cost_factor = module_cost_factor if module_cost_factor is not None else self._default_module_cost_factor
        self.module_cost_exponent = module_cost_exponent if module_cost_exponent is not None else self._default_module_cost_exponent
        self.base_CEPCI = base_CEPCI if base_CEPCI is not None else self._default_base_CEPCI
        self.power_utility = bst.PowerUtility()

    def _run(self):
        feed, = self.ins
        retentate, permeate = self.outs

        retentate.T = permeate.T = feed.T
        retentate.P = permeate.P = feed.P

        permeate.empty()
        retentate.empty()

        feed_water_mass = feed.imass[self.water_ID]
        retentate_water = feed_water_mass * (1.0 - self.FeedWater_Recovery_to_Permeate)
        retentate.imass[self.water_ID] = max(0.0, retentate_water)
        permeate.imass[self.water_ID] = max(0.0, feed_water_mass - retentate.imass[self.water_ID])

        import warnings as _warnings
        available_ids = {chem.ID for chem in self.chemicals}

        def _to_id_list(x):
            if x is None:
                return []
            if isinstance(x, (list, tuple, set)):
                return [i for i in x if i is not None]
            return [x]

        retention_map = {}
        attr_map = (
            (_to_id_list(self.TargetProduct_ID), self.TargetProduct_Retention, "TargetProduct"),
            (_to_id_list(self.OtherLargeMolecules_ID), self.OtherLargeMolecules_Retention, "OtherLargeMolecules"),
            (_to_id_list(self.Salt_ID), self.Salt_Retention, "Salt"),
        )
        for ids, retention, label in attr_map:
            missing = []
            for chem_id in ids:
                if chem_id in available_ids:
                    retention_map[chem_id] = retention
                else:
                    missing.append(chem_id)
            if missing:
                _warnings.warn(
                    f"Ultrafiltration: {label} IDs not found in thermo chemicals and will use default retention: {missing}",
                    RuntimeWarning,
                )

        for chem in self.chemicals:
            ID = chem.ID
            if ID == self.water_ID:
                continue
            mass_in = feed.imass[ID]
            if mass_in < 1e-12:
                continue

            current_retention = retention_map.get(ID, self.DefaultSolutes_Retention)
            retentate_mass = mass_in * current_retention
            permeate_mass = mass_in - retentate_mass

            retentate.imass[ID] = max(0.0, retentate_mass)
            permeate.imass[ID] = max(0.0, permeate_mass)

            if permeate.imass[ID] < 0:
                retentate.imass[ID] += permeate.imass[ID]
                permeate.imass[ID] = 0.0
            if retentate.imass[ID] < 0:
                permeate.imass[ID] += retentate.imass[ID]
                retentate.imass[ID] = 0.0

    def _design(self):
        Design = self.design_results
        permeate_stream = self.outs[1]
        if permeate_stream.isempty() or permeate_stream.rho == 0:
            permeate_vol_L_per_hr = 0.0
        else:
            permeate_vol_L_per_hr = permeate_stream.F_vol * 1000.0

        membrane_area_m2 = permeate_vol_L_per_hr / self.membrane_flux_LMH if (self.membrane_flux_LMH > 0 and permeate_vol_L_per_hr > 0) else 0.0
        Design['Membrane Area'] = membrane_area_m2
        Design['membrane_flux_LMH'] = self.membrane_flux_LMH
        Design['TMP_bar'] = self.TMP_bar
        Design['membrane_cost_USD_per_m2'] = self.membrane_cost_USD_per_m2
        Design['membrane_lifetime_years'] = self.membrane_lifetime_years
        Design['equipment_lifetime_years'] = self.equipment_lifetime_years

        self.feed_pump = bst.Pump(None, None, P=self.TMP_bar * 1e5)
        self.feed_pump.ins[0] = self.ins[0].copy()
        self.feed_pump.simulate()
        self.feed_pump._design()
        Design['pump_efficiency'] = 0.85 * 100.0
        self.power_utility.rate = self.feed_pump.power_utility.rate
    
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
            
        self.baseline_purchase_costs['Pump'] = self.feed_pump.purchase_cost
    

import flexsolve as flx
from biosteam import main_flowsheet as F
class IonExchangeCycle(bst.Unit):
    """
    Simulates a complete Ion Exchange (IEX) chromatography cycle as a
    steady-state equivalent unit operation.

    This unit models the loading of a feed stream, binding of target products,
    and subsequent washing, elution, and regeneration steps. It calculates the
    required volumes of all buffers and chemicals based on the column size,
    which is determined by the feed rate and target product concentration.
    The waste streams are separated by cycle phase for downstream processing.

    Parameters
    ----------
    ins :
        [0] Feed: The process stream containing the target product and impurities.
        [1] Equilibration/Wash Buffer (Buffer A): Low-salt buffer.
        [2] Elution Buffer (Buffer B): High-salt buffer to elute the product.
        [3] Regeneration Solution: e.g., NaOH or high molar NaCl.
    outs :
        [0] Product: Concentrated target product in elution buffer.
        [1] FlowthroughWaste: Unbound components from the feed stream.
        [2] WashWaste: Spent equilibration and wash buffer.
        [3] RegenerationWaste: Spent regeneration solution with stripped impurities.

    cycle_time_hr : float
        Total duration of one complete IEX cycle (loading to regeneration), in hours.
    equilibration_CV : float
        Volume of equilibration buffer used, in multiples of Column Volumes (CV).
    wash_CV : float
        Volume of wash buffer used, in multiples of Column Volumes (CV).
    elution_CV : float
        Volume of elution buffer used, in multiples of Column Volumes (CV).
    regeneration_CV : float
        Volume of regeneration solution used, in multiples of Column Volumes (CV).
    resin_DBC_g_L : float
        Dynamic Binding Capacity of the resin in grams of product per liter of resin.
    load_safety_factor : float
        Safety factor for column loading, typically < 1.0 (e.g., 0.8 means column is
        loaded to 80% of its DBC).
    TargetProduct_IDs : tuple[str]
        IDs of the chemical(s) to be captured and eluted as product.
    TargetProduct_Yield : float
        Fraction of the target product in the feed that is recovered in the product stream.
    BoundImpurity_IDs : tuple[str]
        IDs of impurities that also bind to the resin but are mostly removed.
    BoundImpurity_Removal : float
        Fraction of bound impurities in the feed that are sent to the waste stream.
    NonBinding_Carryover : float
        Fraction of non-binding components that ends up in the product stream.

    # ... other costing parameters ...
    """
    _N_ins = 4
    _N_outs = 4 # Changed from 2 to 4

    _F_BM_default = {
        'IEX Column Hardware': 2.5,
        'IEX Resin': 1.5,
        'Pump': 1.89,
    }
    # Auxiliary units for costing
    _auxiliary_unit_names = ('pump',)

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 # Cycle parameters
                 cycle_time_hr=4.0,
                 equilibration_CV=5.0,
                 wash_CV=5.0,
                 elution_CV=3.0,
                 regeneration_CV=5.0,

                 # Sizing and separation parameters
                 resin_DBC_g_L=50.0,
                 load_safety_factor=0.8,
                 TargetProduct_IDs=('Leghemoglobin',),
                 TargetProduct_Yield=0.95,
                 BoundImpurity_IDs=('Heme_b',),
                 BoundImpurity_Removal=0.93,
                 NonBinding_Carryover=0.04,

                 # Costing parameters
                 resin_cost_USD_per_L=30.0,
                # water treatment
                    # strong acid cation resin: 2~5
                    # weak acid cation resin: 3 ~ 6
                    # strong base anion resin: 7~15
                    # weak base anion resin: 5~10
                # for food tech/beverage/nutrition
                # https://lanxess.com/en/products-and-brands/brands/lewatit/industries/food-and-beverage#:~:text=Ion%20exchange%20resins%20for%20the,and%20food%20and%20drink%20additives.
                # https://www.purolite.com/index/bioprocessing/ion-exchange-chromatography#:~:text=Purolite%20Praesto%20resins%20are%20manufactured,over%20a%20wide%20pH%20range.
                # https://waterfilter.net.au/products/softening-resin-lewatit-s1567#:~:text=Ask%20a%20Question,a%20recurring%20or%20deferred%20purchase.
                    # strong acid cation resin: 10~50
                    # weak acid cation resin: 15 ~ 60
                    # strong base anion resin: 30~100
                    # weak base anion resin: 25~80
                # https://www.cytivalifesciences.com/en/us/shop/chromatography/resins/ion-exchange?sort=NameAsc&chunk=1
                #DuPont (FilmTec™) 
                # Suez Water Technologies & Solutions (formerly GE Water & Process Technologies)
                # Synder Filtration
                # Pall Corporation
                # Sartorius Stedim Biotech
                # Merck Millipore
                # Koch Separation Solutions
                # Toray Industries, Inc.
                # Asahi Kasei
                # Pentair

                 resin_lifetime_years=5.0, # 1~10 years
                 column_hardware_cost_factor=30000.0,
                 column_hardware_cost_exponent=0.6
                ):
        bst.Unit.__init__(self, ID, ins, outs, thermo)

        # Operating parameters
        self.cycle_time_hr = cycle_time_hr
        self.equilibration_CV = equilibration_CV
        self.wash_CV = wash_CV
        self.elution_CV = elution_CV
        self.regeneration_CV = regeneration_CV

        # Separation parameters
        self.resin_DBC_g_L = resin_DBC_g_L
        self.load_safety_factor = load_safety_factor
        if isinstance(TargetProduct_IDs, str): self.TargetProduct_IDs = (TargetProduct_IDs,)
        else: self.TargetProduct_IDs = tuple(TargetProduct_IDs)
        self.TargetProduct_Yield = TargetProduct_Yield
        self.BoundImpurity_IDs = tuple(BoundImpurity_IDs)
        self.BoundImpurity_Removal = BoundImpurity_Removal
        self.NonBinding_Carryover = NonBinding_Carryover

        # Costing
        self.resin_cost_USD_per_L = resin_cost_USD_per_L
        self.resin_lifetime_years = resin_lifetime_years
        self.column_hardware_cost_factor = column_hardware_cost_factor
        self.column_hardware_cost_exponent = column_hardware_cost_exponent

        # Initialize internal pump for power calculation
        self.pump = bst.Pump(None, P=4 * 1e5) # Assume 4 bar pressure drop

    def _run(self):
        # Unpack streams
        feed, buffer_A, buffer_B, regen_sol = self.ins
        product, ft_waste, wash_waste, regen_waste = self.outs

        # --- Initialize Output Streams ---
        # The product stream is primarily elution buffer
        product.copy_like(buffer_B)
        # The flowthrough waste initially contains everything from the feed
        ft_waste.copy_like(feed)
        # The wash waste is simply the spent wash buffer
        wash_waste.copy_like(buffer_A)
        # The regeneration waste is primarily the regeneration solution
        regen_waste.copy_like(regen_sol)

        # Create sets for faster checking
        target_ids = set(self.TargetProduct_IDs)
        bound_impurity_ids = set(self.BoundImpurity_IDs)

        # --- Distribute Solutes from the Feed Stream ---
        # Iterate over all components in the feed and assign them to the correct output
        for chem in self.chemicals:
            solute_in_feed = feed.imass[chem.ID]
            if solute_in_feed < 1e-12:
                continue

            if chem.ID in target_ids:
                # Target product binds, then elutes or is stripped in regeneration
                to_product = solute_in_feed * self.TargetProduct_Yield
                to_regen = solute_in_feed - to_product

                product.imass[chem.ID] += to_product
                regen_waste.imass[chem.ID] += to_regen
                # Remove the entire component from flowthrough, as it binds
                ft_waste.imass[chem.ID] = 0

            elif chem.ID in bound_impurity_ids:
                # Bound impurities bind, then are mostly stripped in regeneration
                to_regen = solute_in_feed * self.BoundImpurity_Removal
                to_product = solute_in_feed - to_regen

                product.imass[chem.ID] += to_product
                regen_waste.imass[chem.ID] += to_regen
                # Remove the entire component from flowthrough, as it binds
                ft_waste.imass[chem.ID] = 0

            else:
                # Non-binding components mostly go to flowthrough, with some carryover
                # This applies to water, salts in the feed, and other impurities
                carryover_to_product = solute_in_feed * self.NonBinding_Carryover
                
                product.imass[chem.ID] += carryover_to_product
                # The rest remains in the flowthrough, so we subtract the carryover part
                ft_waste.imass[chem.ID] -= carryover_to_product

        # --- Set Final Stream Temperatures ---
    
        product.T = buffer_B.T
        ft_waste.T = feed.T
        wash_waste.T = buffer_A.T
        regen_waste.T = regen_sol.T

    def _design(self):
        # This method remains unchanged, as it is based on inputs.
        D = self.design_results

        # Calculate resin volume based on target product loading
        target_mass_per_cycle_kg = sum(self.ins[0].imass[ID] for ID in self.TargetProduct_IDs) * self.cycle_time_hr
        effective_DBC = self.resin_DBC_g_L * self.load_safety_factor
        if effective_DBC > 0 and target_mass_per_cycle_kg > 0:
            resin_volume_L = (target_mass_per_cycle_kg * 1000) / effective_DBC
        else:
            resin_volume_L = 1.0  # Minimum volume for costing

        D['Resin Volume (CV)'] = resin_volume_L
        D['Cycle time (hr)'] = self.cycle_time_hr
        D['Target Product Yield'] = self.TargetProduct_Yield * 100
        D['Bound Impurity Removal'] = self.BoundImpurity_Removal * 100

        # Store CV values in design results
        D['Equilibration CV'] = self.equilibration_CV
        D['Wash CV'] = self.wash_CV
        D['Elution CV'] = self.elution_CV
        D['Regeneration CV'] = self.regeneration_CV

        # Cost the pump based on total liquid processed
        total_flow_stream = self.ins[0].copy()
        for s in self.ins[1:]:
            total_flow_stream.mass += s.mass

        if total_flow_stream.F_mass > 0:
            self.pump.ins[0] = total_flow_stream
            self.pump.simulate()
            self.pump._design()
            self.power_utility.rate = self.pump.power_utility.rate

    def _cost(self):
        # This method remains unchanged.
        C = self.baseline_purchase_costs
        D = self.design_results

        resin_volume_L = D.get('Resin Volume (CV)', 0.0)

        # Cost of the column hardware
        if resin_volume_L > 0:
            hw_cost = self.column_hardware_cost_factor * (resin_volume_L ** self.column_hardware_cost_exponent)
            C['IEX Column Hardware'] = bst.CE / 500 * hw_cost
        else:
            C['IEX Column Hardware'] = 0.0

        # Cost of the resin (initial fill + replacements)
        if self.resin_lifetime_years > 0 and resin_volume_L > 0:
            plant_lifetime = getattr(F.stream, 'plant_life', 20)
            replacements = plant_lifetime / self.resin_lifetime_years
            resin_cost = self.resin_cost_USD_per_L * resin_volume_L * replacements
            C['IEX Resin'] = resin_cost
        else:
            C['IEX Resin'] = 0.0

        # Cost of the pump
        self.pump._cost()
        C['Pump'] = self.pump.purchase_cost


class ResinAdsorption(IonExchangeCycle): pass

class SprayDryer(bst.SprayDryer):pass

@cost('Flow rate', 'Tank', S=1171, units='kg/hr',
        CE=522, cost=196000, n=0.7, BM=2)
class AmmoniaStorageTank(bst.StorageTank): pass


@cost('Flow rate', 'Tank', S=1981, units='kg/hr',
        CE=522, cost=96000, n=0.7, BM=1.5)
@cost('Flow rate', 'Pump', S=1981, units='kg/hr',
        CE=522, cost=7493, n=0.8, BM=2.3, kW=0.5)
class SulfuricAcidStorageTank(bst.StorageTank): pass


@cost('Flow rate', 'Pump', S=43149, units='kg/hr',
        CE=522, cost=8200, n=0.8, BM=2.3, kW=7.457)
@cost('Flow rate', 'Tank', S=40414, units='kg/hr',
        CE=522, cost=439000, n=0.7, BM=1.8)
@cost('Flow rate', 'Agitator', S=40414, units='kg/hr',
        CE=522, cost=31800, n=0.5, BM=1.5, kW=11.3205)
class SeedHoldTank(bst.Mixer): pass



class NeutralizationTank1(bst.Unit):
    _N_ins = 2
    _N_outs = 1
    
    _F_BM_default = {
        'Tank': 2.0,
        'Agitator': 1.5,
        'Cooler': 2.3,
    }
    
    _units = {
        'Flow rate': 'kg/hr',
        'Tank volume': 'm3',
        'Agitator power': 'kW',
        'Cooling duty': 'kJ/hr',
    }
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, T=None, reactions=None, 
                 residence_time=2.0, agitator_kW_per_m3=0.5):
        super().__init__(ID, ins, outs, thermo)
        self.T = T or 298.15  # Default temperature if not specified
        self.residence_time = residence_time  # hours
        self.agitator_kW_per_m3 = agitator_kW_per_m3  # kW per m3 of tank volume
        
        # Initialize reactions after thermo is available
        if reactions is None:
            chemicals = self.chemicals
            self.reaction1 = bst.Rxn(
                'H2SO4 + 2 NaOH -> Na2SO4 + 2 H2O', reactant='NaOH', X=1, chemicals=chemicals
            )
            self.reaction2 = bst.Rxn(
                'H2SO4 + Na2SO4 -> 2 NaHSO4', reactant='H2SO4', X=1, chemicals=chemicals
            )
            self.reactions = bst.SRxn([self.reaction1, self.reaction2])
        else:
            self.reactions = reactions
        
        # Initialize auxiliary units
        self.cooler = bst.HXutility(None, None, T=self.T)

    def _run(self):
        feed1, feed2 = self.ins
        out, = self.outs
        
        # Mix the two inlet streams
        out.mix_from(self.ins)
        
        # Store initial enthalpy before reaction
        H_before = out.H
        
        # Run reactions if specified
        if self.reactions:
            if hasattr(self.reactions, '__iter__'):
                # Multiple reactions
                for reaction in self.reactions:
                    reaction(out)
            else:
                # Single reaction
                self.reactions(out)
        
        # Calculate heat of reaction
        H_after_reaction = out.H
        heat_of_reaction = H_after_reaction - H_before
        
        # Cool to target temperature using auxiliary cooler
        self.cooler.ins[0] = out.copy()
        self.cooler.T = self.T
        self.cooler.simulate()
        
        # Update outlet stream
        out.copy_like(self.cooler.outs[0])
        out.T = self.T

    def _design(self):
        Design = self.design_results
        out = self.outs[0]
        
        # Calculate tank volume based on residence time and water flow rate only
        # to avoid molar volume calculation issues with salts
        water_mass_flow = out.imass['H2O']  # kg/hr
        if water_mass_flow > 0:
            water_density = 1000  # kg/m3 approximate density of water
            water_vol_flow = water_mass_flow / water_density  # m3/hr
            tank_volume_m3 = water_vol_flow * self.residence_time
        else:
            tank_volume_m3 = 1.0  # Minimum volume
            
        # Calculate agitator power based on tank volume
        agitator_power_kW = tank_volume_m3 * self.agitator_kW_per_m3
        
        # Get cooling duty from auxiliary cooler
        cooling_duty_kJ_hr = abs(self.cooler.Hnet) if hasattr(self.cooler, 'Hnet') else 0
        
        # Store design results
        Design['Flow rate'] = sum(s.F_mass for s in self.ins)
        Design['Tank volume'] = tank_volume_m3
        Design['Agitator power'] = agitator_power_kW
        Design['Cooling duty'] = cooling_duty_kJ_hr
        
        # Set power utility for agitator
        self.power_utility.rate = agitator_power_kW
        
        # Set heat utility for cooling (cooling water)
        if cooling_duty_kJ_hr > 0:
            self.add_heat_utility(-cooling_duty_kJ_hr, self.T)

    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        
        tank_volume = D['Tank volume']
        flow_rate = D['Flow rate']
        agitator_power = D['Agitator power']
        cooling_duty = D['Cooling duty']
        
        # Cost tank based on volume
        if tank_volume > 0:
            C['Tank'] = bst.CE / 522 * 96000 * (tank_volume / 10) ** 0.7
        
        # Cost agitator based on power
        if agitator_power > 0:
            C['Agitator'] = bst.CE / 522 * 31800 * (agitator_power / 11.3) ** 0.5
        
        # Cost cooler based on cooling duty
        if cooling_duty > 0:
            # Estimate cooler cost based on heat duty (simplified approach)
            C['Cooler'] = bst.CE / 522 * 50000 * (cooling_duty / 1e6) ** 0.6
        else:
            C['Cooler'] = 0




