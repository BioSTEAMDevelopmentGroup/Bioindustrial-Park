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
    #'PSA',

    ##### Downstream #####
    #'CellDisruption',
    #'ProteinCentrifuge',
    #'Evaporator',
    #'daifiltration',
    #'ion_exchange',
    #'nanofiltration',

    #'SprayDrying',

    #'ScrewPress','CellDisruption',
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

# %%
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

# class CellDisruption(bst.Homogenizer): pass
# pressure vessel?

#@cost('Flow rate', units='kg/hr', CE=CEPCI_by_year[2010], cost=100000, S=100000, n=0.6, kW=100)
#@copy_algorithm(bst.SolidLiquidsSplitCentrifuge, run=False)       
#class Centrifuge(bst.SpliSolidLiquidsSplitCentrifugetter): pass

class ProteinCentrifuge(bst.SolidsCentrifuge): pass

class Evaporator(bst.MultiEffectEvaporator): pass

class Diafiltration(bst.Unit):
    _N_ins = 2  # Feed and Wash Solution
    _N_outs = 2 # Permeate and Retentate

    water_ID='H2O'
    _default_TargetProduct_ID = 'Leghemoglobin'
    _default_Salt_ID = 'Salts'
    _default_OtherLargeMolecules_ID = 'OtherLargeMolecules'
    _default_DefaultSolute_ID = 'DefaultSolute'
    _default_TargetProduct_Retention = 0.98
    _default_Salt_Retention = 0.05
    _default_OtherLargeMolecules_Retention = 0.99
    _default_DefaultSolute_Retention = 0.05
    _default_FeedWater_Recovery_to_Permeate = 0.9

    def __init__(self, ID='', ins=None, outs=None, thermo=None,
            # Membrane properties
            TargetProduct_ID=None,
            Salt_ID=None,
            OtherLargeMolecules_ID=None,
            DefaultSolute_ID=None,
            TargetProduct_Retention=None,
            Salt_Retention=None,
            OtherLargeMolecules_Retention=None,
            DefaultSolute_Retention=None,
            FeedWater_Recovery_to_Permeate=None, **kwargs):
        super().__init__(ID, ins,outs,thermo)

        # set membrane properties
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
        self.OtherLargeMolecules_ID = OtherLargeMolecules_ID if OtherLargeMolecules_ID is not None else self._default_OtherLargeMolecules_ID
        self.DefaultSolute_ID = DefaultSolute_ID if DefaultSolute_ID is not None else self._default_DefaultSolute_ID
        self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
        self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
        self.OtherLargeMolecules_Retention = OtherLargeMolecules_Retention if OtherLargeMolecules_Retention is not None else self._default_OtherLargeMolecules_Retention
        self.DefaultSolute_Retention = DefaultSolute_Retention if DefaultSolute_Retention is not None else self._default_DefaultSolute_Retention
        self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate

    def _run(self):
        feed, wash_solution = self.ins
        retentate, permeate = self.outs

        retentate.T = permeate.T = feed.T
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

        current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
        if abs(current_total_water_out - total_incoming_water) > 1e-9: # Tolerance for floating point
            if total_incoming_water >= retentate.imass[self.water_ID]:
                permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
            else:
                retentate.imass[self.water_ID] = total_incoming_water
                permeate.imass[self.water_ID] = 0.0

        # --- Solute Balance ---
        for chem in self.chemicals:
            ID = chem.ID
            if ID == self.water_ID:
                continue

            mass_in_feed = feed.imass[ID]
            mass_in_wash = wash_solution.imass[ID]
            total_mass_in = mass_in_feed + mass_in_wash

            if total_mass_in <= 1e-12: # Effectively zero mass
                retentate.imass[ID] = 0.0
                permeate.imass[ID] = 0.0
                continue

            current_retention = self.DefaultSolute_Retention

            if ID == self.TargetProduct_ID:
                current_retention = self.TargetProduct_Retention
            elif self.Salt_ID and ID in self.Salt_ID:
                current_retention = self.Salt_Retention
            elif self.OtherLargeMolecules_ID and ID in self.OtherLargeMolecules_ID:
                current_retention = self.OtherLargeMolecules_Retention
            
            retentate_mass_solute = total_mass_in * current_retention
            permeate_mass_solute = total_mass_in * (1.0 - current_retention)
            
            retentate.imass[ID] = max(0.0, retentate_mass_solute)
            # Permeate by difference for better mass balance, after retentate is set
            permeate.imass[ID] = max(0.0, total_mass_in - retentate.imass[ID])


            # Final check to ensure solute mass balance due to max(0,...) or floating point nuances
            current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
            mass_balance_error = current_total_solute_out - total_mass_in
            if abs(mass_balance_error) > 1e-9 * total_mass_in and abs(mass_balance_error) > 1e-12 :
                # If there's a significant discrepancy, adjust one stream (e.g. permeate)
                # This logic can be refined, but aims to conserve mass.
                permeate.imass[ID] -= mass_balance_error 
                if permeate.imass[ID] < 0: # If adjustment makes it negative
                    retentate.imass[ID] += permeate.imass[ID] # Add the negative part to retentate
                    permeate.imass[ID] = 0.0
                    if retentate.imass[ID] < 0: # Should not happen if total_mass_in >=0
                        retentate.imass[ID] = 0.0
                        # At this point, mass is lost if total_mass_in was positive. This indicates an issue.
                        # For robustness, one might assign all to one stream if total_mass_in > 0.
                        # print(f"Warning: Mass balance issue for {ID} in {self.ID}")


    def _design(self):
        """Placeholder for design calculations (e.g., membrane area)."""
        # Example: Calculate membrane area based on permeate flux
        # self.membrane_flux_rate = 50 # L/(m^2*hr), would be a parameter
        # permeate_vol_flow_m3hr = self.outs[0].F_vol / 1000 # Assuming F_vol is in L/hr
        
        # if permeate_vol_flow_m3hr > 0 and hasattr(self, 'membrane_flux_rate') and self.membrane_flux_rate > 0:
        #     # Flux often given in L/m2/h or GFD. Ensure units are consistent.
        #     # Example: flux in L/m2/h, F_vol in L/hr
        #     required_area_m2 = self.outs[0].F_vol / self.membrane_flux_rate
        #     self.design_results['Membrane Area (m^2)'] = required_area_m2
        pass

    def _cost(self):
        """Placeholder for cost calculations."""
        # if 'Membrane Area (m^2)' in self.design_results:
        #     area = self.design_results['Membrane Area (m^2)']
        #     # Example: CEPCI = bst.CE # Chemical Engineering Plant Cost Index
        #     # self.purchase_costs['Membrane Module'] = cost_correlation_function(area) * CEPCI / 500
        pass

class IonExchange(bst.Unit): 
    _N_ins = 2  # Feed (conditioned) and Elution Buffer Profile
    _N_outs = 2 # Product (in elution buffer) and Waste Stream

    # --- Default Values for Parameters ---
    _default_water_ID = 'H2O' # Using your preferred ID for water
    _default_TargetProduct_ID = 'Leghemoglobin' # As per your context

    _default_TargetProduct_Yield = 0.90      # e.g., 90% recovery of loaded target in product pool
    
    # Impurities specifically targeted for removal by binding differently than the product
    _default_BoundImpurity_IDs_tuple = ('HostCellProtein', 'DNA', 'Endotoxin') # Example impurity IDs
    # Fraction of these 'Bound Impurities' (from feed) that are successfully removed from the product path (i.e., go to waste)
    _default_BoundImpurity_Removal_Efficiency = 0.99 # e.g., 99% removal (2 LRV)

    # For other solutes in the feed (not TargetProduct, not BoundImpurity, not ElutionBufferSalt)
    # This fraction of these other solutes (from feed) ends up in the product stream. The rest goes to waste.
    _default_NonBinding_Solutes_Carryover_to_Product = 0.05 
    
    # Key components defining the elution buffer matrix (e.g., the elution salt from ins[1])
    # If these components are also present in the feed (ins[0]), the feed's portion is assumed to go to waste.
    _default_ElutionBuffer_Defining_Component_IDs_tuple = ('NaCl', 'KCl') # Example elution salts

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                water_ID=None,
                TargetProduct_ID=None,
                TargetProduct_Yield=None,
                BoundImpurity_IDs_tuple=None, # Expects a tuple of chemical IDs
                BoundImpurity_Removal_Efficiency=None,
                NonBinding_Solutes_Carryover_to_Product=None,
                ElutionBuffer_Defining_Component_IDs_tuple=None, # Expects a tuple
                 **kwargs):
        super().__init__(ID, ins, outs, thermo, **kwargs)

        # Set operational parameters
        self.water_ID = water_ID if water_ID is not None else self._default_water_ID
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.TargetProduct_Yield = TargetProduct_Yield if TargetProduct_Yield is not None else self._default_TargetProduct_Yield
        
        self.BoundImpurity_IDs_tuple = BoundImpurity_IDs_tuple if BoundImpurity_IDs_tuple is not None else self._default_BoundImpurity_IDs_tuple
        self.BoundImpurity_Removal_Efficiency = BoundImpurity_Removal_Efficiency if BoundImpurity_Removal_Efficiency is not None else self._default_BoundImpurity_Removal_Efficiency
        
        self.NonBinding_Solutes_Carryover_to_Product = NonBinding_Solutes_Carryover_to_Product if NonBinding_Solutes_Carryover_to_Product is not None else self._default_NonBinding_Solutes_Carryover_to_Product
        self.ElutionBuffer_Defining_Component_IDs_tuple = ElutionBuffer_Defining_Component_IDs_tuple if ElutionBuffer_Defining_Component_IDs_tuple is not None else self._default_ElutionBuffer_Defining_Component_IDs_tuple

    def _run(self):
        feed = self.ins[0]
        elution_buffer_profile = self.ins[1] # Defines the matrix of the product stream
        product_stream = self.outs[0]
        waste_stream = self.outs[1]

        # Initialize output streams' temperatures
        product_stream.T = elution_buffer_profile.T # Product is in elution buffer
        waste_stream.T = feed.T                   # Waste derives from feed temperature

        # Product stream starts with the exact composition and mass of the elution_buffer_profile stream.
        # Solutes from the feed will be added to this.
        product_stream.imass = elution_buffer_profile.imass.copy()
        
        # Waste stream starts empty and will collect components from the feed that don't go to product.
        waste_stream.empty()

        # Process components from the original feed stream
        for chem in self.chemicals: # Iterate over all chemicals in the system
            ID = chem.ID
            feed_mass_component = feed.imass[ID] # Mass of this chemical in the input feed

            # If component is negligible in feed, it won't contribute from feed.
            # Its presence in product is already set if it was in elution_buffer_profile.
            # Ensure it's initialized (likely to 0) in waste if not already handled.
            if feed_mass_component <= 1e-12:
                if waste_stream.imass[ID] < 0: waste_stream.imass[ID] = 0.0 # Should be 0 from empty()
                continue

            if ID == self.TargetProduct_ID:
                mass_to_product = feed_mass_component * self.TargetProduct_Yield
                product_stream.imass[ID] += mass_to_product # Add recovered product
                waste_stream.imass[ID] = feed_mass_component - mass_to_product
            
            elif ID == self.water_ID:
                # Water from the original feed stream goes to the waste stream.
                # The product stream's water content is defined by the elution_buffer_profile.
                waste_stream.imass[ID] = feed_mass_component

            elif self.ElutionBuffer_Defining_Component_IDs_tuple and \
                ID in self.ElutionBuffer_Defining_Component_IDs_tuple:
                # If a defining component of the elution buffer (e.g., NaCl) was also in the feed,
                # that portion from the feed is assumed to go to the waste stream.
                # The concentration of this component in the product_stream is dictated solely
                # by the elution_buffer_profile stream.
                waste_stream.imass[ID] = feed_mass_component

            elif self.BoundImpurity_IDs_tuple and ID in self.BoundImpurity_IDs_tuple:
                # These are impurities targeted for removal
                mass_removed_to_waste = feed_mass_component * self.BoundImpurity_Removal_Efficiency
                waste_stream.imass[ID] = mass_removed_to_waste
                # The remainder (not removed) carries over to the product stream
                product_stream.imass[ID] += (feed_mass_component - mass_removed_to_waste)

            else: 
                # For all other solutes present in the feed (not target, not water, 
                # not elution buffer defining components, not specified bound impurities)
                mass_to_product = feed_mass_component * self.NonBinding_Solutes_Carryover_to_Product
                product_stream.imass[ID] += mass_to_product
                waste_stream.imass[ID] = feed_mass_component - mass_to_product
        
        # Ensure all chemical amounts are non-negative after all additions/subtractions
        for chem_obj in self.chemicals:
            idx = chem_obj.ID
            if product_stream.imass[idx] < 0: product_stream.imass[idx] = 0.0
            if waste_stream.imass[idx] < 0: waste_stream.imass[idx] = 0.0


    def _design(self):
        """
        Placeholder for Ion Exchange Column design.
        Key design parameters would include resin volume, column dimensions.
        """
        # Example: Calculate resin volume based on Dynamic Binding Capacity (DBC)
        # Parameters needed for this (to be added to __init__ if implementing):
        # self.resin_DBC_g_L = 50 # g of TargetProduct per L of resin (e.g.)
        # self.DBC_safety_factor = 0.8 # Operate at 80% of DBC
        # self.load_flow_rate_CV_hr = 5 # Column Volumes per hour for loading
        # self.num_cycles_per_year = 300 # For equipment sizing based on annual throughput

        # target_product_in_feed_kg_hr = self.ins[0].imass[self.TargetProduct_ID] # If continuous average
        # if hasattr(self, 'resin_DBC_g_L') and self.resin_DBC_g_L > 0:
        #     effective_DBC_kg_L = (self.resin_DBC_g_L / 1000.0) * self.DBC_safety_factor
        #     if effective_DBC_kg_L > 0:
        #          # This calculation depends on whether flow rate is per hour or per batch
        #          # For a batch process, target_product_in_feed_kg_hr would be kg/batch
        #          # resin_volume_L = (target_product_in_feed_kg_hr_or_kg_batch) / effective_DBC_kg_L
        #          # self.design_results['Resin Volume (L)'] = resin_volume_L
        #          # Height, Diameter from H/D ratio (e.g. 10-20 for process scale)
        #          pass
        pass

    def _cost(self):
        """
        Placeholder for Ion Exchange Column and Resin cost.
        Costs would be based on resin volume, column hardware.
        """
        # if 'Resin Volume (L)' in self.design_results:
        #     resin_volume_L = self.design_results['Resin Volume (L)']
        #     # Example cost factors (these would ideally be class attributes or from a config)
        #     cost_per_L_resin = 1500 # $/L (highly variable based on resin type)
        #     column_hardware_factor = 0.6 # Hardware cost as a fraction of resin cost
            
        #     self.purchase_costs['IEX Resin'] = resin_volume_L * cost_per_L_resin
        #     self.purchase_costs['IEX Column Hardware'] = self.purchase_costs['IEX Resin'] * column_hardware_factor
        pass

class Nanofiltration(bst.Unit):
    _N_ins = 2  # Feed and Wash Solution
    _N_outs = 2 # Permeate and Retentate

    water_ID='H2O'
    _default_TargetProduct_ID = 'Leghemoglobin'
    _default_Salt_ID = 'Salts'
    _default_OtherLargeMolecules_ID = 'OtherLargeMolecules'
    _default_DefaultSolute_ID = 'DefaultSolute'
    _default_TargetProduct_Retention = 0.999
    _default_Salt_Retention = 0.01
    _default_OtherLargeMolecules_Retention = 0.8
    _default_DefaultSolute_Retention = 0.01
    _default_FeedWater_Recovery_to_Permeate = 0.9

    def __init__(self, ID='', ins=None, outs=None, thermo=None,
            # Membrane properties
            TargetProduct_ID=None,
            Salt_ID=None,
            OtherLargeMolecules_ID=None,
            DefaultSolute_ID=None,
            TargetProduct_Retention=None,
            Salt_Retention=None,
            OtherLargeMolecules_Retention=None,
            DefaultSolute_Retention=None,
            FeedWater_Recovery_to_Permeate=None, **kwargs):
        super().__init__(ID, ins,outs,thermo)

        # set membrane properties
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
        self.OtherLargeMolecules_ID = OtherLargeMolecules_ID if OtherLargeMolecules_ID is not None else self._default_OtherLargeMolecules_ID
        self.DefaultSolute_ID = DefaultSolute_ID if DefaultSolute_ID is not None else self._default_DefaultSolute_ID
        self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
        self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
        self.OtherLargeMolecules_Retention = OtherLargeMolecules_Retention if OtherLargeMolecules_Retention is not None else self._default_OtherLargeMolecules_Retention
        self.DefaultSolute_Retention = DefaultSolute_Retention if DefaultSolute_Retention is not None else self._default_DefaultSolute_Retention
        self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate

    def _run(self):
        feed, wash_solution = self.ins
        retentate, permeate = self.outs

        retentate.T = permeate.T = feed.T
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

        current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
        if abs(current_total_water_out - total_incoming_water) > 1e-9: # Tolerance for floating point
            if total_incoming_water >= retentate.imass[self.water_ID]:
                permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
            else:
                retentate.imass[self.water_ID] = total_incoming_water
                permeate.imass[self.water_ID] = 0.0

        # --- Solute Balance ---
        for chem in self.chemicals:
            ID = chem.ID
            if ID == self.water_ID:
                continue

            mass_in_feed = feed.imass[ID]
            mass_in_wash = wash_solution.imass[ID]
            total_mass_in = mass_in_feed + mass_in_wash

            if total_mass_in <= 1e-12: # Effectively zero mass
                retentate.imass[ID] = 0.0
                permeate.imass[ID] = 0.0
                continue

            current_retention = self.DefaultSolute_Retention

            if ID == self.TargetProduct_ID:
                current_retention = self.TargetProduct_Retention
            elif self.Salt_ID and ID in self.Salt_ID:
                current_retention = self.Salt_Retention
            elif self.OtherLargeMolecules_ID and ID in self.OtherLargeMolecules_ID:
                current_retention = self.OtherLargeMolecules_Retention
            
            retentate_mass_solute = total_mass_in * current_retention
            permeate_mass_solute = total_mass_in * (1.0 - current_retention)
            
            retentate.imass[ID] = max(0.0, retentate_mass_solute)
            # Permeate by difference for better mass balance, after retentate is set
            permeate.imass[ID] = max(0.0, total_mass_in - retentate.imass[ID])


            # Final check to ensure solute mass balance due to max(0,...) or floating point nuances
            current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
            mass_balance_error = current_total_solute_out - total_mass_in
            if abs(mass_balance_error) > 1e-9 * total_mass_in and abs(mass_balance_error) > 1e-12 :
                # If there's a significant discrepancy, adjust one stream (e.g. permeate)
                # This logic can be refined, but aims to conserve mass.
                permeate.imass[ID] -= mass_balance_error 
                if permeate.imass[ID] < 0: # If adjustment makes it negative
                    retentate.imass[ID] += permeate.imass[ID] # Add the negative part to retentate
                    permeate.imass[ID] = 0.0
                    if retentate.imass[ID] < 0: # Should not happen if total_mass_in >=0
                        retentate.imass[ID] = 0.0
                        # At this point, mass is lost if total_mass_in was positive. This indicates an issue.
                        # For robustness, one might assign all to one stream if total_mass_in > 0.
                        # print(f"Warning: Mass balance issue for {ID} in {self.ID}")


    def _design(self):
        """Placeholder for design calculations (e.g., membrane area)."""
        # Example: Calculate membrane area based on permeate flux
        # self.membrane_flux_rate = 50 # L/(m^2*hr), would be a parameter
        # permeate_vol_flow_m3hr = self.outs[0].F_vol / 1000 # Assuming F_vol is in L/hr
        
        # if permeate_vol_flow_m3hr > 0 and hasattr(self, 'membrane_flux_rate') and self.membrane_flux_rate > 0:
        #     # Flux often given in L/m2/h or GFD. Ensure units are consistent.
        #     # Example: flux in L/m2/h, F_vol in L/hr
        #     required_area_m2 = self.outs[0].F_vol / self.membrane_flux_rate
        #     self.design_results['Membrane Area (m^2)'] = required_area_m2
        pass

    def _cost(self):
        """Placeholder for cost calculations."""
        # if 'Membrane Area (m^2)' in self.design_results:
        #     area = self.design_results['Membrane Area (m^2)']
        #     # Example: CEPCI = bst.CE # Chemical Engineering Plant Cost Index
        #     # self.purchase_costs['Membrane Module'] = cost_correlation_function(area) * CEPCI / 500
        pass

class Nanofiltration(bst.PressureFilter): pass

class SprayDryer(bst.SprayDryer): pass

# class ScrewPress(bst.ScrewPress): pass
