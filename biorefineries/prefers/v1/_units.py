# -*- coding: utf-8 -*-
"""
Created on 2025-04-18 15:20:45
Last Updated: 2026-01-20 (Parameter Validation & Literature Annotation)

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg

Literature References:
----------------------
[1] Yao et al. (2025). "Recent advances in microbial fermentation for the production
    of leghemoglobin." World J Microbiol Biotechnol 41:404. DOI:10.1007/s11274-025-04609-y
    - LegHb fermentation parameters: titer 7.27 g/L (K. marxianus), yield optimization
    - Heme supply strategies, pH 5.0-5.5 optimal for yeast
    - Temperature 28-30°C optimal for K. phaffii/P. pastoris

[2] Perry's Chemical Engineers' Handbook, 8th Ed.
    - Filtration: Table 18-5, solids loading 10-50 kg/m²/hr for RDVF
    - Centrifugation design parameters

[3] DuPont Water Solutions / Axeon Water Technologies
    - RO water recovery: 75-95% (industrial range)
    - Membrane lifetime: 3-5 years
    - Operating pressure: 15-30 bar standard

[4] Membranes.com / Synder Filtration
    - UF flux: 20-80 LMH (ultrafiltration)
    - NF flux: 10-40 LMH (nanofiltration)
    - TMP: 2-10 bar

[5] SNS Insider (2023) / Samcotech
    - Membrane costs: $80-$350/m² (industrial), $1500-$5000/m² (pharma)

[6] Walch & Jungbauer (2017). "Continuous desalting of refolded protein solution."
    Biotechnol J 12:1700082. DOI:10.1002/biot.201700082
    - Ion exchange chromatography parameters for protein capture

[7] NREL Technical Report TP-5100-47764 (Humbird et al. 2011)
    - Process economics, installation factors, CEPCI correlations
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
from biosteam.exceptions import lb_warning
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
    
    # Gas Separation
    'VacuumPSA',

    ##### Downstream #####
    'CellDisruption',
    'ProteinCentrifuge',
    'Evaporator',
    'Diafiltration',
    'ResinColumn',
    'ReverseOsmosis',
    'NanofiltrationDF',
    'SprayDrying',
    'BoilerTurbogenerator',
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

class BoilerTurbogenerator(bst.BoilerTurbogenerator):
    """BoilerTurbogenerator with guarded emissions enthalpy update."""

    def _design(self):
        B_eff = self.boiler_efficiency
        TG_eff = self.turbogenerator_efficiency
        steam_demand = self.steam_demand
        Design = self.design_results
        chemicals = self.chemicals
        self._load_utility_agents()
        mol_steam = sum([i.flow for i in self.steam_utilities])
        feed_solids, feed_gas, makeup_water, fuel, lime, chems, oxygen_rich_gas = self.ins
        oxygen_rich_gas.empty()
        if self.fuel_source == 'CH4':
            fuel.phase = 'g'
            fuel.set_property('T', 60, 'degF')
            fuel.set_property('P', 14.73, 'psi')
        emissions, blowdown_water, ash_disposal = self.outs
        if not lime.price:
            lime.price = 0.19937504680689402
        if not chems.price:
            chems.price = 4.995862254032183
        H_steam = sum([i.duty for i in self.steam_utilities])
        side_steam = self.side_steam
        if side_steam: 
            H_steam += side_steam.H
            mol_steam += side_steam.F_mol
        steam_demand.imol['7732-18-5'] = mol_steam 
        self.combustion_reactions = combustion_rxns = chemicals.get_combustion_reactions()
        non_empty_feeds = [i for i in (feed_solids, feed_gas) if not i.isempty()]
        boiler_efficiency_basis = self.boiler_efficiency_basis
        fuel_source = self.fuel_source
        def calculate_excess_electricity_at_natual_gas_flow(fuel_flow):
            if fuel_flow:
                fuel_flow = abs(fuel_flow)
                fuel.imol[fuel_source] = fuel_flow
            else:
                fuel.empty()
            if boiler_efficiency_basis == 'LHV':
                H_combustion = fuel.LHV
                for feed in non_empty_feeds: H_combustion += feed.LHV
            elif boiler_efficiency_basis == 'HHV':
                H_combustion = fuel.HHV
                for feed in non_empty_feeds: H_combustion += feed.HHV
            else:
                raise ValueError(
                    f"invalid boiler efficiency basis {boiler_efficiency_basis}; "
                    f"valid values include 'LHV', or 'HHV'"
                )
            self.H_content = H_content = B_eff * H_combustion 
            self.H_loss_to_emissions = H_combustion - H_content
            H_electricity = H_content - H_steam # Heat available for the turbogenerator
            electricity = H_electricity * TG_eff  # Electricity produced
            self.cooling_duty = electricity - H_electricity
            Design['Work'] = work = electricity / 3600
            duty_over_mol = 39000 # kJ / mol-superheated steam 
            self.total_steam = H_content / duty_over_mol #: [float] Total steam produced by the boiler (kmol/hr)
            Design['Flow rate'] = flow_rate = self.total_steam * 18.01528
            
            if self.satisfy_system_electricity_demand:
                boiler = self.cost_items['Boiler']
                rate_boiler = boiler.kW * flow_rate / boiler.S
                return work - self.electricity_demand - rate_boiler
            else:
                return work
        
        self._excess_electricity_without_fuel = excess_electricity = calculate_excess_electricity_at_natual_gas_flow(0)
        if excess_electricity < 0:
            f = calculate_excess_electricity_at_natual_gas_flow
            lb = 0.
            fuel.imol[fuel_source] = 1
            ub = - excess_electricity * 3600 / fuel.LHV
            while f(ub) < 0.: 
                lb = ub
                ub *= 2
            flx.IQ_interpolation(f, lb, ub, xtol=1, ytol=1)
        
        if self.cooling_duty > 0.: 
            # In the event that no electricity is produced and the solver
            # solution for natural gas is slightly below the requirement for steam
            # (this would lead to a positive duty).
            self.cooling_duty = 0.
            Design['Work'] = 0.
        
        emissions.T = 298.15 # Will be updated later with the energy balance
        emissions.P = 101325
        emissions.phase = 'g'
        emissions_mol = emissions.mol
        emissions_mol[:] = fuel.mol
        for feed in non_empty_feeds: emissions_mol[:] += feed.mol
        combustion_rxns.force_reaction(emissions_mol)
        O2_consumption = -emissions.imol['O2']
        oxygen_rich_gas.reset_flow(**self.oxygen_rich_gas_composition)
        z_O2 = oxygen_rich_gas.imol['O2'] / oxygen_rich_gas.F_mol
        oxygen_rich_gas.F_mol = O2_consumption / z_O2
        emissions_mol += oxygen_rich_gas.mol
        F_emissions = emissions.F_mass
        z_CO2 = emissions.imass['CO2'] / F_emissions
        z_CO2_target = self.CO2_emissions_concentration
        if z_CO2 > z_CO2_target:
            F_emissions_new = z_CO2 * F_emissions / z_CO2_target
            dF_emissions = F_emissions_new - F_emissions
            oxygen_rich_gas.F_mass = F_mass_O2_new = oxygen_rich_gas.F_mass + dF_emissions
            emissions_mol += oxygen_rich_gas.mol * (dF_emissions / F_mass_O2_new)
        try:
            H_base = emissions.H
        except Exception:
            H_base = 0.0
        emissions._property_cache['H'] = H_base + self.H_loss_to_emissions
        hu_cooling = bst.HeatUtility()
        hu_cooling(self.cooling_duty, steam_demand.T)
        hus_heating = bst.HeatUtility.sum_by_agent(self.steam_utilities)
        for hu in hus_heating: hu.reverse()
        self.heat_utilities = [*hus_heating, hu_cooling]
        water_index = chemicals.index('7732-18-5')
        makeup_water.mol[water_index] = blowdown_water.mol[water_index] = (
                self.total_steam * self.boiler_blowdown * 1 / (1 - self.RO_rejection)   
        )
        ash_IDs = [i.ID for i in self.chemicals if not i.formula]
        emissions_mol = emissions.mol
        SO2_produced = 0
        if 'SO2' in chemicals:
            SO2_produced += emissions.imol['SO2']
        if 'CaSO4' in chemicals:
            SO2_produced += emissions.imol['CaSO4']
            ash_IDs.append('CaSO4')
        if SO2_produced: 
            rxn, ID_lime = self._get_desulfurization_rxn_and_coreactant()
            # FGD lime scaled based on SO2 generated,	
            # 20% stoichiometric excess based on P52 of ref [1]
            rxn.force_reaction(emissions)
            lime.imol[ID_lime] = lime_mol = SO2_produced * 1.2
            emissions_mol.remove_negatives()
        else:
            lime.empty()
        # About 0.4536 kg/hr of boiler chemicals are needed per 234484 kg/hr steam produced
        chems.imol['Ash'] = boiler_chems = 1.9345e-06 * Design['Flow rate']
        ash_disposal.empty()
        ash_disposal.copy_flow(emissions, IDs=tuple(ash_IDs), remove=True)
        ash_disposal.imol['Ash'] += boiler_chems
        dry_ash = ash_disposal.F_mass
        moisture = min(emissions.imass['Water'], dry_ash * 0.3) # ~20% moisture
        ash_disposal.imass['Water'] = moisture
        emissions.imass['Water'] -= moisture
        Design['Ash disposal'] = dry_ash + moisture
        if SO2_produced:
            if ID_lime == 'Ca(OH)2': # Ca(OH)2
                lime.imol['Water'] = 4 * lime_mol # Its a slurry
            else: # CaO
                lime.imol['Water'] = 5 * lime_mol

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


class VacuumPSA(bst.Unit):
    """
    Vacuum Pressure Swing Adsorption (VPSA) unit for gas separation.
    
    Separates gas mixtures using cyclic adsorption/desorption on adsorbent materials.
    Designed for syngas component separation (H2, CO, C2H4, etc.) using zeolite 13X.
    
    Parameters
    ----------
    ID : str, optional
        Unit identifier.
    ins : Stream, optional
        Mixed gas feed stream.
    outs : Stream, optional
        * [0] product: H2-rich raffinate stream
        * [1] purge: CO/C2H4-rich tail gas (extract)
    split : dict, optional
        Component-specific split factors to product stream.
        Default splits are based on zeolite 13X selectivity.
    P_ads : float, optional
        Adsorption pressure in Pa. Default is 6e5 (6 bar).
    P_des : float, optional
        Desorption pressure in Pa. Default is 1e4 (0.1 bar vacuum).
    cycle_time : float, optional
        Total cycle duration in seconds. Default is 600 (10 min).
    N_beds : int, optional
        Number of adsorbent beds. Default is 2 for continuous operation.
    adsorbent_loading : float, optional
        Adsorbent loading capacity in mol/kg. Default is 2.0.
    adsorbent_bulk_density : float, optional
        Adsorbent bulk density in kg/m³. Default is 650.
    vacuum_efficiency : float, optional
        Vacuum pump efficiency. Default is 0.70.
    adsorbent_cost : float, optional
        Adsorbent cost in USD/kg. Default is 5.0.
    
    Notes
    -----
    This is a steady-state approximation of the cyclic PSA process, suitable for 
    TEA/LCA screening but not for detailed process design.
    
    References
    ----------
    [1] Ruthven, D.M., Farooq, S., Knaebel, K.S. "Pressure Swing Adsorption", 
        VCH Publishers, 1994 (Design basis)
    [2] Yang, R.T. "Gas Separation by Adsorption Processes", Imperial College Press, 
        1997 (Zeolite 13X selectivity data)
    [3] Naquash et al., "Hydrogen purification from syngas", ICAE 2021 (Application)
    [4] Peters, M.S., et al. "Plant Design and Economics for Chemical Engineers", 
        5th Ed., 2003 (Costing power laws)
    [5] Humbird et al., NREL/TP-5100-47764, 2011 (Installation factors)
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['H2', 'CO', 'C2H4', 'CO2', 'N2'])
    >>> feed = bst.Stream('feed', H2=50, CO=30, C2H4=15, CO2=3, N2=2, 
    ...                   units='kmol/hr', phase='g')
    >>> V = VacuumPSA('V1', ins=feed)
    >>> V.simulate()
    >>> V.show()
    """
    
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    
    _units = {
        'Feed gas flow': 'kmol/hr',
        'Adsorbent mass': 'kg',
        'Bed volume': 'm3',
        'Vacuum power': 'kW',
        'H2 recovery': '%',
    }
    
    _F_BM_default = {
        'Pressure vessels': 2.5,
        'Vacuum pump': 1.8,
        'Adsorbent': 1.0,
        'Piping and valves': 1.0,
    }
    
    # Default split factors based on zeolite 13X selectivity
    _default_split = {
        'H2': 0.95,      # Low adsorption affinity -> goes to product
        'CO': 0.30,      # Moderate-high affinity -> mostly to purge
        'C2H4': 0.20,    # High affinity (π-bonding) -> mostly to purge
        'CO2': 0.10,     # Strong quadrupole interaction -> mostly to purge
        'N2': 0.85,      # Low affinity -> goes to product
        'CH4': 0.60,     # Moderate affinity
        'O2': 0.80,      # Low affinity
        'H2O': 0.05,     # Strong affinity -> mostly to purge
        'Ethanol': 0.15, # High affinity
        'Ethylene': 0.20, # Alias for C2H4
    }
    
    def _init(self,
              split=None,
              P_ads=6e5,           # Pa (6 bar)
              P_des=1e4,           # Pa (0.1 bar vacuum)
              cycle_time=600,      # s (10 min)
              N_beds=2,
              adsorbent_loading=2.0,      # mol/kg
              adsorbent_bulk_density=650, # kg/m³
              vacuum_efficiency=0.70,
              adsorbent_cost=5.0):        # USD/kg
        
        # Use default splits if not provided
        if split is None:
            self.split = self._default_split.copy()
        else:
            # Merge with defaults (user values override)
            self.split = self._default_split.copy()
            self.split.update(split)
        
        self.P_ads = P_ads
        self.P_des = P_des
        self.cycle_time = cycle_time
        self.N_beds = N_beds
        self.adsorbent_loading = adsorbent_loading
        self.adsorbent_bulk_density = adsorbent_bulk_density
        self.vacuum_efficiency = vacuum_efficiency
        self.adsorbent_cost = adsorbent_cost
    
    def _run(self):
        """Perform mass balance using split factors."""
        feed = self.ins[0]
        product, purge = self.outs
        
        # Initialize output streams
        product.empty()
        purge.empty()
        
        # Set phase to gas
        product.phase = 'g'
        purge.phase = 'g'
        
        # Apply split factors to each component
        for chem in feed.chemicals:
            chem_id = chem.ID
            feed_mol = feed.imol[chem_id]
            
            if feed_mol <= 0:
                continue
            
            # Get split factor (default to 0.5 for unknown components)
            split_factor = self.split.get(chem_id, 0.5)
            
            # Split to product and purge
            product.imol[chem_id] = feed_mol * split_factor
            purge.imol[chem_id] = feed_mol * (1 - split_factor)
        
        # Set temperature and pressure for outlets
        product.T = feed.T
        product.P = self.P_ads  # Product at adsorption pressure
        purge.T = feed.T
        purge.P = self.P_des    # Purge at desorption pressure
    
    def _design(self):
        """Calculate design requirements and utilities."""
        D = self.design_results
        feed = self.ins[0]
        purge = self.outs[1]
        
        # Feed gas flow
        D['Feed gas flow'] = feed.F_mol
        
        # Adsorbent mass calculation
        # Based on gas loading capacity and cycle time
        # m_ads = F_mol * cycle_time / (3600 * q_loading)
        F_mol_per_s = feed.F_mol / 3600  # kmol/s
        m_ads = F_mol_per_s * self.cycle_time / (self.adsorbent_loading / 1000)  # kg per bed
        D['Adsorbent mass'] = m_ads * self.N_beds
        
        # Bed volume
        V_bed = m_ads / self.adsorbent_bulk_density  # m³ per bed
        D['Bed volume'] = V_bed * self.N_beds
        
        # Vacuum pump power
        # P = F_purge * R * T * ln(P_ads/P_des) / η
        R = 8.314  # J/(mol·K)
        F_purge_mol_s = purge.F_mol / 3600  # kmol/s -> mol/s * 1000
        compression_work = F_purge_mol_s * 1000 * R * feed.T * np.log(self.P_ads / self.P_des)
        power_kW = compression_work / (1000 * self.vacuum_efficiency)
        D['Vacuum power'] = power_kW
        
        # Register power utility
        self.power_utility.consumption = power_kW
        
        # H2 recovery (if H2 present)
        if 'H2' in feed.chemicals.IDs and feed.imol['H2'] > 0:
            H2_recovery = self.outs[0].imol['H2'] / feed.imol['H2'] * 100
            D['H2 recovery'] = H2_recovery
        else:
            D['H2 recovery'] = 0
    
    def _cost(self):
        """Calculate equipment purchase costs."""
        D = self.design_results
        C = self.baseline_purchase_costs
        
        V_bed_total = D.get('Bed volume', 0)
        power_kW = D.get('Vacuum power', 0)
        m_ads_total = D.get('Adsorbent mass', 0)
        
        # Pressure vessels (power-law scaling)
        # Base: $50,000 for 10 m³
        if V_bed_total > 0:
            C['Pressure vessels'] = 50000 * (V_bed_total / 10) ** 0.6
        else:
            C['Pressure vessels'] = 0
        
        # Vacuum pump
        # Base: $30,000 for 100 kW
        if power_kW > 0:
            C['Vacuum pump'] = 30000 * (power_kW / 100) ** 0.7
        else:
            C['Vacuum pump'] = 0
        
        # Adsorbent (linear cost)
        C['Adsorbent'] = self.adsorbent_cost * m_ads_total
        
        # Piping and valves (15% of vessel cost)
        C['Piping and valves'] = 0.15 * C['Pressure vessels']


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
            ID=None,
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
        # Goal: use the least number of pumps while keeping per-stage head within
        # a practical centrifugal range to avoid pump-type warnings.
        # Typical industrial centrifugal pumps handle <= 8,000 ft head per stage.
        # This algorithm maximizes per-stage pressure rise under both head and
        # compression-ratio constraints.
        max_head_ft = 8000.0
        max_pressure_ratio = 4.0
        g = 9.80665  # m/s^2
        rho = max(feed.rho, 1e-6)  # kg/m3, prevent divide-by-zero
        max_deltaP_by_head = max_head_ft * 0.3048 * rho * g  # Pa

        total_power = 0.0
        current_stream = feed.copy()
        current_P = max(feed.P, 101325.0)
        target_P = self.P_high

        pressure_stages = []
        while current_P < target_P:
            max_deltaP_by_ratio = current_P * (max_pressure_ratio - 1.0)
            deltaP = min(max_deltaP_by_head, max_deltaP_by_ratio, target_P - current_P)
            next_P = current_P + max(deltaP, 1.0)
            pressure_stages.append(next_P)
            current_P = next_P

        for target_pressure in pressure_stages:
            temp_pump = bst.Pump(
                ID=None,
                ins=current_stream.copy(),
                P=target_pressure,
            )
            temp_pump.simulate()
            total_power += temp_pump.power_utility.rate
            current_stream = temp_pump.outs[0].copy()
        
        # Store total power consumption
        self.power_utility.rate = total_power
        
        # 2. Simulate the valve to find the thermal effect (e.g., cooling)
        temp_valve = bst.IsenthalpicValve(
            ID=None,
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


class Centrifuge(bst.SolidsCentrifuge):
        """
        Centrifuge with expanded solids-loading range for small-scale operation.

        Notes
        -----
        - BioSTEAM default bounds are 1-20 (reciprocating pusher) and 2-40 (scroll
            solid bowl) ton/hr. For pilot-scale operations (0.1-2 ton/hr), we allow
            extrapolation of the Humbird/Seider correlation to avoid over-warning
            while keeping the same cost exponent.
        - Cost scaling remains consistent with BioSTEAM (Humbird et al. 2011; Seider
            et al. 2017), but the lower bound is extended for small solids loading.
        """
        solids_loading_range = {
                'reciprocating_pusher': (0.1, 20),
                'scroll_solid_bowl': (0.1, 40),
        }

        #: Minimum solids loading used for cost extrapolation (ton/hr)
        min_solids_loading_cost = 0.1

        def _design(self):
                solids, centrifuge_type = self._solids, self.centrifuge_type
                ts = sum([s.imass[solids].sum() for s in self.ins if not s.isempty()])
                ts *= 0.0011023  # To short tons (2000 lbs/hr)
                self.design_results['Solids loading'] = ts
                lb, ub = self.solids_loading_range[centrifuge_type]
                if ts < lb:
                    lb_warning(self, 'Solids loading', ts, 'ton/hr', lb)
                self.design_results['Number of centrifuges'] = int(np.ceil(ts/ub)) if ub else 1

                # Extrapolate cost down to min_solids_loading_cost
                ts_cost = max(ts, self.min_solids_loading_cost)
                cost = 68040 * (ts_cost**0.5) if centrifuge_type else 170100 * (ts_cost**0.3)
                cost *= bst.CE / 567
                self.baseline_purchase_costs['Centrifuges'] = cost
                self.F_BM['Centrifuges'] = 2.03
                self.design_results['Flow rate'] = F_vol_in = self.F_vol_in
                self.power_utility(F_vol_in * self.kWhr_per_m3)

class ReverseOsmosis(bst.Unit):
    """
    Reverse Osmosis unit with literature-validated parameters.
    
    An upgraded version of ReverseOsmosis with validated default parameters
    based on industrial references and proper citations.
    
    Parameters
    ----------
    ins : 
        Inlet fluid to be treated.
    outs : 
        * [0] Permeate (treated water)
        * [1] Retentate (concentrated brine)
    water_recovery : float, optional
        Fraction of water recovered to permeate. Defaults to 0.85.
        Ref: Axeon Water, DuPont - typical industrial range 75-95%.
    membrane_flux : float, optional
        Membrane flux in L/m²/hr (LMH). Defaults to 40.
        Ref: membranes.com - typical wastewater RO range 17-40 LMH.
    membrane_cost_per_m2 : float, optional
        Membrane cost in USD/m². Defaults to 200.
        Ref: SNS Insider 2023, Samcotech - industrial range $80-$350/m².
    membrane_lifetime_years : float, optional
        Membrane replacement interval in years. Defaults to 3.
        Ref: DuPont/DOW membrane guidelines - typical 3-5 years.
    plant_lifetime_years : float, optional
        Plant operational lifetime in years. Defaults to 20.
        Ref: NREL TEA guidelines, Peters & Timmerhaus.
    operating_pressure_bar : float, optional
        Operating pressure in bar. Defaults to 25.
        Ref: Axeon Water - low-pressure RO 10-15 bar, standard 15-30 bar.
    specific_energy_consumption : float, optional
        Specific energy consumption in kWh/m³ permeate. Defaults to 3.0.
        Ref: VSEP, DuPont technical bulletins - typical range 2-6 kWh/m³.
    
    References
    ----------
    .. [1] Axeon Water Technologies. "Understanding RO Water Recovery Rates."
    .. [2] DuPont Water Solutions. "FilmTec Reverse Osmosis Technical Manual."
    .. [3] SNS Insider. "Reverse Osmosis Membrane Market Report 2023."
    .. [4] NREL. "Process Design and Economics for Biochemical Conversion." (2011)
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['H2O', 'NaCl'])
    >>> feed = bst.Stream('feed', H2O=1000, NaCl=10, units='kg/hr')
    >>> RO = RO2('RO1', ins=feed)
    >>> RO.simulate()
    >>> RO.show()
    """
    _N_ins = 1
    _N_outs = 2
    
    _F_BM_default = {
        'Pump': 2.3,                    # Ref: Peters & Timmerhaus, Table 16.11
        'Membrane replacement': 1.0,    # Direct replacement cost, no installation
        'Hardware': 1.8,                # (assumed) Pressure vessels, piping
    }

    _units = {
        'Flow rate': 'm3/hr',
        'Permeate flow': 'm3/hr',
        'Power': 'kW',
        'Pump Pressure': 'bar',
        'Water recovery': '%',
        'Membrane Area': 'm2',
        'Membrane flux': 'L/m2/hr',
        'Specific energy': 'kWh/m3',
    }

    def _init(self, 
              water_recovery=0.85,           # Ref: Axeon - typical 75-95%
              membrane_flux=40,               # Ref: membranes.com - 17-40 LMH
              membrane_cost_per_m2=200,       # Ref: SNS 2023 - $80-350/m²
              membrane_lifetime_years=3,      # Ref: DuPont - 3-5 years
              plant_lifetime_years=20,        # Ref: NREL TEA
              operating_pressure_bar=25,      # Ref: Axeon - 15-30 bar standard
              specific_energy_consumption=3.0 # Ref: typical 2-6 kWh/m³
              ):
        self.water_recovery = water_recovery
        self.membrane_flux = membrane_flux  # L/m²/hr
        self.membrane_cost_per_m2 = membrane_cost_per_m2  # USD/m²
        self.membrane_lifetime_years = membrane_lifetime_years
        self.plant_lifetime_years = plant_lifetime_years
        self.operating_pressure_bar = operating_pressure_bar
        self.specific_energy_consumption = specific_energy_consumption
    
    @property
    def permeate(self):
        """Return the permeate (treated water) stream."""
        return self.outs[0]
    
    @property
    def retentate(self):
        """Return the retentate (concentrated brine) stream."""
        return self.outs[1]
    
    def _run(self):
        feed, = self.ins
        permeate, retentate = self.outs
        
        # Copy thermal conditions
        permeate.copy_thermal_condition(feed)
        retentate.copy_like(feed)
        
        # Water balance based on recovery
        water_index = self.chemicals.index('H2O')
        feed_water_flow = feed.mol[water_index]
        water_recovered = self.water_recovery * feed_water_flow
        
        # Permeate is mostly pure water (high rejection of solutes)
        permeate.empty()
        permeate.mol[water_index] = water_recovered
        
        # Retentate contains remaining water and all solutes
        retentate.mol[water_index] = feed_water_flow - water_recovered
        
        # For non-water species, assume high rejection (99%+ stays in retentate)
        # This is characteristic of RO membranes
        for i, chem in enumerate(self.chemicals):
            if chem.ID != 'H2O':
                rejection = 0.99  # Ref: Typical RO rejection >99% for most ions
                permeate.mol[i] = feed.mol[i] * (1 - rejection)
                retentate.mol[i] = feed.mol[i] * rejection

    def _design(self):
        Design = self.design_results
        feed = self.ins[0]
        permeate = self.outs[0]
        
        # Basic flow measurements
        Design['Flow rate'] = feed.F_vol  # m³/hr
        Design['Permeate flow'] = permeate.F_vol  # m³/hr
        Design['Water recovery'] = self.water_recovery * 100  # %
        Design['Pump Pressure'] = self.operating_pressure_bar  # bar
        Design['Membrane flux'] = self.membrane_flux  # L/m²/hr
        Design['Specific energy'] = self.specific_energy_consumption  # kWh/m³

        # Calculate membrane area based on flux and permeate flow rate
        permeate_flow_L_hr = permeate.F_vol * 1000  # m³/hr to L/hr
        if self.membrane_flux > 0 and permeate_flow_L_hr > 0:
            membrane_area_m2 = permeate_flow_L_hr / self.membrane_flux
        else:
            membrane_area_m2 = 0.0
        Design['Membrane Area'] = membrane_area_m2  # m²

        # Calculate power two ways and use the higher value for conservatism:
        # Method 1: Based on specific energy consumption
        permeate_vol_m3_hr = permeate.F_vol
        power_from_SEC = self.specific_energy_consumption * permeate_vol_m3_hr  # kW
        
        # Method 2: Based on pump simulation
        feed_copy = bst.Stream()
        feed_copy.copy_like(feed)
        self.pump = bst.Pump(
            ins=feed_copy,
            outs=(),
            P=self.operating_pressure_bar * 1e5  # bar to Pa
        )
        self.pump.simulate()
        self.pump._design()
        power_from_pump = self.pump.power_utility.rate
        
        # Use higher of the two estimates
        power_kW = max(power_from_SEC, power_from_pump)
        Design['Power'] = power_kW  # kW
        self.power_utility.rate = power_kW

    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        
        # 1. Pump cost
        if hasattr(self, 'pump') and self.pump:
            self.pump._cost()
            C['Pump'] = self.pump.purchase_cost
        else:
            C['Pump'] = 0.0
        
        # 2. Membrane replacement cost over plant lifetime
        membrane_area = D.get('Membrane Area', 0.0)
        if membrane_area > 0 and self.membrane_lifetime_years > 0:
            num_replacements = self.plant_lifetime_years / self.membrane_lifetime_years
            total_membrane_cost = membrane_area * self.membrane_cost_per_m2 * num_replacements
            C['Membrane replacement'] = total_membrane_cost
        else:
            C['Membrane replacement'] = 0.0
        
        # 3. Hardware cost (pressure vessels, piping, instrumentation)
        # Ref: Power-law scaling based on membrane area
        # Base: $50,000 for 100 m² membrane area, exponent 0.7
        if membrane_area > 0:
            base_cost = 50000  # USD
            base_area = 100    # m²
            exponent = 0.7     # (assumed) typical for process equipment
            C['Hardware'] = base_cost * (membrane_area / base_area) ** exponent
        else:
            C['Hardware'] = 0.0


class Diafiltration(bst.Unit):
    """
    Verified Diafiltration unit for separation of solutes based on size.
    
    Validated against industrial literature for protein purification applications.
    Detailed Specification: docs/units/Diafiltration_spec.md
    
    Parameters
    ----------
    TargetProduct_Retention : float
        Retention of target product (0-1). Default: 0.99.
        Ref: Protein-specific, typical >95% for UF with 10kDa MWCO.
    membrane_flux_LMH : float
        Permeate flux [L/m²/hr]. Default: 40.0.
        Ref: Membranes.com, Synder Filtration.
        Valid range: UF 20-80 LMH, NF 10-40 LMH [4].
    membrane_cost_USD_per_m2 : float
        Replacement cost. Default: 150.0.
        Ref: SNS Insider 2023 [5].
        Valid range: $80-$350/m² (industrial), $1500-$5000/m² (pharma).
    membrane_lifetime_years : float
        Frequency of replacement. Default: 2.0.
        Ref: DuPont Water Solutions [3].
        Valid range: 1-3 years (process dependent).
    TMP_bar1, TMP_bar2 : float
        Transmembrane pressure. Default: 2.0 bar.
        Ref: Synder Filtration [4].
        Valid range: UF 2-4 bar, NF 10-25 bar.
    recirculation_ratio : float
        Recirculation to reduce fouling. Default: 10.0.
        Valid range: 5-20 (application dependent).
    
    Notes
    -----
    For LegHb purification, UF with 3-10 kDa MWCO is typical since LegHb monomer
    is ~16 kDa [1]. Diafiltration with 5-7 diavolumes achieves good buffer exchange
    and impurity removal.
    """
    _N_ins = 2
    _N_outs = 2

    # =====================================================================
    # PRESETS: Factory configurations for different membrane types
    # =====================================================================
    PRESETS = {
        'UF': {  # Ultrafiltration - for protein concentration
            'membrane_flux_LMH': 50.0,        # LMH [appliedmembranes.com]
            'TMP_bar1': 2.0,                  # bar [moruiwater.com]
            'TMP_bar2': 1.5,                  # bar
            'membrane_cost_USD_per_m2': 150.0,# $/m² [industry avg]
            'membrane_lifetime_years': 3.0,   # years [jx-purification.com]
            'TargetProduct_Retention': 0.99,  # High protein retention
            'Salt_Retention': 0.05,           # Salts pass through
            'OtherLargeMolecules_Retention': 0.98,
            'DefaultSolutes_Retention': 0.08,
        },
        'NF': {  # Nanofiltration - for buffer exchange / fine separation
            'membrane_flux_LMH': 25.0,        # LMH [lower flux]
            'TMP_bar1': 8.0,                  # bar [dupont.com]
            'TMP_bar2': 5.0,                  # bar
            'membrane_cost_USD_per_m2': 250.0,# $/m² [higher cost]
            'membrane_lifetime_years': 2.0,   # years [more fouling]
            'TargetProduct_Retention': 0.995, # Very high retention
            'Salt_Retention': 0.30,           # Partial salt rejection
            'OtherLargeMolecules_Retention': 0.995,
            'DefaultSolutes_Retention': 0.15,
        },
    }

    @classmethod
    def from_preset(cls, preset, ID='', ins=None, outs=None, thermo=None, **kwargs):
        """
        Create a Diafiltration unit from a named preset.
        
        Parameters
        ----------
        preset : str
            Preset name: 'UF' (Ultrafiltration) or 'NF' (Nanofiltration).
        ID : str
            Unit ID.
        ins, outs : streams
            Input and output streams.
        thermo : Thermo, optional
            Thermodynamic property package.
        **kwargs : 
            Additional parameters to override preset values.
        
        Returns
        -------
        Diafiltration
            Configured unit instance.
        
        Examples
        --------
        >>> U401 = Diafiltration.from_preset('UF', ID='U401', ins=(...), outs=(...))
        >>> U403 = Diafiltration.from_preset('NF', ID='U403', ins=(...), outs=(...))
        """
        if preset not in cls.PRESETS:
            raise ValueError(f"Unknown preset '{preset}'. Choose from: {list(cls.PRESETS.keys())}")
        params = cls.PRESETS[preset].copy()
        params.update(kwargs)
        unit = cls(ID, ins, outs, thermo, **params)
        unit.preset = preset
        return unit


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

import flexsolve as flx
from biosteam import main_flowsheet as F


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


# %%
####################################
##### Verified Unit Upgrades #######
####################################



class Filtration(bst.Unit):
    """
    Continuous Rotary Drum Vacuum Filter (RDVF) for biomass separation.
    
    Replaces generic Centrifuge models with a rigorous filtration model based on
    surface area loading and cake formation dynamics.
    Detailed Specification: docs/units/Filtration_spec.md
    
    Parameters
    ----------
    solids_loading : float
        Specific solids throughput [kg-solids / m² / hr]. 
        Default: 20.0 (Typical for fermentation broth/mycelia).
        Ref: Perry's Handbook 8th Ed, Table 18-5 [2].
        Valid range: 10-50 kg/m²/hr (yeast/mycelia filtration).
    cake_moisture_content : float
        Residual moisture in discharged cake [wt fraction water].
        Default: 0.20 (20%).
        Valid range: 0.20-0.80 depending on cake compressibility.
    solid_capture_efficiency : float
        Fraction of solids retained in the cake.
        Default: 0.99.
        Valid range: 0.95-0.999 (typical industrial).
    power_per_m2 : float
        Specific power consumption (vacuum pump + drive) [kW/m²].
        Default: 1.0.
        Ref: Perry's Handbook [2].
        Valid range: 0.5-2.0 kW/m².
    TMP_bar : float
        Transmembrane/vacuum pressure differential [bar].
        Default: 2.0.
        Valid range: 1.5-3.0 bar for RDVF.
    membrane_cost_USD_per_m2 : float
        Filter cloth/membrane cost. Default: 100.0.
        Valid range: $50-$200/m² for industrial cloth.
    membrane_lifetime_years : float
        Filter medium replacement interval. Default: 3.0.
        Valid range: 2-5 years.
    
    Notes
    -----
    For P. pastoris/K. marxianus biomass, RDVF is suitable for primary
    clarification. Cell disruption releases intracellular products (LegHb)
    and creates smaller debris particles requiring downstream polishing [1].
    """
    _N_ins = 1
    _N_outs = 2
    
    # =====================================================================
    # PRESETS: Factory configurations for different membrane types
    # =====================================================================
    PRESETS = {
        'MF': {  # Microfiltration - for cell separation / primary clarification
            'solids_loading': 30.0,             # kg/m²/hr [higher throughput]
            'TMP_bar': 1.5,                     # bar [low pressure]
            'power_per_m2': 0.8,                # kW/m² [lower power]
            'membrane_cost_USD_per_m2': 80.0,   # $/m² [cheaper]
            'membrane_lifetime_years': 4.0,     # years [longer life]
            'cake_moisture_content': 0.25,      # wetter cake
            'solid_capture_efficiency': 0.98,   # good capture
        },
        'UF': {  # Ultrafiltration - for finer separation / virus removal
            'solids_loading': 15.0,             # kg/m²/hr [lower throughput]
            'TMP_bar': 3.0,                     # bar [higher pressure]
            'power_per_m2': 1.2,                # kW/m² [higher power]
            'membrane_cost_USD_per_m2': 150.0,  # $/m² [more expensive]
            'membrane_lifetime_years': 2.5,     # years [shorter life]
            'cake_moisture_content': 0.20,      # drier cake
            'solid_capture_efficiency': 0.995,  # better capture
        },
    }

    @classmethod
    def from_preset(cls, preset, ID='', ins=None, outs=(), thermo=None, **kwargs):
        """
        Create a Filtration unit from a named preset.
        
        Parameters
        ----------
        preset : str
            Preset name: 'MF' (Microfiltration) or 'UF' (Ultrafiltration).
        ID : str
            Unit ID.
        ins, outs : streams
            Input and output streams.
        thermo : Thermo, optional
            Thermodynamic property package.
        **kwargs : 
            Additional parameters to override preset values.
        
        Returns
        -------
        Filtration
            Configured unit instance.
        
        Examples
        --------
        >>> S402 = Filtration.from_preset('MF', ID='S402', ins=(...), outs=(...))
        >>> S403 = Filtration.from_preset('UF', ID='S403', ins=(...), outs=(...))
        """
        if preset not in cls.PRESETS:
            raise ValueError(f"Unknown preset '{preset}'. Choose from: {list(cls.PRESETS.keys())}")
        params = cls.PRESETS[preset].copy()
        params.update(kwargs)
        unit = cls(ID, ins, outs, thermo, **params)
        unit.preset = preset
        return unit
    
    _units = {
        'Filter Area': 'm2',
        'Power': 'kW',
        'Cost': 'USD'
    }
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 solids_loading=20.0,
                 cake_moisture_content=0.20,
                 solid_capture_efficiency=0.99,
                 power_per_m2=1.0,
                 TMP_bar=2.0,
                 membrane_cost_USD_per_m2=100.0,
                 membrane_lifetime_years=3.0):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.solids_loading = solids_loading
        self.cake_moisture_content = cake_moisture_content
        self.solid_capture_efficiency = solid_capture_efficiency
        self.power_per_m2 = power_per_m2
        self.TMP_bar = TMP_bar
        self.membrane_cost_USD_per_m2 = membrane_cost_USD_per_m2
        self.membrane_lifetime_years = membrane_lifetime_years
        self.preset = None  # Will be set by from_preset
        
    def _run(self):
        feed = self.ins[0]
        cake, filtrate = self.outs
        
        # Initialize
        cake.empty()
        filtrate.empty()
        cake.T = filtrate.T = feed.T
        cake.P = filtrate.P = feed.P
        
        # Calculate total solids in feed
        solids_mass = 0.0
        
        # Standard splitting logic
        for chem in self.chemicals:
            mass_in = feed.imass[chem.ID]
            if mass_in <= 0: continue
            
            # Identify solids: Explicit list + _In suffix checks
            is_solid = False
            if chem.ID in ('cellmass', 'Biomass', 'Pichia_pastoris', 'Cellulose', 'Glucan', 'Lignin'):
                is_solid = True
            elif chem.ID.endswith('_In'):
                # Intracellular components trapped in cells are effectively solids
                is_solid = True
            elif chem.ID in ('Precipitate',): 
                is_solid = True
                
            if is_solid:
                captured = mass_in * self.solid_capture_efficiency
                cake.imass[chem.ID] += captured
                filtrate.imass[chem.ID] += mass_in - captured
                solids_mass += captured
            else:
                # Water and solubles go to filtrate initially
                filtrate.imass[chem.ID] += mass_in
                
        # Adjust cake moisture
        # Cake total mass = Solids + Water
        # Water / (Solids + Water) = moisture_content
        # Water = Solids * moisture / (1 - moisture)
        if self.cake_moisture_content < 1.0:
            req_water = solids_mass * self.cake_moisture_content / (1.0 - self.cake_moisture_content)
        else:
            req_water = solids_mass # Error case
            
        # Move required water from filtrate to cake
        # Only move water (H2O)
        available_water = filtrate.imass['H2O']
        water_transfer = min(available_water, req_water)
        
        filtrate.imass['H2O'] -= water_transfer
        cake.imass['H2O'] += water_transfer
        
        self.design_results['Solids Capture'] = solids_mass
        
    def _design(self):
        D = self.design_results
        
        # Calculate Area
        solids_rate = D.get('Solids Capture', 0.0) # kg/hr
        if self.solids_loading > 0:
            area = solids_rate / self.solids_loading
        else:
            area = 0.0
            
        D['Filter Area'] = area
        
        # Calculate Power
        power = area * self.power_per_m2
        D['Power'] = power
        
        self.power_utility.rate = power
        
    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        area = D.get('Filter Area', 0.0)
        
        # Cost correlation for Rotary Vacuum Filter
        # Source: Peters & Timmerhaus (general equipment)
        # Base cost ~$65,000 for 10 m2
        if area > 0:
            C['Rotary Vacuum Filter'] = 65000 * (area / 10)**0.6
        else:
            C['Rotary Vacuum Filter'] = 0.0

# %%
##########################
##### New Units ##########
##########################

class ResinColumn(bst.Unit):
    """
    Polymorphic unit supporting both Ion Exchange (chromatography) and Adsorption 
    (e.g., activated carbon, resin adsorption) modes.
    
    Standardizes inputs/outputs for both operations:
    - Ion Exchange: Binds target, washes, elutes.
    - Adsorption: Flow-through cleaning or capture.

    Detailed Specification: docs/units/ResinColumn2_spec.md
    
    Parameters
    ----------
    preset : str
        'IonExchange' or 'Adsorption'. Determines the physical model used.
    resin_DBC_g_L : float (IonExchange mode)
        Dynamic Binding Capacity [g protein/L resin]. Default: 50.0.
        Ref: Typical industrial values 30-100 g/L [6].
        Valid range: 30-100 g/L (resin and target dependent).
    cycle_time_hr : float
        Total cycle duration. Default: 4.0 hr.
        Valid range: 2-8 hr (throughput vs. binding optimization).
    resin_cost_USD_per_L : float
        Resin cost. Default: 30.0.
        Valid range: Strong cation $10-50/L, Strong anion $30-100/L.
    resin_lifetime_years : float
        Resin replacement interval. Default: 5.0.
        Ref: Vendor data [6].
        Valid range: 1-10 years.
    EBCT_min : float (Adsorption mode)
        Empty Bed Contact Time [min]. Default: 5.0.
        Ref: urbansaqua.com.
        Valid range: 2-10 min for liquid treatment.
    superficial_velocity_m_h : float (Adsorption mode)
        Liquid flow velocity [m/hr]. Default: 10.0.
        Ref: aquaenergyexpo.com.
        Valid range: 5-20 m/hr.
    adsorbent_bulk_density : float (Adsorption mode)
        Adsorbent packing density [kg/m³]. Default: 450.
        Ref: calgoncarbon.com (GAC: 350-550 kg/m³).
    adsorbent_cost_USD_per_kg : float (Adsorption mode)
        Adsorbent cost. Default: 5.0.
        Valid range: Activated carbon ~$2-8/kg.
    
    Notes
    -----
    For LegHb purification, ion exchange chromatography can be used to remove
    bound impurities like Heme_b while retaining the target protein. See Walch &
    Jungbauer (2017) for continuous desalting approaches [6].
    """
    _N_ins = 4
    _N_outs = 4
    
    _F_BM_default = {
        'Column Hardware': 2.5,
        'Resin': 1.0,
        'Adsorbent': 1.0,
        'Pump': 2.3,
    }

    PRESETS = {
        'IonExchange': {
            'mode': 'IonExchange',
            'resin_DBC_g_L': 50.0,
            'cycle_time_hr': 4.0,
            'regeneration_CV': 5.0,
            'resin_cost_USD_per_L': 30.0,
            'resin_lifetime_years': 5.0,
        },
        'Adsorption': {
            'mode': 'Adsorption',
            'EBCT_min': 5.0, # Ref: urbansaqua.com (5-30 min)
            'superficial_velocity_m_h': 10.0, # Ref: aquaenergyexpo.com (5-20 m/h)
            'adsorbent_bulk_density': 450.0, # Ref: calgoncarbon.com (350-550 kg/m3)
            'adsorbent_cost_USD_per_kg': 5.0,
            'adsorbent_lifetime_years': 3.0,
        }
    }
    
    @classmethod
    def from_preset(cls, preset, ID='', ins=None, outs=(), thermo=None, **kwargs):
        """
        Create a ResinColumn2 unit from a named preset.
        """
        if preset not in cls.PRESETS:
             # Fallback or strict? Let's check if it matches a mode key directly
             pass
        
        # Start with defaults from the preset
        if preset in cls.PRESETS:
            defaults = cls.PRESETS[preset].copy()
            # Remove 'mode' from defaults if it clashes or use it to set preset
            mode_val = defaults.pop('mode', preset)
            # Update with kwargs
            defaults.update(kwargs)
            return cls(ID, ins, outs, thermo, preset=mode_val, **defaults)
        else:
            # If not a standard preset, assume user knows what they are doing or it's a variance
            return cls(ID, ins, outs, thermo, preset=preset, **kwargs)

    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 preset='IonExchange',
                 # Common
                 TargetProduct_IDs=('Leghemoglobin',),
                 TargetProduct_Yield=0.95,
                 BoundImpurity_IDs=('Heme_b',),
                 BoundImpurity_Removal=0.93,
                 NonBinding_Carryover=0.04,
                 
                 # Ion Exchange Specific
                 cycle_time_hr=4.0,
                 equilibration_CV=5.0,
                 wash_CV=5.0,
                 elution_CV=3.0,
                 regeneration_CV=5.0,
                 resin_DBC_g_L=50.0,
                 load_safety_factor=0.8,
                 resin_cost_USD_per_L=30.0,
                 resin_lifetime_years=5.0,
                 column_hardware_cost_factor=30000.0,
                 column_hardware_cost_exponent=0.6,
                 
                 # Adsorption Specific
                 EBCT_min=5.0, # Ref: urbansaqua.com
                 superficial_velocity_m_h=10.0, # Ref: aquaenergyexpo.com
                 adsorbent_bulk_density=450.0, # Ref: calgoncarbon.com
                 adsorbent_cost_USD_per_kg=5.0,
                 adsorbent_lifetime_years=3.0,
                 # Adsorption separation parameters (realistic impurity distribution)
                 NonTarget_Removal=0.99,  # 99% of non-targets to flowthrough
                 Wash_Impurity_Carryover=0.02,  # 2% of organics leak to wash
                 Regen_Impurity_Carryover=0.01,  # 1% of organics leak to regen
                 ):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.preset = preset
        
        # Separation
        self.TargetProduct_IDs = tuple(TargetProduct_IDs) if not isinstance(TargetProduct_IDs, str) else (TargetProduct_IDs,)
        self.TargetProduct_Yield = TargetProduct_Yield
        self.BoundImpurity_IDs = tuple(BoundImpurity_IDs)
        self.BoundImpurity_Removal = BoundImpurity_Removal
        self.NonBinding_Carryover = NonBinding_Carryover
        
        # IEX Params
        self.cycle_time_hr = cycle_time_hr
        self.equilibration_CV = equilibration_CV
        self.wash_CV = wash_CV
        self.elution_CV = elution_CV
        self.regeneration_CV = regeneration_CV
        self.resin_DBC_g_L = resin_DBC_g_L
        self.load_safety_factor = load_safety_factor
        self.resin_cost_USD_per_L = resin_cost_USD_per_L
        self.resin_lifetime_years = resin_lifetime_years
        self.column_hardware_cost_factor = column_hardware_cost_factor
        self.column_hardware_cost_exponent = column_hardware_cost_exponent
        
        # Adsorption Params
        self.EBCT_min = EBCT_min
        self.superficial_velocity_m_h = superficial_velocity_m_h
        self.adsorbent_bulk_density = adsorbent_bulk_density
        self.adsorbent_cost_USD_per_kg = adsorbent_cost_USD_per_kg
        self.adsorbent_lifetime_years = adsorbent_lifetime_years
        # Adsorption separation parameters
        self.NonTarget_Removal = NonTarget_Removal
        self.Wash_Impurity_Carryover = Wash_Impurity_Carryover
        self.Regen_Impurity_Carryover = Regen_Impurity_Carryover
        
        # Components
        # Initialize with dummy streams to avoid MissingStream errors during costing
        self.pump = bst.Pump(None, ins=bst.Stream(None), outs=bst.Stream(None), P=4e5) 
        self._auxiliary_unit_names = ('pump',)

    def _run(self):
        if self.preset == 'IonExchange':
            self._run_ion_exchange()
        elif self.preset == 'Adsorption':
            self._run_adsorption()
        else:
            # Fallback to Ion Exchange logic if unknown? Or Error.
            # Let's try to infer or error.
            # If parameters look like IEX, maybe allow? But plan says safer to raise.
            raise ValueError(f"Unsupported preset mode: {self.preset}. Choose 'IonExchange' or 'Adsorption'.")

    def _run_ion_exchange(self):
        # Implementation copied and adapted from IonExchangeCycle
        feed, buffer_A, buffer_B, regen_sol = self.ins
        product, ft_waste, wash_waste, regen_waste = self.outs
        
        # Setup outputs
        product.copy_like(buffer_B)
        ft_waste.copy_like(feed)
        wash_waste.copy_like(buffer_A)
        regen_waste.copy_like(regen_sol)
        
        target_ids = set(self.TargetProduct_IDs)
        bound_impurity_ids = set(self.BoundImpurity_IDs)
        
        # Solute balance
        for chem in self.chemicals:
            solute_in_feed = feed.imass[chem.ID]
            if solute_in_feed < 1e-12: continue
            
            if chem.ID in target_ids:
                to_product = solute_in_feed * self.TargetProduct_Yield
                to_regen = solute_in_feed - to_product
                product.imass[chem.ID] += to_product
                regen_waste.imass[chem.ID] += to_regen
                ft_waste.imass[chem.ID] = 0
            elif chem.ID in bound_impurity_ids:
                to_regen = solute_in_feed * self.BoundImpurity_Removal
                to_product = solute_in_feed - to_regen
                product.imass[chem.ID] += to_product
                regen_waste.imass[chem.ID] += to_regen
                ft_waste.imass[chem.ID] = 0
            else:
                carryover = solute_in_feed * self.NonBinding_Carryover
                product.imass[chem.ID] += carryover
                ft_waste.imass[chem.ID] -= carryover

        # Temperatures
        product.T = buffer_B.T
        ft_waste.T = feed.T
        wash_waste.T = buffer_A.T
        regen_waste.T = regen_sol.T

    def _run_adsorption(self):
        """
        Adsorption Mode (Hydrophobic Capture + Elution):
        
        Process Description:
        - Feed passes through resin; hydrophobic targets bind (adsorb).
        - Hydrophilic/non-binding solutes (salts, sugars) mostly flow through to waste.
        - Small amounts of organics leak to wash and regeneration streams (realistic behavior).
        - Elution releases adsorbed targets into the eluate (product stream).
        - Regeneration is a separate cleaning step.
        
        Stream Mapping:
        - ins[0] (Feed)  -> outs[0] (Flowthrough / Waste)
        - ins[1] (Wash)  -> outs[3] (ResinWash)
        - ins[2] (Elute) -> outs[1] (Eluate / Product)
        - ins[3] (Regen) -> outs[2] (Regen waste)
        
        Mass Balance (updated for realistic impurity distribution):
        - Flowthrough = Feed * NonTarget_Removal (for non-targets)
        - ResinWash = Wash_In + Feed * Wash_Impurity_Carryover
        - Eluate = Elute_In + Adsorbed_Targets + trace impurities
        - RegenWaste = Regen_In + Feed * Regen_Impurity_Carryover
        """
        feed = self.ins[0]
        wash_in = self.ins[1]
        elute_in = self.ins[2]
        regen_in = self.ins[3]
        
        treated = self.outs[0]      # Flowthrough / Waste
        eluate = self.outs[1]       # Elution stream / Product
        regen_out = self.outs[2]    # Regeneration waste (separate)
        wash_out = self.outs[3]     # Spent Wash Buffer
        
        # 1. Initialize streams
        treated.copy_like(feed)
        eluate.copy_like(elute_in)
        regen_out.copy_like(regen_in)
        wash_out.copy_like(wash_in)
        
        treated.T = feed.T
        treated.P = feed.P
        eluate.T = elute_in.T
        eluate.P = elute_in.P
        regen_out.T = regen_in.T
        regen_out.P = regen_in.P
        wash_out.T = wash_in.T
        wash_out.P = wash_in.P
        
        # 2. Perform Adsorption (Separation)
        target_ids = set(self.TargetProduct_IDs)
        
        for chem in self.chemicals:
            mass_in = feed.imass[chem.ID]
            if mass_in < 1e-12:
                continue
            
            if chem.ID in target_ids:
                # Target is adsorbed, then eluted
                removed = mass_in * self.TargetProduct_Yield
                remaining = mass_in - removed
                treated.imass[chem.ID] = remaining
                eluate.imass[chem.ID] += removed
            else:
                # Non-targets: Distribute across streams (realistic behavior like IonExchange)
                # Most go to flowthrough, small amounts leak to wash/regen/eluate
                to_flowthrough = mass_in * self.NonTarget_Removal
                to_wash = mass_in * self.Wash_Impurity_Carryover
                to_regen = mass_in * self.Regen_Impurity_Carryover
                # Remaining trace goes to eluate (ensures mass balance closure)
                to_eluate = mass_in - to_flowthrough - to_wash - to_regen
                
                treated.imass[chem.ID] = to_flowthrough
                wash_out.imass[chem.ID] += max(0, to_wash)
                regen_out.imass[chem.ID] += max(0, to_regen)
                eluate.imass[chem.ID] += max(0, to_eluate)
                
    def _design(self):

        if self.preset == 'IonExchange':
            self._design_ion_exchange()
        elif self.preset == 'Adsorption':
            self._design_adsorption()

    def _design_ion_exchange(self):
        # Existing logic
        D = self.design_results
        target_mass_kg = sum(self.ins[0].imass[ID] for ID in self.TargetProduct_IDs) * self.cycle_time_hr
        effective_DBC = self.resin_DBC_g_L * self.load_safety_factor
        
        if effective_DBC > 0 and target_mass_kg > 0:
            resin_vol = (target_mass_kg * 1000) / effective_DBC
        else:
            resin_vol = 1.0
            
        D['Resin Volume (L)'] = resin_vol
        D['Cycle time (hr)'] = self.cycle_time_hr
        self._cost_pump(sum(s.F_mass for s in self.ins))

    def _design_adsorption(self):
        # EBCT design
        # V_bed = Q_vol * EBCT
        D = self.design_results
        feed = self.ins[0]
        Q_m3_h = feed.F_vol # m3/hr
        
        # Bed Volume
        bed_vol_m3 = Q_m3_h * (self.EBCT_min / 60.0)
        
        # Check Velocity constraint
        # Area = Q / v
        if self.superficial_velocity_m_h > 0:
            min_area_m2 = Q_m3_h / self.superficial_velocity_m_h
            # Height = V / A
            height_m = bed_vol_m3 / min_area_m2 if min_area_m2 > 0 else 0
        else:
            min_area_m2 = 0
            height_m = 0
            
        D['Bed Volume (m3)'] = bed_vol_m3
        D['Adsorbent Mass (kg)'] = bed_vol_m3 * self.adsorbent_bulk_density
        D['Target Area (m2)'] = min_area_m2
        D['Bed Height (m)'] = height_m
        
        self._cost_pump(feed.F_mass) # Adsorption usually just feed pump

    def _cost_pump(self, flow_mass):
        # Common pump sizing
        if flow_mass > 0 and self.ins[0].F_mass > 0:
            # Use feed properties as approximation for pump sizing
            # This ensures F_mass assignment works because composition is defined
            self.pump.ins[0].copy_like(self.ins[0])
            self.pump.ins[0].F_mass = flow_mass
            self.pump.simulate()
            self.pump._design()
            self.pump._cost()
            self.power_utility.rate = self.pump.power_utility.rate
        
    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        
        if self.preset == 'IonExchange':
            # Hardware
            vol_L = D.get('Resin Volume (L)', 0)
            if vol_L > 0:
                C['Column Hardware'] = bst.CE/500 * self.column_hardware_cost_factor * (vol_L ** self.column_hardware_cost_exponent)
                
                # Resin
                if self.resin_lifetime_years > 0:
                    life = getattr(F.stream, 'plant_life', 20)
                    replacements = life / self.resin_lifetime_years
                    C['Resin'] = self.resin_cost_USD_per_L * vol_L * replacements
            
        elif self.preset == 'Adsorption':
            # Adsorbent
            mass_kg = D.get('Adsorbent Mass (kg)', 0)
            if mass_kg > 0:
                 # Initial + Replacements
                 life = getattr(F.stream, 'plant_life', 20)
                 if self.adsorbent_lifetime_years > 0:
                     replacements = life / self.adsorbent_lifetime_years
                     total_mass = mass_kg * replacements
                     C['Adsorbent'] = total_mass * self.adsorbent_cost_USD_per_kg
            
            # Hardware (Vessel) - can approximate as pressure vessel or tank
            vol_m3 = D.get('Bed Volume (m3)', 0)
            if vol_m3 > 0:
                 # Using similar hardware cost logic but converting m3 to L if using same factor
                 # Or use specific adsorption vessel cost.
                 # Let's Assume similar complexity to IEX column for now, scaling by volume L
                 vol_L = vol_m3 * 1000
                 C['Column Hardware'] = bst.CE/500 * self.column_hardware_cost_factor * (vol_L ** self.column_hardware_cost_exponent)

        C['Pump'] = self.pump.purchase_cost





