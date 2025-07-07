
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040
"""


# %% Setup

import numpy as np
import thermosteam as tmo
from math import exp, pi, log, ceil
from flexsolve import aitken_secant
from biosteam import Unit, BatchCrystallizer
from biosteam.units import Flash, HXutility, Mixer, MixTank, Pump, \
    SolidsSeparator, StorageTank, LiquidsSplitSettler, \
    BatchBioreactor, StirredTankReactor, LiquidsMixingTank
from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year as CEPCI
from thermosteam import Stream, MultiStream
from biorefineries.TAL.process_settings import price
# from biorefineries.TAL.utils import CEPCI, baseline_feedflow, compute_extra_chemical, adjust_recycle
from biorefineries.TAL.utils import baseline_feedflow, compute_extra_chemical, adjust_recycle

from thermosteam import SeriesReaction
from biorefineries.TAL._general_utils import get_pH_polyprotic_acid_mixture, get_molarity
from flexsolve import IQ_interpolation

_kg_per_ton = 907.18474
_Gcal_2_kJ = 4.184 * 1e6 # (also MMkcal/hr)
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

compute_TAL_titer = lambda effluent: effluent.imass['TAL'] / effluent.F_vol

compute_TAL_mass = lambda effluent: effluent.imass['TAL']

#%% pH utils
def get_pH_stream(stream):
    if stream.imol['CitricAcid', 'H3PO4', 'AceticAcid'].sum() > 0.:
        return get_pH_polyprotic_acid_mixture(stream,
                                ['CitricAcid', 'H3PO4', 'AceticAcid'], 
                                [[10**-3.13, 10**-4.76, 10**-6.40], 
                                 [10**-2.16, 10**-7.21, 10**-12.32],
                                 [10**-4.76]],
                                'ideal')
    else:
        molarity_NaOH = get_molarity('NaOH', stream)
        if molarity_NaOH == 0.: return 7.
        else: return 14. + log(molarity_NaOH, 10.) # assume strong base completely dissociates in aqueous solution


def get_pH_given_base_addition(mol_base_per_m3_broth, base_mixer):
    base_mixer.mol_base_per_m3_broth = mol_base_per_m3_broth
    base_mixer.simulate_base_addition_and_acids_neutralization()
    return get_pH_stream(base_mixer.outs[0])

def load_pH(pH, base_mixer):
    obj_f_pH = lambda mol_base_per_m3_broth: get_pH_given_base_addition(mol_base_per_m3_broth=mol_base_per_m3_broth, 
                                                                        base_mixer=base_mixer)\
                                        - pH
    IQ_interpolation(obj_f_pH, 0., 5., ytol=0.001)
    
#%% Reactor
from biosteam.units.design_tools import PressureVessel
from biosteam.exceptions import DesignError
class Reactor(Unit, PressureVessel, isabstract=True):
    '''    
    Create an abstract class for reactor unit, purchase cost of the reactor
    is based on volume calculated by residence time.

    Parameters
    ----------
    ins : stream
        Inlet.        
    outs : stream
        Outlet.
    tau=0.5 : float
        Residence time [hr].        
    V_wf=0.8 : float
        Fraction of working volume over total volume.        
    kW_per_m3=0.985: float
        Power usage of agitator
        (0.985 converted from 5 hp/1000 gal as in [1], for liquid–liquid reaction or extraction).
    wall_thickness_factor=1: float
        A safety factor to scale up the calculated minimum wall thickness.
    vessel_material : str, optional
        Vessel material. Default to 'Stainless steel 316'.
    vessel_type : str, optional
        Vessel type. Can only be 'Horizontal' or 'Vertical'.
        
    References
    ----------
    .. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.; 
        Ng, M. K. Cost Accounting and Capital Cost Estimation. In Product 
        and Process Design Principles; Wiley, 2017; pp 470.
    '''
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    _units = {**PressureVessel._units,
              'Residence time': 'hr',
              'Total volume': 'm3',
              'Reactor volume': 'm3'}
    
    # For a single reactor, based on diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _V_max = pi/4*(20**2)*40/35.3147 
    
    def _init(self, 
                 #  *, 
                  P=101325, tau=0.5, V_wf=0.8,
                  length_to_diameter=2, 
                  kW_per_m3=0.985,
                  wall_thickness_factor=1,
                  vessel_material='Stainless steel 316',
                  vessel_type='Vertical'):
        
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    def _design(self):
        Design = self.design_results
        ins_F_vol = self.F_vol_in
        V_total = ins_F_vol * self.tau / self.V_wf
        P = self.P * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        wall_thickness_factor = self.wall_thickness_factor
        
        N = ceil(V_total/self._V_max)
        if N == 0:
            V_reactor = 0
            D = 0
            L = 0
        else:
            V_reactor = V_total / N
            D = (4*V_reactor/pi/length_to_diameter)**(1/3)
            D *= 3.28084 # convert from m to ft
            L = D * length_to_diameter

        Design['Residence time'] = self.tau
        Design['Total volume'] = V_total
        Design['Single reactor volume'] = V_reactor
        Design['Number of reactors'] = N
        Design.update(self._vessel_design(P, D, L))
        if wall_thickness_factor == 1: pass
        elif wall_thickness_factor < 1:
            raise DesignError('wall_thickness_factor must be larger than 1')
        else:
              Design['Wall thickness'] *= wall_thickness_factor
              # Weight is proportional to wall thickness in PressureVessel design
              Design['Weight'] = round(Design['Weight']*wall_thickness_factor,2)
            
    def _cost(self):
        Design = self.design_results
        purchase_costs = self.purchase_costs
        
        if Design['Total volume'] == 0:
            for i, j in purchase_costs.items():
                purchase_costs[i] = 0
        
        else:
            purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in purchase_costs.items():
                purchase_costs[i] *= Design['Number of reactors']
            
            self.power_utility(self.kW_per_m3 * Design['Total volume'])
    # def _run(self):
    #     PressureVessel._run()
    @property
    def BM(self):
        vessel_type = self.vessel_type
        if not vessel_type:
            raise AttributeError('vessel_type not defined')
        elif vessel_type == 'Vertical':
            return self.BM_vertical
        elif vessel_type == 'Horizontal':
            return self.BM_horizontal 
        else:
            raise RuntimeError("invalid vessel type")
    
    @property
    def F_vol_in(self): # exclude gases
        return sum([i.F_vol for i in self.ins if i.phase=='l' and i.F_mol and not (i.imol['CO2']/i.F_mol==1. 
                                                                       or i.imol['O2']/i.F_mol>0.1
                                                                       or i.imol['H2']/i.F_mol>0.1)])
    
    @property
    def F_vol_out(self): # exclude gases
        return sum([i.F_vol for i in self.outs if i.phase=='l' and i.F_mol and not (i.imol['CO2']/i.F_mol==1. 
                                                                       or i.imol['O2']/i.F_mol>0.1
                                                                       or i.imol['H2']/i.F_mol>0.1)])
    
    
# %% 

# =============================================================================
# Feedstock preprocessing
# =============================================================================

# The system as a whole, capital and operating costs already considered in 
# the cost of feedstock cost
@cost(basis='Flow rate', ID='System', units='kg/hr',
      kW=511.3205, cost=13329690, S=94697, CE=CEPCI[2009], n=0.6, BM=1.7)
class FeedstockPreprocessing(Unit):
    # 2205 U.S. ton/day (2000 metric tonne/day) as in Humbird et al.
    # _baseline_flow_rate = baseline_feedflow.sum()
    _cached_flow_rate = 2205

#%% Conversion

class BatchCoFermentation(BatchBioreactor):
    # Co-Fermentation time (hr)
    
    auxiliary_unit_names = ('base_mix_tank',)
    
    _N_ins = 7
    _N_outs = 2
    
    tau_cofermentation = 120 # initial value; updated by spec.load_productivity
    
    ferm_ratio = 1.
    
    ##### BASELINE AND RANGE FOR CSL LOADING AND CELL MASS YIELD #####
    
    # For N sources, the media in Markham et al. 2018 (https://doi.org/10.1073/pnas.172120311) 
    # contained:
    # 40 g/L peptone
    # 20 g/L yeast extract
    
    # Replacing yeast extract with CSL can actually enhance the pyruvate carboxylation pathway
    # in Y. lipolytica (Liu et al. 2015; https://doi.org/10.1016/j.indcrop.2015.10.029)
    
    # From Tan et al. 2016 (DOI 10.1088/1755-1315/36/1/012058):
    # Yeast extract total N = 10 - 11.8 dry wt%
    # CSL total N = 7.7 - 8.2 dry wt%
    
    # Peptone is rich in tryptophan; assumed 100% tryptophan
    # Tryptophan (C11H12N2O2) molar mass is 204.229 g
    # Therefore peptone total N = (2*14.0067)/204.229 = 13.717 dry wt%
    
    # Therefore, total N in media in Markham et al. 2018 is:
    # (a) Using low end of total N for yeast extract
    # = 40*0.13717 + 20*0.10
    # = 7.4868 g/L
    # (b) Using high end of total N for yeast extract
    # = 40*0.13717 + 20*0.118
    # = 7.8468 g/L
    
    # !!! METHOD I: Therefore, CSL concentration required for same total N
    # (a) Using low end of total N for yeast extract and high end for CSL
    # = 7.8468 / 0.077 = 101.9 g/L
    # (a) Using high end of total N for yeast extract and low end for CSL
    # = 7.4868 / 0.082 = 91.30 g/L
    
    # Dry cell mass in Markham et al. 2018 cannot easily be determined from viable cell count
    # From Cordova et al. 2020 (https://doi.org/10.1016/j.biortech.2020.123639):
    # Reached Yarrowia lipolitica dry cell mass of 47.8 g/L 
    # with a TAL titer of 21.6 g/L using corn stover hydrolysate
    
    # From Niehus et al. 2018 (DOI: 10.1155/2018/6393749)
    # Empirical formula of Yarrowia is CH1.82O0.51N0.12
    # Empirical molar mass of Yarrowia = (12*1 + 1*1.8 + 16*0.5 + 14*0.12) = 23.48 g
    # Therefore, total N content = 14*0.12/23.48
    # = 0.07155 g-N/g-cell
    # Therefore, total cellular N in media 
    # = 0.07155 * 47.8 = 3.42 g/L
    # !!! METHOD II: Therefore, CSL concentration required for same total N 
    # (a) Using low end of total N for CSL
    # = 3.42 / 0.077 = 44.416 g/L
    # (b) Using high end of total N for CSL
    # = 3.42 / 0.082 = 41.707 g/L
    
    # !!! Therefore, assumed range for CSL loading is 41.707 g/L - 112.1 g/L
    # Assumed baseline for CSL loading is mean of bounds = 76.903 g/L
    
    # !!! Assuming the same dry cell mass here, yield of cell mass on glucose
    # = (47.8/180) / (6 * 23.48 / 180.156) = 0.339
    
    ####### BASELINE AND RANGE FOR DAP LOADING #######
    # Humbird et al. 2011 DAP:CSL: 
    # 1. (0.33 g-DAP/L / 2.5 g-CSL/L)
    # 2. (0.67 g-DAP/L / 5 g-CSL/L)
    # ~= 0.133 g-DAP/g-CSL
    # therefore baseline is: 0.133 * 76.903 = 10.228 g/L
    # therefore range is:
        # low end: 0.133 * 41.707 = 5.547 g/L
        # high end 0.133 * 112.1 = 14.909 g/L
    
    ###############################################
    
    CSL_loading = 76.903 # g/L
    # CSL_loading = 32.5 # g/L
    
    DAP_loading = 10.228 # g/L
    
    regular_microbe_conversion = 0.339
    regular_citric_acid_conversion = 0.08856 # from Markham et al.; 16 g/L citrate from 180 g/L glucose
    
    acetate_target_loading = 13.667 * (60.05196/82.033789) # g-AceticAcid-eq / L # 13.667 g-sodium acetate /L as in Markham et al. 2018
    
    acetate_target_loading_default = 13.667 * (60.05196/82.033789) # g-AceticAcid-eq / L # 13.667 g-sodium acetate /L as in Markham et al. 2018 # default to this acetate loading if target TAL titer-yield combination is feasible
    
    air_m3_per_h_per_m3_reactor = 3.5*60/3 # 3.5 slpm for a 3L bioreactor; Markham et al. 2018 # used when aeration_rate_basis=='fixed rate basis'
    
    DO_saturation_concentration_kg_per_m3 = 7.8e-3 # 7.8 mg/L at 28 degrees C # used when aeration_rate_basis=='DO saturation basis'
    
    DO_saturation_target_level = 0.5 # 50% saturation from Markham et al. 2018 # used when aeration_rate_basis=='DO saturation basis'
    
    air_flow_rate_safety_factor_for_DO_saturation_basis = 1.
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 3
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.5 # !!!
    
    # autoselect_N  = True

    def _init(self, 
                 #  
                 T=28+273.15,
                 P=101325., 
                 tau=120, # initial value; updated by spec.load_productivity
                 V=3785.,
                 acetate_ID='SodiumAcetate',
                 aeration_rate_basis='fixed rate basis', # 'fixed rate basis' or 'DO saturation basis'
                 pH_to_load=6.5, # 'unmodified' or float value
                 ):
        BatchBioreactor._init(self,  T=T, P=P, tau=tau, V=V)
        
        self.pH_to_load = pH_to_load
        
        self.aeration_rate_basis = aeration_rate_basis
        
        self.acetate_ID = acetate_ID
        self.acetate_hydrolysis_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('SodiumAcetate + H2O -> AceticAcid + NaOH',        'SodiumAcetate',   1.-1e-4), 
        ])
        
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 0.666667 TAL + 2 CO2',        'Glucose',   0.19), 
        Rxn('Glucose -> 0.3 VitaminA',               'Glucose',   0.), # retinol
        Rxn('Glucose + O2 -> CitricAcid + H2O',               'Glucose',  self.regular_citric_acid_conversion), # 2H+ excluded # from Markham et al.; 16 g/L citrate from 180 g/L glucose
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   self.regular_microbe_conversion),
        
        Rxn('Xylose -> 0.555583 TAL + 1.3334 CO2',       'Xylose',    0.19),
        Rxn('Xylose -> 0.25 VitaminA',       'Xylose',    0.),
        Rxn('Xylose -> 0.8333 CitricAcid + H2O',               'Xylose',  self.regular_citric_acid_conversion),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    self.regular_microbe_conversion),
        
        Rxn('AceticAcid -> 0.333333 TAL + 3H2O',       'AceticAcid',    0.19),
        Rxn('AceticAcid -> 0.1 VitaminA',       'AceticAcid',    0.),
        Rxn('AceticAcid -> 0.0714 VitaminD2',       'AceticAcid',    0.),
        Rxn('AceticAcid -> 2 FermMicrobe',        'AceticAcid',    0.),
        
        ])
        
        self.CO2_generation_rxns = ParallelRxn([
            Rxn('Glucose + 6O2 -> 6CO2 + 6H2O', 'Glucose', 1.-1e-2),
            Rxn('Xylose + 5O2 -> 5CO2 + 5H2O', 'Xylose', 1.-1e-2),
            Rxn('AceticAcid + 2O2 -> 2CO2 + 2H2O', 'AceticAcid', 1.-1e-2),
            ])
        
        self.sucrose_hydrolysis_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Sucrose + H2O -> 2 Glucose',        'Sucrose',   1.-1e-4), 
        ])
        
        
        self.glucose_to_TAL_rxn = self.cofermentation_rxns[0]
        self.xylose_to_TAL_rxn = self.cofermentation_rxns[4]
        self.acetate_to_TAL_rxn = self.cofermentation_rxns[8]
        
        self.glucose_to_VitaminA_rxn = self.cofermentation_rxns[1]
        self.xylose_to_VitaminA_rxn = self.cofermentation_rxns[5]
        self.acetate_to_VitaminA_rxn = self.cofermentation_rxns[9]
        
        self.glucose_to_CitricAcid_rxn = self.cofermentation_rxns[2]
        self.xylose_to_CitricAcid_rxn = self.cofermentation_rxns[6]
        self.acetate_to_VitaminD2_rxn = self.cofermentation_rxns[10]
        
        self.glucose_to_microbe_rxn = self.cofermentation_rxns[3]
        self.xylose_to_microbe_rxn = self.cofermentation_rxns[7]
        self.acetate_to_microbe_rxn = self.cofermentation_rxns[11]
        
        self.glucose_to_CO2_rxn = self.CO2_generation_rxns[0]
        self.xylose_to_CO2_rxn = self.CO2_generation_rxns[1]
        self.acetate_to_CO2_rxn = self.CO2_generation_rxns[2]
        
        self._initialize_base_mix_tank()
    
    def _initialize_base_mix_tank(self):
        base_mix_tank =\
        self.auxiliary('base_mix_tank',
                       LiquidsMixingTank,
                       ins=('broth_before_pH_control', # effluent broth
                            # acetylacetone_decarboxylation_equilibrium, 
                            # recycled_nonevaporated_supernatant,
                            # 'recyled_top_prod_from_evaporating_supernatant',
                            # '', '', '',
                            self.ins[6], # base for pH control
                            ), 
                       outs=('pH_controlled_broth'),
                       tau=1.,
                       )
        
        base_mix_tank.mol_base_per_m3_broth = 0. # actually kmol-bsase/m3-broth, or mol-base/L-broth
        base_mix_tank.base_neutralizes_acids = True
        base_mix_tank.base_ID = 'NaOH'
        
        base_mix_tank.neutralization_rxns = SeriesReaction([
            Rxn('H3PO4 + 3NaOH -> SodiumPhosphate + H2O', 'H3PO4',   1.-1e-5),
            Rxn('CitricAcid + 3NaOH -> SodiumCitrate + H2O', 'CitricAcid',   1.-1e-5),
            Rxn('AceticAcid + NaOH -> SodiumAcetate + H2O', 'AceticAcid',   1.-1e-5),
            ])
        
        base_mix_tank.mol_base_per_m3_broth_needed_to_completely_neutralize_acids = 0.
        
        base_mix_tank.pH_to_load = self.pH_to_load
        base_mix_tank.load_pH = lambda pH: load_pH(pH, base_mix_tank)
        base_mix_tank.get_pH_maintained = lambda: get_pH_stream(base_mix_tank.outs[0])
        
        # @base_mix_tank.add_specification(run=False)
        def base_mix_tank_simulate_base_addition_and_acids_neutralization():
            base_mix_tank_ins_0 = base_mix_tank.ins[0]
            base_mix_tank_in_base = base_mix_tank.ins[1]
            base_mix_tank_in_base.empty()
            base_mix_tank_in_base.imol[base_mix_tank.base_ID] = base_mix_tank.mol_base_per_m3_broth * base_mix_tank_ins_0.F_vol
            
            base_mix_tank.mol_base_per_m3_broth_needed_to_completely_neutralize_acids =\
            min_base_req_to_completely_neutralize =\
                (base_mix_tank_ins_0.imol['AceticAcid']
                +3.* base_mix_tank_ins_0.imol['CitricAcid']
                + 3.*base_mix_tank_ins_0.imol['H3PO4'])/base_mix_tank_ins_0.F_vol
                                        
            base_mix_tank._run()
            
            base_mix_tank_outs_0_l = base_mix_tank.outs[0]['l']
            
            
            if base_mix_tank.base_neutralizes_acids:
                if base_mix_tank.mol_base_per_m3_broth < min_base_req_to_completely_neutralize:
                    
                    # base_mix_tank.neutralization_rxns = SeriesReaction([
                    #     Rxn('0.3333H3PO4 + NaOH -> 0.3333SodiumPhosphate + 0.3333H2O', 'NaOH',   1.-1e-5),
                    #     Rxn('0.3333CitricAcid + NaOH -> 0.3333SodiumCitrate + 0.3333H2O', 'NaOH',   1.-1e-5),
                    #     Rxn('AceticAcid + NaOH -> SodiumAcetate + H2O', 'NaOH',   1.-1e-5),
                    #     ])
                    
                    if base_mix_tank_outs_0_l.imol['NaOH'] > 3.*base_mix_tank_outs_0_l.imol['H3PO4']:
                        mol_H3PO4 = base_mix_tank_outs_0_l.imol['H3PO4']
                        base_mix_tank_outs_0_l.imol['H3PO4'] = 0.
                        base_mix_tank_outs_0_l.imol['NaOH'] -= 3.*mol_H3PO4
                        base_mix_tank_outs_0_l.imol['SodiumPhosphate'] += mol_H3PO4
                        base_mix_tank_outs_0_l.imol['H2O'] += mol_H3PO4
                    else:
                        mol_NaOH = base_mix_tank_outs_0_l.imol['NaOH']
                        base_mix_tank_outs_0_l.imol['NaOH'] = 0.
                        base_mix_tank_outs_0_l.imol['H3PO4'] -= 0.3333*mol_NaOH
                        base_mix_tank_outs_0_l.imol['SodiumPhosphate'] += 0.3333*mol_NaOH
                        base_mix_tank_outs_0_l.imol['H2O'] += 0.3333*mol_NaOH
                        
                        
                    if base_mix_tank_outs_0_l.imol['NaOH'] > 3.*base_mix_tank_outs_0_l.imol['CitricAcid']:
                        mol_CitricAcid = base_mix_tank_outs_0_l.imol['CitricAcid']
                        base_mix_tank_outs_0_l.imol['CitricAcid'] = 0.
                        base_mix_tank_outs_0_l.imol['NaOH'] -= 3.*mol_CitricAcid
                        base_mix_tank_outs_0_l.imol['SodiumCitrate'] += mol_CitricAcid
                        base_mix_tank_outs_0_l.imol['H2O'] += mol_CitricAcid
                    else:
                        mol_NaOH = base_mix_tank_outs_0_l.imol['NaOH']
                        base_mix_tank_outs_0_l.imol['NaOH'] = 0.
                        base_mix_tank_outs_0_l.imol['CitricAcid'] -= 0.3333*mol_NaOH
                        base_mix_tank_outs_0_l.imol['SodiumCitrate'] += 0.3333*mol_NaOH
                        base_mix_tank_outs_0_l.imol['H2O'] += 0.3333*mol_NaOH
                        
                        
                    if base_mix_tank_outs_0_l.imol['NaOH'] > base_mix_tank_outs_0_l.imol['AceticAcid']:
                        mol_AceticAcid = base_mix_tank_outs_0_l.imol['AceticAcid']
                        base_mix_tank_outs_0_l.imol['AceticAcid'] = 0.
                        base_mix_tank_outs_0_l.imol['NaOH'] -= mol_AceticAcid
                        base_mix_tank_outs_0_l.imol['SodiumAcetate'] += mol_AceticAcid
                        base_mix_tank_outs_0_l.imol['H2O'] += mol_AceticAcid
                    else:
                        mol_NaOH = base_mix_tank_outs_0_l.imol['NaOH']
                        base_mix_tank_outs_0_l.imol['NaOH'] = 0.
                        base_mix_tank_outs_0_l.imol['AceticAcid'] -= mol_NaOH
                        base_mix_tank_outs_0_l.imol['SodiumAcetate'] += mol_NaOH
                        base_mix_tank_outs_0_l.imol['H2O'] += mol_NaOH
                    
                else:
                    base_mix_tank.neutralization_rxns = SeriesReaction([
                        Rxn('H3PO4 + 3NaOH -> SodiumPhosphate + H2O', 'H3PO4',   1.),
                        Rxn('CitricAcid + 3NaOH -> SodiumCitrate + H2O', 'CitricAcid',   1.),
                        Rxn('AceticAcid + NaOH -> SodiumAcetate + H2O', 'AceticAcid',   1.),
                        ])
                    
                    base_mix_tank.outs[0].phase='l'
                    # base_mix_tank.neutralization_rxns.adiabatic_reaction(base_mix_tank.outs[0])
                    base_mix_tank.neutralization_rxns(base_mix_tank.outs[0])
            
            base_mix_tank.outs[0].phase='l'
        
        base_mix_tank.simulate_base_addition_and_acids_neutralization = base_mix_tank_simulate_base_addition_and_acids_neutralization
        
        @base_mix_tank.add_specification(run=False)
        def base_mix_tank_pH_loading_spec():
            base_mix_tank.ins[0].copy_like(self.outs[1])
            if base_mix_tank.pH_to_load == 'unmodified':
                base_mix_tank.mol_base_per_m3_broth = 0.
                base_mix_tank.simulate_base_addition_and_acids_neutralization()
            else:
                base_mix_tank.load_pH(base_mix_tank.pH_to_load)
    
    
    def _run(self):
        feed, seed, CSL, Acetate_spiking, DAP, air, base_for_pH_control = self.ins
        for i in [CSL, Acetate_spiking, DAP, air]: i.empty()
        
        vapor, effluent = self.outs
        
        vapor.empty()
        effluent.empty()
        
        effluent.mix_from([feed, seed])
        
        current_acetate_loading = effluent.imass[self.acetate_ID] / effluent.F_vol
        required_acetate_spiking = max(0, self.acetate_target_loading - current_acetate_loading)
        Acetate_spiking.imass[self.acetate_ID] = required_acetate_spiking * effluent.F_vol
        
        # air.P = self.air_pressure
        air.imol['N2'] = 0.79
        air.imol['O2'] = 0.21
        
        if self.aeration_rate_basis == 'DO saturation basis':
            self.air_exit_F_mol_needed = (1./0.21) * (1/32.) * self.DO_saturation_concentration_kg_per_m3 * self.DO_saturation_target_level\
                *(seed.F_vol+feed.F_vol)
            
            air.F_mol = 1e8 # initial value; updated after reactions
        
        elif self.aeration_rate_basis == 'fixed rate basis':
            air.F_vol = self.air_m3_per_h_per_m3_reactor * (seed.F_vol+feed.F_vol) * self.tau
        
        else: raise RuntimeError(f"Unsupported aeration_rate_basis ({self.aeration_rate_basis}); must be 'fixed rate basis' or 'DO saturation basis'.")
        
        CSL.imass['CSL'] = (seed.F_vol+feed.F_vol) * self.CSL_loading 
        
        DAP.imass['DAP'] = (seed.F_vol+feed.F_vol) * self.DAP_loading 
        
        effluent.mix_from([effluent, Acetate_spiking, air, CSL, DAP])
        effluent.T = vapor.T = self.T
        
        self.acetate_hydrolysis_rxns(effluent.mol)
        self.sucrose_hydrolysis_rxns(effluent.mol)
        self.cofermentation_rxns(effluent.mol)
        self.CO2_generation_rxns(effluent.mol)
        
        
        if self.aeration_rate_basis == 'DO saturation basis':
            O2_mol_remaining, N2_mol_remaining = effluent.imol['O2', 'N2']
            O2_mol_excess = O2_mol_remaining - 0.21*self.air_exit_F_mol_needed
            # N2_mol_excess = N2_mol_remaining - 0.79*air_exit_F_mol_needed
            N2_mol_excess = O2_mol_excess * 0.79/0.21
            air_mol_excess = O2_mol_excess + N2_mol_excess
            
            self.air_mol_excess = air_mol_excess
        
            air.F_mol -= air_mol_excess
            effluent.imol['O2'] -= O2_mol_excess
            effluent.imol['N2'] -= N2_mol_excess
        
        aeration_safety_factor = self.air_flow_rate_safety_factor_for_DO_saturation_basis
        effluent.imol['O2'] += (aeration_safety_factor-1) * air.imol['O2']
        effluent.imol['N2'] += (aeration_safety_factor-1) * air.imol['N2']
        air.F_mol += (aeration_safety_factor-1) * air.F_mol
        
        vapor.imol['CO2', 'O2', 'N2'] = effluent.imol['CO2', 'O2', 'N2']
        vapor.phase = 'g'
        effluent.imol['CO2', 'O2', 'N2'] = 0, 0, 0
        effluent.imass['CSL', 'DAP'] = 0, 0
        
        vapor.imol['CO2'] += CSL.get_atomic_flow('C')
        
        # self.effluent_titer = compute_TAL_titer(effluent)
        # self.effluent = effluent
        base_mix_tank = self.base_mix_tank
        base_mix_tank.pH_to_load = self.pH_to_load
        base_mix_tank.specifications[0]()
        # base_mix_tank._run()
        effluent.copy_like(base_mix_tank.outs[0])
        
    @property
    def effluent_titer(self):
        return compute_TAL_titer(self.effluent)
    
# class AeratedTALCoFermentation(AeratedBioreactor):
    
# Seed train, 5 stages, 2 trains
@cost(basis='Seed fermenter size', ID='Stage #1 fermenters', units='kg',
      # 44339, 211, and 26 are streams 303, 309, and 310 in Humbird et al.
      cost=75400, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #2 fermenters', units='kg',
      cost=116600, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #3 fermenters', units='kg',
      cost=157600, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #4 fermenters', units='kg',
      cost=352000, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Seed fermenter size', ID='Stage #4 agitators', units='kg',
      kW=11.1855, cost=26000, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Seed fermenter size', ID='Stage #5 fermenters', units='kg',
      cost=1180000, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Seed fermenter size', ID='Stage #5 agitators', units='kg',
      kW=14.914, cost=43000, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pumps', units='kg/hr',
      kW=59.656, cost=24300, S=43149, CE=CEPCI[2009], n=0.8, BM=2.3)
class SeedTrain(Unit):
    _N_ins = 1
    _N_outs = 2
    _units= {'Seed fermenter size': 'kg',
             'Flow rate': 'kg/hr'}
    
    # Cycle time for each batch (hr), including 12 hr turnaround time 
    tau_batch = 36
    
    
    regular_microbe_conversion = 0.339
    regular_citric_acid_conversion = 0.08856 # from Markham et al.; 16 g/L citrate from 180 g/L glucose
    
    def _init(self,  T=30+273.15, 
                 ferm_ratio=0.95, # ferm_ratio is the ratio of conversion relative to the fermenter
                 ):
        self.T = T
        self.ferm_ratio = ferm_ratio
        self.heat_exchanger = HXutility(None, None, None, T=T) 
        
        
        self.acetate_hydrolysis_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('SodiumAcetate + H2O -> AceticAcid + NaOH',        'SodiumAcetate',   1.-1e-4), 
        ])
        
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 0.666667 TAL + 2 CO2',        'Glucose',   0.19), 
        Rxn('Glucose -> 0.3 VitaminA',               'Glucose',   0.), # retinol
        Rxn('Glucose -> CitricAcid + H2O',               'Glucose',  self.regular_citric_acid_conversion), # 2H+ excluded # from Markham et al.; 16 g/L citrate from 180 g/L glucose
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   self.regular_microbe_conversion),
        
        Rxn('Xylose -> 0.555583 TAL + 1.3334 CO2',       'Xylose',    0.19),
        Rxn('Xylose -> 0.25 VitaminA',       'Xylose',    0.),
        Rxn('Xylose -> 0.8333 CitricAcid + H2O',               'Xylose',  self.regular_citric_acid_conversion),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    self.regular_microbe_conversion),
        
        Rxn('AceticAcid -> 0.333333 TAL + 2CO2 + 6H2O',       'AceticAcid',    0.19),
        Rxn('AceticAcid -> 0.1 VitaminA',       'AceticAcid',    0.),
        Rxn('AceticAcid -> 0.0714 VitaminD2',       'AceticAcid',    0.),
        Rxn('AceticAcid -> 2 FermMicrobe',        'AceticAcid',    0.),
        
        ])
        
        for i in self.cofermentation_rxns:
            i.X *= ferm_ratio
        
        # oxygen currently excluded from balance
        self.CO2_generation_rxns = ParallelRxn([
            Rxn('Glucose -> 6CO2 + 6H2O', 'Glucose', 1.),
            Rxn('Xylose -> 5CO2 + 5H2O', 'Xylose', 1.),
            Rxn('AceticAcid -> 2CO2 + 2H2O', 'AceticAcid', 1.),
            ])
        
        self.sucrose_hydrolysis_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Sucrose + H2O -> 2 Glucose',        'Sucrose',   1.-1e-4), 
        ])
        
        
        self.glucose_to_TAL_rxn = self.cofermentation_rxns[0]
        self.xylose_to_TAL_rxn = self.cofermentation_rxns[4]
        self.acetate_to_TAL_rxn = self.cofermentation_rxns[8]
        
        self.glucose_to_VitaminA_rxn = self.cofermentation_rxns[1]
        self.xylose_to_VitaminA_rxn = self.cofermentation_rxns[5]
        self.acetate_to_VitaminA_rxn = self.cofermentation_rxns[9]
        
        self.glucose_to_CitricAcid_rxn = self.cofermentation_rxns[2]
        self.xylose_to_CitricAcid_rxn = self.cofermentation_rxns[6]
        self.acetate_to_CitricAcid_rxn = self.cofermentation_rxns[10]
        
        self.glucose_to_microbe_rxn = self.cofermentation_rxns[3]
        self.xylose_to_microbe_rxn = self.cofermentation_rxns[7]
        self.acetate_to_microbe_rxn = self.cofermentation_rxns[11]
        
        self.glucose_to_CO2_rxn = self.CO2_generation_rxns[0]
        self.xylose_to_CO2_rxn = self.CO2_generation_rxns[1]
        self.acetate_to_CO2_rxn = self.CO2_generation_rxns[2]

    def _run(self):
        feed, = self.ins
        
        effluent, vapor = self.outs
        effluent.mix_from([feed,])
        
        effluent.T = vapor.T = self.T
        
        self.acetate_hydrolysis_rxns(effluent.mol)
        self.sucrose_hydrolysis_rxns(effluent.mol)
        self.cofermentation_rxns(effluent.mol)
        self.CO2_generation_rxns(effluent.mol)
        
        vapor.imol['CO2'] = effluent.imol['CO2']
        vapor.phase = 'g'
        effluent.imol['CO2'] = 0
        effluent.imass['CSL'] = 0
        

    def _design(self):
        Design = self.design_results
        Design['Flow rate'] = self.outs[0].F_mass
        Design['Seed fermenter size'] = self.outs[0].F_mass * self.tau_batch
        duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
        # duty = self.Hnet
        mixed_feed = tmo.Stream()
        mixed_feed.mix_from(self.outs)
        mixed_feed.T=self.ins[0].T
        # mixed_feed.vle(T=mixed_feed.T, P=mixed_feed.P)
        mixed_feed-0-self.heat_exchanger
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(mixed_feed,), 
                                                            duty=duty,
                                                            vle=False)
        
    def _cost(self):
        Unit._cost()
        self._decorated_cost()
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities += hx.heat_utilities
    
    def mol_atom_in(self, atom):
        return sum([stream.get_atomic_flow(atom) for stream in self.ins])
    
    def mol_atom_out(self, atom):
        return sum([stream.get_atomic_flow(atom) for stream in self.outs])
    
# Seed hold tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=439000, S=40414, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=11.1855, cost=31800, S=40414, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=7.457, cost=8200, S=43149, CE=CEPCI[2009], n=0.8, BM=2.3)
class SeedHoldTank(Unit): pass


# %% 

# =============================================================================
# Separation
# =============================================================================

class TALCrystallizer(BatchCrystallizer):
    TAL_solubility_multiplier = 1.
    _SA_vol_per_mass = 0.0008252419812169215
    def _init(self,  
                 target_recovery=0.6,
                 thermo=None,
                 tau=8, # assumed; uncertainty range is 2-14 h
                 N=5, 
                 V=None, T=273.15 + 1.,
                 basis='water solubility',
                 get_mol_TAL_dissolved_given_T_and_mol_water=None,
                 fixed_operating_T=273.15 + 1.,
                 Nmin=2, Nmax=36, 
                 T_range=(274., 372.5),
                 vessel_material='Stainless steel 316',
                 output_conc_multiplier=1.,
                 kW=0.00746):
        
        BatchCrystallizer._init(self, 
        # BatchCrystallizer._init(self,
        
                                     # thermo,
                     tau, N, V, T,
                     Nmin, Nmax, vessel_material,
                     kW)
        # self._ins = ins 
        # self._outs = outs
        
        self.target_recovery = target_recovery
        self.T_range = T_range
        self.tau = tau
        self.basis = basis
        self.get_mol_TAL_dissolved_given_T_and_mol_water = get_mol_TAL_dissolved_given_T_and_mol_water
        self.fixed_operating_T = fixed_operating_T
        self.output_conc_multiplier = output_conc_multiplier
        
    
    def _run(self):
        in_stream, = self.ins
        out_stream, = self.outs
        target_recovery = self.target_recovery
        
        out_stream.copy_like(in_stream)
        # out_stream.sle(T=self.T, solute='SuccinicAcid')
        
        tot_TAL = in_stream.imol['TAL']
        
        if self.fixed_operating_T:
            self.T = self.fixed_operating_T
            
        out_stream.T = self.T
        
        TAL_solubility = self.TAL_solubility_multiplier *\
            self.get_mol_TAL_dissolved_given_T_and_mol_water(out_stream.T, out_stream.imol['Water'])
        out_stream.phases = ('s', 'l')
        TAL_dissolved = min(TAL_solubility, tot_TAL)
        out_stream.imol['l', 'TAL'] = TAL_dissolved
        out_stream.imol['s', 'TAL'] = max(0, tot_TAL - TAL_dissolved)




# %% 

# =============================================================================
# Storage
# =============================================================================

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=96000, S=1981, CE=CEPCI[2010], n=0.7, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=7493, S=1981, CE=CEPCI[2010], n=0.8, BM=2.3)
class SulfuricAcidStorageTank(Unit): pass

#%% Catalytic upgrading

class HydrogenationEstersReactor(Reactor):
    """
    A hydrogenation reactor to produce esters from TAL.
    """
    _N_ins = 5
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default,
            'PdC catalyst': 1}
    mcat_frac = 0.5 # fraction of catalyst by weight in relation to the reactant (TAL)
    mono_di_hydroxy_esters_conversion_ratio = 0.5
    TAL_to_esters_and_DHL_conversion = 0.9
    esters_DHL_conversion_ratio = 0.5
    TAL_to_esters_conversion = esters_DHL_conversion_ratio * TAL_to_esters_and_DHL_conversion
    # catalyst_deactivation_k = 2e-3/10
    catalyst_replacements_per_year = 0. # deprecated variable
    hydrogenation_rxns = ParallelRxn([
            #   Reaction definition   Reactant   Conversion
            Rxn('TAL + 2H2 + Ethanol -> Ethyl_5_hydroxyhexanoate',         'TAL',   TAL_to_esters_conversion*mono_di_hydroxy_esters_conversion_ratio),
            Rxn('TAL + H2 + Ethanol + H2O -> Ethyl_3_5_dihydroxyhexanoate',         'TAL',   TAL_to_esters_conversion*(1.-mono_di_hydroxy_esters_conversion_ratio)),
            Rxn('TAL + 2H2 -> HMTHP',         'TAL',   1e-6), # conversion from Huber group experimental data
            # Rxn('HMDHP + H2 -> HMTHP',         'HMDHP',   1.-1e-5),
            # ])
    # byproduct_formation_rxns  = ParallelRxn([
            #   Reaction definition   Reactant   Conversion
            Rxn('TAL + 3H2 -> DHL + H2O',         'TAL',   0.5-1e-3), # conversion from Huber group experimental data
            Rxn('TAL + H2 -> HMDHP',         'TAL',   1e-3),  # conversion from Huber group experimental data
            Rxn('HMDHP + H2 -> HMTHP',         'HMDHP',   1.),  # conversion from Huber group experimental data
            Rxn('HMTHP + Ethanol -> Ethyl_5_hydroxyhexanoate',         'HMTHP',   TAL_to_esters_conversion*mono_di_hydroxy_esters_conversion_ratio),  # conversion from Huber group experimental data
            Rxn('DHL + Ethanol -> Ethyl_3_5_dihydroxyhexanoate',         'DHL',   TAL_to_esters_conversion*mono_di_hydroxy_esters_conversion_ratio),  # conversion from Huber group experimental data
            ])
    
    TAL_to_esters_rxns = hydrogenation_rxns[:2]
    TAL_to_DHL_rxn = hydrogenation_rxns[3]
    byproduct_formation_rxns  = ParallelRxn([
            #   Reaction definition   Reactant   Conversion
            Rxn('TAL + H2O -> Pentenone + CO2',         'TAL',   1.-1e-5),
           ])
    TAL_conversion_rxns = [i for i in hydrogenation_rxns if i.reactant=='TAL']
    
    def _init(self,  
                 tau = (7./3.) * 2., # from Huber group
                 T=323., # from Huber group
                 vessel_material='Stainless steel 316', **args):
        Reactor._init(self,  tau=tau, vessel_material=vessel_material)
        self.T = T
        self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        
    def _run(self):
        feed, recycle, reagent, fresh_catalyst, recovered_catalyst = self.ins
        effluent, vented_gas, spent_catalyst = self.outs
        
        # effluent.imol['HMDHP'] += 1e-10
        
        H2_mol_needed_for_HMTHP_formation = (sum(feed.imol['TAL', 'HMDHP']) + sum(recycle.imol['TAL', 'HMDHP']))*2*self.hydrogenation_rxns[0].X
        reagent.imol['H2'] = 20*H2_mol_needed_for_HMTHP_formation
        reagent.phase = 'g'
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle, reagent])
        
        # effluent.imol['HMTHP'] -= 1e-10
        
        effluent.phase = 'l'
        
        
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.hydrogenation_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
        
        extra_H2_mol_to_exclude = effluent.imol['H2']
        
        vented_gas.phase = 'g'
        vented_gas.imol['CO2'] = effluent.imol['CO2']
        # vented_gas.imol['H2'] = effluent.imol['H2']
        
        effluent.imol['CO2'] = 0.
        # effluent.imol['H2'] = 0.
        
        reagent.imol['H2'] -= extra_H2_mol_to_exclude
        effluent.imol['H2'] = 0.
        
        for i in self.outs:
            i.T = self.T
        
        spent_catalyst.imass['Pd'] =\
            self.mcat_frac * self.ins[0].imass['TAL'] * self.tau * self.catalyst_replacements_per_year/ 7884.0
        # spent_catalyst.imass['Pd'] = fresh_catalyst.imass['Pd'] =\
            
            # self.catalyst_deactivation_k * (feed.imol['TAL'] - effluent.imol['TAL'])
            # self.catalyst_deactivation_k * reagent.imol['H2']
            
            # effluent.imol['Ethyl_5_hydroxyhexanoate'])
            # effluent.imol['Ethyl_3_5_dihydroxyhexanoate'])
            # effluent.imol['DHL'])
            
        
            # Alternative 
            # reagent.imol['H2']
        
        # catalyst leaves in the effluent (recovered downstream)
        effluent.imass['Pd'] = self.mcat_frac * self.ins[0].imass['TAL']
        fresh_catalyst.imass['Pd'] = max(0, effluent.imass['Pd'] - recovered_catalyst.imass['Pd'])
        
    def _design(self):
        Reactor._design(self)
        duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
        # mixed_feed = tmo.Stream()
        # mixed_feed.mix_from(self.ins)
        # self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, mixed_feed)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=self.ins, duty=duty)
        
    def _cost(self):
        Reactor._cost(self)
        self.purchase_costs['PdC catalyst'] =\
            self.mcat_frac * self.ins[0].imass['TAL'] * self.tau * price['PdC']
            



class HydrogenationReactor(StirredTankReactor):
    """
    A hydrogenation reactor.
    """
    _N_ins = 5
    _N_outs = 3
    
    TEA_operating_hours = 4320
    
    # auxiliary_unit_names = ('heat_exchanger')
    _F_BM_default = {**StirredTankReactor._F_BM_default,
            'Ni-SiO2 catalyst': 1.}
    
    mcat_frac = 0.2 # fraction of catalyst by weight in relation to the reactant (TAL) # from Huber group
    
    spent_catalyst_replacements_per_year = 1. # number of times the entire catalyst_weight is replaced per year
    
    catalyst_weight = 0. # updated in _run
    
    def _init(self,  
                 tau = 9.4, # from Huber group
                 T=100. + 273.15, # from Huber group
                 P=3.5e6, # 35 bar # from Huber group
                 vessel_material='Stainless steel 316',
                 NiSiO2_catalyst_price=price['Ni-SiO2'],
                 batch=True,
                 rigorous_hx=False,
                 **args):
        super()._init(tau=tau, T=T, P=P, batch=batch, vessel_material=vessel_material)
        # self.T = T
        # self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        self.NiSiO2_catalyst_price = NiSiO2_catalyst_price
        self.heat_exchanger.rigorous = rigorous_hx
        self.hydrogenation_rxns = hydrogenation_rxns = ParallelRxn([
                #   Reaction definition   Reactant   Conversion
                Rxn('TAL + 2H2 -> HMTHP',         'TAL',   0.969), # conversion from Huber group experimental data
                # Rxn('HMDHP + H2 -> HMTHP',         'HMDHP',   1.-1e-5)
                ])
        self.byproduct_formation_rxns  = ParallelRxn([
                #   Reaction definition   Reactant   Conversion
                Rxn('TAL + 3H2 -> DHL + H2O',         'TAL',   1.-1e-5), # conversion from Huber group experimental data
                # Rxn('TAL + H2 -> HMDHP',         'TAL',   1-0.2125),  # conversion from Huber group experimental data
                ])
        
        self.TAL_to_HMP_rxn = hydrogenation_rxns[0]
        
    def _run(self):
        feed, recycle, reagent, recovered_catalyst, fresh_catalyst = self.ins
        vent, spent_catalyst, effluent = self.outs
        effluent.empty()
        effluent.phase = 'l'
        vent.phase = 'g'
        # effluent.imol['HMDHP'] += 1e-10
        tau = self.tau
        
        self.catalyst_weight = cat_weight = self.mcat_frac * sum(i.imass['TAL'] for i in self.ins) * tau
        
        H2_mol_needed_for_HMTHP_formation = 2 * (sum(feed.imol['TAL', 'HMDHP']) + sum(recycle.imol['TAL', 'HMDHP']))
        # *self.hydrogenation_rxns[0].X
        
        reagent.imol['H2'] = 5*H2_mol_needed_for_HMTHP_formation
        reagent.phase = 'g'
        
        effluent.mix_from([feed, recycle, reagent,])
        
        effluent.phase = 'l'
        
        self.hydrogenation_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
        
        extra_H2_mol_to_exclude = effluent.imol['H2']
        reagent.imol['H2'] -= extra_H2_mol_to_exclude
        effluent.imol['H2'] = 0.
        
        # 
        current_cat_mass_flow = recovered_catalyst.imass['NiSiO2']
        spent_catalyst_mass_flow = self.spent_catalyst_replacements_per_year*cat_weight/self.TEA_operating_hours # kg/h
        current_cat_mass_flow-=spent_catalyst_mass_flow
        req_cat_mass_flow = min(spent_catalyst_mass_flow, cat_weight/tau)
        
        spent_catalyst.phase = 's'
        spent_catalyst.imass['NiSiO2'] = spent_catalyst_mass_flow
        
        fresh_catalyst.phase = 's'
        fresh_catalyst.imass['NiSiO2'] = fresh_cat_mass_flow = max(0, req_cat_mass_flow - current_cat_mass_flow)
        
        effluent.T = self.T
        effluent.P = self.P
        effluent.phases = ('l', 's')
        effluent.imass['s', 'NiSiO2'] = fresh_cat_mass_flow+current_cat_mass_flow
        

        
    def _design(self):
        super()._design()
        # duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
        # mixed_feed = tmo.Stream()
        
        # for i in self.outs: i.phase = 'l'
        # mixed_feed.mix_from(self.outs[0])
        # # for i in self.outs:
        #     # if i.imol['NiSiO2']:
        # spent_catalyst = self.outs[1]
        # spent_catalyst.phases = ('l', 's')
        # spent_catalyst.imol['s', 'NiSiO2'] = i.imol['l', 'NiSiO2']
        # spent_catalyst.imol['l', 'NiSiO2'] = 0.
            
        # mixed_feed.T=self.ins[0].T
        # # mixed_feed.vle(T=mixed_feed.T, P=mixed_feed.P)
        # self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(mixed_feed,), 
        #                                                     duty=duty,
        #                                                     vle=False)
    
    def _cost(self):
        super()._cost()
        
        # divide catalyst cost by number of reactors as bst.StirredTankReactor currently multiplies
        # all purchase_cost values by number of reactors
        self.purchase_costs['Ni-SiO2 catalyst'] = self.catalyst_weight * self.NiSiO2_catalyst_price / self.N_reactors
        
        # hx = self.heat_exchanger
        # self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        # self.heat_utilities += hx.heat_utilities
    
class DehydrationReactor(StirredTankReactor):
    """
    A dehydration reactor.
    """
    _N_ins = 4
    _N_outs = 3
    
    TEA_operating_hours = 4320
    
    # auxiliary_unit_names = ('heat_exchanger')
    
    _F_BM_default = {**StirredTankReactor._F_BM_default,
            'Amberlyst-70 catalyst': 1}
    mcat_frac = 0.5 # fraction of catalyst by weight in relation to the reactant (HMTHP)
    
    spent_catalyst_replacements_per_year = 1. # number of times the entire catalyst_weight is replaced per year
    
    def _init(self,  
                 tau = 6.1, # from Huber group
                 T=160. + 273.15, # from Huber group
                 P=2e6, # 20 bar # from Huber group
                 vessel_material='Stainless steel 316',
                 Amberlyst70_catalyst_price=price['Amberlyst-70'],
                 batch=True,
                 rigorous_hx=False,
                 **args):
        super()._init(T = T, tau=tau, P=P, batch=batch, vessel_material=vessel_material)
        # self.T = T
        # self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        self.Amberlyst70_catalyst_price = Amberlyst70_catalyst_price
        self.heat_exchanger.rigorous = rigorous_hx
        self.dehydration_rxns = dehydration_rxns = ParallelRxn([
                #   Reaction definition                                       Reactant   Conversion
                Rxn('HMTHP -> PSA + H2O',         'HMTHP',   0.871) # conversion from Huber group experimental data
                    ])
        self.byproduct_formation_rxns  = ParallelRxn([
                #   Reaction definition   Reactant   Conversion
                Rxn('HMTHP -> DHL + H2O',         'HMTHP',   1.-1e-5), # conversion from Huber group experimental data
                # Rxn('TAL + H2 -> HMDHP',         'TAL',   1-0.2125),  # conversion from Huber group experimental data
                ])
        self.HMP_to_PSA_rxn = dehydration_rxns[0]
        
    def _run(self):
        feed, recycle, recovered_catalyst, fresh_catalyst = self.ins
        vent, spent_catalyst, effluent = self.outs
        effluent.empty()
        effluent.phase = 'l'
        vent.phase = 'g'
        tau = self.tau
        self.catalyst_weight = cat_weight = self.mcat_frac * sum(i.imass['HMTHP'] for i in self.ins) * tau
        
        effluent.mix_from([feed, recycle])
        
        self.dehydration_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
        
        req_cat_mass_flow = cat_weight/tau
        current_cat_mass_flow = recovered_catalyst.imass['Amberlyst70_']
        spent_catalyst_mass_flow = self.spent_catalyst_replacements_per_year*cat_weight/self.TEA_operating_hours # kg/h
        current_cat_mass_flow-=spent_catalyst_mass_flow
        
        spent_catalyst.phase = 's'
        spent_catalyst.imass['Amberlyst70_'] = spent_catalyst_mass_flow
        
        fresh_catalyst.phase = 's'
        fresh_catalyst.imass['Amberlyst70_'] = fresh_cat_mass_flow = max(0, req_cat_mass_flow - current_cat_mass_flow)
        
        effluent.T = spent_catalyst.T = self.T
        effluent.P = self.P
        effluent.phases = ('l', 's')
        effluent.imass['s', 'Amberlyst70_'] = fresh_cat_mass_flow+current_cat_mass_flow
        
        
    def _design(self):
        super()._design()
        # duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
        # mixed_feed = tmo.Stream()
        
        # for i in self.outs: i.phase = 'l'
        # mixed_feed.mix_from(self.outs[0])
        # # for i in self.outs:
        # #     if i.imol['Amberlyst70_']:
        # spent_catalyst = self.outs[1]
        # spent_catalyst.phases = ('l', 's')
        # spent_catalyst.imol['s', 'Amberlyst70_'] = i.imol['l', 'Amberlyst70_']
        # spent_catalyst.imol['l', 'Amberlyst70_'] = 0.
                
        # mixed_feed.T=self.ins[0].T
        # # mixed_feed.vle(T=mixed_feed.T, P=mixed_feed.P)
        # self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(mixed_feed,), 
        #                                                     duty=duty,
        #                                                     vle=False)
    
    def _cost(self):
        super()._cost()
        # divide catalyst cost by number of reactors as bst.StirredTankReactor currently multiplies
        # all purchase_cost values by number of reactors
        self.purchase_costs['Amberlyst-70 catalyst'] = self.catalyst_weight * self.Amberlyst70_catalyst_price / self.N_reactors

        # hx = self.heat_exchanger
        # self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        # self.heat_utilities += hx.heat_utilities

            

class RingOpeningHydrolysisReactor(StirredTankReactor):
    """
    A ring-opening and hydrolysis reactor.
    """
    _N_ins = 3
    _N_outs = 2
    # auxiliary_unit_names = ('heat_exchanger')
    
    # _F_BM_default = {**StirredTankReactor._F_BM_default,
    #         'Heat exchangers': 3.,}
    

    def _init(self,  
                 tau = 19., # from Huber group
                 T=130. + 273.15, # from Huber group
                 P=2e6, # 20 bar # from Huber group
                 vessel_material='Stainless steel 316', 
                 batch=True,
                 rigorous_hx=False,
                 **args):
        super()._init(T=T, tau=tau, P=P, batch=batch, vessel_material=vessel_material)
        self.heat_exchanger.rigorous = rigorous_hx
            
        self.ring_opening_rxns = ring_opening_rxns = ParallelRxn([
                #   Reaction definition                                       Reactant   Conversion
                Rxn('PSA -> SA',         'PSA',   0.999) # conversion from Huber group
                    ])
        
        self.PSA_to_SA_rxn = ring_opening_rxns[0]
        self.hydrolysis_rxns = hydrolysis_rxns= ParallelRxn([
                #   Reaction definition                                       Reactant   Conversion
                Rxn('SA + KOH -> KSA + H2O',         'SA',   1.-1e-5) # assumed 
                    ])
        self.SA_to_KSA_rxn = hydrolysis_rxns[0]
        self.byproduct_formation_rxns = ParallelRxn([
                #   Reaction definition                                       Reactant   Conversion
                Rxn('PSA -> 0.2PolyPSA',         'PSA',   1.-1e-5) # assumed
                    ])
        
    def _run(self):
        feed, recycle, reagent = self.ins
        vent, effluent = self.outs
        vent.phase = 'g'
        # effluent = feed.copy()
        
        
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle])
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.ring_opening_rxns(effluent.mol)
        
        reagent.imol['KOH'] = max(0, effluent.imol['SA'] - effluent.imol['KOH'])
        
        effluent.mix_from([effluent, reagent])
        
        self.hydrolysis_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
        
        effluent.T = self.T
        effluent.P = self.P

        
        
    def _design(self):
        super()._design()
    
    def _cost(self):
        super()._cost()
            

            
class RingOpeningReactor(Reactor):
    """
    A ring-opening reactor.
    """
    _N_ins = 2
    _N_outs = 1
    _F_BM_default = {**Reactor._F_BM_default,
            'RaneyNi catalyst': 1}
    mcat_frac = 0.01 # fraction of catalyst by weight in relation to the reactant (TAL)
    ring_opening_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('PSA -> SA',         'PSA',   1.) # conversion from Chia et al. 2012
                ])
    byproduct_formation_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('PSA -> 0.2PolyPSA',         'PSA',   1.-0.089) # conversion from Chia et al. 2012
                ])
    PSA_to_SA_rxn = ring_opening_rxns[0]
    
    def _run(self):
        feed, recycle = self.ins
        effluent = self.outs[0]
        
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle])
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.ring_opening_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
    def _cost(self):
        super()._cost()
        self.purchase_costs['RaneyNi catalyst'] =\
            self.mcat_frac * self.ins[0].imass['HMTHP'] * price['RaneyNi']
            
class DehydrationRingOpeningReactor(Reactor):
    """
    A dehydration reactor.
    """
    _N_ins = 2
    _N_outs = 1
    _F_BM_default = {**Reactor._F_BM_default,
            'RaneyNi catalyst': 1}
    mcat_frac = 0.01 # fraction of catalyst by weight in relation to the reactant (TAL)
    dehydration_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('HMTHP -> SA',         'HMTHP',   0.667) # conversion from Chia et al. 2012
                ])
    TAL_to_SA_rxn = dehydration_rxns[0]
    
    
    def _run(self):
        feed, recycle = self.ins
        effluent = self.outs[0]
        
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle])
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.dehydration_rxns(effluent.mol)
        
    def _cost(self):
        super()._cost()
        self.purchase_costs['RaneyNi catalyst'] =\
            self.mcat_frac * self.ins[0].imass['HMTHP'] * price['RaneyNi']


class HydrolysisReactor(Reactor):
    """
    A hydrogenation reactor.
    """
    _N_ins = 3
    _N_outs = 1
    _F_BM_default = {**Reactor._F_BM_default,
            'Amberlyst15 catalyst': 1}
    mcat_frac = 0.01 # fraction of catalyst by weight in relation to the reactant (TAL)
    hydrolysis_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('SA + KOH -> KSA + H2O',         'SA',   1.) # not mentioned in Viswanathan et al. 2020
                ])
    TAL_to_SA_rxn = hydrolysis_rxns[0]
    
    def _run(self):
        feed, recycle, reagent = self.ins
        effluent = self.outs[0]
        reagent.imol['KOH'] = max(0, feed.imol['SA'] + recycle.imol['SA'] - recycle.imol['KOH'])
        
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle, reagent])
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.hydrolysis_rxns(effluent.mol)
        
    def _cost(self):
        super()._cost()
        self.purchase_costs['Amberlyst15 catalyst'] =\
            self.mcat_frac * self.ins[0].imass['SA'] * price['Amberlyst15']
            
class Crystallization(Reactor):
    N_ins = 3
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    dehydration_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('KSA + HCl -> SA + KCl',         'KSA',   0.98)
                ])
    TAL_to_SA_rxn = dehydration_rxns[0]
    
    def _init(self,  tau=1., T=293.15, reagent_fraction=1.):
        Reactor._init(self,  tau=tau)
        self.reagent_fraction = reagent_fraction
    

    def _run(self):
        feed, reagent, recycle_reagent, recycle_feed = ins = self.ins
        effluent, KCl = outs = self.outs
        reagent.empty()
        reagent.imol['HCl'] = max(0, sum([i.imol['KSA'] for i in ins]) - sum([i.imol['HCl'] for i in ins]))
        # effluent = feed.copy()
        effluent.mix_from([feed, 
                           reagent, 
                           recycle_reagent, 
                           recycle_feed])
        
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.dehydration_rxns(effluent.mol)
        
        KCl.imol['KCl'] = effluent.imol['KCl']
        KCl.T = effluent.T
        KCl.P = effluent.P
        
        # effluent.imol['H2O'] = effluent.imol['H2O']
        
        effluent.imol['HCl'] = 0.
        effluent.imol['KCl'] = 0
        # effluent.imol['H2O'] = 0
        
    def _cost(self):
        super()._cost()


class HClKOHRecovery(Reactor):
    N_ins = 2
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    dehydration_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('KCl + H2O -> HCl + KOH',         'KCl',   1.)
                ])
    TAL_to_SA_rxn = dehydration_rxns[0]
    
    def _run(self):
        feed, reagent = self.ins
        HCl, KOH = self.outs
        
        reagent.imol['H2O'] = feed.imol['KCl']
        # effluent = feed.copy()
        KOH.mix_from([feed, reagent])
        
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.dehydration_rxns(KOH.mol)
        
        HCl.imol['HCl'] = KOH.imol['HCl']
        HCl.T = KOH.T 
        HCl.P = KOH.P
        
        KOH.imol['HCl'] = 0
        
    def _cost(self):
        super()._cost()