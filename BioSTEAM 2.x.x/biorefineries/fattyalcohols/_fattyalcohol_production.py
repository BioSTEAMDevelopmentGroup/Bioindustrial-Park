# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import os
import biosteam as bst
import thermosteam as tmo
from biosteam import units
from thermosteam.reaction import Reaction as Rxn
from biorefineries.fattyalcohols import units as fa_units

__all__ = ('create_fattyalcohol_production_sys',)


# TODO: 
# - Add CSL and DAP as nutrients (remove tryptone)
# - Include stoichiometry similar to Humbird et. al.
# - Include salt price.

# Nitrogen requirement:
# - 0.0444 g N / g CSL (0.0725 mol N / mol CSL)
# - 0.212 g N / g DAP (2 mol N / mol DAP)
# - 0.0902 g N / g Cells (0.16 mol N / mol Cell)
# - 0.7362 g CSL / g DAP (2.611 mol CSL / mol DAP)

# Growth equation: 
# Let x be # mol DAP and y be # mol Cells
# - Mass balance:
#   MW_glucose + x * (2.611 * MW_csl + MW_dap) = y * MW_cells + 2 * MW_co2 + x * 1.3055 * MW_h2o + x * MW_h3po4
#   y = 4.3398 * x  + 3.71
# - Nitrogen balance:
#   x * (2.611 * 0.0725 + 2) = y * 0.16
#   y = 13.68 * x
# - Solve:
#   4.3398 * x  + 3.71 = 13.68 * x
#   x = 0.3972
#   y = 5.4338
# - Equation:
#   Glucose + 1.037 CSL + 0.3972 DAP ->  5.4338 Cells + 2 CO2 + 0.51858 H2O + 0.3972 H3PO4

# Bench scale stoichiometry: 
# Glucose -> 1.5 Hexanol + 1.5 Octanol + 4.5 Decanol + 
#            6 Dodecanol + 4 Tetradecanol + H2O + CO2
# Remove heavier alcohol products to ease separations
def create_fattyalcohol_production_sys(ID='fattyalcohol_production_sys',
    fa_stoichiometry='Glucose -> 1.5 Hexanol + 1.5 Octanol + 4.5 Decanol + H2O + CO2',
    yg_stoichiometry='Glucose ->  3.72 Cells + 2 CO2'):
    
    chemicals = tmo.settings.get_chemicals()
    
    ### Reactions ###

    fatty_alcohol_production = Rxn(fa_stoichiometry, 'Glucose', 0.5, basis='wt')
    fatty_alcohol_production.basis = 'mol'
    fatty_alcohol_production.correct_atomic_balance(
        constants = ['Hexanol', 'Octanol', 'Decanol', 'Dodecanol', 'Tetradecanol'])
    yeast_growth = Rxn(yg_stoichiometry, 'Glucose', .999, correct_mass_balance=True)
    
    ### Streams ###
    filepath = os.path.dirname(__file__)
    stream_data_path = os.path.join(filepath, 'streams.yaml')   
    stream_data = tmo.ThermoData.from_yaml(stream_data_path)
    glucose, process_water, tridecane, CSL, DAP, salt, \
    cell_recycle, mixed_bioreactor_feed = stream_data.create_streams(
        ['glucose', 'process_water', 'tridecane', 'CSL', 'DAP',
         'salt', 'cell_recycle', 'mixed_bioreactor_feed']
    )
    process_water.T = 310.15 # Assume heat integration so that this stream comes in heated

    ### Unit operations ###
    
    M101 = units.Mixer('M101', ins=(tridecane, 'solvent_recycle'))
    M101.ins[1].T = 310.15 # Assume heat integration so that this stream comes in heated
    
    one_week = 7*24
    T101 = units.StorageTank('T101', glucose, tau=one_week)
    T102 = units.StorageTank('T102', M101-0, tau=0.25)
    P102 = units.Pump('P102', T102-0)
    T103 = fa_units.CSLTank('T103', CSL)
    T104 = units.StorageTank('T104', salt, tau=one_week)
    T106 = fa_units.DAPTank('T106', DAP)
        
    def balance_feeds_to_bioreactor():
        feed_imol = mixed_bioreactor_feed.imol
        glucose.imol['Glucose'] = feed_imol['Glucose']
        process_water.imol['Water'] = feed_imol['Water'] - M101.ins[1].imol['Water']
        tridecane.imol['Tridecane'] = feed_imol['Tridecane'] - M101.ins[1].imol['Tridecane']
        salt.imol['NaCl'] = feed_imol['NaCl']
        cell_recycle.imol['Cells'] = feed_imol['Cells']
        CSL.imol['CSL'] = feed_imol['CSL']
        DAP.imol['DAP'] = feed_imol['DAP']
        M101._run()
    
    M101.specification = balance_feeds_to_bioreactor
    M102 = units.Mixer('M102',
                       ins=(T101-0, P102-0, T103-0, T104-0, T106-0,
                            process_water, cell_recycle))
    H101 = units.HXutility('H101',
                            ins=M102-0,
                            outs=mixed_bioreactor_feed,
                            material='Stainless steel/stainless steel',
                            T=310.15, V=0)
    def make_sure_enough_CSL_and_DAP_are_fed():
        # tmo.reaction.CHECK_FEASIBILITY = False
        feed = R101.ins[0]
        CSL.imass['CSL'] = feed.imass['CSL'] = 0.0025 * feed.F_mass
        DAP.imass['DAP'] = feed.imass['DAP'] = 0.33 * feed.F_vol
        R101._run()
        # effluent = R101.outs[1]
        # if effluent.imol['CSL'] < 0. or effluent.imol['DAP'] < 0.:
        #     CSL.imol['CSL'] = feed.imol['CSL'] = CSL_flow - effluent.imol['CSL']
        #     DAP.imol['DAP'] = feed.imol['DAP'] = DAP_flow - effluent.imol['DAP']
        #     effluent.imol['CSL'] = effluent.imol['DAP'] = 0.
        # tmo.reaction.CHECK_FEASIBILITY = True    
    R101 = fa_units.FattyAlcoholBioreactor('R101', outs=('CO2', ''), ins=H101-0,
                                     fermentation_reaction=fatty_alcohol_production,
                                     postfermentation_reaction=yeast_growth, T=310.15, tau=32, 
                                     V=3785, Nmin=2, Nmax=32)
    R101.specification = make_sure_enough_CSL_and_DAP_are_fed
    T105 = units.StorageTank('T105', ins=R101-1, tau=2)
    solids_split = chemicals.kwarray({'Cells': 0.999})
    solids_split[solids_split == 0] = 0.001
    solids_split[chemicals.index('Water')] = 0
    P104 = units.Pump('P104', ins=T105-0)
    C101 = fa_units.SolidLiquidsSplitCentrifuge('C101',
            ins=P104-0, outs=('', 'aqueous_fraction', 'cell_mass'),
            solids_split=solids_split,
            liquids_split=dict(
                Water=0.0004,
                Tridecane=0.9966,
                Hexanol=0.6781,
                Octanol=0.9608,
                Decanol=0.9966,
                Dodecanol=0.9999,
                Tetradecanol=1.0),
            moisture_content=0.45)
    P107 = units.Pump('P107', C101-0, 'oil_fraction', P=5*101325)
    
    ### System ###
    
    return bst.System(ID, (M101, T102, P102, T101, T104, T103, T106, M102, H101,
                           R101, T105, P104, C101, P107))

