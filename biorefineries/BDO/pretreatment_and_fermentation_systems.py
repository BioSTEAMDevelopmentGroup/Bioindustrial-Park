1# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 03:54:05 2021

@author: yrc2
"""
import biosteam as bst
from biosteam import Stream
from biosteam.process_tools import SystemFactory
from biorefineries import HP as hp
import biorefineries.BDO as bdo 
from biorefineries.BDO.process_settings import price
from biorefineries.BDO.utils import find_split, splits_df, baseline_feedflow
from biorefineries.BDO.chemicals_data import BDO_chemicals, chemical_groups, \
                                soluble_organics, combustibles
import numpy as np


@SystemFactory(
    ID='BDO_pretreatment_and_fermentation_sys',
    ins=[dict(ID='feedstock', phase='l', T=298.15, P=101325,
              price=0.0628405632131775, H2O=2.084e+04, Acetate=1509,
              Extract=1.221e+04, Sucrose=641.8, Protein=2584, 
              Glucan=2.921e+04, Mannan=500.1, Galactan=1192, 
              Xylan=1.628e+04, Arabinan=1984, Lignin=1.314e+04, 
              Ash=4109, units='kg/hr'),
         dict(ID='recycle_acetoin')],
    outs=[dict(ID='filtered_fermentation_effluent'),
          dict(ID='vent'),
          dict(ID='solids'),
          dict(ID='wastewater')],
)
def create_pretreatment_and_fermentation_system(ins, outs):
    """
    Create a pretreatment and fermentation system for the production of BDO 
    from cellulosic biomass.

    Parameters
    ----------
    ins : stream
        Feedstock
    outs : stream sequence
        [0] Filtered fermentation effluent.
        [1] Solids
        [3] Wastewater

    Examples
    --------
    >>> from biorefineries import BDO as bdo
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(bdo.BDO_chemicals)
    >>> BDO_pretreatment_and_fermentation_sys = bdo.create_pretreatment_and_fermentation_system()
    >>> BDO_pretreatment_and_fermentation_sys.simulate()
    >>> BDO_pretreatment_and_fermentation_sys.show()
    System: BDO_pretreatment_and_fermentation_sys
    ins...
    [0] feedstock
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O       1.16e+03
                        Acetate   25.1
                        Extract   67.8
                        Sucrose   1.87
                        Protein   113
                        Glucan    180
                        Mannan    3.08
                        ...
    [1] recycle_acetoin
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    outs...
    [0] filtered_fermentation_effluent
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                3.67e+03
                        AceticAcid         0.586
                        Glucose            27.6
                        2,3-Butanediol     190
                        GlucoseOligomer    6.5
                        Extract            62.8
                        Xylose             34
                        ...
    [1] cell_mass
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                59.9
                        AceticAcid         0.586
                        Glucose            1.04
                        GlucoseOligomer    0.251
                        Extract            2.43
                        Xylose             1.63
                        XyloseOligomer     0.106
                        ...
    [2] solids
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                704
                        AceticAcid         10.1
                        Glucose            5.69
                        GlucoseOligomer    0.261
                        Extract            2.53
                        Xylose             4.84
                        XyloseOligomer     0.11
                        ...
    [3] wastewater
        phase: 'l', T: 387.97 K, P: 101325 Pa
        flow (kmol/hr): H2O         1.68e+03
                        AceticAcid  4.89
                        Furfural    3.54
                        HMF         1.51
                        
    """

    feedstock, recycle_acetoin = ins
    filtered_fermentation_effluent, vent, solids, wastewater = outs

    # =============================================================================
    # Feedstock
    # =============================================================================
    
    U101 = hp.units.FeedstockPreprocessing('U101', ins=feedstock)
    
    # Handling costs/utilities included in feedstock cost thus not considered here
    U101.cost_items['System'].cost = 0
    U101.cost_items['System'].kW = 0

    # =============================================================================
    # Pretreatment streams
    # =============================================================================
    
    # For pretreatment, 93% purity
    sulfuric_acid_T201 = Stream('sulfuric_acid_T201', price=price['Sulfuric acid'],units='kg/hr')
    # To be mixed with sulfuric acid, flow updated in SulfuricAcidMixer
    water_M201 = Stream('water_M201', T=300, units='kg/hr')
        
    # To be used for feedstock conditioning
    water_M202 = Stream('water_M202', T=300, units='kg/hr')
    
    # To be added to the feedstock/sulfuric acid mixture, flow updated by the SteamMixer
    water_M203 = Stream('water_M203', phase='l', T=300, P=13.*101325, units='kg/hr')
    
    
    # =============================================================================
    # Pretreatment units
    # =============================================================================
    H_M201 = bst.HXutility('H_M201', ins=water_M201,
                                     outs='steam_M201',
                                     T=99.+273.15, rigorous=True)
    
    H_M201.heat_utilities[0].heat_transfer_efficiency = 1.
    def H_M201_specification():
        T201._run()
        acid_imass = T201.outs[0].imass['SulfuricAcid']
        H_M201.ins[0].imass['Water'] = acid_imass / 0.05
        # H_M201.ins[0].imass['H2SO4'] = H_M201.ins[0].imass['Water']/1000.
        H_M201._run()
    H_M201.specification = H_M201_specification
    # H_M201._cost = lambda: None
    # H_M201._design = lambda: None
    # H_M201.heat_utilities[0].heat_exchanger = None
    H_M202 = bst.HXutility('H_M202', ins=water_M202,
                                     outs='hot_water_M202',
                                     T=99.+273.15, rigorous=True)
    H_M202.heat_utilities[0].heat_transfer_efficiency = 1.
    def H_M202_specification():
        U101._run()
        H_M201.run()
        M201._run()
        feedstock, acid = U101.outs[0], M201.outs[0]
        recycled_water = H201.outs[0]
        mixture_F_mass = feedstock.F_mass + acid.F_mass
        mixture_imass_water = feedstock.imass['Water'] + acid.imass['Water'] + \
            recycled_water.imass['Water']
        total_mass = (mixture_F_mass - mixture_imass_water)/M202.solid_loading
        H_M202.ins[0].imass['Water'] = total_mass - mixture_F_mass
        # H_M202.ins[0].imass['H2SO4'] = H_M202.ins[0].imass['Water']/1000.
        H_M202._run()
    H_M202.specification = H_M202_specification
    
    
    # Prepare sulfuric acid
    get_feedstock_dry_mass = lambda: feedstock.F_mass - feedstock.imass['H2O']
    T201 = hp.units.SulfuricAcidAdditionTank('T201', ins=sulfuric_acid_T201,
                                          feedstock_dry_mass=get_feedstock_dry_mass())
    
    M201 = hp.units.SulfuricAcidMixer('M201', ins=(T201-0, H_M201-0))
        
    # Mix sulfuric acid and feedstock, adjust water loading for pretreatment
    M202 = hp.units.PretreatmentMixer('M202', ins=(U101-0, M201-0, H_M202-0, ''))
    
    # Mix feedstock/sulfuric acid mixture and steam
    M203 = bst.SteamMixer('M203', ins=(M202-0, water_M203), P=5.5*101325)
    M203.heat_utilities[0].heat_transfer_efficiency = 1.
    R201 = hp.units.PretreatmentReactorSystem('R201', ins=M203-0, outs=('R201_g', 'R201_l'))
    
    # Pump bottom of the pretreatment products to the oligomer conversion tank
    T202 = hp.units.BlowdownTank('T202', ins=R201-1)
    T203 = hp.units.OligomerConversionTank('T203', ins=T202-0)
    F201 = hp.units.PretreatmentFlash('F201', ins=T203-0,
                                   outs=('F201_waste_vapor', 'F201_to_fermentation'),
                                   P=101325, Q=0)
    H201 = bst.HXutility('H201', ins=F201-0,
                              V=0, rigorous=True)
    H202 = bst.HXutility('H202', ins=R201-0,
                              V=0, rigorous=True)
    
    P201 = hp.units.HydrolysatePump('P201', ins=F201-1)
    
    # =============================================================================
    # Conversion streams
    # =============================================================================
    
    # Flow and price will be updated in EnzymeHydrolysateMixer
    enzyme = Stream('enzyme', units='kg/hr', price=price['Enzyme'])
    # Used to adjust enzymatic hydrolysis solid loading, will be updated in EnzymeHydrolysateMixer
    enzyme_water = Stream('enzyme_water', units='kg/hr')
    
    # Corn steep liquor as nitrogen nutrient for microbes, flow updated in R301
    CSL = Stream('CSL', units='kg/hr')
    
    # For diluting concentrated, inhibitor-reduced hydrolysate
    dilution_water = Stream('dilution_water', units='kg/hr')
    
    
    # =============================================================================
    # Conversion units
    # =============================================================================
    
    # Cool hydrolysate down to fermentation temperature at 50°C
    H301 = bst.HXutility('H301', ins=P201-0, T=45+273.15, rigorous=True)
    
    T204 = bdo.units.LimeAdditionTank('T204', ins=(H301-0, 'lime_fresh'))
    
    @T204.add_specification(run=True)
    def update_lime_and_mix():
        hydrolysate = F201.outs[1]
        T204.ins[1].imol['CaO'] = hydrolysate.imol['H2SO4'] * 1.
    
    gypsum_split={'CaSO4':1.}
    S201 = hp.units.GypsumFilter('S201', ins=T204-0, outs=('gypsum', ''), split = gypsum_split,
                              moisture_content=0.2)
    
    
    # Mix enzyme with the cooled pretreatment hydrolysate
    M301 = hp.units.EnzymeHydrolysateMixer('M301', ins=(S201-1, enzyme, enzyme_water))
    
    # Saccharification
    R301 = hp.units.Saccharification('R301', 
                                    ins=M301-0,
                                    outs='saccharification_effluent')
    M303_P = bst.Pump('M303_P', ins=R301-0)
    # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S301_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S301_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S301_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S301 = hp.units.CellMassFilter('S301', ins=M303_P-0, 
                                moisture_content=0.35,
                                split=find_split(S301_index,
                                                  S301_cell_mass_split,
                                                  S301_filtrate_split,
                                                  chemical_groups))
    F301 = bst.MultiEffectEvaporator('F301', ins=S301-1, outs=('F301_l', 'F301_g'),
                                     P = (101325, 73581, 50892, 32777, 20000), V = 0.793)
    M204 = bst.Mixer('M204', ins=(H202-0, H201-0, F301-1), outs=wastewater)
    F301_P = bst.Pump('F301_P', ins=F301-0)
    
    def adjust_M304_water():
        M304.ins[1].imol['Water'] = (M304.water_multiplier - 1) * M304.ins[0].imol['Water']
        M304._run()
        
    M304 = bst.Mixer('M304', ins=(F301_P-0, dilution_water, ''))
    M304.water_multiplier = 1.
    M304.specification = adjust_M304_water
    M304_H = bst.HXutility('M304_H', ins=M304-0, T=30+273.15, rigorous=True)
    M304_H_P = bst.Pump('M304_H_P', ins=M304_H-0)
    S302 = bst.Splitter('S302', 
                        ins=M304_H_P-0, 
                        outs=('to_seedtrain', 'to_cofermentation'),
                        split=0.07)
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    R303 = bdo.units.SeedTrain('R303', ins=S302-0, outs=('seed', ''), ferm_ratio=0.9)
    
    T301 = hp.units.SeedHoldTank('T301', ins=R303-0)
    
    # Cofermentation
    R302 = bdo.units.CoFermentation('R302', 
                                    ins=(recycle_acetoin, S302-1, T301-0, CSL),
                                    outs=('fermentation_effluent', vent))
    
    
    # =============================================================================
    # Separation units
    # =============================================================================
    
    # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S401 = bst.SolidsCentrifuge('S401', ins=R302-0, 
                                      outs=('cell_mass', filtered_fermentation_effluent),
                                    # moisture_content=0.50,
                                    split=find_split(S401_index,
                                                      S401_cell_mass_split,
                                                      S401_filtrate_split,
                                                      chemical_groups), solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan'])
    M305 = bst.Mixer('M305', [S401-0, S301-0], solids)