1# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 03:54:05 2021

@author: sarangbhagwat, yrc2
"""
import thermosteam as tmo
import biosteam as bst
from biosteam.process_tools import SystemFactory
from biorefineries.BDO import units, facilities
from biorefineries.BDO.process_settings import price
import numpy as np
from biosteam import BoundedNumericalSpecification
from thermosteam.separations import partition
from flexsolve import IQ_interpolation
__all__ = (
    'create_separation_system_DPHP',
    'create_separation_system_oleyl_alcohol',
    'create_separation_system_oleyl_alcohol_3'
)

@SystemFactory(
    ID='BDO_separation_sys',
    ins=[dict(ID='filtered_fermentation_effluent', 
              phase='l', T=323.15, P=101325, 
              H2O=3663, AceticAcid=0.5861, 
              Glucose=27.61, BDO=191.1,
              GlucoseOligomer=6.5, Extract=62.82, 
              Xylose=33.96, XyloseOligomer=2.741, 
              Cellobiose=0.8368, Mannose=2.573, 
              MannoseOligomer=0.06861, Galactose=6.132, 
              GalactoseOligomer=0.1635, Arabinose=12.53, 
              ArabinoseOligomer=0.334, SolubleLignin=0.4022, 
              Protein=0.04538, Enzyme=23.95, 
              FermMicrobe=0.07988, Furfural=0.1058, 
              Acetoin=1.678, HMF=0.04494, 
              Glucan=0.003135, Mannan=3.481e-05, 
              Galactan=8.295e-05, Xylan=0.00139, 
              Arabinan=0.0001694, Lignin=0.00331, 
              Ash=0.02955, units='kmol/hr')],
    outs=[dict(ID='BDO'),
          dict(ID='unreacted_acetoin'),
          dict(ID='wastewater')],
)
def create_separation_system_DPHP(ins, outs):
    """
    Create a separation system for BDO using ethanol and DPHP for the 
    "salting out" effect.

    Parameters
    ----------
    ins : stream
        Fermentation effluent.
    outs : stream sequence
        [0] BDO
        [1] Unreacted acetoin
        [2] Wastewater

    Examples
    --------
    >>> from biorefineries import BDO
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(BDO.BDO_chemicals)
    >>> BDO_separation_sys = create_separation_system_DPHP()
    >>> BDO_separation_sys.simulate()
    >>> BDO_separation_sys.show()
    System: BDO_separation_sys
    ins...
    [0] filtered_fermentation_effluent
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                3.66e+03
                        AceticAcid         0.586
                        Glucose            27.6
                        2,3-Butanediol     191
                        GlucoseOligomer    6.5
                        Extract            62.8
                        Xylose             34
                        ...
    outs...
    [0] BDO
        phase: 'l', T: 454.81 K, P: 101325 Pa
        flow (kmol/hr): H2O                9.06e-06
                        Ethanol            0.109
                        Glucose            0.277
                        2,3-Butanediol     189
                        GlucoseOligomer    0.0652
                        Extract            0.63
                        Xylose             0.341
                        ...
    [1] unreacted_acetoin
        phase: 'l', T: 435.47 K, P: 101325 Pa
        flow (kmol/hr): 2,3-Butanediol     0.0945
                        3-Hydroxybutanone  1.59
    [2] wastewater
        phase: 'l', T: 374.32 K, P: 101325 Pa
        flow (kmol/hr): H2O                             3.66e+03
                        Ethanol                         0.454
                        AceticAcid                      0.586
                        Glucose                         27.3
                        2,3-Butanediol                  2.03
                        Dipotassium hydrogen phosphate  4.12
                        GlucoseOligomer                 6.43
                        ...

    """
    filtered_fermentation_effluent, = ins
    BDO, unreacted_acetoin, wastewater = outs
    
    ethanol_fresh = bst.Stream('ethanol_fresh', price=price['Ethanol'])
    DPHP_fresh = bst.Stream('DPHP_fresh', price=price['DPHP'])
    
    # DPHP storage
    #!!! Yalin suggests to use BioSTEAM's storage tank, and maybe we don't need the ConveryingBelt
    # (Yalin removed that from lactic acid biorefinery)
    T604 = units.DPHPStorageTank('T604', ins=DPHP_fresh)
    T604.line = 'DPHP storage tank'
    T604_P = bst.units.ConveyingBelt('T604_P', ins=T604-0)
    
    # 7-day storage time, similar to ethanol's in Humbird et al.
    T605 = bst.StorageTank('T605', ins=ethanol_fresh,
                           tau=7*24, V_wf=0.9,
                           vessel_type='Floating roof',
                           vessel_material='Carbon steel')
    T605.line = 'Ethanol storage tank'
    T605_P = bst.Pump('T605_P', ins=T605-0, P=101325)
    
    M401 = bst.units.Mixer('M401', ins=(filtered_fermentation_effluent, T604_P-0, None))
    
    @M401.add_specification(run=True)
    def adjust_M401_DPHP():
        feed, feed_DPHP, recycle_DPHP = all_ins = M401.ins[0:3]
        tot_mass = sum([i.F_mass - i.imass['DPHP'] - i.imass['Ethanol'] for i in all_ins])
        M401.ins[1].imass['Dipotassium hydrogen phosphate'] = max(0, 0.24 * tot_mass - recycle_DPHP.imass['DPHP'] - feed.imass['DPHP'])
    
    M401_P = bst.Pump('M401_P', ins=M401-0, outs='mixed_stream', P=101325)
    M401_P_H = bst.HXutility('M401_P_H', ins = M401_P-0, V = 0, rigorous = True)
    M401_b = bst.Mixer('M401_b', ins=(T605_P-0, ''))
    
    @M401_b.add_specification(run=True)
    def adjust_M401_b_ethanol():
        M401.run()
        feed_ethanol, recycle_ethanol  = M401_b.ins 
        feed, feed_DPHP, recycle_DPHP = all_ins = M401.ins[0:3] # yes, M401
        tot_mass = sum([i.F_mass - i.imass['DPHP'] - i.imass['Ethanol'] for i in all_ins])
        M401_b.ins[0].imass['Ethanol'] = max(0, 0.25 * tot_mass - recycle_ethanol.imass['Ethanol'])
    
    S402 = bst.units.MultiStageMixerSettlers('S402',
        ins = (M401_P_H-0, M401_b-0),
        partition_data={
            'K': np.array([1/28.34, 1/0.046, 1/10000, 10000, 1/28.34, 1/0.046,
                           1/0.046, 1/0.046, 1/0.046, 1/0.046, 1/0.046, 1/0.046, 1/0.046]),
            'IDs': ('2,3-Butanediol', 'Glucose', 'Ethanol', 'Water', 'Acetoin', 'Xylose',
                    'GlucoseOligomer', 'Extract', 'XyloseOligomer', 'Arabinose', 'ArabinoseOligomer', 'SolubleLignin',
                    'Enzyme'),
            'phi' : 0.5,
        },
        N_stages = 5
    )
    
    @S402.add_specification(run=True)
    def adjust_S402_split():
        M401_b.run()
    
    F401 = bst.units.Flash('F401', ins=S402-0, outs=('F401_g', 'F401_l'),
                            T = 379, P = 101325,
                            vessel_material = 'Stainless steel 316')
    
    # # Condense waste vapor for recycling
    F401_H = bst.units.HXutility('F401_H', ins=F401-0, V=0, rigorous=True)
    F401_P = bst.Pump('F401_P', ins=F401-1, P=101325)
    
    S403 = bst.Splitter('S403', ins = F401_P-0,
                        outs = ('DPHP_recycle', ''),
                        split = 0.98)
    
    M403 = bst.Mixer('M403', ins=(S403-1, F401_H-0), outs=wastewater)
    
    @S403.add_specification
    def adjust_S403_split():
        S403._run()
        separables = ('Extract', 'Arabinose', 'Xylose', 'XyloseOligomer',
                      'Mannose', 'Glucose', 'GlucoseOligomer', 'Galactose') # DPHP is insoluble in ethanol
        extract_mol = S403.ins[0].imol[separables]
        S403.outs[0].imol[separables] = 0
        S403.outs[1].imol[separables] = extract_mol
        
    # DPHP recycle to mixer-settler
    S403-0-2-M401
    
    F402 = bst.units.Flash('F402', ins=S402-1,
                                        outs=('F402_g', 'F402_l'),
                                        T = 368, P = 101325,
                                        vessel_material = 'Stainless steel 316')
    
    F402_H = bst.units.HXutility('F402_H', ins=F402-0, V=0, rigorous=True)
    F402_P = bst.Pump('F402_P', ins=F402-1, P=101325)
    
    D401 = bst.units.ShortcutColumn('D401', ins=F402_P-0,
                                        outs=('D401_g', 'D401_l'),
                                        LHK=('Ethanol', 'BDO'),
                                        partial_condenser=False,
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.9995, Hr=0.9995, k=1.2,
                                        vessel_material = 'Stainless steel 316')
    D401_P = bst.Pump('D401_P', ins=D401-1, P=101325)

    # ethanol recycle
    M402 = bst.units.Mixer('M402', ins = (F402_H-0, D401-0), outs = 'ethanol_mixed')
    M402_P = bst.Pump('M402_P', ins=M402-0, outs='ethanol_recycle', P=101325)
    M402_P-0-1-M401_b

    D402 = bst.units.ShortcutColumn('D402', ins=D401_P-0,
                                        outs=(unreacted_acetoin, 'D402_l'),
                                        LHK=('Acetoin', 'BDO'),
                                        partial_condenser=False,
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.9995, Hr=0.9995, k=1.2,
                                        vessel_material = 'Stainless steel 316')
    
    D402_P = bst.Pump('D402_P', ins=D402-1, outs=BDO, P=101325)


@SystemFactory(
    ID='BDO_separation_sys',
    ins=[dict(ID='filtered_fermentation_effluent', 
              phase='l', T=323.15, P=101325, 
              H2O=3663, AceticAcid=0.5861, 
              Glucose=27.61, BDO=191.1,
              GlucoseOligomer=6.5, Extract=62.82, 
              Xylose=33.96, XyloseOligomer=2.741, 
              Cellobiose=0.8368, Mannose=2.573, 
              MannoseOligomer=0.06861, Galactose=6.132, 
              GalactoseOligomer=0.1635, Arabinose=12.53, 
              ArabinoseOligomer=0.334, SolubleLignin=0.4022, 
              Protein=0.04538, Enzyme=23.95, 
              FermMicrobe=0.07988, Furfural=0.1058, 
              Acetoin=1.678, HMF=0.04494, 
              Glucan=0.003135, Mannan=3.481e-05, 
              Galactan=8.295e-05, Xylan=0.00139, 
              Arabinan=0.0001694, Lignin=0.00331, 
              Ash=0.02955, units='kmol/hr')],
    outs=[dict(ID='BDO'),
          dict(ID='unreacted_acetoin'),
          dict(ID='wastewater')],
)
def create_separation_system_oleyl_alcohol(ins, outs):
    """
    Create a separation system for BDO using ethanol and DPHP for the 
    "salting out" effect.

    Parameters
    ----------
    ins : stream
        Fermentation effluent.
    outs : stream sequence
        [0] BDO
        [1] Unreacted acetoin
        [2] Wastewater

    Examples
    --------
    >>> from biorefineries import BDO as bdo
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(bdo.BDO_chemicals)
    >>> BDO_separation_sys = bdo.create_separation_system_oleyl_alcohol()
    >>> BDO_separation_sys.simulate()
    >>> BDO_separation_sys.show()
    System: BDO_separation_sys
    ins...
    [0] filtered_fermentation_effluent
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                3.66e+03
                        AceticAcid         0.586
                        Glucose            27.6
                        2,3-Butanediol     191
                        GlucoseOligomer    6.5
                        Extract            62.8
                        Xylose             34
                        ...
    outs...
    [0] BDO
        phase: 'l', T: 407.71 K, P: 20265 Pa
        flow (kmol/hr): H2O                0.108
                        2,3-Butanediol     163
                        OleylAlcohol       0.0246
                        3-Hydroxybutanone  0.000715
    [1] unreacted_acetoin
        phase: 'g', T: 385.71 K, P: 20265 Pa
        flow (kmol/hr): 2,3-Butanediol     0.0815
                        3-Hydroxybutanone  1.43
    [2] wastewater
        phase: 'l', T: 343.2 K, P: 101325 Pa
        flow (kmol/hr): H2O                3.66e+03
                        AceticAcid         0.586
                        Glucose            27.6
                        2,3-Butanediol     28.1
                        OleylAlcohol       0.22
                        GlucoseOligomer    6.5
                        Extract            62.8
                        ...
    
    """
    filtered_fermentation_effluent, = ins
    BDO, unreacted_acetoin, wastewater = outs
    
    oleyl_alcohol = bst.Stream('oleyl_alcohol', price=price.get('OleylAlcohol', 0.))
    solvent_recycle = bst.Stream('')
    
    # 7-day storage time, similar to ethanol's in Humbird et al.
    T605 = bst.StorageTank('T605', ins=oleyl_alcohol,
                           tau=7*24, V_wf=0.9,
                           vessel_type='Floating roof',
                           vessel_material='Carbon steel')
    T605.line = 'Oleyl alcohol storage tank'
    T605_P = bst.Pump('T605_P', ins=T605-0, P=101325)
    
    M402 = bst.Mixer('M402', ins=(T605_P-0, solvent_recycle))
    
    preheated_stream = bst.Stream()
    D407 = bst.BinaryDistillation('D407', 
        ins=preheated_stream,
        LHK=('Water', 'BDO'),
        partial_condenser=False,
        k=1.1,
        # product_specification_format='Composition',
        # y_top=0.99999, x_bot=0.93)
        product_specification_format='Recovery',
        Lr=0.5, Hr=0.999)
    D407.target_BDO_x = 0.07
    def get_x(chem_ID, stream):
        return stream.imol[chem_ID]/sum(stream.imol[tuple([i.ID for i in stream.vle_chemicals])])
    def D407_f(Lr):
        D407.Hr = 0.999
        D407.Lr = Lr
        D407._run()
        BDO_x = get_x('BDO', D407.outs[1])
        return  get_x('BDO', D407.outs[1]) - max(get_x('BDO', D407.ins[0]), D407.target_BDO_x)
    D407.specification = bst.BoundedNumericalSpecification(D407_f, 0.001, 0.999)
    
    D407_Pb = bst.Pump('D407_Pb', D407-1, P=101325.)
    
    H407_b = bst.HXprocess('H407_b', 
        ins=[filtered_fermentation_effluent, D407_Pb-0],
        outs=[preheated_stream, ''],
    )
    
    S402 = bst.units.MultiStageMixerSettlers('S402',
        ins = (H407_b-1, M402-0),
        partition_data={
            'K': np.array([1./1.940224889932903, 1./0.16864361850509718, 1./0.37, 1./1.940224889932903, 1./10000.,
                            10000., 10000., 10000., 10000., 
                            10000., 10000., 10000., 10000.,
                            1./0.07147546181140499, 1./1.8386053937013191]),
            # 'K': np.array([1./1.5073253099767108, 1./0.16864361850509718, 1./0.37, 1./1.5073253099767108, 1./10000.,
            #                10000., 10000., 10000., 10000., 
            #                10000., 10000., 10000., 10000.,
            #                1./0.10987236183668542, 1./2.4303113434090213]),
            'IDs': ('2,3-Butanediol', 'Water', 'Ethanol', 'Acetoin', 'OleylAlcohol',
                    'Xylose', 'GlucoseOligomer', 'Extract', 'XyloseOligomer', 
                    'Arabinose', 'ArabinoseOligomer', 'SolubleLignin', 'Enzyme',
                    'Glycerol', 'Ethanol',),
            'phi' : 0.5,
        },
        N_stages = 20,
    )
    

    @S402.add_specification(run=True)
    def adjust_S402_split():
        feed = S402.ins[0]
        Water = feed.imass['Water']
        required_solvent = 8 * Water
        oleyl_alcohol, recycle = M402.ins
        oleyl_alcohol.imass['OleylAlcohol'] = max(0, required_solvent- recycle.imass['OleylAlcohol'])
        if recycle.imass['OleylAlcohol'] > required_solvent:
            recycle.imass['OleylAlcohol'] = required_solvent
        M402._run()
    
    S402_Pr = bst.Pump('S402_Pr', ins=S402-0, P=101325)
    S402_Pe = bst.Pump('S402_Pe', ins=S402-1, P=101325)
    
    D401_H = bst.HXprocess('D401_H', ins=[S402_Pe-0, None], outs=['', solvent_recycle], dT=15.)
    D401 = bst.units.BinaryDistillation('D401', ins=D401_H-0,
                                    outs=('D401_g', 'D401_l'),
                                    LHK=('BDO', 'OleylAlcohol'),
                                    partial_condenser=True,
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9999, k=1.1,
                                    P=0.06 * 101325,
                                    vessel_material = 'Stainless steel 316')
    D401_Pb = bst.Pump('D401_Pb', ins=D401-1, P=101325)
    D401_Pb-0-1-D401_H
    D402 = bst.units.ShortcutColumn('D402', ins=D401-0,
                                    outs=('D402_g', 'D402_l'),
                                    LHK=('Water', 'Acetoin'),
                                    partial_condenser=False,
                                    is_divided=True,
                                    P=0.2 * 101325,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9995, k=1.1,
                                    vessel_material = 'Stainless steel 316')
    
    D402_Pd = bst.Pump('D402_Pd', D402-0, P=101325)
    D402_Pb = bst.Pump('D402_Pb', D402-1, P=101325)
    
    M403 = bst.Mixer('M403', [S402_Pr-0, D402_Pd-0, D407-0], wastewater)
    D403x = bst.units.BinaryDistillation('D403x', ins=D402_Pb-0,
                                    outs=('D403x_g', 'D403x_l'),
                                    LHK=('Acetoin', 'BDO'),
                                    partial_condenser=False,
                                    is_divided=True,
                                    P=0.2 * 101325,
                                    product_specification_format='Recovery',
                                    Lr=0.995, Hr=0.999, k=1.2,
                                    vessel_material = 'Stainless steel 316')
    D403x_H = bst.HXutility('D403x_H', D403x-0, T=305.15, rigorous=True)
    D403x_Pd = bst.Pump('D403x_Pd', D403x_H-0, unreacted_acetoin, P=101325)
    D403x_Pb = bst.Pump('D403x_Pb', D403x-1, BDO, P=101325)


@SystemFactory(
    ID='BDO_separation_sys',
    ins=[dict(ID='filtered_fermentation_effluent', 
              phase='l', T=323.15, P=101325, 
              H2O=3663, AceticAcid=0.5861, 
              Glucose=27.61, BDO=191.1,
              GlucoseOligomer=6.5, Extract=62.82, 
              Xylose=33.96, XyloseOligomer=2.741, 
              Cellobiose=0.8368, Mannose=2.573, 
              MannoseOligomer=0.06861, Galactose=6.132, 
              GalactoseOligomer=0.1635, Arabinose=12.53, 
              ArabinoseOligomer=0.334, SolubleLignin=0.4022, 
              Protein=0.04538, Enzyme=23.95, 
              FermMicrobe=0.07988, Furfural=0.1058, 
              Acetoin=1.678, HMF=0.04494, 
              Glucan=0.003135, Mannan=3.481e-05, 
              Galactan=8.295e-05, Xylan=0.00139, 
              Arabinan=0.0001694, Lignin=0.00331, 
              Ash=0.02955, units='kmol/hr')],
    outs=[dict(ID='BDO'),
          dict(ID='unreacted_acetoin'),
          dict(ID='wastewater')],
)
def create_separation_system_oleyl_alcohol_3(ins, outs):
    """
    Create a separation system for BDO using ethanol and DPHP for the 
    "salting out" effect.

    Parameters
    ----------
    ins : stream
        Fermentation effluent.
    outs : stream sequence
        [0] BDO
        [1] Unreacted acetoin
        [2] Wastewater

    Examples
    --------
    >>> from biorefineries import BDO as bdo
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(bdo.BDO_chemicals)
    >>> BDO_separation_sys = bdo.create_separation_system_oleyl_alcohol()
    >>> BDO_separation_sys.simulate()
    >>> BDO_separation_sys.show()
    System: BDO_separation_sys
    ins...
    [0] filtered_fermentation_effluent
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): H2O                3.66e+03
                        AceticAcid         0.586
                        Glucose            27.6
                        2,3-Butanediol     191
                        GlucoseOligomer    6.5
                        Extract            62.8
                        Xylose             34
                        ...
    outs...
    [0] BDO
        phase: 'l', T: 407.71 K, P: 20265 Pa
        flow (kmol/hr): H2O                0.108
                        2,3-Butanediol     163
                        OleylAlcohol       0.0246
                        3-Hydroxybutanone  0.000715
    [1] unreacted_acetoin
        phase: 'g', T: 385.71 K, P: 20265 Pa
        flow (kmol/hr): 2,3-Butanediol     0.0815
                        3-Hydroxybutanone  1.43
    [2] wastewater
        phase: 'l', T: 343.2 K, P: 101325 Pa
        flow (kmol/hr): H2O                3.66e+03
                        AceticAcid         0.586
                        Glucose            27.6
                        2,3-Butanediol     28.1
                        OleylAlcohol       0.22
                        GlucoseOligomer    6.5
                        Extract            62.8
                        ...
    
    """
    
    
    filtered_fermentation_effluent, = ins
    BDO, unreacted_acetoin, wastewater = outs
    
    oleyl_alcohol = bst.Stream('oleyl_alcohol', price=price.get('OleylAlcohol', 0.))
    solvent_recycle = bst.Stream('')
    
    # 7-day storage time, similar to ethanol's in Humbird et al.
    T605 = bst.StorageTank('T605', ins=oleyl_alcohol,
                           tau=7*24, V_wf=0.9,
                           vessel_type='Floating roof',
                           vessel_material='Carbon steel')
    T605.line = 'Oleyl alcohol storage tank'
    T605_P = bst.Pump('T605_P', ins=T605-0, P=101325)
    
    M402 = bst.Mixer('M402', ins=(T605_P-0, solvent_recycle))
    
    preheated_stream = bst.Stream()
    D407 = bst.BinaryDistillation('D407', 
        ins=preheated_stream,
        LHK=('Water', 'BDO'),
        partial_condenser=False,
        k=1.1,
        # product_specification_format='Composition',
        # y_top=0.99999, x_bot=0.93)
        product_specification_format='Recovery',
        Lr=0.5, Hr=0.999)
    D407.target_BDO_x = 0.07
    def get_x(chem_ID, stream):
        return stream.imol[chem_ID]/sum(stream.imol[tuple([i.ID for i in stream.vle_chemicals])])
    def D407_f(Lr):
        D407.Hr = 0.999
        D407.Lr = Lr
        D407._run()
        BDO_x = get_x('BDO', D407.outs[1])
        return  get_x('BDO', D407.outs[1]) - max(get_x('BDO', D407.ins[0]), D407.target_BDO_x)
    D407.specification = bst.BoundedNumericalSpecification(D407_f, 0.001, 0.999)
    
    D407_Pb = bst.Pump('D407_Pb', D407-1, P=101325.)
    
    H407_b = bst.HXprocess('H407_b', 
        ins=[filtered_fermentation_effluent, D407_Pb-0],
        outs=[preheated_stream, ''],
    )
    
    S402 = bst.units.MultiStageMixerSettlers('S402',
        ins = (H407_b-1, M402-0),
        partition_data={
            'K': np.array([1./1.940224889932903, 1./0.16864361850509718, 1./0.37, 1./1.940224889932903, 1./10000.,
                            10000., 10000., 10000., 10000., 
                            10000., 10000., 10000., 10000.,
                            1./0.07147546181140499, 1./1.8386053937013191]),
            # 'K': np.array([1./1.5073253099767108, 1./0.16864361850509718, 1./0.37, 1./1.5073253099767108, 1./10000.,
            #                10000., 10000., 10000., 10000., 
            #                10000., 10000., 10000., 10000.,
            #                1./0.10987236183668542, 1./2.4303113434090213]),
            'IDs': ('2,3-Butanediol', 'Water', 'Ethanol', 'Acetoin', 'OleylAlcohol',
                    'Xylose', 'GlucoseOligomer', 'Extract', 'XyloseOligomer', 
                    'Arabinose', 'ArabinoseOligomer', 'SolubleLignin', 'Enzyme',
                    'Glycerol', 'Ethanol',),
            'phi' : 0.5,
        },
        N_stages = 20,
    )
    
    def get_K(chem_ID, stream, phase_1, phase_2):
        return (stream[phase_1].imol[chem_ID]/stream[phase_1].F_mol)/(stream[phase_2].imol[chem_ID]/stream[phase_2].F_mol)
    
    from biorefineries.BDO.BDO_interp import rbf
    def get_S402_Ks_interp(required_solvent_factor, feed, mixedstream=None):
        return rbf([[required_solvent_factor, feed.imass['Glycerol']/feed.imass[S402.relevant_IDs].sum()]])[0]
    
    def get_S402_Ks_rigor(required_solvent_factor, feed, mixedstream):
        try:
            mixedstream['l'].imass['OleylAlcohol'] = 0.5 * required_solvent_factor * feed.imass['Water']
            mixedstream['L'].imass['OleylAlcohol'] = 0.5 * required_solvent_factor * feed.imass['Water']
        except:
            mixedstream.imass['OleylAlcohol'] = required_solvent_factor * feed.imass['Water']
        mixedstream.lle(mixedstream.T, top_chemical='H2O')
        # return [get_K(i, mixedstream, 'l', 'L') if not i=='OleylAlcohol' else 1/4e5 for i in S402.relevant_IDs]
        # if mixedstream['l'].imol['OleylAlcohol'] <= mixedstream['L'].imol['OleylAlcohol']:
        #     import pdb
        #     pdb.set_trace()
        l = 'l' if mixedstream['l'].imol['OleylAlcohol'] >= mixedstream['L'].imol['OleylAlcohol'] else 'L'
        L = 'L' if l=='l' else 'l'
        return [get_K('Water', mixedstream, l, L), get_K('Glycerol', mixedstream, l, L),\
                get_K('BDO', mixedstream, l, L), get_K('OleylAlcohol', mixedstream, l, L)]
            
    def S402_f(required_solvent_factor):
        # D401._run()
        D401_Pb._run()
        D401_H._run()
        M402._run()
        
        relevant_IDs = S402.relevant_IDs
        feed = S402.ins[0]
        Water = feed.imass['Water']
        required_solvent = required_solvent_factor * Water
        oleyl_alcohol, recycle = M402.ins
        oleyl_alcohol.imass['OleylAlcohol'] = max(0, required_solvent- recycle.imass['OleylAlcohol'])
        if recycle.imass['OleylAlcohol'] > required_solvent:
            recycle.imass['OleylAlcohol'] = required_solvent
        M402._run()
        # new_Ks = get_S402_Ks_rigor(required_solvent_factor, feed, mixedstream)
        new_Ks = get_S402_Ks_interp(required_solvent_factor, feed)
        S402_Ks = S402.partition_data['K']
        # S402_IDs = S402.partition_data['IDs']
        
        # for i in range(len(new_Ks)):
        #     # print(S402_Ks[S402_IDs.index(relevant_IDs[i])], 1./new_Ks[i])
        #     # import pdb
        #     # pdb.set_trace()
        #     S402_Ks[S402_IDs.index(relevant_IDs[i])] = 1./new_Ks[i]
        
        S402_Ks[1], S402_Ks[13], S402_Ks[0], S402_Ks[4] = 1./new_Ks # ('Water', 'Glycerol', '2,3-Butanediol', 'OleylAlcohol')
        
        mixedstream.empty() # not needed but not time-consuming
        mixedstream.mix_from(S402.ins)
        # mixedstream.lle(mixedstream.T)
        phi = partition(mixedstream, mixedstream['l'], mixedstream['L'], IDs = relevant_IDs, 
                        K=np.array([S402_Ks[1], S402_Ks[13], S402_Ks[0], S402_Ks[4]]))
        extract_phase = 'l' if mixedstream['l'].imol['OleylAlcohol'] >= mixedstream['L'].imol['OleylAlcohol'] else 'L'
        raffinate_phase = 'L' if extract_phase=='l' else 'l'
        
        # S402._run()
        return mixedstream[extract_phase].imol['BDO']/(mixedstream[extract_phase].imol['BDO'] +\
            mixedstream[raffinate_phase].imol['BDO']) -\
            0.588
        
    def S402_spec():
        IQ_interpolation(S402_f, 7., 16., xtol=0.001, ytol=1e-4)
        S402._run()
    
    # S402._has_run = False
    # S402.specification = BoundedNumericalSpecification(S402_f, 7., 16., rtol)
    S402.specification = S402_spec
    
    mixedstream = tmo.Stream('mixedstream')
    mixedstream.phases = ('l', 'L')
    S402.relevant_IDs = ('Water', 'Glycerol', '2,3-Butanediol', 'OleylAlcohol')
    # @S402.add_specification(run=True)
    # def adjust_S402_split():
    #     relevant_IDs = S402.relevant_IDs
    #     feed = S402.ins[0]
    #     Water = feed.imass['Water']
    #     required_solvent = 8. * Water
    #     oleyl_alcohol, recycle = M402.ins
    #     oleyl_alcohol.imass['OleylAlcohol'] = max(0, required_solvent- recycle.imass['OleylAlcohol'])
    #     if recycle.imass['OleylAlcohol'] > required_solvent:
    #         # print(recycle.imass['OleylAlcohol'], required_solvent)
    #         recycle.imass['OleylAlcohol'] = required_solvent
    #     M402._run()
    
    S402_Pr = bst.Pump('S402_Pr', ins=S402-0, P=101325)
    S402_Pe = bst.Pump('S402_Pe', ins=S402-1, P=101325)
    ideal_thermo = S402.thermo.ideal()
    D401_H = bst.HXprocess('D401_H', ins=[S402_Pe-0, None], outs=['', solvent_recycle], dT=15.)
    D401 = bst.units.BinaryDistillation('D401', ins=D401_H-0,
                                    outs=('D401_g', 'D401_l'),
                                    LHK=('BDO', 'OleylAlcohol'),
                                    partial_condenser=True,
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9999, k=1.1,
                                    P=0.06 * 101325,
                                    vessel_material = 'Stainless steel 316',
                                    condenser_thermo = ideal_thermo)
    D401_Pb = bst.Pump('D401_Pb', ins=D401-1, P=101325)
    D401_Pb-0-1-D401_H
    D402 = bst.units.ShortcutColumn('D402', ins=D401-0,
                                    outs=('D402_g', 'D402_l'),
                                    LHK=('Water', 'Acetoin'),
                                    partial_condenser=False,
                                    is_divided=True,
                                    P=0.2 * 101325,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9995, k=1.1,
                                    vessel_material = 'Stainless steel 316')
    
    D402_Pd = bst.Pump('D402_Pd', D402-0, P=101325)
    D402_Pb = bst.Pump('D402_Pb', D402-1, P=101325)
    
    M403 = bst.Mixer('M403', [S402_Pr-0, D402_Pd-0, D407-0], wastewater)
    D403x = bst.units.BinaryDistillation('D403x', ins=D402_Pb-0,
                                    outs=('D403x_g', 'D403x_l'),
                                    LHK=('Acetoin', 'BDO'),
                                    partial_condenser=False,
                                    is_divided=True,
                                    P=0.2 * 101325,
                                    product_specification_format='Recovery',
                                    Lr=0.995, Hr=0.999, k=1.2,
                                    vessel_material = 'Stainless steel 316')
    D403x_H = bst.HXutility('D403x_H', D403x-0, T=305.15, rigorous=True)
    D403x_Pd = bst.Pump('D403x_Pd', D403x_H-0, unreacted_acetoin, P=101325)
    D403x_Pb = bst.Pump('D403x_Pb', D403x-1, BDO, P=101325)
    