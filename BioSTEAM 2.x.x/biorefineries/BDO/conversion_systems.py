1# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 03:54:05 2021

@author: yrc2
"""
import biosteam as bst
from biosteam import Stream
from biosteam.process_tools import SystemFactory
from biorefineries.BDO import units, facilities
from biorefineries.BDO.process_settings import price
from biorefineries.BDO.utils import find_split, splits_df, baseline_feedflow
from biorefineries.BDO.chemicals_data import BDO_chemicals, chemical_groups, \
                                soluble_organics, combustibles
import numpy as np


__all__ = (
    'create_conversion_system',    
)

@SystemFactory(
    ID='BDO_conversion_sys',
    ins=[dict(ID='BDO', phase='l', T=454.62, P=101325,
              H2O=9.469e-06, Ethanol=0.1061, Glucose=0.2162,
              BDO=184.3, GlucoseOligomer=0.0509, Extract=0.4919, 
              Xylose=0.2659, XyloseOligomer=0.02146, Arabinose=0.09811, 
              ArabinoseOligomer=0.08989, SolubleLignin=0.1082, 
              Enzyme=6.445, Acetoin=1.41, units='kmol/hr')],
    outs=[dict(ID='MEK', price=price['MEK']),
          dict(ID='isobutanol', price=price['Isobutanol']),
          dict(ID='wastewater')],
)
def create_conversion_system(ins, outs):
    """
    Create a conversion system for the production of MEK from renewable BDO 
    (from cellulosic biomass).

    Parameters
    ----------
    ins : stream
        Feedstock
    outs : stream sequence
        [0] Filtered fermentation effluent.
        [1] Cell mass
        [2] Solids
        [3] Wastewater

    Examples
    --------
    >>> from biorefineries import BDO
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(BDO.BDO_chemicals)
    >>> BDO_conversion_sys = create_conversion_system()
    >>> BDO_conversion_sys.simulate()
    >>> BDO_conversion_sys.show()
    System: BDO_conversion_sys
    ins...
    [0] BDO
        phase: 'l', T: 454.62 K, P: 101325 Pa
        flow (kmol/hr): H2O                9.47e-06
                        Ethanol            0.106
                        Glucose            0.216
                        2,3-Butanediol     184
                        GlucoseOligomer    0.0509
                        Extract            0.492
                        Xylose             0.266
                        ...
    outs...
    [0] MEK
        phase: 'l', T: 352.72 K, P: 101325 Pa
        flow (kmol/hr): H2O               0.0914
                        Ethanol           0.106
                        MEK               164
                        Isobutyraldehyde  0.00913
    [1] isobutanol
        phase: 'l', T: 372.42 K, P: 101325 Pa
        flow (kmol/hr): Ethanol           2.81e-05
                        MEK               0.0617
                        Isobutyraldehyde  0.0146
                        Isobutanol        17.9
    [2] wastewater
        phase: 'l', T: 858.35 K, P: 101325 Pa
        flow (kmol/hr): H2O                183
                        Ethanol            0.000315
                        Glucose            0.216
                        2,3-Butanediol     1.63
                        MEK                0.0822
                        Isobutyraldehyde   1.58e-12
                        GlucoseOligomer    0.0509
                        ...

    """
    BDO, = ins
    MEK, isobutanol, wastewater = outs
    
    H2_fresh = Stream('H2_fresh', price = price['H2'])
    
    R401 = units.DehydrationReactor('R401', ins = (BDO, ''),
                                    tau = 5,
                                    vessel_material='Stainless steel 316')
    
    D403 = bst.units.ShortcutColumn('D403', ins=R401-0,
                                        outs=('IBA', 'D403_l'),
                                        LHK=('IBA', 'MEK'),
                                        is_divided=True,
                                        partial_condenser=False,
                                        product_specification_format='Recovery',
                                        Lr=0.9995, Hr=0.9995, k=1.2,
                                        vessel_material = 'Stainless steel 316')
    D403_P = bst.Pump('D403_P', ins=D403-1, P=101325)
    ideal_thermo = D403.thermo.ideal()
    D404 = bst.units.ShortcutColumn('D404', ins=D403_P-0,
                                        outs=('D404_g', 'D404_l'),
                                        LHK=('MEK', 'H2O'),
                                        partial_condenser=False,
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        thermo=ideal_thermo,
                                        Lr=0.9995, Hr=0.9995, k=1.2,
                                        vessel_material = 'Stainless steel 316')
    D404_P = bst.Pump('D404_P', ins=D404-1, P=101325)
    
    D405 = bst.units.ShortcutColumn('D405', ins=D404_P-0,
                                        outs=('D405_g', 'D405_l'),
                                        LHK=('H2O', 'BDO'),
                                        partial_condenser=False,
                                        is_divided=True,
                                        P=0.2 * 101325,
                                        thermo=ideal_thermo,
                                        product_specification_format='Recovery',
                                        Lr=0.9995, Hr=0.9995, k=1.2,
                                        vessel_material = 'Stainless steel 316')
    
    D405_P = bst.Pump('D405_P', ins=D405-1, P=101325)
    
    S404 = bst.units.Splitter('S404', ins = D405_P-0, split = 0.92)
    S404-0-1-R401
    
    M404 = bst.Mixer('M404', ins=(D405-0, S404-1), outs=wastewater)
    
    R402 = units.HydrogenationReactor('R402', ins = (D403-0, '', ''),
                                    tau = 2,
                                    vessel_material = 'Stainless steel 316')
    
    D406 = bst.units.BinaryDistillation('D406', ins=R402-0,
                                        outs=('IBA', 'IBO'),
                                        LHK=('IBA', 'Isobutanol'),
                                        is_divided=True,
                                        partial_condenser=False,
                                        product_specification_format='Recovery',
                                        Lr=0.995, Hr=0.995, k=1.2,
                                        vessel_material = 'Stainless steel 316')
    
    D406_P = bst.Pump('D406_P', ins=D406-1, P=101325)
    
    D406-0-2-R402
    
    # 7-day storage time, similar to ethanol's in Humbird et al.
    T606 = units.BDOStorageTank('T606', ins=D404-0, tau=7*24, V_wf=0.9,
                                         vessel_type='Floating roof',
                                         vessel_material='Stainless steel')
    
    
    
    T606.line = 'MEKStorageTank'
    T606_P = bst.Pump('T606_P', ins=T606-0, outs=MEK, P=101325)
    
    T607 = units.HydrogenGasStorageTank('T607', ins=H2_fresh)
    T607.line = 'H2 storage tank'
    T607-0-1-R402
    
    T608 = units.BDOStorageTank('T608', ins=D406_P-0, tau=7*24, V_wf=0.9,
                                          vessel_type='Floating roof',
                                          vessel_material='Stainless steel')
    
    
    
    T608.line = 'IsobutanolStorageTank'
    T608_P = bst.Pump('T608_P', ins=T608-0, outs=isobutanol, P=101325)

