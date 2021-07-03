1# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 03:54:05 2021

@author: yrc2
"""
import biosteam as bst
from biosteam import Stream
from biosteam.process_tools import SystemFactory
from biorefineries import BDO as bdo
import biorefineries.cornstover as cs
import numpy as np

__all__ = (
    'create_system_DPHP',
    'create_system_oleyl_alcohol',
)

@bst.SystemFactory(
    ID='BDO_sys',
    ins=[dict(ID='feedstock', phase='l', T=298.15, P=101325,
              price=0.0628405632131775, H2O=2.084e+04, Acetate=1509,
              Extract=1.221e+04, Sucrose=641.8, Protein=2584, 
              Glucan=2.921e+04, Mannan=500.1, Galactan=1192, 
              Xylan=1.628e+04, Arabinan=1984, Lignin=1.314e+04, 
              Ash=4109, units='kg/hr')],
    outs=[dict(ID='MEK', price=bdo.price['MEK']),
          dict(ID='isobutanol', price=bdo.price['Isobutanol']),],
)
def create_system_oleyl_alcohol(ins, outs):
    """
    Create a system for the production of BDO from cellulosic biomass.

    Parameters
    ----------
    ins : stream
        Feedstock
    outs : stream sequence
        [0] MEK.
        [1] isobutanol

    Examples
    --------
    >>> from biorefineries import BDO as bdo
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(bdo.BDO_chemicals)
    >>> BDO_sys = bdo.create_system_oleyl_alcohol()
    >>> BDO_sys.simulate()
    >>> BDO_sys.show()
    System: BDO_sys
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
    outs...
    [0] MEK
        phase: 'l', T: 352.74 K, P: 101325 Pa
        flow (kmol/hr): H2O               0.0937
                        MEK               168
                        Isobutyraldehyde  0.00936
    [1] isobutanol
        phase: 'l', T: 372.42 K, P: 101325 Pa
        flow (kmol/hr): MEK               0.0632
                        Isobutyraldehyde  0.0149
                        Isobutanol        18.3
    """
    feedstock,  = ins
    MEK, isobutanol = outs
    f = bst.main_flowsheet
    s = f.stream
    u = f.unit
    recycle_acetoin = bst.Stream()
    pretreatment_and_fermentation_sys = bdo.create_pretreatment_and_fermentation_system(
        ins=[feedstock, recycle_acetoin],
        mockup=True,
    )
    filtered_fermentation_effluent, solids, wastewater_a = pretreatment_and_fermentation_sys.outs
    separation_sys = bdo.create_separation_system_oleyl_alcohol(
        ins=filtered_fermentation_effluent,
        outs=['', recycle_acetoin, ''],
        mockup=True,
    )
    BDO, recycle_acetoin, wastewater_b = separation_sys.outs
    conversion_sys = bdo.create_conversion_system(
        ins=BDO,
        outs=[MEK, isobutanol, ''],
        mockup=True,
    )
    MEK, isobutanol, wastewater_c = conversion_sys.outs
    # u.S404.split = 0.6
    # M405 = bst.Mixer('M405', ins=u.S404-1)
    # M405.insert(u.S402_Pe-0)
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater_a, wastewater_b, wastewater_c],
        mockup=True,
    )
    methane, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
    solids.source.ins.append(sludge)
    cs.create_facilities(
        solids_to_boiler=solids,
        gas_to_boiler=methane,
        process_water_streams=(s.caustic, 
                               s.water_M201, 
                               s.water_M202, 
                               s.water_M203, 
                               s.water_M205, 
                               s.enzyme_water),
        feedstock=feedstock,
        RO_water=treated_water,
    )
    # u.M305.ins[0] = None
    # u.M305.ins[1] = None
    # u.M601.ins[1] = None
    HXN = bst.facilities.HeatExchangerNetwork('HXN', 
        ignored=[u.D404.boiler, u.D407.boiler]
    )
    # def HXN_no_run_cost():
    #     HXN.heat_utilities = tuple()
    #     HXN._installed_cost = 0.
    
    # # To simulate without HXN, uncomment the following 3 lines:
    # HXN._cost = HXN_no_run_cost
    # HXN.energy_balance_percent_error = 0.
    # HXN.new_HXs = HXN.new_HX_utils = []

@bst.SystemFactory(
    ID='BDO_sys',
    ins=[dict(ID='feedstock', phase='l', T=298.15, P=101325,
              price=0.0628405632131775, H2O=2.084e+04, Acetate=1509,
              Extract=1.221e+04, Sucrose=641.8, Protein=2584, 
              Glucan=2.921e+04, Mannan=500.1, Galactan=1192, 
              Xylan=1.628e+04, Arabinan=1984, Lignin=1.314e+04, 
              Ash=4109, units='kg/hr')],
    outs=[dict(ID='MEK', price=bdo.price['MEK']),
          dict(ID='isobutanol', price=bdo.price['Isobutanol']),],
)
def create_system_DPHP(ins, outs):
    """
    Create a system for the production of BDO from cellulosic biomass.

    Parameters
    ----------
    ins : stream
        Feedstock
    outs : stream sequence
        [0] MEK.
        [1] isobutanol

    Examples
    --------
    >>> from biorefineries import BDO as bdo
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(bdo.BDO_chemicals)
    >>> BDO_sys = bdo.create_system_DPHP()
    >>> BDO_sys.simulate()
    >>> BDO_sys.show()
    System: BDO_sys
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
    outs...
    [0] MEK
        phase: 'l', T: 352.74 K, P: 101325 Pa
        flow (kmol/hr): H2O               0.0937
                        MEK               168
                        Isobutyraldehyde  0.00936
    [1] isobutanol
        phase: 'l', T: 372.42 K, P: 101325 Pa
        flow (kmol/hr): MEK               0.0632
                        Isobutyraldehyde  0.0149
                        Isobutanol        18.3
    
    """
    feedstock,  = ins
    MEK, isobutanol = outs
    f = bst.main_flowsheet
    s = f.stream
    u = f.unit
    recycle_acetoin = bst.Stream()
    pretreatment_and_fermentation_sys = bdo.create_pretreatment_and_fermentation_system(
        ins=[feedstock, recycle_acetoin],
        mockup=True,
    )
    filtered_fermentation_effluent, solids, wastewater_a = pretreatment_and_fermentation_sys.outs
    separation_sys = bdo.create_separation_system_DPHP(
        ins=filtered_fermentation_effluent,
        outs=['', recycle_acetoin, ''],
        mockup=True,
    )
    BDO, recycle_acetoin, wastewater_b = separation_sys.outs
    conversion_sys = bdo.create_conversion_system(
        ins=BDO,
        outs=[MEK, isobutanol, ''],
        mockup=True,
    )
    MEK, isobutanol, wastewater_c = conversion_sys.outs
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater_a, wastewater_b, wastewater_c],
        mockup=True,
    )
    methane, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
    solids.source.ins.append(sludge)
    cs.create_facilities(
        solids_to_boiler=solids,
        gas_to_boiler=methane,
        process_water_streams=(s.caustic, 
                               s.water_M201, 
                               s.water_M202, 
                               s.water_M203, 
                               s.water_M205, 
                               s.enzyme_water),
        feedstock=feedstock,
        RO_water=treated_water,
    )
    # u.M305.ins[0] = None
    # u.M305.ins[1] = None
    HXN = bst.facilities.HeatExchangerNetwork('HXN')
    # def HXN_no_run_cost():
    #     HXN.heat_utilities = tuple()
    #     HXN._installed_cost = 0.
    
    # # To simulate without HXN, uncomment the following 3 lines:
    # HXN._cost = HXN_no_run_cost
    # HXN.energy_balance_percent_error = 0.
    # HXN.new_HXs = HXN.new_HX_utils = []
