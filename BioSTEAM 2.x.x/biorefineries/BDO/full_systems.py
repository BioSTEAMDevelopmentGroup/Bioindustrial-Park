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
    'create_system_oleyl_alcohol',
    'create_system_broth',
)

def rename_storage_units(units, number):
    
    def is_storage_unit(unit):
        return (
            ('storage' in unit.line.lower() 
             or isinstance(unit, bst.StorageTank)
             or 'storage' in unit.__class__.__name__.lower()) 
        )
    
    bst.rename_units([i for i in units if is_storage_unit(i)], number)


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
    filtered_fermentation_effluent, vent, solids, wastewater_a = pretreatment_and_fermentation_sys.outs
    separation_sys = bdo.create_separation_system_oleyl_alcohol(
        ins=filtered_fermentation_effluent,
        outs=['', recycle_acetoin, ''],
        mockup=True,
    )
    BDO, recycle_acetoin, wastewater_b = separation_sys.outs
    conversion_sys = bdo.create_conversion_system(
        ins=BDO,
        outs=[MEK, isobutanol, '',],
        mockup=True,
    )
    MEK, isobutanol, wastewater_c, spent_TCP_catalyst, spent_KieCNi_catalyst = conversion_sys.outs
    # u.S404.split = 0.6
    # M405 = bst.Mixer('M405', ins=u.S404-1)
    # M405.insert(u.S402_Pe-0)
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater_a, wastewater_b, wastewater_c],
        mockup=True,
        area=500,
    )
    methane, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
    MX = bst.Mixer(700, [vent, methane])
    solids.source.ins.append(sludge)
    cs.create_facilities(
        solids_to_boiler=solids,
        gas_to_boiler=MX-0,
        process_water_streams=(s.caustic, 
                               s.water_M201, 
                               s.water_M202, 
                               s.water_M203, 
                               s.enzyme_water),
        feedstock=feedstock,
        RO_water=treated_water,
    )
    # u.M305.ins[0] = None
    # u.M305.ins[1] = None
    # u.M601.ins[1] = None
    HXN = bst.facilities.HeatExchangerNetwork(1000, 
        ignored=[u.D404.boiler, u.D406, u.D407.boiler]
    )
    u.FT.ID = 900
    BT = u.BT
    CT = u.CT
    BT.ID = 700
    CT.ID = 800
    excluded = (HXN, BT, CT)
    for i in list(u):
        if isinstance(i, bst.Facility) and i not in excluded:
            i.ID = 900
    rename_storage_units(u, 600)
    # def HXN_no_run_cost():
    #     HXN.heat_utilities = tuple()
    #     HXN._installed_cost = 0.
    
    # # To simulate without HXN, uncomment the following 3 lines:
    # HXN._cost = HXN_no_run_cost
    # HXN.energy_balance_percent_error = 0.
    # HXN.new_HXs = HXN.new_HX_utils = []

# @bst.SystemFactory(
#     ID='BDO_sys',
#     ins=[dict(ID='feedstock', phase='l', T=298.15, P=101325,
#               price=0.0628405632131775, H2O=2.084e+04, Acetate=1509,
#               Extract=1.221e+04, Sucrose=641.8, Protein=2584, 
#               Glucan=2.921e+04, Mannan=500.1, Galactan=1192, 
#               Xylan=1.628e+04, Arabinan=1984, Lignin=1.314e+04, 
#               Ash=4109, units='kg/hr')],
#     outs=[dict(ID='MEK', price=bdo.price['MEK']),
#           dict(ID='isobutanol', price=bdo.price['Isobutanol']),],
# )
# def create_system_DPHP(ins, outs):
#     """
#     Create a system for the production of BDO from cellulosic biomass.

#     Parameters
#     ----------
#     ins : stream
#         Feedstock
#     outs : stream sequence
#         [0] MEK.
#         [1] isobutanol

#     Examples
#     --------
#     >>> from biorefineries import BDO as bdo
#     >>> import biosteam as bst
#     >>> bst.settings.set_thermo(bdo.BDO_chemicals)
#     >>> BDO_sys = bdo.create_system_DPHP()
#     >>> BDO_sys.simulate()
#     >>> BDO_sys.show()
#     System: BDO_sys
#     ins...
#     [0] feedstock
#         phase: 'l', T: 298.15 K, P: 101325 Pa
#         flow (kmol/hr): H2O       1.16e+03
#                         Acetate   25.1
#                         Extract   67.8
#                         Sucrose   1.87
#                         Protein   113
#                         Glucan    180
#                         Mannan    3.08
#                         ...
#     outs...
#     [0] MEK
#         phase: 'l', T: 352.74 K, P: 101325 Pa
#         flow (kmol/hr): H2O               0.0937
#                         MEK               168
#                         Isobutyraldehyde  0.00936
#     [1] isobutanol
#         phase: 'l', T: 372.42 K, P: 101325 Pa
#         flow (kmol/hr): MEK               0.0632
#                         Isobutyraldehyde  0.0149
#                         Isobutanol        18.3
    
#     """
#     feedstock,  = ins
#     MEK, isobutanol = outs
#     f = bst.main_flowsheet
#     s = f.stream
#     u = f.unit
#     recycle_acetoin = bst.Stream()
#     pretreatment_and_fermentation_sys = bdo.create_pretreatment_and_fermentation_system(
#         ins=[feedstock, recycle_acetoin],
#         mockup=True,
#     )
#     filtered_fermentation_effluent, solids, wastewater_a = pretreatment_and_fermentation_sys.outs
#     separation_sys = bdo.create_separation_system_DPHP(
#         ins=filtered_fermentation_effluent,
#         outs=['', recycle_acetoin, ''],
#         mockup=True,
#     )
#     BDO, recycle_acetoin, wastewater_b = separation_sys.outs
#     conversion_sys = bdo.create_conversion_system(
#         ins=BDO,
#         outs=[MEK, isobutanol, ''],
#         mockup=True,
#     )
#     MEK, isobutanol, wastewater_c = conversion_sys.outs
    
#     wastewater_treatment_sys = bst.create_wastewater_treatment_system(
#         ins=[wastewater_a, wastewater_b, wastewater_c],
#         area=500,
#         mockup=True,
#     )
#     methane, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
#     solids.source.ins.append(sludge)
#     cs.create_facilities(
#         solids_to_boiler=solids,
#         gas_to_boiler=methane,
#         process_water_streams=(s.caustic, 
#                                s.water_M201, 
#                                s.water_M202, 
#                                s.water_M203, 
#                                s.water_M205, 
#                                s.enzyme_water),
#         feedstock=feedstock,
#         RO_water=treated_water,
#     )
#     # u.M305.ins[0] = None
#     # u.M305.ins[1] = None
#     HXN = bst.facilities.HeatExchangerNetwork(1000)
#     u.FT.ID = 900
#     u.BT.ID = 700
#     u.CT.ID = 800
#     excluded = (HXN, u.BT, u.CT)
#     for i in u:
#         if isinstance(i, bst.Facility) and i not in excluded:
#             i.ID = 900
#     rename_storage_units(u, 600)
#     # def HXN_no_run_cost():
#     #     HXN.heat_utilities = tuple()
#     #     HXN._installed_cost = 0.
    
#     # # To simulate without HXN, uncomment the following 3 lines:
#     # HXN._cost = HXN_no_run_cost
#     # HXN.energy_balance_percent_error = 0.
#     # HXN.new_HXs = HXN.new_HX_utils = []

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
    outs=[dict(ID='conc_aqueous_broth'),
          dict(ID='wastewater_b')],
)
def create_concentration_evaporator_sys(ins, outs):
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
    conc_aqueous_broth, wastewater_b = outs
    
    F401 = bst.units.MultiEffectEvaporator('F401', ins=filtered_fermentation_effluent, outs=('F401_0', wastewater_b),
        P = (101325, 73581, 50892, 32777, 20000), V = 0.5)
    F401.target_BDO_x = 250e-6
    def get_x(chem_ID, stream):
        return stream.imol[chem_ID]/sum(stream.imol[tuple([i.ID for i in stream.vle_chemicals])])
    
    def F401_specification():
        instream = F401.ins[0]
        # ratio = target_water_x/get_x('Water', instream)
        ratio = get_x('BDO', instream)/F401.target_BDO_x
        # no need to check for ratio>1 because our target_water_x is consistently lower than the max possible titer
        F401.V = max(1e-5, 1. - ratio)
        F401._run()
    
    def F401_cost():
        bst.MultiEffectEvaporator._cost()
        if F401.V == 1e-5:
            F401_hu = F401.heat_utilities
            for hu in F401_hu:
                hu.flow, hu.duty, hu.cost = 0., 0., 0.
                hu.heat_exchanger = None
            F401_bpc = F401.baseline_purchase_costs
            F401_bpc_keys = F401_bpc.keys()
            for k in F401_bpc_keys:
                F401_bpc[k] = 0.
    F401._cost = F401_cost
    F401.specification = F401_specification
    
    M401 = bst.Mixer('M401', ins=(F401-0, 'dilution_water'), outs=('diluted_broth',))
    M401.target_BDO_x = 250e-6
    def M401_f(dilution_water_mol):
        if get_x('BDO', M401.ins[0]) >= M401.target_BDO_x:
            M401.ins[1].imol['Water'] = dilution_water_mol
            M401._run()
            return get_x('BDO', M401.outs[0]) -  M401.target_BDO_x
        else:
            M401.ins[1].imol['Water'] = 1e-6
            M401._run()
            return 0.
    M401.specification = bst.BoundedNumericalSpecification(M401_f, 1e-6, 1e8)
    
    M401_H = bst.HXutility('M401_H', ins=M401-0, outs='cooled_broth', T=25.+273.15, rigorous=True)
    # wastewater_b = F401.outs[1]
    T601 = bst.StorageTank('T601', ins=M401_H-0,
                           tau=7*24, V_wf=0.9,
                           vessel_type='Floating roof',
                           vessel_material='Carbon steel')
    T601.line = 'Conc. aq. broth storage tank'
    T601_P = bst.Pump('T601_P', ins=T601-0, outs=(conc_aqueous_broth,), P=101325)
    # conc_aqueous_broth = T601_P-0
    
@bst.SystemFactory(
    ID='BDO_broth_sys',
    ins=[dict(ID='feedstock', phase='l', T=298.15, P=101325,
              price=0.0628405632131775, H2O=2.084e+04, Acetate=1509,
              Extract=1.221e+04, Sucrose=641.8, Protein=2584, 
              Glucan=2.921e+04, Mannan=500.1, Galactan=1192, 
              Xylan=1.628e+04, Arabinan=1984, Lignin=1.314e+04, 
              Ash=4109, units='kg/hr')],
    outs=[dict(ID='conc_aqueous_broth', price=bdo.price['MEK'])],
)
def create_system_broth(ins, outs):
    feedstock,  = ins
    conc_aqueous_broth, = outs
    f = bst.main_flowsheet
    s = f.stream
    u = f.unit
    recycle_acetoin = bst.Stream()
    pretreatment_and_fermentation_sys = bdo.create_pretreatment_and_fermentation_system(
        ins=[feedstock, recycle_acetoin],
        mockup=True,
    )
    filtered_fermentation_effluent, vent, solids, wastewater_a = pretreatment_and_fermentation_sys.outs
    
    concentration_sys = create_concentration_evaporator_sys(
        ins=filtered_fermentation_effluent,
        outs=(conc_aqueous_broth, ''),
        mockup=True,
    )
   
    recycle_acetoin.imol['Acetoin'] = 1e-3
    # conc_aqueous_broth, wastewater_b = concentration_sys.outs
   
    # conc_aqueous_broth_2 = conc_aqueous_broth
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater_a, concentration_sys-1],
        mockup=True,
        area=500,
    )
    methane, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
    MX = bst.Mixer(700, [vent, methane])
    solids.source.ins.append(sludge)
    cs.create_facilities(
        solids_to_boiler=solids,
        gas_to_boiler=MX-0,
        process_water_streams=(s.caustic, 
                               s.water_M201, 
                               s.water_M202, 
                               s.water_M203, 
                               s.enzyme_water),
        feedstock=feedstock,
        RO_water=treated_water,
    )
    # u.M305.ins[0] = None
    # u.M305.ins[1] = None
    # u.M601.ins[1] = None
    HXN = bst.facilities.HeatExchangerNetwork(1000, 
        # ignored=[u.D404.boiler, u.D406, u.D407.boiler]
    )
    u.FT.ID = 900
    BT = u.BT
    CT = u.CT
    BT.ID = 700
    CT.ID = 800
    excluded = (HXN, BT, CT)
    for i in list(u):
        if isinstance(i, bst.Facility) and i not in excluded:
            i.ID = 900
    rename_storage_units(u, 600)



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
def create_system_oleyl_alcohol_3(ins, outs):
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
    filtered_fermentation_effluent, vent, solids, wastewater_a = pretreatment_and_fermentation_sys.outs
    separation_sys = bdo.create_separation_system_oleyl_alcohol_3(
        ins=filtered_fermentation_effluent,
        outs=['', recycle_acetoin, ''],
        mockup=True,
    )
    BDO, recycle_acetoin, wastewater_b = separation_sys.outs
    conversion_sys = bdo.create_conversion_system(
        ins=BDO,
        outs=[MEK, isobutanol, '',],
        mockup=True,
    )
    MEK, isobutanol, wastewater_c, spent_TCP_catalyst, spent_KieCNi_catalyst = conversion_sys.outs
    # u.S404.split = 0.6
    # M405 = bst.Mixer('M405', ins=u.S404-1)
    # M405.insert(u.S402_Pe-0)
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater_a, wastewater_b, wastewater_c],
        mockup=True,
        area=500,
    )
    methane, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
    MX = bst.Mixer(700, [vent, methane])
    solids.source.ins.append(sludge)
    cs.create_facilities(
        solids_to_boiler=solids,
        gas_to_boiler=MX-0,
        process_water_streams=(s.caustic, 
                               s.water_M201, 
                               s.water_M202, 
                               s.water_M203, 
                               s.enzyme_water),
        feedstock=feedstock,
        RO_water=treated_water,
    )
    # u.M305.ins[0] = None
    # u.M305.ins[1] = None
    # u.M601.ins[1] = None
    HXN = bst.facilities.HeatExchangerNetwork(1000, 
        ignored=[u.D404.boiler, u.D406, u.D407.boiler]
    )
    u.FT.ID = 900
    BT = u.BT
    CT = u.CT
    BT.ID = 700
    CT.ID = 800
    excluded = (HXN, BT, CT)
    for i in list(u):
        if isinstance(i, bst.Facility) and i not in excluded:
            i.ID = 900
    rename_storage_units(u, 600)
    # def HXN_no_run_cost():
    #     HXN.heat_utilities = tuple()
    #     HXN._installed_cost = 0.
    
    # # To simulate without HXN, uncomment the following 3 lines:
    # HXN._cost = HXN_no_run_cost
    # HXN.energy_balance_percent_error = 0.
    # HXN.new_HXs = HXN.new_HX_utils = []

# @bst.SystemFactory(
#     ID='BDO_sys',
#     ins=[dict(ID='feedstock', phase='l', T=298.15, P=101325,
#               price=0.0628405632131775, H2O=2.084e+04, Acetate=1509,
#               Extract=1.221e+04, Sucrose=641.8, Protein=2584, 
#               Glucan=2.921e+04, Mannan=500.1, Galactan=1192, 
#               Xylan=1.628e+04, Arabinan=1984, Lignin=1.314e+04, 
#               Ash=4109, units='kg/hr')],
#     outs=[dict(ID='MEK', price=bdo.price['MEK']),
#           dict(ID='isobutanol', price=bdo.price['Isobutanol']),],
# )
# def create_system_DPHP(ins, outs):
#     """
#     Create a system for the production of BDO from cellulosic biomass.

#     Parameters
#     ----------
#     ins : stream
#         Feedstock
#     outs : stream sequence
#         [0] MEK.
#         [1] isobutanol

#     Examples
#     --------
#     >>> from biorefineries import BDO as bdo
#     >>> import biosteam as bst
#     >>> bst.settings.set_thermo(bdo.BDO_chemicals)
#     >>> BDO_sys = bdo.create_system_DPHP()
#     >>> BDO_sys.simulate()
#     >>> BDO_sys.show()
#     System: BDO_sys
#     ins...
#     [0] feedstock
#         phase: 'l', T: 298.15 K, P: 101325 Pa
#         flow (kmol/hr): H2O       1.16e+03
#                         Acetate   25.1
#                         Extract   67.8
#                         Sucrose   1.87
#                         Protein   113
#                         Glucan    180
#                         Mannan    3.08
#                         ...
#     outs...
#     [0] MEK
#         phase: 'l', T: 352.74 K, P: 101325 Pa
#         flow (kmol/hr): H2O               0.0937
#                         MEK               168
#                         Isobutyraldehyde  0.00936
#     [1] isobutanol
#         phase: 'l', T: 372.42 K, P: 101325 Pa
#         flow (kmol/hr): MEK               0.0632
#                         Isobutyraldehyde  0.0149
#                         Isobutanol        18.3
    
#     """
#     feedstock,  = ins
#     MEK, isobutanol = outs
#     f = bst.main_flowsheet
#     s = f.stream
#     u = f.unit
#     recycle_acetoin = bst.Stream()
#     pretreatment_and_fermentation_sys = bdo.create_pretreatment_and_fermentation_system(
#         ins=[feedstock, recycle_acetoin],
#         mockup=True,
#     )
#     filtered_fermentation_effluent, solids, wastewater_a = pretreatment_and_fermentation_sys.outs
#     separation_sys = bdo.create_separation_system_DPHP(
#         ins=filtered_fermentation_effluent,
#         outs=['', recycle_acetoin, ''],
#         mockup=True,
#     )
#     BDO, recycle_acetoin, wastewater_b = separation_sys.outs
#     conversion_sys = bdo.create_conversion_system(
#         ins=BDO,
#         outs=[MEK, isobutanol, ''],
#         mockup=True,
#     )
#     MEK, isobutanol, wastewater_c = conversion_sys.outs
    
#     wastewater_treatment_sys = bst.create_wastewater_treatment_system(
#         ins=[wastewater_a, wastewater_b, wastewater_c],
#         area=500,
#         mockup=True,
#     )
#     methane, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
#     solids.source.ins.append(sludge)
#     cs.create_facilities(
#         solids_to_boiler=solids,
#         gas_to_boiler=methane,
#         process_water_streams=(s.caustic, 
#                                s.water_M201, 
#                                s.water_M202, 
#                                s.water_M203, 
#                                s.water_M205, 
#                                s.enzyme_water),
#         feedstock=feedstock,
#         RO_water=treated_water,
#     )
#     # u.M305.ins[0] = None
#     # u.M305.ins[1] = None
#     HXN = bst.facilities.HeatExchangerNetwork(1000)
#     u.FT.ID = 900
#     u.BT.ID = 700
#     u.CT.ID = 800
#     excluded = (HXN, u.BT, u.CT)
#     for i in u:
#         if isinstance(i, bst.Facility) and i not in excluded:
#             i.ID = 900
#     rename_storage_units(u, 600)
#     # def HXN_no_run_cost():
#     #     HXN.heat_utilities = tuple()
#     #     HXN._installed_cost = 0.
    
#     # # To simulate without HXN, uncomment the following 3 lines:
#     # HXN._cost = HXN_no_run_cost
#     # HXN.energy_balance_percent_error = 0.
#     # HXN.new_HXs = HXN.new_HX_utils = []
