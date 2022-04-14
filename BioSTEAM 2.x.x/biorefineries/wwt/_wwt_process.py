#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
Unit construction and functions for creating wastewater treatment system.

References
----------
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of
Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
National Renewable Energy Lab (NREL), 2011.
https://www.nrel.gov/docs/fy11osti/47764.pdf

[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
https://doi.org/10.2172/1483234
'''


# %%

import thermosteam as tmo, biosteam as bst
from biosteam import Unit
from biosteam.units.decorators import cost
from . import (
    default_insolubles,
    InternalCirculationRx, AnMBR, PolishingFilter, BeltThickener, SludgeCentrifuge,
    get_combustion_energy, prices,
    )

_mgd_to_cmh = 157.7255 # auom('gallon').conversion_factor('m3')*1e6/24
_gpm_to_cmh = 0.2271 # auom('gallon').conversion_factor('m3')*60
_Gcal_to_kJ = 4184000 # auom('kcal').conversion_factor('kJ')*1e6 # (also MMkcal/hr)
_kW_to_kJhr = 3600 # auom('kW').conversion_factor('kJ/hr')

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
CEPCI = bst.units.design_tools.CEPCI_by_year

__all__ = ('create_wastewater_process', 'CHP', 'Skipped',)


# %%

# =============================================================================
# Other units
# =============================================================================

# # The scaling basis of BeltThickener and Centrifuge changed significantly
# # from previous report to this current one (ref [2]),
# # therefore using alternative designs
# @cost(basis='COD flow', ID='Thickeners', units='kg-O2/hr',
#       kW=107.3808, cost=750000, S=5600, CE=CEPCI[2012], n=0.6, BM=1.6)
# class BeltThickener(Unit):
#     _ins_size_is_fixed = False
#     _N_outs = 2
#     _units= {'COD flow': 'kg-O2/hr'}

#     def __init__(self, ID='', ins=None, outs=(), thermo=None,
#                  insolubles=default_insolubles):
#         Unit.__init__(self, ID, ins, outs, thermo)
#         self.insolubles = get_insoluble_IDs(self.chemicals, insolubles)

#     def _run(self):
#         centrate, solids = self.outs
#         insolubles = self.insolubles
#         solubles = get_soluble_IDs(self.chemicals, insolubles)

#         influent = self.ins[0].copy()
#         influent.mix_from(self.ins)

#         solids.copy_flow(influent, insolubles)
#         # Concentrate sludge to 4% solids
#         solids.imass['Water'] = 0.96/0.04 * influent.imass[insolubles].sum()
#         if solids.imass['Water'] < influent.imass['Water']:
#             ratio = solids.imass['Water'] / influent.imass['Water']
#             solids.imass[solubles] = ratio * influent.imass[solubles]
#             solids.T = influent.T

#             centrate.mol = influent.mol - solids.mol
#             centrate.T = influent.T
#         else:
#             centrate.empty()
#             solids.copy_like(influent)

#         self._inf = influent


#     def _design(self):
#         self.design_results['COD flow'] = compute_stream_COD(self._inf)


# @cost(basis='COD flow', ID='Centrifuge', units='kg-O2/hr',
#       # power usage includes feed pumping and centrifuge
#       kW=22.371+123.0405, cost=686800, S=5600, CE=CEPCI[2012], n=0.6, BM=2.7)
# class SludgeCentrifuge(Unit):
#     _N_ins = 1
#     _N_outs = 2
#     _units= {'COD flow': 'kg-O2/hr'}

#     __init__ = BeltThickener.__init__

#     def _run(self):
#         influent = self.ins[0]
#         centrate, solids = self.outs
#         centrate.T = solids.T = influent.T
#         insolubles = self.insolubles
#         solubles = get_soluble_IDs(self.chemicals, insolubles)

#         # Centrifuge captures 95% of the solids at 20% solids
#         solids.imass[insolubles] = 0.95 * influent.imass[insolubles]
#         solids.imass['Water'] = 0.8/0.2 * (influent.imass[insolubles].sum())
#         if solids.imass['Water'] < influent.imass['Water']:
#             ratio = solids.imass['Water'] / influent.imass['Water']
#             solids.imass[solubles] = ratio * influent.imass[solubles]

#             centrate.mol = influent.mol - solids.mol
#         else:
#             centrate.empty()
#             solids.copy_like(influent)

#         self._inf = influent


#     _design = BeltThickener._design


@cost(basis='Volumetric flow', ID='Reactor', units='m3/hr',
      # 2.7 in million gallons per day (MGD)
      cost=2450000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=1, BM=1.8)
@cost(basis='Volumetric flow', ID='Evaporator', units='m3/hr',
      # 2.7 in million gallons per day (MGD)
       kW=1103.636, cost=5000000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=0.6, BM=1.6)
class ReverseOsmosis(Unit):
    _N_ins = 1
    _N_outs = 2
    _units = {'Volumetric flow': 'm3/hr'}

    def _run(self):
        influent = self.ins[0]
        water, brine = self.outs

        self.design_results['Volumetric flow'] = self.F_vol_in

        # Based on stream 626 and 627 in ref [1]
        water.imass['Water'] = 376324/(376324+4967) * influent.imass['Water']
        brine.mol = influent.mol - water.mol
        water.T = brine.T = influent.T


class Skipped(Unit):
    '''
    Copy ins[`main_in`] as ins[`main_out`].
    Can be also used to calculate the cost of wastewater treatment
    by clearing all WWT-related units and streams.
    '''
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    cache_dct = {} # to save some computation effort

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 main_in=0, main_out=0, wwt_units=[], clear_cost=False):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.main_in = main_in
        self.main_out = main_out
        self.wwt_units = wwt_units
        self.clear_cost = clear_cost

    def _run(self):
        self.outs[self.main_out].copy_like(self.ins[self.main_in])

    def _cost(self):
        if not self.clear_cost: return
        for u in self.wwt_units:
            u.baseline_purchase_costs.clear()
            u.installed_costs.clear()
            for hu in u.heat_utilities: hu.empty()
            u.power_utility(0)
            for s in u.ins+u.outs:
                if s.price: s.empty()


class CHP(Unit):
    '''
    Used to estimate the cost of producing electricity.

    Parameters
    ----------
    eff : float
        Combined efficiency for combustion and power generation.
    unit_CAPEX : float
        Capital cost of the CHP per kW of power generated, $/kW.

    References
    ----------
    [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.
    '''
    _ins_size_is_fixed = False
    _F_BM_default = {'CHP': 1.}

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 eff=0.3375, # average of 40.5%, 27%, 36%, and 31.5% for the four types in ref [1]
                 unit_CAPEX=1225):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.eff = eff
        self.unit_CAPEX = unit_CAPEX

    def _run(self):
        mixed, = self.outs
        mixed.mix_from(self.ins)
        H_net = get_combustion_energy(mixed, 1)
        self.H_for_power = H_net * self.eff


    def _design(self):
        kW = self.H_for_power / 3600 # kJ/hr to kW
        self.baseline_purchase_costs['CHP'] = kW*self.unit_CAPEX
        self.power_utility(-kW)


# %%

# =============================================================================
# System function
# =============================================================================

@bst.SystemFactory(
    ID='wastewater_sys',
    outs=[dict(ID='biogas'),
          dict(ID='sludge'),
          dict(ID='recycled_water'),
          dict(ID='brine')],
    fixed_ins_size=False,
)
def create_wastewater_process(ins, outs, process_ID='6', flowsheet=None,
                            skip_IC=False, IC_kwargs={},
                            skip_AnMBR=False, AnMBR_kwargs={},
                            skip_AeF=False, AF_kwargs={}):
    if flowsheet:
        bst.main_flowsheet.set_flowsheet(flowsheet)
    wwt_streams = ins
    biogas, sludge, recycled_water, brine = outs

    ######################## Units ########################
    # Mix waste liquids for treatment
    X = str(process_ID)
    MX01 = bst.units.Mixer(f'M{X}01', ins=wwt_streams)

    RX01_outs = (f'biogas_R{X}01', 'IC_eff', 'IC_sludge')
    if skip_IC:
        RX01 = Skipped(f'R{X}01', ins=MX01-0, outs=RX01_outs, main_in=0, main_out=1) # outs[0] is the gas phase 
    else:
        RX01 = InternalCirculationRx(f'R{X}01', ins=MX01-0, outs=RX01_outs,
                                     T=35+273.15, **IC_kwargs)

    RX02_outs = (f'biogas_R{X}02', f'permeate_R{X}02', f'sludge_R{X}02', f'vent_R{X}02')
    if skip_AnMBR:
        RX02 = Skipped(f'R{X}02', ins=RX01-1, outs=RX02_outs, main_in=0, main_out=1)
    else:
        # Just setting the prices, flows will be updated upon simulation
        naocl_RX02 = tmo.Stream(f'naocl_R{X}02', NaOCl=0.125, Water=1-0.125, units='kg/hr')
        naocl_RX02.price = (naocl_RX02.F_mass/naocl_RX02.F_vol/1000)*prices['naocl'] # $/L to $/kg
        citric_RX02 = tmo.Stream(f'citric_R{X}02', CitricAcid=1, units='kg/hr')
        citric_RX02.price = (citric_RX02.F_mass/citric_RX02.F_vol/1000)*prices['citric_acid'] # $/L to $/kg
        bisulfite_RX02 = tmo.Stream(f'bisulfite_R{X}02', Bisulfite=0.38, Water=1-0.38, units='kg/hr')
        bisulfite_RX02.price = (bisulfite_RX02.F_mass/bisulfite_RX02.F_vol/1000)*prices['bisulfite'] # $/L to $/kg

        RX02 = AnMBR(f'R{X}02', ins=(RX01-1, '', naocl_RX02, citric_RX02,
                                     bisulfite_RX02, f'air_R{X}02'),
                     outs=RX02_outs,
                     reactor_type='CSTR',
                     membrane_configuration='cross-flow',
                     membrane_type='multi-tube',
                     membrane_material='ceramic',
                     include_aerobic_filter=False,
                     add_GAC=False,
                     include_degassing_membrane=True,
                     T=None, # heat loss will be adjusted later
                     # Below include in the TEA
                     include_pump_building_cost=False,
                     include_excavation_cost=False, **AnMBR_kwargs)

    RX03_outs = (f'biogas_R{X}03', f'treated_R{X}03', f'sludge_R{X}03', f'vent_R{X}03')
    if skip_AeF:
        RX03 = Skipped(f'R{X}03', ins=(RX02-1, ''), outs=RX03_outs, main_in=0, main_out=1)
    else:
        RX03 = PolishingFilter(f'R{X}03', ins=(RX02-1, '', f'air_R{X}03'), outs=RX03_outs,
                              filter_type='aerobic',
                              include_degassing_membrane=False,
                              T=None, # heat loss will be adjusted later
                              # Below include in the TEA
                              include_pump_building_cost=False,
                              include_excavation_cost=False,
                              **AF_kwargs)
   
    # # This isn't working, think of a better way to deal with it
    # _RX03_cost = RX03._cost
    # def adjust_heat_loss():
    #     _RX03_cost()
    #     loss_kW = RX02._heat_loss + RX03._heat_loss
    #     # Assume the heat loss in RX02/RX03 can be compensated by heat exchange
    #     # with RX01 with an 80% heat transfer efficiency
    #     RX01.heat_utilities[0].duty += loss_kW * _kW_to_kJhr / 0.8
    #     RX02.power_utility.rate -= RX02._heat_loss
    #     RX03.power_utility.rate -= RX03._heat_loss
    # RX03._cost = adjust_heat_loss

    bst.units.Mixer(f'M{X}02', ins=(RX01-0, RX02-0, RX03-0), outs=biogas)

    # Recycled the majority of sludge (96%) to the aerobic filter,
    # 96% from the membrane bioreactor in ref [2]
    SX01 = bst.units.Splitter(f'S{X}01', ins=RX03-2, outs=(f'recycled_S{X}01', f'wasted_S{X}01'),
                              split=0.96)

    solubles = [i.ID for i in SX01.chemicals if not i.ID in default_insolubles]
    SX02 = BeltThickener(f'S{X}02', ins=(RX01-2, RX02-2, SX01-1),
                         outs=(f'eff_S{X}02', f'sludge_S{X}02'),
                         sludge_moisture=0.96, solubles=solubles)

    SX03 = SludgeCentrifuge(f'S{X}03', ins=SX02-1,
                            outs=(f'centrate_S{X}03', sludge),
                            sludge_moisture=0.8, solubles=solubles,
                            centrifuge_type='reciprocating_pusher')

    # Mix recycles to aerobic digestion
    bst.units.Mixer(f'M{X}03', ins=(SX01-0, SX02-0, SX03-0), outs=1-RX03)

    # Reverse osmosis to treat aerobically polished water
    SX04 = ReverseOsmosis(f'S{X}04', ins=RX03-1, outs=(recycled_water, ''))

    # A process specification for wastewater treatment cost calculation
    Skipped('Caching', ins=SX04-1, outs=brine, main_in=0, main_out=0,
            wwt_units=[u for u in bst.main_flowsheet.unit if (u.ID[1:1+len(X)]==X)])