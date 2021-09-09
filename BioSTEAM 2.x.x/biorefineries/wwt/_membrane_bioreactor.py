#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


'''
TODO:
    - Add algorithms for other configurations
    (AF, submerged, sparging, GAC, flat sheet, hollow fiber)
    - Add OLR and volume calculation
    - Maybe add AeMBR as well (make an MBR superclass)
        - AeMBR can use higher flux and allows for lower transmembrane pressure

References
----------
[1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
Valorization of Dilute Organic Carbon Waste Streams.
Energy Environ. Sci. 2016, 9 (3), 1102–1112.
https://doi.org/10.1039/C5EE03715H.
'''

import math
import biosteam as bst
import thermosteam as tmo
# from collections.abc import Iterable
from biosteam.exceptions import DesignError

# from biorefineries.wwt import (
#     default_insolubles,
#     InternalCirculationRx,
#     WWTpump,
#     PolishingFilter,
#     new_price
#     )
# from biorefineries.utils import (
#     auom,
#     compute_stream_COD,
#     format_str,
#     get_BD_dct,
#     get_split_dct,
#     )
from _chemicals import default_insolubles
from _internal_circulation_rx import InternalCirculationRx
from _wwt_pump import WWTpump
from _polishing_filter import PolishingFilter
from _settings import new_price
from utils import (
    auom,
    compute_stream_COD,
    format_str,
    get_BD_dct,
    get_split_dct,
    )

__all__ = ('AnMBR',)

_ft_to_m = auom('ft').conversion_factor('m')
_ft2_to_m2 = auom('ft2').conversion_factor('m2')
_ft3_to_m3 = auom('ft3').conversion_factor('m3')
_ft3_to_gal = auom('ft3').conversion_factor('gallon')
_m3_to_gal = auom('m3').conversion_factor('gal')
_cmh_to_mgd = _m3_to_gal * 24 / 1e6 # cubic meter per hour to million gallon per day
_lb_to_kg = auom('lb').conversion_factor('kg')



# %%

class AnMBR(bst.Unit):
    '''
    Anaerobic membrane bioreactor (AnMBR) for wastewater treatment as in
    Shoener et al. [1]_ Some assumptions adopted from Humbird et al. [2]_

    In addition to the anaerobic treatment, an optional second stage can be added,
    which can be aerobic filter or granular activated carbon (GAC).

    Parameters
    ----------
    reactor_type : str
        Can either be "CSTR" for continuous stirred tank reactor
        or "AF" for anaerobic filter.
    membrane_configuration : str
        Can either be "cross-flow" or "submerged".
    membrane_type : str
        Can be "hollow fiber" ("submerged" configuration only),
        "flat sheet" (either "cross-flow" or "submerged" configuration),
        or "multi-tube" ("cross-flow" configuration only).
    membrane_material : str
        Can be any of the plastics ("PES", "PVDF", "PET", "PTFE")
        for any of the membrane types ("hollow fiber", "flat sheet", "multi-tube"),
        or "sintered steel" for "flat sheet",
        or "ceramic" for "multi-tube".
    include_aerobic_filter : bool
        Whether to include an aerobic filtration process in this AnMBR,
        can only be True in "AF" (not "CSTR") reactor.
    add_GAC : bool
        If to add granual activated carbon to enhance biomass retention,
        can only be True for the "submerged" configuration.
    include_degassing_membrane : bool
        If to include a degassing membrane to enhance methane
        (generated through the digestion reaction) recovery.
    biodegradability : float or dict
        Biodegradability of chemicals,
        when shown as a float, all biodegradable chemicals are assumped to have
        the same degradability.
    Y : float
        Biomass yield, [kg biomass/kg consumed COD].
    split : dict
        Component-wise split to the treated water.
        E.g., {'Water':1, 'WWTsludge':0} indicates all of the water goes to
        the treated water and all of the WWTsludge goes to the wasted sludge.
        Default splits (based on the membrane bioreactor in [2]_) will be used
        if not provided.
    T : float
        Temperature of the reactor.
        Will not control temperature if provided as None.
    include_pump_building_cost : bool
        Whether to include the construction cost of pump building.
    include_excavation_cost : bool
        Whether to include the construction cost of excavation.
    kwargs : dict
        Other keyword arguments (e.g., J_max, SGD).

    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.

    .. [2] Humbird et al., Process Design and Economics for Biochemical Conversion of
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf
    '''
    _N_ins = 6 # influent, recycle (optional), naocl, citric acid, bisulfite, air (optional)
    _N_outs = 4 # biogas, effluent, waste sludge, air (optional)

    # Equipment-related parameters
    _N_train_min = 2
    _cas_per_tank_spare = 2

    _mod_surface_area = {
        'hollow fiber': 370,
        'flat sheet': 1.45/_ft2_to_m2,
        'multi-tube': 32/_ft2_to_m2
        }

    _mod_per_cas = None
    _mod_per_cas_range = {
        'hollow fiber': (30, 48), # min, max
        'flat sheet': (150, 200),
        'multi-tube': (44, 48)
        }

    _cas_per_tank = None
    _cas_per_tank_range = (16, 22)

    _N_blower = 0

    _W_tank = 21
    _D_tank = 12

    _W_dist = 4.5
    _W_eff = 4.5

    _L_well = 8
    _W_well = 8
    _D_well = 12

    _excav_slope = 1.5
    _constr_access = 3

    # Operation-related parameters
    _HRT = 24
    _J_max = 12
    _TMP_dct = {
        'cross-flow': 2.5,
        'submerged': 2.5,
        }
    _TMP_aerobic = None
    _recir_ratio = 2.25
    _v_cross_flow = 1.2
    _v_GAC = 8
    _SGD = 0.625
    _AFF = 3.33
    _sludge_conc = 10.5

    _refresh_rxns = InternalCirculationRx._refresh_rxns

    # Other equipment
    auxiliary_unit_names = ('heat_exchanger',)
    _pumps =  ('perm', 'retent', 'recir', 'sludge', 'naocl', 'citric', 'bisulfite',
               'AF', 'AeF')


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 reactor_type='CSTR',
                 membrane_configuration='cross-flow',
                 membrane_type='multi-tube',
                 membrane_material='ceramic',
                 include_aerobic_filter=False,
                 add_GAC=False,
                 include_degassing_membrane=True,
                 biodegradability={}, Y=0.05, split=None,
                 T=35+273.15,
                 include_pump_building_cost=False,
                 include_excavation_cost=False,
                 **kwargs):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.reactor_type = reactor_type
        self.include_aerobic_filter = include_aerobic_filter
        self.membrane_configuration = membrane_configuration
        self.membrane_type = membrane_type
        self.membrane_material = membrane_material
        self.add_GAC = add_GAC
        self.include_degassing_membrane = include_degassing_membrane
        self.biodegradability = \
            biodegradability if biodegradability else get_BD_dct(self.chemicals)
        self.Y = Y
        self.split = split if split else get_split_dct(self.chemicals)
        self.T = T
        self.include_pump_building_cost = include_pump_building_cost
        self.include_excavation_cost = include_excavation_cost

        # Initiate the attributes
        self.heat_exchanger = hx = bst.HXutility(None, None, None, T=T)
        self.heat_utilities = hx.heat_utilities
        self._refresh_rxns()

        for k, v in kwargs.items():
            setattr(self, k, v)

        self._check_design()


    def _check_design(self):
        reactor_type = self.reactor_type
        m_config = self.membrane_configuration
        m_type = self.membrane_type
        m_material = self.membrane_material

        if reactor_type == 'CSTR':
            if self.include_aerobic_filter:
                raise DesignError('Aerobic filteration cannot be used in CSTR.')

        if m_config == 'submerged':
            if not m_type in ('hollow fiber', 'flat sheet'):
                raise DesignError('Only "hollow fiber" or "flat sheet" is allowed '
                                  'for "submerged" membrane, not {m_type}.')
        else: # cross-flow
            if not m_type in ('flat sheet', 'multi-tube'):
                raise DesignError('Only "flat sheet" or "multi-tube" is allowed '
                                  'for "cross-flow" membrane, not {m_type}.')
            if self.add_GAC:
                raise DesignError('No GAC should be added '
                                  '(i.e., `add_GAC` can only be False'
                                  'for "cross-flow" membrane.')

        plastics = ('PES', 'PVDF', 'PET', 'PTFE')
        if m_type == 'hollow fiber':
            if not m_material in plastics:
                raise DesignError(f'Only plastic materials {plastics} '
                                  'allowed for "hollow fiber" membrane',
                                  f'not "{m_material}".')
        elif m_type == 'flat sheet':
            if not m_material in (*plastics, 'sintered steel'):
                raise DesignError(f'Only plastic materials {plastics} and "sintered steel"'
                                  'allowed for "flat sheet" membrane',
                                  f'not "{m_material}".')
        else: # multi-tube
            if not m_material in (*plastics, 'ceramic'):
                raise DesignError(f'Only plastic materials {plastics} and "ceramic"'
                                  'allowed for "multi-tube" membrane',
                                  f'not "{m_material}".')


    @staticmethod
    def _degassing(receiving_stream, original_stream):
        InternalCirculationRx._degassing(receiving_stream, original_stream)


    @staticmethod
    def compute_COD(stream):
        return compute_stream_COD(stream)


    # =========================================================================
    # _run
    # =========================================================================
    def _run(self):
        raw, recycled, naocl, citric, bisulfite, air_in = self.ins
        biogas, perm, sludge, air_out = self.outs
        degassing = self._degassing

        # Initiate the streams
        biogas.phase = 'g'
        biogas.empty()

        inf = raw.copy()
        inf.mix_from((raw, recycled))
        self._inf = inf.copy() # this stream will be preserved (i.e., no reaction)

        # Chemicals for cleaning
        # Assume both chemicals are used up
        # 2.2 L/yr/cmd of 12.5 wt% solution (15% vol)
        naocl_solution = (2.2/1e3/365/24) * (inf.F_vol*24) # m3/hr solution
        naocl.empty()
        naocl.imass['NaOCl', 'Water'] = [0.125, 1-0.125]
        naocl.F_vol = naocl_solution
        # Convert the price from $/L in the `new_price` dict to $/kg
        naocl.price = (naocl.F_mass/naocl.F_vol) * new_price['NaOCl']

        # 0.6 L/yr/cmd of 100 wt% solution, 13.8 lb/kg
        citric_solution = (0.6/1e3/365/24) * (inf.F_vol*24) # m3/hr pure
        citric.empty()
        citric.ivol['CitricAcid'] = citric_solution
        citric.price = (citric.F_mass/citric.F_vol) * new_price['CitricAcid']

        # 0.35 L/yr/cmd of 38% solution, 3.5 lb/gal
        bisulfite_solution = (0.35/1e3/365/24) * (inf.F_vol*24) # m3/hr solution
        bisulfite.empty()
        bisulfite.imass['Bisulfite', 'Water'] = [0.38, 1-0.38]
        bisulfite.F_vol = bisulfite_solution
        bisulfite.price = (bisulfite.F_mass/bisulfite.F_vol) * new_price['Bisulfite']

        # For pump design
        self._compute_mod_case_tank_N()
        Q_R_mgd, Q_IR_mgd = self._compute_liq_flows()
        retent, recir = inf.copy(), inf.copy()
        retent.F_mass *= Q_R_mgd / self.Q_mgd
        recir.F_mass *= Q_IR_mgd / self.Q_mgd
        self._retent, self._recir = retent, recir

        # Effluents
        self.growth_rxns(inf.mol)
        self.biogas_rxns(inf.mol)
        tmo.separations.split(inf, perm, sludge, self._isplit.data)

        sludge_conc = self._sludge_conc
        insolubles = tuple(i.ID for i in self.chemicals if i.ID in default_insolubles)
        m_insolubles = sludge.imass[insolubles].sum()
        if m_insolubles/sludge.F_vol <= sludge_conc:
            diff = sludge.ivol['Water'] - m_insolubles/sludge_conc
            sludge.ivol['Water'] = m_insolubles/sludge_conc
            perm.ivol['Water'] += diff

        degassing(biogas, perm)
        degassing(biogas, sludge)

        # Gas for sparging, no sparging needed if submerged or using GAC
        air_out.link_with(air_in)
        air_in.T = 17 + 273.15
        self._design_blower()

        if self.T is not None:
            perm.T = sludge.T = biogas.T = air_out.T = self.T


    # Called by _run
    def _compute_mod_case_tank_N(self):
        N_mod_min, N_mod_max = self.mod_per_cas_range[self.membrane_type]
        N_cas_min, N_cas_max = self.cas_per_tank_range

        mod_per_cas, cas_per_tank = N_mod_min, N_cas_min

        J, J_max, N_train = self.J, self.J_max, self._N_train_min
        while J > J_max:
            mod_per_cas += 1
            if mod_per_cas == N_mod_max + 1:
                if cas_per_tank == N_cas_max + 1:
                    N_train += 1
                    mod_per_cas, cas_per_tank = N_mod_min, N_cas_min
                else:
                    cas_per_tank += 1
                    mod_per_cas = N_mod_min

        self._N_train, self._mod_per_cas, self._cas_per_tank = \
            N_train, mod_per_cas, cas_per_tank


    # Called by _run
    def _compute_liq_flows(self):
        m_type = self.membrane_type

        if m_type == 'multi-tube':
            # Cross-flow flow rate per module,
            # based on manufacture specifications for compact 33, [m3/hr]
            Q_cross_flow = 53.5 * self.v_cross_flow
            Q_R_cmh = self.N_mod_tot * Q_cross_flow # total retentate flow rate, [m3/hr]
            Q_R_mgd = Q_R_cmh * _cmh_to_mgd # [mgd]

        Q_mgd, recir_ratio = self.Q_mgd, self.recir_ratio
        if Q_mgd*recir_ratio >= Q_R_mgd:
            Q_IR_mgd = Q_mgd*recir_ratio - Q_R_mgd
        else:
            Q_IR_mgd = 0
            # # Gives really large recir_ration,
            # # probably shouldn't back-calculate this way
            # self._recir_ratio = Q_R_mgd / Q_mgd

        if self.add_GAC:
            Q_upflow_req = (self.v_GAC/_ft_to_m) * \
                self.L_membrane_tank*self.W_tank*self.N_train * 24 * _ft3_to_gal / 1e6
            Q_IR_add_mgd = max(0, (Q_mgd+Q_IR_mgd)-Q_upflow_req)
            Q_IR_mgd += Q_IR_add_mgd
            self._recir_ratio = Q_IR_mgd / Q_mgd

        return Q_R_mgd, Q_IR_mgd


    # Called by _run
    def _design_blower(self):
        air = self.ins[-1]

        if (not self.add_GAC) and (self.membrane_configuration=='submerged'):
            gas = self.SGD * self.mod_surface_area*_ft2_to_m2 # [m3/h]
            gas /= (_ft3_to_m3 * 60) # [ft3/min]
            gas_train = gas * self.N_train*self.cas_per_tank*self.mod_per_cas

            TCFM = math.ceil(gas_train) # total cubic ft per min
            N = 1
            if TCFM <= 30000:
                CFMB = TCFM / N # cubic ft per min per blower
                while CFMB > 7500:
                    N += 1
                    CFMB = TCFM / N
            elif 30000 < TCFM <= 72000:
                CFMB = TCFM / N
                while CFMB > 18000:
                    N += 1
                    CFMB = TCFM / N
            else:
                CFMB = TCFM / N
                while CFMB > 100000:
                    N += 1
                    CFMB = TCFM / N

            gas_m3_hr = TCFM * _ft3_to_m3 * 60 # ft3/min to m3/hr

            air.ivol['N2'] = 0.79
            air.ivol['O2'] = 0.21
            air.F_vol = gas_m3_hr
        else: # no sparging/blower needed
            TCFM = CFMB = 0.
            N = -1 # to account for the spare
            air.empty()

        D = self.design_results
        D['Total air flow [CFM]'] = TCFM
        D['Blower capacity [CFM]'] = CFMB
        D['Blowers'] = self._N_blower = N + 1 # add a spare


    # =========================================================================
    # _design
    # =========================================================================
    def _design(self):
        D = self.design_results
        D['Treatment train'] = self.N_train
        D['Cassette per train'] = self.cas_per_tank
        D['Module per cassette'] = self.mod_per_cas
        D['Total membrane modules'] = self.N_mod_tot

        # Step A: Reactor and membrane tanks
        # Call the corresponding design function
        # (_design_CSTR or _design_AF)
        func = getattr(self, f'_design_{self.reactor_type}')

        wall, slab, excavation = func()
        D['Wall concrete [ft3]'] = wall
        D['Slab concrete [ft3]'] = slab
        D['Excavation [ft3]'] = excavation

        # Optional addition of packing media (used in filters)
        ldpe, hdpe = 0., 0.
        for i in (self.AF, self.AeF):
            if i is None:
                continue
            ldpe += i.design_results['Packing LDPE [m3]']
            hdpe += i.design_results['Packing LDPE [m3]']

        # Optional addition of GAC
        D['GAC [kg]'] = self._design_GAC()

        # Step B: Membrane
        # Call the corresponding design function
        # (_design_hollow_fiber, _design_flat_sheet, or _design_multi_tube)
        m_type = format_str(self.membrane_type)
        func = getattr(self, f'_design_{m_type}')
        D['Membrane [m3]'] = func()

        # Step C: Pumps
        pipe, pumps, hdpe = self._design_pump()
        D['Pipe stainless steel [kg]'] = pipe
        D['Pump stainless steel [kg]'] = pumps
        D['Pump chemical storage HDPE [m3]'] = hdpe

        # Step D: Degassing membrane
        D['Degassing membrane'] = self.N_degasser


    ### Step A functions ###
    # Called by _design
    def _design_CSTR(self):
        N_train = self.N_train
        W_dist, L_CSTR , W_eff, L_membrane_tank = \
            self.W_dist, self.L_CSTR, self.W_eff, self.L_membrane_tank
        W_PB, W_BB, D_tank = self.W_PB, self.W_BB, self.D_tank
        t_wall, t_slab = self.t_wall, self.t_slab
        SL, CA = self.excav_slope, self.constr_access

        ### Concrete calculation ###
        W_N_trains = (self.W_tank+2*t_wall)*N_train - t_wall*(N_train-1)

        D = D_tank + 2 # add 2 ft of freeboard
        t = t_wall + t_slab

        get_VWC = lambda L1, N: N * t_wall * L1 * D
        get_VSC = lambda L2: t * L2 * W_N_trains

        # Concrete for distribution channel, [ft3]
        VWC_dist = get_VWC(L1=(W_N_trains+W_dist), N=2)
        VSC_dist = get_VSC(L2=(W_dist+2*t_wall))

        # Concrete for CSTR tanks, [ft3]
        VWC_CSTR = get_VWC(L1=L_CSTR, N=(N_train+1))
        VSC_CSTR = get_VSC(L2=L_CSTR)

        # Concrete for effluent channel, [ft3]
        VWC_eff = get_VWC(L1=(W_N_trains+W_eff), N=2)
        VSC_eff = get_VSC(L2=(W_eff+2*t_wall))

        # Concrete for the pump/blower building, [ft3]
        VWC_PBB = get_VWC(L1=(W_N_trains+W_PB+W_BB), N=2)
        VSC_PBB = get_VSC(L2=(W_PB+t_wall+W_BB))

        if self.membrane_configuration == 'submerged':
            VWC_membrane_tank, VSC_membrane_tank, VWC_well, VSC_well = \
                self._design_membrane_tank(self, D, N_train, W_N_trains,
                                           self.L_membrane_tank)
        else:
            VWC_membrane_tank, VSC_membrane_tank, VWC_well, VSC_well = 0., 0., 0., 0.

        # Total volume of wall concrete, [ft3]
        VWC = VWC_dist + VWC_CSTR + VWC_eff + VWC_PBB + VWC_membrane_tank + VWC_well

        # Total volume of slab concrete [ft3]
        VSC = VSC_dist + VSC_CSTR + VSC_eff + VSC_PBB + VSC_membrane_tank + VSC_well

        ### Excavation calculation ###
        get_VEX = lambda L_bttom, W_bottom, diff: \
            0.5 * D_tank * (L_bottom*W_bottom+(L_bottom+diff)*(W_bottom+diff)) # bottom+top

        L_bottom = W_dist + L_CSTR + W_eff + L_membrane_tank + 2*CA
        W_bottom = W_N_trains + 2*CA
        diff = D_tank * SL

        # Excavation volume for membrane tanks, [ft3]
        VEX_membrane_tank = get_VEX(L_bottom, W_bottom, diff)

        # Excavation volume for pump/blower building, [ft3]
        VEX_PBB = get_VEX((W_PB+W_BB+2*CA), W_bottom, diff)

        VEX = VEX_membrane_tank + VEX_PBB

        return VWC, VSC, VEX


    # Called by _design
    def _design_AF(self):
        '''NOT READY YET.'''
        # Use FilterTank


    # Called by _design_CSTR/_design_AF
    def _design_membrane_tank(self, D, N_train, W_N_trains, L_membrane_tank,
                              t_wall, t_slab):
        L_well, W_well, D_well = self.L_well, self.W_well, self.D_well

        # Concrete for membrane tanks, [ft3]
        t = t_wall + t_slab
        VWC_membrane_tank = (N_train+1) * t_wall * L_membrane_tank * D
        VSC_membrane_tank = t * L_membrane_tank * W_N_trains

        # Concrete for wet well (mixed liquor storage), [ft3]
        L = L_well + 2*t_wall
        W = W_well + 2*t_wall

        VWC_well = 2 * t_wall * (L_well+W) * D_well
        VSC_well = (t_slab+t_wall) * L * W

        return VWC_membrane_tank, VSC_membrane_tank, VWC_well, VSC_well


    # Called by _design
    def _design_GAC(self):
        '''NOT READY YET.'''
        if not self.add_GAC:
            return 0

        M_GAC = 1
        return M_GAC


    ### Step B functions ###
    # Called by _design
    def _design_hollow_fiber(self):
        '''NOT READY YET.'''


    # Called by _design
    def _design_flat_sheet(self):
        '''NOT READY YET.'''


    # Called by _design
    def _design_multi_tube(self):
        # # 0.01478 is volume of material for each membrane tube, [m3]
        # #      L_tube               OD        ID
        # V_tube = 3 * math.pi/4 * ((6e-3)**2-(5.2e-3)**2)
        # V_SU = 700 * V_tube # V for each small unit [m3]
        # M_SU = 1.78*1e3 * V_SU # mass = density*volume, [kg/m3], not used
        return self.N_mod_tot*0.01478


    ### Step C function ###
    # Called by _design
    def _design_pump(self):
        rx_type, m_config, pumps = \
            self.reactor_type, self.membrane_configuration, self._pumps

        # IDs = ['perm', 'retent', 'recir', 'sludge', 'naocl', 'citric', 'bisulfite']

        ins_dct = {
            'perm': self.outs[1].proxy(),
            'retent': self._retent,
            'recir': self._recir,
            'sludge': self.outs[2].proxy(),
            'naocl': self.ins[2].proxy(),
            'citric': self.ins[3].proxy(),
            'bisulfite': self.ins[4].proxy(),
            }

        type_dct = {
            'perm': f'permeate_{m_config}',
            'retent': f'retentate_{rx_type}',
            'recir': f'recirculation_{rx_type}',
            'sludge': 'sludge',
            'naocl': 'chemical',
            'citric': 'chemical',
            'bisulfite': 'chemical',
            }

        inputs_dct = {
            'perm': (self.cas_per_tank, self.D_tank, self.TMP_anaerobic,
                     self.include_aerobic_filter),
            'retent': (self.cas_per_tank,),
            'recir': (self.L_CSTR,),
            'sludge': (),
            'naocl': (),
            'citric': (),
            'bisulfite': (),
            }

        WWTpump._batch_adding_pump(self, pumps[:-2], ins_dct, type_dct, inputs_dct)

        self.AF_pump = self.AF.lift_pump if self.AF else None
        self.AeF_pump = self.AeF.lift_pump if self.AeF else None

        pipe_ss, pump_ss, hdpe = 0., 0., 0.
        for i in pumps:
            p = getattr(self, f'{i}_pump')
            if p == None:
                continue

            p.simulate()
            pipe_ss += p.design_results['Pipe stainless steel [kg]']
            pump_ss += p.design_results['Pump stainless steel [kg]']
            hdpe += p.design_results['Chemical storage HDPE [m3]']

        return pipe_ss, pump_ss, hdpe


    # =========================================================================
    # _cost
    # =========================================================================
    def _cost(self):
        D, C, F_BM, lifetime = self.design_results, self.baseline_purchase_costs, \
            self.F_BM, self._default_equipment_lifetime

        ### Capital ###
        # Concrete and excavaction
        VEX, VWC, VSC = \
            D['Excavation [ft3]'], D['Wall concrete [ft3]'], D['Slab concrete [ft3]']
        # 27 is to convert the VEX from ft3 to yard3
        C['Reactor excavation'] = VEX/27*8 if self.include_excavation_cost else 0.
        C['Wall concrete'] = VWC / 27 * 650
        C['Slab concrete'] = VSC / 27 * 350

        # Membrane
        # TODO: now assume $8/ft2 for all membrane materials,
        # maybe check the price for different materials
        C['Membrane'] = 8 * D['Membrane [m3]'] / _ft2_to_m2
        F_BM['Membrane'] = 1 + 0.15 # assume 15% for replacement labor
        lifetime['Membrane'] = 10

        # GAC
        # $13.78/kg
        C['GAC'] = 13.78 * D['GAC [kg]']

        # Packing material
        ldpe, hdpe = 0., 0.
        for i in (self.AF, self.AeF):
            if i is None:
                continue
            ldpe += i.baseline_purchase_costs['Packing LDPE [m3]']
            hdpe += i.baseline_purchase_costs['Packing HDPE [m3]']

        # Pump
        # Note that maintenance and operating costs are included as a lumped
        # number in the biorefinery thus not included here
        # TODO: considering adding the O&M and letting user choose if to include
        pumps, building = self._cost_pump()
        C['Pumps'] = pumps
        C['Pump building'] = building if self.include_pump_building_cost else 0.
        C['Pump excavation'] = VEX/27*0.3 if self.include_excavation_cost else 0.

        F_BM['Pumps'] = F_BM['Pump building'] = F_BM['Pump excavation'] = \
            1.18 * (1+0.007) # 0.007 is for  miscellaneous costs
        lifetime['Pumps'] = 15

        # Blower and air pipe
        TCFM, CFMB = D['Total air flow [CFM]'], D['Blower capacity [CFM]']
        C['Air pipes'], C['Blowers'], C['Blower building'] = self._cost_blower(TCFM, CFMB)
        F_BM['Blowers'] = 2 * 1.11
        F_BM['Blower building'] = 1.11
        lifetime['Blowers'] = 15

        # Degassing membrame
        C['Degassing membrane'] = 10000 * D['Degassing membrane']

        # Set bare module factor to 1 if not otherwise provided
        for k in C.keys():
            F_BM[k] = 1 if not F_BM.get(k) else F_BM.get(k)

        ### Heat and power ###
        # Heat loss, assume air is 17°C, ground is 10°C
        T = self.T
        if T is None:
            loss = 0.
        else:
            N_train, L_CSTR, W_tank, D_tank = \
                self.N_train, self.L_CSTR, self.W_tank, self.D_tank
            A_W = 2 * (L_CSTR+W_tank) * D_tank
            A_F = L_CSTR * W_tank
            A_W *= N_train * _ft2_to_m2
            A_F *= N_train * _ft2_to_m2

            loss = 0.7 * (T-(17+273.15)) * A_W / 1e3 # 0.7 W/m2/°C for wall
            loss += 1.7 * (T-(10+273.15)) * A_F / 1e3 # 1.7 W/m2/°C for floor
            loss += 0.95 * (T-(17+273.15)) * A_F / 1e3 # 0.95 W/m2/°C for floating cover
        self._heat_loss = loss

        # Fluid heating
        inf = self._inf
        if T:
            H_at_T = inf.thermo.mixture.H(mol=inf.mol, phase='l', T=T, P=101325)
            duty = -(inf.H - H_at_T)
        else:
            duty = 0
        self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, inf)

        # Pumping
        pumping = 0.
        for ID in self._pumps:
            p = getattr(self, f'{ID}_pump')
            if p is None:
                continue
            pumping += p.power_utility.rate

        # Gas
        sparging = 0. #!!! output from submerge design
        degassing = 3 * self.N_degasser # assume each uses 3 kW

        self.power_utility.rate = sparging + degassing + pumping + loss


    # Called by _cost
    _cost_pump = PolishingFilter._cost_pump


    # Called by _cost
    def _cost_blower(self, TCFM, CFMB):
        AFF = self.AFF

        # Air pipes
        # Note that the original codes use CFMD instead of TCFM for air pipes,
        # but based on the coding they are equivalent
        if TCFM <= 1000:
            air_pipes = 617.2 * AFF * (TCFM**0.2553)
        elif 1000 < TCFM <= 10000:
            air_pipes = 1.43 * AFF * (TCFM**1.1337)
        else:
            air_pipes = 28.59 * AFF * (TCFM**0.8085)

        # Blowers
        if TCFM <= 30000:
            ratio = 0.7 * (CFMB**0.6169)
            blowers = 58000*ratio / 100
        elif 30000 < TCFM <= 72000:
            ratio = 0.377 * (CFMB**0.5928)
            blowers = 218000*ratio / 100
        else:
            ratio = 0.964 * (CFMB**0.4286)
            blowers  = 480000*ratio / 100

        # Blower building
        area = 128 * (TCFM**0.256) # building area, [ft2]
        building = area * 90 # 90 is the unit price, [$/ft]

        return air_pipes, blowers, building


    ### Reactor configuration ###
    @property
    def reactor_type(self):
        '''
        [str] Can either be "CSTR" for continuous stirred tank reactor
        or "AF" for anaerobic filter.
        '''
        return self._reactor_type
    @reactor_type.setter
    def reactor_type(self, i):
        if not i.upper() in ('CSTR', 'AF'):
            raise ValueError('`reactor_type` can only be "CSTR", or "AF", '
                             f'not "{i}".')
        self._reactor_type = i.upper()

    @property
    def membrane_configuration(self):
        '''[str] Can either be "cross-flow" or "submerged".'''
        return self._membrane_configuration
    @membrane_configuration.setter
    def membrane_configuration(self, i):
        i = 'cross-flow' if i.lower() in ('cross flow', 'crossflow') else i
        if not i.lower() in ('cross-flow', 'submerged'):
            raise ValueError('`membrane_configuration` can only be "cross-flow", '
                             f'or "submerged", not "{i}".')
        self._membrane_configuration = i.lower()

    @property
    def membrane_type(self):
        '''
        [str] Can be "hollow fiber" ("submerged" configuration only),
        "flat sheet" (either "cross-flow" or "submerged" configuration),
        or "multi-tube" ("cross-flow" configuration only).
        '''
        return self._membrane_type
    @membrane_type.setter
    def membrane_type(self, i):
        i = 'multi-tube' if i.lower() in ('multi tube', 'multitube') else i
        if not i.lower() in ('hollow fiber', 'flat sheet', 'multi-tube'):
            raise ValueError('`membrane_type` can only be "hollow fiber", '
                             f'"flat sheet", or "multi-tube", not "{i}".')
        self._membrane_type = i.lower()

    @property
    def membrane_material(self):
        '''
        [str] Can be any of the plastics ("PES", "PVDF", "PET", "PTFE")
        for any of the membrane types ("hollow fiber", "flat sheet", "multi-tube"),
        or "sintered steel" for "flat sheet",
        or "ceramic" for "multi-tube".
        '''
        return self._membrane_material
    @membrane_material.setter
    def membrane_material(self, i):
        plastics = ('PES', 'PVDF', 'PET', 'PTFE')
        if i.upper() in plastics:
            self._membrane_material = i.upper()
        elif i.lower() in ('sintered steel', 'ceramic'):
            self._membrane_material = i.lower()
        else:
            raise ValueError(f'`membrane_material` can only be plastics materials '
                             f'{plastics}, "sintered steel", or "ceramic", not {i}.')


    ### Reactor/membrane tank ###
    @property
    def AF(self):
        '''[:class:`~.FilterTank`] Anaerobic filter tank.'''
        if self.reactor_type == 'CSTR':
            return None
        return self._AF

    @property
    def AeF(self):
        '''[:class:`~.FilterTank`] Aerobic filter tank.'''
        if not self.include_aerobic_filter:
            return None
        return self._AeF

    @property
    def N_train(self):
        '''[int] Number of treatment train.'''
        if not hasattr(self, '_N_train'):
            return self._N_train_min
        return self._N_train
    @N_train.setter
    def N_train(self, i):
        self._N_train = math.ceil(i)

    @property
    def cas_per_tank_spare(self):
        '''[int] Number of spare cassettes per train.'''
        return self._cas_per_tank_spare
    @cas_per_tank_spare.setter
    def cas_per_tank_spare(self, i):
        self._cas_per_tank_spare = math.ceil(i)

    @property
    def mod_per_cas_range(self):
        '''
        [tuple] Range (min, max) of the number of membrane modules per cassette
        for the current membrane type.
        '''
        return self._mod_per_cas_range
    @mod_per_cas_range.setter
    def mod_per_cas_range(self, i):
        self._mod_per_cas_range[self.membrane_type] = \
            tuple(math.floor(i[0]), math.floor(i[1]))

    @property
    def mod_per_cas(self):
        '''
        [float] Number of membrane modules per cassette for the current membrane type.
        '''
        return self._mod_per_cas or self._mod_per_cas_range[self.membrane_type][0]

    @property
    def cas_per_tank_range(self):
        '''
        [tuple] Range (min, max) of the number of membrane cassette per tank
        (same for all membrane types).
        '''
        return self._cas_per_tank_range
    @cas_per_tank_range.setter
    def cas_per_tank_range(self, i):
        self._cas_per_tank_range = tuple(math.floor(i[0]), math.floor(i[1]))

    @property
    def cas_per_tank(self):
        '''
        [float] Number of membrane cassettes per tank for the current membrane type.
        '''
        return self._cas_per_tank or self._cas_per_tank_range[0]

    @property
    def N_mod_tot(self):
        '''[int] Total number of memberane modules.'''
        return self.N_train * self.cas_per_tank * self.mod_per_cas

    @property
    def mod_surface_area(self):
        '''
        [float] Surface area of the membrane for the current membrane type, [m2/module].
        Note that one module is one sheet for plat sheet and one tube for multi-tube.
        '''
        return self._mod_surface_area[self.membrane_type]
    @mod_surface_area.setter
    def mod_surface_area(self, i):
        self._mod_surface_area[self.membrane_type] = float(i)

    @property
    def L_CSTR(self):
        '''[float] Length of the CSTR tank, [ft].'''
        if self.reactor_type == 'AF':
            return 0
        return self._inf.F_vol/_ft3_to_m3*self.HRT/(self.N_train*self.W_tank*self.D_tank)

    @property
    def L_membrane_tank(self):
        '''[float] Length of the membrane tank, [ft].'''
        #!!! Maybe should set this to 0 for cross-flow?
        return math.ceil((self.cas_per_tank+self.cas_per_tank_spare)*3.4)

    @property
    def W_tank(self):
        '''[float] Width of the reactor/membrane tank (same value), [ft].'''
        return self._W_tank
    @W_tank.setter
    def W_tank(self, i):
        self._W_tank = float(i)

    @property
    def D_tank(self):
        '''[float] Depth of the reactor/membrane tank (same value), [ft].'''
        return self._D_tank
    @D_tank.setter
    def D_tank(self, i):
        self._D_tank = float(i)

    @property
    def W_dist(self):
        '''[float] Width of the distribution channel, [ft].'''
        return self._W_dist
    @W_dist.setter
    def W_dist(self, i):
        self._W_dist = float(i)

    @property
    def W_eff(self):
        '''[float] Width of the effluent channel, [ft].'''
        return self._W_eff
    @W_eff.setter
    def W_eff(self, i):
        self._W_eff = float(i)

    @property
    def t_wall(self):
        '''
        [float] Concrete wall thickness, [ft].
        Minimum of 1 ft with 1 in added for every ft of depth over 12 ft.
        '''
        return 1 + max(self.D_tank-12, 0)/12

    @property
    def t_slab(self):
        '''
        [float] Concrete slab thickness, [ft].
        2 in thicker than the wall thickness.
        '''
        return self.t_wall+2/12


    ### Pump/blower ###
    # @property
    # def pump_dct(self):
    #     '''
    #     [dict] All pumps included in this unit, will be automatically updated
    #     during simulation.
    #     Keys are "permeate", "retentate", "recirculation", "lift",
    #     and "chemical" (chemical storage included),
    #     values are :class:`WWTpump` objects (or None if not applicable).
    #     Note that the value for "chemical" is a tuple of two :class:`WWTpump`,
    #     the first one is for NaOCl and the second one is for citric acid
    #     (both used for membrane cleaning).
    #     '''
    #     return self._pump_dct

    @property
    def N_blower(self):
        '''
        [int] Number of blowers needed for gas sparging
        (not needed for some designs).
        Note that this is not used in costing
        (the cost is estimated based on the total sparging gas need).
        '''
        if not self.add_GAC and self.membrane_configuration=='submerged':
            return self._N_blower
        return 0

    @property
    def N_degasser(self):
        '''
        [int] Number of degassing membrane needed for dissolved biogas removal
        (optional).
        '''
        if self.include_degassing_membrane:
            return math.ceil(self.Q_cmd/24/30) # assume each can hand 30 m3/d of influent
        return 0

    @property
    def W_PB(self):
        '''[float] Width of the pump building, [ft].'''
        if self.membrane_configuration == 'submerged':
            N = self.cas_per_tank
        else: # cross-flow
            N = math.ceil(self.L_CSTR/((1+8/12)+(3+4/12)))

        if 0 <= N <= 10:
            W_PB = 27 + 4/12
        elif 11 <= N <= 16:
            W_PB = 29 + 6/12
        elif 17 <= N <= 22:
            W_PB = 31 + 8/12
        elif 23 <= N <= 28:
            W_PB = 35
        elif N >= 29:
            W_PB = 38 + 4/12
        else:
            W_PB = 0

        return W_PB

    @property
    def L_BB(self):
        '''[float] Length of the blower building, [ft].'''
        if self.membrane_configuration == 'submerged':
            return (69+6/12) if self.cas_per_tank<=18 else (76 + 8/12)
        return 0

    @property
    def W_BB(self):
        '''[float] Width of the blower building, [ft].'''
        if self.membrane_configuration == 'submerged':
            return (18+8/12) if self.cas_per_tank<=18 else 22
        return 0


    ### Wet well (submerged only) ###
    @property
    def L_well(self):
        '''
        [float] Length of the wet well, [ft].
        Only needed for submerged configuration.
        '''
        return self._L_well if self.membrane_configuration == 'submerged' else 0
    @L_well.setter
    def L_well(self, i):
        self._L_well = float(i)

    @property
    def W_well(self):
        '''
        [float] Width of the wet well, [ft].
        Only needed for submerged configuration.
        '''
        return self._W_well if self.membrane_configuration == 'submerged' else 0
    @W_well.setter
    def W_well(self, i):
        self._W_well = float(i)

    @property
    def D_well(self):
        '''
        [float] Depth of the wet well, [ft].
        Only needed for submerged configuration.
        '''
        return self._D_well if self.membrane_configuration == 'submerged' else 0
    @D_well.setter
    def D_well(self, i):
        self._D_well = float(i)


    ### Excavation ###
    @property
    def excav_slope(self):
        '''[float] Slope for excavation (horizontal/vertical).'''
        return self._excav_slope
    @excav_slope.setter
    def excav_slope(self, i):
        self._excav_slope = float(i)

    @property
    def constr_access(self):
        '''[float] Extra room for construction access, [ft].'''
        return self._constr_access
    @constr_access.setter
    def constr_access(self, i):
        self._constr_access = float(i)


    ### Operation-related parameters ###
    @property
    def Q_mgd(self):
        '''
        [float] Influent volumetric flow rate in million gallon per day, [mgd].
        '''
        return self.ins[0].F_vol*_m3_to_gal*24/1e6

    @property
    def Q_gpm(self):
        '''[float] Influent volumetric flow rate in gallon per minute, [gpm].'''
        return self.Q_mgd*1e6/24/60

    @property
    def Q_cmd(self):
        '''
        [float] Influent volumetric flow rate in cubic meter per day, [cmd].
        '''
        return self.Q_mgd *1e6/_m3_to_gal # [m3/day]

    @property
    def Q_cfs(self):
        '''[float] Influent volumetric flow rate in cubic feet per second, [cfs].'''
        return self.Q_mgd*1e6/24/60/60/_ft3_to_gal

    @property
    def HRT(self):
        '''
        [float] Hydraulic retention time, [hr].
        '''
        return self._HRT
    @HRT.setter
    def HRT(self, i):
        self._HRT = float(i)

    @property
    def recir_ratio(self):
        '''
        [float] Internal recirculation ratio, will be updated in simulation
        if the originally set ratio is not adequate for the desired flow
        required by GAC (if applicable).
        '''
        return self._recir_ratio
    @recir_ratio.setter
    def recir_ratio(self, i):
        self._recir_ratio = float(i)

    @property
    def J_max(self):
        '''[float] Maximum membrane flux, [L/m2/hr].'''
        return self._J_max
    @J_max.setter
    def J_max(self, i):
        self._J_max = float(i)

    @property
    def J(self):
        '''[float] Membrane flux, [L/m2/hr].'''
        # Based on the flux of one train being offline
        SA = (self.N_train-1) * self.cas_per_tank * self.mod_per_cas * self.mod_surface_area
        return self._inf.F_vol*1e3/SA # 1e3 is conversion from m3 to L

    @property
    def TMP_anaerobic(self):
        '''[float] Transmembrane pressure in the anaerobic reactor, [psi].'''
        return self._TMP_dct[self.membrane_configuration]
    @TMP_anaerobic.setter
    def TMP_anaerobic(self, i):
        self._TMP_dct[self.membrane_configuration] = float(i)

    @property
    def TMP_aerobic(self):
        '''
        [float] Transmembrane pressure in the aerobic filter, [psi].
        Defaulted to half of the reactor TMP.
        '''
        if not self._include_aerobic_filter:
            return 0.
        else:
            return self._TMP_aerobic or self._TMP_dct[self.membrane_configuration]/2
    @TMP_aerobic.setter
    def TMP_aerobic(self, i):
        self._TMP_aerobic = float(i)

    @property
    def SGD(self):
        '''[float] Specific gas demand, [m3 gas/m2 membrane area/h].'''
        return self._SGD
    @SGD.setter
    def SGD(self, i):
        self._SGD = float(i)

    @property
    def AFF(self):
        '''
        [float] Air flow fraction, used in air pipe costing.
        The default value is calculated as STE/6
        (STE stands for standard oxygen transfer efficiency, and default STE is 20).
        If using different STE value, AFF should be 1 if STE/6<1
        and 3.33 if STE/6>1.
        '''
        return self._AFF
    @AFF.setter
    def AFF(self, i):
        self._AFF = float(i)

    @property
    def v_cross_flow(self):
        '''
        [float] Cross-flow velocity, [m/s].
        '''
        return self._v_cross_flow if self.membrane_configuration=='cross-flow' else 0
    @v_cross_flow.setter
    def v_cross_flow(self, i):
        self._v_cross_flow = float(i)

    @property
    def v_GAC(self):
        '''
        [float] Upflow velocity for GAC bed expansion, [m/hr].
        '''
        return self._v_GAC if self.add_GAC==True else 0
    @v_GAC.setter
    def v_GAC(self, i):
        self._v_GAC = float(i)

    @property
    def biodegradability(self):
        '''
        [float of dict] Biodegradability of chemicals,
        when shown as a float, all biodegradable chemicals are assumped to have
        the same degradability.
        '''
        return self._biodegradability
    @biodegradability.setter
    def biodegradability(self, i):
        if isinstance(i, float):
            if not 0<=i<=1:
                raise ValueError('`biodegradability` should be within [0, 1], '
                                 f'the input value {i} is outside the range.')
            self._biodegradability = i
            return

        for k, v in i.items():
            if not 0<=v<=1:
                raise ValueError('`biodegradability` should be within [0, 1], '
                                 f'the input value for chemical "{k}" is '
                                 'outside the range.')
        self._biodegradability = i

    @property
    def sludge_conc(self):
        '''Concentration of biomass ("WWTsludge") in the waste sludge, [g/L].'''
        return self._sludge_conc
    @sludge_conc.setter
    def sludge_conc(self, i):
        self._sludge_conc = i

    @property
    def Y(self):
        '''[float] Biomass yield, [kg biomass/kg consumed COD].'''
        return self._Y
    @Y.setter
    def Y(self, i):
        if not 0 <= i <= 1:
            raise ValueError('`Y` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Y = i

    @property
    def i_rm(self):
        '''[:class:`np.array`] Removal of each chemical in this reactor.'''
        return self._i_rm

    @property
    def split(self):
        '''Component-wise split to the treated water.'''
        return self._split
    @split.setter
    def split(self, i):
        self._split = i
        self._isplit = self.chemicals.isplit(i, order=None)

    @property
    def biogas_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Biogas production reactions.
        '''
        return self._biogas_rxns

    @property
    def growth_rxns(self):
        '''
        [:class:`tmo.ParallelReaction`] Biomass (WWTsludge) growth reactions.
        '''
        return self._growth_rxns

    @property
    def organic_rm(self):
        '''[float] Overall organic (COD) removal rate.'''
        Qi, Qe = self._inf.F_vol, self.outs[1].F_vol
        Si, Se = self.compute_COD(self._inf), self.compute_COD(self.outs[1])
        return 1 - Qe*Se/(Qi*Si)