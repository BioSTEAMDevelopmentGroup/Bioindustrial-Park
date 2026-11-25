#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FT unit operations (Ostadi/Hillestad PBtL), BioSTEAM."""
import math
import os
import sys
import biosteam as bst
import flexsolve as flx
from biosteam import Unit
from biosteam.units.decorators import cost
import thermosteam as tmo

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from ft_chemicals import (
    FT_PARAFFINS, FT_OLEFINS, FT_LIGHTS, FT_LIQUIDS,
)

BAR = 1e5  # Pa


@cost(basis='Water removal rate', ID='Dryer', units='kg/hr',
      cost=3294700, S=12233, CE=522, n=0.8, BM=1.7)
class BiomassDryer(Unit):
    """Steam-heated dryer; duty is the enthalpy rise of solids and evaporated water
    (~2.8 GJ/t H2O, the Hillestad et al. 2018 basis). Vapor and dried solids leave at
    T_out, heated by steam at T_steam."""
    _N_ins = 1
    _N_outs = 2
    _units = {'Water removal rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(),
                 target_moisture=0.05, T_steam=150 + 273.15, T_out=105 + 273.15):
        Unit.__init__(self, ID, ins, outs)
        self.target_moisture = target_moisture
        self.T_steam = T_steam
        self.T_out = T_out

    def _run(self):
        feed = self.ins[0]
        dried, vapor = self.outs
        dry_mass = feed.F_mass - feed.imass['Water']
        final_water = dry_mass * self.target_moisture / (1.0 - self.target_moisture)
        water_to_remove = max(0.0, feed.imass['Water'] - final_water)
        dried.copy_like(feed)
        dried.phases = ('l', 's')
        dried.imass['s', 'Water'] = 0.0
        dried.imass['l', 'Water'] = final_water
        vapor.empty()
        vapor.imass['Water'] = water_to_remove
        vapor.phase = 'g'
        dried.T = vapor.T = self.T_out

    def _design(self):
        self.design_results['Water removal rate'] = self.outs[1].F_mass
        self.add_heat_utility(self.H_out - self.H_in, self.T_out)

    @property
    def specific_heat(self):
        """Dryer heat per tonne of water removed [GJ/t]."""
        water = self.outs[1].F_mass
        return (self.H_out - self.H_in) / water / 1e3 if water else 0.0


@cost(basis='Grinding flow', ID='Grinder', units='kg/hr',
      cost=13329690, S=94697, CE=522, n=0.6, BM=1.7)
class BiomassGrinder(Unit):
    _N_ins = 1
    _N_outs = 1
    _units = {'Grinding flow': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), power_demand=50):
        Unit.__init__(self, ID, ins, outs)
        self.power_demand = power_demand  # kWh/ton

    @property
    def power(self):
        """Grinding power [kW]."""
        return (self.power_demand / 1000.0) * self.ins[0].F_mass

    def _run(self):
        # Grinding work leaves as heat with the mill sweep air (not tracked)
        self.outs[0].copy_like(self.ins[0])

    def _design(self):
        self.power_utility.consumption = self.power
        self.design_results['Grinding flow'] = self.ins[0].F_mass


# Gibbs species (C–H–O). N2/H2S held fixed (H2S in BST Gibbs breaks S balance).
_GIBBS_IDS = ['H2', 'CO', 'CO2', 'Water', 'CH4', 'O2']
_RWGS_GIBBS_IDS = ['H2', 'CO', 'CO2', 'Water', 'O2']  # no CH4 — blocks methanation
_FROZEN_GAS_IDS = ['N2', 'H2S']
_gas_thermo = None
_rwgs_thermo = None


def _get_gas_thermo():
    global _gas_thermo
    if _gas_thermo is None:
        _gas_thermo = tmo.Thermo(tmo.Chemicals(_GIBBS_IDS))
    return _gas_thermo


def _get_rwgs_thermo():
    global _rwgs_thermo
    if _rwgs_thermo is None:
        _rwgs_thermo = tmo.Thermo(tmo.Chemicals(_RWGS_GIBBS_IDS))
    return _rwgs_thermo


def _gibbs_equilibrate(stream, T, P, species=None, thermo=None):
    """Minimize G for C–H–O gases via biosteam.minimize_Gibbs_free_energy (COBYLA)."""
    from biosteam.units.design_tools.Gibbs_equilibrium_reaction import (
        minimize_Gibbs_free_energy,
    )
    ids = species or _GIBBS_IDS
    gas = tmo.Stream(None, thermo=thermo or _get_gas_thermo())
    for ID in ids:
        if ID in stream.chemicals:
            gas.imol[ID] = max(0.0, float(stream.imol[ID]))
    if gas.F_mass <= 0:
        return
    gas.T = T
    gas.P = P
    gas.phase = 'g'
    minimize_Gibbs_free_energy(gas, method='COBYLA')
    for ID in ids:
        if ID in stream.chemicals:
            stream.imol[ID] = max(0.0, float(gas.imol[ID]))


def _reform_hydrocarbons(stream, min_water=0.0, allow_pox=True, min_co2=0.0):
    """Convert paraffins + olefins to syngas at EF ~1400 C.

    min_water / min_co2: leave at least this much H2O / CO2 (kmol/h).
    Real EF: steam reforming is preferred; flame-zone CO2 is not fully dry-reformed
    before quench. allow_pox=False reserves O2 for a dedicated oxidation step.
    """
    def _water_for_sr():
        return max(0.0, float(stream.imol['Water']) - float(min_water))

    def _co2_for_dr():
        return max(0.0, float(stream.imol['CO2']) - float(min_co2))

    # Paraffins C_nH_{2n+2}
    for n, name in enumerate(FT_PARAFFINS, start=1):
        if name not in stream.chemicals:
            continue
        remaining = float(stream.imol[name])
        if remaining <= 0:
            continue
        if allow_pox and float(stream.imol['O2']) > 0:
            extent = min(remaining, 2.0 * float(stream.imol['O2']) / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['O2'] -= 0.5 * n * extent
                stream.imol['CO'] += n * extent
                stream.imol['H2'] += (n + 1) * extent
        remaining = float(stream.imol[name])
        water_sr = _water_for_sr()
        if remaining > 0 and water_sr > 0:
            extent = min(remaining, water_sr / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['Water'] -= n * extent
                stream.imol['CO'] += n * extent
                stream.imol['H2'] += (2 * n + 1) * extent
        remaining = float(stream.imol[name])
        co2_dr = _co2_for_dr()
        if remaining > 0 and co2_dr > 0:
            extent = min(remaining, co2_dr / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['CO2'] -= n * extent
                stream.imol['CO'] += 2 * n * extent
                stream.imol['H2'] += (n + 1) * extent
        remaining = float(stream.imol[name])
        if remaining > 0:
            _reform_with_reserves(stream, name, n, 2 * n + 2)
    # Olefins C_nH_{2n}
    for i, name in enumerate(FT_OLEFINS):
        n = i + 2
        if name not in stream.chemicals:
            continue
        remaining = float(stream.imol[name])
        if remaining <= 0:
            continue
        if allow_pox and float(stream.imol['O2']) > 0:
            extent = min(remaining, 2.0 * float(stream.imol['O2']) / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['O2'] -= 0.5 * n * extent
                stream.imol['CO'] += n * extent
                stream.imol['H2'] += n * extent
        remaining = float(stream.imol[name])
        water_sr = _water_for_sr()
        if remaining > 0 and water_sr > 0:
            extent = min(remaining, water_sr / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['Water'] -= n * extent
                stream.imol['CO'] += n * extent
                stream.imol['H2'] += 2 * n * extent
        remaining = float(stream.imol[name])
        co2_dr = _co2_for_dr()
        if remaining > 0 and co2_dr > 0:
            extent = min(remaining, co2_dr / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['CO2'] -= n * extent
                stream.imol['CO'] += 2 * n * extent
                stream.imol['H2'] += n * extent
        remaining = float(stream.imol[name])
        if remaining > 0:
            _reform_with_reserves(stream, name, n, 2 * n)


def _reform_with_reserves(stream, name, n, h):
    """Reform what is left of C_nH_h with the reserved steam, then CO2.
    Anything that still lacks an oxygen source stays unreformed."""
    remaining = float(stream.imol[name])
    extent = min(remaining, float(stream.imol['Water']) / n)
    if extent > 0:
        stream.imol[name] -= extent
        stream.imol['Water'] -= n * extent
        stream.imol['CO'] += n * extent
        stream.imol['H2'] += (n + 0.5 * h) * extent
    remaining = float(stream.imol[name])
    extent = min(remaining, float(stream.imol['CO2']) / n)
    if extent > 0:
        stream.imol[name] -= extent
        stream.imol['CO2'] -= n * extent
        stream.imol['CO'] += 2 * n * extent
        stream.imol['H2'] += 0.5 * h * extent


def _burn_oxygen(stream):
    """Consume leftover O2: CO→CO2 then H2→H2O (forms CO2/steam for WGS)."""
    o2 = float(stream.imol['O2'])
    if o2 <= 1e-12:
        return
    ext_co = min(float(stream.imol['CO']), 2.0 * o2)
    if ext_co > 0:
        stream.imol['CO'] -= ext_co
        stream.imol['O2'] -= 0.5 * ext_co
        stream.imol['CO2'] += ext_co
        o2 = float(stream.imol['O2'])
    ext_h2 = min(float(stream.imol['H2']), 2.0 * o2)
    if ext_h2 > 0:
        stream.imol['H2'] -= ext_h2
        stream.imol['O2'] -= 0.5 * ext_h2
        stream.imol['Water'] += ext_h2
    # Any residual O2: leave it for Gibbs (rare)


class H2RatioAdjuster(Unit):
    """Add H2 so outlet H2/CO equals target (Ostadi FOCAPD / PFD setpoints).

    Set allow_import=False for conventional BtL (no green/imported H2 between stages).
    """
    _N_ins = 2
    _N_outs = 1

    def _init(self, target_H2_CO=1.92, T_h2=100 + 273.15, allow_import=True):
        self.target_H2_CO = target_H2_CO
        self.T_h2 = T_h2
        self.allow_import = bool(allow_import)

    def _run(self):
        syngas, h2_src = self.ins
        out = self.outs[0]
        co = float(syngas.imol['CO'])
        h2_have = float(syngas.imol['H2'])
        if self.allow_import and co > 0:
            h2_need = max(0.0, self.target_H2_CO * co - h2_have)
        else:
            h2_need = 0.0
        h2_src.empty()
        h2_src.imol['H2'] = h2_need
        h2_src.phase = 'g'
        h2_src.T = self.T_h2
        h2_src.P = syngas.P
        out.mix_from([syngas, h2_src], energy_balance=False)
        out.phase = 'g'
        out.P = syngas.P
        out.H = syngas.H + h2_src.H


_H2_MW = 2.016
_CEPCI_2020 = 596.2
_CEPCI_2022 = 816.0

# Technology presets (both from DOE analyses of 50 t/d central plants). `kWh_per_kg` is
# lifetime-average system electricity at the stated H2 delivery pressure; `cost` is total
# installed capital (direct + indirect) at size `S`, scaled with exponent `n` on
# `size_basis` ('power' in kW or 'H2' in kg/h) and escalated with CEPCI.
ELECTROLYZERS = {
    # DOE Hydrogen Program Record 24005 (2024), current PEM: 57.5 kWh/kg average system
    # use, H2 at 30 bar, ~2,000 $/kW installed (2022$). Only the cathode is pressurized,
    # so O2 leaves near atmospheric pressure.
    'PEM': dict(
        name='PEM', kWh_per_kg=57.5,
        T=80 + 273.15, P_H2=30 * BAR, P_O2=1.3 * BAR,
        cost=2000e3, S=1000.0, size_basis='power', n=1.0, CE=_CEPCI_2022,
    ),
    # Acevedo et al. 2023 (DOE/Strategic Analysis, OSTI 2203367), Table 12, current
    # near-atmospheric alkaline: 52.4 kWh/kg average incl. balance of plant and H2
    # compression to 30 bar; total installed capital 162.4 M$ at 50 t/d and 926.2 M$ at
    # 500 t/d (2020$). Stack at 90 °C, 1.3 bar (Table 2).
    'AEL': dict(
        name='Alkaline', kWh_per_kg=52.4,
        T=90 + 273.15, P_H2=30 * BAR, P_O2=1.3 * BAR,
        cost=926.193e6, S=500e3 / 24, size_basis='H2',
        n=math.log(926.193 / 162.446) / math.log(500 / 50), CE=_CEPCI_2020,
    ),
}


class WaterElectrolyzer(bst.Facility):
    """Water electrolysis supplying all imported H2 (2 H2O -> 2 H2 + O2).

    Runs after the recycle loop converges and sizes itself on the H2 drawn by
    the H2 feed streams. O2 up to the gasifier demand (`gasifier_oxygen`, set by the
    gasifier) is written to `electrolyzer_oxygen`, the rest of the demand to
    `asu_oxygen`; any excess is vented. Auxiliary isentropic compressors raise H2 and
    gasifier O2 from the electrolyzer outlet pressure to P_delivery. Electric input not
    stored in the products (stack and aftercooler heat) goes to cooling water.
    `technology` selects a preset in ELECTROLYZERS ('AEL' or 'PEM').
    Ins: [0] water. Outs: [0] O2 vent.
    """
    network_priority = 0  # before the ASU, which is sized on the O2 left after this unit
    _N_ins = 1
    _N_outs = 1
    _units = {'Electrolyzer power': 'kW', 'H2 production': 'kg/hr'}
    _F_BM_default = {'Electrolyzer': 1.0}
    auxiliary_unit_names = ('h2_compressor', 'o2_compressor')

    def _init(self, technology='PEM', h2_streams=(), electrolyzer_oxygen=None,
              asu_oxygen=None, gasifier_oxygen=None, P_delivery=40 * BAR,
              compressor_eta=0.75):
        self.technology = technology
        self.h2_streams = tuple(h2_streams)
        self.electrolyzer_oxygen = electrolyzer_oxygen
        self.asu_oxygen = asu_oxygen
        self.gasifier_oxygen = gasifier_oxygen
        self.P_delivery = P_delivery        # FT/shift/gasifier pressure
        self.compressor_eta = compressor_eta

    @property
    def O2_demand(self):
        """[float] O2 required by the gasifier [kmol/h]."""
        return float(self.gasifier_oxygen.imol['O2']) if self.gasifier_oxygen is not None else 0.0

    @property
    def technology(self):
        return self._technology

    @technology.setter
    def technology(self, technology):
        if technology not in ELECTROLYZERS:
            raise ValueError(f'technology must be one of {list(ELECTROLYZERS)}, not {technology!r}')
        self._technology = technology

    @property
    def spec(self):
        return ELECTROLYZERS[self._technology]

    @property
    def H2_production(self):
        return sum(float(s.imol['H2']) for s in self.h2_streams)

    def _run(self):
        spec = self.spec
        water = self.ins[0]
        o2_vent = self.outs[0]
        h2 = self.H2_production
        o2 = 0.5 * h2
        o2_to_gasifier = min(o2, float(self.O2_demand))

        water.empty()
        water.imol['Water'] = h2
        water.phase = 'l'
        water.T = 25 + 273.15
        water.P = max(spec['P_H2'], spec['P_O2'])

        o2_vent.empty()
        o2_vent.imol['O2'] = o2 - o2_to_gasifier
        o2_vent.phase = 'g'
        o2_vent.T = spec['T']
        o2_vent.P = 1 * BAR

        if self.electrolyzer_oxygen is not None:
            self.electrolyzer_oxygen.imol['O2'] = o2_to_gasifier
        if self.asu_oxygen is not None:
            self.asu_oxygen.imol['O2'] = max(0.0, float(self.O2_demand) - o2_to_gasifier)

    def _compressor(self, name, chemical, kmol, P_in):
        """Isentropic compression to P_delivery; intercooled stages of at most
        `max_stage_ratio` (cooled back to the electrolyzer temperature)."""
        feed = tmo.Stream(None, thermo=self.thermo, phase='g',
                          T=self.spec['T'], P=P_in, **{chemical: kmol})
        ratio = self.P_delivery / P_in
        n_stages = max(1, math.ceil(math.log(ratio) / math.log(self.max_stage_ratio) - 1e-9))
        # Electric-motor drivers: BioSTEAM's default steam-turbine driver books the
        # power as both consumed and produced
        if n_stages == 1:
            unit = self.auxiliary(name, bst.IsentropicCompressor, ins=feed,
                                  P=self.P_delivery, eta=self.compressor_eta,
                                  driver='Electric motor')
            unit.simulate()
        else:
            unit = self.auxiliary(name, bst.MultistageCompressor, ins=feed,
                                  pr=ratio ** (1 / n_stages), n_stages=n_stages,
                                  eta=self.compressor_eta)
            unit.simulate()
            for stage in unit.compressors:
                stage.driver = 'Electric motor'
            unit.simulate()
        return unit

    max_stage_ratio = 3.5

    @property
    def compression_power(self):
        return sum(c.power_utility.rate
                   for c in (self.h2_compressor, self.o2_compressor) if c)

    def _design(self):
        spec = self.spec
        self.h2_compressor = self.o2_compressor = None
        h2 = self.H2_production
        h2_kg = h2 * _H2_MW
        power_kW = h2_kg * spec['kWh_per_kg']
        self.design_results['Electrolyzer power'] = power_kW
        self.design_results['H2 production'] = h2_kg
        o2 = (float(self.electrolyzer_oxygen.imol['O2'])
              if self.electrolyzer_oxygen is not None else 0.0)
        if h2 > 0 and spec['P_H2'] < self.P_delivery:
            self._compressor('h2_compressor', 'H2', h2, spec['P_H2'])
        if o2 > 0 and spec['P_O2'] < self.P_delivery:
            self._compressor('o2_compressor', 'O2', o2, spec['P_O2'])
        self.power_utility.consumption = power_kW  # compressor power is added as auxiliary
        power = power_kW + self.compression_power
        products = [*self.h2_streams, self.outs[0]]
        if self.electrolyzer_oxygen is not None:
            products.append(self.electrolyzer_oxygen)
        aux_duty = sum(hu.duty for c in (self.h2_compressor, self.o2_compressor) if c
                       for hu in c.heat_utilities)
        self.cooling_duty = (sum(s.Hnet for s in products) - self.ins[0].Hnet
                             - power * 3600.0 - aux_duty)
        if h2 > 0:
            self.add_heat_utility(self.cooling_duty, spec['T'])

    def _cost(self):
        spec = self.spec
        D = self.design_results
        size = D['Electrolyzer power'] if spec['size_basis'] == 'power' else D['H2 production']
        self.baseline_purchase_costs['Electrolyzer'] = (
            spec['cost'] * (size / spec['S']) ** spec['n'] * bst.CE / spec['CE']
            if size > 0 else 0.0
        )


@cost(basis='O2 production', ID='Air separation unit', units='kg/hr',
      cost=36.8e6 * 1.326, S=76.6e3, CE=550.8, n=0.5, BM=1.0)
class AirSeparationUnit(bst.Facility):
    """Cryogenic ASU supplying the gasifier O2 not covered by the electrolyzer.

    Sized on `oxygen` (the ASU O2 feed to the gasifier, set by the electrolyzer) after
    the recycle loop converges. Installed cost incl. O2 compression from Hannula 2016
    (Table 2: 36.8 M EUR-2010 at 76.6 t/h O2, exponent 0.5; 1.326 USD/EUR).
    Power: `kWh_per_kg` delivered at gasifier pressure. Ins: [0] air. Outs: [0] N2 vent.
    """
    network_priority = 1
    _N_ins = 1
    _N_outs = 1
    _units = {'O2 production': 'kg/hr'}

    def _init(self, oxygen=None, kWh_per_kg=0.35, O2_recovery=0.95):
        self.oxygen = oxygen
        self.kWh_per_kg = kWh_per_kg
        self.O2_recovery = O2_recovery

    def _run(self):
        air, vent = self.ins[0], self.outs[0]
        O2 = float(self.oxygen.imol['O2']) if self.oxygen is not None else 0.0
        O2_air = O2 / self.O2_recovery
        air.empty()
        air.imol['O2'] = O2_air
        air.imol['N2'] = O2_air * 79 / 21
        air.phase = 'g'
        vent.empty()
        vent.imol['O2'] = O2_air - O2
        vent.imol['N2'] = air.imol['N2']
        vent.phase = 'g'

    def _design(self):
        O2_kg = self.oxygen.F_mass if self.oxygen is not None else 0.0
        self.design_results['O2 production'] = O2_kg
        self.power_utility.consumption = power = O2_kg * self.kWh_per_kg
        if O2_kg:  # compressor heat to cooling water
            H_out = self.oxygen.Hnet + self.outs[0].Hnet
            self.add_heat_utility(H_out - self.ins[0].Hnet - power * 3600.0, 40 + 273.15)

    def _cost(self):
        if self.design_results['O2 production'] > 0:
            self._decorated_cost()
        else:
            self.baseline_purchase_costs.clear()


@cost(basis='Flow rate', ID='Gasifier', units='kg/hr',
      cost=45000000, S=80000, CE=522, n=0.7, BM=2.7)
class EntrainedFlowGasifier(Unit):
    """Ins: [0] biomass, [1] steam, [2] O2, [3] recycle. Outs: [0] syngas, [1] slag.

    Real EF behavior:
      - Recycle CH4/lights/CO2 are reformed; injected steam is partly reserved for WGS
      - Leftover O2 forms CO2/H2O by oxidation
      - Hot C–H–O equilibrium is only partially approached (finite residence time /
        quench freeze) so outlet CO2 is nonzero — full Gibbs would over-convert
        CO2 via reverse WGS when recycle H2 is high

    Adiabatic (Hillestad et al. 2018, Sec. 3.1): the O2 flow in ins[2] is solved so
    that the gas leaves at T_op with no external heat; heat losses are neglected.
    With adiabatic=False, ins[2] is used as given and the outlet is still set to T_op.
    """
    _N_ins = 4
    _N_outs = 2
    _units = {'Flow rate': 'kg/hr'}

    def _init(self, T_op=1400 + 273.15, P_op=40 * BAR,
              steam_reserve_frac=0.55, eq_approach=0.25, adiabatic=True):
        self.T_op = T_op
        self.P_op = P_op
        self.steam_reserve_frac = steam_reserve_frac
        # 0 = freeze after oxidation/reform; 1 = full hot Gibbs
        self.eq_approach = eq_approach
        self.adiabatic = adiabatic
        self.biomass_rxn = tmo.Reaction(
            'DryBiomass + 0.1965 O2 -> CO + 0.695 H2 + 0.0014 N2 + 0.0007 H2S',
            reactant='DryBiomass', X=1.0,
        )

    def _partial_gibbs(self, syngas):
        """Move a fraction of the way from current state to hot Gibbs equilibrium."""
        ids = ['H2', 'CO', 'CO2', 'Water', 'CH4', 'O2']
        before = {ID: float(syngas.imol[ID]) for ID in ids if ID in syngas.chemicals}
        _gibbs_equilibrate(syngas, self.T_op, self.P_op)
        a = min(max(float(self.eq_approach), 0.0), 1.0)
        for ID, n0 in before.items():
            n_eq = float(syngas.imol[ID])
            syngas.imol[ID] = n0 + a * (n_eq - n0)

    def _heat_deficit(self, O2):
        """Hnet(out) - Hnet(in) [kJ/h] at O2 [kmol/h]; positive means heat is missing."""
        self.ins[2].imol['O2'] = O2
        self._gasify()
        return sum(s.Hnet for s in self.outs) - sum(s.Hnet for s in self.ins)

    def _run(self):
        oxygen = self.ins[2]
        oxygen.phase = 'g'
        if not self.adiabatic:
            self._gasify()
            return
        # Lower bound: O2 consumed by the biomass decomposition reaction itself
        O2_per_biomass = -float(self.biomass_rxn.istoichiometry['O2'])
        lb = 1.001 * O2_per_biomass * float(self.ins[0].imol['DryBiomass'])
        x0 = max(float(oxygen.imol['O2']), lb + 1.0)
        f0 = self._heat_deficit(x0)
        if f0 > 0:
            lb, y_lb = x0, f0
            ub = 1.25 * x0
            y_ub = self._heat_deficit(ub)
            while y_ub > 0:
                lb, y_lb = ub, y_ub
                ub *= 1.25
                y_ub = self._heat_deficit(ub)
        else:
            ub, y_ub = x0, f0
            y_lb = self._heat_deficit(lb)
            if y_lb <= 0:  # excess heat even at minimum O2
                self._heat_deficit(lb)
                return
        O2 = flx.IQ_interpolation(self._heat_deficit, lb, ub, y_lb, y_ub, xtol=0.01, ytol=1e3)
        self._heat_deficit(O2)  # the solver may stop without evaluating at its answer

    @property
    def O2_equivalence_ratio(self):
        """Lambda: O2 fed / O2 for complete combustion of the fuel fed (biomass and
        recycled tail gas); Hillestad et al. 2018 use 0.2-0.4."""
        O2_stoich = 0.0
        for fuel in (self.ins[0], self.ins[3]):  # H2O and CO2 cancel out
            O2_stoich += (fuel.get_atomic_flow('C') + 0.25 * fuel.get_atomic_flow('H')
                          - 0.5 * fuel.get_atomic_flow('O'))
        return float(self.ins[2].imol['O2']) / O2_stoich if O2_stoich else 0.0

    def _gasify(self):
        syngas, slag = self.outs
        syngas.mix_from(self.ins, energy_balance=False)
        syngas.phase = 'g'
        slag.empty()
        if 'Ash' in self.chemicals:
            slag.imass['Ash'] = syngas.imass['Ash']
            syngas.imass['Ash'] = 0.0
        if syngas.imol['DryBiomass'] > 0:
            self.biomass_rxn.force_reaction(syngas)
            syngas.imol['DryBiomass'] = 0.0

        steam_feed = float(self.ins[1].imol['Water']) if self.ins[1] is not None else 0.0
        min_water = float(self.steam_reserve_frac) * steam_feed

        # Flame/oxidation zone: O2 → CO2. Then steam-preferring reform of recycle.
        _burn_oxygen(syngas)
        co2_flame = float(syngas.imol['CO2'])
        min_co2 = 0.60 * co2_flame  # keep most flame-zone CO2 (not fully dry-reformed)
        _reform_hydrocarbons(
            syngas, min_water=min_water, allow_pox=False, min_co2=min_co2,
        )
        syngas.phase = 'g'
        syngas.T = self.T_op
        syngas.P = self.P_op
        self._partial_gibbs(syngas)
        # Freeze: hot reverse-WGS must not erase flame-zone CO2 before quench
        if float(syngas.imol['CO2']) < min_co2:
            _forward_wgs(syngas, min_co2 - float(syngas.imol['CO2']))
        _reform_hydrocarbons(
            syngas, min_water=min_water, allow_pox=False, min_co2=min_co2,
        )
        self._partial_gibbs(syngas)
        if float(syngas.imol['CO2']) < min_co2:
            _forward_wgs(syngas, min_co2 - float(syngas.imol['CO2']))

        slag.phase = 'l'
        slag.T = self.T_op
        slag.P = self.P_op

    def _design(self):
        self.design_results['Flow rate'] = sum(s.F_mass for s in self.ins)


def _forward_wgs(stream, extent):
    """CO + H2O -> CO2 + H2 (cooling / forward WGS), limited by the CO and steam present."""
    extent = min(max(0.0, float(extent)), float(stream.imol['CO']), float(stream.imol['Water']))
    if extent <= 0:
        return 0.0
    stream.imol['CO'] -= extent
    stream.imol['Water'] -= extent
    stream.imol['CO2'] += extent
    stream.imol['H2'] += extent
    return extent


def _reverse_wgs(stream, extent):
    """CO2 + H2 -> CO + H2O (RWGS)."""
    extent = max(0.0, float(extent))
    if extent <= 0:
        return 0.0
    extent = min(extent, float(stream.imol['CO2']), float(stream.imol['H2']))
    if extent <= 0:
        return 0.0
    stream.imol['CO2'] -= extent
    stream.imol['H2'] -= extent
    stream.imol['CO'] += extent
    stream.imol['Water'] += extent
    return extent


@cost(basis='Flow rate', ID='RWGS', units='kg/hr',
      cost=5000000, S=80000, CE=522, n=0.7, BM=2.5)
class RWGSReactor(Unit):
    """
    Ins: [0] raw syngas from gasifier, [1] H2 make-up, [2] shift steam.
    Outs: [0] conditioned syngas.

    Adiabatic: the outlet temperature follows from the energy balance. In pbtl mode
    this is the chemical quench of Hillestad et al. 2018 (Sec. 2.8; ~1400 -> ~1000 C);
    in btl_wgs mode the shift heat leaves with the gas and is recovered in the
    waste-heat boiler.

    Modes
    -----
    btl_wgs (default, real-world conventional BtL/FT):
        Forward WGS: CO + H2O -> CO2 + H2 to approach target H2/CO (makes AGR CO2).
        Steam from the boiler (ins[2]) tops the inlet water up to
        steam_to_shifted_CO mol per mol CO shifted (1.5 keeps the outlet below WGS
        equilibrium at a ~400 C sour-shift exit). Imported H2 may finish the ratio
        trim. No reverse-WGS CO2 recovery. Extra CO2 is stripped in AGR
        (Kreutz/Leibbrandt/NETL-type).

    conventional:
        No CO2→CO recovery; H2 imported only to hit H2/CO (hybrid).

    pbtl:
        Imported-H2 reverse WGS converts CO2 → CO down to residual_CO2 for AGR
        (Hillestad/Ostadi power-and-biomass-to-liquids).
    """
    _N_ins = 3
    _N_outs = 1
    _units = {'Flow rate': 'kg/hr'}

    def _init(self, T_op=1000 + 273.15, P_op=40 * BAR, target_H2_CO=1.92,
              residual_CO2=305.0, mode='btl_wgs', steam_to_shifted_CO=1.5, adiabatic=True):
        self.T_op = T_op  # outlet temperature only if adiabatic=False
        self.P_op = P_op
        self.target_H2_CO = target_H2_CO
        self.residual_CO2 = residual_CO2
        self.steam_to_shifted_CO = steam_to_shifted_CO
        self.adiabatic = adiabatic
        self.mode = str(mode).lower().strip()
        self.CO2_conversion_achieved = 0.0
        self.CO2_before_RWGS = 0.0
        self.CO2_after_quench = 0.0
        self.cooling_WGS_extent = 0.0
        self.RWGS_extent = 0.0
        self.WGS_extent = 0.0

    def _run(self):
        syngas, h2, steam = self.ins
        out = self.outs[0]
        steam.empty()
        steam.phase = 'g'
        out.copy_like(syngas)
        out.phase = 'g'
        out.T = self.T_op
        out.P = self.P_op
        self._condition(out, h2, steam)
        if self.adiabatic:
            out.H = sum(s.Hnet for s in self.ins) - out.Hf

    def _condition(self, out, h2, steam):
        self.cooling_WGS_extent = 0.0
        self.WGS_extent = 0.0
        self.CO2_after_quench = float(out.imol['CO2'])
        co2_0 = float(out.imol['CO2'])
        self.CO2_before_RWGS = co2_0

        h2.empty()
        h2.phase = 'g'
        h2.T = 100 + 273.15
        h2.P = self.P_op
        h2.imol['H2'] = 0.0

        mode = self.mode
        if mode == 'pbtl':
            convert = max(0.0, co2_0 - float(self.residual_CO2))
            co_after = float(out.imol['CO']) + convert
            h2_after_no_add = float(out.imol['H2']) - convert
            h2_want = self.target_H2_CO * co_after if co_after > 1e-12 else 0.0
            h2_add = max(0.0, h2_want - h2_after_no_add, convert)
            h2.imol['H2'] = h2_add
            out.imol['H2'] = float(out.imol['H2']) + h2_add
            self.RWGS_extent = _reverse_wgs(out, convert)
            self.CO2_conversion_achieved = (
                self.RWGS_extent / co2_0 if co2_0 > 1e-12 else 0.0
            )
            co = float(out.imol['CO'])
            if co > 1e-12:
                h2_have = float(out.imol['H2'])
                h2_want = self.target_H2_CO * co
                delta = h2_want - h2_have
                out.imol['H2'] = h2_want
                h2.imol['H2'] = max(0.0, float(h2.imol['H2']) + delta)

        elif mode in ('btl_wgs', 'btl', 'real'):
            # Forward WGS to raise H2/CO: CO + H2O -> CO2 + H2 (adds steam if needed)
            # Then trim with imported H2 if still short of target (keeps FT PFD ratios).
            self.RWGS_extent = 0.0
            self.CO2_conversion_achieved = 0.0
            R = float(self.target_H2_CO)
            co = float(out.imol['CO'])
            h2_have = float(out.imol['H2'])
            if co > 1e-12 and R > 0:
                x = (R * co - h2_have) / (1.0 + R)
                x = max(0.0, min(x, co * 0.95))
            else:
                x = 0.0
            steam.imol['Water'] = max(0.0, self.steam_to_shifted_CO * x - float(out.imol['Water']))
            out.imol['Water'] += float(steam.imol['Water'])
            self.WGS_extent = _forward_wgs(out, x)
            self.cooling_WGS_extent = self.WGS_extent
            co = float(out.imol['CO'])
            if co > 1e-12:
                h2_have = float(out.imol['H2'])
                h2_want = R * co
                delta = h2_want - h2_have
                if delta > 0:
                    out.imol['H2'] = h2_want
                    h2.imol['H2'] = delta
                else:
                    h2.imol['H2'] = 0.0

        else:
            # hybrid conventional: leave CO2; import H2 for ratio
            self.RWGS_extent = 0.0
            self.CO2_conversion_achieved = 0.0
            co = float(out.imol['CO'])
            if co > 1e-12:
                h2_have = float(out.imol['H2'])
                h2_want = self.target_H2_CO * co
                delta = h2_want - h2_have
                out.imol['H2'] = h2_want
                h2.imol['H2'] = max(0.0, delta)

    def _design(self):
        self.design_results['Flow rate'] = sum(s.F_mass for s in self.ins)


@cost(basis='Flow rate', ID='Cyclone', units='kg/hr',
      cost=150000, S=50000, CE=522, n=0.6, BM=2.0)
class SyngasCyclone(bst.SolidsSeparator):
    _N_ins = 1
    _N_outs = 2
    _units = {'Flow rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), split=0.999):
        bst.SolidsSeparator.__init__(self, ID, ins, outs, split=split)

    def _run(self):
        super()._run()
        self.outs[0].T = self.outs[1].T = self.ins[0].T
        self.outs[0].P = self.outs[1].P = self.ins[0].P

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass


@cost(basis='Heat duty', ID='Waste Heat Boiler', units='kJ/hr',
      cost=2500000, S=1e7, CE=522, n=0.7, BM=2.2)
class SyngasWasteHeatBoiler(Unit):
    _N_ins = 2
    _N_outs = 2
    _units = {'Heat duty': 'kJ/hr'}

    def __init__(self, ID='', ins=None, outs=(),
                 T_syngas_out=400 + 273.15,
                 T_steam_out=700 + 273.15,
                 P_steam=117 * BAR):
        Unit.__init__(self, ID, ins, outs)
        self.T_syngas_out = T_syngas_out
        self.T_steam_out = T_steam_out
        self.P_steam = P_steam

    def _run(self):
        hot_syngas, bfw = self.ins
        cooled_syngas, superheated_steam = self.outs
        cooled_syngas.copy_like(hot_syngas)
        cooled_syngas.T = self.T_syngas_out
        heat_available = hot_syngas.H - cooled_syngas.H
        basis = tmo.Stream(None, thermo=self.thermo)
        basis.imol['Water'] = 1.0
        basis.phase = 'l'
        basis.T = bfw.T
        basis.P = bfw.P
        h_water = basis.H
        basis.phase = 'g'
        basis.T = self.T_steam_out
        basis.P = self.P_steam
        heat_needed_per_kmol = basis.H - h_water
        moles_of_steam = (heat_available / heat_needed_per_kmol
                          if heat_needed_per_kmol > 0 else 0.0)
        bfw.empty()
        bfw.imol['Water'] = moles_of_steam
        bfw.phase = 'l'
        bfw.T = 32 + 273.15
        bfw.P = self.P_steam
        superheated_steam.empty()
        superheated_steam.imol['Water'] = moles_of_steam
        superheated_steam.phase = 'g'
        superheated_steam.T = self.T_steam_out
        superheated_steam.P = self.P_steam

    def _design(self):
        self.design_results['Heat duty'] = self.ins[0].H - self.outs[0].H


class SyngasCooler(bst.HXutility):
    def __init__(self, ID='', ins=None, outs=(), T_out=180 + 273.15):
        bst.HXutility.__init__(self, ID, ins, outs, T=T_out)


@cost(basis='Flow rate', ID='Syngas Filter', units='kg/hr',
      cost=200000, S=50000, CE=522, n=0.6, BM=2.3)
class SyngasFilter(bst.SolidsSeparator):
    _units = {'Flow rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), split=0.999):
        bst.SolidsSeparator.__init__(self, ID, ins, outs, split=split)

    def _run(self):
        super()._run()
        self.outs[0].T = self.outs[1].T = self.ins[0].T
        self.outs[0].P = self.outs[1].P = self.ins[0].P

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass


@cost(basis='Flow rate', ID='Water Wash', units='kg/hr',
      cost=600000, S=50000, CE=522, n=0.6, BM=2.0)
class SyngasWaterWash(Unit):
    _N_ins = 2
    _N_outs = 2
    _units = {'Flow rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)

    def _run(self):
        syngas, wash_water = self.ins
        clean_gas, wastewater = self.outs
        clean_gas.copy_like(syngas)
        clean_gas.T = 80 + 273.15
        wastewater.copy_like(wash_water)
        wastewater.T = 80 + 273.15

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass
        # Circulating wash water is cooled against cooling water
        self.add_heat_utility(self.H_out - self.H_in, 80 + 273.15)


class SyngasChiller(bst.HXutility):
    def __init__(self, ID='', ins=None, outs=(), T_out=40 + 273.15):
        bst.HXutility.__init__(self, ID, ins, outs, T=T_out)


@cost(basis='Flow rate', ID='Knockout Drum', units='kg/hr',
      cost=150000, S=50000, CE=522, n=0.6, BM=2.0)
class KnockOutDrum(bst.Flash):
    _units = {'Flow rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=()):
        bst.Flash.__init__(self, ID, ins, outs, T=40 + 273.15, P=40 * BAR)

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass
        self.add_heat_utility(self.H_out - self.H_in, self.T)


@cost(basis='Flow rate', ID='Acid Gas Removal', units='kg/hr',
      cost=4810000, S=24123, CE=522, n=0.6, BM=4.3)
class AcidGasRemoval(bst.Splitter):
    """H2S fully removed; CO2 removed by fraction CO2_removal."""
    _N_ins = 1
    _N_outs = 2
    _units = {'Flow rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), CO2_removal=0.96):
        bst.Splitter.__init__(self, ID, ins, outs, split=0.5)
        self.CO2_removal = CO2_removal

    def _run(self):
        self.split = [1.0] * len(self.chemicals)

        def idx(chem_id):
            try:
                return self.chemicals.index(chem_id)
            except ValueError:
                return None

        co2_idx = idx('CO2')
        if co2_idx is not None:
            self.split[co2_idx] = 1.0 - self.CO2_removal
        h2s_idx = idx('H2S')
        if h2s_idx is not None:
            self.split[h2s_idx] = 0.0
        water_idx = idx('Water')
        if water_idx is not None:
            self.split[water_idx] = 0.0
        super()._run()
        self.outs[0].P = self.outs[1].P = self.ins[0].P
        self.outs[0].phase = 'g'

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass


class SyngasPreheater(bst.HXutility):
    def __init__(self, ID='', ins=None, outs=(), T_out=210 + 273.15):
        bst.HXutility.__init__(self, ID, ins, outs, T=T_out)


@cost(basis='Flow rate', ID='FT Reactor', units='kg/hr',
      cost=28000000, S=65000, CE=522, n=0.6, BM=3.1)
class FischerTropschReactor(Unit):
    """
    Ostadi/Hillestad FT: dual ASF (paraffins α1, olefins α2=0.7·α1).
    α1 from Ostadi CERD 2016 Eq. (17); >90% paraffins (Fuel 2018).
    Ins: [0] syngas, [1] H2, [2] BFW. Outs: [0] effluent, [1] steam.
    """
    _N_ins = 3
    _N_outs = 2
    _units = {'Flow rate': 'kg/hr'}

    def _init(self,
              CO_conversion=0.60,
              T_op=210 + 273.15,
              P_op=37 * BAR,
              paraffin_frac=0.92,
              methane_extra=0.02):
        self.CO_conversion = CO_conversion
        self.T_op = T_op
        self.P_op = P_op
        self.paraffin_frac = paraffin_frac
        self.methane_extra = methane_extra  # extra CH4 beyond ASF (paper)

    @staticmethod
    def _alpha1(stream, T, P):
        """
        Chain growth α1: Song et al. (2004) / Yermakova–Anikeev form
        as cited in Ostadi CERD 2016 Eq. (10)–(11).
        """
        F = float(stream.F_mol)
        if F <= 0:
            return 0.90
        y_h2 = float(stream.imol['H2']) / F
        y_co = float(stream.imol['CO']) / F
        den = y_h2 + y_co
        if den <= 0:
            return 0.90
        A, B = 0.2332, 0.6330
        alpha = (A * (y_co / den) + B) * (1.0 - 0.0039 * (T - 533.0))
        # Understoichiometric cobalt LTFT: keep in paper-relevant band
        return min(0.95, max(0.85, alpha))

    @staticmethod
    def _asf_weights(alpha, n_min, n_max):
        # Carbon (mass) fraction W_n = n(1-α)²α^(n-1): each weight is the share of
        # converted CO that ends up in C_n, which is what the reactions below need.
        return [n * (1.0 - alpha) ** 2 * (alpha ** (n - 1)) for n in range(n_min, n_max + 1)]

    def _build_reactions(self, X_co, alpha1, paraffin_frac):
        alpha2 = 0.70 * alpha1  # Fuel 2018: α2 ≈ 70% of α1
        reaction_list = []
        # Extra methane (beyond ASF) — Fuel/Ostadi note higher CH4 selectivity
        x_ch4_extra = X_co * self.methane_extra
        x_main = X_co * (1.0 - self.methane_extra)
        x_par = x_main * paraffin_frac
        x_ole = x_main * (1.0 - paraffin_frac)
        if x_ch4_extra > 0:
            reaction_list.append(
                tmo.Reaction('CO + 3 H2 -> CH4 + Water', reactant='CO', X=x_ch4_extra)
            )
        # Paraffins C1–C20
        w_p = self._asf_weights(alpha1, 1, 20)
        tot_p = sum(w_p)
        for n, (product, w) in enumerate(zip(FT_PARAFFINS, w_p), start=1):
            X_n = x_par * (w / tot_p)
            if X_n <= 0:
                continue
            reaction_list.append(tmo.Reaction(
                f'{n} CO + {2 * n + 1} H2 -> {product} + {n} Water',
                reactant='CO', X=X_n,
            ))
        # Olefins C2–C20
        w_o = self._asf_weights(alpha2, 2, 20)
        tot_o = sum(w_o)
        for n, (product, w) in enumerate(zip(FT_OLEFINS, w_o), start=2):
            X_n = x_ole * (w / tot_o)
            if X_n <= 0:
                continue
            reaction_list.append(tmo.Reaction(
                f'{n} CO + {2 * n} H2 -> {product} + {n} Water',
                reactant='CO', X=X_n,
            ))
        return tmo.ParallelReaction(reaction_list), alpha1, alpha2

    def _run(self):
        syngas, injected_h2, bfw = self.ins
        effluent, steam = self.outs
        effluent.mix_from([syngas, injected_h2], energy_balance=False)
        # Hnet includes heats of formation, so Q_released captures the FT reaction heat
        H_reactants = syngas.Hnet + injected_h2.Hnet
        alpha1 = self._alpha1(effluent, self.T_op, self.P_op)
        self.alpha1 = alpha1
        self.alpha2 = 0.70 * alpha1
        co0 = float(effluent.imol['CO'])
        h20 = float(effluent.imol['H2'])
        # Mean H2 demand (paraffin-dominated)
        w_p = self._asf_weights(alpha1, 1, 20)
        tot_p = sum(w_p)
        h2_par = sum((w / tot_p) * (2 * n + 1) / n for n, w in enumerate(w_p, start=1))
        w_o = self._asf_weights(0.70 * alpha1, 2, 20)
        tot_o = sum(w_o)
        h2_ole = sum((w / tot_o) * 2.0 for w in w_o)  # 2n H2 / n CO = 2
        h2_per_co = (self.paraffin_frac * h2_par
                     + (1.0 - self.paraffin_frac) * h2_ole
                     + self.methane_extra * 3.0)
        if co0 > 0 and h2_per_co > 0:
            x_eff = min(self.CO_conversion, max(0.0, 0.999 * h20 / (h2_per_co * co0)))
        else:
            x_eff = 0.0
        self.CO_conversion_achieved = x_eff
        rxn, _, _ = self._build_reactions(x_eff, alpha1, self.paraffin_frac)
        rxn.force_reaction(effluent)
        for chem in ('H2', 'CO', 'Water'):
            if effluent.imol[chem] < 0:
                effluent.imol[chem] = 0.0
        effluent.phase = 'g'
        effluent.T = self.T_op
        effluent.P = self.P_op
        Q_released = H_reactants - effluent.Hnet
        basis = tmo.Stream(None, thermo=self.thermo)
        basis.imol['Water'] = 1.0
        basis.phase = 'l'
        basis.T = bfw.T
        basis.P = bfw.P
        h_liq = basis.H
        basis.phase = 'g'
        basis.T = self.T_op
        basis.P = self.P_op
        hvap = basis.H - h_liq
        moles_steam = Q_released / hvap if (Q_released > 0 and hvap > 0) else 0.0
        bfw.empty()
        bfw.imol['Water'] = moles_steam
        steam.empty()
        steam.imol['Water'] = moles_steam
        steam.phase = 'g'
        steam.T = self.T_op
        steam.P = self.P_op

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass + self.ins[1].F_mass


@cost(basis='Flow rate', ID='FT Separator', units='kg/hr',
      cost=450000, S=65000, CE=522, n=0.6, BM=3.0)
class FTSeparator(Unit):
    """Outs: [0] C1–C4 + syngas, [1] C5–C20 liquids, [2] water."""
    _N_ins = 1
    _N_outs = 3
    _units = {'Flow rate': 'kg/hr'}

    def _init(self, T_op=50 + 273.15, P_op=35 * BAR):
        self.T_op = T_op
        self.P_op = P_op

    def _run(self):
        feed = self.ins[0]
        gas, organics, water = self.outs
        gas.empty()
        organics.empty()
        water.empty()
        gas_species = ['CO', 'H2', 'CO2', 'N2', 'H2S', 'O2'] + FT_LIGHTS
        for chem in gas_species:
            if chem in self.chemicals:
                gas.imol[chem] = feed.imol[chem]
        for chem in FT_LIQUIDS:
            if chem in self.chemicals:
                organics.imol[chem] = feed.imol[chem]
        if 'Water' in self.chemicals:
            water.imol['Water'] = feed.imol['Water']
        gas.T = organics.T = water.T = self.T_op
        gas.P = organics.P = water.P = self.P_op
        gas.phase = 'g'
        organics.phase = 'l'
        water.phase = 'l'

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass
        feed = self.ins[0]
        self.add_heat_utility(sum(s.H for s in self.outs) - feed.H, self.T_op)


class BoilerTurbogenerator(bst.BoilerTurbogenerator):
    """BioSTEAM boiler-turbogenerator that also expands steam raised in process units.

    Burns ins[1] (tail-gas purge) and supplies all process steam utilities.
    `recovered_steam` holds (steam, boiler feed water) stream pairs, e.g. the syngas
    waste-heat boiler and the FT reactor jackets; their enthalpy rise first displaces
    boiler steam and the surplus goes to the turbogenerator. Natural gas is bought
    only if burned feeds plus recovered steam fall short of the steam demand; the
    electricity deficit (electrolyzer) is imported, since
    satisfy_system_electricity_demand defaults to False here.
    """

    def __init__(self, ID='', ins=None,
                 outs=('emissions', 'rejected_water_and_blowdown', 'ash_disposal'),
                 thermo=None, *, recovered_steam=(), side_steams=(), **kwargs):
        kwargs.setdefault('satisfy_system_electricity_demand', False)
        super().__init__(ID, ins, outs, thermo, **kwargs)
        self.recovered_steam = tuple(recovered_steam)
        self._recovered_steam_utility = None
        # Live steam drawn by process units as feeds (gasifier, shift)
        self.side_steams = tuple(side_steams)
        if self.side_steams:
            self.side_steam = tmo.Stream(None, thermo=self.thermo)

    @property
    def recovered_steam_duty(self):
        """[float] Enthalpy rise of steam raised in process units [kJ/hr]."""
        return sum(max(0.0, steam.H - bfw.H) for steam, bfw in self.recovered_steam)

    def _load_utility_agents(self):
        super()._load_utility_agents()
        self._recovered_steam_utility = None
        duty = self.recovered_steam_duty
        if duty > 0:
            hu = bst.HeatUtility()
            hu.load_agent(self.agent)
            hu.duty = hu.unit_duty = -duty
            hu.flow = hu.cost = 0.0
            self.steam_utilities.append(hu)
            self._recovered_steam_utility = hu

    def _design(self):
        if self.side_steams:
            self.side_steam.mix_from(self.side_steams)
        super()._design()
        credit = self._recovered_steam_utility
        if credit is None:
            return
        process = [hu for hu in self.steam_utilities if hu is not credit]
        self.steam_utilities.remove(credit)
        heating = bst.HeatUtility.sum_by_agent(process)
        for hu in heating:
            hu.reverse()
        self.heat_utilities = [*heating, self.heat_utilities[-1]]
        self.steam_demand.imol['Water'] = sum(hu.flow for hu in process)


class FTHeatExchangerNetwork(bst.HeatExchangerNetwork):
    """HXN that skips costing if pinch synthesis fails (BioSTEAM SAF pattern)."""

    def _cost(self):
        try:
            super()._cost()
        except (AssertionError, ValueError, RuntimeError) as error:
            import warnings
            warnings.warn(
                f'HXN synthesis failed: {error}. Skipping heat-exchanger network costing.',
                stacklevel=2,
            )
            self.baseline_purchase_costs.clear()
            self.installed_costs.clear()
