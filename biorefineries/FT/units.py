#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FT unit operations (Ostadi/Hillestad PBtL), BioSTEAM."""
import os
import sys
import biosteam as bst
from biosteam import Unit
from biosteam.units.decorators import cost
import thermosteam as tmo

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from ft_chemicals import (
    FT_PARAFFINS, FT_OLEFINS, FT_LIGHTS, FT_LIQUIDS,
)

_GJ2kJ = 1e6
BAR = 1e5  # Pa


@cost(basis='Water removal rate', ID='Dryer', units='kg/hr',
      cost=3294700, S=12233, CE=522, n=0.8, BM=1.7)
class BiomassDryer(Unit):
    _N_ins = 1
    _N_outs = 2
    _units = {'Water removal rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(),
                 target_moisture=0.05, T_steam=150 + 273.15):
        Unit.__init__(self, ID, ins, outs)
        self.target_moisture = target_moisture
        self.T_steam = T_steam

    def _run(self):
        feed = self.ins[0]
        dried, vapor = self.outs
        dry_mass = feed.F_mass - feed.imass['Water']
        final_water = dry_mass * self.target_moisture / (1.0 - self.target_moisture)
        water_to_remove = max(0.0, feed.imass['Water'] - final_water)
        dried.copy_like(feed)
        dried.imass['Water'] = final_water
        vapor.empty()
        vapor.imass['Water'] = water_to_remove
        vapor.phase = 'g'
        dried.T = vapor.T = self.T_steam

    def _design(self):
        water_removed_kg = self.outs[1].F_mass
        heat_duty = (2.8 / 1000.0) * water_removed_kg * _GJ2kJ  # 2.8 GJ/t H2O
        self.design_results['Water removal rate'] = water_removed_kg
        self.add_heat_utility(heat_duty, self.T_steam)


@cost(basis='Grinding flow', ID='Grinder', units='kg/hr',
      cost=13329690, S=94697, CE=522, n=0.6, BM=1.7)
class BiomassGrinder(Unit):
    _N_ins = 1
    _N_outs = 1
    _units = {'Grinding flow': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), power_demand=50):
        Unit.__init__(self, ID, ins, outs)
        self.power_demand = power_demand  # kWh/ton

    def _run(self):
        self.outs[0].copy_like(self.ins[0])

    def _design(self):
        feed_mass_kg = self.ins[0].F_mass
        self.power_utility.consumption = (self.power_demand / 1000.0) * feed_mass_kg
        self.design_results['Grinding flow'] = feed_mass_kg


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


def _reform_hydrocarbons(stream):
    """Convert paraffins + olefins to syngas at EF ~1400 C."""
    # Paraffins C_nH_{2n+2}: POX H2=(n+1); SR H2=(2n+1); DR H2=(n+1)
    for n, name in enumerate(FT_PARAFFINS, start=1):
        if name not in stream.chemicals:
            continue
        remaining = float(stream.imol[name])
        if remaining <= 0:
            continue
        if float(stream.imol['O2']) > 0:
            extent = min(remaining, 2.0 * float(stream.imol['O2']) / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['O2'] -= 0.5 * n * extent
                stream.imol['CO'] += n * extent
                stream.imol['H2'] += (n + 1) * extent
        remaining = float(stream.imol[name])
        if remaining > 0 and float(stream.imol['Water']) > 0:
            extent = min(remaining, float(stream.imol['Water']) / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['Water'] -= n * extent
                stream.imol['CO'] += n * extent
                stream.imol['H2'] += (2 * n + 1) * extent
        remaining = float(stream.imol[name])
        if remaining > 0 and float(stream.imol['CO2']) > 0:
            extent = min(remaining, float(stream.imol['CO2']) / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['CO2'] -= n * extent
                stream.imol['CO'] += 2 * n * extent
                stream.imol['H2'] += (n + 1) * extent
        remaining = float(stream.imol[name])
        if remaining > 0:
            stream.imol[name] -= remaining
            stream.imol['CO'] += n * remaining
            stream.imol['H2'] += (2 * n + 1) * remaining
    # Olefins C_nH_{2n} (n=2..20): POX H2=n; SR H2=2n; DR H2=n
    for i, name in enumerate(FT_OLEFINS):
        n = i + 2
        if name not in stream.chemicals:
            continue
        remaining = float(stream.imol[name])
        if remaining <= 0:
            continue
        if float(stream.imol['O2']) > 0:
            extent = min(remaining, 2.0 * float(stream.imol['O2']) / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['O2'] -= 0.5 * n * extent
                stream.imol['CO'] += n * extent
                stream.imol['H2'] += n * extent
        remaining = float(stream.imol[name])
        if remaining > 0 and float(stream.imol['Water']) > 0:
            extent = min(remaining, float(stream.imol['Water']) / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['Water'] -= n * extent
                stream.imol['CO'] += n * extent
                stream.imol['H2'] += 2 * n * extent
        remaining = float(stream.imol[name])
        if remaining > 0 and float(stream.imol['CO2']) > 0:
            extent = min(remaining, float(stream.imol['CO2']) / n)
            if extent > 0:
                stream.imol[name] -= extent
                stream.imol['CO2'] -= n * extent
                stream.imol['CO'] += 2 * n * extent
                stream.imol['H2'] += n * extent
        remaining = float(stream.imol[name])
        if remaining > 0:
            stream.imol[name] -= remaining
            stream.imol['CO'] += n * remaining
            stream.imol['H2'] += 2 * n * remaining


class H2RatioAdjuster(Unit):
    """Add H2 so outlet H2/CO equals target (Ostadi FOCAPD / PFD setpoints)."""
    _N_ins = 2
    _N_outs = 1

    def _init(self, target_H2_CO=1.92, T_h2=100 + 273.15):
        self.target_H2_CO = target_H2_CO
        self.T_h2 = T_h2

    def _run(self):
        syngas, h2_src = self.ins
        out = self.outs[0]
        co = float(syngas.imol['CO'])
        h2_have = float(syngas.imol['H2'])
        h2_need = max(0.0, self.target_H2_CO * co - h2_have) if co > 0 else 0.0
        h2_src.empty()
        h2_src.imol['H2'] = h2_need
        h2_src.phase = 'g'
        h2_src.T = self.T_h2
        h2_src.P = syngas.P
        out.mix_from([syngas, h2_src], energy_balance=False)
        out.phase = 'g'
        out.T = syngas.T
        out.P = syngas.P


@cost(basis='Flow rate', ID='Gasifier', units='kg/hr',
      cost=45000000, S=80000, CE=522, n=0.7, BM=2.7)
class EntrainedFlowGasifier(Unit):
    """Ins: [0] biomass, [1] steam, [2] O2, [3] recycle. Outs: [0] syngas, [1] slag."""
    _N_ins = 4
    _N_outs = 2
    _units = {'Flow rate': 'kg/hr'}

    def _init(self, T_op=1400 + 273.15, P_op=40 * BAR):
        self.T_op = T_op
        self.P_op = P_op
        self.biomass_rxn = tmo.Reaction(
            'DryBiomass + 0.1965 O2 -> CO + 0.695 H2 + 0.0014 N2 + 0.0007 H2S',
            reactant='DryBiomass', X=1.0,
        )

    def _run(self):
        syngas, slag = self.outs
        syngas.mix_from(self.ins, energy_balance=False)
        slag.empty()
        if 'Ash' in self.chemicals:
            slag.imass['Ash'] = syngas.imass['Ash']
            syngas.imass['Ash'] = 0.0
        if syngas.imol['DryBiomass'] > 0:
            self.biomass_rxn.force_reaction(syngas)
            syngas.imol['DryBiomass'] = 0.0
        _reform_hydrocarbons(syngas)
        syngas.phase = 'g'
        syngas.T = self.T_op
        syngas.P = self.P_op
        _gibbs_equilibrate(syngas, self.T_op, self.P_op)
        # Second pass: reform any CH4/C2+ left after Gibbs water/CO2 shift
        _reform_hydrocarbons(syngas)
        _gibbs_equilibrate(syngas, self.T_op, self.P_op)
        slag.phase = 'l'
        slag.T = self.T_op
        slag.P = self.P_op

    def _design(self):
        self.design_results['Flow rate'] = sum(s.F_mass for s in self.ins)


@cost(basis='Flow rate', ID='RWGS', units='kg/hr',
      cost=5000000, S=80000, CE=522, n=0.7, BM=2.5)
class RWGSReactor(Unit):
    """
    Ins: [0] syngas, [1] H2 make-up. Outs: [0] shifted syngas.
    H2 make-up is set so post-Gibbs H2/CO equals target_H2_CO (FOCAPD).
    CH4 is excluded from the Gibbs pool to suppress methanation.
    """
    _N_ins = 2
    _N_outs = 1
    _units = {'Flow rate': 'kg/hr'}

    def _init(self, T_op=1000 + 273.15, P_op=40 * BAR, target_H2_CO=1.92):
        self.T_op = T_op
        self.P_op = P_op
        self.target_H2_CO = target_H2_CO

    def _equilibrate(self, out):
        # Freeze CH4: pull out, Gibbs without it, put back
        ch4 = float(out.imol['CH4']) if 'CH4' in out.chemicals else 0.0
        if 'CH4' in out.chemicals:
            out.imol['CH4'] = 0.0
        _gibbs_equilibrate(
            out, self.T_op, self.P_op,
            species=_RWGS_GIBBS_IDS, thermo=_get_rwgs_thermo(),
        )
        if 'CH4' in out.chemicals:
            out.imol['CH4'] = ch4

    def _run(self):
        syngas, h2 = self.ins
        out = self.outs[0]
        # Analytical H2 estimate, then one Gibbs (fast, stable for recycle)
        co0 = float(syngas.imol['CO'])
        co2_0 = float(syngas.imol['CO2'])
        h2_0 = float(syngas.imol['H2'])
        h2_add = max(0.0, self.target_H2_CO * (co0 + 0.5 * co2_0) - h2_0)
        # Refine with a few bracket steps
        lo, hi = max(0.0, 0.5 * h2_add), max(h2_add * 2.0 + 100.0, 1000.0)
        best = h2_add
        for _ in range(8):
            mid = 0.5 * (lo + hi)
            h2.empty()
            h2.imol['H2'] = mid
            h2.phase = 'g'
            h2.T = 100 + 273.15
            h2.P = self.P_op
            out.mix_from([syngas, h2], energy_balance=False)
            out.phase = 'g'
            out.T = self.T_op
            out.P = self.P_op
            self._equilibrate(out)
            co = float(out.imol['CO'])
            ratio = float(out.imol['H2']) / co if co > 1e-12 else 0.0
            best = mid
            if ratio < self.target_H2_CO:
                lo = mid
            else:
                hi = mid
        h2.empty()
        h2.imol['H2'] = best
        h2.phase = 'g'
        h2.T = 100 + 273.15
        h2.P = self.P_op
        out.mix_from([syngas, h2], energy_balance=False)
        out.phase = 'g'
        out.T = self.T_op
        out.P = self.P_op
        self._equilibrate(out)

    def _design(self):
        self.design_results['Flow rate'] = sum(s.F_mass for s in self.ins)


@cost(basis='Flow rate', ID='Cyclone', units='kg/hr',
      cost=150000, S=50000, CE=522, n=0.6, BM=2.0)
class SyngasCyclone(bst.SolidsSeparator):
    _N_ins = 1
    _N_outs = 2
    _units = {'Flow rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), split=0.995):
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
        return [(1.0 - alpha) * (alpha ** (n - 1)) for n in range(n_min, n_max + 1)]

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
        H_reactants = syngas.H + injected_h2.H
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
        Q_released = H_reactants - effluent.H
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
