#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Mechanism probe for the high-glycolysis x low-ethanol-inhibition cliff of
the glycolysis x inhib_ethanol sweep (analyses/evaluate_EtOH_glycolysis_inhib_ethanol.py,
scenario A, enzyme burden ON, A's fed-batch strategy).

In the sweep, below a diagonal (inhib_ethanol ~ 0.5 x glycolysis) the ethanol
yield collapses and the IRR falls from +12.6 % to -4 % in one grid step (e.g.
glycolysis 2.564x: inhib 1.026x -> 0.923x -> 0.821x gives yield 0.475 ->
0.396 -> 0.324 g/g). This script re-simulates points across that cliff and a
set of counterfactuals at the collapsed point, and records for each:

* the fermentation outcomes (yield, titer, tau, spikes, cells, IRR);
* the node-by-node carbon routing at harvest from the cumulative reaction
  fluxes (nskinetics compute_flux_summary, integrated 0..tau): glucose ->
  glycolysis (r1) vs growth (r7); pyruvate -> Pdc (r3) vs TCA (r2) vs ALS
  (r13) vs accumulation; acetaldehyde -> Adh (r6) vs acetate (r4) vs
  accumulation; acetate -> TCA (r5) vs growth (r8) vs accumulation;
* trajectory diagnostics: peak pyruvate / acetaldehyde, the time fraction r3
  (Pdc) runs > 90 % saturated, and the share of r1 flux carried by the
  acetaldehyde-activated k_1e term;
* the fraction of r1 flux removed by ethanol (compute_flux_summary
  fraction_lost, an along-trajectory counterfactual, not a re-simulation);
* the enzyme-burden factor d (growth derating) at the point.

Hypotheses and their decisive counterfactuals (all at glycolysis 2.564x,
inhib_ethanol 0.821x):
  H1 Pdc (r3) bottleneck            -> 'k_3 x G' restores yield
  H2 Adh (r6) bottleneck            -> 'k_6 x G' restores yield
     + acetaldehyde-activated r1    -> 'k_1e unscaled' restores yield
  H3 the protective brake is ethanol inhibition OF GLYCOLYSIS
                                    -> 'k_1ie at 1x' restores, 'k_7ie at 1x'
                                       (growth-inhibition control) does not
Raising k_3 / k_6 enlarges the modeled proteome pool, which derates growth
more (a yield-favouring confound), so those counterfactuals and the cliff
point are also run with the burden OFF.

Writes analyses/results/glycolysis_cliff_probe_<stamp>.csv (one row per case)
and ..._trajectories.png (diagnostic). One load; ~20 simulations.
"""

import os
import traceback
from datetime import datetime

import numpy as np
import pandas as pd

from biorefineries import isobutanol
isobutanol.load()

from biorefineries.isobutanol import scenarios
from biorefineries.isobutanol import kinetic_optimization as ko
from biorefineries.isobutanol import system as isys

from nskinetics import compute_flux_summary
from nskinetics.engine import flux_analysis as fa
from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import FLUX_MAP_SPEC

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

#%% Set-up (as in the sweep)

model = isobutanol.models.models_EtOH_IBO_corn.model
fbs_spec = isobutanol.models.models_EtOH_IBO_corn.fbs_spec
solve_TEA = isys.solve_TEA
model_specification = model.specification
f = model.system.flowsheet
V406 = f.V406
km = V406.nsk_kinetic_model
r = km._te
product = f.ethanol
IBO_product = f.isobutanol

bundle = scenarios.load_scenario('A', burden=True)
burden_model = bundle['burden_model']

_available = ko.discover_kinetic_parameters(r)
GLYCOLYSIS = [m for m in ko.METABOLIC_14D_RATE_GROUPS['glycolysis'] if m in _available]
INHIB_ETHANOL = [m for m in ko.METABOLIC_MINIMAL_SUBSET_GROUPS['inhib_ethanol'] if m in _available]
EXTRA = ['k_3', 'k_6']
TOUCHED = GLYCOLYSIS + INHIB_ETHANOL + EXTRA
BASELINE = {m: float(getattr(r, m)) for m in TOUCHED}
print('glycolysis members:', {m: BASELINE[m] for m in GLYCOLYSIS})
print('inhib_ethanol members:', {m: BASELINE[m] for m in INHIB_ETHANOL})
print('k_3, k_6:', BASELINE['k_3'], BASELINE['k_6'])

# the sweep's grid values
X = np.linspace(1e-3, 4.0, 40)
Y = np.linspace(1e-3, 2.0, 40)
G = X[25]  # 2.5643
I_CLIFF = Y[16]  # 0.8207

#%% Cases: (name, glycolysis mult, inhib mult, {param: multiplier override}, burden_on)
# Overrides are multipliers of the scenario-A baseline, applied AFTER the
# group multipliers (so 'k_1ie': 1.0 restores that member alone).
CASES = [
    ('baseline',                   1.0,   1.0,    {}, True),
    ('IRR optimum',                X[13], Y[15],  {}, True),
    ('cliff col, inhib 1.026',     G,     Y[20],  {}, True),
    ('cliff col, inhib 0.923',     G,     Y[18],  {}, True),
    ('cliff col, inhib 0.821',     G,     I_CLIFF, {}, True),
    ('cliff col, inhib 0.513',     G,     Y[10],  {}, True),
    ('CF k_3 x G',                 G,     I_CLIFF, {'k_3': G}, True),
    ('CF k_6 x G',                 G,     I_CLIFF, {'k_6': G}, True),
    ('CF k_3,k_6 x G',             G,     I_CLIFF, {'k_3': G, 'k_6': G}, True),
    ('CF k_1e unscaled',           G,     I_CLIFF, {'k_1e': 1.0}, True),
    ('CF k_1ie at 1x',             G,     I_CLIFF, {'k_1ie': 1.0}, True),
    ('CF k_7ie at 1x (control)',   G,     I_CLIFF, {'k_7ie': 1.0}, True),
    ('CF k_10ie at 1x (control)',  G,     I_CLIFF, {'k_10ie': 1.0}, True),
    ('[burden off] cliff 0.821',   G,     I_CLIFF, {}, False),
    ('[burden off] CF k_3 x G',    G,     I_CLIFF, {'k_3': G}, False),
    ('[burden off] CF k_6 x G',    G,     I_CLIFF, {'k_6': G}, False),
    ('[burden off] CF k_3,k_6 x G', G,    I_CLIFF, {'k_3': G, 'k_6': G}, False),
    ('[burden off] cliff 1.026',   G,     Y[20],  {}, False),
]


def apply_case(g, i, overrides):
    for m in TOUCHED:
        setattr(r, m, BASELINE[m])
    for m in GLYCOLYSIS:
        setattr(r, m, BASELINE[m]*g)
    for m in INHIB_ETHANOL:
        setattr(r, m, BASELINE[m]*i)
    for m, mult in overrides.items():
        if m not in BASELINE:
            raise KeyError(f'override {m} not a touched parameter')
        setattr(r, m, BASELINE[m]*mult)


#%% Diagnostics of the run just simulated

def _col(df, name):
    for c in (f'[{name}]', name):
        if c in df.columns:
            return df[c].to_numpy()
    raise KeyError(name)


def stoich(species, reaction):
    m = r.getFullStoichiometryMatrix()
    rows = list(getattr(m, 'rownames', None) or r.getFloatingSpeciesIds())
    cols = list(getattr(m, 'colnames', None) or r.getReactionIds())
    return float(np.asarray(m)[rows.index(species), cols.index(reaction)])


def trajectory_rates(reactions):
    """Rates of `reactions` [g/h, extensive] along the harvested trajectory
    (rows 0..tau), replayed exactly as compute_flux_summary does; the model
    state is restored afterwards."""
    df = km.results_df
    df = df.iloc[:fa._end_index(df, V406.tau, V406) + 1]
    ordered = fa._write_order(km.state_selections())
    write_cols = [c for c in ordered if c != 'time']
    arrs = {c: df[c].to_numpy() for c in write_cols}
    rxn_ids = list(r.getReactionIds())
    idx_of = {rid: rxn_ids.index(rid) for rid in reactions}
    snap = {c: r[c] for c in ordered}
    try:
        rates = fa._rates_along(r, arrs, write_cols, idx_of, len(df))
    finally:
        fa._restore(r, snap, ordered)
    return df, rates


def diagnostics():
    reactions = ['r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r13']
    s = compute_flux_summary(V406, FLUX_MAP_SPEC.inhibition_map,
                             reactions=FLUX_MAP_SPEC.reactions)
    cm = s.cumulative_mass
    df, rates = trajectory_rates(['r1', 'r3', 'r6'])
    t = df['time'].to_numpy()
    env = _col(df, 'env') if ('env' in df.columns) else np.ones_like(t)
    pyr, ald = _col(df, 's_pyr'), _col(df, 's_acetald')
    ace = _col(df, 's_acetate')
    glu = _col(df, 's_glu')

    def accum(sp):
        a = _col(df, sp)
        return a[-1]*env[-1] - a[0]*env[0]

    # node balances [g of that node's species]
    glu_r1, glu_r7 = cm['r1'], cm['r7']
    pyr_in = stoich('s_pyr', 'r1')*cm['r1']
    ald_in = stoich('s_acetald', 'r3')*cm['r3']
    ace_in = stoich('s_acetate', 'r4')*cm['r4']
    out = {
        'glu -> glycolysis r1': glu_r1/(glu_r1 + glu_r7),
        'glu -> growth r7': glu_r7/(glu_r1 + glu_r7),
        'pyr -> Pdc r3': cm['r3']/pyr_in,
        'pyr -> TCA r2': cm['r2']/pyr_in,
        'pyr -> ALS r13': cm['r13']/pyr_in,
        'pyr accumulated': accum('s_pyr')/pyr_in,
        'ald -> Adh r6 (EtOH)': cm['r6']/ald_in,
        'ald -> acetate r4': cm['r4']/ald_in,
        'ald accumulated': accum('s_acetald')/ald_in,
        'ace -> TCA r5': (cm['r5']/ace_in) if ace_in > 0 else np.nan,
        'ace -> growth r8': (cm['r8']/ace_in) if ace_in > 0 else np.nan,
        'ace accumulated': (accum('s_acetate')/ace_in) if ace_in > 0 else np.nan,
        'cum r1 [g glu]': glu_r1, 'cum r7 [g glu]': glu_r7,
        'cum r3 [g pyr]': cm['r3'], 'cum r2 [g pyr]': cm['r2'],
        'cum r6 [g ald]': cm['r6'], 'cum r4 [g ald]': cm['r4'],
        'final pyr [g/L]': pyr[-1], 'peak pyr [g/L]': pyr.max(),
        'final ald [g/L]': ald[-1], 'peak ald [g/L]': ald.max(),
        'final acetate [g/L]': ace[-1], 'final glu [g/L]': glu[-1],
        'r1 lost to ethanol': s.fraction_lost.get('r1', {}).get('ethanol', np.nan),
        'r7 lost to ethanol': s.fraction_lost.get('r7', {}).get('ethanol', np.nan),
        'r4 lost to ethanol': s.fraction_lost.get('r4', {}).get('ethanol', np.nan),
        'r6 lost to ethanol': s.fraction_lost.get('r6', {}).get('ethanol', np.nan),
    }
    # Pdc saturation (r3's own Hill term) along the trajectory
    sat3 = pyr**4/(pyr**4 + r.K_3)
    dt = np.diff(t)
    out['frac time r3 > 90 % sat'] = float(np.sum(dt*(sat3[1:] > 0.9))/max(t[-1] - t[0], 1e-12))
    r1 = rates['r1']
    out['flux-wtd r3 saturation'] = float(np.trapezoid(rates['r3']*sat3, t)/max(np.trapezoid(rates['r3'], t), 1e-12))
    # share of r1 carried by the acetaldehyde-activated k_1e term (shared
    # inhibition exponentials cancel in the ratio)
    k1l, k1h, k1e = r.k_1l, r.k_1h, r.k_1e
    K1l, K1h, K1e, K1i = r.K_1l, r.K_1h, r.K_1e, r.K_1i
    term_l = k1l*glu/(glu + K1l) + k1h*glu/(glu + K1h)
    term_e = k1e*ald*glu/(glu*(K1i*ald + 1) + K1e)
    share_e = np.where(term_l + term_e > 0, term_e/(term_l + term_e), 0.)
    out['r1 share from k_1e term'] = float(np.trapezoid(r1*share_e, t)/max(np.trapezoid(r1, t), 1e-12))
    traj = {'t': t, 'pyr': pyr, 'ald': ald, 'ace': ace, 'glu': glu,
            'EtOH': _col(df, 's_EtOH'), 'x': _col(df, 'x'),
            'r1': r1/env, 'r3': rates['r3']/env, 'r6': rates['r6']/env}
    return out, traj


#%% Run

rows, trajs = [], {}
for name, g, i, overrides, burden_on in CASES:
    print(f'\n=== {name}: glycolysis {g:.4f}x, inhib_ethanol {i:.4f}x, '
          f'overrides {overrides}, burden {"ON" if burden_on else "OFF"}')
    row = {'case': name, 'glycolysis': g, 'inhib_ethanol': i,
           'overrides': str(overrides), 'burden': burden_on}
    try:
        apply_case(g, i, overrides)
        values = {n: float(getattr(r, n)) for n in burden_model.required_capacities()}
        b = burden_model.evaluate(values)
        row.update({'Phi_M': b.Phi_M, 'burden factor d': b.burden_factor if burden_on else 1.0})
        isys.set_active_burden(burden_model if burden_on else None)
        model_specification(**fbs_spec.current_specifications)
        sol = solve_TEA(stream_IDs=(product.ID, IBO_product.ID))
        d = V406.nsk_results_specific_tau_dict
        row.update({'IRR': sol['IRR'], 'EtOH MPSP': sol['MPSPs'][product.ID],
                    'yield': d['y_EtOH_IBO_glu_added'], 'titer': d['[s_EtOH]'],
                    'productivity': d['prod_EtOH'], 'tau': V406.tau,
                    'spikes': d['curr_n_glu_spikes'], 'cells': d['[x]'],
                    'TCI': model.system.TEA.TCI/1e6})
        diag, traj = diagnostics()
        row.update(diag)
        trajs[name] = traj
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})
    except Exception as e:
        row['error'] = f'{type(e).__name__}: {e}'
        print('ERROR:', row['error'])
        traceback.print_exc()
    rows.append(row)

# leave the model at the scenario-A baseline with the burden on
apply_case(1.0, 1.0, {})
isys.set_active_burden(burden_model)

#%% Save

stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
csv_path = os.path.join(out_dir, f'glycolysis_cliff_probe_{stamp}.csv')
pd.DataFrame(rows).to_csv(csv_path, index=False)
print('\nWrote', csv_path)

# Sweep reproduction check (kinetic outcomes do not depend on the hensmith
# version; IRR does, by ~+0.2 pp since the sweep ran)
P = os.path.join(out_dir, 'ibo_(40, 40, 1)_glyco_inhib_Spike_rb0.001-4_ib0.001-2_opt=False_max_n=16__')
try:
    sweep = {k: pd.read_csv(P + k + '.csv').iloc[:, 1:].to_numpy(float)
             for k in ('Combined Yield', 'EtOH Titer', 'Fermentation time', 'IRR')}
    print('\nSweep reproduction (probe vs sweep CSV):')
    for row in rows:
        if row['overrides'] != '{}' or not row['burden'] or 'yield' not in row:
            continue
        jj = int(np.argmin(abs(X - row['glycolysis'])))
        ii = int(np.argmin(abs(Y - row['inhib_ethanol'])))
        if abs(X[jj] - row['glycolysis']) > 1e-9 or abs(Y[ii] - row['inhib_ethanol']) > 1e-9:
            continue
        print(f"  {row['case']:<28} yield {row['yield']:.4f} vs {sweep['Combined Yield'][ii, jj]:.4f} | "
              f"titer {row['titer']:.2f} vs {sweep['EtOH Titer'][ii, jj]:.2f} | "
              f"tau {row['tau']:.2f} vs {sweep['Fermentation time'][ii, jj]:.2f} | "
              f"IRR {row['IRR']:.4f} vs {sweep['IRR'][ii, jj]:.4f}")
except FileNotFoundError:
    print('sweep CSVs not found; reproduction check skipped')

# Diagnostic trajectories
plot_cases = [c for c in ('cliff col, inhib 1.026', 'cliff col, inhib 0.923',
                          'cliff col, inhib 0.821', 'CF k_3 x G', 'CF k_6 x G',
                          'CF k_1e unscaled', 'CF k_1ie at 1x',
                          'CF k_7ie at 1x (control)') if c in trajs]
panels = [('pyr', 'pyruvate [g/L]'), ('ald', 'acetaldehyde [g/L]'),
          ('ace', 'acetate [g/L]'), ('EtOH', 'ethanol [g/L]'),
          ('x', 'cells [g/L]'), ('r1', 'r1 glycolysis [g glu/L/h]'),
          ('r3', 'r3 Pdc [g pyr/L/h]'), ('r6', 'r6 Adh [g ald/L/h]')]
fig, axes = plt.subplots(2, 4, figsize=(18, 8))
for ax, (key, title) in zip(axes.flat, panels):
    for c in plot_cases:
        ax.plot(trajs[c]['t'], trajs[c][key], label=c)
    ax.set_title(title)
    ax.set_xlabel('time [h]')
axes.flat[0].legend(fontsize=7)
fig.tight_layout()
png_path = csv_path.replace('.csv', '_trajectories.png')
fig.savefig(png_path, dpi=150)
print('Wrote', png_path)
