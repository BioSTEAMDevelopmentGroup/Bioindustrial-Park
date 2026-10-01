#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Main figure of the kinetic-BO campaigns, "Scouts explore, profit
selects" -> analyses/results/publication/Optimization-figures/
kinBO_main_scouts_<stamp>.{png,pdf} (+ _latest copies) and the caption draft
kinBO_main_scouts_caption.md next to this script.

Panels (figwork/figure_spec.md section 3; take-aways T1-T5):
    a  when: IRR of the best design so far vs simulated trials (log x),
       uninformed vs TRY-informed, with the 1,000 seeds as a strip on the
       left and the refinement inset (T2, T4)
    b  where: IRR of every simulated design vs isobutanol share of the
       alcohol titer (T2, T3, T4)
    c  visited vs returned: every trial of each campaign as a strip, the
       best visit (open) and the returned design (filled), plus a two-column
       table (trials above the plateau, trials losing) (T3)
    d  titers of each returned design (T1, T4)
    e  metabolic proteome of each returned design (T1)
    f  exploration (trials making isobutanol) vs IRR of the best visit and
       the returned design (T5)

SIM-SAFE: pandas / numpy / matplotlib only, through the sibling modules
_common (data, definitions, facts) and _style (palette, axes, checks);
enzyme_burden.py is read BY FILE PATH for the caption's housekeeping sector.
Every annotation number is built from `_common.compute_facts()` (asserted
against `_common.EXPECTED` by `check_facts()` before drawing); the script also
asserts the figure-level claims it prints (e.g. "ethanol only", "all
co-production"), counts the loss dots against the facts, and checks that no
annotation touches a key marker, connector, trajectory line, arrow or the
inset (`mark_clashes`), on top of `_style.check_figure`.

Run:  "C:/Users/saran/anaconda3/envs/IBO_2026/python.exe" fig_main_scouts.py
      [--out-dir DIR] [--no-latest] [--no-caption]
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C                                          # noqa: E402
import _style as S                                           # noqa: E402

from matplotlib import transforms as mtransforms             # noqa: E402
from matplotlib.collections import PathCollection            # noqa: E402
from matplotlib.lines import Line2D                          # noqa: E402
from matplotlib.patches import (ConnectionPatch,             # noqa: E402
                                FancyArrowPatch)
from matplotlib.text import Text                             # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator   # noqa: E402

STEM = 'kinBO_main_scouts'
CAPTION_PATH = os.path.join(C.FIGDIR, 'kinBO_main_scouts_caption.md')
EB_PATH = os.path.join(C.PKG, 'enzyme_burden.py')

# cross-process reproducibility margin of the (deterministic) PI objective:
# the relay campaign's reproduction of its best trial (project memory
# relay-campaign-12d); differences below it are never called real
REPRO_MARGIN_PI = 0.0105
# scenario-A IRR under the CURRENT model version (CLAUDE.md "Current
# baseline", hensmith 6d4776f); the campaigns ran under the 2b5b27d export
START_IRR_CURRENT_MODEL = 0.1279
SEED = 350                     # the campaigns' optuna seed (_rs350)
HURDLE_PCT = 15                # PI = NPV at this hurdle rate / TCI

TEXT, NOTE = S.TEXT, S.NOTE
P = S.PALETTE
GREY_ARROW = '0.35'
GREY_LABEL = '#555555'
GREY_FOOT = '#888888'
GREY_TAG = '#444444'
CONNECTOR = dict(color='0.25', lw=1.1, solid_capstyle='butt')

# %% Layout (inches, origin bottom-left; spec 3.1) ------------------------------
L = {
    'a_strip': (0.72, 4.55, 0.40, 2.85),
    'a': (1.17, 4.55, 4.20, 2.85),
    'a_inset': (0.585, 0.075, 0.385, 0.33),     # axes fraction of a
    'b': (6.02, 4.55, 3.86, 2.85),
    'c': (1.52, 0.72, 2.02, 2.58),
    'd': (4.44, 0.72, 0.84, 2.58),
    'e': (5.40, 0.72, 1.50, 2.58),
    'f': (7.52, 0.72, 2.36, 2.58),
}
TABLE_X = (3.75, 4.23)            # centres of the two table columns
GROUP_X = 0.12                    # group headers, left-aligned
ROW_LABEL_X = 1.30                # row labels, right-aligned
ROW_MARK_X = 1.40                 # role marker after the row label
SEP_X = (0.12, 6.90)              # separators: labels .. e
TOP_LETTER_Y, BOT_LETTER_Y = 7.66, 3.86
LETTER_X = {'a': 0.06, 'b': 5.42, 'c': 0.06, 'd': 4.30, 'e': 5.30, 'f': 7.06}
BAND_Y = (3.62, 3.44)             # the two key lines of the header band
KEY_DE_X = 4.44                   # d/e key start (= d's left edge)
KEY_GAP_IN = 0.22                 # min gap between the d/e and f keys

# bottom-row rows (shared y of c, the table, d, e; y grows downward)
ROW_Y = {'base': 0.0, 'unin': 1.8, 'relay': 2.8, 'ey': 4.6, 'et': 5.6,
         'ep': 6.6, 'iy': 8.4, 'it': 9.4, 'ip': 10.4}
HEADER_Y = {'profit': 1.0, 'etoh': 3.8, 'ibo': 7.6}
ROW_YLIM = (10.9, -0.5)
BAR_H = 0.62

# role-marker sizes (scatter s, pt^2)
S_RING, S_FILLED = 44, 38
S_DIAMOND_C, S_STAR_C, S_BASE_C = 60, 120, 40

ANNOT = []        # annotation Text artists checked against the marks
MARKS = []        # key marks / lines / arrows annotations must not touch


# %% Small helpers ---------------------------------------------------------------
def _ann(ax, s, xy, xytext=(0, 0), **kw):
    """Annotation anchored at data `xy`, offset `xytext` points; registered
    for the mark-clash check."""
    kw.setdefault('fontsize', S.FS['annot'])
    kw.setdefault('color', TEXT)
    t = ax.annotate(s, xy=xy, xytext=xytext, textcoords='offset points',
                    annotation_clip=False, **kw)
    ANNOT.append(t)
    return t


def _mark(ax, x, y, marker, s, face, edge, lw=0.0, z=5, check=True):
    sc = ax.scatter([x], [y], s=s, marker=marker, facecolors=face,
                    edgecolors=edge, linewidths=lw, zorder=z, clip_on=False)
    if check:
        MARKS.append(sc)
    return sc


def _ring(ax, x, y, marker, edge, s=S_RING, lw=1.4, z=5):
    return _mark(ax, x, y, marker, s, 'white', edge, lw, z)


def _filled(ax, x, y, marker, color, s=S_FILLED, z=5.5, edge='none',
             lw=0.0):
    return _mark(ax, x, y, marker, s, color, edge, lw, z)


def _step_xy(sims, vals, x_end):
    """Vertices of a post-step line through (sims, vals), extended to x_end
    (drawn with plot so the vertices ARE the step path)."""
    x = np.repeat(np.asarray(sims, float), 2)[1:]
    y = np.repeat(np.asarray(vals, float), 2)[:-1]
    return np.r_[x, x_end], np.r_[y, y[-1]]


def _text_w_in(fig, s, fontsize, **kw):
    t = Text(0, 0, s, fontsize=fontsize, **kw)
    t.set_figure(fig)
    return t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi


def _rows_transform(fig, ax):
    """x in figure inches, y in the bottom-row data coordinates."""
    return mtransforms.blended_transform_factory(fig.dpi_scale_trans,
                                                 ax.transData)


def _irr(v_pct):
    """IRR in % (facts) -> plot coordinate (losses at the band centre)."""
    return S.irr_plot(float(v_pct) / 100.0)


def _pct0(v):
    return f'{float(v):.0f} %'


def _g(v, nd=1):
    return C.fmt_num(v, nd)


WORDS = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six',
         7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten'}


def _ordinal(n):
    n = int(n)
    suf = 'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(
        n % 10, 'th')
    return f'{n}{suf}'


# %% Claims the figure text makes (asserted, not assumed) ---------------------------
def assert_claims(F):
    camp = F['campaigns']
    U = F['U_pct']
    rel = F['relay']
    # T2: the uninformed plateau is an ethanol-only design, flat to the end
    assert camp['unin']['best_visit']['class'] == 'etoh'
    assert camp['unin']['best_visit']['ibo'] < C.IBO_THRESHOLD
    assert F['unin']['plateau_sim'] == F['U_sim']
    assert F['unin']['bsf_pct'][2000] == F['unin']['bsf_pct'][F['U_sim']]
    assert F['unin']['ibo_ge5']['max_irr_pct'] < U          # all worse
    assert (F['unin']['ibo_ge5']['pct_losing']
            > F['unin']['ibo_lt5']['pct_losing'])          # and riskier
    # b title: every design above the plateau makes isobutanol
    assert F['global']['n_gt_U_ibo_lt5'] == 0
    assert F['global']['n_alcohol_lt1_with_irr'] == 0
    # 296 relay designs above the best seed, all co-production, no other
    assert F['global']['n_nonrelay_gt_best_seed'] == 0
    assert (F['global']['n_relay_gt_best_seed']
            == F['global']['n_relay_gt_best_seed_coprod'])
    # seeds above U all from isobutanol scouts
    assert F['seeds']['n_gt_U'] == F['seeds']['n_gt_U_ibo']
    assert F['seeds']['max_etoh_irr_pct'] < U
    # a: first trial lost money, trial 2 co-production above U, trial 13
    # first above every seed
    assert C.is_loss(rel['sims'][1]['irr_pct'])
    assert rel['first_gt_U'] == 2 and rel['sims'][2]['class'] == 'coprod'
    assert rel['sims'][2]['irr_pct'] == rel['sims'][2]['bsf_pct']
    assert rel['first_gt_best_seed'] in rel['sims']
    assert rel['sims'][25]['irr_pct'] == rel['bsf_pct'][25]
    # c title: every scout's best visit is above its returned design
    for k in C.SCOUT_KEYS:
        st = camp[k]
        assert st['best_visit']['irr_pct'] > st['returned']['irr_pct'], k
    # ethanol scouts stayed below the plateau
    for k in C.ETOH_SCOUTS:
        assert camp[k]['n_gt_U'] == 0
        assert camp[k]['best_visit']['irr_pct'] < U
    # T5 synthesis: no non-relay campaign returned a co-production design,
    # the relay's best is the highest IRR of the eight
    for k in ('unin',) + C.SCOUT_KEYS:
        assert camp[k]['returned']['class'] != 'coprod', k
    assert rel['best']['class'] == 'coprod'
    assert all(rel['best']['irr_pct'] >= camp[k]['best_visit']['irr_pct']
               for k in C.MAIN_KEYS)
    # the profitability campaigns: returned == best visit
    for k in C.PROFIT_KEYS:
        assert (camp[k]['returned']['trial']
                == camp[k]['best_visit']['trial']), k
    # isobutanol scouts' best: pure isobutanol (yield scout)
    assert camp['iy']['best_visit']['share_pct'] > 99.9
    assert (camp['iy']['best_visit']['irr_pct']
            == max(camp[k]['best_visit']['irr_pct'] for k in C.IBO_SCOUTS))
    assert abs(camp['iy']['best_visit']['irr_pct']
               - F['seeds']['best_irr_pct']) < 1e-9
    # proteome: idle dehydrogenases named in the caption
    pr = F['proteome']
    assert pr['et']['ehr'] < 1e-3 < pr['et']['adh6']
    assert pr['iy']['pdc'] < 1e-3 < pr['iy']['adh1']
    assert pr['relay']['Phi_M'] <= F['budget'] and pr['relay']['growth'] > .99


# %% Panel a: when ------------------------------------------------------------------
def panel_a(fig, F):
    U, start = F['U_pct'], F['start_irr_pct']
    rel, unin = F['relay'], F['unin']
    ax = S.inch_axes(fig, *L['a'])
    strip = S.inch_axes(fig, *L['a_strip'], sharey=ax)

    # --- seed strip (shares y with a; carries the y title and labels)
    S.irr_axis(strip, 'y', step=5)
    S.plateau_line(strip, U)
    strip.set_xlim(0, 1)
    man = C.load_manifest('relay')
    rng = S.panel_rng()
    best_idx = man['IRR'].idxmax()
    n_loss = 0
    for fam, (lo, hi), col in (('etoh', (0.08, 0.42), P['etoh_light']),
                               ('ibo', (0.58, 0.92), P['ibo_light'])):
        m = man[man['family'] == fam]
        x = rng.uniform(lo, hi, len(m))
        y = S.irr_plot(m['IRR'].to_numpy(float), rng)
        x[m.index.to_numpy() == best_idx] = 0.75      # ring sits on its dot
        sc = strip.scatter(x, y, s=3, c=col, alpha=0.55, lw=0,
                           rasterized=True, zorder=2)
        n_loss += int((sc.get_offsets()[:, 1] < 0).sum())
    assert n_loss == F['seeds']['n_losing'], 'seed loss dots != facts'
    assert man.loc[best_idx, 'family'] == 'ibo'
    _ring(strip, 0.75, F['seeds']['best_irr_pct'], 'o', P['ibo_dark'], s=40,
          lw=1.2)
    S.style_ticks(strip, minor_x=False)
    strip.set_xticks([0.25, 0.75])
    strip.set_xticklabels(['EtOH', 'IBO'], fontsize=S.FS['key'],
                          rotation=90)
    for lab, col in zip(strip.get_xticklabels(),
                        (P['etoh_dark'], P['ibo_dark'])):
        lab.set_color(col)
    strip.tick_params(axis='x', which='both', length=0, top=False, pad=3)
    strip.set_ylabel(S.bold_axis_title('IRR [%]'), labelpad=2)
    S.fig_text(fig, L['a_strip'][0] + L['a_strip'][2] / 2,
               L['a_strip'][1] + L['a_strip'][3] + 0.05,
               f"{C.fmt_int(F['seeds']['n'])} seeds", fontsize=S.FS['key'],
               fontweight='bold', ha='center', va='bottom')

    # --- main axes
    S.log_trial_axis(ax)
    S.loss_band(ax)
    n_start = F['checks']['shared_startup_rows']
    ax.axvspan(1, n_start + 1, color=P['startup_band'], lw=0, zorder=0.05)
    MARKS.append(S.plateau_line(ax, U))
    MARKS.append(S.start_line(ax, start))
    lines = {}
    for key, z in (('unin', 3.0), ('relay', 3.5)):          # TRY-informed last
        s, b = C.best_so_far(C.complete(key))
        x, y = _step_xy(s, S.irr_plot(b), F['campaigns'][key]['max_sim'])
        col = P['relay'] if key == 'relay' else P['unin']
        ln, = ax.plot(x, y, color=col, lw=2.6, zorder=z,
                      solid_joinstyle='miter', solid_capstyle='butt')
        lines[key] = (s, b)
        MARKS.append(ln)
    # spec 6.B.7: cyan flat from the plateau trial to 2,000; teal >= the
    # trial-2 value from trial 2 on, a loss at trial 1
    s_u, b_u = lines['unin']
    assert np.all(b_u[s_u >= unin['plateau_sim']] == b_u[-1])
    s_r, b_r = lines['relay']
    assert np.isneginf(b_r[0]) or b_r[0] < 0
    assert np.all(b_r[s_r >= 2] >= rel['sims'][2]['bsf_pct'] / 100 - 1e-12)

    # end markers and milestones
    rb = rel['best']
    _mark(ax, F['U_sim'], U, 'D', 55, P['unin'], 'white', 0.8, z=6)
    for sim in (rel['first_gt_U'], rel['first_ge25'], rb['sim']):
        v = rel['sims'][sim]['bsf_pct'] if sim in rel['sims'] else \
            rel['bsf_pct'][sim]
        _mark(ax, sim, v, 'o', 22, P['relay'], 'white', 0.6, z=6)
    _mark(ax, rb['sim'], rb['irr_pct'], '*', 140, P['relay'], 'white', 0.8,
          z=7)
    v2 = rel['sims'][rel['first_gt_U']]['bsf_pct']
    v18 = rel['sims'][rel['first_ge25']]['bsf_pct']
    fs = S.FS['note']
    _ann(ax, C.fmt_irr(v2), (rel['first_gt_U'], v2), (5, -3), ha='left',
         va='top', fontsize=fs, color=P['relay'])
    _ann(ax, f"{C.fmt_irr(v18)} · {C.fmt_trial(rel['first_ge25'])}",
         (rel['first_ge25'], v18), (5, -4), ha='left', va='top',
         fontsize=fs, color=P['relay'])
    _ann(ax, C.fmt_irr(rb['irr_pct']), (rb['sim'], rb['irr_pct']), (0, 7),
         ha='center', va='bottom', fontsize=fs, color=P['relay'])

    # start-up band foot label, starting strain, direct labels
    _ann(ax, 'space-filling start-up', (np.sqrt(n_start + 1), 1.0), ha='center',
         va='center', fontsize=fs, color=GREY_FOOT)
    _ann(ax, f'starting strain {C.fmt_irr(start)}', (2000, start), (0, -3),
         ha='right', va='top', fontsize=fs, color=NOTE)
    _ann(ax, 'TRY-informed', (1.5, 28.6), ha='left', va='bottom',
         fontweight='bold', color=P['relay'])
    _ann(ax, f"uninformed: {C.fmt_irr(U)} at {C.fmt_trial(unin['plateau_sim'])},"
             f"\nethanol only; flat for {C.fmt_int(unin['flat_trials'])} trials",
         (2000, U), (0, 7), ha='right', va='bottom', color=P['unin_text'],
         linespacing=1.15)

    # seeds -> GP arrow
    best_seed = F['seeds']['best_irr_pct']
    con = ConnectionPatch(
        xyA=(1.0, best_seed), coordsA=strip.get_yaxis_transform(),
        xyB=(1.85, 20.6), coordsB=ax.transData, arrowstyle='-|>',
        connectionstyle='arc3,rad=-0.25', color=GREY_ARROW, lw=0.9,
        mutation_scale=8, shrinkA=4, shrinkB=1, zorder=8)
    fig.add_artist(con)
    MARKS.append(con)
    t = S.fig_text(fig, L['a'][0] + 0.08, 6.86, 'preloaded into the GP',
                   fontsize=fs, color=GREY_LABEL, ha='left', va='bottom')
    ANNOT.append(t)

    # inset: refinement of the TRY-informed best so far
    ins = ax.inset_axes(L['a_inset'])
    ins.set_facecolor('white')
    for sp in ins.spines.values():
        sp.set_color('0.6')
        sp.set_linewidth(0.6)
    x, y = _step_xy(s_r, 100 * b_r, F['campaigns']['relay']['max_sim'])
    ins.plot(x, y, color=P['relay'], lw=1.8, solid_joinstyle='miter')
    ins_sims = (25, 100, rel['first_ge27'], rb['sim'])
    for sim in ins_sims:
        ins.scatter([sim], [rel['bsf_pct'][sim]], s=16, color=P['relay'],
                    edgecolors='white', linewidths=0.5, zorder=5)
    ins.set_xlim(20, 1000)
    ins.set_ylim(25.6, 27.5)
    ins.xaxis.set_major_locator(FixedLocator([25, 500, 1000]))
    ins.xaxis.set_major_formatter(FixedFormatter(['25', '500', '1,000']))
    ins.yaxis.set_major_locator(FixedLocator([26, 27]))
    ins.yaxis.set_major_formatter(FixedFormatter(['26', '27']))
    ins.tick_params(labelsize=S.FS['note'], pad=1.5)
    S.style_ticks(ins)
    pi25 = rel['sims'][25]['PI']
    pi_best = rel['bsf_PI'][rb['sim']]
    t = ins.text(0.04, 0.93, f'refined: PI {pi25:.2f} → {pi_best:.2f}',
                 transform=ins.transAxes, ha='left', va='top',
                 fontsize=S.FS['note'], color=P['relay'])
    S.style_ticks(ax)
    ax.tick_params(axis='y', which='both', labelleft=False)
    return ax, strip, ins


# %% Panel b: where -----------------------------------------------------------------
def panel_b(fig, F):
    U, start = F['U_pct'], F['start_irr_pct']
    camp = F['campaigns']
    ax = S.inch_axes(fig, *L['b'])
    S.irr_axis(ax, 'y', step=10)
    ax.set_xlim(-3, 103)
    ax.xaxis.set_major_locator(FixedLocator([0, 25, 50, 75, 100]))
    ax.set_xlabel(S.bold_axis_title('Isobutanol share of alcohol titer [%]'))
    S.tint_above(ax, U)
    MARKS.append(S.plateau_line(ax, U))
    MARKS.append(S.start_line(ax, start))
    rng = S.panel_rng()
    n_drawn = n_loss = n_loss_exp = 0
    for keys, col, s, a in ((C.ETOH_SCOUTS, P['etoh_light'], 3, 0.45),
                            (C.IBO_SCOUTS, P['ibo_light'], 3, 0.45),
                            (('relay',), P['relay'], 5, 0.55),
                            (('unin',), P['unin'], 5, 0.7)):
        d = pd.concat([C.complete(k) for k in keys], ignore_index=True)
        d = d[d['alcohol'] >= C.SHARE_MIN_ALCOHOL]
        y = S.irr_plot(d['IRR'].to_numpy(float), rng)
        sc = ax.scatter(d['ibo_share_pct'], y, s=s, c=col, alpha=a, lw=0,
                        rasterized=True, zorder=2)
        n_drawn += len(d)
        n_loss += int((sc.get_offsets()[:, 1] < 0).sum())
        n_loss_exp += int(d['loss'].sum())
    assert n_drawn == F['global']['n_alcohol_ge1'], n_drawn
    assert n_loss == n_loss_exp

    # marks: isobutanol scouts' best visits, the two profitability bests
    for k in C.IBO_SCOUTS:
        bv = camp[k]['best_visit']
        _ring(ax, bv['share_pct'], bv['irr_pct'], C.campaign(k).marker,
              P['ibo_dark'], s=42, lw=1.3, z=6)
    rb = F['relay']['best']
    ub = camp['unin']['best_visit']
    _mark(ax, rb['share_pct'], rb['irr_pct'], '*', 140, P['relay'], 'white',
          0.8, z=7)
    _mark(ax, ub['share_pct'], ub['irr_pct'], 'D', 50, P['unin'], 'white',
          0.6, z=7)

    box = dict(boxstyle='square,pad=0.15', fc='white', ec='none', alpha=0.8)
    _ann(ax, f'uninformed plateau {C.fmt_irr(U)}', (3, U), (0, 2.5),
         ha='left', va='bottom', color=P['unin_text'], bbox=box)
    ui = F['unin']
    _ann(ax, f"uninformed designs with ≥ {C.IBO_THRESHOLD:g} "
             f"{S.UNIT_TITER} isobutanol: {ui['ibo_ge5']['n']},\n"
             f"best {C.fmt_irr(ui['ibo_ge5']['max_irr_pct'])}; "
             f"{_pct0(ui['ibo_ge5']['pct_losing'])} lose money "
             f"({_pct0(ui['ibo_lt5']['pct_losing'])} below "
             f"{C.IBO_THRESHOLD:g} {S.UNIT_TITER})",
         (50, 6.5), ha='center', va='center', fontsize=S.FS['note'],
         color=P['unin_text'], bbox=box, linespacing=1.2)
    best_seed = F['seeds']['best_irr_pct']
    iy_bv = camp['iy']['best_visit']
    _ann(ax, f"isobutanol scouts' best: {C.fmt_irr(best_seed)}",
         (iy_bv['share_pct'], iy_bv['irr_pct']), (-8, 0), ha='right',
         va='center',
         color=P['ibo_dark'], bbox=box)
    g = F['global']
    _ann(ax, f"{g['n_relay_gt_best_seed']} designs above "
             f"{C.fmt_irr(best_seed)}: all TRY-informed co-production",
         (50, 31.0), (0, -4), ha='center', va='top', fontsize=9.5,
         color=P['relay'])
    _ann(ax, f"{C.fmt_irr(rb['irr_pct'])}: {rb['etoh']:.0f} {S.UNIT_TITER} "
             f"ethanol\n+ {rb['ibo']:.0f} {S.UNIT_TITER} isobutanol",
         (rb['share_pct'], rb['irr_pct']), (-14, -9), ha='right',
         va='center', fontsize=S.FS['note'], color=P['relay'], bbox=box,
         linespacing=1.15)
    S.style_ticks(ax)
    return ax


# %% Panels c, d, e: the bottom rows ----------------------------------------------------
def _row_axes(fig, key):
    ax = S.inch_axes(fig, *L[key])
    ax.set_ylim(*ROW_YLIM)
    return ax


def panel_c(fig, F):
    U, start = F['U_pct'], F['start_irr_pct']
    camp = F['campaigns']
    ax = _row_axes(fig, 'c')
    S.irr_axis(ax, 'x', step=10, loss_fs=S.FS['note'])
    S.tint_above(ax, U, which='x')
    S.plateau_line(ax, U, which='x')
    S.start_line(ax, start, which='x')
    rng = S.panel_rng()
    for k in C.MAIN_KEYS:
        c = C.campaign(k)
        d = C.complete(k)
        x = S.irr_plot(d['IRR'].to_numpy(float), rng)
        y = ROW_Y[k] + rng.uniform(-0.32, 0.32, len(d))
        if c.family == 'profit':
            col, a = (P['relay'] if c.is_relay else P['unin']), 0.35
        else:
            col, a = P[c.family + '_light'], 0.5
        sc = ax.scatter(x, y, s=2.5, c=col, alpha=a, lw=0, rasterized=True,
                        zorder=2)
        n_loss = int((sc.get_offsets()[:, 0] < 0).sum())
        exp = round(camp[k]['pct_losing'] * camp[k]['n_complete'] / 100)
        assert n_loss == exp, (k, n_loss, exp)       # spec 6.B.6
    # markers
    for k in C.SCOUT_KEYS:
        c = C.campaign(k)
        st = camp[k]
        dark = P[c.family + '_dark']
        xb, xr = _irr(st['best_visit']['irr_pct']), _irr(st['returned'][
            'irr_pct'])
        ln, = ax.plot([xr, xb], [ROW_Y[k]] * 2, zorder=4, **CONNECTOR)
        _ring(ax, xb, ROW_Y[k], c.marker, dark, s=S_RING, lw=1.4, z=5)
        _filled(ax, xr, ROW_Y[k], c.marker, dark, s=S_FILLED, z=5.5)
    _mark(ax, U, ROW_Y['unin'], 'D', S_DIAMOND_C, P['unin'], 'white', 0.7,
          z=6)
    _mark(ax, F['relay']['best']['irr_pct'], ROW_Y['relay'], '*', S_STAR_C,
          P['relay'], 'white', 0.8, z=6)
    _mark(ax, start, ROW_Y['base'], 'D', S_BASE_C, 'white', P['base'], 1.2,
          z=6)
    S.style_ticks(ax, y=False)
    return ax


def row_labels(fig, F, ax_c):
    tr = _rows_transform(fig, ax_c)
    fs = S.FS['row']
    fig.text(ROW_LABEL_X, ROW_Y['base'], 'Starting strain', transform=tr,
             ha='right', va='center', fontsize=fs, color=TEXT)
    for fam, y in HEADER_Y.items():
        col = {'profit': TEXT, 'etoh': P['etoh_dark'], 'ibo': P['ibo_dark']}[
            fam]
        fig.text(GROUP_X, y, C.FAMILY_LABEL[fam], transform=tr, ha='left',
                 va='center', fontsize=fs, fontweight='bold', color=col)
        sep_y = y - 0.45
        fig.add_artist(Line2D(SEP_X, [sep_y, sep_y], transform=tr,
                              color='0.85', lw=0.6, zorder=0.5))
    marks = [('base', 'D', 'white', P['base'], 6.0, 1.1)]
    for k in C.MAIN_KEYS:
        c = C.campaign(k)
        if c.family == 'profit':
            col = P['relay'] if c.is_relay else P['unin_text']
            fig.text(ROW_LABEL_X, ROW_Y[k], c.label, transform=tr, ha='right',
                     va='center', fontsize=fs, fontweight='bold', color=col)
            mcol = P['relay'] if c.is_relay else P['unin']
            ms = 9.0 if c.is_relay else 6.0
            marks.append((k, c.marker, mcol, 'white', ms, 0.5))
        else:
            fig.text(ROW_LABEL_X, ROW_Y[k], c.label, transform=tr, ha='right',
                     va='center', fontsize=fs, color=TEXT)
            marks.append((k, c.marker, 'white', P[c.family + '_dark'], 6.2,
                          1.2))
    for k, m, fc, ec, ms, mew in marks:
        fig.add_artist(Line2D([ROW_MARK_X], [ROW_Y[k]], transform=tr,
                              ls='none', marker=m, ms=ms, mfc=fc, mec=ec,
                              mew=mew))


def table_c(fig, F, ax_c):
    tr = _rows_transform(fig, ax_c)
    camp = F['campaigns']
    U = F['U_pct']
    hy = BAND_Y[1] + 0.06
    for x, h in zip(TABLE_X, (f'trials\n> {C.fmt_irr(U)}', 'trials\nlosing')):
        S.fig_text(fig, x, hy, h, fontsize=S.FS['note'], ha='center',
                   va='center', color=TEXT, linespacing=1.1)
    fs = S.FS['table']
    for k in C.MAIN_KEYS:
        st = camp[k]
        y = ROW_Y[k]
        if k == 'unin':
            n_txt = '–'
        else:
            n_txt = C.fmt_int(st['n_gt_U'])
        bold = k == 'iy'
        if C.campaign(k).is_relay:
            fig.text(TABLE_X[0], y - 0.13, n_txt, transform=tr, ha='center',
                     va='center', fontsize=fs, color=TEXT)
            fig.text(TABLE_X[0], y + 0.47,
                     f"of {C.fmt_int(st['n_complete'])}", transform=tr,
                     ha='center', va='center', fontsize=S.FS['note'],
                     color=NOTE, zorder=3,
                     bbox=dict(boxstyle='square,pad=0.1', fc='white',
                               ec='none'))
        else:
            fig.text(TABLE_X[0], y, n_txt, transform=tr, ha='center',
                     va='center', fontsize=fs, color=TEXT,
                     fontweight='bold' if bold else 'normal')
        fig.text(TABLE_X[1], y, _pct0(st['pct_losing']), transform=tr,
                 ha='center', va='center', fontsize=fs, color=TEXT)


def panel_d(fig, F):
    ax = _row_axes(fig, 'd')
    pr = F['proteome']
    for k in C.ROW_KEYS:
        y = ROW_Y[k]
        et, ib = pr[k]['etoh'], pr[k]['ibo']
        ax.barh(y, et, height=BAR_H, color=P['etoh'], ec='white', lw=0.4,
                zorder=2)
        ax.barh(y, ib, left=et, height=BAR_H, color=P['ibo'], ec='white',
                lw=0.4, zorder=2)
    ax.set_xlim(0, 175)
    ax.xaxis.set_major_locator(FixedLocator([0, 100]))
    ax.xaxis.set_minor_locator(FixedLocator([50, 150]))
    ax.set_xlabel(S.bold_axis_title(f'Titer\n[{S.UNIT_TITER}]'))
    S.style_ticks(ax, y=False)
    return ax


E_XMAX = 0.27                     # e's x limit (room for the x0.30 tag)


def panel_e(fig, F):
    ax = _row_axes(fig, 'e')
    pr = F['proteome']
    for k in C.ROW_KEYS:
        y = ROW_Y[k]
        left = 0.0
        for sec in C.SECTORS:
            v = pr[k][sec.key]
            ax.barh(y, v, left=left, height=BAR_H, color=P[sec.color_key],
                    ec='white', lw=0.4, zorder=2)
            left += v
        assert abs(left - pr[k]['Phi_M']) < 1e-6, k
        g = pr[k]['growth']
        if g < 0.99:
            _ann(ax, f'×{g:.2f}', (left + 0.004, y), ha='left', va='center',
                 fontsize=S.FS['note'], fontstyle='italic', color=GREY_TAG)
    bud = F['budget']
    ax.axvline(bud, color='k', ls=':', lw=0.9, zorder=3)
    assert pr['base']['Phi_M'] < bud
    _ann(ax, 'penalty →', (bud + 0.0033, ROW_Y['base']), ha='left',
         va='center', fontsize=S.FS['note'], color=TEXT)
    ax.set_xlim(0, E_XMAX)
    ax.xaxis.set_major_locator(FixedLocator([0, 0.1, 0.2]))
    ax.xaxis.set_major_formatter(FixedFormatter(['0', '0.1', '0.2']))
    ax.set_xlabel(S.bold_axis_title(
        f'Metabolic proteome\n[{S.UNIT_PROTEOME}]'))
    S.style_ticks(ax, y=False)
    return ax


# %% Panel f: exploration vs selection --------------------------------------------------
L_F_XLIM = (0.008, 0.93)          # f's logit x limits (fractions)


def panel_f(fig, F):
    U, start = F['U_pct'], F['start_irr_pct']
    camp = F['campaigns']
    ax = S.inch_axes(fig, *L['f'])
    S.logit_pct_axis(ax, lim_pct=(100 * L_F_XLIM[0], 100 * L_F_XLIM[1]))
    S.irr_axis(ax, 'y', step=10)
    ax.yaxis.labelpad = 1.5
    S.tint_above(ax, U)
    MARKS.append(S.plateau_line(ax, U))
    MARKS.append(S.start_line(ax, start))
    pos = {}
    for k in C.SCOUT_KEYS:
        c = C.campaign(k)
        st = camp[k]
        x = st['exploration_pct'] / 100
        yb, yr = _irr(st['best_visit']['irr_pct']), _irr(st['returned'][
            'irr_pct'])
        dark = P[c.family + '_dark']
        ln, = ax.plot([x, x], [yr, yb], zorder=4, **CONNECTOR)
        MARKS.append(ln)
        _ring(ax, x, yb, c.marker, dark, s=S_RING, lw=1.4, z=5)
        _filled(ax, x, yr, c.marker, dark, s=S_FILLED, z=5.5)
        pos[k] = (x, yb, yr)
    xu = camp['unin']['exploration_pct'] / 100
    _mark(ax, xu, U, 'D', S_DIAMOND_C, P['unin'], 'white', 0.7, z=6)
    xs = F['seeds']['exploration_pct'] / 100
    ys = F['seeds']['best_irr_pct']
    _mark(ax, xs, ys, 'D', 50, 'white', P['seed'], 1.3, z=6)
    xr = camp['relay']['exploration_pct'] / 100
    yr = F['relay']['best']['irr_pct']
    _mark(ax, xr, yr, '*', 150, P['relay'], 'white', 0.8, z=7)

    arr = FancyArrowPatch((xs, ys), (xr, yr), arrowstyle='-|>',
                          connectionstyle='arc3,rad=0.0', color=GREY_ARROW,
                          lw=1.0, ls=(0, (3, 2)), mutation_scale=9,
                          shrinkA=6, shrinkB=8, zorder=6.5)
    ax.add_patch(arr)
    MARKS.append(arr)

    fs, nfs = S.FS['annot'], S.FS['note']
    _ann(ax, f"{C.fmt_int(F['seeds']['n'])} seeds", (xs, ys), (-7, 0),
         ha='right', va='center', fontsize=nfs, color=GREY_LABEL)
    # group labels (data-anchored; offsets in points)
    x_left = L_F_XLIM[0] * 1.06            # just inside the left spine
    _ann(ax, 'ethanol scouts', (x_left, 7.0), (3, 0), ha='left', va='top',
         fontsize=fs, color=P['etoh_dark'])
    _ann(ax, 'uninformed', (xu, U), (5, -5), ha='left', va='top',
         fontsize=fs, color=P['unin_text'])
    _ann(ax, 'isobutanol scouts', (pos['ip'][0], pos['ip'][2]), (-7, 0),
         ha='right', va='center', fontsize=fs, color=P['ibo_dark'])
    note = dict(fontsize=nfs, fontstyle='italic', color=NOTE,
                linespacing=1.1)
    _ann(ax, 'explore + select', (xr, yr), (-20, -9), ha='right',
         va='bottom', **note)
    _ann(ax, 'TRY-informed', (xr, yr), (-20, 2), ha='right', va='bottom',
         fontsize=fs, color=P['relay'])
    # corner notes
    _ann(ax, 'profit alone: little\nexploration, plateau', (x_left, U),
         (3, 9), ha='left', va='bottom', **note)
    _ann(ax, 'TRY alone: explores,\nreturns low IRR',
         (camp['iy']['exploration_pct'] / 100, pos['iy'][2]), (-6, 0),
         ha='right', va='center', **note)
    S.style_ticks(ax)
    return ax


# %% Header-band keys ------------------------------------------------------------------
def key_c(fig):
    """'○ best visit  ● returned design' + '(argmax of its own objective)'."""
    x0 = GROUP_X
    ring = dict(marker='o', ms=6.4, mfc='white', mec=TEXT, mew=1.2)
    dot = dict(marker='o', ms=6.0, mfc=TEXT, mec=TEXT, mew=0.0)
    x1, _ = S.inline_key(fig, x0, BAND_Y[0], [{**ring, 'text': 'best visit'}])
    x2 = x1 + 0.22
    S.inline_key(fig, x2, BAND_Y[0], [{**dot, 'text': 'returned design'}])
    S.fig_text(fig, x2 + 6.0 / 72 + 0.05, BAND_Y[1],
               '= argmax of its own objective', fontsize=S.FS['key'],
               color=NOTE, va='center', ha='left')


def key_de(fig, x0):
    """Shared d/e colour key: aligned columns of swatch + text."""
    rows = [[('etoh', 'ethanol / Pdc'), ('adh1', 'Adh1'), ('gly', 'glycolysis')],
            [('ibo', 'isobutanol / ALS→Aro10'), ('adh6', 'Adh6'),
             ('tca', 'TCA/acetate')]]
    fs, sw, gap, cgap = S.FS['key'], 0.11, 0.04, 0.12
    widths = [max(sw + gap + _text_w_in(fig, r[j][1], fs) for r in rows)
              for j in range(3)]
    xs = np.r_[x0, x0 + np.cumsum(np.array(widths) + cgap)[:-1]]
    for y, r in zip(BAND_Y, rows):
        for x, (ck, txt) in zip(xs, r):
            S.inline_key(fig, x, y, [{'swatch': P[ck], 'text': txt,
                                      'size_in': sw}], fontsize=fs,
                         gap_in=gap)
    return xs[-1] + widths[-1]


def key_f(fig, x0):
    """Role-marker key for f: ○ yield □ titer △ productivity, then the
    open / filled rule."""
    it = dict(ms=6.0, mfc='white', mec=TEXT, mew=1.1)
    x1, _ = S.inline_key(fig, x0, BAND_Y[0], [
        {'marker': 'o', **it, 'text': 'yield'},
        {'marker': 's', **it, 'ms': 5.6, 'text': 'titer'},
        {'marker': '^', **it, 'text': 'productivity'}])
    x2, _ = S.inline_key(fig, x0, BAND_Y[1], [
        {'text': 'open: best visit · filled: returned', 'color': NOTE}])
    return max(x1, x2)


# %% Mark-clash check (text vs markers / lines / arrows / inset) -----------------------
def _mark_points(fig, m, renderer):
    """(N x 2 display points, radius px) of a mark artist."""
    dpi = fig.dpi
    if isinstance(m, PathCollection):
        pts = m.get_offset_transform().transform(m.get_offsets())
        r = float(np.sqrt(np.max(m.get_sizes())) / 2 * dpi / 72)
        return pts, r
    if isinstance(m, Line2D):
        xy = m.get_transform().transform(np.column_stack(m.get_data()))
        r = m.get_linewidth() / 2 * dpi / 72
    elif isinstance(m, FancyArrowPatch):
        paths = m._get_path_in_displaycoord()[0]
        paths = paths if isinstance(paths, list) else [paths]
        polys = [p for path in paths
                 for p in path.to_polygons(closed_only=False)]
        xy = np.vstack(polys)
        r = m.get_linewidth() / 2 * dpi / 72
    else:
        raise TypeError(type(m))
    # densify segments (2 px)
    out = [xy[:1]]
    for a, b in zip(xy[:-1], xy[1:]):
        if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
            continue
        n = max(int(np.hypot(*(b - a)) / 2), 1)
        out.append(a + (b - a) * np.linspace(0, 1, n + 1)[1:, None])
    return np.vstack(out), r


def mark_clashes(fig, texts, marks, regions=(), pad_px=1.0):
    """Annotation texts touching a mark (marker, connector, step line,
    reference line, arrow) or a region (an axes bbox, e.g. the inset).
    Returns a list of messages; [] = pass."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    msgs = []
    pts = [(m, *_mark_points(fig, m, r)) for m in marks]
    for t in texts:
        if not t.get_visible() or not t.get_text().strip():
            continue
        bb = t.get_window_extent(r)
        for m, xy, rad in pts:
            if getattr(t, 'arrow_patch', None) is m:
                continue
            p = pad_px + rad
            hit = ((xy[:, 0] > bb.x0 - p) & (xy[:, 0] < bb.x1 + p)
                   & (xy[:, 1] > bb.y0 - p) & (xy[:, 1] < bb.y1 + p))
            if hit.any():
                msgs.append(f'text {t.get_text()!r} touches '
                            f'{type(m).__name__} {m.get_label()!r}')
        for name, rb in regions:
            if bb.overlaps(rb):
                msgs.append(f'text {t.get_text()!r} overlaps {name}')
    for name, rb in regions:              # lines must not run under the inset
        for m, xy, rad in pts:
            if isinstance(m, Line2D) and m.axes is not None and \
                    m.get_linestyle() == '-':
                inside = ((xy[:, 0] > rb.x0) & (xy[:, 0] < rb.x1)
                          & (xy[:, 1] > rb.y0) & (xy[:, 1] < rb.y1))
                if inside.any():
                    msgs.append(f'line {m.get_label()!r} runs under {name}')
    return msgs


# %% Figure ------------------------------------------------------------------------------
def build(F):
    ANNOT.clear()
    MARKS.clear()
    fig = S.new_figure(*S.MAIN_SIZE)
    ax_a, strip, ins = panel_a(fig, F)
    ax_b = panel_b(fig, F)
    ax_c = panel_c(fig, F)
    row_labels(fig, F, ax_c)
    table_c(fig, F, ax_c)
    ax_d = panel_d(fig, F)
    ax_e = panel_e(fig, F)
    ax_f = panel_f(fig, F)
    key_c(fig)
    x_de = key_de(fig, KEY_DE_X)
    x_f0 = max(L['f'][0], x_de + KEY_GAP_IN)
    x_fk = key_f(fig, x_f0)
    print(f'header band: d/e key {KEY_DE_X:.2f}-{x_de:.2f} in, '
          f'f key {x_f0:.2f}-{x_fk:.2f} in')
    assert x_fk <= S.MAIN_SIZE[0] - 0.05, f'f key runs to {x_fk:.2f} in'

    rel = F['relay']
    titles = {
        'a': f"Seeded search beats the plateau at {C.fmt_trial(rel['first_gt_U'])}",
        'b': f"Every design above {C.fmt_irr(F['U_pct'])} makes isobutanol",
        'c': 'Scouts visit higher IRR than they return',
        'd': 'Titers', 'e': 'Proteome', 'f': 'Explore, then select'}
    for p, ttl in titles.items():
        y = TOP_LETTER_Y if p in 'ab' else BOT_LETTER_Y
        S.panel_letter(fig, LETTER_X[p], y, p)
        S.panel_title(fig, LETTER_X[p], y, ttl)
    return fig, dict(a=ax_a, strip=strip, inset=ins, b=ax_b, c=ax_c, d=ax_d,
                     e=ax_e, f=ax_f)


def run_checks(fig, axes, raise_on_fail=True):
    res = S.check_figure(fig, size=S.MAIN_SIZE, raise_on_fail=False)
    fig.canvas.draw()
    inset_bb = axes['inset'].get_tightbbox(fig.canvas.get_renderer())
    clashes = mark_clashes(fig, ANNOT, MARKS, regions=[('inset', inset_bb)])
    for m in clashes:
        S._safe_print(f'[marks] {m}')
    print(f'mark checks: {len(clashes)} problem(s)')
    n = sum(len(v) for v in res.values()) + len(clashes)
    if n and raise_on_fail:
        raise S.FigureCheckError(f'{n} figure-check problem(s)')
    return n


# %% Caption ------------------------------------------------------------------------------
def _eb():
    spec = importlib.util.spec_from_file_location('eb_by_path', EB_PATH)
    eb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(eb)
    return eb


def caption(F):
    camp, rel, seeds, g, ui = (F['campaigns'], F['relay'], F['seeds'],
                               F['global'], F['unin'])
    U = C.fmt_irr(F['U_pct'])
    start = C.fmt_irr(F['start_irr_pct'])
    rb = rel['best']
    best_seed = C.fmt_irr(seeds['best_irr_pct'])
    n_dec = len(C.DECISION_VARS)
    thr = f'{C.IBO_THRESHOLD:g}'
    unit = 'g·L⁻¹'
    pr = F['proteome']
    eb = _eb()
    housekeeping = eb.PROTEIN_CONTENT * eb.HOUSEKEEPING_FRACTION
    pi25, pi_best = rel['sims'][25]['PI'], rel['bsf_PI'][rb['sim']]
    ratio = (pi_best - pi25) / REPRO_MARGIN_PI
    ratio_word = WORDS.get(int(round(ratio)), f'{ratio:.0f}')
    n_scout = [camp[k]['n_complete'] for k in ('unin',) + C.SCOUT_KEYS]
    ibo_ret = [camp[k]['returned']['irr_pct'] for k in C.IBO_SCOUTS]
    ibo_fin = [v for v in ibo_ret if np.isfinite(v) and v >= 0]
    n_noirr = sum(1 for v in ibo_ret if not np.isfinite(v))
    etoh_ret = [camp[k]['returned']['irr_pct'] for k in C.ETOH_SCOUTS]
    etoh_bv = max(camp[k]['best_visit']['irr_pct'] for k in C.ETOH_SCOUTS)
    growth = [f"×{pr[k]['growth']:.2f}" for k in C.ROW_KEYS
              if pr[k]['growth'] < 0.99]
    rep = F['replicates']['progress']
    seeded_25 = max(rep['relay']['first_ge25'], rep['relay_rep']['first_ge25'])
    s1 = rel['sims'][1]
    s2 = rel['sims'][2]
    s13 = rel['sims'][rel['first_gt_best_seed']]
    first_word = ('lost money' if C.is_loss(s1['irr_pct']) else
                  f"reached {C.fmt_irr(s1['irr_pct'])}")
    n_hits = rel['n_best_bound_hits']
    ibo_counts = ', '.join(C.fmt_int(camp[k]['n_gt_U'])
                           for k in C.IBO_SCOUTS[:-1])
    ibo_counts += f" and {C.fmt_int(camp[C.IBO_SCOUTS[-1]]['n_gt_U'])}"
    ibo_max_bv = max(camp[k]['best_visit']['irr_pct'] for k in C.IBO_SCOUTS)
    ret_txt = (', '.join(C.fmt_irr(v) for v in ibo_fin)
               + (f' and {WORDS[n_noirr]} with no IRR' if n_noirr else ''))
    lines = [
        f"**Figure N | Process-level (TRY) campaigns scout the isobutanol "
        f"designs that a profitability-only search rarely explores, and "
        f"seeding the profitability search with their trials finds a "
        f"profitable ethanol + isobutanol co-production strain.** Eight "
        f"Gaussian-process campaigns searched the same {n_dec}-dimensional "
        f"strain-design and feeding space (scenario-A start, enzyme-burden "
        f"constraint, seed {SEED}). Two maximized profitability "
        f"(profitability index, PI = net present value at a {HURDLE_PCT} % "
        f"hurdle rate / total capital investment; shown as internal rate of "
        f"return, IRR): **uninformed** (cyan; "
        f"{C.fmt_int(camp['unin']['max_sim'])} trials, the first "
        f"{F['checks']['shared_startup_rows']} a space-filling design shared "
        f"with the TRY campaigns) and **TRY-informed** (teal; "
        f"{C.fmt_int(camp['relay']['max_sim'])} trials; its Gaussian process "
        f"was preloaded with {C.fmt_int(seeds['n'])} trials of the six TRY "
        f"campaigns). Six TRY campaigns (\"scouts\") each maximized one "
        f"titer, rate (productivity) or yield of ethanol (amber) or "
        f"isobutanol (violet) (○ yield, □ titer, △ productivity). Dashed "
        f"cyan: the uninformed plateau, {U} IRR; tinted: IRR above it; "
        f"dotted grey: the starting strain (scenario A), {start}; \"loss\": "
        f"IRR < 0 or no IRR (an outright money-loser). Titers are per litre "
        f"of water; \"making isobutanol\" means ≥ {thr} {unit}; "
        f"\"co-production\" means ≥ {thr} {unit} of each alcohol.",
        '',
        f"**(a)** IRR of the best design so far, which is also each "
        f"campaign's incumbent, vs simulated trials. The uninformed campaign "
        f"reached {U} with an ethanol-only design at "
        f"{C.fmt_trial(ui['plateau_sim'])} and did not improve in the "
        f"remaining {C.fmt_int(ui['flat_trials'])} trials. Left strip: the "
        f"{C.fmt_int(seeds['n'])} seed trials by scout family. All "
        f"{seeds['n_gt_U']} seeds above {U} come from isobutanol scouts; the "
        f"best is {best_seed}. The TRY-informed campaign's first trial "
        f"{first_word}; its second ({C.fmt_irr(s2['irr_pct'])}) already "
        f"co-produced both alcohols; its "
        f"{_ordinal(rel['first_gt_best_seed'])} "
        f"({C.fmt_irr(s13['irr_pct'])}) beat every seed; it passed 25 % at "
        f"{C.fmt_trial(rel['first_ge25'])} and "
        f"{C.fmt_irr(rel['bsf_pct'][25])} at trial 25. Inset: refinement to "
        f"{C.fmt_irr(rb['irr_pct'])} at {C.fmt_trial(rb['sim'])} (PI "
        f"{pi25:.2f} → {pi_best:.2f}, {ratio_word} times the cross-process "
        f"reproducibility margin of ΔPI ≈ {REPRO_MARGIN_PI:.2f}).",
        '',
        f"**(b)** IRR of every simulated design vs isobutanol's share of the "
        f"alcohol titer. {C.fmt_int(g['n_alcohol_ge1'])} designs; the "
        f"{C.fmt_int(g['n_alcohol_lt1'])} designs making < "
        f"{C.SHARE_MIN_ALCOHOL:g} {unit} alcohol all have no IRR and are "
        f"omitted. None of the {C.fmt_int(g['n_ibo_lt5'])} designs making < "
        f"{thr} {unit} isobutanol, in any campaign, exceeds {U}. The "
        f"uninformed campaign's {ui['ibo_ge5']['n']} isobutanol-making "
        f"designs were all worse (best "
        f"{C.fmt_irr(ui['ibo_ge5']['max_irr_pct'])}) and lost money more "
        f"often ({_pct0(ui['ibo_ge5']['pct_losing'])} vs "
        f"{_pct0(ui['ibo_lt5']['pct_losing'])}). The isobutanol scouts' best "
        f"visits (open violet markers) reached {C.fmt_irr(ibo_max_bv)} with "
        f"pure isobutanol. All {g['n_relay_gt_best_seed']} designs above "
        f"{best_seed} are TRY-informed co-production designs (★, "
        f"{C.fmt_irr(rb['irr_pct'])}: {rb['etoh']:.1f} {unit} ethanol + "
        f"{rb['ibo']:.1f} {unit} isobutanol).",
        '',
        f"**(c)** Every completed trial of each campaign (strip), its "
        f"highest-IRR visit (open) and the design it returned, the argmax of "
        f"its own objective (filled). Right: trials above {U} and the share "
        f"of trials losing money ({C.fmt_int(min(n_scout))}–"
        f"{C.fmt_int(max(n_scout))} completed trials per campaign; "
        f"TRY-informed {C.fmt_int(camp['relay']['n_complete'])}). The "
        f"isobutanol scouts visited {ibo_counts} designs above {U} (up to "
        f"{C.fmt_irr(ibo_max_bv)}) but returned designs at {ret_txt}. The "
        f"ethanol scouts never exceeded {C.fmt_irr(etoh_bv)} and returned "
        f"{min(etoh_ret):.1f}–{C.fmt_irr(max(etoh_ret))}.",
        '',
        f"**(d, e)** Titers and metabolic proteome of each returned design. "
        f"Pdc and the ALS→Aro10 enzymes (ALS, Ilv5, Ilv3, Aro10) commit "
        f"carbon to ethanol and to isobutanol. The light shades are the "
        f"terminal alcohol dehydrogenases (Adh1, Adh6), which carry flux only "
        f"when their upstream branch is expressed: e.g. Adh6 in the "
        f"ethanol-titer design and Adh1 in the isobutanol-yield design are "
        f"idle. Translation ({F['phi_T']:.3f} g·(g DCW)⁻¹) and housekeeping "
        f"({housekeeping:.3f}) sectors are the same in every design and are "
        f"omitted. Dotted line: the penalty-free budget ({F['budget']:.3f}). "
        f"Designs beyond it grow more slowly (growth factors "
        f"{', '.join(growth)}).",
        '',
        f"**(f)** Share of each campaign's trials making isobutanol (logit "
        f"scale) vs the IRR of its best visit (open) and of its returned "
        f"design (filled). ◇: the {C.fmt_int(seeds['n'])} seeds "
        f"({seeds['exploration_pct']:.1f} % making isobutanol; best "
        f"{best_seed}).",
        '',
        "In these campaigns neither objective alone returned a profitable "
        "co-production design. The TRY objectives explored the isobutanol "
        "region but returned low-IRR designs. The unseeded profitability "
        "objective stayed at an ethanol-only plateau. Seeded with the "
        "scouts' trials, profitability search reached the highest IRR found "
        "by any of the eight campaigns.",
        '',
        '*Notes.*',
        f"* Single seed. A replicate unseeded campaign (different start-up "
        f"design) left the plateau at "
        f"{C.fmt_trial(rep['unin_rep']['first_gt_U'])}, reached 25 % at "
        f"{C.fmt_trial(rep['unin_rep']['first_ge25'])} and "
        f"{C.fmt_irr(rep['unin_rep']['best_pct'])} at "
        f"{C.fmt_trial(rep['unin_rep']['best_sim'])}. Both seeded campaigns "
        f"reached 25 % within {seeded_25} trials (Fig. S1).",
        f"* The seeds are all {seeds['n_keep_above']} scout trials at least "
        f"as profitable as the starting strain (PI), plus "
        f"{seeds['n_maximin']} maximin space-filling scout trials. Their "
        f"recorded values were preloaded without re-simulation; the scouts' "
        f"{C.fmt_int(g['n_try'])} simulations are by-products of the TRY "
        f"campaigns.",
        f"* The best design lies at {WORDS.get(n_hits, n_hits)} search "
        f"bounds, so higher IRR may exist outside the searched ranges.",
        f"* IRRs are at default prices under the model version used for the "
        f"campaigns (starting strain {start}; "
        f"{C.fmt_irr(100 * START_IRR_CURRENT_MODEL)} under the current model "
        f"version).",
    ]
    # the caption's quoted numbers must agree with the figure's facts
    assert ratio > 1, ratio
    return '\n'.join(lines) + '\n'


# %% Main ----------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out-dir', default=C.OUT_DIR,
                    help='output folder (default: the publication folder)')
    ap.add_argument('--no-latest', action='store_true',
                    help='do not write the _latest copies')
    ap.add_argument('--no-caption', action='store_true',
                    help='do not (re)write the caption draft')
    ap.add_argument('--draft', action='store_true',
                    help='save even if a render check fails (layout '
                         'iteration only; exits 1 after saving)')
    args = ap.parse_args(argv)

    F = C.check_facts()                 # raises FactsMismatch on any mismatch
    S.cvd_check()
    hits = C.literal_scan()
    assert not hits, f'hard-coded annotation literals: {hits}'
    assert_claims(F)

    fig, axes = build(F)
    n_bad = run_checks(fig, axes, raise_on_fail=not args.draft)
    out = S.save(fig, STEM, out_dir=args.out_dir, latest=not args.no_latest)
    if n_bad:
        print(f'DRAFT: {n_bad} render-check problem(s); not a release render')
        sys.exit(1)
    if not args.no_caption:
        txt = caption(F)
        with open(CAPTION_PATH, 'w', encoding='utf-8') as fh:
            fh.write(txt)
        print(f'caption -> {CAPTION_PATH} ({len(txt.split())} words)')
    C.assert_sim_safe()
    print('sim-safety: OK')
    return out


if __name__ == '__main__':
    main()
