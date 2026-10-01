#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Supplementary Fig. S1 of the kinetic-BO campaign figures: the replicate /
honesty check (spec section 5.1, validated per sections 6 and 7).

(a) Best IRR so far vs simulated trials for the two unseeded (uninformed)
    profitability campaigns -- this work's seed-350 campaign (solid) and the
    2026-09-16 replicate with a different space-filling start-up design
    (dashed) -- and the two seeded (TRY-informed) campaigns: this work's
    `_rl15c111dc` relay (solid) and the 2026-09-23 `_rlba1b2315` relay seeded
    from the seven replicate scouts (dashed). A summary table above the axes
    gives the best preloaded seed and the trials to beat it (seeded only:
    the like-for-like speed, since the replicate's seeds already held a
    > 25 % design), trials to >= 25 %, the IRR after 25 trials and the best
    IRR. Every best-design marker sits at its true position; the replicate
    TRY-informed best (trial 824, next to this work's 913) is a small open
    star drawn on top of this work's larger filled one.
(b) The seven replicate scouts in main panel c's family order (ethanol TRY,
    isobutanol TRY, then price-weighted yield): every COMPLETE trial's IRR as
    a strip, the best visit (open role marker), the returned design (filled;
    argmax of the campaign's own objective) and a connector; table of trials
    above this work's uninformed plateau and trials losing (main c's two
    columns) and the returned IRR ('loss' in the loss band, as main c / S2).

Writes `kinBO_S1_replicates_<stamp>.{png,pdf}` (+ `_latest`) to
`analyses/results/publication/Optimization-figures/` and the caption draft
`kinBO_S1_replicates_caption.md` next to this script; every number of the
figure text and the caption is built from `_common.compute_facts()` (asserted
by `check_facts()`) or from the figure-local facts below (asserted against
`EXPECTED_S1`), so caption and figure never disagree.

SIM-SAFE: pandas / numpy / matplotlib only, via `_common` / `_style`; never
imports `biorefineries.*`, nskinetics, biosteam, thermosteam or optuna and
never calls `load()`.

Run::

    "C:/Users/saran/anaconda3/envs/IBO_2026/python.exe" fig_s1_replicates.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                                            # noqa: E402
from matplotlib.lines import Line2D                           # noqa: E402
from matplotlib.text import Text                              # noqa: E402

import _common as C                                           # noqa: E402
import _style as S                                            # noqa: E402

STEM = 'kinBO_S1_replicates'
CAPTION_PATH = os.path.join(C.FIGDIR, 'kinBO_S1_replicates_caption.md')
SIZE = S.S1_SIZE                                  # 10.0 x 3.9 in

# The replicate unseeded campaign's best design (trial 1,603 = CSV trial 1602)
# re-simulated to IRR 0.2729 vs the recorded 0.2726
# (analyses/reproduce_split12d_trial.py, 2026-09-19; CLAUDE.md "Structure &
# Key Files"). External to the CSVs, so a documented constant.
REPRO_IRR_TOL_MD = '3 × 10⁻⁴'
REPLICATE_RUN_DATE = '2026-09-16'                 # replicate uninformed launch
RELAY_REP_RUN_DATE = '2026-09-23'                 # _rlba1b2315 relay launch
# The TEA's hurdle rate: PI = NPV at this rate / TCI, so PI > 0 <=> IRR above
# it; 'profitable' is kept for that, 'positive IRR' for 0 < IRR < hurdle
# (CLAUDE.md, Objectives: 'PI' = NPV(15 % hurdle) / TCI)
HURDLE_PCT = 15

# %% Layout (inches, origin bottom-left) ------------------------------------
AX_Y, AX_H = 0.62, 2.04
A_X, A_W = 0.75, 4.10                             # panel a axes
B_HDR_X = 5.17                                    # b group headers (left)
B_LAB_R = 6.10                                    # b row labels (right edge)
B_MRK_X = 6.22                                    # b role marker centre
B_X = 6.36                                        # b strip axes (left)
B_TABLE_R = 9.97                                  # b table right edge
B_COL_GAP = 0.13                                  # b gap between columns
B_AX_GAP = 0.14                                   # b strip axes -> table
B_MIN_W = 1.70                                    # narrowest b strip axes
LETTER_Y = 3.68                                   # panel letters / titles
LETTER_A_X, LETTER_B_X = 0.06, 5.02
# header-band text rows (centres): r0 = column heads, r1..r5 = table rows
BAND_ROWS = (3.50, 3.355, 3.21, 3.065, 2.92, 2.775)
LINE_LW = 2.2                                     # step-line width (a)
DASH = (0, (4, 2))                                # replicate line style
# 'replicate: co-producing' (a): no leader; the label hangs directly under
# the dashed line's first >= 25 % segment, its left end just right of that
# segment's open circle. Offsets in points from the circle centre
CIRC25_S, CIRC25_LW = 30, 1.3                     # first-25 % circles
CP_GAP_X_PT = 1.5                                 # circle edge -> text left
CP_GAP_Y_PT = 1.2                                 # line edge -> text top
# Best-design markers sit at their TRUE positions. The two TRY-informed bests
# (trials 913 / 824, IRR 27.3 / 27.1 %) are ~4 pt apart on the log axis.
# Fresh round 1: the small open star on top of the large filled one read as
# one smudge; now the replicate's open star is drawn LARGE (a ring-like
# outline on a white face) and this work's smaller filled star sits inside
# it, on top: two concentric-looking stars, neither moved
STAR_S_FILLED = 100                               # this work's star (s, pt^2)
STAR_S_OPEN = 260                                 # replicate's open star
OPEN_STAR_LW = 1.2                                # its outline (pt)
OPEN_STAR_HALO_LW = 2.4                           # its white halo (pt)
DIAMOND_S = 42                                    # unseeded best diamonds
# bold 'TRY-informed (seeded)' label: top edge this far below the axes top
# (IRR points), in the empty upper-left corner above the seeded lines
SEEDED_LABEL_TOP_GAP = 1.2

# panel b rows (y increases downward); (kind, key or label, y). Family order
# as in main panel c (ethanol TRY, then isobutanol TRY), price-weighted last
B_ROWS = (
    ('header', 'etoh', 0.0),
    ('row', 'rep_ey', 0.8), ('row', 'rep_et', 1.8), ('row', 'rep_ep', 2.8),
    ('header', 'ibo', 3.8),
    ('row', 'rep_iy', 4.6), ('row', 'rep_it', 5.6), ('row', 'rep_ip', 6.6),
    ('header', 'pw', 7.6),
    ('row', 'rep_pw', 8.4),
)
B_YLIM = (8.9, -0.5)
B_HEADER_TEXT = {'ibo': 'Isobutanol TRY', 'pw': 'Price-weighted',
                 'etoh': 'Ethanol TRY'}
STRIP_HALF = 0.32

# Price-weighted-yield dots: _style's pw_light (#C4C4C4) is too close to
# ibo_light under tritanopia (CIE76 dE 7.4 < 12); #959595 passes (13.5) and
# stays a neutral grey. Figure-local override, checked by the cvd_check below.
PW_DOT = '#959595'
S1_CVD_PAIRS = S.CVD_PAIRS + (('pw', 'ibo_dark'), ('pw', 'etoh_dark'),
                              ('pw_dot', 'ibo_light'),
                              ('pw_dot', 'etoh_light'))

# %% Figure-local facts (beyond _common.EXPECTED) ----------------------------
# Values quoted only by S1 (caption / annotations) and not covered by
# _common.EXPECTED: asserted here with the printed-decimals tolerance rule.
EXPECTED_S1 = {
    'rep_seeds_n': C.E(1000),
    'rep_seeds_n_donors': C.E(7),
    'rep_seeds_by_family': C.E({'ibo': 451, 'etoh': 348, 'pw': 201}),
    'rep_seeds_n_gt_U': C.E(142),
    'rep_best_seed_pct': C.E('25.13'),
    'rep_best_seed_donor': C.E('rep_iy'),
    'relay_rep_first_gt_best_seed': C.E(19),
    'rep_startup_rows': C.E(50),
    'main_startup_rows': C.E(51),
    'rep_max_returned_pct': C.E('15.45'),
    'rep_etoh_max_visit_pct': C.E('15.31'),
    'best_classes': C.E({'relay': 'coprod', 'relay_rep': 'coprod',
                         'unin_rep': 'coprod'}),
    'unin_rep_best_ibo': C.E('41.7'),
    'unin_rep_best_etoh': C.E('29.2'),
    'relay_rep_best_ibo': C.E('42.9'),
    'relay_rep_best_etoh': C.E('11.3'),
    'n_ge25_markers': C.E(3),
    # product class of every finite best-so-far step (incumbent improvement)
    'inc_classes': C.E({'unin': {'etoh': 19}, 'unin_rep': {'coprod': 21}}),
    'unin_first_profitable_sim': C.E(53),
    'unin_rep_first_profitable_sim': C.E(12),
    'unin_rep_first_profitable_pct': C.E('9.63'),
    'unin_rep_first_profitable_ibo': C.E('28.0'),
    'unin_rep_first_profitable_etoh': C.E('62.1'),
    'unin_rep_second_step_sim': C.E(52),
}


def _incumbent_steps(key):
    """COMPLETE rows that raise the finite best-so-far IRR (sim order)."""
    c = C.complete(key).sort_values('sim')
    irr = c['IRR'].to_numpy(float)
    run = np.maximum.accumulate(irr)
    new = np.r_[True, irr[1:] > run[:-1]] & np.isfinite(irr)
    return c[new]


def _close(computed, exp):
    v = exp.value
    if isinstance(v, (bool, str, dict, list)) or v is None:
        return computed == v
    c = float(computed)
    if not np.isfinite(v):
        return c == v
    return abs(c - float(v)) <= exp.tol + 1e-12


def s1_facts(facts):
    """Numbers quoted by S1 beyond facts['replicates'], recomputed from the
    CSVs (relay manifests, start-up designs, best designs) and asserted
    against EXPECTED_S1."""
    U = facts['U']
    man = C.load_manifest('relay_rep')
    mi = man['IRR'].to_numpy(float)
    bs = man.loc[man['IRR'].idxmax()]
    rel_rep = C.complete('relay_rep')
    dec = list(C.DECISION_VARS)
    rep = facts['replicates']['scouts']
    best = {k: C.row_record(C.best_visit_row(k))
            for k in ('relay', 'relay_rep', 'unin_rep')}
    prog = facts['replicates']['progress']
    F = {
        'rep_seeds_n': int(len(man)),
        'rep_seeds_n_donors': int(man['donor_key'].nunique()),
        'rep_seeds_by_family': {str(k): int(v) for k, v in
                                man['family'].value_counts().items()},
        'rep_seeds_n_gt_U': int((mi > U + C.ABOVE_EPS).sum()),
        'rep_best_seed_pct': 100.0 * float(bs['IRR']),
        'rep_best_seed_donor': str(bs['donor_key']),
        'relay_rep_first_gt_best_seed': C.first_sim(
            rel_rep, float(bs['IRR']) + C.ABOVE_EPS, strict=True),
        'rep_startup_rows': C._leading_identical_rows(
            ('unin_rep',) + C.REP_SCOUT_KEYS, dec),
        'main_startup_rows': C._leading_identical_rows(
            ('unin',) + C.SCOUT_KEYS, dec),
        'rep_max_returned_pct': max(v['returned_irr_pct']
                                    for v in rep.values()),
        'rep_etoh_max_visit_pct': max(rep[k]['best_visit_irr_pct'] for k in
                                      ('rep_ey', 'rep_et', 'rep_ep')),
        'best_classes': {k: v['class'] for k, v in best.items()},
        'unin_rep_best_ibo': best['unin_rep']['ibo'],
        'unin_rep_best_etoh': best['unin_rep']['etoh'],
        'relay_rep_best_ibo': best['relay_rep']['ibo'],
        'relay_rep_best_etoh': best['relay_rep']['etoh'],
        'n_ge25_markers': sum(prog[k]['first_ge25'] is not None for k in
                              ('unin', 'unin_rep', 'relay', 'relay_rep')),
        'best': best,
    }
    steps = {k: _incumbent_steps(k) for k in ('unin', 'unin_rep')}
    F['inc_classes'] = {k: {str(a): int(b) for a, b in
                            v['product_class'].value_counts().items()}
                        for k, v in steps.items()}
    prof = {k: v[v['IRR'] >= 0.0].iloc[0] for k, v in steps.items()}
    F['unin_first_profitable_sim'] = int(prof['unin']['sim'])
    F['unin_rep_first_profitable_sim'] = int(prof['unin_rep']['sim'])
    F['unin_rep_first_profitable_pct'] = 100.0 * float(prof['unin_rep']['IRR'])
    F['unin_rep_first_profitable_ibo'] = float(prof['unin_rep']['ibo'])
    F['unin_rep_first_profitable_etoh'] = float(prof['unin_rep']['etoh'])
    after = steps['unin_rep'][steps['unin_rep']['sim']
                              > F['unin_rep_first_profitable_sim']]
    F['unin_rep_second_step_sim'] = int(after['sim'].iloc[0])
    F['n_steps'] = {k: int(len(v)) for k, v in steps.items()}
    # "co-producing from trial 12" / "ethanol only" / "from its start-up"
    if not (F['unin_rep_first_profitable_sim']
            <= F['rep_startup_rows'] < F['unin_rep_second_step_sim']):
        raise C.FactsMismatch('replicate first profitable step is not the '
                              'only one inside its start-up design')
    bad = [f'{k}: computed {F[k]!r} vs expected {e.value!r} (tol {e.tol})'
           for k, e in EXPECTED_S1.items() if not _close(F[k], e)]
    # the claims the figure / caption make in words
    if not F['rep_max_returned_pct'] < 100.0 * U:
        bad.append('a replicate scout returned a design above the plateau')
    if any(rep[k]['n_gt_U'] for k in ('rep_ey', 'rep_et', 'rep_ep')):
        bad.append('an ethanol replicate scout visited above the plateau')
    if not all(rep[k]['n_gt_U'] > 0 for k in ('rep_iy', 'rep_it', 'rep_ip',
                                               'rep_pw')):
        bad.append('an isobutanol / price-weighted scout never beat U')
    if any(rep[k]['best_visit_irr_pct'] < rep[k]['returned_irr_pct']
           for k in rep):
        bad.append('a best visit below its returned design')
    if facts['replicates']['unin_rep_startup_identical']:
        bad.append('replicate start-up identical to this work')
    # caption: the replicate's best seed was already above 25 %, so its
    # 'trials to 25 %' and 'trials to beat every seed' are the same event
    prr = facts['replicates']['progress']['relay_rep']
    if not (F['rep_best_seed_pct'] >= 25.0
            and prr['first_ge25'] == F['relay_rep_first_gt_best_seed']):
        bad.append('replicate best seed / first >= 25 % wording is stale')
    if facts['seeds']['best_irr_pct'] >= 25.0:
        bad.append("this work's best seed reached 25 %")
    if bad:
        raise C.FactsMismatch('S1 facts: ' + '; '.join(bad))
    return F


# %% Formatting helpers ----------------------------------------------------------
MINUS = '−'
LOSS_WORD = C.fmt_irr(float('-inf'))              # 'loss', as main c / S2


def fmt_signed_irr(v, nd=1):
    """IRR in % for a table cell: '15.5 %', '−6.0 %'; -inf -> 'no IRR'."""
    v = float(v)
    if np.isneginf(v):
        return 'no IRR'
    s = f'{abs(v):.{nd}f} %'
    return (MINUS + s) if v < 0 else s


def _irr_words(v):
    """Caption wording of an IRR in %: '9.6 %', '−11.2 % (a loss)',
    'no IRR (a loss)'."""
    s = fmt_signed_irr(v)
    return s if C.fmt_irr(v) != 'loss' else f'{s} (a loss)'


def fmt_trials(sim):
    return 'never' if sim is None else C.fmt_int(sim)


def pair(a, b, sep=' / '):
    return f'{a}{sep}{b}'


def _text_w_in(fig, s, **kw):
    t = Text(0, 0, s, **kw)
    t.set_figure(fig)
    return t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi


def _fig_marker(fig, x, y, marker, ms, mfc, mec, mew=1.0, zorder=5):
    ln = Line2D([x], [y], transform=fig.dpi_scale_trans, ls='none',
                marker=marker, ms=ms, mfc=mfc, mec=mec, mew=mew,
                zorder=zorder)
    fig.add_artist(ln)
    return ln


def _fig_line(fig, x0, x1, y, **kw):
    ln = Line2D([x0, x1], [y, y], transform=fig.dpi_scale_trans, **kw)
    fig.add_artist(ln)
    return ln


def _segment_hits_box(p0, p1, box):
    """Liang-Barsky: does segment p0-p1 (display px) enter box (x0, y0, x1,
    y1)?"""
    (x0, y0), (x1, y1) = p0, p1
    dx, dy = x1 - x0, y1 - y0
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, x0 - box[0]), (dx, box[2] - x0),
                 (-dy, y0 - box[1]), (dy, box[3] - y0)):
        if p == 0:
            if q < 0:
                return False
            continue
        r = q / p
        if p < 0:
            t0 = max(t0, r)
        else:
            t1 = min(t1, r)
        if t0 > t1:
            return False
    return True


def text_line_clashes(fig, ax, pad_px=0.5):
    """Texts of `ax` whose window extent meets one of the axes' drawn lines
    (step lines incl. their risers, reference lines), the line's half-width
    added as padding. text_overlaps only checks text against text; this is
    the text-against-data check for panel a. Returns a list of messages."""
    from matplotlib.cbook import STEP_LOOKUP_MAP
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    msgs = []
    for line in ax.lines:
        xy = line.get_xydata()
        ds = line.get_drawstyle()
        if ds in STEP_LOOKUP_MAP and ds != 'default':
            xy = np.column_stack(STEP_LOOKUP_MAP[ds](xy[:, 0], xy[:, 1]))
        pts = line.get_transform().transform(xy)
        half = line.get_linewidth() * fig.dpi / 72.0 / 2.0
        for t in ax.texts:
            bb = t.get_window_extent(r)
            box = (bb.x0 - half - pad_px, bb.y0 - half - pad_px,
                   bb.x1 + half + pad_px, bb.y1 + half + pad_px)
            if any(_segment_hits_box(pts[i], pts[i + 1], box)
                   for i in range(len(pts) - 1)
                   if np.all(np.isfinite(pts[i:i + 2]))):
                msgs.append(f'text {t.get_text()!r} meets line '
                            f'{line.get_label()!r}')
    return msgs


def text_marker_clashes(fig, ax, pad_px=1.0):
    """Texts of `ax` whose window extent meets a scatter marker (its
    bounding circle, radius sqrt(s)/2 pt + edge); the marker counterpart of
    text_line_clashes. Returns a list of messages."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    msgs = []
    for coll in ax.collections:
        offs = coll.get_offsets()
        if len(offs) == 0:
            continue
        sizes = np.broadcast_to(coll.get_sizes(), (len(offs),))
        lws = np.broadcast_to(coll.get_linewidths(), (len(offs),))
        pts = coll.get_offset_transform().transform(offs)
        rad = (np.sqrt(sizes) / 2.0 + lws / 2.0) * fig.dpi / 72.0
        for t in ax.texts:
            bb = t.get_window_extent(r)
            for (px, py), rr in zip(pts, rad):
                dx = max(bb.x0 - px, 0.0, px - bb.x1)
                dy = max(bb.y0 - py, 0.0, py - bb.y1)
                if np.hypot(dx, dy) < rr + pad_px:
                    msgs.append(f'text {t.get_text()!r} meets a marker at '
                                f'({px:.0f}, {py:.0f}) px')
    return msgs


def text_spine_clashes(fig, ax, min_pt=2.0):
    """Texts of `ax` closer than `min_pt` to (or across) a spine."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ab = ax.get_window_extent(r)
    m = min_pt * fig.dpi / 72.0
    return [f'text {t.get_text()!r} within {min_pt} pt of a spine'
            for t in ax.texts
            for bb in (t.get_window_extent(r),)
            if (bb.x0 - ab.x0 < m or ab.x1 - bb.x1 < m
                or bb.y0 - ab.y0 < m or ab.y1 - bb.y1 < m)]


# %% Panel a -------------------------------------------------------------------
A_LINES = (   # draw order: TRY-informed (this work) last, on top
    ('unin', 'unin', '-'), ('unin_rep', 'unin', DASH),
    ('relay_rep', 'relay', DASH), ('relay', 'relay', '-'),
)


def draw_panel_a(fig, facts, F):
    P = S.PALETTE
    prog = facts['replicates']['progress']
    ax = S.inch_axes(fig, A_X, AX_Y, A_W, AX_H)
    S.log_trial_axis(ax)
    S.irr_axis(ax, 'y', step=5)
    S.start_line(ax, facts['start_irr_pct'])
    S.style_ticks(ax)

    series = {}
    for key, ckey, ls in A_LINES:
        s, b = C.best_so_far(C.complete(key))
        y = S.irr_plot(b)                         # losses at -3, no jitter
        series[key] = (s, b)
        ax.step(s, y, where='post', color=P[ckey], lw=LINE_LW, ls=ls,
                solid_capstyle='butt', dash_capstyle='butt', zorder=3,
                label=key)

    # open circle at the first trial >= 25 % (none for this work's unseeded)
    for key, ckey, _ in A_LINES:
        f25 = prog[key]['first_ge25']
        if f25 is None:
            continue
        s, b = series[key]
        ax.scatter([f25], [100.0 * C.bsf_at(s, b, f25)], s=CIRC25_S,
                   marker='o', facecolor='white', edgecolor=P[ckey],
                   lw=CIRC25_LW, zorder=6)
    # best design: filled marker for this work, white-faced for the
    # replicate, every one at its true (trial, IRR). The two TRY-informed
    # stars nearly coincide: this work's is larger, the replicate's smaller
    # open star is drawn on top of it (no marker sits at a false position)
    for key, ckey, ls in A_LINES:
        mk = '*' if C.campaign(key).is_relay else 'D'
        filled = not C.campaign(key).replicate
        if mk == '*':
            size = STAR_S_FILLED if filled else STAR_S_OPEN
        else:
            size = DIAMOND_S
        x_b, y_b = prog[key]['best_sim'], prog[key]['best_pct']
        if mk == '*' and not filled:
            # the replicate's LARGE open star: white face (hides the lines
            # under it), outline on a white halo, all UNDER this work's
            # smaller filled star (z 7.6), which sits inside it
            ax.scatter([x_b], [y_b], s=size, marker=mk, facecolor='white',
                       edgecolor='none', lw=0, zorder=6.8)
            for ec, lw, z in (('white', OPEN_STAR_HALO_LW, 6.9),
                              (P[ckey], OPEN_STAR_LW, 7.0)):
                ax.scatter([x_b], [y_b], s=size, marker=mk,
                           facecolor='none', edgecolor=ec, lw=lw, zorder=z)
            continue
        ax.scatter([x_b], [y_b], s=size,
                   marker=mk, facecolor=P[ckey] if filled else 'white',
                   edgecolor='white' if filled else P[ckey],
                   lw=0.8 if filled else 1.2, zorder=7.6 if filled else 7.5)

    # direct labels; the seeded label hangs from just below the axes top,
    # in the empty corner left of the seeded lines' climb
    ax.text(1.18, S.IRR_LIM[1] - SEEDED_LABEL_TOP_GAP,
            'TRY-informed (seeded)', color=P['relay_text'],
            fontsize=S.FS['annot'], fontweight='bold', ha='left',
            va='top', zorder=8)
    U_pct = facts['U_pct']
    x_end = facts['campaigns']['unin']['max_sim']
    # product labels of the two unseeded lines, each naming its line
    # ("this work" / "replicate") and attached to it: the solid line's label
    # sits on the line; the dashed line's hangs directly under it, starting
    # just right of the line's first-25 % circle, where the dashed line is
    # the nearest line above the text (the TRY-informed lines run higher)
    ax.text(x_end, U_pct + 0.45, 'this work: ethanol only',
            color=P['unin_text'], fontsize=S.FS['note'], style='italic',
            ha='right', va='bottom', zorder=8)
    ax.text(x_end, U_pct + 2.9, 'uninformed (unseeded)',
            color=P['unin_text'], fontsize=S.FS['annot'], fontweight='bold',
            ha='right', va='bottom', zorder=8)
    s_r, b_r = series['unin_rep']
    x_cp = prog['unin_rep']['first_ge25']            # its open circle
    if x_cp is None:
        raise C.FactsMismatch('replicate uninformed never reached 25 %')
    y_cp = 100.0 * C.bsf_at(s_r, b_r, x_cp)
    # the best-so-far line only rises, so the segment at the label's left
    # end is the lowest stretch of line above the whole label
    circ_r_pt = np.sqrt(CIRC25_S) / 2.0 + CIRC25_LW / 2.0   # its radius
    ax.annotate('replicate: co-producing', xy=(x_cp, y_cp),
                xycoords='data', textcoords='offset points',
                xytext=(circ_r_pt + CP_GAP_X_PT,
                        -(LINE_LW / 2.0 + CP_GAP_Y_PT)),
                color=P['unin_text'], fontsize=S.FS['note'],
                style='italic', ha='left', va='top', zorder=8)
    ax.text(x_end, facts['start_irr_pct'] - 0.7,
            f"starting strain {C.fmt_pct(facts['start_irr_pct'])}",
            color=S.NOTE, fontsize=S.FS['note'], ha='right', va='top',
            zorder=8)

    # --- header band: summary table (values: this work / replicate). The
    # seed rows make the seeded speed comparable: the replicate's seed pool
    # already held a design at 25 %, so 'trials to beat it' (its best seed)
    # is the like-for-like speed of the two seeded campaigns
    fs = S.FS['note']
    r0 = BAND_ROWS[0]
    lab_kw = dict(fontsize=fs, color=S.TEXT)
    seeds = facts['seeds']
    rows = (
        ('best seed IRR',
         pair(C.fmt_pct(seeds['best_irr_pct']).replace(' %', ''),
              C.fmt_pct(F['rep_best_seed_pct'])),
         'no seeds'),
        ('trials to beat it',
         pair(fmt_trials(facts['relay']['first_gt_best_seed']),
              fmt_trials(F['relay_rep_first_gt_best_seed'])),
         '–'),
        ('trials to ≥ 25 %',
         pair(fmt_trials(prog['relay']['first_ge25']),
              fmt_trials(prog['relay_rep']['first_ge25'])),
         pair(fmt_trials(prog['unin']['first_ge25']),
              fmt_trials(prog['unin_rep']['first_ge25']))),
        ('IRR after 25 trials',
         pair(C.fmt_irr(prog['relay']['at25_pct']).replace(' %', ''),
              C.fmt_irr(prog['relay_rep']['at25_pct'])),
         pair(C.fmt_irr(prog['unin']['at25_pct']).replace(' %', ''),
              C.fmt_irr(prog['unin_rep']['at25_pct']))),
        ('best IRR',
         pair(C.fmt_irr(prog['relay']['best_pct']).replace(' %', ''),
              C.fmt_irr(prog['relay_rep']['best_pct'])),
         pair(C.fmt_irr(prog['unin']['best_pct']).replace(' %', ''),
              C.fmt_irr(prog['unin_rep']['best_pct']))),
    )
    glyph_gap, glyph_w = 0.06, 0.10
    # glyphs tying rows to plot markers: o (first >= 25 %), star + diamond
    n_glyphs = tuple({'trials to ≥ 25 %': 1, 'best IRR': 2}.get(r[0], 0)
                     for r in rows)
    lab_w = max(_text_w_in(fig, r[0], **lab_kw)
                + (glyph_gap + n * glyph_w if n else 0.0)
                for r, n in zip(rows, n_glyphs))
    x_c1 = A_X + lab_w + 0.14
    c1_w = max(max(_text_w_in(fig, r[1], **lab_kw) for r in rows),
               _text_w_in(fig, 'TRY-informed', fontsize=fs,
                          fontweight='bold'))
    x_c2 = x_c1 + c1_w + 0.20
    S.fig_text(fig, A_X, r0, 'this work / replicate', fontsize=fs,
               color=S.NOTE, style='italic', va='center', ha='left')
    S.fig_text(fig, x_c1, r0, 'TRY-informed', fontsize=fs, fontweight='bold',
               color=P['relay_text'], va='center', ha='left')
    S.fig_text(fig, x_c2, r0, 'uninformed', fontsize=fs, fontweight='bold',
               color=P['unin_text'], va='center', ha='left')
    row_y = BAND_ROWS[1:1 + len(rows)]
    for (lab, v1, v2), y in zip(rows, row_y):
        S.fig_text(fig, A_X, y, lab, va='center', ha='left', **lab_kw)
        S.fig_text(fig, x_c1, y, v1, fontsize=fs, color=P['relay_text'],
                   va='center', ha='left')
        no_seed = v2 in ('no seeds', '–')
        S.fig_text(fig, x_c2, y, v2, fontsize=fs,
                   color=S.NOTE if no_seed else P['unin_text'],
                   style='italic' if v2 == 'no seeds' else 'normal',
                   va='center', ha='left')
    g = '#444444'
    for (lab, _, _), y, n in zip(rows, row_y, n_glyphs):
        if not n:
            continue
        x_g = A_X + _text_w_in(fig, lab, **lab_kw) + glyph_gap + 0.05
        if n == 1:
            _fig_marker(fig, x_g - 0.005, y, 'o', 5.0, 'white', g, mew=1.1)
        else:
            _fig_marker(fig, x_g, y, '*', 8.0, g, g, mew=0.4)
            _fig_marker(fig, x_g + glyph_w, y, 'D', 4.2, g, g, mew=0.4)

    # line-style key, right-aligned to the axes, on the two seed rows (their
    # uninformed cells are short)
    x_r = A_X + A_W
    seg = 0.34
    key_kw = dict(fontsize=fs, color=S.TEXT, va='center', ha='left')
    key_rows = (('this work', '-', row_y[0]), ('replicate', DASH, row_y[1]))
    tw = max(_text_w_in(fig, lab, fontsize=fs) for lab, _, _ in key_rows)
    x0 = x_r - tw - 0.06 - seg                    # one column of segments
    key_cells_right = x_c2 + max(_text_w_in(fig, r[2], fontsize=fs)
                                 for r in rows[:2])
    if x0 < key_cells_right + 0.15:
        raise S.FigureCheckError('S1a line-style key meets the table')
    for lab, ls, y in key_rows:
        _fig_line(fig, x0, x0 + seg, y, color='#444444', lw=LINE_LW, ls=ls,
                  solid_capstyle='butt', dash_capstyle='butt')
        S.fig_text(fig, x0 + seg + 0.06, y, lab, **key_kw)
    return ax, (x_c1, x_c2)


# %% Panel b -------------------------------------------------------------------
def _family_colors(fam):
    P = S.PALETTE
    if fam == 'pw':      # grey dots pile up darker: lower alpha
        return {'dot': PW_DOT, 'dark': P['pw'], 'alpha': 0.3}
    return {'dot': P[fam + '_light'], 'dark': P[fam + '_dark'], 'alpha': 0.5}


def _pct0(v):
    """Share in % for a table cell (as main panel c): '76 %'."""
    return f'{float(v):.0f} %'


def b_table_columns(fig, facts):
    """The panel-b table: main panel c's two columns ('trials > U',
    'trials losing') then S1's 'returned IRR', laid out right to left from
    B_TABLE_R by measured widths. Returns [(header, {key: (text, bold)},
    centre_in, value_right_in)] and the left edge of the first column."""
    rep = facts['replicates']['scouts']
    keys = [k for kind, k, _ in B_ROWS if kind == 'row']
    max_gt = max(rep[k]['n_gt_U'] for k in keys)
    cols = [
        (f'trials\n> {C.fmt_pct(facts["U_pct"])}',
         {k: (C.fmt_int(rep[k]['n_gt_U']), rep[k]['n_gt_U'] == max_gt)
          for k in keys}),
        ('trials\nlosing',
         {k: (_pct0(rep[k]['pct_losing']), False) for k in keys}),
        # 'loss' for every design drawn in the loss band (IRR < 0 or no
        # IRR), as main panel c / f and Fig. S2 print it; the caption
        # gives the values behind each 'loss'
        ('returned\nIRR',
         {k: (C.fmt_irr(rep[k]['returned_irr_pct']), False)
          for k in keys}),
    ]
    out, right = [], B_TABLE_R
    for head, vals in reversed(cols):
        wh = max(_text_w_in(fig, ln, fontsize=S.FS['note'])
                 for ln in head.split('\n'))
        wv = max(_text_w_in(fig, t, fontsize=S.FS['table'],
                            fontweight='bold' if b else 'normal')
                 for t, b in vals.values())
        w = max(wh, wv)
        cx = right - w / 2.0
        out.append((head, vals, cx, cx + wv / 2.0))
        right -= w + B_COL_GAP
    return out[::-1], right + B_COL_GAP


def draw_panel_b(fig, facts, F):
    P = S.PALETTE
    U_pct = facts['U_pct']
    rep = facts['replicates']['scouts']
    cols, table_left = b_table_columns(fig, facts)
    b_w = table_left - B_AX_GAP - B_X
    if b_w < B_MIN_W:
        raise S.FigureCheckError(f'panel b strip axes only {b_w:.2f} in')
    ax = S.inch_axes(fig, B_X, AX_Y, b_w, AX_H)
    # fresh round 1: 'loss' at the tick size (it was 9 pt) and '0'
    # labelled, as in main panel c: 'loss' moves left inside the band
    S.irr_axis(ax, 'x', step=10, loss_at=S.loss_label_x(b_w))
    ax.set_ylim(*B_YLIM)
    S.tint_above(ax, U_pct, which='x')
    S.plateau_line(ax, U_pct, which='x')
    S.style_ticks(ax, y=False)

    rng = S.panel_rng()
    n_n = {}
    for kind, key, y in B_ROWS:
        if kind != 'row':
            continue
        c = C.campaign(key)
        cols_f = _family_colors(c.family)
        d = C.complete(key)
        xs = S.irr_plot(d['IRR'].to_numpy(float), rng)
        ys = y + rng.uniform(-STRIP_HALF, STRIP_HALF, len(xs))
        ax.scatter(xs, ys, s=2.5, color=cols_f['dot'], alpha=cols_f['alpha'],
                   lw=0, rasterized=True, zorder=2)
        n_n[key] = (int(np.sum((xs >= S.LOSS_BAND[0])
                               & (xs <= S.LOSS_BAND[1]))),
                    int(C.is_loss(d['IRR'].to_numpy(float)).sum()),
                    int(round(rep[key]['pct_losing'] * len(d) / 100.0)))
        bv = S.irr_plot(rep[key]['best_visit_irr_pct'] / 100.0)
        rt = S.irr_plot(rep[key]['returned_irr_pct'] / 100.0)
        ax.plot([rt, bv], [y, y], color='0.25', lw=1.1, zorder=3,
                solid_capstyle='butt')
        ax.scatter([bv], [y], s=44, marker=c.marker, facecolor='white',
                   edgecolor=cols_f['dark'], lw=1.4, zorder=5)
        # white edge lifts the filled marker off a dense strip
        ax.scatter([rt], [y], s=40, marker=c.marker, facecolor=cols_f['dark'],
                   edgecolor='white', lw=0.7, zorder=6)
    # every loss is drawn in the band, and nothing else is (spec 6.B.6); the
    # 'trials losing' column counts exactly those dots
    for key, (in_band, n_loss, n_col) in n_n.items():
        if not in_band == n_loss == n_col:
            raise C.FactsMismatch(f'{key}: {in_band} dots in the loss band '
                                  f'vs {n_loss} losses vs {n_col} in table')

    # row labels, headers, markers, separators, table
    y2fig = lambda yy: (AX_Y + AX_H * (B_YLIM[0] - yy)  # noqa: E731
                        / (B_YLIM[0] - B_YLIM[1]))
    for kind, key, y in B_ROWS:
        yf = y2fig(y)
        if kind == 'header':
            col = _family_colors(key)['dark']
            S.fig_text(fig, B_HDR_X, yf, B_HEADER_TEXT[key], color=col,
                       fontsize=S.FS['row'], fontweight='bold', va='center',
                       ha='left')
            if y > 0:
                ys = y2fig(y - 0.45)
                _fig_line(fig, B_HDR_X, B_TABLE_R + 0.01, ys, color='0.85',
                          lw=0.6)
            continue
        c = C.campaign(key)
        cols_f = _family_colors(c.family)
        S.fig_text(fig, B_LAB_R, yf, c.role, fontsize=S.FS['row'],
                   color=S.TEXT, va='center', ha='right')
        ms = 7.0 if c.marker != '^' else 7.6
        _fig_marker(fig, B_MRK_X, yf, c.marker, ms, 'white', cols_f['dark'],
                    mew=1.3)
        for _, vals, _, x_r in cols:
            txt, bold = vals[key]
            S.fig_text(fig, x_r, yf, txt, fontsize=S.FS['table'],
                       color=S.NOTE if txt == LOSS_WORD else S.TEXT,
                       va='center', ha='right',
                       fontweight='bold' if bold else 'normal')

    # header band: marker key, table headers, plateau label
    fs = S.FS['note']
    r0 = BAND_ROWS[0]
    g = '#444444'
    S.inline_key(fig, B_HDR_X, r0, [
        {'marker': 'o', 'ms': 6.5, 'mfc': 'white', 'mec': g, 'mew': 1.3,
         'text': 'best visit'},
        {'marker': 'o', 'ms': 6.0, 'mfc': g, 'mec': g, 'mew': 0.6,
         'text': 'returned design (argmax of own objective)'},
    ], fontsize=fs)
    y_head = (BAND_ROWS[2] + BAND_ROWS[3]) / 2 + 0.005
    for head, _, cx, _ in cols:
        S.fig_text(fig, cx, y_head, head, fontsize=fs, color=S.NOTE,
                   va='center', ha='center', linespacing=1.15)
    # the reference is THIS work's plateau (S1a: the plateau did not
    # replicate), named as such over its line
    x_u = B_X + b_w * (U_pct - S.IRR_LIM[0]) / (S.IRR_LIM[1] - S.IRR_LIM[0])
    S.fig_text(fig, x_u, BAND_ROWS[-1],
               f"this work's uninformed plateau {C.fmt_pct(U_pct)}",
               fontsize=fs, color=P['unin_text'], va='center', ha='center')
    return ax


# %% Caption -----------------------------------------------------------------------
def build_caption(facts, F):
    prog = facts['replicates']['progress']
    rep = facts['replicates']['scouts']
    walk = facts['replicates']['unin_rep_walk']
    seeds = facts['seeds']
    U = C.fmt_pct(facts['U_pct'])
    pr, prr = prog['relay'], prog['relay_rep']
    pu, pur = prog['unin'], prog['unin_rep']
    w_esc = walk[min(walk)]              # first walk step = leaves plateau
    w_25 = walk[max(walk)]               # last walk step = first >= 25 %
    if w_esc['sim'] != pur['first_gt_U'] or w_25['sim'] != pur['first_ge25']:
        raise C.FactsMismatch('unin_rep walk does not match its progress')
    fam = F['rep_seeds_by_family']
    best = F['best']

    def g(v):
        return f'{v:.1f} g·L⁻¹'

    n3 = [C.fmt_int(rep[k]['n_gt_U']) for k in ('rep_iy', 'rep_it', 'rep_ip')]
    ibo_n = f'{n3[0]}, {n3[1]} and {n3[2]}'

    # the values behind each 'loss' of S1b's returned-IRR column
    fam_word = {'ibo': 'isobutanol', 'etoh': 'ethanol',
                'pw': 'price-weighted'}
    b_keys = [k for kind, k, _ in B_ROWS if kind == 'row']

    def _and(words):
        return words[0] if len(words) == 1 else (', '.join(words[:-1])
                                                 + ' and ' + words[-1])

    def _scouts(keys):
        """'price-weighted-yield scout', 'isobutanol titer and
        productivity scouts' (grouped by family, in B_ROWS order)."""
        fams = {}
        for k in keys:
            fams.setdefault(C.campaign(k).family, []).append(
                C.campaign(k).role)
        names = [f'{fam_word[f]}-{r[0]} scout' if len(r) == 1 else
                 f'{fam_word[f]} {_and(r)} scouts' for f, r in fams.items()]
        return _and(names)

    ret = {k: float(rep[k]['returned_irr_pct']) for k in b_keys}
    ret_none = [k for k in b_keys if np.isneginf(ret[k])]
    ret_neg = [k for k in b_keys if np.isfinite(ret[k]) and ret[k] < 0]
    if not ret_none and not ret_neg:
        raise C.FactsMismatch('S1b prints no returned loss: caption stale')
    parts = []
    if ret_none:
        parts.append(f'the {_scouts(ret_none)} returned '
                     + ('an outright money-loser' if len(ret_none) == 1
                        else 'outright money-losers') + ' (no IRR)')
    if ret_neg:
        parts.append(f'the {_scouts(ret_neg)} returned '
                     + ' and '.join(fmt_signed_irr(ret[k]) for k in ret_neg))
    loss_ret = ' and '.join(parts)
    # "first positive-IRR design after this work's start-up design"
    if not 0.0 < F['unin_rep_first_profitable_pct'] < HURDLE_PCT:
        raise C.FactsMismatch('replicate first positive-IRR design is not '
                              'below the hurdle: caption stale')
    if not F['unin_first_profitable_sim'] > F['main_startup_rows']:
        raise C.FactsMismatch("this work's first positive-IRR step is "
                              "inside its start-up design: caption stale")
    lines = [
        '# Figure S1 | Replicate check',
        '',
        '<!-- Generated by fig_s1_replicates.py from the campaign CSVs; '
        'do not edit by hand (re-run the script). -->',
        '',
        '**Figure S1 | Replicate check: seeding is fast in both runs; the '
        'unseeded plateau did not replicate.**',
        '',
        f'**(a)** Best IRR found so far against simulated trials (log '
        f'scale) for the two unseeded (uninformed) profitability campaigns '
        f'and the two seeded (TRY-informed) ones. Solid lines are this '
        f'work\'s campaigns: the uninformed campaign and the TRY-informed '
        f'campaign seeded with {C.fmt_int(seeds["n"])} trials of the six '
        f'TRY scouts ({C.fmt_int(seeds["n_ibo"])} isobutanol-scout and '
        f'{C.fmt_int(seeds["n_etoh"])} ethanol-scout trials: '
        f'{C.fmt_int(seeds["n_keep_above"])} of the '
        f'{C.fmt_int(seeds["n_eligible_above"])} scout trials with a '
        f'profitability index at least the starting strain\'s, '
        f'{seeds["n_quarantined_above"]} convergence-quarantined trials and '
        f'{seeds["n_duplicate_above"]} duplicate dropped, plus '
        f'{C.fmt_int(seeds["n_maximin"])} space-filling ones). Dashed lines '
        f'are replicates: an uninformed campaign (run {REPLICATE_RUN_DATE}) '
        f'whose {F["rep_startup_rows"]}-trial space-filling start-up design, '
        f'shared with seven replicate scouts, differs from this work\'s '
        f'({F["main_startup_rows"]} trials, shared with the six scouts), and '
        f'a TRY-informed campaign (run {RELAY_REP_RUN_DATE}) seeded with '
        f'{C.fmt_int(F["rep_seeds_n"])} trials of those seven scouts '
        f'({C.fmt_int(fam["ibo"])} isobutanol-scout, '
        f'{C.fmt_int(fam["etoh"])} ethanol-scout and {C.fmt_int(fam["pw"])} '
        f'price-weighted-yield trials). Open circles: first trial with IRR '
        f'≥ 25 %; stars and diamonds: each campaign\'s best design (filled '
        f'for this work, open for the replicate), every marker at its true '
        f'position. The two seeded campaigns found their best designs at '
        f'trials {C.fmt_int(pr["best_sim"])} and '
        f'{C.fmt_int(prr["best_sim"])} ({C.fmt_pct(pr["best_pct"])} and '
        f'{C.fmt_pct(prr["best_pct"])}), so their stars nearly coincide: '
        f'this work\'s smaller filled star is drawn inside the replicate\'s '
        f'larger open star. Italic labels give the products of every '
        f'best-so-far '
        f'design of the two unseeded campaigns (co-producing = at least '
        f'{C.IBO_THRESHOLD:g} g·L⁻¹ each of isobutanol and ethanol; ethanol '
        f'only = less than {C.IBO_THRESHOLD:g} g·L⁻¹ isobutanol). The '
        f'table above the axes gives, for each pair, this work / replicate. '
        f'Its first two rows apply to the seeded campaigns only: the IRR of '
        f'the best preloaded seed and the first simulated trial that beat '
        f'it, which compares the two seeded campaigns like for like (the '
        f'replicate\'s seeds already held a design above 25 %).',
        '',
        f'**(b)** The seven replicate scouts, drawn as in panel c of the '
        f'main figure and in its family order (ethanol TRY, then isobutanol '
        f'TRY; the price-weighted-yield campaign, which has no counterpart '
        f'in this work, last): every '
        f'completed trial (dots), the best-IRR trial visited (open marker), '
        f'the design returned (filled marker; argmax of the campaign\'s own '
        f'objective) and their connector. Marker shape gives the objective: '
        f'circle yield, square titer, triangle productivity. The dashed line '
        f'and the tint mark this work\'s uninformed plateau ({U}), kept as '
        f'the reference although the replicate did not stay on it (a). The '
        f'table gives, as in main panel c, the trials above that plateau and '
        f'the share of trials losing money, and then the returned design\'s '
        f'IRR, printed "loss" (as in the main figure and Fig. S2) for every '
        f'returned design drawn in the loss band: {loss_ret}.',
        '',
        '**What replicates.**',
        '',
        f'* Seeded search beat every seed within '
        f'{facts["relay"]["first_gt_best_seed"]} and '
        f'{F["relay_rep_first_gt_best_seed"]} trials (best seeds '
        f'{C.fmt_pct(seeds["best_irr_pct"])} and '
        f'{C.fmt_pct(F["rep_best_seed_pct"])}; the replicate\'s was its '
        f'isobutanol-yield scout\'s best visit, already above 25 %, so its '
        f'{prr["first_ge25"]} trials to 25 % are the same event and not a '
        f'second measure of speed). This work\'s seeded campaign passed '
        f'25 % at trial {pr["first_ge25"]}. The two reached '
        f'{C.fmt_pct(pr["best_pct"])} and {C.fmt_pct(prr["best_pct"])} '
        f'(trials {C.fmt_int(pr["best_sim"])} and '
        f'{C.fmt_int(prr["best_sim"])}). After 25 trials the seeded '
        f'campaigns stood at {C.fmt_pct(pr["at25_pct"])} and '
        f'{C.fmt_pct(prr["at25_pct"])}; the unseeded ones at '
        f'{_irr_words(pu["at25_pct"])} and {_irr_words(pur["at25_pct"])}.',
        f'* The isobutanol scouts visited designs above {U} '
        f'({ibo_n} trials for yield, titer and productivity; price-weighted '
        f'yield {C.fmt_int(rep["rep_pw"]["n_gt_U"])}), while the ethanol '
        f'scouts '
        f'never did (best {C.fmt_pct(F["rep_etoh_max_visit_pct"])}).',
        f'* Every replicate scout returned a design below this work\'s '
        f'uninformed plateau '
        f'(highest: {C.fmt_pct(F["rep_max_returned_pct"])}).',
        f'* All three campaigns that passed 27 % returned co-production '
        f'designs: this work\'s TRY-informed '
        f'{g(best["relay"]["ibo"])} isobutanol + {g(best["relay"]["etoh"])} '
        f'ethanol; replicate TRY-informed {g(best["relay_rep"]["ibo"])} + '
        f'{g(best["relay_rep"]["etoh"])}; replicate uninformed '
        f'{g(best["unin_rep"]["ibo"])} + {g(best["unin_rep"]["etoh"])}.',
        '',
        '**What does not replicate.**',
        '',
        f'* The unseeded plateau: this work\'s uninformed campaign stayed at '
        f'{U} (ethanol only) from trial {facts["U_sim"]} to trial '
        f'{C.fmt_int(facts["campaigns"]["unin"]["max_sim"])}; every '
        f'improvement of the replicate co-produced, from a start-up design '
        f'on, and it passed 25 % at trial {w_25["sim"]}. With two runs this '
        f'is an observation, not '
        f'a test: 1 of 2 unseeded runs stayed on the ethanol-only plateau, '
        f'while seeded search was fast in 2 of 2.',
        '',
        '*Methods notes (not part of the legend).*',
        f'* Every improvement of this work\'s uninformed campaign was an '
        f'ethanol-only design (all {F["n_steps"]["unin"]} best-so-far steps; '
        f'first positive-IRR design at trial '
        f'{F["unin_first_profitable_sim"]}, after its '
        f'{F["main_startup_rows"]}-trial space-filling start-up design). The '
        f'replicate\'s {F["rep_startup_rows"]}-trial start-up design already '
        f'held a co-production design with a positive IRR below the '
        f'{HURDLE_PCT} % hurdle (trial {F["unin_rep_first_profitable_sim"]}: '
        f'{C.fmt_pct(F["unin_rep_first_profitable_pct"])} with '
        f'{g(F["unin_rep_first_profitable_ibo"])} isobutanol + '
        f'{g(F["unin_rep_first_profitable_etoh"])} ethanol), and every later '
        f'improvement kept co-producing (all {F["n_steps"]["unin_rep"]} '
        f'steps): it passed {U} at trial {w_esc["sim"]} '
        f'({C.fmt_pct(w_esc["irr_pct"])}; {g(w_esc["ibo"])} + '
        f'{g(w_esc["etoh"])}), 25 % at trial {w_25["sim"]} '
        f'({C.fmt_pct(w_25["irr_pct"])}; {g(w_25["ibo"])} + '
        f'{g(w_25["etoh"])}) and {C.fmt_pct(pur["best_pct"])} only at trial '
        f'{C.fmt_int(pur["best_sim"])}. In both unseeded runs every '
        f'improvement kept the product class of the campaign\'s first '
        f'positive-IRR design.',
        '* "Trial" is the simulated index, 1-based: CSV trial_number + 1, '
        'or, for a TRY-informed campaign, trial_number minus its preload '
        'size + 1 (preloaded seeds excluded); only completed trials are '
        'used. IRRs are at default prices; a loss '
        '(IRR < 0, or no IRR for an outright money-loser) is drawn in the '
        'grey band. The replicate uninformed and scout campaigns ran before '
        'the model version was pinned; the replicate uninformed campaign\'s '
        f'best design re-simulated to within {REPRO_IRR_TOL_MD} IRR of its '
        'recorded value.',
        '',
    ]
    return '\n'.join(lines)


# %% Main ------------------------------------------------------------------------
def main():
    hits = C.literal_scan([os.path.abspath(__file__)])
    if hits:
        raise C.FactsMismatch(f'hard-coded annotation literals: {hits}')
    facts = C.check_facts()
    F = s1_facts(facts)
    S.cvd_check(pairs=S1_CVD_PAIRS,
                palette={**S.PALETTE, 'pw_dot': PW_DOT}, dataviz=False)

    fig = S.new_figure(*SIZE)
    ax_a, _ = draw_panel_a(fig, facts, F)
    draw_panel_b(fig, facts, F)
    S.panel_letter(fig, LETTER_A_X, LETTER_Y, 'a')
    S.panel_title(fig, LETTER_A_X, LETTER_Y,
                  'Seeding is fast in both runs; the plateau did not replicate')
    S.panel_letter(fig, LETTER_B_X, LETTER_Y, 'b')
    S.panel_title(fig, LETTER_B_X, LETTER_Y,
                  'Replicate scouts again visit higher IRR than they return')

    S.check_figure(fig, size=SIZE)
    clashes = (text_line_clashes(fig, ax_a) + text_marker_clashes(fig, ax_a)
               + text_spine_clashes(fig, ax_a))
    for m in clashes:
        S._safe_print(f'[lines] {m}')
    if clashes:
        raise S.FigureCheckError(f'{len(clashes)} label / line clash(es)')
    out = S.save(fig, STEM)
    with open(CAPTION_PATH, 'w', encoding='utf-8') as fh:
        fh.write(build_caption(facts, F))
    print(f'caption: {CAPTION_PATH}')
    C.assert_sim_safe()
    return out


if __name__ == '__main__':
    main()
