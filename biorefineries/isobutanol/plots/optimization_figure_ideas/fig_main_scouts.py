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
       left (an arrow from the strip as a whole into the TRY-informed
       label; the co-producing seeds above the plateau ringed), the
       best-seed level, the uninformed start-up strip and the refinement
       inset (T2, T4)
    b  where: IRR of every simulated design vs isobutanol share of the
       alcohol titer; the uninformed isobutanol designs outlined (T2, T3,
       T4)
    c  visited vs returned: every trial of each campaign as a strip, the
       best visit (open) and the returned design (filled) (T3)
    d  titers of each returned design (T1, T4)
    e  metabolic proteome of each returned design (T1)
    f  exploration (trials making isobutanol) vs IRR of the best visit and
       the returned design, with the TRY-informed campaign's first trial
       >= 25 % (T5: speed, not only the endpoint; also marked in a)

Round 3: b's two free labels are placed by search (`_place_label`: no key
mark / reference line / label under them, fewest outlined uninformed
designs under the glyph lines, fewest other dots); the TRY families carry
one name set ('ethanol TRY' / 'isobutanol TRY', = C.FAMILY_LABEL), 'scouts'
being the collective noun for the six TRY campaigns.

Round 4: a names its seed columns in full inside the strip, rings the
co-producing seeds above the plateau (the preload already held them), says
what the TRY-informed line produces and puts the inset's linear trial
labels on top; b's labels are drawn as a sandwich (white under-copy below
the outlined uninformed designs, text above them; SANDWICH_Z) and placed
per word; c's title states the finding; f's IRR axis stops at the shared
limit, with every label next to its own marker; the caption states that
the seeds held co-producers and is trimmed (detail moved to the methods
notes).

Round 5: no data marker under any letter. b's labels are placed by an
exhaustive search (`_spot` / `_best_spot`, replacing `_place_label` and the
round-4 sandwich, whose outlined designs ran through the letters) on spots
no marker disk touches (asserted): '27.3 %' beside the star, the
best-seed label above its ring with a leader, the 296-designs finding top
left, the plateau label above the plateau's empty left end with a leader;
the outlined-designs block, for which no spot below the starting strain is
free, is the one label on an opaque white backing box, placed to hide the
fewest outlined designs (never the best one; counts in the caption's
methods notes). a's strip label moves off its two seeds, and the inset
moves into the empty region between IRR 0 and the starting-strain label
(`INSET_DATA`, clearances asserted by `inset_clearance`).

Fresh round 1: the d/e key in three rows ending inside e (its third column
ran under f's letter), the bottom-row axes 0.11 in shorter for it; the
TRY-informed text in its own darker teal (`relay_text`); a's seed columns
named as scout families in a header key (no in-strip labels), the
first trial >= 25 % marked, the start-up strip kept below the lines and the
starting-strain label above its line; b's outlined-designs block on a
translucent backing (dims, does not hide; the one marker-check exemption)
and the number of overplotted isobutanol-only designs beside the 100 %
column; c's '0' labelled and its filled markers above the open ones; f's
seeds arrow ending on the TRY-informed star; the caption's two-sentence
headline, trial / FAIL definitions, the seed-selection arithmetic and a
compute caveat (`facts['sobol']`).

SIM-SAFE: pandas / numpy / matplotlib only, through the sibling modules
_common (data, definitions, facts) and _style (palette, axes, checks);
enzyme_burden.py is read BY FILE PATH for the caption's housekeeping sector.
Every annotation number is built from `_common.compute_facts()` (asserted
against `_common.EXPECTED` by `check_facts()` before drawing); the script also
asserts the figure-level claims it prints (e.g. "ethanol only", "all
co-production"), counts the loss dots against the facts, and checks that no
annotation touches a key marker, connector, trajectory line, arrow or the
inset (`mark_clashes`) and that the inset sits clear of the 0 line and a's
labels (`inset_clearance`), on top of `_style.check_figure` (which includes
`marker_text_hits`: no scatter / marker centre under any annotation unless
an opaque backing box hides it).

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
from matplotlib.patches import FancyArrowPatch               # noqa: E402
from matplotlib.text import Text                             # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator   # noqa: E402

STEM = 'kinBO_main_scouts'
CAPTION_PATH = os.path.join(C.FIGDIR, 'kinBO_main_scouts_caption.md')
EB_PATH = os.path.join(C.PKG, 'enzyme_burden.py')

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
    'b': (6.02, 4.55, 3.86, 2.85),
    # fresh round 1: 0.11 in shorter (2.58 -> 2.47) so the header band
    # holds the three-row d/e key
    # 2026-10-01: the two-column table between c and d (trials above the
    # plateau, trials losing) removed at the user's request; c widened into
    # its space, ending 0.30 in short of d's left spine
    'c': (1.52, 0.72, 2.62, 2.47),
    'd': (4.44, 0.72, 0.72, 2.47),
    'e': (5.26, 0.72, 1.64, 2.47),
    'f': (7.52, 0.72, 2.36, 2.47),
}
# a's inset in a's DATA coordinates ((trials), (IRR %)), round 5: inside
# the empty region above IRR = 0 and below the starting-strain label (its
# frame dipped into the loss band and its top labels crowded the label
# before); the right edge keeps its '1,000' label clear of a's right-spine
# ticks (round 2). The clearances are asserted (inset_clearance)
INSET_DATA = ((125.0, 1400.0), (0.9, 7.6))
INSET_GAP_PT = 2.0                # min gap: inset (labels incl.) <-> the 0
                                  # line, a's labels, a's right-spine ticks
GROUP_X = 0.12                    # group headers, left-aligned
ROW_LABEL_X = 1.30                # row labels, right-aligned
ROW_MARK_X = 1.40                 # role marker after the row label
# row separators: one run across the labels, c, d and e (no table between
# c and d to break around since 2026-10-01)
SEP_SEGMENTS = ((0.12, 6.90),)
TOP_LETTER_Y, BOT_LETTER_Y = 7.66, 3.86
# bottom-row letters: >= 0.3 in between a title's end and the next letter
# (round 3: 'Titers' ran into 'e'); f's title is longer, so its letter sits
# further left with a smaller letter-title gap (TITLE_DX)
# fresh round 2: 'd' sits over d's left spine (at 4.20 in, 'd Titers' read
# as the title of the former c table's 'trials losing' column); e moves
# right to keep MIN_TITLE_GAP_IN, with tighter letter-title gaps (TITLE_DX)
LETTER_X = {'a': 0.06, 'b': 5.42, 'c': 0.06, 'd': 4.40, 'e': 5.36, 'f': 6.94}
TITLE_DX = {'d': 0.20, 'e': 0.22, 'f': 0.22}   # default 0.28 in
MIN_TITLE_GAP_IN = 0.30           # title end -> next letter (asserted)
BAND_Y = (3.645, 3.475)           # the two key lines of the header band
# fresh round 1: the d/e key in three rows (two columns + an inline third
# row), so it ends inside e's column (its third column ran under f's letter)
BAND_Y_DE = BAND_Y + (3.305,)
HEAD_Y = 7.45                     # top-row sub-header line (a strip, b key)
KEY_DE_X = 4.52                   # d/e key start (0.08 in right of d's
                                  # left edge)
KEY_GAP_IN = 0.30                 # min gap between the d/e and f keys
KEY_DE_MAX_X = 6.90               # the d/e key ends inside e (e's right)

# bottom-row rows (shared y of c, the table, d, e; y grows downward)
ROW_Y = {'base': 0.0, 'unin': 1.8, 'relay': 2.8, 'ey': 4.6, 'et': 5.6,
         'ep': 6.6, 'iy': 8.4, 'it': 9.4, 'ip': 10.4}
HEADER_Y = {'profit': 1.0, 'etoh': 3.8, 'ibo': 7.6}
ROW_YLIM = (10.9, -0.5)
BAR_H = 0.62

STARTUP_TOP = S.LOSS_CENTER - 0.75   # a: top of the start-up strip (IRR)

# role-marker sizes (scatter s, pt^2)
S_RING, S_FILLED = 44, 38
S_DIAMOND_C, S_STAR_C, S_BASE_C = 60, 120, 40

# scout dots in the a seed strip and in b: the light family shades darkened
# by 8 L* for DOTS only (they sat near the visibility floor at 180 mm);
# strips in c and every marker keep the palette shades. Fresh round 2: every
# dot cloud of b (and the seed strip) is drawn from S.B_DOTS, whose
# rendered pairs cvd_check verifies under deutan / protan / tritan vision
DOT_ETOH = S.B_DOTS['etoh']['face']
DOT_IBO = S.B_DOTS['ibo']['face']
BEST_SEED_LINE = dict(color=P['ibo_dark'], lw=0.9, ls=(0, (1.2, 1.6)),
                      zorder=1.15)

RIGHT_PAD_PT = 5.0                # a: right-aligned labels, clear of ticks
# b labels (round 5): every label on a spot no marker touches (LABEL_PAD_PT
# clearance), inside the axes clear of the inward ticks (B_MARGINS_PT:
# left, bottom, right, top); the outlined-designs block alone on a tight
# white backing box above every dot layer (BLOCK_Z) and below the key marks
LABEL_PAD_PT = 1.0
LINE_CLEAR_PT = 1.0               # extra clearance from a reference line
LABEL_LS = (1.15, 1.1, 1.0)       # multi-line labels: the most generous
                                  # line spacing a free spot holds
BLOCK_LS = (1.1,)                 # the block's line spacing
INSET_LS = 1.05                   # the inset's two-line note
B_MARGINS_PT = (3.0, 3.0, 5.0, 4.5)
BLOCK_BOX_PAD = 0.25              # backing-box pad, fraction of the font size
# fresh round 1: TRANSLUCENT (alpha 0.7 < S.BACKING_MIN_ALPHA), so the
# designs behind the block stay visible (dimmed to 30 %) instead of hidden
# (round 5's opaque box hid 15 designs, disclosed only in the methods
# notes). The block is the one label exempt from the marker check
# (BLOCK_TEXT -> check_figure(marker_exempt=...)); it is still placed to
# cover the fewest outlined designs, never the best one
# fresh round 2: OPAQUE again (alpha >= S.BACKING_MIN_ALPHA): through the
# translucent box the dimmed designs read as specks inside the text at
# print size. The box is placed to hide the fewest large dots, never the
# best one; what it hides is counted in the methods notes (BLOCK_HIDDEN)
BLOCK_BOX = dict(boxstyle=f'square,pad={BLOCK_BOX_PAD}', fc='white',
                 ec='none', alpha=0.92)
BLOCK_TEXT = []                   # the block's Text (marker-check exempt)
BLOCK_Z = 2.9                     # dots z 2 / outlined 2.6 < box < marks 5-7
BLOCK_MAX_HIDDEN_OUTLINED = 4     # outlined designs the box may hide
BLOCK_HIDDEN = {}                 # what the box dims (printed)
COL_SHARE_MIN = 99.9              # b: the 100 % column (isobutanol share)
COL_JITTER = 2.0                  # fresh round 2: that column's dots spread
                                  # over [100 - COL_JITTER, 100] % (stated in
                                  # the caption), so its two campaigns show
RING_B_S = 42                     # b's isobutanol-TRY best-visit rings (pt^2)
LEADER_LW = 0.7
LEADER_MIN_GAP_PT = 3.0           # a label further above its line gets a
                                  # leader
# in-figure family names (lower case of C.FAMILY_LABEL: 'ethanol TRY')
FAM_TEXT = {f: C.FAMILY_LABEL[f][0].lower() + C.FAMILY_LABEL[f][1:]
            for f in ('etoh', 'ibo')}

ANNOT = []        # annotation Text artists checked against the marks
COLUMN_COUNT = {}  # b: designs in the 100 % column above the plateau
MARKS = []        # key marks / lines / arrows annotations must not touch
A_TEXTS = {}      # a's labels the inset check refers to


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


def _zero_line(ax):
    """The IRR = 0 line S.irr_axis just drew (registered for the mark
    check; round 3: f's note sat on it)."""
    ln = ax.lines[-1]
    assert np.allclose(ln.get_ydata(), 0.0), 'not the zero line'
    return ln


def _dots(ax, x, y, spec, z=2.0):
    """A dot cloud in the S.B_DOTS style `spec` (rasterized)."""
    kw = (dict(edgecolors=spec['edge'], linewidths=spec['lw'])
          if spec.get('edge') else dict(lw=0))
    return ax.scatter(x, y, s=spec['s'], c=spec['face'], alpha=spec['alpha'],
                      rasterized=True, zorder=z, **kw)


def _key_dot(name, ms=5.0):
    """An inline_key marker item showing the B_DOTS cloud `name` as it
    renders (alpha pre-blended on white; the rim scaled to the key dot so
    it covers the same share of the disc)."""
    spec = S.B_DOTS[name]
    a = spec['alpha']
    item = dict(marker='o', ms=ms, mfc=S.blend_white(spec['face'], a),
                mew=0.0)
    if spec.get('edge'):
        r = np.sqrt(spec['s']) / 2.0
        frac = 1.0 - ((r - spec['lw'] / 2) / (r + spec['lw'] / 2)) ** 2
        q = np.sqrt(1.0 - frac)                  # ri / ro on the key dot
        h = ms / 2.0 * (1.0 - q) / (1.0 + q)     # half the key rim width
        item.update(mec=S.blend_white(spec['edge'], a), mew=2.0 * h)
    else:
        item['mec'] = item['mfc']
    return item


def _first_coprod_run(rc):
    """The first simulated trial from which EVERY best-so-far improvement
    of the (sim-ordered, COMPLETE) TRY-informed rows co-produces."""
    ri = rc['IRR'].to_numpy(float)
    step = np.r_[True, ri[1:] > np.maximum.accumulate(ri)[:-1]]
    st = rc[step]
    cls = st['product_class'].to_numpy()
    sims = st['sim'].to_numpy(int)
    k = len(cls)
    while k > 0 and cls[k - 1] == 'coprod':
        k -= 1
    assert k < len(cls), 'the last improvement does not co-produce'
    return int(sims[k])


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
    # c: every scout's best visit is above its returned design; c title
    # (round 4): the isobutanol scouts visit above the plateau ('high IRR')
    # and return below the starting strain or no IRR ('low')
    for k in C.SCOUT_KEYS:
        st = camp[k]
        assert st['best_visit']['irr_pct'] > st['returned']['irr_pct'], k
    for k in C.IBO_SCOUTS:
        assert camp[k]['n_gt_U'] >= 1, k
        assert camp[k]['returned']['irr_pct'] < F['start_irr_pct'], k
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
    # round 2 claims: the title's 'passes every seed by trial 13' is the
    # FIRST trial above the best seed (and it is a newly simulated design)
    s13 = rel['sims'][rel['first_gt_best_seed']]
    assert s13['irr_pct'] == s13['bsf_pct'] > F['seeds']['best_irr_pct']
    assert rel['sims'][2]['irr_pct'] < F['seeds']['best_irr_pct']
    # round 4, a: 'ethanol + isobutanol' under the TRY-informed line -- every
    # best-so-far step from the first trial above every seed on co-produces
    rc = C.complete('relay').sort_values('sim')
    ri = rc['IRR'].to_numpy(float)
    step = np.r_[True, ri[1:] > np.maximum.accumulate(ri)[:-1]]
    steps = rc[step & (rc['sim'].to_numpy(int) >= rel['first_gt_best_seed'])]
    assert len(steps) and (steps['product_class'] == 'coprod').all()
    assert 25 in rel['sims'] and rel['sims'][25]['class'] == 'coprod'
    # fresh round 2: the line label says 'from trial N': every improvement
    # from N on co-produces, the incumbent just before N does not (trial 3's
    # isobutanol-only design), and N is a pinned fact
    n_cop = _first_coprod_run(rc)
    assert n_cop in rel['sims'] and rel['sims'][n_cop]['class'] == 'coprod'
    assert rel['sims'][3]['class'] == 'ibo' and 3 < n_cop
    assert n_cop <= rel['first_gt_best_seed']
    # round 4, caption: the preload held co-producers above the plateau,
    # all below the best seed, and the first trial above every seed beats
    # them all
    sd_ = F['seeds']
    assert 0 < sd_['n_coprod_gt_U'] <= sd_['n_gt_U']
    assert U < sd_['best_coprod_irr_pct'] < sd_['best_irr_pct']
    assert rel['sims'][2]['irr_pct'] < sd_['best_coprod_irr_pct']
    assert camp['ip']['best_visit']['class'] == 'coprod'
    assert abs(camp['ip']['best_visit']['irr_pct']
               - sd_['best_coprod_irr_pct']) < 1e-9
    # the uninformed isobutanol designs: start-up + post-plateau = all
    ib = F['unin']['ibo_ge5']
    assert ib['n_startup'] + ib['n_after_plateau'] == ib['n']
    # e tags: every design beyond the budget is tagged (and only those)
    for k in C.ROW_KEYS:
        assert (pr[k]['Phi_M'] > F['budget']) == (pr[k]['growth'] < 1.0), k
    # round 3 caption claims: the scouts' visits above the plateau split
    # into pure-isobutanol + co-production; the title's 'visit designs more
    # profitable than any the uninformed campaign reached' is about
    # isobutanol TRY campaigns only
    g = F['global']
    assert (g['n_try_gt_U_ibo_only'] + g['n_try_gt_U_coprod']
            == g['n_try_gt_U'] == sum(camp[k]['n_gt_U']
                                      for k in C.IBO_SCOUTS))
    assert camp['unin']['n_gt_U'] == 0
    # f / caption: the uninformed campaign never reached 25 %; the
    # TRY-informed campaign did (its trial is f's speed label)
    prg = F['progress']
    assert prg['unin']['first_ge25'] is None
    assert prg['relay']['first_ge25'] == rel['first_ge25'] is not None


# %% Panel a: when ------------------------------------------------------------------
def panel_a(fig, F):
    U, start = F['U_pct'], F['start_irr_pct']
    rel, unin = F['relay'], F['unin']
    best_seed = F['seeds']['best_irr_pct']
    ax = S.inch_axes(fig, *L['a'])
    strip = S.inch_axes(fig, *L['a_strip'], sharey=ax)

    # --- seed strip (shares y with a; carries the y title and labels)
    S.irr_axis(strip, 'y', step=5)
    S.plateau_line(strip, U)
    strip.axhline(best_seed, **BEST_SEED_LINE)
    strip.set_xlim(0, 1)
    man = C.load_manifest('relay')
    rng = S.panel_rng()
    best_idx = man['IRR'].idxmax()
    n_loss = 0
    seed_xy = {}
    for fam, (lo, hi), col in (('etoh', (0.08, 0.42), DOT_ETOH),
                               ('ibo', (0.58, 0.92), DOT_IBO)):
        m = man[man['family'] == fam]
        x = rng.uniform(lo, hi, len(m))
        y = S.irr_plot(m['IRR'].to_numpy(float), rng)
        x[m.index.to_numpy() == best_idx] = 0.75      # ring sits on its dot
        sc = _dots(strip, x, y, S.B_DOTS[fam])
        n_loss += int((sc.get_offsets()[:, 1] < 0).sum())
        seed_xy.update(zip(m.index, zip(x, y)))
    assert n_loss == F['seeds']['n_losing'], 'seed loss dots != facts'
    assert man.loc[best_idx, 'family'] == 'ibo'
    _ring(strip, 0.75, best_seed, 'o', P['ibo_dark'], s=40, lw=1.2)
    # round 4: the co-producing seeds above the plateau ringed in the
    # TRY-informed teal (the preload already held profitable co-producers;
    # caption (a)), so trial 2's co-producer is not read as new by itself
    cop = man.index[man['coprod'].to_numpy(bool)
                    & C.above_plateau(man['IRR'], F['U'])]
    assert len(cop) == F['seeds']['n_coprod_gt_U']
    assert best_idx not in cop
    cx, cy = np.array([seed_xy[i] for i in cop]).T
    assert np.isclose(cy.max(), F['seeds']['best_coprod_irr_pct'])
    strip.scatter(cx, cy, s=20, facecolors='none', edgecolors=P['relay'],
                  linewidths=0.9, zorder=3)
    S.style_ticks(strip, minor_x=False)
    strip.set_xticks([0.25, 0.75])
    strip.tick_params(axis='x', which='both', length=0, top=False,
                      labelbottom=False)
    assert F['seeds']['max_etoh_irr_pct'] < F['U_pct']
    strip.set_ylabel(S.bold_axis_title('IRR [%]'), labelpad=2)
    # header: the seeds, the scout trials they were picked from (so the
    # trial axis is not read as the whole compute) and, fresh round 1, the
    # two columns named as SCOUT FAMILIES in a colour key ('ethanol' /
    # 'isobutanol' inside the strip read as product classes, and the
    # isobutanol label sat on seeds: no free stretch of that column holds
    # 'isobutanol TRY')
    hx = L['a_strip'][0] - 0.08
    h1 = f"{C.fmt_int(F['seeds']['n'])} seeds"
    S.fig_text(fig, hx, HEAD_Y, h1, fontsize=S.FS['key'], fontweight='bold',
               ha='left', va='bottom')
    h2 = f" from {C.fmt_int(F['global']['n_try'])} scout trials:"
    x2 = hx + _text_w_in(fig, h1, S.FS['key'], fontweight='bold')
    t2 = S.fig_text(fig, x2, HEAD_Y, h2, fontsize=S.FS['key'],
                    color=GREY_LABEL, ha='left', va='bottom')
    fig.canvas.draw()
    bb = t2.get_window_extent(fig.canvas.get_renderer())
    yc = (bb.y0 + bb.y1) / 2 / fig.dpi
    xk, _ = S.inline_key(
        fig, x2 + _text_w_in(fig, h2, S.FS['key']) + 0.08, yc,
        [{**_key_dot('etoh'), 'text': FAM_TEXT['etoh']},
         {**_key_dot('ibo'), 'text': FAM_TEXT['ibo']}],
        fontsize=S.FS['key'], item_gap_in=0.13)
    assert xk <= LETTER_X['b'] - 0.25, f'a header runs to {xk:.2f} in'

    # --- main axes
    S.log_trial_axis(ax)
    S.loss_band(ax)
    # the UNINFORMED campaign's space-filling start-up: a cyan strip along
    # its flat loss segment only (the TRY-informed campaign had none).
    # Fresh round 1: the strip stops below the two lines' loss level
    # (STARTUP_TOP), so it holds the label and not the TRY-informed line's
    # trial-1 loss
    n_start = F['checks']['shared_startup_rows']
    ax.fill_between([1, n_start + 1], S.LOSS_BAND[0], STARTUP_TOP,
                    color=P['unin'], alpha=0.10, lw=0, zorder=0.12)
    MARKS.append(S.plateau_line(ax, U))
    MARKS.append(S.start_line(ax, start))
    MARKS.append(ax.axhline(best_seed, **BEST_SEED_LINE))
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

    # end markers and milestones: trial 2 (first above the plateau), the
    # first trial above every seed, the best
    rb = rel['best']
    sim2, sim_s = rel['first_gt_U'], rel['first_gt_best_seed']
    sim25 = rel['first_ge25']
    _mark(ax, F['U_sim'], U, 'D', 55, P['unin'], 'white', 0.8, z=6)
    for sim in (sim2, sim_s, sim25, rb['sim']):
        v = rel['sims'][sim]['bsf_pct'] if sim in rel['sims'] else \
            rel['bsf_pct'][sim]
        _mark(ax, sim, v, 'o', 22, P['relay'], 'white', 0.6, z=6)
    _mark(ax, rb['sim'], rb['irr_pct'], '*', 140, P['relay'], 'white', 0.8,
          z=7)
    v2 = rel['sims'][sim2]['bsf_pct']
    fs = S.FS['note']
    _ann(ax, f"{C.fmt_irr(v2)} · {C.fmt_trial(sim2)}", (sim2, v2), (5, -3),
         ha='left', va='top', fontsize=fs, color=P['relay_text'])
    # fresh round 1: the first trial >= 25 % (also quoted in f). Fresh
    # round 2: the trial-13 point keeps its dot
    # but loses its label (the title states it; one of the labels the
    # density review asked to drop), and this label, above-left of its dot,
    # gets a leader to it (it sat over the trial-13 dot, unattached)
    v25 = rel['sims'][sim25]['bsf_pct']
    assert v25 >= 25.0 > rel['sims'][sim_s]['bsf_pct']
    t25 = _ann(ax, f"{C.fmt_irr(v25)} · {C.fmt_trial(sim25)}", (sim25, v25),
               (-9, 8.5), ha='right', va='bottom', fontsize=fs,
               color=P['relay_text'],
               arrowprops=dict(arrowstyle='-', color=P['relay_text'],
                               lw=LEADER_LW, relpos=(1.0, 0.0), shrinkA=0.5,
                               shrinkB=3.0))
    MARKS.append(t25.arrow_patch)
    _ann(ax, C.fmt_irr(rb['irr_pct']), (rb['sim'], rb['irr_pct']), (0, 7),
         ha='center', va='bottom', fontsize=fs, color=P['relay_text'])
    # round 4: what the TRY-informed line produces, as the uninformed line
    # says 'ethanol only', under the line from trial 25. Fresh round 2:
    # 'from trial N' (trials 3-6 held an isobutanol-only incumbent; every
    # improvement from N on co-produces: _first_coprod_run, asserted)
    s25 = 25
    n_cop = _first_coprod_run(C.complete('relay').sort_values('sim'))
    _ann(ax, f'ethanol + isobutanol from {C.fmt_trial(n_cop)}',
         (s25, rel['sims'][s25]['bsf_pct']), (4, -4), ha='left', va='top',
         fontsize=fs, color=P['relay_text'])

    # start-up label (in the loss band, under the flat cyan segment),
    # starting strain, best seed, direct labels
    _ann(ax, 'uninformed space-filling start-up',
         (np.sqrt(n_start + 1), S.LOSS_BAND[0] + 1.15), ha='center',
         va='center', fontsize=fs, color=P['unin_text'])
    # the right-aligned labels end RIGHT_PAD_PT short of x = 2,000 (round 3:
    # at 0 pt 'trial 104,' and 'trials' touched the right-spine ticks)
    rp = -RIGHT_PAD_PT
    # fresh round 1: above its line, like the best-seed and plateau labels
    # (below it, the inset's top tick labels sat tight under it)
    A_TEXTS['start'] = _ann(
        ax, f'starting strain {C.fmt_irr(start)}', (2000, start), (rp, 3),
         ha='right', va='bottom', fontsize=fs, color=NOTE)
    # fresh round 2: the best-seed label at the LEFT end of its line, where
    # the line leaves the seed strip (on the right it took the room the
    # 'ethanol + isobutanol from trial N' label now needs)
    _ann(ax, f'best seed {C.fmt_irr(best_seed)}', (1.0, best_seed), (4, 3),
         ha='left', va='bottom', fontsize=fs, color=P['ibo_dark'])
    # fresh round 2: the lines are BEST-SO-FAR curves (the strip's dots are
    # single seeds): said in the empty right end of the loss band
    _ann(ax, 'lines: best IRR so far', (2000, S.LOSS_CENTER), (rp, 0),
         ha='right', va='center', fontsize=fs, fontstyle='italic',
         color=NOTE)
    _ann(ax, f"uninformed: {C.fmt_irr(U)} at {C.fmt_trial(unin['plateau_sim'])},"
             f"\nethanol only; flat for {C.fmt_int(unin['flat_trials'])} trials",
         (2000, U), (rp, 7), ha='right', va='bottom', color=P['unin_text'],
         linespacing=1.15)
    # 'TRY-informed: search starts from these seeds', fed by an arrow out of
    # the seed strip as a whole (from above every seed, not from one seed)
    # fresh round 2: 'TRY-informed (seeded)' defines the name once (the
    # titles say 'seeded search')
    y_lab, x_lab = 29.4, 1.5
    t1 = _ann(ax, 'TRY-informed', (x_lab, y_lab), ha='left', va='center',
              fontweight='bold', color=P['relay_text'])
    fig.canvas.draw()
    bb = t1.get_window_extent(fig.canvas.get_renderer())
    x_end = ax.transData.inverted().transform((bb.x1, bb.y0))[0]
    _ann(ax, ' (seeded): starts from these seeds', (x_end, y_lab), (0, 0),
         ha='left', va='center', fontsize=fs, color=GREY_LABEL)
    # a figure-level arrow in inches (the layout is fixed in inches): from
    # the middle of the strip to the label
    to_in = fig.dpi_scale_trans.inverted()
    pA = to_in.transform(strip.get_yaxis_transform().transform((0.5, y_lab)))
    pB = to_in.transform(ax.transData.transform((x_lab, y_lab)))
    con = FancyArrowPatch(tuple(pA), tuple(pB), transform=fig.dpi_scale_trans,
                          arrowstyle='-|>', color=GREY_ARROW, lw=0.9,
                          mutation_scale=8, shrinkA=0, shrinkB=2.5, zorder=8)
    fig.add_artist(con)
    MARKS.append(con)

    # inset: refinement of the TRY-informed best so far
    (ix0, ix1), (iy0, iy1) = INSET_DATA
    f0, f1 = ax.transAxes.inverted().transform(
        ax.transData.transform([(ix0, iy0), (ix1, iy1)]))
    ins = ax.inset_axes([f0[0], f0[1], f1[0] - f0[0], f1[1] - f0[1]])
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
    # round 4: the inset's (linear) trial labels on top, away from a's own
    # log trial axis
    ins.tick_params(axis='x', which='both', labeltop=True, labelbottom=False)
    # bottom right, under the flat line (top left touched the trial-871 dot
    # once the inset was narrowed). In IRR, the inset's own axis (round 3:
    # the PI values, a second metric, live in caption (a)). Round 5: two
    # lines, right of the trial 25-100 steps (one line now spans the
    # shorter inset's width and ran into the trial-25 dot)
    v25 = rel['sims'][25]['bsf_pct']
    ins.text(0.95, 0.07, f"refined:\n{C.fmt_num(v25, 1)} → "
                         f"{C.fmt_irr(rb['irr_pct'])}",
             transform=ins.transAxes, ha='right', va='bottom',
             fontsize=S.FS['note'], color=P['relay_text'],
             linespacing=INSET_LS)
    S.style_ticks(ax)
    ax.tick_params(axis='y', which='both', labelleft=False)
    return ax, strip, ins


# %% Panel b: where -----------------------------------------------------------------
def panel_b(fig, F):
    U, start = F['U_pct'], F['start_irr_pct']
    camp = F['campaigns']
    best_seed = F['seeds']['best_irr_pct']
    ax = S.inch_axes(fig, *L['b'])
    S.irr_axis(ax, 'y', step=10)
    MARKS.append(_zero_line(ax))
    ax.set_xlim(-3, 103)
    ax.xaxis.set_major_locator(FixedLocator([0, 25, 50, 75, 100]))
    ax.set_xlabel(S.bold_axis_title('Isobutanol share of alcohol titer [%]'))
    S.tint_above(ax, U)
    MARKS.append(S.plateau_line(ax, U))
    MARKS.append(S.start_line(ax, start))
    # the best seed (= the isobutanol scouts' best visit): lets the reader
    # check '296 designs above 22.9 %' by eye
    MARKS.append(ax.axhline(best_seed, **BEST_SEED_LINE))
    rng = S.panel_rng()
    n_drawn = n_loss = n_loss_exp = 0
    unin_ibo = None
    col_n = {}
    # fresh round 2: the clouds in the S.B_DOTS styles (cyan and the scout
    # violet separated by lightness, not hue: they collapsed under deutan /
    # protan vision); the 100 % column spread over COL_JITTER points
    for name, keys in (('etoh', C.ETOH_SCOUTS), ('ibo', C.IBO_SCOUTS),
                       ('relay', ('relay',)), ('unin', ('unin',))):
        d = pd.concat([C.complete(k) for k in keys], ignore_index=True)
        d = d[d['alcohol'] >= C.SHARE_MIN_ALCOHOL]
        y = S.irr_plot(d['IRR'].to_numpy(float), rng)
        x = d['ibo_share_pct'].to_numpy(float).copy()
        incol = x >= COL_SHARE_MIN
        x[incol] = 100.0 - rng.uniform(0.0, COL_JITTER, int(incol.sum()))
        col_n[name] = int((incol & np.asarray(
            C.above_plateau(d['IRR'], F['U']), bool)).sum())
        scs = []
        if name == 'unin':
            # the uninformed campaign's isobutanol-making designs: larger,
            # drawn on top (T2: it tried them)
            ib = d['makes_ibo'].to_numpy(bool)
            scs.append(_dots(ax, x[~ib], y[~ib], S.B_DOTS['unin']))
            sc2 = _dots(ax, x[ib], y[ib], S.B_DOTS['unin_ibo'], z=2.6)
            scs.append(sc2)
            unin_ibo = int(ib.sum())
        elif name == 'ibo':
            # fresh round 2: the column's scout dots ABOVE the TRY-informed
            # dots (T3: the scout designs above the plateau sit there)
            scs.append(_dots(ax, x[~incol], y[~incol], S.B_DOTS['ibo']))
            scs.append(_dots(ax, x[incol], y[incol], S.B_DOTS['ibo'],
                             z=2.2))
        else:
            scs.append(_dots(ax, x, y, S.B_DOTS[name]))
        n_drawn += len(d)
        n_loss += sum(int((c.get_offsets()[:, 1] < 0).sum()) for c in scs)
        n_loss_exp += int(d['loss'].sum())
    assert n_drawn == F['global']['n_alcohol_ge1'], n_drawn
    assert n_loss == n_loss_exp
    assert unin_ibo == F['unin']['ibo_ge5']['n'], unin_ibo

    # marks: isobutanol scouts' best visits, the two profitability bests
    for k in C.IBO_SCOUTS:
        bv = camp[k]['best_visit']
        _ring(ax, bv['share_pct'], bv['irr_pct'], C.campaign(k).marker,
              P['ibo_dark'], s=RING_B_S, lw=1.3, z=6)
    rb = F['relay']['best']
    ub = camp['unin']['best_visit']
    _mark(ax, rb['share_pct'], rb['irr_pct'], '*', 140, P['relay'], 'white',
          0.8, z=7)
    _mark(ax, ub['share_pct'], ub['irr_pct'], 'D', 50, P['unin'], 'white',
          0.6, z=7)

    # labels (round 5): no data marker under any letter. Every label but
    # the outlined-designs block sits where no marker disk touches its text
    # box (_spot; asserted): '27.3 %' beside the star, the TRY-informed
    # finding top left, the best-seed label top right with a leader down to
    # its ring, the plateau label over the plateau's empty left end with a
    # leader down to the line. No position between the 0 and
    # starting-strain lines is free for the block at 9 pt (dots everywhere
    # over 0-100 %), so it is the one label on a tight white backing box
    # (BLOCK_BOX, above every dot layer, below the key marks), placed where
    # the box hides the fewest outlined designs, then the fewest dots; the
    # counts go to the caption's methods notes (BLOCK_HIDDEN)
    strict = [m for m in MARKS if isinstance(m, PathCollection)
              and m.axes is ax]
    data = [c for c in ax.collections if isinstance(c, PathCollection)
            and c not in strict]
    pts, outl, strict_pts = (_disks(fig, data), _disks(fig, [sc2]),
                             _disks(fig, strict))
    common = dict(pts=pts, strict_pts=strict_pts,
                  hlines=[0.0, start, U, best_seed], margins_pt=B_MARGINS_PT)
    fs = S.FS['note']
    pt = fig.dpi / 72.0
    to_px = ax.transData.transform
    abb = ax.bbox
    placed = []                       # boxes (px) later labels keep clear of

    # '27.3 %' beside the star: right of it, else above it
    sx, sy = to_px((rb['share_pct'], rb['irr_pct']))
    o = np.array([(dx, dy) for dx in np.arange(6.0, 20.01, 0.5)
                  for dy in np.arange(-8.0, 10.01, 0.5)])
    sp_r = _spot(ax, C.fmt_irr(rb['irr_pct']),
                 dict(fontsize=fs, color=P['relay_text'], ha='left',
                      va='center'),
                 np.c_[sx + o[:, 0] * pt, sy + o[:, 1] * pt],
                 cost=np.hypot(o[:, 0] - 6.0, o[:, 1]) * pt, **common)
    o = np.array([(dx, dy) for dx in np.arange(-6.0, 6.01, 0.5)
                  for dy in np.arange(7.0, 16.01, 0.5)])
    sp_a = _spot(ax, C.fmt_irr(rb['irr_pct']),
                 dict(fontsize=fs, color=P['relay_text'], ha='center',
                      va='bottom'),
                 np.c_[sx + o[:, 0] * pt, sy + o[:, 1] * pt],
                 cost=(np.hypot(o[:, 0], o[:, 1] - 7.0) + 6.0) * pt, **common)
    sp = min(sp_r, sp_a, key=lambda q: q['score'])
    _put(ax, sp, placed)

    # fresh round 1: the isobutanol-only scout designs above the plateau
    # overplot into one stroke at 100 % share. Fresh round 2: the column is
    # jittered and its scout dots drawn above the TRY-informed ones, so they
    # show; their count (and the TRY-informed designs' in the same column)
    # moves to the caption -- the bare violet number was unreadable, and no
    # free spot near the column holds a label with a noun (a leader from
    # the nearest free spot ran 80 pt across the ring's label)
    g = F['global']
    d = pd.concat([C.complete(k) for k in C.SCOUT_KEYS], ignore_index=True)
    a = d[C.above_plateau(d['IRR'], F['U'])]
    col = a[a['ibo_share_pct'] >= COL_SHARE_MIN]
    n_col = len(col)
    assert n_col and (col['etoh'] < C.ETOH_THRESHOLD).all()
    assert n_col <= g['n_try_gt_U_ibo_only']
    assert col_n['ibo'] == n_col and col_n['etoh'] == col_n['unin'] == 0
    COLUMN_COUNT.update(n=n_col, n_ibo_only=g['n_try_gt_U_ibo_only'],
                        n_relay=col_n['relay'])

    # the ring (isobutanol TRY best visit = best seed): the lowest free
    # spot above it that spans its x, with a vertical leader down to it
    # (placed before the TRY-informed finding, which can go anywhere)
    iy_bv = camp['iy']['best_visit']
    rx, ry = to_px((iy_bv['share_pct'], iy_bv['irr_pct']))
    xs = np.arange(rx - 1.0, abb.x1 + 0.5, 1.0)
    ys = np.arange(abb.y1, ry, -1.0)
    A = np.array([(x, y) for x in xs for y in ys])
    sp = _best_spot(ax, [f"{FAM_TEXT['ibo']} best\n= best seed, "
                         f"{C.fmt_irr(best_seed)}"],
                    dict(fontsize=fs, color=P['ibo_dark'], ha='right',
                         va='top'), A,
                    lambda box: (box[:, 1] - ry) + 0.3 * (box[:, 2] - rx),
                    within_x=rx, avoid=placed, **common)
    _put(ax, sp, placed, leader_to=(rx, ry),
         shrink_pt=np.sqrt(RING_B_S) / 2 + 1.5, color=P['ibo_dark'],
         pts_all=pts)

    # the TRY-informed finding: the free spot nearest the top-left corner,
    # under the top ticks; two lines if a free spot holds them, else four
    # short ones (the cloud's top reaches the two-line box's bottom edge)
    g = F['global']
    n296, bs = g['n_relay_gt_best_seed'], C.fmt_irr(best_seed)
    xs = np.arange(abb.x0, abb.x0 + 0.6 * abb.width, 1.0)
    ys = np.arange(abb.y1, to_px((0, best_seed))[1], -1.0)
    A = np.array([(x, y) for x in xs for y in ys])
    sp = _best_spot(ax, [f"{n296} designs above {bs}:\n"
                         f"all TRY-informed co-production",
                         f"{n296} designs\nabove {bs}:\nall TRY-informed"
                         f"\nco-production"],
                    dict(fontsize=fs, color=P['relay_text'], ha='left',
                         va='top'),
                    A, (abb.y1 - A[:, 1]) + 0.3 * (A[:, 0] - abb.x0),
                    avoid=placed, **common)
    _put(ax, sp, placed)

    # the plateau label: above the dashed line, as close as an empty spot
    # allows (one or two lines); a vertical leader down to the line when it
    # sits more than LEADER_MIN_GAP_PT above it
    _, uy = to_px((0, U))
    xs = np.arange(abb.x0, abb.x0 + 0.75 * abb.width, 1.5)
    ys = np.arange(uy + 2.0 * pt, to_px((0, best_seed))[1], 1.0)
    A = np.array([(x, y) for x in xs for y in ys])
    # (cost = the gap to the line: the nearest free spot wins, then the
    # preferred wrap / spacing)
    sp = _best_spot(ax, [f'uninformed plateau {C.fmt_irr(U)}',
                         f'uninformed\nplateau {C.fmt_irr(U)}'],
                    dict(fontsize=fs, color=P['unin_text'], ha='left',
                         va='bottom'), A,
                    lambda box: (box[:, 1] - uy) + 0.2 * (box[:, 0] - abb.x0),
                    steps=(1e-3, 1e-2), avoid=placed, **common)
    lead = sp['box'].y0 - uy > LEADER_MIN_GAP_PT * pt
    _put(ax, sp, placed, color=P['unin_text'],
         leader_to=(None, uy) if lead else None,
         pts_all=np.vstack([pts, strict_pts]))

    # the outlined uninformed designs' block, on its backing box: between
    # the 0 and starting-strain lines, wrapped three ways
    ui = F['unin']['ibo_ge5']
    lt = F['unin']['ibo_lt5']
    thr = f'{C.IBO_THRESHOLD:g} {S.UNIT_TITER}'
    n_ui, best_ui = ui['n'], C.fmt_irr(ui['max_irr_pct'])
    lose_ui, lose_lt = _pct0(ui['pct_losing']), _pct0(lt['pct_losing'])
    # fresh round 2: '(large)' -- every uninformed dot now has a dark rim;
    # these are the large ones; narrower wraps let the box reach the empty
    # stretch between 40 and 70 %
    wraps = (f"{n_ui} uninformed designs with ≥ {thr}\nisobutanol "
             f"(large): best {best_ui},\n{lose_ui} lose money (vs "
             f"{lose_lt} of the rest)",
             f"{n_ui} uninformed designs with\n≥ {thr} isobutanol "
             f"(large):\nbest {best_ui}, {lose_ui} lose money\n(vs "
             f"{lose_lt} of the rest)",
             f"{n_ui} uninformed designs\nwith ≥ {thr} isobutanol\n"
             f"(large): best {best_ui},\n{lose_ui} lose money\n(vs "
             f"{lose_lt} of the rest)",
             f"{n_ui} uninformed\ndesigns with\n≥ {thr} isobutanol\n"
             f"(large): best {best_ui},\n{lose_ui} lose money\n(vs "
             f"{lose_lt} of the rest)",
             f"{n_ui} uninformed designs\nmaking isobutanol (large):\n"
             f"best {best_ui}, {lose_ui} lose\nmoney (vs {lose_lt} of the "
             f"rest)",
             f"{n_ui} uninformed designs making\nisobutanol (large): best "
             f"{best_ui},\n{lose_ui} lose money (vs {lose_lt} of the rest)",
             f"{n_ui} uninformed\ndesigns making\nisobutanol (large):\n"
             f"best {best_ui}, {lose_ui}\nlose money (vs\n{lose_lt} of the "
             f"rest)")
    y_lo, y_hi = to_px((0, 0.0))[1], to_px((0, start))[1]
    xs = np.arange(abb.x0, abb.x1, 1.5)
    ys = np.arange(y_lo, y_hi, 1.0)
    A = np.array([(x, y) for x in xs for y in ys])
    blk_kw = dict(fontsize=fs, color=P['unin_text'], ha='center',
                  va='center', bbox=BLOCK_BOX, zorder=BLOCK_Z)
    # fresh round 2: never over the dense 0 % / 100 % columns (a narrow
    # wrap there hid fewer large dots but hundreds of others)
    cols_px = [mtransforms.Bbox([[to_px((xa, 0))[0], abb.y0],
                                 [to_px((xb, 0))[0], abb.y1]])
               for xa, xb in ((-3.0, 4.0), (100.0 - COL_JITTER - 2.0, 103.0))]
    sp = _best_spot(ax, wraps, blk_kw, A,
                    np.abs(A[:, 0] - 0.5 * (abb.x0 + abb.x1)) * 0.01,
                    lss=BLOCK_LS, key_pts=outl,
                    grow_px=BLOCK_BOX_PAD * fs * pt,
                    avoid=placed + cols_px, **common)
    assert sp['score'][0] == 0, f"block: {sp['score']}"
    BLOCK_TEXT.append(_put(ax, sp, placed, boxed=True))
    hb = sp['box'].padded(0.0)
    BLOCK_HIDDEN.update(
        outlined=int(_n_touch(outl, hb.extents[None], 0.0)[0]),
        dots=int(_n_touch(pts, hb.extents[None], 0.0)[0]),
        alpha=BLOCK_BOX['alpha'])
    BLOCK_HIDDEN['other'] = BLOCK_HIDDEN['dots'] - BLOCK_HIDDEN['outlined']
    assert BLOCK_HIDDEN['outlined'] <= BLOCK_MAX_HIDDEN_OUTLINED, BLOCK_HIDDEN
    # the design the block quotes ('best ...') stays visible
    top = outl[[int(np.argmax(outl[:, 1]))]]
    assert _n_touch(top, hb.extents[None], LABEL_PAD_PT * pt)[0] == 0,         'the block hides the best outlined design'
    print(f"b block: opaque box (alpha {BLOCK_BOX['alpha']}) over "
          f"{BLOCK_HIDDEN['outlined']} of {n_ui} large dots + "
          f"{BLOCK_HIDDEN['other']} other dots (score {sp['score'][:4]})")

    S.style_ticks(ax)
    return ax


# %% Label placement (round 5) --------------------------------------------------
def _disks(fig, cols):
    """(N x 3) [x, y, radius] display px of the finite points of the
    scatter collections `cols`."""
    out = [np.zeros((0, 3))]
    for c in cols:
        p = np.asarray(c.get_offset_transform().transform(c.get_offsets()),
                       float).reshape(-1, 2)
        p = p[np.all(np.isfinite(p), axis=1)]
        r = float(np.sqrt(np.max(c.get_sizes())) / 2 * fig.dpi / 72)
        out.append(np.column_stack([p, np.full(len(p), r)]))
    return np.vstack(out)


def _n_touch(pts, boxes, pad, chunk=256):
    """(M,) number of disks of `pts` (N x 3 [x, y, r] px) touching each box
    of `boxes` (M x 4 [x0, y0, x1, y1] px) grown by `pad` px."""
    boxes = np.atleast_2d(np.asarray(boxes, float))
    out = np.zeros(len(boxes), int)
    if not len(pts) or not len(boxes):
        return out
    e = pts[:, 2] + pad
    lo, hi = boxes.min(0), boxes.max(0)
    keep = ((pts[:, 0] > lo[0] - e) & (pts[:, 0] < hi[2] + e)
            & (pts[:, 1] > lo[1] - e) & (pts[:, 1] < hi[3] + e))
    x, y, e = pts[keep, 0], pts[keep, 1], e[keep]
    for i in range(0, len(boxes), chunk):
        b = boxes[i:i + chunk]
        hit = ((x > b[:, :1] - e) & (x < b[:, 2:3] + e)
               & (y > b[:, 1:2] - e) & (y < b[:, 3:4] + e))
        out[i:i + chunk] = hit.sum(1)
    return out


def _extent(ax, s, kw):
    """Window extent of label `s` (annotate kwargs `kw`, without any bbox
    patch) relative to its anchor: [dx0, dy0, dx1, dy1] px."""
    r = ax.figure.canvas.get_renderer()
    probe = ax.annotate(s, xy=(0.5, 0.5), xycoords='axes fraction',
                        annotation_clip=False, **kw)
    bb = probe.get_window_extent(r)
    probe.remove()
    a = ax.transAxes.transform((0.5, 0.5))
    return np.array([bb.x0 - a[0], bb.y0 - a[1], bb.x1 - a[0], bb.y1 - a[1]])


def _spot(ax, s, kw, anchors, *, pts, strict_pts, hlines=(), avoid=(),
          key_pts=None, grow_px=0.0, cost=None, within_x=None,
          margins_pt=(1.0, 1.0, 1.0, 1.0), pad_pt=LABEL_PAD_PT):
    """The best anchor (display px, M x 2 `anchors`) for label `s`.

    The label's box (its text extent grown by `grow_px`, e.g. a backing
    patch) must lie inside the axes by `margins_pt` (left, bottom, right,
    top: clear of the inward ticks) and, with `within_x`, reach that x
    (from 3 pt right of its left edge to 1 pt past its right edge).
    Candidates are ranked lexicographically by (STRICT: key-mark disks
    touching + reference lines `hlines` (data y) crossing + `avoid` boxes
    overlapped; `key_pts` disks touching; `pts` disks touching; `cost`),
    every touch within `pad_pt`. `cost`: an (M,) array, or a function of
    the (M x 4) boxes. Returns {'xy' (data), 'box' (Bbox px), 'score',
    's', 'kw'}."""
    fig = ax.figure
    pt = fig.dpi / 72.0
    pad = pad_pt * pt
    ext = _extent(ax, s, {k: v for k, v in kw.items() if k != 'bbox'})
    A = np.asarray(anchors, float).reshape(-1, 2)
    boxes = np.column_stack([A[:, 0] + ext[0] - grow_px,
                             A[:, 1] + ext[1] - grow_px,
                             A[:, 0] + ext[2] + grow_px,
                             A[:, 1] + ext[3] + grow_px])
    abb = ax.bbox
    ml, mb, mr, mt = (m * pt for m in margins_pt)
    ok = ((boxes[:, 0] > abb.x0 + ml) & (boxes[:, 2] < abb.x1 - mr)
          & (boxes[:, 1] > abb.y0 + mb) & (boxes[:, 3] < abb.y1 - mt))
    if within_x is not None:
        ok &= (boxes[:, 0] + 3.0 * pt < within_x) & \
              (boxes[:, 2] + 1.0 * pt > within_x)
    c = (np.zeros(len(A)) if cost is None
         else cost(boxes) if callable(cost) else np.asarray(cost, float))
    A, boxes, c = A[ok], boxes[ok], c[ok]
    assert len(A), f'no candidate inside the axes for {s!r}'
    n_strict = _n_touch(strict_pts, boxes, pad)
    lpad = pad + LINE_CLEAR_PT * pt       # lines: + their half width
    for h in hlines:
        y = ax.transData.transform((1.0, h))[1]
        n_strict += ((boxes[:, 1] - lpad < y)
                     & (boxes[:, 3] + lpad > y)).astype(int)
    for b in avoid:
        n_strict += ((boxes[:, 0] < b.x1 + pad) & (boxes[:, 2] > b.x0 - pad)
                     & (boxes[:, 1] < b.y1 + pad)
                     & (boxes[:, 3] > b.y0 - pad)).astype(int)
    n_key = (_n_touch(key_pts, boxes, pad) if key_pts is not None
             else np.zeros(len(A), int))
    n_all = _n_touch(pts, boxes, pad)
    i = np.lexsort((c, n_all, n_key, n_strict))[0]
    return {'xy': ax.transData.inverted().transform(A[i]),
            'box': mtransforms.Bbox(boxes[i].reshape(2, 2)),
            'score': (int(n_strict[i]), int(n_key[i]), int(n_all[i]),
                      float(c[i])),
            's': s, 'kw': kw}


def _best_spot(ax, texts, kw, anchors, cost, lss=LABEL_LS, steps=(1e4, 1e5),
               **spot_kw):
    """The best `_spot` over a label's wraps `texts` (preferred first) x line
    spacings `lss` (most generous first): each rank down costs steps[1]
    (wrap) / steps[0] (spacing) on top of `cost` (array or function of the
    boxes), so a free spot always wins, then -- with the default steps --
    the preferred wrap, then the most generous spacing, then `cost`."""
    best = None
    for k, s in enumerate(texts):
        for i, ls in enumerate(lss):
            extra = steps[1] * k + steps[0] * i
            c = ((lambda box, f=cost, e=extra: f(box) + e) if callable(cost)
                 else np.asarray(cost, float) + extra)
            sp = _spot(ax, s, {**kw, 'linespacing': ls}, anchors, cost=c,
                       **spot_kw)
            if best is None or sp['score'] < best['score']:
                best = sp
    return best


def _put(ax, sp, placed, boxed=False, leader_to=None,
         shrink_pt=0.5, color=None, pts_all=None):
    """Draw the label `_spot` chose; unboxed labels must be free (no key
    mark, line, label or dot touching: asserted). `leader_to` = (x, y)
    display px of a vertical leader's end below the label (x None = pick
    the x along the label's bottom whose segment touches no dot of
    `pts_all`, nearest its left quarter; a given x is checked against
    `pts_all` too). Appends the label's box (and the
    leader's) to `placed`; returns the annotation."""
    s, kw = sp['s'], dict(sp['kw'])
    print(f"label {s.splitlines()[0]!r} (+{s.count(chr(10))} lines, spacing "
          f"{kw.get('linespacing', 1.2):g}): score {sp['score'][:3]}"
          f"{' boxed' if boxed else ''}{' + leader' if leader_to else ''}")
    if not boxed:
        assert sp['score'][:3] == (0, 0, 0), \
            f'{s!r} is not on a free spot: {sp["score"]}'
    box = sp['box']
    fig = ax.figure
    pt = fig.dpi / 72.0
    if leader_to is None:
        t = _ann(ax, s, tuple(sp['xy']), **kw)
        placed.append(box)
        return t
    lx, ly = leader_to
    y_end = ly + shrink_pt * pt * (1 if box.y0 > ly else -1)
    cand = (np.arange(box.x0 + 4 * pt, box.x1 - 4 * pt, 1.0) if lx is None
            else np.array([lx]))
    segs = np.column_stack([cand - pt, np.full(len(cand), min(y_end, box.y0)),
                            cand + pt, np.full(len(cand), max(y_end, box.y0))])
    n = _n_touch(pts_all, segs, LABEL_PAD_PT * pt)
    pref = box.x0 + 0.25 * box.width
    j = np.lexsort((np.abs(cand - pref), n))[0]
    assert n[j] == 0, f'the leader of {s!r} runs through {n[j]} dot(s)'
    lx = cand[j]
    tx, ty = ax.transData.inverted().transform((lx, ly))
    rel = ((lx - box.x0) / box.width, 0.0)
    t = ax.annotate(s, xy=(tx, ty), xytext=tuple(sp['xy']),
                    textcoords='data', annotation_clip=False,
                    arrowprops=dict(arrowstyle='-', color=color or TEXT,
                                    lw=LEADER_LW, relpos=rel, shrinkA=1.0,
                                    shrinkB=shrink_pt), **kw)
    ANNOT.append(t)
    MARKS.append(t.arrow_patch)
    placed += [box, mtransforms.Bbox([[lx - pt, min(ly, box.y0)],
                                      [lx + pt, max(ly, box.y0)]])]
    return t


def key_b(fig):
    """One-line colour key of b's dot clouds under b's title (dots as
    markers, same colours as the clouds)."""
    # round 3: the TRY families carry ONE name set everywhere ('ethanol TRY'
    # / 'isobutanol TRY', as C.FAMILY_LABEL in c and S2); 'scouts' is
    # only the collective noun for the six TRY campaigns
    # starts under the title text (round 3: the longer family names ran
    # past b's right spine from b's left edge)
    # fresh round 2: the key dots as the clouds render (B_DOTS: rim,
    # alpha), so they differ under deutan / protan vision as the clouds do
    x_end, _ = S.inline_key(fig, LETTER_X['b'] + 0.28, HEAD_Y + 0.07, [
        {**_key_dot('unin'), 'text': 'uninformed'},
        {**_key_dot('relay'), 'text': 'TRY-informed'},
        {**_key_dot('etoh'), 'text': FAM_TEXT['etoh']},
        {**_key_dot('ibo'), 'text': FAM_TEXT['ibo']}],
        fontsize=S.FS['key'], item_gap_in=0.13)
    assert x_end <= L['b'][0] + L['b'][2], f'b key runs to {x_end:.2f} in'


# %% Panels c, d, e: the bottom rows ----------------------------------------------------
def _row_axes(fig, key):
    ax = S.inch_axes(fig, *L[key])
    ax.set_ylim(*ROW_YLIM)
    return ax


def panel_c(fig, F):
    U, start = F['U_pct'], F['start_irr_pct']
    camp = F['campaigns']
    ax = _row_axes(fig, 'c')
    # fresh round 1: '0' labelled like a, b and f; 'loss' (full size) moves
    # left inside the band, clear of it
    S.irr_axis(ax, 'x', step=10,
               loss_at=S.loss_label_x(L['c'][2]))
    S.tint_above(ax, U, which='x')
    # fresh round 2: the plateau line ABOVE the role markers (z 5-5.8, under
    # the diamond at 6), so the ethanol-yield ring 0.44 points short of it
    # shows the line crossing it right of its centre, not touching it
    S.plateau_line(ax, U, which='x', zorder=5.9)
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
        # fresh round 1: the returned (filled) marker drawn ABOVE the open
        # one with a white edge, so a close pair (ethanol-TRY titer: 12.7
        # vs 14.4 %) shows both outlines
        _filled(ax, xr, ROW_Y[k], c.marker, dark, s=S_FILLED, z=5.8,
                edge='white', lw=0.8)
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
        for seg in SEP_SEGMENTS:
            fig.add_artist(Line2D(seg, [sep_y, sep_y], transform=tr,
                                  color='0.85', lw=0.6, zorder=0.5))
    marks = [('base', 'D', 'white', P['base'], 6.0, 1.1)]
    for k in C.MAIN_KEYS:
        c = C.campaign(k)
        if c.family == 'profit':
            col = P['relay_text'] if c.is_relay else P['unin_text']
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


# fresh round 2: each d bar ends in its row's returned-design marker (as in
# c), so d and e need no tracking back to c's row labels; d's x limit makes
# room for the longest bar's marker
D_XMAX = 200.0
D_MARK_GAP_PT = 1.5               # bar end -> marker edge


def _row_marker(k):
    """(marker, s, face, edge, lw) of row k's returned design, as in c."""
    if k == 'base':
        return 'D', 26, 'white', P['base'], 1.1
    c = C.campaign(k)
    if c.family == 'profit':
        return (('*', 70, P['relay'], 'white', 0.6) if c.is_relay
                else ('D', 30, P['unin'], 'white', 0.6))
    return c.marker, 26, P[c.family + '_dark'], 'white', 0.5


def panel_d(fig, F):
    ax = _row_axes(fig, 'd')
    pr = F['proteome']
    per_pt = D_XMAX / (L['d'][2] * 72.0)          # g/L per point
    for k in C.ROW_KEYS:
        y = ROW_Y[k]
        et, ib = pr[k]['etoh'], pr[k]['ibo']
        ax.barh(y, et, height=BAR_H, color=P['etoh'], ec='white', lw=0.4,
                zorder=2)
        ax.barh(y, ib, left=et, height=BAR_H, color=P['ibo'], ec='white',
                lw=0.4, zorder=2)
        m, sz, fc, ec, lw = _row_marker(k)
        r = np.sqrt(sz) / 2.0
        xm = et + ib + (r + D_MARK_GAP_PT) * per_pt
        assert xm + r * per_pt < D_XMAX, (k, xm)
        _mark(ax, xm, y, m, sz, fc, ec, lw, z=5)
    ax.set_xlim(0, D_XMAX)
    # fresh round 2: a major tick every 50 g/L (labels at 0 and 100 only:
    # '50' / '150' do not fit at the tick size on 0.72 in)
    ax.xaxis.set_major_locator(FixedLocator([0, 50, 100, 150]))
    ax.xaxis.set_major_formatter(FixedFormatter(['0', '', '100', '']))
    ax.xaxis.set_minor_locator(FixedLocator([25, 75, 125, 175]))
    ax.set_xlabel(S.bold_axis_title(f'Titer\n[{S.UNIT_TITER}]'))
    S.style_ticks(ax, y=False)
    return ax


E_XMAX = 0.30                     # e's x limit (room for the growth tags:
                                  # the first, 'growth x0.84', ends 0.03 in
                                  # short of the spine; asserted)
ADH_DARK = {'adh1': 'etoh_dark', 'adh6': 'ibo_dark'}


def _growth_tag(g):
    """Growth-factor tag: two decimals, three when two would read '1.00'
    (the isobutanol-productivity design, x0.997, is just over the budget)."""
    s = f'{g:.2f}'
    return f'×{g:.3f}' if s == '1.00' else f'×{s}'


E_TAGS = []                       # e's growth tags (their fit is asserted)


def panel_e(fig, F):
    E_TAGS.clear()
    ax = _row_axes(fig, 'e')
    pr = F['proteome']
    bud = F['budget']
    first_tag = True
    for k in C.ROW_KEYS:
        y = ROW_Y[k]
        left = 0.0
        for sec in C.SECTORS:
            v = pr[k][sec.key]
            # fresh round 2: Adh1 / Adh6 hatched in their branch's dark
            # shade (their light fills are the scout families' dot colours)
            kw = ({} if sec.key not in ADH_DARK else
                  dict(hatch=S.ADH_HATCH, hatchcolor=P[ADH_DARK[sec.key]],
                       hatch_linewidth=0.45))
            ax.barh(y, v, left=left, height=BAR_H, color=P[sec.color_key],
                    ec='white', lw=0.4, zorder=2, **kw)
            left += v
        assert abs(left - pr[k]['Phi_M']) < 1e-6, k
        g = pr[k]['growth']
        # tag every design beyond the budget (growth < 1); the first tag
        # names the quantity
        if g < 1.0:
            tag = _growth_tag(g)
            t = _ann(ax, f'growth {tag}' if first_tag else tag,
                     (left + 0.003, y), ha='left', va='center',
                     fontsize=S.FS['note'], fontstyle='italic',
                     color=GREY_TAG)
            E_TAGS.append(t)
            first_tag = False
    ax.axvline(bud, color='k', ls=':', lw=0.9, zorder=3)
    assert pr['base']['Phi_M'] < bud
    # fresh round 2: 4 decimals, as in Fig. S2 and both captions
    _ann(ax, f'budget {C.fmt_budget(bud)}', (bud + 0.0033, ROW_Y['base']),
         ha='left', va='center', fontsize=S.FS['note'], color=TEXT)
    ax.set_xlim(0, E_XMAX)
    ax.xaxis.set_major_locator(FixedLocator([0, 0.1, 0.2]))
    ax.xaxis.set_major_formatter(FixedFormatter(['0', '0.1', '0.2']))
    ax.set_xlabel(S.bold_axis_title(
        f'Metabolic proteome\n[{S.UNIT_PROTEOME}]'))
    S.style_ticks(ax, y=False)
    return ax


# %% Panel f: exploration vs selection --------------------------------------------------
L_F_XLIM = (0.008, 0.93)          # f's logit x limits (fractions)
F_TICKS_PCT = (1, 2, 5, 10, 20, 50, 80)    # 90 dropped: '80 90' crowded
F_ETOH_LABEL_Y = 10.0             # 'ethanol TRY' label centre (IRR, %)
F_LABELS = []                     # f labels whose fit inside f is asserted
# Fresh round 1: the seeds -> star arrow ENDS ON the star; fresh round 2:
# a flatter arc (-0.20 -> -0.08), bowing up-left (F_ARROW_RAD)
F_ARROW_RAD = -0.08
F_ARROW_SHRINK_PT = 7.5           # the star's radius + 1.5 pt
F_MINOR_PCT = (3, 4, 6, 7, 8, 9, 30, 40, 60, 70)   # logit minors (round 2)
# fresh round 2: f alone gets headroom above 30 % (F_YLIM) for the
# TRY-informed speed label above its star (below or left of the star the
# seeds arrow crosses it; right of it the spine is too close), so the
# speed half of T5 is read in f itself
F_YLIM = (S.IRR_LIM[0], 35.5)


def panel_f(fig, F):
    U, start = F['U_pct'], F['start_irr_pct']
    camp = F['campaigns']
    ax = S.inch_axes(fig, *L['f'])
    # fresh round 2: the logit scale named on the axis and shown by its
    # unlabelled minors (only the caption said so)
    S.logit_pct_axis(ax, ticks_pct=F_TICKS_PCT,
                     lim_pct=(100 * L_F_XLIM[0], 100 * L_F_XLIM[1]),
                     minor_pct=F_MINOR_PCT,
                     title='Trials making isobutanol\n[%, logit scale]')
    S.irr_axis(ax, 'y', step=10)
    MARKS.append(_zero_line(ax))
    assert ax.get_ylim() == S.IRR_LIM
    ax.set_ylim(*F_YLIM)
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
    # the seeds: a grey hexagon (the open grey diamond is the starting
    # strain in c)
    xs = F['seeds']['exploration_pct'] / 100
    ys = F['seeds']['best_irr_pct']
    _mark(ax, xs, ys, 'h', 62, 'white', P['seed'], 1.3, z=6)
    xr = camp['relay']['exploration_pct'] / 100
    yr = F['relay']['best']['irr_pct']
    _mark(ax, xr, yr, '*', 150, P['relay'], 'white', 0.8, z=7)

    arr = FancyArrowPatch((xs, ys), (xr, yr), arrowstyle='-|>',
                          color=GREY_ARROW, lw=1.0, ls=(0, (3, 2)),
                          mutation_scale=9, shrinkA=6,
                          shrinkB=F_ARROW_SHRINK_PT,
                          connectionstyle=f'arc3,rad={F_ARROW_RAD}',
                          zorder=6.5)
    ax.add_patch(arr)
    MARKS.append(arr)

    fs, nfs = S.FS['annot'], S.FS['note']
    note = dict(fontsize=nfs, fontstyle='italic', color=NOTE,
                linespacing=1.1)
    # the seeds: label under the hexagon, right-aligned short of the
    # isobutanol-yield stem
    _ann(ax, f"{C.fmt_int(F['seeds']['n'])} seeds", (xs, ys), (13, -7),
         ha='right', va='top', fontsize=nfs, color=GREY_LABEL)
    # group labels (data-anchored; offsets in points). The ethanol label
    # sits right of the ethanol markers
    _ann(ax, FAM_TEXT['etoh'], (pos['ep'][0], F_ETOH_LABEL_Y), (6, 0),
         ha='left', va='center', fontsize=fs, color=P['etoh_dark'])
    # round 4: the uninformed label straddles its own plateau line right
    # of the diamond: the name just above the dashed line, the qualifier
    # (one campaign run: 'run', not 'seed', which here means a preloaded
    # trial) just below it
    _ann(ax, 'uninformed', (xu, U), (6, 2), ha='left', va='bottom',
         fontsize=fs, color=P['unin_text'])
    # fresh round 2: the 'profitability alone' half of T5, parallel to the
    # isobutanol TRY note ('alone: visits, returns low IRR')
    _ann(ax, 'alone: plateau', (xu, U), (6, -2.5), ha='left', va='top',
         **note)
    # round 4: one isobutanol-TRY label in the loss band, beside the
    # isobutanol markers, carrying the note (it sat under the ethanol label)
    _ann(ax, FAM_TEXT['ibo'], (pos['ip'][0], S.LOSS_CENTER), (-7, 1.5),
         ha='right', va='bottom', fontsize=fs, color=P['ibo_dark'])
    _ann(ax, 'alone: visits, returns low IRR', (pos['ip'][0], S.LOSS_CENTER),
         (-7, -0.5), ha='right', va='top', **note)
    # the TRY-informed label centred above its star (F_YLIM headroom), with
    # its speed: the first trial >= 25 % (also marked in a), worded as a
    # threshold so it is not read as the star's own IRR; the uninformed
    # campaign never reached 25 % (its 'alone: plateau' note). Its fit
    # inside f is asserted in build()
    r_star = np.sqrt(150) / 2
    F_LABELS.clear()
    F_LABELS.append(_ann(
        ax, f"TRY-informed:\n≥ 25 % by {C.fmt_trial(F['relay']['first_ge25'])}",
        (xr, yr), (0, r_star + 1.5), ha='center', va='bottom', fontsize=nfs,
        color=P['relay_text'], linespacing=1.1))
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
    """Shared d/e colour key: three rows of a dark / light swatch pair and
    one label (fresh round 2: a strict grid -- 'TCA/acetate' sat mid-row,
    out of line with the Adh column; a two-column grid ran past e). The
    light Adh swatches are hatched as in e. Returns the right end (in)."""
    rows = [('etoh', 'adh1', 'ethanol / Pdc, Adh1'),
            ('ibo', 'adh6', 'isobutanol / ALS→Aro10, Adh6'),
            ('gly', 'tca', 'glycolysis, TCA/acetate')]
    fs, sw, gap = S.FS['key'], 0.11, 0.04
    x_end = x0
    for y, (dk, lk, txt) in zip(BAND_Y_DE, rows):
        light = {'swatch': P[lk], 'size_in': sw, 'text': txt}
        if lk in ADH_DARK:
            light.update(hatch=S.ADH_HATCH, hatchcolor=P[ADH_DARK[lk]])
        xe, _ = S.inline_key(fig, x0, y, [{'swatch': P[dk], 'size_in': sw},
                                          light],
                             fontsize=fs, gap_in=gap, item_gap_in=0.02)
        x_end = max(x_end, xe)
    return x_end


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
        bb = Text.get_window_extent(t, r)         # the text, not its arrow
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
    BLOCK_HIDDEN.clear()
    BLOCK_TEXT.clear()
    COLUMN_COUNT.clear()
    A_TEXTS.clear()
    fig = S.new_figure(*S.MAIN_SIZE)
    ax_a, strip, ins = panel_a(fig, F)
    ax_b = panel_b(fig, F)
    ax_c = panel_c(fig, F)
    row_labels(fig, F, ax_c)
    ax_d = panel_d(fig, F)
    ax_e = panel_e(fig, F)
    ax_f = panel_f(fig, F)
    key_b(fig)
    key_c(fig)
    x_de = key_de(fig, KEY_DE_X)
    assert x_de <= KEY_DE_MAX_X, f'd/e key runs to {x_de:.2f} in'
    x_f0 = max(L['f'][0], x_de + KEY_GAP_IN)
    x_fk = key_f(fig, x_f0)
    print(f'header band: d/e key {KEY_DE_X:.2f}-{x_de:.2f} in, '
          f'f key {x_f0:.2f}-{x_fk:.2f} in')
    assert x_fk <= S.MAIN_SIZE[0] - 0.05, f'f key runs to {x_fk:.2f} in'

    rel = F['relay']
    # titles state findings; a = the seeded search's OWN first gain (the
    # first trial above every seed, not trial 2, which the seeds already
    # beat); e = the profitable co-producer's proteome (asserted below)
    pr = F['proteome']
    assert pr['relay']['pdc'] < 0.5 * pr['base']['pdc']
    assert pr['relay']['Phi_M'] <= F['budget']
    # fresh round 2: a states the co-production (the metabolic-engineering
    # message) before the seed comparison; f no longer generalizes 'visit'
    # to every scout (the ethanol scouts never passed the plateau)
    titles = {
        'a': ("Seeded search co-produces; beats every seed by "
              f"{C.fmt_trial(rel['first_gt_best_seed'])}"),
        'b': f"Every design above {C.fmt_irr(F['U_pct'])} makes isobutanol",
        # round 4: the finding, not the by-construction 'visit > return'
        'c': 'Isobutanol scouts visit high IRR, return low',
        # round 3: e states the panel-wide finding (T1), not the co-producer
        # row's (that is in caption (d, e)); f names who visits and who
        # selects
        'd': 'Titers', 'e': 'Strains differ',
        'f': 'Scouts explore, seeded run selects'}
    ends = {}
    for p, ttl in titles.items():
        y = TOP_LETTER_Y if p in 'ab' else BOT_LETTER_Y
        S.panel_letter(fig, LETTER_X[p], y, p)
        dx = TITLE_DX.get(p, 0.28)
        ends[p] = LETTER_X[p] + dx + _text_w_in(
            fig, ttl, S.FS['title'], fontweight='bold')
        S.panel_title(fig, LETTER_X[p], y, ttl, dx=dx)
    for p, q in (('c', 'd'), ('d', 'e'), ('e', 'f')):
        gap = LETTER_X[q] - ends[p]
        assert gap >= MIN_TITLE_GAP_IN, f'title {p} ends {gap:.2f} in from {q}'
    assert ends['f'] <= S.MAIN_SIZE[0] - 0.05, f"f title ends {ends['f']:.2f}"
    assert ends['a'] <= LETTER_X['b'] - MIN_TITLE_GAP_IN, ends['a']
    # e's growth tags end inside e (E_XMAX)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for t in E_TAGS:
        assert t.get_window_extent(r).x1 < ax_e.bbox.x1 - 2.0, t.get_text()
    # f's star label stays inside f (clear of the right spine's ticks and
    # below the top spine)
    for t in F_LABELS:
        bb = t.get_window_extent(r)
        assert bb.x1 < ax_f.bbox.x1 - 6.0, t.get_text()
        assert bb.y1 < ax_f.bbox.y1 - 6.0, t.get_text()
    return fig, dict(a=ax_a, strip=strip, inset=ins, b=ax_b, c=ax_c, d=ax_d,
                     e=ax_e, f=ax_f)


def inset_clearance(fig, axes, gap_pt=INSET_GAP_PT):
    """a's inset, its tick labels included, must sit wholly in a's empty
    region: >= `gap_pt` above the IRR = 0 line, >= `gap_pt` below every a
    label it lies under, and inside a clear of the right-spine ticks.
    Returns a list of messages; [] = pass."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ax, ins = axes['a'], axes['inset']
    tb = ins.get_tightbbox(r)
    gap = gap_pt * fig.dpi / 72.0
    msgs = []
    y0 = ax.transData.transform((1.0, 0.0))[1]
    if tb.y0 < y0 + gap:
        msgs.append(f'inset reaches {(tb.y0 - y0) * 72 / fig.dpi:.1f} pt '
                    f'from the IRR = 0 line (< {gap_pt} pt)')
    for t in ANNOT:
        if t.axes is not ax:
            continue
        lb = t.get_window_extent(r)
        if lb.x1 > tb.x0 and lb.x0 < tb.x1 and lb.y0 >= tb.y0 \
                and lb.y0 < tb.y1 + gap:
            msgs.append(f'inset top {(lb.y0 - tb.y1) * 72 / fig.dpi:.1f} pt '
                        f'below {t.get_text()!r} (< {gap_pt} pt)')
    tick_pt = 4.0                         # a's inward right-spine ticks
    if (tb.x1 > ax.bbox.x1 - (tick_pt * fig.dpi / 72.0 + gap)
            or tb.x0 < ax.bbox.x0 or tb.y1 > ax.bbox.y1):
        msgs.append(f'inset {tb} not inside a clear of its spine ticks')
    a = A_TEXTS['start'].get_window_extent(r)
    print(f'inset: {(tb.y0 - y0) * 72 / fig.dpi:.1f} pt above IRR = 0, '
          f'{(a.y0 - tb.y1) * 72 / fig.dpi:.1f} pt under the starting-'
          f'strain label, {(ax.bbox.x1 - tb.x1) * 72 / fig.dpi:.1f} pt from '
          f"a's right spine; frame {ins.bbox.width / fig.dpi:.2f} x "
          f'{ins.bbox.height / fig.dpi:.2f} in')
    return msgs


def run_checks(fig, axes, raise_on_fail=True):
    res = S.check_figure(fig, size=S.MAIN_SIZE, raise_on_fail=False,
                         marker_exempt=tuple(BLOCK_TEXT))
    fig.canvas.draw()
    inset_bb = axes['inset'].get_tightbbox(fig.canvas.get_renderer())
    clashes = mark_clashes(fig, ANNOT, MARKS, regions=[('inset', inset_bb)])
    clashes += inset_clearance(fig, axes)
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


METHODS_MARK = '*Methods notes (not part of the legend).*'
LEGEND_MAX_WORDS = 350            # fresh round 2 (journal legend limits)


def caption(F):
    """The legend (title + panels + synthesis, every number from the facts)
    followed by a short block of methods notes that belong in Methods, not
    in the legend (round 2: the legend was 788 words)."""
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
    pi25, pi_best = rel['sims'][25]['PI'], rel['bsf_PI'][rb['sim']]
    assert pi25 == rel['bsf_PI'][25]
    # the reproducibility margin, from the eight campaigns' own duplicate
    # designs (facts['repro']; round-1 fixer)
    rp = F['repro']
    margin = rp['max_dPI']
    ratio = (pi_best - pi25) / margin
    ratio_word = WORDS.get(int(round(ratio)), f'{ratio:.0f}')
    ibo_ret = [camp[k]['returned']['irr_pct'] for k in C.IBO_SCOUTS]
    ibo_fin = [v for v in ibo_ret if np.isfinite(v) and v >= 0]
    n_noirr = sum(1 for v in ibo_ret if not np.isfinite(v))
    etoh_ret = [camp[k]['returned']['irr_pct'] for k in C.ETOH_SCOUTS]
    etoh_bv = max(camp[k]['best_visit']['irr_pct'] for k in C.ETOH_SCOUTS)
    tagged = [k for k in C.ROW_KEYS if pr[k]['growth'] < 1.0]
    growth = ', '.join(_growth_tag(pr[k]['growth']) for k in tagged)
    prg = F['progress']
    s1, s2 = rel['sims'][1], rel['sims'][2]
    s_seed = rel['sims'][rel['first_gt_best_seed']]
    first_word = ('lost money' if C.is_loss(s1['irr_pct']) else
                  f"reached {C.fmt_irr(s1['irr_pct'])}")
    n_hits = rel['n_best_bound_hits']
    ibo_counts = ', '.join(C.fmt_int(camp[k]['n_gt_U'])
                           for k in C.IBO_SCOUTS[:-1])
    ibo_counts += f" and {C.fmt_int(camp[C.IBO_SCOUTS[-1]]['n_gt_U'])}"
    ret_txt = (', '.join(C.fmt_irr(v) for v in ibo_fin)
               + (f' and {WORDS[n_noirr]} with no IRR' if n_noirr else ''))
    pdc_ratio = pr['relay']['pdc'] / pr['base']['pdc']
    ib5 = ui['ibo_ge5']
    fam = FAM_TEXT
    # round 4: the seeds already held co-producers above the plateau (the
    # headline and (a) say so; skeptic-referee, round 3)
    n_cop_word = WORDS.get(seeds['n_coprod_gt_U'], str(seeds['n_coprod_gt_U']))
    best_cop = C.fmt_irr(seeds['best_coprod_irr_pct'])
    # fresh round 1: the headline in two short sentences, one per finding
    # (the one-sentence, 50-word headline nested 'four of them ... above its
    # plateau'); the seeded co-producers and their best move to (a)
    sobol = F['sobol']
    # one basis for every simulation count: simulated trials, FAIL rows
    # included (max_sim; round-1 fixer)
    n_seeded_sims = g['n_try_sims'] + camp['relay']['max_sim']
    fails = {k: camp[k]['n_fail'] for k in C.MAIN_KEYS}
    role_order = ' / '.join(C.campaign(k).role for k in C.ETOH_SCOUTS)
    v25 = rel['sims'][rel['first_ge25']]['bsf_pct']
    col = COLUMN_COUNT
    pi_gain = pi_best - pi25
    # fresh round 2: the legend cut to <= LEGEND_MAX_WORDS (it was 653):
    # definitions, the seed selection and the compute stay; the per-panel
    # numeric recaps the panels already show move to the methods notes.
    # The synthesis speaks of these single runs only
    n_cop = _first_coprod_run(C.complete('relay').sort_values('sim'))
    s3 = rel['sims'][3]
    bh = BLOCK_HIDDEN
    legend = [
        f"**Figure N | Isobutanol TRY scouts visit designs above the {U} "
        f"ethanol-only plateau of the uninformed profitability campaign but "
        f"return low-IRR designs; seeded with their trials, "
        f"profitability search co-produces above that plateau on "
        f"{C.fmt_trial(rel['first_gt_U'])}, beats every seed by "
        f"{C.fmt_trial(rel['first_gt_best_seed'])} and reaches "
        f"{C.fmt_irr(rb['irr_pct'])}.** "
        f"Eight Gaussian-process campaigns searched one "
        f"{n_dec}-dimensional strain-design and feeding space. Two maximized "
        f"profitability (profitability index, shown as internal rate of "
        f"return, IRR): **uninformed** (cyan, "
        f"{C.fmt_int(camp['unin']['max_sim'])} trials) and **TRY-informed** "
        f"(teal, {C.fmt_int(camp['relay']['max_sim'])} trials after "
        f"{C.fmt_int(seeds['n'])} preloaded seeds). Six TRY campaigns (the "
        f"scouts) each maximized one titer, rate or yield of ethanol (amber) "
        f"or isobutanol (violet); their {C.fmt_int(g['n_try_sims'])} "
        f"simulations, by-products of those campaigns, supplied the "
        f"seeds: "
        f"{seeds['n_keep_above']} scout trials at least as profitable as "
        f"the starting strain and {seeds['n_maximin']} space-filling "
        f"(maximin) ones. Trial: simulated index. Dashed cyan: plateau "
        f"(tinted above); dotted grey: starting strain, {start}; dotted "
        f"violet: best seed, {best_seed}; loss: IRR < 0 or none. Making "
        f"isobutanol: ≥ {thr} {unit}; co-production: ≥ {thr} {unit} of each "
        f"alcohol.",
        '',
        f"**(a)** Best IRR so far; strip: the seeds (teal rings: "
        f"co-producers above the plateau). Inset: refinement, linear trial "
        f"axis.",
        '',
        f"**(b)** IRR vs isobutanol share of the alcohol titer, every design "
        f"making ≥ {C.SHARE_MIN_ALCOHOL:g} {unit} alcohol (the 100 % column "
        f"spread over {100 - COL_JITTER:g}–100 %). Large cyan dots: "
        f"uninformed designs making isobutanol. Open violet ○ □ △: "
        f"{fam['ibo']} best visits; ◆ ★: profitability bests.",
        '',
        f"**(c)** Each campaign's trials (strip), highest-IRR visit (open) "
        f"and returned design (filled: the argmax of its own objective).",
        '',
        f"**(d, e)** Titers and metabolic proteome of the returned designs "
        f"(markers as in c). Dotted: penalty-free budget, "
        f"{C.fmt_budget(F['budget'])} g·(g DCW)⁻¹; ×: growth factor beyond "
        f"it. Hatched: Adh1, Adh6.",
        '',
        f"**(f)** Share of trials making isobutanol vs best-visit (open) and "
        f"returned (filled) IRR. Hexagon: the seeds.",
        '',
        f"In these single runs, neither objective alone returned a "
        f"profitable co-producer: the unseeded profitability campaign "
        f"plateaued at {U} on ethanol, and the scouts returned ethanol-only "
        f"or isobutanol-only designs. Seeded with the scouts' trials, the "
        f"profitability campaign beat every seed by "
        f"{C.fmt_trial(rel['first_gt_best_seed'])} and reached "
        f"{C.fmt_irr(rb['irr_pct'])}. The seeded route used "
        f"{C.fmt_int(n_seeded_sims)} simulations "
        f"({C.fmt_int(g['n_try_sims'])} scout + "
        f"{C.fmt_int(camp['relay']['max_sim'])}), the uninformed campaign "
        f"{C.fmt_int(camp['unin']['max_sim'])}.",
    ]
    methods = [
        METHODS_MARK,
        f"* Random seed {SEED} for all eight campaigns. Every campaign starts "
        f"from scenario A and is subject to the enzyme-burden constraint; "
        f"titers are per litre of water. The TRY-informed campaign's preload "
        f"({C.fmt_int(seeds['n'])} trials: {C.fmt_int(seeds['n_ibo'])} from "
        f"the {fam['ibo']} and {C.fmt_int(seeds['n_etoh'])} from the "
        f"{fam['etoh']} campaigns) exceeded its start-up count, so it had no "
        f"space-filling start-up: its GP proposed every trial. The seeds "
        f"were preloaded with their recorded objective values, without "
        f"re-simulation. The {seeds['n_keep_above']} are the "
        f"{seeds['n_eligible_above']} scout trials at least as profitable as "
        f"the starting strain minus {seeds['n_quarantined_above']} "
        f"convergence-quarantined and {WORDS[seeds['n_duplicate_above']]} "
        f"duplicate; picked by PI, so the seed set is not blind to "
        f"profitability. The {g['n_try_gt_U']} scout trials above {U} are "
        f"among them ({g['n_try_gt_U_ibo_only']} isobutanol-only, < {thr} "
        f"{unit} ethanol; {g['n_try_gt_U_coprod']} co-producers, best "
        f"{best_cop}).",
        f"* Trial = CSV trial_number + 1 (TRY-informed: trial_number − "
        f"{rel['sim_offset']}; its best, {C.fmt_trial(rb['sim'])}, is "
        f"trial_number {C.fmt_int(rb['sim'] + rel['sim_offset'])}), 1-based, "
        f"preloaded seeds excluded; only completed trials are shown. FAIL "
        f"rows, excluded everywhere: uninformed {fails['unin']}, "
        f"TRY-informed {fails['relay']}, {fam['etoh']} "
        f"{'/'.join(str(fails[k]) for k in C.ETOH_SCOUTS)} and {fam['ibo']} "
        f"{'/'.join(str(fails[k]) for k in C.IBO_SCOUTS)} ({role_order}); "
        f"simulation counts include them.",
        f"* PI = net present value at a {HURDLE_PCT} % hurdle rate / total "
        f"capital investment. PI is deterministic but depends slightly on "
        f"the converged state each simulation starts from (the previous "
        f"trial's). The reproducibility margin, ΔPI {margin:.4f}, is the "
        f"largest PI difference between two campaigns' simulations of the "
        f"same design after different trial histories ({rp['n_pairs']} "
        f"pairs; {rp['max_dIRR_pts']:.2f} IRR points at that pair); smaller "
        f"differences are not interpreted.",
        f"* (a) The uninformed campaign plateaued at {U} (ethanol only) from "
        f"{C.fmt_trial(ui['plateau_sim'])} to "
        f"{C.fmt_int(camp['unin']['max_sim'])}; the cyan strip is its "
        f"space-filling start-up. The TRY-informed campaign {first_word} on "
        f"{C.fmt_trial(1)}, co-produced at {C.fmt_irr(s2['irr_pct'])} on "
        f"{C.fmt_trial(2)}, held an isobutanol-only incumbent "
        f"({C.fmt_irr(s3['irr_pct'])}) from {C.fmt_trial(3)} to "
        f"{C.fmt_trial(n_cop - 1)}, and every improvement from "
        f"{C.fmt_trial(n_cop)} on co-produced; it beat every seed on "
        f"{C.fmt_trial(rel['first_gt_best_seed'])} "
        f"({C.fmt_irr(s_seed['irr_pct'])}), reached {C.fmt_irr(v25)} on "
        f"{C.fmt_trial(rel['first_ge25'])}, and refined "
        f"{C.fmt_irr(rel['sims'][25]['bsf_pct'])} ({C.fmt_trial(25)}) to "
        f"{C.fmt_irr(rb['irr_pct'])} ({C.fmt_trial(rb['sim'])}), a PI gain "
        f"of {pi_gain:.3f}, {ratio_word} times the reproducibility margin. "
        f"All {seeds['n_gt_U']} seeds above {U} are from {fam['ibo']} "
        f"campaigns, {n_cop_word} of them co-producers.",
        f"* (b) Shown: {C.fmt_int(g['n_alcohol_ge1'])} designs; the "
        f"{C.fmt_int(g['n_alcohol_lt1'])} making < "
        f"{C.SHARE_MIN_ALCOHOL:g} {unit} alcohol (none has an IRR) are "
        f"omitted; {fam['etoh']} designs lie mostly under the uninformed "
        f"column. The uninformed campaign's {ib5['n']} isobutanol-making "
        f"designs were all worse than its plateau (best "
        f"{C.fmt_irr(ib5['max_irr_pct'])}) and more often lost money "
        f"({_pct0(ib5['pct_losing'])} vs "
        f"{_pct0(ui['ibo_lt5']['pct_losing'])}); it proposed "
        f"{ib5['n_after_plateau']} of them after "
        f"{C.fmt_trial(ui['plateau_sim'])} "
        f"({ib5['pct_of_post_plateau']:.0f} % of its later trials). Their "
        f"label's opaque backing hides {bh['outlined']} of them and "
        f"{bh['other']} other designs (never the best one). The 100 % column "
        f"holds {col['n']} isobutanol-only scout designs above {U} (drawn "
        f"over the TRY-informed dots) and {col['n_relay']} TRY-informed "
        f"ones; the {WORDS[col['n_ibo_only'] - col['n']]} other "
        f"isobutanol-only scout design above {U} is the titer best visit "
        f"(□). No design making < {thr} {unit} isobutanol exceeds {U} (b's "
        f"title); all {g['n_relay_gt_best_seed']} designs above {best_seed} "
        f"are TRY-informed co-producers; ★ makes {rb['etoh']:.0f} {unit} "
        f"ethanol + {rb['ibo']:.0f} {unit} isobutanol.",
        f"* (c) The {fam['ibo']} campaigns visited {ibo_counts} designs above "
        f"{U} but returned {ret_txt}; the {fam['etoh']} campaigns never "
        f"exceeded {C.fmt_irr(etoh_bv)} and returned {min(etoh_ret):.1f}–"
        f"{C.fmt_irr(max(etoh_ret))}.",
        f"* (d, e) The returned strains differ in products and proteome. The "
        f"TRY-informed co-producer keeps a reduced ethanol branch (Pdc "
        f"{pdc_ratio:.2f}× the starting strain's) and stays within the "
        f"budget. Pdc and ALS→Aro10 commit carbon to ethanol and isobutanol. "
        f"The translation (φ_T = {F['phi_T']:.3f}) and housekeeping sectors "
        f"are omitted; the {WORDS[len(tagged)]} designs beyond the budget "
        f"grow at {growth} of their penalty-free rate.",
        f"* (f) The uninformed campaign made isobutanol in "
        f"{camp['unin']['exploration_pct']:.1f} % of its trials and never "
        f"reached 25 %; the TRY-informed campaign made isobutanol in "
        f"{camp['relay']['exploration_pct']:.1f} % and reached 25 % on "
        f"{C.fmt_trial(prg['relay']['first_ge25'])}. Hexagon: x, share of "
        f"seeds making isobutanol ({seeds['exploration_pct']:.1f} %); y, "
        f"best seed.",
        f"* No compute-matched control (e.g. "
        f"{C.fmt_int(n_seeded_sims)} random or space-filling seeds) was run. "
        f"Two {C.fmt_int(sobol['campaign']['n'])}-design Sobol' random "
        f"samples of the same space bracket it: under the campaigns' own "
        f"sampling measure no design exceeded the {HURDLE_PCT} % hurdle "
        f"(best {C.fmt_irr(sobol['campaign']['max_irr_pct'])}); under a "
        f"linear engineering prior {sobol['screening']['n_gt_hurdle']} did "
        f"({sobol['screening']['n_gt_U']} above the {U} plateau; best "
        f"{C.fmt_irr(sobol['screening']['max_irr_pct'])}, below the best "
        f"seed). These single runs therefore show seeding buying speed "
        f"here, not unique access to designs above the plateau.",
        f"* The best design lies at {WORDS.get(n_hits, n_hits)} search "
        f"bounds, so higher IRR may exist outside the searched ranges.",
        f"* IRRs are at default prices under the model version used for the "
        f"campaigns (starting strain {start}; "
        f"{C.fmt_irr(100 * START_IRR_CURRENT_MODEL)} under the current model "
        f"version).",
    ]
    n_leg = len('\n'.join(legend).split())
    assert n_leg <= LEGEND_MAX_WORDS, f'legend {n_leg} words'
    assert 0 < sobol['screening']['n_gt_U'] <= sobol['screening'][
        'n_gt_hurdle']
    assert col['n_relay'] > 0
    # fresh round 1 claims: the seed-selection arithmetic, the compute
    # caveat's bracket, the column count, the PI gain vs the margin
    assert (seeds['n_eligible_above'] - seeds['n_quarantined_above']
            - seeds['n_duplicate_above'] == seeds['n_keep_above'])
    assert seeds['n_quarantined_kept'] == 0
    assert seeds['duplicate_max_rel_diff'] < 1e-9
    assert sobol['campaign']['n_gt_hurdle'] == 0
    assert (0 < sobol['screening']['n_gt_hurdle']
            and sobol['screening']['max_irr_pct'] < seeds['best_irr_pct'])
    assert col['n_ibo_only'] - col['n'] == 1
    assert camp['it']['best_visit']['share_pct'] < COL_SHARE_MIN
    assert camp['it']['best_visit']['etoh'] < C.ETOH_THRESHOLD
    assert abs(pi_gain / margin - ratio) < 1e-12
    # the margin's basis: identical histories reproduce PI bitwise, so
    # the margin is load-path drift, not objective noise
    assert rp['startup_max_dPI'] < 1e-12 and rp['n_pairs'] > 0
    assert all(rel['sims'][k]['class'] == 'coprod'
               for k in (rel['first_gt_best_seed'], rel['first_ge25']))
    # the caption's quoted claims must agree with the figure's facts
    assert ratio > 1, ratio
    # the synthesis: the unseeded campaign's best = its plateau, ethanol
    # only; no scout returned a co-producer; (f): the uninformed campaign
    # never reached 25 %, the TRY-informed one did
    assert abs(prg['unin']['best_pct'] - F['U_pct']) < 1e-9
    assert camp['unin']['returned']['class'] == 'etoh'
    assert all(camp[k]['returned']['class'] in ('etoh', 'ibo')
               for k in C.SCOUT_KEYS)
    assert prg['unin']['first_ge25'] is None
    assert prg['relay']['first_ge25'] == rel['first_ge25'] is not None
    # round 4: the headline's seeded claims
    assert seeds['n_coprod_gt_U'] == g['n_try_gt_U_coprod']
    assert rel['first_gt_best_seed'] < 100
    return '\n'.join(legend) + '\n\n' + '\n'.join(methods) + '\n'


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
        n_leg = len(txt.split(METHODS_MARK)[0].split())
        print(f'caption -> {CAPTION_PATH} (legend {n_leg} words, '
              f'{len(txt.split())} with the methods notes)')
    C.assert_sim_safe()
    print('sim-safety: OK')
    return out


if __name__ == '__main__':
    main()
