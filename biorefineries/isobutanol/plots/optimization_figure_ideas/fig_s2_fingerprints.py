#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Fig. S2 of the kinetic-BO campaign figures: design fingerprints.

One heatmap of the 12 decision variables of the `metabolic_split_12d`
campaigns for the starting strain (scenario A), every campaign's RETURNED
design (argmax of its own objective) and the three isobutanol TRY scouts'
BEST VISITS (highest-IRR trial, italic rows). It replaces the 15-panel
`parameters` variant of plot_kin_opt_parameter_sets.py and carries take-away
T1 (the returned strains differ in their enzyme levers and feeding) at the
parameter level, plus the T3 detail that the isobutanol scouts' best visits
are different designs from the ones they returned.

Encoding (spec section 5.2, figwork/figure_spec.md)
    Native capacities  k_3 (Pdc), k_6 (Adh1), glycolysis, k_17 (Adh6):
                       log2 fold vs the starting strain, PuOr_r, clipped +-4
    Isobutanol pathway k_13 (ALS), ehrlich_downstream (Ilv5/Ilv3/Aro10):
                       absolute g/L/h, Purples 0-4
    Inhibition         the three product-inhibition multipliers: BrBG on
                       log2(multiplier), 1 = the starting strain
    Feeding            threshold, target delta, spikes (actual / cap):
                       greys within the search range
    Outlined cells sit at a search bound (_common.at_bound, 1e-6 relative).
    Right-hand table: IRR, own objective, Phi_M and the growth factor.

Every number in the figure and in the caption file it writes
(kinBO_S2_fingerprints_caption.md, next to this script) is computed from the
campaign CSVs through _common; the caption's claims are asserted against
_common.EXPECTED and the figure-local S2_EXPECTED before anything is drawn.

SIM-SAFE: pandas / numpy / matplotlib only; pk is loaded by file path inside
_common.baseline_record(). Run with
    "C:/Users/saran/anaconda3/envs/IBO_2026/python.exe" fig_s2_fingerprints.py
Output: analyses/results/publication/Optimization-figures/
    kinBO_S2_fingerprints_<YYYY.MM.DD-HH.MM>.{png,pdf} (+ _latest copies).
"""
import os
import sys
from dataclasses import dataclass

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C                                           # noqa: E402
import _style as S                                            # noqa: E402
from matplotlib import colormaps                              # noqa: E402
from matplotlib.colors import ListedColormap, Normalize, to_rgb  # noqa: E402
from matplotlib.lines import Line2D                           # noqa: E402
from matplotlib.patches import Rectangle                      # noqa: E402

STEM = 'kinBO_S2_fingerprints'
CAPTION_PATH = os.path.join(C.FIGDIR, STEM + '_caption.md')

# %% Layout (inches on the 10 x 5.6 in canvas; origin bottom-left) -----------------
W, H = S.S2_SIZE
X_GROUP = 0.10            # group header (left-aligned)
X_MARK = 0.26             # role marker centre
X_LABEL = 0.37            # row label (left-aligned)
X_TRIAL = 1.53            # trial number (right-aligned)
X0 = 1.61                 # heatmap left edge
CW = 0.46                 # cell width
BLOCK_GAP = 0.10          # gap between column blocks
CELL_PAD = 0.016          # white separator between cells
RH = 0.27                 # design-row pitch
HEADER_GAP = 0.25         # group-header row pitch
Y_BOT = 0.70              # heatmap bottom edge
CBAR_Y, CBAR_H = 0.33, 0.075   # colour bars under the heatmap
TABLE_GAP = 0.10          # heatmap -> right table (IRR, own objective,
RIGHT_MARGIN = 0.05       #   Phi_M, growth; placed by place_table)
KEY_PITCH = 0.19          # line pitch of the keys / notes above the table
LABEL_ROT = 45.0          # column-label rotation [deg]
OUTLINE_LW = 1.4          # bound-hit outline [pt]

FS_ROW = 10               # row labels / right-table values
FS_CELL = S.FS['note']    # 9-pt cell text (the floor)
FS_HEAD = S.FS['note']    # 9-pt table headers / column labels / notes
FS_BLOCK = S.FS['annot']  # 10-pt bold block titles


# %% Figure-local expected values (the caption's claims; spec section 5.2) ----------
E = C.E
S2_EXPECTED = {
    # TRY-informed returned design: reduced Pdc, baseline-like glycolysis,
    # ALS and Adh6 at their upper bounds, within the proteome budget
    'relay_k3_fold': E('0.33'), 'relay_k6_fold': E('2.6'),
    'relay_gly_fold': E('0.98'),
    'relay_hits': {'k_13': E('hi'), 'k_17': E('hi')},
    'relay_Phi_M': E('0.130'), 'relay_within_budget': E(True),
    # uninformed: Pdc at its 4x bound, isobutanol pathway at its floors
    'unin_k3_fold': E('4.0'),
    'unin_hits': {'k_3': E('hi'), 'k_13': E('lo'),
                  'ehrlich_downstream': E('lo')},
    'unin_Phi_M': E('0.153'), 'unin_growth': E('0.84'),
    # isobutanol-yield pair: consecutive trials, returned vs best visit
    'iy_sims': {'returned': E(842), 'best_visit': E(843)},
    'iy_yield_gain_pct': E('7.6'), 'iy_irr_drop_pts': E('21.2'),
    # isobutanol titer / productivity: returned designs knock Adh1 down,
    # best visits keep it above the starting strain
    'it_ip_ret_k6_fold_max': E('0.002'),
    'it_ip_ret_etoh_max': E('0.2'),
    'it_ip_bv_k6_fold': {'it': E('1.4'), 'ip': E('1.4')},
    'it_ip_bv_etoh': {'it': E('4.0'), 'ip': E('50.1')},
    'it_ip_bv_ibo': {'it': E('34.8'), 'ip': E('29.9')},
    'n_rows': E(12),
}


def _check(exp, got, path='S2'):
    """Assert `got` (nested dict) against `exp` (nested dict of Expect)."""
    bad = []
    if isinstance(exp, dict):
        for k, e in exp.items():
            bad += _check(e, got[k], f'{path}.{k}')
        return bad
    v, ok = got, False
    if isinstance(exp.value, (str, bool)) or exp.value is None:
        ok = v == exp.value
    else:
        ok = abs(float(v) - float(exp.value)) <= exp.tol
    if not ok:
        bad.append(f'{path}: computed {v!r}, expected {exp.value!r} '
                   f'(tol {exp.tol:.2g})')
    return bad


# %% Columns and blocks ----------------------------------------------------------------
@dataclass(frozen=True)
class Column:
    var: str               # decision variable (_common.DECISION_VARS)
    label: str             # rotated column label


@dataclass(frozen=True)
class Block:
    key: str
    title: str
    descriptor: str        # 9-pt line above the colour bar
    columns: tuple


BLOCKS = (
    Block('native', 'Native capacities', 'fold vs starting strain', (
        Column('k_3', 'Pdc (k$_3$)'), Column('k_6', 'Adh1 (k$_6$)'),
        Column('glycolysis', 'glycolysis'),
        Column('k_17', 'Adh6 (k$_{17}$)'))),
    Block('path', 'Isobutanol pathway', 'absolute, ' + S.UNIT_PROD, (
        Column('k_13', 'ALS (k$_{13}$)'),
        Column('ehrlich_downstream', 'Ilv5–Aro10'))),
    Block('inhib', 'Inhibition', 'multiplier', (
        Column('inhib_ethanol', 'by ethanol'),
        Column('inhib_isobutanol', 'by isobutanol'),
        Column('inhib_acetate', 'by acetate'))),
    Block('feed', 'Feeding', 'within search range', (
        Column('threshold_conc', 'threshold'),
        Column('target_delta', 'target Δ'),
        Column('max_n_spikes', 'spikes / cap'))),
)
FOLD_CLIP = 4.0                                   # log2 fold, +-4 (1/16x-16x)
INHIB_LIM = max(abs(np.log2(C.BANDS['inhib_ethanol'][0])),
                abs(np.log2(C.BANDS['inhib_ethanol'][1])))
PATH_LIM = C.BANDS['k_13'][1]                     # 0-4 g/L/h
GREYS = ListedColormap(colormaps['Greys'](np.linspace(0.03, 0.80, 256)))
CMAPS = {'native': colormaps['PuOr_r'], 'path': colormaps['Purples'],
         'inhib': colormaps['BrBG'], 'feed': GREYS}
NORMS = {'native': Normalize(-FOLD_CLIP, FOLD_CLIP),
         'path': Normalize(0.0, PATH_LIM),
         'inhib': Normalize(-INHIB_LIM, INHIB_LIM),
         'feed': Normalize(0.0, 1.0)}


def scale_value(block, var, value, base):
    """Decision value -> the block's colour coordinate."""
    if block == 'native':
        return float(np.clip(np.log2(value / base[var]), -FOLD_CLIP,
                             FOLD_CLIP))
    if block == 'path':
        return float(np.clip(value, 0.0, PATH_LIM))
    if block == 'inhib':
        return float(np.log2(value))
    lo, hi = C.BANDS[var]
    return float(np.clip((value - lo) / (hi - lo), 0.0, 1.0))


def fmt_fold(f):
    if f < 0.01:
        return '<0.01×'
    if f < 1.0:
        return f'{f:.2f}×'
    if f < 10.0:
        return f'{f:.1f}×'
    return f'{f:.0f}×'


def fmt_abs(v):
    if v == 0:
        return '0'
    if v < 0.01:
        return f'{v:.3f}'
    if v < 1.0:
        return f'{v:.2f}'
    return f'{v:.1f}'


def cell_text(block, var, value, base):
    if block == 'native':
        return fmt_fold(value / base[var])
    if block == 'path':
        return fmt_abs(value)
    if block == 'inhib':
        return f'{value:.2f}'
    return f'{value:.0f}'


def fmt_growth(g):
    """'×0.84'; a slight derate that would round to ×1.00 keeps 3 decimals
    ('×0.997'), so a Phi_M just over the budget never reads as not
    derated."""
    return f'×{g:.3f}' if 0.99 <= g < 0.9995 else f'×{g:.2f}'


def luma(rgb):
    """Rec. 601 luma of gamma-encoded sRGB (white text below 0.45)."""
    r, g, b = rgb[:3]
    return 0.299 * r + 0.587 * g + 0.114 * b


# %% Rows ----------------------------------------------------------------------------------
OBJ_FMT = {   # campaign role -> (own-objective formatter on the row record)
    'yield': lambda v: f'{v:.3f} g·g$^{{-1}}$',
    'titer': lambda v: f'{v:.1f} {S.UNIT_TITER}',
    'productivity': lambda v: f'{v:.2f} {S.UNIT_PROD}',
}


def design_rows(facts):
    """The 12 heatmap rows in drawing order (top to bottom), with group
    headers: list of ('header', family) / ('row', dict) items."""
    base = C.baseline_record()
    bA = facts['baseline_A']
    items = [('row', {
        'kind': 'base', 'label': 'Starting strain', 'style': 'normal',
        'color': S.TEXT, 'marker': dict(marker='D', mfc='white',
                                        mec=S.PALETTE['base'], mew=1.3),
        'sim': None, 'decision': dict(base['decision']),
        'n_spikes': float(bA['n_glu_spikes']),
        'irr_pct': facts['start_irr_pct'], 'obj': '–',
        'Phi_M': float(base['Phi_M']), 'growth': float(base['burden_factor']),
        'hits': {}})]
    groups = (('profit', C.PROFIT_KEYS), ('etoh', C.ETOH_SCOUTS),
              ('ibo', C.IBO_SCOUTS))
    for fam, keys in groups:
        items.append(('header', fam))
        for k in keys:
            c = C.campaign(k)
            st = S.campaign_style(c)
            kinds = ('ret', 'bv') if fam == 'ibo' else ('ret',)
            for kind in kinds:
                row = C.returned_row(k) if kind == 'ret' else \
                    C.best_visit_row(k)
                rec = C.row_record(row)
                camp = facts['campaigns'][k]['returned' if kind == 'ret'
                                             else 'best_visit']
                assert rec['trial'] == camp['trial'], (k, kind)
                if fam == 'profit':
                    obj = f'PI {rec["PI"]:.3f}'
                else:
                    obj = OBJ_FMT[c.role](rec['objective'])
                ms = 9.5 if c.role == 'profit_relay' else 6.0
                if kind == 'ret':
                    mk = dict(marker=c.marker, mfc=st['dark'], mec='white'
                              if fam == 'profit' else st['dark'], mew=0.6,
                              ms=ms)
                else:
                    mk = dict(marker=c.marker, mfc='white', mec=st['dark'],
                              mew=1.3, ms=ms)
                items.append(('row', {
                    'kind': kind, 'key': k,
                    'label': c.label if kind == 'ret' else 'best visit',
                    'style': 'italic' if kind == 'bv' else 'normal',
                    'color': st['text'], 'marker': mk, 'sim': rec['sim'],
                    'decision': dict(rec['decision']),
                    'n_spikes': rec['n_spikes'], 'irr_pct': rec['irr_pct'],
                    'obj': obj, 'Phi_M': rec['Phi_M'],
                    'growth': rec['growth'], 'hits': C.bound_hits(row),
                    'rec': rec}))
    return items


# %% Facts behind the caption -----------------------------------------------------------------
def s2_facts(facts, items):
    base = facts['baseline_decision']
    rows = {(r['key'], r['kind']): r for t, r in items
            if t == 'row' and r['kind'] != 'base'}
    rel, uni = rows[('relay', 'ret')], rows[('unin', 'ret')]
    iy_r, iy_b = rows[('iy', 'ret')]['rec'], rows[('iy', 'bv')]['rec']
    fold = lambda r, v: r['decision'][v] / base[v]        # noqa: E731
    out = {
        'relay_k3_fold': fold(rel, 'k_3'),
        'relay_k6_fold': fold(rel, 'k_6'),
        'relay_gly_fold': fold(rel, 'glycolysis'),
        'relay_hits': {v: rel['hits'].get(v) for v in ('k_13', 'k_17')},
        'relay_Phi_M': rel['Phi_M'],
        'relay_within_budget': bool(rel['Phi_M'] < facts['budget']),
        'unin_k3_fold': fold(uni, 'k_3'),
        'unin_hits': {v: uni['hits'].get(v)
                      for v in ('k_3', 'k_13', 'ehrlich_downstream')},
        'unin_Phi_M': uni['Phi_M'], 'unin_growth': uni['growth'],
        'iy_sims': {'returned': iy_r['sim'], 'best_visit': iy_b['sim']},
        'iy_yield_gain_pct': 100.0 * (iy_r['ibo_yield'] / iy_b['ibo_yield']
                                      - 1.0),
        'iy_irr_drop_pts': iy_b['irr_pct'] - iy_r['irr_pct'],
        'it_ip_ret_k6_fold_max': max(fold(rows[(k, 'ret')], 'k_6')
                                     for k in ('it', 'ip')),
        'it_ip_ret_etoh_max': max(rows[(k, 'ret')]['rec']['etoh']
                                  for k in ('it', 'ip')),
        'it_ip_bv_k6_fold': {k: fold(rows[(k, 'bv')], 'k_6')
                             for k in ('it', 'ip')},
        'n_rows': sum(1 for t, _ in items if t == 'row'),
        # quoted in the caption, asserted through EXPECTED / above
        'iy_r': iy_r, 'iy_b': iy_b, 'relay': rel, 'unin': uni,
        'it_ip_bv_ibo': {k: rows[(k, 'bv')]['rec']['ibo']
                         for k in ('it', 'ip')},
        'it_ip_bv_etoh': {k: rows[(k, 'bv')]['rec']['etoh']
                          for k in ('it', 'ip')},
    }
    bad = _check(S2_EXPECTED, out)
    if bad:
        raise C.FactsMismatch('Fig. S2 facts differ from S2_EXPECTED:\n  '
                              + '\n  '.join(bad))
    return out


# %% Rotated-label overlap check -------------------------------------------------------------
def _rotated_corners(t, renderer):
    """Display-space corners of a rotated Text (rotation_mode 'anchor')."""
    rot = t.get_rotation()
    t.set_rotation(0.0)
    bb = t.get_window_extent(renderer)
    t.set_rotation(rot)
    ax, ay = t.get_transform().transform(t.get_position())
    th = np.deg2rad(rot)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    pts = np.array([[bb.x0, bb.y0], [bb.x1, bb.y0], [bb.x1, bb.y1],
                    [bb.x0, bb.y1]]) - (ax, ay)
    return pts @ R.T + (ax, ay)


def _sat_overlap(p, q, tol=0.5):
    for poly in (p, q):
        for i in range(4):
            e = poly[(i + 1) % 4] - poly[i]
            n = np.array([-e[1], e[0]])
            n = n / (np.hypot(*n) or 1.0)
            a, b = p @ n, q @ n
            if a.max() - b.min() <= tol or b.max() - a.min() <= tol:
                return False
    return True


def rotated_label_overlaps(fig, rotated):
    """Overlaps of rotated texts (as rotated rectangles) with each other and
    with every other drawn text (axis-aligned boxes), plus canvas bounds.
    _style.text_overlaps compares axis-aligned boxes, which flags parallel
    rotated labels falsely; these labels are exempted there and checked
    here instead."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    polys = [(t, _rotated_corners(t, renderer)) for t in rotated]
    ids = {id(t) for t in rotated}
    others = []
    for t, bb, _, w in S._drawn_texts(fig):
        if id(t) not in ids:
            others.append((t, np.array([[bb.x0, bb.y0], [bb.x1, bb.y0],
                                        [bb.x1, bb.y1], [bb.x0, bb.y1]])))
    msgs = []
    Wp, Hp = fig.bbox.width, fig.bbox.height
    for i, (t, p) in enumerate(polys):
        if (p[:, 0].min() < 0 or p[:, 1].min() < 0 or p[:, 0].max() > Wp
                or p[:, 1].max() > Hp):
            msgs.append(f'rotated label outside canvas: {t.get_text()!r}')
        for t2, q in polys[i + 1:] + others:
            if _sat_overlap(p, q):
                msgs.append(f'rotated overlap: {t.get_text()!r} x '
                            f'{t2.get_text()!r}')
    return msgs


def cell_text_fit(fig, cells, margin_px=0.0):
    """Every cell's text must lie inside its cell's free area (inside the
    bound-hit outline where there is one), with `margin_px` to spare."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    msgs = []
    for txts, (x0, y0, x1, y1) in cells:
        tr = txts[0].get_transform()
        (X0, Y0), (X1, Y1) = tr.transform([(x0, y0), (x1, y1)])
        bbs = [t.get_window_extent(renderer) for t in txts]
        bx0, bx1 = min(b.x0 for b in bbs), max(b.x1 for b in bbs)
        by0, by1 = min(b.y0 for b in bbs), max(b.y1 for b in bbs)
        if (bx0 < X0 + margin_px or bx1 > X1 - margin_px
                or by0 < Y0 + margin_px or by1 > Y1 - margin_px):
            msgs.append('cell text does not fit: '
                        + ''.join(t.get_text() for t in txts)
                        + f' (text {bx1 - bx0:.0f} x {by1 - by0:.0f} px, '
                        f'cell {X1 - X0:.0f} x {Y1 - Y0:.0f} px)')
    return msgs


# %% Drawing --------------------------------------------------------------------------------------
def layout_rows(items):
    """[(item, y_bottom, height)] top to bottom; returns (rows, top)."""
    n_rows = sum(1 for t, _ in items if t == 'row')
    n_head = len(items) - n_rows
    top = Y_BOT + n_rows * RH + n_head * HEADER_GAP
    y, out = top, []
    for t, r in items:
        h = RH if t == 'row' else HEADER_GAP
        y -= h
        out.append((t, r, y, h))
    return out, top


def column_x():
    """{var: (x_left, block key)} and {block key: (x0, x1)}."""
    xs, spans, x = {}, {}, X0
    for b in BLOCKS:
        x_start = x
        for col in b.columns:
            xs[col.var] = (x, b.key)
            x += CW
        spans[b.key] = (x_start, x)
        x += BLOCK_GAP
    return xs, spans


def draw(facts, items):
    fig = S.new_figure(W, H)
    ax = S.inch_axes(fig, 0.0, 0.0, W, H)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis('off')
    base = facts['baseline_decision']
    rows, top = layout_rows(items)
    xs, spans = column_x()
    x_end = spans['feed'][1]
    fam_color = {'profit': S.TEXT, 'etoh': S.PALETTE['etoh_dark'],
                 'ibo': S.PALETTE['ibo_dark']}
    cells = []                  # [(cell texts, free area in inches)]
    table = {'irr': [], 'obj': [], 'phi': [], 'growth': []}

    # --- rows: labels, markers, trial numbers, cells, right-hand table
    for t, r, y, h in rows:
        yc = y + h / 2
        if t == 'header':
            ax.text(X_GROUP, yc, C.FAMILY_LABEL[r], fontsize=FS_ROW,
                    fontweight='bold', color=fam_color[r], ha='left',
                    va='center')
            continue
        indent = 0.10 if r['kind'] == 'bv' else 0.0
        mk = dict(r['marker'])
        mk.setdefault('ms', 6.0)
        ax.add_line(Line2D([X_MARK + indent], [yc], ls='none', **mk))
        ax.text(X_LABEL + indent, yc, r['label'], fontsize=FS_ROW,
                fontstyle=r['style'], color=r['color'], ha='left',
                va='center')
        if r['sim'] is not None:
            ax.text(X_TRIAL, yc, C.fmt_int(r['sim']), fontsize=FS_HEAD,
                    color=S.NOTE, ha='right', va='center')
        for b in BLOCKS:
            for col in b.columns:
                x, _ = xs[col.var]
                v = r['decision'][col.var]
                shown = r['n_spikes'] if col.var == 'max_n_spikes' else v
                z = scale_value(b.key, col.var, shown, base)
                rgb = CMAPS[b.key](NORMS[b.key](z))
                ax.add_patch(Rectangle(
                    (x + CELL_PAD / 2, y + CELL_PAD / 2), CW - CELL_PAD,
                    h - CELL_PAD, facecolor=rgb, edgecolor='none', zorder=1))
                tc = 'white' if luma(rgb) < 0.45 else S.TEXT
                if col.var == 'max_n_spikes':
                    sub = S.NOTE if tc == S.TEXT else '0.85'
                    txts = [ax.text(x + CW / 2 + 0.01, yc, f'{shown:.0f}',
                                    fontsize=FS_CELL, color=tc, ha='right',
                                    va='center', zorder=3),
                            ax.text(x + CW / 2 + 0.01, yc, f'/{v:.0f}',
                                    fontsize=FS_CELL, color=sub, ha='left',
                                    va='center', zorder=3)]
                else:
                    txts = [ax.text(x + CW / 2, yc,
                                    cell_text(b.key, col.var, v, base),
                                    fontsize=FS_CELL, color=tc, ha='center',
                                    va='center', zorder=3)]
                d = CELL_PAD / 2     # outline stroke centred on the cell edge
                inner = d
                if r['kind'] != 'base' and col.var in r['hits']:
                    ax.add_patch(Rectangle(
                        (x + d, y + d), CW - 2 * d, h - 2 * d,
                        facecolor='none', edgecolor='black', lw=OUTLINE_LW,
                        zorder=4))
                    inner = d + OUTLINE_LW / 72.0 / 2
                cells.append((txts, (x + inner, y + inner, x + CW - inner,
                                     y + h - inner)))
        # right-hand table (x positions set by place_table once measured)
        irr = C.fmt_irr(r['irr_pct'])
        above = bool(np.any(C.above_plateau(r['irr_pct'] / 100.0,
                                            facts['U'])))
        table['irr'].append(ax.text(
            0, yc, irr, fontsize=FS_ROW, ha='right', va='center',
            color=S.NOTE if irr == 'loss' else S.TEXT,
            fontweight='bold' if above else 'normal'))
        table['obj'].append(ax.text(0, yc, r['obj'], fontsize=FS_ROW,
                                    ha='left', va='center', color=S.TEXT))
        table['phi'].append(ax.text(
            0, yc, f'{r["Phi_M"]:.3f}', fontsize=FS_ROW, ha='right',
            va='center', color=S.TEXT if r['Phi_M'] > facts['budget']
            else S.NOTE))
        derated = r['growth'] < 0.99         # main-figure tag rule
        table['growth'].append(ax.text(
            0, yc, fmt_growth(r['growth']), fontsize=FS_ROW, ha='right',
            va='center', fontweight='bold' if derated else 'normal',
            color=S.TEXT if r['growth'] < 0.9995 else S.NOTE))

    # --- column labels (rotated), block titles, right-table headers
    rotated = []
    y_lab = top + 0.05
    for b in BLOCKS:
        for col in b.columns:
            x, _ = xs[col.var]
            rotated.append(ax.text(
                x + CW / 2 - 0.03, y_lab, col.label, fontsize=FS_HEAD,
                rotation=LABEL_ROT, rotation_mode='anchor', ha='left',
                va='center', color=S.TEXT))
    y_head = top + 0.06
    ax.text(X_TRIAL, y_head, 'trial', fontsize=FS_HEAD, color=S.NOTE,
            ha='right', va='bottom')
    for key, head, ha in (('irr', 'IRR', 'right'),
                          ('obj', 'own objective', 'left'),
                          ('phi', 'Φ$_\\mathrm{M}$', 'right'),
                          ('growth', 'growth', 'right')):
        table[key].append(ax.text(0, y_head, head, fontsize=FS_HEAD,
                                  fontweight='bold', ha=ha, va='bottom'))
    place_table(fig, table, x_end + TABLE_GAP, W - RIGHT_MARGIN)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lab_top = max(_rotated_corners(t, renderer)[:, 1].max() for t in rotated)
    y_title = ax.transData.inverted().transform((0, lab_top))[1] + 0.05
    for b in BLOCKS:
        x0, x1 = spans[b.key]
        ax.text((x0 + x1) / 2, y_title, b.title, fontsize=FS_BLOCK,
                fontweight='bold', ha='center', va='bottom', color=S.TEXT)

    # --- key (top left): marker fill, bound outline
    ky = y_title + 0.09
    S.inline_key(fig, X_GROUP, ky, [
        dict(marker='o', ms=6.0, mfc=S.TEXT, mec=S.TEXT,
             text='returned design')], fontsize=FS_HEAD)
    S.inline_key(fig, X_GROUP, ky - KEY_PITCH, [
        dict(marker='o', ms=6.0, mfc='white', mec=S.TEXT, mew=1.3,
             text='best visit (highest IRR)', style='italic')],
        fontsize=FS_HEAD)
    S.inline_key(fig, X_GROUP - 0.005, ky - 2 * KEY_PITCH, [
        dict(swatch='white', edgecolor='black', lw=OUTLINE_LW, size_in=0.11,
             text='at a search bound')], fontsize=FS_HEAD, gap_in=0.04)

    # --- table notes (top right, right-aligned)
    xr = W - RIGHT_MARGIN
    ax.text(xr, ky, f'bold IRR: above the plateau '
            f'({C.fmt_pct(facts["U_pct"])})', fontsize=FS_HEAD,
            color=S.NOTE, ha='right', va='center')
    ax.text(xr, ky - KEY_PITCH, f'growth derated at Φ$_\\mathrm{{M}}$ > '
            f'{facts["budget"]:.4f}', fontsize=FS_HEAD, color=S.NOTE,
            ha='right', va='center')

    # --- colour bars under the heatmap, one per block
    for b in BLOCKS:
        x0, x1 = spans[b.key]
        ax.text((x0 + x1) / 2, Y_BOT - 0.10, b.descriptor, fontsize=FS_HEAD,
                color=S.NOTE, ha='center', va='center')
        cax = S.inch_axes(fig, x0, CBAR_Y, x1 - x0, CBAR_H)
        draw_colorbar(cax, b.key)
    return fig, rotated, cells


def place_table(fig, table, x0, x1, min_gap=0.10):
    """Place the right-hand table's columns (dict key -> Texts, in order)
    between x0 and x1 inches: each column is as wide as its widest text and
    the slack is shared equally between the gaps. Raises if a gap would be
    narrower than `min_gap` inches."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    widths = [max(t.get_window_extent(renderer).width for t in txts)
              / fig.dpi for txts in table.values()]
    gap = (x1 - x0 - sum(widths)) / (len(widths) - 1)
    if gap < min_gap:
        raise S.FigureCheckError(f'right table needs {sum(widths):.2f} in + '
                                 f'gaps; only {x1 - x0:.2f} in available')
    cursor = x0
    for txts, w in zip(table.values(), widths):
        for t in txts:
            t.set_x(cursor + w if t.get_ha() == 'right' else cursor)
        cursor += w + gap
    print(f'right table: column widths {[round(w, 2) for w in widths]} in, '
          f'gap {gap:.3f} in')
    return gap


def draw_colorbar(cax, key):
    """A thin horizontal colour bar of a block's scale, 9-pt ticks."""
    cmap, norm = CMAPS[key], NORMS[key]
    if key == 'native':
        lo, hi = -FOLD_CLIP, FOLD_CLIP
        ticks = [-4, -2, 0, 2, 4]
        labels = ['1/16', '1/4', '1', '4', '16']
    elif key == 'path':
        lo, hi = 0.0, PATH_LIM
        ticks = [0, 2, 4]
        labels = ['0', '2', '4']
    elif key == 'inhib':
        blo, bhi = C.BANDS['inhib_ethanol']
        lo, hi = np.log2(blo), np.log2(bhi)
        ticks = [lo, 0.0, hi]
        labels = [f'{blo:g}', '1', f'{bhi:g}']
    else:
        lo, hi = 0.0, 1.0
        ticks = [0.0, 1.0]
        labels = ['min', 'max']
    z = np.linspace(lo, hi, 256)
    cax.imshow(cmap(norm(z))[None, :, :], aspect='auto',
               extent=(lo, hi, 0, 1), interpolation='bilinear')
    cax.set_xlim(lo, hi)
    cax.set_ylim(0, 1)
    cax.set_xticks(ticks)
    cax.set_xticklabels(labels)
    cax.tick_params(axis='x', labelsize=FS_HEAD, pad=2)
    tl = cax.get_xticklabels()          # end labels aligned inward, so the
    tl[0].set_ha('left')                # labels of adjacent bars never meet
    tl[-1].set_ha('right')
    cax.set_yticks([])
    S.style_ticks(cax, y=False, minor_x=False)
    cax.tick_params(axis='x', which='major', top=False, length=3)
    for sp in cax.spines.values():
        sp.set_linewidth(0.6)


# %% Caption ------------------------------------------------------------------------------------
def write_caption(facts, f2, path=CAPTION_PATH):
    """Write the caption draft (Markdown) with every number taken from the
    facts / S2 facts asserted above, so caption and figure never disagree."""
    base = facts['baseline_decision']
    rel, uni = f2['relay'], f2['unin']
    iy_r, iy_b = f2['iy_r'], f2['iy_b']
    budget = facts['budget']
    bv_k6 = sorted({f'{v:.1f}' for v in f2['it_ip_bv_k6_fold'].values()})
    bv_k6 = bv_k6[0] if len(bv_k6) == 1 else '–'.join(bv_k6)
    lo_i, hi_i = C.BANDS['inhib_ethanol']
    et, ib = f2['it_ip_bv_etoh'], f2['it_ip_bv_ibo']
    ph = 'Φ<sub>M</sub>'
    gl = 'g·L<sup>−1</sup>'
    glh = 'g·L<sup>−1</sup>·h<sup>−1</sup>'
    lines = [
        '**Figure S2 | Design fingerprints of the returned designs, and of '
        "the isobutanol scouts' best visits.**",
        '',
        '* **Rows.** The starting strain (scenario A); the design each '
        'campaign returned (argmax of its own objective; filled markers); '
        'and, in italics with open markers, the highest-IRR trial ("best '
        'visit") of each isobutanol TRY scout. Trial numbers are simulated '
        'indices.',
        '* **Columns.** Native capacities are fold change vs the starting '
        f'strain (Pdc k<sub>3</sub> {base["k_3"]:.2f}, Adh1 k<sub>6</sub> '
        f'{base["k_6"]:.2f} and Adh6 k<sub>17</sub> {base["k_17"]:.4f} {glh}; '
        f'glycolysis multiplier {base["glycolysis"]:.0f}); the colour is '
        'clipped at 1/16× and 16×. The isobutanol-pathway capacities (ALS '
        'k<sub>13</sub> and the Ilv5/Ilv3/Aro10 group, "Ilv5–Aro10") are '
        f'absolute, in {glh}. Inhibition multipliers scale the ethanol, '
        'isobutanol and acetate product-inhibition coefficients (search range '
        f'{lo_i:g}–{hi_i:g}; 1 = starting strain). Feeding variables are '
        f'shaded within their search ranges (threshold and target Δ in {gl}; '
        'spikes = actual / cap, shaded by the actual count). Outlined cells '
        'sit at a search bound.',
        '* **Right.** IRR ("loss" = negative or no IRR; bold = above the '
        f'uninformed campaign\'s plateau, {C.fmt_pct(facts["U_pct"])}); the '
        "campaign's own objective; the metabolic proteome "
        f'{ph} [g·(g DCW)<sup>−1</sup>]; and the growth factor (growth is '
        f'derated above {ph} = {budget:.4f}, the penalty-free budget). The '
        f'starting strain\'s IRR ({C.fmt_pct(facts["start_irr_pct"])}) is '
        'under the model version the campaigns ran with.',
        '* **Returned strains differ in products and proteome.** The '
        'profitable TRY-informed co-production design throttles the entry to '
        f'the ethanol branch (Pdc {fmt_fold(f2["relay_k3_fold"])}, Adh1 '
        f'{fmt_fold(f2["relay_k6_fold"])}) at baseline-like glycolysis '
        f'({fmt_fold(f2["relay_gly_fold"])}), puts ALS and Adh6 at their '
        f'upper bounds, and stays within the proteome '
        f'budget ({ph} {rel["Phi_M"]:.3f}, growth ×{rel["growth"]:.2f}). The '
        f'uninformed design puts Pdc at its {fmt_fold(f2["unin_k3_fold"])} '
        'bound and switches the isobutanol pathway off (ALS and '
        f'Ilv5/Ilv3/Aro10 at their floors); its {ph} ({uni["Phi_M"]:.3f}) '
        f'exceeds the budget (growth ×{uni["growth"]:.2f}).',
        "* **Returned vs best visit.** The isobutanol-yield scout's returned "
        f'design (trial {C.fmt_int(iy_r["sim"])}) and its best visit (trial '
        f'{C.fmt_int(iy_b["sim"])}) are consecutive trials. The '
        f'{f2["iy_yield_gain_pct"]:.1f} % higher yield it was rewarded for '
        f'({iy_r["ibo_yield"]:.3f} vs {iy_b["ibo_yield"]:.3f} '
        f'g·g<sup>−1</sup>) cost {f2["iy_irr_drop_pts"]:.0f} IRR points '
        f'({C.fmt_irr(iy_r["irr_pct"], 2)} vs '
        f'{C.fmt_irr(iy_b["irr_pct"], 2)}): TCI {iy_r["TCI"]:.0f} vs '
        f'{iy_b["TCI"]:.0f} MM$, batch {iy_r["tau"]:.0f} vs '
        f'{iy_b["tau"]:.0f} h, growth ×{iy_r["growth"]:.2f} vs '
        f'×{iy_b["growth"]:.2f}. The isobutanol titer and productivity scouts '
        'returned designs with Adh1 knocked down (≤ '
        f'{f2["it_ip_ret_k6_fold_max"]:.3f}× the starting strain) that make ≤ '
        f'{f2["it_ip_ret_etoh_max"]:.1f} {gl} ethanol and have no IRR; their '
        f'best visits keep Adh1 at {bv_k6}× and make {et["it"]:.1f} and '
        f'{et["ip"]:.1f} {gl} ethanol alongside {ib["it"]:.1f} and '
        f'{ib["ip"]:.1f} {gl} isobutanol.',
        '* The starting strain is shown for reference: its zero '
        'isobutanol-pathway capacities and its target Δ '
        f'({base["target_delta"]:.1f} {gl}) lie below the search floors '
        f'({C.BANDS["k_13"][0]:g} {glh} and '
        f'{C.BANDS["target_delta"][0]:g} {gl}), so its cells carry no '
        'outline.',
    ]
    text = '\n'.join(lines) + '\n'
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(text)
    print(f'caption written: {path}')
    return text


# %% Main --------------------------------------------------------------------------------------
def main():
    facts = C.check_facts()
    S.cvd_check()
    items = design_rows(facts)
    f2 = s2_facts(facts, items)
    fig, rotated, cells = draw(facts, items)
    probs = rotated_label_overlaps(fig, rotated) + cell_text_fit(fig, cells)
    for m in probs:
        S._safe_print(f'[s2] {m}')
    S.check_figure(fig, size=S.S2_SIZE, exempt=rotated)
    if probs:
        raise S.FigureCheckError(f'{len(probs)} figure-local problem(s)')
    hits = C.literal_scan([os.path.abspath(__file__)])
    if hits:
        raise C.FactsMismatch(f'hard-coded annotation literals: {hits}')
    write_caption(facts, f2)
    S.save(fig, STEM)
    C.assert_sim_safe()
    print('Fig. S2 done: ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
