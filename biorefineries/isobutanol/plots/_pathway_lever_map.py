#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Enzyme lever-map figure -- the fourth variant of the kinetic-optimization
parameter-set comparison figure (companion to plot_kin_opt_parameter_sets.py).

Re-presents the metabolic_split_12d decision variables ENZYME-FIRST, laid on
the glucose -> ethanol / isobutanol carbon backbone, each lever shown as a
fold-change from wild type (or an absolute newly-expressed rate), and closes
to profitability with an outcome ribbon (max-product != max-profit).

Sim-safe: this module imports only matplotlib and numpy. The parent passes in
its already-file-path-loaded ko / eb modules and a bundle of its own
callables/tables (`helpers`); nothing here loads a module by path or imports
the biorefineries package. See plot_pathway_lever_map."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

# --- pathway node -> lever mapping (metabolic_split_12d layout) --------------
# step id -> (record key, kind). 'fold' = a x-baseline fold-change lever
# (baseline 1.0: glycolysis group + the k_*_rel fermentation/Aro10/Adh6 cells);
# 'absolute' = a newly-expressed Ehrlich enzyme, off (0) in wild type, drawn in
# g/L/h. Keys match what load_set / baseline_set expose under use_split_12d_layout
# (SPLIT_12D_REL_RATE_VARS -> k_*_rel; SPLIT_12D_RATE_VARS -> k_13/k_14/k_15).
NODE_LEVERS = {
    'r1':  ('glycolysis', 'fold'),
    'r3':  ('k_3_rel', 'fold'),
    'r6':  ('k_6_rel', 'fold'),
    'r13': ('k_13', 'absolute'),
    'r14': ('k_14', 'absolute'),
    'r15': ('k_15', 'absolute'),
    'r16': ('k_16_rel', 'fold'),
    'r17': ('k_17_rel', 'fold'),
}
# product-inhibition multipliers (baseline 1.0, fold-change). Larger multiplier
# => larger exp(-k*[product]) exponent in the nskinetics rate laws (r1/r4/r6/
# r7/r10/r17) => MORE inhibition => LESS tolerant, so a bar above 1x reads as
# less tolerant and engineered tolerance rises as the multiplier drops below 1x.
TOLERANCE_LEVERS = ('inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate')
# outcome-ribbon columns (the punchline): financial + the two product titers.
OUTCOME_COLS = ('IRR', 'IBO titer', 'EtOH titer')


# --- pure-logic value helpers (unit-tested offline) -------------------------
def lever_value(rec, key):
    """The lever's plotted value for set record `rec`: rec[key] (a fold-change
    for a 'fold' lever, g/L/h for an 'absolute' lever), or nan when the key is
    absent or the value is non-finite."""
    v = rec.get(key)
    try:
        v = float(v)
    except (TypeError, ValueError):
        return np.nan
    return v if np.isfinite(v) else np.nan


def is_off_in_wildtype(baseline_rec, key):
    """True if an absolute-lever enzyme is genuinely off (0 or absent) in the
    wild type, so its 'off in wild type' marker is drawn and its baseline bar
    suppressed."""
    v = lever_value(baseline_rec, key)
    return not (np.isfinite(v) and v > 0.0)


def clamp_outcome(col, value, floor, clamp_neg):
    """Outcome-ribbon value with the parent's CLAMP_NEG_TO_ZERO treatment: for
    a clamp column (IRR) a finite loss or -inf is drawn at `floor`; nan (never
    solved) stays nan so the caller omits that bar."""
    v = float(value)
    if np.isnan(v):
        return np.nan
    if col in clamp_neg:
        if not np.isfinite(v) or v < floor:   # -inf or a finite loss -> floor
            return floor
    if not np.isfinite(v):
        return np.nan
    return v


def proteome_segments(rec, eb, categories):
    """Reduce a record's burden pools to the slim stacked-strip segments:
    the isobutanol-pathway pool, the rest of the modeled metabolic pool, the
    flexible slack, the derated translation sector, and housekeeping. Sums to
    eb.PROTEIN_CONTENT. `categories` is helpers['BURDEN_CATEGORIES'] (the parent's
    ((name, [steps]), ...) partition of Phi_M)."""
    PC = float(eb.PROTEIN_CONTENT)
    housekeeping = PC * float(eb.HOUSEKEEPING_FRACTION)
    ibo_steps = next(steps for name, steps in categories
                     if name == 'Isobutanol production')
    ibo = sum(float(rec[f'pool_{st}']) for st in ibo_steps)
    Phi_M = float(rec['Phi_M'])
    rest = max(0.0, Phi_M - ibo)
    phi_T_built = float(rec['burden_factor']) * float(rec['phi_T'])
    slack = max(0.0, PC - housekeeping - Phi_M - phi_T_built)
    return {'isobutanol': ibo, 'rest': rest, 'slack': slack,
            'translation': phi_T_built, 'housekeeping': housekeeping}


# --- drawing ----------------------------------------------------------------
# pathway node coordinates in the pathway axis' 0..1 data space (top-down; the
# pyruvate split feeds the ethanol arm on the left and the Ehrlich/isobutanol
# arm on the right). Tunable in the Task-5 visual pass.
_NODE_POS = {
    'r1':  (0.30, 0.90),   # glycolysis lump
    'r3':  (0.30, 0.66),   # Pdc1
    'r6':  (0.30, 0.42),   # Adh1
    'r13': (0.72, 0.78),   # Ilv2+Ilv6
    'r14': (0.72, 0.62),   # Ilv5
    'r15': (0.72, 0.46),   # Ilv3
    'r16': (0.72, 0.30),   # Aro10
    'r17': (0.72, 0.14),   # Adh6
}
# metabolite marker positions (label, x, y, emphasise). The isobutanol-arm
# intermediates sit in the gap BELOW each enzyme's glyph cluster and ABOVE the
# next box, clear of the upward-growing absolute bars.
_METABOLITES = [
    ('Glucose',            0.30, 0.99, False),
    ('Pyruvate',           0.51, 0.79, False),
    ('Acetaldehyde',       0.30, 0.54, False),
    ('ETHANOL',            0.30, 0.19, True),
    ('2-acetolactate',     0.72, 0.685, False),
    ('2,3-dihydroxy-\nisovalerate', 0.72, 0.520, False),
    ('KIV',                0.72, 0.360, False),
    ('Isobutyraldehyde',   0.72, 0.200, False),
    ('ISOBUTANOL',         0.72, 0.03, True),
]
_NODE_W, _NODE_H = 0.19, 0.066   # enzyme box size in pathway-axis data coords


def _cluster(ax, cx, cy, w, h, values, cols, kind, off_mask=None):
    """A mini bar cluster in axis data coords, centred at (cx, cy), width w,
    half-height h. For kind 'fold' the bars deflect from a 1x reference line at
    cy (up = over-expression, down = knockdown) on a symmetric log2 scale
    clipped to +-2 (0.25x..4x); for kind 'absolute' the bars grow up from cy-h
    scaled to the cluster's own max. `off_mask[j]` True suppresses bar j (a
    baseline 'off in wild type' cell)."""
    n = len(values)
    bw = w / max(n, 1) * 0.72
    xs = np.linspace(cx - w / 2 + bw / 2, cx + w / 2 - bw / 2, n)
    if kind == 'fold':
        ax.plot([cx - w / 2, cx + w / 2], [cy, cy], color='0.35', lw=0.6,
                zorder=3)   # 1x reference
        for x, v, c, j in zip(xs, values, cols, range(n)):
            if off_mask is not None and off_mask[j]:
                continue
            if not np.isfinite(v) or v <= 0:
                continue
            s = np.clip(np.log2(v), -2.0, 2.0) / 2.0     # -1..1
            ax.add_patch(Rectangle((x - bw / 2, cy), bw, s * h,
                                   facecolor=c, edgecolor='none', zorder=4))
    else:  # absolute
        finite = [v for v in values if np.isfinite(v) and v > 0]
        vmax = max(finite) if finite else 1.0
        base = cy - h
        for x, v, c, j in zip(xs, values, cols, range(n)):
            if off_mask is not None and off_mask[j]:
                ax.plot(x, base, marker='x', ms=4, color=c, zorder=4)  # off-in-WT
                continue
            if not np.isfinite(v) or v <= 0:
                continue
            ax.add_patch(Rectangle((x - bw / 2, base), bw, (v / vmax) * 2 * h,
                                   facecolor=c, edgecolor='none', zorder=4))


def _enzyme_node(ax, step, cx, cy, sets, cols, helpers):
    """A rounded enzyme box (name / reaction id / parameter symbol) with its
    mini fold-change or absolute glyph cluster below the labels."""
    _mathify = helpers['_mathify']
    name = helpers['STEP_ENZYME'][step]
    param = helpers['STEP_PARAMS'][step]
    key, kind = NODE_LEVERS[step]
    ax.add_patch(FancyBboxPatch(
        (cx - _NODE_W / 2, cy - _NODE_H / 2), _NODE_W, _NODE_H,
        boxstyle='round,pad=0.006,rounding_size=0.012',
        facecolor='white', edgecolor='0.25', lw=0.9, zorder=3))
    ax.text(cx, cy + _NODE_H * 0.28, name, ha='center', va='center',
            fontsize=helpers['FONTS']['tick'], fontweight='bold', zorder=5)
    ax.text(cx, cy - _NODE_H * 0.02,
            f'{step}  ' + _mathify(param.split(",")[0].strip()),
            ha='center', va='center', fontsize=helpers['FONTS']['tick'] - 1,
            color='0.25', zorder=5)
    base = sets[0]
    values = [lever_value(s, key) for s in sets]
    off = [is_off_in_wildtype(base, key) and s.get('is_baseline') and
           kind == 'absolute' for s in sets] if kind == 'absolute' else None
    _cluster(ax, cx, cy - _NODE_H * 0.62, _NODE_W * 0.9, _NODE_H * 0.42,
             values, [cols[id(s)] for s in sets], kind, off_mask=off)


def _arrow(ax, p0, p1, **kw):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=9,
                                 lw=1.0, color=kw.pop('color', '0.4'),
                                 shrinkA=6, shrinkB=6, zorder=2, **kw))


def _draw_pathway(ax, sets, cols, helpers):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    # metabolite markers
    for label, x, y, emph in _METABOLITES:
        ax.text(x, y, label, ha='center', va='center',
                fontsize=helpers['FONTS']['tick'] - (0 if emph else 1),
                fontweight='bold' if emph else 'normal',
                color='0.1' if emph else '0.35', zorder=5)
    # ethanol arm arrows: glucose -> r1 -> pyruvate -> r3 -> acetald -> r6 -> EtOH
    _arrow(ax, (0.30, 0.965), _NODE_POS['r1'])
    _arrow(ax, _NODE_POS['r1'], (0.44, 0.80))       # r1 -> pyruvate
    _arrow(ax, (0.51, 0.74), _NODE_POS['r3'])       # pyruvate -> r3
    _arrow(ax, _NODE_POS['r3'], (0.30, 0.575))
    _arrow(ax, (0.30, 0.515), _NODE_POS['r6'])
    _arrow(ax, _NODE_POS['r6'], (0.30, 0.24))
    # isobutanol arm: pyruvate -> r13 -> ... -> r17 -> ISOBUTANOL
    _arrow(ax, (0.58, 0.78), _NODE_POS['r13'])
    for a, b in (('r13', 'r14'), ('r14', 'r15'), ('r15', 'r16'), ('r16', 'r17')):
        _arrow(ax, _NODE_POS[a], _NODE_POS[b])
    _arrow(ax, _NODE_POS['r17'], (0.72, 0.06))
    for step, (cx, cy) in _NODE_POS.items():
        _enzyme_node(ax, step, cx, cy, sets, cols, helpers)
    # faint dashed feedback arrows from the product pools back to the pathway
    for prod in ((0.30, 0.20), (0.72, 0.03)):
        ax.add_patch(FancyArrowPatch(prod, (0.50, 0.50), arrowstyle='-|>',
                                     mutation_scale=7, lw=0.7, color='0.7',
                                     ls=(0, (2, 2)), shrinkA=8, shrinkB=8,
                                     zorder=1))


def _draw_process_box(ax, sets, cols, helpers):
    """Fermentor-side box: the three fed-batch feeding levers as ABSOLUTE bars
    (threshold + target sugar on a shared 0-300 axis; spike count on its own)."""
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    ax.add_patch(FancyBboxPatch((0.05, 0.08), 0.9, 0.84,
                                boxstyle='round,pad=0.01,rounding_size=0.03',
                                facecolor='0.96', edgecolor='0.4', lw=1.0))
    ax.text(0.5, 0.965, 'Process levers\n(fed-batch feeding)', ha='center',
            va='top', fontsize=helpers['FONTS']['cell'], fontweight='bold')
    rows = [('threshold_conc', 300.0), ('target_conc', 300.0),
            ('n_glu_spikes', 50.0)]
    for r, (key, vmax) in enumerate(rows):
        y = 0.63 - r * 0.25
        title = helpers['FEED_LABELS'][key][0].replace('\n', ' ')
        ax.text(0.10, y + 0.155, title, ha='left', va='bottom',
                fontsize=helpers['FONTS']['tick'] - 1, color='0.2')
        vals = [lever_value(s, key) for s in sets]
        n = len(sets)
        xs = np.linspace(0.14, 0.86, n)
        bw = (xs[1] - xs[0]) * 0.7 if n > 1 else 0.1
        for x, v, s in zip(xs, vals, sets):
            if not np.isfinite(v) or v <= 0:
                continue
            ax.add_patch(Rectangle((x - bw / 2, y), bw, (v / vmax) * 0.13,
                                   facecolor=cols[id(s)], edgecolor='none'))


def _draw_tolerance(ax, sets, cols, helpers):
    """Cellular-tolerance lever group: the three product-inhibition multipliers
    as x-baseline clusters (baseline 1.0). A larger multiplier = MORE inhibition
    = LESS tolerant, so annotate engineered tolerance toward < 1x."""
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    ax.text(0.5, 0.95, 'Cellular tolerance levers (product inhibition)',
            ha='center', va='top', fontsize=helpers['FONTS']['cell'],
            fontweight='bold')
    labels = ('Ethanol', 'Isobutanol', 'Acetate')
    xs = np.linspace(0.2, 0.8, len(TOLERANCE_LEVERS))
    for x, key, lab in zip(xs, TOLERANCE_LEVERS, labels):
        vals = [lever_value(s, key) for s in sets]
        _cluster(ax, x, 0.5, 0.22, 0.24, vals, [cols[id(s)] for s in sets],
                 'fold')
        ax.text(x, 0.16, lab, ha='center', va='top',
                fontsize=helpers['FONTS']['tick'] - 1, color='0.2')
    ax.annotate('more tolerant', xy=(0.06, 0.34), xytext=(0.06, 0.66),
                ha='center', fontsize=helpers['FONTS']['tick'] - 1, color='0.4',
                arrowprops=dict(arrowstyle='-|>', color='0.4', lw=0.8))


def _draw_proteome_strip(ax, sets, cols, eb, helpers):
    """One slim horizontal stacked bar per design summing to eb.PROTEIN_CONTENT:
    isobutanol-pathway pool | rest of the modeled pool | slack | translation |
    housekeeping -- the causal 'why' between levers and outcome."""
    PC = float(eb.PROTEIN_CONTENT)
    order = ['isobutanol', 'rest', 'slack', 'translation', 'housekeeping']
    hatch = {'isobutanol': 'xxx', 'rest': None, 'slack': None,
             'translation': 'oo', 'housekeeping': '++'}
    fill = {'isobutanol': True, 'rest': True, 'slack': False,
            'translation': False, 'housekeeping': False}
    n = len(sets)
    for i, s in enumerate(sets):
        y = n - i
        seg = proteome_segments(s, eb, helpers['BURDEN_CATEGORIES'])
        c = cols[id(s)]
        x = 0.0
        for name in order:
            w = seg[name]
            if w <= 0:
                x += w; continue
            ax.barh(y, w, left=x, height=0.6,
                    facecolor=(c if fill[name] else 'white'),
                    edgecolor=(c if not fill[name] else '0.15'),
                    lw=0.5, hatch=hatch[name], zorder=2)
            x += w
    ax.axvline(PC, color='0.15', ls='--', lw=1.0)
    ax.set_xlim(0, PC * 1.02); ax.set_ylim(0.4, n + 0.6)
    ax.set_yticks([])
    ax.set_xlabel(helpers['_bold_axis_title'](
        'Proteome allocation [g protein·(g DCW)$^{-1}$]'),
        fontsize=helpers['FONTS']['cell'])
    ax.tick_params(axis='x', direction='inout', top=False, length=4)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def _draw_ribbon(fig, gs, sets, cols, helpers):
    """Outcome ribbon: IRR (%), isobutanol titer, ethanol titer -- one grouped
    bar cluster per metric, the four designs side by side; the financial
    design starred."""
    specs = [('IRR', 'Financial\n(IRR %)', 0.30, True),
             ('IBO titer', 'Isobutanol\ntiter [g·L$^{-1}$]', None, False),
             ('EtOH titer', 'Ethanol\ntiter [g·L$^{-1}$]', None, False)]
    sub = gs.subgridspec(1, 3, wspace=0.5)
    cn = helpers['CLAMP_NEG_TO_ZERO']
    n = len(sets)
    for k, (col, title, cap, is_irr) in enumerate(specs):
        ax = fig.add_subplot(sub[0, k])
        raw = [clamp_outcome(col, s.get(col, np.nan), 0.0, cn) for s in sets]
        vals = [v for v in raw if np.isfinite(v)]
        vmax = cap if cap else (max(vals) * 1.15 if vals else 1.0)
        for j, (s, v) in enumerate(zip(sets, raw)):
            if not np.isfinite(v):
                continue
            ax.bar(j, v, width=0.7, color=cols[id(s)], zorder=2)
            if s.get('objective') in ('PI', 'PI (log-tail)'):
                ax.text(j, v, '$\\bigstar$', ha='center', va='bottom',
                        fontsize=helpers['FONTS']['cell'], color='0.1')
        ax.set_ylim(0, vmax)
        ax.set_title(helpers['_bold_axis_title'](title),
                     fontsize=helpers['FONTS']['tick'])
        ax.set_xticks([])
        if is_irr:
            from matplotlib.ticker import PercentFormatter
            ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0,
                                                          symbol=''))
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)


def plot_pathway_lever_map(sets, colors, out_stem, *, ko, eb, helpers, dpi=300):
    """Render the enzyme lever-map figure and write <out_stem>.png and .pdf.

    `sets` is the baseline-first list of set records already built by the parent
    under the split_12d layout; `colors` is the objective-keyed {id(set): hex}
    map; `helpers` bundles the parent callables/tables the drawing reuses. `ko`
    is accepted for symmetry with the other companions (unused here). Returns
    out_stem."""
    helpers['apply_fonts']()
    fig = plt.figure(figsize=(12.0, 12.0))
    gs = fig.add_gridspec(
        4, 2, height_ratios=[0.8, 4.9, 1.3, 1.7], width_ratios=[1.0, 2.1],
        left=0.05, right=0.97, top=0.94, bottom=0.06, hspace=0.35, wspace=0.12)
    # title + design legend (top strip spanning both columns)
    axT = fig.add_subplot(gs[0, :]); axT.axis('off')
    axT.text(0.0, 0.7, 'Levers for isobutanol co-production: strain, process, '
             'and the price of pushing a product to its maximum',
             ha='left', va='center', fontsize=helpers['FONTS']['panel'],
             fontweight='bold')
    handles = [Rectangle((0, 0), 1, 1, fc=colors[id(s)],
                         label=s['label'] + (' $\\bigstar$'
                         if s.get('objective') in ('PI', 'PI (log-tail)')
                         else '')) for s in sets]
    axT.legend(handles=handles, loc='lower left', ncol=len(sets), frameon=False,
               fontsize=helpers['FONTS']['legend'], bbox_to_anchor=(0.0, -0.15))
    # process levers (left of the pathway) + the pathway (right)
    _draw_process_box(fig.add_subplot(gs[1, 0]), sets, colors, helpers)
    _draw_pathway(fig.add_subplot(gs[1, 1]), sets, colors, helpers)
    # tolerance group (spanning) then proteome strip then ribbon
    _draw_tolerance(fig.add_subplot(gs[2, 0]), sets, colors, helpers)
    _draw_proteome_strip(fig.add_subplot(gs[2, 1]), sets, colors, eb, helpers)
    _draw_ribbon(fig, gs[3, :], sets, colors, helpers)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out_stem}.{ext}', dpi=dpi)
    plt.close(fig)
    return out_stem
