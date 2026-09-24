#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Stage 2 of the TEA-breakdown figure of the metabolic_split_12d campaign
optima (spec docs/superpowers/specs/2026-09-24-tea-breakdowns-split12d-design.md;
stage 1, which simulates, is analyses/collect_tea_breakdowns_split12d.py).

A 3x3 grid of cost & utility breakdowns, one panel per scenario (LAYOUT):
the scenario-A baseline, the profitability and flagship (relay) campaigns,
then the isobutanol and the ethanol yield / titer / productivity campaigns,
each at its objective-optimum trial. Every panel has five 100 %-stacked bars
(the unit groups' installed equipment cost, cooling duty, heating duty,
electricity consumption and operating cost), each group drawn as its share of
the metric's POSITIVE total, so credits (the heat-exchanger-network savings,
sold excess electricity) sit below zero. Each bar's NET total -- operating
cost converted to MM$/y, i.e. the AOC -- is printed above it.

Sim-safe: reads the stage-1 JSON only (json / numpy / matplotlib); never
imports the biorefineries package, never load()s. Run:

    python plots/plot_tea_breakdowns_split12d.py [--data <json>] [--out-dir DIR]

--data defaults to the newest analyses/results/tea_breakdowns_split12d_*.json.
Writes <stem>_<stamp>.png and .pdf to --out-dir (default analyses/results).
"""
import os
import glob
import json
import math
import argparse
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import TICKDOWN, TICKLEFT
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, MultipleLocator, NullLocator

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')
DATA_GLOB = 'tea_breakdowns_split12d_*.json'
DEFAULT_STEM = 'tea_breakdowns_split12d_figure'

#: panel grid: rows by theme (baseline + the two profitability campaigns;
#: isobutanol; ethanol), columns yield / titer / productivity below row 1
LAYOUT = (('baseline', 'profitability', 'flagship'),
          ('ibo_yield', 'ibo_titer', 'ibo_productivity'),
          ('etoh_yield', 'etoh_titer', 'etoh_productivity'))

#: metric -> (bottom-row tick label, units the stage-1 JSON must carry,
#: factor to the displayed unit; operating cost USD/hr -> MM$/y is applied
#: separately because it needs the operating hours)
OPERATING = 'Operating cost'
METRIC_SPECS = {
    'Installed equipment cost': ('Installed\nequipment\ncost\n[MM\\$]', 'MM$'),
    'Cooling duty': ('Cooling\nduty\n[$\\mathrm{GJ·h}^{-1}$]', 'GJ/hr'),
    'Heating duty': ('Heating\nduty\n[$\\mathrm{GJ·h}^{-1}$]', 'GJ/hr'),
    'Electricity consumption': ('Electricity\nconsumption\n[MW]', 'MW'),
    OPERATING: ('Operating\ncost\n[$\\mathrm{MM\\$·y}^{-1}$]', 'USD/hr'),
}

#: unit group -> (face colour, hatch). Process areas take the validated
#: eight-slot categorical palette in order (adjacent stack pairs clear the CVD
#: and normal-vision floors); the two past eight repeat slots 3 and 1 with a
#: texture (never a generated ninth hue). Facilities are neutral greys, the
#: turbogenerator sharing the boiler's grey (both are BT801) under a texture.
#: The four cost pseudo-groups (operating cost only) are white + texture.
GROUP_STYLES = {
    'feedstock acquisition':              ('#2a78d6', ''),
    'feedstock saccharification':         ('#eb6834', ''),
    'sugar solution preparation':         ('#1baf7a', ''),
    'fermentation':                       ('#eda100', ''),
    'alcohol recovery':                   ('#e87ba4', ''),
    'ethanol purification':               ('#008300', ''),
    'isobutanol purification':            ('#4a3aa7', ''),
    'storage and handling':               ('#e34948', ''),
    'DDGS recovery':                      ('#1baf7a', '////'),
    'wastewater treatment':               ('#2a78d6', '////'),
    'heat exchanger network':             ('#52514e', ''),
    'boiler':                             ('#8a8984', ''),
    'turbogenerator':                     ('#8a8984', '////'),
    'cooling utility facilities':         ('#c3c2b7', ''),
    'other facilities':                   ('#c3c2b7', '....'),
    'natural gas (for steam generation)': ('#ffffff', '\\\\\\\\'),
    'natural gas (for product drying)':   ('#ffffff', '////'),
    'fixed operating cost':               ('#ffffff', 'xxxx'),
    'excess electricity':                 ('#ffffff', '....'),
}
#: legend columns: each family starts a new column (process areas fill two),
#: columns padded to LEGEND_ROWS with blank entries so families line up
GROUP_FAMILIES = (
    tuple(list(GROUP_STYLES)[:10]),     # process areas
    tuple(list(GROUP_STYLES)[10:15]),   # facilities
    tuple(list(GROUP_STYLES)[15:]),     # cost pseudo-groups
)
LEGEND_ROWS = 5
EDGE_COLOR = '0.15'
EDGE_WIDTH = 0.4
BAR_WIDTH = 0.64
#: a segment below this |share| [%] is not drawn (a zero-height bar would
#: still draw its outline as a stray line at the stack top)
MIN_DRAWN_SHARE = 1e-6
#: a group whose |share| stays below this [%] in every bar of every panel is
#: left out of the legend (it is invisible in the figure)
LEGEND_MIN_SHARE = 0.1
TOTAL_SIG_FIGS = 3
Y_TOP = 116.0          # room above 100 % for the totals
TOTAL_PAD = 2.0        # gap between the stack top (100 %) and a total [%]
Y_MAJOR, Y_MINOR = 25.0, 5.0

FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis': 12, 'title': 12, 'subtitle': 10,
         'category': 10, 'total': 9, 'legend': 9}
TICK_LEN = {'major': 4.0, 'minor': 2.0}


#%% Data
def newest_data(results_dir=RESULTS_DIR):
    paths = sorted(glob.glob(os.path.join(results_dir, DATA_GLOB)),
                   key=os.path.getmtime)
    if not paths:
        raise FileNotFoundError(f'no {DATA_GLOB} in {results_dir}: run '
                                'analyses/collect_tea_breakdowns_split12d.py')
    return paths[-1]


def load_breakdowns(path):
    """The stage-1 document, validated against what this figure assumes."""
    with open(path) as f:
        doc = json.load(f)
    meta = doc['meta']
    metrics = meta['metrics']
    if set(metrics) != set(METRIC_SPECS):
        raise ValueError(f'{path}: metrics {metrics} != {list(METRIC_SPECS)}')
    for name, (_, units) in METRIC_SPECS.items():
        if meta['metric_units'][name] != units:
            raise ValueError(f'{path}: {name!r} is in '
                             f'{meta["metric_units"][name]!r}, the figure '
                             f'labels assume {units!r}')
    unknown = [g for g in meta['groups'] if g not in GROUP_STYLES]
    if unknown:
        raise ValueError(f'{path}: unit groups without a style {unknown}')
    missing = [k for row in LAYOUT for k in row if k not in doc['scenarios']]
    if missing:
        raise ValueError(f'{path}: no scenario record for {missing}')
    return doc


def breakdown_shares(breakdown, groups, metric):
    """(shares [%] per group in `groups` order, positive total, net total)
    of one metric: each group's value over the sum of the POSITIVE values
    (all zeros when nothing is positive), and the plain sum of all values."""
    values = np.array([breakdown[g][metric] for g in groups], dtype=float)
    positive = float(values[values > 0].sum())
    shares = 100*values/positive if positive > 0 else np.zeros_like(values)
    return shares, positive, float(values.sum())


def display_total(net, metric, operating_hours):
    """A metric's net total in its displayed unit (operating cost USD/hr ->
    MM$/y; the others unchanged)."""
    return net*operating_hours/1e6 if metric == OPERATING else net


def format_total(x, sig=TOTAL_SIG_FIGS):
    """`x` to `sig` significant figures in plain notation (never 1.37e+03),
    with a typographic minus."""
    if not math.isfinite(x):
        return '—'
    if x == 0:
        return '0'
    digits = sig - 1 - math.floor(math.log10(abs(x)))
    r = round(x, digits)
    if r != 0:   # rounding can carry into the next decade (9.996 -> 10.0)
        digits = sig - 1 - math.floor(math.log10(abs(r)))
        r = round(x, digits)
    return f'{r:.{max(digits, 0)}f}'.replace('-', '−')


def panel_subtitle(record):
    """'#1912 · IRR 27.3 %' ('Scenario A' for the baseline; 'IRR —' when the
    IRR is not finite, e.g. -inf for an outright money-loser)."""
    irr = record['IRR']
    irr_text = (f'IRR {100*irr:.1f} %'.replace('-', '−')
                if irr is not None and math.isfinite(irr) else 'IRR —')
    head = ('Scenario A' if record['trial_number'] is None
            else f'#{record["trial_number"]}')
    return f'{head} · {irr_text}'


def y_bottom(doc):
    """The shared y-axis floor: the most negative share of any bar, rounded
    down to 10 % (0 when nothing is negative)."""
    groups, low = doc['meta']['groups'], 0.0
    for key in (k for row in LAYOUT for k in row):
        for metric in doc['meta']['metrics']:
            shares, _, _ = breakdown_shares(doc['scenarios'][key]['breakdown'],
                                            groups, metric)
            low = min(low, float(shares[shares < 0].sum()))
    return 10.0*math.floor(low/10.0)


def legend_groups(doc):
    """Groups visible somewhere in the figure (|share| >= LEGEND_MIN_SHARE in
    at least one bar), in stack order; and the omitted ones."""
    groups, shown = doc['meta']['groups'], set()
    for key in (k for row in LAYOUT for k in row):
        for metric in doc['meta']['metrics']:
            shares, _, _ = breakdown_shares(doc['scenarios'][key]['breakdown'],
                                            groups, metric)
            shown.update(g for g, s in zip(groups, shares)
                         if abs(s) >= LEGEND_MIN_SHARE)
    return ([g for g in groups if g in shown],
            [g for g in groups if g not in shown])


def legend_columns(shown, rows=LEGEND_ROWS):
    """Column-major legend slots: each GROUP_FAMILIES family starts a new
    column, split into columns of `rows` and padded with None (a blank
    entry) so every column has `rows` slots."""
    columns = []
    for family in GROUP_FAMILIES:
        members = [g for g in family if g in shown]
        for i in range(0, len(members), rows):
            chunk = members[i:i + rows]
            columns.append(chunk + [None]*(rows - len(chunk)))
    stray = [g for g in shown if not any(g in f for f in GROUP_FAMILIES)]
    if stray:
        raise ValueError(f'unit groups in no legend family: {stray}')
    return columns


def legend_label(group):
    """Sentence-case display name ('feedstock acquisition' -> 'Feedstock
    acquisition'; 'DDGS recovery' unchanged)."""
    return group[:1].upper() + group[1:]


#%% Figure
def apply_font_rcparams():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = [FONT_FAMILY, 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = FONT_FAMILY
    plt.rcParams['mathtext.it'] = f'{FONT_FAMILY}:italic'
    plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'
    plt.rcParams['hatch.color'] = EDGE_COLOR
    plt.rcParams['hatch.linewidth'] = 0.5
    plt.rcParams['pdf.fonttype'] = 42


def style_ticks(ax):
    """Ticks on all four sides: top/right inward only, left/bottom in and out
    the same length each way (call after a canvas draw)."""
    for which, L in TICK_LEN.items():
        ax.tick_params(axis='both', which=which, direction='inout',
                       length=2*L, top=True, right=True)
        get = 'get_major_ticks' if which == 'major' else 'get_minor_ticks'
        for tick in getattr(ax.xaxis, get)():
            tick.tick2line.set_marker(TICKDOWN)
            tick.tick2line.set_markersize(L)
        for tick in getattr(ax.yaxis, get)():
            tick.tick2line.set_marker(TICKLEFT)
            tick.tick2line.set_markersize(L)


def draw_panel(ax, record, groups, metrics, operating_hours):
    """Five stacked bars of one scenario + their totals above the stacks."""
    for j, metric in enumerate(metrics):
        shares, _, net = breakdown_shares(record['breakdown'], groups, metric)
        up = down = 0.0
        for group, share in zip(groups, shares):
            if abs(share) < MIN_DRAWN_SHARE:
                continue
            color, hatch = GROUP_STYLES[group]
            bottom = up if share > 0 else down
            ax.bar(j, share, bottom=bottom, width=BAR_WIDTH, color=color,
                   hatch=hatch, edgecolor=EDGE_COLOR, linewidth=EDGE_WIDTH,
                   zorder=2)
            if share > 0:
                up += share
            else:
                down += share
        total = display_total(net, metric, operating_hours)
        ax.text(j, 100.0 + TOTAL_PAD, format_total(total), ha='center',
                va='bottom', fontsize=FONTS['total'])
    ax.axhline(0.0, color='k', linewidth=0.8, zorder=3)
    ax.set_title(record['label'], fontsize=FONTS['title'], fontweight='bold',
                 pad=18)
    ax.annotate(panel_subtitle(record), xy=(0.5, 1.0), xycoords='axes fraction',
                xytext=(0, 4), textcoords='offset points', ha='center',
                va='bottom', fontsize=FONTS['subtitle'])


def make_figure(doc):
    apply_font_rcparams()
    meta = doc['meta']
    groups, metrics = meta['groups'], meta['metrics']
    hours = meta['operating_hours']
    bottom = y_bottom(doc)
    fig = plt.figure(figsize=(11.5, 12.8), layout='constrained')
    axs = fig.subplots(len(LAYOUT), len(LAYOUT[0]), sharex=True, sharey=True)
    for row, keys in zip(axs, LAYOUT):
        for ax, key in zip(row, keys):
            draw_panel(ax, doc['scenarios'][key], groups, metrics, hours)
    ax0 = axs[0, 0]
    ax0.set_xlim(-0.6, len(metrics) - 0.4)
    ax0.set_ylim(bottom, Y_TOP)
    ax0.xaxis.set_major_locator(FixedLocator(range(len(metrics))))
    ax0.xaxis.set_minor_locator(NullLocator())
    ax0.yaxis.set_major_locator(FixedLocator(
        np.arange(Y_MAJOR*math.ceil(bottom/Y_MAJOR), 100.0 + 1e-9, Y_MAJOR)))
    ax0.yaxis.set_minor_locator(MultipleLocator(Y_MINOR))
    for ax in axs[-1]:
        ax.set_xticks(range(len(metrics)),
                      [METRIC_SPECS[m][0] for m in metrics],
                      fontsize=FONTS['category'])
    for ax in axs.flat:
        ax.label_outer()
        ax.tick_params(labelsize=FONTS['tick'])
        ax.tick_params(axis='x', labelsize=FONTS['category'])
    fig.supylabel('Cost and utility breakdown [%]', fontsize=FONTS['axis'])
    shown, _ = legend_groups(doc)
    columns = legend_columns(shown)
    handles = [Patch(facecolor=GROUP_STYLES[g][0], hatch=GROUP_STYLES[g][1],
                     edgecolor=EDGE_COLOR, linewidth=EDGE_WIDTH,
                     label=legend_label(g))
               if g is not None else
               Patch(facecolor='none', edgecolor='none', label=' ')
               for column in columns for g in column]
    fig.legend(handles=handles, loc='outside lower center', ncol=len(columns),
               fontsize=FONTS['legend'], frameon=False, handlelength=2.2,
               handleheight=1.3, columnspacing=1.6)
    fig.canvas.draw()
    for ax in axs.flat:
        style_ticks(ax)
        ax.tick_params(axis='x', which='minor', length=0)
    return fig


def print_summary(doc):
    meta = doc['meta']
    groups, metrics, hours = meta['groups'], meta['metrics'], meta['operating_hours']
    head = ''.join(f'{m.split()[0][:11]:>13}' for m in metrics)
    print(f'{"scenario":<34}{head}')
    for key in (k for row in LAYOUT for k in row):
        r = doc['scenarios'][key]
        cells = []
        for m in metrics:
            _, _, net = breakdown_shares(r['breakdown'], groups, m)
            cells.append(f'{format_total(display_total(net, m, hours)):>13}')
        print(f'{r["label"] + " (" + panel_subtitle(r) + ")":<34}'[:34]
              + ''.join(cells))
    _, omitted = legend_groups(doc)
    if omitted:
        print(f'not in the legend (|share| < {LEGEND_MIN_SHARE} % in every '
              f'bar): {omitted}')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--data', default=None,
                    help=f'stage-1 JSON (default: newest {DATA_GLOB})')
    ap.add_argument('--out-dir', default=RESULTS_DIR)
    ap.add_argument('--stem', default=DEFAULT_STEM)
    ap.add_argument('--dpi', type=int, default=300)
    args = ap.parse_args(argv)
    path = args.data or newest_data()
    print(f'data: {path}')
    doc = load_breakdowns(path)
    print_summary(doc)
    fig = make_figure(doc)
    stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
    base = os.path.join(args.out_dir, f'{args.stem}_{stamp}')
    fig.savefig(base + '.png', dpi=args.dpi)
    fig.savefig(base + '.pdf')
    plt.close(fig)
    print(f'wrote {base}.png / .pdf')
    return base


if __name__ == '__main__':
    main()
