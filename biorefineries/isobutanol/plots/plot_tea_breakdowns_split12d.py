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
sold excess electricity) sit below zero; stage 1's two natural-gas groups
(steam generation, product drying) are drawn and tabulated as one 'natural
gas' category (MERGED_GROUPS). The operating-cost bar also carries
the sales revenue of the four products (ethanol, isobutanol, DDGS, corn oil;
spec docs/superpowers/specs/2026-09-24-tea-breakdowns-revenue-design.md) as
credits under the unit groups, on the same scale: a revenue stack deeper than
-100 % means sales exceed the operating cost. Each bar's NET unit-group total
-- operating cost converted to MM$/y, i.e. the AOC -- is printed above it,
and the total revenue [MM$/y] under the operating-cost bar's revenue stack.

Sim-safe: reads the stage-1 JSON only (json / numpy / matplotlib); never
imports the biorefineries package, never load()s. Run:

    python plots/plot_tea_breakdowns_split12d.py [--data <json>] [--out-dir DIR]

--data defaults to the newest analyses/results/tea_breakdowns_split12d_*.json.
Writes <stem>_<stamp>.png and .pdf to --out-dir (default analyses/results),
and one breakdown CSV per scenario named after the data file,
<data stem>_<key>[_trial<N>].csv (rows = the unit groups in stack order,
the four '<product> revenue' credits (negative, operating-cost column
only), 'Total (net)' (the unit-group sum) and 'Total revenue' (the credits'
sum); columns = each metric in its displayed unit, then each row's share of
the metric's positive total [%]); re-running on the same data overwrites
them.
"""
import os
import csv
import copy
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
#: displayed unit per metric (plain ASCII, for the CSV headers)
DISPLAY_UNITS = {'Installed equipment cost': 'MM$', 'Cooling duty': 'GJ/hr',
                 'Heating duty': 'GJ/hr', 'Electricity consumption': 'MW',
                 OPERATING: 'MM$/yr'}
METRIC_SPECS = {
    'Installed equipment cost': ('Installed\nequipment\ncost\n[MM\\$]', 'MM$'),
    'Cooling duty': ('Cooling\nduty\n[$\\mathrm{GJ·h}^{-1}$]', 'GJ/hr'),
    'Heating duty': ('Heating\nduty\n[$\\mathrm{GJ·h}^{-1}$]', 'GJ/hr'),
    'Electricity consumption': ('Electricity\nconsumption\n[MW]', 'MW'),
    # 'and revenue' on line 3, beside the short '[MW]' ('cost and' on line 2
    # ran into 'consumption')
    OPERATING: ('Operating\ncost\nand revenue\n[$\\mathrm{MM\\$·y}^{-1}$]',
                'USD/hr'),
}

#: unit group -> (face colour, hatch). Process areas take the validated
#: eight-slot categorical palette in order (adjacent stack pairs clear the CVD
#: and normal-vision floors); the two past eight repeat slots 3 and 1 with a
#: texture (never a generated ninth hue). Facilities are neutral greys, the
#: turbogenerator sharing the boiler's grey (both are BT801) under a texture.
#: The cost pseudo-groups (operating cost only) are white + texture; stage
#: 1's two natural-gas groups are drawn as ONE (MERGED_GROUPS), under the
#: texture the steam-generation share (the bulk of it) had.
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
    'natural gas':                        ('#ffffff', '\\\\\\\\'),
    'fixed operating cost':               ('#ffffff', 'xxxx'),
    'excess electricity':                 ('#ffffff', '....'),
}
#: drawn category -> the stage-1 unit groups summed into it (per metric, at
#: the first member's place in the stack)
MERGED_GROUPS = {
    'natural gas': ('natural gas (for steam generation)',
                    'natural gas (for product drying)'),
}
#: product (stage-1 `revenue` key, USD/hr) -> (face colour, hatch), in the
#: order the revenue credits stack down the operating-cost bar. Each product
#: takes its own process area's hue (ethanol / isobutanol purification, DDGS
#: recovery; corn oil the amber slot) under a small-circle texture that no
#: unit group uses, drawn in REVENUE_HATCH_COLOR so it shows on the dark
#: hues, so a product never reads as its unit group. Adjacent pairs pass the
#: dataviz validator (worst CVD dE 9.1, DDGS / corn oil).
REVENUE_STYLES = {
    'ethanol':    ('#008300', 'oo'),
    'isobutanol': ('#4a3aa7', 'oo'),
    'DDGS':       ('#1baf7a', 'oo'),
    'corn oil':   ('#eda100', 'oo'),
}
REVENUE_HATCH_COLOR = 'white'
#: legend / CSV category name of a product's revenue credit
REVENUE_CATEGORIES = {f'{p} revenue': p for p in REVENUE_STYLES}
#: legend columns: each family starts a new column (process areas fill two),
#: columns padded to LEGEND_ROWS with blank entries so families line up
GROUP_FAMILIES = (
    tuple(list(GROUP_STYLES)[:10]),     # process areas
    tuple(list(GROUP_STYLES)[10:15]),   # facilities
    tuple(list(GROUP_STYLES)[15:]),     # cost pseudo-groups
    tuple(REVENUE_CATEGORIES),          # product revenue (credits)
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
TOTAL_PAD = 2.0        # gap between a stack end and its total [%]
#: room kept under the deepest revenue stack for its total [%]
REVENUE_LABEL_ROOM = 12.0
Y_MAJOR, Y_MINOR = 25.0, 5.0
FIGSIZE = (11.5, 15.5)

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


def merge_groups(doc):
    """A copy of a stage-1 document with each MERGED_GROUPS set of unit
    groups summed into one category (per metric, in every scenario) at its
    first member's place in meta['groups']. An already-merged document
    passes through; a partial set raises ValueError."""
    doc = copy.deepcopy(doc)
    for name, members in MERGED_GROUPS.items():
        groups = doc['meta']['groups']
        present = [g for g in members if g in groups]
        if not present and name in groups:
            continue
        if len(present) != len(members):
            raise ValueError(f'cannot merge {name!r}: unit groups '
                             f'{[g for g in members if g not in groups]} '
                             'missing')
        position = min(groups.index(g) for g in members)
        merged = [g for g in groups if g not in members]
        merged.insert(position, name)
        doc['meta']['groups'] = merged
        for record in doc['scenarios'].values():
            parts = [record['breakdown'].pop(g) for g in members]
            record['breakdown'][name] = {m: sum(p[m] for p in parts)
                                         for m in parts[0]}
    return doc


def load_breakdowns(path):
    """The stage-1 document with MERGED_GROUPS merged, validated against
    what this figure assumes."""
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
    try:
        doc = merge_groups(doc)
    except ValueError as e:
        raise ValueError(f'{path}: {e}') from None
    meta = doc['meta']
    unknown = [g for g in meta['groups'] if g not in GROUP_STYLES]
    if unknown:
        raise ValueError(f'{path}: unit groups without a style {unknown}')
    missing = [k for row in LAYOUT for k in row if k not in doc['scenarios']]
    if missing:
        raise ValueError(f'{path}: no scenario record for {missing}')
    no_revenue = [k for row in LAYOUT for k in row
                  if 'revenue' not in doc['scenarios'][k]]
    if 'revenue_products' not in meta or no_revenue:
        raise ValueError(f'{path} has no product revenue '
                         f'({no_revenue or "meta"}): re-run '
                         'analyses/collect_tea_breakdowns_split12d.py')
    products = list(REVENUE_STYLES)
    wrong = [k for row in LAYOUT for k in row
             if sorted(doc['scenarios'][k]['revenue']) != sorted(products)]
    if meta['revenue_products'] != products or wrong:
        raise ValueError(f'{path}: revenue products '
                         f'{meta["revenue_products"]} (records {wrong}) != '
                         f'{products}')
    return doc


def breakdown_shares(breakdown, groups, metric):
    """(shares [%] per group in `groups` order, positive total, net total)
    of one metric: each group's value over the sum of the POSITIVE values
    (all zeros when nothing is positive), and the plain sum of all values."""
    values = np.array([breakdown[g][metric] for g in groups], dtype=float)
    positive = float(values[values > 0].sum())
    shares = 100*values/positive if positive > 0 else np.zeros_like(values)
    return shares, positive, float(values.sum())


def revenue_shares(record, positive):
    """{product: share [%]} of the product revenues drawn as credits in the
    operating-cost bar: -revenue over the bar's POSITIVE cost total, the
    rule of every segment (zeros when nothing is positive)."""
    return {p: (-100*record['revenue'][p]/positive if positive > 0 else 0.0)
            for p in REVENUE_STYLES}


def bar_segments(record, groups, metric):
    """One bar in stack order: ([(category, share [%]), ...], the unit
    groups' net total, the total revenue [USD/hr]). The categories are the
    unit groups, then -- operating cost only -- the '<product> revenue'
    credits; the revenue is 0 for every other metric."""
    shares, positive, net = breakdown_shares(record['breakdown'], groups,
                                             metric)
    segments = list(zip(groups, (float(s) for s in shares)))
    if metric != OPERATING:
        return segments, net, 0.0
    segments += [(f'{p} revenue', s)
                 for p, s in revenue_shares(record, positive).items()]
    return segments, net, float(sum(record['revenue'].values()))


def category_style(category):
    """(face colour, hatch, hatch colour) of a unit group or revenue
    category."""
    if category in GROUP_STYLES:
        return (*GROUP_STYLES[category], EDGE_COLOR)
    return (*REVENUE_STYLES[REVENUE_CATEGORIES[category]], REVENUE_HATCH_COLOR)


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
    """The shared y-axis floor: the deepest negative stack of any bar -- less
    REVENUE_LABEL_ROOM under a revenue stack, for its total -- rounded down
    to 10 % (0 when nothing is negative)."""
    groups, low = doc['meta']['groups'], 0.0
    for key in (k for row in LAYOUT for k in row):
        for metric in doc['meta']['metrics']:
            segments, _, revenue = bar_segments(doc['scenarios'][key], groups,
                                                metric)
            down = sum(s for _, s in segments if s < 0)
            low = min(low, down - (REVENUE_LABEL_ROOM if revenue else 0.0))
    return 10.0*math.floor(low/10.0)


def legend_groups(doc):
    """Categories (unit groups, then revenue credits) visible somewhere in
    the figure (|share| >= LEGEND_MIN_SHARE in at least one bar), in stack
    order; and the omitted ones."""
    groups, shown = doc['meta']['groups'], set()
    for key in (k for row in LAYOUT for k in row):
        for metric in doc['meta']['metrics']:
            segments, _, _ = bar_segments(doc['scenarios'][key], groups,
                                          metric)
            shown.update(c for c, s in segments if abs(s) >= LEGEND_MIN_SHARE)
    categories = list(groups) + list(REVENUE_CATEGORIES)
    return ([c for c in categories if c in shown],
            [c for c in categories if c not in shown])


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
        raise ValueError(f'categories in no legend family: {stray}')
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
    """Five stacked bars of one scenario, their unit-group totals above the
    stacks and the total revenue under the operating-cost bar."""
    for j, metric in enumerate(metrics):
        segments, net, revenue = bar_segments(record, groups, metric)
        up = down = 0.0
        for category, share in segments:
            if abs(share) < MIN_DRAWN_SHARE:
                continue
            color, hatch, hatch_color = category_style(category)
            bottom = up if share > 0 else down
            ax.bar(j, share, bottom=bottom, width=BAR_WIDTH, color=color,
                   hatch=hatch, hatchcolor=hatch_color, edgecolor=EDGE_COLOR,
                   linewidth=EDGE_WIDTH, zorder=2)
            if share > 0:
                up += share
            else:
                down += share
        total = display_total(net, metric, operating_hours)
        ax.text(j, 100.0 + TOTAL_PAD, format_total(total), ha='center',
                va='bottom', fontsize=FONTS['total'])
        if revenue:
            ax.text(j, down - TOTAL_PAD,
                    format_total(display_total(revenue, metric,
                                               operating_hours)),
                    ha='center', va='top', fontsize=FONTS['total'])
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
    fig = plt.figure(figsize=FIGSIZE, layout='constrained')
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
    fig.supylabel('Cost, utility and revenue breakdown [%]',
                  fontsize=FONTS['axis'])
    shown, _ = legend_groups(doc)
    columns = legend_columns(shown)
    handles = []
    for column in columns:
        for c in column:
            if c is None:
                handles.append(Patch(facecolor='none', edgecolor='none',
                                     label=' '))
                continue
            color, hatch, hatch_color = category_style(c)
            handles.append(Patch(facecolor=color, hatch=hatch,
                                 hatchcolor=hatch_color, edgecolor=EDGE_COLOR,
                                 linewidth=EDGE_WIDTH, label=legend_label(c)))
    fig.legend(handles=handles, loc='outside lower center', ncol=len(columns),
               fontsize=FONTS['legend'], frameon=False, handlelength=2.2,
               handleheight=1.3, columnspacing=1.6)
    fig.canvas.draw()
    for ax in axs.flat:
        style_ticks(ax)
        ax.tick_params(axis='x', which='minor', length=0)
    return fig


def csv_path(data_path, out_dir, record):
    """<data stem>_<key>[_trial<N>].csv in `out_dir`."""
    stem = os.path.splitext(os.path.basename(data_path))[0]
    trial = record['trial_number']
    suffix = record['key'] + ('' if trial is None else f'_trial{trial}')
    return os.path.join(out_dir, f'{stem}_{suffix}.csv')


def breakdown_rows(doc, key):
    """(header, rows) of one scenario's breakdown table: a row per unit group
    (stack order) with each metric in its displayed unit and its share of
    the metric's positive total [%]; a row per '<product> revenue' credit
    (SIGNED as drawn: -revenue in the operating-cost columns, the other
    metrics blank); 'Total (net)' (the unit-group sum; shares blank); and
    'Total revenue' (the credits' sum, operating cost only)."""
    meta = doc['meta']
    groups, metrics = meta['groups'], meta['metrics']
    hours = meta['operating_hours']
    record = doc['scenarios'][key]
    header = (['Unit group']
              + [f'{m} [{DISPLAY_UNITS[m]}]' for m in metrics]
              + [f'{m} share [% of positive total]' for m in metrics])
    values, shares, nets = {}, {}, {}
    for m in metrics:
        s, positive, net = breakdown_shares(record['breakdown'], groups, m)
        shares[m] = s
        nets[m] = display_total(net, m, hours)
        values[m] = [display_total(record['breakdown'][g][m], m, hours)
                     for g in groups]
        if m == OPERATING:
            revenue_share = revenue_shares(record, positive)
    rows = [[g] + [values[m][i] for m in metrics]
            + [float(shares[m][i]) for m in metrics]
            for i, g in enumerate(groups)]

    def operating_only(value):
        return [value if m == OPERATING else '' for m in metrics]
    for product, share in revenue_share.items():
        credit = display_total(-record['revenue'][product], OPERATING, hours)
        rows.append([f'{product} revenue'] + operating_only(credit)
                    + operating_only(share))
    rows.append(['Total (net)'] + [nets[m] for m in metrics]
                + ['']*len(metrics))
    revenue = display_total(-sum(record['revenue'].values()), OPERATING, hours)
    rows.append(['Total revenue'] + operating_only(revenue)
                + ['']*len(metrics))
    return header, rows


def write_breakdown_csvs(doc, data_path, out_dir):
    """One CSV per scenario (LAYOUT order); returns the paths."""
    paths = []
    for key in (k for row in LAYOUT for k in row):
        path = csv_path(data_path, out_dir, doc['scenarios'][key])
        header, rows = breakdown_rows(doc, key)
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)
        paths.append(path)
    return paths


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
    products = list(REVENUE_STYLES)
    print(f'\n{"revenue [MM$/y], % of it":<34}{"total":>13}'
          + ''.join(f'{p:>12}' for p in products))
    for key in (k for row in LAYOUT for k in row):
        r = doc['scenarios'][key]
        total = sum(r['revenue'].values())
        cells = [f'{format_total(100*r["revenue"][p]/total)} %' if total
                 else '—' for p in products]
        print(f'{r["label"]:<34}'[:34]
              + f'{format_total(display_total(total, OPERATING, hours)):>13}'
              + ''.join(f'{c:>12}' for c in cells))
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
    csvs = write_breakdown_csvs(doc, path, args.out_dir)
    print(f'wrote {len(csvs)} breakdown CSVs: '
          f'{os.path.join(args.out_dir, os.path.basename(csvs[0]))} ...')
    return base


if __name__ == '__main__':
    main()
