#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Offline logic + render test of plots/plot_tea_breakdowns_split12d.py
(loaded BY FILE PATH; no load(), no simulation, never imports the package),
on a synthetic stage-1 document. Checks: (1) positive shares sum to 100 and
credits stay negative; (2) the net total is the plain sum and operating cost
converts USD/hr -> MM$/y; (3) an all-zero metric gives zero shares, not NaN;
(4) format_total's plain 3-significant-figure notation; (5) panel subtitles;
(6) the shared y floor (revenue stacks + their label room) and the legend
omission of an invisible category; (7) legend columns keep each family in
its own padded columns; (8) load_breakdowns refuses wrong units / an
unstyled group / a missing panel / a document without (or with other)
product revenue / half of a merged group set; (9) main() renders PNG + PDF
from a stage-1 (unmerged) document; (10) the per-scenario breakdown CSVs
(names, shape, displayed units, revenue rows, net-total and revenue-total
rows, shares); (11) the operating-cost bar's revenue credits (order, shares
of the positive cost total, AOC and the other bars unchanged); (12) the two
natural-gas groups merged into one category (per-metric sums, stack place,
input untouched, idempotent, closure). Exit 0 + ALL 12 CHECKS PASSED =
clean."""
import os
import csv
import sys
import json
import math
import tempfile
import importlib.util

import numpy as np

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location(
    'ptb', os.path.join(PKG_DIR, 'plots', 'plot_tea_breakdowns_split12d.py'))
ptb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ptb)

HOURS = 8000.0
#: drawn categories (natural gas merged) and the stage-1 unit groups
GROUPS = list(ptb.GROUP_STYLES)
STAGE1_GROUPS = [m for g in GROUPS for m in ptb.MERGED_GROUPS.get(g, (g,))]
METRICS = list(ptb.METRIC_SPECS)
UNITS = {m: u for m, (_, u) in ptb.METRIC_SPECS.items()}
PRODUCTS = list(ptb.REVENUE_STYLES)
OPERATING = 'Operating cost'
#: USD/hr; about 1.3x the synthetic positive operating cost (19 x ~10)
REVENUE = {'ethanol': 150.0, 'isobutanol': 40.0, 'DDGS': 60.0,
           'corn oil': 8.0}


def synthetic_doc(seed=0):
    """A stage-1 document (natural gas NOT merged)."""
    rng = np.random.default_rng(seed)
    scenarios = {}
    keys = [k for row in ptb.LAYOUT for k in row]
    for n, key in enumerate(keys):
        breakdown = {}
        for g in STAGE1_GROUPS:
            breakdown[g] = {m: float(rng.uniform(0.5, 20.0)) for m in METRICS}
        # credits: HXN savings on the duties; excess electricity ~0 everywhere
        breakdown['heat exchanger network']['Heating duty'] = -40.0
        breakdown['heat exchanger network']['Cooling duty'] = -10.0
        for m in METRICS:
            breakdown['excess electricity'][m] = -1e-9
        # the baseline sells no isobutanol; the flagship sells the most
        revenue = {p: (0.0 if key == 'baseline' and p == 'isobutanol'
                       else v*(1.3 if key == 'flagship' else 1.0))
                   for p, v in REVENUE.items()}
        scenarios[key] = dict(
            key=key, label=key.replace('_', ' '),
            trial_number=None if key == 'baseline' else 100 + n,
            IRR=-math.inf if key == 'ibo_titer' else 0.1 + 0.01*n,
            breakdown=breakdown, revenue=revenue)
    # copies: check 8 mutates its documents
    meta = dict(metrics=list(METRICS), metric_units=dict(UNITS),
                groups=list(STAGE1_GROUPS), operating_hours=HOURS,
                revenue_products=list(PRODUCTS))
    return dict(meta=meta, order=keys, scenarios=scenarios)


n_pass, failures = 0, []
def CHECK(label, fn):
    global n_pass
    try:
        detail = fn()
    except Exception as e:
        failures.append(label)
        print(f'FAIL: {label}\n      {type(e).__name__}: {e}', flush=True)
    else:
        n_pass += 1
        print(f'PASS {n_pass}: {label}' + (f' ({detail})' if detail else ''),
              flush=True)


raw_doc = synthetic_doc()
doc = ptb.merge_groups(raw_doc)      # what load_breakdowns hands the figure
bd = doc['scenarios']['flagship']['breakdown']


def check_1():
    shares, positive, _ = ptb.breakdown_shares(bd, GROUPS, 'Heating duty')
    assert abs(shares[shares > 0].sum() - 100.0) < 1e-9, shares[shares > 0].sum()
    i = GROUPS.index('heat exchanger network')
    assert shares[i] < 0 and abs(shares[i] + 100*40.0/positive) < 1e-9, shares[i]
    return f'HXN heating share {shares[i]:.1f} %'
CHECK('positive shares sum to 100; credits negative', check_1)


def check_2():
    _, _, net = ptb.breakdown_shares(bd, GROUPS, 'Operating cost')
    ref = sum(bd[g]['Operating cost'] for g in GROUPS)
    assert abs(net - ref) < 1e-9, (net, ref)
    shown = ptb.display_total(net, 'Operating cost', HOURS)
    assert abs(shown - ref*HOURS/1e6) < 1e-12, shown
    assert ptb.display_total(7.0, 'Cooling duty', HOURS) == 7.0
    return f'AOC {shown:.4f} MM$/y'
CHECK('net total = plain sum; operating cost USD/hr -> MM$/y', check_2)


def check_3():
    zero = {g: {'Cooling duty': 0.0} for g in GROUPS}
    zero['heat exchanger network']['Cooling duty'] = -3.0
    shares, positive, net = ptb.breakdown_shares(zero, GROUPS, 'Cooling duty')
    assert positive == 0.0 and net == -3.0, (positive, net)
    assert np.all(shares == 0.0) and np.all(np.isfinite(shares)), shares
CHECK('a metric with nothing positive gives zero shares, not NaN', check_3)


def check_4():
    cases = {1370.4: '1370', 137.04: '137', 13.704: '13.7', 9.404: '9.40',
             9.996: '10.0', 999.6: '1000', 0.012345: '0.0123',
             -2.5: '−2.50', 0.0: '0', math.nan: '—'}
    for x, want in cases.items():
        got = ptb.format_total(x)
        assert got == want, f'format_total({x!r}) = {got!r}, want {want!r}'
    return f'{len(cases)} cases'
CHECK('format_total: plain 3-sig-fig notation, typographic minus', check_4)


def check_5():
    s = doc['scenarios']
    assert ptb.panel_subtitle(s['baseline']).startswith('Scenario A · IRR ')
    assert ptb.panel_subtitle(s['ibo_titer']) == '#104 · IRR —', \
        ptb.panel_subtitle(s['ibo_titer'])
    got = ptb.panel_subtitle(dict(trial_number=7, IRR=-0.052))
    assert got == '#7 · IRR −5.2 %', got
CHECK('panel subtitles (baseline, -inf IRR, negative IRR)', check_5)


def check_6():
    # the deepest stack is an operating-cost bar's revenue credits
    lows = {}
    for row in ptb.LAYOUT:
        for k in row:
            r = doc['scenarios'][k]
            for m in METRICS:
                shares, positive, _ = ptb.breakdown_shares(r['breakdown'],
                                                           GROUPS, m)
                low = float(shares.clip(max=0).sum())
                if m == OPERATING:
                    low -= 100*sum(r['revenue'].values())/positive
                    low -= ptb.REVENUE_LABEL_ROOM
                lows[k, m] = low
    low = min(lows.values())
    assert min(lows, key=lows.get)[1] == OPERATING, lows
    bottom = ptb.y_bottom(doc)
    assert bottom <= low and bottom > low - 10 and bottom % 10 == 0, (bottom, low)
    shown, omitted = ptb.legend_groups(doc)
    assert omitted == ['excess electricity'], omitted
    assert shown == ([g for g in GROUPS if g != 'excess electricity']
                     + list(ptb.REVENUE_CATEGORIES)), shown
    return f'y floor {bottom:g} % (deepest stack + label room {low:.1f} %)'
CHECK('shared y floor incl. revenue; invisible group left out of the legend',
      check_6)


def check_7():
    shown, _ = ptb.legend_groups(doc)
    cols = ptb.legend_columns(shown)
    assert all(len(c) == ptb.LEGEND_ROWS for c in cols), [len(c) for c in cols]
    fams = ptb.GROUP_FAMILIES
    assert cols[0] + cols[1] == list(fams[0]), 'process areas not in columns 1-2'
    assert cols[2] == list(fams[1]), 'facilities not in column 3'
    assert [g for g in cols[3] if g] == [g for g in fams[2] if g in shown]
    assert cols[3][-1] is None, 'pseudo-group column not padded'
    assert cols[4] == list(ptb.REVENUE_CATEGORIES) + [None], cols[4]
    assert len(cols) == 5, len(cols)
    return f'{len(cols)} columns'
CHECK('legend columns: one family per column set, padded', check_7)


def check_8():
    with tempfile.TemporaryDirectory() as tmp:
        def refused(mutate, needle=''):
            bad = synthetic_doc()
            mutate(bad)
            path = os.path.join(tmp, 'bad.json')
            with open(path, 'w') as f:
                json.dump(bad, f)
            try:
                ptb.load_breakdowns(path)
            except ValueError as e:
                return needle in str(e)
            return False
        assert refused(lambda d: d['meta']['metric_units'].update(
            {'Operating cost': 'MM$/yr'})), 'wrong units accepted'
        assert refused(lambda d: d['meta']['groups'].append('mystery group')), \
            'unstyled group accepted'
        assert refused(lambda d: d['scenarios'].pop('flagship')), \
            'missing panel accepted'
        # a pre-revenue stage-1 document asks for a re-run
        assert refused(lambda d: d['meta'].pop('revenue_products'),
                       're-run'), 'document without revenue_products accepted'
        assert refused(lambda d: d['scenarios']['ibo_yield'].pop('revenue'),
                       're-run'), 'record without revenue accepted'
        assert refused(lambda d: d['scenarios']['etoh_titer']['revenue'].pop(
            'corn oil')), 'record missing a product accepted'
        assert refused(lambda d: d['meta'].update(
            revenue_products=PRODUCTS[::-1])), 'other product order accepted'
        # half of a merged set cannot be merged
        assert refused(lambda d: d['meta']['groups'].remove(
            'natural gas (for product drying)'), 'cannot merge'), \
            'partial natural-gas set accepted'
CHECK('load_breakdowns refuses wrong units / unstyled group / missing panel / '
      'missing or other revenue / a partial merge set', check_8)


def check_9():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'tea_breakdowns_split12d_test.json')
        with open(path, 'w') as f:
            json.dump(raw_doc, f)            # stage-1 form: main() merges
        base = ptb.main(['--data', path, '--out-dir', tmp, '--dpi', '60'])
        for ext in ('.png', '.pdf'):
            assert os.path.getsize(base + ext) > 0, base + ext
CHECK('main() renders PNG + PDF from a synthetic document', check_9)


def check_10():
    with tempfile.TemporaryDirectory() as tmp:
        data = os.path.join(tmp, 'tea_breakdowns_split12d_X.json')
        paths = ptb.write_breakdown_csvs(doc, data, tmp)
        names = sorted(os.path.basename(p) for p in paths)
        assert len(paths) == 9, names
        assert 'tea_breakdowns_split12d_X_baseline.csv' in names, names
        assert 'tea_breakdowns_split12d_X_flagship_trial102.csv' in names, names
        with open(os.path.join(tmp, 'tea_breakdowns_split12d_X_flagship_trial102.csv'),
                  newline='') as f:
            rows = list(csv.reader(f))
        header, body = rows[0], rows[1:]
        assert header[0] == 'Unit group' and len(header) == 1 + 2*len(METRICS)
        assert 'Operating cost [MM$/yr]' in header, header
        names = [f'{p} revenue' for p in PRODUCTS]
        assert [r[0] for r in body] == (GROUPS + names
                                        + ['Total (net)', 'Total revenue']), \
            [r[0] for r in body]
        n = len(GROUPS)
        groups_rows, revenue_rows = body[:n], body[n:n + len(names)]
        net_row, revenue_total = body[-2], body[-1]
        j = header.index('Operating cost [MM$/yr]')
        i = GROUPS.index('boiler')
        want = bd['boiler']['Operating cost']*HOURS/1e6
        assert abs(float(body[i][j]) - want) < 1e-12, (body[i][j], want)
        for m in METRICS:
            col = header.index(f'{m} [{ptb.DISPLAY_UNITS[m]}]')
            scol = header.index(f'{m} share [% of positive total]')
            net = sum(float(r[col]) for r in groups_rows)
            assert abs(float(net_row[col]) - net) <= 1e-9*abs(net), (m, net)
            s = [float(r[scol]) for r in groups_rows]
            assert abs(sum(x for x in s if x > 0) - 100.0) < 1e-9, (m, sum(s))
            if m == OPERATING:
                continue
            # the revenue rows carry operating cost only
            assert all(r[col] == '' and r[scol] == '' for r in revenue_rows)
            assert revenue_total[col] == ''
        flagship = doc['scenarios']['flagship']
        _, positive, _ = ptb.breakdown_shares(bd, GROUPS, OPERATING)
        for p, r in zip(PRODUCTS, revenue_rows):
            credit = -flagship['revenue'][p]*HOURS/1e6
            assert abs(float(r[j]) - credit) < 1e-12, (p, r[j], credit)
            share = float(r[header.index(
                f'{OPERATING} share [% of positive total]')])
            assert abs(share + 100*flagship['revenue'][p]/positive) < 1e-9
        total = sum(float(r[j]) for r in revenue_rows)
        assert total < 0 and abs(float(revenue_total[j]) - total) < 1e-12
        assert net_row[1 + len(METRICS):] == ['']*len(METRICS)
        assert revenue_total[1 + len(METRICS):] == ['']*len(METRICS)
    return f'{len(paths)} CSVs x {len(body)} rows'
CHECK('per-scenario breakdown CSVs: names, shape, units, revenue rows, totals, '
      'shares', check_10)


def check_11():
    flagship = doc['scenarios']['flagship']
    segments, net, revenue = ptb.bar_segments(flagship, GROUPS, OPERATING)
    # the unit groups first, then the credits in REVENUE_STYLES order
    assert [c for c, _ in segments] == GROUPS + list(ptb.REVENUE_CATEGORIES)
    shares, positive, net_ref = ptb.breakdown_shares(bd, GROUPS, OPERATING)
    assert net == net_ref, (net, net_ref)       # the printed AOC is unchanged
    assert abs(revenue - sum(flagship['revenue'].values())) < 1e-12
    got = dict(segments[len(GROUPS):])
    for p in PRODUCTS:
        want = -100*flagship['revenue'][p]/positive
        assert abs(got[f'{p} revenue'] - want) < 1e-12, (p, got, want)
    depth = sum(got.values())
    assert depth < -100, depth    # sales (1.3x positive cost) exceed the AOC
    # every other bar is exactly the unit groups
    for m in METRICS:
        if m == OPERATING:
            continue
        seg, _, rev = ptb.bar_segments(flagship, GROUPS, m)
        assert rev == 0.0 and [c for c, _ in seg] == GROUPS, m
    assert ptb.category_style('corn oil revenue') == (
        *ptb.REVENUE_STYLES['corn oil'], ptb.REVENUE_HATCH_COLOR)
    assert ptb.category_style('boiler')[2] == ptb.EDGE_COLOR
    return f'flagship revenue stack {depth:.1f} %'
CHECK('operating-cost bar: revenue credits after the unit groups, shares of '
      'the positive cost, AOC unchanged', check_11)


def check_12():
    members = ptb.MERGED_GROUPS['natural gas']
    groups = doc['meta']['groups']
    assert groups == GROUPS, groups              # the drawn order
    assert groups.index('natural gas') == STAGE1_GROUPS.index(members[0])
    for key, record in doc['scenarios'].items():
        raw = raw_doc['scenarios'][key]['breakdown']
        assert not any(m in record['breakdown'] for m in members), key
        for m in METRICS:
            want = sum(raw[g][m] for g in members)
            assert record['breakdown']['natural gas'][m] == want, (key, m)
    # the input is untouched; merging a merged document is a no-op
    assert all(m in raw_doc['meta']['groups'] for m in members)
    assert all(m in raw_doc['scenarios']['baseline']['breakdown']
               for m in members)
    assert ptb.merge_groups(doc) == doc
    # the merged column still closes: net = the stage-1 unit-group sum
    _, _, net = ptb.breakdown_shares(bd, GROUPS, OPERATING)
    raw_bd = raw_doc['scenarios']['flagship']['breakdown']
    ref = sum(raw_bd[g][OPERATING] for g in STAGE1_GROUPS)
    assert abs(net - ref) < 1e-9, (net, ref)
    assert ptb.GROUP_STYLES['natural gas'] == ('#ffffff', '\\\\\\\\')
    return f'{len(STAGE1_GROUPS)} stage-1 groups -> {len(GROUPS)} categories'
CHECK('natural gas merged: summed per metric at the first member\'s place, '
      'input untouched, idempotent', check_12)


if failures:
    print(f'\n{len(failures)} CHECK(S) FAILED: {failures}')
    sys.exit(1)
print(f'ALL {n_pass} CHECKS PASSED')
