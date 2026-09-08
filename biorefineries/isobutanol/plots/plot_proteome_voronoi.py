#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Proteome-allocation Voronoi treemap -- static publication figure.

A small-multiples grid of amCharts5 Voronoi treemaps, one tile per
parameter set (scenario-A baseline + a handful of kinetic-optimization
campaign trials). Each tile partitions the modeled proteome (cap
eb.PROTEIN_CONTENT = 0.49 g protein/gDCW) into: Housekeeping (fixed),
Metabolic Phi_M (subdivided into the five pathway categories of
plot_kin_opt_parameter_sets.py panel C), Translation phi_T (the model's
active, growth-scaled ribosomal sector -- the "reduced" translation
value), and Unallocated flexible slack. Every tile's cells sum to 0.49.

Two stages: this script (Stage 1, IBO_2026, sim-safe) computes the
allocation and writes JSON, then invokes the voronoi-treemaps conda env's
node on plots/voronoi/render_treemap.mjs (Stage 2) to render PNG/SVG/PDF.
Use --no-render to write only the JSON.

Sim-safe: plot_kin_opt_parameter_sets.py (ps) -- and through it ko and eb
-- are loaded by file path (no biosteam import, no load()); campaign CSVs
and workbooks are plain pandas reads. Runnable while a campaign is in
flight. Run:

    python plots/plot_proteome_voronoi.py \
        --set "Best IRR" <campaign> best \
        --set "Best isobutanol titer" <campaign> "best:IBO titer"

With no --set arguments it plots the same five-study default as the
parameter-sets figure (one optimum per objective) against the baseline
-> 6 tiles.
"""
import os
import json
import argparse
import importlib.util
import subprocess
from datetime import datetime

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLOTS_DIR = os.path.join(PKG_DIR, 'plots')
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')  # matches ps
VORONOI_DIR = os.path.join(PLOTS_DIR, 'voronoi')
# the voronoi-treemaps conda env's node; override with --node or $VORONOI_NODE
DEFAULT_VORONOI_NODE = (r'C:\Users\saran\anaconda3\envs\voronoi-treemaps'
                        r'\node.exe')
RENDER_SCRIPT = os.path.join(VORONOI_DIR, 'render_treemap.mjs')

TOL = 1e-6


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# The existing sim-safe plotter -- single source of truth for set
# selection and burden pools. Loading it by file path triggers its own
# by-file-path loads of ko and eb (still sim-safe).
ps = _load('plot_kin_opt_parameter_sets',
           os.path.join(PLOTS_DIR, 'plot_kin_opt_parameter_sets.py'))
eb = ps.eb

# map each of ps.BURDEN_CATEGORIES' names to a stable palette key (piece).
# Asserted against ps.BURDEN_CATEGORIES below so a category rename/drift in
# the reused plotter raises here at import.
CATEGORY_PIECES = {
    'Glycolysis': 'cat_glycolysis',
    'TCA cycle': 'cat_tca',
    'Acetate / acetyl-CoA\nproduction': 'cat_acetate',
    'Ethanol production': 'cat_ethanol',
    'Isobutanol production': 'cat_isobutanol',
}
_ps_cat_names = [name for name, _steps in ps.BURDEN_CATEGORIES]
if set(_ps_cat_names) != set(CATEGORY_PIECES):
    raise AssertionError(
        'BURDEN_CATEGORIES drift: ps names %r != CATEGORY_PIECES keys %r'
        % (sorted(_ps_cat_names), sorted(CATEGORY_PIECES)))

# legend / color order (top-level pieces with the five metabolic cats
# nested where "metabolic" would sit)
PIECE_ORDER = ('housekeeping', 'cat_glycolysis', 'cat_tca', 'cat_acetate',
               'cat_ethanol', 'cat_isobutanol', 'translation', 'slack')


def tile_from_record(rec):
    """One set record -> one JSON tile (spec schema).

    Housekeeping and F_flex come from eb constants / the record; the five
    metabolic children are the ps.BURDEN_CATEGORIES pool sums; Phi_M is
    their total (checked against rec['Phi_M']). The translation cell is the
    translation the cell can actually build after metabolism -- the derated
    allocation min(phi_T,demand, F_flex - Phi_M) = d * phi_T,demand, not the
    full demand rec['phi_T'] -- so slack = PROTEIN_CONTENT - housekeeping -
    Phi_M - phi_T stays >= 0 and a feasible tile's cells sum to exactly
    PROTEIN_CONTENT. A derated-growth trial (burden_factor d < 1) shows a
    smaller translation cell and zero slack; growth is throttled, not
    infeasible. Only a truly infeasible point (Phi_M >= F_flex, d = 0:
    metabolism alone overruns the flexible sector) clamps slack to 0 and
    records a `warning` rather than emitting a negative cell.
    """
    PC = float(eb.PROTEIN_CONTENT)
    housekeeping = PC * float(eb.HOUSEKEEPING_FRACTION)

    children_metabolic = []
    Phi_M = 0.0
    for name, steps in ps.BURDEN_CATEGORIES:
        w = sum(float(rec[f'pool_{st}']) for st in steps)
        children_metabolic.append(
            {'name': name.replace('\n', ' '),
             'piece': CATEGORY_PIECES[name],
             'value': w})
        Phi_M += w
    if abs(Phi_M - float(rec['Phi_M'])) > 1e-4:
        raise AssertionError(
            f'{rec.get("label")!r}: category-sum Phi_M {Phi_M:.6f} != '
            f'record Phi_M {float(rec["Phi_M"]):.6f}')

    # Translation actually built: the burden model derates growth linearly as
    # the flexible sector fills, so a trial with Phi_M + phi_T,demand > F_flex
    # builds only what it can afford, d * phi_T,demand = F_flex - Phi_M, and
    # grows slower. Draw that allocated translation, not the full demand, so a
    # feasible tile's cells sum to exactly PROTEIN_CONTENT with slack >= 0.
    F_flex = float(rec['F_flex'])
    phi_T_demand = float(rec['phi_T'])
    phi_T = max(0.0, min(phi_T_demand, F_flex - Phi_M))
    slack = PC - housekeeping - Phi_M - phi_T
    warning = None
    if slack < -TOL:
        warning = (f'burden-infeasible: metabolic pool Phi_M {Phi_M:.5f} '
                   f'>= F_flex {F_flex:.5f} (d '
                   f'{float(rec["burden_factor"]):.3f}); no room for '
                   f'translation (demand phi_T {phi_T_demand:.5f})')
        slack = 0.0
    else:
        total = housekeeping + Phi_M + phi_T + slack
        if abs(total - PC) > 1e-4:
            raise AssertionError(
                f'{rec.get("label")!r}: cells sum to {total:.6f}, not '
                f'PROTEIN_CONTENT {PC:.6f}')

    tile = {
        'label': rec['label'],
        'is_baseline': bool(rec.get('is_baseline')),
        'campaign': rec.get('campaign'),
        'trial_number': (None if rec.get('trial_number') is None
                         else int(rec['trial_number'])),
        'Phi_M': Phi_M,
        'phi_T': phi_T,
        'burden_factor': float(rec['burden_factor']),
        'children': [
            {'name': 'Housekeeping', 'piece': 'housekeeping',
             'value': housekeeping},
            {'name': 'Metabolic', 'piece': 'metabolic', 'value': Phi_M,
             'children': children_metabolic},
            {'name': 'Translation', 'piece': 'translation', 'value': phi_T},
            {'name': 'Unallocated flexible', 'piece': 'slack',
             'value': slack},
        ],
    }
    if warning:
        tile['warning'] = warning
    return tile


def build_document(sets, band_campaign):
    """All sets -> the render document (meta + tiles)."""
    return {
        'meta': {
            'protein_content': float(eb.PROTEIN_CONTENT),
            'housekeeping': float(eb.PROTEIN_CONTENT
                                  * eb.HOUSEKEEPING_FRACTION),
            'F_flex': float(eb.F_FLEX),
            'phi_T_wt': float(eb.PHI_T_WT),
            'generated': datetime.now().isoformat(timespec='seconds'),
            'band_campaign': band_campaign,
            'piece_order': list(PIECE_ORDER),
        },
        'tiles': [tile_from_record(s) for s in sets],
    }


# --- default studies (mirror the parameter-sets figure) ---------------------
DEFAULT_SPECS = (
    ('Financial attractiveness', ps.DEFAULT_STUDY, 'best'),
    ('Isobutanol titer', ps.IBO_TITER_STUDY, 'best'),
    ('Ethanol titer', ps.ETOH_TITER_STUDY, 'best'),
    ('Isobutanol yield', ps.IBO_YIELD_STUDY, 'best'),
    ('Ethanol yield', ps.ETOH_YIELD_STUDY, 'best'),
)


def _norm_trial(t):
    if isinstance(t, str) and t.startswith('best'):
        return t
    return int(t)


def resolve_specs(args):
    """(specs, include_baseline) from parsed CLI args, mirroring
    plot_kin_opt_parameter_sets.main."""
    if args.sets:
        specs = [(lab, camp, _norm_trial(tr)) for lab, camp, tr in args.sets]
    else:
        specs = list(DEFAULT_SPECS)
    include_baseline = not args.no_baseline
    if len(specs) + (1 if include_baseline else 0) > ps.MAX_SETS:
        raise ValueError(
            f'at most {ps.MAX_SETS} tiles (baseline + {len(ps.HUE_COLORS)} '
            'campaign sets)')
    return specs, include_baseline


def write_document(doc, out_dir, stem, stamp):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f'{stem}_{stamp}.json')
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(doc, fh, indent=2)
    return path


def console_report(doc):
    m = doc['meta']
    print(f'band source (campaign): {m["band_campaign"]}')
    print(f'proteome cap {m["protein_content"]:.3f}  housekeeping '
          f'{m["housekeeping"]:.3f}  F_flex {m["F_flex"]:.3f}')
    for t in doc['tiles']:
        kind = 'baseline' if t['is_baseline'] \
            else f'{t["campaign"]} trial {t["trial_number"]}'
        slack = next(c['value'] for c in t['children'] if c['piece'] == 'slack')
        print(f'[{t["label"]}] {kind}: Phi_M {t["Phi_M"]:.4f}  phi_T '
              f'{t["phi_T"]:.4f}  d {t["burden_factor"]:.2f}  slack {slack:.4f}'
              + ('  !! ' + t['warning'] if 'warning' in t else ''))


def build_parser():
    ap = argparse.ArgumentParser(
        description='Proteome-allocation Voronoi treemap (small multiples).')
    ap.add_argument('--set', dest='sets', action='append', nargs=3,
                    metavar=('LABEL', 'CAMPAIGN', 'TRIAL'), default=None,
                    help='a tile to plot; repeatable, in grid order after the '
                         'baseline. TRIAL is an int, "best", or "best:COL".')
    ap.add_argument('--no-baseline', action='store_true',
                    help='drop the scenario-A baseline tile')
    ap.add_argument('--out-dir', default=RESULTS_DIR)
    ap.add_argument('--stem', default=None)
    ap.add_argument('--cols', type=int, default=None,
                    help='grid columns (default: chosen from the tile count)')
    ap.add_argument('--scale', type=float, default=2.0,
                    help='PNG device scale factor (default 2)')
    ap.add_argument('--node', default=None,
                    help='path to the voronoi-treemaps node.exe '
                         '(overrides VORONOI_NODE)')
    ap.add_argument('--no-render', action='store_true',
                    help='write only the JSON (skip the Node/Puppeteer render)')
    return ap


def resolve_node(cli_node=None):
    """Path to the voronoi-treemaps node.exe: --node, else $VORONOI_NODE,
    else DEFAULT_VORONOI_NODE. Errors clearly if it does not exist."""
    node = cli_node or os.environ.get('VORONOI_NODE') or DEFAULT_VORONOI_NODE
    if not os.path.isfile(node):
        raise FileNotFoundError(
            f'voronoi-treemaps node not found at {node!r}. Create the env '
            '(conda create -n voronoi-treemaps -c conda-forge nodejs; '
            'npm --prefix plots/voronoi install) or pass --node / set '
            '$VORONOI_NODE. Use --no-render to write only the JSON.')
    return node


def render(json_path, args):
    """Invoke Stage 2 (Node/Puppeteer) on a written allocation JSON."""
    node = resolve_node(args.node)
    out_stem = os.path.splitext(json_path)[0]   # <out-dir>/<stem>_<stamp>
    cmd = [node, RENDER_SCRIPT, '--doc', json_path, '--out', out_stem,
           '--scale', str(args.scale)]
    if args.cols:
        cmd += ['--cols', str(args.cols)]
    print('rendering: ' + ' '.join(f'"{c}"' if ' ' in c else c for c in cmd))
    subprocess.run(cmd, check=True)
    print(f'wrote {out_stem}.png / .pdf')
    return out_stem


def main(argv=None):
    args = build_parser().parse_args(argv)
    specs, include_baseline = resolve_specs(args)
    sets, _band, band_campaign = ps.build_sets(specs, include_baseline)
    doc = build_document(sets, band_campaign)
    stem = args.stem or (
        f'{os.path.basename(specs[0][1]).replace(".csv", "")}_proteome_voronoi')
    stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
    json_path = write_document(doc, args.out_dir, stem, stamp)
    console_report(doc)
    print(f'wrote {json_path}')
    if args.no_render:
        return json_path
    render(json_path, args)
    return json_path


if __name__ == '__main__':
    main()
