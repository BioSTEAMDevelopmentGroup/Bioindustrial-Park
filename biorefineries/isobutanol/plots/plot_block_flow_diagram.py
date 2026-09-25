#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Simplified block flow diagram of the corn ethanol-isobutanol biorefinery.

Two panels plus a legend, in the layout of the group's earlier corn
ethanol-isobutanol block flow diagram (fermentation process on top,
separation process below, cream / light-blue panel tints, italic
underlined boundary streams, bold products, dashed co-production units):

* **A** — fermentation process: the saccharified slurry from the corn
  dry-grind process is split between the initial-feed and spike-feed
  trains (multi-effect evaporation ``F301`` / ``F302``, dilution-water
  mixing ``M301`` / ``M302``), fermented in the aerated fermentor
  (``V406``; fed-batch, or batch when no spikes are fed, e.g. scenario B;
  compressed air ``K330``, seed yeast), and the fermentation off-gas is
  scrubbed (``V409``) so the stripped alcohols rejoin the broth (``MX8``)
  ahead of the separation process.
* **B** — separation process of the default build at its baseline gate
  (``S201.split = 1.0``: all broth to the solvent-free IBO/EtOH train of
  ``separations.create_IBO_EtOH_separation_system``): beer column
  (``D101``) -> ethanol rectifier (``D102``) -> molecular sieves
  (``MS201``) -> denaturant mixing -> ethanol. The rectifier bottoms carry
  the isobutanol to the stripper (``D103``; bottoms water to process-water
  reuse), the heteroazeotrope is crossed by the decanter (aqueous phase
  back to the rectifier feed) and the drying column (``D104``; azeotrope
  back to the decanter) yields isobutanol (> 99.9 wt%). Units drawn dashed
  run only when isobutanol is produced (at zero isobutanol ``D103`` passes
  its feed to its bottoms and the decanter / ``D104`` idle; none of the
  three is designed or costed). The beer-column stillage goes to corn's
  DDGS train (centrifugation, a thin-stillage purge to wastewater
  treatment, multi-effect evaporation, oil separation, drum drying with a
  thermal-oxidizer exhaust).

Deliberately left out as block-level detail: the idle ethanol-primary
train (zero flow at the baseline gate), heat exchangers, pumps and tanks,
the molecular-sieve regenerant recycle, and the wastewater-treatment /
boiler-turbogenerator / utility facilities.

Suggested caption abbreviation list: DDGS, distillers' dried grains with
solubles; WWT, wastewater treatment (also defined in the legend).

Standalone on purpose: imports only matplotlib (NOT
``biorefineries.isobutanol``), so it runs no simulation and touches no
numba cache. Sized to the 180 mm double-column width; Arial embedded as
TrueType (fonttype 42); text 6-7 pt, panel letters 8 pt bold (capitals,
as in the reference).

Run::

    python plots/plot_block_flow_diagram.py

Writes ``block_flow_diagram.{pdf,svg,png}`` under ``plots/results/``.
"""

import os

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.textpath import TextPath, text_to_path

__all__ = ('make_figure',)

# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------

MM = 1 / 25.4            # inch per mm
PT_PER_MM = 72 / 25.4
FIG_W = 180.0            # mm (double-column width)
FIG_H = 158.0            # mm

# Palette sampled from the reference block flow diagram.
C_FERM_BOX = '#F8E1A9'   # fermentation process / unit fill
C_SEP_BOX = '#ABDEE5'    # separation process / unit fill
C_FERM_BG = '#FFFBF3'    # panel a tint
C_SEP_BG = '#F4FBFB'     # panel b tint
C_INK = '#000000'

FS_BOX = 7.0             # unit names (Nature: 5-7 pt)
FS_STREAM = 6.5          # stream labels
FS_SMALL = 6.0           # internal-stream labels in tight gaps
FS_PRODUCT = 7.0         # product names (bold)
FS_PANEL = 8.0           # panel letters and titles (bold)
LW_BOX = 0.9             # box outline (pt)
LW_STREAM = 0.8          # stream line (pt)
LW_RULE = 0.45           # underline rule (pt)
HEAD_L, HEAD_W = 1.6, 1.1   # arrowhead length / half-width (mm)
LABEL_OFF = 1.3             # label offset from a bare stream line (mm)
HEAD_OFF = HEAD_W + 0.9     # ... from a line beside an arrowhead (mm)
DASH_BOX = (0, (3.0, 1.6))

# Arial vertical metrics (fractions of the em) used to anchor underlined
# lines on their baselines.
ASCENT, DESCENT = 0.80, 0.20
LINE_PITCH = 1.2
UNDERLINE_DROP = 0.13


def _set_style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.default': 'regular',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
        'text.color': C_INK,
    })


def _mm(pt):
    return pt / PT_PER_MM


# --------------------------------------------------------------------------
# Drawing helpers (all coordinates in mm, origin bottom-left)
# --------------------------------------------------------------------------

class Diagram:
    """Draws blocks, streams and labels; underlines are fitted after layout."""

    def __init__(self, ax):
        self.ax = ax
        self._underlined = []    # (text, x anchor, ha, baseline, size)

    # ---- boxes -----------------------------------------------------------
    def box(self, cx, cy, w, h, text, fill, dashed=False, bold=False):
        """Process / unit block centred on (cx, cy); returns its edges."""
        self.ax.add_patch(Rectangle(
            (cx - w / 2, cy - h / 2), w, h, facecolor=fill, edgecolor=C_INK,
            linewidth=LW_BOX, linestyle=DASH_BOX if dashed else 'solid',
            zorder=3))
        self.ax.text(cx, cy, text, ha='center', va='center',
                     fontsize=FS_BOX, fontweight='bold' if bold else 'normal',
                     linespacing=1.15, zorder=4)
        return dict(l=cx - w / 2, r=cx + w / 2, b=cy - h / 2, t=cy + h / 2,
                    cx=cx, cy=cy)

    # ---- streams ---------------------------------------------------------
    def stream(self, pts):
        """Orthogonal polyline through ``pts`` ending in a filled head."""
        ax = self.ax
        if len(pts) > 2:
            xs, ys = zip(*pts[:-1])
            ax.add_line(Line2D(xs, ys, color=C_INK, lw=LW_STREAM,
                               solid_joinstyle='miter',
                               solid_capstyle='butt', zorder=2))
        ax.add_patch(FancyArrowPatch(
            pts[-2], pts[-1], zorder=2, color=C_INK, lw=LW_STREAM,
            shrinkA=0, shrinkB=0, mutation_scale=1.0,
            arrowstyle=(f'-|>,head_length={HEAD_L * PT_PER_MM:.3f},'
                        f'head_width={HEAD_W * PT_PER_MM:.3f}'),
            joinstyle='miter', capstyle='butt'))

    # ---- labels ----------------------------------------------------------
    def label(self, x, y, text, ha='center', va='center', fs=FS_STREAM):
        """Plain label of an intermediate stream."""
        return self.ax.text(x, y, text, ha=ha, va=va, fontsize=fs,
                            linespacing=1.15, zorder=5)

    def boundary_label(self, x, y, text, ha='center', va='center',
                       fs=FS_STREAM):
        """Italic, underlined input / outlet stream, one artist per line.

        ``va`` anchors the whole block (top / center / bottom of the ink
        envelope); each line is placed on its own baseline so that its rule
        can sit a fixed distance below it.
        """
        lines = text.split('\n')
        em = _mm(fs)
        pitch = LINE_PITCH * em
        height = (len(lines) - 1) * pitch + (ASCENT + DESCENT) * em
        top = {'top': y, 'center': y + height / 2, 'bottom': y + height}[va]
        for i, line in enumerate(lines):
            baseline = top - ASCENT * em - i * pitch
            self.ax.text(x, baseline, line, ha=ha, va='baseline',
                         fontsize=fs, fontstyle='italic', zorder=5)
            self._underlined.append((line, x, ha, baseline, fs))

    def product_label(self, x, y, name, spec=None):
        self.ax.text(x, y + (1.5 if spec else 0), name, ha='left',
                     va='center', fontsize=FS_PRODUCT, fontweight='bold',
                     zorder=5)
        if spec:
            self.ax.text(x, y - 1.9, spec, ha='left', va='center',
                         fontsize=FS_STREAM, zorder=5)

    def finalize(self, fig):
        """Draw each underline rule across the ink extent of its line.

        Extents come from the unhinted font outlines (the vector glyphs
        written to the PDF), not from a screen-resolution text layout.
        """
        family = plt.rcParams['font.sans-serif'][0]
        for line, x, ha, baseline, fs in self._underlined:
            prop = FontProperties(family=family, style='italic', size=fs)
            advance, _, _ = text_to_path.get_text_width_height_descent(
                line, prop, ismath=False)
            left = x - _mm(advance) * {'left': 0.0, 'center': 0.5,
                                       'right': 1.0}[ha]
            ink = TextPath((0, 0), line, size=fs, prop=prop).get_extents()
            # a slanted glyph's rightmost ink is at cap height; end the rule
            # under the glyph's foot instead
            x0 = left + _mm(ink.x0)
            x1 = left + _mm(ink.x1) - 0.15 * _mm(fs)
            y_rule = baseline - UNDERLINE_DROP * _mm(fs)
            self.ax.add_line(Line2D([x0, x1], [y_rule, y_rule], color=C_INK,
                                    lw=LW_RULE, zorder=5,
                                    solid_capstyle='butt'))


def panel_background(ax, x0, y0, x1, y1, color):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=color,
                           edgecolor='none', zorder=0))


def panel_title(ax, y1, letter, title):
    ax.text(2.5, y1 - 3.2, letter, fontsize=FS_PANEL, fontweight='bold',
            va='center')
    ax.text(6.0, y1 - 3.2, title, fontsize=FS_PANEL, fontweight='bold',
            va='center')


# --------------------------------------------------------------------------
# Panel a — fermentation process
# --------------------------------------------------------------------------

A_Y0, A_Y1 = 93.5, 158.0     # panel a vertical extent (mm)


def draw_panel_a(d):
    ax = d.ax
    panel_background(ax, 0, A_Y0, FIG_W, A_Y1, C_FERM_BG)
    panel_title(ax, A_Y1, 'A', 'Fermentation process')

    fill = C_FERM_BOX
    y_top, y_bot = 132.0, 110.0          # the two feed trains
    y_mid = (y_top + y_bot) / 2
    bh = 10.0

    d.label(12.0, y_mid, 'saccharified\nslurry (from\ncorn dry-grind\nprocess)')
    spl = d.box(32, y_mid, 14, bh, 'splitter', fill)
    d.stream([(22.0, y_mid), (spl['l'], y_mid)])

    ev1 = d.box(55, y_top, 22, bh, 'multi-effect\nevaporation', fill)
    ev2 = d.box(55, y_bot, 22, bh, 'multi-effect\nevaporation', fill)
    x_split = spl['cx'] + 3
    d.stream([(x_split, spl['t']), (x_split, y_top), (ev1['l'], y_top)])
    d.stream([(x_split, spl['b']), (x_split, y_bot), (ev2['l'], y_bot)])

    mx1 = d.box(78, y_top, 14, bh, 'mixing', fill)
    mx2 = d.box(78, y_bot, 14, bh, 'mixing', fill)
    d.stream([(ev1['r'], y_top), (mx1['l'], y_top)])
    d.stream([(ev2['r'], y_bot), (mx2['l'], y_bot)])

    # evaporator vapours and dilution water (boundary streams)
    gap = 4.5
    d.stream([(ev1['cx'], ev1['t']), (ev1['cx'], ev1['t'] + gap)])
    d.boundary_label(ev1['cx'], ev1['t'] + gap + 0.8,
                     'water & other\nvolatiles (to WWT)', va='bottom')
    d.stream([(ev2['cx'], ev2['b']), (ev2['cx'], ev2['b'] - gap)])
    d.boundary_label(ev2['cx'], ev2['b'] - gap - 0.8,
                     'water & other\nvolatiles (to WWT)', va='top')
    d.stream([(mx1['cx'], mx1['t'] + gap), (mx1['cx'], mx1['t'])])
    d.boundary_label(mx1['cx'], mx1['t'] + gap + 0.8, 'dilution water',
                     va='bottom')
    d.stream([(mx2['cx'], mx2['b'] - gap), (mx2['cx'], mx2['b'])])
    d.boundary_label(mx2['cx'], mx2['b'] - gap - 0.8, 'dilution water',
                     va='top')

    # fermentation: initial feed straight into the left side, spike feed
    # from below
    fe = d.box(116, 126.0, 25, 16, 'fermentation\n(fed-batch\nor batch)', fill)
    d.stream([(mx1['r'], y_top), (fe['l'], y_top)])
    d.label((mx1['r'] + fe['l']) / 2, y_top + 1.9, 'initial feed')
    x_spike = fe['l'] + 5.0
    d.stream([(mx2['r'], y_bot), (x_spike, y_bot), (x_spike, fe['b'])])
    d.label((mx2['r'] + x_spike) / 2, y_bot + 1.9, 'spike feed')
    x_yeast = fe['cx'] + 5.0
    d.stream([(x_yeast, fe['b'] - 5.5), (x_yeast, fe['b'])])
    d.boundary_label(x_yeast, fe['b'] - 6.3, 'yeast', va='top')

    # aeration
    y_aux = 143.5
    co = d.box(fe['l'] + 6.5, y_aux, 19, 8, 'compression', fill)
    d.stream([(co['cx'], co['t'] + 4.0), (co['cx'], co['t'])])
    d.boundary_label(co['cx'], co['t'] + 4.8, 'air', va='bottom')
    d.stream([(co['cx'], co['b']), (co['cx'], fe['t'])])

    # vent scrubbing: fermentation off-gas -> scrubber; the scrubber bottoms
    # (recovered alcohols) rejoin the broth ahead of the separation process
    vs = d.box(147, y_aux, 18, 8, 'vent\nscrubbing', fill)
    x_vent = fe['r'] - 4.0
    d.stream([(x_vent, fe['t']), (x_vent, vs['cy']), (vs['l'], vs['cy'])])
    d.label(x_vent + LABEL_OFF, (fe['t'] + vs['cy']) / 2 - 0.4, 'vented\ngases',
            ha='left')
    d.stream([(vs['cx'], vs['t'] + 4.0), (vs['cx'], vs['t'])])
    d.boundary_label(vs['cx'], vs['t'] + 4.8, 'wash water', va='bottom')
    d.stream([(vs['r'], vs['cy']), (vs['r'] + 4.5, vs['cy'])])
    d.boundary_label(vs['r'] + 5.3, vs['cy'], 'scrubbed\ngases', ha='left')

    sep = d.box(167.0, fe['cy'], 19, 12, 'Separation\nprocess\n(to B)',
                C_SEP_BOX, bold=True)
    d.stream([(fe['r'], fe['cy']), (sep['l'], fe['cy'])])
    d.label((fe['r'] + vs['cx']) / 2, fe['cy'] - 5.0,
            'broth\n(crude ethanol\n& isobutanol)')
    d.stream([(vs['cx'], vs['b']), (vs['cx'], fe['cy'])])
    d.label(vs['cx'] + LABEL_OFF, (vs['b'] + sep['t']) / 2 + 0.3,
            'scrubber bottoms\n(recovered alcohols)', ha='left')


# --------------------------------------------------------------------------
# Panel b — separation process
# --------------------------------------------------------------------------

B_Y0, B_Y1 = 13.0, 92.0


def draw_panel_b(d):
    ax = d.ax
    panel_background(ax, 0, B_Y0, FIG_W, B_Y1, C_SEP_BG)
    panel_title(ax, B_Y1, 'B', 'Separation process')

    fill = C_SEP_BOX
    bh = 10.0
    r1, r2, r3, r4 = 78.0, 60.0, 41.0, 25.0    # row centres
    x1, x2, x3, x4 = 40.0, 66.0, 94.0, 122.0   # column centres
    x_prod = 144.0                             # product-label column
    gap = 4.5

    def to_product(unit, y):
        d.stream([(unit['r'], y), (x_prod - 1.5, y)])

    # ---- row 1: ethanol ---------------------------------------------------
    fp = d.box(14, r1, 19, bh + 1, 'Fermentation\nprocess\n(from A)',
               C_FERM_BOX, bold=True)
    bc = d.box(x1, r1, 18, bh, 'distillation\n(beer column)', fill)
    rc = d.box(x2, r1, 18, bh, 'distillation\n(rectifier)', fill)
    ms = d.box(x3, r1, 18, bh, 'molecular\nsieves', fill)
    mx = d.box(x4, r1, 14, bh, 'mixing', fill)
    d.stream([(fp['r'], r1), (bc['l'], r1)])
    d.stream([(bc['r'], r1), (rc['l'], r1)])
    d.stream([(rc['r'], r1), (ms['l'], r1)])
    d.stream([(ms['r'], r1), (mx['l'], r1)])
    d.stream([(mx['cx'], mx['t'] + gap), (mx['cx'], mx['t'])])
    d.boundary_label(mx['cx'], mx['t'] + gap + 0.8, 'denaturant',
                     va='bottom')
    to_product(mx, r1)
    d.product_label(x_prod, r1, 'ethanol', '(denatured)')

    # ---- row 2: isobutanol (co-production only) ---------------------------
    st = d.box(x2, r2, 18, bh, 'distillation\n(stripper)', fill, dashed=True)
    dc = d.box(x3, r2, 16, bh, 'decanting', fill, dashed=True)
    dr = d.box(x4, r2, 20, bh, 'distillation\n(drying column)', fill,
               dashed=True)
    x_rb = rc['cx'] - 4.0
    d.stream([(x_rb, rc['b']), (x_rb, st['t'])])
    d.label(x_rb - HEAD_OFF, (rc['b'] + st['t']) / 2, 'isobutanol-\nrich bottoms',
            ha='right', fs=FS_SMALL)
    d.stream([(st['r'], r2), (dc['l'], r2)])
    d.stream([(dc['r'], r2), (dr['l'], r2)])
    d.label((dc['r'] + dr['l'] - HEAD_L) / 2, r2 + 2.2, 'organic',
            fs=FS_SMALL)
    to_product(dr, r2)
    d.product_label(x_prod, r2, 'isobutanol', '(>99.9 wt%)')
    # aqueous phase back to the rectifier
    y_aq = (rc['b'] + st['t']) / 2
    x_aq_up = rc['cx'] + 4.0
    d.stream([(dc['cx'], dc['t']), (dc['cx'], y_aq), (x_aq_up, y_aq),
              (x_aq_up, rc['b'])])
    d.label((x_aq_up + dc['cx']) / 2, y_aq + 1.6, 'aqueous phase',
            fs=FS_SMALL)
    # azeotrope back to the decanter
    y_az = dr['b'] - 3.5
    d.stream([(dr['cx'], dr['b']), (dr['cx'], y_az), (dc['cx'], y_az),
              (dc['cx'], dc['b'])])
    d.label((dr['cx'] + dc['cx']) / 2, y_az - 1.7, 'azeotrope', fs=FS_SMALL)
    # stripper bottoms: water reused as process water
    d.stream([(st['cx'], st['b']), (st['cx'], st['b'] - 4.0)])
    d.boundary_label(st['cx'], st['b'] - 4.8, 'water (reused as process water)',
                     va='top')

    # ---- rows 3-4: stillage to DDGS and corn oil ------------------------
    ce = d.box(x1, r3, 18, bh - 1, 'centrifugation', fill)
    d.stream([(bc['cx'], bc['b']), (bc['cx'], ce['t'])])
    d.label(bc['cx'] - LABEL_OFF, (bc['b'] + ce['t']) / 2, 'whole\nstillage',
            ha='right')
    mx6 = d.box(x3, r3, 14, bh - 1, 'mixing', fill)
    dy = d.box(x4, r3, 16, bh - 1, 'drying', fill)
    d.stream([(ce['r'], r3), (mx6['l'], r3)])
    d.label((ce['r'] + mx6['l']) / 2, r3 - 1.9, 'wet cake')
    d.stream([(mx6['r'], r3), (dy['l'], r3)])
    to_product(dy, r3)
    d.product_label(x_prod, r3, 'DDGS')
    # dryer: natural gas and air in, exhaust out
    x_ng, x_ex = dy['l'] + 4.0, dy['r'] - 4.0
    d.stream([(x_ng, dy['b'] - 6.5), (x_ng, dy['b'])])
    d.boundary_label(x_ng - HEAD_OFF, dy['b'] - 3.6, 'natural gas\n& air',
                     ha='right')
    d.stream([(x_ex, dy['b']), (x_ex, dy['b'] - 6.5)])
    d.boundary_label(x_ex + HEAD_OFF, dy['b'] - 3.6,
                     'exhaust (to thermal\noxidation)', ha='left')

    ev = d.box(x2 - 2.0, r4, 22, bh - 1, 'multi-effect\nevaporation', fill)
    os_ = d.box(x3, r4, 18, bh - 1, 'oil\nseparation', fill)
    d.stream([(ce['cx'], ce['b']), (ce['cx'], r4), (ev['l'], r4)])
    d.label(ce['cx'] + LABEL_OFF, (ce['b'] + r4) / 2 + 0.6, 'thin\nstillage',
            ha='left')
    d.stream([(ev['r'], r4), (os_['l'], r4)])
    d.label((ev['r'] + os_['l'] - HEAD_L) / 2, r4 + 2.2, 'syrup',
            fs=FS_SMALL)
    d.stream([(os_['cx'], os_['t']), (os_['cx'], mx6['b'])])
    d.label(os_['cx'] - HEAD_OFF, (os_['t'] + mx6['b']) / 2, 'de-oiled\nsyrup',
            ha='right', fs=FS_SMALL)
    to_product(os_, r4)
    d.product_label(x_prod, r4, 'corn oil')
    # wastes to wastewater treatment
    d.stream([(ev['cx'], ev['b']), (ev['cx'], ev['b'] - 3.5)])
    d.boundary_label(ev['cx'], ev['b'] - 4.3, 'condensate (to WWT)',
                     va='top')
    d.stream([(ce['cx'], r4), (ce['cx'], r4 - 7.0)])
    d.boundary_label(ce['cx'] - HEAD_OFF, r4 - 4.8, 'thin-stillage\npurge (to WWT)',
                     ha='right')


# --------------------------------------------------------------------------
# Legend
# --------------------------------------------------------------------------

def draw_legend(d):
    ax = d.ax
    y0, y1 = 1.0, 11.5
    ax.add_patch(Rectangle((1.5, y0), FIG_W - 3.0, y1 - y0, facecolor='white',
                           edgecolor=C_INK, lw=LW_BOX, zorder=1))
    yc = (y0 + y1) / 2
    ax.text(4.5, yc, 'Legend', fontsize=FS_PANEL, fontweight='bold',
            va='center')
    # two-tone process swatch
    x, w, h = 19.0, 11.0, 6.0
    ax.add_patch(Rectangle((x, yc - h / 2), w / 2, h, facecolor=C_FERM_BOX,
                           edgecolor='none', zorder=3))
    ax.add_patch(Rectangle((x + w / 2, yc - h / 2), w / 2, h,
                           facecolor=C_SEP_BOX, edgecolor='none', zorder=3))
    ax.add_patch(Rectangle((x, yc - h / 2), w, h, facecolor='none',
                           edgecolor=C_INK, lw=LW_BOX, zorder=4))
    ax.text(x + w + 1.5, yc, 'process or unit', ha='left', va='center',
            fontsize=FS_STREAM, zorder=5)
    # intermediate stream
    d.stream([(53.0, yc - 1.6), (71.0, yc - 1.6)])
    d.label(62.0, yc + 1.4, 'intermediate stream')
    # boundary stream
    d.boundary_label(84.5, yc, 'input or\noutlet stream')
    # product
    ax.text(100.5, yc, 'product', fontsize=FS_PRODUCT, fontweight='bold',
            va='center', ha='center')
    # dashed unit
    x, w, h = 110.0, 38.0, 7.0
    ax.add_patch(Rectangle((x, yc - h / 2), w, h, facecolor='none',
                           edgecolor=C_INK, lw=LW_BOX, linestyle=DASH_BOX,
                           zorder=3))
    ax.text(x + w / 2, yc, 'operates only with\nisobutanol co-production',
            ha='center', va='center', fontsize=FS_STREAM, linespacing=1.15,
            zorder=5)
    # abbreviation
    ax.text(152.0, yc, 'WWT, wastewater\ntreatment', ha='left', va='center',
            fontsize=FS_STREAM, linespacing=1.15)


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

def make_figure():
    _set_style()
    fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, FIG_H)
    ax.set_aspect('equal')
    ax.axis('off')
    d = Diagram(ax)
    draw_panel_a(d)
    draw_panel_b(d)
    draw_legend(d)
    d.finalize(fig)
    return fig


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(here, 'results')
    os.makedirs(out_dir, exist_ok=True)
    fig = make_figure()
    stem = os.path.join(out_dir, 'block_flow_diagram')
    fig.savefig(stem + '.pdf')
    fig.savefig(stem + '.svg')
    fig.savefig(stem + '.png', dpi=600)
    plt.close(fig)
    print('wrote', stem + '.{pdf,svg,png}')


if __name__ == '__main__':
    main()
