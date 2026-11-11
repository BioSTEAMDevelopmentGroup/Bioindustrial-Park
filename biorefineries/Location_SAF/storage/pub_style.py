"""
pub_style.py
============
Reusable publication figure style for biomass storage diagnostics.
A4-ready, Arial font, white/transparent background.

USAGE:
    from pub_style import apply_pub_style, style_axes, save_fig

    fig, axes = plt.subplots(...)
    apply_pub_style()                  # call once before creating figures
    for ax in axes.flat:
        style_axes(ax)                 # call on every axes object

TICK ARCHITECTURE:
    Bottom/left (primary)  → inout ticks, labeled
    Top/right              → twin axes owned independently, in-only ticks, no labels
    Minor ticks            → off everywhere

    Twin axes are the correct approach: setting top/right via tick_params on the
    primary axes can be silently overwritten by plotting calls (violinplot, hist,
    etc.).  The twins are created with no data, no labels, no independent scale —
    they exist only to own their spines' tick appearance.
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.font_manager as fm
import numpy as np


# =============================================================================
# FONT RESOLUTION
# =============================================================================
# Prefer Arial (standard on Windows/macOS); fall back to Liberation Sans
# (metric-identical, available on most Linux systems) then DejaVu Sans.

def _resolve_font():
    available = {f.name for f in fm.fontManager.ttflist}
    for candidate in ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]:
        if candidate in available:
            return candidate
    return "sans-serif"

FONT_FAMILY = _resolve_font()


# =============================================================================
# STYLE CONSTANTS
# =============================================================================
FONT_SIZE_TITLE    = 13    # panel titles  e.g. "(a) Total Cost"
FONT_SIZE_LABEL    = 12    # axis labels
FONT_SIZE_TICK     = 11    # tick labels
FONT_SIZE_LEGEND   = 10
FONT_SIZE_ANNOT    = 10    # in-panel annotations
FONT_SIZE_SUPTITLE = 14

TICK_LENGTH_MAJOR  = 5     # pts — applies to both primary (inout) and twin (in)
TICK_LENGTH_MINOR  = 0     # minor ticks disabled everywhere
LINE_WIDTH         = 1.2   # spine / tick line width


# =============================================================================
# APPLY GLOBAL RC PARAMS  (call once per script)
# =============================================================================

def apply_pub_style():
    """
    Sets global rcParams: Arial font, white background, publication font sizes.
    Tick direction is NOT set here — it is applied per-axes via style_axes()
    to avoid rcParam bleed between primary and twin axes.
    """
    plt.rcParams.update({
        # --- font ---
        "font.family":         FONT_FAMILY,
        "font.size":           FONT_SIZE_TICK,
        "axes.titlesize":      FONT_SIZE_TITLE,
        "axes.labelsize":      FONT_SIZE_LABEL,
        "xtick.labelsize":     FONT_SIZE_TICK,
        "ytick.labelsize":     FONT_SIZE_TICK,
        "legend.fontsize":     FONT_SIZE_LEGEND,

        # --- background ---
        "figure.facecolor":    "white",
        "axes.facecolor":      "white",
        "savefig.facecolor":   "white",
        "savefig.transparent": False,

        # --- spines ---
        "axes.spines.top":     True,
        "axes.spines.right":   True,
        "axes.edgecolor":      "black",
        "axes.linewidth":      LINE_WIDTH,

        # --- ticks: leave direction unset globally; style_axes handles it ---
        "xtick.major.size":    TICK_LENGTH_MAJOR,
        "ytick.major.size":    TICK_LENGTH_MAJOR,
        "xtick.minor.size":    TICK_LENGTH_MINOR,
        "ytick.minor.size":    TICK_LENGTH_MINOR,
        "xtick.major.width":   LINE_WIDTH,
        "ytick.major.width":   LINE_WIDTH,
        "xtick.color":         "black",
        "ytick.color":         "black",

        # --- grid (off by default; enable per-panel via style_axes grid=True) ---
        "axes.grid":           False,

        # --- legend ---
        "legend.frameon":      True,
        "legend.framealpha":   0.9,
        "legend.edgecolor":    "0.8",

        # --- lines ---
        "lines.linewidth":     1.8,
    })


# =============================================================================
# PER-AXES STYLE
# =============================================================================

def style_axes(ax, grid=False):
    """
    Apply publication tick styling to one axes object.

    Strategy
    --------
    Primary axes  → bottom + left ticks only, direction='inout', labeled.
    Twin axes     → top (twiny) and right (twinx) ticks, direction='in',
                    no labels, no independent scale, transparent background.

    This keeps the primary tick direction completely isolated from top/right,
    so plotting calls on the primary axes can never overwrite the twin ticks.

    Parameters
    ----------
    ax   : matplotlib Axes
    grid : bool  — draw a light dashed grid on the primary axes (default False)

    Returns
    -------
    ax_top, ax_right : the two twin axes (store them if you need to adjust
                       their limits to match the primary after the fact, though
                       that is rarely necessary since they carry no data).
    """

    # ── 1. PRIMARY axes: bottom + left only, inout ───────────────────────────
    ax.tick_params(
        axis="both", which="both",
        top=False, right=False,        # twins own those sides
        bottom=True, left=True,
        direction="inout",
        length=TICK_LENGTH_MAJOR,
        width=LINE_WIDTH,
        labelsize=FONT_SIZE_TICK,
    )
    ax.minorticks_off()

    for spine in ax.spines.values():
        spine.set_linewidth(LINE_WIDTH)

    # ── 2. TWIN-X (shares x-axis) → owns top spine ───────────────────────────
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())     # match primary x limits

    # hide everything except the top tick marks
    ax_top.tick_params(
        axis="x", which="major",
        top=True, bottom=False,
        direction="in",
        length=TICK_LENGTH_MAJOR,
        width=LINE_WIDTH,
        labeltop=False,
    )
    ax_top.tick_params(axis="x", which="minor", top=False)
    ax_top.set_xticklabels([])
    ax_top.xaxis.set_ticks_position("top")

    # make twin background transparent so the primary shows through
    ax_top.patch.set_visible(False)

    for spine_name, spine in ax_top.spines.items():
        if spine_name == "top":
            spine.set_linewidth(LINE_WIDTH)
            spine.set_visible(True)
        else:
            spine.set_visible(False)

    # ── 3. TWIN-Y (shares y-axis) → owns right spine ─────────────────────────
    ax_right = ax.twinx()
    ax_right.set_ylim(ax.get_ylim())   # match primary y limits

    ax_right.tick_params(
        axis="y", which="major",
        right=True, left=False,
        direction="in",
        length=TICK_LENGTH_MAJOR,
        width=LINE_WIDTH,
        labelright=False,
    )
    ax_right.tick_params(axis="y", which="minor", right=False)
    ax_right.set_yticklabels([])
    ax_right.yaxis.set_ticks_position("right")

    ax_right.patch.set_visible(False)

    for spine_name, spine in ax_right.spines.items():
        if spine_name == "right":
            spine.set_linewidth(LINE_WIDTH)
            spine.set_visible(True)
        else:
            spine.set_visible(False)

    # ── 4. Optional grid on primary ──────────────────────────────────────────
    if grid:
        ax.grid(True, which="major", color="0.88", linewidth=0.6, linestyle="--")
        ax.set_axisbelow(True)

    return ax_top, ax_right


# =============================================================================
# CONVENIENCE: SYNC TWIN LIMITS  (call after all data is plotted)
# =============================================================================

def sync_twins(ax, ax_top, ax_right):
    """
    Re-sync twin axes limits to the primary after data has been plotted
    and matplotlib has auto-scaled the primary axes.

    Usage:
        ax_top, ax_right = style_axes(ax)
        ax.violinplot(...)          # primary auto-scales
        sync_twins(ax, ax_top, ax_right)
    """
    ax_top.set_xlim(ax.get_xlim())
    ax_right.set_ylim(ax.get_ylim())


# =============================================================================
# CONVENIENCE: SAVE FIGURE
# =============================================================================

def save_fig(fig, path, dpi=300, transparent=False):
    """Save with consistent settings."""
    fig.savefig(
        path, dpi=dpi, bbox_inches="tight",
        facecolor="white" if not transparent else "none",
        transparent=transparent,
    )
    print(f"Saved: {path}")
