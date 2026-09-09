#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 4 — Revenue split Sankey
Caption: Revenue at SAF MPSP; coproducts at market.
"""

from pathlib import Path
import sys

# --- paths: script lives in biorefineries/covercress/ ---
BIOREFINERIES = Path(__file__).resolve().parents[1]
if str(BIOREFINERIES) not in sys.path:
    sys.path.insert(0, str(BIOREFINERIES))

import pandas as pd
import plotly.graph_objects as go

# ---------------------------------------------------------------------------
# Load covercress system (skip auto run_tea at bottom of system file)
# ---------------------------------------------------------------------------
SYSTEM_FILE = BIOREFINERIES / "covercress" / "sys_protein_ext_with_combined_ext_centri.py"

with open(SYSTEM_FILE) as f:
    code = f.read()
code = code.replace("MSP_SAF, MFPP = run_tea()", "pass  # skipped for figure script")
code = code.replace("MSP_SAF = run_tea()", "pass")

_ns = {"__name__": "__figure4__", "__file__": str(SYSTEM_FILE)}
exec(compile(code, str(SYSTEM_FILE), "exec"), _ns)

import biosteam as bst
from biorefineries.covercress.process_settings import price

oil_extraction_sys = _ns["oil_extraction_sys"]
saf_tea = _ns["saf_tea"]
set_market_prices = _ns["set_market_prices"]
num_sims = _ns.get("num_sims", 3)
num_solve_tea = _ns.get("num_solve_tea", 3)
s = _ns["s"]

PRODUCTS = {
    "SAF": s.SAF,
    "Protein": s.protein_isolate,
    "Naphtha": s.naphtha_product,
    "Green diesel": s.green_diesel,
    "Propane": s.propane,
}

# CABBI-style colors (order: SAF, protein, naphtha, GD, propane)
COLORS = ["#63C6CE", "#7BBD84", "#F7C652", "#B97A57", "#94948C"]
TOTAL_COLOR = "#734A8C"


# ---------------------------------------------------------------------------
# TEA: MPSP on SAF, coproducts stay at market
# ---------------------------------------------------------------------------
def solve_saf_mpsp():
    set_market_prices()
    for _ in range(num_sims):
        oil_extraction_sys.simulate()
    for _ in range(num_solve_tea):
        s.SAF.price = saf_tea.solve_price(s.SAF)
    print(f"SAF MPSP: ${s.SAF.price:.4f}/kg  (market SAF: ${price['SAF']:.4f}/kg)")
    print(f"NPV @ MPSP: ${saf_tea.NPV:,.0f}")


def revenue_table(tea):
    """Annual product revenue [MM$/yr] and share [%]."""
    op_hr = tea.operating_hours
    rows = []
    for name, stream in PRODUCTS.items():
        mm_yr = stream.cost * op_hr / 1e6
        rows.append({
            "Product": name,
            "Price [$/kg]": stream.price,
            "Flow [kg/hr]": stream.F_mass,
            "Revenue [MM$/yr]": mm_yr,
        })
    df = pd.DataFrame(rows).sort_values("Revenue [MM$/yr]", ascending=False)
    total = df["Revenue [MM$/yr]"].sum()
    df["Share [%]"] = 100.0 * df["Revenue [MM$/yr]"] / total
    summary = pd.DataFrame([{
        "Product": "TOTAL",
        "Price [$/kg]": float("nan"),
        "Flow [kg/hr]": float("nan"),
        "Revenue [MM$/yr]": total,
        "Share [%]": 100.0,
    }])
    return pd.concat([df, summary], ignore_index=True), total


# ---------------------------------------------------------------------------
# Sankey: each product → Total revenue
# ---------------------------------------------------------------------------
def plot_revenue_sankey(df, total_mm, out_html, out_png=None):
    df_plot = df[df["Product"] != "TOTAL"].copy()
    labels = list(df_plot["Product"]) + ["Total revenue"]
    n = len(df_plot)
    total_idx = n

    # hover text with $ and %
    link_labels = [
        f"{row['Product']}: ${row['Revenue [MM$/yr]']:.2f} MM$/yr ({row['Share [%]']:.1f}%)"
        for _, row in df_plot.iterrows()
    ]

    fig = go.Figure(data=[go.Sankey(
        arrangement="snap",
        node=dict(
            label=labels,
            pad=22,
            thickness=20,
            color=COLORS[:n] + [TOTAL_COLOR],
            line=dict(color="white", width=0.5),
        ),
        link=dict(
            source=list(range(n)),
            target=[total_idx] * n,
            value=df_plot["Revenue [MM$/yr]"].tolist(),
            label=link_labels,
            color=[c.replace(")", ", 0.45)").replace("rgb", "rgba")
                   if c.startswith("rgb") else f"rgba(100,100,100,0.35)"
                   for c in COLORS[:n]],
        ),
    )])

    fig.update_layout(
        title=dict(
            text="<b>Figure 4 — Revenue split by product</b><br>"
                 "<sup>Revenue at SAF MPSP; coproducts at market</sup>",
            x=0.5,
            xanchor="center",
        ),
        font=dict(size=12),
        height=480,
        width=760,
        margin=dict(t=80, b=30, l=30, r=30),
    )

    fig.write_html(str(out_html))
    print(f"Saved Sankey (HTML): {out_html}")

    if out_png:
        try:
            fig.write_image(str(out_png), scale=2)
            print(f"Saved Sankey (PNG):  {out_png}")
        except Exception as e:
            print(f"PNG export skipped (install kaleido: pip install kaleido): {e}")

    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    out_dir = Path(__file__).resolve().parent
    out_html = out_dir / "figure4_revenue_sankey.html"
    out_png = out_dir / "figure4_revenue_sankey.png"
    out_xlsx = out_dir / "figure4_revenue_table.xlsx"

    solve_saf_mpsp()
    df, total = revenue_table(saf_tea)

    print("\n========== Revenue table (MPSP scenario) ==========")
    print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\ntea.sales check: ${saf_tea.sales/1e6:.3f} MM$/yr")

    df.to_excel(out_xlsx, index=False)
    print(f"Saved table: {out_xlsx}")

    plot_revenue_sankey(df, total, out_html, out_png)
    return df, fig


if __name__ == "__main__":
    df, fig = main()
    # fig.show()  # uncomment in Spyder/Jupyter