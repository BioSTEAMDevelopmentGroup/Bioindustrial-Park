#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HEFA Mass Balance Sankey
"""

from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey

import biosteam as bst
from biorefineries.covercress.process_settings import price

# ---------------------------------------------------------------------------
# Load covercress system (skip auto run_tea at bottom of system file)
# ---------------------------------------------------------------------------
BIOREFINERIES = Path(__file__).resolve().parents[1]

if str(BIOREFINERIES) not in sys.path:
    sys.path.insert(0, str(BIOREFINERIES))

SYSTEM_FILE = BIOREFINERIES / "covercress" / "sys_protein_ext_with_combined_ext_centri.py"

with open(SYSTEM_FILE) as f:
    code = f.read()

# Disable TEA auto-run
code = code.replace("MSP_SAF, MFPP = run_tea()", "pass")
code = code.replace("MSP_SAF = run_tea()", "pass")

_ns = {"__name__": "__covercress__", "__file__": str(SYSTEM_FILE)}
exec(compile(code, str(SYSTEM_FILE), "exec"), _ns)

# ---------------------------------------------------------------------------
# 1. RUN SYSTEM
# ---------------------------------------------------------------------------
sys = _ns['oil_extraction_sys']
COLORS = [
    "#B97A57", "#7BBD84", "#F7C652", "#63C6CE", "#94948C",
    "#734A8C", "#D1C0E1", "#648496", "#9C7FB8", "#F8858A",
    "#C94C68", "#5B8FA8",
]

sys.simulate()

# ---------------------------------------------------------------------------
# 2. EXTRACT MASS FLOWS
# ---------------------------------------------------------------------------
SAF = sys.outs[0].F_mass
protein = sys.outs[1].F_mass
naphtha = sys.outs[2].F_mass
green_diesel = sys.outs[3].F_mass
propane = sys.outs[4].F_mass

# Hydrogen input
H2_in = sys.flowsheet.stream.hydrogen_fresh.F_mass

# Oil input (sum of TAG species)
feed = sys.flowsheet.stream.flake_feed
oil_in = (
    feed.imass['OOO'] +
    feed.imass['LLL'] +
    feed.imass['LnLnLn'] +
    feed.imass['SSS']
)

# CO2 from HDO
CO2 = sys.flowsheet.stream.hdo_gases.imass['CO2']

# Water from HDO + desolventizer
water = (
    sys.flowsheet.stream.hdo_gases.imass['Water'] +
    sys.flowsheet.stream.desolventizer_vapor.imass['Water']
)

# ---------------------------------------------------------------------------
# 3. NORMALIZE FLOWS
# ---------------------------------------------------------------------------
total_in = oil_in + H2_in
total_out = SAF + green_diesel + naphtha + propane + protein + CO2 + water

flows = [
    oil_in / total_in,
    H2_in / total_in,
    -SAF / total_out,
    -green_diesel / total_out,
    -naphtha / total_out,
    -propane / total_out,
    -protein / total_out,
    -CO2 / total_out,
    -water / total_out,
]

labels = [
    'Oil',
    'Hydrogen',
    'SAF',
    'Green Diesel',
    'Naphtha',
    'Propane',
    'Protein',
    'CO2',
    'Water'
]
# ---------------------------------------------------------------------------
# 4. CREATE SHAPED SANKEY DIAGRAM (MATCH FIRST FIGURE)
# ---------------------------------------------------------------------------

fig = plt.figure(figsize=(13, 6))
ax = plt.gca()

s = Sankey(ax=ax, unit=None)

# ---------------------------------------------------------
# SEGMENT 1 — Inputs → HEFA process (left side)
# ---------------------------------------------------------
d1 = s.add(
    flows=[
        flows[0],   # Oil
        flows[1],   # Hydrogen
        -(flows[2] + flows[3] + flows[4] + flows[5] + flows[6] + flows[7] + flows[8])
    ],
    labels=['Oil substrate', 'Hydrogen', 'Decarboxylation / Hydrocracking'],
    orientations=[1, 1, 0],     # Inputs flow right, process stays centered
    facecolor=COLORS[0],
    alpha=0.85,
    trunklength=1.4,            # Long horizontal trunk (matches first figure)
    pathlengths=[0.6, 0.6, 0.8]
)

# ---------------------------------------------------------
# SEGMENT 2 — HEFA process → Outputs (right side)
# ---------------------------------------------------------
d2 = s.add(
    flows=[
        (flows[2] + flows[3] + flows[4] + flows[5] + flows[6] + flows[7] + flows[8]),
        flows[5],   # LPG
        flows[4],   # Gasoline
        flows[2],   # Jet
        flows[3],   # Diesel
        flows[8],   # Water
        flows[7],   # CO2
        flows[6],   # Protein
    ],
    labels=['', 'LPG', 'Gasoline', 'Jet', 'Diesel', 'Water', 'CO2', 'Protein'],
    orientations=[0, -1, -1, -1, -1, -1, -1, -1],   # Vertical right-side outputs
    facecolor=COLORS[1],
    alpha=0.85,
    trunklength=1.4,
    pathlengths=[0.8] * 8,
    prior=0,
    connect=(2, 0)   # Connect process → process
)

# ---------------------------------------------------------
# DRAW
# ---------------------------------------------------------
s.finish()

# ---------------------------------------------------------
# APPLY CABBI COLORS — collect patches from Axes
# ---------------------------------------------------------
patches = [p for p in ax.patches]   # <-- FIXED: works on all Matplotlib versions

for i, p in enumerate(patches):
    p.set_facecolor(COLORS[i % len(COLORS)])
    p.set_edgecolor("black")
    p.set_linewidth(0.7)

plt.title("HEFA Mass Balance Sankey Diagram", fontsize=15)
plt.show()
# ---------------------------------------------------------------------------
# 3B. PRINT FRACTIONS (normalized mass balance)
# ---------------------------------------------------------------------------

print("\n=== HEFA Mass Balance Fractions ===")
print(f"Oil substrate fraction:     {flows[0]:.4f}")
print(f"Hydrogen fraction:          {flows[1]:.4f}")
print(f"SAF fraction:               {flows[2]:.4f}")
print(f"Green diesel fraction:      {flows[3]:.4f}")
print(f"Naphtha fraction:           {flows[4]:.4f}")
print(f"Propane (LPG) fraction:     {flows[5]:.4f}")
print(f"Protein fraction:           {flows[6]:.4f}")
print(f"CO2 fraction:               {flows[7]:.4f}")
print(f"Water fraction:             {flows[8]:.4f}")

print("\nTotal input fraction:       {:.4f}".format(flows[0] + flows[1]))
print("Total output fraction:      {:.4f}".format(
    flows[2] + flows[3] + flows[4] + flows[5] + flows[6] + flows[7] + flows[8]
))
print("====================================\n")
