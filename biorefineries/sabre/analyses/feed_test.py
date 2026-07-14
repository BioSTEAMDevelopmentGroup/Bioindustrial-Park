# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Feed testing

1) Making sure the feed stream is created correctly with the expected mass flow and composition based on the YAML assumptions for "pelagic_high_quality" Sargassum.
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorefineries.sabre._chemicals import set_thermo
from biorefineries.sabre._process_settings import load_assumptions, get_quality_params, get_scale_feed_kgph
from biorefineries.sabre.streams import make_sargassum_feed

set_thermo()
A = load_assumptions()
q = get_quality_params(A, "pelagic_high_quality")
fresh_feed_kgph = get_scale_feed_kgph(A)

feed = make_sargassum_feed(
    fresh_feed_kgph=fresh_feed_kgph,
    moisture_frac=q["moisture_frac"],
    quality="pelagic_high_quality",
)

print(f"\nFeed total: {feed.F_mass:.2f} kg/hr")
print(f"Moisture frac from YAML: {q['moisture_frac']:.4f}")
print()

water = float(feed.imass["Water"])
total = float(feed.F_mass)
dry = total - water

print(f"{'Component':<20} {'kg/hr':>12} {'wet wt%':>10} {'dry wt%':>10}")
print("-" * 55)

for cid in feed.chemicals.IDs:
    m = float(feed.imass[cid])
    if m > 0:
        wet_pct = m / total * 100
        dry_pct = m / dry * 100 if cid != "Water" else float("nan")
        dry_str = f"{dry_pct:>9.2f}%" if cid != "Water" else "       —"
        print(f"  {cid:<18} {m:>12.2f} {wet_pct:>9.2f}% {dry_str}")

print("-" * 55)
print(f"  {'TOTAL':<18} {total:>12.2f} {100.0:>9.1f}%")
print(f"\n  Dry mass:   {dry:.2f} kg/hr ({dry/total*100:.2f}% wet)")

ash = float(feed.imass["Water"]) if "Ash" not in feed.chemicals.IDs else float(feed.imass["Ash"])
ash = float(feed.imass["Ash"])
vs_dry = (dry - ash) / dry * 100
ts_dry = dry / total * 100

print(f"  Ash (dry):  {ash/dry*100:.2f}%")
print(f"  VS/TS:      {vs_dry/100 * dry / dry * 100:.2f}%  (i.e. VS = {dry-ash:.2f} kg/hr)")
print(f"  Moisture:   {water/total*100:.2f}%")
print(f"\n  YAML targets:")
print(f"    moisture_frac:    {q['moisture_frac']:.4f} ({q['moisture_frac']*100:.2f}%)")
print(f"    ash_wt_frac_dry:  {q.get('ash_wt_frac_dry', 'n/a')}")
print(f"    vs_ts:            {q.get('vs_ts', 'n/a')}")