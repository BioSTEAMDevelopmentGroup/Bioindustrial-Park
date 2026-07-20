# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import os
import sys
from pathlib import Path

import biosteam as bst

os.environ["PATH"] = "/opt/homebrew/bin:" + os.environ["PATH"]

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.systems import create_ad_integrated_system

OUTDIR = SCRIPT_DIR.parent / "results" / "flowsheets"
OUTDIR.mkdir(parents=True, exist_ok=True)

create_chemicals()
if __name__ == "__main__":

    # Build at alpha=0.5 so both pathways are populated --> just to have it
    sys, streams, units, _ = create_ad_integrated_system(
        alpha=0.5,
        pretreatment_case="combined_PE",
    )
    sys.simulate()

    # Full integrated system diagram
    sys.diagram(file=str(OUTDIR / "integrated_biorefinery_full"))
    print(f"Saved: {OUTDIR / 'integrated_biorefinery_full'}")

    # Preprocessing only (PR → PC → ML)
    pre_units = [u for u in sys.units if u.ID in ("PR", "PC", "ML")]
    pre_sys = bst.System("preprocessing", path=pre_units)
    pre_sys.diagram(file=str(OUTDIR / "preprocessing"))
    print(f"Saved: {OUTDIR / 'preprocessing'}")

    # Methanogenic pathway only (PX, EZ, AD → H2SR → UP → SP_AD)
    meth_ids = ("PX", "EZ", "AD", "H2SR", "UP", "SP_AD")
    meth_units = [u for u in sys.units if u.ID in meth_ids]
    if meth_units:
        meth_sys = bst.System("methanogenic_path", path=meth_units)
        meth_sys.diagram(file=str(OUTDIR / "methanogenic_pathway"))
        print(f"Saved: {OUTDIR / 'methanogenic_pathway'}")

    # VFA fermentation pathway only
    exclude_ids = {
        "PR", "PC", "ML", "SPL",      # shared front end
        "PX", "EZ", "HT", "MX",       # methane-side pretreatment
        "AD", "H2SR", "UP", "SP_AD",  # methane pathway
    }

    vfa_units = [u for u in sys.units if u.ID not in exclude_ids]

    if vfa_units:
        vfa_sys = bst.System("vfa_ferm_path", path=vfa_units)
        vfa_sys.diagram(file=str(OUTDIR / "vfa_fermentation_pathway"))
        print(f"Saved: {OUTDIR / 'vfa_fermentation_pathway'}")

    from biorefineries.sabre.systems import create_ad_vfa_system

    sys = create_ad_vfa_system(feedstock="pelagic")
    sys.simulate()
    sys.diagram(file=str(OUTDIR / "vfa_subset"))

    print(f"\nAll flowsheets saved to: {OUTDIR.resolve()}")