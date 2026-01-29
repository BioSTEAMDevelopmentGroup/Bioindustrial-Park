
import biosteam as bst
from biorefineries.prefers.v1.LegHb.system import create_LegHb_system, optimize_NH3_loading
from biorefineries.prefers.v1._process_settings import load_process_settings
import numpy as np
import traceback

print("1. Loading Process Settings...")
load_process_settings()

print(f"PowerUtility GWP (from settings): {bst.PowerUtility.characterization_factors.get('GWP')}")

# Manual Sanitization
raw_ = bst.PowerUtility.characterization_factors.get('GWP')
if isinstance(raw_, (tuple, list)):
    bst.PowerUtility.characterization_factors['GWP'] = raw_[0]
    print(f"Sanitized PowerUtility GWP to: {bst.PowerUtility.characterization_factors['GWP']}")

print("2. Creating System...")
sys = create_LegHb_system()
# Optimization needed to get flows?
print("3. Optimizing NH3...")
optimize_NH3_loading(sys, verbose=False)
sys.simulate()

print("4. Checking Stream GWP factors...")
for s in sys.feeds:
    cf = s.characterization_factors.get('GWP')
    print(f"Stream {s.ID}: {cf} (Type: {type(cf)})")

print("5. Running LCA Table...")
try:
    lca_table = bst.report.lca_displacement_allocation_table(
        systems=[sys],
        key='GWP',
        items=[sys.products[0]],
    )
    print("LCA Table Success!")
    print(lca_table)
except Exception:
    traceback.print_exc()

print("6. Manual LCA Calculation check...")
# Replicate logic roughly
key = 'GWP'
total = 0
for s in sys.feeds:
    if s.price: # Only economic feeds usually count? No, LCA counts all with CF
        cf = s.characterization_factors.get(key, 0)
        impact = s.F_mass * cf
        if cf:
            print(f"  {s.ID}: Mass={s.F_mass:.2f} * CF={cf} = {impact}")
            total += impact
        # Check if any operation fails
        try:
            val = s.F_mass * cf
        except Exception as e:
            print(f"ERROR calculating impact for {s.ID}: {e}")

elec = sys.power_utility.rate # kW
elec_cf = bst.PowerUtility.characterization_factors.get(key, 0)
print(f"  Electricity: Rate={elec:.2f} * CF={elec_cf} = {elec * elec_cf}")
try:
    val = elec * elec_cf
except Exception as e:
    print(f"ERROR calculating electricity impact: {e}")

