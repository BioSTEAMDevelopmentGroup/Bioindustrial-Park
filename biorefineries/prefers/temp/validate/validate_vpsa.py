# Simple VacuumPSA unit validation script
import sys
import biosteam as bst
import sys
import biosteam as bst
import thermosteam as tmo
import functools
print = functools.partial(print, flush=True)

# 1. Import the unit class FIRST (before setting thermo)
# This prevents submodule imports (like LegH) from overwriting our thermo settings later
try:
    from biorefineries.prefers.v1._units import VacuumPSA
except ImportError:
    # If package structure is an issue (e.g. running from root), try appending path
    sys.path.append(r'c:\Programming\PreFerS\Bioindustrial-Park')
    from biorefineries.prefers.v1._units import VacuumPSA

print("Step 1: Imports successful")

# 2. Define chemical property package
chemicals = tmo.Chemicals(['H2', 'CO', 'C2H4', 'CO2', 'N2'])
bst.settings.set_thermo(chemicals)
print("Step 2: Thermo set to: H2, CO, C2H4, CO2, N2")

# 3. Create test stream
feed = bst.Stream('feed', H2=50, CO=30, C2H4=15, CO2=3, N2=2, phase='g')
print(f"Step 3: Feed created: {feed.F_mol:.2f} kmol/hr")

# 4. Instantiate Unit
# Explicitly passing thermo ensures compatibility, though set_thermo should handle it
print("Step 4: Instantiating VacuumPSA...")
V = VacuumPSA('V1', ins=feed, thermo=chemicals)
print("  OK")

# 5. Simulate
print("Step 5: Simulating...")
V.simulate()
print("  OK")

print("\n=== Results ===")
print(f"Feed:    {V.ins[0].F_mol:.2f} kmol/hr")
print(f"Product: {V.outs[0].F_mol:.2f} kmol/hr")
print(f"Purge:   {V.outs[1].F_mol:.2f} kmol/hr")

print("\n=== Design Results ===")
for k, v in V.design_results.items():
    print(f"  {k}: {v:.4g}")

print("\n=== Purchase Costs ===")
for k, v in V.purchase_costs.items():
    print(f"  {k}: ${v:,.0f}")

print(f"\nInstalled cost: ${V.installed_cost:,.0f}")
print("\n✓ Validation complete!")
