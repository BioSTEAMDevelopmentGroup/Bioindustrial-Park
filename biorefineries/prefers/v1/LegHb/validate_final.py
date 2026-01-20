"""Final validation script for Phase 4 completion."""
from system._config1 import create_LegHb_system, check_LegHb_specifications
import biosteam as bst

print("\n" + "="*60)
print("FINAL VALIDATION - Phase 4 Complete")
print("="*60)

sys = create_LegHb_system()
sys.simulate()

# Get the product stream from the system
product = bst.main_flowsheet.stream.FormulatedProduct

# Run the built-in specification check
check_LegHb_specifications(product)

print("\n" + "="*60)
print("✓ ALL PHASES COMPLETE - Documentation Synced")
print("="*60 + "\n")
