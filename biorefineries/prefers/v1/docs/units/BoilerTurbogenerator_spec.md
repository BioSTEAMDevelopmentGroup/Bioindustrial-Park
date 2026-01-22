# BoilerTurbogenerator (Prefers-local)

**Source:** `v1/_units.py`  
**Base Class:** `bst.BoilerTurbogenerator`

## Purpose

Provides steam and electricity from cell-mass and residue streams while safely handling incomplete thermodynamic data for large biomolecules (e.g., proteins).

## Prefers-local Enhancements

1. **Guarded emissions enthalpy update**
   - Wraps the emissions enthalpy calculation and caches the enthalpy directly when property evaluation fails.
   - Prevents failures from missing or non-physical vaporization properties.

2. **Thermo fallbacks for large biomolecules**
   - Glucose model copy for robust fallback properties (e.g., Tb, Hvap, Psat, Cp, V).
   - Protein-derived heating values for globin/leghemoglobin.
   - Nonvolatile solid and vaporization fallbacks for boiler feed components.

## Usage

In process configs, use the prefers-local class:

```
from biorefineries.prefers.v1 import _units as u
BT = u.BoilerTurbogenerator(...)
```

## Notes

- This class is required for robust BT simulation in the LegHb configs where protein-rich debris is routed to the boiler.
- The base BioSTEAM implementation remains unchanged.
