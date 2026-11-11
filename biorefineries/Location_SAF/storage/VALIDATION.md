# Validation of the publication copy

Validated with the author's Storage3 environment (Python 3.9.21).

- All Python sources parse successfully.
- Default coupled simulation: N=2, seed=42, four configurations, 8 rows and 65 columns.
- Two runs with identical seeds produced identical DataFrames.
- Import and simulation succeeded from a different working directory.
- All numeric outputs finite; service level in 0-100 percent; nonnegative storage and transport costs.
- Example adjusted ethanol MSP: approximately 2.92-2.99 USD/gal. These two draws are a smoke check, not a study estimate.
- Model equations and parameters are preserved; input paths and stale summary/plot column references were corrected, and trailing whitespace was removed in the publication copy.

This is not validation of manuscript results or a full 1,000-draw reproduction.

- Summary and diagnostic plotting completed after correcting stale column names.
- All 28 reused repository input files pass checksum verification and preparation.
