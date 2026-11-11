# Input data provenance and preparation

Run `python prepare_inputs.py` before importing the model. The script copies 28
existing repository files, after verifying their SHA-256 checksums, into the
paths used by the storage model. It does not download data or overwrite the
shared source datasets. `INPUTS.json` records each exact source path, byte count,
and checksum. Prepared copies are ignored by Git.

## Existing repository inputs

- GIS boundaries, rainfed region, state transport layer, miscanthus/switchgrass
  feedstock tables, and tortuosity factors: `../GIS_data/`, the inputs already
  distributed with the parent SAF location study.
- Refinery coefficient arrays: existing files in `biorefineries/Location_multi-product/`.
  Most are in `Notebooks/bc_coef/`. The succinic-acid and potassium-sorbate
  input-coordinate arrays are byte-identical to the existing
  `Biorefinery_simulations/Acrylic_acid/outputs/X_data_AA.npy` and reuse that file;
  their product-specific output coefficients remain separate. See `INPUTS.json`.

All source files were compared by content with the author's working storage inputs
and found identical at parent commit `30c66715e5007519a7a5af68862a11695cb0ab53`.
Refer to the parent SAF and multiproduct study documentation for scientific
provenance and references. No new license for third-party datasets is asserted.

## Newly included input

`Truck_cost_per_state.csv` is the author-supplied transport-cost table from the
storage model working directory, copied without numerical changes. It has columns
`State` and `Cost_per_km` (USD/km as used by the transport model). A source
publication and cost-year citation were not supplied with this release.

## Data structure

Keep shapefile sidecars together. Feedstock CSV headers are
`siteid, lat, long, y_data, Net_CI, Farm_price`; see the source models and prior
study documentation for units and transformations. Each refinery coefficient
array has shape `(1000000, 2)` and dtype `float64`. The code uses an 80/20 split
with `random_state=42` and linear interpolation; outside-support fallback uses
minima of test predictions. All five products are loaded at import even though
the default coupled model uses ethanol.

The prepared inputs total approximately 214 MB. This release reuses existing
repository datasets instead of adding duplicate binary files. Manuscript-specific
results and plotting notebooks are not included in the example workflow.
