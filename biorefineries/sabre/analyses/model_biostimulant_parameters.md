# Biostimulant model parameters

Full parameter list for `model_biostimulant.create_biostimulant_model()` --
26 parameters, all uniform in this model (the two `tea.yaml` entries with
`triangular` distributions -- `vfa`, `biomethane_mmbtu` -- aren't part of
the biostimulant system).

| Element | Parameter | Units | Baseline | Distribution | Bounds |
|---|---|---|---|---|---|
| Press-PR | Press solids capture fraction | | 0.98 | Uniform | [0.49, 1.0] |
| Press-PR | Press cake solids wt fraction | | 0.15 | Uniform | [0.075, 0.225] |
| Press-PR | Press solubles to pressate fraction | | 1.0 | Uniform | [0.5, 1.0] |
| Press-PR | Press power intensity | kWh/dry ton | 5.0 | Uniform | [2.5, 7.5] |
| Pressate concentrator-PC | PC water recovery to permeate | | 0.93 | Uniform | [0.465, 1.0] |
| Pressate concentrator-PC | PC retained solute recovery to concentrate | | 0.95 | Uniform | [0.475, 1.0] |
| Pressate concentrator-PC | PC nontarget solute recovery to permeate | | 0.70 | Uniform | [0.35, 1.0] |
| Pressate concentrator-PC | PC design flux | L/m2-h | 35.0 | Uniform | [17.5, 52.5] |
| Pressate concentrator-PC | PC electricity intensity | kWh/m3 | 2.5 | Uniform | [1.25, 3.75] |
| Pressate concentrator-PC | PC membrane capex per m2 | USD/m2 | 500.0 | Uniform | [250, 750] |
| Pressate concentrator-PC | PC maintenance fraction of capex | | 0.035 | Uniform | [0.0175, 0.0525] |
| TEA | IRR | | 0.10 | Uniform | [0.05, 0.15] |
| TEA | Income tax | | 0.21 | Uniform | [0.105, 0.315] |
| TEA | WC over FCI | | 0.05 | Uniform | [0.025, 0.075] |
| TEA | Finance interest | | 0.08 | Uniform | [0.04, 0.12] |
| TEA | Finance fraction | | 0.60 | Uniform | [0.3, 0.9] |
| TEA | Startup FOC fraction | | 1.0 | Uniform | [0.5, 1.0] |
| TEA | Startup VOC fraction | | 0.50 | Uniform | [0.25, 0.75] |
| TEA | Startup sales fraction | | 0.50 | Uniform | [0.25, 0.75] |
| TEA | FOC fraction of FCI | | 0.04 | Uniform | [0.02, 0.06] |
| TEA | Finance years | | 10 | Uniform | [5, 15] |
| TEA | Startup months | | 3 | Uniform | [1.5, 4.5] |
| TEA | Operating days | | 330 | Uniform | [165, 365] (capped, not 495) |
| Stream-Sargassum feed | Sargassum feed price | USD/kg | 0.0 | Uniform | [-0.02, 0.1] (yaml's own range) |
| Stream-Pressed cake | Pressed cake disposal price | USD/kg | -0.04 | Uniform | [-0.08, 0.0] (yaml's own range) |
| Stream-Permeate | Permeate disposal price | USD/kg | -0.004 | Uniform | [-0.008, 0.0] (yaml's own range) |

All rows except the last three are baseline-only +/-50% (fraction-type rows
clipped to [0,1] where the upper bound would otherwise exceed 1); the last
three use the explicit `range`/`uniform` already written in `data/tea.yaml`.

See `model_utils.py` (`distribution_from_yaml`, `add_tea_parameters`) and
`model_biostimulant.py` (`add_press_parameters`,
`add_pressate_concentrator_parameters`, `create_biostimulant_model`) for the
code that builds these, and
`docs/superpowers/specs/2026-07-31-sabre-uncertainty-model-design.md` for
the design rationale (distribution rule, exclusions).
