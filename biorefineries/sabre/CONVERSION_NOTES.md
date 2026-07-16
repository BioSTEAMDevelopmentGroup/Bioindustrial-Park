# SaBRe conversion notes

Bugs found while converting `biorefineries/sabre` from its pre-conversion form
(`BioSTEAM-platform/sabre_ad`) that **caused reported metric values (TCI, VOC,
MSP, NPV, etc.) to deviate from the `sabre_ad` baseline**. Updated as the
conversion proceeds subsystem by subsystem.

**`assumptions.yaml` is off-limits for edits until the user has reviewed its
citations themselves.** It was reverted to its exact original (pre-conversion)
content — the `vfa_fermentation:` section that had been deleted earlier in this
conversion (on the reasoning that its values were now baked into code) is back.
Restoring it is inert from the code's perspective: no Python file reads
`A["vfa_fermentation"]` anymore, those call sites were already removed from
`analyses/vfa_fermentation_tea.py` and `systems/_integrated_system.py`
independent of whether the yaml section exists. This reversal was requested by
the user directly because of the citation-transfer failures below — do not
re-delete or otherwise edit `assumptions.yaml` without explicit approval.

Scope note: several stale unit-level defaults were found and corrected during the
VFA fermentation conversion (e.g. `YarrowiaLipidFermenter`'s own
`product_yield_kg_per_kg_vfa_consumed` default was `0.45` instead of `0.144`).
None of those changed any reported number, because every real code path always
passed the correct value explicitly, overriding the wrong default — so the stale
defaults were dead code, never actually exercised. Those don't belong here; this
file is only for cases where the actual output changed or would change.

Two apparent discrepancies were found and deliberately **not** changed, to keep
output matching the `sabre_ad` baseline exactly:
`target_oil_and_solids_content` (60 vs 70) and `ferm_mgso4_dose` (0.0 vs 0.49,
integrated pathway only) — see git history / conversation log if these need
revisiting later. Neither is logged below since neither changed the baseline.

---

## Metric-deviation bugs

*(none yet — VFA fermentation subsystem converted with zero baseline deviation,
confirmed by `tests.py` matching `sabre_ad` output exactly)*

---

## Citation-accuracy issues (no metric impact — flagging separately)

These don't belong in the section above (they never changed any reported number),
but are real documentation problems worth fixing before anyone treats these
numbers/citations as trustworthy.

**General warning, not just about the three items below**: the AI assistant
doing this conversion is unreliable at transferring citations and has produced
wrong ones repeatedly in this subsystem alone — not just miscopying an existing
citation onto the wrong number, but in one case inventing an attribution with
no basis anywhere in either codebase, and in another carrying a citation over
from a completely different unit. All three were only caught because the user
explicitly pushed back and asked for a self-audit; none were self-reported
proactively. **Do not trust any citation added to a docstring/comment by this
conversion process without independently re-checking it against
`assumptions.yaml` (and, ideally, the actual primary source) yourself.** This
applies to every subsystem still to be converted (AD-biogas, preprocessing,
VFA-AD, pretreatment units), not only the ones listed below.

### `OilExtractionPlaceholder` electricity citation — NEEDS VERIFICATION

`units/_oil_extraction.py`'s docstring currently reads:

> Electricity: high-pressure homogenization 0.203 kWh/kg dry biomass,
> harmonized to the NREL algal-lipid basis (Doucha & Livansky 2008; Postma et
> al. 2017)

This citation attribution is **not verified** and was introduced during the
conversion by carelessly merging two unrelated pieces of pre-existing text:

1. `assumptions.yaml`'s `vfa_fermentation.downstream_recovery` notes say only
   `"Homogenization electricity aligned to the harmonized NREL algal-lipid basis
   rather than the earlier 1.5 kWh/kg placeholder."` — no specific paper named,
   no `sources:` block for this value.
2. The pre-conversion docstring (`sabre_ad/sabre/units/oil_extraction.py`,
   `biorefineries/sabre` initial commit `13e2aa6`) cited `(Doucha & Livansky
   2008; Postma et al. 2017)` — but attached to the **old, superseded** value of
   `~1.5 kWh/kg`, not the new `0.203`.

When the unit's default was updated from `1.5` to `0.203` (to match what
`create_vfa_fermentation_system()` was already actually using — see below), the
old citation names were carried over onto the new number without checking
whether those two papers actually report `0.203 kWh/kg`. They may not.

**Also unverified**: the capital cost anchor's "back-calculated from NREL 55431"
derivation. `assumptions.yaml` states the result ($7.848M) but doesn't show the
calculation. A separate pre-existing docstring in this codebase claimed the
report itself gives ~$8.5M for the same throughput — nobody involved in this
conversion has read NREL/TP-5100-55431 to check that claim either. Don't repeat
"~$8.5M" as if it's a confirmed figure from the report; it's just what an
earlier, unverified docstring asserted.

**STATUS: OPEN — needs primary-source verification** (find and check Doucha &
Livansky 2008, Postma et al. 2017, and NREL/TP-5100-55431 directly, or find
whoever wrote the `assumptions.yaml` notes and ask what they actually derived
from). Flagged in the unit's docstring; see `units/_oil_extraction.py`.

**Why the underlying number swap (1.5→0.203, $8.5M→$7.848M) itself is not
listed as a metric-deviation bug**: `OilExtractionPlaceholder`'s own `__init__`
default was always overridden by `create_vfa_fermentation_system()`'s own
matching default (`0.203`, `$7,848,000`) at every real call site, in both the
pre-conversion and post-conversion code — see
`systems/_vfa_fermentation_system.py`, where `oil_extraction_homogenization_kWh_per_kg`
and `oil_extraction_ref_installed_cost_usd` are declared as that function's own
parameters and explicitly threaded into the `OilExtractionPlaceholder(...)`
call. So the unit's own default was dead text; updating it to match reality
caused zero simulation output change. What's unverified is whether "reality"
(0.203, $7.848M) itself is correct.

### `YarrowiaLipidFermenter.conversion` — fabricated citation, now fixed

While writing the docstring "Sources" section, `conversion=0.85` was grouped
with `target_pH=8.0` under a shared "Gao et al. 2020" citation. Checked
directly (grep across both `sabre_ad` and `biorefineries/sabre`, all `.py` and
`.yaml`): there is **no comment, yaml note, or `sources:` entry anywhere that
ties `conversion=0.85` to Gao et al. 2020 or any other source** — only
`target_pH` (via `assumptions.yaml`'s `pH_control` source note) and
`product_yield_kg_per_kg_vfa_consumed` (via an inline comment in
`vfa_fermentation_tea.py`, `"Gao et al. 2020, Y_L/S at pH 8"`) are actually
grounded. The `conversion` attribution was invented outright, not
mis-transcribed from existing text. **FIXED**: docstring and inline comment in
`units/_vfa_fermentation.py` now correctly show `conversion` as having no
identified literature source (grouped with the other unsourced scoping
assumptions instead).

### `VFAMicrofilter` design flux — cross-contaminated citation from a different unit, now fixed

The docstring attributed `design_flux_L_m2_h=20.0` to "seaweed/algal membrane
concentration literature." That citation (Sievers et al. 2017; Díaz-Reinoso
and Domínguez 2020) actually belongs to a **different unit** —
`assumptions.yaml`'s `pressate_biostimulant.concentrator.sources.flux`, which
documents `PressateConcentrator`'s flux value of **35 LMH**, not
`VFAMicrofilter`'s 20 LMH. `vfa_microfilter`'s own yaml notes give no citation
at all, just descriptive text ("20 LMH selected as a conservative design flux
for low-pressure clarification"). **FIXED**: docstring in
`units/_vfa_fermentation.py` now states no source is given for this value and
explicitly flags the earlier mix-up so it doesn't recur when
`PressateConcentrator` gets converted later.
