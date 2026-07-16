# SaBRe conversion notes

Bugs found while converting `biorefineries/sabre` from its pre-conversion form
(`BioSTEAM-platform/sabre_ad`) that **caused reported metric values (TCI, VOC,
MSP, NPV, etc.) to deviate from the `sabre_ad` baseline**. Updated as the
conversion proceeds subsystem by subsystem.

**Editing `assumptions.yaml` is out of scope for this conversion, permanently
— not a temporary pause.** The original plan called for slimming yaml down
subsystem by subsystem; that's no longer the plan. The `vfa_fermentation:`
section was deleted (125 lines) during this conversion, then the citation
errors below were found, and the user reverted it to its exact original
content and removed yaml editing from scope entirely, for every subsystem,
not just this one. Restoring it was inert from the code's perspective — no
Python file reads `A["vfa_fermentation"]` anymore regardless of whether the
yaml section exists. Going forward: bake assumptions into unit/system-builder
code as usual; leave `assumptions.yaml` itself completely untouched. Yaml
curation is the user's own task, done separately.

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

*(none yet — VFA fermentation, AD-biogas, VFA-AD, and digestate handling
subsystems converted with zero baseline deviation, confirmed by `tests.py`'s
full 20-test suite matching pre-conversion output exactly, `rtol=1e-6`)*

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

### AD-biogas subsystem (`units/_ad.py`, `units/_h2s_removal.py`, `units/_biogas_upgrading.py`)

Converted with zero metric deviation (`tests.py`'s AD-related tests —
`test_ad_tea_final_*`, `test_ad_biostimulant_price_*`, `test_ad_feed_tea_*`,
`test_ad_heatmap_*`, `test_export_ad_report_fixed_*`,
`test_hauke_stream_tables_*`, `test_plot_two_pretreatment_figures_*`,
`test_integrated_tea_*`, `test_integrated_system_sensitivity_*` — all
unchanged at `rtol=1e-6`). Scope was limited to these three unit files plus
the `H2SR` construction call in `systems/_ad_biogas_system.py`; the rest of
that system-builder file (Press, Mill, pretreatment units, digestate units,
pressate concentration) belongs to other not-yet-converted subsystems and
was left untouched, mirroring how `_integrated_system.py` was scoped to
"fermentation block only" during the VFA conversion.

**`AnaerobicDigester.ADBC_VOL_M3`/`ADBC_CAPEX` — unverified provenance, not
fixed, flagged only.** The six-point interpolation table actually used by
`interp_capex()` (the real installed-cost calculation) has no corresponding
entry anywhere in `assumptions.yaml`. A *separate*, structurally unused
single-point anchor (`ad_costing.base_volume_m3`/`base_capex_usd`, stored on
the unit's own `base_volume_m3`/`base_capex_usd` `__init__` params but never
read by `_design`/`_cost`) is cited in yaml to the "ADBC" spreadsheet tool
(`ADBCv2-48 (2).xlsm`, sheet "AnaerobicDigester"). Whether the six-point
table was pulled from the same spreadsheet run as that single-point anchor
has not been verified by anyone in this conversion. **Do not** state or
imply that `ADBC_VOL_M3`/`ADBC_CAPEX` come from that spreadsheet without
opening the actual file and checking. Flagged in the unit's docstring;
status: OPEN, same as the `OilExtractionPlaceholder` NREL/TP-5100-55431 item
above.

**`H2SRemoval.h2s_removal_efficiency` (0.99) is a deliberately conservative
modeling choice, not a transcription of its own citation.** yaml's
`h2s_removal.sources.efficiency` note (Choudhury et al. 2019) states the
cited iron-sponge technology achieves >99.9% removal in practice, but the
actual modeled parameter (`h2s_removal.h2s_removal_efficiency` in yaml, and
this unit's own default) is `0.99`, not `0.999`. The unit's docstring now
states this distinction explicitly rather than implying the citation
supports the exact modeled number.

**`BiogasUpgrading` intentionally kept custom `_design`/`_cost` instead of
switching to the `@cost` decorator**, unlike `H2SRemoval`. Its formula is
single-train and linear in flow (`n=1`), which would otherwise be a clean
`@cost` candidate — but its annual maintenance OPEX is
`maintenance_frac_of_capex_per_yr * <this unit's own installed cost>`, and
`@cost`'s decorator-owned `_cost()` only computes that installed cost
*after* `_design()` has already run (BioSTEAM lifecycle: `_run` → `_design`
→ `_cost`). There's no way to read the not-yet-computed CAPEX from
`_design()` without duplicating the cost formula there too — a latent-bug
risk if the two copies of the formula ever drift. `H2SRemoval` didn't have
this problem because its only OPEX (reagent replacement) depends on
`_design()`-computed raw biogas flow, not on its own CAPEX. This is a
deliberate, documented exception to the general "single-train power-law
units use `@cost`" convention, not an inconsistency to "fix" later.

**`H2SRemoval` zero-flow edge case changed behavior, non-triggering.** The
pre-conversion custom `_cost()` had a `Q_Nm3ph <= 0` guard that returned the
full reference installed cost (`$450,000`) even at zero raw-biogas flow. The
`@cost` decorator's formula naturally gives `0**0.7 = 0` instead. This never
triggers in any of the 20 regression tests or any real call site (the AD
unit always produces nonzero biogas in every case exercised in this
codebase), so it wasn't preserved — flagging it here in case a future
zero-throughput sensitivity run surfaces it.

**Several stale unit-level `__init__` defaults were corrected** (dead code —
always overridden by `systems/_ad_biogas_system.py`, zero output impact,
verified by the unchanged test suite): `AnaerobicDigester.vs_destruction`
(0.50→0.20), `ch4_kg_per_kg_vs_fed` (0.0555→0.100, the `press_mill_only`
baseline case value — not universal, see unit docstring),
`raw_biogas_molfrac` (CarbonDioxide/HydrogenSulfide 0.43/0.02→0.42/0.03),
`headspace_frac` (0.20→0.25), `maintenance_usd_per_m3_yr` (None→10.0),
`digestible_IDs`/`biodegradability` (dropped `Xylan`/`Mannan`/`Galactan`,
absent from yaml's list — these are real thermo chemicals, not typos, just
not part of what `assumptions.yaml`'s AD sections track as biodegradable);
`H2SRemoval.reagent_cost_usd_per_Nm3_raw` (0.005→0.002);
`BiogasUpgrading.ch4_recovery` (0.98→0.99). Per the scope note at the top of
this file, none of these belong in the metric-deviation section above.

### VFA-AD pathway (`units/_ad_vfa.py`, `AcidogenicDigester`)

Converted with zero metric deviation (full `tests.py` suite unchanged at
`rtol=1e-6`, including the VFA-pathway-dependent tests --
`test_vfa_pathway_figures_vfa_composition`,
`test_vfa_price_scenarios_biostimulant_price`,
`test_vfa_yield_sensitivities_*`, `test_hauke_stream_tables_*` -- since all
of these route through `create_vfa_ad_system()` → `AcidogenicDigester`).

`AcidogenicDigester` is structurally already correct (multi-train, custom
`_design`/`_cost`, reuses `AnaerobicDigester`'s `interp_capex`/`ADBC_VOL_M3`/
`ADBC_CAPEX` table via `from biorefineries.sabre.units._ad import
GAL_TO_M3, interp_capex`) -- no structural change, only a Sources docstring
and stale-default corrections, same treatment as the AD-biogas subsystem.
The same **unverified ADBC_VOL_M3/ADBC_CAPEX provenance** caveat logged
above for `AnaerobicDigester` applies here too (shared table) and is not
re-verified separately.

No edit was needed to `systems/_vfa_ad_system.py`: unlike `H2SRemoval` in
the AD-biogas subsystem, `AcidogenicDigester`'s `__init__` signature did not
change (no `@cost`-decorator conversion here — it's multi-train, same
reasoning as `AnaerobicDigester`), so the existing construction call in
`create_vfa_ad_system()` needed no changes. Confirmed via grep that
`AcidogenicDigester` is constructed in exactly one place in this codebase.

**Stale `__init__` defaults corrected** (dead code — always overridden by
`create_vfa_ad_system()`, zero output impact): `vfa_kg_per_kg_vs`
(0.80→0.55, yaml `vfa_ad_performance.cases.seaweed_arrested_fitted.
vfa_kg_per_kg_vs`), `headspace_frac` (0.2→0.25, yaml `vfa_ad.
gas_storage_frac_of_total_volume`), `digestible_IDs` (dropped `Xylan`,
`Mannan`, `Galactan` — real thermo chemicals absent from yaml's list, same
pattern as `AnaerobicDigester` — and `Cellulose`, which is not even a
chemical defined in `_chemicals.py` at all, a stale leftover reference with
no yaml counterpart to check against). `vs_destruction` (0.50) already
matched yaml exactly, no change needed.

**`_resolve_vfa_split()`'s internal fallback split left untouched,
deliberately** (AceticAcid 0.40/PropionicAcid 0.10/ButyricAcid
0.30/ValericAcid 0.05/HexanoicAcid 0.15) — while dead in the current
codebase (the one real call site always supplies yaml's fitted split,
AceticAcid 0.648/PropionicAcid 0.186/ButyricAcid 0.084/ValericAcid
0.082/HexanoicAcid 0.000), this fallback is a generic placeholder for any
thermo with the five common VFA IDs, not a stale copy of the Sargassum
case's fitted split. Overwriting it with the case-specific numbers would
conflate a general-purpose fallback with a case-specific fitted result —
a different situation from the other stale-default corrections in this
file, so not "corrected" the same way.

**`_get_vfa_case()`/`_get_vfa_split()` shallow-merge split, reviewed and
confirmed non-fragile as currently used.** `_get_vfa_case()` does a
one-level `dict.update()` of the selected case onto the top-level
`vfa_ad_performance` dict; `_get_vfa_split()` is kept as a separate function
specifically to avoid relying on that shallow merge for the nested
`vfa_split` dict. As currently used, this isn't an active bug — top-level
`vfa_ad_performance` has no `vfa_split` key of its own for the case's key to
shallow-clobber, so `.update()` behaves correctly either way. The separate
accessor is defensive, not a workaround for an active bug. Not changed;
already the correct level of caution given the current yaml shape.

### Digestate handling (`units/_screwpress.py`, `units/_centrifuge.py`)

Scope note: this subsystem is `DigestateScrewPress` and
`DigestateDecanterCentrifuge` only. `PressateConcentrator` and
`BiostimulantEvaporator` (also under `units/`, also digestate-adjacent by
proximity in `_ad_biogas_system.py`) are the *pressate-biostimulant*
side-stream, a separate subsystem, not touched here.

Converted with zero metric deviation (full `tests.py` suite unchanged at
`rtol=1e-6` — `DigestateScrewPress` is exercised by every AD-biogas and
VFA-AD test; `DigestateDecanterCentrifuge` is not constructed by any system
builder in this codebase at all, confirmed by grep, so it carries zero test
coverage either way).

**`DigestateScrewPress.solids_IDs` is a dead `__init__` parameter — accepted
and stored on `self`, but never read by `_run()`, `_design()`, or
`_cost()`.** `_run()` instead computes TS dynamically as "everything except
Water and `dissolved_IDs`." This is a different, deeper kind of staleness
than the usual "always-overridden by the system builder" pattern seen
elsewhere in this conversion — the parameter has no effect at any value, not
just its default. Left unchanged (default still references `"Cellulose"`,
not a chemical defined anywhere in `_chemicals.py`) rather than "corrected"
to `assumptions.yaml`'s `digestate_screw_press.solids_IDs` list, since doing
so would misleadingly imply the parameter does something it doesn't.
Flagging here in case a future pass wants to either wire it in or remove it
outright — out of scope for a pure conversion pass either way.

**`DigestateScrewPress` stale defaults corrected** (dead code — always
overridden by both real call sites, `systems/_ad_biogas_system.py` and
`systems/_vfa_ad_system.py` — zero output impact): `ts_capture_frac`
(0.33→0.40), `cake_moisture_frac` (0.77→0.45), `eur_to_usd` (1.0→1.19), all
against `assumptions.yaml`'s `digestate_screw_press` section. `capacity_tph_each`
(6.0) and `kWh_per_m3` (0.67) already matched. `capex_eur_table` was *not*
stale — yaml has no `capex_eur_table` key at all, so the unit's own
hardcoded SYSTEMIC Table 2-11 default is genuinely exercised at every real
call site, unlike the other params in this class.

**`DigestateDecanterCentrifuge` — different situation from every other unit
converted so far: it's not constructed anywhere in this codebase**, so its
`__init__` defaults aren't "dead code always overridden" — they're the only
values that would ever apply if someone instantiated it. `solids_IDs`,
`ts_capture_frac`, `cake_moisture_frac`, and `capacity_tph_each` already
matched `assumptions.yaml`'s `digestate_decanter_centrifuge` section
exactly, no change needed. `centrifuge_purchase_cost_usd_each` was `None`
(silently giving $0 installed cost in `_cost()`) — corrected to yaml's
297,500.0, since leaving it `None` would make the unit silently free to
build if anyone does wire it in later. Neither yaml digestate section
(`digestate_screw_press` or `digestate_decanter_centrifuge`) has a
`sources:` sub-block, so no yaml-level literature citations exist for
either unit — all citations present are this codebase's own pre-existing
module-docstring text (SYSTEMIC D3.2 report 2021, SYSTEMIC database,
SYSTEMIC Table 2-11, Tchobanoglous et al. Wastewater Engineering 5th ed.),
carried into the class docstrings as-is, not independently re-verified
against the primary SYSTEMIC/Tchobanoglous sources in this conversion.
