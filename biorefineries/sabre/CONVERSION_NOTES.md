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

All 7 subsystems converted so far (VFA fermentation, AD-biogas, VFA-AD,
digestate handling, preprocessing, pressate-biostimulant, pretreatment
units) preserve `tests.py`'s full 20-test suite, `rtol=1e-6`, **except for
one deliberate fix** — see the `EnzymaticPretreatment`/
`PeroxidePretreatment` entry immediately below, whose corrected values are
now what `tests.py` asserts against. Two substantive discrepancies were
found in
the preprocessing subsystem and deliberately NOT fixed (see "Preprocessing"
under Citation-accuracy issues below; fixing either would change simulation
output). One functional fix was made in the pressate-biostimulant subsystem
(`BiostimulantEvaporator`'s missing `add_OPEX` wiring) but not logged here
since that unit carries zero test coverage and is never actually simulated,
so the fix changed no measured output — see "Pressate-biostimulant" below.

One **real, verified, currently-live bug** was found in the pretreatment
units subsystem. Initially left unfixed (see the superseded writeup that
used to be here, in git history); **fixed on explicit user request** after
a systematic check confirmed no other unit in the codebase had the same
gap. Details below.

### `EnzymaticPretreatment`/`PeroxidePretreatment` — maintenance OPEX never reached the TEA — FIXED

`_cost()` in both units computed `annual_maintenance =
maintenance_frac_of_capex_per_yr * capex_usd` and stored it in
`design_results`, but — unlike every other capex-dependent-maintenance unit
in this codebase (`HeatingPretreatment`, `PressateConcentrator`,
`BiogasUpgrading`, `AnaerobicDigester`'s volume-based maintenance,
`H2SRemoval`'s reagent cost) — never assigned it to `self.add_OPEX`. Since
`capex_usd` is genuinely nonzero at these units' one real call site
(`assumptions.yaml`: `enzymatic.capex_usd: 7,280,000`,
`peroxide.capex_usd: 1,000,000`), this silently omitted real maintenance
cost from VOC/MSP for every pretreatment case that includes an enzymatic or
peroxide stage — `enzymatic`, `peroxide`, `combined_PE`, `combined_PTE`
(everything except `press_mill_only`, which has no pretreatment stage at
all).

**Systematic check for the same gap elsewhere, before fixing**: grepped
every file in `units/` for `maintenance_frac_of_capex_per_yr`/
`Annual maintenance`/`maintenance_usd_per_m3_yr` (7 files matched) and
inspected each `_cost()`. Confirmed clean (already wired to `add_OPEX`
correctly): `AnaerobicDigester`, `BiogasUpgrading`, `PressateConcentrator`,
`HeatingPretreatment`, and `BiostimulantEvaporator` (fixed earlier this
conversion, see the Pressate-biostimulant entry above). No other unit in
`units/` computes a capex-dependent maintenance figure at all (confirmed
by the same grep returning only those 7 files) — so `EnzymaticPretreatment`
and `PeroxidePretreatment` were the only two live instances of this bug.

**Fix**: both units' `_cost()` now build a combined `add_OPEX` dict
(reagent cost + maintenance, when either is nonzero) and assign it once,
instead of only ever assigning the reagent-cost half.

**Metric impact, verified by direct recomputation against
`tests.py`'s 20-test suite** (`press_mill_only`-only tests — 16 of the
20 — are completely unaffected, since that case has no pretreatment
stage; TCI/FCI/material_cost/utility_cost/mass-balance outputs are
unaffected in every test, since only OPEX/VOC and everything downstream of
it (MSP, NPV) changed):

| Test | Metric | Before (buggy) | After (fixed) | Change |
|---|---|---|---|---|
| `test_ad_tea_final_combined_pe_near_zero` | VOC | $22,433,166.02 | $22,722,966.02 | +$289,800.00/yr |
| | MSP ($/MMBtu) | 13.5431 | 13.5922 | +0.0491 |
| | MSP ($/kg CH4) | 0.71237 | 0.71495 | +0.00258 |
| | NPV @ $3/MMBtu | -$515,738,681.60 | -$518,172,980.55 | -$2,434,298.95 |
| `test_ad_biostimulant_price_build_case` (combined_PE/$1.00) | MSP ($/MMBtu) | 21.2140 | 21.2631 | +0.0491 |
| `test_ad_feed_tea_build_case` (enzymatic/-0.02) | MSP ($/MMBtu) | -12.3030 | -12.2057 | +0.0973 |
| `test_ad_heatmap_build_case` (combined_PE/0.02/0.50) | MSP ($/MMBtu) | 25.7729 | 25.8221 | +0.0491 |

The combined_PE MSP shift (+0.0491 $/MMBtu, consistent across all 3
combined_PE-based tests) reflects the *combined* enzymatic + peroxide
maintenance now flowing through; the enzymatic-only case shifts more
(+0.0973 $/MMBtu) since that test uses a `-$0.02/kg` tipping-fee feed
price, which changes the MSP-to-cost sensitivity. `tests.py`'s expected
values were updated to match the fixed (correct) behavior — anyone running
this suite against a pre-fix checkout will see the values above as the old
baseline.

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

### Preprocessing (`units/_press.py`, `units/_mill.py`)

Converted with zero metric deviation (full `tests.py` suite unchanged at
`rtol=1e-6` — both units are exercised by every pathway: AD-biogas, VFA-AD,
and integrated). `capex_model` string-dispatch removed from both classes
(matching the `AnaerobicDigester`/`AcidogenicDigester` "hardcoded anchor,
no string dispatch" convention) since each dispatched string was always
set to the same one value at every real call site (`"scaled_anchor"` for
`Press`, `"inl_hammermill_anchor"` for `Mill`), confirmed by grep across
this codebase and `assumptions.yaml` — `Press`'s second branch
(`"pca_screwpress_curve"`) was never set anywhere and is now deleted
entirely as genuinely dead code, not merely a stale default. Removing the
parameter required editing the 6 real construction call sites across 3
system-builder files (`systems/_ad_biogas_system.py`,
`systems/_vfa_ad_system.py`, `systems/_integrated_system.py` — Press and
Mill are each constructed once per file) to drop the now-nonexistent
`capex_model=...` kwarg. No other part of any of those 3 files was
touched.

**Two substantive discrepancies found, deliberately NOT fixed** (fixing
either would change simulation output — out of scope for a conversion
pass; flagging per the same precedent as the VFA fermentation subsystem's
`target_oil_and_solids_content`/`ferm_mgso4_dose` items):

1. **`Press` capex reference-capacity basis mismatch.**
   `assumptions.yaml`'s `preprocessing.press` section declares
   `basis: dry_tph` and `ref_capacity_tph_dry: 0.953` — but `Press._cost()`
   (both before and after this conversion) only ever implements a
   *wet*-tph reference capacity (`ref_capacity_tph_wet`, default 50.0, no
   yaml counterpart at all). Nothing in this codebase reads
   `ref_capacity_tph_dry` or `basis` from yaml. The $3.1M installed-cost
   anchor is therefore being scaled against a 50 wet-tph reference
   capacity in the actual simulation, not the 0.953 dry-tph reference
   capacity yaml's own structure suggests it should be. Not changed here.
2. **`Press.cake_solids_wt_frac` vs. its own citation note.** yaml's
   `preprocessing.press.cake_solids_wt_frac` is `0.15`, but the adjacent
   `sources.performance` note says "Sargassum-specific pressing data
   indicate cake solids around 26-29 wt% TS; 27 wt% TS selected as central
   baseline" — the modeled value (15%) doesn't match the citation's stated
   baseline (27%). Not changed here.

**Stale defaults corrected** (dead code — always overridden at every real
call site, zero output impact): `Press.cake_solids_wt_frac` (0.35→0.15,
*not* to be confused with discrepancy #2 above — 0.15 is what yaml
actually sets, the mismatch is between yaml's own value and yaml's own
citation note, not between the unit default and yaml), `capex_installed_ref_usd`
($5,000,000→$3,100,000), `power_kWh_per_dry_ton_TS` (None→5.0),
`solids_IDs` (dropped `Xylan`/`Mannan`/`Galactan`, absent from yaml's list
— same pattern as `AnaerobicDigester`/`AcidogenicDigester`, but note
`Press.solids_IDs` *is* actively read by `_run()`, unlike
`DigestateScrewPress.solids_IDs`); `Mill.loss_frac` (0.15→0.03),
`power_kWh_per_dry_ton_dry` (None→25.0). `ref_capacity_tph_wet` (50.0) was
*not* changed — see discrepancy #1 above, there is no yaml value to align
it to.

**Citation attribution note**: `assumptions.yaml`'s `preprocessing.sources`
block is indented as a sibling of `press:`/`mill:` (not nested inside
`mill:` specifically), but its three citations (Oyedeji et al. 2020/NREL
FY11 for grinding power, INL hammer mill cost anchor for capex, INL
preprocessing/logistics accounting for loss) are all specific to `Mill`'s
own parameters, verified by reading the citation text itself rather than
inferring from yaml indentation — attributed to `Mill` in its docstring,
not to "preprocessing in general," to avoid the kind of citation
cross-contamination flagged earlier in this file for
`VFAMicrofilter`/`PressateConcentrator`.

### Pressate-biostimulant side-stream (`units/_pressate_concentrator.py`, `units/_biostimulant_evaporator.py`)

Converted with zero metric deviation (full `tests.py` suite unchanged at
`rtol=1e-6` — `PressateConcentrator` is exercised by every pathway that
enables the biostimulant side-stream, including `test_ad_biostimulant_price_*`
and `test_vfa_price_scenarios_*`; `BiostimulantEvaporator` carries zero test
coverage, see below).

**`VFAMicrofilter`'s cross-contaminated citation, now correctly attributed
here.** The VFA fermentation subsystem's conversion notes (above) already
documented that `VFAMicrofilter`'s design flux citation (Sievers et al.
2017; Diaz-Reinoso and Dominguez 2020) actually belonged to
`PressateConcentrator`. Confirmed and closed out here: `PressateConcentrator
.design_flux_L_m2_h` (35.0) matches `assumptions.yaml`'s
`pressate_biostimulant.concentrator.sources.flux` citation exactly — that
citation is now correctly attached to this unit's docstring.

**`PressateConcentrator` stale defaults corrected** (dead code — always
overridden at every one of the 3 real call sites, zero output impact):
`water_recovery_to_permeate` (0.70→0.93), `electricity_kWh_per_m3_feed`
(0.8→2.5), `capex_usd_per_m2` (120.0→500.0). Kept as custom `_design`/`_cost`
rather than converting to `@cost`, same reasoning documented for
`BiogasUpgrading`: its maintenance OPEX depends on its own computed CAPEX.

**New citation-accuracy finding: a duplicate-key yaml bug, verified by
actually parsing the file, not just reading it.**
`assumptions.yaml`'s `pressate_biostimulant.concentrator.sources.capex`
block contains two `citation:`/`note:` key pairs under the same mapping —
the real capex citation (Sethi and Wiesner 2000; Nguyen and Yoshikawa
2020; Cheryan 1998, "$200-500/m2; $800/m2 for conservative assumption")
immediately followed by a second, unrelated citation/note pair about
maintenance fraction, using the identical key names. Standard YAML
semantics keep only the *last* of duplicate keys within a mapping, so
`yaml.safe_load()` on this file actually returns the maintenance-fraction
text for `sources.capex.citation`/`.note` — confirmed directly by running
`yaml.safe_load()` on the file during this conversion, not inferred from
indentation or assumed from the raw text. The real Sethi/Wiesner/Nguyen/
Cheryan citation is present in the raw file text but is invisible to any
code that actually parses the yaml. Not fixable here — editing
`assumptions.yaml` is out of scope for this conversion regardless of
whether the edit would be a bug fix. Documented in
`units/_pressate_concentrator.py`'s docstring so the correct (raw-text)
citation is at least discoverable, with an explicit warning that it is not
what `yaml.safe_load()` actually returns.

**`BiostimulantEvaporator` is never instantiated anywhere in this
codebase** (same situation as `DigestateDecanterCentrifuge`, logged
above): both real call sites gate construction behind
`evA.get("enabled", False)`, and `assumptions.yaml` has **no**
`pressate_biostimulant.evaporator` section at all (confirmed by grepping
the whole yaml file for "evaporator" — zero matches), so `evA` is always
`{}` and the gate is always `False`. Unlike `DigestateDecanterCentrifuge`,
there isn't even an orphaned/unused yaml block to cross-reference this
unit's defaults against — every parameter is a pure code-level scoping
assumption with no yaml or literature grounding whatsoever, documented as
such rather than left implying a citation exists.

**Functional fix (not a metric-deviation bug, since the unit is never
simulated by any code path): `BiostimulantEvaporator._cost()` computed
annual maintenance into `design_results` but never assigned it to
`self.add_OPEX`**, unlike every other capex-dependent-maintenance unit in
this codebase (`PressateConcentrator`, `BiogasUpgrading`). Fixed to match
that pattern. Zero regression risk (confirmed by the unchanged 20-test
suite) since this unit has no test coverage and isn't constructed by any
system builder — flagged here rather than silently fixed, since it is a
genuine behavioral change to the unit's own logic, not merely a docstring
or default-value correction, even though it currently affects no measured
output.

### Pretreatment units (`units/_enzymatic_pretreatment.py`, `units/_peroxide_pretreatment.py`, `units/_heating_pretreatment.py`)

Converted with zero metric deviation from the docstring/stale-default work
itself. The live, verified maintenance-OPEX bug in
`EnzymaticPretreatment`/`PeroxidePretreatment` — found during this same
pass and initially left unfixed, then fixed on explicit user request, with
`tests.py`'s baseline updated to match — is logged above under
Metric-deviation bugs with full before/after values, not here, since it's
a real correctness defect rather than a citation/documentation issue.

All three units are constructed in exactly one place in this codebase —
`systems._ad_biogas_system._build_methanogenic_pathway` — confirmed by
grep. None needed a `capex_model` string-dispatch removal (they never had
one; `capex_usd` was already a flat `__init__` default in all three, just
stale/zero) and none were converted to an `@cost` decorator: `capex_usd`
is a flat, non-scaling placeholder cost in all three, not a function of
any design basis, so the decorator's power-law machinery doesn't apply.

**`HeatingPretreatment` is only ever reachable via the `combined_PTE`
pretreatment case's peroxide→heating→enzymatic sequence.**
`_build_methanogenic_pathway` has a standalone `pt_kind == "heating"`
dispatch branch with its own inline fallback values (338.15 K / 0.25 h,
matching this unit's pre-conversion defaults exactly) — but
`assumptions.yaml` defines no `ad_pretreatment_cases` entry with
`kind: heating` (only `press_mill_only`/`enzymatic`/`peroxide`/
`combined_PE`/`combined_PTE`), so that branch is dead code, unreachable
with the current yaml configuration. The unit's real (and only reachable)
values come from `combined_PTE.heating` instead (393.15 K / 0.25 h /
$1.2M), which is what its corrected defaults now reflect.

**Stale defaults corrected** (dead code — always overridden at the one
real call site, zero output impact, confirmed by the unchanged test
suite): `EnzymaticPretreatment.enzyme_dose_kg_per_kg_dry_feed` (0.02→0.005),
`treated_fraction` (1.0→0.25), `enzyme_recycle_factor` (1.0→2.0),
`capex_usd` (0.0→$7,280,000); `PeroxidePretreatment.capex_usd`
(0.0→$1,000,000); `HeatingPretreatment.target_temperature_K`
(338.15→393.15 K, see reachability note above), `capex_usd`
(0.0→$1,200,000).

**Citation scope note**: `assumptions.yaml`'s `ad_pretreatment_cases.*`
`sources.methane_yield`/`biogas_composition`/`vs_destruction` entries
(Chikani-Cabrera et al. 2022) describe the *downstream AD effect* of each
pretreatment case (i.e. parameters consumed by `AnaerobicDigester` via
`ad_effects`), not these pretreatment reactors' own operating parameters
(temperature, residence time, dose, recycle factor). Not cited in these
three units' docstrings, to avoid attributing an AD-outcome citation to a
pretreatment-reactor parameter it doesn't actually describe.
