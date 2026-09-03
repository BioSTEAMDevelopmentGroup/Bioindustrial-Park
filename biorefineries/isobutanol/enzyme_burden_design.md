# Design: proteome / enzyme-burden coupling for the kinetic Bayesian optimization

**Status:** draft for review · **Author:** Sarang Bhagwat (with Claude) · **Date:** 2026-09-03
**Scope:** enzyme-burden constraint in depth; the other five constraint families as a
prioritized roadmap (§8).

---

## 1. Problem statement

The kinetic BO (`kinetic_optimization.py`) samples every Vmax-type rate constant
`k_*` independently on `[0.1×, 10×]` its scenario baseline, log-scale, with **no
coupling between them and no cost on growth**
([`build_search_space`](kinetic_optimization.py:130)). The r7 biomass yield
(`0.732` in `s_glu => 0.732 x + …`) is a fixed stoichiometric constant, and the
growth rate `k_7 = 1.203` is a free, uncoupled variable. A single trial can
therefore raise heterologous flux (e.g. 10× the Ehrlich branch `k_13…k_16`) **and**
keep growth untouched.

No real strain behaves this way. Enzyme expression draws on a finite proteome and
energy budget, so pushing expression into any pathway steals from growth. The
observed consequence in production runs is that the IRR optimum walks to a
biologically infeasible corner (kin_opt_B_irr: EtOH abandoned, Ehrlich branch
maxed, IRR 0.39). Adding a proteome-burden coupling makes the optimum represent an
**achievable engineering target** rather than a proteome-infeasible one — without
changing the pinned baselines (see §5.3, the inert-at-baseline property).

## 2. Design decisions (locked)

| Decision | Choice | Rationale |
| --- | --- | --- |
| Growth channel | **Growth *rate* (r7), not yield (0.732)** | The proteome-allocation growth law (Scott/Hwa) is fundamentally a *rate* law: heterologous expression steals ribosomes → lower µ. `k_7` is a plain parameter, so no runtime stoichiometry mutation is needed. The yield channel (maintenance-energy cost on 0.732) is deferred to the roadmap (§8). |
| Burdened set | **All raised Vmax `k_*`** (per-reaction), increment above baseline | Total enzyme investment is capped; pouring expression into *any* pathway — native overexpression or heterologous — taxes growth. |
| Coupling locus | **Mechanistic, inside NSKinetics** | Self-consistent: the burden feeds back through the dynamic simulation (lower µ → less biomass → less flux → lower titer), exactly as an enzyme-constrained model behaves. Interpretable and testable. |
| Doc scope | Burden in depth + roadmap | This document. |

## 3. Biological and literature basis

**Proteome allocation / growth laws (primary basis).** Cells allocate a finite
proteome among sectors; the fraction devoted to non-growth (heterologous or
over-expressed) protein reduces the ribosomal/growth sector, and growth rate falls
**linearly** with that fraction:

> λ / λ₀ = 1 − φ / φ_max

- Scott, Gunderson, Mateescu, Zhang & Hwa. *Interdependence of cell growth and gene
  expression: origins and consequences.* **Science** 2010, 330:1099–1102. (The
  empirical linear growth law; φ_max ≈ 0.48 in *E. coli*.)
- Weiße, Oyarzún, Danos & Swain. *Mechanistic links between cellular trade-offs,
  gene expression, and growth.* **PNAS** 2015, 112:E1038–E1047. (Mechanistic
  derivation of the same trade-off.)
- Basan, Hui, Okano, Zhang, Shen, Williamson & Hwa. *Overflow metabolism in E. coli
  results from efficient proteome allocation.* **Nature** 2015, 528:99–104.

**Enzyme-constrained genome-scale models (the same idea, quantified for yeast).**
Total enzyme mass is bounded: Σ_i (v_i / k_cat,i)·MW_i ≤ P_total. Raising one flux's
enzyme demand forces others (including growth-machinery) down.

- Sánchez, Zhang, Nilsson, Lahtvee, Kerkhoven & Nielsen. *Improving the phenotype
  predictions of a yeast GEM by incorporating enzymatic constraints (GECKO).*
  **Mol. Syst. Biol.** 2017, 13:935.
- Nilsson & Nielsen. *Metabolic trade-offs in yeast are caused by F1F0-ATP
  synthase.* **Sci. Rep.** 2016, 6:22264. (Proteome-constrained yeast overflow.)

**The burden is real and measurable, including in yeast.**
- Ceroni, Algar, Stan & Ellis. *Quantifying cellular capacity identifies gene
  expression designs with reduced burden.* **Nat. Methods** 2015, 12:415–418.
- Kafri, Metzl-Raz, Jona & Barkai. *The cost of protein production.* **Cell Rep.**
  2016, 14:22–31. (Yeast; growth cost of unneeded expression.)

**Isobutanol-strain evidence (candidate calibration anchors, §6).**
- Avalos, Fink & Stephanopoulos. *Compartmentalization of metabolic pathways in
  yeast mitochondria improves branched-chain alcohol production.* **Nat.
  Biotechnol.** 2013, 31:335–341.
- Brat & Boles. *Isobutanol production from cellobiose in S. cerevisiae.* **FEMS
  Yeast Res.** 2013 (and related cytosolic-pathway strains) — report growth/biomass
  penalties on Ehrlich-pathway overexpression.

**Assumption — quasi-steady proteome allocation.** The growth law above is a
*steady-state* proteome partition; we apply it pointwise in a dynamic fed-batch.
This is justified because proteome reallocation (minutes to ≲1 h) is fast relative
to the fermentation timescale (40–80 h), and the model already carries fast enzyme
pools (`X_a`, `X_AcDH`) with their own sub-hour dynamics. Stated explicitly so it
can be challenged.

**Why rate, not yield (record of the deliberation).** Both channels are defensible:
proteome allocation lowers *rate* (this design), while the ATP/maintenance cost of
synthesizing protein lowers *yield* (Pirt maintenance framework — Pirt, *Proc. R.
Soc. B* 1965, 163:224). We take the rate channel because (a) it is the primary,
best-quantified proteome-law effect and (b) it needs no stoichiometry mutation. The
yield channel is a clean roadmap extension (§8).

## 4. What must NOT be double-counted

The mechanism already penalizes heterologous flux three ways; the burden term must
add **only** the proteome/expression cost that is *not* already present:

1. **Carbon competition** — `r13` draws from the same pyruvate pool as `r2`/`r3`,
   so diverting carbon to isobutanol already starves TCA/growth. Do not re-charge
   carbon.
2. **Product inhibition of growth** — `k_7ii·s_IBO` (and `k_7ie`, `k_7ia`) already
   throttle r7 as products accumulate.
3. **Product-enhanced death** — `r10` with `P_10i` (IBO), `P_10e` (EtOH), `P_10a`
   (acetate).

The burden term is the **cost of expressing the enzyme itself** (ribosome/proteome
diversion), independent of what the enzyme's product later does. Keeping these
separate is what makes the constraint defensible rather than a second, hidden IBO
penalty.

## 5. Formal specification

### 5.1 Burden index φ

Over the burdened reaction set **V** (§5.4), define a dimensionless expression-burden
index from the *increment above the scenario baseline*:

> **φ = Σ_{i∈V} w_i · max(0, (k_i − k_i,base) / k_i,ref)**

- `k_i,base` — the parameter's **scenario baseline**, snapshotted at load time
  (§7.2), not the antimony literal. This makes φ = 0 at baseline in **both**
  scenarios.
- `k_i,ref` — the normalization scale. `k_i,ref = k_i,base` when `k_i,base > 0` (so
  the term is the fold-overexpression `m_i − 1`). When `k_i,base = 0` (a pathway
  *off* at baseline, e.g. Ehrlich in scenario A), the ratio is undefined; use a
  declared reference expression scale `k_i,ref` for that reaction so "switching on"
  a pathway carries a finite, defensible cost.
- `w_i` — per-reaction proteome weight. **Default `w_i = 1`** (transparent first
  pass). MW/abundance weighting (`w_i ∝ MW_i·[E]_i,base`, from PaxDb yeast
  abundances or ecYeast enzyme usage) is a refinement (§6) — and is ill-defined for
  the lumped native reactions anyway (see §5.4), which is a second reason to start
  equal-weighted.

Only positive increments count (`max(0, …)`): *reducing* a k below baseline frees
proteome but we do not credit growth for it (conservative; downregulation savings
are small and uncertain).

### 5.2 Realized growth-rate multiplier

> **f_burden = clip(1 − φ / φ_max, f_min, 1)**

Applied as a multiplicative factor on **r7's kinetic law** (the growth reaction):

```
r7_rate = f_burden · (anaerobic_growth_mult + …) · (k_7·s_glu·exp(−k_7ie·s_EtOH)·… / (s_glu + K_7)) · a · env
```

`k_7` (the optimizer's setpoint = *intended/expressed* µmax) is left untouched;
`f_burden` scales the *realized* rate. Parameters:
- `φ_max` — total tolerable expression index before growth hits the floor
  (calibration constant, §6).
- `f_min` — growth-rate floor (∈ ~0.1–0.2): burdened cells still grow slowly, and
  it keeps the integrator away from exactly zero growth.

### 5.3 Properties (each is a test in §7.3)

1. **Inert at baseline (critical).** All `m_i = 1` ⇒ φ = 0 ⇒ f_burden = 1 ⇒ every
   baseline simulation and every pinned smoke-test MPSP is unchanged **by
   construction**, in both scenarios (because `k_i,base` is the per-scenario
   snapshot).
2. **Monotone & bounded.** Raising any Vmax above baseline can only lower growth;
   f_burden ∈ [f_min, 1].
3. **Self-consistent on k_7.** `k_7 ∈ V`, so dialing up µmax itself raises φ and is
   partially self-cancelling — realistic diminishing returns on growth-rate
   engineering.
4. **Feeds back dynamically.** Lower µ → less biomass → less total flux → lower
   titer, propagated through the simulation, not bolted onto the objective.

### 5.4 The burdened set V (parameter selection)

**Rule:** the Vmax-type constants — those declared `has g_per_l_per_h` in the
antimony. Explicitly:

```
k_1l, k_1h, k_1e   (r1 glycolysis, lumped)
k_2                (r2 TCA-pyruvate)
k_3                (r3 PDH)
k_4                (r4 acetaldehyde→acetate)
k_5, k_5e          (r5 TCA-acetate)
k_6                (r6 ADH)
k_7                (r7 growth-glucose)      ← also the reaction being scaled
k_8                (r8 growth-acetate)
k_9, k_9e, k_9c    (r9 AcDH production)
k_13, k_14, k_15, k_16   (Ehrlich branch → isobutanol)
```

**Excluded** (not Vmax / not biosynthetic expression):
- All `l_per_g` and dimensionless constants — `k_*ie/ia/ii` (inhibition), `k_6r`,
  `k_16r` (reversibility), `K_*` (affinities). These are enzyme *properties*, not
  expression levels.
- **`k_10` (death) and `k_11` (AcDH degradation)** — degradation/decay, not
  biosynthetic expression. *Recommend excluding; flagged as open decision O1.*

**Lumped reactions (open decision O2).** r1 carries three Vmax terms (`k_1l/h/e`)
for one lumped step; equal-weighting would triple-count glycolysis. Recommended
handling: aggregate to a **per-reaction** contribution — sum the increments within a
reaction, or weight each of a reaction's Vmax terms by `1/n_terms`. Documented so the
choice is explicit.

## 6. Calibration

Land the plumbing first with the burden **disabled** (`burden_on = 0`, §7), then fit
these before enabling:

- **φ_max** — sets how much total expression the proteome tolerates. Scott/Hwa's
  φ_max ≈ 0.48 is a literal proteome *fraction*; our φ is an index, so φ_max is a
  fitting constant, not that number directly.
- **f_min** — growth floor.
- **weights `w_i`** — equal (default) or MW/abundance-weighted.

**Fitting procedure.** Hold `w_i = 1`, choose (φ_max, f_min) so a *defined
engineering scenario* reproduces a measured growth/biomass penalty:
- Anchor point: one well-characterized isobutanol *S. cerevisiae* strain's reported
  (fold-expression of the Ehrlich branch, growth-rate or final-biomass reduction) —
  candidates in §3 (Avalos 2013; Brat & Boles 2013). **Open decision O3: which
  dataset.**
- Sanity band: the general yeast burden curve (Kafri 2016) — a few % heterologous
  proteome should already give a measurable (single-digit-%) µ cost, and the full
  10× pathway corner should give a large (tens-of-%) cost, not a total collapse.

Until O3 is resolved, φ_max/f_min stay unset and `burden_on = 0`.

## 7. Implementation plan

### 7.1 NSKinetics (antimony source of truth; regenerate SBML)

Edit
`nskinetics/models/s_cerevisiae_ferm_fb_inhib_mod_ibo/s_cerevisiae_ferm_fb_inhib_mod_ibo_antimony.txt`:

1. Declare, for each `k_i ∈ V`, a baseline `k_i_base` (default = current literal)
   and a `k_i_ref` (for zero-baseline pathways).
2. Declare burden constants: `phi_max`, `f_min`, weights `w_*` (default 1), and a
   master switch **`burden_on` (default 0)** so the model is byte-for-byte the
   current behavior until explicitly enabled — this protects every existing pinned
   result and lets the plumbing land before calibration.
3. Add assignment rules `phi_burden := Σ …` and
   `f_burden := piecewise(clip(1 − phi_burden/phi_max, f_min, 1), burden_on == 1, 1)`.
4. Multiply **r7**'s rate law by `f_burden`.
   - *Open decision O4:* also multiply **r8** (growth on acetate) for consistency —
     recommended, since it is also a growth reaction; defaulted to r7-only to match
     the "k_7" framing.
5. Expose `phi_burden`, `f_burden` as variables so they can be tracked per trial.
6. Regenerate the SBML (`…_sbml.xml`) from the antimony and bump the nskinetics local
   clone version (0.6.0 → 0.6.1).

### 7.2 isobutanol package

1. In `system.load()`, **after** the scenario parameter distributions are applied
   and **before** the baseline simulation, snapshot each `k_i ∈ V` into `k_i_base`
   on `r_te` (per-scenario baseline ⇒ f_burden = 1 at baseline). Leave
   `burden_on = 0` until calibrated.
2. `kinetic_optimization.py`: add `phi_burden` and `f_burden` to `TRACKED_METRICS`
   so the trajectory CSV records realized burden per trial. **No structural change**
   otherwise — the BO already sets `k_*` on `r_te`, and the mechanistic `f_burden`
   does the rest; `restore_baseline`'s `finally` already resets the `k_*`, so burden
   auto-resets.

### 7.3 Verification

1. **Plumbing inert (`burden_on = 0`):** all 8 smoke tests reproduce current pinned
   MPSPs exactly. Gate for landing the plumbing.
2. **Inert at baseline (`burden_on = 1`):** at the baseline draw, `f_burden == 1` to
   machine precision in both scenarios; smoke tests still pass.
3. **Mechanism check:** manually 10× the Ehrlich branch ⇒ `f_burden < 1`, µ / final
   biomass / titer drop; the penalized IRR surface should pull the optimum back from
   the previous unconstrained corner (kin_opt_B_irr).
4. **Offline unit test:** the φ / f_burden algebra with fakes (extend
   `analyses/test_kinetic_optimization_offline.py` or add an nskinetics-level test).

### 7.4 Staging

**Stage 1** — land §7.1–7.2 with `burden_on = 0`; verify §7.3(1). No baseline moves.
**Stage 2** — resolve O1–O4, calibrate (§6), set `burden_on = 1`; verify §7.3(2–4);
re-review and re-pin any smoke test whose non-baseline gate legitimately moves.

## 8. Roadmap — the other five constraint families

Prioritized; each defensible, several cheap because the model already exposes the
quantity:

1. **Yield channel of the burden** — maintenance-energy cost lowers Y_x/glu (the
   `0.732`) alongside the rate channel here. Requires runtime stoichiometry mutation
   or a parameterized yield. Basis: Pirt 1965; GECKO growth-associated maintenance.
2. **k_cat ceiling** — cap individual-k increases interpreted as catalytic
   improvement; beyond the cap, an increase is *expression* → already burdened.
   Basis: Bar-Even et al., *Biochemistry* 2011, 50:4402 ("the moderately efficient
   enzyme").
3. **Thermodynamic reversibility bound** — `r6`, `r14`, `r15`, reductive `r16` are
   near-equilibrium (`k_6r`, `k_16r`); bound forward/reverse and driving force.
   Basis: Noor et al., *PLoS Comput. Biol.* 2014 (max-min driving force).
4. **Redox/cofactor balance** — `r16` consumes `0.363 $Red`; bound total reductive
   draw against regeneration capacity (a known cytosolic-isobutanol bottleneck).
5. **Oxygen transfer (OTR/kLa)** — the model already computes `qO2`; impose
   OTR_max so a high-flux optimum stays bioreactor-feasible. Cheap to add.

Softer systems-level envelope: the **yield–rate Pareto** trade-off (Molenaar et al.,
*Mol. Syst. Biol.* 2009; Pfeiffer et al., *Science* 2001) as a sanity bound.

## 9. Open decisions (need answers before Stage 2 / enabling burden)

- **O1** — Exclude `k_10` (death) and `k_11` (AcDH degradation) from V? *Rec: yes.*
- **O2** — Lumped-reaction handling for r1's `k_1l/h/e` (per-reaction aggregation vs
  per-term). *Rec: per-reaction (sum increments, or 1/n_terms weights).*
- **O3** — Calibration anchor dataset (which isobutanol strain's expression↔growth
  point). *Needs your input.*
- **O4** — Apply `f_burden` to r8 (acetate growth) as well as r7? *Rec: yes.*
- **O5** — φ_max, f_min numeric values (follows from O3). *Land plumbing with
  `burden_on = 0` first.*

## 10. References

See §3 and §8 (inline). Model provenance: Lei, Rotbøll & Jørgensen, *J. Biotechnol.*
2001, 88:205 (base structured yeast model); this fork
`Bhagwat2026_Yeast_Ethanol_Isobutanol`.
