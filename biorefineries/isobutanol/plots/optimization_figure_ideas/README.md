# Kinetic-BO campaign figures: "Scouts explore, profit selects"

Publication figures for the eight `metabolic_split_12d` Gaussian-process kinetic
(strain-design) optimization campaigns. These are the seed-350 uninformed profitability
campaign, the six seed-350 TRY "scout" campaigns, and the TRY-informed relay campaign
(`_rl15c111dc`), which was preloaded with 1,000 scout trials. The replicate campaigns used
in Fig. S1 are the untagged 09-16/17 runs and the `_rlba1b2315` relay.

Everything here is **sim-safe**. It reads the campaign CSVs with pandas. It never imports
`biorefineries.*`, nskinetics, biosteam, thermosteam or optuna, and it never calls `load()`.
`plots/plot_kin_opt_parameter_sets.py` (pk) is loaded by file path, and only lazily, for
the scenario-A baseline. That means any script here can run on any numba-cache state,
alongside a running simulation.

## Files

| file | role |
|---|---|
| `_common.py` | Data layer. It holds the paths and the campaign registry, the cached loaders, and the one set of definitions. It also holds the proteome sectors, the decision bands, `compute_facts()`, `EXPECTED` and `check_facts()`. The public API is documented at the top of the file. |
| `_style.py` | Style layer. It holds the palette, rcParams, `style_ticks` and inch-based placement. It provides the IRR / log-trial / logit axes, the loss band, plateau and tint helpers, and `inline_key`. It also holds the render checks (`text_overlaps`, `tick_label_collisions`, `min_font_check`, `glyph_check`, `marker_text_hits`, `check_figure`), the palette check (`cvd_check`) and `save()`. The public API is documented at the top of the file. |
| `check_facts.py` | The folder's offline test. It prints every computed-vs-expected fact, the hard-coded-literal scan of the figure scripts, the palette check and the sim-safety check. It exits 1 on any mismatch. |
| `fig_main_scouts.py` | Main figure (panels a–f), stem `kinBO_main_scouts`. |
| `fig_s1_replicates.py` | Fig. S1, the replicate / honesty check, stem `kinBO_S1_replicates`. |
| `fig_s2_fingerprints.py` | Fig. S2, the design-fingerprint heatmap, stem `kinBO_S2_fingerprints`. |

The full figure specification is in `figwork/figure_spec.md` in the session scratchpad. It
is not tracked.

## Running

```powershell
$py = "C:\Users\saran\anaconda3\envs\IBO_2026\python.exe"
& $py plots\optimization_figure_ideas\check_facts.py        # exit 0 = ALL CHECKS PASSED
& $py plots\optimization_figure_ideas\fig_main_scouts.py
```

Optional flags for `check_facts.py`:

* `--failures-only` prints only the facts that differ from `EXPECTED`.
* `--no-literal-scan` skips the literal scan.
* `--no-palette` skips the palette check.

A full check takes about 1–2 s.

Every figure script must:

1. call `_common.check_facts()` before drawing. It raises `FactsMismatch` on any mismatch.
2. build every annotation number from the facts dict with `fmt_irr` / `fmt_int` / `fmt_pct`.
   The literals `15.8`, `27.3`, `22.9`, `563`, `296`, `110` and `97` must not appear in
   the code. `check_facts.py` scans for them.
3. call `_style.check_figure(fig, size=...)` after building. This draws the figure and
   checks text overlaps, same-axis tick collisions, the 9-pt font floor, Arial glyph
   coverage, markers under text and the canvas size. It raises `FigureCheckError` on
   any problem.
4. call `_style.save(fig, stem)`.

Minimal skeleton of a figure script:

```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C
import _style as S

facts = C.check_facts()
S.cvd_check()
fig = S.new_figure(*S.MAIN_SIZE)
ax = S.inch_axes(fig, 1.17, 4.55, 4.20, 2.85)
S.log_trial_axis(ax)
S.irr_axis(ax, 'y')
S.plateau_line(ax, facts['U_pct'])
S.style_ticks(ax)
...
S.check_figure(fig, size=S.MAIN_SIZE)
S.save(fig, 'kinBO_main_scouts')
C.assert_sim_safe()
```

## Outputs

Outputs go to `analyses/results/publication/Optimization-figures/`, which is gitignored with
the rest of `results/`.

* `save()` writes `<stem>_<YYYY.MM.DD-HH.MM>.png` (300 dpi) and `.pdf`. It never uses
  `bbox_inches='tight'`, because the layout is fixed in inches.
* It also writes stable copies `<stem>_latest.png` / `.pdf`, overwritten on every run, so
  reviewers can always find the newest render.
* It asserts the following:
  * every path is at most 259 characters (LongPathsEnabled is off on this machine);
  * the PNG is exactly the canvas size × 300 dpi (main figure: 3000 × 2370 px);
  * the PDF is at most 10 MB;
  * the PDF has no Type-3 fonts (`pdf.fonttype = 42`, so Arial is embedded and the text
    is selectable).

## Data sources

All files are in `analyses/results/`. Study stems have the form
`kin_opt_ethanol_isobutanol_metabolic_split_12d_<slug>_gp_rb0.001-4_ib0.75-1.5_aA<tag>_burden`.

| key | slug / tag | figure label |
|---|---|---|
| `unin` | `pi_log-tail` / `_rs350` | uninformed (profitability) |
| `relay` | `pi_log-tail` / `_rl15c111dc` (+ `_relay_manifest.csv` = the 1,000 seeds) | TRY-informed (profitability) |
| `ey` `et` `ep` | `etoh_{yield,titer,productivity}` / `_rs350` | Ethanol TRY scouts |
| `iy` `it` `ip` | `ibo_{yield,titer,productivity}` / `_rs350` | Isobutanol TRY scouts |
| `unin_rep`, `relay_rep` | `pi_log-tail` / none and `_rlba1b2315` | Fig. S1 replicates |
| `rep_iy` … `rep_ep`, `rep_pw` | the untagged scouts + `price-weighted_yield` | Fig. S1 replicate scouts |

* The scenario-A starting strain comes from `pk.baseline_set()` and `pk.BASELINE_A`. Its
  IRR is 12.60 % under the model version the campaigns ran with (12.79 % under the
  current one).
* The decision values k_3 (5.81) and k_6 (2.82) come from the scenario-A workbook.
* k_17 (0.1077) is `pk.SPLIT_12D_REL_RATE_BASELINES_A`.

## Definitions (identical in every figure, panel and caption)

**Rows and trial indices**
* Only `state == 'COMPLETE'` rows are used. FAIL rows are excluded everywhere.
* `sim` is the simulated index: `trial_number - min(trial_number) + 1`, 1-based. For the
  relay this is `trial_number - 999`.
* In figure text, "trial N" means `sim`.
* The `trial` fields in the facts are CSV trial numbers. For the 2,000-trial campaigns,
  `sim = trial + 1`; for example, the isobutanol-yield best visit is trial 842 = sim 843.

**IRR and loss**
* IRR is shown in % (`100·IRR`).
* A **loss** is `IRR < 0` or `IRR = -inf` (no IRR: an outright money-loser). Losses are
  never clamped to 0.
* Losses are drawn in a grey band from −6 to 0 %, centred at −3 %. Dots are jittered by
  ±2.4 with `np.random.default_rng(0)`, one generator per panel; lines are not jittered.
  The axis tick label is "loss".

**Plateau and best so far**
* **Plateau U** is the uninformed campaign's maximum IRR: 15.8254 % at sim 104 (trial 103).
* "Above the plateau" means `IRR > U + 1e-9`.
* **Best so far** is `np.maximum.accumulate` of IRR over COMPLETE rows ordered by `sim`,
  with `-inf` kept. For both profitability campaigns it equals the IRR of the running
  objective incumbent at every row; this is asserted.

**Designs**
* The **returned design** is the argmax of the campaign's own `objective`.
* The **best visit** is the argmax of IRR. For the relay, only simulated rows count.

**Products** (titers are g/L of water)
* **Makes isobutanol**: IBO titer ≥ 5 g/L.
* **Co-production**: IBO ≥ 5 and EtOH ≥ 5 g/L.
* **Isobutanol share**: `100·IBO/(IBO+EtOH)`, with titers clipped at 0. It is defined only
  when the sum is at least 1 g/L.

**Campaign-level measures**
* **Exploration** is the % of a campaign's COMPLETE trials that make isobutanol.

**Proteome sectors** (g protein·(g DCW)⁻¹)

| sector | pools |
|---|---|
| glycolysis | r1 |
| TCA/acetate | r2 + r4 + r5 |
| Pdc | r3 |
| Adh1 | r6 |
| ALS→Aro10 | r13–r16 |
| Adh6 | r17 |

* The sectors sum to Φ_M.
* φ_T = 0.11025 is constant.
* The penalty-free budget is F_flex − φ_T = 0.13475.

## Palette and render rules

* House hexes are used for the uninformed (cyan `#18C4DC`) and TRY-informed (teal
  `#0B6E7A`) campaigns and for the baseline grey.
* The scouts are coloured by product family: ethanol amber, isobutanol violet. Within a
  family, the role is the marker: ○ yield, □ titer, △ productivity.
* `cvd_check()` requires CIE76 ΔE ≥ 12 for the seven co-occurring pairs. This holds under
  normal vision and under simulated deutan, protan and tritan vision (Machado 2009). The
  worst pair is adh1/adh6 under tritan, at 17.7.
* The grayscale lightness gaps are 30.6 for cyan/teal and 25.5 for amber/violet.
* The dataviz skill's validator is run as an advisory check when it is available.
* All text is in Arial, at 9 pt or larger on the 10-in canvas (6.4 pt at the 180-mm print
  width).
* Arial has no △ ★ ◆ ◇ glyphs, so draw them as markers (see `_style.inline_key`).
  `glyph_check` flags any character that would fall back to another font.
* Ticks follow the house rule: all four sides, left/bottom in-and-out, top/right inward,
  major 4 pt, minor 2 pt.
* On a horizontal IRR axis narrower than about 2.8 in, the "loss" and "0" tick labels
  collide. Use `irr_axis(ax, 'x', loss_fs=9)`.
* No data marker may cross a letter. `marker_text_hits` (run by `check_figure`) flags
  every scatter or marker point whose centre falls inside an annotation's text extent
  (1-pt pad), per axes, insets included, plus figure-level text over any axes. A point
  under an opaque backing box (face alpha ≥ 0.85, drawn above the point's layer) is
  hidden and passes. Put labels on empty spots first; use a backing box only where no
  spot is free, and say in the caption's methods notes what it hides. The main figure's
  panel-b labels are placed this way by `_spot` / `_best_spot`.
