# Litter v3.2 — implementation plan

Follow-up to the v3.1 review (branch `claude/litter-hazard-review-i05ken`,
commit `72cc9c9`). The review fixed the minor bugs directly; this plan
addresses the four **major/design** findings plus one documentation item. It
is written to be executed step by step; every work package (WP) states the
exact behaviour change, the tests to write first, worked oracle numbers, and
which existing pinned numbers must NOT move.

Execute the WPs **in order** (WP2 depends on WP1's constants table; the rest
are independent but ordering avoids merge conflicts in `litter_hazard.R` /
`litter_exposure.R`). Each WP is a separate commit using the repo's
conventional-commit style (`feat(litter): ...` / `fix(litter): ...`).

## Ground rules (read before starting)

1. Read first: `CLAUDE.md`, `R/litter_hazard.R`, `R/litter_exposure.R`,
   `R/constants.R` (the `DUST_CONSTANTS` block shows the constants style),
   `tests/testthat/test-litter-hazard.R`, `tests/testthat/test-litter-exposure.R`.
2. **Write tests before implementation** (repo rule). Add tests to the
   existing litter test files, in the existing `describe()/it()` style.
3. **Pinned oracles that must not change.** All of these are already asserted
   in `test-litter-hazard.R` / `test-litter-exposure.R`; if your change moves
   any of them, the change is wrong (they were verified analytically during
   the v3.1 review):
   - `litter_hazard_vec(16, 12.5, 0, 0.02)` = **86.8421**
   - `litter_hazard_vec(10, 7, 0, 0.10)` = **1.2127** (film AND paper)
   - `litter_hazard_vec(c(16,35), c(13,15), c(0,0), c(0.20,0.30), material="film")` = **c(4.5248, 30)**
   - same call with `material="paper"` = **c(0, 0)**
   - `litter_hazard_vec(16, 15, 0.4, 0.02)` = **33.333**
   - wetness path at `wetness == 1`: film = **30**, paper = **0**
   - `litter_exposure` demo-site regressions (exposure **c(86, 7.5, 2.5)** etc.)
4. New tunables go into a named list in `R/constants.R` (new
   `LITTER_CONSTANTS`, WP1), documented with a roxygen `@keywords internal`
   block mirroring `DUST_CONSTANTS`. No inline literals in physics code.
5. After roxygen edits run `devtools::document()`; after each WP run
   `devtools::test(filter = "litter")`, and `devtools::test()` before the
   final push. Add a NEWS.md bullet per WP under a `# meteoHazard (litter
   v3.2)` heading.
6. Every default introduced here is an **uncalibrated placeholder** — say so
   in its roxygen, matching the package's existing "high calibration
   uncertainty" phrasing (issues #11/#26 own calibration).

---

## WP1 — smooth material-saturation ramp (removes the interior discontinuity)

### Problem

On the wetness path, film entrainment drops by `saturation_penalty` (70%) as
a **step** at `wetness >= paper_veto_wetness` (0.8), and paper drops to zero
as a step at the same mark, while the moisture-raised threshold below the
mark is smooth. An interior cliff in an otherwise smooth index is physically
unmotivated and will fight the planned calibration.

### Design

Replace the binary `saturated` treatment with a **linear penalty ramp** on
normalised wetness `w_norm`, shared by both input paths, with the full
penalty landing exactly at `w_norm == 1`:

```
ramp(w)  = clamp((w - saturation_onset) / (1 - saturation_onset), 0, 1)
E        = E0 * (1 - pen_max * ramp(w_norm))
```

- `saturation_onset` — new formal, default `0.8`, replaces (and deletes) the
  `paper_veto_wetness` formal. Validate `0 <= saturation_onset < 1` (a
  classed `meteoHazard_input_error` when `>= 1`: the ramp denominator must be
  positive).
- `pen_max` is per material: **film 0.7** (the current
  `saturation_penalty`), **paper 1.0** (the hard veto becomes the ramp's
  endpoint). Keep the exported knob `saturation_penalty` as the *film* value;
  paper's 1.0 is a constant. Put both in `LITTER_CONSTANTS$MATERIALS` (see
  WP2 for the table's final shape; in WP1 it only needs
  `film$saturation_penalty = 0.7`, `paper$saturation_penalty = 1.0`) plus
  `LITTER_CONSTANTS$SATURATION_ONSET = 0.8`.
- Delete the `saturated` variable and the `if (material == "paper") E[saturated] <- 0`
  block entirely; both paths now flow through the one formula. On the soil
  path the ramp input is the existing Fecan-clamped `w_norm`; on the wetness
  path it is `wetness` itself.

### Why the pinned oracles survive (verified analytically)

- All dry / damp oracles have `w_norm <= 0.33 < 0.8` → `ramp = 0` → unchanged.
- Every saturated oracle sits at exactly `w_norm == 1` (SM 0.20, 0.25, 0.30
  all clamp to 1; wetness-path tests use `wetness == 1`) → `ramp = 1` → film
  gets the full 0.7 penalty, paper gets 1.0 (zero) → identical numbers.
- The moisture-monotonicity test (SM up to 0.19 → `w_norm = 0.933`) now takes
  a partial penalty — output values change *between* 0.17 < SM < 0.20 but the
  test only asserts monotone decrease, which the ramp strengthens.

### Tests to write first (add to `test-litter-hazard.R`)

```r
describe("smooth material-saturation ramp (v3.2)", {
  it("film mid-ramp: wetness=0.9, G=35, W=15 -> 65", {
    # E0 = 50 (excess 23.49 > denom 9.93), ramp = (0.9-0.8)/0.2 = 0.5,
    # E = 50*(1 - 0.7*0.5) = 32.5, T = 2 -> 65.
    out <- litter_hazard_vec(35, 15, 0, wetness = 0.9, material = "film")
    expect_equal(out, 65, tolerance = 1e-3)
  })
  it("paper mid-ramp: wetness=0.9 -> 50 (partial, no longer a hard veto)", {
    out <- litter_hazard_vec(35, 15, 0, wetness = 0.9, material = "paper")
    expect_equal(out, 50, tolerance = 1e-3)
  })
  it("no step at the onset: hazard is continuous and monotone over wetness [0.75, 1]", {
    w <- seq(0.75, 1, by = 0.01)
    out <- litter_hazard_vec(rep(35, length(w)), rep(15, length(w)),
                             rep(0, length(w)), wetness = w)
    expect_true(all(diff(out) <= 0))
    expect_lt(max(abs(diff(out))), 4)  # 0.01-step bound; a cliff would be ~35
  })
  it("errors (classed) when saturation_onset >= 1", {
    expect_error(
      litter_hazard_vec(16, 15, 0, wetness = 0.5, saturation_onset = 1),
      class = "meteoHazard_input_error"
    )
  })
})
```

### Implementation steps

1. Add the `LITTER_CONSTANTS` list + roxygen block to `R/constants.R`.
2. In `litter_hazard_vec()`: swap the `paper_veto_wetness` formal for
   `saturation_onset = 0.8`; replace the saturated-block with the ramp
   formula; update validation.
3. Grep for `paper_veto_wetness` across `R/`, `tests/` (one comment at
   `test-litter-hazard.R:140` references it), `man/` — update or delete every
   hit. The formal was introduced on the unmerged v3.1 branch, so no
   deprecation shim is needed.
4. Update roxygen: the "Material-aware graded saturation" Model bullet, the
   `@param saturation_penalty` text, the new `@param saturation_onset`.
   `devtools::document()`.

**Done when:** new tests pass, full litter filter passes with zero re-pins.

---

## WP2 — material-resolved mobilization thresholds + `"rigid"` class

### Problem

The docs lean on Mellink et al. (2024) — near-100% mobility for bags vs
near-0% for bottles at the same ~2.3 m/s wind, "material identity dominates"
— yet `material` currently has **no effect on dry mobilization**: film and
paper share one `gust_threshold`, and rigid items aren't representable at
all. The model's headline citation and its structure disagree.

### Design

Resolve `gust_threshold` / `gust_reference` **per material** when the caller
does not supply them:

- Change formals to `gust_threshold = NULL, gust_reference = NULL`. When
  `NULL`, resolve from `LITTER_CONSTANTS$MATERIALS[[material]]`; when
  supplied, the explicit number wins (calibration override). Run the existing
  `gust_reference > gust_threshold` validation *after* resolution.
- Extend the material table (and `material = c("film", "paper", "rigid")`):

| material | gust_threshold (m/s) | gust_reference (m/s) | saturation_penalty | rationale |
|---|---|---|---|---|
| film  | 3.9737355 | 13.9080743 | 0.7 | current defaults, behaviour-preserving |
| paper | 3.9737355 | 13.9080743 | 1.0 | dry paper is as mobile as film (the Mellink film/paper contrast is about water absorption, not dry mobility); keeping it equal to film preserves the pinned film==paper unsaturated oracle |
| rigid | 12.0 | 25.0 | 0.15 | bottles ~0% mobile at winds that fully mobilize bags (Mellink 2024); a wet surface barely holds a rigid item down. **Uncalibrated placeholders** |

- Docs: rewrite the `@param gust_threshold`/`@param gust_reference` text to
  describe NULL-resolution (keep the `0.30/(0.40/ln 200)` derivation note,
  scoped to film/paper), and fix the Model-section claim so it now truthfully
  says material enters both the dry threshold and the wet-saturation
  treatment. Mention `rigid` in `litter_hazard()` / `litter_risk()` docs
  (both forward `material` via `...`; `match.arg` picks up the new level
  automatically once the formal's default vector is extended).

### Tests to write first

```r
describe("material-resolved mobilization thresholds (v3.2)", {
  it("rigid needs a much stronger gust than film: G=10 moves film, not rigid", {
    expect_gt(litter_hazard_vec(10, 12, 0, 0.02, material = "film"), 0)
    expect_equal(litter_hazard_vec(10, 12, 0, 0.02, material = "rigid"), 0)
  })
  it("rigid worked oracle: G=16, W=15, dry -> 9.4675", {
    # threshold 12, reference 25, denom 13; E = 50*((16-12)/13)^2 = 4.7337; T=2.
    out <- litter_hazard_vec(16, 15, 0, 0.02, material = "rigid")
    expect_equal(out, 9.4675, tolerance = 1e-3)
  })
  it("rigid saturates at the rigid reference: G=26 and G=40 both -> 100", {
    out <- litter_hazard_vec(c(26, 40), c(15, 15), c(0, 0), c(0.02, 0.02),
                             material = "rigid")
    expect_equal(out, c(100, 100), tolerance = 1e-3)
  })
  it("an explicit gust_threshold overrides the material default", {
    a <- litter_hazard_vec(16, 15, 0, 0.02, material = "rigid",
                           gust_threshold = 3.9737355, gust_reference = 13.9080743)
    b <- litter_hazard_vec(16, 15, 0, 0.02, material = "film")
    expect_equal(a, b, tolerance = 1e-3)
  })
  it("film and paper share the dry threshold (pinned film==paper oracle intact)", {
    expect_equal(litter_hazard_vec(10, 7, 0, 0.10, material = "paper"),
                 1.2127, tolerance = 1e-3)
  })
})
```

Also confirm (they already exist — just run them): every default-material
pinned oracle in rule 3, since `material = "film"` resolution must reproduce
the old hard-coded defaults bit-for-bit.

**Done when:** new tests pass; full litter filter passes with zero re-pins;
`man/` regenerated.

---

## WP3 — consistent gap-bearing destination semantics

### Problem

When the downwind bearing falls in a sector gap, `exposure` is scaled by
`default_permeability` (litter passes the unconfigured boundary) but
`leaves_site` is hard-`FALSE` (litter never passes). The magnitude and the
flag are computed from contradictory assumptions.

### Design

Compute the gap fall-through destination with the **same** assumptions the
magnitude already uses: permeability = `default_permeability`, sensitivity =
`FALSE`, distance = unknown.

In `litter_exposure()`, after the source/barrier loop (`best_perm` is `-Inf`
exactly on no-hit hours):

```r
gap_hours   <- !is.finite(best_perm)
gap_open    <- default_permeability >= p_open_min
gap_reached <- if (refined) FALSE else hazard >= offsite_threshold
leaves_site <- leaves_site | (gap_hours & gap_open & gap_reached)
```

- `sensitive_receptor` stays `FALSE` for gaps, so `zone` can still never be
  `off_site` through a gap — zone semantics unchanged.
- Refined mode has no distance for a gap, so `gap_reached = FALSE` there;
  document this limitation explicitly in the `@param default_permeability` /
  Method roxygen.
- Reword the `meteoHazard_litter_no_barriers` warning body (keep the class):
  "no hour can be off-site" remains true, but the second bullet should now
  say gaps use `default_permeability` for both exposure and `leaves_site`,
  and are never sensitive.

### Tests to write first / to re-pin (this WP intentionally re-pins two expectations)

1. `test-litter-exposure.R`, test *"applies the default permeability when the
   downwind bearing hits no barrier"* (hazard 86, wind 0, default 0.5): change
   `expect_false(out$leaves_site)` → `expect_true(...)` (86 ≥ 45 and 0.5 ≥
   `p_open_min` 0.5) with a comment citing this WP. `sensitive_receptor`
   stays `FALSE`; exposure stays 43.
2. Zero-barrier warning test: add `expect_true(out$leaves_site)` (hazard 90)
   and keep `zone == "on_site"`.
3. New cases:

```r
it("a gap bearing does NOT leave site when default_permeability < p_open_min", {
  site <- .make_demo_mh_site()
  out <- litter_exposure(86, 0, site, default_permeability = 0.4)  # p_open_min 0.5
  expect_false(out$leaves_site)
})
it("a gap bearing does NOT leave site below offsite_threshold", {
  site <- .make_demo_mh_site()
  out <- litter_exposure(30, 0, site)   # 30 < 45
  expect_false(out$leaves_site)
})
it("refined mode: gap bearings never satisfy the reach test (no distance)", {
  site <- .make_demo_mh_site()
  out <- litter_exposure(86, 0, site, mean_wind = 50, reach_per_ms = 100)
  expect_false(out$leaves_site)
})
```

**Done when:** the two re-pins and three new tests pass; no other exposure
test changes.

---

## WP4 — robust `.barrier_bearing_range()` (boundary densification)

### Problem

The largest-gap heuristic sees only polygon **vertices**. A coarse concave
barrier can put a bigger gap *between consecutive vertex bearings inside the
true arc* than the real outside gap, returning the complement. Concrete
failure: a square-cornered C-shape (below) whose true subtended arc is ~323°;
the heuristic returns a 270° arc and reports a bearing pointing straight into
a wall as a miss. (The v3.1-review wide-arc fix in `.litter_arc_contains()`
does not save this: 270 + 2·15 < 360.)

### Design

Densify the boundary before taking bearings, so consecutive sampled bearings
are at most a few degrees apart and the largest gap can only be the true
outside gap. In `.barrier_bearing_range()`, after the enclosure guard:

```r
src_xy <- sf::st_coordinates(source_pt)[1, c("X", "Y")]
verts0 <- sf::st_coordinates(barrier_poly)[, c("X", "Y")]
d2     <- (verts0[, "X"] - src_xy["X"])^2 + (verts0[, "Y"] - src_xy["Y"])^2
d_min  <- sqrt(min(d2[d2 > 0]))          # nearest non-coincident vertex
seg_len <- max(d_min * 0.05, 0.5)         # <= ~2.9 deg angular step at range d_min
poly_dense <- sf::st_segmentize(sf::st_geometry(barrier_poly), dfMaxLength = seg_len)
verts <- sf::st_coordinates(poly_dense)[, c("X", "Y")]
```

then proceed with the existing coincident-vertex exclusion, dedup, and
largest-gap logic on `verts`. Notes:

- Keep the exclusion of exactly-coincident vertices — segmentize can land a
  sample exactly on the source (the multi-source test's `barrier_open` edge
  passes through source A; its midpoint sample coincides and must be dropped,
  which the existing `dist2 > 0` filter already does).
- `d_min * 0.05` bounds the angular sampling error at ≈ 2.9°, far inside
  `direction_tol` (15°). The `0.5` m floor bounds the vertex count when the
  source is very close to the boundary.
- Wedges from `site_from_sectors()` already have dense arcs; adding points
  along their straight radial edges introduces no new extreme bearings
  (bearings along a straight segment vary monotonically between the endpoint
  bearings), so existing results are unchanged.

### Tests to write first (add to the "wide barrier arcs" describe block in `test-litter-exposure.R`)

```r
it("a coarse square-cornered C-shape is hit through its walls and missed through its opening", {
  cx <- 335000; cy <- 6250000
  # C-shape opening due NORTH; source sits in the cavity. 8 corners only.
  ring <- matrix(c(
    cx - 300, cy - 300,   cx + 300, cy - 300,   cx + 300, cy + 300,
    cx + 100, cy + 300,   cx + 100, cy - 100,   cx - 100, cy - 100,
    cx - 100, cy + 300,   cx - 300, cy + 300,   cx - 300, cy - 300
  ), ncol = 2, byrow = TRUE)
  feats <- sf::st_sf(
    id       = c("source", "cshape"),
    geometry = sf::st_sfc(sf::st_point(c(cx, cy)),
                          sf::st_polygon(list(ring)), crs = 32755)
  )
  roles <- data.frame(
    feature_id = c("source", "cshape"), hazard = "litter",
    role = c("source", "barrier"),
    permeability = c(NA_real_, 1.0), sensitive = c(NA, TRUE),
    stringsAsFactors = FALSE
  )
  site <- mh_site(feats, roles, epsg = 32755)
  # Due EAST (wind 270 -> theta_down 90) points straight at the x = +100 wall:
  east <- litter_exposure(50, 270, site, default_permeability = 0)
  expect_equal(east$exposure, 50, tolerance = 1e-3)
  # South and west walls likewise:
  expect_equal(litter_exposure(50, 0,  site, default_permeability = 0)$exposure, 50, tolerance = 1e-3)
  expect_equal(litter_exposure(50, 90, site, default_permeability = 0)$exposure, 50, tolerance = 1e-3)
  # Due NORTH (wind 180 -> theta_down 0) exits through the opening
  # (inner mouth corners subtend ~[18.4, 341.6]; 0 is outside even with tol 15):
  north <- litter_exposure(50, 180, site, default_permeability = 0)
  expect_equal(north$exposure, 0, tolerance = 1e-6)
})
```

Run this test before implementing to confirm it fails on the east/south/west
assertions (TDD-red), then implement.

**Done when:** the C-shape test passes; the existing horseshoe, offset-source,
multi-source, and demo-site regressions all pass unchanged.

---

## WP5 — document the transport/reach calibration coupling (docs only)

### Problem

In refined mode, mean wind drives the hazard's transport multiplier `T` *and*
the exposure's reach test `c_L · mean_wind >= distance_m`. That is defensible
— `T` amplifies the mobilized **flux**, the reach test is a geometric
**gate** — but it means `transport_max`/`wind_transport_onset/ref` and
`reach_per_ms` cannot be calibrated independently, and nothing currently says
so.

### Change

No code. Add:

1. `litter_hazard_vec()` roxygen, transport bullet: one sentence — in
   combined use with `litter_exposure()`'s refined mode, `T` acts as a flux
   amplifier while the reach test gates the destination; calibrate the
   transport ramp and `reach_per_ms` **jointly** (issues #11/#26).
2. `litter_exposure()` roxygen, `@param reach_per_ms`: mirror sentence.
3. `litter_risk()` roxygen: pointer to both.
4. NEWS.md bullet.

**Done when:** `devtools::document()` output is clean and the wording appears
in the three `.Rd` files.

---

## Final checklist (after all WPs)

- [ ] `devtools::test()` fully green; the only intentional expectation changes
      are the two WP3 re-pins.
- [ ] `devtools::document()` run; `man/` diff reviewed.
- [ ] Every pinned oracle in Ground-rule 3 asserted unchanged.
- [ ] NEWS.md has one bullet per WP.
- [ ] All new numeric defaults appear in `LITTER_CONSTANTS` with
      "uncalibrated placeholder" roxygen flags — no new inline literals.
