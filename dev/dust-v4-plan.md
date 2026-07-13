# Dust v4 — implementation plan

Follow-up to the scientific review of the dust v3 PR (#28, merged as
`7fe349b`). The review confirmed the v3 physics chain (Shao & Lu threshold,
Fecan moisture correction, Owen/White saltation, MB95 sandblasting) is
implemented faithfully to the cited sources — both known-answer flux anchors
were reproduced independently to 8 significant figures — and fixed the minor
consistency issues directly (commit `b5a67a5`). This plan operationalises the
**major** items: the idealisations v3 itself declared as TODOs (T8/T10), the
smooth-surface/no-drag-partition gap, and the clock-driven crust decay.

Execute the WPs **in order**. WP2 intentionally re-pins two known-answer
fixtures (called out explicitly there); nothing else may move a pinned value.
Each WP is a separate conventional commit (`feat(dust): ...`).

## Ground rules

1. Read first: `CLAUDE.md`, `R/dust_hazard.R`, `R/constants.R`
   (`DUST_CONSTANTS`), `tests/testthat/test-dust-hazard.R`, and for WP5
   `R/site-model.R`, `R/odour_exposure.R`, `R/odour-ventilation.R`.
2. **Write tests before implementation** (repo rule), in the existing
   `describe()/it()` style with `skip_if_no_dust_v2()` gating.
3. **Pinned oracles.** Verified analytically during the review; only WP2 may
   re-pin the two marked *:
   - `dust_flux(20L, 10, 0, 20, 0.02, z0 = 0.005)` = **1.35399760e-06** *
   - same at `clay = 50` (capped) = **2.96222398e-05** *
   - smooth-bed (z0 = NULL) threshold: gust 17 → **0**, gust 18 → **> 0**
   - z0 = 0.005 threshold: gust 10.5 → **0**, gust 10.75 → **> 0** *
   - `dust_hazard(met, crust = FALSE)` ≡ `dust_flux()` on the same columns
   - crust cold-start behaviours (CC-4a, CC-4c)
4. New tunables go into `DUST_CONSTANTS` (`R/constants.R`) with roxygen
   `\item{}` entries; no inline literals. Every new default is an
   **uncalibrated placeholder** — say so in its roxygen.
5. Dust conditions are classed (since `b5a67a5`): errors
   `meteoHazard_input_error`, warnings `meteoHazard_dust_*`. New conditions
   follow that pattern.
6. After each WP: `devtools::document()`, `devtools::test(filter = "dust")`;
   before the final push: full `devtools::test()`. NEWS.md bullet per WP.
7. **Not planned, stays documented-only:** supply limitation / surface
   armouring (research-grade; the Idealisations section already records it),
   and any recalibration of `gust_factor`.

---

## WP1 — met-driven air density + direct grain-size interface (v3's T10)

### Problem

`rho_a` is fixed at 1.225 kg/m³. Air density enters the threshold
(u*t ∝ ρ_a^-1/2) *and* the flux prefactor (Q ∝ ρ_a), so cold dense air both
lowers the entrainment threshold and carries more sand: at −5 °C / 1013 hPa
(ρ = 1.3164 kg/m³) the same 20 m/s gust yields **1.375×** the reference-density
flux. Winter mornings are systematically under-predicted. Separately, sites
with a measured d50 must reverse-engineer a Tyler sieve number.

### Design

**(a) Air density.** New optional vector args on `dust_flux()`:
`temperature_2m = NULL` (degC, plain numeric — matches `litter_wetness_vec()`)
and `surface_pressure = NULL` (hPa, through `.drop_to(x, "hPa")`; Open-Meteo
`surface_pressure`). Both supplied → per-hour
`rho_a = 100 * surface_pressure / (R_D * (temperature_2m + 273.15))` (the
`* 100` is hPa→Pa); neither → `DUST_CONSTANTS$RHO_A_REF` (behaviour
preserving); exactly one → classed `meteoHazard_input_error`. Add
`R_D = 287.05` (J kg⁻¹ K⁻¹, dry-air gas constant) and
`KELVIN_OFFSET = 273.15` to `DUST_CONSTANTS`. Sanity-assert the computed
density lies in `[0.8, 1.6]` kg/m³ (classed error otherwise — catches
Pa-vs-hPa mistakes). `rho_a` becomes a length-n vector through the chain
(`u_star_t_dry`, `Q`); verify every use vectorises.

On `dust_hazard()`: new `air_density = c("reference", "met")`. `"met"`
requires the `temperature_2m` and `surface_pressure` columns (via
`.assert_required_cols()`) and forwards them; `"reference"` (default) ignores
them. Do NOT auto-detect columns — silent behaviour changes are worse than an
explicit flag.

**(b) d50 interface.** New `d50 = NULL` arg on `dust_flux()` and
`dust_hazard()` (m; through `.drop_to(x, "m")`; assert in `[1e-5, 0.02]`).
Non-NULL `d50` supersedes the sieve table; if the caller *explicitly* passed
`tyler_sieve_no` too (`!missing(tyler_sieve_no)`), warn
(`meteoHazard_dust_d50_supersedes`) — mirroring the litter
`soil_moisture`/`wetness` seam.

### Tests first (worked oracles, verified during planning)

```r
it("met-driven air density: -5 degC / 1013.25 hPa raises the gust-20 flux by ~1.375x", {
  base <- dust_flux(20L, 10, 0, 20, 0.02)
  cold <- dust_flux(20L, 10, 0, 20, 0.02,
                    temperature_2m = -5, surface_pressure = 1013.25)
  expect_equal(base, 7.86633315e-08, tolerance = 1e-4)   # new smooth-bed KAT
  expect_equal(cold, 1.08192029e-07, tolerance = 1e-4)
  expect_equal(cold / base, 1.37538, tolerance = 1e-3)
})
it("supplying only one of temperature_2m/surface_pressure errors (classed)", {
  expect_error(dust_flux(20L, 10, 0, 20, 0.02, temperature_2m = -5),
               class = "meteoHazard_input_error")
})
it("d50 equal to the sieve-20 opening reproduces the sieve-20 flux exactly", {
  expect_equal(dust_flux(20L, 10, 0, 20, 0.02),
               dust_flux(clay_percent = 10, wind_speed_10m = 0,
                         wind_gusts_10m = 20, soil_moisture = 0.02,
                         d50 = 0.000833))
})
```

Plus: `dust_hazard(air_density = "met")` column requirement; default path
bit-identical to current (all Ground-rule-3 oracles).

**Done when:** new tests pass; every pinned oracle unchanged; `man/` synced.

---

## WP2 — MB95 drag partition (make a caller-supplied z0 physically meaningful)

### Problem

v3 removed the backwards roughness behaviour by defaulting to a smooth bed
and *warning* on rougher z0 — but a warned-and-biased-high path is still the
only way to represent a rough site. MB95's own drag partition is the
literature answer: roughness absorbs part of the shear stress, so only a
fraction reaches the erodible bed.

### Design

Efficient-fraction ratio (Marticorena & Bergametti 1995, Eqs. 18–19;
**lengths in cm**, internal-boundary-layer height X = 10 cm):

```
feff(z0, z0s) = 1 - ln(z0/z0s) / ln(MB95_DP_COEF * (MB95_DP_X_CM/z0s_cm)^MB95_DP_EXP)
```

with new constants `MB95_DP_COEF = 0.7`, `MB95_DP_EXP = 0.8`,
`MB95_DP_X_CM = 10`. Apply as a threshold divisor: `u_star_t <- u_star_t / feff`
(u* itself keeps using the caller z0 in the log law, as now).

- `z0 <= z0_smooth` (including the NULL default) → `feff = 1` exactly:
  smooth-bed behaviour bit-identical.
- `feff <= DUST_CONSTANTS$FEFF_MIN` (new, `0.01`) → the surface is **fully
  sheltered**: return all-zero flux with classed warning
  `meteoHazard_dust_fully_sheltered`. (Verified: feff goes *negative* at
  z0 ≳ 2 cm over a 28 µm bed — the formula is a fit, not global physics.)
- **Delete** the `meteoHazard_dust_z0_rough` warning and its test — rough z0
  is now modelled, not warned about.

Reference numbers (verified during planning, sieve 20, dry, defaults):

| z0 (m) | feff | threshold gust at 10 m (m/s) |
|---|---|---|
| 2.78e-5 (= d/30) | 1.0000 | 17.90 |
| 1e-4 | 0.7932 | 20.31 |
| 5e-4 | 0.5333 | 25.98 |
| 1e-3 | 0.4214 | 30.58 |
| 5e-3 | 0.1616 | 65.80 |
| 5e-2 | −0.21 → fully sheltered | ∞ |

### Sanctioned re-pins (the two *-marked oracles)

The z0 = 0.005 KAT fixtures were pinned to the *no-partition* model, in which
z0 = 0.005 raised u* without sheltering; under the partition the same fixture
is sub-threshold at gust 20 (threshold ≈ 65.8 m/s) and the old numbers are
unreproducible **by design**. Re-point the fixtures:

- CC-5a / CC-2a: change the fixture to `z0 = NULL` (smooth bed). New pinned
  values: **7.86633315e-08** (clay 10) and **1.72096e-06** (clay 50, capped —
  recompute to 6 s.f. as `7.86633315e-08 * 10^(0.134*10)` and pin what you
  compute).
- CC-5b: replace with the partition threshold at z0 = 5e-4:
  gust **25.5 → 0**, gust **26.5 → > 0**.
- CC-1b (z0 warning): replace with the fully-sheltered test:
  `dust_flux(..., z0 = 0.05)` returns all zeros and warns
  `meteoHazard_dust_fully_sheltered`.

New property test: the threshold 10-m gust is strictly increasing across
`z0 = c(2.78e-5, 1e-4, 5e-4, 1e-3, 5e-3)` (bisect the zero/nonzero boundary
or assert flux at a fixed 26 m/s gust is non-increasing in z0 over the first
four entries — do NOT assert monotonicity of supra-threshold flux magnitude
in general; the u*³ prefactor and feff compete). Document in the roxygen
Idealisations that feff is highly sensitive to z0s = d/30 and is an
uncalibrated screening treatment; screening users should keep `z0 = NULL`.

**Done when:** partition tests + re-pins pass; smooth-bed path bit-identical
(all non-starred oracles); `man/` synced.

---

## WP3 — within-hour Weibull intermittency (v3's T8)

### Problem

The steady hourly-gust forcing biases the flux high near threshold (documented
idealisation) and — the mirror image — reports exactly zero for hours whose
*mean* is sub-threshold but whose within-hour distribution tail is not.

### Design

New `forcing = c("gust", "weibull")` on `dust_flux()` and `dust_hazard()`
(default `"gust"`, behaviour preserving). In `"weibull"` mode the within-hour
10-m wind is `U ~ Weibull(k, c)` with shape `k = DUST_CONSTANTS$WEIBULL_SHAPE`
(new, `2.0`, uncalibrated placeholder) and scale `c = wind_speed_10m /
gamma(1 + 1/k)` (so the mean is preserved); `wind_gusts_10m` and
`gust_factor` are unused on this path (document this).

The hourly-expected flux has a closed form. With `a = kappa / log(z/z0)`
(so `u* = a·U`), threshold wind `Ut = u_star_t / a`, and the truncated
Weibull moments

```
E[U^m; U > Ut] = c^m * gamma(1 + m/k) *
                 pgamma((Ut/c)^k, shape = 1 + m/k, rate = 1, lower.tail = FALSE)
```

the expectation of `Q = (rho_a/g) * (a U)^3 * (1 - (u*t/(a U))^2)` over the
active tail is

```
E[Q] = (rho_a/g) * (a^3 * E[U^3; U>Ut]  -  a * u_star_t^2 * E[U; U>Ut])
```

then `F = alpha * E[Q]` as usual. All vectorised; `u_star_t` (moisture/crust
adjusted, WP2 feff included) is per-hour, so `Ut` is per-hour.

### Tests first (worked oracles, verified against numeric integration to 9 s.f. during planning)

```r
it("weibull forcing matches the closed-form oracle (mean 12, k = 2, dry, smooth bed)", {
  f <- dust_flux(20L, 10, wind_speed_10m = 12, wind_gusts_10m = 20,
                 soil_moisture = 0.02, forcing = "weibull")
  expect_equal(f, 1.062634e-07, tolerance = 1e-4)
})
it("weibull forcing emits where the steady mean is sub-threshold (intermittency tail)", {
  steady <- dust_flux(20L, 10, 12, 12 / 0.84, 0.02)          # U_fm = 12 -> 0
  weib   <- dust_flux(20L, 10, 12, 12 / 0.84, 0.02, forcing = "weibull")
  expect_equal(steady, 0)
  expect_gt(weib, 0)
})
it("a near-degenerate Weibull (k large) converges to the steady flux at the mean", {
  # k = 200: c ~ mean; compare against gust-mode with U_fm forced equal to the mean.
  # (Set the constant via a with_mocked/testthat local constant override, or expose
  #  weibull_shape as an argument defaulting to DUST_CONSTANTS$WEIBULL_SHAPE.)
})
```

Expose `weibull_shape = DUST_CONSTANTS$WEIBULL_SHAPE` as a formal (needed for
the convergence test and for calibration). Also assert: `forcing = "gust"`
path bit-identical to current (all pinned oracles).

**Done when:** oracle + tail + convergence tests pass; gust path unchanged.

---

## WP4 — saltation-gated crust decay

### Problem

Crust decay is a pure clock (`exp(-age/72h)`), but real crust breakdown is
mechanical — abrasion by saltating grains (Gillette 1982; Rice & McEwan
2001), already flagged in v3's Idealisations. Under a week of calm weather
the clock erodes the crust while nothing is hitting it.

### Design

New `crust_decay = c("clock", "saltation")` on `dust_hazard()` (default
`"clock"`, behaviour preserving). In `"saltation"` mode the crust age
advances only during hours when saltation is actually occurring *on the
crusted surface*:

1. **Refactor first (pure, no behaviour change):** extract the u* chain from
   `dust_flux()` into internal helpers `.dust_u_star(wind_speed_10m,
   wind_gusts_10m, gust_factor, z0)` and `.dust_u_star_t_dry(d, rho_a)`, and
   have `dust_flux()` call them. Run the full dust filter — every oracle must
   be bit-identical — and commit the refactor separately.
2. New `.dust_crust_factor_saltation(precipitation, u_star, u_star_t_moist,
   threshold, factor_max, decay_hours, age0)`: sequential loop over hours;
   at hour i, compute `mult_i` from the *current* age, form the crusted
   threshold `u_star_t_i = u_star_t_moist_i * max(1, mult_i)`
   (note: combine with the moisture factor via max, as `dust_flux()` does);
   then advance: rain hour (`precip >= threshold`) → `age <- 0`; else if
   `u_star_i > u_star_t_i` → `age <- age + 1` (crust is being sandblasted);
   else age holds (calm hours do not weather the crust). Return the per-hour
   multiplier vector, which `dust_hazard()` injects as
   `threshold_multiplier` exactly as the clock path does.
3. `dust_hazard()` computes `u_star`/`u_star_t_moist` for the gate via the
   step-1 helpers using the same site parameters it forwards to
   `dust_flux()`. (Yes, u* is computed twice; correctness over micro-perf.)

### Tests first

```r
it("saltation mode: crust persists through a calm week, clock mode decays it", {
  # rain hour 1, then 167 dead-calm hours, then a strong hour 169
  met <- data.frame(wind_speed_10m = c(0, rep(0, 167), 0),
                    wind_gusts_10m = c(0, rep(2, 167), 25),
                    soil_moisture_0_to_1cm = 0.02,
                    precipitation = c(5, rep(0, 168)))
  clock <- dust_hazard(met, crust = TRUE, crust_decay_hours = 24)
  salt  <- dust_hazard(met, crust = TRUE, crust_decay_hours = 24,
                       crust_decay = "saltation")
  expect_lt(salt[169], clock[169])   # crust still ~fresh in saltation mode
})
it("saltation mode: sustained supra-threshold wind decays the crust", { ... })
it("crust_decay = 'clock' is bit-identical to the current behaviour", { ... })
```

**Done when:** step-1 refactor lands with zero numeric change; saltation
tests pass; clock path bit-identical.

---

## WP5 — dust_exposure() / dust_risk() (v3's T9; design outline, own PR)

Largest item; this WP is an architecture contract, not line-level
instructions. Mirror the odour split exactly:

- **`dust_exposure(met_data, site, ...)`**: for each `(dust, source)` feature
  in an `mh_site`, compute the emission flux via the WP1–WP4 chain (surface
  parameters `tyler_sieve_no`/`d50`, `clay_percent`, `z0`, `bulk_density`
  carried as role attributes — first check what extra role columns
  `mh_site()` accepts, see `R/site-model.R`, and extend its validation the
  way `emit_height`/`permeability` are handled); disperse to
  `(dust, receptor)` features with the ISC3 area-source Gaussian machinery
  already in `R/odour_exposure.R`, driven by `ventilation_state()`
  (`u_eff`, `h_mix`, stability class). Output shape mirrors
  `odour_exposure()` (per-hour, per-receptor).
- **`dust_risk(met_data, site, ...)`**: thin wrapper, mirrors
  `odour_risk()`/`litter_risk()`.
- Reuse, don't copy: any plume helper needed by both odour and dust moves to
  a shared internal file (the precedent is `R/geometry.R` for the litter/odour
  wind helpers, R-D3).
- Tests: contract tests mirroring `test-odour-exposure.R`'s structure; at
  minimum a worked single-source/single-receptor oracle, a direction
  sensitivity test, and hazard/exposure consistency (calm ⇒ both zero).

**Done when:** a design note + failing contract-test skeleton exist and are
agreed; implementation proceeds as its own PR against that skeleton.

---

## Final checklist

- [ ] `devtools::test()` fully green; the only re-pins are WP2's sanctioned
      ones.
- [ ] `devtools::document()` run after each WP; `man/` diffs reviewed.
- [ ] All new constants in `DUST_CONSTANTS` with roxygen entries and
      uncalibrated-placeholder flags.
- [ ] NEWS.md: one bullet per WP.
- [ ] Idealisations sections updated: WP1 removes the fixed-density bullet,
      WP2 rewrites the smooth-surface bullet, WP3 rewrites the steady-gust
      bullet, WP4 rewrites the crust-clock paragraph. Supply limitation
      remains, explicitly, the one documented-not-modelled item.
