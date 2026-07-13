# meteoHazard (development version)

## Dust hazard v4

Follow-up to the dust v3 review: operationalises the T8/T10 idealisation
TODOs, the smooth-surface/no-drag-partition gap, and clock-driven crust decay.
`dust_exposure()`/`dust_risk()` (v3's T9) remain out of scope for this PR —
deferred to their own follow-up against a design-note/contract-test skeleton.

* **Met-driven air density + direct grain-size interface** (WP1): new optional
  `temperature_2m`/`surface_pressure` args on `dust_flux()` (and
  `air_density = c("reference", "met")` on `dust_hazard()`) compute a per-hour
  air density via the ideal gas law instead of the fixed
  `DUST_CONSTANTS$RHO_A_REF`; cold, dense air both lowers the entrainment
  threshold and raises the flux (-5 degC / 1013.25 hPa raises a 20 m/s-gust
  flux ~1.375x). New `d50` arg (m) bypasses the Tyler sieve table for sites
  with a measured grain size; `tyler_sieve_no` is now `NULL`-default and
  exactly one of `tyler_sieve_no`/`d50` is required (supplying both warns,
  `d50` wins). New constants `R_D`, `KELVIN_OFFSET`.
* **MB95 drag partition** (WP2): a caller-supplied `z0` above the smooth-bed
  value now engages the Marticorena & Bergametti (1995, Eqs 18-19)
  efficient-fraction drag partition -- roughness shelters the erodible bed by
  raising the entrainment threshold (`u_star_t / feff`) instead of just
  raising u* unsheltered. `z0 <= z0_smooth` (including the `NULL` default)
  still gives `feff = 1` exactly (bit-identical smooth-bed path). A `z0` large
  enough to fully shelter the bed (`feff <= DUST_CONSTANTS$FEFF_MIN`) now
  returns all-zero flux with a new classed warning
  (`meteoHazard_dust_fully_sheltered`), replacing the removed
  `meteoHazard_dust_z0_rough` warning. New constants `MB95_DP_COEF`,
  `MB95_DP_EXP`, `MB95_DP_X_CM`, `FEFF_MIN` (uncalibrated screening use).
  **Breaking (re-pin):** the two `z0 = 0.005` known-answer fixtures move to
  the smooth bed (`7.86633315e-08` clay-10, `1.72096e-06` clay-50 capped) and
  the `z0 = 5e-4` partition threshold (gust 25.5 -> 0, 26.5 -> positive)
  respectively -- the old fixtures assumed no partition and are
  unreproducible by design once roughness shelters the bed.
* **Within-hour Weibull intermittency** (WP3): new
  `forcing = c("gust", "weibull")` on `dust_flux()`/`dust_hazard()` (default
  `"gust"`, bit-identical). `"weibull"` replaces the steady hourly-gust
  forcing with a within-hour `Weibull(weibull_shape, c)` wind distribution
  (mean-preserving on `wind_speed_10m`; `wind_gusts_10m`/`gust_factor` unused)
  and a closed-form hourly-expected flux, correcting the steady-forcing bias
  near threshold and emitting where the hourly mean is sub-threshold but the
  within-hour tail is not. New `weibull_shape` argument (default
  `DUST_CONSTANTS$WEIBULL_SHAPE`, 2.0, uncalibrated placeholder).
* **Saltation-gated crust decay** (WP4): new `crust_decay = c("clock",
  "saltation")` on `dust_hazard()` (default `"clock"`, bit-identical). The
  pure elapsed-time clock decays a crust through a calm week exactly as fast
  as a windy one; `"saltation"` instead advances the crust age only during
  hours when the (currently crusted) surface is actually being sandblasted
  (rain resets it, calm/sub-threshold hours leave it untouched), so it
  persists through a calm spell and decays once winds resume. Internally
  refactors the u* chain into `.dust_u_star()`/`.dust_u_star_t_dry()` (pure,
  no-behaviour-change) and adds `.dust_crust_factor_saltation()`; `u*` is
  computed twice for the `"saltation"` gate (once in the gate, once in the
  real flux) -- correctness over micro-perf.

## Litter hazard v3.2

Follow-up to the v3.1 review: four design corrections plus a documentation item.

* **Smooth material-saturation ramp** (WP1): the wet penalty now ramps linearly
  from `saturation_onset` (0.8 normalised wetness) to full at wetness 1, removing
  the interior discontinuity at the saturation mark. Shared by the soil-moisture
  and wetness-state paths.
* **Material-resolved mobilization thresholds** (WP2): `gust_threshold`,
  `gust_reference`, and `saturation_penalty` now default per `material` from the
  new internal `LITTER_CONSTANTS$MATERIALS` (`NULL` formals; explicit values
  still override). A new `material = "rigid"` (bottles/cans) makes `material`
  affect dry mobility, matching the Mellink (2024) bags-vs-bottles evidence; its
  numbers are uncalibrated placeholders. `paper_veto_wetness` is replaced by
  `saturation_onset`.
* **Consistent gap-bearing semantics** (WP3): a downwind bearing in a sector gap
  now derives `leaves_site` from the same `default_permeability` the exposure
  magnitude uses (previously the boundary was treated as both passable and not).
  `sensitive_receptor` stays `FALSE` for gaps, so a gap can never be off-site.
* **Robust concave barriers** (WP4): `litter_exposure()` densifies each barrier
  boundary before extracting bearings, fixing a largest-gap inversion that could
  report a bearing pointing straight at a coarse concave (C-shaped) barrier as a
  miss.
* **Transport/reach coupling documented** (WP5): the refined-mode reach test and
  the hazard transport multiplier are both driven by the mean wind and must be
  calibrated jointly; noted in `litter_hazard_vec()`, `litter_exposure()`, and
  `litter_risk()`.

# meteoHazard 0.3.0

## Odour hazard physics corrections (v3)

Five physical corrections to the odour hazard/exposure/ventilation chain,
following a scientific review of the model.

* **Nocturnal cold-pool cap on `h_mix`** (the biggest fix): on stable nights
  the dilution depth is now `min(boundary_layer_height, pool_top)` instead of
  always the raw forecast boundary-layer height, so the already-computed
  thermal cold-pool depth actually caps ventilation. This corrects the known
  worst case at the primary site — a valley landfill with ridge-top receptors
  under winter temperature inversions. New `pool_cap = TRUE` argument on
  `ventilation_state()`, `odour_hazard()`, `odour_exposure()`, and
  `odour_risk()`; set `FALSE` to restore the old raw-BLH behaviour. Gated to
  night-only hours where the pool is thermally derived from temperature + RH.
* **Breaking (numeric):** the generation modifier `G` is now **multiplicative**
  rather than additive (`(1+dP_mod)*(1+R_mod)*(1+S_seal)*(1+H_mod)*(1+V_mod)`),
  correctly compounding independent fractional emission effects instead of
  under-predicting coincident worst cases. New range ≈ [0.80, 2.33] (was
  [0.80, 1.95]).
* **Breaking (numeric):** `V_mod` (surface volatilisation) is now
  **exponential** in temperature rather than linear, matching the Henry's-law
  / Clausius-Clapeyron basis (vapour pressure roughly doubles every 10 °C).
  Still anchored at 0 at 10 °C and `V_MOD_MAX` (0.30) at 35 °C.
* Rain scavenging `W_rain` is now **solubility-aware**: new
  `odorant_solubility` argument (numeric in `[0, 1]`, default 0.5) on the same
  four functions as `pool_cap`. `0` models poorly soluble reduced-sulfur
  odorants (barely washed out by rain); `1` reproduces the old fixed washout
  tiers and models highly soluble species (ammonia, amines); the default 0.5
  is a mixed profile.
* Removed the internal `.odour_dispersion_state()` backward-compatibility
  shim (pre-v0, unexported, no users); everything now goes through
  `ventilation_state()` directly.

## Hazards expose physical quantities; the fixed 0–100 scale is removed (issue #11)

The package no longer ships a fixed 0–100 operational scale or tiers for odour,
dust, or litter. Each hazard now returns a physical / relative quantity;
mapping it onto a site-specific operational index (potentially 0–100) and tiers
is a calibration step to be delivered by forthcoming calibration tooling
(issues #11/#8). TWL keeps its tiers (`categorise_twl()`/`twl_colour()`) because
its W/m² zones are physically grounded (Brake & Bates 2002).

* **Breaking:** `odour_exposure()` and `odour_risk()` return the unbounded
  **per-receptor relative concentration** as a `n_hours x n_receptors` matrix
  (column names are the receptor `id`s), instead of a worst-case 0–100 vector.
  Reduce over receptors yourself (e.g. `apply(out, 1, max)`). The `map_c50`
  argument is removed.
* **Breaking:** `dust_hazard()` now returns the crust-adjusted **relative dust
  flux** (same units as `dust_flux()`), not a 0–100 index. The `scale_ref_gust`
  argument and the reference normalisation are removed; the crust-persistence
  gate and `met_data` interface are retained.
* **Breaking:** `litter_hazard()` / `litter_hazard_vec()` drop the `min(., 100)`
  cap and return the relative entrainment × transport index directly (numerically
  unchanged — the default `entrainment_max`/`transport_max` already bound it to
  ~0–100; it is no longer presented as a calibrated scale).
* **Breaking:** removed the fixed-cut-point tier helpers
  `categorise_odour()`/`odour_colour()`, `categorise_dust()`/`dust_colour()`,
  `categorise_litter()`/`litter_colour()`, and the interim
  `odour_index_interim()`.
* A new issue tracks the **site calibration tooling** (fit physical → site index
  from complaint/observation records).

## C8 — Upslope rim-venting to elevated rim receptors (issue #24)

* New `rim_venting = FALSE` parameter on `odour_exposure()` and `odour_risk()`.
  When `TRUE` with `terrain_backend = "descriptors"`, the morning pathway-1a
  crosswind factor is gated by a **vertical reach function**: the pulse reaches a
  receptor at height `z_j` only once the morning vented layer
  `h_vent = pool_top + RIM_LIFT_COEF * cbl_cumsum` has risen above it.
* Per-receptor **upslope alignment** (`align_j = max(0, cos(brng_rec_src − aspect_j))`)
  further scales the 1a term; receptors whose downslope aspect faces away from the
  source contribute zero regardless of pool depth.
* **Default-off and uncalibrated.** `rim_venting = FALSE` produces bit-identical
  output to pre-C8. Constants `RIM_LIFT_COEF = 0.2` and `RIM_DELTA = 25` are
  screening defaults flagged for calibration in #8.
* Receptor height `z_j` follows a D1 priority ladder:
  `receptors$elevation − source_elevation` → `receptors$rel_elevation` (C5
  setup-time DEM) → 0 (no-op). All receptors with `z_j = 0` bypass the gate.
* `mh_terrain_from_dem()` now returns a per-receptor `aspect` column (downslope
  direction, degrees) in `$receptor_fields`, consumed by C8 on the run path.
* **Run-path safety contract:** `odour_exposure()` never calls `terra::terrain()`
  or `terra::extract()`; the aspect and relative elevation columns are read from
  stored receptor features (setup-time C5 output). A §5 regression test locks
  this boundary.
* Replaces the removed M2 receptor impaction (#20): M2 had the wrong sign for
  basin geometry; C8 encodes the correct pathway (morning anabatic venting after
  overnight cold-pool accumulation).

## C9 — Shared hazard core: wind dilution active in odour_exposure() (issue #25)

* **Bug fix:** `odour_exposure()` previously normalised relative concentration
  against a constant `H_ref = PM_MAX / (U_CALM_FLOOR * H_MIX_FALLBACK_STABLE)`,
  meaning `u_eff` did not appear in the exposure formula and M3 valley
  sheltering had no effect on `odour_risk()` outputs.
* New internal helper `.odour_hazard_raw(G, vs)` — the ventilation flux
  `G * PM * W_rain / (u_eff * h_mix)` — is now the single source of truth
  shared by `odour_hazard()` and `odour_exposure()`.
* The exposure normaliser is now geometry-based: `PM_MAX / (U_CALM_FLOOR *
  sigma_y_ref * sigma_z_ref)` at `X_REF_EXPOSURE = 250 m`, Briggs class F,
  so the reference denominator matches the `1 / (u_eff * sigma_y * sigma_z)`
  form used in the per-hour concentration.
* `odour_hazard()` gains `terrain`, `shelter`, and `shelter_h_mix` parameters,
  wiring M3 valley sheltering directly into the ventilation index.
* **Consequence:** exposure and risk outputs are lower at short distances (
  < ~1 km, neutral–unstable conditions) where the old constant normaliser
  under-estimated dilution; ratios across hours and sites are now physically
  consistent with the `1/u_eff` advective term.
* Pre-existing bug in the terrain backend fixed: morning-release `r_scale`
  was not suppressed for calm/pool-NA hours, producing a spurious spike for
  calm mornings. Now zeroed via `r_scale_eff`.

## Phase 2 — DEM terrain helper (C5, issue #19)

* New setup-time function `mh_terrain_from_dem()` derives all C1 `mh_terrain`
  descriptors automatically from a DEM using data-driven analysis scales.
  Requires `whitebox` + the WhiteboxTools binary (`whitebox::install_whitebox()`)
  and `terra`. Never on the run path.
* Analysis scales are recorded in `terrain@meta` (`relief_radius`,
  `valley_dev_scale`, `shelter_fetch_L`, `drainage_catchment_radius`).
* `valley_depth` uses multi-scale DEV delineation — **no flow-accumulation
  channel threshold** — matching the plan's §7 threshold-free requirement.
* Caller overrides compose: any descriptor may be pre-set to skip its DEM
  derivation; `meta` flags it as `"user_supplied"`.
* New `Suggests`: `whitebox`, `terra`, `elevatr`, `withr`.

## Phase 3 — Valley sheltering (M3) (C6, issue #20)

* `ventilation_state()` gains `shelter = FALSE` and `shelter_h_mix = FALSE`.
  When `shelter = TRUE`, a wind-speed-regime-resolved valley-sheltering
  reduction is applied to `u_eff` using `terrain@shelter_index`. **Off by
  default until calibrated** (calibration → #8).
* **Terrain-modifier precedence:** M1 drainage confinement (C3b) and M3
  valley sheltering are mutually exclusive on drainage-active hours
  (`DRAINAGE_SHELTER_OVERLAP = 1.0`). The no-stack contract is regression-
  guarded by a behaviour spec in `test-odour-pathways.R`.
* All new coefficients are named constants in `ODOUR_CONSTANTS` (provenance-
  commented; flagged uncalibrated → #8). No inline literals.

# meteoHazard 0.2.0

## Performance

* `odour_exposure()` was rewritten to evaluate the source/receptor dispersion
  as `(hour x receptor)` matrix operations rather than a per-receptor loop, and
  the flat Gaussian dispersion is now computed once and shared with the
  terrain backend instead of being recomputed. `.pool_partition()` and the
  terrain pathway factors (`.cw_venting()`, `.cw_fumigation()`) are fully
  vectorised. On a 3-source x 100-receptor x 168-hour case this is ~6x faster
  for `terrain_backend = "none"` and ~10x faster for `"descriptors"`; output is
  unchanged (guarded by golden characterisation tests).

## Phase 1 — unified geo-referenced site model and corrected dispersion physics

New S7 domain objects, a corrected Gaussian-plume odour architecture,
terrain-aware morning-pulse pathways, a reworked ventilation-state
decomposition, and a geometry-aware litter-exposure rewrite. **Pre-1.0
breaking-change release: no compatibility shims are provided.**

### New domain objects

* `mh_terrain` — S7 value object carrying site terrain descriptors
  (`drainage_bearing`, `flow_convergence`).
* `mh_site` — S7 object holding geo-referenced source, receptor, and barrier
  features (sf geometry) plus a hazard-role table. Automatically reprojects any
  input CRS to the supplied metric `epsg`.

### New functions

* `ventilation_state(met_data, terrain)` — decomposes an hourly met data frame
  into the full dispersion state: effective wind speed (`u_eff`), mixing depth
  (`h_mix`), Pasquill-Gifford stability class (`s`), calm flag (`is_calm`),
  daytime flag (`is_day`), peak-to-mean factor (`PM`), rainfall factor
  (`W_rain`), cold-pool depth (`pool_top`), CBL growth rate (`cbl_growth`), and
  overnight residual-layer wind at 80/120/180/10 m (`residual_wind`).
* `site_from_sectors(sectors, centroid, ring_radius, epsg)` — builds an
  `mh_site` from a compass-sector data frame, constructing wedge-arc barrier
  polygons projected from the site centroid. Used to feed `litter_exposure()`.

### Breaking API changes

* **`odour_exposure(met_data, site, stability, map_c50, terrain_backend)`** —
  `site` is now an `mh_site` object (previously a `receptors` data frame with
  `bearing` and `distance` columns). New `terrain_backend = c("none",
  "descriptors")` argument activates cold-pool / morning-fumigation pathways.
  Area sources (polygon features) use ISC3-derived initial spreads; multi-source
  concentrations are summed, mapped 0–100, then maximised over receptors.
* **`odour_risk()`** — thin wrapper; passes through all new arguments.
* **`litter_exposure(hazard, wind_direction_10m, site)`** — `site` is now an
  `mh_site` (previously a compass-sector data frame). Use `site_from_sectors()`
  to build an `mh_site` from a sector data frame.

### New required met columns

`ventilation_state()` (and therefore `odour_exposure()` / `odour_risk()`) now
accepts six multi-level wind columns. All are optional individually; the code
degrades gracefully along the 180→120→80→10 m fallback ladder.

| Column | Unit | Purpose |
|---|---|---|
| `wind_speed_80m` | m/s | Residual-layer wind; shear-stability estimator |
| `wind_direction_80m` | ° | Residual-layer direction |
| `wind_speed_120m` | m/s | Improves residual-wind accuracy (optional) |
| `wind_direction_120m` | ° | |
| `wind_speed_180m` | m/s | Improves residual-wind accuracy (optional) |
| `wind_direction_180m` | ° | |

Fetch all six from Open-Meteo by adding `&models=best_match` with hourly
variables `wind_speed_80m,wind_direction_80m,wind_speed_120m,wind_direction_120m,
wind_speed_180m,wind_direction_180m` (the existing `&wind_speed_unit=ms`
request parameter applies).

### New dependencies

* `sf` — geo-referenced feature geometries.
* `S7` — the OOP framework for `mh_site` / `mh_terrain`.

# meteoHazard 0.1.0

First release under the name **meteoHazard** (formerly the `TWL` package),
reframed as a collection of meteorological hazard predictors for waste
management: Thermal Work Limit (`generate_twl()`), odour, wind-blown litter,
and dust.

## Explicit units (the `units` package)

Every dimensional quantity in the package is now unit-checked, so a wrong-unit
input can no longer be silently misread (the class of bug behind the earlier
km/h-vs-m/s litter/dust error).

* **Inputs** may be supplied either as bare numerics in the documented unit or
  as `units` objects, which are converted automatically; a dimensionally
  incompatible unit (e.g. a temperature where a speed is expected) is a classed
  error. This applies to the wind, temperature, pressure, radiation, length and
  density inputs of `generate_twl()`, `odour_hazard()`/`odour_exposure()`,
  `litter_hazard()`/`litter_hazard_vec()` and `dust_hazard()`/`dust_flux()`
  (including receptor distances). Percentage and ratio fields (relative
  humidity, cloud cover, soil moisture) and bearings (degrees) are taken as-is.
* **Outputs** that are genuine physical quantities are returned as `units`
  objects: `generate_twl()` now returns W/m^2. `categorise_twl()` and
  `twl_colour()` are units-aware (they accept the typed output or a bare
  numeric). The dimensionless scores -- the 0-100 litter/dust/odour indices and
  the relative odour hazard -- stay plain numeric.
* Internal conversions (km/h<->m/s, hPa<->kPa, etc.) now go through `units`
  rather than hand-rolled factors. `units` is a new hard dependency.

## Cross-function consistency

The four hazard families were swept for naming, unit, and helper consistency
(pre-v1; no backwards-compatibility shims):

* **Units.** All wind inputs are now in **m/s** across every function (matching
  the bundled Open-Meteo fetcher, which requests `&wind_speed_unit=ms`).
  Previously the litter and dust functions expected km/h while odour expected
  m/s — the same column name silently meant different units. Affected defaults
  were converted to preserve their physical thresholds.
* **Naming.** The risk-index families converge on a `*_hazard()` /
  `*_exposure()` scheme: `litter_hazard()` (was `generate_litter_risk_index`),
  `litter_hazard_vec()` (was `litter_risk_index`), `dust_hazard()` (was
  `generate_dust_risk_index`), `dust_flux()` (was `dust_emission_potential`),
  and `odour_risk()` (was `generate_odour_risk_index`). `generate_twl()` is
  unchanged (it returns W/m^2 with physiological zones, not a 0–100 index).
* **Hazard tiers.** New `categorise_litter()`/`litter_colour()`,
  `categorise_dust()`/`dust_colour()`, and `categorise_odour()`/`odour_colour()`
  apply the documented operational bands on a shared green/amber/orange/red
  palette, mirroring `categorise_twl()`/`twl_colour()`.
* Shared internal validation helpers (`.assert_required_cols()`,
  `.assert_numeric_cols()`, `.na_fill()`) replace the duplicated input-check
  blocks; the litter/dust missing-column errors are now classed
  `meteoHazard_input_error` like the rest.

## Odour model: hazard / exposure split

* The monolithic odour risk index is split into two layers, mirroring the
  litter functions: `odour_hazard()` (receptor-independent, direction-agnostic)
  and `odour_exposure()` (geometry-aware). `odour_risk()` is now a thin
  convenience wrapper over the two.
* `odour_hazard()` returns a relative **ventilation index** (source emission
  divided by wind speed times mixing depth), capturing the dominant
  calm/stable/shallow-boundary-layer signal without baking in receptor
  geometry. Stability now uses the Pasquill-Turner insolation/cloud scheme by
  default (self-consistent with the Briggs dispersion curves), with the legacy
  10 m/80 m shear exponent available via `stability = "shear"`;
  `wind_speed_80m` is therefore now optional. The temperature response
  (`V_mod`) is widened, and a new stability-dependent peak-to-mean factor
  accounts for the sub-minute concentration peaks that drive odour annoyance.
* `odour_exposure()` reconstructs the Pasquill-Gifford Gaussian relative
  concentration per receptor and maps it to a 0–100 worst-case band. Direction
  is now physically coupled to distance and stability (replacing the old fixed
  10-degree gate), with forecast-direction uncertainty added as a separate
  term, and receptors upwind of the source are correctly excluded from the
  plume. The terrain-aware drainage / morning-fumigation refinement is
  relocated here, gated on the optional `drainage_axes`.
* The model defaults are provisional and uncalibrated (#8); the 0–100 mapping
  is tunable via `map_c50`.

## TWL function-chain hardening (#1)

* `generate_twl()` now validates its inputs up front, before any API call,
  raising a classed `meteoHazard_input_error` for a non-`POSIXct` `datetime`,
  length mismatches, and out-of-range coordinates or weather values. `NA`
  values are permitted and propagate.
* Fixed the observation count: results are now always `length(datetime)` long,
  with length-1 inputs recycled, rather than being derived from `temp` alone.
* `convert_pressure` now applies only to **user-supplied** pressure;
  API-fetched pressure (Open-Meteo `surface_pressure`, in hPa) is always
  converted to kPa. This fixes the previously incorrect
  `convert_pressure = FALSE` + API-fetch combination.
* `fetch_openmeteo()` now splits a request across the forecast and archive
  endpoints when the date range spans the ~92-day boundary, instead of
  silently returning `NA` for the half a single endpoint cannot serve.
* `solve_twl_single()` now checks that the heat-balance root is bracketed and
  warns (naming the observation) when it is not, instead of silently returning
  a non-root estimate.
* The globe-temperature and natural-wet-bulb solvers now draw air thermal
  conductivity and kinematic viscosity from shared constants; the globe solver
  uses the sourced value 0.028 W/(m·K) (Brake 2001) it previously did not.
* Internal clean-ups: shared `TWL_FLOOR`/`TWL_CEILING` constants, removal of an
  unreachable clamp branch, corrected progress-bar accounting, and an updated
  Open-Meteo user-agent.
