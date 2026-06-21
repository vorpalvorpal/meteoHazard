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
