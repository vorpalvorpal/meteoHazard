# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development commands

```r
devtools::load_all()          # load package for interactive use
devtools::test()              # run full test suite
devtools::test(filter = "odour-exposure")  # run one test file
devtools::document()          # regenerate man/ from roxygen comments
devtools::check()             # full R CMD check (build + test + lint)
```

Run single test file from the shell:
```bash
Rscript -e "devtools::test(filter = 'odour-exposure')"
```

Accept/review snapshot changes after an intentional algorithmic change:
```r
testthat::snapshot_accept()   # accept all changed snaps
testthat::snapshot_review()   # interactive diff review
```

`bench/` is excluded from the package build (`.Rbuildignore`). Benchmarks live there and are run manually.

## Architecture

### Domain objects (S7)

Two S7 value objects underpin the geometry-aware layers:

- **`mh_terrain`** ([`R/site-model.R`](R/site-model.R)) — scalar terrain descriptors (`drainage_bearing`, `flow_convergence`, `valley_depth`, `shelter_index`, `taf`, …). All fields default to `NA` (flat site). Dimensional inputs accept bare numerics (in documented units) or `units` objects (converted; mismatch = hard error).
- **`mh_site`** ([`R/site-model.R`](R/site-model.R)) — an `sf` feature collection (sources, receptors, barriers) plus a `roles` data frame that maps feature IDs to `(hazard, role)` pairs with role-specific attributes (`emit_height`, `permeability`, `sensitive`). Automatically reprojects to the supplied metric `epsg`.

### Odour pipeline (three layers)

```
odour_hazard(met_data, ...)          # direction-agnostic ventilation index (0–∞, ref=1)
    ↓
odour_exposure(met_data, site, ...)  # Gaussian plume onto mh_site receptors → 0–100
    ↓
odour_risk(met_data, site, ...)      # thin wrapper: hazard + exposure in one call
```

The shared internal helper **`.odour_hazard_raw(G, vs)`** (`R/odour_hazard.R`) computes the raw ventilation flux `G * PM * W_rain / (u_eff * h_mix)`. Both `odour_hazard()` and `odour_exposure()` consume it — this is the single source of truth that ensures wind dilution acts consistently through both layers (C9, issue #25).

### Ventilation state hub

**`ventilation_state(met_data, terrain, stability, shelter, shelter_h_mix)`** ([`R/odour-ventilation.R`](R/odour-ventilation.R)) decomposes met data into the full per-hour dispersion state:

| Field | Meaning |
|---|---|
| `u_eff` | effective wind speed, floored at `U_CALM_FLOOR` |
| `h_mix` | mixing depth (m) |
| `s` | Pasquill-Gifford stability index 0–5 (A=0, F=5) |
| `PM` | peak-to-mean ratio (`PM_MIN`–`PM_MAX`) |
| `W_rain` | below-cloud scavenging factor |
| `pool_top` | nocturnal cold-pool depth (m) |
| `cbl_growth` | convective boundary-layer growth rate |
| `residual_wind` | overnight multi-level wind (80/120/180 m) |

Stability defaults to Pasquill-Turner (`stability = "turner"`); `"shear"` uses the 10 m/80 m log-wind exponent. M3 valley sheltering (`shelter = TRUE`) reduces `u_eff` via `terrain@shelter_index` — off by default (uncalibrated, issue #8).

### Terrain mechanisms (two, mutually exclusive on drainage-active hours)

- **M1 drainage confinement** (`R/odour-pathways.R`) — on nocturnal hours when `drainage_active = TRUE` (`is_channelled & !is_day & !is.na(pool_top) & pool_top > 0`), drainage flow channels along the valley axis. Implemented as morning-pulse pathway helpers: `.pool_partition()`, `.cw_venting()`, `.cw_fumigation()`, `.morning_release()`.
- **M3 valley sheltering** (in `ventilation_state()`) — reduces `u_eff`; scaled by `shelter_index` and wind speed via a two-regime transfer model. Suppressed on `drainage_active` hours (`DRAINAGE_SHELTER_OVERLAP = 1.0`).

Both are active only when `terrain_backend = "descriptors"` and `site@terrain` is non-null.

### Constants

All tunable model parameters live in three named lists in [`R/constants.R`](R/constants.R):

- **`ODOUR_CONSTANTS`** — odour model: `U_CALM_FLOOR`, `PM_MIN/MAX`, stability coefficients, Brunt cooling, pool accumulation, C3b pathway factors, C6 sheltering coefficients, `X_REF_EXPOSURE`.
- **`TWL_CONSTANTS`** — Brake & Bates TWL model.
- **`DUST_CONSTANTS`** — dust model: Shao & Lu threshold (`A_N`, `GAMMA`), densities, `Z0_SMOOTH_RATIO` (smooth-bed `z0 = d/30`), Fécan moisture coefficients, MB95 sandblasting (`MB95_ALPHA_SLOPE/INTERCEPT`, `MB95_CLAY_CAP`).

No inline literals in the physics code.

### Units contract

All dimensional inputs accept bare numerics (assumed in documented unit) or `units` objects (converted; mismatch errors). The helpers `.as_units()` and `.drop_to()` in [`R/units-helpers.R`](R/units-helpers.R) enforce this. Internal arithmetic runs on plain doubles. Physical outputs (`generate_twl()` → W/m²) are returned as `units` objects; dimensionless indices stay plain numeric.

### Internal naming conventions

- **Dot-prefix** (`.function_name()`) marks unexported helpers. Access in tests via `meteoHazard:::`.
- One file per conceptual layer; the file header comment describes its scope.
- `specs/` holds design specs (authoritative: `Odour_v2.md`); retained for provenance.

## Testing

- **Style:** `describe("fn()", { it("behaviour", { ... }) })` throughout.
- **Snapshot tests:** golden values in `tests/testthat/_snaps/`. Re-pin intentionally after algorithmic changes with `snapshot_accept()`.
- **GIS/DEM tests:** gated behind `skip_if_not_installed("whitebox")` / `skip_if_not_installed("terra")` and `whitebox::check_whitebox_binary()`. These skip gracefully without the binary.
- **Unexported helpers:** tested directly via `meteoHazard:::` (e.g. `meteoHazard:::.odour_hazard_raw(G, vs)`).
- **Write all tests before implementation** when adding new behaviour.
