# meteoHazard

An R package that turns meteorological data into **management-relevant predictions and warnings** for waste-management operations. It is part of the [`tidyWaste`](https://github.com/vorpalvorpal) family of packages.

## Functions

| Function | Hazard | Status |
|---|---|---|
| `generate_twl()` | **Thermal Work Limit** — heat stress on workers | ✅ Implemented |
| `ventilation_state()` | **Ventilation state** — full dispersion state including pool-top and residual wind | ✅ Implemented |
| `odour_hazard()` | **Odour nuisance hazard** — site ventilation index (direction-agnostic) | ✅ Implemented |
| `odour_exposure()` | **Odour nuisance exposure** — Gaussian-plume downwind impact on geo-referenced receptors | ✅ Implemented |
| `litter_hazard()` | **Wind-blown litter hazard** — entrainment & transport at the working face | ✅ Implemented |
| `litter_exposure()` | **Wind-blown litter exposure** — where litter goes, given wind direction and site geometry | ✅ Implemented |
| `dust_hazard()` | **Dust hazard** — wind erosion of exposed surfaces | ✅ Implemented |
| `site_from_sectors()` | **Site constructor** — builds an `mh_site` from a compass-sector data frame | ✅ Implemented |

## Units

Every dimensional input is unit-checked via the [`units`](https://r-quantities.github.io/units/) package. You can pass a value either as a **bare numeric in the documented unit** or as a **`units` object**, which is converted automatically — a dimensionally incompatible unit (say, a temperature where a wind speed is expected) is a hard error. This makes the historical "same column name, different unit" bug impossible:

```r
library(units)
# these two calls are equivalent — the tagged wind is converted to m/s
litter_hazard_vec(set_units(57.6, "km/h"), set_units(45, "km/h"), 0, 0.02)
litter_hazard_vec(16, 12.5, 0, 0.02)
```

Outputs that are genuine physical quantities are returned as `units` objects — `generate_twl()` returns W/m² (use `units::drop_units()` for a bare numeric; `categorise_twl()`/`twl_colour()` accept either). The hazard/exposure indices (relative odour concentration, relative dust flux, litter index) are dimensionless relative quantities and stay plain numeric. Percentages (relative humidity, cloud cover), ratios (soil moisture) and bearings (degrees) are taken as-is.

## Installation

Install from GitHub using the `remotes` or `pak` package:

```r
# Using remotes
install.packages("remotes")
remotes::install_github("vorpalvorpal/meteoHazard")

# Using pak
install.packages("pak")
pak::pkg_install("vorpalvorpal/meteoHazard")
```

Or install from a local clone:

```r
# Clone the repository, then:
install.packages(".", repos = NULL, type = "source")
```

### Dependencies

The package requires R (>= 4.1) and the following packages, which are installed automatically:

- `checkmate`
- `cli`
- `dplyr`
- `httr2`
- `purrr`
- `sf`
- `S7`
- `units`

## Thermal Work Limit (TWL)

The **Thermal Work Limit (TWL)** is the maximum sustainable metabolic rate (W/m²) that well-hydrated, acclimatised workers can maintain in a thermally stressful environment. It implements the Brake & Bates (2002) methodology.

### Basic calculation (all parameters supplied)

```r
library(meteoHazard)

twl <- generate_twl(
  datetime     = as.POSIXct("2024-01-15 10:00:00", tz = "UTC"),
  latitude     = -31.95,
  longitude    = 115.86,
  temp         = 36,    # °C dry bulb
  wind_speed   = 0.5,   # m/s at body level
  RH           = 40,    # % relative humidity
  direct_solar = 600,   # W/m²
  diffuse_solar = 100,  # W/m²
  pressure     = 1013   # hPa
)
```

### Auto-fetch missing weather from Open-Meteo

Any weather parameter can be omitted and will be fetched automatically from the [Open-Meteo API](https://open-meteo.com/) (free, no API key required). An internet connection is required.

```r
twl <- generate_twl(
  datetime  = as.POSIXct("2024-01-15 10:00:00", tz = "UTC"),
  latitude  = -31.95,
  longitude = 115.86,
  temp      = 36
  # wind_speed, RH, direct_solar, diffuse_solar, pressure fetched from API
)
```

### Categorise results

```r
categorise_twl(twl)
# "Unrestricted", "Acclimatisation", "Buffer", or "Withdrawal"

twl_colour(twl)
# Hex colour codes for visualisation
```

### TWL Categories

| TWL (W/m²) | Category | Colour |
|---|---|---|
| ≥ 220 | Unrestricted | Green |
| 140 – 220 | Acclimatisation | Amber |
| 115 – 140 | Buffer | Orange |
| < 115 | Withdrawal | Red |

## Site model

The geometry-aware hazard layers (`odour_exposure()`, `litter_exposure()`) share a common geo-referenced site model built from [sf](https://r-spatial.github.io/sf/) and [S7](https://rconsortium.github.io/S7/).

### `mh_site` and `mh_terrain`

An `mh_site` object holds:
- An `sf` feature collection of source, receptor, and barrier geometries (points or polygons), all in a common metric CRS.
- A **roles table** that maps each feature to a hazard type (`"odour"`, `"litter"`, …) and a role (`"source"`, `"receptor"`, `"barrier"`), plus role-specific attributes (e.g. `emit_height`, `permeability`, `sensitive`).

An `mh_terrain` object carries optional site terrain descriptors (`drainage_bearing`, `flow_convergence`) used by the cold-pool / morning-fumigation physics.

```r
library(meteoHazard)
library(sf)

terrain <- mh_terrain(drainage_bearing = 135, flow_convergence = 0.7)
```

### Deriving terrain descriptors from a DEM (setup-time)

`mh_terrain_from_dem()` derives the `mh_terrain` descriptors automatically from a digital elevation model. It is a **setup-time-only** function (runs WhiteboxTools flow-routing; takes seconds to minutes); call it once and save the result.

```r
# One-time setup: derive terrain from a DEM raster.
dem     <- terra::rast("site_dem.tif")
source  <- sf::read_sf("source_location.gpkg")
terrain <- mh_terrain_from_dem(dem, source = source, epsg = 32755L)
# terrain is an mh_terrain object — pass to mh_site() as shown above.
```

Requires the `whitebox` and `terra` packages and the WhiteboxTools binary (`whitebox::install_whitebox()`).

### `site_from_sectors()` — litter sites

For the litter model, `site_from_sectors()` builds an `mh_site` from a compact compass-sector data frame: it places a point source at the site centroid and constructs wedge-arc barrier polygons for each sector.

```r
sectors <- data.frame(
  arc_start    = c("NE", "SW"),
  arc_end      = c("SE", "NW"),
  permeability = c(1.0, 0.3),    # open to the E; tree belt to the W
  sensitive    = c(TRUE, FALSE)
)

site <- site_from_sectors(
  sectors,
  centroid = c(115.86, -31.95),  # lon/lat of the site
  epsg     = 32755L              # metric CRS to project into (UTM 55S here)
)
```

## Odour dispersal risk

The odour model is split into two layers:

* `odour_hazard()` — a **receptor-independent, direction-agnostic** hourly hazard index (source emission divided by ventilation). Answers *how strong is the odour situation around the site this hour* and returns a relative index (baseline = 1.0).
* `odour_exposure()` — the **geometry-aware** layer: maps the hazard through a Pasquill-Gifford Gaussian plume onto each receptor in an `mh_site`, applies distance decay, direction uncertainty, and (optionally) terrain-driven cold-pool and morning-fumigation corrections, then returns the **per-receptor relative concentration** as a `n_hours × n_receptors` matrix (the unbounded physical layer).

Mapping that physical value onto a bounded index (e.g. 0–100) and into site-specific tiers is a **calibration** step left to the consumer (e.g. a dashboard that knows the site's complaint history). The package no longer ships a fixed 0–100 odour scale; site-specific calibration tooling is planned (issue #11). Reduce over receptors for a worst-case summary with e.g. `apply(odour_risk(met_data, site), 1, max)`.

`odour_risk()` is a convenience wrapper that runs both in one call. None of these functions call the weather API — the caller fetches the hourly variables from [Open-Meteo](https://open-meteo.com/) (requesting `&wind_speed_unit=ms`) and passes them as a data frame, one row per consecutive hour.

```r
library(meteoHazard)
library(sf)

# Build the site: a point source and two receptor points in metric CRS
src  <- st_sfc(st_point(c(0, 0)),       crs = 32755)
rec1 <- st_sfc(st_point(c(500, 0)),     crs = 32755)
rec2 <- st_sfc(st_point(c(-800, 200)),  crs = 32755)

feats <- st_sf(
  id       = c("source", "rec_east", "rec_west"),
  geometry = c(src, rec1, rec2)
)
roles <- data.frame(
  feature_id  = c("source", "rec_east", "rec_west"),
  hazard      = "odour",
  role        = c("source", "receptor", "receptor"),
  emit_height = c(5, NA, NA)
)
site <- mh_site(features = feats, roles = roles, epsg = 32755L)

# met_data: one row per hour, required Open-Meteo columns:
# wind_direction_10m, wind_speed_10m, boundary_layer_height, temperature_2m,
# pressure_msl, precipitation, relative_humidity_2m, cloud_cover,
# direct_radiation, soil_moisture_0_to_1cm, soil_moisture_1_to_3cm
# Optional upper-level wind columns improve pool-top and terrain physics:
# wind_speed_80m, wind_direction_80m (and 120m, 180m variants).

odour <- odour_risk(met_data, site)

# With terrain-aware cold-pool / morning-fumigation physics:
terrain <- mh_terrain(drainage_bearing = 135, flow_convergence = 0.7)
site_t  <- mh_site(features = feats, roles = roles, epsg = 32755L,
                   terrain = terrain)
odour_terrain <- odour_risk(met_data, site_t, terrain_backend = "descriptors")
```

### Terrain mechanisms

The terrain backend activates two physical mechanisms, plus an optional rim-venting extension:

- **M1 drainage confinement** — on nocturnal hours when a cold pool is present (`drainage_active = TRUE`), drainage flow channels along the valley axis, concentrating odour downslope.
- **M3 valley sheltering** (`shelter = TRUE`) — reduces effective wind speed on hours when the site is enclosed by surrounding terrain, scaled by the `shelter_index` (topographic openness) and wind speed.
- **Upslope rim-venting** (`rim_venting = TRUE`) — extends the M1 morning pulse to elevated rim receptors: during the morning inversion break-up the vented layer climbs the slopes, so a ridge-top receptor is exposed only once the morning vented layer (cold-pool depth plus convective growth) reaches its elevation *and* it sits up the slope the venting flow is climbing. Requires per-receptor `rel_elevation`/`elevation` and `aspect` (from `mh_terrain_from_dem()`); a strict no-op without them. **Off by default — uncalibrated screening defaults (`RIM_LIFT_COEF`, `RIM_DELTA`); calibration pending.**

**Nocturnal cold-pool cap on `h_mix`** (`pool_cap = TRUE`, default, on `ventilation_state()`/`odour_hazard()`/`odour_exposure()`/`odour_risk()`) — on stable nights the model caps the dilution depth at `min(boundary_layer_height, pool_top)` instead of always using the raw forecast boundary-layer height, so a shallow thermal cold pool (not an unreliable synoptic BLH field) sets the ventilation depth. This is the mechanism that captures the site's known worst case: winter temperature inversions trapping odour in the valley, with the strongest impact at ridge-top receptors. The cap only engages at night, and only when the pool is thermally derived from `temperature_2m` + `relative_humidity_2m`; set `pool_cap = FALSE` to restore the raw boundary-layer-height behaviour.

**Precedence**: M1 takes priority over M3 on any hour when `drainage_active` is `TRUE`. The suppression strength is controlled by `ODOUR_CONSTANTS$DRAINAGE_SHELTER_OVERLAP` (default 1.0 = full mutual exclusion on that hour).

```r
terrain <- mh_terrain(
  drainage_bearing = 135,
  flow_convergence = 0.7,
  shelter_index    = 60    # degrees; ~90 = open plain, <60 = strongly enclosed
)
site_t <- mh_site(features = feats, roles = roles, epsg = 32755L, terrain = terrain)
odour_shelter <- odour_risk(met_data, site_t,
                            terrain_backend = "descriptors", shelter = TRUE)
```

`shelter_index` can be provided directly or derived automatically via `mh_terrain_from_dem()`.

### Rain scavenging and odorant solubility

`odorant_solubility` (default 0.5, on the same four functions as `pool_cap`) controls how strongly rain washes odorant out of the plume. Below-cloud washout scales with Henry's-law solubility: `0` represents poorly soluble reduced-sulfur odorants (the compounds that drive most landfill complaints), which barely wash out; `1` represents highly soluble species (ammonia, amines) and reproduces the old fixed washout tiers; the default `0.5` is a mixed sulfur/soluble profile.

```r
# Model a poorly soluble reduced-sulfur odorant (little rain washout):
odour_dry <- odour_risk(met_data, site, odorant_solubility = 0)
```

`odour_exposure()` returns the physical per-receptor relative concentration; the package no longer ships a fixed 0–100 odour scale or operational tiers. Mapping the physical output onto a site-specific index and tiers (potentially 0–100) is a calibration step delivered by forthcoming calibration tooling (issue #11). Reduce over receptors for a worst-case summary with e.g. `apply(odour_risk(met_data, site), 1, max)`.

Stability defaults to Pasquill-Turner (insolation/cloud and wind), with a legacy 10 m/80 m shear estimator available via `stability = "shear"`. See `?odour_hazard`, `?odour_exposure`, and `?ventilation_state` for the full model and references.

## Wind-blown litter

Litter is split into two composable layers: a **hazard** index (the meteorology at the working face) and an **exposure** layer (where the litter goes, given wind direction and site geometry).

### Hazard index

`litter_hazard()` computes an hourly litter **hazard** as a relative index: how strongly is litter being entrained from the working face and moved. It is grounded in aeolian wind-erosion physics — friction-velocity (`u*`) entrainment with the EPA AP-42 excess-squared form, a moisture-raised threshold (Fécan et al. 1999) plus a saturation veto, a rainfall hard gate, and a mean-wind transport-potential multiplier. It is **direction-agnostic**: wind direction and barriers belong to the exposure layer below. (With the default `entrainment_max`/`transport_max` the index spans ~0–100 by construction, but that scaling is a modelling choice, not a calibrated operational scale.)

The caller supplies pre-fetched Open-Meteo data (one row per hour) with four columns:

```r
# met_data columns: wind_gusts_10m (m/s), wind_speed_10m (m/s),
#                   precipitation (mm), soil_moisture_0_to_1cm (m³/m³).
# Fetch winds from Open-Meteo with &wind_speed_unit=ms.

hazard <- litter_hazard(met_data)
```

The vector API `litter_hazard_vec()` takes the columns directly (handy inside `dplyr::mutate()`) and exposes all calibration parameters (roughness `z0`, threshold/reference friction velocities, moisture-threshold gain/curvature, transport onset/reference, rain threshold). Mapping the index onto operational tiers is a site-specific calibration step (forthcoming tooling, issue #11) — the package no longer ships fixed litter tiers.

### Exposure layer

`litter_exposure()` sits on top of the hazard index and answers the consequence question: given the hazard, the wind direction, and the site geometry, where does the litter end up? The site geometry is an `mh_site` object; use `site_from_sectors()` to build one from a compact compass-sector data frame (see the **Site model** section above).

```r
sectors <- data.frame(
  arc_start    = c("NE", "SW"),
  arc_end      = c("SE", "NW"),
  permeability = c(1.0, 0.3),    # open receptor to the E; tree belt to the W
  sensitive    = c(TRUE, FALSE)
)
site <- site_from_sectors(sectors, centroid = c(115.86, -31.95), epsg = 32755L)

exposure <- litter_exposure(hazard, met_data$wind_direction_10m, site)
# data frame: exposure (hazard attenuated by direction) and a severity zone
#   within_face < on_site < off_site
```

The directional factor `M` is the highest permeability among barrier sectors whose arc (expanded by `direction_tol` on each edge) contains the downwind bearing — worst-case passage probability. When no sector is hit, `default_permeability` applies. The severity zone uses the raw hazard for mobility and the hit-sector `sensitive` flag for destination classification.

See `?litter_hazard_vec`, `?litter_hazard`, `?litter_exposure`, and `?site_from_sectors` for the full parameter lists, model structure, and references.

## Dust hazard

`dust_hazard()` computes an hourly **relative dust flux** from wind erosion of an exposed, erodible surface (a `met_data`-frame convenience wrapper over `dust_flux()` that adds the crust gate). It uses a physical saltation-to-emission chain — a Shao & Lu (2000) threshold friction velocity, a Fécan et al. (1999) soil-moisture correction, an Owen/White saltation flux, and the MB95 sandblasting efficiency — driven by the gust (fastest-mile proxy). The erodible surface is assumed smooth (no non-erodible roughness drag partition). The package no longer ships a fixed 0–100 dust index.

Meteorological inputs are pre-fetched Open-Meteo columns; the surface is described by one-time site-survey parameters (Tyler sieve number for the modal aggregate size, clay %, roughness, bulk density).

```r
# met_data columns: wind_speed_10m (m/s), wind_gusts_10m (m/s),
#                   soil_moisture_0_to_1cm (m³/m³); + precipitation (mm) when
#                   crust = TRUE. Fetch winds with &wind_speed_unit=ms.

dust <- dust_hazard(
  met_data,
  tyler_sieve_no = 20L,   # modal aggregate size of the erodible surface
  clay_percent   = 10
)
```

An optional precipitation **crust gate** (`crust = TRUE`) raises the erosion threshold for days after rain — a memory effect that instantaneous soil moisture cannot capture — and is off by default so it can be enabled per site. The lower-level `dust_flux()` returns the underlying relative dust flux without the `met_data`/crust convenience. See `?dust_hazard` and `?dust_flux` for the full parameter list and references.

Mapping the relative flux onto operational tiers is a site-specific calibration step (forthcoming tooling, issue #11) — the package no longer ships fixed dust tiers.

## Reference

Brake, D.J. and Bates, G.P. (2002) Limiting Metabolic Rate (Thermal Work Limit) as an Index of Thermal Stress. *Applied Occupational and Environmental Hygiene*, 17:3, 176–186. <https://doi.org/10.1080/104732202317070947>
