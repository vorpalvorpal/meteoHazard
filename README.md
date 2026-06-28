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

Outputs that are genuine physical quantities are returned as `units` objects — `generate_twl()` returns W/m² (use `units::drop_units()` for a bare numeric; `categorise_twl()`/`twl_colour()` accept either). The 0–100 indices and the relative odour hazard are dimensionless and stay plain numeric. Percentages (relative humidity, cloud cover), ratios (soil moisture) and bearings (degrees) are taken as-is.

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

Mapping that physical value onto a bounded 0–100 index and into site-specific tiers is a **calibration** step left to the consumer (e.g. a dashboard that knows the site's complaint history). The previous saturating 0–100 map is retained — parked and uncalibrated — as `odour_index_interim()` pending a calibration helper; reduce over receptors and map it yourself with e.g. `odour_index_interim(odour_risk(met_data, site))`. See issue #11 for the rationale.

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

The **interim** 0–100 index from `odour_index_interim()` maps to provisional operational tiers (uncalibrated screening defaults — parked pending the calibration helper; the map is tunable via its `map_c50` argument):

| Exposure | Tier | Response |
|---|---|---|
| < 15 | LOW | Normal operations |
| 15 – 40 | MODERATE | Heightened awareness; check cover integrity |
| 40 – 70 | HIGH | Active mitigation — reduce tipping face, deploy suppression |
| > 70 | VERY HIGH | Maximum response — consider ceasing tipping |

`categorise_odour()` returns these tier labels and `odour_colour()` the matching colours.

Stability defaults to Pasquill-Turner (insolation/cloud and wind), with a legacy 10 m/80 m shear estimator available via `stability = "shear"`. See `?odour_hazard`, `?odour_exposure`, and `?ventilation_state` for the full model and references.

## Wind-blown litter

Litter is split into two composable layers: a **hazard** index (the meteorology at the working face) and an **exposure** layer (where the litter goes, given wind direction and site geometry).

### Hazard index

`litter_hazard()` computes an hourly litter **hazard** in the range 0–100: how strongly is litter being entrained from the working face and moved. It is grounded in aeolian wind-erosion physics — friction-velocity (`u*`) entrainment with the EPA AP-42 excess-squared form, a moisture-raised threshold (Fécan et al. 1999) plus a saturation veto, a rainfall hard gate, and a mean-wind transport-potential multiplier. It is **direction-agnostic**: wind direction and barriers belong to the exposure layer below.

The caller supplies pre-fetched Open-Meteo data (one row per hour) with four columns:

```r
# met_data columns: wind_gusts_10m (m/s), wind_speed_10m (m/s),
#                   precipitation (mm), soil_moisture_0_to_1cm (m³/m³).
# Fetch winds from Open-Meteo with &wind_speed_unit=ms.

hazard <- litter_hazard(met_data)
```

The vector API `litter_hazard_vec()` takes the columns directly (handy inside `dplyr::mutate()`) and exposes all calibration parameters (roughness `z0`, threshold/reference friction velocities, moisture-threshold gain/curvature, transport onset/reference, rain threshold). Suggested tier mapping (high uncertainty — requires site calibration):

| Hazard | Tier | Meaning |
|---|---|---|
| 0 – 19 | LOW | Entrainment-and-transport unlikely |
| 20 – 44 | MODERATE | Enhanced controls warranted |
| 45 – 69 | HIGH | Likely without intervention |
| 70 – 100 | EXTREME | Maximum controls or cessation required |

`categorise_litter()` returns these tier labels and `litter_colour()` the matching colours.

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

`dust_hazard()` computes an hourly dust hazard (0–100) from wind erosion of an exposed, erodible surface. It uses a physical saltation-to-emission chain — a Shao & Lu (2000) threshold friction velocity, a Fécan et al. (1999) soil-moisture correction, an Owen/White saltation flux, and the MB95 sandblasting efficiency — driven by the gust (fastest-mile proxy). The erodible surface is assumed smooth (no non-erodible roughness drag partition). The index is normalised against a reference gust on a dry surface, so it keeps resolution rather than saturating.

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

An optional precipitation **crust gate** (`crust = TRUE`) raises the erosion threshold for days after rain — a memory effect that instantaneous soil moisture cannot capture — and is off by default so it can be enabled per site. The lower-level `dust_flux()` returns the underlying relative dust flux. See `?dust_hazard` and `?dust_flux` for the full parameter list and references.

Suggested tier mapping (pre-calibration); `categorise_dust()` returns these labels and `dust_colour()` the matching colours:

| Hazard | Tier | Meaning |
|---|---|---|
| 0 – 24 | LOW | Erosion unlikely |
| 25 – 49 | MODERATE | Enhanced controls warranted |
| 50 – 74 | HIGH | Likely without intervention |
| 75 – 100 | EXTREME | Maximum controls or cessation required |

## Reference

Brake, D.J. and Bates, G.P. (2002) Limiting Metabolic Rate (Thermal Work Limit) as an Index of Thermal Stress. *Applied Occupational and Environmental Hygiene*, 17:3, 176–186. <https://doi.org/10.1080/104732202317070947>
