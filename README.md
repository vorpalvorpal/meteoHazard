# meteoHazard

An R package that turns meteorological data into **management-relevant predictions and warnings** for waste-management operations. It is part of the [`tidyWaste`](https://github.com/vorpalvorpal) family of packages.

## Functions

| Function | Hazard | Status |
|---|---|---|
| `generate_twl()` | **Thermal Work Limit** — heat stress on workers | ✅ Implemented |
| `generate_odour_risk_index()` | **Odour nuisance** — downwind odour dispersion | ✅ Implemented |
| `generate_litter_risk_index()` | **Wind-blown litter hazard** — entrainment & transport at the working face | ✅ Implemented |
| `litter_exposure()` | **Wind-blown litter exposure** — where litter goes, given direction & site geometry | ✅ Implemented |
| `generate_dust_risk_index()` | **Dust hazard** — wind erosion of exposed surfaces | ✅ Implemented |

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

## Odour dispersal risk

`generate_odour_risk_index()` computes an hourly odour risk score for a landfill, representing the worst-case impact at any of the specified downwind receptors. It couples emission generation with atmospheric dispersion through a physics-based Pasquill-Gifford Gaussian-plume factor (distance-, stability-, and boundary-layer-height-dependent), and optionally accounts for terrain-driven katabatic drainage flow and morning fumigation.

Unlike `generate_twl()`, this function does **not** call the weather API. The caller fetches the required hourly variables from the [Open-Meteo `/v1/forecast` endpoint](https://open-meteo.com/) (requesting `&wind_speed_unit=ms`) and passes them as a data frame, one row per consecutive hour.

```r
library(meteoHazard)

# met_data: one row per hour, with the 12 required Open-Meteo columns
# (wind_direction_10m, wind_speed_10m, wind_speed_80m, boundary_layer_height,
#  temperature_2m, pressure_msl, precipitation, relative_humidity_2m,
#  cloud_cover, direct_radiation, soil_moisture_0_to_1cm, soil_moisture_1_to_3cm)

# receptors: bearing (° from site centroid) and distance (m) to each location
receptors <- data.frame(
  bearing  = c(90, 270),
  distance = c(500, 800)
)

odour <- generate_odour_risk_index(met_data, receptors)
```

The raw score is **not** clamped (theoretical maximum ≈ 1.80). Suggested operational tiers (subject to calibration against complaint records):

| Score | Tier | Response |
|---|---|---|
| < 0.15 | LOW | Normal operations |
| 0.15 – 0.40 | MODERATE | Heightened awareness; check cover integrity |
| 0.40 – 0.80 | HIGH | Active mitigation — reduce tipping face, deploy suppression |
| > 0.80 | VERY HIGH | Maximum response — consider ceasing tipping |

Passing the optional `drainage_axes` argument enables terrain-aware drainage and fumigation handling — see `?generate_odour_risk_index` for the full model description and references.

## Wind-blown litter

Litter is split into two composable layers: a **hazard** index (the meteorology at the working face) and an **exposure** layer (where the litter goes, given wind direction and site geometry).

### Hazard index

`generate_litter_risk_index()` computes an hourly litter **hazard** in the range 0–100: how strongly is litter being entrained from the working face and moved. It is grounded in aeolian wind-erosion physics — friction-velocity (`u*`) entrainment with the EPA AP-42 excess-squared form, a moisture-raised threshold (Fécan et al. 1999) plus a saturation veto, a rainfall hard gate, and a mean-wind transport-potential multiplier. It is **direction-agnostic**: wind direction and barriers belong to the exposure layer below.

The caller supplies pre-fetched Open-Meteo data (one row per hour) with four columns:

```r
# met_data columns: wind_gusts_10m, wind_speed_10m, precipitation,
#                   soil_moisture_0_to_1cm

hazard <- generate_litter_risk_index(met_data)
```

The vector API `litter_risk_index()` takes the columns directly (handy inside `dplyr::mutate()`) and exposes all calibration parameters (roughness `z0`, threshold/reference friction velocities, moisture-threshold gain/curvature, transport onset/reference, rain threshold). Suggested tier mapping (high uncertainty — requires site calibration):

| Hazard | Tier | Meaning |
|---|---|---|
| 0 – 19 | LOW | Entrainment-and-transport unlikely |
| 20 – 44 | MODERATE | Enhanced controls warranted |
| 45 – 69 | HIGH | Likely without intervention |
| 70 – 100 | EXTREME | Maximum controls or cessation required |

### Exposure layer

`litter_exposure()` sits on top of the hazard index and answers the consequence question: given the hazard, the wind direction, and the site geometry, where does the litter end up? The site is described as boundary sectors, each with a `permeability` (1 = open, lower = better barrier) and a `sensitive` flag.

```r
site <- data.frame(
  arc_start    = c("NE", "SW"),
  arc_end      = c("SE", "NW"),
  permeability = c(1.0, 0.3),    # open receptor to the E; tree belt to the W
  sensitive    = c(TRUE, FALSE)
)

exposure <- litter_exposure(hazard, met_data$wind_direction_10m, site)
# data frame: exposure (hazard attenuated by direction) and a severity zone
#   within_face < on_site < off_site
```

See `?litter_risk_index`, `?generate_litter_risk_index`, and `?litter_exposure` for the full parameter lists, model structure, and references.

## Dust hazard

`generate_dust_risk_index()` computes an hourly dust hazard (0–100) from wind erosion of an exposed, erodible surface. It uses a physical saltation-to-emission chain — a Shao & Lu (2000) threshold friction velocity, a Fécan et al. (1999) soil-moisture correction, an optional Marticorena & Bergametti (1995) drag partition, an Owen/White saltation flux, and the MB95 sandblasting efficiency — driven by the gust (fastest-mile proxy). The index is normalised against a reference gust on a dry surface, so it keeps resolution rather than saturating.

Meteorological inputs are pre-fetched Open-Meteo columns; the surface is described by one-time site-survey parameters (Tyler sieve number for the modal aggregate size, clay %, roughness, bulk density).

```r
# met_data columns: wind_speed_10m, wind_gusts_10m, soil_moisture_0_to_1cm
#                   (+ precipitation when crust = TRUE)

dust <- generate_dust_risk_index(
  met_data,
  tyler_sieve_no = 20L,   # modal aggregate size of the erodible surface
  clay_percent   = 10
)
```

An optional precipitation **crust gate** (`crust = TRUE`) raises the erosion threshold for days after rain — a memory effect that instantaneous soil moisture cannot capture — and is off by default so it can be enabled per site. The lower-level `dust_emission_potential()` returns the underlying relative dust flux. See `?generate_dust_risk_index` and `?dust_emission_potential` for the full parameter list and references.

## Reference

Brake, D.J. and Bates, G.P. (2002) Limiting Metabolic Rate (Thermal Work Limit) as an Index of Thermal Stress. *Applied Occupational and Environmental Hygiene*, 17:3, 176–186. <https://doi.org/10.1080/104732202317070947>
