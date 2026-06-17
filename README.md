# meteoHazard

An R package that turns meteorological data into **management-relevant predictions and warnings** for waste-management operations. It is part of the [`tidyWaste`](https://github.com/vorpalvorpal) family of packages.

## Functions

| Function | Hazard | Status |
|---|---|---|
| `generate_twl()` | **Thermal Work Limit** — heat stress on workers | ✅ Implemented |
| `generate_odour_risk_index()` | **Odour nuisance** — downwind odour dispersion | ✅ Implemented |
| `generate_litter_risk_index()` | **Wind-blown litter** — material escaping the site boundary | ✅ Implemented |

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

## Wind-blown litter risk

`generate_litter_risk_index()` computes an hourly **Litter Risk Index (LRI)** in the range 0–100 for windblown litter dispersal at a landfill. It uses a multiplicative model with two wetness gates (active rainfall and soil moisture), driven primarily by gust energy and modulated by sustained wind, atmospheric instability, and wind direction relative to sensitive receptors.

As with odour, the caller supplies pre-fetched Open-Meteo data (one row per hour). Receptor directions are given as compass arcs — wind blowing *from* a direction inside an arc is treated as blowing *toward* that receptor.

```r
# met_data columns: wind_gusts_10m, wind_speed_10m, wind_direction_10m,
#   precipitation, soil_moisture_0_to_1cm, cape, is_day

litter <- generate_litter_risk_index(
  met_data,
  receptor_arcs = list(c("W", "N"))   # one receptor: wind from W clockwise to N
)
```

The vector API `litter_risk_index()` takes the individual columns directly (handy inside `dplyr::mutate()`) and exposes all calibration parameters (gust thresholds, soil-moisture gates, off-wind attenuation, etc.). Suggested tier mapping (high uncertainty — requires site calibration):

| LRI | Tier | Meaning |
|---|---|---|
| 0 – 19 | LOW | Dispersal unlikely |
| 20 – 44 | MODERATE | Enhanced controls warranted |
| 45 – 69 | HIGH | Dispersal likely without intervention |
| 70 – 100 | EXTREME | Maximum controls or cessation required |

See `?generate_litter_risk_index` and `?litter_risk_index` for the full parameter list, model structure, and references.

## Reference

Brake, D.J. and Bates, G.P. (2002) Limiting Metabolic Rate (Thermal Work Limit) as an Index of Thermal Stress. *Applied Occupational and Environmental Hygiene*, 17:3, 176–186. <https://doi.org/10.1080/104732202317070947>
