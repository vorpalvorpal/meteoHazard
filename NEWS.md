# meteoHazard 0.1.0

First release under the name **meteoHazard** (formerly the `TWL` package),
reframed as a collection of meteorological hazard predictors for waste
management. The Thermal Work Limit (`generate_twl()`) is the first implemented
function; `predict_odour()` and `predict_litter()` are stubs.

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
