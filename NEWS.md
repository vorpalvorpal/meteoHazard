# meteoHazard 0.1.0

First release under the name **meteoHazard** (formerly the `TWL` package),
reframed as a collection of meteorological hazard predictors for waste
management. The Thermal Work Limit (`generate_twl()`) is the first implemented
function; `predict_odour()` and `predict_litter()` are stubs.

## Odour model: hazard / exposure split

* The monolithic `generate_odour_risk_index()` is split into two layers,
  mirroring the litter functions: `odour_hazard()` (receptor-independent,
  direction-agnostic) and `odour_exposure()` (geometry-aware).
  `generate_odour_risk_index()` is now a backward-compatible convenience wrapper
  over the two.
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
