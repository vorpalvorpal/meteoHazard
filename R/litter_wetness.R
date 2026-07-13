#' Hourly litter-surface wetness state
#'
#' Computes a sequential (hour-over-hour) surface-wetness state `w` in
#' `[0, 1]` for windblown litter: an hour with enough rain resets the surface
#' to fully wet (`w = 1`); otherwise the surface dries exponentially at a rate
#' driven by vapour-pressure deficit (VPD), wind, and insolation. This is a
#' fast, thin-film state deliberately distinct from the gridded
#' `soil_moisture_0_to_1cm` layer used elsewhere in the package: soil-moisture
#' memory operates on a scale of *days* (McColl et al. 2017,
#' \doi{10.1038/ngeo2868}), whereas a film of water or a wet paper/plastic
#' surface holds well under a millimetre and dries off in a matter of hours
#' under sun/wind (Zhou et al. 2021, \doi{10.3390/ijerph18041790}). Feeding
#' [litter_hazard_vec()] the slow soil-moisture proxy when a fast-drying
#' litter surface is actually what governs entrainment understates the
#' litter hazard in the hours just after a shower; this state is the
#' fast-memory alternative, supplied to [litter_hazard_vec()] via its
#' `wetness` argument.
#'
#' The drying-rate form (VPD x wind x insolation multipliers on a base rate)
#' mirrors the aerodynamic-plus-radiative structure of open-surface evaporation
#' (Penman 1948, \doi{10.1098/rspa.1948.0037}: evaporation rises with the
#' vapour-pressure deficit, wind, and available radiation) and the
#' leaf-wetness-duration literature, where those same drivers govern the drying
#' phase of a surface-wetness event (Sentelhas et al. 2008,
#' \doi{10.1016/j.agrformet.2007.09.011}). All rate coefficients here
#' (`dry_rate_base`, `vpd_coef`, `wind_coef`, `sw_coef`) are litter-specific
#' calibration placeholders, not values taken from a validated litter/leaf
#' study -- flagged high calibration uncertainty throughout.
#'
#' @section Units:
#' `precipitation`, `wind_speed_10m`, and `shortwave_radiation` may be
#' supplied either as bare numerics in the documented unit or as
#' \pkg{units} objects, which are converted automatically (a dimensionally
#' incompatible unit is an error). `temperature_2m` and
#' `relative_humidity_2m` are taken as plain numerics (degC, %). The
#' returned wetness state is dimensionless and is a plain numeric.
#'
#' @param precipitation Numeric vector. Hourly precipitation (mm),
#'   Open-Meteo `precipitation`. Drives the rain-reset trigger.
#' @param temperature_2m Numeric vector. Air temperature at 2 m (degC),
#'   Open-Meteo `temperature_2m`. Feeds the Tetens VPD calculation. Not
#'   bounds-checked -- sub-zero values are physically valid.
#' @param relative_humidity_2m Numeric vector. Relative humidity at 2 m (%),
#'   Open-Meteo `relative_humidity_2m`, in `0`-`100`. Feeds the Tetens VPD
#'   calculation.
#' @param wind_speed_10m Numeric vector. Mean wind speed at 10 m (m/s),
#'   Open-Meteo `wind_speed_10m`. Increases the drying rate (evaporative
#'   ventilation of the thin film).
#' @param shortwave_radiation Numeric vector. Incoming shortwave radiation
#'   (W/m^2), Open-Meteo `shortwave_radiation`. Increases the drying rate
#'   (radiative/convective heating of the surface).
#' @param wetness_set_precip Hourly precipitation (mm) at or above which the
#'   surface is reset to fully wet (`w = 1`), overriding the prior hour's
#'   state. Default `0.5`.
#' @param dry_rate_base Base exponential drying rate (per hour) at zero VPD,
#'   calm wind, and no insolation. Default `0.7`. High calibration
#'   uncertainty.
#' @param vpd_coef Drying-rate gain per kPa of vapour-pressure deficit.
#'   Default `1.0`. High calibration uncertainty.
#' @param wind_coef Drying-rate gain per m/s of mean wind. Default `0.1`.
#'   High calibration uncertainty.
#' @param sw_coef Drying-rate gain per unit of `shortwave_radiation / sw_ref`.
#'   Default `0.5`. High calibration uncertainty.
#' @param sw_ref Reference insolation (W/m^2) that normalises
#'   `shortwave_radiation` in the drying-rate term. Default `500`.
#' @param w0 Initial wetness state, in `[0, 1]`, assumed for the hour just
#'   before row 1 (i.e. what the surface looked like going into the forecast
#'   window). Default `0` (dry start).
#'
#' @return Numeric vector of length equal to the inputs, the litter-surface
#'   wetness state (dimensionless, `[0, 1]`) for each forecast hour.
#'
#' @references
#' McColl, K. A. et al. (2017). Global characterization of surface
#' soil-moisture drydowns. \emph{Nature Geoscience} 10, 100-104.
#' \doi{10.1038/ngeo2868}. Motivates a fast litter-wetness state distinct
#' from the (multi-day-memory) gridded soil-moisture layer.
#'
#' Zhou, Q. et al. (2021). Thin liquid films on solid surfaces: formation,
#' properties and health implications. \emph{International Journal of
#' Environmental Research and Public Health} 18(4), 1790.
#' \doi{10.3390/ijerph18041790}. Basis for treating litter surface wetness as
#' a thin film (sub-mm capacity) that dries off quickly.
#'
#' Tetens, O. (1930). Uber einige meteorologische Begriffe.
#' \emph{Zeitschrift fur Geophysik} 6, 297-309. Basis for the saturation
#' vapour pressure approximation used in `.litter_vpd()`.
#'
#' @seealso [litter_hazard_vec()] for how `wetness` is consumed by the
#'   hazard index (`use_wetness_state = TRUE` in [litter_hazard()]).
#' @export
litter_wetness_vec <- function(
  precipitation,
  temperature_2m,
  relative_humidity_2m,
  wind_speed_10m,
  shortwave_radiation,
  wetness_set_precip = 0.5,
  dry_rate_base       = 0.7,
  vpd_coef            = 1.0,
  wind_coef           = 0.1,
  sw_coef             = 0.5,
  sw_ref              = 500,
  w0                  = 0
) {
  # ---- Normalise dimensional inputs (bare = documented unit; units = converted) #
  # temperature_2m and relative_humidity_2m are taken as plain numerics (degC,
  # %); see @section Units.
  precipitation       <- .drop_to(precipitation, "mm", arg = "precipitation")
  wind_speed_10m      <- .drop_to(wind_speed_10m, "m/s", arg = "wind_speed_10m")
  shortwave_radiation <- .drop_to(shortwave_radiation, "W/m^2", arg = "shortwave_radiation")

  # ---- Validate meteorological inputs (complete, in-range, aligned) -------- #
  # precipitation/wind/shortwave are physically >= 0; RH is a
  # percentage in [0, 100]; temperature_2m is intentionally unrestricted
  # (sub-zero air temperatures are physically valid and litter can still be
  # present on a frozen/damp surface).
  n <- length(precipitation)
  checkmate::assert_numeric(precipitation, lower = 0, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(temperature_2m, any.missing = FALSE, len = n)
  checkmate::assert_numeric(relative_humidity_2m, lower = 0, upper = 100,
                            any.missing = FALSE, len = n)
  checkmate::assert_numeric(wind_speed_10m, lower = 0, any.missing = FALSE, len = n)
  checkmate::assert_numeric(shortwave_radiation, lower = 0, any.missing = FALSE, len = n)

  # ---- Validate parameters -------------------------------------------------- #
  checkmate::assert_number(wetness_set_precip, lower = 0)
  # The drying-rate coefficients are kept non-negative so that dry_rate >= 0
  # always, which (together with w0 in [0, 1]) is what keeps w in [0, 1]
  # without needing to clamp the sequential update below (exp(-dry_rate) is
  # then always in (0, 1]).
  checkmate::assert_number(dry_rate_base, lower = 0)
  checkmate::assert_number(vpd_coef, lower = 0)
  checkmate::assert_number(wind_coef, lower = 0)
  checkmate::assert_number(sw_coef, lower = 0)
  checkmate::assert_number(sw_ref, lower = .Machine$double.eps)
  checkmate::assert_number(w0, lower = 0, upper = 1)

  # ---- Vapour-pressure deficit (Tetens 1930) -------------------------------- #
  vpd <- .litter_vpd(temperature_2m, relative_humidity_2m)

  # ---- Drying rate: VPD x wind x insolation multipliers on a base rate ----- #
  # Aerodynamic + radiative evaporation structure (Penman 1948,
  # doi:10.1098/rspa.1948.0037) as reused in leaf-wetness-duration models
  # (Sentelhas et al. 2008, doi:10.1016/j.agrformet.2007.09.011): drying is
  # driven by atmospheric demand (VPD), ventilation (wind), and radiative/
  # convective heating (insolation).
  # All multipliers are >= 1 (coefficients are non-negative and the driving
  # variables are non-negative), so dry_rate >= dry_rate_base always.
  dry_rate <- dry_rate_base *
    (1 + vpd_coef * vpd) *
    (1 + wind_coef * wind_speed_10m) *
    (1 + sw_coef * shortwave_radiation / sw_ref)

  # ---- Sequential state update (mirrors .dust_crust_factor's hourly loop) -- #
  # A rain hour hard-resets the surface to fully wet; otherwise the previous
  # hour's wetness decays exponentially at this hour's drying rate.
  w <- numeric(n)
  w_prev <- w0
  for (i in seq_len(n)) {
    if (precipitation[i] >= wetness_set_precip) {
      w[i] <- 1
    } else {
      w[i] <- w_prev * exp(-dry_rate[i])
    }
    w_prev <- w[i]
  }

  w
}


#' Litter-surface wetness state for a forecast tibble
#'
#' Computes the hourly litter-surface wetness state ([litter_wetness_vec()])
#' for each row of a meteorological forecast tibble. Wraps the vector API,
#' accepting a tibble with named columns rather than individual vectors.
#'
#' @param met_data A tibble (or data frame) with one row per hourly timestep,
#'   in chronological order (the wetness state is sequential), containing at
#'   least the columns `precipitation` (mm), `temperature_2m` (degC),
#'   `relative_humidity_2m` (%), `wind_speed_10m` (m/s), and
#'   `shortwave_radiation` (W/m^2).
#' @param ... Additional calibration parameters forwarded to
#'   [litter_wetness_vec()] (e.g. `wetness_set_precip`, `dry_rate_base`).
#'
#' @return Numeric vector of length `nrow(met_data)`, the litter-surface
#'   wetness state (`[0, 1]`) for each forecast hour.
#'
#' @seealso [litter_wetness_vec()], [litter_hazard()] (`use_wetness_state =
#'   TRUE`).
#' @export
litter_wetness <- function(met_data, ...) {
  checkmate::assert_data_frame(met_data, min.rows = 1)

  required_cols <- c(
    "precipitation", "temperature_2m", "relative_humidity_2m",
    "wind_speed_10m", "shortwave_radiation"
  )
  .assert_required_cols(
    met_data, required_cols, arg = "met_data",
    info = paste0(
      "Required: precipitation (mm), temperature_2m (degC), ",
      "relative_humidity_2m (%), wind_speed_10m (m/s), ",
      "shortwave_radiation (W/m^2)."
    )
  )

  litter_wetness_vec(
    precipitation         = met_data$precipitation,
    temperature_2m        = met_data$temperature_2m,
    relative_humidity_2m  = met_data$relative_humidity_2m,
    wind_speed_10m        = met_data$wind_speed_10m,
    shortwave_radiation   = met_data$shortwave_radiation,
    ...
  )
}


# ---- Internal helpers ------------------------------------------------------ #

# Vapour-pressure deficit (kPa) from air temperature (degC) and relative
# humidity (%), via the Tetens (1930) saturation-vapour-pressure
# approximation: es = 0.6108 * exp(17.27*T / (T+237.3)). vpd = es*(1-RH/100)
# is floored at zero so a supplied RH slightly above 100 (measurement noise)
# never yields a negative deficit.
.litter_vpd <- function(temperature_2m, relative_humidity_2m) {
  es <- 0.6108 * exp(17.27 * temperature_2m / (temperature_2m + 237.3))
  pmax(0, es * (1 - relative_humidity_2m / 100))
}
