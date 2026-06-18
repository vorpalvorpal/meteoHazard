#' Litter Hazard Index for windblown litter at landfill sites
#'
#' Computes an hourly **meteorological litter hazard** in the range `[0, 100]`:
#' given the weather in a forecast hour, how strong is the propensity for loose
#' litter to be entrained from the working face and moved. This is a
#' point-source hazard index, **not** a dispersion model, and it is
#' direction-agnostic — wind *direction*, site geometry, and barriers are
#' handled separately by [litter_exposure()].
#'
#' The function does not query any API. The caller fetches the four required
#' hourly variables from the Open-Meteo `/v1/forecast` endpoint (default metric
#' units) and passes them as numeric vectors, one element per forecast hour.
#'
#' @section Model:
#' The index is a bounded, multiplicative combination of entrainment, transport
#' potential, and a rainfall gate:
#' \deqn{LRI = \min(E \cdot T \cdot R,\; 100)}
#'
#' \describe{
#'   \item{Entrainment `E`}{Driven by the **gust** through friction velocity
#'     `u*` (the surface-shear-stress scale). Using the neutral logarithmic wind
#'     profile, \eqn{u_{*g} = \kappa G / \ln(z/z_0)} with the gust `G` in m/s.
#'     Entrainment follows the EPA AP-42 wind-erosion shape
#'     \eqn{E = a\,\min(1, ((u_{*g}-u_{*t})/\Delta u_*)^n)} with
#'     \eqn{\Delta u_* = u_{*ref} - u_{*t0}}, bounded above by `entrainment_max`.}
#'   \item{Moisture-raised threshold `u*t`}{A damp surface needs a stronger gust
#'     to release litter (Fécan et al. 1999): the threshold rises with surface
#'     moisture, \eqn{u_{*t} = u_{*t0}(1 + \gamma s^{\beta})} where
#'     \eqn{s = \mathrm{clamp}((SM-SM_{dry})/(SM_{wet}-SM_{dry}), 0, 1)}. At or
#'     above `soil_wet` the surface is treated as saturated and entrainment is
#'     vetoed (`E = 0`).}
#'   \item{Transport potential `T`}{Driven by the **mean wind** directly (not
#'     `u*`; transport is flight-height advection). A linear ramp from 1 to
#'     `transport_max` between `wind_transport_onset` and `wind_transport_ref`
#'     m/s — the "how far is it moved" penalty.}
#'   \item{Rainfall gate `R`}{Binary: `R = 0` when `precipitation >=
#'     rain_threshold`, else 1.}
#' }
#' Any single suppressor (rain, saturated surface, sub-threshold gust) drives
#' the index to zero. The maximum attainable value is exactly 100
#' (`entrainment_max * transport_max = 50 * 2`).
#'
#' @param wind_gusts_10m Numeric vector. Peak wind gust at 10 m (m/s),
#'   Open-Meteo `wind_gusts_10m` (fetch with `&wind_speed_unit=ms`). Drives
#'   entrainment.
#' @param wind_speed_10m Numeric vector. Mean wind speed at 10 m (m/s),
#'   Open-Meteo `wind_speed_10m`. Drives transport potential.
#' @param precipitation Numeric vector. Hourly precipitation (mm), Open-Meteo
#'   `precipitation`. Feeds the rainfall hard gate.
#' @param soil_moisture_0_to_1cm Numeric vector. Volumetric water content of the
#'   0--1 cm layer (m^3/m^3), Open-Meteo `soil_moisture_0_to_1cm`. Used as a
#'   relative surface-wetness surrogate that raises the entrainment threshold and
#'   supplies the saturation veto.
#' @param kappa von Karman constant for the log wind profile. Default 0.40.
#' @param z0 Aerodynamic roughness length of the working face (m). Default 0.05.
#'   Must be less than the 10 m reference height. High calibration uncertainty.
#' @param ustar_t0 Dry threshold friction velocity (m/s). Default 0.30.
#' @param ustar_ref Friction velocity at which entrainment saturates (m/s).
#'   Default 1.05. Must exceed `ustar_t0`.
#' @param entrainment_max Maximum entrainment score. Default 50.
#' @param excess_exponent Power-law exponent on the friction-velocity excess.
#'   Default 2 (AP-42); calibratable in `[2, 3]`.
#' @param moisture_gain Maximum fractional increase of the threshold as the
#'   surface approaches wet. Default 2.0. High calibration uncertainty.
#' @param moisture_curve Curvature of the moisture-threshold rise. Default 0.5
#'   (concave, Fécan-type). High calibration uncertainty.
#' @param soil_dry Soil moisture at or below which the surface is fully dry
#'   (m^3/m^3). Default 0.05.
#' @param soil_wet Soil moisture at or above which the surface is saturated and
#'   entrainment is vetoed (m^3/m^3). Default 0.20. Must exceed `soil_dry`.
#' @param wind_transport_onset Mean wind below which transport adds nothing
#'   (m/s). Default 5.5 (~20 km/h).
#' @param wind_transport_ref Mean wind at which transport saturates (m/s).
#'   Default 15 (~54 km/h). Must exceed `wind_transport_onset`.
#' @param transport_max Maximum transport multiplier. Default 2.0.
#' @param rain_threshold Hourly precipitation at or above which all litter
#'   hazard is suppressed (mm). Default 0.5.
#'
#' @return Numeric vector of length equal to the inputs, the litter hazard index
#'   in `[0, 100]` for each forecast hour.
#'
#' @references
#' EPA (2006). AP-42, Section 13.2.5: Industrial Wind Erosion. Basis for the
#' friction-velocity excess entrainment form.
#'
#' Fecan, F., Marticorena, B. and Bergametti, G. (1999). Parametrization of the
#' increase of the aeolian erosion threshold wind friction velocity due to soil
#' moisture for arid and semi-arid areas. \emph{Annales Geophysicae}, 17,
#' 149--157. Basis for the moisture-raised threshold.
#'
#' Mellink, Y. et al. (2024). Wind- and rain-driven macroplastic mobilization
#' and transport on land. \emph{Scientific Reports} 14, 5006.
#' \doi{10.1038/s41598-024-53971-8}.
#'
#' @seealso [litter_exposure()] for the direction- and geometry-aware exposure
#'   layer that sits on top of this hazard index.
#' @export
litter_risk_index <- function(
  wind_gusts_10m,
  wind_speed_10m,
  precipitation,
  soil_moisture_0_to_1cm,
  kappa                = 0.40,
  z0                   = 0.05,
  ustar_t0             = 0.30,
  ustar_ref            = 1.05,
  entrainment_max      = 50,
  excess_exponent      = 2,
  moisture_gain        = 2.0,
  moisture_curve       = 0.5,
  soil_dry             = 0.05,
  soil_wet             = 0.20,
  wind_transport_onset = 5.5,
  wind_transport_ref   = 15.0,
  transport_max        = 2.0,
  rain_threshold       = 0.5
) {

  # ---- Validate meteorological inputs (complete, non-negative, aligned) ---- #
  n <- length(wind_gusts_10m)
  checkmate::assert_numeric(wind_gusts_10m, lower = 0, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(wind_speed_10m, lower = 0, any.missing = FALSE, len = n)
  checkmate::assert_numeric(precipitation, lower = 0, any.missing = FALSE, len = n)
  checkmate::assert_numeric(soil_moisture_0_to_1cm, lower = 0, upper = 1,
                            any.missing = FALSE, len = n)

  # ---- Validate parameters and the ordering constraints --------------------- #
  checkmate::assert_number(kappa, lower = .Machine$double.eps)
  checkmate::assert_number(z0, lower = .Machine$double.eps, upper = 10 - 1e-9)
  checkmate::assert_number(ustar_t0, lower = 0)
  checkmate::assert_number(ustar_ref)
  if (ustar_ref <= ustar_t0) {
    cli::cli_abort(c(
      "{.arg ustar_ref} ({ustar_ref}) must be greater than {.arg ustar_t0} ({ustar_t0}).",
      "i" = "The entrainment excess range (ustar_ref - ustar_t0) must be positive."
    ))
  }
  checkmate::assert_number(entrainment_max, lower = 0)
  checkmate::assert_number(excess_exponent, lower = 0)
  checkmate::assert_number(moisture_gain, lower = 0)
  checkmate::assert_number(moisture_curve, lower = 0)
  checkmate::assert_number(soil_dry, lower = 0)
  checkmate::assert_number(soil_wet, lower = 0)
  if (soil_dry >= soil_wet) {
    cli::cli_abort(c(
      "{.arg soil_dry} ({soil_dry}) must be less than {.arg soil_wet} ({soil_wet}).",
      "i" = "The moisture ramp denominator (soil_wet - soil_dry) must be positive."
    ))
  }
  checkmate::assert_number(wind_transport_onset, lower = 0)
  checkmate::assert_number(wind_transport_ref)
  if (wind_transport_ref <= wind_transport_onset) {
    cli::cli_abort(c(
      "{.arg wind_transport_ref} ({wind_transport_ref}) must be greater than {.arg wind_transport_onset} ({wind_transport_onset}).",
      "i" = "The transport ramp denominator must be positive."
    ))
  }
  checkmate::assert_number(transport_max, lower = 1)
  checkmate::assert_number(rain_threshold, lower = 0)

  # ---- Entrainment friction velocity from the gust ------------------------- #
  # Neutral logarithmic wind profile: u* = kappa * U / ln(z / z0), with the
  # 10 m gust in m/s. Using the gust (peak wind) as the effective driving wind
  # for entrainment follows AP-42's "fastest-mile" device. (specs/Litter_v3.md
  # S3.1)
  ln_factor <- kappa / log(10 / z0)
  ustar_g   <- ln_factor * wind_gusts_10m

  # ---- Moisture-raised entrainment threshold (Fecan et al. 1999) ----------- #
  # A damp surface raises u*t; a dry surface (SM <= soil_dry) sits at the base
  # threshold. (specs/Litter_v3.md S4.4)
  s       <- pmin(1, pmax(0, (soil_moisture_0_to_1cm - soil_dry) / (soil_wet - soil_dry)))
  ustar_t <- ustar_t0 * (1 + moisture_gain * s^moisture_curve)

  # ---- Entrainment score: bounded excess-power (AP-42 shape) ---------------- #
  # E = a * min(1, (max(0, u*g - u*t) / dustar)^n). The excess is floored at 0
  # before exponentiation so a fractional exponent never meets a negative base.
  delta_ustar <- ustar_ref - ustar_t0
  excess      <- pmax(0, ustar_g - ustar_t)
  entrainment <- entrainment_max * pmin(1, (excess / delta_ustar)^excess_exponent)

  # Saturation veto: a saturated surface releases no litter regardless of gust.
  entrainment[soil_moisture_0_to_1cm >= soil_wet] <- 0

  # ---- Transport potential from the mean wind ------------------------------ #
  # Linear ramp [1, transport_max] over the mean wind (m/s). Transport is
  # flight-height advection, so it is driven by the mean wind, not u*.
  # (specs/Litter_v3.md S4.3)
  transport <- 1 + (transport_max - 1) *
    pmin(1, pmax(0, wind_speed_10m - wind_transport_onset) /
           (wind_transport_ref - wind_transport_onset))

  # ---- Rainfall hard gate -------------------------------------------------- #
  rain_gate <- ifelse(precipitation >= rain_threshold, 0, 1)

  # ---- Composite, capped at 100 -------------------------------------------- #
  pmin(entrainment * transport * rain_gate, 100)
}


#' Litter hazard index for a forecast tibble
#'
#' Computes the hourly litter hazard index ([litter_risk_index()]) for each row
#' of a meteorological forecast tibble. Wraps the vector API, accepting a tibble
#' with named columns rather than individual vectors.
#'
#' @param met_data A tibble (or data frame) with one row per hourly timestep
#'   containing at least the columns `wind_gusts_10m` (m/s), `wind_speed_10m`
#'   (m/s), `precipitation` (mm), and `soil_moisture_0_to_1cm` (m^3/m^3).
#' @param ... Additional calibration parameters forwarded to
#'   [litter_risk_index()] (e.g. `rain_threshold`, `soil_wet`, `z0`).
#'
#' @return Numeric vector of length `nrow(met_data)`, the litter hazard index in
#'   `[0, 100]` for each forecast hour.
#'
#' @seealso [litter_risk_index()], [litter_exposure()].
#' @export
generate_litter_risk_index <- function(met_data, ...) {
  checkmate::assert_data_frame(met_data, min.rows = 1)

  required_cols <- c(
    "wind_gusts_10m", "wind_speed_10m", "precipitation", "soil_moisture_0_to_1cm"
  )
  .assert_required_cols(
    met_data, required_cols, arg = "met_data",
    info = paste0(
      "Required: wind_gusts_10m (m/s), wind_speed_10m (m/s), ",
      "precipitation (mm), soil_moisture_0_to_1cm (m\u00b3/m\u00b3)."
    )
  )

  litter_risk_index(
    wind_gusts_10m         = met_data$wind_gusts_10m,
    wind_speed_10m         = met_data$wind_speed_10m,
    precipitation          = met_data$precipitation,
    soil_moisture_0_to_1cm = met_data$soil_moisture_0_to_1cm,
    ...
  )
}
