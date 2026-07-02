# ===========================================================================
# Odour hazard layer (receptor-independent, direction-agnostic).
#
# See specs/Odour_v2.md for the full model. This file holds the
# odour_hazard() ventilation index (built on ventilation_state(), in
# R/odour-ventilation.R). The geometry-aware exposure layer lives in
# R/odour_exposure.R.
# ===========================================================================


#' Landfill odour hazard index (meteorological, receptor-independent)
#'
#' Computes an hourly **relative odour hazard** for a landfill site: given the
#' weather in each forecast hour, how strong is the potential for odour to
#' accumulate and reach the surroundings. This is a direction-agnostic,
#' receptor-independent point hazard, **not** a dispersion model -- wind
#' direction, receptor distance, and site geometry are handled separately by
#' [odour_exposure()].
#'
#' The function does not query any API. The caller fetches the required hourly
#' variables from the Open-Meteo `/v1/forecast` (or archive) endpoint and passes
#' them as a data frame, one row per consecutive hourly timestep.
#'
#' @section Units:
#' The dimensional columns (`wind_speed_10m`, `wind_speed_80m`,
#' `direct_radiation`, `boundary_layer_height`, `temperature_2m`, `pressure_msl`,
#' `precipitation`) may each be a bare numeric in the documented unit or a
#' \pkg{units} object, which is converted automatically (a dimensionally
#' incompatible unit is an error). The percentage / ratio columns (`cloud_cover`,
#' `relative_humidity_2m`, `soil_moisture_*`) and `wind_direction_10m` (degrees)
#' are taken as-is. The returned hazard is a dimensionless relative index (plain
#' numeric).
#'
#' @section Model:
#' The hazard is the **ventilation index** -- source emission strength divided
#' by the atmosphere's ability to ventilate it (wind speed times mixing depth):
#' \deqn{H = G \cdot PM(s) \cdot W_{rain} / (u_{eff} \cdot h_{mix}) / H_{ref}}
#' normalised so a calm, shallow, stable baseline gives \eqn{H \approx 1}. It
#' captures the dominant signal -- calm, stable, shallow boundary layer = high
#' hazard -- without baking in any receptor distance or direction. See
#' `specs/Odour_v2.md` for the full derivation.
#'
#' \describe{
#'   \item{\eqn{G}}{Source generation modifier
#'     \eqn{(1 + \Delta P_{mod})(1 + R_{mod})(1 + S_{seal})(1 + H_{mod})(1 + V_{mod})}
#'     (barometric pumping, post-rain piston, soil sealing, humidity, surface
#'     volatilisation), combined multiplicatively -- independent fractional
#'     modifiers compound rather than add -- range about 0.80 to 2.33.}
#'   \item{\eqn{PM(s)}}{Peak-to-mean ratio, rising with stability from
#'     `PM_MIN` (unstable) to `PM_MAX` (very stable); odour annoyance is driven
#'     by sub-minute peaks that hourly-average dispersion smooths out.}
#'   \item{\eqn{W_{rain}}}{Below-cloud scavenging, blended between no washout
#'     and the soluble-limit tiers by `odorant_solubility` (Henry's-law
#'     solubility, default 0.5 = mixed sulfur/soluble profile).}
#'   \item{\eqn{u_{eff}, h_{mix}, s}}{Effective wind, mixing depth, and stability
#'     from the shared dispersion state (see Details / `stability`). By default
#'     (`pool_cap = TRUE`) `h_mix` is capped at the nocturnal cold-pool depth on
#'     stable nights -- see [ventilation_state()].}
#' }
#'
#' @param met_data A data frame (or tibble), one row per consecutive hourly
#'   timestep, with numeric columns `wind_speed_10m` (m/s), `direct_radiation`
#'   (W/m^2), `cloud_cover` (\%), `boundary_layer_height` (m), `temperature_2m`
#'   (deg C), `pressure_msl` (hPa), `precipitation` (mm), `relative_humidity_2m`
#'   (\%), `soil_moisture_0_to_1cm` and `soil_moisture_1_to_3cm` (m^3/m^3).
#'   `NA` values are permitted and handled conservatively. Row order must be
#'   consecutive hourly (for the 3-hour pressure tendency and 24-hour rainfall
#'   lookback).
#' @param stability Stability estimator: `"turner"` (default, Pasquill-Turner
#'   from insolation/cloud and wind -- self-consistent with the dispersion
#'   curves) or `"shear"` (legacy 10 m/80 m power-law exponent, which also
#'   requires a `wind_speed_80m` column).
#'
#' @return A numeric vector of length `nrow(met_data)`: the relative odour
#'   hazard index for each hour (reference baseline = 1.0; not clamped).
#'
#'   **Scale note.** Unlike [litter_hazard()] and [dust_hazard()], which return
#'   bounded 0--100 indices, this is an *unbounded relative* index: a ventilation
#'   index has no natural ceiling, and the 0--100 framing is applied downstream
#'   by [odour_exposure()]. Whether to unify all hazards onto one scale is a
#'   deferred modelling/API decision (GitHub issue #11).
#'
#'   **`H approx 1` is no longer the worst baseline.** With the v3 nocturnal
#'   cold-pool cap (`pool_cap = TRUE`, default) an actively-trapped inversion
#'   night can push `h_mix` well below the calm-stable reference depth used to
#'   normalise `H`, so trapped nights now routinely exceed 1.0.
#'
#' @references
#' Turner, D.B. (1970). \emph{Workbook of Atmospheric Dispersion Estimates}.
#' US EPA AP-26.
#' Holzworth, G.C. (1972). \emph{Mixing Heights, Wind Speeds, and Potential for
#' Urban Air Pollution Throughout the Contiguous United States}. US EPA.
#' Czepiel, P.M. et al. (2003). The influence of atmospheric pressure on
#' landfill methane emissions. \emph{Waste Management}, 23(7), 593--598.
#' Zou, S.C. et al. (2003). Volatile organic compound emissions from landfills.
#' \emph{Atmospheric Environment}, 37(16), 2197--2211.
#'
#' @param datetime Optional `POSIXct` vector, one value per row. When supplied,
#'   the rows are checked for consecutive hourly spacing -- the 3-hour pressure
#'   tendency and 24-hour rainfall lookback assume it -- and a warning is issued
#'   if they are not (the computation proceeds on row order regardless).
#'
#' @seealso [odour_exposure()] for the geometry-aware exposure layer, and
#'   [odour_risk()] for the combined convenience wrapper.
#' @param terrain An [mh_terrain()] object, or `NULL` (default). When supplied,
#'   `shelter` and `shelter_h_mix` may modify the effective wind speed and
#'   mixing depth via M3 valley sheltering before the ventilation index is
#'   computed. Ignored when `shelter = FALSE`.
#' @param shelter Logical (default `FALSE`). If `TRUE` and `terrain` is
#'   non-`NULL` with a finite `shelter_index`, apply M3 valley sheltering
#'   to reduce the effective wind speed.
#' @param shelter_h_mix Logical (default `FALSE`). If `TRUE`, also reduce the
#'   mixing depth via M3 shelter.
#' @param pool_cap Logical (default `TRUE`). Passed to [ventilation_state()];
#'   caps `h_mix` at the nocturnal cold-pool depth on stable nights.
#' @param odorant_solubility Number in `[0, 1]` (default
#'   `ODOUR_CONSTANTS$ODORANT_SOLUBILITY_DEFAULT`, 0.5). Passed to
#'   [ventilation_state()]; blends `W_rain` between no washout (0) and the
#'   soluble-limit tiers (1).
#' @export
odour_hazard <- function(met_data, stability = c("turner", "shear"),
                         datetime = NULL,
                         terrain = NULL,
                         shelter = FALSE, shelter_h_mix = FALSE,
                         pool_cap = TRUE,
                         odorant_solubility = ODOUR_CONSTANTS$ODORANT_SOLUBILITY_DEFAULT) {
  stability <- match.arg(stability)
  .assert_hourly(datetime)

  required_cols <- c(
    "wind_speed_10m", "direct_radiation", "cloud_cover",
    "boundary_layer_height", "temperature_2m", "pressure_msl",
    "precipitation", "relative_humidity_2m",
    "soil_moisture_0_to_1cm", "soil_moisture_1_to_3cm"
  )

  checkmate::assert_data_frame(met_data, min.rows = 1)
  .assert_required_cols(met_data, required_cols, arg = "met_data",
                        info = "See {.code ?odour_hazard} for the required Open-Meteo columns.")
  .assert_numeric_cols(met_data, required_cols, arg = "met_data")

  met_data <- .odour_normalise_met(met_data)

  vs <- ventilation_state(met_data, terrain = terrain, stability = stability,
                          shelter = shelter, shelter_h_mix = shelter_h_mix,
                          pool_cap = pool_cap,
                          odorant_solubility = odorant_solubility)
  G  <- .odour_generation(met_data)

  hazard_ref <- ODOUR_CONSTANTS$PM_MAX /
    (ODOUR_CONSTANTS$U_CALM_FLOOR * ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE)

  .odour_hazard_raw(G, vs) / hazard_ref
}


# Raw ventilation flux: G * PM * W_rain / (u_eff * h_mix).
# Single source of truth consumed by odour_hazard() and odour_exposure() alike,
# ensuring both layers share an identical numerator.
.odour_hazard_raw <- function(G, vs) {
  G * vs$PM * vs$W_rain / (vs$u_eff * vs$h_mix)
}


# ---- Dimensional-column normalisation -------------------------------------- #
# Convert the dimensional met_data columns to canonical-unit plain doubles: a
# bare numeric is assumed already in the canonical unit, a units object is
# converted (a dimensional mismatch errors). Percentage / ratio columns
# (cloud_cover, relative_humidity_2m, soil_moisture_*) and wind_direction_10m
# (degrees) are dimensionless-by-convention here and left untouched. Columns
# that are absent are skipped (presence is validated by the caller).
.odour_normalise_met <- function(met_data) {
  canonical <- c(
    wind_speed_10m        = "m/s",
    wind_speed_80m        = "m/s",
    direct_radiation      = "W/m^2",
    boundary_layer_height = "m",
    temperature_2m        = "degree_C",
    pressure_msl          = "hPa",
    precipitation         = "mm"
  )
  for (col in names(canonical)) {
    if (!is.null(met_data[[col]])) {
      met_data[[col]] <- .drop_to(met_data[[col]], canonical[[col]], arg = col)
    }
  }
  met_data
}


# ---- Pasquill-Turner stability index --------------------------------------- #
# Turner (1964/1970) day/night table mapped to a numeric PG class (A = 0 ...
# F = 5, with half steps for the A-B / B-C / C-D pairings). Daytime insolation
# grade comes from the direct-radiation magnitude (strong >= 700, moderate
# >= 350, else slight W/m^2); nighttime branches on cloud cover (>= 50% ~ 4/8
# low cloud). Vectorised over the inputs.
.turner_stability <- function(u, rad, cloud, is_day) {
  grade  <- dplyr::case_when(rad >= 700 ~ 3, rad >= 350 ~ 2, TRUE ~ 1)
  cloudy <- cloud >= 50

  day_s <- dplyr::case_when(
    u < 2 & grade == 3 ~ 0.0,
    u < 2 & grade == 2 ~ 0.5,
    u < 2              ~ 1.0,
    u < 3 & grade == 3 ~ 0.5,
    u < 3 & grade == 2 ~ 1.0,
    u < 3              ~ 2.0,
    u < 5 & grade == 3 ~ 1.0,
    u < 5 & grade == 2 ~ 1.5,
    u < 5              ~ 2.0,
    u < 6 & grade == 3 ~ 2.0,
    u < 6 & grade == 2 ~ 2.5,
    u < 6              ~ 3.0,
    grade == 3         ~ 2.0,
    grade == 2         ~ 3.0,
    TRUE               ~ 3.0
  )

  night_s <- dplyr::case_when(
    u < 2          ~ 5.0,
    u < 3 & cloudy ~ 4.0,
    u < 3          ~ 5.0,
    u < 5 & cloudy ~ 3.0,
    u < 5          ~ 4.0,
    TRUE           ~ 3.0
  )

  ifelse(is_day, day_s, night_s)
}


# ---- Shear power-law exponent -> continuous stability ---------------------- #
# Maps alpha = ln(u80/u10)/ln(8) to s in [0, 5] (Irwin 1979; Counihan 1975).
# Legacy estimator, retained as the optional `stability = "shear"` override.
.alpha_to_s <- function(alpha) {
  dplyr::case_when(
    alpha <= 0.07 ~ 0.0,
    alpha <= 0.10 ~ (alpha - 0.07) / 0.03,
    alpha <= 0.13 ~ 1.0 + (alpha - 0.10) / 0.03,
    alpha <= 0.15 ~ 2.0 + (alpha - 0.13) / 0.02,
    alpha <= 0.22 ~ 3.0 + (alpha - 0.15) / 0.07,
    alpha <= 0.40 ~ 4.0 + (alpha - 0.22) / 0.18,
    TRUE          ~ 5.0
  )
}


# ---- Consecutive-hourly spacing guard -------------------------------------- #
# The pressure-tendency and rainfall lookbacks index by row, assuming each row
# is the next consecutive hour. When the caller supplies a datetime, warn (do
# not abort) if the rows are not regularly spaced one hour apart.
.assert_hourly <- function(datetime) {
  if (is.null(datetime)) return(invisible())
  if (!inherits(datetime, "POSIXct")) {
    cli::cli_abort(
      "{.arg datetime} must be a {.cls POSIXct} vector, not {.cls {class(datetime)}}.",
      class = "meteoHazard_input_error"
    )
  }
  if (length(datetime) > 1) {
    gaps <- as.numeric(diff(datetime), units = "secs")
    if (any(is.na(gaps)) || any(abs(gaps - 3600) > 1)) {
      cli::cli_warn(c(
        "{.arg datetime} is not consecutive hourly; row order is used regardless.",
        "i" = "The 3-hour pressure tendency and 24-hour rainfall lookback assume one row per hour."
      ))
    }
  }
  invisible()
}
