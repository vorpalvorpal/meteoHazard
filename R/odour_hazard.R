# ===========================================================================
# Odour hazard layer (receptor-independent, direction-agnostic).
#
# See specs/Odour_v2.md for the full model. This file holds the shared
# dispersion-state helper and the odour_hazard() ventilation index. The
# geometry-aware exposure layer lives in R/odour_exposure.R.
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
#'     \eqn{1 + \Delta P_{mod} + R_{mod} + S_{seal} + H_{mod} + V_{mod}}
#'     (barometric pumping, post-rain piston, soil sealing, humidity, surface
#'     volatilisation), range ~[0.80, 1.95].}
#'   \item{\eqn{PM(s)}}{Peak-to-mean ratio, rising with stability from
#'     `PM_MIN` (unstable) to `PM_MAX` (very stable); odour annoyance is driven
#'     by sub-minute peaks that hourly-average dispersion smooths out.}
#'   \item{\eqn{W_{rain}}}{Below-cloud scavenging of soluble odorants.}
#'   \item{\eqn{u_{eff}, h_{mix}, s}}{Effective wind, mixing depth, and stability
#'     from the shared dispersion state (see Details / `stability`).}
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
#'   [generate_odour_risk_index()] for the combined convenience wrapper.
#' @export
odour_hazard <- function(met_data, stability = c("turner", "shear"),
                         datetime = NULL) {
  stability <- match.arg(stability)
  .assert_hourly(datetime)

  required_cols <- c(
    "wind_speed_10m", "direct_radiation", "cloud_cover",
    "boundary_layer_height", "temperature_2m", "pressure_msl",
    "precipitation", "relative_humidity_2m",
    "soil_moisture_0_to_1cm", "soil_moisture_1_to_3cm"
  )

  checkmate::assert_data_frame(met_data, min.rows = 1)
  missing_cols <- setdiff(required_cols, names(met_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      c(
        "{.arg met_data} is missing required columns: {.val {missing_cols}}.",
        "i" = "See {.code ?odour_hazard} for the required Open-Meteo columns."
      ),
      class = "meteoHazard_input_error"
    )
  }
  for (col in required_cols) {
    if (!is.numeric(met_data[[col]])) {
      cli::cli_abort(
        "{.arg met_data} column {.val {col}} must be numeric, not {.cls {class(met_data[[col]])}}.",
        class = "meteoHazard_input_error"
      )
    }
  }

  state <- .odour_dispersion_state(met_data, stability)
  n_t   <- nrow(met_data)

  # ---- G: source generation modifier (additive, labelled) ----------------- #
  pressure    <- met_data$pressure_msl
  precip_safe <- ifelse(is.na(met_data$precipitation), 0, met_data$precipitation)
  temp        <- met_data$temperature_2m
  rh_safe     <- ifelse(is.na(met_data$relative_humidity_2m), 0, met_data$relative_humidity_2m)
  sm01_safe   <- ifelse(is.na(met_data$soil_moisture_0_to_1cm), 0, met_data$soil_moisture_0_to_1cm)
  sm13_safe   <- ifelse(is.na(met_data$soil_moisture_1_to_3cm), 0, met_data$soil_moisture_1_to_3cm)

  # Barometric pumping: falling pressure increases advective gas flux.
  dP3 <- pressure - dplyr::lag(pressure, 3)
  dP_mod <- dplyr::case_when(
    is.na(dP3) ~ 0.0,
    dP3 <= -5  ~ 0.30,
    dP3 < 0    ~ -0.06 * dP3,
    TRUE       ~ 0.0
  )

  # Post-rain piston effect; active-rain guard MUST be the first branch
  # (P_24 includes currently-falling rain). P_24[i] sums the 24 preceding rows
  # via a cumulative-sum window: cs[i] - cs[max(1, i-24)].
  cs   <- cumsum(c(0, precip_safe))
  idx  <- seq_len(n_t)
  P_24 <- cs[idx] - cs[pmax(1L, idx - 24L)]
  R_mod <- dplyr::case_when(
    precip_safe > 0.5 ~ 0.0,
    P_24 > 15         ~ 0.20,
    P_24 > 5          ~ 0.10,
    TRUE              ~ 0.0
  )

  # Soil-moisture cover sealing (wettest layer is the diffusion bottleneck).
  sm_seal <- pmax(sm01_safe, sm13_safe)
  S_seal <- dplyr::case_when(
    sm_seal >= 0.40 ~ -0.20,
    sm_seal >= 0.25 ~ -0.20 * (sm_seal - 0.25) / 0.15,
    TRUE            ~ 0.0
  )

  H_mod <- dplyr::case_when(
    rh_safe >= 85 ~ 0.15,
    rh_safe >= 60 ~ 0.15 * (rh_safe - 60) / 25,
    TRUE          ~ 0.0
  )

  # Surface NMOC volatilisation (Henry's law); ceiling widened to V_MOD_MAX.
  vmax <- ODOUR_CONSTANTS$V_MOD_MAX
  V_mod <- dplyr::case_when(
    is.na(temp) ~ 0.0,
    temp <= 10  ~ 0.0,
    temp >= 35  ~ vmax,
    TRUE        ~ vmax * (temp - 10) / 25
  )

  G <- 1.0 + dP_mod + R_mod + S_seal + H_mod + V_mod

  # ---- Peak-to-mean and scavenging overlays -------------------------------- #
  PM <- ODOUR_CONSTANTS$PM_MIN +
    (ODOUR_CONSTANTS$PM_MAX - ODOUR_CONSTANTS$PM_MIN) * (state$s / 5)

  W_rain <- dplyr::case_when(
    precip_safe > 4.0 ~ 0.05,
    precip_safe > 1.0 ~ 0.15,
    precip_safe > 0.2 ~ 0.40,
    TRUE              ~ 1.0
  )

  # ---- Ventilation index, normalised to the calm/stable/shallow baseline --- #
  hazard_raw <- G * PM * W_rain / (state$u_eff * state$h_mix)
  hazard_ref <- ODOUR_CONSTANTS$PM_MAX /
    (ODOUR_CONSTANTS$U_CALM_FLOOR * ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE)

  hazard_raw / hazard_ref
}


# ---- Shared dispersion state ----------------------------------------------- #
# Computes the per-hour atmospheric state used by BOTH odour_hazard() and
# odour_exposure(): the Pasquill-Gifford stability index s in [0, 5] (A = 0,
# F = 5), the effective wind u_eff (calm-floored), the mixing depth h_mix, and
# the is_calm / is_day flags. `stability` selects the estimator:
#   "turner" (default) - Pasquill-Turner from insolation/cloud + wind. Primary,
#       because the Briggs sigma curves are tabulated by the PG class that the
#       Turner scheme defines (self-consistent).
#   "shear"  - legacy 10 m/80 m power-law exponent (needs wind_speed_80m).
.odour_dispersion_state <- function(met_data, stability = c("turner", "shear")) {
  stability <- match.arg(stability)

  u10   <- met_data$wind_speed_10m
  rad   <- met_data$direct_radiation
  cloud <- met_data$cloud_cover
  bl    <- met_data$boundary_layer_height

  u10_safe   <- ifelse(is.na(u10), 0, u10)
  rad_safe   <- ifelse(is.na(rad), 0, rad)
  cloud_safe <- ifelse(is.na(cloud), 50, cloud)

  is_calm <- is.na(u10) | u10 < ODOUR_CONSTANTS$U_CALM_FLOOR
  is_day  <- rad_safe > 10

  if (stability == "turner") {
    s_raw <- .turner_stability(u10_safe, rad_safe, cloud_safe, is_day)
    s <- ifelse(is_calm, 4.25, s_raw)
  } else {
    u80 <- met_data$wind_speed_80m
    use_calm_stab <- is_calm | is.na(u80) | u80 <= 0
    u10_for_alpha <- pmax(u10_safe, 0.001)
    u80_for_alpha <- ifelse(is.na(u80) | u80 <= 0, 0.001, u80)
    alpha <- log(u80_for_alpha / u10_for_alpha) / log(8)
    s <- ifelse(use_calm_stab, 4.25, .alpha_to_s(alpha))
  }

  u_eff <- pmax(u10_safe, ODOUR_CONSTANTS$U_CALM_FLOOR)

  h_mix <- ifelse(
    is.na(bl),
    ifelse(is_calm | s >= 4, ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE,
           ODOUR_CONSTANTS$H_MIX_FALLBACK_UNSTABLE),
    bl
  )

  list(s = s, u_eff = u_eff, h_mix = h_mix, is_calm = is_calm, is_day = is_day)
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
