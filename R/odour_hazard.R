# ===========================================================================
# Odour hazard layer (receptor-independent, direction-agnostic).
#
# See specs/Odour_v2.md for the full model. This file holds the shared
# dispersion-state helper and the odour_hazard() ventilation index. The
# geometry-aware exposure layer lives in R/odour_exposure.R.
# ===========================================================================


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
