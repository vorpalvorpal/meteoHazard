# ===========================================================================
# Ventilation state (C2, issue #15): source-independent atmospheric state
# used by odour_hazard(), odour_exposure(), and future terrain-aware layers.
#
# ventilation_state() replaces .odour_dispersion_state() and extends it with
# peak-to-mean (PM), below-cloud scavenging (W_rain), nocturnal cold-pool
# depth (pool_top), convective boundary-layer growth rate (cbl_growth), and
# overnight residual winds at multiple levels.
#
# .odour_generation() lifts the G computation out of odour_hazard() so the
# two components can be tested and extended independently.
# ===========================================================================


#' Atmospheric ventilation state for the odour dispersion model
#'
#' Computes the per-hour atmospheric state used by [odour_hazard()] and
#' [odour_exposure()], extended with the nocturnal cold-pool depth
#' (`pool_top`), convective boundary-layer growth rate (`cbl_growth`), and
#' overnight residual winds at 80 m / 120 m / 180 m AGL.
#'
#' @param met_data A data frame, one row per consecutive hour. Required columns:
#'   `wind_speed_10m` (m/s), `direct_radiation` (W/m^2), `cloud_cover` (\%),
#'   `boundary_layer_height` (m). Optional columns: `temperature_2m` (°C) and
#'   `relative_humidity_2m` (\%) for the Brunt cooling rate; multi-level winds
#'   `wind_speed_80m`, `wind_direction_80m`, `wind_speed_120m`,
#'   `wind_direction_120m`, `wind_speed_180m`, `wind_direction_180m` for
#'   `residual_wind`. Missing optional columns yield `NA` in the corresponding
#'   output fields (not an error). Malformed optional columns (present but
#'   non-numeric) are an error.
#' @param terrain An [mh_terrain()] object, or `NULL` (default, flat site).
#'   When non-`NULL`, `terrain@valley_depth` caps the cold-pool depth and
#'   `terrain@taf` amplifies the Brunt cooling rate before accumulation.
#' @param stability Stability estimator: `"turner"` (default) or `"shear"`.
#'   `"shear"` requires `wind_speed_80m`.
#' @param shelter Logical (default `FALSE`). When `TRUE` and `terrain` has a
#'   finite `shelter_index`, reduces `u_eff` using the M3 valley-sheltering
#'   transfer model. **Off by default until calibrated** (see #8). `FALSE`
#'   gives bit-identical output to the pre-C6 baseline.
#' @param shelter_h_mix Logical (default `FALSE`). When `TRUE` (and
#'   `shelter = TRUE`), also scales `h_mix` by `(1 - reduction_effective)`.
#'   Undocumented interaction with the forecast `boundary_layer_height`
#'   (potential double-count); leave `FALSE` until reconciled.
#'
#' @return A named list, all vectors of length `nrow(met_data)`:
#' \describe{
#'   \item{`u_eff`}{Effective wind speed (m/s), floored at `U_CALM_FLOOR`.}
#'   \item{`h_mix`}{Mixing depth (m), positive; falls back to stable/unstable
#'     constant when `boundary_layer_height` is `NA` or non-positive.}
#'   \item{`s`}{Pasquill-Gifford stability index in \[0, 5\] (A = 0, F = 5).}
#'   \item{`is_calm`}{Logical; wind below `U_CALM_FLOOR` or `NA`.}
#'   \item{`is_day`}{Logical; direct radiation > 10 W/m^2.}
#'   \item{`PM`}{Peak-to-mean ratio, interpolated between `PM_MIN` (s = 0)
#'     and `PM_MAX` (s = 5).}
#'   \item{`W_rain`}{Below-cloud scavenging factor in (0, 1]; requires
#'     `precipitation` in `met_data` (falls back to 1.0 when absent/NA).}
#'   \item{`pool_top`}{Nocturnal cold-pool depth (m). Accumulated nightly from
#'     Brunt net-longwave cooling then Whiteman (1999) saturation depth, floored
#'     by the Venkatram (1980) mechanical depth. Frozen at morning release and
#'     reset at next nightfall. `NA` when temperature and RH are unavailable.}
#'   \item{`cbl_growth`}{Convective boundary-layer growth rate (m/h); positive
#'     during morning transition, zero otherwise.}
#'   \item{`residual_wind`}{Named list with elements `speed_80m`, `dir_80m`,
#'     `speed_120m`, `dir_120m`, `speed_180m`, `dir_180m`, `speed_10m`,
#'     `dir_10m`: per-hour overnight circular-mean direction (degrees) and
#'     concurrent speed at each level. The 10 m surface-wind level is a last-
#'     resort fallback used when all upper-level columns are absent; it
#'     represents the surface-layer wind inside the boundary layer rather than
#'     the decoupled residual layer above. `NA` for absent or all-daytime
#'     windows.}
#' }
#'
#' @seealso [odour_hazard()], [odour_exposure()]
#' @references
#' Brunt, D. (1932). Notes on radiation in the atmosphere. \emph{Quarterly
#'   Journal of the Royal Meteorological Society}, 58(247), 389--420.
#' Whiteman, C.D. (1999). Wintertime evolution of the temperature inversion
#'   in the Sinbad Valley. \emph{Journal of Applied Meteorology}, 38, 178--188.
#' Venkatram, A. (1980). Estimating the Monin-Obukhov length in the stable
#'   boundary layer for dispersion calculations. \emph{Boundary-Layer
#'   Meteorology}, 19, 481--485.
#'
#' @section Terrain-modifier contract:
#' When multiple terrain mechanisms are active simultaneously, the following
#' precedence applies: **M1 drainage confinement** (C3b) supersedes **M3
#' valley sheltering** on any hour when `drainage_active` is `TRUE`
#' (`is_channelled & !is_day & pool_top > 0`). The suppression strength is
#' controlled by `ODOUR_CONSTANTS$DRAINAGE_SHELTER_OVERLAP` (default 1.0 =
#' full mutual exclusion). **M2 receptor impaction** is orthogonal (a vertical
#' geometry term, not a wind-field modifier) and never participates in
#' precedence. A future external wind field (C7) will supersede all native
#' terms.
#'
#' The M3 hazard/exposure re-entanglement caveat: enabling `shelter = TRUE`
#' makes `ventilation_state()` depend on a terrain descriptor, so the
#' ventilation state is no longer purely meteorological. This is inherent to
#' M3; keep `shelter = FALSE` in all meteorological-analysis contexts.
#'
#' @export
ventilation_state <- function(met_data, terrain = NULL,
                              stability = c("turner", "shear"),
                              shelter      = FALSE,
                              shelter_h_mix = FALSE) {
  stability <- match.arg(stability)

  # ---- Validation ---------------------------------------------------------- #
  checkmate::assert_data_frame(met_data, min.rows = 1)

  required_cols <- c(
    "wind_speed_10m", "direct_radiation", "cloud_cover", "boundary_layer_height"
  )
  .assert_required_cols(met_data, required_cols, arg = "met_data",
                        info = "See {.code ?ventilation_state} for the required Open-Meteo columns.")
  .assert_numeric_cols(met_data, required_cols, arg = "met_data")

  # Optional multi-level wind columns: absent is fine (NA output); present but
  # non-numeric is an error.
  ml_cols <- c(
    "wind_speed_80m", "wind_direction_80m",
    "wind_speed_120m", "wind_direction_120m",
    "wind_speed_180m", "wind_direction_180m"
  )
  ml_present <- ml_cols[ml_cols %in% names(met_data)]
  if (length(ml_present) > 0L) {
    .assert_numeric_cols(met_data, ml_present, arg = "met_data")
  }

  met_data <- .odour_normalise_met(met_data)

  if (!is.logical(shelter) || length(shelter) != 1L || is.na(shelter))
    cli::cli_abort("{.arg shelter} must be logical (TRUE or FALSE).",
                   class = "meteoHazard_input_error")
  if (!is.logical(shelter_h_mix) || length(shelter_h_mix) != 1L || is.na(shelter_h_mix))
    cli::cli_abort("{.arg shelter_h_mix} must be logical (TRUE or FALSE).",
                   class = "meteoHazard_input_error")

  # ---- Core dispersion state (replaces .odour_dispersion_state) ------------ #
  u10   <- met_data$wind_speed_10m
  rad   <- met_data$direct_radiation
  cloud <- met_data$cloud_cover
  bl    <- met_data$boundary_layer_height

  u10_safe   <- .na_fill(u10, 0)
  rad_safe   <- .na_fill(rad, 0)
  cloud_safe <- .na_fill(cloud, 50)

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
    is.na(bl) | bl <= 0,
    ifelse(is_calm | s >= 4, ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE,
           ODOUR_CONSTANTS$H_MIX_FALLBACK_UNSTABLE),
    bl
  )

  # ---- PM and W_rain ------------------------------------------------------- #
  PM <- ODOUR_CONSTANTS$PM_MIN +
    (ODOUR_CONSTANTS$PM_MAX - ODOUR_CONSTANTS$PM_MIN) * (s / 5)

  if (!is.null(met_data$precipitation)) {
    precip_safe <- .na_fill(met_data$precipitation, 0)
    W_rain <- dplyr::case_when(
      precip_safe > 4.0 ~ 0.05,
      precip_safe > 1.0 ~ 0.15,
      precip_safe > 0.2 ~ 0.40,
      TRUE              ~ 1.0
    )
  } else {
    W_rain <- rep(1.0, nrow(met_data))
  }

  # ---- pool_top ------------------------------------------------------------ #
  pool_top <- .pool_top(met_data, u10_safe, is_day, terrain)

  # ---- cbl_growth ---------------------------------------------------------- #
  cbl_growth <- pmax(0, c(0, diff(h_mix)))

  # ---- Residual winds ------------------------------------------------------ #
  residual_wind <- .residual_wind(met_data, is_day)

  # ---- M3 valley sheltering ------------------------------------------------ #
  if (isTRUE(shelter)) {
    sheltered <- .valley_shelter(u_eff, h_mix, u10_safe, is_day, pool_top,
                                 terrain, shelter_h_mix)
    u_eff <- sheltered$u_eff
    h_mix <- sheltered$h_mix
  }

  list(
    u_eff         = u_eff,
    h_mix         = h_mix,
    s             = s,
    is_calm       = is_calm,
    is_day        = is_day,
    PM            = PM,
    W_rain        = W_rain,
    pool_top      = pool_top,
    cbl_growth    = cbl_growth,
    residual_wind = residual_wind
  )
}


# ---- Nocturnal cold-pool depth --------------------------------------------- #
# Returns pool_top (m) for each hour: the Brunt/Whiteman nocturnal heat-deficit
# depth, floored by the Venkatram mechanical depth. Accumulated from night-start
# for each hour, then frozen at daily maximum until next nightfall.
#
# Returns NA when temperature and RH are both absent (Brunt cooling requires
# them; mechanical floor still uses u10).
.pool_top <- function(met_data, u10_safe, is_day, terrain) {
  n_t   <- nrow(met_data)
  temp  <- met_data$temperature_2m
  rh    <- met_data$relative_humidity_2m
  cloud <- .na_fill(met_data$cloud_cover, 50)

  is_dark <- !is_day

  # TAF from terrain (amplifies cooling rate; >= 1).
  taf <- if (!is.null(terrain) && !is.na(terrain@taf)) terrain@taf else 1.0

  # Basin-sill cap.
  valley_cap <- if (!is.null(terrain) && !is.na(terrain@valley_depth)) {
    terrain@valley_depth
  } else {
    Inf
  }

  # Venkatram mechanical floor per hour: h_floor = C * u_star^1.5
  kappa  <- ODOUR_CONSTANTS$POOL_KAPPA
  z0     <- ODOUR_CONSTANTS$POOL_Z0
  C_vent <- ODOUR_CONSTANTS$VENKATRAM_COEF
  u_star <- kappa * u10_safe / log(10 / z0)
  h_floor <- dplyr::case_when(
    u10_safe <= 0 ~ ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE,
    TRUE          ~ C_vent * u_star^1.5
  )

  # Brunt net-longwave cooling rate (W/m^2) — meaningful only in dark hours.
  # Q*[i] = sigma * T_K^4 * (a - b*sqrt(e)) * (1 - k*C/100)
  # where e = RH/100 * 6.112 * exp(17.67*T / (T+243.5)) [hPa]
  sigma <- 5.67e-8
  a_b   <- ODOUR_CONSTANTS$BRUNT_A
  b_b   <- ODOUR_CONSTANTS$BRUNT_B
  k_b   <- ODOUR_CONSTANTS$BRUNT_CLOUD_K

  have_t_rh <- !is.null(temp) && !is.null(rh) &&
    any(!is.na(temp)) && any(!is.na(rh))

  if (!have_t_rh) {
    # No Brunt cooling; fall back to mechanical floor everywhere.
    pool_raw <- pmax(h_floor, ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE)
    pool_raw <- pmin(pool_raw, valley_cap)
    return(pool_raw)
  }

  T_safe  <- .na_fill(temp, 15)
  rh_safe <- .na_fill(rh,   50)
  T_K     <- T_safe + 273.15
  e_vp    <- rh_safe / 100 * 6.112 * exp(17.67 * T_safe / (T_safe + 243.5))
  Qstar   <- sigma * T_K^4 * (a_b - b_b * sqrt(pmax(e_vp, 0))) *
    (1 - k_b * cloud / 100)

  # Q* only accumulates in the dark (is_dark). Outside, force to 0.
  Qstar_night <- ifelse(is_dark, pmax(Qstar, 0), 0)

  # For each hour t, scan back to the start of the current night and sum
  # taf * Q*[i] * 3600 J/m^2.  Then take running max since night start.
  Q_accum    <- numeric(n_t)
  pool_accum <- numeric(n_t)

  for (t in seq_len(n_t)) {
    if (!is_dark[t]) {
      Q_accum[t]    <- 0
      pool_accum[t] <- 0
      next
    }
    # Find start of the current night: scan backward to the most recent
    # is_day hour (night start is the hour after that).
    night_start <- t
    if (t > 1L) {
      for (i in seq(t - 1L, 1L)) {
        if (!is_dark[i]) {
          night_start <- i + 1L
          break
        }
        night_start <- i
      }
    }

    Q_accum[t] <- taf * sum(Qstar_night[night_start:t]) * 3600

    # Whiteman saturation depth.
    H_SAT <- ODOUR_CONSTANTS$POOL_H_SAT
    Q_SAT <- ODOUR_CONSTANTS$POOL_Q_SAT
    h_def  <- H_SAT * (1 - exp(-Q_accum[t] / Q_SAT))

    # Running maximum since night start (pool grows, never shrinks overnight).
    prev_max <- if (t > 1L) pool_accum[t - 1L] else 0
    h_with_floor <- max(h_def, h_floor[t])
    pool_accum[t] <- max(h_with_floor, if (is_dark[t] && t > 1L && is_dark[t - 1L]) prev_max else 0)
  }

  # Freeze overnight maximum into daytime hours (pool stays until next night).
  pool_top <- numeric(n_t)
  last_night_max <- 0
  for (t in seq_len(n_t)) {
    if (is_dark[t]) {
      last_night_max <- pool_accum[t]
      pool_top[t]    <- pool_accum[t]
    } else {
      pool_top[t] <- last_night_max
    }
  }

  pmin(pool_top, valley_cap)
}


# ---- Overnight circular-mean residual winds -------------------------------- #
# Returns a list of speed/direction vectors for each available level (80m, 120m,
# 180m). Direction is the circular mean of overnight (is_day == FALSE) hours
# in the current night window; speed is the concurrent arithmetic mean.
.residual_wind <- function(met_data, is_day) {
  n_t     <- nrow(met_data)
  is_dark <- !is_day

  levels <- list(
    list(speed = "wind_speed_80m",  dir = "wind_direction_80m"),
    list(speed = "wind_speed_120m", dir = "wind_direction_120m"),
    list(speed = "wind_speed_180m", dir = "wind_direction_180m"),
    list(speed = "wind_speed_10m",  dir = "wind_direction_10m")
  )
  out_names <- c("speed_80m", "dir_80m",
                 "speed_120m", "dir_120m",
                 "speed_180m", "dir_180m",
                 "speed_10m",  "dir_10m")

  result <- vector("list", 8L)
  names(result) <- out_names

  for (lev_i in seq_along(levels)) {
    spd_col <- levels[[lev_i]]$speed
    dir_col <- levels[[lev_i]]$dir
    sp_out  <- out_names[2L * lev_i - 1L]
    dr_out  <- out_names[2L * lev_i]

    if (!(spd_col %in% names(met_data)) || !(dir_col %in% names(met_data))) {
      result[[sp_out]] <- rep(NA_real_, n_t)
      result[[dr_out]] <- rep(NA_real_, n_t)
      next
    }

    spd_v <- met_data[[spd_col]]
    dir_v <- met_data[[dir_col]]

    spd_out <- numeric(n_t)
    dir_out <- numeric(n_t)

    for (t in seq_len(n_t)) {
      # Current night window: from start of this night to t.
      if (!is_dark[t]) {
        spd_out[t] <- NA_real_
        dir_out[t] <- NA_real_
        next
      }
      night_start <- t
      if (t > 1L) {
        for (i in seq(t - 1L, 1L)) {
          if (!is_dark[i]) {
            night_start <- i + 1L
            break
          }
          night_start <- i
        }
      }
      idx    <- night_start:t
      d_dark <- dir_v[idx]
      s_dark <- spd_v[idx]
      valid  <- !is.na(d_dark) & !is.na(s_dark)
      if (!any(valid)) {
        spd_out[t] <- NA_real_
        dir_out[t] <- NA_real_
      } else {
        theta    <- d_dark[valid] * pi / 180
        circ_dir <- atan2(mean(sin(theta)), mean(cos(theta))) * 180 / pi
        dir_out[t] <- circ_dir %% 360
        spd_out[t] <- mean(s_dark[valid])
      }
    }

    result[[sp_out]] <- spd_out
    result[[dr_out]] <- dir_out
  }

  result
}


# ---- Source generation modifier -------------------------------------------- #
# Lifts the G computation from odour_hazard() verbatim. Returns a numeric
# vector of length nrow(met_data). Requires: pressure_msl, precipitation,
# temperature_2m, relative_humidity_2m, soil_moisture_0_to_1cm,
# soil_moisture_1_to_3cm.
.odour_generation <- function(met_data) {
  n_t         <- nrow(met_data)
  pressure    <- met_data$pressure_msl
  precip_safe <- .na_fill(met_data$precipitation, 0)
  temp        <- met_data$temperature_2m
  rh_safe     <- .na_fill(met_data$relative_humidity_2m, 0)
  sm01_safe   <- .na_fill(met_data$soil_moisture_0_to_1cm, 0)
  sm13_safe   <- .na_fill(met_data$soil_moisture_1_to_3cm, 0)

  # Barometric pumping: falling pressure increases advective gas flux.
  dP3 <- pressure - dplyr::lag(pressure, 3)
  dP_mod <- dplyr::case_when(
    is.na(dP3) ~ 0.0,
    dP3 <= -5  ~ 0.30,
    dP3 < 0    ~ -0.06 * dP3,
    TRUE       ~ 0.0
  )

  # Post-rain piston effect; active-rain guard MUST be the first branch.
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

  1.0 + dP_mod + R_mod + S_seal + H_mod + V_mod
}


# ---- M3 valley sheltering transfer model ---------------------------------- #
# Modifies u_eff (and optionally h_mix) based on terrain shelter_index.
# All args are length-n_t vectors or scalars. Returns list(u_eff, h_mix).
#
# Transfer model (all constants from ODOUR_CONSTANTS; uncalibrated → #8):
#   s_f = clamp((SHELTER_OPEN_REF - shelter_index) /
#               (SHELTER_OPEN_REF - SHELTER_ENCLOSED_REF), 0, 1)
#   w_r = clamp((SHELTER_U_FLUSH - u10) /
#               (SHELTER_U_FLUSH - SHELTER_U_FULL), 0, 1)
#   reduction = SHELTER_MAX_REDUCTION * s_f * w_r
#
# Precedence: on drainage-active hours (channelled + night + pool > 0),
# reduction is suppressed by DRAINAGE_SHELTER_OVERLAP (default 1 = full).
#   drainage_active = is_channelled & !is_day & !is.na(pool_top) & pool_top > 0
#   reduction_effective = reduction * (1 - DRAINAGE_SHELTER_OVERLAP * drainage_active)
#   u_eff_sheltered = max(u_eff * (1 - reduction_effective), U_CALM_FLOOR)
.valley_shelter <- function(u_eff, h_mix, u10, is_day, pool_top, terrain, shelter_h_mix) {
  K <- ODOUR_CONSTANTS

  # No-op when terrain is NULL or shelter_index is NA.
  si <- if (!is.null(terrain) && !is.na(terrain@shelter_index)) terrain@shelter_index else NA_real_
  if (is.na(si)) return(list(u_eff = u_eff, h_mix = h_mix))

  # Shelter strength (from openness angle).
  s_f <- pmax(0, pmin(1, (K$SHELTER_OPEN_REF - si) /
                          (K$SHELTER_OPEN_REF - K$SHELTER_ENCLOSED_REF)))

  # Wind-regime taper (clamp to [0,1]).
  w_r <- pmax(0, pmin(1, (K$SHELTER_U_FLUSH - u10) /
                          (K$SHELTER_U_FLUSH - K$SHELTER_U_FULL)))

  # Base reduction.
  reduction <- K$SHELTER_MAX_REDUCTION * s_f * w_r

  # Drainage-confinement precedence: suppress shelter on drainage-active hours.
  is_channelled <- !is.null(terrain) &&
    !is.na(terrain@flow_convergence) &&
    !is.na(terrain@drainage_bearing) &&
    terrain@flow_convergence >= 0.5
  drainage_active <- is_channelled &
    !is_day &
    !is.na(pool_top) &
    pool_top > 0
  reduction_effective <- reduction * (1 - K$DRAINAGE_SHELTER_OVERLAP * as.numeric(drainage_active))

  u_eff_new <- pmax(u_eff * (1 - reduction_effective), K$U_CALM_FLOOR)
  h_mix_new <- if (isTRUE(shelter_h_mix)) {
    pmax(h_mix * (1 - reduction_effective), 1)
  } else {
    h_mix
  }

  list(u_eff = u_eff_new, h_mix = h_mix_new)
}
