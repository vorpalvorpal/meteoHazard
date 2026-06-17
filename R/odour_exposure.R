# ===========================================================================
# Odour exposure layer (geometry- and direction-aware).
#
# Maps the receptor-independent odour hazard (R/odour_hazard.R) onto a
# per-receptor operational consequence using a Pasquill-Gifford Gaussian-plume
# geometry, returning a 0-100 worst-case band. The optional terrain-aware
# drainage / morning-fumigation module is relocated here from the monolith;
# see the revisit note on .odour_drainage_state().
#
# See specs/Odour_v2.md for the model.
# ===========================================================================


#' Landfill odour exposure (direction- and geometry-aware)
#'
#' Maps the direction-agnostic odour hazard ([odour_hazard()]) onto an
#' operational consequence given the wind direction and the receptor geometry.
#' Where the hazard answers *how strong is the odour situation around the site*,
#' the exposure layer answers *who is downwind and how bad is it for them*.
#'
#' For each hour and receptor the full Pasquill-Gifford Gaussian-plume relative
#' concentration is reconstructed from the hazard (the `h_mix` in the hazard's
#' ventilation term cancels), giving distance decay and a stability- and
#' distance-aware directional Gaussian. Forecast wind-direction uncertainty is
#' added as a separate convolution term. The relative concentration is mapped to
#' a 0-100 band and the worst-case across receptors is returned.
#'
#' @param hazard Numeric vector, the relative odour hazard from [odour_hazard()],
#'   one value per forecast hour.
#' @param met_data The same data frame passed to [odour_hazard()] (must contain
#'   `wind_direction_10m` plus the columns the dispersion state needs:
#'   `wind_speed_10m`, `direct_radiation`, `cloud_cover`,
#'   `boundary_layer_height`). `nrow(met_data)` must equal `length(hazard)`.
#' @param receptors A data frame, one row per receptor, with `bearing` (degrees
#'   from north, 0-360, landfill to receptor) and `distance` (m, strictly
#'   positive).
#' @param drainage_axes Optional data frame enabling the terrain-aware katabatic
#'   drainage and morning fumigation refinement (`bearing_from` 0-360, `weight`
#'   positive). When `NULL` (default) that logic is skipped.
#' @param stability Stability estimator passed to the shared dispersion state;
#'   see [odour_hazard()].
#' @param map_c50 Relative concentration at which the 0-100 map reaches ~63
#'   (`exposure = 100 * (1 - exp(-C_rel / map_c50))`). Provisional, calibratable
#'   operational knob; default 0.3.
#'
#' @return A numeric vector of length `length(hazard)`: the worst-case 0-100
#'   odour exposure across receptors for each hour.
#'
#' @references
#' Briggs, G.A. (1973). \emph{Diffusion Estimation for Small Emissions}. NOAA.
#' Schauberger, G. et al. (2002). Calculating direction-dependent separation
#' distance by a dispersion model. \emph{Biosystems Engineering}, 82(1), 25--37.
#' Whiteman, C.D. (2000). \emph{Mountain Meteorology}. Oxford University Press.
#'
#' @seealso [odour_hazard()] for the upstream hazard index, and
#'   [generate_odour_risk_index()] for the combined wrapper.
#' @export
odour_exposure <- function(hazard, met_data, receptors, drainage_axes = NULL,
                           stability = c("turner", "shear"), map_c50 = 0.3) {
  stability <- match.arg(stability)

  # ---- Validation ---------------------------------------------------------- #
  checkmate::assert_numeric(hazard, lower = 0, any.missing = FALSE, min.len = 1)
  checkmate::assert_data_frame(met_data, min.rows = 1)
  if (nrow(met_data) != length(hazard)) {
    cli::cli_abort(
      "{.arg met_data} has {nrow(met_data)} row{?s} but {.arg hazard} has length {length(hazard)}.",
      class = "meteoHazard_input_error"
    )
  }
  required_cols <- c(
    "wind_direction_10m", "wind_speed_10m", "direct_radiation",
    "cloud_cover", "boundary_layer_height"
  )
  missing_cols <- setdiff(required_cols, names(met_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg met_data} is missing required columns: {.val {missing_cols}}.",
      class = "meteoHazard_input_error"
    )
  }
  checkmate::assert_data_frame(receptors, min.rows = 1)
  if (!all(c("bearing", "distance") %in% names(receptors))) {
    cli::cli_abort(
      "{.arg receptors} must contain columns {.val bearing} and {.val distance}.",
      class = "meteoHazard_input_error"
    )
  }
  checkmate::assert_numeric(receptors$bearing, lower = 0, upper = 360,
                            any.missing = FALSE, .var.name = "receptors$bearing")
  checkmate::assert_numeric(receptors$distance, lower = .Machine$double.eps,
                            any.missing = FALSE, .var.name = "receptors$distance")
  if (!is.null(drainage_axes)) {
    checkmate::assert_data_frame(drainage_axes, min.rows = 1)
    if (!all(c("bearing_from", "weight") %in% names(drainage_axes))) {
      cli::cli_abort(
        "{.arg drainage_axes} must contain columns {.val bearing_from} and {.val weight}.",
        class = "meteoHazard_input_error"
      )
    }
    checkmate::assert_numeric(drainage_axes$bearing_from, lower = 0, upper = 360,
                              any.missing = FALSE, .var.name = "drainage_axes$bearing_from")
    checkmate::assert_numeric(drainage_axes$weight, lower = .Machine$double.eps,
                              any.missing = FALSE, .var.name = "drainage_axes$weight")
  }
  checkmate::assert_number(map_c50, lower = .Machine$double.eps)

  # ---- Shared dispersion state and drainage/fumigation state --------------- #
  state <- .odour_dispersion_state(met_data, stability)
  n_t   <- nrow(met_data)
  n_r   <- nrow(receptors)

  drn <- if (!is.null(drainage_axes)) {
    .odour_drainage_state(met_data, state, receptors, drainage_axes)
  } else {
    list(is_drainage = rep(FALSE, n_t), is_fumigation = rep(FALSE, n_t),
         effective_severity = rep(0, n_t), max_alignment = rep(0, n_r))
  }

  wind_dir   <- met_data$wind_direction_10m
  theta_down <- (ifelse(is.na(wind_dir), NA_real_, wind_dir) + 180) %% 360
  sigma_fc   <- ODOUR_CONSTANTS$SIGMA_FC_DEG * pi / 180
  c_y        <- ODOUR_CONSTANTS$SIGMA_Y_COEF

  s_low  <- floor(state$s)
  s_high <- pmin(s_low + 1, 5)
  frac   <- state$s - s_low

  CALM_CROSSWIND <- 0.5 # omnidirectional meander factor under calm/NA direction

  # ---- Per-receptor relative concentration -> 0-100 ------------------------ #
  risk_matrix <- matrix(0.0, nrow = n_t, ncol = n_r)

  for (j in seq_len(n_r)) {
    x_j <- receptors$distance[j]

    # Briggs sigma_y / sigma_z at this distance for all 6 PG classes.
    sigma_y_classes <- c_y * x_j / sqrt(1 + 0.0001 * x_j)
    sigma_z_classes <- c(
      0.20  * x_j,
      0.12  * x_j,
      0.08  * x_j / sqrt(1 + 0.0002 * x_j),
      0.06  * x_j / sqrt(1 + 0.0015 * x_j),
      0.03  * x_j / (1 + 0.0003 * x_j),
      0.016 * x_j / (1 + 0.0003 * x_j)
    )

    sigma_y_t <- sigma_y_classes[s_low + 1] +
      frac * (sigma_y_classes[s_high + 1] - sigma_y_classes[s_low + 1])
    sigma_z_t <- sigma_z_classes[s_low + 1] +
      frac * (sigma_z_classes[s_high + 1] - sigma_z_classes[s_low + 1])

    # Effective lateral spread = physical width (+) forecast-direction uncertainty.
    sigma_y_eff <- sqrt(sigma_y_t^2 + (x_j * sigma_fc)^2)

    # Dilution / distance factor (the hazard's h_mix cancels here).
    geom_base <- state$h_mix / (sigma_y_eff * pmin(sigma_z_t, state$h_mix))

    # Crosswind / direction factor: stability- and distance-aware Gaussian on
    # the crosswind offset y = x * sin(dTheta).
    theta_j     <- receptors$bearing[j]
    diff_raw    <- abs(theta_down - theta_j)
    delta_theta <- pmin(diff_raw, 360 - diff_raw)
    y_cross     <- x_j * sin(delta_theta * pi / 180)
    cw <- exp(-0.5 * (y_cross / sigma_y_eff)^2)
    # Upwind hemisphere (delta_theta > 90): the receptor is behind the source,
    # not in the plume. The crosswind offset y = x*sin(theta) is symmetric about
    # 90 deg, so without this gate an upwind receptor would read as fully
    # exposed. At 90 deg the Gaussian is already negligible, so this is smooth.
    cw[delta_theta > 90] <- 0
    cw[state$is_calm | is.na(wind_dir)] <- CALM_CROSSWIND

    # Drainage override: katabatic flow confines emissions in the hollow,
    # directing them away from upslope receptors (0.05 on-axis to 0.10 off).
    if (any(drn$is_drainage)) {
      cw[drn$is_drainage] <- 0.05 + 0.05 * (1 - drn$max_alignment[j])
    }
    # Fumigation floor: anabatic flow lofts the overnight pool toward aligned
    # receptors (ceiling 0.6, weaker than synoptic transport).
    if (any(drn$is_fumigation)) {
      floor_fum <- 0.6 * drn$effective_severity * drn$max_alignment[j]
      cw[drn$is_fumigation] <- pmax(cw[drn$is_fumigation], floor_fum[drn$is_fumigation])
    }

    c_rel <- hazard * geom_base * cw

    # Fumigation source-entrainment boost (was a G multiplier in the monolith;
    # applied here because G now lives in the receptor-independent hazard).
    if (any(drn$is_fumigation)) {
      c_rel[drn$is_fumigation] <- c_rel[drn$is_fumigation] *
        (1 + 0.5 * drn$effective_severity[drn$is_fumigation])
    }

    risk_matrix[, j] <- 100 * (1 - exp(-c_rel / map_c50))
  }

  apply(risk_matrix, 1, max)
}


# ---- Drainage / morning-fumigation state (relocated, optional) ------------- #
# revisit: this is the single-site terrain heuristic carried over unchanged in
# logic from the monolithic generate_odour_risk_index(). It is a bespoke state
# machine with uncalibrated magic numbers, kept to prompt future generalisation
# (a generic nocturnal-accumulation / morning-release pulse) and calibration.
# It now rides the Pasquill-Turner stability index s (so its s > 4.0 thresholds
# behave slightly differently than under the old shear-based s).
#
# Returns, per hour: is_drainage, is_fumigation, effective_severity; and per
# receptor: max_alignment (time-invariant drainage-axis alignment).
.odour_drainage_state <- function(met_data, state, receptors, drainage_axes) {
  n_t <- nrow(met_data)
  n_r <- nrow(receptors)

  is_calm     <- state$is_calm
  s_t         <- state$s
  h_effective <- state$h_mix
  rad_safe    <- ifelse(is.na(met_data$direct_radiation), 0, met_data$direct_radiation)
  cloud_safe  <- ifelse(is.na(met_data$cloud_cover), 50, met_data$cloud_cover)

  is_drainage          <- rep(FALSE, n_t)
  is_fumigation        <- rep(FALSE, n_t)
  effective_severity_v <- rep(0.0, n_t)
  max_alignment        <- rep(0.0, n_r)

  # Per-receptor drainage-axis alignment (time-invariant).
  for (j in seq_len(n_r)) {
    alignments <- vapply(seq_len(nrow(drainage_axes)), function(k) {
      delta_k <- abs(receptors$bearing[j] - drainage_axes$bearing_from[k])
      delta_k <- min(delta_k, 360 - delta_k)
      if (delta_k <= 30) {
        cos(pi * delta_k / 60)^2 * drainage_axes$weight[k]
      } else {
        0.0
      }
    }, numeric(1))
    max_alignment[j] <- min(1.0, max(alignments))
  }

  drainage_hours           <- integer(n_t)
  latched_drainage_dur     <- integer(n_t)
  overnight_stable_hours   <- integer(n_t)
  latched_stable_dur       <- integer(n_t)
  fumigation_pool_consumed <- logical(n_t)

  fumigation_event_start <- NA_integer_
  fum_severity_at_start  <- 0.0
  latch_consumed_flag    <- FALSE

  for (t in seq_len(n_t)) {
    is_dark_t      <- rad_safe[t] < 10
    is_new_night_t <- is_dark_t && (t == 1L || rad_safe[t - 1L] >= 10)

    is_drainage[t] <- is_calm[t] && is_dark_t && cloud_safe[t] < 70

    drainage_hours[t] <- if (t == 1L) {
      if (is_drainage[t]) 1L else 0L
    } else if (is_drainage[t]) {
      drainage_hours[t - 1L] + 1L
    } else {
      0L
    }

    is_stable_t  <- s_t[t] > 4.0
    is_shallow_t <- h_effective[t] < 300
    overnight_stable_hours[t] <- if (is_stable_t && is_shallow_t && is_dark_t) {
      if (t == 1L) 1L else overnight_stable_hours[t - 1L] + 1L
    } else {
      0L
    }

    if (is_new_night_t || latch_consumed_flag) {
      latched_drainage_dur[t] <- 0L
      latched_stable_dur[t]   <- 0L
      latch_consumed_flag     <- FALSE
    } else if (t == 1L) {
      latched_drainage_dur[t] <- 0L
      latched_stable_dur[t]   <- 0L
    } else {
      if (drainage_hours[t - 1L] > 0L && drainage_hours[t] == 0L) {
        latched_drainage_dur[t] <- drainage_hours[t - 1L]
      } else if (drainage_hours[t] > 0L) {
        latched_drainage_dur[t] <- 0L
      } else {
        latched_drainage_dur[t] <- latched_drainage_dur[t - 1L]
      }
      if (overnight_stable_hours[t - 1L] > 0L && overnight_stable_hours[t] == 0L) {
        latched_stable_dur[t] <- overnight_stable_hours[t - 1L]
      } else if (overnight_stable_hours[t] > 0L) {
        latched_stable_dur[t] <- 0L
      } else {
        latched_stable_dur[t] <- latched_stable_dur[t - 1L]
      }
    }

    if (is_dark_t) {
      fumigation_pool_consumed[t] <- FALSE
    } else if (t == 1L) {
      fumigation_pool_consumed[t] <- FALSE
    } else {
      fumigation_pool_consumed[t] <- fumigation_pool_consumed[t - 1L]
    }

    preceding_idx   <- if (t > 1L) max(1L, t - 6L):(t - 1L) else integer(0)
    has_recent_dark <- length(preceding_idx) > 0 && any(rad_safe[preceding_idx] < 10)
    activation_criteria <- rad_safe[t] > 50 && has_recent_dark

    if (!is.na(fumigation_event_start)) {
      if ((t - fumigation_event_start) >= 3L || is_dark_t) {
        fumigation_event_start <- NA_integer_
        fum_severity_at_start  <- 0.0
      }
    }

    if (activation_criteria && !fumigation_pool_consumed[t] &&
        is.na(fumigation_event_start)) {
      ldd <- latched_drainage_dur[t]
      lsd <- latched_stable_dur[t]
      accumulation_hours <- if (ldd > 0L) as.numeric(ldd) else as.numeric(lsd) * 0.5
      fum_severity_at_start       <- min(1.0, accumulation_hours / 8.0)
      fumigation_event_start      <- t
      fumigation_pool_consumed[t] <- TRUE
      latch_consumed_flag         <- TRUE
    }

    is_fumigation[t] <- !is.na(fumigation_event_start)

    if (is_fumigation[t]) {
      hours_since_start       <- t - fumigation_event_start
      decay                   <- max(0.0, 1.0 - hours_since_start / 3.0)
      effective_severity_v[t] <- fum_severity_at_start * decay
    }
  }

  list(is_drainage = is_drainage, is_fumigation = is_fumigation,
       effective_severity = effective_severity_v, max_alignment = max_alignment)
}
