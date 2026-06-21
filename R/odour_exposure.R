# ===========================================================================
# Odour exposure layer (C3a): mh_site + area-source ISC3 + multi-source sum.
#
# odour_exposure() takes an mh_site (sources + receptors), met_data, and
# ventilation_state() output.  For each hour it sums relative concentrations
# over all sources, maps the summed value to 0-100, then returns the
# worst-case receptor for the hour.
#
# terrain_backend = "none"  : flat Gaussian (implemented here, C3a).
# terrain_backend = "descriptors": wired in C3b; C3a just returns the flat
#   result unchanged.
#
# See specs/Odour_C3a.md for the model.
# ===========================================================================


#' Landfill odour exposure (direction- and geometry-aware, mh_site API)
#'
#' Maps the atmospheric ventilation state onto an operational consequence given
#' the wind direction, the source geometry, and the receptor layout described
#' by an [mh_site()] object. For each hour the relative concentrations from all
#' `(odour, source)` features are summed at each `(odour, receptor)`, then
#' mapped to a 0-100 band and reduced to the worst-case receptor.
#'
#' Area sources use ISC3 initial spreads: `sigma_y0 = crosswind_halfwidth /
#' 4.3` and `sigma_z0 = emit_height / 2.15` (both added in quadrature with the
#' Briggs physical spread). Point sources have `sigma_y0 = sigma_z0 = 0`.
#'
#' @param met_data A data frame, one row per consecutive hourly timestep. Must
#'   contain at minimum: `wind_direction_10m`, `wind_speed_10m`,
#'   `direct_radiation`, `cloud_cover`, `boundary_layer_height`. Columns for
#'   the generation modifier `G` (`pressure_msl`, `precipitation`,
#'   `temperature_2m`, `relative_humidity_2m`, `soil_moisture_0_to_1cm`,
#'   `soil_moisture_1_to_3cm`) are optional; if absent `G` defaults to 1.
#' @param site An [mh_site()] object with at least one `(odour, source)` role
#'   and one `(odour, receptor)` role.
#' @param stability Stability estimator: `"turner"` (default) or `"shear"`.
#' @param map_c50 Relative concentration at which the 0-100 map reaches ~63
#'   (`100 * (1 - exp(-C_rel / map_c50))`). Default 0.3.
#' @param terrain_backend `"none"` (flat Gaussian, C3a) or `"descriptors"`
#'   (terrain-aware morning pulse, wired in C3b; C3a returns the flat result).
#'
#' @return A plain numeric vector of length `nrow(met_data)`: the worst-case
#'   0-100 odour exposure across receptors for each hour.
#'
#' @section Units:
#' Dimensional `met_data` columns (see [odour_hazard()]) may be supplied as
#' bare numerics (documented unit) or \pkg{units} objects (converted
#' automatically; mismatch errors). Bearings and angles are in degrees and
#' taken as-is.
#'
#' @references
#' Briggs, G.A. (1973). \emph{Diffusion Estimation for Small Emissions}. NOAA.
#' US EPA (1995). \emph{User's Guide for the Industrial Source Complex (ISC3)
#'   Dispersion Models}.
#'
#' @seealso [odour_hazard()], [odour_risk()], [ventilation_state()], [mh_site()]
#' @export
odour_exposure <- function(met_data, site,
                           stability = c("turner", "shear"),
                           map_c50 = 0.3,
                           terrain_backend = c("none", "descriptors")) {
  stability       <- match.arg(stability)
  terrain_backend <- match.arg(terrain_backend)

  # ---- Validate met_data -------------------------------------------------- #
  checkmate::assert_data_frame(met_data, min.rows = 1)
  required_cols <- c(
    "wind_direction_10m", "wind_speed_10m", "direct_radiation",
    "cloud_cover", "boundary_layer_height"
  )
  .assert_required_cols(met_data, required_cols, arg = "met_data")
  met_data <- .odour_normalise_met(met_data)

  # ---- Validate site ------------------------------------------------------ #
  if (!S7::S7_inherits(site, mh_site)) {
    cli::cli_abort(
      "{.arg site} must be an {.cls mh_site} object.",
      class = "meteoHazard_input_error"
    )
  }
  checkmate::assert_number(map_c50, lower = .Machine$double.eps)

  # ---- Extract sources and receptors ------------------------------------- #
  sources   <- .role_features(site, "odour", "source")
  receptors <- .role_features(site, "odour", "receptor")

  if (nrow(sources) == 0L) {
    cli::cli_abort(
      "{.arg site} has no {.val (odour, source)} role.",
      class = "meteoHazard_input_error"
    )
  }
  if (nrow(receptors) == 0L) {
    cli::cli_abort(
      "{.arg site} has no {.val (odour, receptor)} role.",
      class = "meteoHazard_input_error"
    )
  }

  # ---- Ventilation state -------------------------------------------------- #
  vs <- ventilation_state(met_data, terrain = site@terrain, stability = stability)

  # ---- Generation modifier G ---------------------------------------------- #
  # Falls back to 1.0 per row if required columns are absent.
  G <- tryCatch(
    .odour_generation(met_data),
    error = function(e) rep(1.0, nrow(met_data))
  )

  # ---- Setup ------------------------------------------------------------- #
  n_t      <- nrow(met_data)
  n_s      <- nrow(sources)
  n_r      <- nrow(receptors)
  wind_dir <- met_data$wind_direction_10m

  hazard_ref <- ODOUR_CONSTANTS$PM_MAX /
    (ODOUR_CONSTANTS$U_CALM_FLOOR * ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE)

  sigma_fc <- ODOUR_CONSTANTS$SIGMA_FC_DEG * pi / 180
  c_y      <- ODOUR_CONSTANTS$SIGMA_Y_COEF

  s_low  <- floor(vs$s)
  s_high <- pmin(s_low + 1L, 5L)
  frac   <- vs$s - s_low

  CALM_CROSSWIND <- 0.5

  # Pre-compute downwind bearing for all hours.
  theta_down <- .downwind_bearing(wind_dir)

  # ---- Accumulate summed concentration over sources ----------------------- #
  # c_sum_matrix[t, j]: sum over sources k of c_rel[t, j, k]
  c_sum_matrix <- matrix(0.0, nrow = n_t, ncol = n_r)

  for (k in seq_len(n_s)) {
    source_k <- sources[k, ]

    # Use centroid for bearing/distance from source to receptor.
    source_pt <- sf::st_centroid(sf::st_geometry(source_k))

    # Area-source initial spread sigma_y0 (time-varying, depends on wind dir).
    geom_type_k <- as.character(
      sf::st_geometry_type(source_k, by_geometry = FALSE)
    )
    if (grepl("POLYGON", geom_type_k, ignore.case = TRUE)) {
      sigma_y0_t <- vapply(wind_dir, function(wd) {
        if (is.na(wd)) 0.0
        else .crosswind_halfwidth(source_k, wd) / ODOUR_CONSTANTS$ISC3_SIGMA_Y0_COEF
      }, numeric(1))
    } else {
      sigma_y0_t <- rep(0.0, n_t)
    }

    # Area-source initial spread sigma_z0 from emit_height (scalar, time-const).
    emit_ht <- if ("emit_height" %in% names(source_k)) source_k$emit_height[[1]] else 0.0
    if (is.null(emit_ht) || is.na(emit_ht)) emit_ht <- 0.0
    sigma_z0 <- emit_ht / ODOUR_CONSTANTS$ISC3_SIGMA_Z0_COEF

    for (j in seq_len(n_r)) {
      receptor_j <- receptors[j, ]

      bd     <- .bearing_distance(source_pt, receptor_j)
      x_jk   <- bd$distance
      theta_jk <- bd$bearing

      # Skip coincident source/receptor: contribute 0 (avoid division by zero).
      if (is.na(x_jk) || x_jk == 0) next

      # Briggs sigma_y / sigma_z for all 6 PG classes at distance x_jk.
      sigma_y_classes <- c_y * x_jk / sqrt(1 + 0.0001 * x_jk)
      sigma_z_classes <- c(
        0.20  * x_jk,
        0.12  * x_jk,
        0.08  * x_jk / sqrt(1 + 0.0002 * x_jk),
        0.06  * x_jk / sqrt(1 + 0.0015 * x_jk),
        0.03  * x_jk / (1 + 0.0003 * x_jk),
        0.016 * x_jk / (1 + 0.0003 * x_jk)
      )

      # Interpolate between bounding stability classes (time-varying s).
      sigma_y_t <- sigma_y_classes[s_low + 1L] +
        frac * (sigma_y_classes[s_high + 1L] - sigma_y_classes[s_low + 1L])
      sigma_z_t <- sigma_z_classes[s_low + 1L] +
        frac * (sigma_z_classes[s_high + 1L] - sigma_z_classes[s_low + 1L])

      # ISC3 area-source: add initial spreads in quadrature.
      sigma_y_eff_t <- sqrt(sigma_y_t^2 + sigma_y0_t^2 + (x_jk * sigma_fc)^2)
      sigma_z_eff_t <- sqrt(sigma_z_t^2 + sigma_z0^2)

      # Geometry factor (h_mix cancels the hazard-denominator h_mix).
      geom_base_t <- vs$h_mix /
        (sigma_y_eff_t * pmin(sigma_z_eff_t, vs$h_mix))

      # Crosswind / directional factor.
      delta_theta <- .angular_diff(theta_down, theta_jk)
      y_cross     <- x_jk * sin(delta_theta * pi / 180)
      cw <- exp(-0.5 * (y_cross / sigma_y_eff_t)^2)

      # Upwind hemisphere: receptor is behind the source.
      cw[delta_theta > 90] <- 0
      # Calm / missing wind direction: omnidirectional meander.
      cw[vs$is_calm | is.na(wind_dir)] <- CALM_CROSSWIND

      # Relative concentration for this source-receptor pair.
      c_rel_jk <- G * vs$PM * vs$W_rain * geom_base_t * cw / hazard_ref

      c_sum_matrix[, j] <- c_sum_matrix[, j] + c_rel_jk
    }
  }

  # ---- Reduction: sum_k -> map -> max_j ----------------------------------- #
  risk_matrix <- 100 * (1 - exp(-c_sum_matrix / map_c50))
  apply(risk_matrix, 1, max)
}


# ---- Internal helpers (kept here; not exported) ---------------------------- #

# Returns the bearing the wind is blowing TOWARD (downwind bearing).
# wind_direction_10m is the bearing FROM which wind comes (met convention).
.downwind_bearing <- function(wind_dir) {
  (wind_dir + 180) %% 360
}

# Smallest unsigned angular difference between two bearings (0-180 deg).
.angular_diff <- function(theta1, theta2) {
  d <- abs(theta1 - theta2) %% 360
  ifelse(d > 180, 360 - d, d)
}
