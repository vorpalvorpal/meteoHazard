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

  # Receptor coordinates and the per-class sigma_z formula coefficients are
  # source-independent, so extract them once.
  rec_coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(receptors)))

  # Terrain-backend setup (loop-invariant). The descriptors path adds a
  # per-source terrain delta reusing the shared flat dispersion matrices.
  descriptors <- terrain_backend == "descriptors" && !is.null(site@terrain)
  if (descriptors) {
    terrain           <- site@terrain
    pool_top_for_part <- ifelse(is.na(vs$pool_top), 0, vs$pool_top)
    delta_t           <- pmax(ODOUR_CONSTANTS$DELTA_FLOOR,
                              ODOUR_CONSTANTS$DELTA_FRAC * pool_top_for_part)
    vs_aug            <- vs
    vs_aug$wind_dir_surface <- wind_dir
    pool_na   <- is.na(vs$pool_top)   # length n_t
    is_calm_t <- vs$is_calm           # length n_t
    # Morning release is source-independent (depends only on the pool history).
    r <- .morning_release(pool_top_for_part, vs$cbl_growth, vs$is_day)
  }

  # Per-hour generation/ventilation scalar (length n_t), broadcast down columns.
  # Kept separate from hazard_ref so the final division order matches the
  # original per-receptor arithmetic exactly (G*PM*W_rain*geom*cw / hazard_ref).
  ghw <- G * vs$PM * vs$W_rain

  # ---- Accumulate summed concentration over sources ----------------------- #
  # All per-(hour, receptor) quantities are held as (n_t x n_r) matrices so the
  # arithmetic dispatches once per source instead of once per receptor.
  c_sum_matrix     <- matrix(0.0, nrow = n_t, ncol = n_r)
  c_terrain_matrix <- if (descriptors) matrix(0.0, n_t, n_r) else NULL

  for (k in seq_len(n_s)) {
    source_k  <- sources[k, ]
    source_pt <- sf::st_centroid(sf::st_geometry(source_k))

    # Area-source initial spread sigma_y0 (time-varying, depends on wind dir).
    geom_type_k <- as.character(
      sf::st_geometry_type(source_k, by_geometry = FALSE)
    )
    if (grepl("POLYGON", geom_type_k, ignore.case = TRUE)) {
      wd_safe    <- ifelse(is.na(wind_dir), 0.0, wind_dir)
      sigma_y0_t <- .crosswind_halfwidth(source_k, wd_safe) /
                      ODOUR_CONSTANTS$ISC3_SIGMA_Y0_COEF
      sigma_y0_t[is.na(wind_dir)] <- 0.0
    } else {
      sigma_y0_t <- rep(0.0, n_t)
    }

    emit_ht <- if ("emit_height" %in% names(source_k)) source_k$emit_height[[1]] else 0.0
    if (is.null(emit_ht) || is.na(emit_ht)) emit_ht <- 0.0
    sigma_z0 <- emit_ht / ODOUR_CONSTANTS$ISC3_SIGMA_Z0_COEF

    # Bearing and distance from this source to every receptor (length n_r).
    src_xy <- sf::st_coordinates(source_pt)
    dE     <- rec_coords[, "X"] - src_xy[1L, "X"]
    dN     <- rec_coords[, "Y"] - src_xy[1L, "Y"]
    x      <- sqrt(dE^2 + dN^2)
    brng   <- (atan2(dE, dN) * 180 / pi) %% 360
    dead   <- !is.finite(x) | x == 0   # coincident receptors → 0 contribution
    brng[dead] <- NA_real_

    # ---- Shared flat dispersion as (n_t x n_r) matrices ------------------- #
    # Briggs class spreads as (6 x n_r): one column per receptor distance.
    SY6 <- c_y %o% (x / sqrt(1 + 0.0001 * x))                 # sigma_y, 6 x n_r
    SZ6 <- rbind(                                             # sigma_z, 6 x n_r
      0.20  * x,
      0.12  * x,
      0.08  * x / sqrt(1 + 0.0002 * x),
      0.06  * x / sqrt(1 + 0.0015 * x),
      0.03  * x / (1 + 0.0003 * x),
      0.016 * x / (1 + 0.0003 * x)
    )
    # Interpolate between bounding stability classes: row-index by the per-hour
    # class (length n_t) replicates rows → (n_t x n_r).
    SY_lo <- SY6[s_low + 1L, , drop = FALSE]
    SY_hi <- SY6[s_high + 1L, , drop = FALSE]
    SZ_lo <- SZ6[s_low + 1L, , drop = FALSE]
    SZ_hi <- SZ6[s_high + 1L, , drop = FALSE]
    sigma_y_t <- SY_lo + frac * (SY_hi - SY_lo)               # frac down cols
    sigma_z_t <- SZ_lo + frac * (SZ_hi - SZ_lo)

    Xmat        <- matrix(x, n_t, n_r, byrow = TRUE)          # distance per col
    sigma_y_eff <- sqrt(sigma_y_t^2 + sigma_y0_t^2 + (Xmat * sigma_fc)^2)
    sigma_z_eff <- sqrt(sigma_z_t^2 + sigma_z0^2)
    geom_base   <- vs$h_mix / (sigma_y_eff * pmin(sigma_z_eff, vs$h_mix))

    delta_theta <- .angular_diff(matrix(theta_down, n_t, n_r),
                                 matrix(brng, n_t, n_r, byrow = TRUE))
    cw_flat     <- exp(-0.5 * (Xmat * sin(delta_theta * pi / 180) / sigma_y_eff)^2)
    cw_flat[delta_theta > 90] <- 0
    cw_flat[vs$is_calm | is.na(wind_dir), ] <- CALM_CROSSWIND

    # Flat relative concentration; coincident receptors contribute 0.
    c_rel_flat <- ghw * geom_base * cw_flat / hazard_ref
    c_rel_flat[, dead] <- 0
    c_sum_matrix <- c_sum_matrix + c_rel_flat

    # ---- Descriptors terrain delta (reuses geom_base / cw_flat) ----------- #
    if (descriptors) {
      part <- .pool_partition(emit_ht, pool_top_for_part, delta_t)
      f_1a <- part$f_1a
      f_1b <- part$f_1b
      r_scale <- ifelse(f_1b > 0, r / pmax(pool_top_for_part, 1), 0)

      above_pool_ht <- max(0, emit_ht -
                             stats::median(pool_top_for_part, na.rm = TRUE) / 2)
      fumic_prep <- .cw_fumigation_prep(vs_aug, above_pool_ht = above_pool_ht)

      cw_1a <- .cw_venting_matrix(brng, vs_aug, terrain)              # n_t x n_r
      cw_1b <- .cw_fumigation_matrix(brng, vs_aug, terrain,
                                     prep = fumic_prep)               # n_t x n_r

      # cw_1a_safe: fall back to flat cw when calm, pool NA, or pathway NA.
      cond_a <- matrix(is_calm_t | pool_na, n_t, n_r) | is.na(cw_1a)
      cw_1a_safe <- cw_1a
      cw_1a_safe[cond_a] <- cw_flat[cond_a]

      # cw_1b_safe: NA → 0 (no fumigation outside morning), then flat when calm
      # or pool NA.
      cw_1b_safe <- cw_1b
      cw_1b_safe[is.na(cw_1b_safe)] <- 0
      rows_fb <- is_calm_t | pool_na
      cw_1b_safe[rows_fb, ] <- cw_flat[rows_fb, ]

      # Blended crosswind (f_1a, f_1b, r_scale are per-hour, broadcast down cols).
      cw_blended <- f_1a * cw_1a_safe + f_1b * cw_1b_safe * (1 + r_scale)

      c_rel_terrain <- ghw * geom_base * (cw_blended - cw_flat) / hazard_ref
      c_rel_terrain[, dead] <- 0
      c_terrain_matrix <- c_terrain_matrix + c_rel_terrain
    }
  }

  # ---- Reduction: sum_k -> map -> max_j ----------------------------------- #
  c_total     <- if (descriptors) c_sum_matrix + c_terrain_matrix else c_sum_matrix
  risk_matrix <- 100 * (1 - exp(-c_total / map_c50))
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
