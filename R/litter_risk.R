#' Landfill litter risk (combined hazard + exposure)
#'
#' Convenience wrapper that computes the windblown-litter hazard
#' ([litter_hazard()]) and maps it through the geometry-aware
#' [litter_exposure()] layer in one call, so the caller does not have to keep the
#' hazard vector and `wind_direction_10m` aligned by hand. Mirrors
#' [odour_risk()].
#'
#' This is the entry point for the common case; use [litter_exposure()] directly
#' if you already have a hazard vector. Mapping the relative exposure onto an
#' operational index is a site-specific calibration step (issues #11/#8).
#'
#' @param met_data A tibble (or data frame), one row per hourly timestep. Must
#'   contain `wind_direction_10m` (degrees) plus the columns [litter_hazard()]
#'   requires for the chosen `use_wetness_state`. When `reach_per_ms` is
#'   supplied it must also contain `wind_speed_10m` (m/s), forwarded as the
#'   refined-reach mean wind.
#' @param site An [`mh_site`] S7 object; see [litter_exposure()].
#' @param use_wetness_state Logical (default `FALSE`). Passed to
#'   [litter_hazard()]: when `TRUE`, the entrainment wetness signal is the
#'   drying litter-wetness state ([litter_wetness()]) instead of soil moisture.
#' @param direction_tol,p_open_min,move_threshold,offsite_threshold,default_permeability
#'   Exposure-layer parameters, passed to [litter_exposure()].
#' @param reach_per_ms Optional positive scalar (metres per m/s). When supplied,
#'   [litter_exposure()] runs in refined distance-reach mode with
#'   `mean_wind = met_data$wind_speed_10m`.
#' @param ... Additional hazard calibration parameters forwarded to
#'   [litter_hazard()] (e.g. `rain_threshold`, `material`, `gust_threshold`).
#'
#' @return The [litter_exposure()] data frame (`exposure`, `zone`,
#'   `directional_factor`, `leaves_site`, `sensitive_receptor`), one row per
#'   forecast hour.
#'
#' @seealso [litter_hazard()], [litter_exposure()], [odour_risk()].
#' @export
litter_risk <- function(met_data, site,
                        use_wetness_state    = FALSE,
                        direction_tol        = 15,
                        p_open_min           = 0.5,
                        move_threshold       = 20,
                        offsite_threshold    = 45,
                        default_permeability = 0.5,
                        reach_per_ms         = NULL,
                        ...) {
  checkmate::assert_data_frame(met_data, min.rows = 1)
  .assert_required_cols(
    met_data, "wind_direction_10m", arg = "met_data",
    info = "Required for the exposure layer: wind_direction_10m (degrees)."
  )

  # `...` carries hazard calibration params only (the exposure params are named
  # formals here), so it reaches litter_hazard() -> litter_hazard_vec() and not
  # litter_exposure().
  hazard <- litter_hazard(met_data, use_wetness_state = use_wetness_state, ...)

  exp_args <- list(
    hazard               = hazard,
    wind_direction_10m   = met_data$wind_direction_10m,
    site                 = site,
    direction_tol        = direction_tol,
    p_open_min           = p_open_min,
    move_threshold       = move_threshold,
    offsite_threshold    = offsite_threshold,
    default_permeability = default_permeability
  )
  # Refined distance-reach mode: forward the mean wind as the reach driver.
  if (!is.null(reach_per_ms)) {
    exp_args$mean_wind    <- met_data$wind_speed_10m
    exp_args$reach_per_ms <- reach_per_ms
  }

  do.call(litter_exposure, exp_args)
}
