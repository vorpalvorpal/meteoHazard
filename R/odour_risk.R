#' Landfill odour risk index (combined hazard + exposure)
#'
#' Convenience wrapper that computes the meteorological odour [odour_hazard()]
#' and maps it through the geometry-aware [odour_exposure()] layer in one call,
#' returning the worst-case 0-100 odour exposure across the receptors for each
#' forecast hour.
#'
#' This is the entry point for the common case; use [odour_hazard()] alone for a
#' receptor-independent site index, or [odour_exposure()] directly to reuse a
#' precomputed hazard. See `specs/Odour_v2.md` for the model.
#'
#' @inheritParams odour_exposure
#' @param met_data A data frame, one row per consecutive hourly timestep, with
#'   the columns required by both layers (see [odour_hazard()] and
#'   [odour_exposure()]). `NA` values are permitted and handled conservatively.
#' @param stability Stability estimator (`"turner"` or `"shear"`); see
#'   [odour_hazard()].
#' @param datetime Optional `POSIXct` vector for the consecutive-hourly spacing
#'   check; see [odour_hazard()].
#'
#' @return A numeric vector of length `nrow(met_data)`: the worst-case 0-100
#'   odour exposure across receptors for each hour.
#'
#' @seealso [odour_hazard()], [odour_exposure()].
#' @export
generate_odour_risk_index <- function(met_data, receptors, drainage_axes = NULL,
                                       stability = c("turner", "shear"),
                                       datetime = NULL) {
  stability <- match.arg(stability)
  hazard <- odour_hazard(met_data, stability = stability, datetime = datetime)
  odour_exposure(hazard, met_data, receptors, drainage_axes = drainage_axes,
                 stability = stability)
}
