#' Predict odour nuisance (stub)
#'
#' Will predict odour-nuisance potential from a waste-management site given
#' meteorological conditions, returning a management-relevant warning level.
#' The intended approach couples emission characteristics with atmospheric
#' dispersion (wind speed and direction, atmospheric stability) to estimate
#' downwind odour concentration at sensitive receptors.
#'
#' @details
#' This function is a placeholder and is **not yet implemented**. The planned
#' signature and return value may change.
#'
#' @param datetime POSIXct datetime vector (required).
#' @param latitude Latitude in decimal degrees (required).
#' @param longitude Longitude in decimal degrees (required).
#' @param ... Reserved for future meteorological and source parameters.
#'
#' @return Not yet defined. Intended to return a numeric odour-nuisance index
#'   or an ordered warning category.
#'
#' @keywords internal
#' @export
predict_odour <- function(datetime, latitude, longitude, ...) {
  cli::cli_abort(
    "{.fn predict_odour} is not yet implemented (planned for a future release)."
  )
}
