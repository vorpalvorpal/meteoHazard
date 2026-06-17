#' Predict wind-blown litter risk (stub)
#'
#' Will predict the risk of wind-blown litter escaping a waste-management site
#' given meteorological conditions, returning a management-relevant warning
#' level. The intended approach relates wind speed (and gust) to the threshold
#' at which lightweight material becomes airborne and is transported beyond the
#' site boundary.
#'
#' @details
#' This function is a placeholder and is **not yet implemented**. The planned
#' signature and return value may change.
#'
#' @param datetime POSIXct datetime vector (required).
#' @param latitude Latitude in decimal degrees (required).
#' @param longitude Longitude in decimal degrees (required).
#' @param ... Reserved for future meteorological and site parameters.
#'
#' @return Not yet defined. Intended to return a numeric litter-risk index or
#'   an ordered warning category.
#'
#' @keywords internal
#' @export
predict_litter <- function(datetime, latitude, longitude, ...) {
  cli::cli_abort(
    "{.fn predict_litter} is not yet implemented (planned for a future release)."
  )
}
