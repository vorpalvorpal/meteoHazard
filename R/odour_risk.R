#' Landfill odour risk index (combined hazard + exposure)
#'
#' Convenience wrapper that computes the atmospheric ventilation state and maps
#' it through the geometry-aware [odour_exposure()] layer in one call,
#' returning the worst-case 0-100 odour exposure across the receptors for each
#' forecast hour.
#'
#' This is the entry point for the common case; use [odour_exposure()] directly
#' when you need to inspect the per-source or per-receptor concentrations.
#'
#' @inheritParams odour_exposure
#' @param datetime Optional `POSIXct` vector for the consecutive-hourly spacing
#'   check; see [odour_hazard()].
#'
#' @return A numeric vector of length `nrow(met_data)`: the worst-case 0-100
#'   odour exposure across receptors for each hour.
#'
#' @seealso [odour_exposure()], [odour_hazard()], [ventilation_state()].
#' @export
odour_risk <- function(met_data, site,
                       stability = c("turner", "shear"),
                       map_c50 = 0.3,
                       terrain_backend = c("none", "descriptors"),
                       shelter      = FALSE,
                       shelter_h_mix = FALSE,
                       impaction    = FALSE,
                       datetime = NULL) {
  stability <- match.arg(stability)
  .assert_hourly(datetime)
  odour_exposure(met_data, site,
                 stability       = stability,
                 map_c50         = map_c50,
                 terrain_backend = terrain_backend,
                 shelter         = shelter,
                 shelter_h_mix   = shelter_h_mix,
                 impaction       = impaction)
}
