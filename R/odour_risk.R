#' Landfill odour risk (combined hazard + exposure)
#'
#' Convenience wrapper that computes the atmospheric ventilation state and maps
#' it through the geometry-aware [odour_exposure()] layer in one call,
#' returning the per-receptor **relative concentration** matrix for each
#' forecast hour (the physical layer; see [odour_exposure()] for the return
#' shape and the rationale for not baking in a 0-100 scale).
#'
#' This is the entry point for the common case; use [odour_exposure()] directly
#' for the same result, and [odour_index_interim()] for the parked, uncalibrated
#' 0-100 screening index.
#'
#' @inheritParams odour_exposure
#' @param datetime Optional `POSIXct` vector for the consecutive-hourly spacing
#'   check; see [odour_hazard()].
#' @param rim_venting Logical. Passed to [odour_exposure()]; activates C8
#'   upslope rim-venting when `TRUE`. Default `FALSE`.
#'
#' @return A numeric matrix (`nrow(met_data)` x n_receptors) of relative odour
#'   concentration; see [odour_exposure()].
#'
#' @seealso [odour_exposure()], [odour_index_interim()], [odour_hazard()],
#'   [ventilation_state()].
#' @export
odour_risk <- function(met_data, site,
                       stability = c("turner", "shear"),
                       terrain_backend = c("none", "descriptors"),
                       shelter      = FALSE,
                       shelter_h_mix = FALSE,
                       datetime = NULL,
                       rim_venting = FALSE) {
  stability <- match.arg(stability)
  .assert_hourly(datetime)
  odour_exposure(met_data, site,
                 stability       = stability,
                 terrain_backend = terrain_backend,
                 shelter         = shelter,
                 shelter_h_mix   = shelter_h_mix,
                 rim_venting     = rim_venting)
}
