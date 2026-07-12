#' Landfill odour risk (combined hazard + exposure)
#'
#' Convenience wrapper that computes the atmospheric ventilation state and maps
#' it through the geometry-aware [odour_exposure()] layer in one call,
#' returning the per-receptor **relative concentration** matrix for each
#' forecast hour (the physical layer; see [odour_exposure()] for the return
#' shape and the rationale for not baking in a 0-100 scale).
#'
#' This is the entry point for the common case; use [odour_exposure()] directly
#' for the same result. Mapping the relative concentration onto an operational
#' index is a site-specific calibration step (see issues #11/#8).
#'
#' @inheritParams odour_exposure
#' @param datetime Optional `POSIXct` vector for the consecutive-hourly spacing
#'   check; see [odour_hazard()].
#' @param pool_cap Logical (default `TRUE`). Passed to [odour_exposure()];
#'   caps `h_mix` at the nocturnal cold-pool depth on stable nights.
#' @param odorant_solubility Number in `[0, 1]` (default
#'   `ODOUR_CONSTANTS$ODORANT_SOLUBILITY_DEFAULT`, 0.5). Passed to
#'   [odour_exposure()]; blends `W_rain` between no washout (0) and the
#'   soluble-limit tiers (1).
#' @param rim_venting Logical. Passed to [odour_exposure()]; activates C8
#'   upslope rim-venting when `TRUE`. Default `FALSE`.
#'
#' @return A numeric matrix (`nrow(met_data)` x n_receptors) of relative odour
#'   concentration; see [odour_exposure()].
#'
#' @seealso [odour_exposure()], [odour_hazard()], [ventilation_state()].
#' @export
odour_risk <- function(met_data, site,
                       stability = c("turner", "shear"),
                       terrain_backend = c("none", "descriptors"),
                       shelter      = FALSE,
                       shelter_h_mix = FALSE,
                       pool_cap = TRUE,
                       odorant_solubility = ODOUR_CONSTANTS$ODORANT_SOLUBILITY_DEFAULT,
                       datetime = NULL,
                       rim_venting = FALSE) {
  stability <- match.arg(stability)
  .assert_hourly(datetime)
  odour_exposure(met_data, site,
                 stability       = stability,
                 terrain_backend = terrain_backend,
                 shelter         = shelter,
                 shelter_h_mix   = shelter_h_mix,
                 pool_cap        = pool_cap,
                 odorant_solubility = odorant_solubility,
                 rim_venting     = rim_venting)
}
