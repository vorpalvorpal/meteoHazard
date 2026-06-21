# ===========================================================================
# S7 site-model classes: mh_terrain() and mh_site().
#
# mh_terrain holds pinned topographic descriptors (all optional; NA = flat).
# mh_site bundles an sf feature set, a roles data.frame, optional terrain,
# a projected CRS (epsg), and a vertical datum string.
#
# Units contract (see R/units-helpers.R): dimensional inputs accept a bare
# numeric (assumed in the documented unit) OR a `units` object (converted;
# mismatch errors). Stored internally as plain doubles in canonical units.
# ===========================================================================

.VALID_HAZARDS <- c("odour", "dust", "litter")
.VALID_ROLES   <- c("source", "receptor", "barrier", "sink")
.VALID_DATUMS  <- "AGL"

# ---------------------------------------------------------------------------
# mh_terrain
# ---------------------------------------------------------------------------

#' Pinned topographic terrain descriptors
#'
#' Constructs a validated `mh_terrain` object holding the scalar terrain
#' descriptors used by the meteoHazard dispersion and hazard models. All
#' descriptors are optional; pass `NA` (the default) for a flat site with no
#' terrain correction.
#'
#' Dimensional arguments (`relief`, `valley_depth`, `basin_capacity`) accept a
#' bare numeric (assumed already in the documented unit) or a [`units`][units::units]
#' object (converted; dimensional mismatch errors). They are stored internally
#' as plain doubles in canonical units.
#'
#' @param relief Non-negative numeric or units object. Height above local base
#'   in metres (NOT a standardised deviation index).
#' @param valley_depth Non-negative numeric or units object. Depth of the
#'   containing valley in metres.
#' @param basin_capacity Non-negative numeric or units object. Capacity of the
#'   containing basin in m^3.
#' @param drainage_bearing Numeric. Dominant drainage-flow bearing in degrees,
#'   in `[0, 360)`.
#' @param flow_convergence Numeric. Dimensionless flow-convergence index; must
#'   be finite (not Inf/NaN).
#' @param slope Numeric. Mean slope angle in degrees.
#' @param aspect Numeric. Aspect angle in degrees.
#' @param taf Numeric. Topographic amplification factor (dimensionless, >= 1).
#' @param shelter_index Numeric. Openness-based shelter index in degrees; must
#'   be finite (not Inf/NaN).
#' @param meta Named list. Scale fields: `relief_radius`, `channel_threshold`,
#'   `fetch_L`, `dem_resolution`, `datum`.
#'
#' @return An `mh_terrain` S7 object.
#' @export
mh_terrain <- S7::new_class(
  "mh_terrain",
  properties = list(
    relief          = S7::class_double,
    valley_depth    = S7::class_double,
    basin_capacity  = S7::class_double,
    drainage_bearing = S7::class_double,
    flow_convergence = S7::class_double,
    slope           = S7::class_double,
    aspect          = S7::class_double,
    taf             = S7::class_double,
    shelter_index   = S7::class_double,
    meta            = S7::class_list
  ),
  constructor = function(relief = NA_real_,
                         valley_depth = NA_real_,
                         basin_capacity = NA_real_,
                         drainage_bearing = NA_real_,
                         flow_convergence = NA_real_,
                         slope = NA_real_,
                         aspect = NA_real_,
                         taf = NA_real_,
                         shelter_index = NA_real_,
                         meta = list()) {
    relief           <- .drop_to(relief,           "m",   "relief")
    valley_depth     <- .drop_to(valley_depth,     "m",   "valley_depth")
    basin_capacity   <- .drop_to(basin_capacity,   "m^3", "basin_capacity")

    # Store plain numeric NA when NULL (bare NA passes through .drop_to as-is
    # because .as_units returns NULL for NULL, but bare NA_real_ is numeric)
    if (is.null(relief))           relief           <- NA_real_
    if (is.null(valley_depth))     valley_depth     <- NA_real_
    if (is.null(basin_capacity))   basin_capacity   <- NA_real_

    S7::new_object(
      S7::S7_object(),
      relief           = as.double(relief),
      valley_depth     = as.double(valley_depth),
      basin_capacity   = as.double(basin_capacity),
      drainage_bearing = as.double(drainage_bearing),
      flow_convergence = as.double(flow_convergence),
      slope            = as.double(slope),
      aspect           = as.double(aspect),
      taf              = as.double(taf),
      shelter_index    = as.double(shelter_index),
      meta             = as.list(meta)
    )
  },
  validator = function(self) {
    errs <- character(0)

    check_nonneg <- function(val, nm) {
      if (!is.na(val) && val < 0)
        errs <<- c(errs, paste0(nm, " must be non-negative (got ", val, ")"))
    }
    check_nonneg(self@relief,          "relief")
    check_nonneg(self@valley_depth,    "valley_depth")
    check_nonneg(self@basin_capacity,  "basin_capacity")

    if (!is.na(self@taf) && self@taf < 1)
      errs <- c(errs, paste0("taf must be >= 1 (got ", self@taf, ")"))

    if (!is.na(self@drainage_bearing) &&
        (self@drainage_bearing < 0 || self@drainage_bearing >= 360))
      errs <- c(errs,
                paste0("drainage_bearing must be in [0, 360) (got ",
                       self@drainage_bearing, ")"))

    if ((!is.na(self@flow_convergence) || is.nan(self@flow_convergence)) &&
        !is.finite(self@flow_convergence))
      errs <- c(errs, "flow_convergence must be finite")

    if ((!is.na(self@shelter_index) || is.nan(self@shelter_index)) &&
        !is.finite(self@shelter_index))
      errs <- c(errs, "shelter_index must be finite")

    if (length(errs) > 0) {
      cli::cli_abort(
        c("Invalid {.cls mh_terrain}:", setNames(errs, rep("x", length(errs)))),
        class = "meteoHazard_input_error"
      )
    }
  }
)

# ---------------------------------------------------------------------------
# mh_site
# ---------------------------------------------------------------------------

#' Geo-referenced site model
#'
#' Bundles an `sf` feature collection, a roles table, optional terrain
#' descriptors, a projected CRS, and a vertical datum into a single validated
#' domain object.
#'
#' @param features An `sf` object with at least an `id` column and a geometry
#'   column. May optionally carry `elevation` and `emit_height` columns.
#' @param roles A data.frame with columns `feature_id`, `hazard`, and `role`.
#'   Optional columns: `weight`, `sensitive`, `permeability`. Every
#'   `feature_id` must exist in `features$id`.
#' @param terrain An `mh_terrain` object or `NULL` (default).
#' @param epsg Integer. EPSG code for the projected metric CRS to use
#'   internally. Must identify a projected (non-geographic) CRS. Features in any
#'   other CRS (geographic or a different projection) are reprojected to `epsg`
#'   automatically.
#' @param datum Character. Vertical datum string. Currently only `"AGL"` is
#'   recognised.
#'
#' @return An `mh_site` S7 object.
#' @export
mh_site <- S7::new_class(
  "mh_site",
  properties = list(
    features = S7::class_any,
    roles    = S7::class_any,
    terrain  = S7::class_any,
    epsg     = S7::class_integer,
    datum    = S7::class_character
  ),
  constructor = function(features,
                         roles,
                         terrain = NULL,
                         epsg,
                         datum = "AGL") {
    # --- datum check -------------------------------------------------------
    if (!datum %in% .VALID_DATUMS) {
      valid_datums <- .VALID_DATUMS
      cli::cli_abort(
        c("{.arg datum} must be one of {.val {valid_datums}}; got {.val {datum}}."),
        class = "meteoHazard_input_error"
      )
    }

    # --- features must be sf -----------------------------------------------
    if (!inherits(features, "sf")) {
      cli::cli_abort(
        "{.arg features} must be an {.cls sf} object.",
        class = "meteoHazard_input_error"
      )
    }

    # --- features must have id column --------------------------------------
    if (!"id" %in% names(features)) {
      cli::cli_abort(
        "{.arg features} must have an {.val id} column.",
        class = "meteoHazard_input_error"
      )
    }

    # --- roles required columns --------------------------------------------
    .assert_required_cols(
      roles,
      c("feature_id", "hazard", "role"),
      arg = "roles"
    )

    # --- every roles$feature_id must exist in features$id ------------------
    bad_ids <- setdiff(roles$feature_id, features$id)
    if (length(bad_ids) > 0) {
      cli::cli_abort(
        c("{.arg roles} contains {.val feature_id} values not found in {.arg features$id}:",
          "x" = "{.val {bad_ids}}"),
        class = "meteoHazard_input_error"
      )
    }

    # --- hazard / role values ----------------------------------------------
    bad_hazards <- setdiff(roles$hazard, .VALID_HAZARDS)
    if (length(bad_hazards) > 0) {
      valid_hazards <- .VALID_HAZARDS
      cli::cli_abort(
        c("Unknown {.arg roles$hazard} value(s): {.val {bad_hazards}}.",
          "i" = "Recognised hazards: {.val {valid_hazards}}."),
        class = "meteoHazard_input_error"
      )
    }

    bad_roles <- setdiff(roles$role, .VALID_ROLES)
    if (length(bad_roles) > 0) {
      valid_roles <- .VALID_ROLES
      cli::cli_abort(
        c("Unknown {.arg roles$role} value(s): {.val {bad_roles}}.",
          "i" = "Recognised roles: {.val {valid_roles}}."),
        class = "meteoHazard_input_error"
      )
    }

    # --- terrain must be mh_terrain or NULL --------------------------------
    if (!is.null(terrain) && !S7::S7_inherits(terrain, mh_terrain)) {
      cli::cli_abort(
        "{.arg terrain} must be an {.cls mh_terrain} object or {.val NULL}.",
        class = "meteoHazard_input_error"
      )
    }

    # --- epsg: validate the projection is metric/projected -----------------
    epsg_int <- as.integer(epsg)

    target_crs <- tryCatch(
      sf::st_crs(epsg_int),
      error = function(e) {
        cli::cli_abort(
          "{.arg epsg} {.val {epsg_int}} is not a recognised EPSG code.",
          class = "meteoHazard_input_error"
        )
      }
    )

    if (isTRUE(sf::st_is_longlat(target_crs))) {
      cli::cli_abort(
        c("{.arg epsg} {.val {epsg_int}} is a geographic (lon/lat) CRS.",
          "i" = "Supply a projected metric EPSG (e.g. a UTM zone)."),
        class = "meteoHazard_input_error"
      )
    }

    # --- CRS handling for features -----------------------------------------
    feat_crs <- sf::st_crs(features)

    if (is.na(feat_crs)) {
      cli::cli_abort(
        "{.arg features} has no CRS set. Assign one before constructing {.fn mh_site}.",
        class = "meteoHazard_input_error"
      )
    }

    if (sf::st_crs(features) != target_crs) {
      # reproject any differing CRS (geographic or projected) to the metric epsg
      features <- sf::st_transform(features, epsg_int)
    }

    S7::new_object(
      S7::S7_object(),
      features = features,
      roles    = as.data.frame(roles),
      terrain  = terrain,
      epsg     = epsg_int,
      datum    = as.character(datum)
    )
  }
)

# ---------------------------------------------------------------------------
# print / format method for mh_site
# ---------------------------------------------------------------------------

S7::method(print, mh_site) <- function(x, ...) {
  n_feat  <- nrow(x@features)
  n_roles <- nrow(x@roles)
  haz_tbl <- if (n_roles > 0) table(x@roles$hazard) else integer(0)
  crs_nm  <- sf::st_crs(x@features)$Name
  if (is.na(crs_nm) || is.null(crs_nm)) crs_nm <- paste0("EPSG:", x@epsg)

  cli::cli_h2("mh_site")
  cli::cli_alert_info("{n_feat} feature{?s}, {n_roles} role{?s}")
  if (length(haz_tbl) > 0) {
    role_str <- paste(names(haz_tbl), haz_tbl, sep = ": ", collapse = ", ")
    cli::cli_alert_info("Roles by hazard: {role_str}")
  }
  cli::cli_alert_info("CRS: {crs_nm}  |  datum: {x@datum}")
  if (!is.null(x@terrain))
    cli::cli_alert_info("Terrain: attached")
  else
    cli::cli_alert_info("Terrain: none")

  invisible(x)
}
