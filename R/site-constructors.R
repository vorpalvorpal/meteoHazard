# ===========================================================================
# Convenience constructors that build mh_site objects from non-GIS input
# formats. These are the "lightweight" entry points: they accept simple
# tabular descriptions and construct the corresponding sf geometry + roles,
# then delegate to mh_site() for all validation.
#
# C4 Stage 1 (issue #18).
# ===========================================================================

# ---------------------------------------------------------------------------
# site_from_sectors()
# ---------------------------------------------------------------------------

#' Build an mh_site from compass-sector barrier descriptions
#'
#' Converts a compass-sector data frame (the format accepted by the legacy
#' `litter_exposure()` `site` argument) into an [`mh_site`] S7 object.
#' A wedge-shaped polygon is constructed for each sector row on a ring of
#' radius `ring_radius` centred on `centroid`. The returned site has:
#'
#' * A `(litter, source)` point feature at `centroid` (id `"source"`).
#' * One `(litter, barrier)` polygon per sector row (ids `"barrier_1"`,
#'   `"barrier_2"`, ...), carrying numeric columns `bearing_start` and
#'   `bearing_end` so that downstream functions can use bearing containment
#'   tests directly.
#'
#' @param sectors data.frame with columns:
#'   \describe{
#'     \item{`arc_start`, `arc_end`}{Compass labels (`N`, `NE`, `E`, `SE`,
#'       `S`, `SW`, `W`, `NW`) giving the clockwise start and end of each
#'       boundary sector.}
#'     \item{`permeability`}{Numeric `[0, 1]`.}
#'     \item{`sensitive`}{Logical.}
#'   }
#' @param centroid An `sf` object containing a single POINT feature in any
#'   CRS.  Must not be `NULL`.  If the CRS is geographic it will be
#'   reprojected to `epsg` before use.
#' @param ring_radius Numeric (metres). Radius of the ring on which barrier
#'   arc polygons are placed. Default `1000`.
#' @param epsg Integer. Projected metric EPSG code passed to [mh_site()].
#'
#' @return An [`mh_site`] object.
#' @export
site_from_sectors <- function(sectors, centroid, ring_radius = 1000, epsg) {

  # ---- Guard: centroid must not be NULL ------------------------------------
  if (is.null(centroid)) {
    cli::cli_abort(
      "{.arg centroid} must be an {.cls sf} POINT object, not {.val NULL}.",
      class = "meteoHazard_input_error"
    )
  }

  # ---- Validate sectors ----------------------------------------------------
  checkmate::assert_data_frame(sectors, min.rows = 1)
  required_cols <- c("arc_start", "arc_end", "permeability", "sensitive")
  missing_cols  <- setdiff(required_cols, names(sectors))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg sectors} is missing required columns: {.val {missing_cols}}.",
      class = "meteoHazard_input_error"
    )
  }
  valid_labels <- names(LITTER_COMPASS_DEGREES)
  bad_labels   <- setdiff(
    c(sectors$arc_start, sectors$arc_end), valid_labels
  )
  if (length(bad_labels) > 0) {
    cli::cli_abort(
      c("{.arg sectors} contains invalid compass label(s): {.val {bad_labels}}.",
        "i" = "Valid labels are: {.val {valid_labels}}."),
      class = "meteoHazard_input_error"
    )
  }
  checkmate::assert_numeric(sectors$permeability, lower = 0, upper = 1,
                            any.missing = FALSE,
                            .var.name = "sectors$permeability")
  checkmate::assert_logical(sectors$sensitive, any.missing = FALSE,
                            .var.name = "sectors$sensitive")

  # ---- Validate centroid is an sf object -----------------------------------
  if (!inherits(centroid, "sf")) {
    cli::cli_abort(
      "{.arg centroid} must be an {.cls sf} object.",
      class = "meteoHazard_input_error"
    )
  }
  checkmate::assert_number(ring_radius, lower = 0, finite = TRUE)
  checkmate::assert_integerish(epsg, len = 1)

  # ---- Reproject centroid to the target EPSG if needed --------------------
  target_crs <- sf::st_crs(as.integer(epsg))
  if (sf::st_crs(centroid) != target_crs) {
    centroid <- sf::st_transform(centroid, as.integer(epsg))
  }
  ctr_coords <- sf::st_coordinates(centroid)
  cx <- ctr_coords[1, "X"]
  cy <- ctr_coords[1, "Y"]

  # ---- Build geometries ----------------------------------------------------
  n_sectors <- nrow(sectors)

  # Source point
  source_geom <- sf::st_sfc(
    sf::st_point(c(cx, cy)),
    crs = as.integer(epsg)
  )

  # Barrier arc polygons (wedge shape)
  barrier_geoms      <- vector("list", n_sectors)
  bearing_start_vals <- numeric(n_sectors)
  bearing_end_vals   <- numeric(n_sectors)

  N_arc <- 30L  # number of points along each outer arc

  for (k in seq_len(n_sectors)) {
    alpha <- unname(LITTER_COMPASS_DEGREES[sectors$arc_start[k]])
    beta  <- unname(LITTER_COMPASS_DEGREES[sectors$arc_end[k]])

    bearing_start_vals[k] <- alpha
    bearing_end_vals[k]   <- beta

    # Generate arc angles (clockwise from north).
    # When alpha == beta (full 360 circle) or alpha < beta: simple range.
    # When alpha > beta: wraps through north.
    if (alpha <= beta) {
      thetas <- seq(alpha, beta, length.out = N_arc)
    } else {
      # wraps through north: e.g. 315 -> 0 -> 45
      half <- N_arc %/% 2L
      thetas <- c(
        seq(alpha, 360, length.out = half),
        seq(0,     beta, length.out = N_arc - half)
      )
    }

    # Convert bearing (deg from N, clockwise) to (x, y) offsets
    # bearing θ → x = r·sin(θ), y = r·cos(θ)
    theta_rad <- thetas * pi / 180
    arc_x <- cx + ring_radius * sin(theta_rad)
    arc_y <- cy + ring_radius * cos(theta_rad)

    # Polygon: outer arc + centroid + close (repeat first point)
    poly_x <- c(arc_x, cx, arc_x[1])
    poly_y <- c(arc_y, cy, arc_y[1])

    barrier_geoms[[k]] <- sf::st_polygon(list(cbind(poly_x, poly_y)))
  }

  barrier_sfc <- sf::st_sfc(barrier_geoms, crs = as.integer(epsg))

  # ---- Assemble sf feature collection -------------------------------------
  # Source row
  source_sf <- sf::st_sf(
    id            = "source",
    bearing_start = NA_real_,
    bearing_end   = NA_real_,
    geometry      = source_geom,
    stringsAsFactors = FALSE
  )

  # Barrier rows
  barrier_sf <- sf::st_sf(
    id            = paste0("barrier_", seq_len(n_sectors)),
    bearing_start = bearing_start_vals,
    bearing_end   = bearing_end_vals,
    geometry      = barrier_sfc,
    stringsAsFactors = FALSE
  )

  features <- rbind(source_sf, barrier_sf)

  # ---- Roles data.frame ---------------------------------------------------
  roles <- data.frame(
    feature_id   = c("source", paste0("barrier_", seq_len(n_sectors))),
    hazard       = "litter",
    role         = c("source", rep("barrier", n_sectors)),
    permeability = c(NA_real_, sectors$permeability),
    sensitive    = c(NA, sectors$sensitive),
    stringsAsFactors = FALSE
  )

  # ---- Construct and return mh_site ----------------------------------------
  mh_site(
    features = features,
    roles    = roles,
    epsg     = as.integer(epsg)
  )
}
