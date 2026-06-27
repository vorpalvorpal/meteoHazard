# helper-synthetic-dem.R
# ===========================================================================
# Builders for analytically-known synthetic DEMs used in test-terrain-from-dem.R.
# Each returns a terra::SpatRaster in EPSG:32755 (or a specified metric CRS),
# centred at local coordinates (0, 0) so that a source point placed at
# sf::st_point(c(0, 0)) sits at the DEM origin.
#
# All helpers call skip_if_not_installed("terra") so the entire helper file
# degrades gracefully in environments without the GIS stack.
#
# .wbt_ok() / .wbt_skip() are defined here (not in the test file) so every
# describe() in test-terrain-from-dem.R picks them up automatically.
# ===========================================================================

# Guard helpers -----------------------------------------------------------

#' TRUE when whitebox (package + binary) and terra are both present.
.wbt_ok <- function() {
  requireNamespace("whitebox", quietly = TRUE) &&
    requireNamespace("terra",    quietly = TRUE) &&
    isTRUE(whitebox::check_whitebox_binary())
}

#' Skip the current test when the full GIS stack is absent.
.wbt_skip <- function() {
  skip_if_not_installed("whitebox")
  skip_if_not_installed("terra")
  if (!isTRUE(whitebox::check_whitebox_binary()))
    skip("WhiteboxTools binary absent; run whitebox::install_whitebox()")
}

# Core raster builder -----------------------------------------------------

# Populate a raster by evaluating fun(x, y) at each cell centre, where x and
# y are the easting and northing in the raster's CRS.
.make_dem <- function(fun, ncell = 64, extent_m = 2000, crs_epsg = 32755) {
  skip_if_not_installed("terra")
  r <- terra::rast(
    nrows = ncell, ncols = ncell,
    xmin = -extent_m / 2, xmax = extent_m / 2,
    ymin = -extent_m / 2, ymax = extent_m / 2,
    crs  = paste0("EPSG:", crs_epsg)
  )
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- fun(xy[, "x"], xy[, "y"])
  r
}

# Analytic DEM shapes -----------------------------------------------------

#' Gaussian hill: z = H * exp(-r^2 / (2 * sigma^2)).
#'
#' Source at peak (0, 0): `relief ≈ H`; `relief_radius` on the order of
#' `sigma`. All elevations are ≥ 0.
.dem_gaussian_hill <- function(H = 100, sigma = 200, ncell = 64,
                                extent_m = 2000) {
  .make_dem(
    function(x, y) H * exp(-(x^2 + y^2) / (2 * sigma^2)),
    ncell = ncell, extent_m = extent_m
  )
}

#' Gaussian hill embedded in a broad basin.
#'
#' A hill of height `H_hill` sits inside a basin of depth `H_basin` below a
#' flat plateau. The DEVmax-selected radius should be shorter than for the
#' isolated hill, yielding a smaller `relief`.
.dem_hill_in_basin <- function(H_hill = 80, H_basin = 200, sigma = 150,
                                R_basin = 600, ncell = 64, extent_m = 2000) {
  .make_dem(
    function(x, y) {
      r     <- sqrt(x^2 + y^2)
      hill  <- H_hill * exp(-(x^2 + y^2) / (2 * sigma^2))
      basin <- H_basin * pmin(r / R_basin, 1)  # conical basin around origin
      hill + basin
    },
    ncell = ncell, extent_m = extent_m
  )
}

#' V-valley running E-W: floor at y = 0, rims at y = ±extent_m/2 at height D.
#'
#' Source at valley floor (0, 0): `valley_depth ≈ D`.
.dem_v_valley <- function(D = 80, ncell = 64, extent_m = 2000) {
  .make_dem(
    function(x, y) D * abs(y) / (extent_m / 2),
    ncell = ncell, extent_m = extent_m
  )
}

#' Flat plain: z = 0 everywhere.
#'
#' `relief ≈ 0`, `valley_depth ≈ 0`, `taf ≈ 1`, `drainage_bearing = NA`.
.dem_flat <- function(ncell = 32, extent_m = 1000) {
  .make_dem(function(x, y) rep(0, length(x)), ncell = ncell, extent_m = extent_m)
}

#' Conical basin: z = D * pmin(r/R, 1) where r = sqrt(x^2 + y^2).
#'
#' Source at floor (0, 0). Sill at height D, radius R.
#' Analytic volume to sill: `pi * R^2 * D / 3` (cone).
.dem_basin <- function(D = 50, R = 400, ncell = 64, extent_m = 1200) {
  .make_dem(
    function(x, y) {
      r <- sqrt(x^2 + y^2)
      D * pmin(r / R, 1)
    },
    ncell = ncell, extent_m = extent_m
  )
}

#' Eastward-tilted plane.
#'
#' z = slope * (extent_m/2 - x): higher to west, lower to east, all ≥ 0.
#' Source at (0, 0) is at mid-height. Dominant downslope direction = east ≈ 90°.
.dem_tilted_east <- function(slope_m_per_m = 0.05, ncell = 64, extent_m = 2000) {
  .make_dem(
    function(x, y) slope_m_per_m * (extent_m / 2 - x),
    ncell = ncell, extent_m = extent_m
  )
}

#' Converging hollow (concave): z = -H * exp(-r^2/(2*sigma^2)) + H (inverted Gaussian).
#'
#' Source at the lowest point (0, 0). Flow converges toward origin.
#' `flow_convergence` here should exceed that of the diverging spur below.
.dem_converging_hollow <- function(H = 50, sigma = 300, ncell = 64,
                                    extent_m = 2000) {
  .make_dem(
    function(x, y) H * (1 - exp(-(x^2 + y^2) / (2 * sigma^2))),
    ncell = ncell, extent_m = extent_m
  )
}

#' Diverging spur (convex): z = H * exp(-r^2/(2*sigma^2)).
#'
#' Source at the highest point (0, 0). Flow diverges away from origin.
#' `flow_convergence` here should be lower than for the converging hollow.
.dem_diverging_spur <- function(H = 50, sigma = 300, ncell = 64,
                                 extent_m = 2000) {
  .make_dem(
    function(x, y) H * exp(-(x^2 + y^2) / (2 * sigma^2)),
    ncell = ncell, extent_m = extent_m
  )
}

# Convenience point constructors ------------------------------------------

#' Source point at the DEM centre (0, 0) in EPSG:32755.
.src_centre <- function() {
  sf::st_sf(id = "src",
            geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 32755))
}

#' Receptor point at offset (dx_m, dy_m) from origin in EPSG:32755.
.rec_offset <- function(dx_m, dy_m, id = "rec") {
  sf::st_sf(id = id,
            geometry = sf::st_sfc(sf::st_point(c(dx_m, dy_m)), crs = 32755))
}
