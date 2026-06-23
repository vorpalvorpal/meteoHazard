# ===========================================================================
# DEM infrastructure helpers for mh_terrain_from_dem() (C5, issue #19).
# All functions are internal (dot-prefixed, not exported).
# Behind Suggests: whitebox, terra, elevatr, withr.
# ===========================================================================


# ---------------------------------------------------------------------------
# .require_dem_stack(need_download = FALSE)
# ---------------------------------------------------------------------------
# Checks that the GIS Suggests stack is available. Throws a classed
# meteoHazard_dependency_error when any required piece is missing, with a
# message that names the missing piece and suggests the install command.
# Set need_download = TRUE when an elevatr download is needed.
.require_dem_stack <- function(need_download = FALSE) {
  missing_pkgs <- character(0)
  missing_binary <- FALSE

  if (!requireNamespace("whitebox", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "whitebox")
  } else {
    if (!isTRUE(whitebox::check_whitebox_binary())) {
      missing_binary <- TRUE
    }
  }

  if (!requireNamespace("terra", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "terra")
  }

  if (need_download && !requireNamespace("elevatr", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "elevatr")
  }

  if (length(missing_pkgs) > 0L) {
    cli::cli_abort(
      c("mh_terrain_from_dem() requires packages that are not installed:",
        "i" = "Missing: {.pkg {missing_pkgs}}",
        "i" = "Install with: {.code install.packages(c({paste0('\"', missing_pkgs, '\"', collapse=', ')}))}"),
      class = "meteoHazard_dependency_error"
    )
  }

  if (missing_binary) {
    cli::cli_abort(
      c("mh_terrain_from_dem() requires the {.pkg whitebox} binary (WhiteboxTools).",
        "i" = "The binary is not installed. Run: {.code whitebox::install_whitebox()}"),
      class = "meteoHazard_dependency_error"
    )
  }

  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# .dem_ingest(dem, epsg)
# ---------------------------------------------------------------------------
# Accepts a terra::SpatRaster, a file path string, or a bounding-box object
# (sf bbox or similar), reprojects to the metric epsg CRS if needed, and
# returns a terra::SpatRaster in that CRS.
# Throws a meteoHazard_input_error when:
#   - dem is a string but the file does not exist
#   - dem is a SpatRaster with no CRS
#   - the epsg resolves to a geographic (lon/lat) CRS
.dem_ingest <- function(dem, epsg) {
  # Validate epsg first (must be a projected metric CRS).
  target_crs <- terra::crs(paste0("EPSG:", as.integer(epsg)))
  if (nchar(target_crs) == 0L) {
    cli::cli_abort(
      "epsg {.val {epsg}} is not a recognised CRS.",
      class = "meteoHazard_input_error"
    )
  }
  test_rast <- terra::rast(nrows = 1L, ncols = 1L,
                             xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                             crs = paste0("EPSG:", as.integer(epsg)))
  if (isTRUE(terra::is.lonlat(test_rast))) {
    cli::cli_abort(
      "epsg {.val {epsg}} is a geographic (lon/lat) CRS; a projected metric CRS is required.",
      class = "meteoHazard_input_error"
    )
  }

  # Load the DEM.
  if (is.character(dem)) {
    if (!file.exists(dem)) {
      cli::cli_abort(
        "DEM file not found: {.path {dem}}",
        class = "meteoHazard_input_error"
      )
    }
    dem <- terra::rast(dem)
  }

  if (!inherits(dem, "SpatRaster")) {
    cli::cli_abort(
      "{.arg dem} must be a {.cls SpatRaster}, a file path, or a bounding box.",
      class = "meteoHazard_input_error"
    )
  }

  # Check the DEM has a CRS.
  if (terra::crs(dem) == "") {
    cli::cli_abort(
      "The DEM has no CRS. Supply a CRS-tagged {.cls SpatRaster} or a file with embedded CRS.",
      class = "meteoHazard_input_error"
    )
  }

  # Reproject to target if needed.
  dem_crs_epsg <- tryCatch(
    as.integer(gsub(".*EPSG:(\\d+).*", "\\1", terra::crs(dem, describe = TRUE)$code)),
    error = function(e) NA_integer_
  )
  if (!is.na(dem_crs_epsg) && dem_crs_epsg == as.integer(epsg)) {
    return(dem)
  }
  terra::project(dem, paste0("EPSG:", as.integer(epsg)), method = "bilinear")
}


# ---------------------------------------------------------------------------
# .dem_scan_range(dem)
# ---------------------------------------------------------------------------
# Returns a named list(min_scale, max_scale) in metres, bounding the multi-
# scale DEV scan. min_scale is one cell width; max_scale is half the DEM's
# shorter dimension. Both are in the DEM's native (metric) units.
.dem_scan_range <- function(dem) {
  res_xy    <- terra::res(dem)      # (xres, yres) in CRS units (metres)
  cell_size <- max(res_xy)          # conservative: larger cell as floor
  ext       <- terra::ext(dem)
  width_x   <- ext$xmax - ext$xmin
  width_y   <- ext$ymax - ext$ymin
  max_scale <- min(width_x, width_y) / 2

  list(min_scale = cell_size, max_scale = max_scale)
}

