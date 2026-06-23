# ===========================================================================
# DEM terrain helper (C5, issue #19): mh_terrain_from_dem().
# Setup-time only — NEVER on the run path.
# Requires: whitebox + binary, terra. Behind Suggests.
# ===========================================================================


#' Derive terrain descriptors from a DEM
#'
#' A setup-time helper that derives the C1 `mh_terrain` descriptors from a
#' digital elevation model using data-driven analysis scales. The result is a
#' frozen descriptor table for the run path. **Not for use on the run path** —
#' deriving descriptors is expensive (minutes); the run path reads the saved
#' `mh_terrain` object.
#'
#' @param dem A `terra::SpatRaster`, an absolute file path to a raster file, or
#'   an `sf` bounding box (triggers an `elevatr` download at setup time).
#' @param source An `sf` object containing the source location(s). Used to
#'   determine the descriptor reference point.
#' @param receptors An `sf` object of receptor points, or `NULL` (default).
#'   When supplied, per-receptor `rel_elevation` is derived and returned as
#'   `$receptor_fields`.
#' @param epsg Integer. EPSG code for the projected metric CRS used for all
#'   geometry. Must be a projected (non-geographic) CRS.
#' @param conditioning Breach/fill conditioning applied to the DEM before flow
#'   routing: `"breach"` (default — tries `wbt_breach_depressions`, then
#'   `wbt_fill_depressions_wang_and_liu` as fallback, then the raw DEM),
#'   `"fill"` (runs `wbt_fill_depressions_wang_and_liu` directly), or
#'   `"none"`. The method actually applied is recorded in `@meta$conditioning`.
#' @param shelter_fetch_m Numeric. Lookup radius in metres for the topographic
#'   openness (shelter_index) calculation. `NULL` (default) derives a
#'   physically-motivated default from the DEM scan range.
#' @param ... Additional arguments reserved for future use (caller descriptor
#'   overrides are passed here; see Details).
#'
#' @details
#' **Caller overrides:** any descriptor name (`relief`, `valley_depth`,
#' `basin_capacity`, `taf`, `shelter_index`, `drainage_bearing`,
#' `flow_convergence`, `slope`, `aspect`) may be passed as a named scalar in
#' `...` to skip its DEM derivation and use the supplied value directly. The
#' corresponding `meta` entry is set to `"user_supplied"`.
#'
#' @return When `receptors = NULL`: a validated `mh_terrain` S7 object with
#'   every applicable scalar descriptor and `@meta` recording the analysis
#'   scale for each derived descriptor plus `dem_resolution` and `dem_source`.
#'   When `receptors` is an `sf` object: a named list with
#'   `$terrain` (the `mh_terrain`) and `$receptor_fields` (a data frame with
#'   columns `feature_id`, `rel_elevation`).
#'
#' @section Setup-time only:
#' This function runs WhiteboxTools flow-routing and openness calculations that
#' take seconds to minutes. Call it once at site setup, save the returned
#' `mh_terrain`, and pass that to [odour_exposure()], [odour_risk()], etc.
#'
#' @section Dependencies:
#' Requires the `whitebox` R package and the WhiteboxTools binary
#' (`whitebox::install_whitebox()`), and the `terra` package.
#'
#' @references
#' Lindsay, J.B. (2015). Geomorphons: a new automated landform classification
#'   method for characterizing the topographic position of land. *Geomorphology*,
#'   238, 18--29.
#' Yokoyama, R., Shirasawa, M., & Pike, R.J. (2002). Visualizing topography by
#'   openness: a new application of image processing to digital elevation models.
#'   *Photogrammetric Engineering and Remote Sensing*, 68(3), 257--265.
#'
#' @export
mh_terrain_from_dem <- function(dem,
                                  source,
                                  receptors  = NULL,
                                  epsg,
                                  conditioning = c("breach", "fill", "none"),
                                  shelter_fetch_m = NULL,
                                  ...) {
  conditioning <- match.arg(conditioning)
  overrides    <- list(...)

  # ---- 1. Input validation (before binary check) --------------------------- #
  # Validate source is an sf object
  if (!inherits(source, "sf")) {
    cli::cli_abort(
      "{.arg source} must be an {.cls sf} object.",
      class = "meteoHazard_input_error"
    )
  }

  # Load and validate DEM (checks: file exists, has CRS, epsg is projected).
  # .dem_ingest() is terra-only; no binary needed.
  dem_r <- .dem_ingest(dem, epsg)

  # ---- 2. Binary check ----------------------------------------------------- #
  .require_dem_stack(need_download = FALSE)

  # ---- 3. Setup ------------------------------------------------------------ #
  # Resolution and source string for meta.
  res_xy       <- terra::res(dem_r)
  dem_res_m    <- max(res_xy)
  dem_source_s <- if (is.character(dem)) dem else "SpatRaster"

  # Reproject source to the DEM CRS.
  src_reproj <- sf::st_transform(source, crs = sf::st_crs(as.integer(epsg)))
  src_pt     <- sf::st_centroid(sf::st_geometry(src_reproj))

  # Scan range for multi-scale derivations.
  scan_range <- .dem_scan_range(dem_r)
  # Warn when the DEM extent is too small for a meaningful multi-scale scan
  # (max scan radius < 15 cells → fewer than 15 scale steps).
  if (scan_range$max_scale / dem_res_m < 15) {
    cli::cli_warn(
      "DEM extent is small; the multi-scale scan range is narrow and results may be unreliable.",
      class = "meteoHazard_scan_clamp_warning"
    )
  }

  # Temp directory for all WBT file I/O (cleaned up on function exit).
  workdir <- withr::local_tempdir(clean = TRUE)

  # Write DEM to disk (WBT is file-based).
  dem_path <- file.path(workdir, "dem.tif")
  terra::writeRaster(dem_r, dem_path, overwrite = TRUE)

  # ---- 4. DEM conditioning (breach/fill for flow routing) ------------------ #
  # Track the method actually applied so it can be recorded in meta.
  cond_path <- file.path(workdir, "dem_cond.tif")
  cond_used <- "raw"   # fallback when all WBT calls fail

  if (conditioning == "breach") {
    tryCatch({
      whitebox::wbt_breach_depressions(
        dem = dem_path, output = cond_path, fill_pits = TRUE
      )
      if (file.exists(cond_path)) cond_used <- "breach"
    }, error = function(e) NULL)

    if (cond_used == "raw") {
      tryCatch({
        whitebox::wbt_fill_depressions_wang_and_liu(
          dem = dem_path, output = cond_path
        )
        if (file.exists(cond_path)) cond_used <- "fill"
      }, error = function(e) NULL)
    }

    if (cond_used == "raw") cond_path <- dem_path

  } else if (conditioning == "fill") {
    tryCatch({
      whitebox::wbt_fill_depressions_wang_and_liu(
        dem = dem_path, output = cond_path
      )
      if (file.exists(cond_path)) cond_used <- "fill"
    }, error = function(e) NULL)
    if (cond_used == "raw") cond_path <- dem_path

  } else {
    cond_path <- dem_path
    cond_used <- "none"
  }

  # ---- 5. Derive scalar descriptors ---------------------------------------- #
  meta <- list(
    dem_resolution = dem_res_m,
    dem_source     = dem_source_s,
    conditioning   = cond_used
  )

  # Apply overrides: if the caller supplied a descriptor, skip derivation.
  descriptor_names <- c("relief", "valley_depth", "basin_capacity",
                         "taf", "shelter_index", "drainage_bearing",
                         "flow_convergence", "slope", "aspect")
  override_applied <- intersect(names(overrides), descriptor_names)
  for (nm in override_applied) {
    meta[[paste0(nm, "_radius")]] <- "user_supplied"
  }

  if ("relief" %in% override_applied) {
    relief <- as.double(overrides$relief)
  } else {
    r_rel  <- .derive_relief(dem_r, src_pt, scan_range, workdir)
    relief <- r_rel$relief
    meta[["relief_radius"]] <- r_rel$relief_radius
  }

  if ("valley_depth" %in% override_applied) {
    valley_depth <- as.double(overrides$valley_depth)
  } else {
    r_vd         <- .derive_valley_depth(dem_r, src_pt, scan_range, workdir)
    valley_depth <- r_vd$valley_depth
    meta[["valley_dev_scale"]] <- r_vd$valley_dev_scale
  }

  if ("flow_convergence" %in% override_applied) {
    flow_convergence <- as.double(overrides$flow_convergence)
    drainage_bearing <- if ("drainage_bearing" %in% override_applied)
                          as.double(overrides$drainage_bearing) else NA_real_
  } else {
    r_flow           <- .derive_flow(dem_r, cond_path, src_pt, scan_range, workdir)
    flow_convergence <- r_flow$flow_convergence
    drainage_bearing <- if ("drainage_bearing" %in% override_applied)
                          as.double(overrides$drainage_bearing)
                        else r_flow$drainage_bearing
    meta[["drainage_catchment_radius"]] <- r_flow$drainage_catchment_radius
    meta[["flow_method"]]               <- r_flow$flow_method
  }

  basin_capacity <- if ("basin_capacity" %in% override_applied)
                      as.double(overrides$basin_capacity)
                    else .derive_basin_capacity(dem_r, src_pt, workdir)

  taf_val <- if ("taf" %in% override_applied) as.double(overrides$taf)
             else .derive_taf(valley_depth, dem_r, src_pt, dem_res_m)

  if ("shelter_index" %in% override_applied) {
    si_val <- as.double(overrides$shelter_index)
  } else {
    r_si   <- .derive_shelter(dem_r, src_pt, scan_range, shelter_fetch_m, workdir)
    si_val <- r_si$shelter_index
    meta[["shelter_fetch_L"]] <- r_si$shelter_fetch_L
  }

  if (all(c("slope", "aspect") %in% override_applied)) {
    slope_val  <- as.double(overrides$slope)
    aspect_val <- as.double(overrides$aspect)
  } else {
    sa         <- .derive_slope_aspect(dem_r, src_pt)
    slope_val  <- if ("slope"  %in% override_applied) as.double(overrides$slope)  else sa$slope
    aspect_val <- if ("aspect" %in% override_applied) as.double(overrides$aspect) else sa$aspect
  }

  terrain <- mh_terrain(
    relief           = relief,
    valley_depth     = valley_depth,
    basin_capacity   = basin_capacity,
    drainage_bearing = drainage_bearing,
    flow_convergence = flow_convergence,
    slope            = slope_val,
    aspect           = aspect_val,
    taf              = taf_val,
    shelter_index    = si_val,
    meta             = meta
  )

  # ---- 6. Per-receptor fields (Stage 3; stub if receptors = NULL) ---------- #
  if (is.null(receptors)) return(terrain)

  rec_fields <- .derive_receptor_fields(dem_r, src_pt, receptors, epsg)
  list(terrain = terrain, receptor_fields = rec_fields)
}


# ---------------------------------------------------------------------------
# Internal descriptor derivers
# ---------------------------------------------------------------------------

# .derive_relief(): DEVmax-selected scale; magnitude = source_elev - regional_min
# Returns list(relief, relief_radius).
.derive_relief <- function(dem_r, src_pt, scan_range, workdir) {
  # Bound the scan range; warn if DEM is too small.
  min_s <- ceiling(scan_range$min_scale)
  max_s <- floor(scan_range$max_scale)
  if (max_s <= min_s) {
    cli::cli_warn(
      "DEM extent is smaller than the minimum scan range; clamping to DEM extent.",
      class = "meteoHazard_scan_clamp_warning"
    )
    max_s <- min_s * 2L
  }

  dem_path  <- file.path(workdir, "dem.tif")
  mag_path  <- file.path(workdir, "dev_mag.tif")
  scl_path  <- file.path(workdir, "dev_scale.tif")

  whitebox::wbt_max_elevation_deviation(
    dem       = dem_path,
    out_mag   = mag_path,
    out_scale = scl_path,
    min_scale = min_s,
    max_scale = max_s,
    step      = max(1L, as.integer((max_s - min_s) / 10L))
  )

  # Extract the DEV-selected scale at the source point.
  scl_r   <- terra::rast(scl_path)
  sel_scl <- terra::extract(scl_r, terra::vect(src_pt))[[1, 2]]
  # Flat DEMs produce NA/0 scale from DEV; clamp to minimum (1 cell).
  if (is.na(sel_scl) || sel_scl <= 0) sel_scl <- ceiling(scan_range$min_scale)
  relief_radius <- sel_scl * terra::res(dem_r)[1]

  # Source elevation and regional minimum within sel_scl radius.
  src_elev  <- terra::extract(dem_r, terra::vect(src_pt))[[1, 2]]
  buf_m     <- relief_radius
  buf       <- sf::st_buffer(src_pt, dist = buf_m)
  regional  <- terra::extract(dem_r, terra::vect(buf))
  reg_min   <- min(regional[[2]], na.rm = TRUE)

  list(relief = max(0, src_elev - reg_min), relief_radius = relief_radius)
}


# .derive_valley_depth(): threshold-light multi-scale DEV valley delineation.
# Returns list(valley_depth, valley_dev_scale).
.derive_valley_depth <- function(dem_r, src_pt, scan_range, workdir) {
  min_s <- ceiling(scan_range$min_scale)
  max_s <- floor(scan_range$max_scale)
  if (max_s <= min_s) max_s <- min_s * 2L

  dem_path  <- file.path(workdir, "dem.tif")
  mag_path  <- file.path(workdir, "vd_mag.tif")
  scl_path  <- file.path(workdir, "vd_scale.tif")

  whitebox::wbt_max_elevation_deviation(
    dem       = dem_path,
    out_mag   = mag_path,
    out_scale = scl_path,
    min_scale = min_s,
    max_scale = max_s,
    step      = max(1L, as.integer((max_s - min_s) / 10L))
  )

  # Extract the DEV-selected scale at the source point.
  scl_r   <- terra::rast(scl_path)
  sel_scl <- terra::extract(scl_r, terra::vect(src_pt))[[1, 2]]
  # Flat DEMs produce NA/0 scale from DEV; clamp to minimum (1 cell).
  if (is.na(sel_scl) || sel_scl <= 0) sel_scl <- ceiling(scan_range$min_scale)
  valley_dev_scale <- sel_scl * terra::res(dem_r)[1]

  # Extract DEV magnitude at source.
  mag_r    <- terra::rast(mag_path)
  dev_val  <- terra::extract(mag_r, terra::vect(src_pt))[[1, 2]]

  buf_m    <- valley_dev_scale
  buf      <- sf::st_buffer(src_pt, dist = buf_m)
  regional <- terra::extract(dem_r, terra::vect(buf))
  reg_max  <- max(regional[[2]], na.rm = TRUE)
  src_elev <- terra::extract(dem_r, terra::vect(src_pt))[[1, 2]]

  # Valley depth = how high the surrounding rim is above the source.
  depth <- max(0, reg_max - src_elev)
  list(valley_depth = depth, valley_dev_scale = valley_dev_scale)
}


# .derive_flow(): D-infinity flow direction + plan curvature at source.
# Returns list(flow_convergence, drainage_bearing, drainage_catchment_radius, flow_method).
.derive_flow <- function(dem_r, cond_path, src_pt, scan_range, workdir) {
  catchment_r  <- scan_range$max_scale / 2
  flat_default <- list(flow_convergence = 0, drainage_bearing = NA_real_,
                       drainage_catchment_radius = catchment_r, flow_method = "flat")

  # On effectively flat DEMs, fill conditioning creates spurious flow paths.
  # Return physical defaults (no convergence, no meaningful bearing).
  dem_range <- diff(range(terra::values(dem_r, na.rm = TRUE)))
  if (is.na(dem_range) || dem_range < 1e-6) return(flat_default)

  pointer_path <- file.path(workdir, "dinf_ptr.tif")
  accum_path   <- file.path(workdir, "dinf_acc.tif")

  whitebox::wbt_d_inf_pointer(dem = cond_path, output = pointer_path)
  if (!file.exists(pointer_path)) return(flat_default)
  whitebox::wbt_d_inf_flow_accumulation(
    input = pointer_path, output = accum_path,
    out_type = "Cells", log = FALSE, clip = FALSE
  )

  ptr_r   <- terra::rast(pointer_path)
  # Single-cell read at the source: approximation. A more defensible estimator
  # would weight by upstream accumulation across the catchment, but the single-
  # cell D-inf pointer is adequate for a screening tool.
  dir_val <- terra::extract(ptr_r, terra::vect(src_pt))[[1, 2]]   # D-inf flow direction, radians

  # Convert D-inf bearing (radians, from east CCW) to met bearing (from north CW).
  drainage_bearing <- if (is.na(dir_val)) NA_real_
                      else (90 - dir_val * 180 / pi) %% 360

  # Plan curvature at source for flow_convergence.
  plan_path <- file.path(workdir, "plan_curv.tif")
  whitebox::wbt_plan_curvature(dem = cond_path, output = plan_path)
  plan_r  <- terra::rast(plan_path)
  curv    <- terra::extract(plan_r, terra::vect(src_pt))[[1, 2]]

  # Accumulation at source, normalised to [0, 1] (log scale, capped).
  acc_r   <- terra::rast(accum_path)
  acc_val <- terra::extract(acc_r, terra::vect(src_pt))[[1, 2]]
  max_acc <- max(terra::values(acc_r, na.rm = TRUE))
  flow_convergence <- if (is.na(acc_val) || max_acc <= 1) 0
                      else log1p(acc_val) / log1p(max_acc)

  list(flow_convergence          = pmax(0, pmin(1, flow_convergence)),
       drainage_bearing          = drainage_bearing,
       drainage_catchment_radius = catchment_r,
       flow_method               = "d_inf")
}


# .derive_basin_capacity(): volume of the depression enclosing the source.
# Fills the original DEM with wbt_fill_depressions_wang_and_liu (pour-point
# exact), computes depth = fill - dem, then uses terra::patches() to localise
# to the source's connected depression. Returns 0 when the source sits on a
# hillslope, open valley, or flat plain (depth at source < 0.01 m).
.derive_basin_capacity <- function(dem_r, src_pt, workdir) {
  dem_path  <- file.path(workdir, "dem.tif")
  fill_path <- file.path(workdir, "fill_bc.tif")

  ok <- tryCatch({
    whitebox::wbt_fill_depressions_wang_and_liu(dem = dem_path, output = fill_path)
    file.exists(fill_path)
  }, error = function(e) FALSE)
  if (!ok) return(0)

  fill_r <- terra::rast(fill_path)
  dem_disk <- terra::rast(dem_path)   # same projection/extent as fill
  depth  <- fill_r - dem_disk

  src_depth <- terra::extract(depth, terra::vect(src_pt))[[1, 2]]
  if (is.na(src_depth) || src_depth < 0.01) return(0)

  # Label connected depression patches; find which one contains the source.
  dep_mask <- depth > 0.01
  patches  <- terra::patches(dep_mask, directions = 8L, zeroAsNA = TRUE)
  src_patch_id <- terra::extract(patches, terra::vect(src_pt))[[1, 2]]
  if (is.na(src_patch_id)) return(0)

  patch_mask <- terra::values(patches) == src_patch_id
  patch_mask[is.na(patch_mask)] <- FALSE
  cell_area <- prod(terra::res(dem_r))
  sum(terra::values(depth)[patch_mask], na.rm = TRUE) * cell_area
}


# .derive_taf(): topographic amplification factor from valley geometry.
# TAF = (flat equivalent volume) / (valley volume) = 1 + valley_depth / h_mix
# Simplified: ratio of atmospheric volume over flat vs valley (Steinacker 1984).
.derive_taf <- function(valley_depth, dem_r, src_pt, dem_res_m) {
  if (is.na(valley_depth) || valley_depth <= 0) return(1.0)
  # Simple Steinacker TAF: ratio of volume above flat vs volume above valley floor.
  h_ref <- ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE
  taf   <- (h_ref + valley_depth) / h_ref
  pmax(1.0, taf)
}


# .derive_shelter(): positive topographic openness via horizon-scan (terra).
# `wbt_openness` is a paid WBT Pro extension tool; this function provides a
# free reimplementation of the same method: Yokoyama et al. (2002),
# 16-azimuth mean of maximum horizon angles.
# Returns list(shelter_index, shelter_fetch_L).
.derive_shelter <- function(dem_r, src_pt, scan_range, shelter_fetch_m, workdir) {
  L <- if (is.null(shelter_fetch_m)) scan_range$max_scale / 4 else shelter_fetch_m
  L <- pmin(L, scan_range$max_scale)
  list(shelter_index = .openness_terra(dem_r, src_pt, L), shelter_fetch_L = L)
}


# .openness_terra(): positive topographic openness in degrees.
# Scans 16 equally-spaced azimuths to distance L, finds the maximum horizon
# elevation angle per azimuth, returns 90 - mean(max_angles).
# A flat site → ~90°; a valley bottom → ~50-65°.
.openness_terra <- function(dem_r, src_pt, L_m) {
  src_coords <- sf::st_coordinates(src_pt)
  src_elev   <- terra::extract(dem_r, terra::vect(src_pt))[[1, 2]]
  if (is.na(src_elev)) return(NA_real_)

  res_m  <- terra::res(dem_r)[1]
  n_pts  <- max(1L, as.integer(ceiling(L_m / res_m)))
  n_az   <- 16L
  azimuths <- (seq_len(n_az) - 1L) * 2 * pi / n_az
  steps  <- seq_len(n_pts) * res_m

  az_rep   <- rep(azimuths, each = n_pts)
  step_rep <- rep(steps,    times = n_az)
  pts_xy   <- cbind(
    src_coords[1] + step_rep * sin(az_rep),
    src_coords[2] + step_rep * cos(az_rep)
  )
  pts_v  <- terra::vect(pts_xy, crs = terra::crs(dem_r))
  elev   <- terra::extract(dem_r, pts_v)[[2]]
  angles <- atan2(elev - src_elev, step_rep) * 180 / pi

  az_idx     <- rep(seq_len(n_az), each = n_pts)
  max_angles <- as.numeric(tapply(angles, az_idx, max, na.rm = TRUE))
  90 - mean(max_angles, na.rm = TRUE)
}


# .derive_slope_aspect(): standard terra terrain at source point.
.derive_slope_aspect <- function(dem_r, src_pt) {
  slope_r  <- terra::terrain(dem_r, v = "slope",  unit = "degrees")
  aspect_r <- terra::terrain(dem_r, v = "aspect", unit = "degrees")
  slope    <- terra::extract(slope_r,  terra::vect(src_pt))[[1, 2]]
  aspect   <- terra::extract(aspect_r, terra::vect(src_pt))[[1, 2]]
  list(slope = slope, aspect = aspect)
}


# .derive_receptor_fields(): rel_elevation per receptor (elevation relative to source).
.derive_receptor_fields <- function(dem_r, src_pt, receptors, epsg) {
  rec_reproj <- sf::st_transform(receptors, crs = sf::st_crs(as.integer(epsg)))
  rec_elev   <- terra::extract(dem_r, terra::vect(rec_reproj))[[2]]
  src_elev   <- terra::extract(dem_r, terra::vect(src_pt))[[1, 2]]
  n_r        <- nrow(rec_reproj)

  data.frame(
    feature_id    = if ("id" %in% names(receptors)) receptors$id else seq_len(n_r),
    rel_elevation = rec_elev - src_elev,
    stringsAsFactors = FALSE
  )
}
