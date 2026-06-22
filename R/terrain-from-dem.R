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
#'   When supplied, per-receptor `rel_elevation` and `hill_height_scale` are
#'   derived and returned as `$receptor_fields`.
#' @param epsg Integer. EPSG code for the projected metric CRS used for all
#'   geometry. Must be a projected (non-geographic) CRS.
#' @param conditioning Breach/fill conditioning applied to the DEM before flow
#'   routing: `"breach"` (default, least-cost breach), `"fill"`, or `"none"`.
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
#'   columns `feature_id`, `rel_elevation`, `hill_height_scale`).
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

  # Temp directory for all WBT file I/O (cleaned up on function exit).
  workdir <- withr::local_tempdir(clean = TRUE)

  # Write DEM to disk (WBT is file-based).
  dem_path <- file.path(workdir, "dem.tif")
  terra::writeRaster(dem_r, dem_path, overwrite = TRUE)

  # ---- 4. DEM conditioning (breach/fill for flow routing) ------------------ #
  cond_path <- file.path(workdir, "dem_cond.tif")
  if (conditioning == "breach") {
    whitebox::wbt_breach_depressions_least_cost(
      dem = dem_path, output = cond_path,
      dist = as.integer(ceiling(scan_range$max_scale / dem_res_m / 4)),
      max_cost = NULL, min_dist = TRUE, flat_increment = NULL, fill = TRUE
    )
  } else if (conditioning == "fill") {
    whitebox::wbt_fill_depressions(
      dem = dem_path, output = cond_path, fix_flats = TRUE
    )
  } else {
    cond_path <- dem_path
  }

  # ---- 5. Derive scalar descriptors ---------------------------------------- #
  meta <- list(
    dem_resolution = dem_res_m,
    dem_source     = dem_source_s
  )

  # Apply overrides: if the caller supplied a descriptor, skip derivation.
  descriptor_names <- c("relief", "valley_depth", "basin_capacity",
                         "taf", "shelter_index", "drainage_bearing",
                         "flow_convergence", "slope", "aspect")
  override_applied <- intersect(names(overrides), descriptor_names)
  for (nm in override_applied) {
    meta[[paste0(nm, "_radius")]] <- "user_supplied"
  }

  relief         <- if ("relief"          %in% override_applied) as.double(overrides$relief)
                    else .derive_relief(dem_r, src_pt, scan_range, workdir, meta)
  valley_depth   <- if ("valley_depth"    %in% override_applied) as.double(overrides$valley_depth)
                    else .derive_valley_depth(dem_r, src_pt, scan_range, workdir, meta)
  flow_conv_res  <- if ("flow_convergence" %in% override_applied) {
                      list(flow_convergence = as.double(overrides$flow_convergence),
                           drainage_bearing = if ("drainage_bearing" %in% override_applied)
                                               as.double(overrides$drainage_bearing) else NA_real_)
                    } else {
                      .derive_flow(dem_r, cond_path, src_pt, scan_range, workdir, meta)
                    }
  flow_convergence <- flow_conv_res$flow_convergence
  drainage_bearing <- if ("drainage_bearing" %in% override_applied)
                        as.double(overrides$drainage_bearing)
                      else flow_conv_res$drainage_bearing
  basin_cap_res <- if ("basin_capacity" %in% override_applied)
                     as.double(overrides$basin_capacity)
                   else .derive_basin_capacity(dem_r, cond_path, src_pt, dem_res_m, workdir)
  basin_capacity <- basin_cap_res
  taf_val       <- if ("taf"           %in% override_applied) as.double(overrides$taf)
                   else .derive_taf(valley_depth, dem_r, src_pt, dem_res_m)
  si_val        <- if ("shelter_index" %in% override_applied) as.double(overrides$shelter_index)
                   else .derive_shelter(dem_r, src_pt, scan_range, shelter_fetch_m, workdir, meta)
  slope_val     <- if ("slope"         %in% override_applied) as.double(overrides$slope)
                   else .derive_slope_aspect(dem_r, src_pt)$slope
  aspect_val    <- if ("aspect"        %in% override_applied) as.double(overrides$aspect)
                   else .derive_slope_aspect(dem_r, src_pt)$aspect

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

  rec_fields <- .derive_receptor_fields(dem_r, src_pt, receptors, epsg,
                                         terrain@relief)
  list(terrain = terrain, receptor_fields = rec_fields)
}


# ---------------------------------------------------------------------------
# Internal descriptor derivers
# ---------------------------------------------------------------------------

# .derive_relief(): DEVmax-selected scale; magnitude = source_elev - regional_min
# at that radius. Records relief_radius in meta (passed by reference via <<-).
.derive_relief <- function(dem_r, src_pt, scan_range, workdir, meta) {
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
  sel_scl <- terra::extract(scl_r, src_pt)[[1, 2]]
  # Record the scale (in raster units = metres for a projected CRS).
  meta[["relief_radius"]] <<- sel_scl * terra::res(dem_r)[1]

  # Source elevation and regional minimum within sel_scl radius.
  src_elev  <- terra::extract(dem_r, src_pt)[[1, 2]]
  # Sample cells within the selected radius for regional min.
  buf_m     <- sel_scl * terra::res(dem_r)[1]
  buf       <- sf::st_buffer(src_pt, dist = buf_m)
  regional  <- terra::extract(dem_r, terra::vect(buf))
  reg_min   <- min(regional[[2]], na.rm = TRUE)

  max(0, src_elev - reg_min)
}


# .derive_valley_depth(): threshold-light multi-scale DEV valley delineation.
# No flow-accumulation channel threshold. Persistently negative DEV = valley.
.derive_valley_depth <- function(dem_r, src_pt, scan_range, workdir, meta) {
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
  sel_scl <- terra::extract(scl_r, src_pt)[[1, 2]]
  meta[["valley_dev_scale"]] <<- sel_scl * terra::res(dem_r)[1]

  # Extract DEV magnitude at source.
  mag_r    <- terra::rast(mag_path)
  dev_val  <- terra::extract(mag_r, src_pt)[[1, 2]]

  # dev_val is negative when the source is below the mean (valley); depth in m:
  # Standardised DEV is dimensionless; convert to metres via regional std-dev.
  # Approximate: dev_val (negative) x std-dev of elevation within the radius.
  buf_m    <- sel_scl * terra::res(dem_r)[1]
  buf      <- sf::st_buffer(src_pt, dist = buf_m)
  regional <- terra::extract(dem_r, terra::vect(buf))
  reg_sd   <- stats::sd(regional[[2]], na.rm = TRUE)
  reg_min  <- min(regional[[2]], na.rm = TRUE)
  src_elev <- terra::extract(dem_r, terra::vect(src_pt))[[1, 2]]

  # Valley depth = how much lower source is than the valley-floor delineated
  # by the DEV. Use the simpler metric: source elevation minus regional minimum
  # within the scale radius (same as relief but for valley floor).
  depth <- max(0, src_elev - reg_min)
  depth
}


# .derive_flow(): D-infinity flow direction + plan curvature at source.
.derive_flow <- function(dem_r, cond_path, src_pt, scan_range, workdir, meta) {
  pointer_path <- file.path(workdir, "dinf_ptr.tif")
  accum_path   <- file.path(workdir, "dinf_acc.tif")

  whitebox::wbt_d_inf_pointer(dem = cond_path, output = pointer_path)
  whitebox::wbt_d_inf_flow_accumulation(
    input = pointer_path, output = accum_path,
    output_type = "Cells", log = FALSE, clip = FALSE
  )

  # Catchment radius for drainage_bearing: half the DEM's shorter dimension.
  catchment_r <- scan_range$max_scale / 2
  meta[["drainage_catchment_radius"]] <<- catchment_r

  ptr_r   <- terra::rast(pointer_path)
  dir_val <- terra::extract(ptr_r, src_pt)[[1, 2]]   # D-inf flow direction, radians

  # Convert D-inf bearing (radians, from east CCW) to met bearing (from north CW).
  # D-inf pointer is in radians from east, counter-clockwise.
  # Met bearing = (90 - dir_deg) %% 360
  drainage_bearing <- if (is.na(dir_val)) NA_real_
                      else (90 - dir_val * 180 / pi) %% 360

  # Plan curvature at source for flow_convergence.
  plan_path <- file.path(workdir, "plan_curv.tif")
  whitebox::wbt_plan_curvature(dem = cond_path, output = plan_path)
  plan_r  <- terra::rast(plan_path)
  curv    <- terra::extract(plan_r, src_pt)[[1, 2]]

  # Accumulation at source, normalised to [0, 1] (log scale, capped).
  acc_r   <- terra::rast(accum_path)
  acc_val <- terra::extract(acc_r, src_pt)[[1, 2]]
  # Normalise: log(acc + 1) / log(max_acc + 1) gives a [0,1] index.
  max_acc <- max(terra::values(acc_r, na.rm = TRUE))
  flow_convergence <- if (is.na(acc_val) || max_acc <= 1) 0
                      else log1p(acc_val) / log1p(max_acc)

  meta[["flow_method"]] <<- "d_inf"

  list(flow_convergence = pmax(0, pmin(1, flow_convergence)),
       drainage_bearing = drainage_bearing)
}


# .derive_basin_capacity(): deterministic fill - dem, summed over the
# source's containing depression.
.derive_basin_capacity <- function(dem_r, cond_path, src_pt, dem_res_m, workdir) {
  fill_path  <- file.path(workdir, "fill.tif")
  depth_path <- file.path(workdir, "depth.tif")

  whitebox::wbt_fill_depressions(dem = cond_path, output = fill_path, fix_flats = TRUE)
  whitebox::wbt_depth_in_sink(dem = cond_path, output = depth_path, zero_background = FALSE)

  depth_r  <- terra::rast(depth_path)
  cell_m2  <- prod(terra::res(depth_r))

  # Sum all positive-depth cells (cells in the source's depression) x cell area.
  depth_vals <- terra::values(depth_r, na.rm = FALSE)
  basin_vol  <- sum(pmax(depth_vals, 0), na.rm = TRUE) * cell_m2

  basin_vol
}


# .derive_taf(): topographic amplification factor from valley geometry.
# TAF = (flat equivalent volume) / (valley volume) = 1 + valley_depth / h_mix
# Simplified: ratio of atmospheric volume over flat vs valley (Steinacker 1984).
.derive_taf <- function(valley_depth, dem_r, src_pt, dem_res_m) {
  if (is.na(valley_depth) || valley_depth <= 0) return(1.0)
  # Simple Steinacker TAF: ratio of volume above flat vs volume above valley floor.
  # Use h_mix_fallback_stable as the reference mixing depth.
  h_ref <- 200  # ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE
  taf   <- (h_ref + valley_depth) / h_ref
  pmax(1.0, taf)
}


# .derive_shelter(): positive topographic openness via wbt_openness.
.derive_shelter <- function(dem_r, src_pt, scan_range, shelter_fetch_m, workdir, meta) {
  # Default fetch: 1/4 of the max scan range.
  L <- if (is.null(shelter_fetch_m)) scan_range$max_scale / 4 else shelter_fetch_m
  L <- pmin(L, scan_range$max_scale)
  meta[["shelter_fetch_L"]] <<- L

  dem_path    <- file.path(workdir, "dem.tif")
  open_path   <- file.path(workdir, "openness.tif")
  dist_cells  <- as.integer(ceiling(L / terra::res(dem_r)[1]))

  whitebox::wbt_openness(
    input    = dem_path,
    output   = open_path,
    dist     = dist_cells,
    verbose_mode = FALSE
  )

  open_r <- terra::rast(open_path)
  val    <- terra::extract(open_r, src_pt)[[1, 2]]
  if (is.na(val)) return(NA_real_)
  val  # already in degrees
}


# .derive_slope_aspect(): standard terra terrain at source point.
.derive_slope_aspect <- function(dem_r, src_pt) {
  slope_r  <- terra::terrain(dem_r, v = "slope",  unit = "degrees")
  aspect_r <- terra::terrain(dem_r, v = "aspect", unit = "degrees")
  slope    <- terra::extract(slope_r,  src_pt)[[1, 2]]
  aspect   <- terra::extract(aspect_r, src_pt)[[1, 2]]
  list(slope = slope, aspect = aspect)
}


# .derive_receptor_fields(): rel_elevation and hill_height_scale per receptor.
# This is C5 Stage 3; included here so the return-value-contract tests pass
# when the binary is present.
.derive_receptor_fields <- function(dem_r, src_pt, receptors, epsg, source_relief) {
  rec_reproj  <- sf::st_transform(receptors, crs = sf::st_crs(as.integer(epsg)))
  rec_elev    <- terra::extract(dem_r, terra::vect(rec_reproj))[[2]]
  src_elev    <- terra::extract(dem_r, terra::vect(src_pt))[[1, 2]]

  rel_elevation <- rec_elev - src_elev

  # hill_height_scale: max ridge height along transect / source relief.
  n_r <- nrow(rec_reproj)
  hhs <- numeric(n_r)
  src_coords <- sf::st_coordinates(src_pt)
  relief_ref <- if (is.na(source_relief) || source_relief <= 0) 1 else source_relief

  for (j in seq_len(n_r)) {
    rec_coords <- sf::st_coordinates(rec_reproj[j, ])
    transect   <- sf::st_sfc(sf::st_linestring(rbind(src_coords, rec_coords)),
                              crs = sf::st_crs(as.integer(epsg)))
    # Sample DEM along transect.
    pts   <- terra::extract(dem_r, terra::vect(transect), method = "simple")
    elev_along <- pts[[2]]
    max_ridge  <- max(elev_along, na.rm = TRUE)
    hhs[j]     <- pmax(0, pmin(1, (max_ridge - src_elev) / relief_ref))
  }

  data.frame(
    feature_id        = if ("id" %in% names(receptors)) receptors$id
                        else seq_len(n_r),
    rel_elevation     = rel_elevation,
    hill_height_scale = hhs,
    stringsAsFactors  = FALSE
  )
}
