# test-terrain-from-dem.R
# ===========================================================================
# Behaviour spec for mh_terrain_from_dem() — C5 / Phase 2 (#19).
#
# GIS-dependent specs require the whitebox package AND the WhiteboxTools
# binary. All such specs begin with .wbt_skip() (defined in
# helper-synthetic-dem.R). They will SKIP in environments without the full
# GIS stack and PASS once it is present and the function is implemented.
#
# The graceful-skip and input-validation groups can run without the binary.
# Correctness oracles use the analytically-known synthetic DEMs from
# helper-synthetic-dem.R (Gaussian hill, V-valley, flat plain, conical basin,
# eastward-tilted plane).
# ===========================================================================


# ---------------------------------------------------------------------------
# 1. Graceful skip when the GIS stack is absent
#    These tests run in environments where whitebox IS installed but the
#    binary is absent (common CI before whitebox::install_whitebox()).
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — graceful skip when binary is absent", {

  it("raises a classed meteoHazard_dependency_error when the binary is missing", {
    skip_if_not_installed("whitebox")
    skip_if_not_installed("terra")
    skip_if(.wbt_ok(), "Binary present; skip the binary-absent oracle")

    src <- .src_centre()
    dem <- .dem_flat(ncell = 10, extent_m = 500)

    expect_error(
      mh_terrain_from_dem(dem, source = src, epsg = 32755L),
      class = "meteoHazard_dependency_error"
    )
  })

  it("names 'whitebox' (the binary) in the error message", {
    skip_if_not_installed("whitebox")
    skip_if_not_installed("terra")
    skip_if(.wbt_ok(), "Binary present; skip the binary-absent oracle")

    src <- .src_centre()
    dem <- .dem_flat(ncell = 10, extent_m = 500)

    expect_error(
      mh_terrain_from_dem(dem, source = src, epsg = 32755L),
      regexp = "whitebox",
      class  = "meteoHazard_dependency_error"
    )
  })

  it("suggests whitebox::install_whitebox() in the error message", {
    skip_if_not_installed("whitebox")
    skip_if_not_installed("terra")
    skip_if(.wbt_ok(), "Binary present; skip the binary-absent oracle")

    src <- .src_centre()
    dem <- .dem_flat(ncell = 10, extent_m = 500)

    expect_error(
      mh_terrain_from_dem(dem, source = src, epsg = 32755L),
      regexp = "install_whitebox",
      class  = "meteoHazard_dependency_error"
    )
  })

})


# ---------------------------------------------------------------------------
# 2. Input validation (does not require binary; function must exist)
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — input validation", {

  it("raises a classed error for a DEM SpatRaster with no CRS", {
    # A DEM without a CRS cannot be placed in space; the function must reject
    # it before attempting any GIS computation.
    skip_if_not_installed("terra")

    src    <- .src_centre()
    dem_nc <- terra::rast(nrows = 10, ncols = 10,
                          xmin = -500, xmax = 500,
                          ymin = -500, ymax = 500)   # no CRS
    terra::values(dem_nc) <- runif(100, 0, 10)

    expect_error(
      mh_terrain_from_dem(dem_nc, source = src, epsg = 32755L),
      class = "meteoHazard_input_error"
    )
  })

  it("raises a classed error for a non-existent file path", {
    skip_if_not_installed("terra")
    src <- .src_centre()
    expect_error(
      mh_terrain_from_dem("/no/such/file.tif", source = src, epsg = 32755L),
      class = "meteoHazard_input_error"
    )
  })

  it("raises a classed error for a non-projected (geographic) epsg", {
    skip_if_not_installed("terra")
    src <- .src_centre()
    dem <- .dem_flat(ncell = 10, extent_m = 500)
    expect_error(
      mh_terrain_from_dem(dem, source = src, epsg = 4326L),  # WGS84 lon/lat
      class = "meteoHazard_input_error"
    )
  })

  it("raises a classed error when source is not an sf object", {
    skip_if_not_installed("terra")
    dem <- .dem_flat(ncell = 10, extent_m = 500)
    expect_error(
      mh_terrain_from_dem(dem, source = c(0, 0), epsg = 32755L),
      class = "meteoHazard_input_error"
    )
  })

})


# ---------------------------------------------------------------------------
# 3. Return value contract
#    All tests in this group require the full GIS stack.
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — return value contract", {

  it("returns a valid mh_terrain object (passes the C1 validator)", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_gaussian_hill(H = 80, sigma = 200)
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_true(S7::S7_inherits(ter, mh_terrain))
    # C1 validator fires inside mh_terrain(); no separate call needed —
    # construction success is the assertion.
  })

  it("populates meta with dem_resolution and dem_source", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_flat()
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_true("dem_resolution" %in% names(ter@meta))
    expect_true("dem_source"     %in% names(ter@meta))
    expect_true(is.numeric(ter@meta$dem_resolution))
    expect_gt(ter@meta$dem_resolution, 0)
  })

  it("records a scale entry in meta for every scanned descriptor", {
    # Plan §2: every descriptor derived by a multi-scale scan must record its
    # chosen scale in meta (e.g. relief_radius, valley_dev_scale, shelter_fetch_L,
    # drainage_catchment_radius, flow_method).
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_gaussian_hill()
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    # At minimum these keys must be present:
    expected_keys <- c("relief_radius", "valley_dev_scale", "shelter_fetch_L",
                       "drainage_catchment_radius", "dem_resolution", "dem_source")
    present <- expected_keys %in% names(ter@meta)
    expect_true(
      all(present),
      label = paste("Missing meta keys:", paste(expected_keys[!present], collapse = ", "))
    )
  })

  it("returns scalar mh_terrain (no per-receptor frame) when receptors = NULL", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_flat()
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    # Should be a plain mh_terrain, not a list with a receptor frame.
    expect_true(S7::S7_inherits(ter, mh_terrain))
  })

  it("attaches a per-receptor data frame keyed to receptor id when receptors are supplied", {
    .wbt_skip()
    src  <- .src_centre()
    recs <- rbind(.rec_offset(200, 0, "r1"), .rec_offset(0, 300, "r2"))
    dem  <- .dem_gaussian_hill()
    out  <- mh_terrain_from_dem(dem, source = src, receptors = recs, epsg = 32755L)

    # Function returns a list with $terrain (mh_terrain) and $receptor_fields (df)
    expect_true(is.list(out))
    expect_true(S7::S7_inherits(out$terrain, mh_terrain))
    expect_true(is.data.frame(out$receptor_fields))
    expect_true("feature_id"        %in% names(out$receptor_fields))
    expect_true("rel_elevation" %in% names(out$receptor_fields))
    expect_setequal(out$receptor_fields$feature_id, recs$id)
  })

})


# ---------------------------------------------------------------------------
# 4. Descriptor oracles: relief
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — relief descriptor", {

  it("returns relief ≈ H for a Gaussian hill of peak height H", {
    # DEVmax selects the scale at max |DEV|; the magnitude is source_elev −
    # regional_min within that radius. For an isolated Gaussian hill of height H
    # on a flat base, regional_min ≈ 0 at the right radius → relief ≈ H.
    # Tolerance: 20% to accommodate DEV scale selection and smoothing.
    .wbt_skip()
    H   <- 100
    src <- .src_centre()
    dem <- .dem_gaussian_hill(H = H, sigma = 200)
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_equal(ter@relief, H, tolerance = 0.20)   # 20% of H
  })

  it("records relief_radius in meta (positive, in metres)", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_gaussian_hill(H = 80, sigma = 200)
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_true(is.numeric(ter@meta$relief_radius))
    expect_gt(ter@meta$relief_radius, 0)
  })

  it("relief is scale-aware: hill embedded in broad basin yields smaller relief", {
    # The same hill inside a broader basin makes the DEVmax-selected radius
    # smaller (the broad basin dominates at large scales), so the regional_min
    # within the selected short radius is higher → relief decreases.
    .wbt_skip()
    src          <- .src_centre()
    dem_isolated <- .dem_gaussian_hill(H = 80, sigma = 200)
    dem_in_basin <- .dem_hill_in_basin(H_hill = 80, H_basin = 200, sigma = 200)
    ter_isolated <- mh_terrain_from_dem(dem_isolated, source = src, epsg = 32755L)
    ter_in_basin <- mh_terrain_from_dem(dem_in_basin, source = src, epsg = 32755L)
    expect_lt(ter_in_basin@relief, ter_isolated@relief)
  })

  it("relief ≈ 0 on a flat plain", {
    # Flat DEM: DEV is identically 0 everywhere; relief should be near 0.
    # Tolerance: 5m (rounding/floating-point in the raster operations).
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_flat()
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_equal(ter@relief, 0, tolerance = 5)
  })

})


# ---------------------------------------------------------------------------
# 5. Descriptor oracles: drainage_bearing
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — drainage_bearing descriptor", {

  it("returns drainage_bearing ≈ 90° (east) for an eastward-tilted plane", {
    # A plane tilted toward the east: dominant flow direction = east = 90°.
    # Circular-mean of D-inf flow directions at source → ≈ 90°.
    # Tolerance: ±20° to account for boundary effects and D-inf discretisation.
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_tilted_east(slope_m_per_m = 0.05)
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_equal(ter@drainage_bearing, 90, tolerance = 20)
  })

  it("returns drainage_bearing = NA on a flat plain (no preferred direction)", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_flat()
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_true(is.na(ter@drainage_bearing))
  })

  it("records drainage_catchment_radius in meta (positive, in metres)", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_tilted_east()
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_true(is.numeric(ter@meta$drainage_catchment_radius))
    expect_gt(ter@meta$drainage_catchment_radius, 0)
  })

})


# ---------------------------------------------------------------------------
# 6. Descriptor oracles: flow_convergence
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — flow_convergence descriptor", {

  it("flow_convergence is higher at a converging valley than at a diverging spur", {
    # Physics: water flows toward a valley floor (many upstream cells accumulate
    # there) and away from a hill peak (no upstream cells, accumulation = 1).
    #
    # .dem_converging_hollow() is a CLOSED depression: fill conditioning sets the
    # entire bowl to a flat plateau at rim elevation, destroying convergent flow.
    # .dem_v_valley() is an open form — the floor drains to the DEM boundary —
    # so conditioning leaves its convergent structure intact.
    .wbt_skip()
    src       <- .src_centre()
    dem_conv  <- .dem_v_valley(D = 80)
    dem_spur  <- .dem_diverging_spur(H = 50, sigma = 300)
    ter_conv  <- mh_terrain_from_dem(dem_conv, source = src, epsg = 32755L)
    ter_spur  <- mh_terrain_from_dem(dem_spur, source = src, epsg = 32755L)
    expect_gt(ter_conv@flow_convergence, ter_spur@flow_convergence)
  })

  it("flow_convergence is ≥ 0 on a flat plain (no divergence)", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_flat()
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_gte(ter@flow_convergence, 0)
  })

})


# ---------------------------------------------------------------------------
# 7. Descriptor oracles: valley_depth (multi-scale DEV, no channel threshold)
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — valley_depth (threshold-free multi-scale DEV)", {

  it("valley_depth ≈ D for a V-valley of known depth D", {
    # Source at the valley floor (0, 0). The valley rim is at height D.
    # Multi-scale DEV delineates the floor; valley_depth = source_elev −
    # valley-floor_elev ≈ D.
    # Tolerance: 20% to accommodate DEV scale choice and DEM resolution.
    .wbt_skip()
    D   <- 80
    src <- .src_centre()
    dem <- .dem_v_valley(D = D)
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_equal(ter@valley_depth, D, tolerance = 0.20)
  })

  it("records valley_dev_scale (in metres) in meta", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_v_valley(D = 80)
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_true(is.numeric(ter@meta$valley_dev_scale))
    expect_gt(ter@meta$valley_dev_scale, 0)
  })

  it("valley_depth is stable regardless of DEM conditioning (no channel threshold)", {
    # Re-run with different DEM conditioning (pre-breached vs. raw). Since
    # multi-scale DEV has no flow-accumulation threshold, the delineated valley
    # floor should change by < 10% across conditioning methods.
    # The implementation exposes conditioning via `conditioning = c("breach","fill","none")`.
    .wbt_skip()
    D   <- 80
    src <- .src_centre()
    dem <- .dem_v_valley(D = D)
    ter_breach <- mh_terrain_from_dem(dem, source = src, epsg = 32755L,
                                      conditioning = "breach")
    ter_fill   <- mh_terrain_from_dem(dem, source = src, epsg = 32755L,
                                      conditioning = "fill")
    # valley_depth should be within 10% of D regardless of conditioning choice
    expect_equal(ter_breach@valley_depth, ter_fill@valley_depth, tolerance = 0.10)
  })

  it("valley_depth ≈ 0 on a flat plain", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_flat()
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_equal(ter@valley_depth, 0, tolerance = 5)  # within 5m rounding
  })

})


# ---------------------------------------------------------------------------
# 8. Descriptor oracles: basin_capacity
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — basin_capacity", {

  it("basin_capacity ≈ pi * R^2 * D / 3 for a conical basin of depth D, radius R", {
    # Analytic volume to sill for a conical depression (linear z vs r up to R):
    #   V = integral_0^R 2*pi*r*(D - D*r/R) dr = 2*pi*D * (R^2/2 - R^2/3)
    #     = 2*pi*D * R^2/6 = pi * R^2 * D / 3.
    # Tolerance: 20% to account for DEM cell discretisation and fill algorithm.
    .wbt_skip()
    D   <- 50; R <- 400
    V_analytic <- pi * R^2 * D / 3   # ≈ 8,377,580 m^3
    src <- .src_centre()
    dem <- .dem_basin(D = D, R = R)
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_equal(ter@basin_capacity, V_analytic, tolerance = 0.20)
  })

  it("basin_capacity is non-negative on a flat plain (degenerate depression)", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_flat()
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_gte(ter@basin_capacity, 0)
  })

})


# ---------------------------------------------------------------------------
# 9. Descriptor oracles: taf and shelter_index
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — taf", {

  it("taf >= 1 always (C1 invariant)", {
    .wbt_skip()
    for (dem_fn in list(
      .dem_flat(), .dem_gaussian_hill(), .dem_v_valley(), .dem_basin()
    )) {
      ter <- mh_terrain_from_dem(dem_fn, source = .src_centre(), epsg = 32755L)
      if (!is.na(ter@taf))
        expect_gte(ter@taf, 1,
                   label = paste("taf ≥ 1 violated for DEM of class",
                                 class(dem_fn)))
    }
  })

  it("taf ≈ 1 on a flat plain", {
    .wbt_skip()
    ter <- mh_terrain_from_dem(.dem_flat(), source = .src_centre(), epsg = 32755L)
    if (!is.na(ter@taf)) expect_equal(ter@taf, 1, tolerance = 0.05)
  })

  it("taf is larger for a deeper valley than a shallow one", {
    .wbt_skip()
    ter_deep    <- mh_terrain_from_dem(.dem_v_valley(D = 100),
                                        source = .src_centre(), epsg = 32755L)
    ter_shallow <- mh_terrain_from_dem(.dem_v_valley(D = 20),
                                        source = .src_centre(), epsg = 32755L)
    if (!is.na(ter_deep@taf) && !is.na(ter_shallow@taf))
      expect_gte(ter_deep@taf, ter_shallow@taf)
  })

})

describe("mh_terrain_from_dem() — shelter_index", {

  it("shelter_index(enclosed basin) < shelter_index(open plain)", {
    # Topographic openness: ~90° for an open plain, smaller in an enclosed basin.
    .wbt_skip()
    ter_basin <- mh_terrain_from_dem(.dem_basin(D = 50, R = 400),
                                      source = .src_centre(), epsg = 32755L)
    ter_plain <- mh_terrain_from_dem(.dem_flat(),
                                      source = .src_centre(), epsg = 32755L)
    if (!is.na(ter_basin@shelter_index) && !is.na(ter_plain@shelter_index))
      expect_lt(ter_basin@shelter_index, ter_plain@shelter_index)
  })

  it("records shelter_fetch_L in meta", {
    .wbt_skip()
    ter <- mh_terrain_from_dem(.dem_gaussian_hill(), source = .src_centre(),
                                epsg = 32755L)
    expect_true("shelter_fetch_L" %in% names(ter@meta))
  })

  it("shelter_fetch_L in meta matches a requested fetch radius", {
    # If the caller supplies a fetch radius via shelter_fetch_m, meta should
    # record the same value.
    .wbt_skip()
    requested_L <- 500
    ter <- mh_terrain_from_dem(.dem_gaussian_hill(), source = .src_centre(),
                                epsg = 32755L, shelter_fetch_m = requested_L)
    expect_equal(ter@meta$shelter_fetch_L, requested_L, tolerance = 1)
  })

})


# ---------------------------------------------------------------------------
# 10. Flat-DEM edge case (degenerate input)
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — flat DEM edge case", {

  it("returns valid (non-error) mh_terrain with near-zero descriptors on a flat DEM", {
    .wbt_skip()
    src <- .src_centre()
    ter <- mh_terrain_from_dem(.dem_flat(), source = src, epsg = 32755L)
    expect_true(S7::S7_inherits(ter, mh_terrain))

    # Numeric descriptors that exist must be non-negative and near zero.
    if (!is.na(ter@relief))           expect_equal(ter@relief,        0, tolerance = 5)
    if (!is.na(ter@valley_depth))     expect_equal(ter@valley_depth,  0, tolerance = 5)
    if (!is.na(ter@basin_capacity))   expect_gte(ter@basin_capacity,  0)
    if (!is.na(ter@flow_convergence)) expect_equal(ter@flow_convergence, 0, tolerance = 0.1)
    if (!is.na(ter@taf))              expect_equal(ter@taf,           1, tolerance = 0.05)
    # drainage_bearing should be NA on a flat plain (no preferred direction).
    expect_true(is.na(ter@drainage_bearing))
  })

})


# ---------------------------------------------------------------------------
# 11. Reproducibility
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — reproducibility", {

  it("two identical calls are bit-identical on all numeric descriptor fields", {
    .wbt_skip()
    src  <- .src_centre()
    dem  <- .dem_gaussian_hill()
    ter1 <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    ter2 <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    expect_identical(ter1@relief,            ter2@relief)
    expect_identical(ter1@valley_depth,      ter2@valley_depth)
    expect_identical(ter1@basin_capacity,    ter2@basin_capacity)
    expect_identical(ter1@drainage_bearing,  ter2@drainage_bearing)
    expect_identical(ter1@flow_convergence,  ter2@flow_convergence)
    expect_identical(ter1@taf,               ter2@taf)
    expect_identical(ter1@shelter_index,     ter2@shelter_index)
  })

})


# ---------------------------------------------------------------------------
# 12. Caller overrides
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — caller overrides compose", {

  it("a caller-supplied relief is used as-is and meta flags it as user_supplied", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_gaussian_hill(H = 80)
    ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L, relief = 42)
    expect_equal(ter@relief, 42)
    expect_equal(ter@meta$relief_radius, "user_supplied")
  })

  it("other descriptors are still derived when one is overridden", {
    .wbt_skip()
    src <- .src_centre()
    dem <- .dem_gaussian_hill(H = 80)
    # Override only relief; valley_depth should still be DEM-derived (not 0 or NA).
    ter_override <- mh_terrain_from_dem(dem, source = src, epsg = 32755L, relief = 42)
    ter_derived  <- mh_terrain_from_dem(dem, source = src, epsg = 32755L)
    # Other descriptors should match the fully-derived version.
    expect_equal(ter_override@valley_depth,     ter_derived@valley_depth,    tolerance = 1e-6)
    expect_equal(ter_override@drainage_bearing, ter_derived@drainage_bearing, tolerance = 1e-6)
    expect_equal(ter_override@shelter_index,    ter_derived@shelter_index,    tolerance = 1e-6)
  })

})


# ---------------------------------------------------------------------------
# 13. CRS handling
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — CRS reprojection", {

  it("a geographic-CRS DEM yields the same descriptors as a pre-projected DEM", {
    # When the DEM is in WGS84 (EPSG:4326) and epsg = 32755, the function
    # should reproject internally. Descriptors should match those from a DEM
    # already in EPSG:32755, within 10% (reprojection resampling tolerance).
    #
    # NOTE: this test builds its own DEM at a valid EPSG:32755 coordinate
    # (near Sydney, ~335000 E / 6240000 N) rather than using .src_centre() /
    # .dem_gaussian_hill(), because the (0, 0) UTM origin is 10 000 km south
    # of the equator (invalid) and the round-trip WGS84 projection there
    # introduces a ~125 m systematic shift that is not a code bug.
    .wbt_skip()
    cx <- 335000L; cy <- 6240000L          # valid zone-55S coordinate, near Sydney
    src <- sf::st_sf(
      id       = "src",
      geometry = sf::st_sfc(sf::st_point(c(cx, cy)), crs = 32755L)
    )
    dem_metric <- terra::rast(
      nrows = 64L, ncols = 64L,
      xmin  = cx - 1000L, xmax = cx + 1000L,
      ymin  = cy - 1000L, ymax = cy + 1000L,
      crs   = "EPSG:32755"
    )
    xy <- terra::xyFromCell(dem_metric, seq_len(terra::ncell(dem_metric)))
    terra::values(dem_metric) <-
      80 * exp(-((xy[, "x"] - cx)^2 + (xy[, "y"] - cy)^2) / (2 * 200^2))
    dem_geographic <- terra::project(dem_metric, "EPSG:4326")

    ter_metric     <- mh_terrain_from_dem(dem_metric,     source = src, epsg = 32755L)
    ter_geographic <- mh_terrain_from_dem(dem_geographic, source = src, epsg = 32755L)

    # Bilinear resampling through a valid UTM location should agree within 10%.
    expect_equal(ter_geographic@relief, ter_metric@relief, tolerance = 0.10)
  })

})


# ---------------------------------------------------------------------------
# 14. Per-receptor fields
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — per-receptor fields", {

  it("rel_elevation > 0 for a receptor above the source on a Gaussian hill flank", {
    # Source at the peak (0, 0): DEM elevation = H.
    # Receptor at base of hill (500 m away): DEM elevation ≈ 0.
    # The source is HIGHER, so rel_elevation = rec_elev - src_elev < 0.
    # But a receptor placed ABOVE the source should give rel_elevation > 0.
    # On a Gaussian hill the source is at the peak, so all receptors are LOWER.
    # Instead, use a valley where source is at the floor and receptor is upslope.
    .wbt_skip()
    src  <- .src_centre()          # valley floor: lowest point
    rec  <- .rec_offset(0, 300)   # 300m north, up the hillside
    dem  <- .dem_v_valley(D = 80)
    out  <- mh_terrain_from_dem(dem, source = src, receptors = rec, epsg = 32755L)
    # rec at y=300 on V-valley: z = 80 * 300/1000 = 24 m; source at y=0: z = 0.
    # rel_elevation = rec_elev - src_elev = 24 - 0 = 24 > 0
    expect_gt(out$receptor_fields$rel_elevation[1], 0)
  })

  it("rel_elevation < 0 for a receptor in a valley below the source", {
    # Source on a hill flank at positive elevation; receptor at valley floor.
    # Use the V-valley: place source 500 m north (upslope), receptor at origin.
    .wbt_skip()
    src_upslope <- sf::st_sf(id = "src",
                              geometry = sf::st_sfc(sf::st_point(c(0, 500)), crs = 32755))
    rec_floor   <- sf::st_sf(id = "rec",
                              geometry = sf::st_sfc(sf::st_point(c(0, 0)),   crs = 32755))
    dem <- .dem_v_valley(D = 80)
    out <- mh_terrain_from_dem(dem, source = src_upslope, receptors = rec_floor,
                                epsg = 32755L)
    # src at y=500: z = 80*500/1000 = 40m; rec at y=0: z = 0. rel = 0 - 40 = -40 < 0.
    expect_lt(out$receptor_fields$rel_elevation[1], 0)
  })

  it("receptor data frame is keyed to receptor feature id", {
    .wbt_skip()
    src  <- .src_centre()
    recs <- rbind(.rec_offset(200, 0, "alpha"), .rec_offset(-200, 0, "beta"))
    dem  <- .dem_flat()
    out  <- mh_terrain_from_dem(dem, source = src, receptors = recs, epsg = 32755L)
    expect_setequal(out$receptor_fields$feature_id, c("alpha", "beta"))
  })

})


# ---------------------------------------------------------------------------
# 15. DEM extent edge case: scan range exceeds DEM
# ---------------------------------------------------------------------------

describe("mh_terrain_from_dem() — DEM smaller than scan range", {

  it("clamps to DEM extent and records achieved scale with a cli warning", {
    # A very small DEM (200m x 200m) cannot support a scan up to the default
    # maximum scale. The function must clamp, warn, and record the achieved
    # (not requested) scale in meta.
    .wbt_skip()
    src <- sf::st_sf(id = "src",
                     geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 32755))
    dem <- .make_dem(function(x, y) 10 * exp(-(x^2 + y^2) / (2 * 40^2)),
                     ncell = 20, extent_m = 200)  # tiny DEM
    expect_warning(
      ter <- mh_terrain_from_dem(dem, source = src, epsg = 32755L),
      regexp = "scan"  # the warning mentions 'scan range' being clamped
    )
    # Output is still a valid mh_terrain, not an error.
    expect_true(S7::S7_inherits(ter, mh_terrain))
  })

})
