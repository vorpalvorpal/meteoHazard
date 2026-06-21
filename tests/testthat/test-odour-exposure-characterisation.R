# Characterisation (golden) tests for odour_exposure() and the vectorised
# geometry/pathway helpers.
#
# PURPOSE: these tests pin the EXACT numeric output of the hot path so that a
# performance refactor (vectorising the source/receptor loop into matrix ops)
# can be proven behaviour-preserving. They are deliberately value-exact, unlike
# the relational specs in test-odour-exposure.R / test-odour-pathways.R.
#
# The snapshots in _snaps/ were generated from the pre-restructure
# implementation and are the regression oracle. A genuine algorithmic change
# (not a refactor) is expected to update them deliberately.

# ---------------------------------------------------------------------------
# A rich, fully deterministic scenario that exercises the whole hot path in
# two dimensions: 3 sources (2 point at different emit heights + 1 polygon
# area source), 6 receptors at varied bearings/distances, and 24 hours of met
# spanning night/day/calm/windy with multi-level winds and a morning CBL jump.
# ---------------------------------------------------------------------------

.char_site <- function(with_terrain = FALSE) {
  ox <- 335000
  oy <- 6250000

  # Polygon area source: 120 m square centred at the origin.
  hw <- 60
  poly <- sf::st_polygon(list(matrix(
    c(ox - hw, oy - hw,
      ox + hw, oy - hw,
      ox + hw, oy + hw,
      ox - hw, oy + hw,
      ox - hw, oy - hw),
    ncol = 2, byrow = TRUE
  )))

  # Two point sources, offset from the origin.
  src_p1 <- sf::st_point(c(ox + 40, oy - 20))
  src_p2 <- sf::st_point(c(ox - 30, oy + 35))

  # Six receptors at varied bearings (deg from N, clockwise) and distances.
  rec_specs <- list(
    c(b =   0, d =  300),
    c(b =  60, d =  600),
    c(b = 120, d =  900),
    c(b = 180, d =  450),
    c(b = 240, d = 1200),
    c(b = 300, d =  750)
  )
  rec_geoms <- lapply(rec_specs, function(s) {
    sf::st_point(c(ox + s["d"] * sin(s["b"] * pi / 180),
                   oy + s["d"] * cos(s["b"] * pi / 180)))
  })

  ids   <- c("poly", "pt1", "pt2",
             "r1", "r2", "r3", "r4", "r5", "r6")
  geoms <- sf::st_sfc(c(list(poly, src_p1, src_p2), rec_geoms), crs = 32755)
  feats <- sf::st_sf(id = ids, geometry = geoms)

  roles <- data.frame(
    feature_id  = ids,
    hazard      = "odour",
    role        = c("source", "source", "source",
                    rep("receptor", 6)),
    emit_height = c(8, 3, 15, rep(NA_real_, 6)),
    stringsAsFactors = FALSE
  )

  ter <- if (with_terrain)
    mh_terrain(drainage_bearing = 30, flow_convergence = 0.75) else NULL
  mh_site(feats, roles, terrain = ter, epsg = 32755L)
}

.char_met <- function() {
  n <- 24
  # Night (1-8), day (9-20), night (21-24); morning onset at hour 9.
  direct <- c(rep(0, 8),
              50, 150, 300, 450, 600, 650, 600, 450, 300, 150, 50, 0,
              0, 0, 0, 0)
  # BLH: shallow overnight, jumps at the morning onset (hour 9) to make
  # cbl_growth[9] > 0, then grows and collapses again at dusk.
  blh <- c(rep(150, 8),
           400, 600, 800, 1000, 1300, 1500, 1500, 1300, 1000, 700, 400, 200,
           150, 150, 150, 150)
  # Wind direction sweeps the full circle so the downwind axis crosses every
  # receptor over the day.
  wdir <- (seq(0, 345, length.out = n)) %% 360
  # Mostly breezy, with one calm hour (12) to exercise the calm branch.
  wspd <- rep(3, n); wspd[12] <- 0.3
  # Multi-level winds present overnight so residual_wind is non-NA and the
  # fumigation pathway fires in the morning.
  data.frame(
    wind_direction_10m     = wdir,
    wind_speed_10m         = wspd,
    direct_radiation       = direct,
    cloud_cover            = rep(40, n),
    boundary_layer_height  = blh,
    temperature_2m         = 8 + 10 * (direct / 650),
    pressure_msl           = rep(1011, n),
    precipitation          = rep(0, n),
    relative_humidity_2m   = 80 - 20 * (direct / 650),
    soil_moisture_0_to_1cm = rep(0.12, n),
    soil_moisture_1_to_3cm = rep(0.12, n),
    wind_speed_80m         = rep(5, n),
    wind_direction_80m     = rep(250, n),
    wind_speed_120m        = rep(6, n),
    wind_direction_120m    = rep(255, n),
    wind_speed_180m        = rep(7, n),
    wind_direction_180m    = rep(260, n)
  )
}

describe("odour_exposure() golden output (regression oracle for the refactor)", {

  it("matches the pinned flat-backend output on the 3-source x 6-receptor x 24-hour grid", {
    site <- .char_site(with_terrain = FALSE)
    out  <- odour_exposure(.char_met(), site, terrain_backend = "none")
    expect_length(out, 24)
    expect_snapshot_value(round(out, 8), style = "json2", tolerance = 1e-6)
  })

  it("matches the pinned descriptors-backend output on the same grid", {
    site <- .char_site(with_terrain = TRUE)
    out  <- odour_exposure(.char_met(), site, terrain_backend = "descriptors")
    expect_length(out, 24)
    expect_snapshot_value(round(out, 8), style = "json2", tolerance = 1e-6)
  })

  it("matches the pinned flat output for a single polygon source over swept wind", {
    # Isolates the vectorised .crosswind_halfwidth() path inside the exposure
    # loop: one polygon source, swept wind direction over 24 hours.
    ox <- 335000; oy <- 6250000; hw <- 75
    poly <- sf::st_polygon(list(matrix(
      c(ox-hw, oy-hw, ox+hw, oy-hw, ox+hw, oy+hw, ox-hw, oy+hw, ox-hw, oy-hw),
      ncol = 2, byrow = TRUE
    )))
    rec <- sf::st_point(c(ox, oy + 500))
    feats <- sf::st_sf(id = c("poly", "rec"),
                       geometry = sf::st_sfc(poly, rec, crs = 32755))
    roles <- data.frame(feature_id = c("poly", "rec"), hazard = "odour",
                        role = c("source", "receptor"),
                        emit_height = c(10, NA_real_),
                        stringsAsFactors = FALSE)
    site <- mh_site(feats, roles, epsg = 32755L)
    out  <- odour_exposure(.char_met(), site, terrain_backend = "none")
    expect_snapshot_value(round(out, 8), style = "json2", tolerance = 1e-6)
  })
})


# ---------------------------------------------------------------------------
# Helper-level value guards (localised failure messages for the vectorised
# primitives the refactor touches).
# ---------------------------------------------------------------------------

describe(".crosswind_halfwidth() vector path equals the scalar path", {

  it("agrees elementwise with a per-direction scalar evaluation", {
    ox <- 335000; oy <- 6250000; hw <- 90
    poly <- sf::st_sfc(sf::st_polygon(list(matrix(
      c(ox-hw, oy-hw, ox+hw, oy-2*hw, ox+2*hw, oy+hw, ox-hw, oy+hw, ox-hw, oy-hw),
      ncol = 2, byrow = TRUE
    ))), crs = 32755)
    dirs <- seq(0, 359, by = 1)
    vec_path <- .crosswind_halfwidth(poly, dirs)
    scl_path <- vapply(dirs, function(wd) .crosswind_halfwidth(poly, wd),
                       numeric(1))
    expect_equal(vec_path, scl_path, tolerance = 1e-12)
  })
})

describe(".pool_partition() golden output with a varying emission band", {

  it("matches the pinned fractions for vector emit_extent / pool_top / delta", {
    # All three arguments vary across the vector — the case the current tests
    # never exercise (they broadcast a scalar emit_extent).
    emit  <- c(0, 2, 5, 10, 20, 40, 5, 30, 12, 1)
    ptop  <- c(50, 50, 30, 30, 30, 30, 0.5, 0.5, 25, 100)
    delta <- c(20, 20, 10, 10, 10, 10, 0.1, 0.1, 8, 25)
    p <- .pool_partition(emit, ptop, delta)
    # Budget invariant holds elementwise.
    expect_equal(p$f_1a + p$f_1b, rep(1, length(emit)), tolerance = 1e-12)
    expect_snapshot_value(round(p$f_1b, 10), style = "json2", tolerance = 1e-8)
  })
})

describe(".cw_venting() / .cw_fumigation() golden vectors", {

  # Reuse the 24-h pathway fixture style locally.
  .vmet <- function(...) {
    n <- 24
    base <- list(
      wind_speed_10m         = rep(2, n),
      direct_radiation       = c(rep(0, 12), rep(200, 12)),
      cloud_cover            = rep(20, n),
      boundary_layer_height  = c(rep(150, 11), 150, 300, 500, 700, 900, 1100,
                                  1300, 1500, 1600, 1700, 1800, 1900, 1950),
      temperature_2m         = rep(10, n),
      relative_humidity_2m   = rep(70, n),
      precipitation          = rep(0, n),
      wind_direction_10m     = rep(180, n),
      pressure_msl           = rep(1013, n),
      soil_moisture_0_to_1cm = rep(0.1, n),
      soil_moisture_1_to_3cm = rep(0.1, n),
      wind_direction_80m     = rep(250, n),
      wind_speed_80m         = rep(3, n),
      wind_direction_120m    = rep(255, n),
      wind_speed_120m        = rep(3, n),
      wind_direction_180m    = rep(260, n),
      wind_speed_180m        = rep(3, n)
    )
    ov <- list(...); base[names(ov)] <- ov
    as.data.frame(lapply(base, rep_len, n))
  }

  it("pins .cw_venting() for a channelled terrain at a 60-degree receptor", {
    ter <- mh_terrain(drainage_bearing = 30, flow_convergence = 0.75)
    vs  <- ventilation_state(.vmet(), terrain = ter)
    cw  <- .cw_venting(bearing_to_receptor = 60, vs = vs, terrain = ter)
    expect_snapshot_value(round(cw, 10), style = "json2", tolerance = 1e-8)
  })

  it("pins .cw_fumigation() for a channelled terrain at a 90-degree receptor", {
    ter <- mh_terrain(drainage_bearing = 30, flow_convergence = 0.75)
    vs  <- ventilation_state(.vmet(), terrain = ter)
    vs$wind_dir_surface <- .vmet()$wind_direction_10m
    cw  <- .cw_fumigation(bearing_to_receptor = 90, vs = vs, terrain = ter,
                          above_pool_ht = 10)
    expect_snapshot_value(round(cw, 10), style = "json2", tolerance = 1e-8)
  })
})
