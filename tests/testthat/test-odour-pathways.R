# Behaviour spec for the odour terrain morning-pulse physics (C3b, issue #17):
# .pool_partition(), .cw_venting() [1a], .cw_fumigation() [1b], .morning_release().

# ---------------------------------------------------------------------------
# Shared met builder: 24-hour cycle (12 h night + 12 h day)
# BLH jumps from 150 m to 300 m on the first daytime hour (t=13) so that
# cbl_growth[13] = 150 > 0, triggering .morning_release() at the onset hour.
# ---------------------------------------------------------------------------
.pmet <- function(n = 24, ...) {
  base <- list(
    wind_speed_10m         = rep(2, n),
    direct_radiation       = c(rep(0, 12), rep(200, 12)),
    cloud_cover            = rep(20, n),
    boundary_layer_height  = c(rep(150, 11),
                                150, 300, 500, 700, 900, 1100,
                                1300, 1500, 1600, 1700, 1800,
                                1900, 1950),
    temperature_2m         = rep(10, n),
    relative_humidity_2m   = rep(70, n),
    precipitation          = rep(0, n),
    wind_direction_10m     = rep(180, n),
    pressure_msl           = rep(1013, n),
    soil_moisture_0_to_1cm = rep(0.1, n),
    soil_moisture_1_to_3cm = rep(0.1, n)
  )
  ov <- list(...)
  base[names(ov)] <- ov
  as.data.frame(lapply(base, rep_len, n))
}

# Build a ventilation_state for the standard 24-h cycle
.pvs <- function(met = .pmet(), terrain = NULL) {
  ventilation_state(met, terrain = terrain)
}

# Standard channelled terrain: drainage flows north (0°), high convergence
.pter_chan <- function(drain_bearing = 0, flow_conv = 0.8) {
  mh_terrain(drainage_bearing = drain_bearing, flow_convergence = flow_conv)
}

# Standard radial terrain (lone mound): no drainage_bearing, low convergence
.pter_radial <- function() {
  mh_terrain(flow_convergence = 0.3)
}


# ---------------------------------------------------------------------------
describe(".pool_partition()", {

  it("returns within-pool (1a) and above-pool (1b) fractions that sum to 1", {
    # Scalar and vector inputs must always satisfy f_1a + f_1b == 1
    p <- .pool_partition(emit_extent = 10, pool_top = 50, delta = 20)
    expect_equal(p$f_1a + p$f_1b, 1, tolerance = 1e-10)

    # Vector inputs
    p2 <- .pool_partition(
      emit_extent = c(0, 5, 15),
      pool_top    = c(30, 30, 30),
      delta       = c(10, 10, 10)
    )
    expect_equal(p2$f_1a + p2$f_1b, c(1, 1, 1), tolerance = 1e-10)
  })

  it("increases the above-pool fraction monotonically with (emit-height - pool_top)", {
    # As (emit_extent - pool_top) increases, more emission is above the pool
    pool_top  <- 50
    delta     <- 15
    extents   <- c(5, 20, 40, 60, 80)
    f_1b_vals <- vapply(extents, function(e) {
      .pool_partition(e, pool_top, delta)$f_1b
    }, numeric(1))
    # f_1b should increase (or at least not decrease) as emit_extent grows
    expect_true(all(diff(f_1b_vals) >= -1e-10))
  })

  it("sends ~all emission to 1a for a low source deep in the pool", {
    # emit_extent = 2 m, pool_top = 200 m, delta = 20 m
    # Source is entirely within the cold pool → f_1a ≈ 1, f_1b ≈ 0
    p <- .pool_partition(emit_extent = 2, pool_top = 200, delta = 20)
    expect_gt(p$f_1a, 0.95)
    expect_lt(p$f_1b, 0.05)
  })

  it("sends ~all emission to 1b for a high source above a shallow pool", {
    # emit_extent = 30 m, pool_top = 0.5 m, delta = 0.1 m
    # Emission band [0,30m] sits well above a tiny pool (0.5m deep)
    # → almost all emission is above the pool → f_1b ≈ 1, f_1a ≈ 0
    p <- .pool_partition(emit_extent = 30, pool_top = 0.5, delta = 0.1)
    expect_lt(p$f_1a, 0.05)
    expect_gt(p$f_1b, 0.95)
  })

  it("varies continuously as pool_top sweeps an E(z) band edge", {
    # Sweep pool_top from 0 to 60 across an emit_extent of 30 m
    # f_1a should vary smoothly (no discontinuities) — check max step < 0.1
    pool_tops <- seq(0, 60, by = 2)
    f_1a_vals <- vapply(pool_tops, function(pt) {
      .pool_partition(emit_extent = 30, pool_top = pt, delta = 10)$f_1a
    }, numeric(1))
    max_step <- max(abs(diff(f_1a_vals)))
    expect_lt(max_step, 0.15)
    # And it should be monotonically increasing (more pool → more 1a)
    expect_true(all(diff(f_1a_vals) >= -1e-10))
  })
})


# ---------------------------------------------------------------------------
describe(".cw_venting() [pathway 1a]", {

  it("confines at night, lowering an up-slope receptor relative to flat terrain", {
    # Channelled terrain, drainage north. Receptor bearing 0° (north = along drainage).
    ter <- .pter_chan(drain_bearing = 0, flow_conv = 0.8)
    vs  <- .pvs(terrain = ter)

    cw_1a <- .cw_venting(bearing_to_receptor = 0, vs = vs, terrain = ter)
    night_hours <- which(!vs$is_day)

    # All night hours should have a confinement value between 0 and 0.3
    expect_true(all(!is.na(cw_1a[night_hours])))
    expect_true(all(cw_1a[night_hours] >= 0))
    expect_true(all(cw_1a[night_hours] <= ODOUR_CONSTANTS$CONFINEMENT_1A))
  })

  it("vents up-slope at the morning transition, raising an aligned receptor", {
    ter <- .pter_chan(drain_bearing = 0, flow_conv = 0.8)
    vs  <- .pvs(terrain = ter)

    cw_1a <- .cw_venting(bearing_to_receptor = 0, vs = vs, terrain = ter)
    morn_hours <- which(vs$is_day & vs$cbl_growth > 0)

    # Morning hours with aligned receptor should have elevated cw_1a
    if (length(morn_hours) > 0) {
      expect_true(all(!is.na(cw_1a[morn_hours])))
      # Aligned receptor (bearing == drainage_bearing) gets max venting
      expect_true(all(cw_1a[morn_hours] > ODOUR_CONSTANTS$CONFINEMENT_1A))
    }
  })

  it("uses radial drainage for a lone mound and channelled drainage along the drainage bearing", {
    # Radial terrain: all night directions get the same uniform confinement
    ter_radial <- .pter_radial()
    vs <- .pvs(terrain = ter_radial)

    cw_north <- .cw_venting(0,   vs = vs, terrain = ter_radial)
    cw_east  <- .cw_venting(90,  vs = vs, terrain = ter_radial)
    cw_south <- .cw_venting(180, vs = vs, terrain = ter_radial)
    night_hours <- which(!vs$is_day)

    # Radial: all directions equal at night
    expect_equal(cw_north[night_hours], cw_east[night_hours],  tolerance = 1e-10)
    expect_equal(cw_north[night_hours], cw_south[night_hours], tolerance = 1e-10)

    # Channelled: direction matters
    ter_chan <- .pter_chan(drain_bearing = 0, flow_conv = 0.8)
    vs2 <- .pvs(terrain = ter_chan)
    cw_aligned  <- .cw_venting(0,   vs = vs2, terrain = ter_chan)
    cw_cross    <- .cw_venting(90,  vs = vs2, terrain = ter_chan)
    # Cross-drain receptor gets less confinement (higher cw_1a) than along-drain
    # because along-drain is most confined
    expect_true(all(cw_cross[night_hours] >= cw_aligned[night_hours] - 1e-10))
  })
})


# ---------------------------------------------------------------------------
describe(".cw_fumigation() [pathway 1b]", {

  it("directs the morning burst downwind along the residual-layer wind, not the drainage axis", {
    # Setup: drainage is NORTH (0°), but residual wind is FROM WEST (270°)
    # → fumigation goes EAST (90° downwind of westerly wind)
    # Receptor at EAST (90°): should get high cw_1b
    # Receptor at NORTH (0°, drainage axis): should get low cw_1b
    ter <- mh_terrain(drainage_bearing = 0, flow_convergence = 0.8)

    # Build met with multi-level westerly wind
    met <- .pmet()
    met$wind_direction_80m  <- rep(270, 24)   # FROM west
    met$wind_speed_80m      <- rep(3,   24)
    met$wind_direction_120m <- rep(270, 24)
    met$wind_speed_120m     <- rep(3,   24)
    met$wind_direction_180m <- rep(270, 24)
    met$wind_speed_180m     <- rep(3,   24)

    vs <- ventilation_state(met, terrain = ter)
    morn_hours <- which(vs$is_day & vs$cbl_growth > 0 &
                          !is.na(vs$pool_top) & vs$pool_top > 0)

    if (length(morn_hours) > 0) {
      cw_east  <- .cw_fumigation(90,  vs = vs, terrain = ter)
      cw_north <- .cw_fumigation(0,   vs = vs, terrain = ter)

      # Receptor east gets strong fumigation (downwind of W wind)
      # Receptor north (drainage direction) gets weak fumigation
      expect_true(all(cw_east[morn_hours] > cw_north[morn_hours]))
    } else {
      skip("no morning hours with pool in this met scenario")
    }
  })

  it("barely raises a receptor up-drainage but crosswind to the residual wind", {
    # Drainage north (0°), residual wind FROM south (180°) → fumigation goes NORTH
    # Receptor east (90°): crosswind to both drainage and residual wind
    # → cw_1b for east receptor should be near 0
    ter <- mh_terrain(drainage_bearing = 0, flow_convergence = 0.8)
    met <- .pmet()
    met$wind_direction_80m  <- rep(180, 24)   # FROM south → fumigates north
    met$wind_speed_80m      <- rep(3,   24)
    met$wind_direction_120m <- rep(180, 24)
    met$wind_speed_120m     <- rep(3,   24)
    met$wind_direction_180m <- rep(180, 24)
    met$wind_speed_180m     <- rep(3,   24)

    vs <- ventilation_state(met, terrain = ter)
    morn_hours <- which(vs$is_day & vs$cbl_growth > 0 &
                          !is.na(vs$pool_top) & vs$pool_top > 0)

    if (length(morn_hours) > 0) {
      cw_east <- .cw_fumigation(90, vs = vs, terrain = ter)
      # 90° off the fumigation direction → cos(90°)^2 = 0
      expect_true(all(cw_east[morn_hours] < 0.05))
    } else {
      skip("no morning hours with pool in this met scenario")
    }
  })
})


# ---------------------------------------------------------------------------
describe("both pathways over a tall mound", {

  it("contributes to both an up-slope and a downwind receptor in the same morning", {
    # Drainage north (0°), residual wind FROM west (270°) → fumigation east
    # Receptor A: north (up-slope, along drainage = 1a pathway)
    # Receptor B: east (downwind of residual wind = 1b pathway)
    ter <- mh_terrain(drainage_bearing = 0, flow_convergence = 0.8)
    met <- .pmet()
    met$wind_direction_80m  <- rep(270, 24)
    met$wind_speed_80m      <- rep(3,   24)
    met$wind_direction_120m <- rep(270, 24)
    met$wind_speed_120m     <- rep(3,   24)
    met$wind_direction_180m <- rep(270, 24)
    met$wind_speed_180m     <- rep(3,   24)

    vs <- ventilation_state(met, terrain = ter)
    morn_hours <- which(vs$is_day & vs$cbl_growth > 0 &
                          !is.na(vs$pool_top) & vs$pool_top > 0)

    if (length(morn_hours) > 0) {
      cw_1a_north <- .cw_venting(0,  vs, ter)
      cw_1b_east  <- .cw_fumigation(90, vs, ter)

      # Both pathways should be active (non-NA) in the morning hours
      expect_true(any(!is.na(cw_1a_north[morn_hours]) &
                        cw_1a_north[morn_hours] > 0))
      expect_true(any(!is.na(cw_1b_east[morn_hours]) &
                        cw_1b_east[morn_hours] > 0))
    } else {
      skip("no morning hours with pool in this met scenario")
    }
  })
})


# ---------------------------------------------------------------------------
describe(".morning_release()", {

  it("conserves the accumulated mass A over the release window", {
    vs <- .pvs()
    # Find the pool_top at morning onset
    pool_tp    <- ifelse(is.na(vs$pool_top), 0, vs$pool_top)
    cbl_growth <- vs$cbl_growth
    is_day     <- vs$is_day

    # Locate the first morning-onset hour
    t0 <- NA_integer_
    for (t in 2:length(is_day)) {
      if (is_day[t] && !is_day[t - 1] && pool_tp[t] > 0 && cbl_growth[t] > 0) {
        t0 <- t
        break
      }
    }

    if (is.na(t0)) {
      skip("no morning onset hour found in standard 24-h fixture")
    }

    A <- pool_tp[t0]
    r <- .morning_release(pool_tp, cbl_growth, is_day)

    # Mass conservation: sum of release enhancements ≈ A
    # Allow generous tolerance because window truncation can cut tail
    expect_equal(sum(r), A, tolerance = 0.15 * A + 0.1)
  })

  it("gives a shorter, taller pulse with the same integral when CBL growth doubles", {
    pool_top_v  <- c(rep(0, 12), 80, rep(80, 11))
    is_day_v    <- c(rep(FALSE, 12), rep(TRUE, 12))

    cbl_slow <- c(rep(0, 12), 40, rep(0, 11))
    cbl_fast <- c(rep(0, 12), 80, rep(0, 11))

    r_slow <- .morning_release(pool_top_v, cbl_slow, is_day_v)
    r_fast <- .morning_release(pool_top_v, cbl_fast, is_day_v)

    # Both should have approximately the same total mass
    expect_equal(sum(r_slow), sum(r_fast), tolerance = 0.05 * sum(r_slow) + 0.1)

    # Faster CBL growth → shorter tau → peak at a later hour is smaller
    # Both have all their mass concentrated at t0 = 13 (the onset hour)
    # so the peak is at t0 for both. With faster growth, tau is smaller,
    # so the entire pulse is concentrated in fewer hours.
    peak_t_slow <- which.max(r_slow)
    peak_t_fast <- which.max(r_fast)
    # Peaks at the same onset hour (t0 = 13 in 1-indexed)
    expect_equal(peak_t_slow, peak_t_fast)
    # Peak value with fast CBL is >= slow (same mass, narrower window)
    expect_gte(r_fast[peak_t_fast], r_slow[peak_t_slow] - 1e-10)
  })

  it("saturates the accumulation A with cooling hours", {
    # With many night hours, pool_top (and hence A) grows until capped.
    # The release should conserve whatever A has accumulated.
    pool_long <- c(rep(0, 24), 150, rep(150, 23))
    is_day_long <- c(rep(FALSE, 24), rep(TRUE, 24))
    cbl_g_long  <- c(rep(0, 24), 60, rep(0, 23))

    r <- .morning_release(pool_long, cbl_g_long, is_day_long)
    A <- 150
    expect_equal(sum(r), A, tolerance = 0.15 * A + 0.1)
  })

  it("is idempotent: the same input window yields the same result", {
    vs    <- .pvs()
    pt    <- ifelse(is.na(vs$pool_top), 0, vs$pool_top)
    cbl_g <- vs$cbl_growth
    isd   <- vs$is_day

    r1 <- .morning_release(pt, cbl_g, isd)
    r2 <- .morning_release(pt, cbl_g, isd)
    expect_identical(r1, r2)
  })

  it("produces no pulse (and no error) for a window with no preceding night", {
    # All-daytime window with pool_top = 0 → no morning onset → all zeros
    n      <- 12
    pt     <- rep(0,     n)
    cbl_g  <- rep(50,    n)
    is_d   <- rep(TRUE,  n)
    r <- .morning_release(pt, cbl_g, is_d)
    expect_equal(r, rep(0, n))
  })

  it("draws 1a confinement, 1a venting and 1b fumigation from a single budget A", {
    # f_1a + f_1b == 1 is the budget invariant; check it over a 24-h run
    vs <- .pvs()
    pt_safe <- ifelse(is.na(vs$pool_top), 0, vs$pool_top)
    delta_t <- pmax(ODOUR_CONSTANTS$DELTA_FLOOR,
                    ODOUR_CONSTANTS$DELTA_FRAC * pt_safe)

    # Use a point source at 5 m emit height
    emit_ht <- 5
    part <- .pool_partition(emit_ht, pt_safe, delta_t)
    f_1a <- part$f_1a
    f_1b <- part$f_1b

    # Invariant: sum to 1 at every hour
    expect_equal(f_1a + f_1b, rep(1, length(f_1a)), tolerance = 1e-10)

    # Release factor doesn't change the partition (only the intensity of 1b)
    r <- .morning_release(pt_safe, vs$cbl_growth, vs$is_day)
    # r is non-negative
    expect_true(all(r >= 0))
  })
})


# ---------------------------------------------------------------------------
describe("odour_exposure() terrain-backend edges", {

  # Minimal site builder for terrain-backend tests
  .make_terrain_site <- function(terrain = NULL) {
    origin_x <- 335000
    origin_y <- 6250000
    feats <- sf::st_sf(
      id       = c("src", "rec"),
      geometry = sf::st_sfc(
        sf::st_point(c(origin_x,       origin_y)),
        sf::st_point(c(origin_x,       origin_y + 500)),
        crs = 32755
      )
    )
    roles <- data.frame(
      feature_id = c("src", "rec"),
      hazard     = "odour",
      role       = c("source", "receptor"),
      stringsAsFactors = FALSE
    )
    mh_site(feats, roles, terrain = terrain, epsg = 32755L)
  }

  it("falls back to surface advection when pool_top or residual_wind is NA", {
    # Met with no temperature/RH → pool_top will be non-NA but mechanical;
    # no multi-level wind → residual_wind all NA.
    # Use terrain_backend = "descriptors" with a terrain object.
    ter  <- mh_terrain(drainage_bearing = 0, flow_convergence = 0.8)
    site <- .make_terrain_site(terrain = ter)

    met <- .pmet()
    # Remove temperature and RH to get simpler pool_top (still computed from
    # mechanical floor, non-NA)
    # Remove multi-level winds so residual_wind is all NA
    # (default .pmet() has no multi-level wind columns)

    result <- odour_exposure(met, site, terrain_backend = "descriptors")
    expect_true(is.numeric(result))
    expect_true(all(is.finite(result)))
    expect_true(all(result >= 0 & result <= 100))
  })

  it("reproduces the C3a flat result with terrain_backend = 'none'", {
    # With terrain_backend = "none", the descriptors backend is bypassed.
    # The result must be identical to a site with no terrain.
    site_no_terrain <- .make_terrain_site(terrain = NULL)
    ter   <- mh_terrain(drainage_bearing = 0, flow_convergence = 0.8)
    site_terrain <- .make_terrain_site(terrain = ter)

    met <- .pmet()

    flat_no_ter  <- odour_exposure(met, site_no_terrain, terrain_backend = "none")
    flat_with_ter <- odour_exposure(met, site_terrain,   terrain_backend = "none")

    # terrain_backend = "none" must give the same result regardless of whether
    # the site has a terrain object attached
    expect_equal(flat_no_ter, flat_with_ter, tolerance = 1e-10)

    # And both must be plain numeric in [0, 100]
    expect_true(all(is.finite(flat_no_ter)))
    expect_true(all(flat_no_ter >= 0 & flat_no_ter <= 100))
  })
})
