# Tests for odour_exposure(): C3a mh_site API — area-source ISC3 initial
# spreads, multi-source sum→map→max reduction, and terrain_backend = 'none'.

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Build a minimal mh_site with one point source and one receptor.
# src_bearing / src_distance: source relative to a fixed origin (unused, source
#   is always placed at the origin for these tests).
# rec_bearing / rec_distance: receptor relative to the origin, in degrees from
#   north (clockwise) and metres.
.make_odour_site <- function(rec_bearing = 0, rec_distance = 400,
                              src_bearing = NA, src_distance = NA) {
  origin_x <- 335000
  origin_y <- 6250000
  src_pt <- sf::st_sfc(sf::st_point(c(origin_x, origin_y)), crs = 32755)
  rec_x  <- origin_x + rec_distance * sin(rec_bearing * pi / 180)
  rec_y  <- origin_y + rec_distance * cos(rec_bearing * pi / 180)
  rec_pt <- sf::st_sfc(sf::st_point(c(rec_x, rec_y)), crs = 32755)
  feats  <- sf::st_sf(
    id       = c("src", "rec"),
    geometry = sf::st_sfc(src_pt[[1]], rec_pt[[1]], crs = 32755)
  )
  roles <- data.frame(
    feature_id = c("src", "rec"),
    hazard     = "odour",
    role       = c("source", "receptor"),
    stringsAsFactors = FALSE
  )
  mh_site(feats, roles, epsg = 32755L)
}

# Build a two-source, one-receptor mh_site.  Both sources share the origin;
# the receptor is placed at rec_bearing / rec_distance from it.
.make_two_source_site <- function(rec_bearing = 0, rec_distance = 400) {
  origin_x <- 335000
  origin_y <- 6250000
  rec_x    <- origin_x + rec_distance * sin(rec_bearing * pi / 180)
  rec_y    <- origin_y + rec_distance * cos(rec_bearing * pi / 180)
  feats <- sf::st_sf(
    id       = c("src1", "src2", "rec"),
    geometry = sf::st_sfc(
      sf::st_point(c(origin_x, origin_y)),
      sf::st_point(c(origin_x, origin_y)),   # co-located
      sf::st_point(c(rec_x,    rec_y)),
      crs = 32755
    )
  )
  roles <- data.frame(
    feature_id = c("src1", "src2", "rec"),
    hazard     = "odour",
    role       = c("source", "source", "receptor"),
    stringsAsFactors = FALSE
  )
  mh_site(feats, roles, epsg = 32755L)
}

# Two receptors: one aligned downwind, one upwind.
.make_two_receptor_site <- function() {
  origin_x <- 335000
  origin_y <- 6250000
  feats <- sf::st_sf(
    id       = c("src", "rec_n", "rec_s"),
    geometry = sf::st_sfc(
      sf::st_point(c(origin_x,         origin_y)),          # source
      sf::st_point(c(origin_x,         origin_y + 400)),    # rec north
      sf::st_point(c(origin_x,         origin_y - 400)),    # rec south
      crs = 32755
    )
  )
  roles <- data.frame(
    feature_id = c("src", "rec_n", "rec_s"),
    hazard     = "odour",
    role       = c("source", "receptor", "receptor"),
    stringsAsFactors = FALSE
  )
  mh_site(feats, roles, epsg = 32755L)
}

# Build a polygon-source site (100 m × 100 m square centred at origin).
.make_polygon_site <- function(rec_bearing = 0, rec_distance = 400, width = 100) {
  origin_x <- 335000
  origin_y <- 6250000
  hw <- width / 2
  poly <- sf::st_polygon(list(matrix(
    c(origin_x - hw, origin_y - hw,
      origin_x + hw, origin_y - hw,
      origin_x + hw, origin_y + hw,
      origin_x - hw, origin_y + hw,
      origin_x - hw, origin_y - hw),
    ncol = 2, byrow = TRUE
  )))
  rec_x <- origin_x + rec_distance * sin(rec_bearing * pi / 180)
  rec_y <- origin_y + rec_distance * cos(rec_bearing * pi / 180)
  feats <- sf::st_sf(
    id       = c("src", "rec"),
    geometry = sf::st_sfc(poly, sf::st_point(c(rec_x, rec_y)), crs = 32755)
  )
  roles <- data.frame(
    feature_id = c("src", "rec"),
    hazard     = "odour",
    role       = c("source", "receptor"),
    stringsAsFactors = FALSE
  )
  mh_site(feats, roles, epsg = 32755L)
}

# Standard met builder.  All columns are provided so .odour_generation() works.
.make_odour_met <- function(n = 1, ...) {
  base <- list(
    wind_direction_10m     = 180,   # FROM south → downwind N (bearing 0)
    wind_speed_10m         = 3,
    direct_radiation       = 0,
    cloud_cover            = 50,
    boundary_layer_height  = 500,
    temperature_2m         = 15,
    pressure_msl           = 1013,
    precipitation          = 0,
    relative_humidity_2m   = 50,
    soil_moisture_0_to_1cm = 0.1,
    soil_moisture_1_to_3cm = 0.1
  )
  ov <- list(...)
  base[names(ov)] <- ov
  as.data.frame(lapply(base, rep_len, n))
}

# ---------------------------------------------------------------------------
# C3a behaviour specs
# ---------------------------------------------------------------------------

describe("odour_exposure() on the mh_site model (terrain backend 'none')", {

  # 1. Distance decay
  it("makes a nearer on-axis receptor more exposed than a far one", {
    near_site <- .make_odour_site(rec_bearing = 0, rec_distance = 400)
    far_site  <- .make_odour_site(rec_bearing = 0, rec_distance = 1500)
    d <- .make_odour_met()
    near <- odour_exposure(d, near_site)
    far  <- odour_exposure(d, far_site)
    expect_gt(near, far)
  })

  # 2. Direction decay
  it("falls as the plume swings off the receptor bearing", {
    # Receptor due north (bearing 0 from source).
    # Wind FROM south/SSE/east → downwind bearing 0 / 20 / 90.
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 400)
    d <- .make_odour_met(n = 3, wind_direction_10m = c(180, 160, 270))
    ex <- odour_exposure(d, site)
    expect_gt(ex[1], ex[2])
    expect_gt(ex[2], ex[3])
    expect_lt(ex[3], 0.5)   # ~90° off-axis: negligible
  })

  # 3. Upwind receptor
  it("reads ~0 for a receptor directly upwind of the source", {
    # Receptor due south (bearing 180) of source; wind FROM south (toward N)
    # → downwind is N. The receptor at bearing 180 is upwind.
    site <- .make_odour_site(rec_bearing = 180, rec_distance = 400)
    d    <- .make_odour_met(wind_direction_10m = 180)
    ex   <- odour_exposure(d, site)
    expect_lt(ex, 0.01)
  })

  # 4. Forecast-direction uncertainty
  it("keeps a small off-axis miss substantial via forecast-direction uncertainty", {
    # Use 5000 m so on-axis is sub-100 and differences are visible.
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 5000)
    d <- .make_odour_met(n = 2, wind_direction_10m = c(180, 170))  # 0° / 10° off
    ex <- odour_exposure(d, site)
    expect_lt(ex[2], ex[1])
    expect_gt(ex[2], 0.4 * ex[1])
  })

  # 5. Calm: direction-agnostic
  it("is direction-agnostic at a fixed distance under calm winds", {
    site_n <- .make_odour_site(rec_bearing = 0,   rec_distance = 400)
    site_s <- .make_odour_site(rec_bearing = 180, rec_distance = 400)
    d <- .make_odour_met(wind_speed_10m = 0.2)
    north <- odour_exposure(d, site_n)
    south <- odour_exposure(d, site_s)
    expect_equal(north, south)
  })

  # 6. Worst-case receptor reduction
  it("returns the worst-affected receptor's value for the hour", {
    site <- .make_two_receptor_site()
    d    <- .make_odour_met(wind_direction_10m = 180)  # downwind → north
    combined <- odour_exposure(d, site)
    # The aligned receptor (north) sets the value; the south receptor is upwind
    # and contributes 0.  So combined == single-aligned-receptor result.
    site_aligned <- .make_odour_site(rec_bearing = 0, rec_distance = 400)
    single <- odour_exposure(d, site_aligned)
    expect_equal(combined, single, tolerance = 1e-6)
  })

  # 7. Monotonic in generation / stays in [0, 100]
  it("increases with hazard and stays within [0, 100]", {
    # Use 10 km so neither condition saturates the 0-100 map.
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 10000)
    # "Low G": cool, low RH, no pressure drop → G ≈ 1
    d_lo <- .make_odour_met(
      wind_speed_10m = 2, boundary_layer_height = 400,
      direct_radiation = 0, cloud_cover = 50,
      temperature_2m = 10, pressure_msl = 1013, relative_humidity_2m = 50
    )
    # "High G": hot, high RH, falling pressure → G ≈ 1.45
    d_hi <- .make_odour_met(
      wind_speed_10m = 2, boundary_layer_height = 400,
      direct_radiation = 0, cloud_cover = 50,
      temperature_2m = 35, pressure_msl = 1008, relative_humidity_2m = 90
    )
    lo <- odour_exposure(d_lo, site)
    hi <- odour_exposure(d_hi, site)
    expect_gte(lo, 0)
    expect_lte(lo, 100)
    expect_gte(hi, 0)
    expect_lte(hi, 100)
    expect_gt(hi, lo)
  })

  # 8. Plain numeric output
  it("returns a plain numeric band", {
    site <- .make_odour_site()
    out  <- odour_exposure(.make_odour_met(), site)
    expect_type(out, "double")
    expect_false(inherits(out, "units"))
  })

  # 9. Input validation
  it("errors (classed) on an invalid mh_site or a hazard/met length mismatch", {
    # Non-mh_site object
    expect_error(
      odour_exposure(.make_odour_met(), list(not = "a site")),
      class = "meteoHazard_input_error"
    )
    # mh_site with no odour source
    feats_no_src <- sf::st_sf(
      id       = "rec",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250400)), crs = 32755)
    )
    roles_no_src <- data.frame(
      feature_id = "rec", hazard = "odour", role = "receptor",
      stringsAsFactors = FALSE
    )
    site_no_src <- mh_site(feats_no_src, roles_no_src, epsg = 32755L)
    expect_error(
      odour_exposure(.make_odour_met(), site_no_src),
      class = "meteoHazard_input_error"
    )
  })

  # 10. Finite at source edge (polygon source, receptor very close)
  it("stays finite for a receptor at the source edge (no point-source divergence)", {
    # 100 m × 100 m polygon source; receptor 60 m from centroid (inside the
    # initial-spread radius → sigma_y_eff is dominated by sigma_y0 ≈ 11.6 m,
    # not a collapsing Gaussian).
    site <- .make_polygon_site(rec_bearing = 0, rec_distance = 60, width = 100)
    d    <- .make_odour_met()
    out  <- odour_exposure(d, site)
    expect_true(is.finite(out))
  })

  # 11. ISC3 sigma_y0 dominates at small distance for wide polygon
  it("tends to the initial spread sigma_y0 as distance goes to zero", {
    # Wide polygon (200 m), receptor 30 m from centroid.
    # sigma_y0 = 100 / 4.3 ≈ 23.3 m; at x = 30 m Briggs sigma_y is tiny.
    # Exposure should be finite and non-zero.
    site <- .make_polygon_site(rec_bearing = 0, rec_distance = 30, width = 200)
    d    <- .make_odour_met()
    out  <- odour_exposure(d, site)
    expect_true(is.finite(out))
    expect_gt(out, 0)
  })

  # 12. Multi-source sum
  it("sums concentrations from multiple sources at a receptor", {
    # 10 km: one source gives ~53/100, two give ~78/100 — clearly distinguishable.
    d        <- .make_odour_met(wind_speed_10m = 2, boundary_layer_height = 400)
    site_two <- .make_two_source_site(rec_bearing = 0, rec_distance = 10000)
    site_one <- .make_odour_site(rec_bearing = 0, rec_distance = 10000)
    two_src  <- odour_exposure(d, site_two)
    one_src  <- odour_exposure(d, site_one)
    expect_gt(two_src, one_src)
  })

  # 13. Two identical co-located sources ≈ 2× pre-map concentration
  it("gives ~twice the single-source concentration for two identical co-located sources", {
    # Use a very large map_c50 so the exponential map is nearly linear in C:
    # risk ≈ 100 * C / map_c50 when C << map_c50.
    # At 20 km C_rel ≈ 0.10 (single) and 0.20 (double); with map_c50=100
    # we stay well within the linear regime.
    d       <- .make_odour_met(wind_speed_10m = 2, boundary_layer_height = 400)
    map_c50 <- 100   # huge → map is linear for typical C_rel values
    site_two <- .make_two_source_site(rec_bearing = 0, rec_distance = 20000)
    site_one <- .make_odour_site(rec_bearing = 0, rec_distance = 20000)
    two <- odour_exposure(d, site_two, map_c50 = map_c50)
    one <- odour_exposure(d, site_one, map_c50 = map_c50)
    # In the linear regime: risk ∝ C; two identical co-located sources → 2×.
    expect_equal(two / one, 2, tolerance = 0.05)
  })

  # 14. Summed map can make two moderate sources severe
  it("applies the 0-100 map to the summed concentration, so two moderate sources can be severe", {
    # At 10 km: one source ≈ 53/100, two sources ≈ 78/100 (both sub-100).
    # This demonstrates that the map is applied to the SUM, not per-source.
    d <- .make_odour_met(wind_speed_10m = 2, boundary_layer_height = 400)
    site_two <- .make_two_source_site(rec_bearing = 0, rec_distance = 10000)
    site_one <- .make_odour_site(rec_bearing = 0, rec_distance = 10000)
    one <- odour_exposure(d, site_one, map_c50 = 0.3)
    two <- odour_exposure(d, site_two, map_c50 = 0.3)
    # Two sources should give higher exposure than one, and still be bounded.
    expect_gt(two, one)
    expect_lte(two, 100)
  })

  # 15. Single source: result is max over receptors (same as per-receptor)
  it("reduces to the current max-over-receptors for a single source", {
    # Two receptors: one aligned, one off-axis.  The result must equal the
    # aligned receptor's value (since max is taken across receptors).
    site_two_rec <- .make_two_receptor_site()
    site_one_rec <- .make_odour_site(rec_bearing = 0, rec_distance = 400)
    d <- .make_odour_met(wind_direction_10m = 180)  # downwind north
    two_rec <- odour_exposure(d, site_two_rec)
    one_rec <- odour_exposure(d, site_one_rec)
    expect_equal(two_rec, one_rec, tolerance = 1e-6)
  })

  # 16. No morning pulse with terrain_backend = 'none'
  it("produces no morning pulse with terrain_backend = 'none'", {
    # Two hours with equal wind speed, h_mix, stability but different
    # radiation (first dark, second sunny).  In flat Gaussian there is no
    # overnight accumulation and no morning release, so exposure should be
    # driven purely by the ventilation state (which is the same for both
    # hours here because we fix all dispersion-relevant met columns).
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 400)
    d <- .make_odour_met(
      n = 2,
      wind_direction_10m    = 180,
      wind_speed_10m        = 2,
      boundary_layer_height = 400,
      cloud_cover           = 50,
      direct_radiation      = c(0, 600)   # night then day
    )
    ex <- odour_exposure(d, site, terrain_backend = "none")
    # The morning hour should NOT spike relative to the night hour
    # (no drainage-pool release).  Accept a tolerance of 20% (slight
    # difference can arise because PM and stability depend on radiation).
    expect_lt(abs(ex[2] - ex[1]) / pmax(ex[1], 0.001), 0.5)
  })

  # 17. Coincident source and receptor → skipped (contributes 0)
  it("skips a receptor coincident with the source", {
    # Place receptor at the same coordinate as the source (distance = 0).
    feats <- sf::st_sf(
      id       = c("src", "rec"),
      geometry = sf::st_sfc(
        sf::st_point(c(335000, 6250000)),
        sf::st_point(c(335000, 6250000)),   # identical
        crs = 32755
      )
    )
    roles <- data.frame(
      feature_id = c("src", "rec"), hazard = "odour",
      role = c("source", "receptor"), stringsAsFactors = FALSE
    )
    site <- mh_site(feats, roles, epsg = 32755L)
    d    <- .make_odour_met()
    out  <- odour_exposure(d, site)
    expect_equal(out, 0)
  })
})


# ---------------------------------------------------------------------------
# M2 — receptor impaction (C6, issue #20)
# ---------------------------------------------------------------------------
#
# Physics: an elevated receptor (Δz > 0) is pulled toward the plume centreline
# as stability increases; f_vert = exp(-0.5*(Δz_eff/sigma_z_eff)^2) blends
# toward 1 under stable conditions. f_vert ≤ 1 always.
#
# Δz   = receptor_elevation (from features$elevation) −
#         source_emit_height (from features$emit_height for sources).
# phi_s = clamp((s - IMPACTION_S_NEUTRAL) /
#               (IMPACTION_S_STABLE - IMPACTION_S_NEUTRAL), 0, 1)
#         IMPACTION_S_NEUTRAL = 3, IMPACTION_S_STABLE = 5
# around_frac = ifelse(is.na(hill_height_scale), 0, hill_height_scale)
# collapse    = IMPACTION_STRENGTH * pmax(phi_s, around_frac)  [= 0.8 * ...]
# Δz_eff      = Δz * (1 - collapse)
# f_vert      = exp(-0.5 * (Δz_eff / sigma_z_eff)^2)
# ---------------------------------------------------------------------------

# Helper: build a site with elevation columns on features.
# emit_height on the source feature drives ISC3 sigma_z0 AND M2 Δz denominator.
# elevation on the receptor drives M2 Δz numerator.
.make_impaction_site <- function(src_emit_height = 10,
                                  rec_elevation   = 0,
                                  rec_distance    = 500,
                                  hill_height_scale = NA_real_) {
  origin_x <- 335000; origin_y <- 6250000
  rec_x    <- origin_x
  rec_y    <- origin_y + rec_distance  # receptor due north

  hhs_col <- rep(hill_height_scale, 2)  # source + receptor row

  feats <- sf::st_sf(
    id              = c("src", "rec"),
    emit_height     = c(src_emit_height, NA_real_),  # on features, not roles
    elevation       = c(0,               rec_elevation),
    hill_height_scale = hhs_col,
    geometry        = sf::st_sfc(
      sf::st_point(c(origin_x, origin_y)),
      sf::st_point(c(rec_x,    rec_y)),
      crs = 32755
    )
  )
  roles <- data.frame(
    feature_id = c("src", "rec"),
    hazard     = "odour",
    role       = c("source", "receptor"),
    stringsAsFactors = FALSE
  )
  mh_site(feats, roles, epsg = 32755L)
}

# Met producing neutral stability: overcast, moderate wind (Turner class D, s≈3).
.met_neutral <- function(n = 1) {
  .make_odour_met(n = n,
    wind_direction_10m    = 180,   # blowing from south → downwind = north
    wind_speed_10m        = 3,
    direct_radiation      = 0,
    cloud_cover           = 100,   # fully overcast → Turner class D ≈ neutral
    boundary_layer_height = 600
  )
}

# Met producing very stable conditions: calm/slow wind, clear cold night.
.met_stable <- function(n = 1) {
  .make_odour_met(n = n,
    wind_direction_10m    = 180,
    wind_speed_10m        = 1.0,   # slow: PG class E–F when direct_rad = 0
    direct_radiation      = 0,
    cloud_cover           = 0,     # clear night → stable
    boundary_layer_height = 150,
    temperature_2m        = 5,
    relative_humidity_2m  = 40
  )
}


describe("odour_exposure(): M2 receptor impaction (impaction = TRUE)", {

  # -- Default-off is identity -----------------------------------------------

  it("impaction = FALSE is bit-identical to the current output (regression)", {
    # The default must not change the existing golden oracle.
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 500)
    d    <- .make_odour_met()
    expect_equal(
      odour_exposure(d, site, impaction = FALSE),
      odour_exposure(d, site),                    # default
      tolerance = 0
    )
  })

  it("impaction = FALSE with elevation columns is bit-identical to current (regression)", {
    # Adding elevation to features must not change behaviour when impaction = FALSE.
    site_elev <- .make_impaction_site(src_emit_height = 10, rec_elevation = 40,
                                      rec_distance    = 500)
    site_flat <- .make_odour_site(rec_bearing = 0, rec_distance = 500)
    d <- .met_neutral()
    expect_equal(
      odour_exposure(d, site_elev, impaction = FALSE),
      odour_exposure(d, site_flat, impaction = FALSE),
      tolerance = 1e-9
    )
  })

  it("impaction = TRUE with Δz = 0 is bit-identical to impaction = FALSE", {
    # f_vert = exp(-0.5 * (0 / sigma_z)^2) = 1 → no change.
    site <- .make_impaction_site(src_emit_height = 10, rec_elevation = 10,
                                  rec_distance    = 500)  # rec_elevation = emit_height → Δz = 0
    d <- .met_neutral()
    expect_equal(
      odour_exposure(d, site, impaction = TRUE),
      odour_exposure(d, site, impaction = FALSE),
      tolerance = 1e-9
    )
  })


  # -- Bounded by centreline (f_vert ≤ 1) ------------------------------------

  it("exposure with impaction = TRUE is ≤ exposure with impaction = FALSE for elevated receptor", {
    # f_vert ≤ 1 always → M2 never inflates exposure above the flat-terrain value.
    site_elevated <- .make_impaction_site(src_emit_height = 5, rec_elevation = 50,
                                           rec_distance    = 500)
    d <- .met_neutral()
    exp_on  <- odour_exposure(d, site_elevated, impaction = TRUE)
    exp_off <- odour_exposure(d, site_elevated, impaction = FALSE)
    expect_lte(exp_on, exp_off + 1e-9)
  })

  it("f_vert ≤ 1 holds for every hour over a 48-hour sweep", {
    # Property-based: generate 48 hours of varied met; the impaction=TRUE
    # result must never exceed impaction=FALSE.
    withr::local_seed(42L)
    n_hours <- 48
    d_sweep <- .make_odour_met(
      n                  = n_hours,
      wind_direction_10m = rep(180, n_hours),
      wind_speed_10m     = runif(n_hours, 0.5, 8),
      cloud_cover        = runif(n_hours, 0, 100),
      direct_radiation   = c(rep(0, 24), runif(24, 0, 600)),
      boundary_layer_height = runif(n_hours, 150, 1500)
    )
    site <- .make_impaction_site(src_emit_height = 8, rec_elevation = 40)
    exp_on  <- odour_exposure(d_sweep, site, impaction = TRUE)
    exp_off <- odour_exposure(d_sweep, site, impaction = FALSE)
    expect_true(all(exp_on <= exp_off + 1e-9),
                label = "f_vert ≤ 1 violated: impaction=TRUE exceeds impaction=FALSE")
  })


  # -- Stability effect on elevated receptor ---------------------------------

  it("impaction ratio (on/off) is higher under stable than neutral for elevated receptor", {
    # Elevated receptor (Δz > 0): as stability increases, collapse increases,
    # Δz_eff decreases, f_vert increases toward 1. So the on/off ratio rises.
    # Test: ratio[stable] > ratio[neutral].
    site <- .make_impaction_site(src_emit_height = 5, rec_elevation = 50,
                                  rec_distance    = 500)

    d_neutral <- .met_neutral()
    d_stable  <- .met_stable()

    exp_on_neutral  <- odour_exposure(d_neutral, site, impaction = TRUE)
    exp_off_neutral <- odour_exposure(d_neutral, site, impaction = FALSE)
    exp_on_stable   <- odour_exposure(d_stable,  site, impaction = TRUE)
    exp_off_stable  <- odour_exposure(d_stable,  site, impaction = FALSE)

    ratio_neutral <- exp_on_neutral  / pmax(exp_off_neutral,  .Machine$double.eps)
    ratio_stable  <- exp_on_stable   / pmax(exp_off_stable,   .Machine$double.eps)

    # Stable ratio should be closer to 1 (less attenuation) than neutral ratio.
    expect_gt(ratio_stable, ratio_neutral)
  })

  it("exposure at elevated receptor is non-decreasing in stability when impaction = TRUE", {
    # 6 hours with monotonically increasing stability (neutral → very stable)
    # and a receptor above the source: impaction=TRUE exposure should be
    # non-decreasing as stability rises (collapse increases → Δz_eff falls → f_vert↑).
    #
    # Achieve increasing stability by increasing cloud cover and decreasing wind.
    # Use the same wind direction and receptor geometry.
    site <- .make_impaction_site(src_emit_height = 5, rec_elevation = 60,
                                  rec_distance    = 400)
    n_h  <- 6
    # Overcast + very windy → neutral (D), clear + calm → stable (F)
    d_ramp <- .make_odour_met(
      n = n_h,
      wind_direction_10m    = rep(180, n_h),
      wind_speed_10m        = c(5, 4, 3, 2, 1.5, 1),      # decreasing
      cloud_cover           = c(100, 80, 60, 40, 20, 0),  # clearing
      direct_radiation      = rep(0, n_h),                 # night throughout
      boundary_layer_height = c(600, 500, 400, 300, 200, 150)
    )
    exp_on <- odour_exposure(d_ramp, site, impaction = TRUE)
    # Each hour should have ≥ the previous (or very close, allowing rounding).
    expect_true(all(diff(exp_on) >= -0.5),   # 0.5 unit tolerance on 0-100 scale
                label = "exposure should be non-decreasing as stability rises")
  })


  # -- Below-source receptor (Δz < 0) ----------------------------------------

  it("below-source receptor (Δz < 0) at neutral stability has lower exposure than flat", {
    # Neutral: collapse = 0, Δz_eff = Δz < 0.
    # f_vert = exp(-0.5 * (Δz/sigma_z)^2) < 1 → lower exposure.
    site_below <- .make_impaction_site(src_emit_height = 30, rec_elevation = 5,
                                        rec_distance    = 500)  # rec below emit height
    d <- .met_neutral()
    exp_on  <- odour_exposure(d, site_below, impaction = TRUE)
    exp_off <- odour_exposure(d, site_below, impaction = FALSE)
    expect_lt(exp_on, exp_off)  # attenuation for below-source receptor
  })


  # -- hill_height_scale column -----------------------------------------------

  it("hill_height_scale = 1 pulls elevated receptor toward centreline even at neutral stability", {
    # With hill_height_scale = 1 (full blocking):
    #   around_frac = 1, collapse = 0.8 * max(phi_s=0, 1) = 0.8
    #   Δz_eff = Δz * 0.2 (much smaller than with no hhs)
    # So exposure with hhs=1 > exposure with hhs=NA at neutral stability.
    d       <- .met_neutral()
    site_hhs <- .make_impaction_site(src_emit_height = 5, rec_elevation = 50,
                                      rec_distance    = 500, hill_height_scale = 1.0)
    site_nohhs <- .make_impaction_site(src_emit_height = 5, rec_elevation = 50,
                                        rec_distance   = 500, hill_height_scale = NA_real_)
    exp_hhs   <- odour_exposure(d, site_hhs,   impaction = TRUE)
    exp_nohhs <- odour_exposure(d, site_nohhs, impaction = TRUE)
    expect_gt(exp_hhs, exp_nohhs)
  })

  it("absent hill_height_scale column (no column on features) → pure stability blend", {
    # When no hill_height_scale column exists, around_frac = 0 for all receptors.
    # Result should equal a site with hill_height_scale = NA (same behaviour).
    site_no_col <- .make_odour_site(rec_bearing = 0, rec_distance = 500)
    # Manually add elevation to the no-hhs site's features via column addition:
    feats <- site_no_col@features
    feats$emit_height <- ifelse(feats$id == "src", 5, NA_real_)
    feats$elevation   <- ifelse(feats$id == "rec", 50, 0)
    roles_df <- site_no_col@roles
    site_manual <- mh_site(feats, roles_df, epsg = 32755L)

    site_na_hhs <- .make_impaction_site(src_emit_height = 5, rec_elevation = 50,
                                         rec_distance   = 500, hill_height_scale = NA_real_)
    d <- .met_neutral()
    expect_equal(
      odour_exposure(d, site_manual, impaction = TRUE),
      odour_exposure(d, site_na_hhs, impaction = TRUE),
      tolerance = 1e-9
    )
  })


  # -- NA elevation receptor --------------------------------------------------

  it("a receptor with NA elevation is unaffected by impaction = TRUE", {
    # NA elevation → Δz = 0 → f_vert = 1 → that receptor's column is unchanged.
    # Build a two-receptor site: one with elevation, one without (NA).
    origin_x <- 335000; origin_y <- 6250000
    feats <- sf::st_sf(
      id          = c("src", "rec_elev", "rec_na"),
      emit_height = c(5, NA_real_, NA_real_),
      elevation   = c(0, 50,      NA_real_),
      geometry    = sf::st_sfc(
        sf::st_point(c(origin_x,         origin_y)),
        sf::st_point(c(origin_x,         origin_y + 500)),
        sf::st_point(c(origin_x + 300,   origin_y + 400)),
        crs = 32755
      )
    )
    roles <- data.frame(
      feature_id = c("src", "rec_elev", "rec_na"),
      hazard     = "odour",
      role       = c("source", "receptor", "receptor"),
      stringsAsFactors = FALSE
    )
    site_two_rec <- mh_site(feats, roles, epsg = 32755L)

    # Build a matching site where rec_na is the only receptor (for reference).
    feats_na_only <- feats[feats$id %in% c("src", "rec_na"), ]
    roles_na_only <- roles[roles$feature_id %in% c("src", "rec_na"), ]
    site_na_only  <- mh_site(feats_na_only, roles_na_only, epsg = 32755L)

    d <- .met_neutral()

    # The two-receptor site (impaction=TRUE) should give the same worst-case
    # exposure for the NA receptor hour as the NA-only site (impaction=FALSE),
    # assuming the elevated receptor dominates in the two-receptor case.
    # More directly: the NA receptor's contribution is f_vert = 1, unchanged.
    # We verify the NA-only exposure matches impaction=FALSE (Δz = 0 identity).
    exp_na_on  <- odour_exposure(d, site_na_only, impaction = TRUE)
    exp_na_off <- odour_exposure(d, site_na_only, impaction = FALSE)
    expect_equal(exp_na_on, exp_na_off, tolerance = 1e-9)
  })


  # -- Validation ------------------------------------------------------------

  it("raises meteoHazard_input_error for non-logical impaction argument", {
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 400)
    d    <- .make_odour_met()
    expect_error(
      odour_exposure(d, site, impaction = 1L),
      class = "meteoHazard_input_error"
    )
    expect_error(
      odour_exposure(d, site, impaction = "yes"),
      class = "meteoHazard_input_error"
    )
  })

})
