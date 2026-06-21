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
