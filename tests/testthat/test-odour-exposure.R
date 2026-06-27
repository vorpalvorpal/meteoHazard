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
    # Without terrain_backend, there is no drainage-pool release.
    # The morning (day) hour may be LOWER than the night hour (wider plume,
    # lower PM at unstable class A vs stable class E), but must not EXCEED
    # the night by more than 20% (no artificial upward spike from pool release).
    expect_lte(ex[2], ex[1] * 1.2)
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
# C9 — wind dilution (1/u_eff advective term) and M3 end-to-end
# ---------------------------------------------------------------------------

# Site with terrain for M3 shelter tests (no drainage → M1 never fires).
.make_sheltered_site <- function(rec_bearing = 0, rec_distance = 5000,
                                  shelter_index = 50) {
  ox  <- 335000; oy <- 6250000
  rec_x <- ox + rec_distance * sin(rec_bearing * pi / 180)
  rec_y <- oy + rec_distance * cos(rec_bearing * pi / 180)
  feats <- sf::st_sf(
    id       = c("src", "rec"),
    geometry = sf::st_sfc(sf::st_point(c(ox, oy)),
                           sf::st_point(c(rec_x, rec_y)), crs = 32755)
  )
  roles <- data.frame(feature_id = c("src", "rec"), hazard = "odour",
                       role = c("source", "receptor"), stringsAsFactors = FALSE)
  mh_site(feats, roles,
          terrain = mh_terrain(shelter_index = shelter_index),
          epsg = 32755L)
}

describe("odour_exposure(): wind dilution — 1/(u_eff) advective term (C9)", {

  it("exposure is strictly lower at higher wind speed (stability fixed via shear ratio)", {
    # u80/u10 = 2 in both rows → same alpha = log(2)/log(8) ≈ 0.33 → same PG class.
    # Only u_eff differs: 2 vs 4 m/s → c_rel ∝ 1/u_eff → e2 > e4.
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 5000)
    met2 <- .make_odour_met(wind_speed_10m = 2, wind_speed_80m = 4,
                             wind_direction_10m = 180, boundary_layer_height = 500)
    met4 <- .make_odour_met(wind_speed_10m = 4, wind_speed_80m = 8,
                             wind_direction_10m = 180, boundary_layer_height = 500)
    e2 <- odour_exposure(met2, site, stability = "shear")
    e4 <- odour_exposure(met4, site, stability = "shear")
    expect_gt(e2, e4)
  })

  it("exposure ratio ≈ u_eff ratio (2:1) at large distance in the linear map regime", {
    # At 5 km, c_rel << map_c50 → linear regime → ratio ≈ u_eff4/u_eff2 = 4/2 = 2.
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 5000)
    met2 <- .make_odour_met(wind_speed_10m = 2, wind_speed_80m = 4,
                             wind_direction_10m = 180, boundary_layer_height = 500)
    met4 <- .make_odour_met(wind_speed_10m = 4, wind_speed_80m = 8,
                             wind_direction_10m = 180, boundary_layer_height = 500)
    e2 <- odour_exposure(met2, site, stability = "shear")
    e4 <- odour_exposure(met4, site, stability = "shear")
    expect_equal(e2 / e4, 2, tolerance = 0.15)
  })

  it("odour_hazard and odour_risk respond with the same ratio to a pure u_eff change", {
    # Both layers share .odour_hazard_raw() → same ventilation core → same ratio.
    # Use stability = "shear" with constant u80/u10 = 2 so stability class is
    # identical for both rows — only u_eff changes.  Large distance stays in
    # the linear map regime.
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 5000)
    d2 <- .make_odour_met(wind_speed_10m = 2, wind_speed_80m = 4,
                           wind_direction_10m = 180, boundary_layer_height = 500)
    d6 <- .make_odour_met(wind_speed_10m = 6, wind_speed_80m = 12,
                           wind_direction_10m = 180, boundary_layer_height = 500)
    hz_ratio  <- odour_hazard(d2, stability = "shear") / odour_hazard(d6, stability = "shear")
    exp_ratio <- odour_risk(d2, site, stability = "shear") /
                 odour_risk(d6, site, stability = "shear")
    expect_equal(exp_ratio, hz_ratio, tolerance = 0.10)
  })
})

describe("odour_exposure(): M3 valley shelter raises exposure end-to-end (C9)", {

  it("shelter = TRUE raises exposure on a sheltered, low-wind, non-drainage night", {
    # shelter_index=50 (max enclosure), no flow_convergence (M1 inactive).
    # u10=1.5 (=SHELTER_U_FULL) → u_eff_on=0.5, u_eff_off=1.5.
    # c_rel ∝ 1/u_eff → r_on > r_off.
    site <- .make_sheltered_site(rec_bearing = 0, rec_distance = 5000)
    d    <- .make_odour_met(wind_speed_10m = 1.5, wind_direction_10m = 180,
                             direct_radiation = 0, cloud_cover = 5,
                             boundary_layer_height = 150)
    r_off <- odour_exposure(d, site, shelter = FALSE)
    r_on  <- odour_exposure(d, site, shelter = TRUE)
    expect_gt(r_on, r_off)
  })

  it("shelter exposure ratio ≈ u_eff ratio (3:1) in the linear map regime", {
    # u_eff_off = 1.5, u_eff_on = max(1.5*(1-0.7), 0.5) = 0.5 → ratio = 1.5/0.5 = 3.
    site <- .make_sheltered_site(rec_bearing = 0, rec_distance = 5000)
    d    <- .make_odour_met(wind_speed_10m = 1.5, wind_direction_10m = 180,
                             direct_radiation = 0, cloud_cover = 5,
                             boundary_layer_height = 150)
    r_off <- odour_exposure(d, site, shelter = FALSE)
    r_on  <- odour_exposure(d, site, shelter = TRUE)
    expect_equal(r_on / r_off, 3, tolerance = 0.15)
  })
})


# ---------------------------------------------------------------------------
# C8 — Upslope rim-venting: odour_exposure() integration specs
# ---------------------------------------------------------------------------
# Written BEFORE implementation (TDD).  These fail until odour_exposure() gains
# a `rim_venting` parameter and the gate is wired through the 1a venting term.
# ---------------------------------------------------------------------------

# Site builder: receptor with explicit rel_elevation and downslope aspect.
.make_rim_site <- function(rec_bearing     = 0,
                            rec_distance    = 1000,
                            rel_elevation   = 80,
                            aspect          = NA_real_,
                            drainage_bearing = 0,
                            flow_convergence = 0.8) {
  ox <- 335000; oy <- 6250000
  rec_x <- ox + rec_distance * sin(rec_bearing * pi / 180)
  rec_y <- oy + rec_distance * cos(rec_bearing * pi / 180)
  feats <- sf::st_sf(
    id            = c("src", "rec"),
    rel_elevation = c(NA_real_, rel_elevation),
    aspect        = c(NA_real_, aspect),
    geometry      = sf::st_sfc(sf::st_point(c(ox, oy)),
                                sf::st_point(c(rec_x, rec_y)), crs = 32755)
  )
  roles <- data.frame(
    feature_id = c("src", "rec"), hazard = "odour",
    role = c("source", "receptor"), stringsAsFactors = FALSE
  )
  ter <- mh_terrain(
    drainage_bearing = drainage_bearing,
    flow_convergence = flow_convergence,
    valley_depth     = 80,
    taf              = 1.5
  )
  mh_site(feats, roles, terrain = ter, epsg = 32755L)
}

# Met: n_night cold clear nights building a deep pool, then n_day morning hours.
.make_rim_met <- function(n_night              = 12,
                           n_day               = 3,
                           wind_direction_night = 0,
                           wind_direction_day   = 180) {
  n   <- n_night + n_day
  blh <- c(rep(150, n_night), 400, 800, 1500)[seq_len(n)]
  data.frame(
    wind_direction_10m     = c(rep(wind_direction_night, n_night),
                                rep(wind_direction_day,   n_day)),
    wind_speed_10m         = rep(1.5, n),
    direct_radiation       = c(rep(0,   n_night), rep(400, n_day)),
    cloud_cover            = c(rep(0,   n_night), rep(20,  n_day)),
    boundary_layer_height  = blh,
    temperature_2m         = rep(5,    n),
    pressure_msl           = rep(1013, n),
    precipitation          = rep(0,    n),
    relative_humidity_2m   = rep(60,   n),
    soil_moisture_0_to_1cm = rep(0.1,  n),
    soil_moisture_1_to_3cm = rep(0.1,  n)
  )
}


describe("odour_exposure(): rim_venting default-off identity (C8)", {

  it("rim_venting = FALSE reproduces the pre-C8 output bit-exactly", {
    # The new parameter must be a strict no-op when FALSE.
    site <- .make_odour_site(rec_bearing = 0, rec_distance = 400)
    d    <- .make_odour_met()
    baseline <- odour_exposure(d, site)
    flagged  <- odour_exposure(d, site, rim_venting = FALSE)
    expect_equal(flagged, baseline, tolerance = .Machine$double.eps)
  })

  it("rim_venting = TRUE with no elevation/rel_elevation column is a no-op", {
    # z_j defaults to 0 → reach = 1 everywhere → bit-identical to FALSE.
    # MUST use terrain_backend = "descriptors" with actual terrain so the descriptors
    # path fires and the z_j=0 fallback is exercised (not just bypassed by terrain=NULL).
    ter <- mh_terrain(drainage_bearing = 0, flow_convergence = 0.8,
                       valley_depth = 60, taf = 1.5)
    ox  <- 335000; oy <- 6250000
    # No rel_elevation or elevation column on the receptor.
    feats_no_elev <- sf::st_sf(
      id       = c("src", "rec"),
      geometry = sf::st_sfc(sf::st_point(c(ox, oy)),
                             sf::st_point(c(ox, oy + 1000L)), crs = 32755)
    )
    roles <- data.frame(feature_id = c("src", "rec"), hazard = "odour",
                         role = c("source", "receptor"), stringsAsFactors = FALSE)
    site <- mh_site(feats_no_elev, roles, terrain = ter, epsg = 32755L)
    d    <- .make_rim_met()
    off  <- odour_exposure(d, site, terrain_backend = "descriptors", rim_venting = FALSE)
    on   <- odour_exposure(d, site, terrain_backend = "descriptors", rim_venting = TRUE)
    expect_equal(off, on, tolerance = .Machine$double.eps)
  })
})


describe("odour_exposure(): rim-venting behaviour (C8)", {

  it("aligned-slope rim receptor is more exposed on morning hours after deep pool", {
    # Receptor at bearing 60° (ENE) from source; drainage_bearing = 0° (north).
    # Current alignment from drainage = cos(60°) = 0.5.
    # With C8: aspect_j = 240° (downslope toward source at bearing 60+180 = 240°);
    # brng_rec→src = 240°; align_j = cos(0°) = 1 > 0.5.
    # Deep pool (12 cold nights): reach → 1 at morning.
    # ⇒ morning-hour exposure with rim_venting=TRUE > rim_venting=FALSE.
    site_on  <- .make_rim_site(rec_bearing = 60, rec_distance = 2000,
                                rel_elevation = 60, aspect = 240,
                                drainage_bearing = 0, flow_convergence = 0.8)
    site_off <- .make_rim_site(rec_bearing = 60, rec_distance = 2000,
                                rel_elevation = 60, aspect = 240,
                                drainage_bearing = 0, flow_convergence = 0.8)
    d <- .make_rim_met(wind_direction_day = 240)  # wind from SW → downwind ENE
    r_off <- odour_exposure(d, site_off,
                             terrain_backend = "descriptors", rim_venting = FALSE)
    r_on  <- odour_exposure(d, site_on,
                             terrain_backend = "descriptors", rim_venting = TRUE)
    morning_hours <- (nrow(d) - 2L):nrow(d)
    expect_gt(max(r_on[morning_hours]), max(r_off[morning_hours]),
              label = "morning exposure rises when align_j > alignment and reach ≈ 1")
  })

  it("high summit + weak morning: reach stays low, venting suppressed vs rim_venting=FALSE", {
    # z_j = 300 m.  BLH barely rises (150 → 155 → 160 m), giving
    # cbl_growth ≈ [5, 5] and cbl_cumsum_max ≈ 10.
    # h_vent = pool_top + RIM_LIFT_COEF * 10.  Even with RIM_LIFT_COEF = 5:
    # h_vent ≈ pool_top + 50 ≪ 300 → reach ≈ 0 → 1a venting suppressed.
    # On morning hours: exposure with rim_venting=TRUE ≤ rim_venting=FALSE.
    site <- .make_rim_site(rec_bearing = 0, rec_distance = 2000,
                            rel_elevation = 300, aspect = 180)
    d_weak <- data.frame(
      wind_direction_10m     = c(0, 180, 180),    # night then morning from S
      wind_speed_10m         = rep(1.5, 3),
      direct_radiation       = c(0, 400, 400),
      cloud_cover            = c(0,  20,  20),
      boundary_layer_height  = c(150, 155, 160),  # tiny morning rise
      temperature_2m         = rep(5, 3),
      pressure_msl           = rep(1013, 3),
      precipitation          = rep(0, 3),
      relative_humidity_2m   = rep(60, 3),
      soil_moisture_0_to_1cm = rep(0.1, 3),
      soil_moisture_1_to_3cm = rep(0.1, 3)
    )
    r_off <- odour_exposure(d_weak, site, terrain_backend = "descriptors",
                             rim_venting = FALSE)
    r_on  <- odour_exposure(d_weak, site, terrain_backend = "descriptors",
                             rim_venting = TRUE)
    morning_hours <- 2:3
    expect_true(
      all(r_on[morning_hours] <= r_off[morning_hours] + 1e-6),
      label = "weak morning / high summit: rim_venting=TRUE ≤ FALSE (reach ≈ 0)"
    )
  })

  it("facing receptor (align≈1) more exposed than away-facing (align≈0) on morning hours", {
    # Both receptors: same bearing (0°), distance (2000 m), rel_elevation (60 m).
    # rec_facing: aspect = 180° (downslope toward source) → align_j = 1.
    # rec_away:   aspect = 0°   (downslope away from src) → align_j = 0.
    # With deep pool and rim_venting=TRUE: facing receptor >> away-facing.
    ox <- 335000; oy <- 6250000
    rec_x <- ox; rec_y <- oy + 2000L
    .rim_site_1rec <- function(asp) {
      mh_site(
        sf::st_sf(id = c("src", "rec"),
                  rel_elevation = c(NA_real_, 60),
                  aspect        = c(NA_real_, asp),
                  geometry = sf::st_sfc(sf::st_point(c(ox, oy)),
                                         sf::st_point(c(rec_x, rec_y)), crs = 32755)),
        data.frame(feature_id = c("src", "rec"), hazard = "odour",
                   role = c("source", "receptor"), stringsAsFactors = FALSE),
        terrain = mh_terrain(drainage_bearing = 0, flow_convergence = 0.8,
                              valley_depth = 80, taf = 1.5),
        epsg = 32755L
      )
    }
    d  <- .make_rim_met(wind_direction_day = 180)   # wind from S → plume to N
    r_facing <- odour_exposure(d, .rim_site_1rec(180),
                                terrain_backend = "descriptors", rim_venting = TRUE)
    r_away   <- odour_exposure(d, .rim_site_1rec(0),
                                terrain_backend = "descriptors", rim_venting = TRUE)
    morning  <- nrow(d)   # last hour (peak morning)
    expect_gt(r_facing[morning], r_away[morning],
              label = "facing slope (align≈1) more exposed than away-facing (align≈0)")
  })

  it("produces finite, non-NaN values under near-calm venting morning (composition with #25)", {
    # Calm floor (u_eff ≈ U_CALM_FLOOR) × reach ≈ 1 (low z_j, deep pool):
    # hz = G*PM*W_rain/(u_eff*h_mix) is finite at the floor; result must be finite.
    site <- .make_rim_site(rel_elevation = 10, aspect = 180)
    d    <- .make_rim_met()
    r    <- odour_exposure(d, site, terrain_backend = "descriptors", rim_venting = TRUE)
    expect_true(all(is.finite(r)), label = "no NaN/Inf under calm + reach≈1")
    expect_true(all(r >= 0 & r <= 100))
  })
})
