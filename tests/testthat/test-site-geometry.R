# Behaviour spec for the shared site geometry helpers (C1, issue #14):
# .bearing_distance(), .relative_elevation(), .role_features(),
# .crosswind_halfwidth(). All computed in a projected metric CRS.

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# Origin point and directional offsets in EPSG:32755 (UTM zone 55S)
# We construct points manually in the projected CRS to get exact offsets.
.origin_pt <- function() {
  sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
}

.point_east <- function(d = 500) {
  sf::st_sfc(sf::st_point(c(335000 + d, 6250000)), crs = 32755)
}

.point_north <- function(d = 300) {
  sf::st_sfc(sf::st_point(c(335000, 6250000 + d)), crs = 32755)
}

.point_south <- function(d = 200) {
  sf::st_sfc(sf::st_point(c(335000, 6250000 - d)), crs = 32755)
}

# A unit-square polygon (side = 100 m) centred on the origin, aligned N-S.
# Vertices at corners going clockwise.
.square_poly <- function(side = 100) {
  h <- side / 2
  coords <- matrix(
    c(-h,  h,
       h,  h,
       h, -h,
      -h, -h,
      -h,  h),          # close the ring
    ncol = 2, byrow = TRUE
  )
  sf::st_sfc(sf::st_polygon(list(coords)), crs = 32755)
}

# A degenerate polygon collapsed to a point.
.point_poly <- function() {
  coords <- matrix(c(0, 0, 0, 0, 0, 0), ncol = 2, byrow = TRUE)
  sf::st_sfc(sf::st_polygon(list(coords)), crs = 32755)
}

# An mh_site with two features and a known roles table.
.make_site_for_roles <- function() {
  feats <- sf::st_sf(
    id = c("src1", "rec1", "barrier1"),
    elevation = c(10, 25, 15),
    geometry = sf::st_sfc(
      sf::st_point(c(335000, 6250000)),
      sf::st_point(c(335500, 6250000)),
      sf::st_point(c(335250, 6250000)),
      crs = 32755
    )
  )
  roles <- data.frame(
    feature_id = c("src1",   "rec1",    "barrier1"),
    hazard     = c("odour",  "odour",   "odour"),
    role       = c("source", "receptor","barrier"),
    stringsAsFactors = FALSE
  )
  mh_site(feats, roles, epsg = 32755)
}

# ---------------------------------------------------------------------------
# .bearing_distance()
# ---------------------------------------------------------------------------

describe(".bearing_distance()", {

  it("returns a bearing in degrees from north (clockwise) and a distance in metres", {
    result <- .bearing_distance(.origin_pt(), .point_east(500))
    expect_type(result, "list")
    expect_named(result, c("bearing", "distance"))
    expect_type(result$bearing,  "double")
    expect_type(result$distance, "double")
  })

  it("gives bearing 90 and the coordinate separation for a receptor due east", {
    result <- .bearing_distance(.origin_pt(), .point_east(500))
    expect_equal(result$bearing,  90, tolerance = 1e-6)
    expect_equal(result$distance, 500, tolerance = 1e-6)
  })

  it("gives bearing 0 due north and 180 due south", {
    north <- .bearing_distance(.origin_pt(), .point_north(300))
    south <- .bearing_distance(.origin_pt(), .point_south(200))
    expect_equal(north$bearing, 0,   tolerance = 1e-6)
    expect_equal(south$bearing, 180, tolerance = 1e-6)
  })

  it("is symmetric in distance between two features", {
    a <- .origin_pt()
    b <- .point_east(500)
    d_ab <- .bearing_distance(a, b)$distance
    d_ba <- .bearing_distance(b, a)$distance
    expect_equal(d_ab, d_ba, tolerance = 1e-9)
  })

  it("gives a reverse bearing that differs by 180 degrees", {
    a <- .origin_pt()
    b <- .point_east(500)
    brg_ab <- .bearing_distance(a, b)$bearing
    brg_ba <- .bearing_distance(b, a)$bearing
    diff_mod <- abs(brg_ab - brg_ba) %% 360
    expect_equal(min(diff_mod, 360 - diff_mod), 180, tolerance = 1e-6)
  })

  it("returns distance 0 and bearing NA for coincident points", {
    result <- .bearing_distance(.origin_pt(), .origin_pt())
    expect_equal(result$distance, 0)
    expect_true(is.na(result$bearing))
  })

})

# ---------------------------------------------------------------------------
# .relative_elevation()
# ---------------------------------------------------------------------------

describe(".relative_elevation()", {

  it("returns receptor elevation minus the source base on the AGL datum", {
    src <- data.frame(elevation = 10)
    rec <- data.frame(elevation = 25)
    expect_equal(.relative_elevation(src, rec), 15)
  })

  it("propagates NA when an elevation is missing", {
    src_na  <- data.frame(elevation = NA_real_)
    rec_ok  <- data.frame(elevation = 25)
    src_ok  <- data.frame(elevation = 10)
    rec_na  <- data.frame(elevation = NA_real_)

    expect_true(is.na(.relative_elevation(src_na, rec_ok)))
    expect_true(is.na(.relative_elevation(src_ok, rec_na)))
  })

})

# ---------------------------------------------------------------------------
# .role_features()
# ---------------------------------------------------------------------------

describe(".role_features()", {

  it("returns exactly the features whose roles match the hazard and role", {
    site    <- .make_site_for_roles()
    sources <- .role_features(site, "odour", "source")
    expect_true(inherits(sources, "sf"))
    expect_equal(nrow(sources), 1L)
    expect_equal(sources$id, "src1")

    receptors <- .role_features(site, "odour", "receptor")
    expect_equal(nrow(receptors), 1L)
    expect_equal(receptors$id, "rec1")

    barriers <- .role_features(site, "odour", "barrier")
    expect_equal(nrow(barriers), 1L)
    expect_equal(barriers$id, "barrier1")
  })

  it("returns a zero-row sf for an absent hazard or role", {
    site <- .make_site_for_roles()

    absent_hazard <- .role_features(site, "dust", "source")
    expect_true(inherits(absent_hazard, "sf"))
    expect_equal(nrow(absent_hazard), 0L)

    absent_role <- .role_features(site, "odour", "sink")
    expect_true(inherits(absent_role, "sf"))
    expect_equal(nrow(absent_role), 0L)
  })

})

# ---------------------------------------------------------------------------
# .crosswind_halfwidth()
# ---------------------------------------------------------------------------

describe(".crosswind_halfwidth()", {

  it("returns half the footprint extent perpendicular to the wind", {
    poly   <- .square_poly(side = 100)
    result <- .crosswind_halfwidth(poly, wind_dir = 0)   # wind from north
    expect_type(result, "double")
    expect_length(result, 1L)
  })

  it("gives half the side for a square aligned N-S with wind from the north", {
    # Wind from north (0 deg): crosswind is east-west → width = side
    poly <- .square_poly(side = 100)
    hw   <- .crosswind_halfwidth(poly, wind_dir = 0)
    expect_equal(hw, 50, tolerance = 1e-9)
  })

  it("gives half of side*sqrt(2) for that square with wind from the north-east", {
    # Wind from NE (45 deg): crosswind width = side * sqrt(2)
    poly <- .square_poly(side = 100)
    hw   <- .crosswind_halfwidth(poly, wind_dir = 45)
    expect_equal(hw, 100 * sqrt(2) / 2, tolerance = 1e-9)
  })

  it("is non-negative", {
    poly <- .square_poly(side = 100)
    for (wd in c(0, 45, 90, 135, 180, 270, 315)) {
      expect_gte(.crosswind_halfwidth(poly, wd), 0)
    }
  })

  it("is symmetric under a 180-degree wind reversal", {
    poly <- .square_poly(side = 100)
    for (wd in c(0, 30, 45, 90, 135)) {
      hw1 <- .crosswind_halfwidth(poly, wd)
      hw2 <- .crosswind_halfwidth(poly, wd + 180)
      expect_equal(hw1, hw2, tolerance = 1e-9)
    }
  })

  it("returns 0 for a degenerate (point) geometry", {
    degen <- .point_poly()
    expect_equal(.crosswind_halfwidth(degen, wind_dir = 45), 0)
  })

})
