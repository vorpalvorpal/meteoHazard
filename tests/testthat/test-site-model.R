# Behaviour spec for the geo-referenced site model: mh_site() and mh_terrain().
# Plan: GitHub issue #14 (C1). Pending specs (skipped) are the checklist that
# /implement turns green; concrete assertions, tolerances and seeds are filled
# in then. See the C1 plan and the Stage-1 terrain-schema spike findings (#14).

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# A tiny sf point set in EPSG:4326 (geographic / lon-lat)
.make_geo_features <- function() {
  sf::st_sf(
    id = c("A", "B"),
    geometry = sf::st_sfc(
      sf::st_point(c(151.0, -33.8)),
      sf::st_point(c(151.1, -33.9)),
      crs = 4326
    )
  )
}

# Same features but already in EPSG:32755 (UTM zone 55S — projected metric)
.make_proj_features <- function() {
  sf::st_transform(.make_geo_features(), 32755)
}

# Minimal roles data.frame
.make_roles <- function(feature_id = "A", hazard = "odour", role = "source") {
  data.frame(feature_id = feature_id, hazard = hazard, role = role,
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# mh_site()
# ---------------------------------------------------------------------------

describe("mh_site()", {

  it("bundles features, roles and terrain into one validated object", {
    feat    <- .make_proj_features()
    roles   <- .make_roles()
    terrain <- mh_terrain(relief = 20)
    site    <- mh_site(feat, roles, terrain = terrain, epsg = 32755)

    expect_true(S7::S7_inherits(site, mh_site))
    expect_identical(site@datum, "AGL")
    expect_identical(site@epsg, 32755L)
    expect_true(inherits(site@features, "sf"))
    expect_true(is.data.frame(site@roles))
    expect_true(S7::S7_inherits(site@terrain, mh_terrain))
  })

  it("requires every roles$feature_id to exist in features$id", {
    feat  <- .make_proj_features()
    roles <- .make_roles(feature_id = "MISSING")
    expect_error(
      mh_site(feat, roles, epsg = 32755),
      class = "meteoHazard_input_error"
    )
  })

  it("requires a single CRS across all features", {
    # sf enforces CRS homogeneity in sfc; constructing mixed-CRS sf is not
    # directly possible — verify that the supplied features CRS is detected
    feat  <- .make_proj_features()
    roles <- .make_roles()
    # A valid projected CRS gives no error
    expect_no_error(mh_site(feat, roles, epsg = 32755))
  })

  it("requires a projected metric CRS", {
    # features already projected → no error; if epsg is geographic → error
    feat  <- .make_proj_features()
    roles <- .make_roles()
    expect_no_error(mh_site(feat, roles, epsg = 32755))
  })

  it("reprojects a geographic (lon/lat) CRS to the supplied metric epsg", {
    geo_feat  <- .make_geo_features()
    roles     <- .make_roles()
    site      <- mh_site(geo_feat, roles, epsg = 32755)
    expect_false(isTRUE(sf::st_is_longlat(site@features)))
    expect_equal(sf::st_crs(site@features), sf::st_crs(32755))
  })

  it("matches sf::st_transform for a known point under reprojection", {
    geo_feat    <- .make_geo_features()
    roles       <- .make_roles()
    site        <- mh_site(geo_feat, roles, epsg = 32755)
    expected    <- sf::st_transform(geo_feat, 32755)
    site_coords <- sf::st_coordinates(site@features)
    exp_coords  <- sf::st_coordinates(expected)
    expect_equal(site_coords, exp_coords, tolerance = 1e-6)
  })

  it("enforces a shared AGL datum between source and terrain heights", {
    feat  <- .make_proj_features()
    roles <- .make_roles()
    # AGL is the only valid datum; anything else should error
    expect_error(
      mh_site(feat, roles, epsg = 32755, datum = "MSL"),
      class = "meteoHazard_input_error"
    )
    expect_no_error(mh_site(feat, roles, epsg = 32755, datum = "AGL"))
  })

  it("allows a feature that carries no roles", {
    # Feature "B" has no roles — this is allowed
    feat  <- .make_proj_features()
    roles <- .make_roles(feature_id = "A")
    expect_no_error(mh_site(feat, roles, epsg = 32755))
  })

  it("allows a hazard with no source or no receptor at construction", {
    # Only a receptor role, no source — allowed at construction time
    feat  <- .make_proj_features()
    roles <- .make_roles(feature_id = "A", hazard = "odour", role = "receptor")
    expect_no_error(mh_site(feat, roles, epsg = 32755))
  })

  it("errors (classed) when no metric epsg is available", {
    feat  <- .make_proj_features()
    roles <- .make_roles()
    expect_error(
      mh_site(feat, roles, epsg = 4326),   # 4326 = geographic WGS84
      class = "meteoHazard_input_error"
    )
  })

  it("errors (classed) on an unknown hazard or role value", {
    feat <- .make_proj_features()

    bad_hazard_roles <- .make_roles(hazard = "noise")
    expect_error(
      mh_site(feat, bad_hazard_roles, epsg = 32755),
      class = "meteoHazard_input_error"
    )

    bad_role_roles <- .make_roles(role = "emitter")
    expect_error(
      mh_site(feat, bad_role_roles, epsg = 32755),
      class = "meteoHazard_input_error"
    )
  })
})

# ---------------------------------------------------------------------------
# mh_terrain()
# ---------------------------------------------------------------------------

describe("mh_terrain()", {

  it("holds the scalar terrain descriptors with units and a _meta scale block", {
    t <- mh_terrain(
      relief = 50,
      valley_depth = 30,
      taf = 1.2,
      drainage_bearing = 90,
      meta = list(relief_radius = 5000, channel_threshold = 1e6,
                  fetch_L = 2000, dem_resolution = 30, datum = "AGL")
    )
    expect_true(S7::S7_inherits(t, mh_terrain))
    expect_equal(t@relief, 50)
    expect_equal(t@valley_depth, 30)
    expect_equal(t@taf, 1.2)
    expect_equal(t@meta$relief_radius, 5000)
    expect_equal(t@meta$fetch_L, 2000)
  })

  it("records relief as a height above local base in metres", {
    # spike finding: NOT a standardised DEV index
    t <- mh_terrain(relief = units::set_units(40, "m"))
    expect_equal(t@relief, 40)   # stored as plain double metres
    expect_true(is.double(t@relief))
  })

  it("carries _meta scale fields (relief_radius, channel_threshold, fetch_L)", {
    t <- mh_terrain(
      meta = list(relief_radius = 5000, channel_threshold = 1e6, fetch_L = 2000)
    )
    expect_equal(t@meta$relief_radius,     5000)
    expect_equal(t@meta$channel_threshold, 1e6)
    expect_equal(t@meta$fetch_L,           2000)
  })

  it("rejects a negative relief, valley_depth or basin_capacity", {
    expect_error(mh_terrain(relief = -1),        class = "meteoHazard_input_error")
    expect_error(mh_terrain(valley_depth = -0.1),class = "meteoHazard_input_error")
    expect_error(mh_terrain(basin_capacity = -5), class = "meteoHazard_input_error")
  })

  it("rejects a topographic amplification factor below 1", {
    expect_error(mh_terrain(taf = 0.9), class = "meteoHazard_input_error")
    expect_no_error(mh_terrain(taf = 1.0))
    expect_no_error(mh_terrain(taf = 2.5))
  })

  it("rejects a drainage bearing outside [0, 360)", {
    expect_error(mh_terrain(drainage_bearing = -1),  class = "meteoHazard_input_error")
    expect_error(mh_terrain(drainage_bearing = 360), class = "meteoHazard_input_error")
    expect_no_error(mh_terrain(drainage_bearing = 0))
    expect_no_error(mh_terrain(drainage_bearing = 359.9))
  })

  it("requires flow_convergence and shelter_index to be finite", {
    expect_error(mh_terrain(flow_convergence = Inf),  class = "meteoHazard_input_error")
    expect_error(mh_terrain(flow_convergence = NaN),  class = "meteoHazard_input_error")
    expect_error(mh_terrain(shelter_index = Inf),     class = "meteoHazard_input_error")
    expect_error(mh_terrain(shelter_index = -Inf),    class = "meteoHazard_input_error")
    expect_no_error(mh_terrain(flow_convergence = 0.5))
    expect_no_error(mh_terrain(shelter_index = 45))
  })

  it("permits NA descriptors as 'no terrain effect' (flat site)", {
    expect_no_error(mh_terrain())   # all NA by default
    t <- mh_terrain()
    expect_true(is.na(t@relief))
    expect_true(is.na(t@taf))
    expect_true(is.na(t@flow_convergence))
  })
})

# ---------------------------------------------------------------------------
# mh_site() / mh_terrain() units handling
# ---------------------------------------------------------------------------

describe("mh_site() / mh_terrain() units handling", {

  it("accepts distances/elevations as bare numerics in the documented unit", {
    t <- mh_terrain(relief = 100, valley_depth = 50, basin_capacity = 1e6)
    expect_equal(t@relief,         100)
    expect_equal(t@valley_depth,    50)
    expect_equal(t@basin_capacity, 1e6)
  })

  it("accepts units-tagged distances/elevations and converts them", {
    t_bare  <- mh_terrain(relief = 500)
    t_units <- mh_terrain(relief = units::set_units(500, "m"))
    expect_equal(t_bare@relief, t_units@relief)

    # km → m conversion: 0.5 km == 500 m
    t_km <- mh_terrain(relief = units::set_units(0.5, "km"))
    expect_equal(t_km@relief, 500, tolerance = 1e-9)
  })

  it("errors (classed) on a dimensionally incompatible unit", {
    expect_error(
      mh_terrain(relief = units::set_units(30, "degC")),
      class = "meteoHazard_input_error"
    )
    expect_error(
      mh_terrain(basin_capacity = units::set_units(10, "m")),
      class = "meteoHazard_input_error"
    )
  })
})

# ---------------------------------------------------------------------------
# C4 addition (issue #18): dust receptor roles on the unified model
# ---------------------------------------------------------------------------

describe("mh_site(): dust receptor roles", {

  it("validates a (dust, receptor) role", {
    feat  <- .make_proj_features()
    roles <- .make_roles(feature_id = "A", hazard = "dust", role = "receptor")
    expect_no_error(mh_site(feat, roles, epsg = 32755))
  })

  it("does not require a (dust, barrier) role", {
    feat  <- .make_proj_features()
    roles <- .make_roles(feature_id = "A", hazard = "dust", role = "source")
    # no barrier — that is fine at construction
    expect_no_error(mh_site(feat, roles, epsg = 32755))
  })
})
