# Behaviour spec for the site convenience constructors (C4, issue #18):
# site_from_sectors() — compass-sector descriptions -> features geometry.
# Pending specs (skipped) are the checklist /implement turns green.
#
# Extended per scratchpad/tdd/PLAN.md (Litter pipeline v3.1) with:
#   R-A4 - full-circle sector (arc_start == arc_end yields a real ring)
#   R-A7 - sector-coverage-gap warning
#   distance_m - optional per-sector reach distance (defaults to ring_radius),
#                carried onto the barrier features/roles (spec S3.4)

describe("site_from_sectors()", {

  it("places a barrier from arc [NE, SE] spanning bearings 45-135 from the centroid", {
    ctr <- sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    sectors <- data.frame(
      arc_start    = "NE",
      arc_end      = "SE",
      permeability = 0.5,
      sensitive    = FALSE,
      stringsAsFactors = FALSE
    )
    site <- site_from_sectors(sectors, ctr, ring_radius = 500, epsg = 32755)

    barrier <- site@features[site@features$id == "barrier_1", ]
    coords  <- sf::st_coordinates(barrier)
    cx <- 335000; cy <- 6250000
    bearings <- atan2(coords[, "X"] - cx, coords[, "Y"] - cy) * 180 / pi %% 360
    expect_true(min(bearings) <= 45 + 5)   # starts near 45
    expect_true(max(bearings) >= 135 - 5)  # ends near 135
  })

  it("carries permeability and sensitive onto the feature roles", {
    ctr <- sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    sectors <- data.frame(
      arc_start    = "N",
      arc_end      = "E",
      permeability = 0.3,
      sensitive    = TRUE,
      stringsAsFactors = FALSE
    )
    site <- site_from_sectors(sectors, ctr, ring_radius = 500, epsg = 32755)

    barrier_role <- site@roles[site@roles$feature_id == "barrier_1", ]
    expect_equal(barrier_role$permeability, 0.3)
    expect_true(barrier_role$sensitive)
  })

  it("wraps an arc that crosses north (e.g. [NW, NE]) correctly", {
    ctr <- sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    sectors <- data.frame(
      arc_start    = "NW",
      arc_end      = "NE",
      permeability = 0.5,
      sensitive    = FALSE,
      stringsAsFactors = FALSE
    )
    site <- site_from_sectors(sectors, ctr, ring_radius = 500, epsg = 32755)

    barrier <- site@features[site@features$id == "barrier_1", ]
    coords   <- sf::st_coordinates(barrier)
    cx <- 335000; cy <- 6250000
    northmost <- coords[, "Y"] - cy   # max should be close to ring_radius = 500
    expect_gt(max(northmost), 400)    # the barrier extends toward north
  })

  it("requires an explicit origin (centroid / working face)", {
    sectors <- data.frame(
      arc_start    = "N",
      arc_end      = "E",
      permeability = 0.3,
      sensitive    = FALSE,
      stringsAsFactors = FALSE
    )
    expect_error(
      site_from_sectors(sectors, centroid = NULL, epsg = 32755),
      class = "meteoHazard_input_error"
    )
  })
})


describe("site_from_sectors() full-circle sector (R-A4)", {

  .full_circle_ctr <- function() {
    sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
  }

  it("arc_start == arc_end yields a real ring, area ~ pi * ring_radius^2", {
    sectors <- data.frame(
      arc_start = "N", arc_end = "N",
      permeability = 1.0, sensitive = TRUE,
      stringsAsFactors = FALSE
    )
    site <- site_from_sectors(sectors, .full_circle_ctr(), ring_radius = 500, epsg = 32755)

    barrier <- site@features[site@features$id == "barrier_1", ]
    area <- as.numeric(sf::st_area(barrier))
    expect_gt(area, 0)
    expect_equal(area, pi * 500^2, tolerance = 0.05)
  })

  it("the full-circle barrier is hit from all 8 principal wind directions", {
    sectors <- data.frame(
      arc_start = "N", arc_end = "N",
      permeability = 1.0, sensitive = TRUE,
      stringsAsFactors = FALSE
    )
    site <- site_from_sectors(sectors, .full_circle_ctr(), ring_radius = 500, epsg = 32755)

    for (d in seq(0, 315, 45)) {
      out <- litter_exposure(50, d, site, default_permeability = 0)
      expect_equal(out$exposure, 50 * 1.0, tolerance = 1e-3, info = paste("wind_direction_10m", d))
    }
  })
})


describe("site_from_sectors() sector-coverage-gap warning (R-A7)", {

  .gap_ctr <- function() {
    sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
  }

  it("a single NE-SE sector leaves ~270 degrees uncovered and warns meteoHazard_litter_sector_gap", {
    sectors <- data.frame(
      arc_start = "NE", arc_end = "SE",
      permeability = 1.0, sensitive = FALSE,
      stringsAsFactors = FALSE
    )
    expect_warning(
      site_from_sectors(sectors, .gap_ctr(), ring_radius = 500, epsg = 32755),
      class = "meteoHazard_litter_sector_gap"
    )
  })

  it("a full tiling of sectors (N-E, E-S, S-W, W-N) covers 360 degrees and does not warn", {
    sectors <- data.frame(
      arc_start    = c("N", "E", "S", "W"),
      arc_end      = c("E", "S", "W", "N"),
      permeability = c(0.5, 0.5, 0.5, 0.5),
      sensitive    = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    )
    expect_no_warning(
      site_from_sectors(sectors, .gap_ctr(), ring_radius = 500, epsg = 32755),
      class = "meteoHazard_litter_sector_gap"
    )
  })
})


describe("site_from_sectors() distance_m column (spec S3.4)", {

  .distance_ctr <- function() {
    sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
  }

  it("defaults distance_m to ring_radius when sectors$distance_m is absent", {
    sectors <- data.frame(
      arc_start = "NE", arc_end = "SE",
      permeability = 1.0, sensitive = FALSE,
      stringsAsFactors = FALSE
    )
    site <- site_from_sectors(sectors, .distance_ctr(), ring_radius = 750, epsg = 32755)

    expect_true("distance_m" %in% names(site@features))
    brow <- site@features[site@features$id == "barrier_1", ]
    expect_equal(brow$distance_m, 750)

    expect_true("distance_m" %in% names(site@roles))
    role_row <- site@roles[site@roles$feature_id == "barrier_1", ]
    expect_equal(role_row$distance_m, 750)
  })

  it("carries an explicit sectors$distance_m through to the barrier feature/role", {
    sectors <- data.frame(
      arc_start = "NE", arc_end = "SE",
      permeability = 1.0, sensitive = FALSE, distance_m = 250,
      stringsAsFactors = FALSE
    )
    site <- site_from_sectors(sectors, .distance_ctr(), ring_radius = 1000, epsg = 32755)

    brow <- site@features[site@features$id == "barrier_1", ]
    expect_equal(brow$distance_m, 250)
  })

  it("rejects a non-positive distance_m", {
    sectors <- data.frame(
      arc_start = "NE", arc_end = "SE",
      permeability = 1.0, sensitive = FALSE, distance_m = 0,
      stringsAsFactors = FALSE
    )
    expect_error(
      site_from_sectors(sectors, .distance_ctr(), ring_radius = 1000, epsg = 32755),
      class = "meteoHazard_input_error"
    )
  })
})
