# Behaviour spec for the site convenience constructors (C4, issue #18):
# site_from_sectors() — compass-sector descriptions -> features geometry.
# Pending specs (skipped) are the checklist /implement turns green.

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
