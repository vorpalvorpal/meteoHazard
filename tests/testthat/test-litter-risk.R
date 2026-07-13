# Behaviour specification for litter_risk() (R-D1).
#
# Encodes scratchpad/tdd/PLAN.md section 0.2/1 (R-D1). litter_risk() is the
# convenience wrapper (mirroring odour_risk()) that computes the hazard and
# maps it through litter_exposure() in one call:
#   litter_risk(met_data, site) == litter_exposure(litter_hazard(met_data),
#                                                   met_data$wind_direction_10m,
#                                                   site)
#
# R/litter_risk.R does not exist yet (new file): every test below is expected
# to fail with "could not find function 'litter_risk'" until it is created --
# that is the intended TDD-red state; no skip-gating is applied.

.risk_site <- function(epsg = 32755, ring_radius = 1000) {
  ctr <- sf::st_sf(
    id       = "ctr",
    geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = epsg)
  )
  sectors <- data.frame(
    arc_start    = c("NE", "SW"),
    arc_end      = c("SE", "NW"),
    permeability = c(1.0, 0.3),
    sensitive    = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  # suppressWarnings: partial-coverage fixture (sector-gap warning is exercised
  # in test-site-constructors.R, not here).
  suppressWarnings(site_from_sectors(sectors, ctr, ring_radius = ring_radius, epsg = epsg))
}

.risk_met <- function() {
  data.frame(
    wind_gusts_10m         = c(6, 12, 16),
    wind_speed_10m         = c(3, 9, 15),
    precipitation          = c(0, 0, 0),
    soil_moisture_0_to_1cm = c(0.02, 0.02, 0.02),
    wind_direction_10m     = c(270, 90, 0)
  )
}

describe("litter_risk()", {

  it("R-D1: equals the two-step litter_hazard() + litter_exposure() call, column for column", {
    met  <- .risk_met()
    site <- .risk_site()

    combined <- litter_risk(met, site)
    two_step <- litter_exposure(litter_hazard(met), met$wind_direction_10m, site)

    expect_equal(combined, two_step)
  })

  it("R-D1: forwards litter_exposure()'s named parameters (direction_tol etc.) identically to the two-step call", {
    met  <- .risk_met()
    site <- .risk_site()

    combined <- litter_risk(met, site, direction_tol = 20, p_open_min = 0.6,
                             move_threshold = 15, offsite_threshold = 40,
                             default_permeability = 0.2)
    two_step <- litter_exposure(
      litter_hazard(met), met$wind_direction_10m, site,
      direction_tol = 20, p_open_min = 0.6,
      move_threshold = 15, offsite_threshold = 40,
      default_permeability = 0.2
    )

    expect_equal(combined, two_step)
  })

  it("R-D1: errors, naming wind_direction_10m, when that column is absent", {
    met <- .risk_met()
    met$wind_direction_10m <- NULL
    site <- .risk_site()

    expect_error(litter_risk(met, site), regexp = "wind_direction_10m")
  })

  it("R-D1: use_wetness_state = TRUE runs the wetness path and matches the manual two-step call", {
    met <- data.frame(
      wind_gusts_10m        = 16,
      wind_speed_10m        = 15,
      precipitation         = 0,
      temperature_2m        = 25,
      relative_humidity_2m  = 40,
      shortwave_radiation   = 600,
      wind_direction_10m    = 270
    )
    site <- .risk_site()

    combined <- litter_risk(met, site, use_wetness_state = TRUE)
    two_step <- litter_exposure(
      litter_hazard(met, use_wetness_state = TRUE), met$wind_direction_10m, site
    )

    expect_equal(combined, two_step)
  })

  it("R-D1: reach_per_ms forwards mean_wind = met_data$wind_speed_10m to litter_exposure()", {
    met <- data.frame(
      wind_gusts_10m         = 16,
      wind_speed_10m         = 12,
      precipitation          = 0,
      soil_moisture_0_to_1cm = 0.02,
      wind_direction_10m     = 270
    )
    site <- .risk_site()  # NE-SE barrier distance_m == ring_radius == 1000

    out <- litter_risk(met, site, reach_per_ms = 100)
    hazard <- litter_hazard(met)
    two_step <- litter_exposure(
      hazard, met$wind_direction_10m, site,
      mean_wind = met$wind_speed_10m, reach_per_ms = 100
    )

    expect_equal(out, two_step)
    # mean_wind 12 * reach_per_ms 100 = reach 1200 >= distance_m 1000 -> reaches -> off_site
    expect_equal(as.character(out$zone), "off_site")
    expect_true(out$leaves_site)
  })

  it("without reach_per_ms, basic mode is used (mean_wind is not forwarded)", {
    met <- data.frame(
      wind_gusts_10m         = 16,
      wind_speed_10m         = 12,
      precipitation          = 0,
      soil_moisture_0_to_1cm = 0.02,
      wind_direction_10m     = 270
    )
    site <- .risk_site()

    out <- litter_risk(met, site)
    two_step <- litter_exposure(litter_hazard(met), met$wind_direction_10m, site)
    expect_equal(out, two_step)
  })
})
