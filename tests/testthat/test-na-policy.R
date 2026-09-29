# Missing-input (NA) policy shared by dust_hazard(), litter_risk(),
# odour_risk() and generate_twl() (and the layers under them):
#   * a row with a missing required input returns NA for that row only;
#   * the other rows are exactly what a complete-data call gives them;
#   * one classed summary warning (meteoHazard_missing_input) per call;
#   * genuinely invalid (non-missing) values still error.
# Fixtures (.count_missing, .na_met, ...) live in helper-na-policy.R.

# ---- dust ------------------------------------------------------------------ #

describe("dust NA policy", {
  it("dust_hazard() returns NA only for rows with a missing input", {
    met <- .na_met()
    full <- dust_hazard(met, tyler_sieve_no = 48L, gust_factor = 1, air_density = "met")
    for (col in c("wind_speed_10m", "wind_gusts_10m", "soil_moisture_0_to_1cm",
                  "temperature_2m", "surface_pressure")) {
      m <- met
      m[[col]][3] <- NA
      res <- .count_missing(
        dust_hazard(m, tyler_sieve_no = 48L, gust_factor = 1, air_density = "met")
      )
      expect_equal(res$n, 1L, info = col)
      expect_true(is.na(res$value[3]), info = col)
      expect_identical(res$value[-3], full[-3], info = col)
    }
  })

  it("dust_flux() returns NA for missing rows with one warning", {
    expect_warning(
      out <- dust_flux(20L, 10, c(0, NA), c(20, 20), c(0.02, 0.02)),
      class = "meteoHazard_missing_input"
    )
    expect_true(is.na(out[2]))
    expect_identical(out[1], dust_flux(20L, 10, 0, 20, 0.02))
  })

  it("a fully-sheltered surface still returns NA (not 0) for missing rows", {
    out <- suppressWarnings(
      dust_flux(20L, 10, c(5, NA), c(20, 20), c(0.02, 0.02), z0 = 0.5)
    )
    expect_equal(out[1], 0)
    expect_true(is.na(out[2]))
  })

  it("crust = TRUE treats a missing precipitation hour as missing", {
    met <- .na_met()
    met$precipitation[2] <- NA
    res <- .count_missing(dust_hazard(met, crust = TRUE))
    expect_equal(res$n, 1L)
    expect_true(is.na(res$value[2]))
    expect_false(anyNA(res$value[-2]))
    res_s <- .count_missing(dust_hazard(met, crust = TRUE, crust_decay = "saltation"))
    expect_true(is.na(res_s$value[2]))
    expect_false(anyNA(res_s$value[-2]))
  })

  it("still errors on genuinely invalid values", {
    met <- .na_met()
    met$wind_speed_10m[2] <- -1
    expect_error(dust_hazard(met))
    met <- .na_met()
    met$soil_moisture_0_to_1cm[2] <- 1.5
    expect_error(dust_hazard(met))
  })

  it("is silent when nothing is missing", {
    expect_no_warning(dust_hazard(.na_met()))
  })
})

# ---- litter ---------------------------------------------------------------- #

describe("litter NA policy", {
  it("litter_risk() returns NA only for rows with a missing input, one warning", {
    met  <- .na_met()
    site <- .na_litter_site()
    full <- litter_risk(met, site)
    for (col in c("wind_speed_10m", "wind_gusts_10m", "precipitation",
                  "soil_moisture_0_to_1cm", "wind_direction_10m")) {
      m <- met
      m[[col]][4] <- NA
      res <- .count_missing(litter_risk(m, site))
      expect_equal(res$n, 1L, info = col)
      expect_true(all(is.na(unlist(res$value[4, ]))), info = col)
      expect_identical(res$value[-4, ], full[-4, ], info = col)
    }
  })

  it("litter_hazard_vec() propagates NA per row", {
    expect_warning(
      out <- litter_hazard_vec(c(12, NA), c(8, 8), c(0, 0), c(0.1, 0.1)),
      class = "meteoHazard_missing_input"
    )
    expect_true(is.na(out[2]))
    expect_identical(out[1], litter_hazard_vec(12, 8, 0, 0.1))
  })

  it("litter_exposure() returns an all-NA row for a missing hazard or direction", {
    site <- .na_litter_site()
    full <- litter_exposure(c(30, 60, 10), c(270, 90, 0), site)
    res <- .count_missing(litter_exposure(c(30, NA, 10), c(270, 90, NA), site))
    expect_equal(res$n, 1L)
    expect_true(all(is.na(unlist(res$value[2:3, ]))))
    expect_identical(res$value[1, ], full[1, ])
    expect_s3_class(res$value$zone, "ordered")
  })

  it("litter_wetness_vec() gives NA for a missing hour and carries the state on", {
    args <- list(
      precipitation = c(1, 0, 0, 0), temperature_2m = c(20, 20, NA, 20),
      relative_humidity_2m = rep(50, 4), wind_speed_10m = rep(2, 4),
      shortwave_radiation = rep(0, 4)
    )
    res <- .count_missing(do.call(litter_wetness_vec, args))
    expect_equal(res$n, 1L)
    expect_true(is.na(res$value[3]))
    expect_equal(res$value[1:2], do.call(litter_wetness_vec,
                                         lapply(args, `[`, 1:2)))
    # The unknown hour holds the state: hour 4 dries once from hour 2's state.
    expect_equal(res$value[4], do.call(litter_wetness_vec,
                                       lapply(args, `[`, c(1, 2, 4)))[3])
  })

  it("litter_hazard_vec() still errors on invalid values", {
    expect_error(litter_hazard_vec(c(12, -1), c(8, 8), c(0, 0), c(0.1, 0.1)))
    expect_error(litter_hazard_vec(c(12, 12), c(8, 8), c(0, 0), c(0.1, 1.5)))
  })
})

# ---- odour ----------------------------------------------------------------- #

describe("odour NA policy", {
  it("odour_risk() returns an NA row for a missing required input", {
    met  <- .na_met()
    site <- .na_odour_site()
    for (col in c("wind_speed_10m", "wind_direction_10m", "direct_radiation",
                  "cloud_cover", "boundary_layer_height", "temperature_2m",
                  "relative_humidity_2m", "pressure_msl", "precipitation",
                  "soil_moisture_0_to_1cm", "soil_moisture_1_to_3cm")) {
      m <- met
      m[[col]][5] <- NA
      res <- .count_missing(odour_risk(m, site))
      expect_equal(res$n, 1L, info = col)
      expect_true(all(is.na(res$value[5, ])), info = col)
      expect_false(anyNA(res$value[-5, ]), info = col)
    }
  })

  it("odour_hazard() returns NA for a missing input row", {
    met <- .na_met()
    met$boundary_layer_height[2] <- NA
    res <- .count_missing(odour_hazard(met))
    expect_equal(res$n, 1L)
    expect_true(is.na(res$value[2]))
    expect_false(anyNA(res$value[-2]))
  })

  it("is silent when nothing is missing", {
    expect_no_warning(odour_risk(.na_met(), .na_odour_site()))
  })
})

# ---- TWL ------------------------------------------------------------------- #

describe("generate_twl() NA policy", {
  twl <- function(met, wind_height = 10) {
    generate_twl(
      datetime = met$datetime, latitude = -33.7, longitude = 150.3,
      temp = met$temperature_2m, wind_speed = met$wind_speed_10m / 3,
      wind_height = wind_height, RH = met$relative_humidity_2m,
      direct_solar = met$direct_radiation, diffuse_solar = met$diffuse_radiation,
      pressure = met$surface_pressure, verbose = FALSE
    )
  }

  it("returns NA only for rows with a missing input, with one warning", {
    met  <- .na_met()
    full <- twl(met)
    for (col in c("temperature_2m", "wind_speed_10m", "relative_humidity_2m",
                  "direct_radiation", "diffuse_radiation", "surface_pressure")) {
      m <- met
      m[[col]][2] <- NA
      res <- .count_missing(twl(m))
      expect_equal(res$n, 1L, info = col)
      expect_true(is.na(res$value[2]), info = col)
      expect_identical(res$value[-2], full[-2], info = col)
    }
  })

  it("treats an NA wind_height beside a known wind as a missing input", {
    met <- .na_met()
    res <- .count_missing(twl(met, wind_height = c(10, NA, 10, 10, 10, 10)))
    expect_equal(res$n, 1L)
    expect_true(is.na(res$value[2]))
    expect_identical(res$value[-2], twl(met)[-2])
  })

  it("treats an NA datetime or location as a missing input, not a crash", {
    met <- .na_met()
    full <- twl(met)
    m <- met
    m$datetime[2] <- NA
    res <- .count_missing(generate_twl(
      datetime = m$datetime, latitude = c(-33.7, -33.7, NA, -33.7, -33.7, -33.7),
      longitude = 150.3, temp = m$temperature_2m,
      wind_speed = m$wind_speed_10m / 3, wind_height = 10,
      RH = m$relative_humidity_2m, direct_solar = m$direct_radiation,
      diffuse_solar = m$diffuse_radiation, pressure = m$surface_pressure,
      verbose = FALSE
    ))
    expect_equal(res$n, 1L)
    expect_true(all(is.na(res$value[2:3])))
    expect_identical(res$value[-(2:3)], full[-(2:3)])
  })

  it("is silent when nothing is missing", {
    expect_no_warning(twl(.na_met()))
  })
})

