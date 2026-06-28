# Behaviour specification for the v2 dust hazard model (specs/Dust.md, issue #5).
#
# v2 replaces the v1 AP-42 erosion-potential model with a physical
# saltation->dust framework. Specs skip until the v2 API is present (marker: the
# `gust_factor` formal on dust_flux(), absent in v1).
#
# Wind inputs are in m/s (Open-Meteo fetched with &wind_speed_unit=ms). For
# sieve 20 / clay 10 / z0 0.005 the dry entrainment threshold sits at a gust of
# ~10.6 m/s. dust_hazard() returns a crust-adjusted RELATIVE dust flux (issue
# #11 removed the fixed 0-100 index); it shares dust_flux()'s relative units.

skip_if_no_dust_v2 <- function() {
  testthat::skip_if_not(
    exists("dust_flux", mode = "function") &&
      "gust_factor" %in% names(formals(dust_flux)),
    "dust v2 API not yet implemented"
  )
}

DUST_TOL <- 1e-3

# A dry soil moisture (gravimetric w < w' for clay 10, bulk 1.6): w = sm/1.6*100,
# w' = 1.84, so sm < 0.0294 is "dry". 0.02 is safely dry.


describe("dust_flux() [v2]", {

  it("returns zero flux below the entrainment threshold", {
    skip_if_no_dust_v2()
    # Light gust on dry coarse-sand surface: u* below u*t (threshold ~10.6 m/s).
    f <- dust_flux(20L, 10, wind_speed_10m = 1.5, wind_gusts_10m = 8,
                                 soil_moisture = 0.02)
    expect_equal(f, 0)
  })

  it("gives positive flux above threshold and increases monotonically with gust", {
    skip_if_no_dust_v2()
    f <- dust_flux(20L, 10, wind_speed_10m = rep(0, 5),
                                 wind_gusts_10m = c(12, 14, 16, 20, 25),
                                 soil_moisture = rep(0.02, 5))
    expect_true(all(f >= 0))
    expect_true(all(diff(f) >= 0))
    expect_gt(f[5], 0)
  })

  it("is vectorised — one flux per input hour", {
    skip_if_no_dust_v2()
    f <- dust_flux(20L, 10, c(0, 0, 0), c(12, 16, 22), c(0.02, 0.02, 0.02))
    expect_length(f, 3)
  })

  it("falls as soil moisture rises (Fécan threshold increase)", {
    skip_if_no_dust_v2()
    f <- dust_flux(20L, 10, wind_speed_10m = rep(0, 3),
                                 wind_gusts_10m = rep(20, 3),
                                 soil_moisture = c(0.02, 0.06, 0.10))
    expect_true(all(diff(f) <= 0))
    expect_gt(f[1], f[3])
  })

  it("raising threshold_multiplier reduces the flux (crust hook)", {
    skip_if_no_dust_v2()
    f_base  <- dust_flux(20L, 10, 0, 16, 0.02, threshold_multiplier = 1)
    f_crust <- dust_flux(20L, 10, 0, 16, 0.02, threshold_multiplier = 2)
    expect_lt(f_crust, f_base)
  })

  it("rejects an invalid Tyler sieve number", {
    skip_if_no_dust_v2()
    expect_error(dust_flux(7L, 10, 0, 16, 0.02))
  })

  it("rejects out-of-range clay, soil moisture, and missing values", {
    skip_if_no_dust_v2()
    expect_error(dust_flux(20L, 150, 0, 16, 0.02))   # clay > 100
    expect_error(dust_flux(20L, 10, 0, 16, 1.5))     # soil moisture > 1
    expect_error(dust_flux(20L, 10, c(0, NA), c(16, 16), c(0.02, 0.02)))
  })
})


describe("dust_hazard() [v2]", {

  dust_met <- function(gust = c(12, 16, 20), wind = 0, sm = 0.02, precip = NULL, n = NULL) {
    if (is.null(n)) n <- length(gust)
    df <- data.frame(
      wind_speed_10m         = rep_len(wind, n),
      wind_gusts_10m         = rep_len(gust, n),
      soil_moisture_0_to_1cm = rep_len(sm, n)
    )
    if (!is.null(precip)) df$precipitation <- rep_len(precip, n)
    df
  }

  it("returns one relative flux per row, non-negative", {
    skip_if_no_dust_v2()
    out <- dust_hazard(dust_met(gust = c(8, 14, 18, 25)))
    expect_length(out, 4)
    expect_true(all(out >= 0))
  })

  it("equals dust_flux() with crust off (a met-frame convenience wrapper)", {
    skip_if_no_dust_v2()
    met <- dust_met(gust = c(8, 14, 18, 25), wind = 0, sm = 0.02)
    out <- dust_hazard(met, crust = FALSE)
    ref <- dust_flux(20L, 10, wind_speed_10m = met$wind_speed_10m,
                     wind_gusts_10m = met$wind_gusts_10m,
                     soil_moisture = met$soil_moisture_0_to_1cm)
    expect_equal(out, ref, tolerance = DUST_TOL)
  })

  it("is zero for sub-threshold winds and positive above threshold", {
    skip_if_no_dust_v2()
    out <- dust_hazard(dust_met(gust = c(8, 25)))
    expect_equal(out[1], 0)
    expect_gt(out[2], 0)
  })

  it("increases monotonically with gust", {
    skip_if_no_dust_v2()
    out <- dust_hazard(dust_met(gust = c(12, 14, 16, 18, 25)))
    expect_true(all(diff(out) >= 0))
  })

  it("with crust enabled, a rain event suppresses dry hours that recover with time", {
    skip_if_no_dust_v2()
    # Hour 1 rains (>= threshold); the rest are dry. Crust raises the threshold,
    # decaying over crust_decay_hours. Use a long series so the late hours fully
    # recover, making both the suppression and the recovery observable.
    met <- dust_met(gust = 20, wind = 0, sm = 0.02,
                    precip = c(5, rep(0, 120)), n = 121)
    on  <- dust_hazard(met, crust = TRUE,
                                    crust_factor_max = 3, crust_decay_hours = 24)
    off <- dust_hazard(met, crust = FALSE)
    # Just after rain: crust suppresses relative to no-crust.
    expect_lt(on[2], off[2])
    # Long after rain: crust has decayed, so the hour recovers toward no-crust.
    expect_gt(on[121], on[2])
    expect_equal(on[121], off[121], tolerance = DUST_TOL)
  })

  it("ignores precipitation when crust is disabled", {
    skip_if_no_dust_v2()
    with_p    <- dust_hazard(dust_met(gust = 20, precip = 10), crust = FALSE)
    without_p <- dust_hazard(dust_met(gust = 20), crust = FALSE)
    expect_equal(with_p, without_p)
  })

  it("errors if crust is enabled but precipitation is absent", {
    skip_if_no_dust_v2()
    expect_error(dust_hazard(dust_met(gust = 20), crust = TRUE))
  })

  it("errors, naming the missing column, when a required column is absent", {
    skip_if_no_dust_v2()
    met <- dust_met()
    met$soil_moisture_0_to_1cm <- NULL
    expect_error(dust_hazard(met), regexp = "soil_moisture_0_to_1cm")
  })

  it("forwards site parameters to the engine", {
    skip_if_no_dust_v2()
    # A finer modal aggregate (higher sieve number, smaller d) lowers the
    # threshold, raising the index at fixed wind.
    coarse <- dust_hazard(dust_met(gust = 14), tyler_sieve_no = 20L)
    fine   <- dust_hazard(dust_met(gust = 14), tyler_sieve_no = 60L)
    expect_gte(fine, coarse)
  })
})


describe("dust units handling", {

  it("dust_flux accepts units-tagged winds and converts them to m/s", {
    skip_if_no_dust_v2()
    # 16 m/s = 57.6 km/h.
    bare   <- dust_flux(20L, 10, 0, 16, 0.02)
    tagged <- dust_flux(20L, 10, units::set_units(0, "m/s"),
                        units::set_units(57.6, "km/h"), 0.02)
    expect_equal(tagged, bare, tolerance = DUST_TOL)
  })

  it("dust_flux rejects a wind tagged with incompatible units", {
    skip_if_no_dust_v2()
    expect_error(
      dust_flux(20L, 10, 0, units::set_units(16, "degree_C"), 0.02),
      class = "meteoHazard_input_error"
    )
  })

  it("dust_hazard accepts a met tibble carrying units wind columns", {
    skip_if_no_dust_v2()
    met <- data.frame(soil_moisture_0_to_1cm = 0.02)
    met$wind_speed_10m <- units::set_units(0, "m/s")
    met$wind_gusts_10m <- units::set_units(72, "km/h")  # 20 m/s
    bare <- dust_hazard(data.frame(wind_speed_10m = 0, wind_gusts_10m = 20,
                                   soil_moisture_0_to_1cm = 0.02))
    expect_equal(dust_hazard(met), bare, tolerance = DUST_TOL)
  })

  it("dust_hazard returns a plain numeric relative flux", {
    skip_if_no_dust_v2()
    out <- dust_hazard(data.frame(wind_speed_10m = 0, wind_gusts_10m = 18,
                                  soil_moisture_0_to_1cm = 0.02))
    expect_false(inherits(out, "units"))
    expect_type(out, "double")
  })
})
