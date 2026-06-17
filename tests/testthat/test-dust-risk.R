# Behaviour specification for the v2 dust hazard model (specs/Dust.md, issue #5).
#
# v2 replaces the v1 AP-42 erosion-potential model with a physical
# saltation->dust framework. Specs skip until the v2 API is present (marker: the
# `gust_factor` formal on dust_emission_potential(), absent in v1).

skip_if_no_dust_v2 <- function() {
  testthat::skip_if_not(
    exists("dust_emission_potential", mode = "function") &&
      "gust_factor" %in% names(formals(dust_emission_potential)),
    "dust v2 API not yet implemented"
  )
}

DUST_TOL <- 1e-3

# A dry soil moisture (gravimetric w < w' for clay 10, bulk 1.6): w = sm/1.6*100,
# w' = 1.84, so sm < 0.0294 is "dry". 0.02 is safely dry.


describe("dust_emission_potential() [v2]", {

  it("returns zero flux below the entrainment threshold", {
    skip_if_no_dust_v2()
    # Light gust on dry coarse-sand surface: u* below u*t.
    f <- dust_emission_potential(20L, 10, wind_speed_10m = 5, wind_gusts_10m = 20,
                                 soil_moisture = 0.02)
    expect_equal(f, 0)
  })

  it("gives positive flux above threshold and increases monotonically with gust", {
    skip_if_no_dust_v2()
    f <- dust_emission_potential(20L, 10, wind_speed_10m = rep(0, 5),
                                 wind_gusts_10m = c(30, 45, 55, 70, 90),
                                 soil_moisture = rep(0.02, 5))
    expect_true(all(f >= 0))
    expect_true(all(diff(f) >= 0))
    expect_gt(f[5], 0)
  })

  it("is vectorised — one flux per input hour", {
    skip_if_no_dust_v2()
    f <- dust_emission_potential(20L, 10, c(0, 0, 0), c(40, 60, 80), c(0.02, 0.02, 0.02))
    expect_length(f, 3)
  })

  it("falls as soil moisture rises (Fécan threshold increase)", {
    skip_if_no_dust_v2()
    f <- dust_emission_potential(20L, 10, wind_speed_10m = rep(0, 3),
                                 wind_gusts_10m = rep(70, 3),
                                 soil_moisture = c(0.02, 0.06, 0.10))
    expect_true(all(diff(f) <= 0))
    expect_gt(f[1], f[3])
  })

  it("raising threshold_multiplier reduces the flux (crust hook)", {
    skip_if_no_dust_v2()
    f_base  <- dust_emission_potential(20L, 10, 0, 60, 0.02, threshold_multiplier = 1)
    f_crust <- dust_emission_potential(20L, 10, 0, 60, 0.02, threshold_multiplier = 2)
    expect_lt(f_crust, f_base)
  })

  it("drag partition (roughness_z0) lowers the flux versus the default no-partition", {
    skip_if_no_dust_v2()
    # roughness_z0 > z0s = d/30 raises the effective threshold (f_eff < 1).
    f_none <- dust_emission_potential(20L, 10, 0, 70, 0.02, roughness_z0 = NULL)
    f_rough <- dust_emission_potential(20L, 10, 0, 70, 0.02, roughness_z0 = 0.001)
    expect_lt(f_rough, f_none)
  })

  it("rejects an invalid Tyler sieve number", {
    skip_if_no_dust_v2()
    expect_error(dust_emission_potential(7L, 10, 0, 60, 0.02))
  })

  it("rejects out-of-range clay, soil moisture, and missing values", {
    skip_if_no_dust_v2()
    expect_error(dust_emission_potential(20L, 150, 0, 60, 0.02))   # clay > 100
    expect_error(dust_emission_potential(20L, 10, 0, 60, 1.5))     # soil moisture > 1
    expect_error(dust_emission_potential(20L, 10, c(0, NA), c(60, 60), c(0.02, 0.02)))
  })
})


describe("generate_dust_risk_index() [v2]", {

  dust_met <- function(gust = c(30, 50, 70), wind = 0, sm = 0.02, precip = NULL, n = NULL) {
    if (is.null(n)) n <- length(gust)
    df <- data.frame(
      wind_speed_10m         = rep_len(wind, n),
      wind_gusts_10m         = rep_len(gust, n),
      soil_moisture_0_to_1cm = rep_len(sm, n)
    )
    if (!is.null(precip)) df$precipitation <- rep_len(precip, n)
    df
  }

  it("returns one index per row, bounded to [0, 100]", {
    skip_if_no_dust_v2()
    out <- generate_dust_risk_index(dust_met(gust = c(20, 40, 60, 90)))
    expect_length(out, 4)
    expect_true(all(out >= 0 & out <= 100))
  })

  it("anchors the index to 100 at the reference gust on a dry surface", {
    skip_if_no_dust_v2()
    # By construction: an hour at scale_ref_gust (65), zero mean wind, dry,
    # crust off, equals the normalisation reference -> 100.
    out <- generate_dust_risk_index(dust_met(gust = 65, wind = 0, sm = 0.02))
    expect_equal(out, 100, tolerance = DUST_TOL)
  })

  it("is zero for sub-threshold winds and saturates at 100 above the reference", {
    skip_if_no_dust_v2()
    out <- generate_dust_risk_index(dust_met(gust = c(30, 80)))
    expect_equal(out[1], 0)
    expect_equal(out[2], 100, tolerance = DUST_TOL)
  })

  it("increases monotonically with gust", {
    skip_if_no_dust_v2()
    out <- generate_dust_risk_index(dust_met(gust = c(30, 40, 50, 65, 80)))
    expect_true(all(diff(out) >= 0))
  })

  it("with crust enabled, a rain event suppresses later dry hours that decay back", {
    skip_if_no_dust_v2()
    # Hour 1 rains (5 mm >= threshold), hours 2.. are dry but recently crusted.
    met <- dust_met(gust = 70, wind = 0, sm = 0.02,
                    precip = c(5, 0, 0, 0), n = 4)
    on  <- generate_dust_risk_index(met, crust = TRUE,
                                    crust_factor_max = 3, crust_decay_hours = 72)
    off <- generate_dust_risk_index(met, crust = FALSE)
    # Crust suppresses the post-rain hours relative to no-crust.
    expect_lt(on[2], off[2])
    # Recovery: suppression eases as time since rain grows.
    expect_lt(on[2], on[4])
  })

  it("ignores precipitation when crust is disabled", {
    skip_if_no_dust_v2()
    with_p    <- generate_dust_risk_index(dust_met(gust = 70, precip = 10), crust = FALSE)
    without_p <- generate_dust_risk_index(dust_met(gust = 70), crust = FALSE)
    expect_equal(with_p, without_p)
  })

  it("errors if crust is enabled but precipitation is absent", {
    skip_if_no_dust_v2()
    expect_error(generate_dust_risk_index(dust_met(gust = 70), crust = TRUE))
  })

  it("errors, naming the missing column, when a required column is absent", {
    skip_if_no_dust_v2()
    met <- dust_met()
    met$soil_moisture_0_to_1cm <- NULL
    expect_error(generate_dust_risk_index(met), regexp = "soil_moisture_0_to_1cm")
  })

  it("forwards site parameters to the engine", {
    skip_if_no_dust_v2()
    # A finer modal aggregate (higher sieve number, smaller d) lowers the
    # threshold, raising the index at fixed wind.
    coarse <- generate_dust_risk_index(dust_met(gust = 45), tyler_sieve_no = 20L)
    fine   <- generate_dust_risk_index(dust_met(gust = 45), tyler_sieve_no = 60L)
    expect_gte(fine, coarse)
  })
})
