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


describe("DUST_CONSTANTS [v3]", {

  it("exposes the pinned physical constants used by dust_flux()/dust_hazard() [CC-0a]", {
    skip_if_no_dust_v2()
    # Guards against a silent coefficient typo in the v3 constants refactor.
    expect_true(is.list(meteoHazard:::DUST_CONSTANTS))
    expect_equal(meteoHazard:::DUST_CONSTANTS$A_N, 0.0123)
    expect_equal(meteoHazard:::DUST_CONSTANTS$GAMMA, 3.0e-4)
    expect_equal(meteoHazard:::DUST_CONSTANTS$Z0_SMOOTH_RATIO, 1 / 30)
    expect_equal(meteoHazard:::DUST_CONSTANTS$MB95_CLAY_CAP, 20)
    expect_equal(meteoHazard:::DUST_CONSTANTS$FECAN_B, 0.68)
  })
})


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
    # gust bumped from 16 -> 20 (CONTRACT.md): 16 m/s sits below the v3 default
    # smooth-bed threshold (~17.9 m/s), which would make both fluxes 0.
    f_base  <- dust_flux(20L, 10, 0, 20, 0.02, threshold_multiplier = 1)
    f_crust <- dust_flux(20L, 10, 0, 20, 0.02, threshold_multiplier = 2)
    expect_lt(f_crust, f_base)
  })

  it("rejects an invalid Tyler sieve number", {
    skip_if_no_dust_v2()
    expect_error(dust_flux(7L, 10, 0, 16, 0.02),
                 class = "meteoHazard_input_error")
  })

  it("rejects z0 at or above the 10 m reference height (log-law singularity)", {
    skip_if_no_dust_v2()
    expect_error(dust_flux(20L, 10, 0, 20, 0.02, z0 = 10))
    expect_error(dust_flux(20L, 10, 0, 20, 0.02, z0 = 12))
  })

  it("rejects out-of-range clay, soil moisture, and missing values", {
    skip_if_no_dust_v2()
    expect_error(dust_flux(20L, 150, 0, 16, 0.02))   # clay > 100
    expect_error(dust_flux(20L, 10, 0, 16, 1.5))     # soil moisture > 1
    expect_error(dust_flux(20L, 10, c(0, NA), c(16, 16), c(0.02, 0.02)))
  })

  # ---- T1: smooth-bed roughness (z0 default NULL -> d/30) ------------------ #

  it("uses the smooth-bed roughness threshold by default (z0 = NULL) [CC-1a]", {
    skip_if_no_dust_v2()
    # Sieve 20 / clay 10 dry entrainment threshold at z0 = d/30 is ~17.9 m/s.
    expect_equal(
      dust_flux(20L, 10, wind_speed_10m = 0, wind_gusts_10m = 17, soil_moisture = 0.02),
      0
    )
    expect_gt(
      dust_flux(20L, 10, wind_speed_10m = 0, wind_gusts_10m = 18, soil_moisture = 0.02),
      0
    )
  })

  it("a z0 that fully shelters the bed under the MB95 drag partition returns all-zero flux and warns [CC-1b, v4 re-pin]", {
    skip_if_no_dust_v2()
    # z0 = 0.05 -> feff ~ -0.21 (negative, fully sheltered) at sieve 20/clay 10.
    expect_warning(
      dust_flux(20L, 10, 0, 20, 0.02, z0 = 0.05),
      regexp = "shelter"
    )
    expect_warning(
      dust_flux(20L, 10, 0, 20, 0.02, z0 = 0.05),
      class = "meteoHazard_dust_fully_sheltered"
    )
    expect_equal(suppressWarnings(dust_flux(20L, 10, c(0, 0), c(20, 40), c(0.02, 0.02), z0 = 0.05)),
                 c(0, 0))
  })

  it("z0 = NULL matches the explicit smooth-bed z0 exactly, without a warning [CC-1c]", {
    skip_if_no_dust_v2()
    z0_smooth <- meteoHazard:::TYLER_SIEVE_DIAMETERS_M[["20"]] * (1 / 30)
    f_default  <- dust_flux(20L, 10, 0, 20, 0.02, z0 = NULL)
    f_explicit <- dust_flux(20L, 10, 0, 20, 0.02, z0 = z0_smooth)
    expect_equal(f_default, f_explicit)
    # Guard is `>`, not `>=`: the exact smooth-bed value must not warn.
    expect_warning(dust_flux(20L, 10, 0, 20, 0.02, z0 = z0_smooth), NA)
  })

  # ---- T2: clamp MB95 sandblasting alpha at the clay validity ceiling ------ #

  it("caps the MB95 sandblasting alpha at the clay validity ceiling [CC-2a, v4 re-pin: z0 = NULL]", {
    skip_if_no_dust_v2()
    # Re-pinned under WP2 (MB95 drag partition): the old z0 = 0.005 fixture is
    # now sub-threshold at gust 20 (partition threshold ~65.8 m/s), so this
    # moves to the smooth bed (z0 = NULL), matching the WP1 oracle values.
    f_over_cap <- suppressWarnings(dust_flux(20L, 50, 0, 20, 0.02))
    f_at_cap   <- suppressWarnings(dust_flux(20L, 20, 0, 20, 0.02))
    expect_equal(f_over_cap, f_at_cap)
    expect_equal(f_at_cap, 7.86633315e-08 * 10^(0.134 * 10), tolerance = 1e-4)
    expect_equal(f_over_cap, 1.72096e-06, tolerance = 1e-4)
  })

  it("warns when clay_percent exceeds the MB95 validity ceiling [CC-2b]", {
    skip_if_no_dust_v2()
    # Use the default z0 = NULL so ONLY the clay warning fires (an explicit
    # z0 > smooth-bed would add a second, order-dependent warning). The alpha
    # clamp/warning is independent of the entrainment threshold.
    expect_warning(dust_flux(20L, 50, 0, 20, 0.02), regexp = "clay")
    expect_warning(dust_flux(20L, 50, 0, 20, 0.02),
                   class = "meteoHazard_dust_clay_capped")
  })

  # ---- T3: input validation (threshold_multiplier length; gust >= wind) ---- #

  it("rejects a threshold_multiplier whose length matches neither 1 nor n [CC-3a]", {
    skip_if_no_dust_v2()
    expect_error(
      dust_flux(20L, 10, c(0, 0, 0), c(20, 20, 20), c(0.02, 0.02, 0.02),
                threshold_multiplier = c(1, 2)),
      class = "meteoHazard_input_error"
    )
  })

  it("accepts length-1 and length-n threshold_multiplier without error [CC-3b]", {
    skip_if_no_dust_v2()
    expect_no_error(
      dust_flux(20L, 10, c(0, 0, 0), c(20, 20, 20), c(0.02, 0.02, 0.02),
                threshold_multiplier = 1)
    )
    expect_no_error(
      dust_flux(20L, 10, c(0, 0, 0), c(20, 20, 20), c(0.02, 0.02, 0.02),
                threshold_multiplier = c(1, 1, 2))
    )
  })

  it("rejects a gust below the mean wind speed [CC-3c]", {
    skip_if_no_dust_v2()
    expect_error(
      dust_flux(20L, 10, wind_speed_10m = 5, wind_gusts_10m = 3, soil_moisture = 0.02),
      regexp = "gust"
    )
    expect_error(
      dust_flux(20L, 10, wind_speed_10m = 5, wind_gusts_10m = 3, soil_moisture = 0.02),
      class = "meteoHazard_input_error"
    )
  })

  # ---- T5: known-answer values (KAT) ---------------------------------------- #
  # v4 re-pin (WP2): the pre-partition z0 = 0.005 fixtures assumed z0 raised u*
  # without sheltering the bed; under the MB95 drag partition the same z0
  # instead raises the threshold (~65.8 m/s at gust 20), making the old
  # fixtures sub-threshold and unreproducible by design. CC-5a moves to the
  # smooth bed; CC-5b moves to the z0 = 5e-4 partition threshold (~25.98 m/s).

  it("matches the hand-computed known-answer flux at the smooth bed (KAT) [CC-5a, v4 re-pin]", {
    skip_if_no_dust_v2()
    expect_equal(
      dust_flux(20L, 10, 0, 20, 0.02),
      7.86633315e-08,
      tolerance = 1e-4
    )
  })

  it("known-answer: the MB95 partition threshold at z0 = 5e-4 sits between gust 25.5 (zero) and 26.5 (positive) [CC-5b, v4 re-pin]", {
    skip_if_no_dust_v2()
    expect_equal(dust_flux(20L, 10, 0, 25.5, 0.02, z0 = 5e-4), 0)
    expect_gt(dust_flux(20L, 10, 0, 26.5, 0.02, z0 = 5e-4), 0)
  })

  # ---- WP2: MB95 drag partition -------------------------------------------- #

  it("feff = 1 (no drag-partition threshold effect) for any z0 <= the smooth-bed value", {
    skip_if_no_dust_v2()
    # At z0 == z0_smooth (or NULL) the flux is bit-identical to the smooth-bed
    # KAT. Below z0_smooth, feff is still 1 (no partition), but u* itself still
    # depends on the caller's smaller z0 through the log law -- so the overall
    # flux legitimately differs from the z0_smooth value; verify against a
    # feff = 1 (unpartitioned) reference computation instead.
    z0_smooth <- meteoHazard:::TYLER_SIEVE_DIAMETERS_M[["20"]] * (1 / 30)
    smooth  <- dust_flux(20L, 10, 0, 20, 0.02)
    at_val  <- dust_flux(20L, 10, 0, 20, 0.02, z0 = z0_smooth)
    expect_equal(at_val, smooth)

    z0_below <- z0_smooth / 2
    DC <- meteoHazard:::DUST_CONSTANTS
    d  <- meteoHazard:::TYLER_SIEVE_DIAMETERS_M[["20"]]
    u_star_t_dry <- sqrt(DC$A_N * (DC$RHO_P / DC$RHO_A_REF * DC$G * d + DC$GAMMA / (DC$RHO_A_REF * d)))
    u_star       <- DC$KAPPA * (0.84 * 20) / log(DC$Z_REF / z0_below)
    excess       <- 1 - (u_star_t_dry / u_star)^2
    alpha        <- 10^(DC$MB95_ALPHA_SLOPE * 10 + DC$MB95_ALPHA_INTERCEPT)
    expected     <- alpha * (DC$RHO_A_REF / DC$G) * u_star^3 * excess
    expect_equal(dust_flux(20L, 10, 0, 20, 0.02, z0 = z0_below), expected, tolerance = 1e-8)
  })

  it("the partition threshold gust at 10 m is strictly increasing across z0", {
    skip_if_no_dust_v2()
    # Reference table (sieve 20, clay 10, dry, verified during planning):
    # z0 = 2.78e-5 -> 17.90, 1e-4 -> 20.31, 5e-4 -> 25.98, 1e-3 -> 30.58,
    # 5e-3 -> 65.80 m/s. Bisect the zero/nonzero boundary for each.
    threshold_gust <- function(z0) {
      lo <- 0; hi <- 200
      for (i in 1:40) {
        mid <- (lo + hi) / 2
        if (dust_flux(20L, 10, 0, mid, 0.02, z0 = z0) > 0) hi <- mid else lo <- mid
      }
      (lo + hi) / 2
    }
    z0s <- c(2.78e-5, 1e-4, 5e-4, 1e-3, 5e-3)
    gusts <- vapply(z0s, threshold_gust, numeric(1))
    expect_true(all(diff(gusts) > 0))
    expect_equal(gusts, c(17.90, 20.31, 25.98, 30.58, 65.80), tolerance = 1e-2)
  })

})


describe("dust_flux() [v4 WP1: met-driven air density + d50 interface]", {

  it("met-driven air density: -5 degC / 1013.25 hPa raises the gust-20 flux by ~1.375x", {
    skip_if_no_dust_v2()
    base <- dust_flux(20L, 10, 0, 20, 0.02)
    cold <- dust_flux(20L, 10, 0, 20, 0.02,
                      temperature_2m = -5, surface_pressure = 1013.25)
    expect_equal(base, 7.86633315e-08, tolerance = 1e-4)   # new smooth-bed KAT
    expect_equal(cold, 1.08192029e-07, tolerance = 1e-4)
    expect_equal(cold / base, 1.37538, tolerance = 1e-3)
  })

  it("supplying only one of temperature_2m/surface_pressure errors (classed)", {
    skip_if_no_dust_v2()
    expect_error(dust_flux(20L, 10, 0, 20, 0.02, temperature_2m = -5),
                 class = "meteoHazard_input_error")
    expect_error(dust_flux(20L, 10, 0, 20, 0.02, surface_pressure = 1013.25),
                 class = "meteoHazard_input_error")
  })

  it("temperature_2m/surface_pressure vectorise per-hour (rho_a length-n through the chain)", {
    skip_if_no_dust_v2()
    f <- dust_flux(20L, 10, wind_speed_10m = rep(0, 2), wind_gusts_10m = rep(20, 2),
                   soil_moisture = rep(0.02, 2),
                   temperature_2m = c(15, -5), surface_pressure = c(1013.25, 1013.25))
    expect_length(f, 2)
    # Colder hour (index 2) is denser air -> larger flux at the same gust.
    expect_gt(f[2], f[1])
  })

  it("an implausible computed air density is a classed error (catches Pa-vs-hPa mistakes)", {
    skip_if_no_dust_v2()
    # surface_pressure supplied in Pa instead of hPa -> rho_a ~100x too high.
    expect_error(
      dust_flux(20L, 10, 0, 20, 0.02, temperature_2m = 15, surface_pressure = 101325),
      class = "meteoHazard_input_error"
    )
  })

  it("d50 equal to the sieve-20 opening reproduces the sieve-20 flux exactly", {
    skip_if_no_dust_v2()
    expect_equal(dust_flux(20L, 10, 0, 20, 0.02),
                 dust_flux(clay_percent = 10, wind_speed_10m = 0,
                           wind_gusts_10m = 20, soil_moisture = 0.02,
                           d50 = 0.000833))
  })

  it("rejects an out-of-range d50", {
    skip_if_no_dust_v2()
    expect_error(dust_flux(clay_percent = 10, wind_speed_10m = 0, wind_gusts_10m = 20,
                           soil_moisture = 0.02, d50 = 1e-6))
    expect_error(dust_flux(clay_percent = 10, wind_speed_10m = 0, wind_gusts_10m = 20,
                           soil_moisture = 0.02, d50 = 0.03))
  })

  it("requires exactly one of tyler_sieve_no or d50 (classed error when neither supplied)", {
    skip_if_no_dust_v2()
    expect_error(
      dust_flux(clay_percent = 10, wind_speed_10m = 0, wind_gusts_10m = 20, soil_moisture = 0.02),
      class = "meteoHazard_input_error"
    )
  })

  it("warns (classed) when both tyler_sieve_no and d50 are supplied; d50 supersedes", {
    skip_if_no_dust_v2()
    expect_warning(
      dust_flux(20L, 10, 0, 20, 0.02, d50 = 0.000701),   # sieve-24 diameter
      class = "meteoHazard_dust_d50_supersedes"
    )
    with_both <- suppressWarnings(dust_flux(20L, 10, 0, 20, 0.02, d50 = 0.000701))
    d50_only  <- dust_flux(clay_percent = 10, wind_speed_10m = 0, wind_gusts_10m = 20,
                           soil_moisture = 0.02, d50 = 0.000701)
    expect_equal(with_both, d50_only)
  })
})


describe("dust_flux() [v4 WP3: within-hour Weibull intermittency]", {

  it("weibull forcing matches the closed-form oracle (mean 12, k = 2, dry, smooth bed)", {
    skip_if_no_dust_v2()
    f <- dust_flux(20L, 10, wind_speed_10m = 12, wind_gusts_10m = 20,
                   soil_moisture = 0.02, forcing = "weibull")
    expect_equal(f, 1.062634e-07, tolerance = 1e-4)
  })

  it("weibull forcing emits where the steady mean is sub-threshold (intermittency tail)", {
    skip_if_no_dust_v2()
    steady <- dust_flux(20L, 10, 12, 12 / 0.84, 0.02)          # U_fm = 12 -> 0
    weib   <- dust_flux(20L, 10, 12, 12 / 0.84, 0.02, forcing = "weibull")
    expect_equal(steady, 0)
    expect_gt(weib, 0)
  })

  it("a near-degenerate Weibull (k large) converges to the steady flux at the mean", {
    skip_if_no_dust_v2()
    # k = 200: c ~ mean (gamma(1 + 1/k) -> 1), so weibull mode should converge
    # to the deterministic gust-mode flux with U_fm forced exactly equal to
    # the mean (wind_speed_10m = wind_gusts_10m * gust_factor = mean).
    mean_u <- 20
    steady <- dust_flux(20L, 10, wind_speed_10m = mean_u,
                        wind_gusts_10m = mean_u / 0.84, soil_moisture = 0.02,
                        forcing = "gust")
    weib_k200 <- dust_flux(20L, 10, wind_speed_10m = mean_u,
                           wind_gusts_10m = mean_u / 0.84, soil_moisture = 0.02,
                           forcing = "weibull", weibull_shape = 200)
    expect_equal(weib_k200, steady, tolerance = 1e-3)
  })

  it("forcing = 'gust' (default) is bit-identical to the pre-WP3 gust path", {
    skip_if_no_dust_v2()
    explicit <- dust_flux(20L, 10, 0, 20, 0.02, forcing = "gust")
    default  <- dust_flux(20L, 10, 0, 20, 0.02)
    expect_equal(explicit, default)
    expect_equal(default, 7.86633315e-08, tolerance = 1e-4)
  })

  it("weibull forcing vectorises per-hour (one mean per hour)", {
    skip_if_no_dust_v2()
    f <- dust_flux(20L, 10, wind_speed_10m = c(5, 12, 20), wind_gusts_10m = rep(30, 3),
                   soil_moisture = rep(0.02, 3), forcing = "weibull")
    expect_length(f, 3)
    expect_true(all(diff(f) >= 0))   # monotone in the mean
  })

  it("rejects a non-positive weibull_shape", {
    skip_if_no_dust_v2()
    expect_error(
      dust_flux(20L, 10, 0, 20, 0.02, forcing = "weibull", weibull_shape = 0)
    )
    expect_error(
      dust_flux(20L, 10, 0, 20, 0.02, forcing = "weibull", weibull_shape = -1)
    )
  })

  it("rejects an invalid forcing value", {
    skip_if_no_dust_v2()
    expect_error(dust_flux(20L, 10, 0, 20, 0.02, forcing = "bogus"))
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

  it("hours_since_last_rain = 0 suppresses hour 1 relative to Inf (crust cold-start) [CC-4a]", {
    skip_if_no_dust_v2()
    # No rain at all in the series: the only difference between the two calls
    # is the assumed crust age going into hour 1.
    met  <- dust_met(gust = 20, wind = 0, sm = 0.02, precip = 0, n = 5)
    cold <- dust_hazard(met, crust = TRUE, hours_since_last_rain = 0)
    warm <- dust_hazard(met, crust = TRUE, hours_since_last_rain = Inf)
    expect_lt(cold[1], warm[1])
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

  # ---- WP1: air_density = "met" ------------------------------------------- #

  it("air_density = 'met' requires temperature_2m/surface_pressure columns", {
    skip_if_no_dust_v2()
    expect_error(
      dust_hazard(dust_met(gust = 20), air_density = "met"),
      regexp = "temperature_2m|surface_pressure"
    )
  })

  it("air_density = 'reference' (default) ignores temperature_2m/surface_pressure columns even if present", {
    skip_if_no_dust_v2()
    met <- dust_met(gust = 20)
    met$temperature_2m   <- -5
    met$surface_pressure <- 1013.25
    with_cols    <- dust_hazard(met, air_density = "reference")
    without_cols <- dust_hazard(dust_met(gust = 20))
    expect_equal(with_cols, without_cols)
  })

  it("air_density = 'met' forwards temperature_2m/surface_pressure and raises the flux for cold dense air", {
    skip_if_no_dust_v2()
    met <- dust_met(gust = 20)
    met$temperature_2m   <- -5
    met$surface_pressure <- 1013.25
    met_out <- dust_hazard(met, air_density = "met")
    ref_out <- dust_hazard(dust_met(gust = 20), air_density = "reference")
    expect_true(all(met_out >= ref_out))
    expect_true(any(met_out > ref_out))
  })

  it("dust_hazard() default path (air_density = 'reference') is bit-identical to pre-WP1 behaviour", {
    skip_if_no_dust_v2()
    out <- dust_hazard(dust_met(gust = c(8, 14, 18, 25)))
    ref <- dust_flux(20L, 10, wind_speed_10m = rep(0, 4), wind_gusts_10m = c(8, 14, 18, 25),
                     soil_moisture = rep(0.02, 4))
    expect_equal(out, ref)
  })

  # ---- WP3: forcing = "weibull" -------------------------------------------- #

  it("forwards forcing/weibull_shape to the engine", {
    skip_if_no_dust_v2()
    met <- dust_met(gust = 20, wind = 12, sm = 0.02)
    steady  <- dust_hazard(met, forcing = "gust")
    weibull <- dust_hazard(met, forcing = "weibull")
    ref <- dust_flux(20L, 10, wind_speed_10m = 12, wind_gusts_10m = 20,
                     soil_moisture = 0.02, forcing = "weibull")
    expect_equal(weibull, ref, tolerance = DUST_TOL)
    expect_false(isTRUE(all.equal(weibull, steady)))
  })

  it("default forcing = 'gust' is bit-identical to pre-WP3 dust_hazard() output", {
    skip_if_no_dust_v2()
    out <- dust_hazard(dust_met(gust = c(8, 14, 18, 25)))
    explicit <- dust_hazard(dust_met(gust = c(8, 14, 18, 25)), forcing = "gust")
    expect_equal(out, explicit)
  })
})


describe(".dust_u_star() / .dust_u_star_t_dry() [v4 WP4 refactor]", {

  it(".dust_u_star_t_dry() reproduces the dry threshold used inside dust_flux()", {
    skip_if_no_dust_v2()
    DC <- meteoHazard:::DUST_CONSTANTS
    d  <- meteoHazard:::TYLER_SIEVE_DIAMETERS_M[["20"]]
    expected <- sqrt(DC$A_N * (DC$RHO_P / DC$RHO_A_REF * DC$G * d + DC$GAMMA / (DC$RHO_A_REF * d)))
    expect_equal(meteoHazard:::.dust_u_star_t_dry(d, DC$RHO_A_REF), expected)
  })

  it(".dust_u_star() reproduces the gust-driven friction velocity used inside dust_flux()", {
    skip_if_no_dust_v2()
    DC <- meteoHazard:::DUST_CONSTANTS
    d  <- meteoHazard:::TYLER_SIEVE_DIAMETERS_M[["20"]]
    z0 <- d * DC$Z0_SMOOTH_RATIO
    expected <- (DC$KAPPA / log(DC$Z_REF / z0)) * max(0, 0.84 * 20)
    expect_equal(meteoHazard:::.dust_u_star(0, 20, 0.84, z0), expected)
  })

  it("dust_flux() is bit-identical after the WP4 refactor (smooth-bed KAT)", {
    skip_if_no_dust_v2()
    expect_equal(dust_flux(20L, 10, 0, 20, 0.02), 7.86633315e-08, tolerance = 1e-4)
  })
})


describe(".dust_crust_factor() [v3 cold-start]", {

  it("age0 sets the initial crust age for the cold-start row [CC-4c]", {
    skip_if_no_dust_v2()
    # No rain in the series (all precip 0): with age0 = 0, hour 1 is treated as
    # freshly rained-on (age 0 -> full factor_max). With age0 = Inf, hour 1 has
    # no crust memory at all (age Inf -> factor 1, i.e. old default behaviour).
    fresh <- meteoHazard:::.dust_crust_factor(rep(0, 3), 2, 3, 24, age0 = 0)
    stale <- meteoHazard:::.dust_crust_factor(rep(0, 3), 2, 3, 24, age0 = Inf)
    expect_equal(fresh[1], 3, tolerance = 1e-6)
    expect_equal(stale[1], 1)
  })
})


describe(".dust_crust_factor_saltation() [v4 WP4]", {

  it("holds the crust age (and multiplier) through calm/sub-threshold hours", {
    skip_if_no_dust_v2()
    precip <- rep(0, 5)
    u_star <- rep(0.1, 5)          # always below the crusted threshold
    u_star_t_moist <- rep(1, 5)
    mult <- meteoHazard:::.dust_crust_factor_saltation(
      precip, u_star, u_star_t_moist, threshold = 2, factor_max = 3,
      decay_hours = 24, age0 = 0
    )
    expect_equal(mult, rep(3, 5))   # age stays 0 throughout -> full factor_max
  })

  it("advances the crust age only on hours where u_star exceeds the crusted threshold", {
    skip_if_no_dust_v2()
    precip <- rep(0, 5)
    u_star <- rep(10, 5)           # always exceeds the threshold, even fully crusted
    u_star_t_moist <- rep(1, 5)
    mult <- meteoHazard:::.dust_crust_factor_saltation(
      precip, u_star, u_star_t_moist, threshold = 2, factor_max = 3,
      decay_hours = 24, age0 = 0
    )
    expect_true(all(diff(mult) < 0))   # decays every hour once saltation occurs
  })

  it("a rain hour resets the age for the FOLLOWING hour, not retroactively", {
    skip_if_no_dust_v2()
    precip <- c(0, 5, 0, 0)
    u_star <- rep(0, 4)             # never exceeds the threshold
    u_star_t_moist <- rep(1, 4)
    mult <- meteoHazard:::.dust_crust_factor_saltation(
      precip, u_star, u_star_t_moist, threshold = 2, factor_max = 3,
      decay_hours = 24, age0 = Inf
    )
    expect_equal(mult[1], 1)                     # age0 = Inf -> no crust memory yet
    expect_equal(mult[2], 1)                     # this hour's rain not yet reflected
    expect_equal(mult[3], 3, tolerance = 1e-6)    # next hour: fresh crust
    expect_equal(mult[4], 3, tolerance = 1e-6)    # held (calm, sub-threshold)
  })
})


describe("dust_hazard() [v4 WP4: crust_decay = 'saltation']", {

  it("saltation mode: crust persists through a calm week, clock mode decays it", {
    skip_if_no_dust_v2()
    # rain hour 1, then 167 dead-calm hours, then a strong hour 169
    met <- data.frame(wind_speed_10m = c(0, rep(0, 167), 0),
                      wind_gusts_10m = c(0, rep(2, 167), 25),
                      soil_moisture_0_to_1cm = 0.02,
                      precipitation = c(5, rep(0, 168)))
    clock <- dust_hazard(met, crust = TRUE, crust_decay_hours = 24)
    salt  <- dust_hazard(met, crust = TRUE, crust_decay_hours = 24,
                         crust_decay = "saltation")
    expect_lt(salt[169], clock[169])   # crust still fully fresh in saltation mode
    expect_equal(salt[169], 0)         # still-elevated crust threshold -> zero flux
  })

  it("saltation mode: sustained supra-threshold wind decays the crust", {
    skip_if_no_dust_v2()
    n <- 200
    met <- data.frame(wind_speed_10m = rep(0, n),
                      wind_gusts_10m = c(0, rep(60, n - 1)),   # exceeds even the
                      soil_moisture_0_to_1cm = 0.02,           # fully-crusted threshold
                      precipitation = c(5, rep(0, n - 1)))
    base <- dust_hazard(met, crust = FALSE)
    salt <- dust_hazard(met, crust = TRUE, crust_decay_hours = 24,
                        crust_decay = "saltation")
    expect_lt(salt[10], base[10])                       # still suppressed early on
    expect_equal(salt[n], base[n], tolerance = 1e-3)     # decayed back to baseline
  })

  it("crust_decay = 'clock' is bit-identical to the current (pre-WP4) behaviour", {
    skip_if_no_dust_v2()
    met <- data.frame(wind_speed_10m = 0, wind_gusts_10m = c(20, 22, 25),
                      soil_moisture_0_to_1cm = 0.02, precipitation = c(5, 0, 0))
    implicit <- dust_hazard(met, crust = TRUE, crust_decay_hours = 24)
    explicit <- dust_hazard(met, crust = TRUE, crust_decay_hours = 24, crust_decay = "clock")
    expect_equal(implicit, explicit)
  })

  it("crust_decay is irrelevant when crust = FALSE", {
    skip_if_no_dust_v2()
    met <- data.frame(wind_speed_10m = 0, wind_gusts_10m = c(20, 22, 25),
                      soil_moisture_0_to_1cm = 0.02)
    clock_off <- dust_hazard(met, crust = FALSE, crust_decay = "clock")
    salt_off  <- dust_hazard(met, crust = FALSE, crust_decay = "saltation")
    expect_equal(clock_off, salt_off)
  })

  it("rejects an invalid crust_decay value", {
    skip_if_no_dust_v2()
    met <- data.frame(wind_speed_10m = 0, wind_gusts_10m = 20, soil_moisture_0_to_1cm = 0.02)
    expect_error(dust_hazard(met, crust_decay = "bogus"))
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
