# Behaviour specification for the v3 Litter Hazard Index.
#
# These specs encode the behaviour defined in specs/Litter_v3.md (the approved
# v3 design). v3 is NOT yet implemented — R/litter_risk.R currently implements
# v2 — so every spec below is written against the v3 contract and SKIPS until
# the v3 API is present. This is the BDD "red" checklist for the v3
# reimplementation: the implement step turns these green.
#
# Detection: the v3 `litter_risk_index()` introduces a roughness-length argument
# `z0` (friction-velocity conversion) that v2 does not have. Presence of that
# formal is used as the "v3 is implemented" signal.
#
# Proposed v3 signature these specs assume (names may be adjusted by the
# implementer, but the *behaviour* and *default values* are the contract):
#
#   litter_risk_index(
#     wind_gusts_10m, wind_speed_10m, precipitation, soil_moisture_0_to_1cm,
#     kappa = 0.40, z0 = 0.05,
#     ustar_t0 = 0.30, ustar_ref = 1.05, entrainment_max = 50, excess_exponent = 2,
#     moisture_gain = 2.0, moisture_curve = 0.5, soil_dry = 0.05, soil_wet = 0.20,
#     wind_transport_onset = 20, wind_transport_ref = 55, transport_max = 2.0,
#     rain_threshold = 0.5
#   )
#
# Friction-velocity conversion used by the worked-answer oracles (spec §3.1):
#   u* = (kappa / ln(10 / z0)) * W   with W in m/s
#   ln(10 / 0.05) = ln(200) = 5.2983174
#   per-km/h factor = 0.40 / 5.2983174 / 3.6 = 0.02097219

litter_v3_available <- function() {
  exists("litter_risk_index", mode = "function") &&
    "z0" %in% names(formals(litter_risk_index))
}

skip_if_no_litter_v3 <- function() {
  testthat::skip_if_not(litter_v3_available(), "litter v3 API not yet implemented")
}

# Tolerance for known-answer comparisons: the index is O(100) and computed with
# a handful of double-precision operations, so 1e-3 is far tighter than any
# operational meaning while leaving ample headroom over floating-point error.
LRI_TOL <- 1e-3


describe("litter_risk_index() [v3]", {

  describe("output contract and bounds", {

    it("returns one numeric value per input hour", {
      skip_if_no_litter_v3()
      out <- litter_risk_index(
        wind_gusts_10m         = c(10, 35, 55),
        wind_speed_10m         = c(5, 25, 45),
        precipitation          = c(0, 0, 0),
        soil_moisture_0_to_1cm = c(0.02, 0.10, 0.02)
      )
      expect_type(out, "double")
      expect_length(out, 3)
    })

    it("never returns a value outside [0, 100]", {
      skip_if_no_litter_v3()
      # Sweep a wide grid; the bound must hold everywhere.
      g <- expand.grid(
        gust = c(0, 15, 40, 80, 150),
        wind = c(0, 20, 60, 120),
        rain = c(0, 1),
        sm   = c(0, 0.05, 0.12, 0.30)
      )
      out <- litter_risk_index(g$gust, g$wind, g$rain, g$sm)
      expect_true(all(out >= 0 & out <= 100))
    })

    it("reaches exactly 100 only with saturated entrainment and saturating wind", {
      skip_if_no_litter_v3()
      # E = entrainment_max (50) and T = transport_max (2) -> 50 * 2 = 100 (spec §5.1).
      out <- litter_risk_index(
        wind_gusts_10m = 60, wind_speed_10m = 55,
        precipitation = 0, soil_moisture_0_to_1cm = 0.02
      )
      expect_equal(out, 100, tolerance = LRI_TOL)
    })

    it("is direction-agnostic: takes no wind direction, CAPE, is_day, or receptor arcs", {
      skip_if_no_litter_v3()
      # v3 changes 5 & 6: these v2 inputs are removed from the hazard index.
      fmls <- names(formals(litter_risk_index))
      expect_false("wind_direction_10m" %in% fmls)
      expect_false("cape" %in% fmls)
      expect_false("is_day" %in% fmls)
      expect_false("receptor_arcs" %in% fmls)
    })

    it("requires only the four v3 meteorological inputs", {
      skip_if_no_litter_v3()
      fmls <- names(formals(litter_risk_index))
      expect_true(all(
        c("wind_gusts_10m", "wind_speed_10m", "precipitation",
          "soil_moisture_0_to_1cm") %in% fmls
      ))
    })
  })

  describe("entrainment (gust -> friction velocity)", {

    it("is zero when the gust friction velocity is below the dry threshold", {
      skip_if_no_litter_v3()
      # u*_g(14 km/h) = 0.02097219 * 14 = 0.2936 < ustar_t0 = 0.30 -> E = 0.
      out <- litter_risk_index(
        wind_gusts_10m = c(0, 10, 14), wind_speed_10m = c(45, 45, 45),
        precipitation = c(0, 0, 0), soil_moisture_0_to_1cm = c(0.02, 0.02, 0.02)
      )
      expect_equal(out, c(0, 0, 0))
    })

    it("increases monotonically with gust above the threshold", {
      skip_if_no_litter_v3()
      out <- litter_risk_index(
        wind_gusts_10m = c(10, 20, 30, 40, 55, 70),
        wind_speed_10m = rep(45, 6),
        precipitation = rep(0, 6),
        soil_moisture_0_to_1cm = rep(0.02, 6)
      )
      expect_true(all(diff(out) >= 0))
    })

    it("saturates above the entrainment reference gust", {
      skip_if_no_litter_v3()
      # Once u*_g exceeds ustar_ref the entrainment caps at entrainment_max, so
      # further gust gives no increase (dry surface, identical wind).
      out <- litter_risk_index(
        wind_gusts_10m = c(55, 70, 120), wind_speed_10m = rep(45, 3),
        precipitation = rep(0, 3), soil_moisture_0_to_1cm = rep(0.02, 3)
      )
      expect_equal(out[1], out[2], tolerance = LRI_TOL)
      expect_equal(out[2], out[3], tolerance = LRI_TOL)
    })

    it("follows the excess-squared shape (convex ramp between threshold and reference)", {
      skip_if_no_litter_v3()
      # (u*_g - u*t)^2 is convex, so equal gust steps in the ramp produce
      # increasing increments. Dry surface, calm wind (T = 1) isolates E.
      out <- litter_risk_index(
        wind_gusts_10m = c(20, 30, 40, 50), wind_speed_10m = rep(5, 4),
        precipitation = rep(0, 4), soil_moisture_0_to_1cm = rep(0.02, 4)
      )
      incr <- diff(out)
      expect_true(all(diff(incr) > 0))
    })
  })

  describe("moisture-raised entrainment threshold (Alt-1b)", {

    it("lowers the hazard as the surface dampens, at fixed gust", {
      skip_if_no_litter_v3()
      # Rising soil moisture (below wet) raises u*t, shrinking the excess.
      out <- litter_risk_index(
        wind_gusts_10m = rep(40, 4), wind_speed_10m = rep(30, 4),
        precipitation = rep(0, 4),
        soil_moisture_0_to_1cm = c(0.05, 0.10, 0.15, 0.19)
      )
      expect_true(all(diff(out) <= 0))
    })

    it("treats all moisture at or below the dry threshold identically", {
      skip_if_no_litter_v3()
      # Below soil_dry the threshold is not raised, so SM = 0.00, 0.02, 0.05 match.
      out <- litter_risk_index(
        wind_gusts_10m = rep(40, 3), wind_speed_10m = rep(30, 3),
        precipitation = rep(0, 3), soil_moisture_0_to_1cm = c(0.00, 0.02, 0.05)
      )
      expect_equal(out[1], out[2], tolerance = LRI_TOL)
      expect_equal(out[2], out[3], tolerance = LRI_TOL)
    })

    it("vetoes entirely when the surface is saturated (SM >= soil_wet)", {
      skip_if_no_litter_v3()
      # Saturation veto: hazard is exactly 0 regardless of gust strength.
      out <- litter_risk_index(
        wind_gusts_10m = c(55, 120), wind_speed_10m = c(45, 55),
        precipitation = c(0, 0), soil_moisture_0_to_1cm = c(0.20, 0.30)
      )
      expect_equal(out, c(0, 0))
    })

    it("a damp surface needs a stronger gust than a dry surface to entrain", {
      skip_if_no_litter_v3()
      # A gust that produces hazard on a dry face produces less (here ~zero) on a
      # damp face because the threshold has risen above it.
      dry  <- litter_risk_index(35, 5, 0, 0.02)
      damp <- litter_risk_index(35, 5, 0, 0.14)
      expect_gt(dry, damp)
    })
  })

  describe("rain hard gate", {

    it("forces the hazard to zero when precipitation meets the threshold", {
      skip_if_no_litter_v3()
      # Even extreme dry-surface wind is fully suppressed by rain >= rain_threshold.
      out <- litter_risk_index(
        wind_gusts_10m = c(60, 60), wind_speed_10m = c(55, 55),
        precipitation = c(0.5, 3.0), soil_moisture_0_to_1cm = c(0.02, 0.02)
      )
      expect_equal(out, c(0, 0))
    })

    it("does not gate when precipitation is below the threshold", {
      skip_if_no_litter_v3()
      out <- litter_risk_index(
        wind_gusts_10m = 60, wind_speed_10m = 55,
        precipitation = 0.4, soil_moisture_0_to_1cm = 0.02
      )
      expect_gt(out, 0)
    })
  })

  describe("transport multiplier (mean wind)", {

    it("applies no amplification at or below the transport onset wind", {
      skip_if_no_litter_v3()
      # W <= wind_transport_onset (20 km/h) -> T = 1, so the result equals the
      # entrainment-only value and calm vs onset winds match.
      out <- litter_risk_index(
        wind_gusts_10m = c(55, 55), wind_speed_10m = c(5, 20),
        precipitation = c(0, 0), soil_moisture_0_to_1cm = c(0.02, 0.02)
      )
      expect_equal(out[1], out[2], tolerance = LRI_TOL)
      # With saturated entrainment (E = 50) and T = 1 the hazard is exactly 50.
      expect_equal(out[1], 50, tolerance = LRI_TOL)
    })

    it("increases the hazard with mean wind above the onset (the distance penalty)", {
      skip_if_no_litter_v3()
      out <- litter_risk_index(
        wind_gusts_10m = rep(55, 6),
        wind_speed_10m = c(0, 20, 30, 40, 55, 70),
        precipitation = rep(0, 6),
        soil_moisture_0_to_1cm = rep(0.02, 6)
      )
      expect_true(all(diff(out) >= 0))
    })

    it("saturates the transport amplification at the reference wind", {
      skip_if_no_litter_v3()
      # At/above wind_transport_ref (55 km/h) T = transport_max; more wind adds nothing.
      out <- litter_risk_index(
        wind_gusts_10m = rep(55, 3), wind_speed_10m = c(55, 70, 120),
        precipitation = rep(0, 3), soil_moisture_0_to_1cm = rep(0.02, 3)
      )
      expect_equal(out[1], out[2], tolerance = LRI_TOL)
      expect_equal(out[2], out[3], tolerance = LRI_TOL)
    })

    it("never reduces the hazard below the entrainment-only value", {
      skip_if_no_litter_v3()
      # T >= 1 always, so any wind can only amplify, never suppress.
      calm   <- litter_risk_index(40, 0, 0, 0.02)
      windy  <- litter_risk_index(40, 60, 0, 0.02)
      expect_gte(windy, calm)
    })
  })

  describe("multiplicative veto property", {

    it("returns zero if any single suppressing condition holds", {
      skip_if_no_litter_v3()
      # Each row has strong gust+wind but one veto active: rain / saturated /
      # sub-threshold gust. All must be exactly zero.
      out <- litter_risk_index(
        wind_gusts_10m         = c(60, 60, 10),
        wind_speed_10m         = c(55, 55, 55),
        precipitation          = c(2.0, 0.0, 0.0),
        soil_moisture_0_to_1cm = c(0.02, 0.25, 0.02)
      )
      expect_equal(out, c(0, 0, 0))
    })
  })

  describe("worked examples from spec section 5.2", {

    it("Example 1 — hot, dry, strong gust -> 85.7 (EXTREME)", {
      skip_if_no_litter_v3()
      # G=55, W=45, P=0, SM=0.02 (dry). E saturates at 50; T(45)=1.714286.
      # LRI = 50 * 1.714286 = 85.7143.
      out <- litter_risk_index(55, 45, 0, 0.02)
      expect_equal(out, 85.7143, tolerance = LRI_TOL)
    })

    it("Example 2 — moderate gust on a damp surface -> ~0.78 (LOW)", {
      skip_if_no_litter_v3()
      # G=35, W=25, P=0, SM=0.10. u*t rises to 0.6464; E=0.68237; T(25)=1.142857.
      # LRI = 0.779847.
      out <- litter_risk_index(35, 25, 0, 0.10)
      expect_equal(out, 0.7798, tolerance = LRI_TOL)
    })

    it("Example 3 — rain now -> 0 (LOW)", {
      skip_if_no_litter_v3()
      out <- litter_risk_index(55, 45, 3.0, 0.02)
      expect_equal(out, 0)
    })
  })

  describe("input validation", {

    it("rejects missing values in any input (no NA imputation; spec section 1.1)", {
      skip_if_no_litter_v3()
      expect_error(litter_risk_index(c(40, NA), c(30, 30), c(0, 0), c(0.1, 0.1)))
      expect_error(litter_risk_index(c(40, 40), c(30, NA), c(0, 0), c(0.1, 0.1)))
      expect_error(litter_risk_index(c(40, 40), c(30, 30), c(0, NA), c(0.1, 0.1)))
      expect_error(litter_risk_index(c(40, 40), c(30, 30), c(0, 0), c(0.1, NA)))
    })

    it("rejects inputs of differing length", {
      skip_if_no_litter_v3()
      expect_error(litter_risk_index(c(40, 40), 30, 0, 0.1))
    })

    it("rejects zero-length input", {
      skip_if_no_litter_v3()
      expect_error(litter_risk_index(numeric(0), numeric(0), numeric(0), numeric(0)))
    })

    it("rejects physically impossible meteorology", {
      skip_if_no_litter_v3()
      expect_error(litter_risk_index(-1, 30, 0, 0.1))    # negative gust
      expect_error(litter_risk_index(40, -1, 0, 0.1))    # negative wind
      expect_error(litter_risk_index(40, 30, -1, 0.1))   # negative precip
      expect_error(litter_risk_index(40, 30, 0, 1.5))    # soil moisture > 1
      expect_error(litter_risk_index(40, 30, 0, -0.1))   # soil moisture < 0
    })

    it("rejects parameter sets that violate ordering constraints", {
      skip_if_no_litter_v3()
      # Each denominator/threshold pair must be strictly ordered.
      expect_error(litter_risk_index(40, 30, 0, 0.1, soil_dry = 0.20, soil_wet = 0.10))
      expect_error(litter_risk_index(40, 30, 0, 0.1, ustar_t0 = 1.2, ustar_ref = 1.0))
      expect_error(litter_risk_index(
        40, 30, 0, 0.1, wind_transport_onset = 60, wind_transport_ref = 55
      ))
    })
  })
})


describe("generate_litter_risk_index() [v3]", {

  v3_met <- function(n = 3) {
    data.frame(
      wind_gusts_10m         = seq(10, 60, length.out = n),
      wind_speed_10m         = seq(5, 55, length.out = n),
      precipitation          = rep(0, n),
      soil_moisture_0_to_1cm = rep(0.02, n)
    )
  }

  it("computes one bounded index per row of the met tibble", {
    skip_if_no_litter_v3()
    out <- generate_litter_risk_index(v3_met(4))
    expect_type(out, "double")
    expect_length(out, 4)
    expect_true(all(out >= 0 & out <= 100))
  })

  it("agrees with the vector API called on the same columns", {
    skip_if_no_litter_v3()
    met <- v3_met(5)
    expect_equal(
      generate_litter_risk_index(met),
      litter_risk_index(
        met$wind_gusts_10m, met$wind_speed_10m,
        met$precipitation, met$soil_moisture_0_to_1cm
      )
    )
  })

  it("does not require the removed v2 columns (cape, is_day, wind_direction)", {
    skip_if_no_litter_v3()
    # A tibble carrying only the four v3 columns must be accepted.
    expect_no_error(generate_litter_risk_index(v3_met()))
  })

  it("errors, naming the missing column, when a required column is absent", {
    skip_if_no_litter_v3()
    met <- v3_met()
    met$soil_moisture_0_to_1cm <- NULL
    expect_error(generate_litter_risk_index(met), regexp = "soil_moisture_0_to_1cm")
  })

  it("forwards calibration parameters through to the vector API", {
    skip_if_no_litter_v3()
    # Raising the rain threshold lets a light-rain hour through that the default
    # (0.5 mm) would gate; the two calls must differ for such an hour.
    met <- data.frame(
      wind_gusts_10m = 60, wind_speed_10m = 55,
      precipitation = 0.6, soil_moisture_0_to_1cm = 0.02
    )
    expect_equal(generate_litter_risk_index(met), 0)            # gated by default
    expect_gt(generate_litter_risk_index(met, rain_threshold = 1.0), 0)
  })
})
