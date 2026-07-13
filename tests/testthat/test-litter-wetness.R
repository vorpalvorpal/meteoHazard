# Behaviour specification for the litter-wetness state subsystem.
#
# `litter_wetness_vec()` is a sequential (hourly) surface-wetness model for
# windblown litter: an hour with enough rain resets the surface to fully wet
# (w=1); otherwise the surface dries exponentially at a rate driven by
# vapour-pressure deficit (VPD), wind, and insolation. `litter_wetness()` is
# the met-tibble wrapper (mirrors `litter_hazard()`).
#
# Pinned formulas:
#   .litter_vpd(T, RH): es = 0.6108*exp(17.27*T/(T+237.3)); vpd = pmax(0, es*(1-RH/100))
#   dry_rate = dry_rate_base*(1+vpd_coef*vpd)*(1+wind_coef*U)*(1+sw_coef*SW/sw_ref)
#   w[1] = if precip[1] >= wetness_set_precip 1 else w0*exp(-dry_rate[1])
#   w[i] = if precip[i] >= wetness_set_precip 1 else w[i-1]*exp(-dry_rate[i])
#
# Worked oracle (T=25 degC, RH=40%, U=5 m/s, SW=600 W/m^2, w0=1, no rain):
#   vpd = 0.6108*exp(17.27*25/262.3)*(1-0.40) = 1.900667 kPa
#   dry_rate = 0.7*(1+1.900667)*(1+0.5)*(1+0.6) = 4.87312 /h
#   w = 1*exp(-4.87312) = 0.007649462

WET_TOL <- 1e-4


describe("litter_wetness_vec()", {

  describe("rain reset", {

    it("an hour with precipitation >= wetness_set_precip sets w == 1, regardless of w0", {
      out_dry_start <- litter_wetness_vec(
        precipitation = 0.5, temperature_2m = 25, relative_humidity_2m = 40,
        wind_speed_10m = 5, shortwave_radiation = 600, w0 = 0
      )
      out_wet_start <- litter_wetness_vec(
        precipitation = 0.5, temperature_2m = 25, relative_humidity_2m = 40,
        wind_speed_10m = 5, shortwave_radiation = 600, w0 = 1
      )
      expect_equal(out_dry_start, 1)
      expect_equal(out_wet_start, 1)
    })

    it("precipitation strictly below wetness_set_precip does not force a reset", {
      out <- litter_wetness_vec(
        precipitation = 0.49, temperature_2m = 25, relative_humidity_2m = 40,
        wind_speed_10m = 5, shortwave_radiation = 600, w0 = 1
      )
      expect_lt(out, 1)
    })

    it("uses the previous hour's wetness as sequential state (a prior rain hour wets a later dry hour)", {
      no_rain_before <- litter_wetness_vec(
        precipitation = c(0, 0), temperature_2m = c(20, 20),
        relative_humidity_2m = c(50, 50), wind_speed_10m = c(3, 3),
        shortwave_radiation = c(300, 300), w0 = 0
      )
      rain_before <- litter_wetness_vec(
        precipitation = c(1, 0), temperature_2m = c(20, 20),
        relative_humidity_2m = c(50, 50), wind_speed_10m = c(3, 3),
        shortwave_radiation = c(300, 300), w0 = 0
      )
      expect_equal(rain_before[1], 1)
      expect_gt(rain_before[2], no_rain_before[2])
    })
  })

  describe("dry-down", {

    it("is monotone non-increasing across a dry spell (constant met, no rain)", {
      out <- litter_wetness_vec(
        precipitation = rep(0, 5), temperature_2m = rep(20, 5),
        relative_humidity_2m = rep(50, 5), wind_speed_10m = rep(3, 5),
        shortwave_radiation = rep(300, 5), w0 = 1
      )
      expect_true(all(diff(out) <= 0))
    })

    it("worked oracle -- T=25, RH=40, U=5, SW=600, w0=1 -> w ~ 0.007649", {
      out <- litter_wetness_vec(
        precipitation = 0, temperature_2m = 25, relative_humidity_2m = 40,
        wind_speed_10m = 5, shortwave_radiation = 600, w0 = 1
      )
      expect_equal(out, 0.007649, tolerance = WET_TOL)
    })

    it("stronger VPD/wind/insolation forcing each dry faster (smaller w after one hour, w0=1)", {
      base <- litter_wetness_vec(0, 25, 40, 5, 600, w0 = 1)

      hotter  <- litter_wetness_vec(0, 35, 40, 5, 600, w0 = 1)   # higher VPD
      windier <- litter_wetness_vec(0, 25, 40, 15, 600, w0 = 1)  # higher wind
      sunnier <- litter_wetness_vec(0, 25, 40, 5, 1000, w0 = 1)  # higher insolation

      expect_lt(hotter, base)
      expect_lt(windier, base)
      expect_lt(sunnier, base)
    })
  })

  describe("output contract and validation", {

    it("stays within [0, 1] across an extreme parameter sweep, including rain resets", {
      g <- expand.grid(
        precip = c(0, 0.6, 5),
        temp   = c(-5, 15, 40),
        rh     = c(5, 50, 100),
        wind   = c(0, 10, 30),
        sw     = c(0, 400, 1000)
      )
      out <- litter_wetness_vec(g$precip, g$temp, g$rh, g$wind, g$sw)
      expect_false(anyNA(out))
      expect_true(all(out >= 0 & out <= 1))
    })

    it("rejects missing values in any input", {
      expect_error(litter_wetness_vec(c(0, NA), c(20, 20), c(50, 50), c(3, 3), c(300, 300)))
      expect_error(litter_wetness_vec(c(0, 0), c(20, NA), c(50, 50), c(3, 3), c(300, 300)))
      expect_error(litter_wetness_vec(c(0, 0), c(20, 20), c(50, NA), c(3, 3), c(300, 300)))
      expect_error(litter_wetness_vec(c(0, 0), c(20, 20), c(50, 50), c(3, NA), c(300, 300)))
      expect_error(litter_wetness_vec(c(0, 0), c(20, 20), c(50, 50), c(3, 3), c(300, NA)))
    })

    it("rejects negative precipitation, wind speed, and shortwave radiation, and out-of-range humidity", {
      expect_error(litter_wetness_vec(-1, 20, 50, 3, 300))     # negative precip
      expect_error(litter_wetness_vec(0, 20, 50, -1, 300))     # negative wind
      expect_error(litter_wetness_vec(0, 20, 50, 3, -1))       # negative SW
      expect_error(litter_wetness_vec(0, 20, -1, 3, 300))      # RH < 0
      expect_error(litter_wetness_vec(0, 20, 101, 3, 300))     # RH > 100
    })

    it("allows sub-zero temperature_2m (physically valid; not bounds-checked)", {
      expect_no_error(litter_wetness_vec(0, -10, 50, 3, 300))
    })

    it("rejects w0 outside [0, 1]", {
      expect_error(litter_wetness_vec(0, 20, 50, 3, 300, w0 = 1.5))
      expect_error(litter_wetness_vec(0, 20, 50, 3, 300, w0 = -0.1))
    })

    it("rejects inputs of differing length", {
      expect_error(litter_wetness_vec(c(0, 0), 20, 50, 3, 300))
    })

    it("rejects zero-length input", {
      expect_error(litter_wetness_vec(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0)))
    })
  })
})


describe(".litter_vpd() [internal Tetens VPD helper]", {

  it("matches the Tetens-based worked oracle (T=25, RH=40 -> vpd ~ 1.9007 kPa)", {
    expect_equal(.litter_vpd(25, 40), 1.9007, tolerance = WET_TOL)
  })

  it("is floored at zero at saturation (RH=100) and does not go negative for RH > 100", {
    expect_equal(.litter_vpd(25, 100), 0, tolerance = 1e-9)
    expect_gte(.litter_vpd(25, 150), 0)
  })
})


describe("litter_wetness() [met-tibble wrapper]", {

  .wet_met <- function(n = 3) {
    data.frame(
      precipitation         = rep(0, n),
      temperature_2m        = seq(15, 25, length.out = n),
      relative_humidity_2m  = seq(70, 40, length.out = n),
      wind_speed_10m        = seq(2, 8, length.out = n),
      shortwave_radiation   = seq(100, 600, length.out = n)
    )
  }

  it("agrees with the vector API called on the same columns", {
    met <- .wet_met(4)
    expect_equal(
      litter_wetness(met),
      litter_wetness_vec(
        met$precipitation, met$temperature_2m, met$relative_humidity_2m,
        met$wind_speed_10m, met$shortwave_radiation
      )
    )
  })

  it("returns one value in [0, 1] per row", {
    met <- .wet_met(5)
    out <- litter_wetness(met)
    expect_length(out, 5)
    expect_true(all(out >= 0 & out <= 1))
  })

  it("errors, naming the missing column, when a required column is absent", {
    met <- .wet_met()
    met$shortwave_radiation <- NULL
    expect_error(litter_wetness(met), regexp = "shortwave_radiation")
  })

  it("forwards calibration parameters through to the vector API", {
    met <- data.frame(
      precipitation = 0, temperature_2m = 25, relative_humidity_2m = 40,
      wind_speed_10m = 5, shortwave_radiation = 600
    )
    default_w <- litter_wetness(met, w0 = 1)
    slower_w  <- litter_wetness(met, w0 = 1, dry_rate_base = 0.1)
    expect_gt(slower_w, default_w)
  })
})
