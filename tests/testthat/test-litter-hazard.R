# Behaviour specification for the v3.1 Litter Hazard Index.
#
# v3.1 replaces the friction-velocity parameterisation (kappa, z0, ustar_t0,
# ustar_ref) with an identifiable gust-only parameterisation (gust_threshold,
# gust_reference), algebraically identical at the pinned defaults; adds a
# `material`-aware graded saturation treatment; adds an optional precomputed
# `wetness` input that supersedes soil moisture; smooths the hard rain gate
# into a ramp; and drops the `<=100` claim in favour of an unbounded relative
# index.
#
# Detection: the v3.1 `litter_hazard_vec()` introduces the `gust_threshold`
# formal (absent in v2/v3.0, which used `z0`); presence of that formal is the
# "v3.1 is implemented" signal.
#
# Wind inputs are in m/s (Open-Meteo fetched with &wind_speed_unit=ms),
# matching the rest of the package.
#
# Pinned gust-space constants, derived so the reparameterisation is
# algebraically identical to the old friction-velocity form at defaults:
#   gust_threshold = 0.30 / (0.40 / ln(200)) = 3.9737355 m/s
#   gust_reference  = 1.05 / (0.40 / ln(200)) = 13.9080743 m/s

litter_v3_available <- function() {
  exists("litter_hazard_vec", mode = "function") &&
    "gust_threshold" %in% names(formals(litter_hazard_vec))
}

skip_if_no_litter_v3 <- function() {
  testthat::skip_if_not(litter_v3_available(), "litter v3.1 API (gust_threshold) not yet implemented")
}

# Tolerance for known-answer comparisons: the index is O(100) and computed with
# a handful of double-precision operations, so 1e-3 is far tighter than any
# operational meaning while leaving ample headroom over floating-point error.
LRI_TOL <- 1e-3


describe("litter_hazard_vec() [v3.1]", {

  describe("output contract and bounds", {

    it("returns one numeric value per input hour", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m         = c(3, 10, 15),
        wind_speed_10m         = c(1.5, 7, 12.5),
        precipitation          = c(0, 0, 0),
        soil_moisture_0_to_1cm = c(0.02, 0.10, 0.02)
      )
      expect_type(out, "double")
      expect_length(out, 3)
    })

    it("is always non-negative and finite (the <=100 claim is dropped)", {
      skip_if_no_litter_v3()
      # Sweep a wide grid at default parameters; only the lower bound and
      # finiteness are contractual now (issue #11).
      g <- expand.grid(
        gust = c(0, 4, 11, 22, 42),
        wind = c(0, 6, 17, 33),
        rain = c(0, 0.3, 1),
        sm   = c(0, 0.05, 0.12, 0.30)
      )
      out <- litter_hazard_vec(g$gust, g$wind, g$rain, g$sm)
      expect_true(all(out >= 0 & is.finite(out)))
    })

    it("an inflated entrainment_max pushes the output past the historical 100 cap", {
      skip_if_no_litter_v3()
      # entrainment_max=60 -> E0=60 (saturated), T(15)=2 -> 120 (> 100 allowed).
      out <- litter_hazard_vec(16, 15, 0, 0.02, entrainment_max = 60)
      expect_equal(out, 120, tolerance = LRI_TOL)
      expect_gt(out, 100)
    })

    it("reaches exactly 100 only with saturated entrainment and saturating wind", {
      skip_if_no_litter_v3()
      # E = entrainment_max (50) and T = transport_max (2) -> 50 * 2 = 100.
      out <- litter_hazard_vec(
        wind_gusts_10m = 16, wind_speed_10m = 15,
        precipitation = 0, soil_moisture_0_to_1cm = 0.02
      )
      expect_equal(out, 100, tolerance = LRI_TOL)
    })

    it("is direction-agnostic: takes no wind direction, CAPE, is_day, or receptor arcs", {
      skip_if_no_litter_v3()
      fmls <- names(formals(litter_hazard_vec))
      expect_false("wind_direction_10m" %in% fmls)
      expect_false("cape" %in% fmls)
      expect_false("is_day" %in% fmls)
      expect_false("receptor_arcs" %in% fmls)
    })

    it("carries the four v3 meteorological inputs as named formals", {
      skip_if_no_litter_v3()
      fmls <- names(formals(litter_hazard_vec))
      expect_true(all(
        c("wind_gusts_10m", "wind_speed_10m", "precipitation",
          "soil_moisture_0_to_1cm") %in% fmls
      ))
    })

    it("exposes gust_threshold/gust_reference and no longer exposes z0", {
      skip_if_no_litter_v3()
      fmls <- names(formals(litter_hazard_vec))
      expect_true("gust_threshold" %in% fmls)
      expect_true("gust_reference" %in% fmls)
      expect_false("z0" %in% fmls)
      expect_false("kappa" %in% fmls)
      expect_false("ustar_t0" %in% fmls)
      expect_false("ustar_ref" %in% fmls)
    })
  })

  describe("soil moisture / wetness contract seam", {

    it("soil_moisture_0_to_1cm defaults to NULL (optional when wetness is supplied)", {
      skip_if_no_litter_v3()
      expect_null(eval(formals(litter_hazard_vec)$soil_moisture_0_to_1cm))
    })

    it("accepts wetness in place of soil_moisture_0_to_1cm", {
      skip_if_no_litter_v3()
      expect_no_error(
        litter_hazard_vec(16, 15, 0, wetness = 0.02)
      )
    })

    it("errors when neither soil_moisture_0_to_1cm nor wetness is supplied", {
      skip_if_no_litter_v3()
      expect_error(
        litter_hazard_vec(16, 15, 0),
        class = "meteoHazard_input_error"
      )
    })

    it("wetness supersedes soil_moisture_0_to_1cm when both are supplied, with a warning", {
      skip_if_no_litter_v3()
      # wetness=0.02 is far below paper_veto_wetness (0.8): if wetness wins,
      # the (saturated, SM=0.30) veto path is NOT taken and the result is > 0.
      expect_warning(
        out <- litter_hazard_vec(
          16, 15, 0,
          soil_moisture_0_to_1cm = 0.30, wetness = 0.02, material = "paper"
        ),
        class = "meteoHazard_litter_wetness_supersedes"
      )
      expect_gt(out, 0)
    })

    it("positional callers litter_hazard_vec(g, w, r, sm) are unaffected", {
      skip_if_no_litter_v3()
      # The 4th positional argument still binds to soil_moisture_0_to_1cm.
      out <- litter_hazard_vec(16, 12.5, 0, 0.02)
      expect_equal(out, 86.8421, tolerance = LRI_TOL)
    })
  })

  describe("entrainment (gust -> friction velocity)", {

    it("is zero when the gust is below the dry threshold", {
      skip_if_no_litter_v3()
      # dry threshold gust = 3.9737355 m/s; 3.9 sits below it -> E = 0.
      out <- litter_hazard_vec(
        wind_gusts_10m = c(0, 3, 3.9), wind_speed_10m = c(12, 12, 12),
        precipitation = c(0, 0, 0), soil_moisture_0_to_1cm = c(0.02, 0.02, 0.02)
      )
      expect_equal(out, c(0, 0, 0))
    })

    it("dry threshold gust 3.9 gives 0, gust 4.1 gives > 0", {
      skip_if_no_litter_v3()
      out_below <- litter_hazard_vec(3.9, 12, 0, 0.02)
      out_above <- litter_hazard_vec(4.1, 12, 0, 0.02)
      expect_equal(out_below, 0)
      expect_gt(out_above, 0)
    })

    it("increases monotonically with gust above the threshold", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = c(3, 6, 9, 12, 15, 20),
        wind_speed_10m = rep(12, 6),
        precipitation = rep(0, 6),
        soil_moisture_0_to_1cm = rep(0.02, 6)
      )
      expect_true(all(diff(out) >= 0))
    })

    it("saturates above the reference gust (14/18/35 all equal, dry, W=12)", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = c(14, 18, 35), wind_speed_10m = rep(12, 3),
        precipitation = rep(0, 3), soil_moisture_0_to_1cm = rep(0.02, 3)
      )
      expect_equal(out[1], out[2], tolerance = LRI_TOL)
      expect_equal(out[2], out[3], tolerance = LRI_TOL)
    })

    it("follows the excess-squared shape (convex ramp between threshold and reference)", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = c(6, 8, 10, 12), wind_speed_10m = rep(1.5, 4),
        precipitation = rep(0, 4), soil_moisture_0_to_1cm = rep(0.02, 4)
      )
      incr <- diff(out)
      expect_true(all(diff(incr) > 0))
    })

    it("errors (classed) when gust_reference does not exceed gust_threshold", {
      skip_if_no_litter_v3()
      expect_error(
        litter_hazard_vec(12, 8, 0, 0.1, gust_threshold = 5, gust_reference = 5),
        class = "meteoHazard_input_error"
      )
      expect_error(
        litter_hazard_vec(12, 8, 0, 0.1, gust_threshold = 5, gust_reference = 4),
        class = "meteoHazard_input_error"
      )
    })
  })

  describe("moisture-raised entrainment threshold (below saturation, both materials)", {

    it("lowers the hazard as the surface dampens, at fixed gust", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = rep(12, 4), wind_speed_10m = rep(8, 4),
        precipitation = rep(0, 4),
        soil_moisture_0_to_1cm = c(0.05, 0.10, 0.15, 0.19)
      )
      expect_true(all(diff(out) <= 0))
    })

    it("treats all moisture at or below the dry threshold identically", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = rep(12, 3), wind_speed_10m = rep(8, 3),
        precipitation = rep(0, 3), soil_moisture_0_to_1cm = c(0.00, 0.02, 0.05)
      )
      expect_equal(out[1], out[2], tolerance = LRI_TOL)
      expect_equal(out[2], out[3], tolerance = LRI_TOL)
    })

    it("a damp surface needs a stronger gust than a dry surface to entrain", {
      skip_if_no_litter_v3()
      dry  <- litter_hazard_vec(9, 1.5, 0, 0.02)
      damp <- litter_hazard_vec(9, 1.5, 0, 0.14)
      expect_gt(dry, damp)
    })
  })

  describe("material-aware graded saturation", {

    it("default material is film -- an unspecified material gives the graded residual, not a hard veto", {
      skip_if_no_litter_v3()
      # entrainment_max*(1-saturation_penalty)*T(15) = 50*0.3*2 = 30.
      out <- litter_hazard_vec(35, 15, 0, 0.30)
      expect_gt(out, 0)
      expect_equal(out, 30, tolerance = LRI_TOL)
    })

    it("film: saturated entrainment is reduced by saturation_penalty, not vetoed", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = c(16, 35), wind_speed_10m = c(13, 15),
        precipitation = c(0, 0), soil_moisture_0_to_1cm = c(0.20, 0.30),
        material = "film"
      )
      expect_equal(out, c(4.5248, 30), tolerance = LRI_TOL)
    })

    it("paper: a saturated surface is hard-vetoed to zero (re-pointed from the old default-material test)", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = c(16, 35), wind_speed_10m = c(13, 15),
        precipitation = c(0, 0), soil_moisture_0_to_1cm = c(0.20, 0.30),
        material = "paper"
      )
      expect_equal(out, c(0, 0))
    })

    it("the graded penalty does not touch the unsaturated ramp, for either material", {
      skip_if_no_litter_v3()
      # G=10, W=7, P=0, SM=0.10 (damp, not saturated) -> 1.212686 regardless
      # of material; the film/paper split only applies at s == 1.
      out_film  <- litter_hazard_vec(10, 7, 0, 0.10, material = "film")
      out_paper <- litter_hazard_vec(10, 7, 0, 0.10, material = "paper")
      expect_equal(out_film,  1.2127, tolerance = LRI_TOL)
      expect_equal(out_paper, 1.2127, tolerance = LRI_TOL)
    })

    it("rejects an unrecognised material value", {
      skip_if_no_litter_v3()
      expect_error(litter_hazard_vec(16, 12.5, 0, 0.02, material = "cardboard"))
    })
  })

  describe("smooth rain ramp (replaces the hard gate)", {

    it("P <= rain_onset (0.2 mm) is fully ungated", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(16, 15, 0.2, 0.02)
      expect_equal(out, 100, tolerance = LRI_TOL)
    })

    it("P >= rain_threshold (0.5 mm) is fully gated to zero", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = c(16, 16), wind_speed_10m = c(15, 15),
        precipitation = c(0.5, 3.0), soil_moisture_0_to_1cm = c(0.02, 0.02)
      )
      expect_equal(out, c(0, 0))
    })

    it("mid-ramp P=0.4 sits strictly between 0 and 100, at the pinned ramp value", {
      skip_if_no_litter_v3()
      # 100 * (1 - (0.4-0.2)/(0.5-0.2)) = 100 * (1 - 0.66667) = 33.333.
      out <- litter_hazard_vec(16, 15, 0.4, 0.02)
      expect_gt(out, 0)
      expect_lt(out, 100)
      expect_equal(out, 33.333, tolerance = LRI_TOL)
    })

    it("errors (classed) when rain_threshold does not exceed rain_onset", {
      skip_if_no_litter_v3()
      expect_error(
        litter_hazard_vec(12, 8, 0, 0.1, rain_onset = 0.5, rain_threshold = 0.5),
        class = "meteoHazard_input_error"
      )
      expect_error(
        litter_hazard_vec(12, 8, 0, 0.1, rain_onset = 0.6, rain_threshold = 0.5),
        class = "meteoHazard_input_error"
      )
    })
  })

  describe("transport multiplier (mean wind)", {

    it("applies no amplification at or below the transport onset wind", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = c(16, 16), wind_speed_10m = c(1.5, 5.5),
        precipitation = c(0, 0), soil_moisture_0_to_1cm = c(0.02, 0.02)
      )
      expect_equal(out[1], out[2], tolerance = LRI_TOL)
      expect_equal(out[1], 50, tolerance = LRI_TOL)
    })

    it("increases the hazard with mean wind above the onset", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = rep(16, 6),
        wind_speed_10m = c(0, 5.5, 8, 11, 15, 20),
        precipitation = rep(0, 6),
        soil_moisture_0_to_1cm = rep(0.02, 6)
      )
      expect_true(all(diff(out) >= 0))
    })

    it("saturates the transport amplification at the reference wind", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(
        wind_gusts_10m = rep(16, 3), wind_speed_10m = c(15, 20, 35),
        precipitation = rep(0, 3), soil_moisture_0_to_1cm = rep(0.02, 3)
      )
      expect_equal(out[1], out[2], tolerance = LRI_TOL)
      expect_equal(out[2], out[3], tolerance = LRI_TOL)
    })

    it("never reduces the hazard below the entrainment-only value", {
      skip_if_no_litter_v3()
      calm   <- litter_hazard_vec(12, 0, 0, 0.02)
      windy  <- litter_hazard_vec(12, 16, 0, 0.02)
      expect_gte(windy, calm)
    })

    it("errors (classed) when wind_transport_ref does not exceed wind_transport_onset", {
      skip_if_no_litter_v3()
      expect_error(
        litter_hazard_vec(12, 8, 0, 0.1, wind_transport_onset = 16, wind_transport_ref = 15),
        class = "meteoHazard_input_error"
      )
    })
  })

  describe("multiplicative veto property", {

    it("returns zero if any single suppressing condition holds (paper material isolates the saturation veto)", {
      skip_if_no_litter_v3()
      # row 1: rain veto; row 2: saturation veto (paper); row 3: sub-threshold gust veto.
      out <- litter_hazard_vec(
        wind_gusts_10m         = c(16, 16, 3),
        wind_speed_10m         = c(15, 15, 15),
        precipitation          = c(2.0, 0.0, 0.0),
        soil_moisture_0_to_1cm = c(0.02, 0.25, 0.02),
        material = "paper"
      )
      expect_equal(out, c(0, 0, 0))
    })
  })

  describe("worked examples (preserved unchanged under the gust-space reparameterisation)", {

    it("Example 1 -- strong gust, strong wind, dry -> 86.8421 (EXTREME)", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(16, 12.5, 0, 0.02)
      expect_equal(out, 86.8421, tolerance = LRI_TOL)
    })

    it("Example 2 -- moderate gust on a damp surface -> 1.2127 (LOW)", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(10, 7, 0, 0.10)
      expect_equal(out, 1.2127, tolerance = LRI_TOL)
    })

    it("Example 3 -- rain now -> 0 (LOW)", {
      skip_if_no_litter_v3()
      out <- litter_hazard_vec(16, 12.5, 3.0, 0.02)
      expect_equal(out, 0)
    })
  })

  describe("input validation", {

    it("rejects missing values in any input (no NA imputation)", {
      skip_if_no_litter_v3()
      expect_error(litter_hazard_vec(c(12, NA), c(8, 8), c(0, 0), c(0.1, 0.1)))
      expect_error(litter_hazard_vec(c(12, 12), c(8, NA), c(0, 0), c(0.1, 0.1)))
      expect_error(litter_hazard_vec(c(12, 12), c(8, 8), c(0, NA), c(0.1, 0.1)))
      expect_error(litter_hazard_vec(c(12, 12), c(8, 8), c(0, 0), c(0.1, NA)))
    })

    it("rejects inputs of differing length", {
      skip_if_no_litter_v3()
      expect_error(litter_hazard_vec(c(12, 12), 8, 0, 0.1))
    })

    it("rejects zero-length input", {
      skip_if_no_litter_v3()
      expect_error(litter_hazard_vec(numeric(0), numeric(0), numeric(0), numeric(0)))
    })

    it("rejects physically impossible meteorology", {
      skip_if_no_litter_v3()
      expect_error(litter_hazard_vec(-1, 8, 0, 0.1))    # negative gust
      expect_error(litter_hazard_vec(12, -1, 0, 0.1))   # negative wind
      expect_error(litter_hazard_vec(12, 8, -1, 0.1))   # negative precip
      expect_error(litter_hazard_vec(12, 8, 0, 1.5))    # soil moisture > 1
      expect_error(litter_hazard_vec(12, 8, 0, -0.1))   # soil moisture < 0
    })

    it("rejects parameter sets that violate the soil moisture ordering constraint", {
      skip_if_no_litter_v3()
      expect_error(litter_hazard_vec(12, 8, 0, 0.1, soil_dry = 0.20, soil_wet = 0.10))
    })
  })
})


describe("litter_hazard() [v3.1, basic mode: use_wetness_state = FALSE]", {

  v3_met <- function(n = 3) {
    data.frame(
      wind_gusts_10m         = seq(3, 16, length.out = n),
      wind_speed_10m         = seq(1.5, 15, length.out = n),
      precipitation          = rep(0, n),
      soil_moisture_0_to_1cm = rep(0.02, n)
    )
  }

  it("computes one finite, non-negative index per row of the met tibble", {
    skip_if_no_litter_v3()
    out <- litter_hazard(v3_met(4))
    expect_type(out, "double")
    expect_length(out, 4)
    expect_true(all(out >= 0 & is.finite(out)))
  })

  it("agrees with the vector API called on the same columns", {
    skip_if_no_litter_v3()
    met <- v3_met(5)
    expect_equal(
      litter_hazard(met),
      litter_hazard_vec(
        met$wind_gusts_10m, met$wind_speed_10m,
        met$precipitation, met$soil_moisture_0_to_1cm
      )
    )
  })

  it("does not require the removed v2 columns (cape, is_day, wind_direction)", {
    skip_if_no_litter_v3()
    expect_no_error(litter_hazard(v3_met()))
  })

  it("errors, naming the missing column, when soil_moisture_0_to_1cm is absent (default use_wetness_state=FALSE)", {
    skip_if_no_litter_v3()
    met <- v3_met()
    met$soil_moisture_0_to_1cm <- NULL
    expect_error(litter_hazard(met), regexp = "soil_moisture_0_to_1cm")
  })

  it("forwards calibration parameters through to the vector API", {
    skip_if_no_litter_v3()
    met <- data.frame(
      wind_gusts_10m = 16, wind_speed_10m = 15,
      precipitation = 0.6, soil_moisture_0_to_1cm = 0.02
    )
    expect_equal(litter_hazard(met), 0)            # gated by default rain_threshold
    expect_gt(litter_hazard(met, rain_threshold = 1.0), 0)
  })
})


describe("litter_hazard() [wetness-state integration: use_wetness_state = TRUE]", {

  # precip = 0.5 (>= default wetness_set_precip) forces wetness = 1 on this
  # single hour regardless of T/RH/wind/SW; rain_onset/rain_threshold are
  # raised so the SAME hour's hazard rain_gate stays fully open (1), which
  # decouples the wetness-reset trigger from the hazard rain gate (both
  # default to the coincidental value 0.5mm otherwise).
  wet_met <- function(material_gust = 35, material_wind = 15) {
    data.frame(
      wind_gusts_10m       = material_gust,
      wind_speed_10m       = material_wind,
      precipitation        = 0.5,
      temperature_2m       = 25,
      relative_humidity_2m = 40,
      shortwave_radiation  = 600
    )
  }

  it("runs on the five drying columns without soil_moisture_0_to_1cm", {
    skip_if_no_litter_v3()
    met <- wet_met()
    out <- litter_hazard(met, use_wetness_state = TRUE)
    expect_type(out, "double")
    expect_length(out, 1)
    expect_true(is.finite(out) && out >= 0)
  })

  it("errors, naming the missing column, when a required drying column is absent", {
    skip_if_no_litter_v3()
    met <- wet_met()
    met$shortwave_radiation <- NULL
    expect_error(
      litter_hazard(met, use_wetness_state = TRUE),
      regexp = "shortwave_radiation"
    )
  })

  it("wetness == 1 + material='paper' -> hazard is 0 (paper veto)", {
    skip_if_no_litter_v3()
    met <- wet_met()
    out <- litter_hazard(
      met, use_wetness_state = TRUE, material = "paper",
      rain_onset = 1, rain_threshold = 2
    )
    expect_equal(out, 0)
  })

  it("wetness == 1 + material='film' -> hazard is the graded residual (== 30, matching the SM-saturated oracle)", {
    skip_if_no_litter_v3()
    met <- wet_met()
    out <- litter_hazard(
      met, use_wetness_state = TRUE, material = "film",
      rain_onset = 1, rain_threshold = 2
    )
    expect_equal(out, 30, tolerance = LRI_TOL)
  })
})


describe("citation hygiene", {

  # Grep-style check on the built .Rd source. Reads relative to the working
  # directory testthat sets during
  # devtools::test()/test_check() (tests/testthat/), i.e. two levels up to
  # the package root's man/ directory.
  .litter_hazard_rd_text <- function() {
    rd_path <- "../../man/litter_hazard_vec.Rd"
    testthat::skip_if_not(file.exists(rd_path), "man/litter_hazard_vec.Rd not found")
    paste(readLines(rd_path, warn = FALSE), collapse = "\n")
  }

  it("cites Scientific Reports article number 3898, not the old (wrong) 5006", {
    txt <- .litter_hazard_rd_text()
    expect_true(grepl("3898", txt, fixed = TRUE))
    expect_false(grepl("5006", txt, fixed = TRUE))
  })

  it("credits Mani as the citation's author, not Roebroek", {
    txt <- .litter_hazard_rd_text()
    expect_true(grepl("Mani", txt, fixed = TRUE))
    expect_false(grepl("Roebroek", txt, fixed = TRUE))
  })
})


describe("litter_hazard_vec() [units handling]", {

  it("accepts units-tagged wind inputs and converts them to m/s", {
    skip_if_no_litter_v3()
    bare   <- litter_hazard_vec(16, 12.5, 0, 0.02)
    tagged <- litter_hazard_vec(
      units::set_units(57.6, "km/h"), units::set_units(45, "km/h"), 0, 0.02
    )
    expect_equal(tagged, bare, tolerance = LRI_TOL)
  })

  it("converts a units-tagged precipitation input", {
    skip_if_no_litter_v3()
    tagged <- litter_hazard_vec(16, 15, units::set_units(0.05, "cm"), 0.02)
    expect_equal(tagged, 0)
  })

  it("rejects a wind input tagged with dimensionally incompatible units", {
    skip_if_no_litter_v3()
    expect_error(
      litter_hazard_vec(units::set_units(16, "degree_C"), 15, 0, 0.02),
      class = "meteoHazard_input_error"
    )
  })

  it("returns a plain numeric index, not a units object", {
    skip_if_no_litter_v3()
    out <- litter_hazard_vec(16, 12.5, 0, 0.02)
    expect_false(inherits(out, "units"))
    expect_type(out, "double")
  })

  it("litter_hazard() accepts a met tibble carrying units columns", {
    skip_if_no_litter_v3()
    met <- data.frame(precipitation = 0, soil_moisture_0_to_1cm = 0.02)
    met$wind_gusts_10m <- units::set_units(57.6, "km/h")
    met$wind_speed_10m <- units::set_units(45, "km/h")
    expect_equal(litter_hazard(met), litter_hazard_vec(16, 12.5, 0, 0.02),
                 tolerance = LRI_TOL)
  })
})
