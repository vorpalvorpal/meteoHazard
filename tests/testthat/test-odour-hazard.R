# Tests for odour_hazard(): the receptor-independent ventilation index
# (source / ventilation), its G / peak-to-mean / scavenging overlays, the
# normalisation, and NA handling.

# met_data builder with sensible defaults; override any column by name.
mh <- function(n = 1, ...) {
  base <- list(
    wind_speed_10m         = 3,
    direct_radiation       = 0,
    cloud_cover            = 50,
    boundary_layer_height  = 500,
    temperature_2m         = 15,
    pressure_msl           = 1013,
    precipitation          = 0,
    relative_humidity_2m   = 50,
    soil_moisture_0_to_1cm = 0.1,
    soil_moisture_1_to_3cm = 0.1
  )
  ov <- list(...)
  base[names(ov)] <- ov
  as.data.frame(lapply(base, rep_len, n))
}

# ── Normalisation / wiring (characterisation) ───────────────────────────────
test_that("hazard matches the closed-form ventilation index for a known row", {
  # Defaults, single row: dP_mod = 0 (no lookback), R = 0, S_seal = 0, H = 0,
  # V_mod = 0.30*(15-10)/25 = 0.06 -> G = 1.06; night, wind 3-5, cloud 50%
  # -> s = 3 -> PM = 1 + 2*(3/5) = 2.2; W_rain = 1; u_eff = 3, h_mix = 500.
  # raw = 1.06*2.2/(3*500) = 0.00155467; ref = 3/(0.5*200) = 0.03.
  expect_equal(odour_hazard(mh()), 0.00155467 / 0.03, tolerance = 1e-4)
})

# ── Ventilation dominates ───────────────────────────────────────────────────
test_that("low-ventilation conditions vastly exceed high-ventilation ones", {
  bad  <- odour_hazard(mh(wind_speed_10m = 1.5, direct_radiation = 0,
                          cloud_cover = 10, boundary_layer_height = 100))
  good <- odour_hazard(mh(wind_speed_10m = 8, direct_radiation = 800,
                          cloud_cover = 0, boundary_layer_height = 2000))
  expect_gt(bad, good * 10)
})

# ── Source modifier G ───────────────────────────────────────────────────────
test_that("V_mod widens the temperature response to V_MOD_MAX (~0.30)", {
  h_cold <- odour_hazard(mh(temperature_2m = 5))   # V_mod = 0   -> G = 1.00
  h_hot  <- odour_hazard(mh(temperature_2m = 35))  # V_mod = 0.30 -> G = 1.30
  expect_equal(h_hot / h_cold, 1.30, tolerance = 1e-6)
})

test_that("falling pressure raises hazard once the 3-hour lookback exists", {
  # Pressure falling 3 hPa/h; everything else constant so hazard tracks G.
  d <- mh(n = 5, pressure_msl = c(1013, 1010, 1007, 1004, 1001))
  h <- odour_hazard(d)
  expect_equal(h[1], h[2]) # first 3 rows: no lookback, dP_mod = 0
  expect_gt(h[4], h[1])    # dP3 = -9 -> dP_mod = 0.30
})

test_that("post-rain piston surge raises hazard when the surface has dried", {
  wet  <- odour_hazard(mh(n = 26, precipitation = c(rep(1, 25), 0)))
  dry  <- odour_hazard(mh(n = 26, precipitation = rep(0, 26)))
  expect_gt(wet[26], dry[26]) # P_24 = 25 > 15 -> R_mod = 0.20
})

test_that("active-rain guard suppresses the piston surge during heavy rain", {
  # 25 h of rain (P_24 = 25 > 15) then a heavy active hour: R_mod = 0 (guard),
  # W_rain = 0.05. raw = G(1.06)*PM(s=3,=2.2)*0.05/(3*500); ref = 0.03.
  d <- mh(n = 26, precipitation = c(rep(1, 25), 5))
  expect_equal(odour_hazard(d)[26],
               (1.06 * 2.2 * 0.05 / (3 * 500)) / 0.03, tolerance = 1e-4)
})

# ── Peak-to-mean / stability ────────────────────────────────────────────────
test_that("hazard rises with stability via PM(s) at fixed ventilation", {
  # Same wind (4 m/s) and BL, so u_eff and h_mix are fixed; only s (hence PM)
  # varies: night clear -> E(4); night cloudy -> D(3); day strong -> B(1).
  h_stable  <- odour_hazard(mh(wind_speed_10m = 4, direct_radiation = 0, cloud_cover = 10))
  h_neutral <- odour_hazard(mh(wind_speed_10m = 4, direct_radiation = 0, cloud_cover = 90))
  h_unstab  <- odour_hazard(mh(wind_speed_10m = 4, direct_radiation = 800, cloud_cover = 0))
  expect_gt(h_stable, h_neutral)
  expect_gt(h_neutral, h_unstab)
})

# ── Scavenging ──────────────────────────────────────────────────────────────
test_that("heavy rain suppresses hazard via W_rain", {
  dry   <- odour_hazard(mh(precipitation = 0))
  heavy <- odour_hazard(mh(precipitation = 5))
  expect_equal(heavy / dry, 0.05, tolerance = 1e-6) # W_rain 0.05 vs 1.0
})

# ── NA handling / robustness ────────────────────────────────────────────────
test_that("an all-NA row yields a finite, conservative hazard (no NA out)", {
  d <- mh(
    wind_speed_10m = NA_real_, direct_radiation = NA_real_, cloud_cover = NA_real_,
    boundary_layer_height = NA_real_, temperature_2m = NA_real_,
    pressure_msl = NA_real_, precipitation = NA_real_,
    relative_humidity_2m = NA_real_,
    soil_moisture_0_to_1cm = NA_real_, soil_moisture_1_to_3cm = NA_real_
  )
  h <- odour_hazard(d)
  expect_false(is.na(h))
  expect_gt(h, 0)
})

test_that("a zero boundary-layer height yields a finite hazard, not Inf", {
  h <- odour_hazard(mh(boundary_layer_height = 0))
  expect_true(is.finite(h))
  expect_gt(h, 0)
})

test_that("output length equals the number of rows", {
  expect_length(odour_hazard(mh(n = 24)), 24)
})

# ── Validation ──────────────────────────────────────────────────────────────
test_that("missing required columns raise a classed input error", {
  expect_error(
    odour_hazard(mh()[, -1]),
    class = "meteoHazard_input_error"
  )
})

test_that("a non-numeric required column raises a classed input error", {
  d <- mh()
  d$temperature_2m <- "warm"
  expect_error(odour_hazard(d), class = "meteoHazard_input_error")
})

# ── Units handling ───────────────────────────────────────────────────────────
test_that("odour_hazard accepts units-tagged met columns and converts them", {
  bare   <- odour_hazard(mh())
  tagged <- mh()
  tagged$wind_speed_10m        <- units::set_units(10.8, "km/h")    # 3 m/s
  tagged$temperature_2m        <- units::set_units(59, "degree_F")  # 15 degC
  tagged$boundary_layer_height <- units::set_units(0.5, "km")       # 500 m
  tagged$pressure_msl          <- units::set_units(101300, "Pa")    # 1013 hPa
  expect_equal(odour_hazard(tagged), bare, tolerance = 1e-6)
})

test_that("odour_hazard rejects a met column tagged with incompatible units", {
  bad <- mh()
  bad$wind_speed_10m <- units::set_units(3, "degree_C")
  expect_error(odour_hazard(bad), class = "meteoHazard_input_error")
})

test_that("odour_hazard returns a plain numeric (relative) index", {
  expect_false(inherits(odour_hazard(mh()), "units"))
})

# === C2 behaviour spec (issue #15): odour_hazard() as a summary over the split ===
describe("odour_hazard() after the C2 ventilation/generation split", {

  it("returns numerically identical values to the current implementation for the existing cases", {
    # Run every fixture from the characterisation tests above and verify the
    # values are bit-for-bit identical.  These are the golden-test inputs.
    cases <- list(
      mh(),
      mh(wind_speed_10m = 1.5, direct_radiation = 0, cloud_cover = 10,
         boundary_layer_height = 100),
      mh(wind_speed_10m = 8, direct_radiation = 800, cloud_cover = 0,
         boundary_layer_height = 2000),
      mh(temperature_2m = 5),
      mh(temperature_2m = 35),
      mh(n = 5, pressure_msl = c(1013, 1010, 1007, 1004, 1001)),
      mh(n = 26, precipitation = c(rep(1, 25), 0)),
      mh(n = 26, precipitation = c(rep(1, 25), 5)),
      mh(precipitation = 5),
      mh(wind_speed_10m = 4, direct_radiation = 0, cloud_cover = 10),
      mh(wind_speed_10m = 4, direct_radiation = 800, cloud_cover = 0),
      mh(n = 24)
    )
    for (d in cases) {
      vs  <- ventilation_state(d, terrain = NULL)
      G   <- .odour_generation(d)
      raw <- G * vs$PM * vs$W_rain / (vs$u_eff * vs$h_mix)
      ref <- ODOUR_CONSTANTS$PM_MAX /
        (ODOUR_CONSTANTS$U_CALM_FLOOR * ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE)
      expected <- raw / ref
      expect_equal(odour_hazard(d), expected, tolerance = 1e-12)
    }
  })

  it("is composed from ventilation_state() and .odour_generation()", {
    d  <- mh(n = 4, wind_speed_10m = c(2, 5, 0.3, 8),
             direct_radiation = c(0, 600, 0, 800))
    vs <- ventilation_state(d, terrain = NULL)
    G  <- .odour_generation(d)
    ref <- ODOUR_CONSTANTS$PM_MAX /
      (ODOUR_CONSTANTS$U_CALM_FLOOR * ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE)
    manual <- G * vs$PM * vs$W_rain / (vs$u_eff * vs$h_mix) / ref
    expect_equal(odour_hazard(d), manual, tolerance = 1e-12)
  })

  # gap-fill — characterises EXISTING behaviour the current tests do not assert:
  it("is higher for a calm, shallow, stable hour than a windy, deep, neutral one", {
    h_bad  <- odour_hazard(mh(wind_speed_10m = 0.3, direct_radiation = 0,
                              cloud_cover = 5, boundary_layer_height = 100))
    h_good <- odour_hazard(mh(wind_speed_10m = 8,   direct_radiation = 500,
                              cloud_cover = 80, boundary_layer_height = 1500))
    expect_gt(h_bad, h_good)
  })

})


# ---------------------------------------------------------------------------
# C9 — .odour_hazard_raw() single-source-of-truth and shelter wiring
# ---------------------------------------------------------------------------

describe("odour_hazard(): .odour_hazard_raw() — single source of truth (C9)", {

  it(".odour_hazard_raw(G, vs) returns G*PM*W_rain/(u_eff*h_mix) (unit definition)", {
    vs <- list(PM = 2.5, W_rain = 0.8, u_eff = 3.0, h_mix = 400.0)
    G  <- 1.2
    expect_equal(meteoHazard:::.odour_hazard_raw(G, vs),
                 1.2 * 2.5 * 0.8 / (3.0 * 400.0))
  })

  it("odour_hazard() * H_ref_vent == .odour_hazard_raw(G, vs) for varied inputs", {
    # Structural tie: multiplying the normalized hazard back by the reference
    # constant must recover the raw helper exactly.
    d   <- mh(n = 4, wind_speed_10m = c(1, 3, 6, 0.3),
               direct_radiation = c(0, 600, 0, 0))
    vs  <- ventilation_state(d, terrain = NULL)
    G   <- .odour_generation(d)
    K   <- ODOUR_CONSTANTS
    ref <- K$PM_MAX / (K$U_CALM_FLOOR * K$H_MIX_FALLBACK_STABLE)
    expect_equal(odour_hazard(d) * ref,
                 meteoHazard:::.odour_hazard_raw(G, vs),
                 tolerance = 1e-12)
  })
})


describe("odour_hazard(): terrain / shelter params (C9)", {

  it("shelter = FALSE leaves hazard bit-identical to calling without terrain (regression)", {
    d   <- mh()
    ter <- mh_terrain(shelter_index = 50)
    expect_equal(
      odour_hazard(d, terrain = ter, shelter = FALSE),
      odour_hazard(d),
      tolerance = 0
    )
  })

  it("odour_hazard(shelter = TRUE) raises hazard on a sheltered low-wind night", {
    # shelter_index=50 (= SHELTER_ENCLOSED_REF → s_f=1, maximum shelter strength),
    # u10=1.5 (= SHELTER_U_FULL → w_r=1), no flow_convergence (M1 inactive).
    # reduction = SHELTER_MAX_REDUCTION * 1 * 1 = 0.7
    # u_eff_on  = max(1.5 * (1-0.7), 0.5) = max(0.45, 0.5) = 0.5 = U_CALM_FLOOR
    # u_eff_off = 1.5
    # hazard ∝ 1/u_eff  →  ratio = 1.5 / 0.5 = 3.
    d   <- mh(wind_speed_10m = 1.5, direct_radiation = 0, cloud_cover = 5,
               boundary_layer_height = 150)
    ter <- mh_terrain(shelter_index = 50)   # no flow_convergence → no drainage suppression
    h_off <- odour_hazard(d, terrain = ter, shelter = FALSE)
    h_on  <- odour_hazard(d, terrain = ter, shelter = TRUE)
    expect_gt(h_on, h_off)
    expect_equal(h_on / h_off, 3, tolerance = 0.01)
  })
})
