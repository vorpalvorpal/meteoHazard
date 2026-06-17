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
