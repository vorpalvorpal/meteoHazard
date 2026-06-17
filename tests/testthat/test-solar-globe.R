# Tests for calculate_solar_position(), solve_globe_temp(), calculate_globe_temp()
#
# Solar position tests use well-known astronomical reference values.
# Globe temperature tests check physically expected behaviour and bounds.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
make_dt <- function(datetime_str, tz = "UTC") {
  as.POSIXct(datetime_str, tz = tz)
}

# Perth, Australia (latitude -31.95, longitude 115.86)
PERTH_LAT <- -31.95
PERTH_LON <- 115.86

# ---------------------------------------------------------------------------
# 1. calculate_solar_position
# ---------------------------------------------------------------------------

test_that("solar zenith is near 90 at sunrise for Perth (June solstice)", {
  # Perth sunrise on June 21 is ~07:04 AWST = 23:04 UTC (June 20).
  # At 23:18 UTC the zenith should be very close to 90°.  Allow ±2°.
  dt_sunrise <- make_dt("2024-06-20 23:18:00")
  pos <- calculate_solar_position(dt_sunrise, PERTH_LAT, PERTH_LON)
  expect_true(abs(pos$zenith - 90) < 2,
    label = sprintf("zenith near Perth sunrise = %.1f° (expected ~90°)", pos$zenith))
})

test_that("solar zenith is smallest near solar noon for Perth (December solstice)", {
  # Perth December solar noon is ~04:00 UTC (12:00 AWST = UTC+8).
  # Zenith at noon (~04:00 UTC) must be lower than at 01:00 or 07:00 UTC.
  dt_noon <- make_dt("2024-12-21 04:00:00")
  dt_morn <- make_dt("2024-12-21 01:00:00")
  dt_aftn <- make_dt("2024-12-21 07:00:00")

  noon_z <- calculate_solar_position(dt_noon, PERTH_LAT, PERTH_LON)$zenith
  morn_z <- calculate_solar_position(dt_morn, PERTH_LAT, PERTH_LON)$zenith
  aftn_z <- calculate_solar_position(dt_aftn, PERTH_LAT, PERTH_LON)$zenith

  expect_lt(noon_z, morn_z)
  expect_lt(noon_z, aftn_z)
})

test_that("solar noon zenith is within expected range for Perth in December", {
  # Dec 21 solar noon zenith ≈ |lat - declination| = |(-31.95) - (-23.4)| = 8.55°.
  # The algorithm is a low-precision approximation; allow ±3°.
  dt_noon <- make_dt("2024-12-21 04:00:00")
  pos <- calculate_solar_position(dt_noon, PERTH_LAT, PERTH_LON)
  expect_true(pos$zenith >= 5 && pos$zenith <= 15,
    label = sprintf("Dec solstice noon zenith = %.1f° (expected 5-15°)", pos$zenith))
})

test_that("zenith is > 90 at night for Perth (midnight AWST = 16:00 UTC)", {
  # Midnight AWST = 16:00 UTC; sun is well below horizon.
  dt_night <- make_dt("2024-06-21 16:00:00")
  pos <- calculate_solar_position(dt_night, PERTH_LAT, PERTH_LON)
  expect_gt(pos$zenith, 90)
})

test_that("azimuth is in [0, 360] range", {
  dts <- make_dt(c("2024-06-21 00:00:00", "2024-06-21 02:00:00",
                   "2024-06-21 04:00:00", "2024-12-21 04:00:00"))
  pos <- calculate_solar_position(dts, PERTH_LAT, PERTH_LON)
  expect_true(all(pos$azimuth >= 0 & pos$azimuth <= 360))
})

test_that("vectorised solar position returns same length as input", {
  dts <- make_dt(c("2024-06-15 00:00:00", "2024-06-15 02:00:00",
                   "2024-06-15 04:00:00"))
  pos <- calculate_solar_position(dts, PERTH_LAT, PERTH_LON)
  expect_length(pos$zenith, 3L)
  expect_length(pos$azimuth, 3L)
})

test_that("zenith increases monotonically from noon toward midnight (June 21 Perth)", {
  # After solar noon (~04:00 UTC June 21), zenith should increase monotonically
  # as the sun sets.
  hours_utc <- c(4, 6, 8, 10, 12)
  dts <- make_dt(sprintf("2024-06-21 %02d:00:00", hours_utc))
  zeniths <- vapply(dts, function(dt)
    calculate_solar_position(dt, PERTH_LAT, PERTH_LON)$zenith, numeric(1))
  # Each zenith should be greater than the previous
  expect_true(all(diff(zeniths) > 0),
    label = paste("Zeniths:", paste(round(zeniths, 1), collapse = ", ")))
})

# ---------------------------------------------------------------------------
# 2. solve_globe_temp
# ---------------------------------------------------------------------------

test_that("globe temp is less than air temp with no solar (sky cooling)", {
  # With zero solar, the sky temperature approximation (T_sky = temp - 20 °C)
  # causes net radiative cooling, so globe temp < air temp.
  tg <- solve_globe_temp(30, wind_speed = 1.0, direct_solar = 0,
                          diffuse_solar = 0, zenith = 90, albedo = 0.12)
  expect_lt(tg, 30)
  # But not excessively below air temp (sanity clamp at air - 5)
  expect_gte(tg, 30 - 5)
})

test_that("globe temp increases with increasing direct solar radiation", {
  tg_low  <- solve_globe_temp(30, 1.0, direct_solar = 100,  diffuse_solar = 50,  zenith = 30)
  tg_mid  <- solve_globe_temp(30, 1.0, direct_solar = 500,  diffuse_solar = 100, zenith = 30)
  tg_high <- solve_globe_temp(30, 1.0, direct_solar = 900,  diffuse_solar = 150, zenith = 30)
  expect_gt(tg_mid,  tg_low)
  expect_gt(tg_high, tg_mid)
})

test_that("globe temp decreases with increasing wind speed (for same solar load)", {
  # Higher wind → more convective cooling → lower globe temperature.
  tg_calm  <- solve_globe_temp(30, wind_speed = 0.2, direct_solar = 600,
                                diffuse_solar = 100, zenith = 30)
  tg_windy <- solve_globe_temp(30, wind_speed = 4.0, direct_solar = 600,
                                diffuse_solar = 100, zenith = 30)
  expect_gt(tg_calm, tg_windy)
})

test_that("globe temp is clamped to [air-5, air+30] sanity range", {
  tg <- solve_globe_temp(30, 1.0, direct_solar = 1000, diffuse_solar = 200, zenith = 0)
  expect_gte(tg, 30 - 5)
  expect_lte(tg, 30 + 30)
})

test_that("globe temp is plausible for typical outdoor summer conditions", {
  # Perth summer noon: ~35 °C air, 900 W/m² direct, 150 diffuse, wind 1 m/s.
  # Globe should be warmer than air (positive solar loading) but ≤ 65 °C.
  tg <- solve_globe_temp(35, wind_speed = 1.0, direct_solar = 900,
                          diffuse_solar = 150, zenith = 15, albedo = 0.12)
  expect_gt(tg, 35)
  expect_lt(tg, 65)
})

# ---------------------------------------------------------------------------
# 3. calculate_globe_temp (vectorised wrapper)
# ---------------------------------------------------------------------------

test_that("calculate_globe_temp returns correct length and type", {
  n <- 4L
  tg <- calculate_globe_temp(
    temp          = rep(30, n),
    wind_speed    = rep(1.0, n),
    direct_solar  = rep(500, n),
    diffuse_solar = rep(100, n),
    zenith        = rep(30,  n),
    albedo        = 0.12
  )
  expect_length(tg, n)
  expect_type(tg, "double")
})

test_that("calculate_globe_temp propagates NA", {
  tg <- calculate_globe_temp(
    temp          = c(NA_real_, 30),
    wind_speed    = c(1.0, 1.0),
    direct_solar  = c(500, 500),
    diffuse_solar = c(100, 100),
    zenith        = c(30,  30),
    albedo        = 0.12
  )
  expect_true(is.na(tg[1]))
  expect_false(is.na(tg[2]))
})

test_that("calculate_globe_temp matches scalar solve_globe_temp", {
  tg_vec <- calculate_globe_temp(30, 1.0, 600, 100, 30, albedo = 0.12)
  tg_scl <- solve_globe_temp(30, 1.0, 600, 100, 30, albedo = 0.12)
  expect_equal(tg_vec, tg_scl, tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# 4. Shared physical constants (Stage 0 / Stage B, new behaviour)
# ---------------------------------------------------------------------------
test_that("air thermal conductivity is the sourced 0.028 W/(m.K)", {
  # Brake (2001), Whillier Table 1 / BB p.470 — now shared by globe and wick.
  expect_equal(TWL_CONSTANTS$AIR_THERMAL_CONDUCTIVITY, 0.028)
})

test_that("air kinematic viscosity is a shared constant (1.5e-5 m^2/s)", {
  expect_equal(TWL_CONSTANTS$AIR_KINEMATIC_VISCOSITY, 1.5e-5)
})

test_that("globe temperature at the reference condition matches the recorded value", {
  # Characterisation guard pinning the post-change (0.028) globe output so it
  # cannot drift silently. Reference condition has no published exact value.
  skip(paste("pending: after the globe solver switches to conductivity 0.028,",
             "record solve_globe_temp(35, 1.0, 900, 150, 15, 0.12) and assert",
             "the result within +/-1.5 C of the recorded value"))
})
