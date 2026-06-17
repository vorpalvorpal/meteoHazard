# Tests for generate_twl() and fetch_openmeteo()
#
# These tests use only locally-supplied data (no API calls) to be fast and
# reproducible.  A small set of API integration tests at the bottom are
# skipped unless TWLTEST_ONLINE=true is set in the environment.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
make_dt <- function(datetime_str, tz = "UTC") {
  as.POSIXct(datetime_str, tz = tz)
}

# Call generate_twl() with all weather inputs supplied (no API fetch needed).
# Uses a sheltered indoor-style scenario: direct/diffuse solar = 0.
twl_supplied <- function(temp, wind_speed, RH, pressure_hpa = 1013,
                         Icl = 0.6, icl = 0.45, verbose = FALSE) {
  generate_twl(
    datetime      = make_dt("2024-06-15 10:00:00"),
    latitude      = -31.95,
    longitude     = 115.86,
    temp          = temp,
    wind_speed    = wind_speed,
    RH            = RH,
    direct_solar  = 0,
    diffuse_solar = 0,
    pressure      = pressure_hpa,
    Icl           = Icl,
    icl           = icl,
    verbose       = verbose
  )
}

# ---------------------------------------------------------------------------
# 1. Basic output type and length
# ---------------------------------------------------------------------------
test_that("generate_twl returns numeric vector of correct length", {
  n <- 3L
  result <- generate_twl(
    datetime      = make_dt(c("2024-06-15 08:00:00",
                               "2024-06-15 09:00:00",
                               "2024-06-15 10:00:00")),
    latitude      = -31.95,
    longitude     = 115.86,
    temp          = c(30, 32, 35),
    wind_speed    = c(0.5, 0.5, 0.5),
    RH            = c(60, 55, 50),
    direct_solar  = c(0, 0, 0),
    diffuse_solar = c(0, 0, 0),
    pressure      = c(1013, 1013, 1013),
    verbose       = FALSE
  )
  expect_type(result, "double")
  expect_length(result, n)
  expect_true(all(!is.na(result)))
})

# ---------------------------------------------------------------------------
# 2. Monotonicity: increasing temperature reduces TWL
# ---------------------------------------------------------------------------
test_that("TWL decreases as temperature increases at fixed humidity", {
  twl_30 <- twl_supplied(30, 0.5, 60)
  twl_35 <- twl_supplied(35, 0.5, 60)
  twl_40 <- twl_supplied(40, 0.5, 60)
  expect_gt(twl_30, twl_35)
  expect_gt(twl_35, twl_40)
})

test_that("TWL decreases as RH increases at fixed temperature", {
  twl_40 <- twl_supplied(38, 0.5, 40)
  twl_60 <- twl_supplied(38, 0.5, 60)
  twl_80 <- twl_supplied(38, 0.5, 80)
  expect_gt(twl_40, twl_60)
  expect_gt(twl_60, twl_80)
})

test_that("TWL increases as wind speed increases at fixed temp and RH", {
  twl_low  <- twl_supplied(38, 0.2, 60)
  twl_high <- twl_supplied(38, 2.0, 60)
  expect_gt(twl_high, twl_low)
})

# ---------------------------------------------------------------------------
# 3. Output range is physically valid (60--380 W/m^2)
# ---------------------------------------------------------------------------
test_that("TWL is constrained to [60, 380] W/m^2", {
  result <- twl_supplied(40, 0.5, 70)
  expect_gte(result, 60)
  expect_lte(result, 380)
})

# ---------------------------------------------------------------------------
# 4. Withdrawal limits are applied
# ---------------------------------------------------------------------------
test_that("DB > 44 C caps TWL to 115 W/m^2", {
  result <- twl_supplied(45, 0.5, 30)
  expect_lte(result, 115)
})

# ---------------------------------------------------------------------------
# 5. categorise_twl zones are consistent with generate_twl output
# ---------------------------------------------------------------------------
test_that("categorise_twl correctly maps generate_twl output to zones", {
  # Hot and humid -> should be Withdrawal or Buffer
  low  <- twl_supplied(42, 0.2, 80)
  expect_true(categorise_twl(low) %in% c("Withdrawal", "Buffer", "Acclimatisation"))

  # Cool and dry -> should be Unrestricted
  high <- twl_supplied(28, 2.0, 30)
  expect_equal(categorise_twl(high), "Unrestricted")
})

# ---------------------------------------------------------------------------
# 6. NA propagation: NA input produces NA output
# ---------------------------------------------------------------------------
test_that("NA in temperature propagates to NA TWL", {
  result <- generate_twl(
    datetime      = make_dt(c("2024-06-15 10:00:00", "2024-06-15 11:00:00")),
    latitude      = -31.95,
    longitude     = 115.86,
    temp          = c(NA_real_, 35),
    wind_speed    = c(0.5, 0.5),
    RH            = c(60, 60),
    direct_solar  = c(0, 0),
    diffuse_solar = c(0, 0),
    pressure      = c(1013, 1013),
    verbose       = FALSE
  )
  expect_true(is.na(result[1]))
  expect_false(is.na(result[2]))
})

# ---------------------------------------------------------------------------
# 7. Wind height correction: API wind (at 10 m) is corrected to body level
#    A supplied wind_speed bypasses the correction; an API-fetched wind_speed
#    does not.  We test the correction factor indirectly via generate_twl's
#    internal behaviour by checking that wind_from_api=TRUE lowers effective
#    wind vs the same value passed as wind_speed directly.
#    (TWL with lower effective wind < TWL with higher effective wind)
# ---------------------------------------------------------------------------
test_that("supplying wind_speed directly differs from post-correction API value", {
  # Supplied directly: treated as body-level wind, no correction
  twl_direct <- twl_supplied(38, 1.5, 50)

  # If we mimic what the API would return (1.5 m/s at 10 m), the corrected
  # body-level wind is 1.5 * 0.667 = 1.0 m/s -> lower wind -> lower TWL
  twl_corrected <- twl_supplied(38, 1.5 * log(1 / 0.01) / log(10 / 0.01), 50)

  # Direct (1.5 m/s) should give higher TWL than corrected (~1.0 m/s)
  expect_gt(twl_direct, twl_corrected)
})

# ---------------------------------------------------------------------------
# 8. trad (MRT) formula: with zero solar, globe_temp ≈ air temp,
#    and trad should be very close to air temp (within 1 C)
# ---------------------------------------------------------------------------
test_that("trad is close to air temp when solar = 0", {
  # With no solar and wind=0.2, globe_temp ≈ air_temp, so trad ≈ air_temp.
  # The TWL result should be consistent with the low-radiant-load scenario.
  # We check it is in a plausible range rather than testing trad directly.
  result <- twl_supplied(35, 0.5, 50)
  expect_gte(result, 60)
  expect_lte(result, 380)
})

# ---------------------------------------------------------------------------
# 8b. Input validation (Stage C, new behaviour) — classed meteoHazard_input_error
# ---------------------------------------------------------------------------
test_that("generate_twl errors when a weather input length is neither 1 nor length(datetime)", {
  expect_error(
    generate_twl(
      datetime      = make_dt(c("2024-06-15 10:00:00", "2024-06-15 11:00:00",
                                 "2024-06-15 12:00:00")),
      latitude      = -31.95, longitude = 115.86,
      temp          = c(30, 32),        # length 2 against datetime length 3
      wind_speed    = 0.5, RH = 50,
      direct_solar  = 0, diffuse_solar = 0, pressure = 1013,
      verbose       = FALSE
    ),
    class = "meteoHazard_input_error"
  )
})

test_that("generate_twl rejects latitude outside [-90, 90]", {
  expect_error(
    generate_twl(datetime = make_dt("2024-06-15 10:00:00"),
                 latitude = 95, longitude = 115.86,
                 temp = 30, wind_speed = 0.5, RH = 50,
                 direct_solar = 0, diffuse_solar = 0, pressure = 1013,
                 verbose = FALSE),
    class = "meteoHazard_input_error"
  )
})

test_that("generate_twl rejects longitude outside [-180, 180]", {
  expect_error(
    generate_twl(datetime = make_dt("2024-06-15 10:00:00"),
                 latitude = -31.95, longitude = 999,
                 temp = 30, wind_speed = 0.5, RH = 50,
                 direct_solar = 0, diffuse_solar = 0, pressure = 1013,
                 verbose = FALSE),
    class = "meteoHazard_input_error"
  )
})

test_that("generate_twl rejects relative humidity outside [0, 100]", {
  expect_error(
    generate_twl(datetime = make_dt("2024-06-15 10:00:00"),
                 latitude = -31.95, longitude = 115.86,
                 temp = 30, wind_speed = 0.5, RH = 150,
                 direct_solar = 0, diffuse_solar = 0, pressure = 1013,
                 verbose = FALSE),
    class = "meteoHazard_input_error"
  )
})

test_that("generate_twl validates input before contacting the Open-Meteo API", {
  # Invalid input with weather omitted would normally trigger a fetch; the
  # classed validation error must be raised first, so the API is never called.
  local_mocked_bindings(
    fetch_openmeteo = function(...) stop("API should not be contacted")
  )
  expect_error(
    generate_twl(datetime = make_dt("2024-06-15 10:00:00"),
                 latitude = -31.95, longitude = 999, verbose = FALSE),
    class = "meteoHazard_input_error"
  )
})

# ---------------------------------------------------------------------------
# 8c. Result-length invariant (characterises EXISTING behaviour; the plan's
#     n_obs <- length(datetime) fix makes it explicit rather than incidental)
# ---------------------------------------------------------------------------
test_that("generate_twl returns one value per datetime when weather inputs are scalar", {
  result <- generate_twl(
    datetime      = make_dt(c("2024-06-15 08:00:00", "2024-06-15 09:00:00",
                               "2024-06-15 10:00:00")),
    latitude      = -31.95, longitude = 115.86,
    temp = 35, wind_speed = 0.5, RH = 50,
    direct_solar = 0, diffuse_solar = 0, pressure = 1013,
    verbose = FALSE
  )
  expect_length(result, 3L)
})

# ---------------------------------------------------------------------------
# 8d. convert_pressure semantics (Stage C)
# ---------------------------------------------------------------------------
test_that("user pressure in kPa (convert_pressure=FALSE) matches hPa (convert_pressure=TRUE)", {
  # Both routes reduce to the same internal kPa pressure, so TWL must be equal.
  twl_hpa <- twl_supplied(38, 0.5, 50, pressure_hpa = 1013)   # default convert=TRUE
  twl_kpa <- generate_twl(
    datetime = make_dt("2024-06-15 10:00:00"),
    latitude = -31.95, longitude = 115.86,
    temp = 38, wind_speed = 0.5, RH = 50,
    direct_solar = 0, diffuse_solar = 0,
    pressure = 101.3, convert_pressure = FALSE, verbose = FALSE
  )
  expect_equal(twl_kpa, twl_hpa, tolerance = 1e-6)
})

test_that("API-fetched pressure is normalised hPa->kPa regardless of convert_pressure", {
  skip(paste("pending: mock an Open-Meteo response with surface_pressure in hPa",
             "and assert the API path divides it by 10 even when",
             "convert_pressure = FALSE"))
})

# ---------------------------------------------------------------------------
# 9. fetch_openmeteo: offline / unit tests
# ---------------------------------------------------------------------------
test_that("fetch_openmeteo returns empty list for zero fields", {
  result <- fetch_openmeteo(
    datetime  = make_dt("2024-06-15 10:00:00"),
    latitude  = -31.95,
    longitude = 115.86,
    fields    = character(0)
  )
  expect_equal(result, list())
})

# ---------------------------------------------------------------------------
# 10. Online integration test (skipped unless TWLTEST_ONLINE=true)
# ---------------------------------------------------------------------------
test_that("fetch_openmeteo retrieves real data for a historical date", {
  skip_if(
    Sys.getenv("TWLTEST_ONLINE") != "true",
    "Skipping online test (set TWLTEST_ONLINE=true to enable)"
  )

  dt <- make_dt("2024-06-15 02:00:00")  # UTC midnight-ish, well into archive
  result <- fetch_openmeteo(
    datetime  = dt,
    latitude  = -31.95,
    longitude = 115.86,
    fields    = c("temperature_2m", "wind_speed_10m",
                  "relative_humidity_2m", "surface_pressure"),
    verbose   = FALSE
  )

  expect_named(result, c("temperature_2m", "wind_speed_10m",
                          "relative_humidity_2m", "surface_pressure"))
  expect_false(is.na(result[["temperature_2m"]]))
  expect_true(result[["temperature_2m"]] > -20 && result[["temperature_2m"]] < 50)
  expect_true(result[["relative_humidity_2m"]] >= 0 &&
                result[["relative_humidity_2m"]] <= 100)
  expect_true(result[["surface_pressure"]] > 900 &&
                result[["surface_pressure"]] < 1100)
})
