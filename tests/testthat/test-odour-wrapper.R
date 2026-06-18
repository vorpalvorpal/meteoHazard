# Tests for odour_risk() (the combined wrapper) and the optional
# consecutive-hourly datetime guard on odour_hazard().

mw <- function(n = 6, ...) {
  base <- list(
    wind_direction_10m     = 180,
    wind_speed_10m         = 3,
    direct_radiation       = 0,
    cloud_cover            = 50,
    boundary_layer_height  = 400,
    temperature_2m         = 18,
    pressure_msl           = 1012,
    precipitation          = 0,
    relative_humidity_2m   = 70,
    soil_moisture_0_to_1cm = 0.1,
    soil_moisture_1_to_3cm = 0.1
  )
  ov <- list(...)
  base[names(ov)] <- ov
  as.data.frame(lapply(base, rep_len, n))
}
rcp <- function(bearing, distance) data.frame(bearing = bearing, distance = distance)

test_that("the wrapper equals exposure(hazard(...)) composed by hand", {
  d   <- mw()
  rec <- rcp(c(0, 90), c(300, 700))
  combined <- odour_risk(d, rec)
  manual   <- odour_exposure(odour_hazard(d), d, rec)
  expect_equal(combined, manual)
})

test_that("the wrapper returns one 0-100 value per row", {
  ex <- odour_risk(mw(n = 12), rcp(0, 400))
  expect_length(ex, 12)
  expect_true(all(ex >= 0 & ex <= 100))
})

test_that("a consecutive-hourly datetime passes the spacing guard silently", {
  d  <- mw(n = 6)
  dt <- as.POSIXct("2024-06-01 00:00", tz = "UTC") + 3600 * (0:5)
  expect_no_warning(odour_hazard(d, datetime = dt))
})

test_that("an irregular datetime warns but still computes", {
  d  <- mw(n = 6)
  dt <- as.POSIXct("2024-06-01 00:00", tz = "UTC") + 3600 * c(0, 1, 2, 4, 5, 6) # gap
  expect_warning(h <- odour_hazard(d, datetime = dt), "consecutive hourly")
  expect_length(h, 6)
})

test_that("a non-POSIXct datetime is a classed input error", {
  expect_error(
    odour_hazard(mw(), datetime = 1:6),
    class = "meteoHazard_input_error"
  )
})
