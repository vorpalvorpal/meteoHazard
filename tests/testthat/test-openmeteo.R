# Tests for fetch_openmeteo(): endpoint selection/splitting and request metadata.
# The network is never contacted — the spanning-range and single-request specs
# mock the HTTP layer (owned by the implement step, which adds the seam).

make_dt <- function(s, tz = "UTC") as.POSIXct(s, tz = tz)

test_that("fetch_openmeteo returns an empty list when no fields are requested", {
  expect_equal(
    fetch_openmeteo(make_dt("2024-06-15 10:00:00"), -31.95, 115.86, character(0)),
    list()
  )
})

test_that("fetch_openmeteo identifies itself as the meteoHazard package", {
  # om_request() builds the httr2 request without performing it; the user-agent
  # is stored on the request and must name the package.
  req <- om_request("https://api.open-meteo.com/v1/forecast?x=1")
  expect_match(req$options$useragent, "meteoHazard")
})

test_that("fetch_openmeteo issues a single request for a range within one endpoint", {
  # An all-historical range (well before today-92) is served by the archive
  # endpoint alone -> exactly one perform call.
  calls <- 0L
  local_mocked_bindings(
    om_perform = function(req) {
      calls <<- calls + 1L
      list(hourly = list(time = "2020-06-15T10:00", temperature_2m = 21.5))
    }
  )
  res <- fetch_openmeteo(make_dt("2020-06-15 10:00:00"), -31.95, 115.86,
                         "temperature_2m")
  expect_equal(calls, 1L)
  expect_equal(res$temperature_2m, 21.5)
})

test_that("fetch_openmeteo splits the request across forecast and archive for a spanning range", {
  # A range straddling today-92 needs the archive endpoint for the old hour and
  # the forecast endpoint for the recent hour; both must resolve to non-NA at
  # their input positions.
  calls <- 0L
  old_hour    <- make_dt("2020-06-15 10:00:00")  # archive side
  recent_hour <- as.POSIXct(
    format(Sys.time() - 86400, "%Y-%m-%d 10:00:00", tz = "UTC"), tz = "UTC"
  )                                              # yesterday -> forecast side
  local_mocked_bindings(
    om_perform = function(req) {
      calls <<- calls + 1L
      if (grepl("archive", req$url)) {
        list(hourly = list(time = "2020-06-15T10:00", temperature_2m = 21.5))
      } else {
        list(hourly = list(
          time = format(recent_hour, "%Y-%m-%dT%H:%M", tz = "UTC"),
          temperature_2m = 14.0
        ))
      }
    }
  )
  res <- fetch_openmeteo(c(old_hour, recent_hour), -31.95, 115.86,
                         "temperature_2m")
  expect_equal(calls, 2L)
  expect_false(any(is.na(res$temperature_2m)))
  expect_equal(res$temperature_2m, c(21.5, 14.0))
})
