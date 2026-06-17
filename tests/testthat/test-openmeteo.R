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

test_that("fetch_openmeteo splits the request across forecast and archive for a spanning range", {
  # A range straddling today-92 needs the archive endpoint for the old end and
  # the forecast endpoint for the recent end; both must resolve to non-NA.
  skip(paste("pending: mock two Open-Meteo responses (archive + forecast); for a",
             "datetime vector spanning today-92, assert an archival hour and a",
             "recent hour both return non-NA values aligned to their positions"))
})

test_that("fetch_openmeteo issues a single request for a range within one endpoint", {
  skip(paste("pending: mock the HTTP layer and assert exactly one request is",
             "performed for an all-historical (single-endpoint) range"))
})

test_that("fetch_openmeteo identifies itself as the meteoHazard package", {
  skip("pending: assert the outgoing request User-Agent string contains 'meteoHazard'")
})
