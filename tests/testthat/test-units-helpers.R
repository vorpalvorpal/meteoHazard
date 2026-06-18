# Tests for the shared units helpers (.as_units / .drop_to): the bare-or-units
# input contract, conversion, and the dimensional-mismatch guard.

test_that(".as_units tags a bare numeric as already being in the unit", {
  u <- .as_units(c(5, 10), "m/s")
  expect_s3_class(u, "units")
  expect_equal(as.numeric(u), c(5, 10))
  expect_equal(as.character(units(u)), "m/s")
})

test_that(".as_units converts a units object to the canonical unit", {
  u <- .as_units(units::set_units(c(18, 36), "km/h"), "m/s")
  expect_s3_class(u, "units")
  expect_equal(as.numeric(u), c(5, 10))
})

test_that(".as_units errors (classed) on a dimensionally incompatible units object", {
  expect_error(
    .as_units(units::set_units(5, "m/s"), "degree_C"),
    class = "meteoHazard_input_error"
  )
})

test_that(".as_units errors (classed) on a non-numeric input", {
  expect_error(.as_units("warm", "degree_C"), class = "meteoHazard_input_error")
})

test_that(".as_units passes NULL through unchanged", {
  expect_null(.as_units(NULL, "m/s"))
})

test_that(".drop_to returns a plain double in the canonical unit", {
  expect_equal(.drop_to(5, "m/s"), 5)                                  # bare assumed m/s
  expect_equal(.drop_to(units::set_units(18, "km/h"), "m/s"), 5)       # converted
  expect_false(inherits(.drop_to(5, "m/s"), "units"))
})

test_that(".drop_to preserves NA and vector length", {
  out <- .drop_to(c(1, NA, 3), "m/s")
  expect_equal(out, c(1, NA, 3))
  expect_length(out, 3)
})

test_that(".drop_to converts pressure hPa -> kPa correctly", {
  expect_equal(.drop_to(units::set_units(1013, "hPa"), "kPa"), 101.3)
  expect_equal(.drop_to(1013, "hPa"), 1013)  # bare assumed hPa
})
