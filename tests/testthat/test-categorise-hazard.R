# Tests for the hazard-tier categorise / colour family (litter, dust, odour),
# which mirror the categorise_twl() / twl_colour() template on a shared
# green/amber/orange/red palette.

test_that("categorise_litter maps the documented bands and boundaries", {
  # Cut-points 20/45/70; "< x" means exactly x falls in the higher band.
  expect_equal(
    categorise_litter(c(0, 19.9, 20, 44.9, 45, 69.9, 70, 100)),
    c("LOW", "LOW", "MODERATE", "MODERATE", "HIGH", "HIGH", "EXTREME", "EXTREME")
  )
  expect_equal(categorise_litter(NA_real_), NA_character_)
})

test_that("categorise_dust maps the documented bands and boundaries", {
  # Cut-points 25/50/75.
  expect_equal(
    categorise_dust(c(0, 24.9, 25, 49.9, 50, 74.9, 75, 100)),
    c("LOW", "LOW", "MODERATE", "MODERATE", "HIGH", "HIGH", "EXTREME", "EXTREME")
  )
  expect_equal(categorise_dust(NA_real_), NA_character_)
})

test_that("categorise_odour maps the documented bands and boundaries", {
  # Cut-points 15/40/70; top band is labelled VERY HIGH (not EXTREME).
  expect_equal(
    categorise_odour(c(0, 14.9, 15, 39.9, 40, 69.9, 70, 100)),
    c("LOW", "LOW", "MODERATE", "MODERATE", "HIGH", "HIGH", "VERY HIGH", "VERY HIGH")
  )
  expect_equal(categorise_odour(NA_real_), NA_character_)
})

test_that("the colour helpers share one palette and align with their tiers", {
  vals <- c(5, 30, 55, 95, NA)

  # All three use the same green/amber/orange/red + grey scheme.
  expect_equal(litter_colour(c(10, 30, 55, 85, NA)),
               c("#4CAF50", "#FFC107", "#FF9800", "#D32F2F", "#CCCCCC"))
  expect_equal(dust_colour(c(10, 35, 60, 90, NA)),
               c("#4CAF50", "#FFC107", "#FF9800", "#D32F2F", "#CCCCCC"))
  expect_equal(odour_colour(c(10, 30, 55, 85, NA)),
               c("#4CAF50", "#FFC107", "#FF9800", "#D32F2F", "#CCCCCC"))

  # NA -> grey for every helper.
  expect_equal(litter_colour(vals)[5], "#CCCCCC")
  expect_equal(dust_colour(vals)[5], "#CCCCCC")
  expect_equal(odour_colour(vals)[5], "#CCCCCC")
})

test_that("the categorise / colour helpers are vectorised and length-preserving", {
  x <- c(5, 22, 48, 72, 99)
  expect_length(categorise_litter(x), length(x))
  expect_length(dust_colour(x), length(x))
  expect_length(categorise_odour(x), length(x))
})
