# Tests for the shared odour dispersion-state helper: Pasquill-Turner stability
# (primary), the optional shear override, the calm-wind floor, and the
# mixing-depth fallbacks. The stability index s runs 0 (class A) to 5 (class F).

# Minimal met_data builder; recycles scalars to a common length.
md <- function(wind_speed_10m, direct_radiation, cloud_cover,
               boundary_layer_height = NA_real_, wind_speed_80m = NA_real_) {
  n <- max(lengths(list(wind_speed_10m, direct_radiation, cloud_cover,
                        boundary_layer_height, wind_speed_80m)))
  data.frame(
    wind_speed_10m        = rep_len(wind_speed_10m, n),
    direct_radiation      = rep_len(direct_radiation, n),
    cloud_cover           = rep_len(cloud_cover, n),
    boundary_layer_height = rep_len(boundary_layer_height, n),
    wind_speed_80m        = rep_len(wind_speed_80m, n)
  )
}

# ── Pasquill-Turner: daytime ────────────────────────────────────────────────
test_that("strong insolation + light wind gives a very unstable class (s ~ 0)", {
  st <- .odour_dispersion_state(md(1.0, 800, 0))
  expect_equal(st$s, 0)
  expect_true(st$is_day)
})

test_that("slight insolation + strong wind gives near-neutral (s = 3, class D)", {
  st <- .odour_dispersion_state(md(7.0, 100, 50))
  expect_equal(st$s, 3)
})

# ── Pasquill-Turner: nighttime ──────────────────────────────────────────────
test_that("clear calm-ish night gives a very stable class (s = 5, class F)", {
  st <- .odour_dispersion_state(md(1.5, 0, 10))
  expect_equal(st$s, 5)
  expect_false(st$is_day)
})

test_that("overcast moderate-wind night is near-neutral (s = 3, class D)", {
  st <- .odour_dispersion_state(md(4.0, 0, 90))
  expect_equal(st$s, 3)
})

test_that("night stability increases as wind drops (5-6 -> 3-5 -> 2-3, clear)", {
  s_strong <- .odour_dispersion_state(md(5.5, 0, 10))$s # D = 3
  s_mid    <- .odour_dispersion_state(md(4.0, 0, 10))$s # E = 4
  s_light  <- .odour_dispersion_state(md(2.5, 0, 10))$s # F = 5
  expect_lt(s_strong, s_mid)
  expect_lt(s_mid, s_light)
})

# ── Calm-wind override ──────────────────────────────────────────────────────
test_that("calm wind forces the s = 4.25 stability override", {
  expect_equal(.odour_dispersion_state(md(0.2, 0, 10))$s, 4.25)
  expect_true(.odour_dispersion_state(md(0.2, 0, 10))$is_calm)
})

test_that("NA wind is treated as calm (s = 4.25)", {
  expect_equal(.odour_dispersion_state(md(NA_real_, 0, 10))$s, 4.25)
})

# ── Effective wind floor ────────────────────────────────────────────────────
test_that("u_eff is floored at U_CALM_FLOOR and otherwise passes through", {
  expect_equal(.odour_dispersion_state(md(0.2, 0, 10))$u_eff, 0.5)
  expect_equal(.odour_dispersion_state(md(3.0, 0, 50))$u_eff, 3.0)
  expect_equal(.odour_dispersion_state(md(NA_real_, 0, 10))$u_eff, 0.5)
})

# ── Mixing-depth fallbacks ──────────────────────────────────────────────────
test_that("h_mix passes through when boundary_layer_height is present", {
  expect_equal(.odour_dispersion_state(md(4, 0, 50, boundary_layer_height = 750))$h_mix, 750)
})

test_that("NA boundary layer falls back by stability: 200 stable/calm, 600 otherwise", {
  expect_equal(.odour_dispersion_state(md(0.2, 0, 10))$h_mix, 200) # calm
  expect_equal(.odour_dispersion_state(md(1.5, 0, 10))$h_mix, 200) # stable night (s = 5)
  expect_equal(.odour_dispersion_state(md(1.0, 800, 0))$h_mix, 600) # unstable day (s = 0)
})

# ── Shear override ──────────────────────────────────────────────────────────
test_that("shear override yields a stable class for strong shear", {
  st <- .odour_dispersion_state(md(1.0, 0, 50, wind_speed_80m = 4.0),
                                stability = "shear")
  expect_equal(st$s, 5) # alpha = ln(4)/ln(8) ~ 0.667 > 0.40 -> class F
})

test_that("shear override falls back to s = 4.25 when 80 m wind is missing", {
  st <- .odour_dispersion_state(md(2.0, 0, 50, wind_speed_80m = NA_real_),
                                stability = "shear")
  expect_equal(st$s, 4.25)
})

# ── Vectorisation ───────────────────────────────────────────────────────────
test_that("the helper is vectorised and returns aligned vectors", {
  st <- .odour_dispersion_state(md(c(1.0, 4.0, 0.2), c(800, 0, 0), c(0, 90, 10)))
  expect_length(st$s, 3)
  expect_equal(st$s, c(0, 3, 4.25))
  expect_equal(st$is_day, c(TRUE, FALSE, FALSE))
})
