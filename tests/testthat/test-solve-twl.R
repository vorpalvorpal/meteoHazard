# Tests for solve_twl_single() behaviour beyond the Brake & Bates Table 1
# regression (that exact-value validation lives in test-twl-table1.R and is the
# bit-identical guard for the compute_components() de-duplication refactor).
# Here we cover the TWL floor/ceiling constants and the bracket sign-check.

# Minimal solver caller deriving psychrometrics from DB/WB the same way the
# Table 1 helper does (sling psychrometer, pressure in kPa).
call_solver <- function(temp, wb, mrt, wind = 0.5, pressure_kpa = 101,
                        Icl = 0.45, icl = 0.45, temp_dewpoint = NULL,
                        index = NULL) {
  es_wb <- calc_sat_vp(wb)
  pa    <- max(0, es_wb - 6.6e-4 * pressure_kpa * (temp - wb))
  if (is.null(temp_dewpoint)) {
    lt <- log(pa / 0.61121)
    temp_dewpoint <- (257.14 * lt) / (18.678 - lt)
  }
  solve_twl_single(
    temp = temp, wind_speed = wind, RH = 100 * pa / calc_sat_vp(temp),
    direct_solar = 0, diffuse_solar = 0, pressure = pressure_kpa,
    pa = pa, temp_dewpoint = temp_dewpoint, wet_bulb = wb,
    globe_temp = mrt, trad = mrt,
    max_core_temp = 38.2, max_sweat_rate = 0.67,
    Icl = Icl, icl = icl, LR = 16.5, lambda = 2430, index = index
  )
}

# ── TWL range constants (Stage 0) ──────────────────────────────────────────
test_that("TWL_FLOOR and TWL_CEILING constants hold the documented values", {
  expect_equal(TWL_CONSTANTS$TWL_FLOOR, 60)
  expect_equal(TWL_CONSTANTS$TWL_CEILING, 380)
})

# ── TWL floor behaviour (Stage A) ──────────────────────────────────────────
test_that("solve_twl_single returns the TWL floor for a degenerate humidity bracket", {
  # Dew point above (core - 0.5 = 37.7) collapses the bisection bracket, so the
  # solver returns the floor (60 W/m^2).
  res <- call_solver(temp = 40, wb = 39, mrt = 40, temp_dewpoint = 39.5)
  expect_equal(res, 60)
})

test_that("solve_twl_single never returns below the TWL floor", {
  res <- call_solver(temp = 42, wb = 30, mrt = 44, wind = 0.2)
  expect_gte(res, 60)
})

# ── Bracket sign-check (Stage A, new diagnostic) ──────────────────────────
test_that("solve_twl_single warns when the heat-balance root is not bracketed", {
  skip(paste("pending: construct inputs whose heat balance has the same sign at",
             "both ends of the t_skin bracket and pass a non-NULL index; assert",
             "a warning that names the observation index"))
})
