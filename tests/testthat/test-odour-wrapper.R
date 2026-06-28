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

# Build a minimal mh_site for wrapper tests.
.make_wrapper_site <- function() {
  feats <- sf::st_sf(
    id       = c("src", "rec1", "rec2"),
    geometry = sf::st_sfc(
      sf::st_point(c(335000,       6250000)),
      sf::st_point(c(335000,       6250400)),   # 400 m north
      sf::st_point(c(335000 + 700, 6250000)),   # 700 m east
      crs = 32755
    )
  )
  roles <- data.frame(
    feature_id = c("src", "rec1", "rec2"),
    hazard     = "odour",
    role       = c("source", "receptor", "receptor"),
    stringsAsFactors = FALSE
  )
  mh_site(feats, roles, epsg = 32755L)
}

# ── datetime guard (these call odour_hazard() directly, unchanged) ───────────

test_that("a consecutive-hourly datetime passes the spacing guard silently", {
  d  <- mw(n = 6)
  dt <- as.POSIXct("2024-06-01 00:00", tz = "UTC") + 3600 * (0:5)
  expect_no_warning(odour_hazard(d, datetime = dt))
})

test_that("an irregular datetime warns but still computes", {
  d  <- mw(n = 6)
  dt <- as.POSIXct("2024-06-01 00:00", tz = "UTC") + 3600 * c(0, 1, 2, 4, 5, 6)
  expect_warning(h <- odour_hazard(d, datetime = dt), "consecutive hourly")
  expect_length(h, 6)
})

test_that("a non-POSIXct datetime is a classed input error", {
  expect_error(
    odour_hazard(mw(), datetime = 1:6),
    class = "meteoHazard_input_error"
  )
})

# === C3a behaviour spec (issue #16): odour_risk() on the new signatures ===
describe("odour_risk() on the mh_site model", {

  it("equals odour_exposure() applied to the ventilation/hazard composed by hand", {
    d    <- mw(n = 6)
    site <- .make_wrapper_site()
    risk_out     <- odour_risk(d, site)
    exposure_out <- odour_exposure(d, site)
    expect_equal(risk_out, exposure_out)
  })

  it("returns a per-receptor relative-concentration matrix (one row per hour)", {
    d    <- mw(n = 12)
    site <- .make_wrapper_site()
    out  <- odour_risk(d, site)
    expect_true(is.matrix(out))
    expect_equal(nrow(out), 12)            # one row per forecast hour
    expect_equal(ncol(out), 2)             # two receptors
    expect_true(all(out >= 0))             # relative concentration, unbounded above
  })
})
