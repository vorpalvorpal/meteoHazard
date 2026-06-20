# Tests for odour_exposure(): Gaussian distance/direction geometry, the 0-100
# map, calm handling, and the relocated drainage/fumigation refinement.

# met_data builder (exposure-relevant columns), defaults overridable by name.
me <- function(n = 1, ...) {
  base <- list(
    wind_direction_10m    = 180, # blows from S -> downwind N (bearing 0)
    wind_speed_10m        = 3,
    direct_radiation      = 0,
    cloud_cover           = 50,
    boundary_layer_height = 500
  )
  ov <- list(...)
  base[names(ov)] <- ov
  as.data.frame(lapply(base, rep_len, n))
}
rcp <- function(bearing, distance) data.frame(bearing = bearing, distance = distance)

# ── Distance decay ──────────────────────────────────────────────────────────
test_that("a nearer on-axis receptor is more exposed than a far one", {
  near <- odour_exposure(0.3, me(), rcp(0, 250))
  far  <- odour_exposure(0.3, me(), rcp(0, 1500))
  expect_gt(near, far)
})

# ── Direction ───────────────────────────────────────────────────────────────
test_that("exposure falls as the plume swings off the receptor bearing", {
  # Receptor due north (bearing 0); vary wind so downwind bearing is 0/20/90.
  d <- me(n = 3, wind_direction_10m = c(180, 160, 270))
  ex <- odour_exposure(rep(0.3, 3), d, rcp(0, 400))
  expect_gt(ex[1], ex[2])
  expect_gt(ex[2], ex[3])
  expect_lt(ex[3], 0.5) # ~90 deg off-axis -> negligible
})

test_that("a receptor directly upwind of the source is not exposed", {
  # Wind blows toward N (downwind bearing 180); receptor due N (bearing 0) is
  # upwind of the source and must not be in the plume.
  ex <- odour_exposure(0.3, me(wind_direction_10m = 0), rcp(0, 400))
  expect_lt(ex, 0.01)
})

test_that("forecast-direction uncertainty keeps a small off-axis miss exposed", {
  d <- me(n = 2, wind_direction_10m = c(180, 170)) # on-axis vs 10 deg off
  ex <- odour_exposure(rep(0.3, 2), d, rcp(0, 400))
  expect_lt(ex[2], ex[1])
  expect_gt(ex[2], 0.4 * ex[1]) # sigma_fc ~ 12 deg keeps it substantial
})

# ── Worst case across receptors ─────────────────────────────────────────────
test_that("the worst-affected receptor sets the hourly exposure", {
  recs <- rcp(c(0, 180), c(400, 400)) # one aligned, one opposite
  combined <- odour_exposure(0.3, me(), recs)
  aimed    <- odour_exposure(0.3, me(), rcp(0, 400))
  expect_equal(combined, aimed)
})

# ── Monotonic in hazard, and bounds ─────────────────────────────────────────
test_that("exposure increases with hazard and stays within [0, 100]", {
  lo <- odour_exposure(0.05, me(), rcp(0, 400))
  hi <- odour_exposure(2.0, me(), rcp(0, 400))
  expect_gt(hi, lo)
  expect_gte(lo, 0)
  expect_lte(hi, 100)
})

# ── Calm: direction-agnostic ────────────────────────────────────────────────
test_that("under calm winds exposure is the same for any bearing at one distance", {
  d <- me(wind_speed_10m = 0.2) # calm
  north <- odour_exposure(0.3, d, rcp(0, 400))
  south <- odour_exposure(0.3, d, rcp(180, 400))
  expect_equal(north, south)
})

# ── Drainage refinement ─────────────────────────────────────────────────────
test_that("katabatic drainage confines emissions away from an aligned receptor", {
  # Calm, dark, clear hour -> drainage. Aligned receptor (bearing 0, drainage
  # axis from 0). The directional factor collapses to 0.3 * ~0.05 = ~0.015
  # (the 0.3 carries the monolith's W_spd confinement) versus the generic calm
  # 0.5 -- a >10x suppression, not just "lower".
  d   <- me(wind_speed_10m = 0.2, direct_radiation = 0, cloud_cover = 10)
  rec <- rcp(0, 300)
  ax  <- data.frame(bearing_from = 0, weight = 1)
  with_drn <- odour_exposure(0.3, d, rec, drainage_axes = ax)
  without  <- odour_exposure(0.3, d, rec, drainage_axes = NULL)
  expect_lt(with_drn, 0.1 * without)
})

test_that("morning fumigation lofts the overnight pool toward aligned receptors", {
  # 6 h drainage night, then a sunny morning hour with the wind pointing AWAY
  # from the receptor. Without the terrain module the receptor sees ~0; with it,
  # the fumigation floor + entrainment boost lift the morning exposure.
  d <- me(
    n = 7,
    wind_speed_10m        = c(rep(0.2, 6), 2),
    direct_radiation      = c(rep(0, 6), 100),
    cloud_cover           = c(rep(10, 6), 20),
    wind_direction_10m    = c(rep(0, 6), 0), # morning wind blows toward N -> away from bearing-0 receptor
    boundary_layer_height = c(rep(100, 6), 500)
  )
  rec <- rcp(0, 300)
  ax  <- data.frame(bearing_from = 0, weight = 1)
  with_drn <- odour_exposure(rep(0.3, 7), d, rec, drainage_axes = ax)
  without  <- odour_exposure(rep(0.3, 7), d, rec, drainage_axes = NULL)
  expect_lt(without[7], 1)            # wind away, no terrain module -> ~0
  expect_gt(with_drn[7], without[7])  # fumigation lifts it
})

# ── Validation ──────────────────────────────────────────────────────────────
test_that("a hazard/met_data length mismatch is a classed error", {
  expect_error(
    odour_exposure(rep(0.3, 2), me(n = 3), rcp(0, 400)),
    class = "meteoHazard_input_error"
  )
})

test_that("missing receptor columns raise a classed error", {
  expect_error(
    odour_exposure(0.3, me(), data.frame(bearing = 0)),
    class = "meteoHazard_input_error"
  )
})

test_that("missing met_data columns raise a classed error", {
  expect_error(
    odour_exposure(0.3, me()[, -1], rcp(0, 400)),
    class = "meteoHazard_input_error"
  )
})

# ── Units handling ───────────────────────────────────────────────────────────
test_that("odour_exposure accepts units-tagged met columns and receptor distance", {
  bare <- odour_exposure(0.3, me(), rcp(0, 250))
  d <- me()
  d$wind_speed_10m <- units::set_units(10.8, "km/h") # 3 m/s
  tagged <- odour_exposure(
    0.3, d,
    data.frame(bearing = 0, distance = units::set_units(0.25, "km")) # 250 m
  )
  expect_equal(tagged, bare, tolerance = 1e-6)
})

test_that("odour_exposure rejects a receptor distance tagged with incompatible units", {
  expect_error(
    odour_exposure(0.3, me(),
                   data.frame(bearing = 0, distance = units::set_units(250, "degree_C"))),
    class = "meteoHazard_input_error"
  )
})

test_that("odour_exposure returns a plain numeric 0-100 band", {
  out <- odour_exposure(0.3, me(), rcp(0, 250))
  expect_false(inherits(out, "units"))
  expect_type(out, "double")
})

# === C3a behaviour spec (issue #16): odour_exposure() on mh_site, area + multi-source ===
# Supersedes the test_that() blocks above when the rewrite lands. Pending here.
describe("odour_exposure() on the mh_site model (terrain backend 'none')", {
  # preserved behaviours, on the new API:
  it("makes a nearer on-axis receptor more exposed than a far one")
  it("falls as the plume swings off the receptor bearing")
  it("reads ~0 for a receptor directly upwind of the source")
  it("keeps a small off-axis miss substantial via forecast-direction uncertainty")
  it("is direction-agnostic at a fixed distance under calm winds")
  it("returns the worst-affected receptor's value for the hour")
  it("increases with hazard and stays within [0, 100]")
  it("returns a plain numeric band")
  it("errors (classed) on an invalid mh_site or a hazard/met length mismatch")
  # area source:
  it("stays finite for a receptor at the source edge (no point-source divergence)")
  it("tends to the initial spread sigma_y0 as distance goes to zero")
  # multiple sources:
  it("sums concentrations from multiple sources at a receptor")
  it("gives ~twice the single-source concentration for two identical co-located sources")
  it("applies the 0-100 map to the summed concentration, so two moderate sources can be severe")
  it("reduces to the current max-over-receptors for a single source")
  # edges:
  it("produces no morning pulse with terrain_backend = 'none'")
  it("skips a receptor coincident with the source")
})
