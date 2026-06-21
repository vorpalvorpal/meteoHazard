# ===========================================================================
# C3a benchmarks: odour_exposure() flat backend
# Issue #16, Section 7
#
# Benchmarks:
#   - odour_exposure() over 168 h x {10, 100} receptors x {1, 3} sources
#     (the (t, j, k) array build is the cost centre)
#   - flat-backend cost vs a single-source 1-receptor baseline
#
# Run interactively; results printed to console.
# ===========================================================================

library(meteoHazard)
library(bench)
library(sf)

# ---------------------------------------------------------------------------
# Met data fixture
# ---------------------------------------------------------------------------

make_met <- function(n_hours) {
  set.seed(42)
  data.frame(
    wind_speed_10m         = runif(n_hours, 1, 8),
    wind_direction_10m     = runif(n_hours, 0, 360),
    wind_gusts_10m         = runif(n_hours, 2, 12),
    boundary_layer_height  = runif(n_hours, 200, 1500),
    temperature_2m         = runif(n_hours, 10, 30),
    pressure_msl           = runif(n_hours, 1000, 1020),
    precipitation          = pmax(0, rnorm(n_hours, 0, 0.5)),
    relative_humidity_2m   = runif(n_hours, 40, 90),
    cloud_cover            = runif(n_hours, 0, 100),
    direct_radiation       = pmax(0, rnorm(n_hours, 300, 200)),
    soil_moisture_0_to_1cm = runif(n_hours, 0.05, 0.3),
    soil_moisture_1_to_3cm = runif(n_hours, 0.05, 0.3),
    wind_speed_80m         = runif(n_hours, 1, 10),
    wind_direction_80m     = runif(n_hours, 0, 360),
    wind_speed_120m        = runif(n_hours, 1, 10),
    wind_direction_120m    = runif(n_hours, 0, 360),
    wind_speed_180m        = runif(n_hours, 1, 10),
    wind_direction_180m    = runif(n_hours, 0, 360)
  )
}

met_168 <- make_met(168)


# ---------------------------------------------------------------------------
# Site fixture builder: n_sources sources, n_receptors receptors
# ---------------------------------------------------------------------------

make_odour_site <- function(n_sources = 1, n_receptors = 1) {
  set.seed(42)

  src_ids <- paste0("src", seq_len(n_sources))
  rec_ids <- paste0("rec", seq_len(n_receptors))

  src_geom <- lapply(seq_len(n_sources), function(k) {
    # Spread sources slightly so they are distinct
    sf::st_point(c(50 * (k - 1), 0))
  })
  rec_geom <- lapply(seq_len(n_receptors), function(k) {
    sf::st_point(c(500 * cos(2 * pi * k / n_receptors),
                   500 * sin(2 * pi * k / n_receptors)))
  })

  all_ids  <- c(src_ids, rec_ids)
  all_geom <- sf::st_sfc(c(src_geom, rec_geom), crs = 32755)

  feats <- sf::st_sf(id = all_ids, geometry = all_geom)
  roles <- data.frame(
    feature_id  = all_ids,
    hazard      = "odour",
    role        = c(rep("source", n_sources), rep("receptor", n_receptors)),
    emit_height = c(rep(5, n_sources), rep(NA_real_, n_receptors))
  )
  mh_site(features = feats, roles = roles, epsg = 32755L)
}

# Pre-build all site variants (construction cost excluded from timings).
site_1src_10rec  <- make_odour_site(n_sources = 1, n_receptors = 10)
site_1src_100rec <- make_odour_site(n_sources = 1, n_receptors = 100)
site_3src_10rec  <- make_odour_site(n_sources = 3, n_receptors = 10)
site_3src_100rec <- make_odour_site(n_sources = 3, n_receptors = 100)

# Baseline: single source, single receptor.
site_1src_1rec   <- make_odour_site(n_sources = 1, n_receptors = 1)


# ---------------------------------------------------------------------------
# --- Section 1: 168 h x receptor count x source count --------------------
# ---------------------------------------------------------------------------

message("Benchmarking: odour_exposure() 168 h x {10,100} receptors x {1,3} sources")

bm_exposure <- bench::mark(
  `exp_1src_1rec_168h`   = odour_exposure(met_168, site_1src_1rec,   terrain_backend = "none"),
  `exp_1src_10rec_168h`  = odour_exposure(met_168, site_1src_10rec,  terrain_backend = "none"),
  `exp_1src_100rec_168h` = odour_exposure(met_168, site_1src_100rec, terrain_backend = "none"),
  `exp_3src_10rec_168h`  = odour_exposure(met_168, site_3src_10rec,  terrain_backend = "none"),
  `exp_3src_100rec_168h` = odour_exposure(met_168, site_3src_100rec, terrain_backend = "none"),
  iterations = 20,
  check = FALSE
)
print(bm_exposure[, c("expression", "min", "median", "mem_alloc", "n_itr")])


# ---------------------------------------------------------------------------
# --- Section 2: Scaling summary -------------------------------------------
# ---------------------------------------------------------------------------

message("Scaling summary (median times relative to 1-src 1-rec baseline):")

baseline <- median(bm_exposure$time[[1]])
ratios <- sapply(bm_exposure$time, function(t) median(t) / baseline)
print(data.frame(
  expression = as.character(bm_exposure$expression),
  median_ms  = round(sapply(bm_exposure$time, function(t) median(t) * 1e3), 2),
  ratio_vs_baseline = round(ratios, 2)
))
