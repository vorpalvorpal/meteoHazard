# ===========================================================================
# C3b benchmarks: terrain morning-pulse pathways
# Issue #17, Section 7
#
# Benchmarks:
#   - .pool_partition() + .morning_release() vectorised over 168 h
#   - odour_exposure(terrain_backend="descriptors") vs "none" at
#     168 h x 100 receptors x 3 sources
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
  # Realistic diurnal cycle: night = 8 pm – 8 am (12 h), day = 8 am – 8 pm (12 h).
  hour_in_day <- (seq_len(n_hours) - 1L) %% 24L
  is_night    <- hour_in_day < 8L | hour_in_day >= 20L

  data.frame(
    wind_speed_10m         = ifelse(is_night, runif(n_hours, 0.5, 2.5),
                                              runif(n_hours, 2.0, 7.0)),
    wind_direction_10m     = runif(n_hours, 0, 360),
    wind_gusts_10m         = ifelse(is_night, runif(n_hours, 1, 4),
                                              runif(n_hours, 3, 10)),
    boundary_layer_height  = ifelse(is_night, runif(n_hours,  80, 300),
                                              runif(n_hours, 500, 1500)),
    temperature_2m         = 15 + 5 * sin(2 * pi * hour_in_day / 24 - pi / 2) +
                               rnorm(n_hours, 0, 0.5),
    pressure_msl           = 1013 + rnorm(n_hours, 0, 2),
    precipitation          = pmax(0, rnorm(n_hours, 0, 0.2)),
    relative_humidity_2m   = ifelse(is_night, runif(n_hours, 70, 95),
                                              runif(n_hours, 40, 70)),
    cloud_cover            = runif(n_hours, 0, 60),
    direct_radiation       = ifelse(is_night, 0,
                                    pmax(10, rnorm(n_hours, 350, 80))),
    soil_moisture_0_to_1cm = runif(n_hours, 0.05, 0.3),
    soil_moisture_1_to_3cm = runif(n_hours, 0.05, 0.3),
    wind_speed_80m         = ifelse(is_night, runif(n_hours, 1.0, 4.0),
                                              runif(n_hours, 3.0, 9.0)),
    wind_direction_80m     = runif(n_hours, 0, 360),
    wind_speed_120m        = ifelse(is_night, runif(n_hours, 1.5, 5.0),
                                              runif(n_hours, 3.5, 10.0)),
    wind_direction_120m    = runif(n_hours, 0, 360),
    wind_speed_180m        = ifelse(is_night, runif(n_hours, 2.0, 6.0),
                                              runif(n_hours, 4.0, 11.0)),
    wind_direction_180m    = runif(n_hours, 0, 360)
  )
}

met_168 <- make_met(168)

# Derive the ventilation state (needed to extract pool_top / cbl_growth
# for the pathway-primitive benchmarks).
terrain_valley <- mh_terrain(
  valley_depth     = 80,
  taf              = 1.4,
  relief           = 30,
  drainage_bearing = 45,
  flow_convergence = 0.2,
  shelter_index    = 30
)
vs_168 <- ventilation_state(met_168, terrain = terrain_valley)


# ---------------------------------------------------------------------------
# Site fixture: 3 sources (heights: 3 m, 10 m, 20 m), 100 receptors
# ---------------------------------------------------------------------------

make_odour_site <- function(n_sources = 3, n_receptors = 100,
                            emit_heights = c(3, 10, 20)) {
  set.seed(42)
  src_ids <- paste0("src", seq_len(n_sources))
  rec_ids <- paste0("rec", seq_len(n_receptors))

  src_geom <- lapply(seq_len(n_sources), function(k) {
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
    emit_height = c(emit_heights[seq_len(n_sources)], rep(NA_real_, n_receptors))
  )
  mh_site(features = feats, roles = roles, terrain = terrain_valley, epsg = 32755L)
}

site_terrain_3src_100rec <- make_odour_site()
# Matching flat site for the "none" comparison (no terrain attached).
make_flat_site <- function(n_sources = 3, n_receptors = 100,
                           emit_heights = c(3, 10, 20)) {
  set.seed(42)
  src_ids <- paste0("src", seq_len(n_sources))
  rec_ids <- paste0("rec", seq_len(n_receptors))
  src_geom <- lapply(seq_len(n_sources), function(k) sf::st_point(c(50 * (k - 1), 0)))
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
    emit_height = c(emit_heights[seq_len(n_sources)], rep(NA_real_, n_receptors))
  )
  mh_site(features = feats, roles = roles, epsg = 32755L)
}

site_flat_3src_100rec <- make_flat_site()


# ---------------------------------------------------------------------------
# --- Section 1: pathway primitives over 168 h ------------------------------
# ---------------------------------------------------------------------------

message("Benchmarking: .pool_partition() over 168 h")

# Representative single-source scenario: emit_height = 10 m.
pool_top_vec <- vs_168$pool_top
emit_height  <- 10

bm_partition <- bench::mark(
  # .pool_partition() is already vectorised over hours; call it once for all 168 h.
  `pool_partition_168h` = {
    pt_safe <- ifelse(is.na(pool_top_vec), 0, pool_top_vec)
    meteoHazard:::.pool_partition(
      emit_extent = rep(emit_height, length(pt_safe)),
      pool_top    = pt_safe,
      delta       = rep(20, length(pt_safe))
    )
  },
  iterations = 20,
  check = FALSE
)
print(bm_partition[, c("expression", "min", "median", "mem_alloc", "n_itr")])

message("Benchmarking: .morning_release() over 168 h")

bm_release <- bench::mark(
  `morning_release_168h` = meteoHazard:::.morning_release(
    pool_top   = vs_168$pool_top,
    cbl_growth = vs_168$cbl_growth,
    is_day     = vs_168$is_day
  ),
  iterations = 20,
  check = FALSE
)
print(bm_release[, c("expression", "min", "median", "mem_alloc", "n_itr")])


# ---------------------------------------------------------------------------
# --- Section 2: terrain backend vs none at 168 h x 100 recs x 3 sources ---
# ---------------------------------------------------------------------------

message("Benchmarking: odour_exposure() descriptors vs none (168 h x 100 recs x 3 src)")

bm_backend <- bench::mark(
  `exp_none_168h`        = odour_exposure(met_168, site_flat_3src_100rec,
                                          terrain_backend = "none"),
  `exp_descriptors_168h` = odour_exposure(met_168, site_terrain_3src_100rec,
                                          terrain_backend = "descriptors"),
  iterations = 20,
  check = FALSE
)
print(bm_backend[, c("expression", "min", "median", "mem_alloc", "n_itr")])
cat("Terrain pathway overhead (median delta):",
    format(median(bm_backend$time[[2]]) - median(bm_backend$time[[1]])), "\n")
