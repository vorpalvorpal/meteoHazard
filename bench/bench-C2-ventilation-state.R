# ===========================================================================
# C2 benchmarks: ventilation_state() + pool_top + odour_hazard()
# Issue #15, Section 7
#
# Benchmarks:
#   - ventilation_state() over realistic horizons (48 h, 168 h)
#   - pool_top heat-deficit at 168 h specifically
#   - odour_hazard() at 168 h (run-path cost)
#
# Run interactively; results printed to console.
# ===========================================================================

library(meteoHazard)
library(bench)

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

met_48  <- make_met(48)
met_168 <- make_met(168)

# Terrain with a modest valley (exercises pool_top TAF + basin-sill cap).
terrain_valley <- mh_terrain(
  valley_depth = 80,
  taf          = 1.4,
  relief       = 30
)


# ---------------------------------------------------------------------------
# --- Section 1: ventilation_state() over 48 h and 168 h -------------------
# ---------------------------------------------------------------------------

message("Benchmarking: ventilation_state() at 48 h / 168 h (flat terrain)")

bm_vs_flat <- bench::mark(
  `vent_state_48h`  = ventilation_state(met_48),
  `vent_state_168h` = ventilation_state(met_168),
  iterations = 20,
  check = FALSE
)
print(bm_vs_flat[, c("expression", "min", "median", "mem_alloc", "n_itr")])

message("Benchmarking: ventilation_state() at 48 h / 168 h (valley terrain)")

bm_vs_terrain <- bench::mark(
  `vent_state_terrain_48h`  = ventilation_state(met_48,  terrain = terrain_valley),
  `vent_state_terrain_168h` = ventilation_state(met_168, terrain = terrain_valley),
  iterations = 20,
  check = FALSE
)
print(bm_vs_terrain[, c("expression", "min", "median", "mem_alloc", "n_itr")])


# ---------------------------------------------------------------------------
# --- Section 2: pool_top heat-deficit at 168 h ----------------------------
# ---------------------------------------------------------------------------
# pool_top is computed inside ventilation_state(); isolate its contribution by
# comparing with- vs without-terrain at 168 h (the delta = the terrain path).

message("Benchmarking: pool_top overhead at 168 h (terrain vs flat delta)")

bm_pool_top <- bench::mark(
  `flat_168h`    = ventilation_state(met_168),
  `terrain_168h` = ventilation_state(met_168, terrain = terrain_valley),
  iterations = 20,
  check = FALSE
)
print(bm_pool_top[, c("expression", "min", "median", "mem_alloc", "n_itr")])
cat("pool_top overhead (median delta):",
    format(median(bm_pool_top$time[[2]]) - median(bm_pool_top$time[[1]])), "\n")


# ---------------------------------------------------------------------------
# --- Section 3: odour_hazard() at 168 h -----------------------------------
# ---------------------------------------------------------------------------

message("Benchmarking: odour_hazard() at 168 h")

bm_hazard <- bench::mark(
  `odour_hazard_168h` = odour_hazard(met_168),
  iterations = 20,
  check = FALSE
)
print(bm_hazard[, c("expression", "min", "median", "mem_alloc", "n_itr")])
