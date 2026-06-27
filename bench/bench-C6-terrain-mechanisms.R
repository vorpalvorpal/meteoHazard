# bench-C6-terrain-mechanisms.R
# ===========================================================================
# Phase 3 (C6) benchmarks: M3 valley sheltering on/off.
# Run with: source("bench/bench-C6-terrain-mechanisms.R")
# Requires the bench package: install.packages("bench")
# ===========================================================================

library(meteoHazard)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

n_168 <- 168
n_r10 <- 10
n_r100 <- 100

# Standard 168-hour met (one week).
make_met <- function(n = 168) {
  data.frame(
    wind_direction_10m     = rep(180, n),
    wind_speed_10m         = c(rep(2, n / 2), rep(1, n / 2)),
    direct_radiation       = c(rep(0, n / 2), rep(200, n / 2)),
    cloud_cover            = rep(20, n),
    boundary_layer_height  = c(rep(150, n / 2), rep(600, n / 2)),
    temperature_2m         = rep(10, n),
    relative_humidity_2m   = rep(70, n),
    precipitation          = rep(0, n),
    pressure_msl           = rep(1013, n),
    soil_moisture_0_to_1cm = rep(0.1, n),
    soil_moisture_1_to_3cm = rep(0.1, n)
  )
}

make_site <- function(n_rec = 10) {
  origin_x <- 335000; origin_y <- 6250000
  angles   <- seq(0, 2 * pi * (1 - 1 / n_rec), length.out = n_rec)
  rec_pts  <- lapply(seq_len(n_rec), function(i)
    sf::st_point(c(origin_x + 500 * sin(angles[i]),
                   origin_y + 500 * cos(angles[i]))))
  feats <- sf::st_sf(
    id          = c("src", paste0("r", seq_len(n_rec))),
    emit_height = c(5, rep(NA_real_, n_rec)),
    geometry    = sf::st_sfc(
      c(list(sf::st_point(c(origin_x, origin_y))), rec_pts),
      crs = 32755
    )
  )
  roles <- data.frame(
    feature_id = c("src", paste0("r", seq_len(n_rec))),
    hazard     = "odour",
    role       = c("source", rep("receptor", n_rec)),
    stringsAsFactors = FALSE
  )
  mh_site(feats, roles, epsg = 32755L)
}

ter_valley <- mh_terrain(
  shelter_index    = 50,
  flow_convergence = 0.8,
  drainage_bearing = 30,
  valley_depth     = 60,
  taf              = 1.5
)

d_168   <- make_met(168)
site_10 <- make_site(10)
site_100 <- make_site(100)

# ---------------------------------------------------------------------------
# M3 shelter benchmarks
# ---------------------------------------------------------------------------

cat("--- M3 ventilation_state() shelter on/off ---\n")
bench::mark(
  shelter_off_48 = ventilation_state(make_met(48),  terrain = ter_valley, shelter = FALSE),
  shelter_on_48  = ventilation_state(make_met(48),  terrain = ter_valley, shelter = TRUE),
  shelter_off_168 = ventilation_state(make_met(168), terrain = ter_valley, shelter = FALSE),
  shelter_on_168  = ventilation_state(make_met(168), terrain = ter_valley, shelter = TRUE),
  iterations = 20, check = FALSE
) |> print()
