# bench-C8-rim-venting.R
# ===========================================================================
# Phase 3 (C8) benchmarks: rim_venting on/off.
# Run with: source("bench/bench-C8-rim-venting.R")
# Requires the bench package: install.packages("bench")
# ===========================================================================

library(meteoHazard)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

n_168 <- 168

make_met <- function(n = 168) {
  n_night <- n - 3L
  blh     <- c(rep(150, n_night), 400, 800, 1500)
  data.frame(
    wind_direction_10m     = c(rep(0, n_night), rep(180, 3L)),
    wind_speed_10m         = rep(1.5, n),
    direct_radiation       = c(rep(0, n_night), rep(400, 3L)),
    cloud_cover            = c(rep(0, n_night), rep(20,  3L)),
    boundary_layer_height  = blh,
    temperature_2m         = rep(5,    n),
    pressure_msl           = rep(1013, n),
    precipitation          = rep(0,    n),
    relative_humidity_2m   = rep(60,   n),
    soil_moisture_0_to_1cm = rep(0.1,  n),
    soil_moisture_1_to_3cm = rep(0.1,  n)
  )
}

make_rim_site <- function(n_rec = 10) {
  origin_x <- 335000; origin_y <- 6250000
  angles   <- seq(0, 2 * pi * (1 - 1 / n_rec), length.out = n_rec)
  distances <- rep(1000, n_rec)

  ids    <- c("src", paste0("r", seq_len(n_rec)))
  pts    <- c(
    list(sf::st_point(c(origin_x, origin_y))),
    lapply(seq_len(n_rec), function(i)
      sf::st_point(c(origin_x + distances[i] * sin(angles[i]),
                     origin_y + distances[i] * cos(angles[i]))))
  )
  feats <- sf::st_sf(
    id            = ids,
    rel_elevation = c(NA_real_, rep(80, n_rec)),
    aspect        = c(NA_real_, (angles * 180 / pi + 180) %% 360),
    geometry      = sf::st_sfc(pts, crs = 32755)
  )
  roles <- data.frame(
    feature_id = ids,
    hazard     = "odour",
    role       = c("source", rep("receptor", n_rec)),
    stringsAsFactors = FALSE
  )
  mh_site(feats, roles,
          terrain = mh_terrain(drainage_bearing = 0, flow_convergence = 0.8,
                               valley_depth = 80, taf = 1.5),
          epsg = 32755L)
}

d_168    <- make_met(168)
site_10  <- make_rim_site(10)
site_100 <- make_rim_site(100)

# ---------------------------------------------------------------------------
# C8 rim_venting on/off benchmarks
# ---------------------------------------------------------------------------

cat("--- C8 rim_venting on/off: descriptors backend, 168 h ---\n")
bench::mark(
  off_10  = odour_exposure(d_168, site_10,  terrain_backend = "descriptors",
                           rim_venting = FALSE),
  on_10   = odour_exposure(d_168, site_10,  terrain_backend = "descriptors",
                           rim_venting = TRUE),
  off_100 = odour_exposure(d_168, site_100, terrain_backend = "descriptors",
                           rim_venting = FALSE),
  on_100  = odour_exposure(d_168, site_100, terrain_backend = "descriptors",
                           rim_venting = TRUE),
  iterations = 10, check = FALSE
) |> print()
