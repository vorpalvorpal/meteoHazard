# Shared fixtures for the missing-input (NA) policy, shortwave derivation and
# dust gust-below-mean tests.

# Run `expr`, muffling and counting meteoHazard_missing_input warnings.
.count_missing <- function(expr) {
  n <- 0L
  value <- withCallingHandlers(
    expr,
    meteoHazard_missing_input = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, n = n)
}

.na_met <- function() {
  data.frame(
    datetime               = as.POSIXct("2024-01-15 00:00", tz = "UTC") + 3600 * (0:5),
    wind_speed_10m         = c(6, 9, 12, 15, 4, 8),
    wind_gusts_10m         = c(12, 16, 22, 28, 7, 14),
    wind_direction_10m     = c(270, 90, 45, 100, 200, 120),
    soil_moisture_0_to_1cm = c(0.02, 0.03, 0.02, 0.05, 0.1, 0.02),
    soil_moisture_1_to_3cm = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    precipitation          = c(0, 0, 0.1, 0, 0, 0),
    temperature_2m         = c(20, 22, 25, 30, 28, 26),
    relative_humidity_2m   = c(60, 55, 40, 30, 35, 45),
    surface_pressure       = c(1000, 1001, 1002, 1000, 999, 1000),
    pressure_msl           = c(1013, 1013, 1012, 1011, 1011, 1012),
    direct_radiation       = c(0, 0, 200, 600, 500, 100),
    diffuse_radiation      = c(0, 0, 80, 150, 120, 60),
    cloud_cover            = c(20, 40, 60, 10, 0, 90),
    boundary_layer_height  = c(300, 400, 800, 1500, 1200, 600)
  )
}

.na_litter_site <- function() {
  ctr <- sf::st_sf(
    id       = "ctr",
    geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
  )
  sectors <- data.frame(
    arc_start = c("NE", "SW"), arc_end = c("SE", "NW"),
    permeability = c(1.0, 0.3), sensitive = c(TRUE, FALSE)
  )
  suppressWarnings(site_from_sectors(sectors, ctr, epsg = 32755L))
}

.na_odour_site <- function() {
  feats <- sf::st_sf(
    id       = c("src", "rec1", "rec2"),
    geometry = sf::st_sfc(
      sf::st_point(c(335000, 6250000)),
      sf::st_point(c(335000, 6250400)),
      sf::st_point(c(335700, 6250000)),
      crs = 32755
    )
  )
  roles <- data.frame(
    feature_id = c("src", "rec1", "rec2"), hazard = "odour",
    role = c("source", "receptor", "receptor")
  )
  mh_site(feats, roles, epsg = 32755L)
}

