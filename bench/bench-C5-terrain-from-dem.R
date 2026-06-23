# bench-C5-terrain-from-dem.R
# ===========================================================================
# C5 benchmarks: mh_terrain_from_dem() setup timing + run-path guard.
#
# (1) Run-path guard: terra::extract on a typical-sized DEM must be sub-ms
#     (C5 success criterion). The run path consumes only the saved mh_terrain
#     scalar fields — no terra I/O — so individual extract calls are the
#     tightest relevant lower bound.
# (2) Reference setup timings for mh_terrain_from_dem() at two DEM sizes
#     (informational; these call WhiteboxTools and take seconds).
#
# Run with: source("bench/bench-C5-terrain-from-dem.R")
# ===========================================================================

library(meteoHazard)
library(terra)
library(bench)
library(sf)

# ---------------------------------------------------------------------------
# (1) Run-path guard: terra::extract must be sub-millisecond
# ---------------------------------------------------------------------------

# Represent a typical production DEM: 500 x 500 cells, 5 m resolution
# (2.5 km × 2.5 km footprint in UTM zone 55S).
dem_r <- terra::rast(
  nrows = 500L, ncols = 500L,
  xmin = 334000, xmax = 336500,
  ymin = 6249000, ymax = 6251500,
  crs = "EPSG:32755"
)
set.seed(42)
terra::values(dem_r) <- cumsum(rnorm(terra::ncell(dem_r), sd = 0.05))

pt <- terra::vect(matrix(c(335250, 6250250), ncol = 2), crs = "EPSG:32755")

cat("--- Run-path guard: terra::extract (single point) ---\n")
bm_extract <- bench::mark(
  `extract_500x500` = terra::extract(dem_r, pt)[[1, 2]],
  iterations = 100,
  check = FALSE
)
print(bm_extract[, c("expression", "min", "median", "mem_alloc", "n_itr")])

extract_ms <- as.numeric(median(bm_extract$time[[1]])) * 1000
if (extract_ms >= 1) {
  warning(sprintf(
    "terra::extract median %.3f ms >= 1 ms -- run-path guard FAILED", extract_ms
  ))
} else {
  cat(sprintf("terra::extract: %.3f ms -- sub-ms guard PASSED\n\n", extract_ms))
}

# Also exercise a larger DEM (1000 x 1000) to confirm linear scaling.
dem_lg <- terra::rast(
  nrows = 1000L, ncols = 1000L,
  xmin = 334000, xmax = 339000,
  ymin = 6246000, ymax = 6251000,
  crs = "EPSG:32755"
)
terra::values(dem_lg) <- cumsum(rnorm(terra::ncell(dem_lg), sd = 0.05))

bm_lg <- bench::mark(
  `extract_1000x1000` = terra::extract(dem_lg, pt)[[1, 2]],
  iterations = 100,
  check = FALSE
)
print(bm_lg[, c("expression", "min", "median", "mem_alloc", "n_itr")])
cat(sprintf("1000x1000 DEM extract: %.3f ms\n\n",
            as.numeric(median(bm_lg$time[[1]])) * 1000))

# ---------------------------------------------------------------------------
# (2) Setup timing: mh_terrain_from_dem() (WhiteboxTools required)
# ---------------------------------------------------------------------------

if (!isTRUE(whitebox::check_whitebox_binary())) {
  cat("WhiteboxTools binary not found -- skipping setup timing.\n")
} else {
  cat("--- Setup timing: mh_terrain_from_dem() (runs WhiteboxTools) ---\n")
  cat("Note: setup is designed to be called ONCE; these timings are informational.\n")

  cx <- 335000L; cy <- 6250000L
  src_sf <- sf::st_sf(id = "src",
    geometry = sf::st_sfc(sf::st_point(c(cx, cy)), crs = 32755L))

  dem_setup <- terra::rast(
    nrows = 128L, ncols = 128L,
    xmin = cx - 1000L, xmax = cx + 1000L,
    ymin = cy - 1000L, ymax = cy + 1000L,
    crs = "EPSG:32755"
  )
  xy <- terra::xyFromCell(dem_setup, seq_len(terra::ncell(dem_setup)))
  terra::values(dem_setup) <-
    80 * exp(-((xy[, "x"] - cx)^2 + (xy[, "y"] - cy)^2) / (2 * 300^2))

  t0 <- proc.time()["elapsed"]
  ter <- mh_terrain_from_dem(dem_setup, source = src_sf, epsg = 32755L)
  t1 <- proc.time()["elapsed"]
  cat(sprintf("mh_terrain_from_dem() on 128x128 DEM: %.1f s\n", t1 - t0))
  cat("Derived descriptors:\n")
  print(ter)
}
