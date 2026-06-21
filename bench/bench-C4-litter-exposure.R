# ===========================================================================
# C4 benchmarks: litter_exposure() + site_from_sectors()
# Issue #18, Section 7
#
# Benchmarks:
#   - litter_exposure() over 168 h x {4, 16} barriers
#   - site_from_sectors() construction (setup-time)
#
# Run interactively; results printed to console.
# ===========================================================================

library(meteoHazard)
library(bench)
library(sf)

# ---------------------------------------------------------------------------
# Hazard + wind fixture
# ---------------------------------------------------------------------------

set.seed(42)
n_hours           <- 168
hazard_168        <- runif(n_hours, 0, 100)
wind_dir_168      <- runif(n_hours, 0, 360)


# ---------------------------------------------------------------------------
# site_from_sectors() fixture builder
# ---------------------------------------------------------------------------

# site_from_sectors() requires centroid as an sf object.
centroid_sf <- sf::st_sf(
  id       = "centroid",
  geometry = sf::st_sfc(sf::st_point(c(339000, 6251000)), crs = 32755)
)

make_sector_site <- function(n_barriers = 4) {
  # Build n_barriers evenly-spaced compass sectors around the full circle.
  # Compass labels cycle through the 8 cardinal/intercardinal points.
  compass_labels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  n_pts <- length(compass_labels)

  sectors <- data.frame(
    arc_start    = compass_labels[((seq_len(n_barriers) - 1) %% n_pts) + 1],
    arc_end      = compass_labels[(seq_len(n_barriers)       %% n_pts) + 1],
    permeability = rep_len(c(1.0, 0.3, 0.7, 0.5), n_barriers),
    sensitive    = rep_len(c(TRUE, FALSE), n_barriers),
    stringsAsFactors = FALSE
  )
  site_from_sectors(sectors, centroid = centroid_sf, epsg = 32755L)
}

site_4barrier  <- make_sector_site(4)
site_16barrier <- make_sector_site(16)


# ---------------------------------------------------------------------------
# --- Section 1: site_from_sectors() construction (setup-time) -------------
# ---------------------------------------------------------------------------

message("Benchmarking: site_from_sectors() construction at 4 / 16 barriers")

bm_construction <- bench::mark(
  `sectors_4barrier`  = make_sector_site(4),
  `sectors_16barrier` = make_sector_site(16),
  iterations = 20,
  check = FALSE
)
print(bm_construction[, c("expression", "min", "median", "mem_alloc", "n_itr")])


# ---------------------------------------------------------------------------
# --- Section 2: litter_exposure() over 168 h x {4, 16} barriers ----------
# ---------------------------------------------------------------------------

message("Benchmarking: litter_exposure() 168 h x 4 barriers")

bm_litter_4 <- bench::mark(
  `litter_4barrier_168h` = litter_exposure(
    hazard             = hazard_168,
    wind_direction_10m = wind_dir_168,
    site               = site_4barrier
  ),
  iterations = 20,
  check = FALSE
)
print(bm_litter_4[, c("expression", "min", "median", "mem_alloc", "n_itr")])

message("Benchmarking: litter_exposure() 168 h x 16 barriers")

bm_litter_16 <- bench::mark(
  `litter_16barrier_168h` = litter_exposure(
    hazard             = hazard_168,
    wind_direction_10m = wind_dir_168,
    site               = site_16barrier
  ),
  iterations = 20,
  check = FALSE
)
print(bm_litter_16[, c("expression", "min", "median", "mem_alloc", "n_itr")])

# Combined for a side-by-side comparison.
message("Side-by-side 4 vs 16 barriers:")
bm_litter_combined <- bench::mark(
  `litter_4barrier_168h`  = litter_exposure(hazard_168, wind_dir_168, site_4barrier),
  `litter_16barrier_168h` = litter_exposure(hazard_168, wind_dir_168, site_16barrier),
  iterations = 20,
  check = FALSE
)
print(bm_litter_combined[, c("expression", "min", "median", "mem_alloc", "n_itr")])
