# ===========================================================================
# C1 benchmarks: mh_site construction + geometry helpers
# Issue #14, Section 7
#
# Benchmarks:
#   - .bearing_distance() + .role_features() over 1 source x {10,100,1000} recs
#   - mh_site() construction + validation + reprojection at same feature counts
#   - .crosswind_halfwidth() over 1000 wind directions for a fixed polygon
#
# Run interactively; results printed to console.
# ===========================================================================

library(meteoHazard)
library(bench)
library(sf)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

make_odour_site <- function(n_receptors = 1) {
  set.seed(42)
  src <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 32755)
  recs <- lapply(seq_len(n_receptors), function(k) {
    sf::st_point(c(500 * cos(2 * pi * k / n_receptors),
                   500 * sin(2 * pi * k / n_receptors)))
  })
  rec_sfc <- sf::st_sfc(recs, crs = 32755)

  feats <- sf::st_sf(
    id       = c("src", paste0("rec", seq_len(n_receptors))),
    geometry = c(src, rec_sfc)
  )
  roles <- data.frame(
    feature_id  = c("src", paste0("rec", seq_len(n_receptors))),
    hazard      = "odour",
    role        = c("source", rep("receptor", n_receptors)),
    emit_height = c(5, rep(NA_real_, n_receptors))
  )
  mh_site(features = feats, roles = roles, epsg = 32755L)
}

# Pre-build one large site for the geometry helper benchmarks (construction
# cost is excluded from those timings).
site_1000 <- make_odour_site(1000)
site_100  <- make_odour_site(100)
site_10   <- make_odour_site(10)

# Extract the source and receptor feature subsets used by the geometry helpers.
src_feat   <- meteoHazard:::.role_features(site_1000, "odour", "source")
rec_feat_10   <- meteoHazard:::.role_features(site_10,   "odour", "receptor")
rec_feat_100  <- meteoHazard:::.role_features(site_100,  "odour", "receptor")
rec_feat_1000 <- meteoHazard:::.role_features(site_1000, "odour", "receptor")

# A square polygon (100 m side) for .crosswind_halfwidth() benchmarks.
set.seed(42)
square_poly <- sf::st_sfc(
  sf::st_polygon(list(cbind(
    c(-50, 50, 50, -50, -50),
    c(-50, -50, 50,  50, -50)
  ))),
  crs = 32755
)

wind_dirs_1000 <- seq(0, 359.64, length.out = 1000)


# ---------------------------------------------------------------------------
# --- Section 1: .bearing_distance() + .role_features() run-path cost -------
# ---------------------------------------------------------------------------

message("Benchmarking: geometry helpers over 10 / 100 / 1000 receptors")

bm_geometry <- bench::mark(
  `role_features_10`      = meteoHazard:::.role_features(site_10,   "odour", "receptor"),
  `role_features_100`     = meteoHazard:::.role_features(site_100,  "odour", "receptor"),
  `role_features_1000`    = meteoHazard:::.role_features(site_1000, "odour", "receptor"),
  `bearing_dist_10`       = meteoHazard:::.bearing_distance(src_feat, rec_feat_10),
  `bearing_dist_100`      = meteoHazard:::.bearing_distance(src_feat, rec_feat_100),
  `bearing_dist_1000`     = meteoHazard:::.bearing_distance(src_feat, rec_feat_1000),
  iterations = 20,
  check = FALSE
)
print(bm_geometry[, c("expression", "min", "median", "mem_alloc", "n_itr")])


# ---------------------------------------------------------------------------
# --- Section 2: mh_site() construction + validation + reprojection ---------
# ---------------------------------------------------------------------------

message("Benchmarking: mh_site() construction at 10 / 100 / 1000 features")

# Helper that builds the raw inputs without calling mh_site() (so we time
# mh_site() itself, not make_odour_site()).
make_site_inputs <- function(n_receptors) {
  src <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 32755)
  recs <- lapply(seq_len(n_receptors), function(k) {
    sf::st_point(c(500 * cos(2 * pi * k / n_receptors),
                   500 * sin(2 * pi * k / n_receptors)))
  })
  rec_sfc <- sf::st_sfc(recs, crs = 32755)
  feats <- sf::st_sf(
    id       = c("src", paste0("rec", seq_len(n_receptors))),
    geometry = c(src, rec_sfc)
  )
  roles <- data.frame(
    feature_id  = c("src", paste0("rec", seq_len(n_receptors))),
    hazard      = "odour",
    role        = c("source", rep("receptor", n_receptors)),
    emit_height = c(5, rep(NA_real_, n_receptors))
  )
  list(features = feats, roles = roles)
}

set.seed(42)
inp_10   <- make_site_inputs(10)
inp_100  <- make_site_inputs(100)
inp_1000 <- make_site_inputs(1000)

bm_construction <- bench::mark(
  `mh_site_10`   = mh_site(features = inp_10$features,   roles = inp_10$roles,   epsg = 32755L),
  `mh_site_100`  = mh_site(features = inp_100$features,  roles = inp_100$roles,  epsg = 32755L),
  `mh_site_1000` = mh_site(features = inp_1000$features, roles = inp_1000$roles, epsg = 32755L),
  iterations = 20,
  check = FALSE
)
print(bm_construction[, c("expression", "min", "median", "mem_alloc", "n_itr")])


# ---------------------------------------------------------------------------
# --- Section 3: .crosswind_halfwidth() over 1000 wind directions -----------
# ---------------------------------------------------------------------------

message("Benchmarking: .crosswind_halfwidth() over 1000 wind directions")

bm_crosswind <- bench::mark(
  `crosswind_1000dirs` = vapply(
    wind_dirs_1000,
    function(wd) meteoHazard:::.crosswind_halfwidth(square_poly, wd),
    numeric(1)
  ),
  iterations = 20,
  check = FALSE
)
print(bm_crosswind[, c("expression", "min", "median", "mem_alloc", "n_itr")])
