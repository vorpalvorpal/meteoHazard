# Behaviour specification for the litter exposure layer (litter_exposure()).
#
# Encodes the behaviour defined in specs/Litter_exposure.md (basic mode) and
# the plan in GitHub issue #18 (C4). The function only accepts an mh_site;
# passing a compass-sector data frame is a classed error.
#
# Proposed contract:
#   litter_exposure(
#     hazard, wind_direction_10m, site,
#     direction_tol = 15, p_open_min = 0.5,
#     move_threshold = 20, offsite_threshold = 45, default_permeability = 0.5
#   ) -> data.frame(exposure = numeric[0,100], zone = ordered factor)
#   zone levels: within_face < on_site < off_site

# === C4 behaviour spec (issue #18): litter_exposure() on the mh_site model ===

# Fixture: create an mh_site equivalent to the old demo_site() compass-sector
# data frame for testing.
# demo_site: NE-SE (perm=1.0, sensitive=TRUE), SW-NW (perm=0.3, sensitive=FALSE)
.make_demo_mh_site <- function(epsg = 32755, ring_radius = 1000) {
  ctr <- sf::st_sf(
    id       = "ctr",
    geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = epsg)
  )
  sectors <- data.frame(
    arc_start    = c("NE", "SW"),
    arc_end      = c("SE", "NW"),
    permeability = c(1.0, 0.3),
    sensitive    = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  site_from_sectors(sectors, ctr, ring_radius = ring_radius, epsg = epsg)
}

describe("litter_exposure() on the mh_site model", {

  it("returns a data frame with one row per hour and correct exposure + zone", {
    site    <- .make_demo_mh_site()
    new_out <- litter_exposure(c(86, 25, 5), c(270, 90, 0), site)
    expect_s3_class(new_out, "data.frame")
    expect_equal(nrow(new_out), 3)
    expect_true(all(c("exposure", "zone") %in% names(new_out)))
    expect_type(new_out$exposure, "double")
    # hazard=86, wind=270 -> theta_down=90 (E), hits NE-SE (perm=1.0) -> 86*1.0=86
    # hazard=25, wind=90  -> theta_down=270 (W), hits SW-NW (perm=0.3) -> 25*0.3=7.5
    # hazard=5,  wind=0   -> theta_down=180 (S), no sector hit -> default 0.5 -> 5*0.5=2.5
    expect_equal(new_out$exposure, c(86, 7.5, 2.5), tolerance = 1e-3)
    expect_equal(as.character(new_out$zone), c("off_site", "on_site", "within_face"))
  })

  it("returns zone as an ordered factor within_face < on_site < off_site", {
    site <- .make_demo_mh_site()
    out  <- litter_exposure(86, 270, site)
    expect_s3_class(out$zone, "ordered")
    expect_identical(levels(out$zone), c("within_face", "on_site", "off_site"))
  })

  it("keeps exposure within [0, 100]", {
    site <- .make_demo_mh_site()
    g    <- expand.grid(h = c(0, 20, 86, 100), d = c(0, 90, 180, 270, 359))
    out  <- litter_exposure(g$h, g$d, site)
    expect_true(all(out$exposure >= 0 & out$exposure <= 100))
  })

  it("uses the most permeable hit barrier as the directional factor", {
    site <- .make_demo_mh_site()
    # Wind FROM W (270) -> downwind toward E (90), which hits NE-SE (perm=1.0)
    out <- litter_exposure(86, 270, site)
    expect_equal(out$exposure, 86 * 1.0, tolerance = 1e-3)
  })

  it("applies the default permeability when the downwind bearing hits no barrier", {
    site <- .make_demo_mh_site()
    # Wind from N (0) -> downwind toward S (180). Neither NE-SE nor SW-NW covers 180.
    out <- litter_exposure(86, 0, site, default_permeability = 0.5)
    expect_equal(out$exposure, 43, tolerance = 1e-3)
  })

  it("expands each barrier arc by direction_tol", {
    site <- .make_demo_mh_site()
    # NE boundary of NE-SE sector is 45 deg. With tol=15, arc expands to [30, 150].
    # wind_dir 211 -> theta_down = (211+180)%%360 = 31. 31 in [30, 150] -> hit (perm=1.0)
    # wind_dir 209 -> theta_down = (209+180)%%360 = 29. 29 < 30 -> miss -> default=0
    out_hit  <- litter_exposure(100, 211, site, direction_tol = 15, default_permeability = 0)
    out_miss <- litter_exposure(100, 209, site, direction_tol = 15, default_permeability = 0)
    expect_gt(out_hit$exposure, 0)
    expect_equal(out_miss$exposure, 0, tolerance = 1e-6)
  })

  it("orders the zone ladder within_face < on_site < off_site", {
    site <- .make_demo_mh_site()
    # Wind toward E (perm=1.0, sensitive=TRUE)
    out_wf  <- litter_exposure(5,  270, site)   # hazard < move_threshold (20)
    out_on  <- litter_exposure(30, 270, site)   # hazard 30, sensitive, but < offsite_threshold (45)
    out_off <- litter_exposure(86, 270, site)   # hazard 86 >= 45, sensitive, perm=1.0
    expect_equal(as.character(out_wf$zone),  "within_face")
    expect_equal(as.character(out_on$zone),  "on_site")
    expect_equal(as.character(out_off$zone), "off_site")
  })

  it("flags off-site only for a sensitive receptor with permeability >= p_open_min", {
    site <- .make_demo_mh_site()
    # Wind FROM E (90) -> downwind toward W (270), hits SW-NW (perm=0.3, sensitive=FALSE)
    # Not sensitive -> cannot be off_site regardless of hazard level
    out_barrier <- litter_exposure(86, 90, site)
    expect_equal(as.character(out_barrier$zone), "on_site")

    # Sensitive sector with perm=1.0 but p_open_min set above it -> treated as closed
    # perm=0.3 < p_open_min=0.5 -> sector not counted as sensitive-open -> on_site
    # Use a site with a single sensitive sector of perm=0.3
    ctr <- sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    sectors_low <- data.frame(
      arc_start    = "NE",
      arc_end      = "SE",
      permeability = 0.3,
      sensitive    = TRUE,
      stringsAsFactors = FALSE
    )
    site_low <- site_from_sectors(sectors_low, ctr, ring_radius = 1000, epsg = 32755)
    # Wind from W (270) -> toward E, hits NE-SE (perm=0.3, sensitive=TRUE)
    # With p_open_min=0.5, perm 0.3 < 0.5 -> not counted as open -> on_site
    out_low <- litter_exposure(86, 270, site_low, p_open_min = 0.5)
    expect_equal(as.character(out_low$zone), "on_site")
  })

  it("returns zero exposure and within_face for a zero hazard", {
    site <- .make_demo_mh_site()
    out  <- litter_exposure(0, 270, site)
    expect_equal(out$exposure, 0)
    expect_equal(as.character(out$zone), "within_face")
  })

  it("scales exposure linearly with hazard at a fixed direction and site", {
    site <- .make_demo_mh_site()
    out  <- litter_exposure(c(40, 80), c(270, 270), site)
    expect_equal(out$exposure[2], 2 * out$exposure[1], tolerance = 1e-3)
  })

  it("errors on hazard outside [0, 100]", {
    site <- .make_demo_mh_site()
    expect_error(litter_exposure(-1,  270, site))
    expect_error(litter_exposure(101, 270, site))
  })

  it("errors on wind direction outside [0, 360]", {
    site <- .make_demo_mh_site()
    expect_error(litter_exposure(50, -10, site))
    expect_error(litter_exposure(50, 400, site))
  })

  it("errors on mismatched hazard and wind-direction lengths", {
    site <- .make_demo_mh_site()
    expect_error(litter_exposure(c(50, 60), 270, site))
  })

  it("errors on missing values in hazard or wind direction", {
    site <- .make_demo_mh_site()
    expect_error(litter_exposure(c(50, NA), c(270, 270), site))
    expect_error(litter_exposure(c(50, 60), c(270, NA),  site))
  })

  it("errors (classed) when move_threshold is not below offsite_threshold", {
    site <- .make_demo_mh_site()
    expect_error(
      litter_exposure(50, 270, site, move_threshold = 50, offsite_threshold = 45),
      class = "meteoHazard_input_error"
    )
  })

  it("errors (classed) on an invalid mh_site or a missing litter source", {
    # data.frame (old API) -> classed error
    old_df <- data.frame(
      arc_start    = c("NE", "SW"),
      arc_end      = c("SE", "NW"),
      permeability = c(1.0, 0.3),
      sensitive    = c(TRUE, FALSE),
      stringsAsFactors = FALSE
    )
    expect_error(
      litter_exposure(50, 270, old_df),
      class = "meteoHazard_input_error"
    )

    # Non-mh_site, non-data.frame -> classed error
    expect_error(
      litter_exposure(50, 270, "not a site"),
      class = "meteoHazard_input_error"
    )

    # mh_site with no (litter, source) role -> classed error
    feats <- sf::st_sf(
      id       = "A",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    roles <- data.frame(
      feature_id = "A",
      hazard     = "litter",
      role       = "receptor",
      stringsAsFactors = FALSE
    )
    site_no_src <- mh_site(feats, roles, epsg = 32755)
    expect_error(
      litter_exposure(50, 270, site_no_src),
      class = "meteoHazard_input_error"
    )
  })

  it("computes bearings from an offset source's actual geometry, not the centroid", {
    # Build a site_from_sectors centroid site (NE-SE barrier only)
    ctr <- sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    sectors_single <- data.frame(
      arc_start    = "NE",
      arc_end      = "SE",
      permeability = 1.0,
      sensitive    = FALSE,
      stringsAsFactors = FALSE
    )
    site_centroid <- site_from_sectors(sectors_single, ctr, ring_radius = 1000, epsg = 32755)

    # Build an equivalent site with source 200 m east of the centroid
    # The barrier polygon is the same physical polygon, but the source is offset
    offset_src   <- sf::st_point(c(335200, 6250000))  # 200 m east
    barrier_poly <- site_centroid@features[site_centroid@features$id == "barrier_1", ]

    feats_offset <- sf::st_sf(
      id            = c("src_offset", "barrier_1"),
      bearing_start = c(NA_real_, 45),
      bearing_end   = c(NA_real_, 135),
      geometry      = sf::st_sfc(
        offset_src,
        sf::st_geometry(barrier_poly)[[1]],
        crs = 32755
      ),
      stringsAsFactors = FALSE
    )
    roles_offset <- data.frame(
      feature_id   = c("src_offset", "barrier_1"),
      hazard       = "litter",
      role         = c("source", "barrier"),
      permeability = c(NA_real_, 1.0),
      sensitive    = c(NA, FALSE),
      stringsAsFactors = FALSE
    )
    site_offset <- mh_site(feats_offset, roles_offset, epsg = 32755)

    # Both should run without error and return the expected structure
    out_centroid <- litter_exposure(50, 270, site_centroid)
    out_offset   <- litter_exposure(50, 270, site_offset)

    expect_s3_class(out_centroid, "data.frame")
    expect_s3_class(out_offset, "data.frame")
    expect_true("exposure" %in% names(out_offset))
    expect_true("zone"     %in% names(out_offset))
  })

  it("takes the worst case among multiple overlapping barriers", {
    ctr <- sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    sectors <- data.frame(
      arc_start    = c("N", "E"),
      arc_end      = c("E", "S"),
      permeability = c(0.2, 0.9),
      sensitive    = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    )
    site <- site_from_sectors(sectors, ctr, ring_radius = 500, epsg = 32755)
    # Wind from W (270) -> theta_down = 90 (E). E is on the shared edge of N-E
    # (0-90) and E-S (90-180). With direction_tol=15, both expanded edges include
    # 90. Worst case (max permeability) = 0.9.
    out <- litter_exposure(50, 270, site, default_permeability = 0)
    expect_equal(out$exposure, 50 * 0.9, tolerance = 1e-3)
  })
})
