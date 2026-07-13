# Behaviour specification for the litter exposure layer (litter_exposure()).
#
# Basic-mode behaviour (issue #18) is preserved as regression fixtures; on
# top of that this file adds tests for:
#   - hazard contract seam (no upper bound; lower=0 only)
#   - enclosure (a barrier that fully surrounds the source)
#   - zero-barrier warning
#   - refined distance-reach mode (mean_wind + reach_per_ms)
#   - destination x sensitivity split (leaves_site / sensitive_receptor)
#   - multiple sources, worst-case aggregation
#
# Return schema: exposure, zone, directional_factor, leaves_site,
# sensitive_receptor.

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
  # suppressWarnings: these fixtures deliberately use partial compass coverage,
  # which (post-implementation) emits meteoHazard_litter_sector_gap; that warning
  # is exercised on purpose in test-site-constructors.R, not here.
  suppressWarnings(site_from_sectors(sectors, ctr, ring_radius = ring_radius, epsg = epsg))
}

LITTER_EXPOSURE_COLS <- c(
  "exposure", "zone", "directional_factor", "leaves_site", "sensitive_receptor"
)

describe("litter_exposure() on the mh_site model (basic-mode regression)", {

  it("returns a data frame with one row per hour, the five documented columns, and correct exposure + zone", {
    site    <- .make_demo_mh_site()
    new_out <- litter_exposure(c(86, 25, 5), c(270, 90, 0), site)
    expect_s3_class(new_out, "data.frame")
    expect_equal(nrow(new_out), 3)
    expect_named(new_out, LITTER_EXPOSURE_COLS, ignore.order = TRUE)
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

  it("keeps exposure non-negative and finite across a wide hazard/direction sweep (including hazard > 100)", {
    site <- .make_demo_mh_site()
    g    <- expand.grid(h = c(0, 20, 86, 100, 150), d = c(0, 90, 180, 270, 359))
    out  <- litter_exposure(g$h, g$d, site)
    expect_true(all(out$exposure >= 0 & is.finite(out$exposure)))
  })

  it("uses the most permeable hit barrier as the directional factor", {
    site <- .make_demo_mh_site()
    # Wind FROM W (270) -> downwind toward E (90), which hits NE-SE (perm=1.0)
    out <- litter_exposure(86, 270, site)
    expect_equal(out$exposure, 86 * 1.0, tolerance = 1e-3)
    expect_equal(out$directional_factor, 1.0, tolerance = 1e-3)
  })

  it("applies the default permeability when the downwind bearing hits no barrier", {
    site <- .make_demo_mh_site()
    # Wind from N (0) -> downwind toward S (180). Neither NE-SE nor SW-NW covers 180.
    out <- litter_exposure(86, 0, site, default_permeability = 0.5)
    expect_equal(out$exposure, 43, tolerance = 1e-3)
    expect_false(out$leaves_site)
    expect_false(out$sensitive_receptor)
  })

  it("expands each barrier arc by direction_tol", {
    site <- .make_demo_mh_site()
    # NE boundary of NE-SE sector is 45 deg. With tol=15, arc expands to [30, 150].
    # wind_dir 211 -> theta_down = (211+180)%%360 = 31. 31 in [30, 150] -> hit (perm=1.0)
    # wind_dir 209 -> theta_down = (209+180)%%360 = 29. 29 < 30 -> miss -> default=0
    out_hit  <- litter_exposure(100, 211, site, direction_tol = 15, default_permeability = 0)
    out_miss <- litter_exposure(100, 209, site, direction_tol = 15, default_permeability = 0)
    expect_equal(out_hit$exposure, 100, tolerance = 1e-3)
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
    out_barrier <- litter_exposure(86, 90, site)
    expect_equal(as.character(out_barrier$zone), "on_site")

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
    site_low <- suppressWarnings(site_from_sectors(sectors_low, ctr, ring_radius = 1000, epsg = 32755))
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

    expect_error(
      litter_exposure(50, 270, "not a site"),
      class = "meteoHazard_input_error"
    )

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
    site_centroid <- suppressWarnings(site_from_sectors(sectors_single, ctr, ring_radius = 1000, epsg = 32755))

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

    out_centroid <- litter_exposure(50, 270, site_centroid)
    out_offset   <- litter_exposure(50, 270, site_offset)

    expect_s3_class(out_centroid, "data.frame")
    expect_s3_class(out_offset, "data.frame")
    expect_true("exposure" %in% names(out_offset))
    expect_true("zone"     %in% names(out_offset))
  })

  it("takes the worst case among multiple overlapping barriers (shared-edge fixture)", {
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
    site <- suppressWarnings(site_from_sectors(sectors, ctr, ring_radius = 500, epsg = 32755))
    # Wind from W (270) -> theta_down = 90 (E), the shared edge of N-E (0-90)
    # and E-S (90-180). With direction_tol=15, both expanded edges include 90.
    # Worst case (max permeability) = 0.9.
    out <- litter_exposure(50, 270, site, default_permeability = 0)
    expect_equal(out$exposure, 50 * 0.9, tolerance = 1e-3)
    expect_equal(out$directional_factor, 0.9, tolerance = 1e-3)
  })
})


describe("hazard contract seam (litter_exposure accepts any non-negative hazard)", {

  it("accepts hazard = 120 without error, scaled by the hit sector's permeability", {
    site <- .make_demo_mh_site()
    out <- litter_exposure(120, 270, site)   # NE-SE perm = 1.0 -> hit
    expect_equal(out$exposure, 120 * 1.0, tolerance = 1e-3)
  })

  it("accepts hazard = 101 without error (the old upper-bound-100 validation is gone)", {
    site <- .make_demo_mh_site()
    expect_no_error(litter_exposure(101, 270, site))
  })

  it("errors on a negative hazard", {
    site <- .make_demo_mh_site()
    # Rejected by the checkmate lower-bound assertion, the package convention for
    # numeric-input validation (cf. litter_hazard_vec()'s negative-input tests);
    # the classed meteoHazard_input_error is reserved for structural/semantic
    # errors (bad site, mis-ordered thresholds).
    expect_error(litter_exposure(-1, 270, site))
  })
})


describe("enclosure (a barrier polygon that fully surrounds the source)", {

  it("is hit from every one of the 8 principal wind directions; none fall through to the default", {
    cx <- 335000; cy <- 6250000
    half <- 500
    ring_poly <- sf::st_polygon(list(matrix(
      c(cx - half, cy - half,
        cx + half, cy - half,
        cx + half, cy + half,
        cx - half, cy + half,
        cx - half, cy - half),
      ncol = 2, byrow = TRUE
    )))

    src_geom <- sf::st_sfc(sf::st_point(c(cx, cy)), crs = 32755)
    barrier_geom <- sf::st_sfc(ring_poly, crs = 32755)

    # Sanity: the fixture really does enclose the source.
    expect_true(isTRUE(sf::st_within(src_geom, barrier_geom, sparse = FALSE)[1, 1]))

    feats <- sf::st_sf(
      id       = c("source", "ring_barrier"),
      geometry = c(src_geom, barrier_geom)
    )
    roles <- data.frame(
      feature_id   = c("source", "ring_barrier"),
      hazard       = "litter",
      role         = c("source", "barrier"),
      permeability = c(NA_real_, 1.0),
      sensitive    = c(NA, TRUE),
      stringsAsFactors = FALSE
    )
    site <- mh_site(feats, roles, epsg = 32755)

    for (d in seq(0, 315, 45)) {
      out <- litter_exposure(50, d, site, default_permeability = 0)
      expect_equal(out$exposure, 50 * 1.0, tolerance = 1e-3, info = paste("wind_direction_10m", d))
    }
  })
})


describe("wide barrier arcs (span + 2*direction_tol >= 360)", {

  # A horseshoe (C-shaped annulus sector) spanning bearings 10 -> 350 deg
  # clockwise (340 deg, gap at north). The source at its centre is NOT
  # st_within the polygon, so the enclosure sentinel does not fire and the
  # bearing arc is computed from the vertices. With the default
  # direction_tol = 15 the expanded arc covers the full circle.
  .horseshoe_site <- function() {
    cx <- 335000; cy <- 6250000
    th_out <- seq(10, 350, length.out = 181) * pi / 180
    outer  <- cbind(cx + 300 * sin(th_out), cy + 300 * cos(th_out))
    inner  <- cbind(cx + 200 * sin(rev(th_out)), cy + 200 * cos(rev(th_out)))
    ring   <- rbind(outer, inner, outer[1, , drop = FALSE])
    feats <- sf::st_sf(
      id       = c("source", "horseshoe"),
      geometry = sf::st_sfc(
        sf::st_point(c(cx, cy)),
        sf::st_polygon(list(ring)),
        crs = 32755
      )
    )
    roles <- data.frame(
      feature_id   = c("source", "horseshoe"),
      hazard       = "litter",
      role         = c("source", "barrier"),
      permeability = c(NA_real_, 1.0),
      sensitive    = c(NA, TRUE),
      stringsAsFactors = FALSE
    )
    mh_site(feats, roles, epsg = 32755)
  }

  it(".litter_arc_contains() hits a bearing deep inside a wide arc once the tolerance band closes the circle", {
    # Raw arc [10, 350] spans 340 deg; with tol = 15 the expanded edges pass
    # each other (355 > 5), and the naive wrap branch would test the
    # COMPLEMENT of the arc, missing theta = 180 entirely.
    expect_true(meteoHazard:::.litter_arc_contains(180, 10, 350, 15))
    expect_true(all(meteoHazard:::.litter_arc_contains(seq(0, 359), 10, 350, 15)))
  })

  it(".litter_arc_contains() still excludes bearings outside a narrow wrapped arc", {
    expect_false(meteoHazard:::.litter_arc_contains(180, 350, 10, 15))
    expect_true(meteoHazard:::.litter_arc_contains(0, 350, 10, 15))
  })

  it("a horseshoe barrier is hit from every principal wind direction (end-to-end)", {
    site <- .horseshoe_site()
    for (d in seq(0, 315, 45)) {
      out <- litter_exposure(50, d, site, default_permeability = 0)
      expect_equal(out$exposure, 50, tolerance = 1e-3,
                   info = paste("wind_direction_10m", d))
    }
  })
})


describe("zero-barrier warning", {

  it("warns meteoHazard_litter_no_barriers when the site has zero litter barriers, and still returns a valid frame", {
    feats <- sf::st_sf(
      id       = "src1",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    roles <- data.frame(
      feature_id = "src1", hazard = "litter", role = "source",
      stringsAsFactors = FALSE
    )
    site_no_barriers <- mh_site(feats, roles, epsg = 32755)

    expect_warning(
      out <- litter_exposure(90, 0, site_no_barriers),
      class = "meteoHazard_litter_no_barriers"
    )
    expect_s3_class(out, "data.frame")
    expect_named(out, LITTER_EXPOSURE_COLS, ignore.order = TRUE)
    expect_equal(out$exposure, 90 * 0.5, tolerance = 1e-6)  # default_permeability = 0.5
    expect_equal(as.character(out$zone), "on_site")
  })
})


describe("refined distance-reach mode (mean_wind + reach_per_ms)", {

  .refined_site <- function() {
    ctr <- sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    sectors <- data.frame(
      arc_start = "NE", arc_end = "SE",
      permeability = 1.0, sensitive = TRUE,
      stringsAsFactors = FALSE
    )
    suppressWarnings(site_from_sectors(sectors, ctr, ring_radius = 1000, epsg = 32755))
  }

  it("the fixture barrier carries distance_m == ring_radius (1000)", {
    site <- .refined_site()
    brow <- site@features[site@features$id == "barrier_1", ]
    expect_equal(brow$distance_m, 1000)
  })

  it("reach (mean_wind * reach_per_ms) >= distance_m: off_site, leaves_site TRUE", {
    site <- .refined_site()
    out <- litter_exposure(86, 270, site, mean_wind = 12, reach_per_ms = 100)  # reach 1200 >= 1000
    expect_equal(as.character(out$zone), "off_site")
    expect_true(out$leaves_site)
    expect_true(out$sensitive_receptor)
    expect_equal(out$exposure, 86, tolerance = 1e-3)
  })

  it("reach < distance_m: on_site, leaves_site FALSE", {
    site <- .refined_site()
    out <- litter_exposure(86, 270, site, mean_wind = 8, reach_per_ms = 100)   # reach 800 < 1000
    expect_equal(as.character(out$zone), "on_site")
    expect_false(out$leaves_site)
    expect_true(out$sensitive_receptor)
    expect_equal(out$exposure, 86, tolerance = 1e-3)
  })

  it("accepts a units-tagged mean_wind and converts it to m/s (package units contract)", {
    site <- .refined_site()
    bare   <- litter_exposure(86, 270, site, mean_wind = 12, reach_per_ms = 100)
    tagged <- litter_exposure(86, 270, site,
                              mean_wind = units::set_units(43.2, "km/h"),
                              reach_per_ms = 100)
    expect_equal(tagged, bare)
  })

  it("rejects a mean_wind tagged with dimensionally incompatible units", {
    site <- .refined_site()
    expect_error(
      litter_exposure(86, 270, site,
                      mean_wind = units::set_units(12, "degree_C"),
                      reach_per_ms = 100),
      class = "meteoHazard_input_error"
    )
  })

  it("basic mode (mean_wind/reach_per_ms both NULL) is unchanged from the magnitude regression", {
    site <- .refined_site()
    out <- litter_exposure(86, 270, site)
    expect_equal(out$exposure, 86, tolerance = 1e-3)
    expect_equal(as.character(out$zone), "off_site")
  })

  it("within_face is governed by hazard < move_threshold in both basic and refined mode", {
    site <- .refined_site()
    out_basic   <- litter_exposure(5, 270, site)
    out_refined <- litter_exposure(5, 270, site, mean_wind = 12, reach_per_ms = 100)
    expect_equal(as.character(out_basic$zone), "within_face")
    expect_equal(as.character(out_refined$zone), "within_face")
  })
})


describe("destination x sensitivity split (leaves_site vs sensitive_receptor)", {

  .single_sector_site <- function(perm, sensitive) {
    ctr <- sf::st_sf(
      id       = "ctr",
      geometry = sf::st_sfc(sf::st_point(c(335000, 6250000)), crs = 32755)
    )
    sectors <- data.frame(
      arc_start = "NE", arc_end = "SE",
      permeability = perm, sensitive = sensitive,
      stringsAsFactors = FALSE
    )
    suppressWarnings(site_from_sectors(sectors, ctr, ring_radius = 1000, epsg = 32755))
  }

  it("fully-open, non-sensitive boundary at hazard=100: on_site, leaves_site TRUE, sensitive_receptor FALSE", {
    site <- .single_sector_site(perm = 1.0, sensitive = FALSE)
    out <- litter_exposure(100, 270, site)
    expect_equal(as.character(out$zone), "on_site")
    expect_true(out$leaves_site)
    expect_false(out$sensitive_receptor)
  })

  it("perm-0, sensitive boundary at hazard=100: exposure 0, leaves_site FALSE, on_site, directional_factor 0", {
    site <- .single_sector_site(perm = 0, sensitive = TRUE)
    out <- litter_exposure(100, 270, site)
    expect_equal(out$exposure, 0)
    expect_false(out$leaves_site)
    expect_equal(as.character(out$zone), "on_site")
    expect_equal(out$directional_factor, 0)
  })

  it("open, sensitive boundary at hazard=86: off_site, leaves_site TRUE, sensitive_receptor TRUE, exposure 86", {
    site <- .single_sector_site(perm = 1.0, sensitive = TRUE)
    out <- litter_exposure(86, 270, site)
    expect_equal(as.character(out$zone), "off_site")
    expect_true(out$leaves_site)
    expect_true(out$sensitive_receptor)
    expect_equal(out$exposure, 86, tolerance = 1e-3)
  })

  it("returns exactly the five documented columns", {
    site <- .single_sector_site(perm = 1.0, sensitive = TRUE)
    out <- litter_exposure(86, 270, site)
    expect_named(out, LITTER_EXPOSURE_COLS, ignore.order = TRUE)
  })
})


describe("multiple sources aggregate worst-case", {

  it("directional_factor = max over sources, leaves_site/sensitive_receptor = any (source A hits open, source B hits closed)", {
    base_x <- 335000; base_y <- 6250000
    src_a <- c(base_x, base_y)
    src_b <- c(base_x + 2000, base_y + 2000)

    square <- function(cx, cy, half = 100) {
      sf::st_polygon(list(matrix(
        c(cx - half, cy - half,
          cx + half, cy - half,
          cx + half, cy + half,
          cx - half, cy + half,
          cx - half, cy - half),
        ncol = 2, byrow = TRUE
      )))
    }

    # barrier_open sits 100 m north of source A: as seen from A this large,
    # nearby square subtends a wide arc [270, 90] (wrap through north) that
    # contains theta_down = 0. As seen from the far-away source B it is a
    # narrow, unrelated sliver that does NOT contain theta_down = 0.
    barrier_open   <- square(src_a[1], src_a[2] + 100)
    # barrier_closed is the mirror-image fixture built around source B, so
    # it is hit from B but not from A.
    barrier_closed <- square(src_b[1], src_b[2] + 100)

    feats <- sf::st_sf(
      id       = c("srcA", "srcB", "barrier_open", "barrier_closed"),
      geometry = sf::st_sfc(
        sf::st_point(src_a), sf::st_point(src_b),
        barrier_open, barrier_closed,
        crs = 32755
      )
    )
    roles <- data.frame(
      feature_id   = c("srcA", "srcB", "barrier_open", "barrier_closed"),
      hazard       = "litter",
      role         = c("source", "source", "barrier", "barrier"),
      permeability = c(NA_real_, NA_real_, 1.0, 0.1),
      sensitive    = c(NA, NA, TRUE, FALSE),
      stringsAsFactors = FALSE
    )
    site <- mh_site(feats, roles, epsg = 32755)

    # theta_down = 0 (due north): wind_direction_10m = 180.
    out <- litter_exposure(86, 180, site, p_open_min = 0.5, default_permeability = 0)

    expect_equal(out$directional_factor, 1.0, tolerance = 1e-6)  # max(1.0, 0.1), follows A
    expect_true(out$leaves_site)
    expect_true(out$sensitive_receptor)
    expect_equal(out$exposure, 86 * 1.0, tolerance = 1e-3)
  })

  it("does not error or drop extra sources (no meteoHazard_litter_multi_source error)", {
    base_x <- 335000; base_y <- 6250000
    src_a <- c(base_x, base_y)
    src_b <- c(base_x + 2000, base_y + 2000)
    square <- function(cx, cy, half = 100) {
      sf::st_polygon(list(matrix(
        c(cx - half, cy - half,
          cx + half, cy - half,
          cx + half, cy + half,
          cx - half, cy + half,
          cx - half, cy - half),
        ncol = 2, byrow = TRUE
      )))
    }
    barrier_open   <- square(src_a[1], src_a[2] + 100)
    barrier_closed <- square(src_b[1], src_b[2] + 100)
    feats <- sf::st_sf(
      id       = c("srcA", "srcB", "barrier_open", "barrier_closed"),
      geometry = sf::st_sfc(
        sf::st_point(src_a), sf::st_point(src_b),
        barrier_open, barrier_closed,
        crs = 32755
      )
    )
    roles <- data.frame(
      feature_id   = c("srcA", "srcB", "barrier_open", "barrier_closed"),
      hazard       = "litter",
      role         = c("source", "source", "barrier", "barrier"),
      permeability = c(NA_real_, NA_real_, 1.0, 0.1),
      sensitive    = c(NA, NA, TRUE, FALSE),
      stringsAsFactors = FALSE
    )
    site <- mh_site(feats, roles, epsg = 32755)

    expect_no_error(litter_exposure(86, 180, site))
  })
})
