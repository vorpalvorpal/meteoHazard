# Behaviour specification for the litter exposure layer (litter_exposure()).
#
# Encodes the behaviour defined in specs/Litter_exposure.md (basic mode) and the
# plan in GitHub issue #2. Specs skip until the function exists, then must pass.
#
# Proposed contract (names are the contract):
#   litter_exposure(
#     hazard, wind_direction_10m, site,
#     direction_tol = 15, p_open_min = 0.5,
#     move_threshold = 20, offsite_threshold = 45, default_permeability = 0.5
#   ) -> data.frame(exposure = numeric[0,100], zone = ordered factor)
#   site: data.frame(arc_start, arc_end [compass labels], permeability [0,1],
#                    sensitive [logical], distance_m [optional, unused here])
#   zone levels: within_face < on_site < off_site

skip_if_no_litter_exposure <- function() {
  testthat::skip_if_not(
    exists("litter_exposure", mode = "function"),
    "litter_exposure() not yet implemented"
  )
}

# Demo site (plan section 6): open sensitive boundary to the E, tree belt to the
# W, everything else at the default permeability.
demo_site <- function() {
  data.frame(
    arc_start   = c("NE", "SW"),
    arc_end     = c("SE", "NW"),
    permeability = c(1.0, 0.3),
    sensitive    = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}

EXP_TOL <- 1e-3


describe("litter_exposure()", {

  describe("output contract and bounds", {

    it("returns a data frame with one row per hour and exposure + zone columns", {
      skip_if_no_litter_exposure()
      out <- litter_exposure(
        hazard = c(86, 25, 5), wind_direction_10m = c(270, 90, 0),
        site = demo_site()
      )
      expect_s3_class(out, "data.frame")
      expect_equal(nrow(out), 3)
      expect_true(all(c("exposure", "zone") %in% names(out)))
      expect_type(out$exposure, "double")
    })

    it("returns zone as an ordered factor within_face < on_site < off_site", {
      skip_if_no_litter_exposure()
      out <- litter_exposure(86, 270, demo_site())
      expect_s3_class(out$zone, "ordered")
      expect_identical(levels(out$zone), c("within_face", "on_site", "off_site"))
    })

    it("keeps the exposure within [0, 100]", {
      skip_if_no_litter_exposure()
      g <- expand.grid(h = c(0, 20, 86, 100), d = c(0, 90, 180, 270, 359))
      out <- litter_exposure(g$h, g$d, demo_site())
      expect_true(all(out$exposure >= 0 & out$exposure <= 100))
    })
  })

  describe("directional attenuation (exposure-adjusted hazard)", {

    it("applies the open sensitive sector's permeability when wind blows toward it", {
      skip_if_no_litter_exposure()
      # Wind FROM the west (270) blows downwind toward the E (open, perm 1.0):
      # EH = hazard * 1.0.
      out <- litter_exposure(86, 270, demo_site())
      expect_equal(out$exposure, 86, tolerance = EXP_TOL)
    })

    it("attenuates by the barrier permeability when wind blows toward a barrier", {
      skip_if_no_litter_exposure()
      # Wind FROM the east (90) blows downwind toward the W tree belt (perm 0.3):
      # EH = 86 * 0.3 = 25.8.
      out <- litter_exposure(86, 90, demo_site())
      expect_equal(out$exposure, 25.8, tolerance = EXP_TOL)
    })

    it("falls back to default_permeability when no sector is hit", {
      skip_if_no_litter_exposure()
      # Downwind bearing toward a gap covered by neither demo sector.
      # dir = 0 -> theta_down = 180 (S); demo site covers E (30-150) and
      # W (210-330), so S is uncovered -> default_permeability (0.5).
      out <- litter_exposure(86, 0, demo_site(), default_permeability = 0.5)
      expect_equal(out$exposure, 43, tolerance = EXP_TOL)
    })

    it("takes the most permeable sector when expanded edges overlap two", {
      skip_if_no_litter_exposure()
      # Two sectors meeting at a bearing: worst case (max permeability) wins.
      site <- data.frame(
        arc_start = c("N", "E"), arc_end = c("E", "S"),
        permeability = c(0.2, 0.9), sensitive = c(FALSE, FALSE),
        stringsAsFactors = FALSE
      )
      # theta_down = 90 (E) lies on the shared edge of both sectors.
      out <- litter_exposure(50, 270, site, default_permeability = 0)
      expect_equal(out$exposure, 45, tolerance = EXP_TOL)  # 50 * 0.9
    })
  })

  describe("arc containment with wraparound (corrected expanded-edge rule)", {

    it("admits the tolerance band west of north for a sector starting at N", {
      skip_if_no_litter_exposure()
      # Single sector N->E (admits theta_down in 345-360 and 0-105 with tol 15);
      # non-hit bearings get permeability 0 so exposure reveals hit/miss.
      site <- data.frame(
        arc_start = "N", arc_end = "E", permeability = 1.0, sensitive = FALSE,
        stringsAsFactors = FALSE
      )
      # wind_direction chosen so theta_down = (dir+180)%%360 hits these bearings.
      # theta_down: 345 (dir 165), 0 (dir 180), 105 (dir 285)  -> HIT
      #             340 (dir 160), 110 (dir 290)               -> MISS
      out <- litter_exposure(
        hazard = rep(50, 5),
        wind_direction_10m = c(165, 180, 285, 160, 290),
        site = site, default_permeability = 0
      )
      expect_equal(out$exposure, c(50, 50, 50, 0, 0), tolerance = EXP_TOL)
    })
  })

  describe("severity-zone classification", {

    it("is within_face when the hazard is below the move threshold", {
      skip_if_no_litter_exposure()
      out <- litter_exposure(c(5, 19), c(270, 270), demo_site())
      expect_true(all(out$zone == "within_face"))
    })

    it("is off_site when a strong hazard blows toward a permeable sensitive sector", {
      skip_if_no_litter_exposure()
      out <- litter_exposure(86, 270, demo_site())  # toward open sensitive E
      expect_equal(as.character(out$zone), "off_site")
    })

    it("is on_site when a strong hazard blows toward a non-sensitive barrier", {
      skip_if_no_litter_exposure()
      out <- litter_exposure(86, 90, demo_site())  # toward W tree belt
      expect_equal(as.character(out$zone), "on_site")
    })

    it("does not reach off_site toward a sensitive sector below the offsite threshold", {
      skip_if_no_litter_exposure()
      out <- litter_exposure(30, 270, demo_site())  # sensitive dir, hazard 30 < 45
      expect_equal(as.character(out$zone), "on_site")
    })
  })

  describe("hazard / exposure separation invariant", {

    it("scales exposure linearly with hazard at a fixed direction and site", {
      skip_if_no_litter_exposure()
      # EH = hazard * M, and M depends only on direction+site, so doubling the
      # hazard doubles the exposure.
      out <- litter_exposure(c(40, 80), c(270, 270), demo_site())
      expect_equal(out$exposure[2], 2 * out$exposure[1], tolerance = EXP_TOL)
    })

    it("gives identical exposure-to-hazard ratios for the same direction", {
      skip_if_no_litter_exposure()
      out <- litter_exposure(c(10, 60), c(90, 90), demo_site())
      expect_equal(out$exposure[1] / 10, out$exposure[2] / 60, tolerance = EXP_TOL)
    })
  })

  describe("edge cases and input validation", {

    it("returns zero exposure and within_face for a zero hazard", {
      skip_if_no_litter_exposure()
      out <- litter_exposure(0, 270, demo_site())
      expect_equal(out$exposure, 0)
      expect_equal(as.character(out$zone), "within_face")
    })

    it("rejects hazard outside [0, 100]", {
      skip_if_no_litter_exposure()
      expect_error(litter_exposure(-1, 270, demo_site()))
      expect_error(litter_exposure(101, 270, demo_site()))
    })

    it("rejects wind direction outside [0, 360]", {
      skip_if_no_litter_exposure()
      expect_error(litter_exposure(50, -10, demo_site()))
      expect_error(litter_exposure(50, 400, demo_site()))
    })

    it("rejects mismatched hazard and wind-direction lengths", {
      skip_if_no_litter_exposure()
      expect_error(litter_exposure(c(50, 60), 270, demo_site()))
    })

    it("rejects missing values in hazard or wind direction", {
      skip_if_no_litter_exposure()
      expect_error(litter_exposure(c(50, NA), c(270, 270), demo_site()))
      expect_error(litter_exposure(c(50, 60), c(270, NA), demo_site()))
    })

    it("rejects a malformed site (bad compass label or out-of-range permeability)", {
      skip_if_no_litter_exposure()
      bad_label <- data.frame(
        arc_start = "ENE", arc_end = "E", permeability = 1, sensitive = FALSE,
        stringsAsFactors = FALSE
      )
      bad_perm <- data.frame(
        arc_start = "N", arc_end = "E", permeability = 1.5, sensitive = FALSE,
        stringsAsFactors = FALSE
      )
      expect_error(litter_exposure(50, 270, bad_label))
      expect_error(litter_exposure(50, 270, bad_perm))
    })

    it("rejects a move threshold not below the off-site threshold", {
      skip_if_no_litter_exposure()
      expect_error(litter_exposure(
        50, 270, demo_site(), move_threshold = 50, offsite_threshold = 45
      ))
    })
  })
})

# === C4 behaviour spec (issue #18): litter_exposure() on the mh_site model ===
# Supersedes the compass-sector specs above when the rewrite lands. Pending here.

# Fixture: create an mh_site equivalent to demo_site() for testing.
# demo_site has: NE-SE (perm=1.0, sensitive=TRUE), SW-NW (perm=0.3, sensitive=FALSE)
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

  it("reproduces the current compass-sector outputs via site_from_sectors() with a centroid source", {
    old_out <- litter_exposure(c(86, 25, 5), c(270, 90, 0), demo_site())
    site    <- .make_demo_mh_site()
    new_out <- litter_exposure(c(86, 25, 5), c(270, 90, 0), site)
    expect_equal(new_out$exposure, old_out$exposure, tolerance = 1e-6)
    expect_equal(as.character(new_out$zone), as.character(old_out$zone))
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
  })

  it("errors (classed) when move_threshold is not below offsite_threshold", {
    site <- .make_demo_mh_site()
    expect_error(
      litter_exposure(50, 270, site, move_threshold = 50, offsite_threshold = 45),
      class = "meteoHazard_input_error"
    )
  })

  it("errors (classed) on an invalid mh_site or a missing litter source", {
    # Non-mh_site, non-data.frame -> error
    expect_error(
      litter_exposure(50, 270, "not a site"),
      class = "meteoHazard_input_error"
    )
    # mh_site with no (litter, source) role -> error
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
    offset_src  <- sf::st_point(c(335200, 6250000))  # 200 m east
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
