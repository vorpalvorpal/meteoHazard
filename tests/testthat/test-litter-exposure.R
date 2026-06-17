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
