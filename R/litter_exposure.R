#' Litter exposure layer (direction- and geometry-aware)
#'
#' Maps the direction-agnostic litter hazard ([litter_hazard_vec()]) onto an
#' operational consequence, given the wind direction and the site geometry.
#' Where the hazard index answers *how strongly is litter being entrained and
#' moved*, the exposure layer answers *where does it go and how bad is that*:
#' litter moved within the working face is minor, across the site is bad, and
#' off site is very bad.
#'
#' Basic mode is the coarse arc-containment band.
#' Supplying `mean_wind` and `reach_per_ms` activates the **refined
#' distance-reach mode**: the off-site decision
#' then compares a characteristic downwind reach to the boundary distance rather
#' than thresholding the hazard magnitude.
#'
#' @section Method:
#' For each hour the downwind bearing is the reciprocal of the (meteorological,
#' blows-from) wind direction, \eqn{\theta_{down} = (dir + 180) \bmod 360}. Each
#' boundary sector whose arc (expanded by `direction_tol` on each edge) contains
#' \eqn{\theta_{down}} is "hit"; the directional factor is the largest
#' permeability among hit sectors across **all** litter sources (worst case), or
#' `default_permeability` if the bearing lands in a gap. The exposure-adjusted
#' hazard is `exposure = hazard * directional_factor`.
#'
#' The severity zone and the two boolean destination flags separate the two
#' distinct questions the old single `zone` conflated — *where does litter go*
#' (`leaves_site`) versus *what does it hit* (`sensitive_receptor`):
#' \itemize{
#'   \item `leaves_site` — litter clears a permeable boundary
#'     (permeability >= `p_open_min`) and, in basic mode, the hazard exceeds
#'     `offsite_threshold` (in refined mode, the reach clears the boundary
#'     distance). This is sensitivity-independent.
#'   \item `sensitive_receptor` — a hit sector is `sensitive` **and** open
#'     (permeability >= `p_open_min`).
#'   \item `zone` is `off_site` only when litter both leaves the site **and**
#'     does so toward a sensitive receptor; hence a fully-open but non-sensitive
#'     boundary reports `on_site` with `leaves_site = TRUE`.
#' }
#'
#' @param hazard Numeric vector, non-negative. The relative litter hazard index
#'   from [litter_hazard_vec()], one value per forecast hour. Issue #11 removed
#'   the fixed 0--100 scale, so no upper bound is imposed; `exposure` inherits
#'   the hazard's relative scale.
#' @param wind_direction_10m Numeric vector, degrees `[0, 360]`. Wind direction
#'   at 10 m (meteorological convention: the direction the wind blows *from*).
#'   Same length as `hazard`.
#' @param site An [`mh_site`] S7 object. Must carry a `(litter, source)`
#'   feature and one or more `(litter, barrier)` features with `permeability`
#'   and `sensitive` columns in the roles table (and, for refined mode, a
#'   `distance_m` column on the barrier features). Use [site_from_sectors()] to
#'   build an `mh_site` from a compass-sector data frame.
#'
#' @param direction_tol Tolerance (degrees) added to each arc edge. Default 15.
#' @param p_open_min Minimum permeability for a sector to count as "open" when
#'   deciding `leaves_site`/`sensitive_receptor`. Default 0.5.
#' @param move_threshold Hazard below which litter stays on the working face.
#'   Default 20.
#' @param offsite_threshold Hazard above which litter can clear the boundary in
#'   **basic** mode. Default 45. Must exceed `move_threshold`. (Ignored in
#'   refined mode, which uses the reach test instead.)
#' @param default_permeability Permeability applied when the downwind bearing
#'   matches no configured sector. Default 0.5. A gap is treated consistently for
#'   magnitude and destination: `exposure` uses this permeability, and
#'   `leaves_site` is `TRUE` when it is `>= p_open_min` and the hour is off-site
#'   (basic mode). `sensitive_receptor` is always `FALSE` for a gap, so a gap can
#'   never be `off_site`. In refined mode a gap has no boundary distance, so it is
#'   never "reached" and `leaves_site` is `FALSE` there.
#' @param mean_wind Optional numeric vector (m/s; a bare numeric or a
#'   \pkg{units} object, converted), same length as `hazard`. The hour's mean
#'   wind. Supplying both `mean_wind` and `reach_per_ms` activates the refined
#'   distance-reach mode.
#' @param reach_per_ms Optional positive scalar (metres per m/s): the calibrated
#'   characteristic reach `c_L`. In refined mode a barrier is "reached" when
#'   `reach_per_ms * mean_wind >= distance_m` for that barrier. Note the mean
#'   wind also drives the hazard's transport multiplier `T`
#'   ([litter_hazard_vec()]); `T` amplifies the mobilized flux while this reach
#'   test gates the destination, so the two must be calibrated **jointly**
#'   (issues #11/#26).
#'
#' @return A data frame with `length(hazard)` rows and columns:
#'   \describe{
#'     \item{`exposure`}{Numeric, non-negative: the off-site (permeability-
#'       attenuated) exposure, `hazard * directional_factor`; same relative
#'       scale as `hazard`.}
#'     \item{`zone`}{Ordered factor `within_face` < `on_site` < `off_site`.}
#'     \item{`directional_factor`}{Numeric `[0, 1]`: worst-case permeability
#'       across hit sectors/sources.}
#'     \item{`leaves_site`}{Logical: litter clears a permeable boundary
#'       (destination; sensitivity-independent).}
#'     \item{`sensitive_receptor`}{Logical: a hit sensitive sector is open.}
#'   }
#'
#' @references
#' NSW EPA (2016). \emph{Environmental Guidelines: Solid Waste Landfills}
#' (2nd edn). Barrier-class guidance informing `permeability`.
#'
#' @seealso [litter_hazard_vec()] for the upstream hazard index, [litter_risk()]
#'   for the combined hazard + exposure wrapper.
#' @export
litter_exposure <- function(
  hazard,
  wind_direction_10m,
  site,
  direction_tol        = 15,
  p_open_min           = 0.5,
  move_threshold       = 20,
  offsite_threshold    = 45,
  default_permeability = 0.5,
  mean_wind            = NULL,
  reach_per_ms         = NULL
) {
  if (!S7::S7_inherits(site, mh_site)) {
    cli::cli_abort(
      "{.arg site} must be an {.cls mh_site}.",
      class = "meteoHazard_input_error"
    )
  }

  # ---- Validate hazard / direction vectors --------------------------------- #
  # Contract seam: the hazard is a RELATIVE index (issue #11 removed the
  # fixed 0-100 scale), so only non-negativity is asserted here. Imposing an
  # upper bound of 100 would make the exposure layer reject legitimate hazard
  # values produced with a non-default entrainment_max.
  n <- length(hazard)
  checkmate::assert_numeric(hazard, lower = 0, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(wind_direction_10m, lower = 0, upper = 360,
                            any.missing = FALSE, len = n)

  # ---- Validate scalar parameters ------------------------------------------ #
  # direction_tol plays the role of the odour module's forecast wind-direction
  # uncertainty SIGMA_FC_DEG (constants.R): it widens each barrier arc to absorb
  # NWP direction error plus lateral spread. Litter uses a HARD tolerance band
  # rather than odour's Gaussian smoothing of sigma_y because this is a coarse
  # screening arc-containment test, not a plume-dispersion calculation.
  checkmate::assert_number(direction_tol, lower = 0, upper = 90)
  checkmate::assert_number(p_open_min, lower = 0, upper = 1)
  checkmate::assert_number(move_threshold, lower = 0)
  checkmate::assert_number(offsite_threshold)
  if (move_threshold >= offsite_threshold) {
    cli::cli_abort(
      c(
        "{.arg move_threshold} ({move_threshold}) must be less than {.arg offsite_threshold} ({offsite_threshold}).",
        "i" = "The zone ladder requires within_face < on_site < off_site."
      ),
      class = "meteoHazard_input_error"
    )
  }
  checkmate::assert_number(default_permeability, lower = 0, upper = 1)

  # ---- Refined distance-reach mode ------------------------------------------ #
  # Active only when BOTH mean_wind and reach_per_ms are supplied. It replaces
  # the hazard-magnitude off-site test with a geometric reach test: a
  # characteristic downwind reach L = reach_per_ms * mean_wind is compared to the
  # downwind distance to the boundary. Deliberately crude (a single linear reach
  # scale, not a trajectory); the natural home for the transport/reach question
  # that the hazard layer's transport term only crudely proxies.
  refined <- !is.null(mean_wind) && !is.null(reach_per_ms)
  if (refined) {
    mean_wind <- .drop_to(mean_wind, "m/s", arg = "mean_wind")
    checkmate::assert_numeric(mean_wind, lower = 0, any.missing = FALSE, len = n)
    checkmate::assert_number(reach_per_ms, lower = .Machine$double.eps)
  }

  # ---- Get source and barrier features ------------------------------------- #
  sources  <- .role_features(site, "litter", "source")
  barriers <- .role_features(site, "litter", "barrier")

  if (nrow(sources) == 0) {
    cli::cli_abort(
      "{.arg site} has no {.val (litter, source)} role.",
      class = "meteoHazard_input_error"
    )
  }

  # A site with a source but no barriers is almost always mis-configured: every
  # bearing falls through to default_permeability and no hour can ever be
  # off-site. Warn rather than error (a bare site is a legitimate, if unusual,
  # screening input).
  if (nrow(barriers) == 0) {
    cli::cli_warn(
      c(
        "{.arg site} has a litter source but no {.val (litter, barrier)} features.",
        "i" = "Every bearing uses {.arg default_permeability} for both exposure and {.field leaves_site}, is never a sensitive receptor, and so no hour can be off-site."
      ),
      class = "meteoHazard_litter_no_barriers"
    )
  }

  # Barrier roles (permeability + sensitive live here)
  barrier_roles <- site@roles[
    site@roles$hazard == "litter" & site@roles$role == "barrier",
  ]

  # ---- Downwind bearing ---------------------------------------------------- #
  theta_down <- .downwind_bearing(wind_direction_10m)

  # ---- Aggregate over ALL sources and barriers (worst-case screening) ------ #
  # Worst-case-across-sources is the natural screening semantics for multiple
  # active tipping faces (mirrors odour_exposure()'s loop over sources): take the
  # most permeable hit, and OR the destination/sensitivity flags.
  best_perm     <- rep(-Inf, n)
  sensitive_hit <- logical(n)
  leaves_site   <- logical(n)

  for (s in seq_len(nrow(sources))) {
    source_pt <- sources[s, ]

    for (k in seq_len(nrow(barriers))) {
      barrier_k <- barriers[k, ]

      # Bearing arc from THIS source to the barrier.
      br    <- .barrier_bearing_range(source_pt, barrier_k)
      hit_k <- .litter_arc_contains(theta_down, br["alpha"], br["beta"], direction_tol)

      # Look up permeability / sensitivity for this barrier.
      role_k <- barrier_roles[barrier_roles$feature_id == barrier_k$id, ]
      perm_k <- if (nrow(role_k) > 0 &&
                    "permeability" %in% names(role_k) &&
                    !is.na(role_k$permeability[1])) {
        role_k$permeability[1]
      } else {
        default_permeability
      }
      sens_k <- if (nrow(role_k) > 0 &&
                    "sensitive" %in% names(role_k) &&
                    !is.na(role_k$sensitive[1])) {
        as.logical(role_k$sensitive[1])
      } else {
        FALSE
      }

      # Downwind distance to this barrier (refined mode only). Use the barrier
      # feature's distance_m attribute (site_from_sectors sets it to the ring
      # radius); absent/NA -> Inf, i.e. never reached. (st_distance() to a
      # wedge polygon is 0 because the polygon includes the source vertex, so the
      # configured distance_m attribute is the reliable source of the reach
      # distance.)
      dist_k <- if ("distance_m" %in% names(barrier_k) &&
                    !is.na(barrier_k$distance_m[1])) {
        barrier_k$distance_m[1]
      } else {
        Inf
      }

      open_k    <- perm_k >= p_open_min
      # "reached": refined -> reach clears the boundary distance; basic -> hazard
      # exceeds the off-site threshold.
      reached_k <- if (refined) {
        reach_per_ms * mean_wind >= dist_k
      } else {
        hazard >= offsite_threshold
      }

      best_perm     <- pmax(best_perm, ifelse(hit_k, perm_k, -Inf))
      sensitive_hit <- sensitive_hit | (hit_k & sens_k & open_k)
      leaves_site   <- leaves_site   | (hit_k & open_k & reached_k)
    }
  }

  # ---- Gap fall-through (no sector hit) ------------------------------------ #
  # On a gap hour (best_perm still -Inf) the magnitude already uses
  # default_permeability, so leaves_site must use the SAME assumptions:
  # permeability = default_permeability, sensitivity = FALSE, distance = unknown.
  # Previously the magnitude said "litter passes the unconfigured boundary" while
  # leaves_site said "never" -- a contradiction. sensitive_receptor stays FALSE
  # for gaps (an unconfigured aspect is not a known sensitive receptor), so zone
  # can still never be off_site through a gap. Refined mode has no distance for a
  # gap, so it can never be "reached" there.
  gap_hours   <- !is.finite(best_perm)
  gap_open    <- default_permeability >= p_open_min
  gap_reached <- if (refined) FALSE else hazard >= offsite_threshold
  leaves_site <- leaves_site | (gap_hours & gap_open & gap_reached)

  directional_factor <- ifelse(is.finite(best_perm), best_perm, default_permeability)

  # ---- Exposure-adjusted hazard and severity zone -------------------------- #
  exposure <- hazard * directional_factor

  # within_face: below move_threshold, litter is not mobile enough to leave the
  # face. off_site: clears a permeable boundary AND does so toward a sensitive
  # receptor. on_site: everything else (mobile but contained, or leaving toward a
  # non-sensitive aspect).
  zone <- ifelse(
    hazard < move_threshold,
    "within_face",
    ifelse(leaves_site & sensitive_hit, "off_site", "on_site")
  )
  zone <- factor(zone, levels = c("within_face", "on_site", "off_site"),
                 ordered = TRUE)

  data.frame(
    exposure           = exposure,
    zone               = zone,
    directional_factor = directional_factor,
    leaves_site        = leaves_site,
    sensitive_receptor = sensitive_hit
  )
}


# ---- Internal helpers ------------------------------------------------------ #

# Compass label -> degrees (meteorological, blows-FROM convention). Eight
# principal points only.
LITTER_COMPASS_DEGREES <- c(
  "N"  =   0,
  "NE" =  45,
  "E"  =  90,
  "SE" = 135,
  "S"  = 180,
  "SW" = 225,
  "W"  = 270,
  "NW" = 315
)

# Does bearing theta (degrees) fall within the arc [alpha, beta] (clockwise),
# expanded by tol on each edge? Branches on the *expanded* edges so the
# tolerance band is handled correctly even when it crosses north.
.litter_arc_contains <- function(theta, alpha, beta, tol) {
  # Full-circle sentinel: .barrier_bearing_range() returns (alpha = 0,
  # beta = 360) for a barrier that ENCLOSES the source, meaning every bearing is
  # hit. The expanded-edge logic below would instead collapse (0, 360) to only
  # +/- tol of north (alpha_exp = 345, beta_exp = 15), so short-circuit here.
  if (beta - alpha >= 360) {
    return(rep(TRUE, length(theta)))
  }
  # Wide-arc guard: when the clockwise span of the raw arc plus the tolerance
  # band on both edges meets a full turn (e.g. a horseshoe barrier spanning
  # 340 deg with tol = 15), the expanded edges pass each other and the wrap
  # branch below would test the COMPLEMENT of the arc — a bearing pointing
  # straight into the barrier would be reported as a miss.
  span <- (beta - alpha) %% 360
  if (span + 2 * tol >= 360) {
    return(rep(TRUE, length(theta)))
  }
  alpha_exp <- (alpha - tol) %% 360
  beta_exp  <- (beta  + tol) %% 360
  if (alpha_exp <= beta_exp) {
    theta >= alpha_exp & theta <= beta_exp
  } else {
    theta >= alpha_exp | theta <= beta_exp
  }
}

# Compute the smallest clockwise arc [alpha, beta] (degrees) that encloses all
# bearings from source_pt to the vertices of barrier_poly. Returns a named
# numeric vector c(alpha = <start>, beta = <end>).
#
# The algorithm: compute all vertex bearings, sort them, find the largest
# clockwise gap between consecutive bearings, and declare the complement as the
# containing arc.
#
# Vertices coincident with the source (distance 0) are excluded — they arise
# when the barrier polygon was built with the source centroid as a vertex (as
# in site_from_sectors()).
.barrier_bearing_range <- function(source_pt, barrier_poly) {
  # Enclosure guard: if the source lies INSIDE the barrier polygon, the barrier
  # surrounds it and is hit from every bearing. The largest-gap heuristic below
  # cannot represent that (it always returns a sub-360 arc), so short-circuit to
  # the full-circle sentinel (0, 360). Verified empirically: a site_from_sectors
  # wedge has the source as an apex BOUNDARY vertex (st_within FALSE -> guard
  # stays silent for ordinary directional barriers); a disk / solid enclosing
  # polygon has the source in the INTERIOR (st_within TRUE -> guard fires).
  if (length(sf::st_within(source_pt, barrier_poly)[[1]]) > 0) {
    return(c(alpha = 0, beta = 360))
  }

  src_xy <- sf::st_coordinates(source_pt)[1, c("X", "Y")]

  # Densify the boundary before taking bearings. The largest-gap heuristic sees
  # only polygon VERTICES, so a coarse CONCAVE barrier can put a larger gap
  # between two consecutive vertex bearings INSIDE the true arc than the real
  # outside gap, returning the complement (a bearing pointing straight at a wall
  # reported as a miss). Segmentizing to <= ~3 deg angular steps (relative to the
  # nearest boundary point) guarantees the largest gap is the true outside gap.
  # Straight radial edges of site_from_sectors wedges are unaffected: bearings
  # vary monotonically along a segment, so densifying introduces no new extreme
  # bearings and existing results are unchanged.
  verts0 <- sf::st_coordinates(barrier_poly)[, c("X", "Y")]
  d2     <- (verts0[, "X"] - src_xy["X"])^2 + (verts0[, "Y"] - src_xy["Y"])^2
  d_min  <- sqrt(min(d2[d2 > 0]))    # nearest non-coincident vertex
  seg_len <- max(d_min * 0.05, 0.5)  # <= ~2.9 deg step at range d_min; 0.5 m floor
  poly_dense <- sf::st_segmentize(sf::st_geometry(barrier_poly), dfMaxLength = seg_len)
  verts  <- sf::st_coordinates(poly_dense)[, c("X", "Y")]

  dE <- verts[, "X"] - src_xy["X"]
  dN <- verts[, "Y"] - src_xy["Y"]

  # Exclude coincident vertices (distance == 0); their bearing is undefined
  dist2 <- dE^2 + dN^2
  keep  <- dist2 > 0
  dE <- dE[keep]
  dN <- dN[keep]

  brgs <- (atan2(dE, dN) * 180 / pi) %% 360

  # Deduplicate and sort
  brgs_sorted <- sort(unique(round(brgs, 6)))

  if (length(brgs_sorted) == 1) {
    return(c(alpha = brgs_sorted[1], beta = brgs_sorted[1]))
  }

  # Clockwise gaps between consecutive bearings (with wrap-around)
  n_brgs <- length(brgs_sorted)
  gaps <- diff(c(brgs_sorted, brgs_sorted[1] + 360))

  # The largest gap is the "outside" of the arc; the arc is its complement
  largest_gap_idx <- which.max(gaps)

  # Arc starts at the bearing just after the largest gap, ends at the bearing
  # just before it.
  alpha_idx <- (largest_gap_idx %% n_brgs) + 1
  beta_idx  <- largest_gap_idx

  c(alpha = brgs_sorted[alpha_idx], beta = brgs_sorted[beta_idx])
}
