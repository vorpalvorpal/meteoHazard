#' Litter exposure layer (direction- and geometry-aware)
#'
#' Maps the direction-agnostic litter hazard ([litter_risk_index()]) onto an
#' operational consequence, given the wind direction and the site geometry.
#' Where the hazard index answers *how strongly is litter being entrained and
#' moved*, the exposure layer answers *where does it go and how bad is that*:
#' litter moved within the working face is minor, across the site is bad, and
#' off site is very bad.
#'
#' This is the "basic mode" of `specs/Litter_exposure.md`: a coarse band, not a
#' dispersion or trajectory model. The refined distance-reach mode is not
#' implemented.
#'
#' @section Method:
#' For each hour the downwind bearing is the reciprocal of the (meteorological,
#' blows-from) wind direction, \eqn{\theta_{down} = (dir + 180) \bmod 360}. Each
#' boundary sector whose arc (expanded by `direction_tol` on each edge) contains
#' \eqn{\theta_{down}} is "hit"; the directional factor `M` is the largest
#' permeability among hit sectors (worst case), or `default_permeability` if the
#' bearing lands in a gap. The exposure-adjusted hazard is `exposure = hazard *
#' M`. The severity zone uses the raw hazard for mobility and the hit sectors for
#' destination.
#'
#' @param hazard Numeric vector in `[0, 100]`. The litter hazard index from
#'   [litter_risk_index()], one value per forecast hour.
#' @param wind_direction_10m Numeric vector, degrees `[0, 360]`. Wind direction
#'   at 10 m (meteorological convention: the direction the wind blows *from*).
#'   Same length as `hazard`.
#' @param site A data frame describing the site boundary as sectors, one row per
#'   sector, with columns:
#'   \describe{
#'     \item{`arc_start`, `arc_end`}{Compass labels (one of N, NE, E, SE, S, SW,
#'       W, NW) giving the clockwise start and end of the boundary sector.}
#'     \item{`permeability`}{Numeric `[0, 1]`: 1 = open / no effective barrier,
#'       lower = better containment (e.g. engineered wall ~0.3).}
#'     \item{`sensitive`}{Logical: `TRUE` if off-site impact in that direction
#'       matters (a receptor).}
#'   }
#'   An optional `distance_m` column is accepted but unused in basic mode.
#' @param direction_tol Tolerance (degrees) added to each arc edge. Default 15.
#' @param p_open_min Minimum permeability for a sensitive sector to count as
#'   "open" when deciding off-site risk. Default 0.5.
#' @param move_threshold Hazard below which litter stays on the working face.
#'   Default 20.
#' @param offsite_threshold Hazard above which litter can clear the boundary.
#'   Default 45. Must exceed `move_threshold`.
#' @param default_permeability Permeability applied when the downwind bearing
#'   matches no configured sector. Default 0.5.
#'
#' @return A data frame with `length(hazard)` rows and columns:
#'   \describe{
#'     \item{`exposure`}{Numeric `[0, 100]`: hazard attenuated by the directional
#'       factor.}
#'     \item{`zone`}{Ordered factor `within_face` < `on_site` < `off_site`.}
#'   }
#'
#' @references
#' NSW EPA (2016). \emph{Environmental Guidelines: Solid Waste Landfills}
#' (2nd edn). Barrier-class guidance informing `permeability`.
#'
#' @seealso [litter_risk_index()] for the upstream hazard index.
#' @export
litter_exposure <- function(
  hazard,
  wind_direction_10m,
  site,
  direction_tol        = 15,
  p_open_min           = 0.5,
  move_threshold       = 20,
  offsite_threshold    = 45,
  default_permeability = 0.5
) {

  # ---- Validate hazard / direction vectors --------------------------------- #
  n <- length(hazard)
  checkmate::assert_numeric(hazard, lower = 0, upper = 100,
                            any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(wind_direction_10m, lower = 0, upper = 360,
                            any.missing = FALSE, len = n)

  # ---- Validate the site configuration ------------------------------------- #
  checkmate::assert_data_frame(site, min.rows = 1)
  required_cols <- c("arc_start", "arc_end", "permeability", "sensitive")
  missing_cols <- setdiff(required_cols, names(site))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg site} is missing required columns: {.val {missing_cols}}."
    )
  }
  valid_labels <- names(LITTER_COMPASS_DEGREES)
  bad_labels <- setdiff(c(site$arc_start, site$arc_end), valid_labels)
  if (length(bad_labels) > 0) {
    cli::cli_abort(c(
      "{.arg site} contains invalid compass label(s): {.val {bad_labels}}.",
      "i" = "Valid labels are: {.val {valid_labels}}."
    ))
  }
  checkmate::assert_numeric(site$permeability, lower = 0, upper = 1,
                            any.missing = FALSE, .var.name = "site$permeability")
  checkmate::assert_logical(site$sensitive, any.missing = FALSE,
                            .var.name = "site$sensitive")

  # ---- Validate scalar parameters ------------------------------------------ #
  checkmate::assert_number(direction_tol, lower = 0, upper = 90)
  checkmate::assert_number(p_open_min, lower = 0, upper = 1)
  checkmate::assert_number(move_threshold, lower = 0)
  checkmate::assert_number(offsite_threshold)
  if (move_threshold >= offsite_threshold) {
    cli::cli_abort(c(
      "{.arg move_threshold} ({move_threshold}) must be less than {.arg offsite_threshold} ({offsite_threshold}).",
      "i" = "The zone ladder requires within_face < on_site < off_site."
    ))
  }
  checkmate::assert_number(default_permeability, lower = 0, upper = 1)

  # ---- Downwind bearing (reciprocal of the blows-from direction) ----------- #
  theta_down <- .downwind_bearing(wind_direction_10m)

  alpha <- unname(LITTER_COMPASS_DEGREES[site$arc_start])
  beta  <- unname(LITTER_COMPASS_DEGREES[site$arc_end])

  # ---- Directional factor and sensitive-hit flag --------------------------- #
  # Loop over the (few) sectors, vectorising each arc test across all hours.
  # best_perm holds the most permeable hit sector per hour (worst case across
  # overlaps); -Inf marks an hour no sector covers, which falls back to the
  # default permeability.
  best_perm     <- rep(-Inf, n)
  sensitive_hit <- logical(n)
  for (k in seq_len(nrow(site))) {
    hit_k <- .litter_arc_contains(theta_down, alpha[k], beta[k], direction_tol)
    best_perm <- pmax(best_perm, ifelse(hit_k, site$permeability[k], -Inf))
    sensitive_hit <- sensitive_hit |
      (hit_k & site$sensitive[k] & site$permeability[k] >= p_open_min)
  }
  directional_factor <- ifelse(is.finite(best_perm), best_perm, default_permeability)

  # ---- Exposure-adjusted hazard and severity zone -------------------------- #
  exposure <- hazard * directional_factor

  zone <- ifelse(
    hazard < move_threshold,
    "within_face",
    ifelse(hazard >= offsite_threshold & sensitive_hit, "off_site", "on_site")
  )
  zone <- factor(zone, levels = c("within_face", "on_site", "off_site"),
                 ordered = TRUE)

  data.frame(exposure = exposure, zone = zone)
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
# tolerance band is handled correctly even when it crosses north (the corrected
# rule from specs/Litter.md S3.2).
.litter_arc_contains <- function(theta, alpha, beta, tol) {
  alpha_exp <- (alpha - tol) %% 360
  beta_exp  <- (beta  + tol) %% 360
  if (alpha_exp <= beta_exp) {
    theta >= alpha_exp & theta <= beta_exp
  } else {
    theta >= alpha_exp | theta <= beta_exp
  }
}
