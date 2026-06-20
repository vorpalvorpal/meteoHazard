# ===========================================================================
# Internal site-geometry helpers. All geometry is computed in a projected
# metric CRS (Euclidean planar trigonometry). These functions are internal
# (dot-prefixed) and not exported.
# ===========================================================================

# Bearing convention: degrees from north, clockwise.
# In a projected CRS: bearing = (atan2(dE, dN) * 180/pi) %% 360
# where dE = easting difference, dN = northing difference.

# ---------------------------------------------------------------------------
# .bearing_distance(from, to)
# ---------------------------------------------------------------------------
# Returns a list(bearing = <deg [0,360) or NA>, distance = <m>) between two
# sf POINT geometries. Accepts sf, sfc, or sfg POINT objects.
# Coincident points → distance 0, bearing NA.

.bearing_distance <- function(from, to) {
  coords_from <- sf::st_coordinates(from)
  coords_to   <- sf::st_coordinates(to)

  dE <- unname(coords_to[1, "X"] - coords_from[1, "X"])   # easting difference
  dN <- unname(coords_to[1, "Y"] - coords_from[1, "Y"])   # northing difference

  distance <- sqrt(dE^2 + dN^2)

  if (distance == 0) {
    return(list(bearing = NA_real_, distance = 0))
  }

  bearing <- (atan2(dE, dN) * 180 / pi) %% 360

  list(bearing = bearing, distance = distance)
}

# ---------------------------------------------------------------------------
# .relative_elevation(source, receptor)
# ---------------------------------------------------------------------------
# Returns receptor elevation minus source base elevation (Δz), both on the
# AGL datum. Reads the `elevation` column from the sf features passed in.
# Missing elevation propagates as NA.

.relative_elevation <- function(source, receptor) {
  src_elev  <- if ("elevation" %in% names(source))  source$elevation[[1]]  else NA_real_
  recv_elev <- if ("elevation" %in% names(receptor)) receptor$elevation[[1]] else NA_real_

  recv_elev - src_elev
}

# ---------------------------------------------------------------------------
# .role_features(site, hazard, role)
# ---------------------------------------------------------------------------
# Returns the sf subset of site@features whose site@roles match the given
# hazard and role. An absent hazard/role → zero-row sf (not an error).

.role_features <- function(site, hazard, role) {
  matched_ids <- site@roles$feature_id[
    site@roles$hazard == hazard & site@roles$role == role
  ]
  site@features[site@features$id %in% matched_ids, ]
}

# ---------------------------------------------------------------------------
# .crosswind_halfwidth(polygon, wind_dir)
# ---------------------------------------------------------------------------
# Returns half the footprint extent measured perpendicular to the wind
# direction. Projects the polygon vertices onto the crosswind unit vector
# (perpendicular to wind_dir) and takes half the range.
# >= 0; symmetric under wind_dir +/- 180; degenerate point geometry → 0.

.crosswind_halfwidth <- function(polygon, wind_dir) {
  coords <- sf::st_coordinates(polygon)

  # Crosswind direction is perpendicular to the wind (rotated 90 deg CW).
  # Wind direction is degrees from north, clockwise.
  # Wind unit vector (downwind): (sin(wind_rad), cos(wind_rad)) in (E, N)
  # Crosswind unit vector (90 deg CW of downwind): (cos(wind_rad), -sin(wind_rad))
  wind_rad <- wind_dir * pi / 180

  # Project each vertex onto the crosswind unit vector
  cw_projections <- coords[, "X"] * cos(wind_rad) - coords[, "Y"] * sin(wind_rad)

  extent <- diff(range(cw_projections))
  extent / 2
}
