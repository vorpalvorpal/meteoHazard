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
  coords   <- sf::st_coordinates(polygon)
  wind_rad <- wind_dir * pi / 180

  # Crosswind unit vector (90 deg CW of downwind): (cos(wr), -sin(wr)) in (E, N).
  # Project each vertex: p[v] = X[v]*cos(wr) - Y[v]*sin(wr).
  if (length(wind_rad) == 1L) {
    # Scalar path — single direction; simple vector operation.
    proj <- coords[, "X"] * cos(wind_rad) - coords[, "Y"] * sin(wind_rad)
    diff(range(proj)) / 2
  } else {
    # Vector path — one direction per hour; extract coordinates only once.
    # proj_mat[v, t] = X[v]*cos(wr[t]) - Y[v]*sin(wr[t])
    proj_mat <- outer(coords[, "X"], cos(wind_rad)) -
                outer(coords[, "Y"], sin(wind_rad))
    (apply(proj_mat, 2L, max) - apply(proj_mat, 2L, min)) / 2
  }
}

# ---------------------------------------------------------------------------
# .receptor_delta_z(source_k, receptors, emit_ht)
# ---------------------------------------------------------------------------
# Returns a numeric vector of length n_r: Δz[j] = receptor_elevation[j] - emit_ht.
# NA elevation → Δz = 0 (receptor unaffected by M2). emit_ht is the source's
# emit_height (already extracted and defaulted to 0).

.receptor_delta_z <- function(source_k, receptors, emit_ht) {
  if (!("elevation" %in% names(receptors))) return(rep(0.0, nrow(receptors)))
  recv_elev <- receptors$elevation
  recv_elev <- ifelse(is.na(recv_elev), emit_ht, recv_elev)
  recv_elev - emit_ht
}

# ---------------------------------------------------------------------------
# .receptor_hill_height_scale(receptors)
# ---------------------------------------------------------------------------
# Returns a numeric vector of length n_r: hill_height_scale[j] ∈ [0,1].
# If the column is absent or NA → 0 (pure stability blend).

.receptor_hill_height_scale <- function(receptors) {
  if (!("hill_height_scale" %in% names(receptors))) return(rep(0.0, nrow(receptors)))
  h <- receptors$hill_height_scale
  ifelse(is.na(h), 0.0, h)
}
