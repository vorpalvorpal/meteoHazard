# ===========================================================================
# Shared wind-geometry helpers used by the odour and litter exposure layers.
# Single-sourcing the meteorological (blows-from) direction convention so the
# modules cannot silently diverge.
# ===========================================================================

# Downwind (plume-travel) bearing: the reciprocal of the meteorological
# "blows-from" wind direction. Vectorised; NA propagates.
.downwind_bearing <- function(wind_direction) {
  (wind_direction + 180) %% 360
}

# Shortest angular separation between two bearings (degrees), in [0, 180].
# Vectorised over either argument; NA propagates.
.angular_diff <- function(a, b) {
  d <- abs(a - b) %% 360
  pmin(d, 360 - d)
}
