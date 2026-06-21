# ===========================================================================
# C3b terrain morning-pulse pathway helpers.
#
# Four internal helpers used by odour_exposure() when
# terrain_backend = "descriptors":
#
#   .pool_partition()   — split emission into above-pool (1b) / within-pool (1a)
#   .cw_venting()       — directional crosswind factor for pathway 1a
#   .cw_fumigation()    — directional crosswind factor for pathway 1b
#   .morning_release()  — multiplicative release enhancement for morning pulse
#
# All are internal (dot-prefix, not exported).
# ===========================================================================


# ---------------------------------------------------------------------------
# .pool_partition(emit_extent, pool_top, delta)
# ---------------------------------------------------------------------------
# Partitions a uniform vertical emission profile E(z) over [0, emit_extent]
# into:
#   f_1b = above-pool fraction  (pool_top acts as the split point)
#   f_1a = within-pool fraction = 1 - f_1b
#
# Physics:
#   Phi(u) = 1 / (1 + exp(-u))   logistic smooth step
#   f_1b   = (1 / emit_extent) * integral_0^emit_extent  Phi((z - pool_top)/delta) dz
#           approximated via 50-point quadrature.
#
# All arguments are vectorised (per-hour).
#
# Returns a list(f_1a, f_1b) with f_1a + f_1b == 1 for every element.

.pool_partition <- function(emit_extent, pool_top, delta) {
  n <- max(length(emit_extent), length(pool_top), length(delta))
  emit_extent <- rep_len(emit_extent, n)
  pool_top    <- rep_len(pool_top,    n)
  delta       <- rep_len(delta,       n)

  phi <- function(u) 1 / (1 + exp(-pmin(pmax(u, -30), 30)))

  # NA in any argument → f_1b = 0 (handled after the vectorised quadrature).
  na_mask <- is.na(emit_extent) | is.na(pool_top) | is.na(delta)

  # Substitute safe values so the quadrature never sees NA; emit_extent <= 0 is
  # the point-source case, which the quadrature reproduces exactly (all 50
  # samples collapse to z = 0, so mean phi = phi(-pool_top / delta)).
  e_eff <- pmax(ifelse(na_mask, 0, emit_extent), 0)
  pt    <- ifelse(na_mask, 0, pool_top)
  d     <- pmax(ifelse(na_mask, 1, delta), 1e-6)

  # 50-point quadrature of phi((z - pool_top)/delta) over z in [0, emit_extent].
  # With z = emit_extent * u and u = linspace(0, 1, 50), the sample grid is
  # identical to seq(0, emit_extent, length.out = 50). Build the (50 x n)
  # argument matrix and average down columns.
  #   arg[q, i] = (e_eff[i] * u[q] - pt[i]) / d[i]
  u   <- seq(0, 1, length.out = 50)
  arg <- (outer(u, e_eff) -
            matrix(pt, nrow = 50L, ncol = n, byrow = TRUE)) /
            matrix(d,  nrow = 50L, ncol = n, byrow = TRUE)
  f_1b <- colMeans(phi(arg))

  f_1b[na_mask] <- 0.0
  f_1b <- pmin(pmax(f_1b, 0), 1)
  list(f_1a = 1 - f_1b, f_1b = f_1b)
}


# ---------------------------------------------------------------------------
# .cw_venting(bearing_to_receptor, vs, terrain, CONFINEMENT_1A, VENTING_1A)
# ---------------------------------------------------------------------------
# Directional crosswind factor for pathway 1a (within-pool venting / anabatic).
#
# Night  (is_day == FALSE):
#   Channelled (flow_convergence >= 0.5 AND drainage_bearing non-NA):
#     cw_1a = CONFINEMENT_1A * (1 - 0.9 * alignment)
#   Radial:
#     cw_1a = CONFINEMENT_1A
#
# Morning transition (is_day == TRUE AND cbl_growth > 0):
#   Channelled:
#     cw_1a = VENTING_1A * max(0, cos(diff_to_drainage))
#   Radial:
#     cw_1a = VENTING_1A   (all directions boosted)
#
# Otherwise: NA  (caller uses flat Gaussian)
#
# bearing_to_receptor: scalar degrees (single receptor)
# vs: list from ventilation_state()
# terrain: mh_terrain or NULL

# Matrix form: `bearings` is a length-n_r vector of receptor bearings; returns
# an (n_t x n_r) matrix. The alignment term is the only receptor-dependent
# quantity, so each hour-mask row is filled with the per-receptor vector.
.cw_venting_matrix <- function(bearings, vs, terrain,
                               CONFINEMENT_1A = ODOUR_CONSTANTS$CONFINEMENT_1A,
                               VENTING_1A     = ODOUR_CONSTANTS$VENTING_1A) {
  n_t <- length(vs$is_day)
  n_r <- length(bearings)

  drain_bearing <- if (!is.null(terrain) && !is.na(terrain@drainage_bearing))
    terrain@drainage_bearing else NA_real_
  flow_conv <- if (!is.null(terrain) && !is.na(terrain@flow_convergence))
    terrain@flow_convergence else 0

  is_channelled <- !is.na(drain_bearing) && flow_conv >= 0.5

  is_night   <- !vs$is_day
  is_morning <- vs$is_day & !is.na(vs$cbl_growth) & vs$cbl_growth > 0

  cw <- matrix(NA_real_, n_t, n_r)
  if (is_channelled) {
    alignment   <- pmax(0, cos(.angular_diff(bearings, drain_bearing) * pi / 180))
    night_val   <- CONFINEMENT_1A * (1 - 0.9 * alignment)
    morning_val <- VENTING_1A * alignment
    if (any(is_night))
      cw[is_night, ]   <- matrix(night_val,   sum(is_night),   n_r, byrow = TRUE)
    if (any(is_morning))
      cw[is_morning, ] <- matrix(morning_val, sum(is_morning), n_r, byrow = TRUE)
  } else {
    if (any(is_night))   cw[is_night, ]   <- CONFINEMENT_1A
    if (any(is_morning)) cw[is_morning, ] <- VENTING_1A
  }
  # Daytime hours with no cbl_growth > 0 remain NA (caller uses flat Gaussian).
  cw
}

# Scalar form (single receptor bearing): thin wrapper over the matrix form so
# the two share one implementation. Returns a length-n_t vector.
.cw_venting <- function(bearing_to_receptor, vs, terrain,
                        CONFINEMENT_1A = ODOUR_CONSTANTS$CONFINEMENT_1A,
                        VENTING_1A     = ODOUR_CONSTANTS$VENTING_1A) {
  .cw_venting_matrix(bearing_to_receptor, vs, terrain,
                     CONFINEMENT_1A, VENTING_1A)[, 1L]
}


# ---------------------------------------------------------------------------
# .cw_fumigation(bearing_to_receptor, vs, terrain, FUMIC_1B)
# ---------------------------------------------------------------------------
# Directional crosswind factor for pathway 1b (above-pool fumigation).
#
# Only active during morning transition (is_day == TRUE AND cbl_growth > 0)
# when pool_top > 0.
#
# Fumigation travels DOWNWIND of the overnight residual layer wind.
# Residual-wind direction fallback ladder:
#   dir_180m → dir_120m → dir_80m → dir_10m (overnight circular mean)
#   → surface wind direction (instantaneous, if available in vs)
#   → cw_1b = 0 (truly no information).
#
# cw_1b = FUMIC_1B * max(0, cos(angular_diff(bearing_to_receptor, downwind)))^2
#
# The residual_wind direction is meteorological convention (blows FROM),
# so downwind = (rw_dir + 180) %% 360.

# ---------------------------------------------------------------------------
# .cw_fumigation_prep(vs, above_pool_ht)
# ---------------------------------------------------------------------------
# Expensive setup for .cw_fumigation(): forward-fill residual-wind directions
# and select the level priority order. This part does NOT depend on the
# receptor bearing, so it should be called ONCE per source (not per receptor).
# Pass the returned prep list to .cw_fumigation() via the `prep` argument.

.cw_fumigation_prep <- function(vs, above_pool_ht = 0) {
  n_t <- length(vs$is_day)
  rw  <- vs$residual_wind
  wind_dir_surf <- if (!is.null(vs$wind_dir_surface)) vs$wind_dir_surface
                   else rep(NA_real_, n_t)

  # Forward-fill non-NA values into subsequent daytime hours.
  # Nighttime NAs remain NA (residual wind is only meaningful overnight).
  .freeze_fwd <- function(dir_vec) {
    out  <- dir_vec
    last <- NA_real_
    for (i in seq_along(dir_vec)) {
      if (!is.na(dir_vec[i])) {
        last <- dir_vec[i]
      } else if (vs$is_day[i]) {
        out[i] <- last
      }
    }
    out
  }

  dir_vecs <- list(
    `80`  = if (!is.null(rw$dir_80m))  .freeze_fwd(rw$dir_80m)  else rep(NA_real_, n_t),
    `120` = if (!is.null(rw$dir_120m)) .freeze_fwd(rw$dir_120m) else rep(NA_real_, n_t),
    `180` = if (!is.null(rw$dir_180m)) .freeze_fwd(rw$dir_180m) else rep(NA_real_, n_t)
  )
  dir_10m_frz <- if (!is.null(rw$dir_10m)) .freeze_fwd(rw$dir_10m) else rep(NA_real_, n_t)

  level_heights <- c(80, 120, 180)
  lev_order     <- order(abs(level_heights - above_pool_ht))
  level_names   <- c("80", "120", "180")[lev_order]

  list(dir_vecs = dir_vecs, dir_10m_frz = dir_10m_frz,
       wind_dir_surf = wind_dir_surf, level_names = level_names)
}


# Matrix form: `bearings` is a length-n_r vector of receptor bearings; returns
# an (n_t x n_r) matrix. The residual-wind direction (and hence the per-hour
# downwind bearing) is receptor-independent, so the only outer combination is
# the cosine of each receptor bearing against the per-hour downwind direction.
.cw_fumigation_matrix <- function(bearings, vs, terrain,
                                  above_pool_ht = 0,
                                  prep     = NULL,
                                  FUMIC_1B = ODOUR_CONSTANTS$FUMIC_1B) {
  n_t <- length(vs$is_day)
  n_r <- length(bearings)

  if (is.null(prep)) prep <- .cw_fumigation_prep(vs, above_pool_ht)

  # Active hours: morning transition with a pool present.
  is_active <- vs$is_day &
               !is.na(vs$cbl_growth) & vs$cbl_growth > 0 &
               !is.na(vs$pool_top)   & vs$pool_top   > 0

  # Coalesce direction across levels (priority order) then fallbacks (len n_t).
  rw_dir_vec <- rep(NA_real_, n_t)
  for (lev_nm in prep$level_names) {
    fill <- is.na(rw_dir_vec)
    rw_dir_vec[fill] <- prep$dir_vecs[[lev_nm]][fill]
  }
  fill <- is.na(rw_dir_vec); rw_dir_vec[fill] <- prep$dir_10m_frz[fill]
  fill <- is.na(rw_dir_vec); rw_dir_vec[fill] <- prep$wind_dir_surf[fill]

  # Outer cosine: downwind per hour (down columns) vs bearing per receptor.
  downwind_dir <- (rw_dir_vec + 180) %% 360
  diff_deg <- .angular_diff(matrix(bearings,     n_t, n_r, byrow = TRUE),
                            matrix(downwind_dir, n_t, n_r))
  cw_val   <- FUMIC_1B * pmax(0, cos(diff_deg * pi / 180))^2
  dim(cw_val) <- c(n_t, n_r)   # pmax() drops the matrix dim; restore it

  cw         <- matrix(NA_real_, n_t, n_r)
  known      <- !is.na(rw_dir_vec)               # len n_t
  rows_known <- is_active & known
  rows_zero  <- is_active & !known
  if (any(rows_known)) cw[rows_known, ] <- cw_val[rows_known, ]
  if (any(rows_zero))  cw[rows_zero, ]  <- 0
  cw
}

# Scalar form (single receptor bearing): thin wrapper over the matrix form.
.cw_fumigation <- function(bearing_to_receptor, vs, terrain,
                           above_pool_ht = 0,
                           prep     = NULL,
                           FUMIC_1B = ODOUR_CONSTANTS$FUMIC_1B) {
  .cw_fumigation_matrix(bearing_to_receptor, vs, terrain,
                        above_pool_ht = above_pool_ht, prep = prep,
                        FUMIC_1B = FUMIC_1B)[, 1L]
}


# ---------------------------------------------------------------------------
# .morning_release(pool_top, cbl_growth, is_day)
# ---------------------------------------------------------------------------
# Computes a per-hour release enhancement factor for the morning pulse.
#
# At each morning-onset hour t0 (first daytime hour after a night with a pool):
#   A   = pool_top[t0]  (accumulated mass proxy, m)
#   tau = A / max(cbl_growth[t0], 0.1)  (characteristic release time, hours)
#   r_raw(t) = (A / tau) * exp(-0.5 * ((t - t0) / tau)^2)
#   r is normalised so sum(r) = A  (mass conservation)
#
# Returns a numeric vector of length n_t. Each element r[t] is the release
# enhancement: the total 1b enhancement at that hour is proportional to r[t].
# In the flat baseline (no morning pulse) r[t] = 0 everywhere.
#
# Invariant: sum(r[t]) over the release window ≈ A (mass conservation).

.morning_release <- function(pool_top, cbl_growth, is_day) {
  n_t     <- length(pool_top)
  release <- rep(0.0, n_t)   # 0 = no enhancement (additive above baseline 1b)

  for (t in seq_len(n_t)) {
    # Morning onset: first daytime hour that immediately follows a night hour
    if (!is_day[t])    next
    if (t == 1L)       next
    if (is_day[t - 1L]) next   # previous hour was already day → not a transition

    if (is.na(pool_top[t])   || pool_top[t]   <= 0) next
    if (is.na(cbl_growth[t]) || cbl_growth[t] <= 0) next

    t0  <- t
    A   <- pool_top[t0]
    tau <- A / max(cbl_growth[t0], 0.1)

    # Gaussian release kernel over the morning window
    window_len <- min(n_t - t0 + 1L, max(1L, as.integer(ceiling(4 * tau + 1))))
    t_range    <- t0:(t0 + window_len - 1L)
    r_raw      <- (A / tau) * exp(-0.5 * (((t_range - t0)) / max(tau, 0.01))^2)

    # Normalise so sum = A (mass conservation)
    total <- sum(r_raw)
    if (total > 0) {
      r_norm <- r_raw * A / total
    } else {
      r_norm <- r_raw
    }

    for (i in seq_along(t_range)) {
      ti <- t_range[i]
      if (ti <= n_t) {
        release[ti] <- release[ti] + r_norm[i]
      }
    }
  }

  release
}
