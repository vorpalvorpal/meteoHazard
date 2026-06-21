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

  f_1b <- numeric(n)

  for (i in seq_len(n)) {
    e  <- emit_extent[i]
    pt <- pool_top[i]
    d  <- max(delta[i], 1e-6)

    if (is.na(e) || is.na(pt) || is.na(d)) {
      f_1b[i] <- 0.0
      next
    }

    if (e <= 0) {
      # Point source at z = 0
      f_1b[i] <- phi(-pt / d)
    } else {
      # Uniform E(z) over [0, e]: mean of Phi over that range
      z_seq   <- seq(0, e, length.out = 50)
      f_1b[i] <- mean(phi((z_seq - pt) / d))
    }
  }

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

.cw_venting <- function(bearing_to_receptor, vs, terrain,
                        CONFINEMENT_1A = ODOUR_CONSTANTS$CONFINEMENT_1A,
                        VENTING_1A     = ODOUR_CONSTANTS$VENTING_1A) {
  n_t   <- length(vs$is_day)
  cw_1a <- rep(NA_real_, n_t)

  drain_bearing <- if (!is.null(terrain) && !is.na(terrain@drainage_bearing))
    terrain@drainage_bearing else NA_real_
  flow_conv <- if (!is.null(terrain) && !is.na(terrain@flow_convergence))
    terrain@flow_convergence else 0

  is_channelled <- !is.na(drain_bearing) && flow_conv >= 0.5

  for (t in seq_len(n_t)) {
    if (!vs$is_day[t]) {
      # Night: confinement
      if (is_channelled) {
        diff_deg  <- .angular_diff(bearing_to_receptor, drain_bearing)
        alignment <- max(0, cos(diff_deg * pi / 180))
        cw_1a[t] <- CONFINEMENT_1A * (1 - 0.9 * alignment)
      } else {
        cw_1a[t] <- CONFINEMENT_1A
      }
    } else if (!is.na(vs$cbl_growth[t]) && vs$cbl_growth[t] > 0) {
      # Morning transition: anabatic venting
      if (is_channelled) {
        diff_deg  <- .angular_diff(bearing_to_receptor, drain_bearing)
        alignment <- max(0, cos(diff_deg * pi / 180))
        cw_1a[t] <- VENTING_1A * alignment
      } else {
        cw_1a[t] <- VENTING_1A
      }
    }
    # else daytime without cbl_growth > 0: NA (flat Gaussian used by caller)
  }

  cw_1a
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

.cw_fumigation <- function(bearing_to_receptor, vs, terrain,
                           above_pool_ht = 0,    # above-pool emission height (m)
                           FUMIC_1B = ODOUR_CONSTANTS$FUMIC_1B) {
  n_t   <- length(vs$is_day)
  cw_1b <- rep(NA_real_, n_t)

  rw <- vs$residual_wind

  # Surface wind direction (10 m), used as ultimate fallback if available in vs.
  # Callers can set vs$wind_dir_surface (e.g. odour_exposure() passes met wind_dir).
  wind_dir_surf <- if (!is.null(vs$wind_dir_surface)) vs$wind_dir_surface
                   else rep(NA_real_, n_t)

  # Pre-freeze residual-wind directions into daytime (carry last nighttime
  # direction forward so morning fumigation can access the overnight wind).
  .freeze_fwd <- function(dir_vec) {
    out  <- dir_vec
    last <- NA_real_
    for (i in seq_along(dir_vec)) {
      if (!is.na(dir_vec[i])) {
        last <- dir_vec[i]
      } else if (vs$is_day[i]) {
        out[i] <- last   # carry overnight value into morning
      }
    }
    out
  }

  # Build frozen direction vectors for each available level.
  dir_vecs <- list(
    `80`  = if (!is.null(rw$dir_80m))  .freeze_fwd(rw$dir_80m)  else rep(NA_real_, n_t),
    `120` = if (!is.null(rw$dir_120m)) .freeze_fwd(rw$dir_120m) else rep(NA_real_, n_t),
    `180` = if (!is.null(rw$dir_180m)) .freeze_fwd(rw$dir_180m) else rep(NA_real_, n_t)
  )
  dir_10m_frz <- if (!is.null(rw$dir_10m)) .freeze_fwd(rw$dir_10m) else rep(NA_real_, n_t)

  # Select the level order by closeness of level height to above_pool_ht.
  # Levels available: 80 m, 120 m, 180 m.
  level_heights <- c(80, 120, 180)
  lev_order     <- order(abs(level_heights - above_pool_ht))
  level_names   <- c("80", "120", "180")[lev_order]

  for (t in seq_len(n_t)) {
    if (!vs$is_day[t]) next
    if (is.na(vs$cbl_growth[t]) || vs$cbl_growth[t] <= 0) next
    if (is.na(vs$pool_top[t])   || vs$pool_top[t]   <= 0) next

    # Select the direction at the level closest to above_pool_ht; fall back
    # through remaining levels, then 10m frozen, then instantaneous surface.
    rw_dir <- NA_real_
    for (lev_nm in level_names) {
      dv <- dir_vecs[[lev_nm]]
      if (!is.na(dv[t])) { rw_dir <- dv[t]; break }
    }
    if (is.na(rw_dir) && !is.na(dir_10m_frz[t])) rw_dir <- dir_10m_frz[t]
    if (is.na(rw_dir) && !is.na(wind_dir_surf[t])) rw_dir <- wind_dir_surf[t]

    if (is.na(rw_dir)) {
      # Truly no directional information
      cw_1b[t] <- 0
      next
    }

    # Fumigation goes downwind (rw_dir is blows-FROM convention)
    downwind_dir <- (rw_dir + 180) %% 360
    diff_deg     <- .angular_diff(bearing_to_receptor, downwind_dir)
    cw_1b[t]    <- FUMIC_1B * max(0, cos(diff_deg * pi / 180))^2
  }

  cw_1b
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
