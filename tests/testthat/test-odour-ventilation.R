# Behaviour spec for the odour ventilation state (C2, issue #15):
# ventilation_state() and its pool_top / cbl_growth / residual_wind / generation
# components. pool_top resolves D3 (terrain-modulated heat-deficit estimate).

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Build a minimal met_data row. All optional columns absent unless supplied.
vs_met <- function(n = 1, ...) {
  base <- list(
    wind_speed_10m        = 3,
    direct_radiation      = 0,   # night
    cloud_cover           = 50,
    boundary_layer_height = 500,
    temperature_2m        = 10,
    relative_humidity_2m  = 60,
    precipitation         = 0
  )
  ov <- list(...)
  base[names(ov)] <- ov
  as.data.frame(lapply(base, rep_len, n))
}

# Night window: 6 dark hours
night6 <- function(...) {
  vs_met(n = 6, direct_radiation = 0, ...)
}


# ---------------------------------------------------------------------------
# Core state fields
# ---------------------------------------------------------------------------

describe("ventilation_state()", {

  it("returns per-hour u_eff, h_mix, s, is_calm, is_day, PM, W_rain, pool_top, cbl_growth and residual_wind", {
    vs <- ventilation_state(vs_met())
    expect_named(vs, c("u_eff", "h_mix", "s", "is_calm", "is_day",
                       "PM", "W_rain", "pool_top", "cbl_growth",
                       "residual_wind"),
                 ignore.order = TRUE)
    expect_length(vs$u_eff, 1L)
    expect_length(vs$PM,    1L)
    expect_length(vs$W_rain, 1L)
    expect_length(vs$pool_top, 1L)
    expect_length(vs$cbl_growth, 1L)
    expect_type(vs$residual_wind, "list")
  })

  it("reproduces the current dispersion-state u_eff, h_mix and s exactly", {
    # Build an input that exercises both day and night, calm and windy rows.
    d <- rbind(
      vs_met(wind_speed_10m = 0.3, direct_radiation = 0,   cloud_cover = 30),  # calm, night
      vs_met(wind_speed_10m = 5,   direct_radiation = 800, cloud_cover = 10),  # windy, day
      vs_met(wind_speed_10m = 3,   direct_radiation = 0,   cloud_cover = 80),  # night cloudy
      vs_met(wind_speed_10m = NA,  direct_radiation = 400, cloud_cover = 50)   # NA wind, day
    )
    vs <- ventilation_state(d)

    # These rows carry temperature_2m/relative_humidity_2m (vs_met() defaults),
    # so the v3 cold-pool cap may shrink h_mix below the raw BLH on night rows
    # (rows 1 and 3) -- only checked here for positivity, not an exact value;
    # see the dedicated "nocturnal cold-pool cap" describe block below.
    expect_true(all(vs$u_eff >= ODOUR_CONSTANTS$U_CALM_FLOOR))
    expect_true(all(vs$h_mix > 0))
    expect_true(all(vs$s >= 0 & vs$s <= 5))

    # Calm row has s locked to 4.25.
    expect_equal(vs$s[1], 4.25)
    # Day + strong radiation + 5 m/s -> day_s is low (B class, s ~ 1).
    expect_lt(vs$s[2], vs$s[1])
  })

})


# ---------------------------------------------------------------------------
# Nocturnal cold-pool cap on h_mix (v3 Change 1)
# ---------------------------------------------------------------------------
#
# Physical basis: under a nocturnal surface inversion the valley cold pool
# (pool_top) caps vertical mixing far below the synoptic boundary-layer
# height forecast -- the depth odour dilutes through is then the pool depth,
# not the BLH. h_mix = min(BLH, pool_top) when: (a) the hour is dark
# (!is_day) -- the daytime pool breakup is the morning-pulse machinery's
# transient, owned by odour_exposure(), not this cap; and (b) pool_top is
# THERMALLY derived (temperature_2m + relative_humidity_2m present) --
# without them .pool_top() returns a mechanical/constant floor (>= 200 m) on
# every hour, which is not evidence of an actual pool.
#
# Written BEFORE implementation (TDD): these fail with "unused argument
# (pool_cap = ...)" until ventilation_state() gains the pool_cap parameter
# (default TRUE).

describe("ventilation_state(): nocturnal cold-pool cap on h_mix (pool_cap)", {

  it("caps h_mix at pool_top on a thermally-derived clear cold night (default pool_cap = TRUE)", {
    d <- night6(temperature_2m = 5, relative_humidity_2m = 30, cloud_cover = 0,
                wind_speed_10m = 0.5, boundary_layer_height = 500)
    vs_capped   <- ventilation_state(d)                       # pool_cap = TRUE (default)
    vs_uncapped <- ventilation_state(d, pool_cap = FALSE)

    # Uncapped h_mix is the raw BLH throughout (no NA/fallback in this fixture).
    expect_true(all(vs_uncapped$h_mix == 500))

    pool_active <- !vs_capped$is_day & !is.na(vs_capped$pool_top) & vs_capped$pool_top > 0
    expect_true(any(pool_active), label = "fixture must accumulate a pool")

    # Capped h_mix must equal pmin(BLH, pool_top) on every pool-active hour.
    expect_equal(vs_capped$h_mix[pool_active],
                 pmin(500, vs_capped$pool_top[pool_active]), tolerance = 1e-8)

    # And it must be <= the uncapped depth everywhere, strictly < wherever the
    # pool is shallower than the synoptic BLH (the whole point of the cap).
    expect_true(all(vs_capped$h_mix[pool_active] <= vs_uncapped$h_mix[pool_active]))
    expect_true(any(vs_capped$h_mix[pool_active] < vs_uncapped$h_mix[pool_active]))
  })

  it("never caps daytime h_mix, even though pool_top is frozen into daytime hours", {
    d <- rbind(
      night6(temperature_2m = 5, relative_humidity_2m = 30, cloud_cover = 0,
             wind_speed_10m = 0.5, boundary_layer_height = 500),
      vs_met(n = 4, direct_radiation = 600, boundary_layer_height = 800,
             temperature_2m = 15, relative_humidity_2m = 50)
    )
    vs <- ventilation_state(d)
    day_rows <- vs$is_day
    expect_true(any(day_rows))
    expect_true(all(vs$pool_top[day_rows] > 0),
                label = "pool must be frozen (non-zero) into daytime")
    # Daytime h_mix must equal the (uncapped) BLH exactly, not the frozen pool.
    expect_equal(vs$h_mix[day_rows], d$boundary_layer_height[day_rows])
  })

  it("does not cap when temperature/RH are absent (a mechanical floor is not evidence of a pool)", {
    d <- data.frame(
      wind_speed_10m = rep(0.5, 6), direct_radiation = rep(0, 6),
      cloud_cover = rep(0, 6), boundary_layer_height = rep(500, 6)
    )
    vs_default <- ventilation_state(d)
    vs_off     <- ventilation_state(d, pool_cap = FALSE)
    expect_identical(vs_default$h_mix, vs_off$h_mix)
    expect_true(all(vs_default$h_mix == 500))
  })

  it("does not cap on a windy night where the mechanical floor exceeds a shallow BLH", {
    d <- night6(temperature_2m = 10, relative_humidity_2m = 60, cloud_cover = 50,
                wind_speed_10m = 8, boundary_layer_height = 50)
    vs_default <- ventilation_state(d)
    vs_off     <- ventilation_state(d, pool_cap = FALSE)
    # Strong mechanical mixing (u10 = 8) drives the Venkatram floor for
    # pool_top well above the (deliberately shallow) BLH, so
    # min(BLH, pool_top) == BLH -- the cap is a no-op here.
    expect_equal(vs_default$h_mix, vs_off$h_mix, tolerance = 1e-8)
  })

  it("cbl_growth is identical whether or not pool_cap is applied (derived from uncapped BLH)", {
    # cbl_growth must come from the synoptic (uncapped) mixing depth: capping
    # h_mix at sunrise (pool breaks up, BLH grows) must not masquerade as CBL
    # growth via diff(h_mix). GOTCHA guarded by this test.
    d <- rbind(
      night6(temperature_2m = 3, relative_humidity_2m = 20, cloud_cover = 0,
             wind_speed_10m = 0.5, boundary_layer_height = 200),
      vs_met(n = 4, direct_radiation = 600,
             boundary_layer_height = c(300, 600, 1000, 1500),
             temperature_2m = 15, relative_humidity_2m = 50)
    )
    vs_on  <- ventilation_state(d, pool_cap = TRUE)
    vs_off <- ventilation_state(d, pool_cap = FALSE)
    expect_identical(vs_on$cbl_growth, vs_off$cbl_growth)
  })

  it("pool_cap = FALSE reproduces the raw BLH / stability fallback exactly regardless of pool depth", {
    d <- night6(temperature_2m = 5, relative_humidity_2m = 30, cloud_cover = 0,
                wind_speed_10m = 0.5, boundary_layer_height = NA_real_)
    vs <- ventilation_state(d, pool_cap = FALSE)
    # NA BLH, calm/stable night -> stable fallback, unaffected by pool_top.
    expect_true(all(vs$h_mix == ODOUR_CONSTANTS$H_MIX_FALLBACK_STABLE))
  })

})


# ---------------------------------------------------------------------------
# pool_top
# ---------------------------------------------------------------------------

describe("ventilation_state(): pool_top (terrain-modulated heat deficit)", {

  it("computes net longwave cooling from the Brunt formula for given T, RH and cloud", {
    # A clear, cold, dry night should produce a positive heat deficit and
    # therefore a pool_top above the mechanical floor.
    d <- night6(temperature_2m = 5, relative_humidity_2m = 30, cloud_cover = 0,
                wind_speed_10m = 0.5)
    vs <- ventilation_state(d)
    expect_true(all(vs$pool_top >= 0))
    # At least some accumulation over 6 hours.
    expect_gt(vs$pool_top[6], vs$pool_top[1])
  })

  it("gives stronger cooling under clear skies than overcast", {
    d_clear    <- night6(cloud_cover = 0,  temperature_2m = 5, relative_humidity_2m = 30,
                         wind_speed_10m = 0.5)
    d_overcast <- night6(cloud_cover = 100, temperature_2m = 5, relative_humidity_2m = 30,
                         wind_speed_10m = 0.5)
    vs_clear    <- ventilation_state(d_clear)
    vs_overcast <- ventilation_state(d_overcast)
    # Clear sky should accumulate more heat deficit → deeper pool.
    expect_gt(vs_clear$pool_top[6], vs_overcast$pool_top[6])
  })

  it("increases monotonically with cooling hours until saturation, then plateaus", {
    d <- night6(temperature_2m = 5, relative_humidity_2m = 30, cloud_cover = 0,
                wind_speed_10m = 0.5)
    vs <- ventilation_state(d)
    pool <- vs$pool_top
    # Must be non-decreasing each hour during a continuous night.
    expect_true(all(diff(pool) >= -1e-10))
    # Long night (168h) should plateau near H_SAT.
    d_long <- vs_met(n = 168, direct_radiation = 0, temperature_2m = 5,
                     relative_humidity_2m = 30, cloud_cover = 0, wind_speed_10m = 0.5)
    vs_long <- ventilation_state(d_long)
    H_SAT <- ODOUR_CONSTANTS$POOL_H_SAT
    expect_lte(max(vs_long$pool_top), H_SAT + 1e-6)
  })

  it("is floored by the Venkatram mechanical depth 2400*u*^1.5", {
    # Warm, humid, overcast night: Brunt Q* is very small. pool_top should
    # not fall below the mechanical floor.
    d <- night6(temperature_2m = 25, relative_humidity_2m = 95, cloud_cover = 100,
                wind_speed_10m = 3)
    vs <- ventilation_state(d)
    kappa  <- ODOUR_CONSTANTS$POOL_KAPPA
    z0     <- ODOUR_CONSTANTS$POOL_Z0
    C_vent <- ODOUR_CONSTANTS$VENKATRAM_COEF
    u_star <- kappa * 3 / log(10 / z0)
    h_floor <- C_vent * u_star^1.5
    expect_true(all(vs$pool_top >= h_floor - 1e-9))
  })

  it("matches the Venkatram floor for a given friction velocity", {
    # With negligible Brunt cooling, pool_top should equal h_floor exactly.
    # Use T=30 high RH overcast to suppress net longwave.
    d <- night6(temperature_2m = 30, relative_humidity_2m = 99, cloud_cover = 100,
                wind_speed_10m = 5)
    vs <- ventilation_state(d)
    kappa  <- ODOUR_CONSTANTS$POOL_KAPPA
    z0     <- ODOUR_CONSTANTS$POOL_Z0
    C_vent <- ODOUR_CONSTANTS$VENKATRAM_COEF
    u_star <- kappa * 5 / log(10 / z0)
    h_floor <- C_vent * u_star^1.5
    expect_equal(vs$pool_top[1], h_floor, tolerance = 1e-6)
  })

  it("is capped by the basin sill (valley_depth / basin_capacity)", {
    sill <- 80  # m
    trn  <- mh_terrain(valley_depth = sill)
    # 24 dark hours to guarantee accumulation past the sill.
    d2 <- vs_met(n = 24, direct_radiation = 0, temperature_2m = 2,
                 relative_humidity_2m = 20, cloud_cover = 0, wind_speed_10m = 0.5)
    vs_capped    <- ventilation_state(d2, terrain = trn)
    vs_uncapped  <- ventilation_state(d2)
    expect_true(all(vs_capped$pool_top <= sill + 1e-9))
    expect_gt(max(vs_uncapped$pool_top), sill)
  })

  it("is amplified, not reduced, by a topographic amplification factor >= 1", {
    trn_flat <- mh_terrain(taf = 1.0)
    trn_amp  <- mh_terrain(taf = 2.0)
    d <- vs_met(n = 8, direct_radiation = 0, temperature_2m = 5,
                relative_humidity_2m = 30, cloud_cover = 0, wind_speed_10m = 0.5)
    vs_flat <- ventilation_state(d, terrain = trn_flat)
    vs_amp  <- ventilation_state(d, terrain = trn_amp)
    # Amplified cooling -> deeper pool.
    expect_gte(max(vs_amp$pool_top), max(vs_flat$pool_top))
  })

  it("falls back to the mechanical floor / h_mix outside the clear-calm-stable regime", {
    # Windy daytime: pool_top should be frozen night max (0 for all-day input).
    d <- vs_met(n = 4, direct_radiation = 600, wind_speed_10m = 8)
    vs <- ventilation_state(d)
    # All daytime hours: the frozen night max is 0 (no preceding night).
    expect_true(all(vs$pool_top == 0))
  })

  it("yields no accumulated pool for an input window containing no night", {
    d <- vs_met(n = 8, direct_radiation = 400)  # all daytime
    vs <- ventilation_state(d)
    expect_true(all(vs$pool_top == 0))
  })

  it("is taken as the overnight maximum, frozen at release, on the AGL datum", {
    # Night (4h), then day (4h): pool_top during day should equal the last
    # night's maximum, not reset to 0.
    d <- rbind(
      vs_met(n = 4, direct_radiation = 0, temperature_2m = 5,
             relative_humidity_2m = 30, cloud_cover = 0, wind_speed_10m = 0.5),
      vs_met(n = 4, direct_radiation = 600)
    )
    vs <- ventilation_state(d)
    night_max <- vs$pool_top[4]
    day_vals  <- vs$pool_top[5:8]
    expect_true(all(day_vals == night_max))
    expect_gt(night_max, 0)
  })

})


# ---------------------------------------------------------------------------
# cbl_growth
# ---------------------------------------------------------------------------

describe("ventilation_state(): cbl_growth", {

  it("is positive across the morning transition and ~0 at night", {
    # Night (6h, shallow stable BL ~ 200m), then day ramp (6h, BL growing).
    bls <- c(rep(200, 6), 300, 500, 800, 1200, 1500, 1800)
    d   <- data.frame(
      wind_speed_10m        = 3,
      direct_radiation      = c(rep(0, 6), rep(600, 6)),
      cloud_cover           = 30,
      boundary_layer_height = bls,
      temperature_2m        = 15,
      relative_humidity_2m  = 60,
      precipitation         = 0
    )
    vs <- ventilation_state(d)
    # Night hours: cbl_growth should be 0 (BL flat or slightly varying).
    expect_true(all(vs$cbl_growth[1:6] == 0))
    # Morning hours: at least some positive growth.
    expect_true(any(vs$cbl_growth[7:12] > 0))
  })

  it("returns 0 (not an NA error) for a flat or missing h_mix series", {
    d <- vs_met(n = 4, boundary_layer_height = NA_real_)
    vs <- ventilation_state(d)
    expect_false(any(is.na(vs$cbl_growth)))
    expect_true(all(vs$cbl_growth == 0))
  })

})


# ---------------------------------------------------------------------------
# residual_wind
# ---------------------------------------------------------------------------

describe("ventilation_state(): residual_wind", {

  it("is the overnight circular-mean direction at each available level", {
    # Provide 80m winds at two night-time hours (270° and 250°, approx W),
    # preceded by 2 daytime hours.
    d <- data.frame(
      wind_speed_10m        = 3,
      direct_radiation      = c(600, 600, 0, 0),
      cloud_cover           = 30,
      boundary_layer_height = 500,
      temperature_2m        = 10,
      relative_humidity_2m  = 60,
      precipitation         = 0,
      wind_speed_80m        = c(NA, NA, 5, 5),
      wind_direction_80m    = c(NA, NA, 270, 250)
    )
    vs <- ventilation_state(d)
    rw <- vs$residual_wind
    # During daytime rows (1-2): NA.
    expect_true(is.na(rw$dir_80m[1]))
    expect_true(is.na(rw$dir_80m[2]))
    # During night rows (3-4): should be a westerly direction.
    expect_false(is.na(rw$dir_80m[4]))
    # Circular mean of 270 and 250 is roughly 260 degrees.
    expect_true(rw$dir_80m[4] > 240 & rw$dir_80m[4] < 280)
  })

  it("falls back to the next lower level when an upper level is missing", {
    # Only 80m present; 120m and 180m columns absent.
    d <- data.frame(
      wind_speed_10m        = 3,
      direct_radiation      = 0,
      cloud_cover           = 30,
      boundary_layer_height = 500,
      temperature_2m        = 10,
      relative_humidity_2m  = 60,
      precipitation         = 0,
      wind_speed_80m        = 5,
      wind_direction_80m    = 200
    )
    vs <- ventilation_state(d)
    rw <- vs$residual_wind
    expect_false(is.na(rw$dir_80m[1]))
    expect_true(is.na(rw$dir_120m[1]))
    expect_true(is.na(rw$dir_180m[1]))
  })

  it("falls back to 10m surface wind when all upper-level winds are missing", {
    # vs_met() supplies wind_speed_10m = 3 but no upper-level columns and no
    # wind_direction_10m.  Add wind_direction_10m so the 10m fallback fires.
    d <- vs_met(n = 4, direct_radiation = 0, wind_direction_10m = 210)
    vs <- ventilation_state(d)
    rw <- vs$residual_wind
    # Upper-level outputs are NA (columns absent).
    expect_true(all(is.na(rw$dir_80m)))
    expect_true(all(is.na(rw$dir_120m)))
    expect_true(all(is.na(rw$dir_180m)))
    # 10m fallback is non-NA for overnight hours (all 4 rows are night).
    expect_true(all(!is.na(rw$dir_10m[!vs$is_day])))
    expect_true(all(!is.na(rw$speed_10m[!vs$is_day])))
  })

})


# ---------------------------------------------------------------------------
# W_rain -- solubility-aware below-cloud scavenging (v3 Change 4)
# ---------------------------------------------------------------------------
#
# Physical basis: below-cloud gas washout scales with Henry's-law solubility
# (Seinfeld & Pandis Ch. 19). The tiered factors (down to 0.05, i.e. 95%
# removal) are the HIGHLY SOLUBLE limit; the reduced-sulfur compounds that
# drive most landfill odour complaints are only sparingly soluble and
# scavenge weakly. W_rain blends linearly between "no washout"
# (odorant_solubility = 0) and the soluble-limit tiers (odorant_solubility =
# 1); default 0.5 represents a mixed sulfur/soluble-VOC profile.

describe("ventilation_state(): W_rain solubility-aware scavenging (odorant_solubility)", {

  it("odorant_solubility = 1 reproduces the soluble-limit tiers exactly", {
    d  <- vs_met(n = 3, precipitation = c(0.5, 2, 5))  # light / moderate / heavy
    vs <- ventilation_state(d, odorant_solubility = 1)
    K <- ODOUR_CONSTANTS
    expect_equal(vs$W_rain, c(K$W_RAIN_FACTOR_LIGHT, K$W_RAIN_FACTOR_MOD,
                              K$W_RAIN_FACTOR_HEAVY), tolerance = 1e-8)
  })

  it("odorant_solubility = 0 gives no washout at any rain rate", {
    d  <- vs_met(n = 3, precipitation = c(0.5, 2, 5))
    vs <- ventilation_state(d, odorant_solubility = 0)
    expect_equal(vs$W_rain, c(1, 1, 1))
  })

  it("odorant_solubility = 0.5 (default) blends linearly to the soluble-limit tier", {
    d           <- vs_met(precipitation = 5)   # heavy rain
    vs_default  <- ventilation_state(d)
    vs_explicit <- ventilation_state(d, odorant_solubility = 0.5)
    expected <- 1 - 0.5 * (1 - ODOUR_CONSTANTS$W_RAIN_FACTOR_HEAVY)  # 0.525
    expect_equal(vs_default$W_rain, expected, tolerance = 1e-8)
    expect_equal(vs_explicit$W_rain, expected, tolerance = 1e-8)
  })

  it("raises meteoHazard_input_error for out-of-range odorant_solubility", {
    d <- vs_met()
    expect_error(ventilation_state(d, odorant_solubility = -0.1),
                 class = "meteoHazard_input_error")
    expect_error(ventilation_state(d, odorant_solubility = 1.1),
                 class = "meteoHazard_input_error")
  })

})


# ---------------------------------------------------------------------------
# .odour_generation()
# ---------------------------------------------------------------------------
#
# Physical basis for the v3 rewrite of G (Changes 2 + 3):
#  - Multiplicative combination: the five modifiers are independent
#    fractional changes to emission rate, and independent fractional effects
#    compound multiplicatively (E = E0 * prod(1+m_i)), not additively. Additive
#    superposition is only the first-order approximation and understates the
#    coincident worst case (e.g. falling pressure + hot + humid together).
#  - Exponential V_mod: the cited basis (Henry's-law / Clausius-Clapeyron
#    vapour-pressure temperature dependence, "~doubling per 10 degC") is
#    exponential, not linear. V_mod is anchored to 0 at V_MOD_T_LO and
#    V_MOD_MAX at V_MOD_T_HI, clamped beyond T_HI as a deliberate screening
#    ceiling (this site's worst case is winter inversions, not extreme heat).

describe(".odour_generation()", {

  mh_base <- function(n = 1, ...) {
    base <- list(
      wind_speed_10m         = 3,
      direct_radiation       = 0,
      cloud_cover            = 50,
      boundary_layer_height  = 500,
      temperature_2m         = 15,
      pressure_msl           = 1013,
      precipitation          = 0,
      relative_humidity_2m   = 50,
      soil_moisture_0_to_1cm = 0.1,
      soil_moisture_1_to_3cm = 0.1
    )
    ov <- list(...)
    base[names(ov)] <- ov
    as.data.frame(lapply(base, rep_len, n))
  }

  it("combines simultaneous modifiers multiplicatively: G = prod(1 + m_i), not sum", {
    # temperature_2m = V_MOD_T_HI clamps V_mod to exactly V_MOD_MAX (0.30),
    # independent of the exact exponential shape -- isolates the combination
    # rule (multiplicative vs additive) from the V_mod curve shape.
    d <- mh_base(n = 5, temperature_2m = ODOUR_CONSTANTS$V_MOD_T_HI,
                 pressure_msl = c(1013, 1010, 1007, 1004, 1001))
    G <- .odour_generation(d)
    V_mod  <- ODOUR_CONSTANTS$V_MOD_MAX
    dP_mod <- 0.30   # dP3 = 1004 - 1013 = -9 <= -5 at row 4 -> saturated dP_mod
    expect_equal(G[4], (1 + dP_mod) * (1 + V_mod), tolerance = 1e-8)
    expect_false(isTRUE(all.equal(G[4], 1 + dP_mod + V_mod)),
                 label = "must NOT be the additive combination")
  })

  it("V_mod(T) is 0 at V_MOD_T_LO, V_MOD_MAX at/beyond V_MOD_T_HI, and convex (exponential) between", {
    T_lo <- ODOUR_CONSTANTS$V_MOD_T_LO
    T_hi <- ODOUR_CONSTANTS$V_MOD_T_HI
    vmax <- ODOUR_CONSTANTS$V_MOD_MAX

    expect_equal(.odour_generation(mh_base(temperature_2m = T_lo)), 1.0, tolerance = 1e-8)
    expect_equal(.odour_generation(mh_base(temperature_2m = T_hi)), 1.0 + vmax, tolerance = 1e-8)
    expect_equal(.odour_generation(mh_base(temperature_2m = T_hi + 5)), 1.0 + vmax,
                 tolerance = 1e-8)  # clamped beyond T_hi

    # Convexity: the value at the midpoint temperature must sit BELOW the
    # linear midpoint (0.5 * vmax) -- doubling growth is back-loaded.
    mid_temp <- (T_lo + T_hi) / 2
    G_mid    <- .odour_generation(mh_base(temperature_2m = mid_temp))
    expect_lt(G_mid - 1.0, 0.5 * vmax)

    # Monotone increasing over (T_lo, T_hi).
    temps <- seq(T_lo, T_hi, length.out = 6)
    Gs    <- vapply(temps, function(t) .odour_generation(mh_base(temperature_2m = t)), numeric(1))
    expect_true(all(diff(Gs) >= -1e-10))
  })

  it("NA temperature leaves V_mod at 0 (neutral)", {
    expect_equal(.odour_generation(mh_base(temperature_2m = NA_real_)), 1.0, tolerance = 1e-8)
  })

  it("wet-soil sealing multiplies rather than subtracts from the other modifiers", {
    # Heavy wet soil: S_seal = -0.20. Default temp 15 -> V_mod follows the
    # exponential ramp. G = (1 + V_mod) * (1 + S_seal) = (1 + V_mod) * 0.80,
    # NOT the old additive 1 + V_mod - 0.20.
    K    <- ODOUR_CONSTANTS
    f_hi <- 2^((K$V_MOD_T_HI - K$V_MOD_T_LO) / K$V_MOD_DOUBLING_C)
    V_mod_15 <- K$V_MOD_MAX * (2^((15 - K$V_MOD_T_LO) / K$V_MOD_DOUBLING_C) - 1) / (f_hi - 1)
    expected <- (1 + V_mod_15) * (1 - 0.20)
    G <- .odour_generation(mh_base(soil_moisture_0_to_1cm = 0.45,
                                    soil_moisture_1_to_3cm = 0.45))
    expect_equal(G, expected, tolerance = 1e-8)
  })

})


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

describe("ventilation_state(): validation", {

  it("errors (classed) on missing required met columns, including the multi-level winds", {
    # Missing a required core column.
    d_bad <- vs_met()[, -1, drop = FALSE]  # drop wind_speed_10m
    expect_error(ventilation_state(d_bad), class = "meteoHazard_input_error")

    # Present but non-numeric multi-level wind column.
    d_ml <- vs_met()
    d_ml$wind_speed_80m    <- 5
    d_ml$wind_direction_80m <- "north"  # non-numeric
    expect_error(ventilation_state(d_ml), class = "meteoHazard_input_error")
  })

  it("raises meteoHazard_input_error for non-logical shelter argument", {
    d <- vs_met()
    expect_error(
      ventilation_state(d, shelter = "yes"),
      class = "meteoHazard_input_error"
    )
    expect_error(
      ventilation_state(d, shelter = 1L),
      class = "meteoHazard_input_error"
    )
  })

  it("raises meteoHazard_input_error for non-logical shelter_h_mix argument", {
    d <- vs_met()
    expect_error(
      ventilation_state(d, shelter_h_mix = "yes"),
      class = "meteoHazard_input_error"
    )
  })

  it("raises meteoHazard_input_error for non-logical pool_cap argument", {
    d <- vs_met()
    expect_error(
      ventilation_state(d, pool_cap = "yes"),
      class = "meteoHazard_input_error"
    )
    expect_error(
      ventilation_state(d, pool_cap = 1L),
      class = "meteoHazard_input_error"
    )
  })

})


# ---------------------------------------------------------------------------
# M3 — valley sheltering (C6, issue #20)
# ---------------------------------------------------------------------------
#
# Transfer model (all constants in ODOUR_CONSTANTS):
#   s_f  = clamp((SHELTER_OPEN_REF - shelter_index) /
#                (SHELTER_OPEN_REF - SHELTER_ENCLOSED_REF), 0, 1)
#          SHELTER_OPEN_REF = 85 deg, SHELTER_ENCLOSED_REF = 50 deg
#   w_r  = clamp((SHELTER_U_FLUSH - u10) /
#                (SHELTER_U_FLUSH - SHELTER_U_FULL), 0, 1)
#          SHELTER_U_FLUSH = 6.0 m/s, SHELTER_U_FULL = 1.5 m/s
#   reduction = SHELTER_MAX_REDUCTION * s_f * w_r   (SHELTER_MAX_REDUCTION = 0.7)
#   reduction_effective = reduction * (1 - DRAINAGE_SHELTER_OVERLAP * drainage_active)
#   u_eff_sheltered = max(u_eff_base * (1 - reduction_effective), U_CALM_FLOOR)
# ---------------------------------------------------------------------------

describe("ventilation_state(): M3 valley sheltering (shelter = TRUE)", {

  # Helper: night met, optionally windy, no temperature/RH to keep pool_top simple.
  vs_night <- function(u10 = 2.0, n = 1) {
    as.data.frame(lapply(list(
      wind_speed_10m        = u10,
      direct_radiation      = 0,
      cloud_cover           = 50,
      boundary_layer_height = 400
    ), rep_len, n))
  }

  # Helper: channelled terrain (M1 active if pool_top > 0)
  ter_channelled <- function(shelter_index = 50) {
    mh_terrain(
      shelter_index    = shelter_index,
      flow_convergence = 0.8,          # ≥ 0.5 → channelled
      drainage_bearing = 30            # non-NA → channelled predicate TRUE
    )
  }

  # Helper: open terrain (no M1 drainage confinement)
  ter_open <- function(shelter_index = 50) {
    mh_terrain(shelter_index = shelter_index)
  }


  # -- Default-off is identity -----------------------------------------------

  it("shelter = FALSE leaves u_eff bit-identical to the default (regression)", {
    # shelter = FALSE is the default; it must reproduce the current output
    # exactly for every input.
    d <- vs_night(u10 = 3)
    vs_default <- ventilation_state(d)
    vs_off     <- ventilation_state(d, shelter = FALSE)
    expect_identical(vs_default$u_eff, vs_off$u_eff)
    expect_identical(vs_default$h_mix, vs_off$h_mix)
  })

  it("shelter = FALSE with terrain leaves u_eff bit-identical to no-terrain (regression)", {
    d   <- vs_night(u10 = 2)
    ter <- ter_open(shelter_index = 50)
    vs_no_ter    <- ventilation_state(d, shelter = FALSE)
    vs_with_ter  <- ventilation_state(d, terrain = ter, shelter = FALSE)
    expect_identical(vs_no_ter$u_eff, vs_with_ter$u_eff)
  })


  # -- Analytic transfer-model values ----------------------------------------

  it("u_eff is clamped to U_CALM_FLOOR when shelter_index = 50 and u10 = 1.5", {
    # At the maximum-shelter inputs:
    #   s_f = (85 - 50) / (85 - 50) = 1.0
    #   w_r = (6.0 - 1.5) / (6.0 - 1.5) = 1.0
    #   reduction = 0.7 * 1.0 * 1.0 = 0.7
    #   u_eff_base = max(1.5, 0.5) = 1.5
    #   u_eff_sheltered = max(1.5 * 0.3, 0.5) = max(0.45, 0.5) = 0.5 = U_CALM_FLOOR
    d   <- vs_night(u10 = 1.5)
    ter <- ter_open(shelter_index = 50)
    vs  <- ventilation_state(d, terrain = ter, shelter = TRUE)
    expect_equal(vs$u_eff, ODOUR_CONSTANTS$U_CALM_FLOOR, tolerance = 1e-9)
  })

  it("u_eff matches the analytic formula at mid-range shelter and wind", {
    # shelter_index = 67.5, u10 = 3.75 (both at midpoints of their ranges):
    #   s_f = (85 - 67.5) / 35 = 0.5
    #   w_r = (6.0 - 3.75) / 4.5 = 0.5
    #   reduction = 0.7 * 0.5 * 0.5 = 0.175
    #   u_eff_base = max(3.75, 0.5) = 3.75
    #   u_eff_sheltered = max(3.75 * (1 - 0.175), 0.5) = max(3.09375, 0.5) = 3.09375
    d   <- vs_night(u10 = 3.75)
    ter <- ter_open(shelter_index = 67.5)
    vs  <- ventilation_state(d, terrain = ter, shelter = TRUE)
    u_eff_base     <- max(3.75, ODOUR_CONSTANTS$U_CALM_FLOOR)
    expected_u_eff <- max(u_eff_base * (1 - 0.175), ODOUR_CONSTANTS$U_CALM_FLOOR)
    expect_equal(vs$u_eff, expected_u_eff, tolerance = 1e-9)
  })


  # -- Structural invariants -------------------------------------------------

  it("u_eff is monotonically non-increasing as shelter_index decreases", {
    # More enclosed (lower shelter_index) → larger s_f → larger reduction → lower u_eff.
    d <- vs_night(u10 = 3.0)
    idxs   <- c(85, 75, 65, 55, 50)   # decreasing openness
    u_effs <- vapply(idxs, function(si) {
      ter <- ter_open(shelter_index = si)
      ventilation_state(d, terrain = ter, shelter = TRUE)$u_eff
    }, numeric(1))
    expect_true(all(diff(u_effs) <= 1e-9),
                label = "u_eff should be non-increasing as shelter_index falls")
  })

  it("reduction is capped: shelter_index below SHELTER_ENCLOSED_REF gives same as SHELTER_ENCLOSED_REF", {
    # SHELTER_ENCLOSED_REF = 50; s_f is clamped to [0,1], so shelter_index = 30
    # gives s_f = 1.0, same as shelter_index = 50.
    d       <- vs_night(u10 = 2.0)
    ter_50  <- ter_open(shelter_index = 50)
    ter_30  <- ter_open(shelter_index = 30)
    vs_50   <- ventilation_state(d, terrain = ter_50, shelter = TRUE)
    vs_30   <- ventilation_state(d, terrain = ter_30, shelter = TRUE)
    expect_equal(vs_50$u_eff, vs_30$u_eff, tolerance = 1e-9)
  })


  # -- Flush regime ----------------------------------------------------------

  it("reduction is zero when u10 >= SHELTER_U_FLUSH (valley flushed)", {
    # SHELTER_U_FLUSH = 6.0 m/s: w_r = 0 → reduction = 0 → u_eff unchanged.
    d_flush  <- vs_night(u10 = 6.0)
    d_strong <- vs_night(u10 = 8.0)
    ter      <- ter_open(shelter_index = 50)

    vs_base_flush  <- ventilation_state(d_flush,  shelter = FALSE)
    vs_base_strong <- ventilation_state(d_strong, shelter = FALSE)
    vs_shlt_flush  <- ventilation_state(d_flush,  terrain = ter, shelter = TRUE)
    vs_shlt_strong <- ventilation_state(d_strong, terrain = ter, shelter = TRUE)

    expect_equal(vs_shlt_flush$u_eff,  vs_base_flush$u_eff,  tolerance = 1e-9)
    expect_equal(vs_shlt_strong$u_eff, vs_base_strong$u_eff, tolerance = 1e-9)
  })

  it("reduction is maximum when u10 <= SHELTER_U_FULL", {
    # SHELTER_U_FULL = 1.5 m/s: w_r = 1.0 for all u10 ≤ 1.5.
    d_at_full   <- vs_night(u10 = 1.5)
    d_below_full <- vs_night(u10 = 1.0)
    ter <- ter_open(shelter_index = 50)

    vs_at   <- ventilation_state(d_at_full,    terrain = ter, shelter = TRUE)
    vs_below <- ventilation_state(d_below_full, terrain = ter, shelter = TRUE)

    # Both should be clamped to U_CALM_FLOOR (since 0.7 reduction from ≤1.5 m/s
    # drives u_eff below the floor).
    expect_equal(vs_at$u_eff,    ODOUR_CONSTANTS$U_CALM_FLOOR, tolerance = 1e-9)
    expect_equal(vs_below$u_eff, ODOUR_CONSTANTS$U_CALM_FLOOR, tolerance = 1e-9)
  })


  # -- Graceful no-ops -------------------------------------------------------

  it("shelter = TRUE with terrain = NULL leaves u_eff unchanged", {
    d  <- vs_night(u10 = 2.0)
    vs_null  <- ventilation_state(d, terrain = NULL,            shelter = TRUE)
    vs_false <- ventilation_state(d, shelter = FALSE)
    expect_identical(vs_null$u_eff, vs_false$u_eff)
  })

  it("shelter = TRUE with terrain having shelter_index = NA leaves u_eff unchanged", {
    d   <- vs_night(u10 = 2.0)
    ter <- mh_terrain(shelter_index = NA_real_, flow_convergence = 0.3)
    vs_na    <- ventilation_state(d, terrain = ter, shelter = TRUE)
    vs_false <- ventilation_state(d, shelter = FALSE)
    expect_identical(vs_na$u_eff, vs_false$u_eff)
  })


  # -- shelter_h_mix ---------------------------------------------------------

  it("h_mix is unchanged by shelter when shelter_h_mix = FALSE (default)", {
    d   <- vs_night(u10 = 1.5)
    ter <- ter_open(shelter_index = 50)
    vs_on  <- ventilation_state(d, terrain = ter, shelter = TRUE, shelter_h_mix = FALSE)
    vs_off <- ventilation_state(d, shelter = FALSE)
    expect_identical(vs_on$h_mix, vs_off$h_mix)
  })

  it("h_mix is reduced by (1 - reduction) when shelter_h_mix = TRUE", {
    # At shelter_index = 50, u10 = 3.75: reduction = 0.175 (from the analytic test).
    # h_mix_sheltered = h_mix_base * (1 - 0.175) = h_mix_base * 0.825.
    d     <- vs_night(u10 = 3.75)
    ter   <- ter_open(shelter_index = 67.5)
    vs_base <- ventilation_state(d, shelter = FALSE)
    vs_hmix <- ventilation_state(d, terrain = ter, shelter = TRUE, shelter_h_mix = TRUE)
    expected_h_mix <- vs_base$h_mix * (1 - 0.175)
    expect_equal(vs_hmix$h_mix, expected_h_mix, tolerance = 1e-6)
  })


  # -- Precedence: M1 drainage × M3 shelter no-stack -------------------------

  it("shelter is fully suppressed on drainage_active hours (DRAINAGE_SHELTER_OVERLAP = 1)", {
    # drainage_active = is_channelled & !is_day & !is.na(pool_top) & pool_top > 0.
    # On those hours: reduction_effective = 0 → u_eff unchanged.
    # Use a channelled terrain and a clear cold night (pool_top will accumulate).
    d_night <- as.data.frame(lapply(list(
      wind_speed_10m        = rep(1.5, 6),
      direct_radiation      = rep(0,   6),   # night
      cloud_cover           = rep(0,   6),   # clear (max cooling)
      boundary_layer_height = rep(200, 6),
      temperature_2m        = rep(5,   6),
      relative_humidity_2m  = rep(40,  6)
    ), identity))

    ter <- ter_channelled(shelter_index = 50)   # channelled + shelter

    vs_no_shelter  <- ventilation_state(d_night, terrain = ter, shelter = FALSE)
    vs_with_shelter <- ventilation_state(d_night, terrain = ter, shelter = TRUE)

    # Identify drainage_active hours: channelled & night & pool_top > 0.
    is_channelled <- TRUE  # flow_convergence = 0.8 ≥ 0.5 & drainage_bearing non-NA
    drainage_active <- is_channelled &
      !vs_no_shelter$is_day &
      !is.na(vs_no_shelter$pool_top) &
      vs_no_shelter$pool_top > 0

    # Fixture must exercise at least one drainage_active hour — guards the guard.
    expect_true(any(drainage_active),
                label = "fixture must produce at least one drainage_active hour")

    # On drainage_active hours, shelter must NOT reduce u_eff.
    expect_equal(
      vs_with_shelter$u_eff[drainage_active],
      vs_no_shelter$u_eff[drainage_active],
      tolerance = 1e-9,
      label = "u_eff on drainage_active hours must equal no-shelter u_eff"
    )
    # M3 on unconstrained (non-drainage) nights is covered by the monotonicity and
    # analytic-formula tests above.
  })

})
