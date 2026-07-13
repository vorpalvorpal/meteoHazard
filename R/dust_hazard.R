#' Vertical dust flux for an exposed landfill surface
#'
#' Computes an hourly vertical dust (PM) emission flux for an exposed, erodible
#' surface using a physical saltation-to-emission chain: a Shao & Lu (2000)
#' threshold friction velocity, a Fecan et al. (1999) soil-moisture correction,
#' an Owen (1964) / White (1979) saltation flux, and the Marticorena &
#' Bergametti (1995, "MB95") sandblasting efficiency. The erodible surface is
#' assumed smooth (no non-erodible roughness drag partition; see `z0`). The
#' returned value is proportional to the instantaneous vertical dust
#' mass-flux rate, returned **unscaled** for relative ranking; for an
#' operational, bounded index use [dust_hazard()].
#'
#' The function does not query any API. The caller supplies the per-hour
#' meteorological vectors (from Open-Meteo) and one-time site-survey parameters.
#'
#' @section Units:
#' The dimensional inputs (`wind_speed_10m`, `wind_gusts_10m`, `z0`,
#' `bulk_density`) may be supplied as bare numerics in the documented unit or as
#' \pkg{units} objects, which are converted automatically (a dimensionally
#' incompatible unit is an error). `clay_percent`, `soil_moisture` (m^3/m^3) and
#' `gust_factor` are dimensionless and taken as-is. The returned value is
#' proportional to the instantaneous vertical dust mass-flux rate
#' (`alpha * Q`; physical dimensions kg m^-2 s^-1), but is returned **unscaled**
#' (plain numeric, not a \pkg{units} object) for relative ranking pending
#' calibration (issues #11/#8). Two unit notes for anyone converting to an
#' absolute rate: (1) the MB95 `alpha` coefficient is used in its native
#' cm^-1 units and is NOT converted to SI (m^-1, which would require `* 100`);
#' (2) this is an instantaneous rate, not an hourly total — the per-hour mass
#' is `F * 3600`.
#'
#' @section Idealisations:
#' This engine makes several simplifying assumptions, in decreasing order of
#' expected impact on absolute accuracy (the relative ranking across hours is
#' more robust than the absolute value):
#' \itemize{
#'   \item **Steady hourly gust forcing.** The gust is applied as a constant
#'     hourly-steady friction velocity rather than integrated over a
#'     within-hour wind-speed distribution. Because saltation flux is a
#'     convex (odd-power, threshold-gated) function of wind speed, this biases
#'     the flux high near the threshold and is a placeholder for a future
#'     within-hour intermittency treatment (see the T8 TODO below).
#'   \item **MB95 alpha not SI-converted, and capped at 20% clay.** `alpha` is
#'     kept in its native cm^-1 units (see @section Units) and is clamped at
#'     the clay % where the Marticorena & Bergametti (1995) fit was validated
#'     (Foroutan et al. 2017; Kok et al. 2014); real high-clay soils aggregate
#'     and emit less than a naive extrapolation would predict.
#'   \item **Fixed reference air density.** `rho_a` is a constant
#'     (`DUST_CONSTANTS$RHO_A_REF`), not adjusted for site temperature or
#'     pressure (see the T10 TODO below).
#'   \item **Smooth surface, no drag partition.** `z0` defaults to the
#'     smooth-bed value `d/30`; non-erodible roughness elements (gravel,
#'     vegetation, stockpiles) are not represented, so u* — and hence the flux
#'     — is not sheltered the way a drag-partition model would predict
#'     (Raupach et al. 1993; Okin 2008; Webb et al. 2014).
#'   \item **Transport-limited, not supply-limited.** The Owen/White saltation
#'     flux assumes an unlimited erodible reservoir. Real surfaces can be
#'     supply-limited (surface armouring, prior deflation), which this model
#'     does not represent; sustained high winds are more likely to be
#'     overestimated than gusty, intermittent ones (Gillette et al. 1982;
#'     Kok et al. 2012). A future reservoir/hybrid treatment is a TODO, not
#'     implemented here.
#' }
#'
#' @param tyler_sieve_no Integer Tyler Standard Sieve number for the modal
#'   aggregate size of the erodible surface. Must be one of the tabulated values
#'   (see `TYLER_SIEVE_DIAMETERS_M`).
#' @param clay_percent Clay fraction of the surface material (% by mass), in
#'   `0`–`100`. Values above `DUST_CONSTANTS$MB95_CLAY_CAP` (20%) are capped
#'   for the MB95 sandblasting efficiency only (with a warning); the Fecan
#'   moisture correction still uses the true `clay_percent`.
#' @param wind_speed_10m Numeric vector. Hourly mean wind speed at 10 m (m/s).
#' @param wind_gusts_10m Numeric vector. Peak gust at 10 m (m/s). Must be
#'   `>= wind_speed_10m` at every hour.
#' @param soil_moisture Numeric vector. Volumetric water content of the top
#'   0--1 cm layer (m^3/m^3), in `0`–`1`.
#' @param z0 Aerodynamic roughness length for the wind profile (m). Default
#'   `NULL`, which derives the smooth-bed value `d/30` (Nikuradse
#'   equivalent-sand roughness; the MB95 `z0s` convention) from the modal
#'   grain diameter implied by `tyler_sieve_no`. A caller-supplied `z0` larger
#'   than this smooth-bed value warns, because larger roughness is not
#'   sheltering the bed here (no drag partition is modelled) and so instead
#'   biases u*, and hence the flux, high.
#' @param bulk_density Dry bulk density (Mg/m^3). Default 1.6.
#' @param gust_factor Gust-duration factor converting the 3-second gust to the
#'   fastest-mile driving wind. Default 0.84 (Durst).
#' @param threshold_multiplier Multiplier on the threshold friction velocity,
#'   length 1 or matching the met vectors. Used by [dust_hazard()]
#'   to inject the crust-persistence factor; defaults to 1 (no effect).
#'
#' @return Numeric vector, one value per hour, proportional to the
#'   instantaneous vertical dust mass-flux rate (relative units, unscaled;
#'   see @section Units); `0` where the wind does not exceed the (possibly
#'   moisture/crust/roughness adjusted) threshold.
#'
#' @references
#' Shao, Y. & Lu, H. (2000) \doi{10.1029/2000JD900304} — threshold friction
#' velocity (`A_N`, `gamma`).
#' Fecan, F., Marticorena, B. & Bergametti, G. (1999)
#' \doi{10.1007/s00585-999-0149-7} — soil-moisture threshold correction.
#' Marticorena, B. & Bergametti, G. (1995) \doi{10.1029/95JD00690} — drag
#' partition and sandblasting efficiency (`z0s = d/30`,
#' `alpha = 10^(0.134*clay - 6)` cm^-1).
#' Owen, P.R. (1964) J. Fluid Mech. 20:225; White, B.R. (1979) J. Geophys.
#' Res. 84:4643 — saltation flux forms.
#' Raupach, M.R., Gillette, D.A. & Leys, J.F. (1993) \doi{10.1029/92JD01922}
#' — roughness elements raise the threshold / shelter the bed.
#' Okin, G.S. (2008) J. Geophys. Res.-Earth Surf. 113:F02S10 — roughness
#' suppresses emission (patchy-vegetation model).
#' Webb, N.P. et al. (2014) \doi{10.1002/2014JD021491} — roughness and
#' shear-stress partitioning.
#' Foroutan, H. et al. (2017) Geosci. Model Dev. 10:2591 (PMC6145470) — MB95
#' alpha validity limited to <=20% clay.
#' Kok, J.F. et al. (2012) Rep. Prog. Phys. 75:106901 — saltation flux;
#' transport- vs supply-limited regimes.
#' Kok, J.F. et al. (2014) Atmos. Chem. Phys. 14:13043 — clay fraction
#' commonly capped at 0.2 in dust-emission schemes.
#' Gillette, D.A. et al. (1982) J. Geophys. Res. 87:9003 — supply limitation
#' and crust rupture by saltation impact.
#' Nikuradse, J. (1933); Sherman, D.J. (1992) — smooth-bed `z0 = d/30`.
#'
#' @seealso [dust_hazard()] for the bounded operational index.
#' @export
dust_flux <- function(
  tyler_sieve_no,
  clay_percent,
  wind_speed_10m,
  wind_gusts_10m,
  soil_moisture,
  z0                   = NULL,
  bulk_density         = 1.6,
  gust_factor          = 0.84,
  threshold_multiplier = 1
) {
  # ---- Normalise dimensional inputs (bare = documented unit; units = converted) #
  # clay_percent, soil_moisture (m^3/m^3 ratio) and gust_factor are dimensionless
  # and taken as-is; the returned flux is in relative units (plain numeric).
  wind_speed_10m <- .drop_to(wind_speed_10m, "m/s",    arg = "wind_speed_10m")
  wind_gusts_10m <- .drop_to(wind_gusts_10m, "m/s",    arg = "wind_gusts_10m")
  # z0 = NULL means "derive the smooth-bed value below"; only normalise/validate
  # a caller-supplied value.
  if (!is.null(z0)) {
    z0 <- .drop_to(z0, "m", arg = "z0")
  }
  bulk_density   <- .drop_to(bulk_density,   "Mg/m^3", arg = "bulk_density")

  # ---- Validation ---------------------------------------------------------- #
  if (!as.character(tyler_sieve_no) %in% names(TYLER_SIEVE_DIAMETERS_M)) {
    cli::cli_abort(c(
      "Invalid {.arg tyler_sieve_no}: {tyler_sieve_no}.",
      "i" = "Valid Tyler Sieve numbers are: {.val {as.integer(names(TYLER_SIEVE_DIAMETERS_M))}}."
    ), class = "meteoHazard_input_error")
  }
  checkmate::assert_number(clay_percent, lower = 0, upper = 100)
  n <- length(wind_speed_10m)
  checkmate::assert_numeric(wind_speed_10m, lower = 0, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(wind_gusts_10m, lower = 0, any.missing = FALSE, len = n)
  checkmate::assert_numeric(soil_moisture, lower = 0, upper = 1, any.missing = FALSE, len = n)
  if (!is.null(z0)) {
    # z0 must sit strictly below the log-law reference height Z_REF (10 m),
    # or log(z/z0) is <= 0 and u* is undefined/negative.
    checkmate::assert_number(z0, lower = 1e-9,
                             upper = DUST_CONSTANTS$Z_REF - 1e-9)
  }
  checkmate::assert_number(bulk_density, lower = 1e-9)
  checkmate::assert_number(gust_factor, lower = 0, upper = 1)
  checkmate::assert_numeric(threshold_multiplier, lower = 0, any.missing = FALSE)
  if (!length(threshold_multiplier) %in% c(1L, n)) {
    cli::cli_abort(
      "{.arg threshold_multiplier} must have length 1 or {n}, not {length(threshold_multiplier)}.",
      class = "meteoHazard_input_error"
    )
  }
  if (any(wind_gusts_10m < wind_speed_10m)) {
    cli::cli_abort(c(
      "{.arg wind_gusts_10m} must be >= {.arg wind_speed_10m} at every hour.",
      "x" = "Row(s) {.val {which(wind_gusts_10m < wind_speed_10m)}} have gust < mean wind."
    ), class = "meteoHazard_input_error")
  }

  # ---- Physical constants (DUST_CONSTANTS; see R/constants.R) -------------- #
  rho_p <- DUST_CONSTANTS$RHO_P     # particle density (kg/m^3, quartz)
  # Reference air density; not yet temperature/pressure adjusted (T10 TODO
  # below).
  rho_a <- DUST_CONSTANTS$RHO_A_REF
  g     <- DUST_CONSTANTS$G
  gamma <- DUST_CONSTANTS$GAMMA     # interparticle cohesion parameter (N/m)
  kappa <- DUST_CONSTANTS$KAPPA     # von Karman constant
  z     <- DUST_CONSTANTS$Z_REF     # reference height (m)

  d <- TYLER_SIEVE_DIAMETERS_M[[as.character(tyler_sieve_no)]]

  # ---- Smooth-bed roughness (Nikuradse; MB95 z0s convention) --------------- #
  # z0 = NULL (default) derives the smooth erodible-bed roughness d/30, i.e.
  # no non-erodible-element drag partition. A caller-supplied z0 greater than
  # this raises u* (rather than sheltering the bed via a drag partition, which
  # is not modelled here), biasing the flux high — warn instead of silently
  # trusting an unmodelled sheltering effect (Raupach et al. 1993; Okin 2008;
  # Webb et al. 2014). The guard is strict `>` so the exact smooth-bed value
  # itself never warns.
  z0_smooth <- d * DUST_CONSTANTS$Z0_SMOOTH_RATIO
  if (is.null(z0)) {
    z0 <- z0_smooth
  } else if (z0 > z0_smooth) {
    cli::cli_warn(c(
      "Supplied {.arg z0} ({z0} m) exceeds the smooth-bed roughness d/30 ({signif(z0_smooth, 3)} m).",
      "!" = "Non-erodible roughness is not represented (the drag partition is not modelled), so a larger z0 raises u* rather than sheltering the bed; the flux may be biased high.",
      "i" = "Leave {.arg z0} = NULL to use the smooth-bed value, or characterise roughness for a drag-partition treatment."
    ), class = "meteoHazard_dust_z0_rough")
  }

  # ---- Dry threshold friction velocity (Shao & Lu 2000) -------------------- #
  u_star_t_dry <- sqrt(DUST_CONSTANTS$A_N * (rho_p / rho_a * g * d + gamma / (rho_a * d)))

  # ---- Moisture correction (Fecan et al. 1999) ----------------------------- #
  # Gravimetric moisture (% by mass) from volumetric content and bulk density;
  # `* 100` is a unit conversion (fraction -> percent), not a physics constant.
  w       <- soil_moisture / bulk_density * 100
  w_prime <- DUST_CONSTANTS$FECAN_WP_QUAD * clay_percent^2 + DUST_CONSTANTS$FECAN_WP_LIN * clay_percent
  # w' always uses the TRUE clay_percent (unlike the MB95 alpha below): the
  # moisture-threshold physics has no analogous >20% validity ceiling.
  f_moist <- ifelse(
    w > w_prime,
    sqrt(1 + DUST_CONSTANTS$FECAN_A * (w - w_prime)^DUST_CONSTANTS$FECAN_B),
    1
  )

  # ---- Combined threshold: wetness vs crust hand off via the maximum ------- #
  u_star_t <- u_star_t_dry * pmax(f_moist, threshold_multiplier)

  # ---- Effective friction velocity (gust-driven, fastest-mile proxy) ------- #
  U_fm   <- pmax(wind_speed_10m, gust_factor * wind_gusts_10m)  # m/s
  # Log law over the smooth-bed z0; the gust is applied as an hourly-steady
  # forcing (idealisation — see @section Idealisations: biases high near
  # threshold; T8 TODO below covers within-hour intermittency).
  u_star <- kappa * U_fm / log(z / z0)

  # ---- Saltation flux (Owen 1964 / White 1979) ----------------------------- #
  # Q = (rho_a/g) u*^3 (1 - (u*t/u*)^2), zero below threshold. The excess factor
  # is computed only where u* exceeds the threshold, so it is never negative.
  # This is transport-limited (unlimited erodible reservoir assumed); see
  # @section Idealisations for the supply-limited caveat.
  excess  <- ifelse(u_star > u_star_t, 1 - (u_star_t / u_star)^2, 0)
  Q       <- (rho_a / g) * u_star^3 * excess

  # ---- Vertical dust flux: MB95 sandblasting efficiency -------------------- #
  # alpha depends only on clay, so for a fixed site it scales the flux without
  # re-ranking hours; retained so this engine yields a physically-meaningful
  # flux (the fixed reference-normalised index it once cancelled out of was
  # removed by issue #11). The MB95 fit
  # is only validated to 20% clay (Foroutan et al. 2017; Kok et al. 2014);
  # above that, real high-clay soils aggregate and emit less, so alpha is
  # capped at its 20%-clay value rather than extrapolated. alpha is in cm^-1
  # (MB95's native unit) and is NOT converted to SI (m^-1, which would need
  # `* 100`) — see @section Units.
  clay_eff <- min(clay_percent, DUST_CONSTANTS$MB95_CLAY_CAP)
  alpha <- 10^(DUST_CONSTANTS$MB95_ALPHA_SLOPE * clay_eff + DUST_CONSTANTS$MB95_ALPHA_INTERCEPT)
  if (clay_percent > DUST_CONSTANTS$MB95_CLAY_CAP) {
    cli::cli_warn(c(
      "{.arg clay_percent} ({clay_percent}%) exceeds the MB95 sandblasting validity limit ({DUST_CONSTANTS$MB95_CLAY_CAP}%).",
      "i" = "alpha capped at its {DUST_CONSTANTS$MB95_CLAY_CAP}% value; the MB95 fit is not calibrated for higher clay (real high-clay soils aggregate and emit less)."
    ), class = "meteoHazard_dust_clay_capped")
  }

  # TODO(dust-v3): T10 — add a direct diameter/d50 interface (bypassing
  # tyler_sieve_no) and a temperature-dependent air density (from
  # temperature_2m/pressure_msl via the ideal gas law) rather than the fixed
  # DUST_CONSTANTS$RHO_A_REF used above. Not implemented in this iteration.
  # TODO(dust-v3): T8 — integrate a within-hour Weibull wind-speed
  # distribution (Stout & Zobeck 1997; Cakmur et al. 2004) rather than the
  # steady hourly-gust forcing above, exposed via a future
  # `forcing = c("gust", "weibull")` argument (Comola et al. 2019 quantify the
  # intermittency bias this would correct; Martin & Kok 2017 note flux
  # scales ~u*^2 for the time-averaged saltation flux under intermittency,
  # vs. the instantaneous u*^3 used here). Not implemented in this iteration.

  alpha * Q
}


#' Dust hazard: crust-adjusted relative dust flux
#'
#' Computes an hourly **relative dust flux** for each row of a meteorological
#' forecast tibble: a `met_data`-frame convenience wrapper over [dust_flux()]
#' that adds the precipitation-driven crust-persistence gate. The returned flux
#' is in the same relative units as [dust_flux()] (a strong gust on a dry
#' surface produces a large value; a sub-threshold or wet/crusted hour produces
#' near-zero).
#'
#' @section Status — physical output, no fixed operational scale (issue #11):
#' This used to return a 0-100 index normalised against a reference gust. That
#' fixed 0-100 scale (and the `categorise_dust()` tiers) was removed: the
#' package now emits the physical relative flux, and mapping it onto a
#' site-specific operational index is a calibration step delivered by
#' forthcoming calibration tooling (issues #11/#8), not by fixed cut-points.
#'
#' @section Crust cold-start:
#' The precipitation-driven crust gate has no look-back before row 1 of
#' `met_data`: unless told otherwise, it assumes row 1 has no crust memory
#' (equivalent to "it has not rained in a very long time"). To represent a
#' forecast that starts partway through a crust-decay window, either (a)
#' prepend roughly `3 * crust_decay_hours` of history to `met_data` so the
#' crust state has settled by the row of interest, or (b) set
#' `hours_since_last_rain` to the actual elapsed time since the last
#' crust-forming rain before the forecast window starts.
#'
#' @section Idealisations:
#' The precipitation-driven crust factor is an uncalibrated exponential
#' clock-decay placeholder (`1 + (crust_factor_max - 1) * exp(-age /
#' crust_decay_hours)`): real crust breakdown is primarily mechanical
#' (abrasion/rupture by saltating-particle impact), not a pure function of
#' elapsed time (Gillette et al. 1982; Rice & McEwan 2001). A future version
#' should gate crust decay on cumulative saltation activity (a TODO, not
#' implemented here) rather than the clock alone. Separately, watering and
#' other surface management (the first-order dust-suppression control on an
#' active landfill working face) are invisible to the gridded
#' `soil_moisture_0_to_1cm` input; `threshold_multiplier` (via [dust_flux()])
#' is the intended injection hook — e.g. a caller can locally raise the
#' threshold multiplier during known watering hours before invoking
#' [dust_flux()] directly.
#'
#' @param met_data A tibble (or data frame), one row per hourly timestep, with at
#'   least `wind_speed_10m` (m/s), `wind_gusts_10m` (m/s), and
#'   `soil_moisture_0_to_1cm` (m^3/m^3); plus `precipitation` (mm) when
#'   `crust = TRUE`.
#' @param tyler_sieve_no,clay_percent,z0,bulk_density,gust_factor
#'   Site and model parameters forwarded to [dust_flux()].
#' @param crust Logical. Enable the precipitation-driven crust-persistence gate.
#'   Default `FALSE`. When `TRUE`, `met_data$precipitation` is required.
#' @param rain_crust_threshold Precipitation (mm) at or above which an hour is
#'   treated as a crust-forming rain event. Default 2.
#' @param crust_factor_max Maximum threshold multiplier immediately after rain.
#'   Default 3.
#' @param crust_decay_hours E-folding time (hours) over which the crust
#'   suppression decays. Default 72.
#' @param hours_since_last_rain Assumed crust age (hours) going into row 1 of
#'   `met_data`, i.e. how long since the last crust-forming rain before the
#'   forecast window starts. Default `Inf` (no crust memory at row 1,
#'   reproducing the pre-v3 behaviour). Set to `0` if the forecast is known to
#'   start immediately after a crust-forming rain event; see @section Crust
#'   cold-start.
#'
#' @return Numeric vector of length `nrow(met_data)`, the crust-adjusted
#'   **relative dust flux** for each forecast hour (same units as [dust_flux()];
#'   unbounded, `0` below the entrainment threshold). Issue #11 removed the fixed
#'   0-100 dust index in favour of this physical output.
#'
#' @seealso [dust_flux()].
#' @export
dust_hazard <- function(
  met_data,
  tyler_sieve_no       = 20L,
  clay_percent         = 10,
  z0                   = NULL,
  bulk_density         = 1.6,
  gust_factor          = 0.84,
  crust                = FALSE,
  rain_crust_threshold = 2,
  crust_factor_max     = 3,
  crust_decay_hours    = 72,
  hours_since_last_rain = Inf
) {
  # Normalise the dimensional scalars (bare = documented unit; units = converted).
  # met_data wind columns are normalised inside dust_flux(); crust_factor_max and
  # crust_decay_hours are dimensionless / in hours and taken as-is.
  rain_crust_threshold <- .drop_to(rain_crust_threshold, "mm", arg = "rain_crust_threshold")

  checkmate::assert_data_frame(met_data, min.rows = 1)
  checkmate::assert_flag(crust)
  checkmate::assert_number(rain_crust_threshold, lower = 0)
  checkmate::assert_number(crust_factor_max, lower = 1)
  checkmate::assert_number(crust_decay_hours, lower = 1e-9)
  checkmate::assert_number(hours_since_last_rain, lower = 0)

  required_cols <- c("wind_speed_10m", "wind_gusts_10m", "soil_moisture_0_to_1cm")
  if (crust) required_cols <- c(required_cols, "precipitation")
  .assert_required_cols(
    met_data, required_cols, arg = "met_data",
    info = paste0(
      "Required: wind_speed_10m (m/s), wind_gusts_10m (m/s), ",
      "soil_moisture_0_to_1cm (m³/m³)",
      if (crust) ", precipitation (mm)." else "."
    )
  )

  # ---- Crust factor per hour (threshold multiplier) ------------------------ #
  crust_mult <- if (crust) {
    .dust_crust_factor(.drop_to(met_data$precipitation, "mm", arg = "precipitation"),
                       rain_crust_threshold, crust_factor_max, crust_decay_hours,
                       age0 = hours_since_last_rain)
  } else {
    1
  }

  common <- list(
    tyler_sieve_no = tyler_sieve_no, clay_percent = clay_percent,
    z0 = z0, bulk_density = bulk_density, gust_factor = gust_factor
  )

  do.call(dust_flux, c(common, list(
    wind_speed_10m       = met_data$wind_speed_10m,
    wind_gusts_10m       = met_data$wind_gusts_10m,
    soil_moisture        = met_data$soil_moisture_0_to_1cm,
    threshold_multiplier = crust_mult
  )))
}

# TODO(dust-v3): T9 — add a receptor-aware dust_exposure() layer, mirroring
# the odour_hazard()/odour_exposure() split: reuse ventilation_state() (see
# R/odour-ventilation.R) for u_eff/h_mix/stability and the ISC3 area-source
# Gaussian-plume machinery already implemented for odour in
# R/odour_exposure.R, plus per-source surface parameters (tyler_sieve_no,
# clay_percent, z0, bulk_density) attached via mh_site roles. A thin
# dust_risk() wrapper (mirroring odour_risk()) would combine dust_hazard()
# and dust_exposure() in one call. Not implemented in this iteration — no new
# functions or arguments have been added for this.


# ---- Internal helpers ------------------------------------------------------ #

# Per-hour crust threshold multiplier. age = hours since the most recent hour
# with precipitation >= threshold (Inf before any such hour, or governed by
# age0 for row 1 — see dust_hazard()'s "Crust cold-start" section). The crust
# is strongest right after rain and decays exponentially over decay_hours.
# This exponential clock-decay is an uncalibrated placeholder; real crust
# breakdown is primarily mechanical (saltation-impact abrasion/rupture), not a
# pure function of elapsed time (Gillette et al. 1982; Rice & McEwan 2001). A
# future version should gate decay on cumulative saltation activity instead
# (TODO, not implemented here).
.dust_crust_factor <- function(precipitation, threshold, factor_max, decay_hours, age0 = Inf) {
  checkmate::assert_numeric(precipitation, lower = 0, any.missing = FALSE, min.len = 1)
  n   <- length(precipitation)
  age <- numeric(n)
  current <- age0
  for (i in seq_len(n)) {
    if (precipitation[i] >= threshold) {
      current <- 0
    } else if (i > 1 && is.finite(current)) {
      current <- current + 1
    }
    age[i] <- current
  }
  ifelse(is.finite(age), 1 + (factor_max - 1) * exp(-age / decay_hours), 1)
}


# Tyler Standard Sieve series: sieve number -> nominal opening in metres.
# Source: Tyler Industrial Products sieve series (US EPA AP-42 Table 13.2.5-1).
TYLER_SIEVE_DIAMETERS_M <- c(
  "3" = 0.006680,
  "4" = 0.004699,
  "5" = 0.003962,
  "6" = 0.003327,
  "8" = 0.002362,
  "9" = 0.001981,
  "10" = 0.001651,
  "14" = 0.001168,
  "20" = 0.000833,
  "24" = 0.000701,
  "28" = 0.000589,
  "32" = 0.000495,
  "35" = 0.000417,
  "42" = 0.000351,
  "48" = 0.000295,
  "60" = 0.000246,
  "65" = 0.000208,
  "80" = 0.000175,
  "100" = 0.000147,
  "115" = 0.000124,
  "150" = 0.000104,
  "170" = 0.000088,
  "200" = 0.000074,
  "250" = 0.000063,
  "270" = 0.000053,
  "325" = 0.000044,
  "400" = 0.000037
)
