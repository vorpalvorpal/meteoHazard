#' Vertical dust flux for an exposed landfill surface
#'
#' Computes an hourly vertical dust (PM) emission flux for an exposed, erodible
#' surface using a physical saltation-to-emission chain: a Shao & Lu (2000)
#' threshold friction velocity, a Fecan et al. (1999) soil-moisture correction,
#' an Owen/White saltation flux, and the Marticorena & Bergametti (1995)
#' sandblasting efficiency. The erodible surface is assumed smooth (no
#' non-erodible roughness drag partition). The returned flux is in relative
#' units; for an operational, bounded index use [generate_dust_risk_index()].
#'
#' The function does not query any API. The caller supplies the per-hour
#' meteorological vectors (from Open-Meteo) and one-time site-survey parameters.
#'
#' @param tyler_sieve_no Integer Tyler Standard Sieve number for the modal
#'   aggregate size of the erodible surface. Must be one of the tabulated values
#'   (see `TYLER_SIEVE_DIAMETERS_M`).
#' @param clay_percent Clay fraction of the surface material (% by mass), in
#'   `0`–`100`.
#' @param wind_speed_10m Numeric vector. Hourly mean wind speed at 10 m (km/h).
#' @param wind_gusts_10m Numeric vector. Peak gust at 10 m (km/h).
#' @param soil_moisture Numeric vector. Volumetric water content of the top
#'   0--1 cm layer (m^3/m^3), in `0`–`1`.
#' @param z0 Aerodynamic roughness length for the wind profile (m). Default
#'   0.005.
#' @param bulk_density Dry bulk density (Mg/m^3). Default 1.6.
#' @param gust_factor Gust-duration factor converting the 3-second gust to the
#'   fastest-mile driving wind. Default 0.84 (Durst).
#' @param threshold_multiplier Multiplier on the threshold friction velocity,
#'   length 1 or matching the met vectors. Used by [generate_dust_risk_index()]
#'   to inject the crust-persistence factor; defaults to 1 (no effect).
#'
#' @return Numeric vector of vertical dust flux (relative units), one per hour;
#'   `0` where the wind does not exceed the (possibly moisture/crust/roughness
#'   adjusted) threshold.
#'
#' @references
#' Shao, Y. & Lu, H. (2000) \doi{10.1029/2000JD900304}.
#' Fecan, F., Marticorena, B. & Bergametti, G. (1999)
#' \doi{10.1007/s00585-999-0149-7}.
#' Marticorena, B. & Bergametti, G. (1995) \doi{10.1029/95JD00690}.
#' Owen, P.R. (1964) saltation of uniform grains in air; White, B.R. (1979).
#'
#' @seealso [generate_dust_risk_index()] for the bounded operational index.
#' @export
dust_emission_potential <- function(
  tyler_sieve_no,
  clay_percent,
  wind_speed_10m,
  wind_gusts_10m,
  soil_moisture,
  z0                   = 0.005,
  bulk_density         = 1.6,
  gust_factor          = 0.84,
  threshold_multiplier = 1
) {
  # ---- Validation ---------------------------------------------------------- #
  if (!as.character(tyler_sieve_no) %in% names(TYLER_SIEVE_DIAMETERS_M)) {
    cli::cli_abort(c(
      "Invalid {.arg tyler_sieve_no}: {tyler_sieve_no}.",
      "i" = "Valid Tyler Sieve numbers are: {.val {as.integer(names(TYLER_SIEVE_DIAMETERS_M))}}."
    ))
  }
  checkmate::assert_number(clay_percent, lower = 0, upper = 100)
  n <- length(wind_speed_10m)
  checkmate::assert_numeric(wind_speed_10m, lower = 0, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(wind_gusts_10m, lower = 0, any.missing = FALSE, len = n)
  checkmate::assert_numeric(soil_moisture, lower = 0, upper = 1, any.missing = FALSE, len = n)
  checkmate::assert_number(z0, lower = 1e-9, upper = 10 - 1e-9)
  checkmate::assert_number(bulk_density, lower = 1e-9)
  checkmate::assert_number(gust_factor, lower = 0, upper = 1)
  checkmate::assert_numeric(threshold_multiplier, lower = 0, any.missing = FALSE)

  # ---- Physical constants -------------------------------------------------- #
  A_N   <- 0.0123     # Shao & Lu (2000) dimensionless coefficient
  rho_p <- 2650.0     # particle density (kg/m^3, quartz)
  rho_a <- 1.225      # air density (kg/m^3)
  g     <- 9.81
  gamma <- 3.0e-4     # interparticle cohesion parameter (N/m)
  kappa <- 0.40       # von Karman constant
  z     <- 10.0       # reference height (m)

  d <- TYLER_SIEVE_DIAMETERS_M[[as.character(tyler_sieve_no)]]

  # ---- Dry threshold friction velocity (Shao & Lu 2000) -------------------- #
  u_star_t_dry <- sqrt(A_N * (rho_p / rho_a * g * d + gamma / (rho_a * d)))

  # ---- Moisture correction (Fecan et al. 1999) ----------------------------- #
  # Gravimetric moisture (% by mass) from volumetric content and bulk density.
  w       <- soil_moisture / bulk_density * 100
  w_prime <- 0.0014 * clay_percent^2 + 0.17 * clay_percent
  f_moist <- ifelse(w > w_prime, sqrt(1 + 1.21 * (w - w_prime)^0.68), 1)

  # ---- Combined threshold: wetness vs crust hand off via the maximum ------- #
  u_star_t <- u_star_t_dry * pmax(f_moist, threshold_multiplier)

  # ---- Effective friction velocity (gust-driven, fastest-mile proxy) ------- #
  U_fm   <- pmax(wind_speed_10m, gust_factor * wind_gusts_10m) / 3.6  # km/h -> m/s
  u_star <- kappa * U_fm / log(z / z0)

  # ---- Saltation flux (Owen 1964 / White 1979) ----------------------------- #
  # Q = (rho_a/g) u*^3 (1 - (u*t/u*)^2), zero below threshold. The excess factor
  # is computed only where u* exceeds the threshold, so it is never negative.
  excess  <- ifelse(u_star > u_star_t, 1 - (u_star_t / u_star)^2, 0)
  Q       <- (rho_a / g) * u_star^3 * excess

  # ---- Vertical dust flux: MB95 sandblasting efficiency -------------------- #
  # alpha depends only on clay, so it cancels in the reference-normalised index;
  # retained so this engine yields a physically-meaningful flux.
  alpha <- 10^(0.134 * clay_percent - 6)

  alpha * Q
}


#' Dust Hazard Index for landfill operations
#'
#' Computes an hourly Dust Hazard Index in the range `[0, 100]` for each row of a
#' meteorological forecast tibble. Wraps [dust_emission_potential()] and
#' normalises the relative dust flux against a reference condition (a strong gust
#' on a dry surface), so the index keeps resolution across ordinary winds rather
#' than saturating.
#'
#' @param met_data A tibble (or data frame), one row per hourly timestep, with at
#'   least `wind_speed_10m` (km/h), `wind_gusts_10m` (km/h), and
#'   `soil_moisture_0_to_1cm` (m^3/m^3); plus `precipitation` (mm) when
#'   `crust = TRUE`.
#' @param tyler_sieve_no,clay_percent,z0,bulk_density,gust_factor
#'   Site and model parameters forwarded to [dust_emission_potential()].
#' @param crust Logical. Enable the precipitation-driven crust-persistence gate.
#'   Default `FALSE`. When `TRUE`, `met_data$precipitation` is required.
#' @param rain_crust_threshold Precipitation (mm) at or above which an hour is
#'   treated as a crust-forming rain event. Default 2.
#' @param crust_factor_max Maximum threshold multiplier immediately after rain.
#'   Default 3.
#' @param crust_decay_hours E-folding time (hours) over which the crust
#'   suppression decays. Default 72.
#' @param scale_ref_gust Gust speed (km/h) that, on a dry crust-free surface,
#'   maps to index 100. Default 65. Must exceed the entrainment threshold.
#'
#' @return Numeric vector of length `nrow(met_data)`, the Dust Hazard Index in
#'   `[0, 100]` for each forecast hour.
#'
#' @seealso [dust_emission_potential()].
#' @export
generate_dust_risk_index <- function(
  met_data,
  tyler_sieve_no       = 20L,
  clay_percent         = 10,
  z0                   = 0.005,
  bulk_density         = 1.6,
  gust_factor          = 0.84,
  crust                = FALSE,
  rain_crust_threshold = 2,
  crust_factor_max     = 3,
  crust_decay_hours    = 72,
  scale_ref_gust       = 65
) {
  checkmate::assert_data_frame(met_data, min.rows = 1)
  checkmate::assert_flag(crust)
  checkmate::assert_number(rain_crust_threshold, lower = 0)
  checkmate::assert_number(crust_factor_max, lower = 1)
  checkmate::assert_number(crust_decay_hours, lower = 1e-9)
  checkmate::assert_number(scale_ref_gust, lower = 1e-9)

  required_cols <- c("wind_speed_10m", "wind_gusts_10m", "soil_moisture_0_to_1cm")
  if (crust) required_cols <- c(required_cols, "precipitation")
  missing_cols <- setdiff(required_cols, names(met_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "{.arg met_data} is missing required columns: {.val {missing_cols}}.",
      "i" = "Required: wind_speed_10m (km/h), wind_gusts_10m (km/h), soil_moisture_0_to_1cm (m³/m³){if (crust) ', precipitation (mm)' else ''}."
    ))
  }

  # ---- Crust factor per hour (threshold multiplier) ------------------------ #
  crust_mult <- if (crust) {
    .dust_crust_factor(met_data$precipitation, rain_crust_threshold,
                       crust_factor_max, crust_decay_hours)
  } else {
    1
  }

  common <- list(
    tyler_sieve_no = tyler_sieve_no, clay_percent = clay_percent,
    z0 = z0, bulk_density = bulk_density, gust_factor = gust_factor
  )

  flux <- do.call(dust_emission_potential, c(common, list(
    wind_speed_10m       = met_data$wind_speed_10m,
    wind_gusts_10m       = met_data$wind_gusts_10m,
    soil_moisture        = met_data$soil_moisture_0_to_1cm,
    threshold_multiplier = crust_mult
  )))

  # Reference flux: the scaling gust on a dry, crust-free surface, calm mean wind.
  flux_ref <- do.call(dust_emission_potential, c(common, list(
    wind_speed_10m       = 0,
    wind_gusts_10m       = scale_ref_gust,
    soil_moisture        = 0,
    threshold_multiplier = 1
  )))
  if (flux_ref <= 0) {
    cli::cli_abort(c(
      "{.arg scale_ref_gust} ({scale_ref_gust} km/h) does not exceed the erosion threshold.",
      "i" = "Increase {.arg scale_ref_gust} so the reference condition produces dust."
    ))
  }

  pmin(100, 100 * flux / flux_ref)
}


# ---- Internal helpers ------------------------------------------------------ #

# Per-hour crust threshold multiplier. age = hours since the most recent hour
# with precipitation >= threshold (Inf before any such hour). The crust is
# strongest right after rain and decays exponentially over decay_hours.
.dust_crust_factor <- function(precipitation, threshold, factor_max, decay_hours) {
  checkmate::assert_numeric(precipitation, lower = 0, any.missing = FALSE, min.len = 1)
  n   <- length(precipitation)
  age <- numeric(n)
  current <- Inf
  for (i in seq_len(n)) {
    if (precipitation[i] >= threshold) {
      current <- 0
    } else if (is.finite(current)) {
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
