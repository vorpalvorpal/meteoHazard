#' Solve globe temperature for a single observation
#'
#' Uses Newton's method to find the equilibrium globe temperature by
#' balancing solar absorption, radiative exchange, and convective
#' heat transfer.
#'
#' @param temp Air temperature in degrees Celsius.
#' @param wind_speed Wind speed in m/s.
#' @param direct_solar Direct beam solar radiation in W/m^2.
#' @param diffuse_solar Diffuse sky solar radiation in W/m^2.
#' @param zenith Solar zenith angle in degrees.
#' @param albedo Ground albedo (0--1). Default 0.12.
#' @param index Observation index for warnings (optional).
#'
#' @return Globe temperature in degrees Celsius.
#' @keywords internal
solve_globe_temp <- function(temp, wind_speed, direct_solar, diffuse_solar,
                             zenith, albedo = 0.12, index = NULL) {

  # Physical constants (shared via TWL_CONSTANTS)
  STEFAN_BOLTZMANN <- TWL_CONSTANTS$STEFAN_BOLTZMANN
  ABS_ZERO         <- TWL_CONSTANTS$ABS_ZERO
  GLOBE_DIAMETER   <- TWL_CONSTANTS$GLOBE_DIAMETER
  GLOBE_EMISSIVITY <- TWL_CONSTANTS$GLOBE_EMISSIVITY

  # Ensure minimum wind speed
  wind_speed <- max(wind_speed, 0.2)

  # Convective heat transfer coefficient for sphere (Kuehn & Goldstein)
  # Nu = 2 + 0.6 * Re^0.5 * Pr^0.33
  Re <- wind_speed * GLOBE_DIAMETER / TWL_CONSTANTS$AIR_KINEMATIC_VISCOSITY
  Pr <- 0.71
  Nu <- 2 + 0.6 * sqrt(Re) * Pr^0.33
  # Air thermal conductivity 0.028 W/(m·K): Brake 2001, Whillier Table 1 / BB p.470
  # Shared with the natural wet-bulb solver via TWL_CONSTANTS$AIR_THERMAL_CONDUCTIVITY
  hc_globe <- Nu * TWL_CONSTANTS$AIR_THERMAL_CONDUCTIVITY / GLOBE_DIAMETER

  # Solar radiation absorbed by globe
  cos_zenith <- cos(zenith * pi / 180)
  cos_zenith <- max(0, cos_zenith)

  # Direct beam on sphere (projected area / surface area = 0.25)
  direct_absorbed <- 0.25 * direct_solar * cos_zenith * GLOBE_EMISSIVITY

  # Diffuse from sky (hemisphere view factor ~0.5)
  diffuse_absorbed <- 0.5 * diffuse_solar * GLOBE_EMISSIVITY

  # Reflected from ground (hemisphere view factor ~0.5)
  total_horizontal <- direct_solar * cos_zenith + diffuse_solar
  reflected_absorbed <- 0.5 * albedo * total_horizontal * GLOBE_EMISSIVITY

  total_solar_absorbed <- direct_absorbed + diffuse_absorbed + reflected_absorbed

  # Sky temperature (approximation)
  T_sky <- temp - 20

  # Ground temperature (approximation - slightly warmer than air in sun)
  T_ground <- temp + 5 * (total_horizontal / 1000)

  # Newton's method to solve globe temperature
  T_globe <- temp + 5

  for (iter in 1:50) {
    T_globe_K <- T_globe + ABS_ZERO
    T_air_K <- temp + ABS_ZERO
    T_sky_K <- T_sky + ABS_ZERO
    T_ground_K <- T_ground + ABS_ZERO

    # Radiative exchange with sky and ground
    rad_sky <- 0.5 * GLOBE_EMISSIVITY * STEFAN_BOLTZMANN * (T_sky_K^4 - T_globe_K^4)
    rad_ground <- 0.5 * GLOBE_EMISSIVITY * STEFAN_BOLTZMANN * (T_ground_K^4 - T_globe_K^4)

    # Heat balance: solar + rad_exchange + convection = 0
    Q_solar <- total_solar_absorbed
    Q_rad <- rad_sky + rad_ground
    Q_conv <- hc_globe * (temp - T_globe)

    balance <- Q_solar + Q_rad + Q_conv

    # Derivative for Newton's method
    drad_dT <- -4 * GLOBE_EMISSIVITY * STEFAN_BOLTZMANN * T_globe_K^3
    dconv_dT <- -hc_globe
    dbalance_dT <- drad_dT + dconv_dT

    # Newton step
    if (abs(dbalance_dT) > 1e-10) {
      delta <- -balance / dbalance_dT
      delta <- max(-5, min(5, delta))
      T_globe <- T_globe + delta

      if (abs(delta) < 0.01) break
    } else {
      break
    }
  }

  # Sanity check - globe temp should be reasonable
  T_globe <- max(temp - 5, min(temp + 30, T_globe))

  return(T_globe)
}

#' Calculate globe temperature (vectorised)
#'
#' Vectorised wrapper around `solve_globe_temp()` using [purrr::pmap_dbl()].
#'
#' @param temp Air temperature in degrees Celsius.
#' @param wind_speed Wind speed in m/s.
#' @param direct_solar Direct beam solar radiation in W/m^2.
#' @param diffuse_solar Diffuse sky solar radiation in W/m^2.
#' @param zenith Solar zenith angle in degrees.
#' @param albedo Ground albedo (0--1). Default 0.12 for asphalt.
#' @param verbose If `TRUE`, pass observation index to solver for warnings.
#'
#' @return Numeric vector of globe temperatures in degrees Celsius.
#' @keywords internal
calculate_globe_temp <- function(temp, wind_speed, direct_solar, diffuse_solar,
                                 zenith, albedo = 0.12, verbose = FALSE) {

  n_obs <- length(temp)

  pmap_dbl(
    list(temp, wind_speed, direct_solar, diffuse_solar, zenith, seq_len(n_obs)),
    function(t, ws, ds, df, z, idx) {
      if (is.na(t) || is.na(ws) || is.na(ds) || is.na(df) || is.na(z)) {
        return(NA_real_)
      }
      solve_globe_temp(t, ws, ds, df, z, albedo, index = if (verbose) idx else NULL)
    }
  )
}
