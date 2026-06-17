#' Solve natural wet bulb temperature for a single observation
#'
#' Iteratively solves the heat balance on an unventilated wet bulb
#' thermometer using Newton's method, based on the Brake (2001) methodology.
#'
#' @param temp Air temperature in degrees Celsius.
#' @param RH Relative humidity in percent.
#' @param pressure Atmospheric pressure in kPa.
#' @param wind_speed Wind speed in m/s.
#' @param globe_temp Globe temperature in degrees Celsius.
#' @param index Observation index for warnings (optional).
#'
#' @return Natural wet bulb temperature in degrees Celsius.
#' @keywords internal
solve_natural_wb_single <- function(temp, RH, pressure, wind_speed, globe_temp,
                                    index = NULL) {
  # Physical constants (shared via TWL_CONSTANTS)
  VIEW_EMISSIVITY_FACTOR <- TWL_CONSTANTS$VIEW_EMISSIVITY_FACTOR
  STEFAN_BOLTZMANN <- TWL_CONSTANTS$STEFAN_BOLTZMANN
  ABS_ZERO <- TWL_CONSTANTS$ABS_ZERO
  BULB_DIA <- TWL_CONSTANTS$BULB_DIA
  AIR_THERMAL_CONDUCTIVITY <- TWL_CONSTANTS$AIR_THERMAL_CONDUCTIVITY

  # Latent heat of evaporation in J/kg at ~20 °C (psychrometric wet bulb
  # context).  Intentionally different from TWL_CONSTANTS$LATENT_HEAT_TWL_KJ
  # (2430 kJ/kg at 30 °C skin temperature) — see ?TWL_CONSTANTS.
  LATENT_HEAT_EVAP <- 2455000 # J/kg at ~20 °C

  ACCURACY_REQUIRED <- 0.02
  MAX_ITERATIONS <- 100

  # Handle NA inputs
  if (is.na(temp) || is.na(RH) || is.na(pressure) ||
    is.na(wind_speed) || is.na(globe_temp)) {
    return(NA_real_)
  }

  # Ensure minimum wind speed
  wind_speed <- max(wind_speed, 0.1)

  # Convert pressure to Pa
  pressure_pa <- pressure * 1000

  # Calculate aspirated wet bulb temperature (WBa) using Stull (2011)
  WBa <- calc_aspirated_wb(temp, RH)

  # Check for invalid WBa
  if (WBa > temp) {
    WBa <- temp - 0.5
  }

  # For very small wet bulb depression, natural ~ aspirated
  if ((temp - WBa) < 0.3) {
    return(WBa)
  }

  # Convective heat transfer coefficient for cylinder
  Re <- wind_speed * BULB_DIA / TWL_CONSTANTS$AIR_KINEMATIC_VISCOSITY
  Nu <- 0.281 * Re^0.6
  hc <- Nu * AIR_THERMAL_CONDUCTIVITY / BULB_DIA

  # Mean radiant temperature from globe temp
  T_mrt <- globe_temp

  # Initial guess - between aspirated WB and air temp
  WBn <- WBa + 0.3 * (temp - WBa)

  # Track best solution
  best_WBn <- WBn
  best_error <- Inf
  prev_WBn <- WBn
  oscillation_count <- 0

  for (iter in 1:MAX_ITERATIONS) {
    # Saturation vapour pressure at wet bulb
    es_wb <- calc_sat_vp(WBn)

    # Actual vapour pressure
    es_air <- calc_sat_vp(temp)
    ea <- (RH / 100) * es_air

    # Evaporative heat loss
    # cp = 1005 J/(kg.K); pressure_pa in Pa; es/ea in kPa converted to Pa by *1000
    E <- (0.622 * LATENT_HEAT_EVAP * hc / (pressure_pa * 1005)) * (es_wb * 1000 - ea * 1000)

    # Convective heat gain
    C <- hc * (temp - WBn)

    # Radiative heat gain (linearised)
    T_mrt_K <- T_mrt + ABS_ZERO
    WBn_K <- WBn + ABS_ZERO
    R <- VIEW_EMISSIVITY_FACTOR * STEFAN_BOLTZMANN * (T_mrt_K^4 - WBn_K^4)

    # Heat balance: C + R - E = 0
    heat_balance <- C + R - E

    # Track best solution
    if (abs(heat_balance) < best_error) {
      best_error <- abs(heat_balance)
      best_WBn <- WBn
    }

    # Check convergence
    if (abs(heat_balance) < ACCURACY_REQUIRED) {
      return(WBn)
    }

    # Check for oscillation
    if (abs(WBn - prev_WBn) < 0.001) {
      oscillation_count <- oscillation_count + 1
      if (oscillation_count > 5) {
        return(best_WBn)
      }
    } else {
      oscillation_count <- 0
    }
    prev_WBn <- WBn

    # Derivative for Newton's method
    des_dT <- es_wb * (18.678 - 2 * WBn / 234.5) / (257.14 + WBn) -
      es_wb * (18.678 - WBn / 234.5) * WBn / (257.14 + WBn)^2
    dE_dT <- (0.622 * LATENT_HEAT_EVAP * hc / (pressure_pa * 1005)) * des_dT * 1000
    dC_dT <- -hc
    dR_dT <- -4 * VIEW_EMISSIVITY_FACTOR * STEFAN_BOLTZMANN * WBn_K^3

    dbalance_dT <- dC_dT + dR_dT - dE_dT

    # Newton step with damping
    if (abs(dbalance_dT) > 1e-10) {
      delta <- -heat_balance / dbalance_dT
      # Adaptive damping
      damping <- ifelse(abs(delta) > 1, 0.3, 0.7)
      delta <- damping * delta
      delta <- max(-0.5, min(0.5, delta))
      WBn <- WBn + delta
    } else {
      # Small perturbation if stuck
      WBn <- WBn + 0.1
    }

    # Keep in valid range
    WBn <- max(WBa - 1, min(temp + 2, WBn))
  }

  # Return best solution if didn't converge
  return(best_WBn)
}

#' Calculate natural wet bulb temperature (vectorised)
#'
#' Vectorised wrapper around `solve_natural_wb_single()` using
#' [purrr::pmap_dbl()], with optional progress bar.
#'
#' @param temp Air temperature in degrees Celsius.
#' @param RH Relative humidity in percent.
#' @param pressure Atmospheric pressure in kPa.
#' @param wind_speed Wind speed in m/s.
#' @param globe_temp Globe temperature in degrees Celsius.
#' @param verbose If `TRUE`, show warnings.
#' @param show_progress If `TRUE`, display a progress bar for large inputs.
#'
#' @return Numeric vector of natural wet bulb temperatures in degrees Celsius.
#' @keywords internal
calculate_natural_wet_bulb <- function(temp, RH, pressure, wind_speed, globe_temp,
                                       verbose = FALSE, show_progress = FALSE) {
  n_obs <- length(temp)

  pb_id <- NULL
  if (show_progress && n_obs > 100) {
    pb_id <- cli_progress_bar("Computing natural wet bulb", total = n_obs)
  }

  wet_bulb_natural <- pmap_dbl(
    list(temp, RH, pressure, wind_speed, globe_temp, seq_len(n_obs)),
    function(t, rh, p, ws, gt, idx) {
      if (show_progress && n_obs > 100 && idx %% 10 == 0) {
        cli_progress_update(id = pb_id)
      }

      solve_natural_wb_single(t, rh, p, ws, gt,
        index = if (verbose) idx else NULL
      )
    }
  )

  if (show_progress && n_obs > 100) {
    cli_progress_done(id = pb_id)
  }

  return(wet_bulb_natural)
}
