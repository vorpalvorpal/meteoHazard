#' Three-zone evaporation model
#'
#' Computes actual evaporation rate E (W/m²) from sweat supply and evaporative
#' capacity, implementing the three-zone model from Brake & Bates (2002),
#' Appendix 1 / Figure 3.
#'
#' Zone A (x < 0.46):   no dripping, E = lambda_Sr
#' Zone B (0.46-1.7):   partial dripping, E = lambda_Sr * exp(-0.4127*(x-0.46)^1.168)
#' Zone C (x > 1.7):    fully wet skin, E = Emax
#'
#' Note: the published Zone B formula has "1.8x-0.46" but the correct
#' continuous form (matching boundaries at x=0.46 and x=1.7) is "x-0.46".
#'
#' @param lambda_Sr Sweat evaporation supply in W/m^2 (= lambda * Sr / 3.6).
#' @param Emax Maximum evaporative capacity from fully-wet skin in W/m^2.
#'
#' @return Actual evaporation rate E in W/m^2.
#' @keywords internal
three_zone_evaporation <- function(lambda_Sr, Emax) {
  if (Emax <= 0) {
    return(0)
  }
  x <- lambda_Sr / Emax
  E <- if (x < 0.46) {
    lambda_Sr
  } else if (x <= 1.7) {
    lambda_Sr * exp(-0.4127 * (x - 0.46)^1.168)
  } else {
    Emax
  }
  max(0, E)
}


#' Solve TWL for a single observation
#'
#' Iteratively solves the heat balance equations to find the maximum
#' sustainable metabolic rate using the three-zone evaporation model
#' from Brake & Bates (2002), Appendix 1.
#'
#' The algorithm iterates on mean skin temperature (t_skin) to find the
#' unique value at which heat flow from core to skin (H) equals heat
#' loss from skin to environment (C + R + E).  At that steady-state,
#' the metabolic rate M = H + B (respiratory losses) equals the TWL.
#'
#' **Bracketing assumption**: the heat-balance function
#' `balance(t_skin) = H - (CR + E)` is assumed to be monotone on the
#' interval `[temp_dewpoint + 1, t_core - 0.5]`, so that a unique root
#' exists and bisection is valid.  When this assumption is violated (i.e.
#' `balance` has the same sign at both endpoints), a warning is emitted
#' naming the observation index and the solver falls back to the
#' best-|balance| result from the bisection loop.
#'
#' @param temp Air temperature in degrees Celsius.
#' @param wind_speed Wind speed in m/s.
#' @param RH Relative humidity in percent.
#' @param direct_solar Direct solar radiation in W/m^2 (not used in heat
#'   balance; retained for interface compatibility).
#' @param diffuse_solar Diffuse solar radiation in W/m^2 (not used in heat
#'   balance; retained for interface compatibility).
#' @param pressure Atmospheric pressure in kPa.
#' @param pa Actual vapour pressure in kPa.
#' @param temp_dewpoint Dew point temperature in degrees Celsius.
#' @param wet_bulb Natural wet bulb temperature in degrees Celsius (not used
#'   in heat balance; retained for interface compatibility).
#' @param globe_temp Globe temperature in degrees Celsius (not used in heat
#'   balance; retained for interface compatibility).
#' @param trad Mean radiant temperature in degrees Celsius.
#' @param max_core_temp Maximum core temperature in degrees Celsius.
#' @param max_sweat_rate Maximum sweat rate in kg/(m^2.hr).
#' @param Icl Clothing insulation in clo.
#' @param icl Clothing vapour permeability (0--1).
#' @param LR Lewis relation in K/kPa.
#' @param lambda Latent heat of evaporation of sweat at skin temperature,
#'   in kJ/kg. Brake & Bates (2002) use 2430 kJ/kg (ASHRAE Eq. 14 at 30°C).
#'   Note: `natural_wet_bulb.R` uses 2455 kJ/kg (value near 20°C) for the
#'   wet-bulb psychrometric equation — a different physical context.
#' @param index Observation index for warnings (optional).
#'
#' @return TWL in W/m^2.
#' @keywords internal
solve_twl_single <- function(temp, wind_speed, RH, direct_solar, diffuse_solar,
                             pressure, pa, temp_dewpoint,
                             wet_bulb, globe_temp, trad,
                             max_core_temp, max_sweat_rate,
                             Icl, icl, LR, lambda,
                             index = NULL) {
  # ------------------------------------------------------------------
  # Constants
  # ------------------------------------------------------------------
  t_core <- max_core_temp # limiting deep body core temperature (°C)

  # ------------------------------------------------------------------
  # Clothing parameters  [ASHRAE Fundamentals Ch. 8]
  # ------------------------------------------------------------------
  Rcl <- 0.155 * Icl # thermal resistance, (m²K)/W  [ASHRAE Eq. 41]
  fcl <- 1.0 + 0.31 * Icl # clothing area factor          [ASHRAE Eq. 47]
  Recl <- Rcl / (LR * icl) # evaporative resistance        [ASHRAE Table II]

  # ------------------------------------------------------------------
  # Convective heat transfer coefficient  [EESAM Eq. 13]
  # hc = 0.608 * P^0.6 * V^0.6
  # ------------------------------------------------------------------
  hc_paper <- 0.608 * pressure^0.6 * wind_speed^0.6 # W/(m²K)

  # ------------------------------------------------------------------
  # Evaporative heat transfer coefficient  [EESAM Eq. 17]
  # he = 1587 * hc * P / (P - pa)^2
  # (the LR*hc approximation is only accurate when pa << P; the exact
  # form is needed for low-pressure or high-humidity conditions)
  # ------------------------------------------------------------------
  he <- 1587 * hc_paper * pressure / (pressure - pa)^2 # W/(m²·kPa)

  # ------------------------------------------------------------------
  # Per-t_skin heat-balance closure
  #
  # All loop-invariant quantities (hc_paper, he, fcl, Rcl, Recl,
  # t_core, trad, temp, pa, lambda, max_sweat_rate) are captured from
  # the enclosing scope.  Only the six quantities that change with
  # t_skin are computed here.
  # ------------------------------------------------------------------
  compute_components <- function(t_skin) {
    # Radiant heat transfer coefficient — depends on iterated t_skin
    # hr = 4.61 * [1 + (trad + t_skin)/546]^3   [EESAM Eq. 9]
    hr <- 4.61 * (1 + (trad + t_skin) / 546)^3 # W/(m²K)
    h <- hr + hc_paper # [ASHRAE Eq. 9]

    # Operative temperature  [ASHRAE Eq. 8]
    toper <- (hr * trad + hc_paper * temp) / h

    # Clothing correction factors
    # Fcle: sensible heat  [ASHRAE Table II, Sensible Heat Flow, last eq.]
    # Fpcl: evaporative   [ASHRAE Table II, Evaporative Heat Flow, last eq.]
    Fcle <- fcl / (1 + fcl * h * Rcl)
    Fpcl <- 1 / (1 + fcl * he * Recl)

    # Sensible heat loss from skin: C + R  [ASHRAE Table III, 3rd eq.]
    CR <- Fcle * h * (t_skin - toper)

    # Maximum evaporative capacity from fully-wet skin (w = 1)
    # ps = sat. vapour pressure at skin temperature
    # Esk = w * Fpcl * fcl * he * (ps - pa)  [ASHRAE Table III Evap., 3rd]
    ps <- calc_sat_vp(t_skin)
    Emax <- max(0, Fpcl * fcl * he * (ps - pa))

    # Thermoregulatory signal  [Cabanac]
    # tz = 0.1 * t_skin + 0.9 * t_core
    tz <- 0.1 * t_skin + 0.9 * t_core

    # Physiological conductance from Wyndham's data  [Figure 1 fit]
    # Kcs = 84 + 72 * tanh(1.3 * (tz - 37.9))   W/(m²K)
    Kcs <- 84 + 72 * tanh(1.3 * (tz - 37.9))

    # Heat flow from core to skin  [EESAM Eq. 23]
    # H = Kcs * (t_core - t_skin)
    H <- max(0, Kcs * (t_core - t_skin))

    # Sweat rate from Wyndham's data  [Figure 2 fit]
    # Sr = 0.42 + 0.44 * tanh(1.16 * (tz - 37.4))   kg/(m²hr)
    Sr <- max(0, min(max_sweat_rate, 0.42 + 0.44 * tanh(1.16 * (tz - 37.4))))
    lambda_Sr <- lambda * Sr / 3.6 # convert kJ/kg * kg/(m²hr) → W/m²

    E <- three_zone_evaporation(lambda_Sr, Emax)

    # Heat balance: H should equal C+R+E at steady state
    # balance > 0 means core is supplying more heat than skin loses
    #   → skin temperature should rise
    # balance < 0 means skin losing more than core supplies
    #   → skin temperature should fall
    balance <- H - (CR + E)

    list(balance = balance, H = H, CR = CR, E = E)
  }

  # ------------------------------------------------------------------
  # Iterative solution: find t_skin where H = C + R + E
  # ------------------------------------------------------------------
  # Starting bracket: dew point + 1 °C  to  core - 0.5 °C
  t_skin_min <- temp_dewpoint + 1.0
  t_skin_max <- t_core - 0.5

  # Guard: if dew point is so high the bracket is degenerate (extreme
  # humidity approaching skin temperature), warn and return the floor.
  if (t_skin_min >= t_skin_max) {
    if (!is.null(index)) {
      cli_alert_warning(
        "Observation {index}: degenerate bisection bracket \\
        (dew point {round(temp_dewpoint, 1)} C >= core limit). \\
        Returning minimum TWL."
      )
    }
    return(TWL_CONSTANTS$TWL_FLOOR)
  }

  # ------------------------------------------------------------------
  # Bracket sign-check (diagnostic only — does NOT alter numerics).
  #
  # Evaluate balance at both bracket endpoints to verify the root is
  # bracketed.  These evaluations feed only the warning; they do NOT
  # update best_balance / best_t_skin and the iteration path is
  # unchanged so results remain bit-identical with the previous code.
  # ------------------------------------------------------------------
  if (!is.null(index)) {
    bal_lo <- compute_components(t_skin_min)$balance
    bal_hi <- compute_components(t_skin_max)$balance
    if (sign(bal_lo) == sign(bal_hi)) {
      cli_alert_warning(
        "Observation {index}: heat-balance root not bracketed \\
        (balance at t_skin_min = {round(bal_lo, 2)}, \\
        balance at t_skin_max = {round(bal_hi, 2)} — same sign). \\
        Bisection may not converge; returning best-residual estimate."
      )
    }
  }

  t_skin <- 0.5 * (t_skin_min + t_skin_max)

  best_t_skin <- t_skin
  best_balance <- Inf

  MAX_ITER <- 200
  TOL_CONVERGENCE <- 0.01 # W/m²

  for (iter in 1:MAX_ITER) {
    comp <- compute_components(t_skin)
    balance <- comp$balance

    if (abs(balance) < best_balance) {
      best_balance <- abs(balance)
      best_t_skin <- t_skin
    }

    if (abs(balance) < TOL_CONVERGENCE) break

    # ----------------------------------------------------------------
    # Bisection update: keep bracket tight
    # ----------------------------------------------------------------
    if (balance > 0) {
      t_skin_min <- t_skin
    } else {
      t_skin_max <- t_skin
    }
    t_skin <- 0.5 * (t_skin_min + t_skin_max)

    # Safety: clamp bracket
    t_skin_min <- max(t_skin_min, temp_dewpoint + 1.0)
    t_skin_max <- min(t_skin_max, t_core - 0.1)
    if (t_skin_min >= t_skin_max) break
  }

  # ------------------------------------------------------------------
  # Compute M at the converged skin temperature
  # M = H + B  where  B = 0.0014*M*(34-ta) + 0.0173*M*(5.87-pa)
  # Solving: M*(1 - resp_coeff) = H  →  M = H/(1-resp_coeff)
  # [ASHRAE Eq. 26]
  # ------------------------------------------------------------------
  # Re-evaluate at best_t_skin using the same closure (eliminates the
  # ~25-line duplicated recomputation block that was here previously).
  comp_final <- compute_components(best_t_skin)
  H_f <- comp_final$H
  CR_f <- comp_final$CR
  E_f <- comp_final$E

  # Use average of H and (CR+E) as our best estimate of heat throughput
  # to avoid systematic bias from imperfect convergence
  H_use <- 0.5 * (H_f + CR_f + E_f)

  # Respiratory coefficient
  resp_coeff <- 0.0014 * (34 - temp) + 0.0173 * (5.87 - pa)
  resp_coeff <- max(0, min(0.5, resp_coeff)) # clamp to physically valid range

  M <- H_use / (1 - resp_coeff)

  return(max(TWL_CONSTANTS$TWL_FLOOR, M))
}
