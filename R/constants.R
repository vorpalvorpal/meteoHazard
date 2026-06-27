#' Physical and physiological constants for TWL calculation
#'
#' A list of shared physical constants used by the globe temperature,
#' natural wet bulb, and TWL heat-balance solvers.  Each function accesses
#' these via `TWL_CONSTANTS$<NAME>` rather than duplicating the numeric
#' literal.
#'
#' \describe{
#'   \item{STEFAN_BOLTZMANN}{Stefan-Boltzmann constant, 5.67e-8 W/(m^2·K^4).}
#'   \item{ABS_ZERO}{Celsius-to-Kelvin offset, 273.15 K.}
#'   \item{GLOBE_DIAMETER}{Standard black-globe diameter, 0.15 m (150 mm).}
#'   \item{GLOBE_EMISSIVITY}{Black-globe emissivity, 0.95.}
#'   \item{VIEW_EMISSIVITY_FACTOR}{Combined view factor (0.8) × emissivity
#'     (0.95) for the natural wet bulb cylinder, 0.76.}
#'   \item{BULB_DIA}{Standard wick/bulb diameter for natural wet bulb,
#'     0.004 m.}
#'   \item{AIR_THERMAL_CONDUCTIVITY}{Thermal conductivity of air,
#'     0.028 W/(m·K). **Shared by the globe-temperature and natural-wet-bulb
#'     solvers** (Brake 2001, Whillier Table 1 / BB p. 470).}
#'   \item{AIR_KINEMATIC_VISCOSITY}{Kinematic viscosity of air at ~300 K,
#'     1.5e-5 m^2/s. Used in the Reynolds-number calculation for both the
#'     black-globe and natural-wet-bulb convective heat transfer
#'     coefficients.}
#'   \item{LATENT_HEAT_TWL_KJ}{Latent heat of evaporation of sweat **at skin
#'     temperature** (~30 °C), 2430 kJ/kg. Used in the TWL heat balance
#'     (Brake & Bates 2002, ASHRAE Eq. 14 at 30 °C). Note: the natural wet
#'     bulb solver uses 2455 kJ/kg (value near 20 °C) in a different
#'     psychrometric context — see [solve_natural_wb_single()].}
#'   \item{TWL_FLOOR}{Lower clamp on the returned TWL, 60 W/m^2. Applied
#'     both when the bisection bracket is degenerate and when the converged
#'     metabolic rate falls below physiologically plausible levels.}
#'   \item{TWL_CEILING}{Upper clamp on the returned TWL, 380 W/m^2.
#'     Corresponds to the upper bound of the published Brake & Bates (2002)
#'     TWL regime chart.}
#' }
#' @keywords internal
TWL_CONSTANTS <- list(
  STEFAN_BOLTZMANN          = 5.67e-8, # W/(m^2·K^4)
  ABS_ZERO                  = 273.15, # K
  GLOBE_DIAMETER            = 0.15, # m
  GLOBE_EMISSIVITY          = 0.95,
  VIEW_EMISSIVITY_FACTOR    = 0.76, # 0.8 * 0.95
  BULB_DIA                  = 0.004, # m
  AIR_THERMAL_CONDUCTIVITY  = 0.028, # W/(m·K); shared globe + wick
  AIR_KINEMATIC_VISCOSITY   = 1.5e-5, # m^2/s at ~300 K; globe + wick Re
  LATENT_HEAT_TWL_KJ        = 2430, # kJ/kg at skin temp ~30 °C
  TWL_FLOOR                 = 60, # W/m^2, lower clamp on TWL
  TWL_CEILING               = 380 # W/m^2, upper clamp on TWL
)

#' Constants for the odour hazard and exposure model
#'
#' Shared parameters for [odour_hazard()] and [odour_exposure()]. The odour
#' model is a two-layer screening tool: a receptor-independent **hazard**
#' (source emission divided by atmospheric ventilation) and a geometry-aware
#' **exposure** layer (Gaussian-plume distance and direction). These constants
#' carry the tunable, partly-uncalibrated parameters; see the calibration issue
#' referenced in `NEWS.md`.
#'
#' \describe{
#'   \item{U_CALM_FLOOR}{Lower clamp on wind speed (m/s) used in the `1/u`
#'     advective-dilution term, 0.5. Below this the anemometer is near its
#'     noise floor and Gaussian-plume theory is invalid.}
#'   \item{SIGMA_FC_DEG}{Forecast wind-direction uncertainty (degrees), 12.
#'     Convolved with the physical plume half-width to widen the directional
#'     Gaussian (separate from, not a substitute for, the stability-dependent
#'     physical width).}
#'   \item{V_MOD_MAX}{Maximum surface-volatilisation contribution to the
#'     generation modifier `G`, 0.30. Widened from an earlier 0.15; high
#'     calibration uncertainty (Henry's-law doubling per ~10 C; Zou 2003).}
#'   \item{PM_MIN, PM_MAX}{Peak-to-mean ratio bounds, 1.0 (unstable) to 3.0
#'     (very stable). Accounts for sub-minute concentration peaks that drive
#'     odour annoyance but are averaged out of ~hourly Gaussian sigmas.
#'     Conservative default; high calibration uncertainty.}
#'   \item{H_MIX_FALLBACK_STABLE, H_MIX_FALLBACK_UNSTABLE}{Mixing-depth
#'     fallbacks (m) when `boundary_layer_height` is `NA`: 200 (stable/calm,
#'     shallow nocturnal boundary layer) and 600 (neutral/unstable).}
#'   \item{SIGMA_Y_COEF}{Briggs (1973) rural lateral-dispersion leading
#'     coefficients for Pasquill-Gifford classes A-F (indices 1-6), used as
#'     `sigma_y = c_y * x / sqrt(1 + 0.0001 * x)`.}
#'   \item{SHELTER_OPEN_REF}{Sky-view openness angle (degrees) for flat/open
#'     terrain with no sheltering effect, 85. Above this value `s_f = 0`.}
#'   \item{SHELTER_ENCLOSED_REF}{Sky-view openness angle (degrees) for strongly
#'     enclosed terrain at full shelter effect, 50. Below this value `s_f = 1`.}
#'   \item{SHELTER_U_FULL}{Wind speed (m/s) at or below which the shelter
#'     regime weight `w_r = 1` (full effect), 1.5.}
#'   \item{SHELTER_U_FLUSH}{Wind speed (m/s) at or above which `w_r = 0`
#'     (shelter flushed out by wind), 6.0.}
#'   \item{SHELTER_MAX_REDUCTION}{Maximum fractional reduction in `u_eff` from
#'     valley sheltering, 0.7. Uncalibrated screening default.}
#'   \item{DRAINAGE_SHELTER_OVERLAP}{Mutual-exclusion weight between M1
#'     drainage confinement and M3 valley sheltering, 1.0. A value of 1 means
#'     full suppression of M3 on hours where M1 drainage is active.}
#'   \item{RIM_LIFT_COEF}{Coefficient (dimensionless) scaling cumulative CBL growth
#'     into the morning vented-layer depth: \code{h_vent = pool_top + RIM_LIFT_COEF *
#'     cbl_cumsum}. Uncalibrated screening default 0.2; calibration deferred to
#'     issue #8.}
#'   \item{RIM_DELTA}{Logistic half-width (m) for the vertical reach gate: controls
#'     sharpness of the transition from 0 to 1 as \code{h_vent} crosses the receptor
#'     height \code{z_j}. Uncalibrated screening default 25 m; calibration → #8.}
#' }
#' @keywords internal
ODOUR_CONSTANTS <- list(
  U_CALM_FLOOR             = 0.5, # m/s, lower clamp on wind in 1/u
  SIGMA_FC_DEG             = 12, # deg, forecast wind-direction uncertainty
  V_MOD_MAX                = 0.30, # max surface-volatilisation bump in G
  PM_MIN                   = 1.0, # peak-to-mean at unstable (s = 0)
  PM_MAX                   = 3.0, # peak-to-mean at very stable (s = 5)
  H_MIX_FALLBACK_STABLE    = 200, # m, NA boundary-layer fallback when calm
  H_MIX_FALLBACK_UNSTABLE  = 600, # m, NA boundary-layer fallback otherwise
  SIGMA_Y_COEF             = c(0.22, 0.16, 0.11, 0.08, 0.06, 0.04), # A-F
  BRUNT_A                  = 0.34,   # Brunt (1932) longwave coefficient a
  BRUNT_B                  = 0.14,   # Brunt (1932) longwave coefficient b
  BRUNT_CLOUD_K            = 0.1,    # cloud-cover reduction factor k (calibratable)
  VENKATRAM_COEF           = 2400,   # h_floor = 2400 * u_star^1.5 (Venkatram 1980)
  POOL_Z0                  = 0.1,    # roughness length z0 for u* (m)
  POOL_H_SAT               = 300,    # max cold-pool depth (m, Whiteman 1999)
  POOL_Q_SAT               = 3e6,    # saturation heat deficit scale (J/m²)
  POOL_KAPPA               = 0.4,    # von Karman constant
  ISC3_SIGMA_Y0_COEF       = 4.3,   # sigma_y0 = crosswind_halfwidth / 4.3
  ISC3_SIGMA_Z0_COEF       = 2.15,  # sigma_z0 = emit_extent / 2.15
  # C3b terrain morning-pulse pathway constants
  CONFINEMENT_1A           = 0.3,   # night confinement factor for pathway 1a
  VENTING_1A               = 0.8,   # morning venting boost for pathway 1a
  FUMIC_1B                 = 0.8,   # fumigation directional factor for pathway 1b
  DELTA_FLOOR              = 20,    # minimum delta for smooth pool partition (m)
  DELTA_FRAC               = 0.25,  # delta = max(DELTA_FLOOR, DELTA_FRAC * pool_top)
  # C6 M3 valley sheltering (uncalibrated screening defaults; calibration → #8)
  SHELTER_OPEN_REF          = 85,   # deg openness: flat/open (no shelter)
  SHELTER_ENCLOSED_REF      = 50,   # deg openness: strongly enclosed
  SHELTER_U_FULL            = 1.5,  # m/s: at/below this, full regime weight
  SHELTER_U_FLUSH           = 6.0,  # m/s: at/above this, shelter flushed to 0
  SHELTER_MAX_REDUCTION     = 0.7,  # maximum fractional u_eff reduction
  DRAINAGE_SHELTER_OVERLAP  = 1.0,  # 1 = full mutual exclusion with M1 drainage
  # C9 — reference distance for the exposure normaliser (Briggs class-F worst case)
  X_REF_EXPOSURE            = 250,  # m
  # C8 upslope rim-venting constants (uncalibrated screening defaults; calibration → #8)
  RIM_LIFT_COEF             = 0.2,  # α: pool_top + α·cbl_cumsum = h_vent (m/m)
  RIM_DELTA                 = 25    # δ: logistic reach sharpness (m)
)
