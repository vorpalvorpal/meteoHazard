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
#' **exposure** layer (Gaussian-plume distance and direction). The generation
#' modifier `G` combines its five fractional terms **multiplicatively**
#' (independent fractional effects compound; v3 change) and these constants
#' carry the tunable, partly-uncalibrated parameters; see the calibration
#' issue referenced in `NEWS.md`.
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
#'     generation modifier `G`, 0.30 (approached asymptotically at
#'     `V_MOD_T_HI` by the exponential Henry's-law ramp). High calibration
#'     uncertainty (Henry's-law doubling per ~10 C; Zou 2003).}
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
#'   \item{ODORANT_SOLUBILITY_DEFAULT}{Default value of the `odorant_solubility`
#'     argument, 0.5. 0 = poorly soluble (reduced-sulfur compounds, the typical
#'     landfill-odour driver); 1 = highly soluble (ammonia, amines, soluble
#'     VOCs); 0.5 = mixed profile.}
#'   \item{W_RAIN_RATE_LIGHT, W_RAIN_RATE_MOD, W_RAIN_RATE_HEAVY}{Precipitation
#'     rate tiers (mm/h) for below-cloud scavenging: 0.2, 1.0, 4.0.}
#'   \item{W_RAIN_FACTOR_LIGHT, W_RAIN_FACTOR_MOD, W_RAIN_FACTOR_HEAVY}{Soluble-limit
#'     washout survival factors (Seinfeld & Pandis Ch. 19) at the light/moderate/heavy
#'     precipitation tiers: 0.40, 0.15, 0.05. These are the `odorant_solubility = 1`
#'     endpoint; `W_rain` blends toward 1.0 (no washout) as solubility falls.}
#'   \item{V_MOD_T_LO}{Temperature (deg C) at/below which surface NMOC
#'     volatilisation is negligible (`V_mod = 0`), 10.}
#'   \item{V_MOD_T_HI}{Temperature (deg C) at/above which `V_mod` is clamped
#'     to `V_MOD_MAX`, 35 -- a deliberate screening ceiling (this site's worst
#'     case is winter inversions, not extreme heat).}
#'   \item{V_MOD_DOUBLING_C}{Temperature interval (deg C) over which the
#'     exponential `V_mod` curve doubles, 10 -- Henry's-law / Clausius-Clapeyron
#'     vapour pressure roughly doubles per 10 degC.}
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
  RIM_DELTA                 = 25,   # δ: logistic reach sharpness (m)
  # Solubility-aware below-cloud washout (v3 Change 4)
  ODORANT_SOLUBILITY_DEFAULT = 0.5,   # 0 = poorly soluble (reduced sulfur), 1 = highly soluble
  W_RAIN_RATE_LIGHT   = 0.2,  W_RAIN_RATE_MOD   = 1.0,  W_RAIN_RATE_HEAVY   = 4.0,   # mm/h
  W_RAIN_FACTOR_LIGHT = 0.40, W_RAIN_FACTOR_MOD = 0.15, W_RAIN_FACTOR_HEAVY = 0.05,  # soluble-limit
  # Exponential V_mod (v3 Change 3)
  V_MOD_T_LO       = 10,   # deg C, below which volatilisation negligible
  V_MOD_T_HI       = 35,   # deg C, screening ceiling
  V_MOD_DOUBLING_C = 10    # deg C per doubling (Henry / Clausius-Clapeyron)
)

#' Constants for the dust hazard model
#'
#' Shared physical parameters for [dust_flux()], [dust_hazard()], and the
#' internal `.dust_crust_factor()` helper. The dust model is a physical
#' saltation-to-emission chain: a Shao & Lu (2000) threshold friction
#' velocity, a Fecan et al. (1999) soil-moisture correction, an Owen (1964) /
#' White (1979) saltation flux, and a Marticorena & Bergametti (1995, "MB95")
#' sandblasting efficiency.
#'
#' \describe{
#'   \item{A_N}{Shao & Lu (2000) dimensionless threshold coefficient, 0.0123.}
#'   \item{GAMMA}{Shao & Lu (2000) interparticle cohesion parameter (N/m),
#'     3.0e-4 (their tabulated range is 1.65e-4-5e-4; 3e-4 is their
#'     recommended best estimate).}
#'   \item{RHO_P}{Particle density (kg/m^3), 2650 (quartz).}
#'   \item{RHO_A_REF}{Reference air density (kg/m^3), 1.225 (sea-level,
#'     15 degC standard atmosphere; not yet temperature/pressure adjusted,
#'     see T10 TODO).}
#'   \item{G}{Gravitational acceleration, 9.81 m/s^2.}
#'   \item{KAPPA}{von Karman constant, 0.40.}
#'   \item{Z_REF}{Reference height (m) for the 10 m wind used in the
#'     logarithmic wind-profile law, 10.}
#'   \item{Z0_SMOOTH_RATIO}{Smooth-bed aerodynamic roughness as a fraction of
#'     the modal grain diameter, `z0 = d * Z0_SMOOTH_RATIO`, 1/30 (Nikuradse
#'     equivalent-sand roughness; MB95 z0s convention).}
#'   \item{FECAN_WP_QUAD, FECAN_WP_LIN}{Fecan et al. (1999) quadratic and
#'     linear coefficients of the residual gravimetric moisture w' (as a
#'     function of clay %), 0.0014 and 0.17.}
#'   \item{FECAN_A}{Fecan et al. (1999) moisture-correction coefficient A in
#'     `sqrt(1 + A * (w - w')^b)`, 1.21.}
#'   \item{FECAN_B}{Fecan et al. (1999) moisture-correction exponent b',
#'     0.68.}
#'   \item{MB95_ALPHA_SLOPE}{Marticorena & Bergametti (1995) sandblasting
#'     log-alpha slope on clay % (cm^-1), 0.134.}
#'   \item{MB95_ALPHA_INTERCEPT}{Marticorena & Bergametti (1995) sandblasting
#'     log-alpha intercept, -6.}
#'   \item{MB95_CLAY_CAP}{Clay fraction (%) above which the MB95 alpha fit is
#'     no longer validated, 20 (Foroutan et al. 2017; Kok et al. 2014).}
#' }
#' @keywords internal
DUST_CONSTANTS <- list(
  A_N                   = 0.0123,   # Shao & Lu (2000) threshold coefficient
  GAMMA                 = 3.0e-4,   # N/m; Shao & Lu (2000) interparticle cohesion (best ~3e-4 of 1.65e-4-5e-4)
  RHO_P                 = 2650,     # kg/m^3 particle density (quartz)
  RHO_A_REF             = 1.225,    # kg/m^3 reference air density
  G                     = 9.81,     # m/s^2
  KAPPA                 = 0.40,     # von Karman constant
  Z_REF                 = 10,       # m, wind reference height
  Z0_SMOOTH_RATIO       = 1 / 30,   # z0 = d/30 smooth-bed roughness (Nikuradse; MB95 z0s)
  FECAN_WP_QUAD         = 0.0014,   # Fecan (1999) w' = a*clay^2 + b*clay
  FECAN_WP_LIN          = 0.17,
  FECAN_A               = 1.21,     # Fecan (1999) correction coefficient A
  FECAN_B               = 0.68,     # Fecan (1999) correction exponent b'
  MB95_ALPHA_SLOPE      = 0.134,    # MB95 sandblasting log-alpha slope (cm^-1)
  MB95_ALPHA_INTERCEPT  = -6,       # MB95 sandblasting log-alpha intercept
  MB95_CLAY_CAP         = 20        # % clay; MB95 alpha validity ceiling
)
