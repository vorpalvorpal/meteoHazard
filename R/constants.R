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
