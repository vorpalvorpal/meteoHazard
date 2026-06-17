#' Landfill Odour Dispersal Risk Score
#'
#' Computes an hourly odour risk score for each timestep in a meteorological
#' forecast, representing the worst-case impact at any of the specified receptor
#' locations downwind of a landfill site.
#'
#' The function does not query the API. The calling code is responsible for
#' fetching the required hourly variables from the Open-Meteo
#' \code{/v1/forecast} endpoint (with \code{&wind_speed_unit=ms}) and passing
#' them as a data frame. Row order is assumed to represent consecutive hourly
#' timesteps; timestamps are not used.
#'
#' @section Model architecture:
#'
#' For each timestep \eqn{t} and each receptor \eqn{j}, the risk score is:
#'
#' \deqn{R_{t,j} = W_{\mathrm{dir},j} \times W_{\mathrm{spd}} \times
#'   F_{\mathrm{disp},j} \times W_{\mathrm{rain}} \times G}
#'
#' During fumigation, \eqn{G} is replaced by \eqn{G_{\mathrm{eff}}} and
#' \eqn{W_{\mathrm{dir},j}} has a floor applied. The returned score for
#' timestep \eqn{t} is the maximum across all receptors:
#' \deqn{R_t = \max_j R_{t,j}}
#'
#' No clamping to \eqn{[0, 1]} is applied — the raw product is returned
#' directly. Individual factor bounds: \eqn{W_{\mathrm{dir}} \in [0, 1]},
#' \eqn{W_{\mathrm{spd}} \in [0.05, 1.0]},
#' \eqn{F_{\mathrm{disp}} \in [0.003, 1.0]},
#' \eqn{W_{\mathrm{rain}} \in [0.05, 1.0]},
#' \eqn{G \in [\sim 0.80, \sim 1.80]}. Theoretical maximum raw score
#' \eqn{\approx 1.80}.
#'
#' @section Factor summary:
#' \describe{
#'   \item{\eqn{W_{\mathrm{dir},j}}}{Gaussian wind direction profile \eqn{[0,
#'     1]}. \eqn{\exp(-0.5(\Delta\theta/\sigma_{\mathrm{dir}})^2)} with
#'     \eqn{\sigma_{\mathrm{dir}} = 10^\circ}. Computed per receptor.}
#'   \item{\eqn{W_{\mathrm{spd}}}}{Wind speed factor \eqn{[0.05, 1.0]}.
#'     Peak at light winds (0.5–2 m/s); declines as dilution increases above
#'     2 m/s; collapses to 0.05 above 8 m/s.}
#'   \item{\eqn{F_{\mathrm{disp},j}}}{Unified dispersion factor \eqn{[0.003,
#'     1.0]}. Replaces the separate \eqn{S_{\mathrm{stab}}}, \eqn{M_h}, and
#'     \eqn{D_j} of the previous implementation with a single
#'     physics-based Pasquill-Gifford Gaussian plume factor accounting for
#'     distance-, stability-, and BL-height-dependent dispersion geometry.
#'     Computed per receptor and per timestep.}
#'   \item{\eqn{W_{\mathrm{rain}}}}{Rainfall scavenging factor \eqn{[0.05,
#'     1.0]}. Below-cloud washout of soluble odorants during precipitation.}
#'   \item{\eqn{G}}{Generation modifier \eqn{[\sim 0.80, \sim 1.80]}.
#'     Additive structure separating bulk gas pathway from surface
#'     volatilisation:
#'     \eqn{G = (1 + \Delta P_{\mathrm{mod}} + R_{\mathrm{mod}} +
#'     S_{\mathrm{seal}} + H_{\mathrm{mod}}) + V_{\mathrm{mod}}}.}
#' }
#'
#' @section Required \code{met_data} columns:
#'
#' All columns correspond to Open-Meteo hourly forecast parameters
#' (request with \code{&wind_speed_unit=ms}). All 12 columns are required.
#' \tabular{ll}{
#'   \code{wind_direction_10m}      \tab Wind direction (°, FROM which wind blows) \cr
#'   \code{wind_speed_10m}          \tab Wind speed at 10 m (m/s) \cr
#'   \code{wind_speed_80m}          \tab Wind speed at 80 m (m/s) — for stability \cr
#'   \code{boundary_layer_height}   \tab Planetary BL height (m) \cr
#'   \code{temperature_2m}          \tab Air temperature at 2 m (°C) \cr
#'   \code{pressure_msl}            \tab Mean sea level pressure (hPa) \cr
#'   \code{precipitation}           \tab Hourly precipitation (mm) \cr
#'   \code{relative_humidity_2m}    \tab Relative humidity at 2 m (\%) \cr
#'   \code{cloud_cover}             \tab Total cloud cover (0–100 \%) \cr
#'   \code{direct_radiation}        \tab Direct solar radiation, preceding hour mean (W/m²) \cr
#'   \code{soil_moisture_0_to_1cm}  \tab Volumetric water content, 0–1 cm depth (m³/m³) \cr
#'   \code{soil_moisture_1_to_3cm}  \tab Volumetric water content, 1–3 cm depth (m³/m³)
#' }
#'
#' @param met_data A data frame (or tibble) with one row per hourly timestep,
#'   containing the 12 columns listed above. \code{NA} values are permitted
#'   and handled conservatively. Row order is assumed to represent consecutive
#'   hourly timesteps (required for the 3-hour pressure tendency and 24-hour
#'   rainfall lookback).
#'
#' @param receptors A data frame with one row per receptor location and two
#'   columns:
#'   \describe{
#'     \item{\code{bearing}}{Bearing from the landfill centroid to the
#'       receptor, in degrees from north (0–360).}
#'     \item{\code{distance}}{Distance from the landfill centroid to the
#'       receptor, in metres. Must be strictly positive.}
#'   }
#'
#' @param drainage_axes Optional data frame enabling terrain-aware drainage
#'   flow and morning fumigation parameterisation. When \code{NULL} (default),
#'   drainage and fumigation logic is skipped and generic calm handling
#'   applies. When provided, must contain:
#'   \describe{
#'     \item{\code{bearing_from}}{Bearing (0–360°) representing the direction
#'       FROM which drainage flow arrives at the landfill.}
#'     \item{\code{weight}}{Relative strength of each drainage axis (must be
#'       positive).}
#'   }
#'
#' @return A numeric vector of length \code{nrow(met_data)}. Each element is
#'   the raw risk score for that timestep (maximum across receptors). Scores
#'   are not clamped — values above 1.0 are possible under extreme conditions.
#'
#' @section Alert tier classification:
#'
#' The risk score maps to operational tiers (provisional; subject to
#' calibration against complaint records):
#' \tabular{rll}{
#'   \strong{Score} \tab \strong{Tier} \tab \strong{Response} \cr
#'   < 0.15         \tab LOW       \tab Normal operations \cr
#'   0.15–0.40      \tab MODERATE  \tab Heightened awareness; check cover integrity \cr
#'   0.40–0.80      \tab HIGH      \tab Active mitigation — reduce tipping face, deploy suppression \cr
#'   \> 0.80        \tab VERY HIGH \tab Maximum response — consider ceasing tipping
#' }
#'
#' @section Calm wind handling:
#'
#' When \code{wind_speed_10m} < 0.5 m/s (or is \code{NA}):
#' \itemize{
#'   \item \eqn{W_{\mathrm{dir},j} = 0.5} for all receptors.
#'   \item \eqn{W_{\mathrm{spd}} = 0.7}.
#'   \item Stability index \eqn{s = 4.25} (moderately-to-strongly stable).
#' }
#'
#' When \code{drainage_axes} is provided and drainage mode is active (calm,
#' dark, clear sky), \eqn{W_{\mathrm{dir},j}} and \eqn{W_{\mathrm{spd}}}
#' are overridden to terrain-specific values (see §6 drainage flow section).
#'
#' @section NA handling:
#'
#' \tabular{ll}{
#'   \code{wind_direction_10m} \code{NA} \tab \eqn{W_{\mathrm{dir},j} = 0.5} \cr
#'   \code{wind_speed_10m} \code{NA}     \tab Treat as calm (<0.5 m/s) \cr
#'   \code{wind_speed_80m} \code{NA}     \tab \eqn{s = 4.25} (calm stability) \cr
#'   \code{boundary_layer_height} \code{NA} \tab \eqn{h = 200} m if calm; \eqn{h = 500} m if not calm \cr
#'   \code{temperature_2m} \code{NA}     \tab \eqn{V_{\mathrm{mod}} = 0} \cr
#'   \code{pressure_msl} \code{NA}       \tab \eqn{\Delta P_{\mathrm{mod}} = 0} \cr
#'   \code{precipitation} \code{NA}      \tab Treat as 0 mm \cr
#'   \code{relative_humidity_2m} \code{NA} \tab \eqn{H_{\mathrm{mod}} = 0} \cr
#'   \code{cloud_cover} \code{NA}        \tab 50\% (neutral; does not block drainage) \cr
#'   \code{direct_radiation} \code{NA}   \tab 0 W/m² (allows drainage activation) \cr
#'   \code{soil_moisture_0_to_1cm} \code{NA} \tab 0 (dry) → \eqn{S_{\mathrm{seal}} = 0} \cr
#'   \code{soil_moisture_1_to_3cm} \code{NA} \tab 0 (dry) → \eqn{S_{\mathrm{seal}} = 0}
#' }
#'
#' @section F_disp physical basis:
#'
#' Ground-level concentration from a continuous ground-level area source at
#' downwind distance \eqn{x} is proportional to:
#' \deqn{C(x) \propto 1 / (u \cdot \sigma_y(x) \cdot h_{\mathrm{eff}}(x))}
#' where \eqn{h_{\mathrm{eff}}(x) = \min(\sigma_z(x), h)} is the effective
#' vertical mixing depth. \eqn{F_{\mathrm{disp}}} captures the dispersion
#' geometry component \eqn{\propto 1 / (\sigma_y \cdot h_{\mathrm{eff}})},
#' normalised so that \eqn{F_{\mathrm{disp}} = 1.0} at the reference
#' condition (class F, 250 m). \eqn{W_{\mathrm{spd}}} remains as a separate
#' heuristic wind-speed factor.
#'
#' @section Drainage flow and morning fumigation:
#'
#' When \code{drainage_axes} is provided, the model accounts for katabatic
#' drainage flow in the creek valleys surrounding the site. At night under
#' clear skies, cold air drains toward the landfill (protective). After
#' sunrise, the flow reverses (anabatic) and the overnight odour pool is
#' transported upslope toward receptors during a fumigation window.
#'
#' Drainage mode activates when: \code{direct_radiation < 10} W/m²,
#' \code{wind_speed_10m < 0.5} m/s, and \code{cloud_cover < 70}\%.
#' During drainage: \eqn{W_{\mathrm{spd}} = 0.3}; \eqn{W_{\mathrm{dir},j}}
#' is 0.05 (on drainage axis) to 0.10 (off axis).
#'
#' Fumigation activates when \code{direct_radiation > 50} W/m² and at least
#' one of the preceding 6 hours was dark. During fumigation: G is boosted by
#' \eqn{1 + 0.5 \times \mathrm{severity}} and \eqn{W_{\mathrm{dir},j}} has
#' a floor of \eqn{0.6 \times \mathrm{severity} \times \mathrm{alignment}_j}.
#' Fumigation severity scales linearly with overnight drainage or stable
#' accumulation hours (capped at 1.0 for 8+ hours), decaying over a 3-hour
#' window.
#'
#' @references
#' Briggs, G.A. (1973). \emph{Diffusion Estimation for Small Emissions}.
#' ATDL Contribution File No. 79. NOAA.
#' Basis for the σ_y and σ_z Pasquill-Gifford parameterisations.
#'
#' Turner, D.B. (1970). \emph{Workbook of Atmospheric Dispersion Estimates}.
#' US EPA AP-26.
#'
#' Gifford, F.A. (1976). Turbulent diffusion-typing schemes: A review.
#' \emph{Nuclear Safety}, 17(1), 68–86.
#'
#' Irwin, J.S. (1979). A theoretical variation of the wind profile power-law
#' exponent as a function of surface roughness and stability.
#' \emph{Atmospheric Environment}, 13(1), 191–194.
#' \doi{10.1016/0004-6981(79)90013-7}.
#'
#' Counihan, J. (1975). Adiabatic atmospheric boundary layers: A review and
#' analysis of data from the period 1880–1972.
#' \emph{Atmospheric Environment}, 9(10), 871–905.
#'
#' Schauberger, G. et al. (2002). Calculating direction-dependent separation
#' distance by a dispersion model. \emph{Biosystems Engineering}, 82(1), 25–37.
#' \doi{10.1006/jaer.2001.0943}.
#'
#' Hsu, S.A., Meindl, E.A. and Gilhousen, D.B. (1994). Determining the
#' power-law wind-profile exponent. \emph{J. Appl. Meteorol.}, 33(6), 757–765.
#'
#' Seibert, P. et al. (2000). Review and intercomparison of operational methods
#' for the determination of the mixing height. \emph{Atmos. Environ.},
#' 34(7), 1001–1027.
#'
#' Seinfeld, J.H. and Pandis, S.N. (2016). \emph{Atmospheric Chemistry and
#' Physics}. 3rd edn. Wiley. Ch. 18 (mixing), Ch. 19 (scavenging).
#'
#' Hanna, S.R. (1983). Lateral turbulence intensity and plume meandering
#' during stable conditions. \emph{J. Appl. Meteorol.}, 22(8), 1424–1430.
#'
#' Czepiel, P.M. et al. (2003). The influence of atmospheric pressure on
#' landfill methane emissions. \emph{Waste Management}, 23(7), 593–598.
#'
#' Xu, L. et al. (2014). Surface-atmosphere exchange in a landfill.
#' \emph{Waste Management}, 34(12), 2571–2580.
#'
#' Rees, J.F. (1980). Optimisation of methane production and refuse
#' decomposition in landfills by temperature control.
#' \emph{J. Chem. Technol. Biotechnol.}, 30(1), 458–465.
#'
#' Zou, S.C. et al. (2003). Characterization of ambient volatile organic
#' compounds at a landfill site. \emph{Chemosphere}, 51(9), 1015–1022.
#'
#' Yesiller, N. et al. (2005). Heat generation in municipal solid waste
#' landfills. \emph{Waste Management}, 25(4), 336–348.
#'
#' Hanson, J.L. et al. (2010). Thermal analysis of decomposition in waste.
#' \emph{J. Geotech. Geoenviron. Eng.}, 136(10), 1395–1405.
#'
#' Moldrup, P. et al. (2001). Three-porosity model for predicting the gas
#' diffusion coefficient in undisturbed soil.
#' \emph{Soil Sci. Soc. Am. J.}, 65(3), 613–623.
#'
#' McBain, M.C. et al. (2005). Micrometeorological measurements of N₂O and
#' CH₄ emissions from a municipal solid waste landfill.
#' \emph{Waste Manage. Res.}, 23(5), 409–419.
#'
#' Whiteman, C.D. (2000). \emph{Mountain Meteorology: Fundamentals and
#' Applications}. Oxford University Press.
#'
#' Zardi, D. and Whiteman, C.D. (2013). Diurnal mountain wind systems.
#' In: Chow, F.K. et al. (eds) \emph{Mountain Weather Research and
#' Forecasting}. Springer.
#'
#' Whiteman, C.D. (1982). Breakup of temperature inversions in deep mountain
#' valleys: Part I. \emph{J. Appl. Meteorol.}, 21(3), 270–289.
#'
#' Banta, R.M. (1984). Daytime boundary-layer evolution over mountainous
#' terrain. Part I. \emph{Mon. Weather Rev.}, 112(2), 340–356.
#'
#' @export
generate_odour_risk_index <- function(met_data, receptors, drainage_axes = NULL) {

  # ---- Input validation ---------------------------------------------------- #

  required_met_cols <- c(
    "wind_direction_10m", "wind_speed_10m", "wind_speed_80m",
    "boundary_layer_height", "temperature_2m", "pressure_msl",
    "precipitation", "relative_humidity_2m",
    "cloud_cover", "direct_radiation",
    "soil_moisture_0_to_1cm", "soil_moisture_1_to_3cm"
  )

  checkmate::assert_data_frame(met_data, min.rows = 1)

  missing_cols <- setdiff(required_met_cols, names(met_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "{.arg met_data} is missing required columns: {.val {missing_cols}}.",
      "i" = "See {.code ?generate_odour_risk_index} for the full list of required Open-Meteo column names."
    ))
  }

  for (col in required_met_cols) {
    if (!is.numeric(met_data[[col]])) {
      cli::cli_abort(
        "{.arg met_data} column {.val {col}} must be numeric, not {.cls {class(met_data[[col]])}}."
      )
    }
  }

  checkmate::assert_data_frame(receptors, min.rows = 1)

  if (!all(c("bearing", "distance") %in% names(receptors))) {
    cli::cli_abort(
      "{.arg receptors} must contain columns {.val bearing} and {.val distance}."
    )
  }

  checkmate::assert_numeric(
    receptors$bearing, lower = 0, upper = 360, any.missing = FALSE,
    .var.name = "receptors$bearing"
  )
  checkmate::assert_numeric(
    receptors$distance, lower = .Machine$double.eps, any.missing = FALSE,
    .var.name = "receptors$distance"
  )

  if (!is.null(drainage_axes)) {
    checkmate::assert_data_frame(drainage_axes, min.rows = 1)
    if (!all(c("bearing_from", "weight") %in% names(drainage_axes))) {
      cli::cli_abort(
        "{.arg drainage_axes} must contain columns {.val bearing_from} and {.val weight}."
      )
    }
    checkmate::assert_numeric(
      drainage_axes$bearing_from, lower = 0, upper = 360, any.missing = FALSE,
      .var.name = "drainage_axes$bearing_from"
    )
    checkmate::assert_numeric(
      drainage_axes$weight, lower = .Machine$double.eps, any.missing = FALSE,
      .var.name = "drainage_axes$weight"
    )
  }

  # ---- Extract met vectors ------------------------------------------------- #

  wind_dir <- met_data$wind_direction_10m
  u10      <- met_data$wind_speed_10m
  u80      <- met_data$wind_speed_80m
  pbl_h    <- met_data$boundary_layer_height
  temp     <- met_data$temperature_2m
  pressure <- met_data$pressure_msl
  precip   <- met_data$precipitation
  rh       <- met_data$relative_humidity_2m
  cloud    <- met_data$cloud_cover
  rad      <- met_data$direct_radiation
  sm_0_1   <- met_data$soil_moisture_0_to_1cm
  sm_1_3   <- met_data$soil_moisture_1_to_3cm

  n_t <- nrow(met_data)
  n_r <- nrow(receptors)

  # ---- Constants ----------------------------------------------------------- #

  # Gaussian wind direction profile standard deviation (degrees).
  # Accounts for NWP wind direction uncertainty (±10–15°) and within-hour
  # variability. At 30° off-axis, W_dir ≈ 0.011 (negligible but never
  # exactly zero — use < 1e-10, not == 0, in tests). (Spec §2)
  SIGMA_DIR <- 10.0

  # Lateral dispersion coefficients c_y for PG classes A–F (Briggs 1973 rural).
  # σ_y(x, class) = c_y × x / sqrt(1 + 0.0001 × x)  [metres, all classes]
  # Standard in regulatory dispersion modelling (ISCST3, US EPA AP-42).
  C_Y <- c(0.22, 0.16, 0.11, 0.08, 0.06, 0.04)  # indices 1–6, classes A–F

  # F_ref: normalisation constant for F_disp.
  # Reference condition: class F (s = 5), x = 250 m.
  # σ_y_F(250) = 0.04 × 250 / sqrt(1 + 0.025)  ≈  9.88 m
  # σ_z_F(250) = 0.016 × 250 / (1 + 0.075)     ≈  3.72 m
  # F_ref = 1 / (9.88 × 3.72)                   ≈  0.0272 m⁻²
  # F_disp = 1.0 at this reference by construction.
  sigma_y_F_ref <- 0.04 * 250 / sqrt(1 + 0.0001 * 250)
  sigma_z_F_ref <- 0.016 * 250 / (1 + 0.0003 * 250)
  F_REF <- 1 / (sigma_y_F_ref * sigma_z_F_ref)

  # ---- Calm flag ----------------------------------------------------------- #
  # Wind is treated as calm when u10 < 0.5 m/s or u10 is NA (spec §7, §8).
  # Below 0.5 m/s the anemometer is near its noise floor and direction
  # measurements are unreliable.
  is_calm <- is.na(u10) | u10 < 0.5

  # Safe u10 for arithmetic (replace NA with 0, treated as calm everywhere).
  u10_safe <- ifelse(is.na(u10), 0, u10)

  # Safe precipitation (replace NA with 0 per spec §8).
  precip_safe <- ifelse(is.na(precip), 0, precip)

  # ---- W_spd: Wind speed factor [0.05, 1.0] ------------------------------- #
  # Piecewise function of u10. Peak risk at 0.5–2 m/s (coherent transport,
  # minimal dilution); linear decline as dilution increases above 2 m/s;
  # collapses to 0.05 above 8 m/s. The 0.7 factor for calm reflects plume
  # meander at very low wind speeds (Hanna 1983).
  W_spd <- dplyr::case_when(
    is_calm           ~ 0.7,
    u10_safe <= 2.0   ~ 1.0,
    u10_safe <= 8.0   ~ 1.0 - 0.95 * (u10_safe - 2.0) / 6.0,
    TRUE              ~ 0.05
  )

  # ---- Stability index s [0, 5]: continuous PG class ---------------------- #
  # Maps shear exponent α = ln(u80/u10) / ln(8) to a continuous index s ∈ [0,5]
  # where s = 0 is class A (very unstable) and s = 5 is class F (very stable).
  # Breakpoints from Irwin (1979) and Counihan (1975). Intervals are unequal:
  # the neutral-to-stable transition (D→F) spans a wider α range.
  #
  # α ≤ 0.07              → s = 0            (class A)
  # 0.07 < α ≤ 0.10       → s = (α−0.07)/0.03
  # 0.10 < α ≤ 0.13       → s = 1+(α−0.10)/0.03
  # 0.13 < α ≤ 0.15       → s = 2+(α−0.13)/0.02
  # 0.15 < α ≤ 0.22       → s = 3+(α−0.15)/0.07
  # 0.22 < α ≤ 0.40       → s = 4+(α−0.22)/0.18
  # α > 0.40              → s = 5            (class F)
  #
  # Calm/u80-missing override: s = 4.25 (moderately-to-strongly stable),
  # applied when u10 < 0.5, u80 = NA, or u80 = 0 (shear exponent unreliable).

  use_calm_stab <- is_calm | is.na(u80) | u80 == 0

  u10_for_alpha <- pmax(u10_safe, 0.001)
  u80_for_alpha <- ifelse(is.na(u80) | u80 <= 0, 0.001, u80)
  alpha_raw     <- log(u80_for_alpha / u10_for_alpha) / log(8)

  alpha_to_s <- function(alpha) {
    dplyr::case_when(
      alpha <= 0.07 ~ 0.0,
      alpha <= 0.10 ~ (alpha - 0.07) / 0.03,
      alpha <= 0.13 ~ 1.0 + (alpha - 0.10) / 0.03,
      alpha <= 0.15 ~ 2.0 + (alpha - 0.13) / 0.02,
      alpha <= 0.22 ~ 3.0 + (alpha - 0.15) / 0.07,
      alpha <= 0.40 ~ 4.0 + (alpha - 0.22) / 0.18,
      TRUE          ~ 5.0
    )
  }

  s_t <- ifelse(use_calm_stab, 4.25, alpha_to_s(alpha_raw))

  # ---- h_effective: fallback-applied BL height ----------------------------- #
  # Used for F_disp computation and overnight stable hours accumulation.
  # Fallbacks reflect conservative uncertainty: 200 m (calm, no BL data)
  # assumes a shallow nocturnal BL; 500 m (not calm, no BL data) assumes
  # a moderate daytime BL.
  h_effective <- ifelse(
    is.na(pbl_h),
    ifelse(is_calm, 200, 500),
    pbl_h
  )

  # ---- W_rain: Rainfall scavenging factor [0.05, 1.0] -------------------- #
  # Falling raindrops physically remove gas-phase odorants via below-cloud
  # scavenging (washout). Key landfill odorants (H₂S, methanethiol, DMS,
  # organic acids, amines) are moderately to highly water-soluble — washout
  # coefficients Λ ≈ 10⁻⁴ s⁻¹ for 1–5 mm/h rain imply a plume half-life of
  # ~1–2 hours (Seinfeld & Pandis 2016, Ch. 19). Rain also enhances turbulent
  # mixing and suppresses outdoor exposure (people close windows).
  # The 0.2 mm/h threshold prevents NWP drizzle artefacts from triggering
  # suppression. NA precipitation → precip_safe = 0 → W_rain = 1.0.
  W_rain <- dplyr::case_when(
    precip_safe > 4.0 ~ 0.05,  # heavy rain: near-complete washout + mixing
    precip_safe > 1.0 ~ 0.15,  # moderate rain: substantial scavenging
    precip_safe > 0.2 ~ 0.40,  # light rain/drizzle: partial scavenging
    TRUE              ~ 1.0    # dry: no scavenging
  )

  # ---- G: Generation modifier [~0.80, ~1.80] ------------------------------ #
  # Two physically distinct emission pathways (Spec §3):
  #
  #   G = (1 + dP_mod + R_mod + S_seal + H_mod) + V_mod
  #
  # The bulk gas pathway (1 + dP_mod + R_mod + S_seal + H_mod) represents
  # deep-waste methanogenesis. This is effectively constant on forecast
  # timescales — the diurnal surface temperature signal is undetectable below
  # ~1 m depth (Yesiller et al. 2005; Hanson et al. 2010) — adjusted by
  # pressure tendency, post-rain piston effect, cover sealing, and humidity.
  #
  # V_mod is the surface NMOC volatilisation pathway (Henry's Law,
  # genuinely temperature-dependent on hourly timescales; Zou et al. 2003).
  # It is additive — not multiplicative — because pressure-driven gas migration
  # from 10 m depth is NOT enhanced by surface temperature. The multiplicative
  # coupling in the previous implementation was physically incorrect.

  # ---- dP_mod: Pressure tendency modifier --------------------------------- #
  # Falling atmospheric pressure reduces confining pressure on landfill gas,
  # increasing advective flux through the cover (barometric pumping).
  # Czepiel et al. (2003): r² = 0.95 between pressure and methane emissions.
  # Xu et al. (2014): 35-fold day-to-day variation from this mechanism.
  dP3 <- pressure - dplyr::lag(pressure, 3)

  dP_mod <- dplyr::case_when(
    is.na(dP3) ~ 0.0,   # first 3 rows or NA pressure: no lookback data
    dP3 <= -5  ~ 0.30,  # rapid pressure fall: strong barometric pumping
    dP3 <  0   ~ -0.06 * dP3,  # linear ramp: −0.06 × hPa change
    TRUE       ~ 0.0   # steady or rising: no enhancement
  )

  # ---- R_mod: Post-rain piston effect (source side) ----------------------- #
  # Infiltrating water displaces landfill gas from cover pores (piston effect).
  # McBain et al. (2005): >60% flux increases following rain.
  # GUARD: R_mod = 0 during active rain (piston effect not yet operative;
  # cover sealing during rain is handled separately by S_seal).
  # This branch MUST remain first — P_24 includes currently-falling rain.
  P_24 <- vapply(seq_len(n_t), function(i) {
    if (i == 1L) return(0.0)
    sum(precip_safe[max(1L, i - 24L):(i - 1L)])
  }, numeric(1))

  R_mod <- dplyr::case_when(
    precip_safe > 0.5 ~  0.0,  # active rain: piston effect not yet operative
    P_24 > 15         ~  0.20, # post-heavy rain surge (piston effect)
    P_24 >  5         ~  0.10, # post-moderate rain: mild enhancement
    TRUE              ~  0.0
  )

  # ---- S_seal: Soil moisture sealing (new, source side) ------------------- #
  # Gas diffusivity through soil scales as a power law of air-filled porosity
  # (Moldrup et al. 2001). The wettest layer is the bottleneck for gas-phase
  # diffusion — hence pmax() of both depths. 0.25 m³/m³ onset ≈ field capacity
  # in HTESSEL; 0.40 m³/m³ ≈ near-complete pore filling. −0.20 maximum is
  # conservative: heterogeneous field cover never achieves ideal lab sealing.
  # NA soil moisture → 0 (dry) → S_seal = 0.
  sm_0_1_safe <- ifelse(is.na(sm_0_1), 0, sm_0_1)
  sm_1_3_safe <- ifelse(is.na(sm_1_3), 0, sm_1_3)
  sm_seal <- pmax(sm_0_1_safe, sm_1_3_safe)

  S_seal <- dplyr::case_when(
    sm_seal >= 0.40 ~ -0.20,
    sm_seal >= 0.25 ~ -0.20 * (sm_seal - 0.25) / 0.15,
    TRUE            ~  0.0
  )

  # ---- H_mod: Humidity modifier ------------------------------------------- #
  # High humidity is correlated with increased odour complaints.
  # NA RH → rh_safe = 0 → H_mod = 0.
  rh_safe <- ifelse(is.na(rh), 0, rh)

  H_mod <- dplyr::case_when(
    rh_safe >= 85 ~ 0.15,
    rh_safe >= 60 ~ 0.15 * (rh_safe - 60) / 25,
    TRUE          ~ 0.0
  )

  # ---- V_mod: Surface volatilisation modifier (new) ----------------------- #
  # Henry's Law: NMOC volatilisation (DMS, methanethiol, limonene, acetic acid)
  # depends strongly on surface temperature. Negligible below 10°C (very low
  # Henry's Law constants); saturates near 35°C (boundary-layer turbulence
  # partially offsets further increase); maximum 0.15 (secondary pathway for a
  # well-covered landfill). NA temperature → 0 (no phantom enhancement).
  # (Zou et al. 2003)
  V_mod <- dplyr::case_when(
    is.na(temp) ~ 0.0,
    temp <= 10  ~ 0.0,
    temp >= 35  ~ 0.15,
    TRUE        ~ 0.15 * (temp - 10) / 25
  )

  G <- (1.0 + dP_mod + R_mod + S_seal + H_mod) + V_mod

  # ---- Drainage / fumigation state ---------------------------------------- #
  # Initialise as all-inactive (used when drainage_axes = NULL and inside loop).
  is_drainage          <- rep(FALSE, n_t)
  is_fumigation        <- rep(FALSE, n_t)
  effective_severity_v <- rep(0.0, n_t)
  max_alignment        <- rep(0.0, n_r)  # per receptor, precomputed

  if (!is.null(drainage_axes)) {

    # NA fallbacks for cloud and radiation (spec §6).
    cloud_safe <- ifelse(is.na(cloud), 50, cloud)
    rad_safe   <- ifelse(is.na(rad),   0,  rad)

    # ---- Precompute drainage alignment per receptor ------------------------ #
    # For receptor j and drainage axis k:
    #   delta_k = angular difference between receptor bearing and drainage axis.
    #   If delta_k ≤ 30°: alignment_k = cos(π × delta_k / 60)² × weight_k
    #   max_alignment_j = min(1, max(alignment_k across all k))
    for (j in seq_len(n_r)) {
      alignments <- vapply(seq_len(nrow(drainage_axes)), function(k) {
        delta_k <- abs(receptors$bearing[j] - drainage_axes$bearing_from[k])
        delta_k <- min(delta_k, 360 - delta_k)
        if (delta_k <= 30) {
          cos(pi * delta_k / 60)^2 * drainage_axes$weight[k]
        } else {
          0.0
        }
      }, numeric(1))
      max_alignment[j] <- min(1.0, max(alignments))
    }

    # ---- Sequential state computation -------------------------------------- #
    # Drainage, accumulation, and fumigation all depend on previous timesteps
    # and must be computed in a sequential loop.

    drainage_hours         <- integer(n_t)
    latched_drainage_dur   <- integer(n_t)
    overnight_stable_hours <- integer(n_t)
    latched_stable_dur     <- integer(n_t)
    fumigation_pool_consumed <- logical(n_t)

    fumigation_event_start <- NA_integer_
    fum_severity_at_start  <- 0.0  # stored at event start, used for all 3 hours
    latch_consumed_flag    <- FALSE

    for (t in seq_len(n_t)) {

      is_dark_t      <- rad_safe[t] < 10
      # Sunset reset trigger: first timestep of a new night.
      is_new_night_t <- is_dark_t && (t == 1L || rad_safe[t - 1L] >= 10)

      # ---- Drainage activation ------------------------------------------- #
      # Katabatic (drainage) flow toward landfill: calm + dark + clear sky.
      is_drainage[t] <- is_calm[t] && is_dark_t && cloud_safe[t] < 70

      # ---- Drainage hours counter ----------------------------------------- #
      if (t == 1L) {
        drainage_hours[t] <- if (is_drainage[t]) 1L else 0L
      } else {
        drainage_hours[t] <- if (is_drainage[t]) drainage_hours[t - 1L] + 1L else 0L
      }

      # ---- Overnight stable hours counter ---------------------------------- #
      # Broader accumulation metric: consecutive hours with s > 4.0 (≈ PG E or
      # more stable), h_effective < 300 m, and no solar radiation. Uses
      # fallback-applied h_effective so NA boundary layer height does not
      # silently prevent accumulation on stable nights. (Spec §6)
      is_stable_t  <- s_t[t] > 4.0
      is_shallow_t <- h_effective[t] < 300

      if (t == 1L) {
        overnight_stable_hours[t] <- if (is_stable_t && is_shallow_t && is_dark_t) 1L else 0L
      } else {
        overnight_stable_hours[t] <- if (is_stable_t && is_shallow_t && is_dark_t) {
          overnight_stable_hours[t - 1L] + 1L
        } else {
          0L
        }
      }

      # ---- Latched durations ---------------------------------------------- #
      # Preserves peak counter across the gap between drainage ending and
      # fumigation starting. Sunset reset clears stale latches so Night 1
      # drainage does not contaminate Day 2 fumigation severity if Day 1
      # had no fumigation. (Spec §6)
      if (is_new_night_t || latch_consumed_flag) {
        # Sunset reset or latch consumed by fumigation: clear to 0.
        latched_drainage_dur[t] <- 0L
        latched_stable_dur[t]   <- 0L
        latch_consumed_flag      <- FALSE

      } else if (t == 1L) {
        latched_drainage_dur[t] <- 0L
        latched_stable_dur[t]   <- 0L

      } else {
        # Drainage latch: capture peak when drainage just ended.
        if (drainage_hours[t - 1L] > 0L && drainage_hours[t] == 0L) {
          latched_drainage_dur[t] <- drainage_hours[t - 1L]
        } else if (drainage_hours[t] > 0L) {
          latched_drainage_dur[t] <- 0L  # still active: not yet latched
        } else {
          latched_drainage_dur[t] <- latched_drainage_dur[t - 1L]
        }

        # Stable hours latch: same logic.
        if (overnight_stable_hours[t - 1L] > 0L && overnight_stable_hours[t] == 0L) {
          latched_stable_dur[t] <- overnight_stable_hours[t - 1L]
        } else if (overnight_stable_hours[t] > 0L) {
          latched_stable_dur[t] <- 0L
        } else {
          latched_stable_dur[t] <- latched_stable_dur[t - 1L]
        }
      }

      # ---- Fumigation pool state ------------------------------------------ #
      # Pool resets each night; carried forward during daytime to prevent
      # a second event on the same morning. (Spec §6)
      if (is_dark_t) {
        fumigation_pool_consumed[t] <- FALSE
      } else if (t == 1L) {
        fumigation_pool_consumed[t] <- FALSE
      } else {
        fumigation_pool_consumed[t] <- fumigation_pool_consumed[t - 1L]
      }

      # ---- Fumigation event management ------------------------------------ #
      # Activation: radiation > 50 W/m² AND at least one of the preceding
      # 6 timesteps was dark. IMPORTANT: fumigation operates independently of
      # is_calm — the morning transition is precisely when synoptic wind starts
      # registering. (Spec §6)
      preceding_idx   <- if (t > 1L) max(1L, t - 6L):(t - 1L) else integer(0)
      has_recent_dark <- length(preceding_idx) > 0 &&
        any(rad_safe[preceding_idx] < 10)

      activation_criteria <- rad_safe[t] > 50 && has_recent_dark

      # End active event if 3-hour window elapsed or it gets dark.
      if (!is.na(fumigation_event_start)) {
        if ((t - fumigation_event_start) >= 3L || is_dark_t) {
          fumigation_event_start <- NA_integer_
          fum_severity_at_start  <- 0.0
        }
      }

      # Start new event if criteria met and pool not consumed.
      if (activation_criteria && !fumigation_pool_consumed[t] &&
          is.na(fumigation_event_start)) {
        # Compute severity from latch at the moment fumigation fires.
        ldd <- latched_drainage_dur[t]
        lsd <- latched_stable_dur[t]
        accumulation_hours <- if (ldd > 0L) {
          as.numeric(ldd)
        } else {
          as.numeric(lsd) * 0.5  # non-drainage: less effective trapping
        }
        fum_severity_at_start  <- min(1.0, accumulation_hours / 8.0)
        fumigation_event_start <- t
        fumigation_pool_consumed[t] <- TRUE
        latch_consumed_flag <- TRUE  # latch cleared on next timestep
      }

      is_fumigation[t] <- !is.na(fumigation_event_start)

      # ---- Effective severity: decay over 3-hour window ------------------- #
      if (is_fumigation[t]) {
        hours_since_start    <- t - fumigation_event_start
        decay                <- max(0.0, 1.0 - hours_since_start / 3.0)
        effective_severity_v[t] <- fum_severity_at_start * decay
      }

    }  # end sequential loop
  }  # end if (!is.null(drainage_axes))

  # ---- W_spd_effective: drainage override ---------------------------------- #
  # During drainage, W_spd = 0.3 (katabatic drainage confines emissions within
  # the hollow, reducing effective transport toward receptors; lower than
  # generic calm's 0.7). (Spec §6 precedence hierarchy)
  W_spd_effective <- W_spd
  if (!is.null(drainage_axes) && any(is_drainage)) {
    W_spd_effective[is_drainage] <- 0.3
  }

  # ---- G_effective: fumigation G boost ------------------------------------ #
  # During fumigation, the convective BL entrains accumulated overnight
  # emissions, boosting the effective source term.
  # G_eff = G × (1 + 0.5 × effective_severity). (Spec §6)
  G_effective <- G
  if (!is.null(drainage_axes) && any(is_fumigation)) {
    G_effective[is_fumigation] <- G[is_fumigation] *
      (1.0 + 0.5 * effective_severity_v[is_fumigation])
  }

  # ---- Downwind direction -------------------------------------------------- #
  # Open-Meteo wind_direction_10m is the direction wind blows FROM. The
  # downwind (plume travel) direction is the reciprocal.
  theta_down <- (ifelse(is.na(wind_dir), NA_real_, wind_dir) + 180) %% 360

  # ---- Compute risk matrix [n_t × n_r] ------------------------------------ #
  risk_matrix <- matrix(0.0, nrow = n_t, ncol = n_r)

  for (j in seq_len(n_r)) {

    x_j <- receptors$distance[j]

    # ---- Precompute σ_y and σ_z at x_j for all 6 PG classes -------------- #
    # σ_y(x, class) = c_y × x / sqrt(1 + 0.0001 × x)  [same form all classes]
    # σ_z formulas differ by class (suppressed growth under stable conditions):
    #   A: σ_z = 0.20 × x
    #   B: σ_z = 0.12 × x
    #   C: σ_z = 0.08 × x / sqrt(1 + 0.0002 × x)
    #   D: σ_z = 0.06 × x / sqrt(1 + 0.0015 × x)
    #   E: σ_z = 0.03 × x / (1 + 0.0003 × x)
    #   F: σ_z = 0.016 × x / (1 + 0.0003 × x)
    # For unstable classes (A, B), σ_z grows linearly (rapid mixing).
    # For stable classes (E, F), the denominator suppresses growth at large
    # distances, reflecting stable stratification.
    #
    # F_disp validation table (F_disp = F_raw / F_ref, clamped to [0.003, 1.0]):
    #   Stability  Dist   BL     σ_y     σ_z    h_eff  F_raw   F_disp
    #   F (s=5)   250 m  50 m   9.88   3.72    3.72  0.0272   1.00
    #   F (s=5)   500 m  100m  19.51   6.96    6.96  0.00736  0.271
    #   F (s=5)  1000 m  100m  38.14  12.31   12.31  0.00213  0.078
    #   D (s=3)   500 m  800m  38.83  22.68   22.68  0.00114  0.042
    #   D (s=3)   250 m  800m  19.70  13.04   13.04  0.00389  0.143
    #   A (s=0)   500 m 2000m 107.76 100.00  100.00  0.0000928 0.003
    #   A (s=0)   250 m 2000m  54.39  50.00   50.00  0.000368  0.014
    #
    # (Values are approximate; implementation uses the exact Briggs formulas.)
    sigma_y_classes <- C_Y * x_j / sqrt(1 + 0.0001 * x_j)
    sigma_z_classes <- c(
      0.20  * x_j,
      0.12  * x_j,
      0.08  * x_j / sqrt(1 + 0.0002 * x_j),
      0.06  * x_j / sqrt(1 + 0.0015 * x_j),
      0.03  * x_j / (1 + 0.0003 * x_j),
      0.016 * x_j / (1 + 0.0003 * x_j)
    )

    # ---- Vectorised F_disp over timesteps --------------------------------- #
    # Linear interpolation between adjacent PG class outputs using continuous
    # stability index s_t. This is smooth, monotonic, and preserves the
    # physical behaviour of each bounding class.
    #
    # O(n_r) precomputation (sigma_y_classes, sigma_z_classes above) plus
    # O(n_t) vectorised arithmetic — avoids an O(n_t × n_r) scalar loop.
    s_low  <- floor(s_t)           # integer 0–5
    s_high <- pmin(s_low + 1, 5)  # integer 1–5 (capped at class F)
    frac   <- s_t - s_low         # fractional part [0, 1)

    # Vectorised indexing: R is 1-based, so add 1 to s_low/s_high.
    sigma_y_t <- sigma_y_classes[s_low + 1] +
      frac * (sigma_y_classes[s_high + 1] - sigma_y_classes[s_low + 1])
    sigma_z_t <- sigma_z_classes[s_low + 1] +
      frac * (sigma_z_classes[s_high + 1] - sigma_z_classes[s_low + 1])

    # h_eff: effective vertical mixing depth.
    # When σ_z < h (near-field): stability controls vertical dispersion.
    # When σ_z ≥ h (well-mixed): concentration scales with h, not σ_z.
    h_eff_j <- pmin(sigma_z_t, h_effective)

    F_raw_j  <- 1 / (sigma_y_t * h_eff_j)
    F_disp_j <- pmax(0.003, pmin(1.0, F_raw_j / F_REF))

    # ---- W_dir_j: Gaussian wind direction profile ------------------------- #
    # exp(-0.5 × (Δθ / σ_dir)²) with σ_dir = 10°. More physically accurate
    # than a hard-cutoff squared-cosine: reflects the lateral Gaussian
    # concentration profile of a Gaussian plume. At 90° off-axis,
    # W_dir ≈ 10⁻¹⁸ (use < 1e-10 in tests, not == 0).
    # Calm/NA direction: W_dir = 0.5 (same as previous implementation).
    theta_j     <- receptors$bearing[j]
    diff_raw    <- abs(theta_down - theta_j)
    delta_theta <- pmin(diff_raw, 360 - diff_raw)

    W_dir_j <- dplyr::case_when(
      is_calm | is.na(wind_dir) ~ 0.5,
      TRUE                      ~ exp(-0.5 * (delta_theta / SIGMA_DIR)^2)
    )

    # ---- Drainage override ------------------------------------------------- #
    # During katabatic drainage, odour disperses within the topographic hollow.
    # W_dir is much lower than generic calm's 0.5: 0.05 on the drainage axis
    # (flow directed toward landfill), 0.10 off-axis. (Spec §6)
    if (!is.null(drainage_axes) && any(is_drainage)) {
      drain_W_dir_j <- 0.05 + 0.05 * (1 - max_alignment[j])
      W_dir_j[is_drainage] <- drain_W_dir_j
    }

    # ---- Fumigation W_dir floor ------------------------------------------- #
    # During fumigation, anabatic (upslope) flow transports the accumulated
    # overnight pool toward aligned receptors. 0.6 ceiling: anabatic flow is
    # weaker and less coherent than synoptically-driven transport. (Spec §6)
    if (!is.null(drainage_axes) && any(is_fumigation)) {
      W_dir_fum_floor <- 0.6 * effective_severity_v * max_alignment[j]
      W_dir_j[is_fumigation] <- pmax(
        W_dir_j[is_fumigation],
        W_dir_fum_floor[is_fumigation]
      )
    }

    # ---- Per-receptor risk: product of all factors (no clamping) ---------- #
    # R_{t,j} = W_dir_j × W_spd × F_disp_j × W_rain × G
    # Raw scores > 1.0 are possible under extreme conditions (~1.80 maximum).
    risk_matrix[, j] <- W_dir_j * W_spd_effective * F_disp_j * W_rain * G_effective
  }

  # ---- Return maximum risk across receptors for each timestep -------------- #
  # The worst-affected receptor determines operational risk at each hour.
  apply(risk_matrix, 1, max)
}
