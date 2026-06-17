#' Litter Risk Index (LRI) for windblown litter dispersal at landfill sites
#'
#' Computes an hourly Litter Risk Index (LRI) in the range \[0, 100\] for
#' forecasting windblown litter dispersal conditions at a landfill site. The
#' LRI is computed from pre-fetched Open-Meteo forecast data and a set of
#' site-specific configuration parameters.
#'
#' The function does not make API calls. The calling code is responsible for
#' fetching the required hourly variables from the Open-Meteo \code{/v1/forecast}
#' endpoint (with \code{forecast_days = 7} and \code{timezone = "Australia/Sydney"})
#' and passing them as arguments. All meteorological inputs accept either a
#' scalar (single forecast hour) or a numeric vector (multiple hours), enabling
#' direct use inside \code{\link[dplyr]{mutate}}.
#'
#' @section Model structure:
#'
#' The LRI uses a \strong{multiplicative model with two wetness gates} — one
#' binary (active rainfall) and one continuous (soil moisture):
#'
#' \deqn{LRI_{raw} = f_G(G) \times S_W(W) \times R(P_0) \times V(SM) \times
#'   (1 + B(C,\,d)) \times M(D)}
#' \deqn{LRI = \min(LRI_{raw},\; 100)}
#'
#' Components:
#' \describe{
#'   \item{\eqn{f_G(G)}}{Gust mobilisation function (primary driver). Power-law
#'     ramp above a minimum threshold that saturates at a reference gust speed.
#'     Exponent \eqn{n = 2} follows from the aerodynamic drag equation
#'     \eqn{F_D = \frac{1}{2}\rho C_D A v^2}; the fraction of waste items
#'     mobilised scales with the excess kinetic energy above the detachment
#'     threshold. Range: \[0, \code{gust_max_score}\].}
#'   \item{\eqn{S_W(W)}}{Sustained wind transport modifier. Multiplier on the
#'     gust score reflecting how far mobilised items travel. Structured as a
#'     multiplier (not an additive term) because sustained wind without
#'     sufficient gust energy cannot mobilise items in the first place. Range:
#'     \[1.0, \code{wind_mult_max}\].}
#'   \item{\eqn{R(P_0)}}{Rainfall hard gate. Binary suppressor: when measurable
#'     rain is falling, all litter risk is instantly set to zero. This provides
#'     immediate response that the soil moisture model (which may lag by one or
#'     more time steps) cannot. Range: \{0, 1\}.}
#'   \item{\eqn{V(SM)}}{Soil moisture wetness gate. Linear ramp between fully
#'     dry and fully wet thresholds, using Open-Meteo's land-surface model
#'     output which already integrates antecedent rainfall, evapotranspiration,
#'     temperature, and drainage history. This offloads drying-rate calculation
#'     to the NWP model. Range: \[0.0, 1.0\].}
#'   \item{\eqn{B(C,\,d)}}{Atmospheric instability bonus. On hot days with high
#'     CAPE, thermal updrafts can loft lightweight items (plastic bags, film) to
#'     heights sufficient to clear perimeter fencing. Zero at night because
#'     surface-driven convection ceases. Range: \[0.0, \code{instability_max}\].}
#'   \item{\eqn{M(D)}}{Directional attenuation factor. Reflects the site
#'     geometry constraint that litter can only escape off-site in specific wind
#'     directions. When wind blows toward a sensitive receptor the factor is 1.0;
#'     when blowing toward barriers or closed aspects it is reduced to
#'     \code{offwind_factor}. Range: \[\code{offwind_factor}, 1.0\].}
#' }
#'
#' The multiplicative structure guarantees that any single suppressing condition
#' drives the entire score to zero — a property that additive models cannot
#' achieve.
#'
#' @section Theoretical range analysis:
#'
#' Maximum possible \eqn{LRI_{raw}} (all components at maxima, wind toward receptor):
#' \deqn{50 \times 2.0 \times 1.0 \times 1.0 \times 1.2 \times 1.0 = 120}
#' The \eqn{\min(\cdot, 100)} cap therefore activates under simultaneous extreme
#' gusts, high sustained wind, full dryness, no rain, high instability, and wind
#' toward a receptor.
#'
#' Maximum \eqn{LRI_{raw}} when wind is NOT toward a receptor (using defaults):
#' \deqn{50 \times 2.0 \times 1.0 \times 1.0 \times 1.2 \times 0.5 = 60}
#' This means EXTREME tier (> 70) is only reachable when wind blows toward a
#' sensitive receptor — a deliberate design feature that concentrates the
#' highest-severity warnings on the directions where off-site escape is
#' physically possible.
#'
#' @section Important caveat:
#'
#' No published composite litter risk index exists in the academic or regulatory
#' literature. This formulation is a novel construction synthesised from the best
#' available empirical data on macroplastic mobilisation thresholds (Mellink
#' et al. 2024), general aerodynamic principles, and standard landfill
#' operational guidance. It should be treated as a first-approximation model
#' requiring empirical calibration against on-site litter escape observations
#' over at least 12 months before tier boundaries are considered reliable.
#' The calling code is responsible for mapping the returned numeric value to
#' operational tiers (see the Tier mapping section below for suggested
#' boundaries).
#'
#' @param wind_gusts_10m Numeric vector. Peak wind gust at 10 m above ground
#'   level (km/h), Open-Meteo variable \code{wind_gusts_10m}. The primary
#'   mobilisation driver: gust energy (not sustained wind) determines whether
#'   discrete waste items are dislodged, consistent with the aerodynamic drag
#'   equation. Mellink et al. (2024) found that 70% of macroplastic item types
#'   mobilised from compacted surfaces below 31.7 km/h.
#'
#' @param wind_speed_10m Numeric vector. Hourly mean wind speed at 10 m (km/h),
#'   Open-Meteo variable \code{wind_speed_10m}. Used as the transport distance
#'   modifier: once items are dislodged by gusts, sustained wind governs how
#'   far they travel and whether partially dislodged items continue to move.
#'
#' @param wind_direction_10m Numeric vector. Wind direction at 10 m (degrees
#'   \[0, 360\], meteorological convention: the direction wind blows \emph{from}),
#'   Open-Meteo variable \code{wind_direction_10m}. Used to determine whether
#'   wind is blowing toward a sensitive receptor, modulating the directional
#'   attenuation factor \eqn{M}.
#'
#' @param precipitation Numeric vector. Hourly precipitation sum (mm), Open-Meteo
#'   variable \code{precipitation}. The rainfall hard gate \eqn{R} suppresses
#'   all litter risk instantly when this exceeds \code{rain_threshold}. Even
#'   light rain (0.5 mm/h) wets exposed waste surfaces and increases item mass
#'   through water absorption; Mellink et al. (2024) found rain significantly
#'   reduced macroplastic mobilisation.
#'
#' @param soil_moisture_0_to_1cm Numeric vector. Volumetric water content of
#'   the top 0--1 cm soil layer (m\eqn{^3}/m\eqn{^3}), Open-Meteo variable
#'   \code{soil_moisture_0_to_1cm}. This is computed by Open-Meteo's ERA5-based
#'   land-surface model and already integrates antecedent rainfall,
#'   evapotranspiration, temperature, wind speed, solar radiation, and drainage.
#'   Using it directly offloads the drying-rate calculation to the NWP model,
#'   eliminating the need for client-side antecedent dry period tracking.
#'   \strong{Important:} absolute values represent natural soil at the grid-cell
#'   location, not a landfill surface. The relative wetting/drying signal (rising
#'   after rain, falling during dry periods) should track correctly, but absolute
#'   thresholds (\code{soil_wet} and \code{soil_dry}) must be calibrated against
#'   observed site conditions.
#'
#' @param cape Numeric vector. Convective Available Potential Energy (J/kg),
#'   Open-Meteo variable \code{cape}. On hot days, convective thermals produce
#'   vertical velocities of 1--5 m/s, easily exceeding the terminal velocity of
#'   plastic bags (~1--2 m/s), lofting them to heights that can clear perimeter
#'   fencing. This is modelled as a secondary bonus — it worsens risk but cannot
#'   create it independently. There is no published quantitative model linking
#'   CAPE to litter lofting probability at landfills.
#'
#' @param is_day Integer or numeric vector (0 or 1). Daytime flag; 1 during
#'   daylight hours, 0 at night. Open-Meteo variable \code{is_day}. At night
#'   the atmospheric instability bonus is set to zero because surface-driven
#'   convection ceases; CAPE-driven lofting is a daytime phenomenon.
#'
#' @param receptor_arcs A list of \code{character(2)} vectors specifying the
#'   compass arcs toward sensitive receptors, in clockwise order. Each element
#'   gives the \emph{start} and \emph{end} compass label of the arc; wind
#'   blowing \emph{from} directions within the arc is considered to be blowing
#'   \emph{toward} that receptor.
#'
#'   Valid compass labels and their degree equivalents:
#'   \tabular{ll}{
#'     \code{"N"}  \tab 0°  \cr
#'     \code{"NE"} \tab 45° \cr
#'     \code{"E"}  \tab 90° \cr
#'     \code{"SE"} \tab 135° \cr
#'     \code{"S"}  \tab 180° \cr
#'     \code{"SW"} \tab 225° \cr
#'     \code{"W"}  \tab 270° \cr
#'     \code{"NW"} \tab 315°
#'   }
#'
#'   Example: \code{list(c("W", "N"), c("SE", "S"))} specifies two receptors.
#'   Receptor 1: wind from 270° clockwise to 0° (i.e. westerly through
#'   northerly). Receptor 2: wind from 135° to 180° (southeasterly to southerly).
#'   An arc that crosses north (e.g. \code{c("W", "N")}: start 270° > end 0°)
#'   is handled correctly via modular arithmetic.
#'
#' @param gust_threshold Numeric scalar. Minimum gust speed for litter
#'   mobilisation, \eqn{G_0} (km/h). Default: 15. Below this gust speed the
#'   LRI is zero regardless of all other conditions. The threshold of 15 km/h
#'   (4.2 m/s) is slightly above the lightest-item critical mobilisation
#'   velocity measured by Mellink et al. (2024) for plastic bags on compacted
#'   surfaces (~8.3 km/h), elevated to account for the partial anchoring effect
#'   of overlying waste and compaction on a landfill working face. Moderate
#'   calibration uncertainty; can be refined by logging gust speed at first
#'   observed litter movement.
#'
#' @param gust_ref Numeric scalar. Gust speed at which the gust score reaches
#'   its maximum and is capped, \eqn{G_{ref}} (km/h). Default: 50. Above this
#'   speed, nearly all loose waste is already mobilised and additional gust
#'   energy does not increase the mobilised fraction. Must be strictly greater
#'   than \code{gust_threshold}. Moderate calibration uncertainty; can be
#'   refined by observing whether mobilisation plateaus above ~50 km/h gusts.
#'
#' @param gust_max_score Numeric scalar. Maximum score contribution from the
#'   gust component, \eqn{a}. Default: 50. This occupies the lower half of the
#'   \[0, 100\] range; reaching the upper half additionally requires elevated
#'   sustained wind, dryness, instability, and directional alignment.
#'
#' @param gust_exponent Numeric scalar. Power-law exponent for the gust
#'   mobilisation function, \eqn{n}. Default: 2. This value is a theoretical
#'   prior from the drag equation (\eqn{F_D \propto v^2}) rather than an
#'   empirically fitted value. It is a strong candidate for calibration via
#'   regression of observed litter frequency against gust speed.
#'
#' @param wind_threshold Numeric scalar. Sustained wind speed below which the
#'   transport modifier is 1.0 (no enhancement), \eqn{W_0} (km/h). Default: 20.
#'   Approximately Beaufort Force 4 ("moderate breeze"), above which lightweight
#'   items already in motion begin to travel significant distances.
#'
#' @param wind_ref Numeric scalar. Sustained wind speed at which the transport
#'   modifier reaches \code{wind_mult_max}, \eqn{W_{ref}} (km/h). Default: 55.
#'   Corresponds to near-gale conditions where all mobilised items will travel
#'   beyond the site boundary regardless of transport distance. Must be strictly
#'   greater than \code{wind_threshold}.
#'
#' @param wind_mult_max Numeric scalar. Maximum sustained wind multiplier,
#'   \eqn{S_{max}}. Default: 2.0. At the maximum, sustained high wind can at
#'   most double the gust-driven score. The upper bound of 2.0 is conservative
#'   by design: sustained wind cannot mobilise items without sufficient gust
#'   energy first, so its influence is bounded as a multiplier. There is no
#'   published quantitative relationship between sustained wind speed and litter
#'   transport distance at landfills; the range \[1.0, 2.0\] reflects engineering
#'   judgement.
#'
#' @param rain_threshold Numeric scalar. Hourly precipitation at or above which
#'   the rainfall hard gate fully suppresses litter risk, \eqn{P_{rain}} (mm).
#'   Default: 0.5. This corresponds to light but measurable rain, sufficient to
#'   wet surfaces within minutes. Drizzle below 0.5 mm/h may not reliably
#'   suppress mobilisation of waterproof plastic film. The rainfall gate
#'   provides immediate suppression at rain onset, compensating for the lag that
#'   may occur in the soil moisture model.
#'
#' @param soil_wet Numeric scalar. Soil moisture at or above which the surface
#'   is considered too wet for litter mobilisation, \eqn{SM_{wet}}
#'   (m\eqn{^3}/m\eqn{^3}). Default: 0.20. \strong{High calibration
#'   uncertainty:} this threshold is an initial estimate. To calibrate, record
#'   Open-Meteo \code{soil_moisture_0_to_1cm} values at times when waste
#'   surfaces are visibly wet after rain. The broad default is deliberately set
#'   to avoid premature false negatives during the calibration period.
#'
#' @param soil_dry Numeric scalar. Soil moisture at or below which the surface
#'   is considered fully dry (gate = 1.0), \eqn{SM_{dry}}
#'   (m\eqn{^3}/m\eqn{^3}). Default: 0.05. Must be strictly less than
#'   \code{soil_wet}. \strong{High calibration uncertainty:} record values when
#'   waste is fully dry and litter is observed to mobilise.
#'
#' @param cape_ref Numeric scalar. CAPE value at which the instability bonus
#'   reaches its maximum, \eqn{C_{ref}} (J/kg). Default: 1000. Represents
#'   moderate convective instability in the Australian context. Bureau of
#'   Meteorology severe thunderstorm warnings typically begin at approximately
#'   1500--2000 J/kg.
#'
#' @param instability_max Numeric scalar. Maximum instability bonus,
#'   \eqn{B_{max}}. Default: 0.2 (20% increase). Reflects the judgement that
#'   convective thermals are a secondary risk factor: they worsen an already-bad
#'   situation but cannot generate risk independently. The bonus is bounded at
#'   20% to maintain the gust mobilisation function as the primary driver.
#'
#' @param offwind_factor Numeric scalar. Attenuation factor applied when wind
#'   does not blow toward any sensitive receptor, \eqn{M_{offwind}}. Default:
#'   0.5. Even when wind blows toward barriers or closed aspects, some litter
#'   can recirculate via turbulent eddies, be lofted over obstacles by thermals,
#'   or escape via vehicle movements — so the factor is not zero. Calibration
#'   guidance (NSW EPA, 2016): use 0.3 for solid engineered walls or dense
#'   tree belts; use 0.7 for sparse vegetation or low earth bunds. High
#'   calibration uncertainty.
#'
#' @param direction_tol Numeric scalar. Wind direction tolerance added to each
#'   edge of a receptor arc, \eqn{\Delta\theta} (degrees). Default: 15.
#'   Accounts for forecast wind direction error (typically ±10--20° for NWP
#'   models) and for the lateral spread of litter dispersal plumes. The arc
#'   containment test handles the 360°/0° wraparound correctly after expansion.
#'
#' @return Numeric vector of length equal to the input vectors, giving the
#'   Litter Risk Index in the range \[0, 100\] for each forecast hour. Returns
#'   exactly 0 when gusts are sub-threshold, rainfall is active, or the surface
#'   is saturated.
#'
#' @section Tier mapping (suggested, requires calibration):
#' The raw numeric LRI can be mapped to four operational tiers by the calling
#' code. Suggested boundaries (high uncertainty — see Calibration priorities):
#' \tabular{rllll}{
#'   \strong{LRI} \tab \strong{Tier} \tab \strong{Label} \tab \strong{Colour}
#'     \tab \strong{Meaning} \cr
#'   0--19  \tab 1 \tab LOW      \tab green  \tab Dispersal unlikely \cr
#'   20--44 \tab 2 \tab MODERATE \tab yellow \tab Enhanced controls warranted \cr
#'   45--69 \tab 3 \tab HIGH     \tab orange \tab Dispersal likely without intervention \cr
#'   70--100 \tab 4 \tab EXTREME \tab red    \tab Maximum controls or cessation required
#' }
#' Example tier assignment: \code{cut(lri, breaks = c(-Inf, 20, 45, 70, Inf),}
#' \code{labels = c("LOW", "MODERATE", "HIGH", "EXTREME"), right = FALSE)}.
#' The recommended calibration approach is ROC (receiver operating
#' characteristic) analysis against a 12-month log of observed litter escape
#' events.
#'
#' @section Calibration priorities:
#' Parameters most in need of site-specific calibration, in approximate priority
#' order:
#' \enumerate{
#'   \item \code{soil_wet} / \code{soil_dry} — the NWP soil moisture is fitted
#'     to natural soil, not a landfill surface.
#'   \item Tier boundaries — initial estimates only; ROC analysis required
#'     after 12 months of litter escape event logging (see Tier mapping section).
#'   \item \code{offwind_factor} — depends strongly on actual barrier quality
#'     and maintenance.
#'   \item \code{gust_threshold} / \code{gust_ref} — can be refined by
#'     logging gust speed at times of observed litter movement.
#'   \item \code{gust_exponent} — regression of observed litter frequency
#'     against gust speed once sufficient event data are available.
#' }
#'
#' @section Known limitations:
#' \enumerate{
#'   \item \strong{Soil moisture model mismatch.} The \code{soil_moisture_0_to_1cm}
#'     field is modelled for natural soil at the NWP grid-cell location. A
#'     landfill surface has different infiltration, retention, and evaporation
#'     characteristics. The relative signal (wetting and drying trends) should
#'     track correctly, but absolute thresholds require calibration.
#'   \item \strong{Topographic wind effects.} Open-Meteo's forecast models use
#'     approximately 9--15 km grid resolution. In topographically complex terrain
#'     (such as the Blue Mountains, NSW), local wind acceleration over ridges,
#'     valley channelling, and lee-side turbulence may cause actual gust speeds
#'     to differ from forecast values by 20--50%.
#'   \item \strong{Waste composition variability.} The model treats litter as a
#'     homogeneous category. A working face receiving heavy construction and
#'     demolition waste has far lower litter risk than one receiving municipal
#'     kerbside waste. A waste-type modifier based on scheduled tipping
#'     operations is a candidate for future refinement.
#'   \item \strong{Litter fence state.} The LRI represents the meteorological
#'     hazard, not residual risk after physical controls. Degraded or poorly
#'     positioned fencing will allow litter escape at lower LRI values than the
#'     tier boundaries suggest.
#'   \item \strong{Soil moisture response lag.} At the very onset of light rain,
#'     the soil moisture field may not yet have risen to \code{soil_wet}. The
#'     separate rainfall gate compensates by providing immediate suppression
#'     when measurable precipitation is detected. In practice, the lag is
#'     sub-hourly and rain onset itself physically suppresses mobilisation before
#'     the surface is fully wetted.
#'   \item \strong{Novel composite index.} No peer-reviewed formulation of this
#'     specific index exists. All inter-component weights and tier boundaries
#'     are based on engineering judgement and must be empirically validated.
#' }
#'
#' @references
#' Mellink, Y., Roebroek, C.T.J., van Emmerik, T.H.M. et al. (2024).
#' Wind- and rain-driven macroplastic mobilization and transport on land.
#' \emph{Scientific Reports} 14, 5006. \doi{10.1038/s41598-024-53971-8}.
#' Empirical basis for \eqn{G_0}: on compacted surfaces, 70% of macroplastic
#' item types mobilised below 31.7 km/h; plastic bags mobilised at ~8.3 km/h.
#' Also the basis for the rainfall suppression gate.
#'
#' NSW EPA (2016). \emph{Environmental Guidelines: Solid Waste Landfills}
#' (2nd Edition).
#' \url{https://www.epa.nsw.gov.au/sites/default/files/solid-waste-landfill-guidelines-160259.pdf}
#' Regulatory framework requiring operators to prevent off-site amenity impacts;
#' litter is a specifically identified regulated nuisance. Source of the
#' site-geometry rationale for the directional attenuation factor \eqn{M}.
#'
#' Open-Meteo. Weather Forecast API documentation.
#' \url{https://open-meteo.com/en/docs}
#' Source of all meteorological input variables.
#'
#' The \eqn{v^2} dependence in \eqn{f_G} follows from the standard
#' aerodynamic drag equation \eqn{F_D = \frac{1}{2}\rho C_D A v^2}; see any
#' introductory aerodynamics or boundary-layer meteorology text.
#'
#' @export
litter_risk_index <- function(
  wind_gusts_10m,
  wind_speed_10m,
  wind_direction_10m,
  precipitation,
  soil_moisture_0_to_1cm,
  cape,
  is_day,
  receptor_arcs,
  # --- Gust mobilisation function (f_G) parameters ---
  gust_threshold  = 15,    # G₀ (km/h): mobilisation onset
  gust_ref        = 50,    # G_ref (km/h): saturation speed
  gust_max_score  = 50,    # a: maximum gust component score
  gust_exponent   = 2,     # n: power-law exponent (v² aerodynamic prior)
  # --- Sustained wind transport modifier (S_W) parameters ---
  wind_threshold  = 20,    # W₀ (km/h): onset of transport enhancement
  wind_ref        = 55,    # W_ref (km/h): speed at which multiplier is maximal
  wind_mult_max   = 2.0,   # S_max: maximum multiplier (doubles gust score)
  # --- Rainfall hard gate (R) parameter ---
  rain_threshold  = 0.5,   # P_rain (mm): hourly precipitation for full suppression
  # --- Soil moisture wetness gate (V) parameters ---
  soil_wet        = 0.20,  # SM_wet (m³/m³): too wet for mobilisation
  soil_dry        = 0.05,  # SM_dry (m³/m³): fully dry
  # --- Atmospheric instability bonus (B) parameters ---
  cape_ref        = 1000,  # C_ref (J/kg): CAPE at which bonus is maximal
  instability_max = 0.2,   # B_max: maximum instability bonus (20%)
  # --- Directional attenuation factor (M) and receptor flag (D) parameters ---
  offwind_factor  = 0.5,   # M_offwind: attenuation when wind is off-receptor
  direction_tol   = 15     # Δθ (degrees): tolerance on each arc edge
) {

  # ---- Input validation: meteorological vectors --------------------------- #

  n_hours <- length(wind_gusts_10m)
  checkmate::assert_numeric(wind_gusts_10m, lower = 0, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(wind_speed_10m, lower = 0, any.missing = FALSE, len = n_hours)
  checkmate::assert_numeric(wind_direction_10m, lower = 0, upper = 360,
                            any.missing = FALSE, len = n_hours)
  checkmate::assert_numeric(precipitation, lower = 0, any.missing = FALSE, len = n_hours)
  checkmate::assert_numeric(soil_moisture_0_to_1cm, lower = 0, upper = 1,
                            any.missing = FALSE, len = n_hours)
  checkmate::assert_numeric(cape, lower = 0, any.missing = FALSE, len = n_hours)
  checkmate::assert_numeric(is_day, lower = 0, upper = 1, any.missing = FALSE, len = n_hours)

  # ---- Input validation: receptor arcs ------------------------------------ #

  checkmate::assert_list(receptor_arcs, min.len = 1)
  valid_compass <- names(LRI_COMPASS_DEGREES)
  for (i in seq_along(receptor_arcs)) {
    arc <- receptor_arcs[[i]]
    checkmate::assert_character(
      arc, len = 2, any.missing = FALSE,
      .var.name = paste0("receptor_arcs[[", i, "]]")
    )
    invalid <- arc[!arc %in% valid_compass]
    if (length(invalid) > 0) {
      arg_name <- paste0("receptor_arcs[[", i, "]]")
      cli::cli_abort(c(
        "Invalid compass label in {.arg {arg_name}}: {.val {invalid}}.",
        "i" = "Valid labels are: {.val {valid_compass}}."
      ))
    }
  }

  # ---- Input validation: model parameters --------------------------------- #

  checkmate::assert_number(gust_threshold, lower = 0)
  checkmate::assert_number(gust_ref)
  if (gust_ref <= gust_threshold) {
    cli::cli_abort(c(
      "{.arg gust_ref} ({gust_ref}) must be greater than {.arg gust_threshold} ({gust_threshold}).",
      "i" = "The gust function denominator (gust_ref - gust_threshold) must be positive."
    ))
  }
  checkmate::assert_number(gust_max_score, lower = 0)
  checkmate::assert_number(gust_exponent, lower = 0)
  checkmate::assert_number(wind_threshold, lower = 0)
  checkmate::assert_number(wind_ref)
  if (wind_ref <= wind_threshold) {
    cli::cli_abort(c(
      "{.arg wind_ref} ({wind_ref}) must be greater than {.arg wind_threshold} ({wind_threshold}).",
      "i" = "The sustained wind modifier denominator (wind_ref - wind_threshold) must be positive."
    ))
  }
  checkmate::assert_number(wind_mult_max, lower = 1)
  checkmate::assert_number(rain_threshold, lower = 0)
  checkmate::assert_number(soil_wet, lower = 0, upper = 1)
  checkmate::assert_number(soil_dry, lower = 0)
  if (soil_dry >= soil_wet) {
    cli::cli_abort(c(
      "{.arg soil_dry} ({soil_dry}) must be less than {.arg soil_wet} ({soil_wet}).",
      "i" = "The soil moisture gate denominator (soil_wet - soil_dry) must be positive."
    ))
  }
  checkmate::assert_number(cape_ref, lower = .Machine$double.eps)
  checkmate::assert_number(instability_max, lower = 0, upper = 1)
  checkmate::assert_number(offwind_factor, lower = 0, upper = 1)
  checkmate::assert_number(direction_tol, lower = 0, upper = 90)

  # ---- Sensitive-receptor direction flag D -------------------------------- #

  arc_alpha_deg <- vapply(
    receptor_arcs,
    function(arc) LRI_COMPASS_DEGREES[[arc[1]]],
    numeric(1)
  )
  arc_beta_deg <- vapply(
    receptor_arcs,
    function(arc) LRI_COMPASS_DEGREES[[arc[2]]],
    numeric(1)
  )

  D <- vapply(wind_direction_10m, function(theta) {
    in_any_arc <- any(mapply(
      function(alpha, beta) {
        .lri_direction_in_arc(theta, alpha, beta, direction_tol)
      },
      arc_alpha_deg, arc_beta_deg
    ))
    as.integer(in_any_arc)
  }, integer(1))

  # ---- f_G: Gust mobilisation function [0, gust_max_score] --------------- #

  f_G <- dplyr::case_when(
    wind_gusts_10m <  gust_threshold ~ 0.0,
    wind_gusts_10m >= gust_ref       ~ gust_max_score,
    TRUE ~ gust_max_score *
      ((wind_gusts_10m - gust_threshold) / (gust_ref - gust_threshold))^gust_exponent
  )

  # ---- S_W: Sustained wind transport modifier [1.0, wind_mult_max] -------- #

  S_W <- 1.0 + (wind_mult_max - 1.0) *
    pmin(1.0, pmax(0.0, wind_speed_10m - wind_threshold) /
           (wind_ref - wind_threshold))

  # ---- R: Rainfall hard gate {0, 1} --------------------------------------- #

  R_gate <- dplyr::if_else(precipitation >= rain_threshold, 0.0, 1.0)

  # ---- V: Soil moisture wetness gate [0.0, 1.0] -------------------------- #

  V <- dplyr::case_when(
    soil_moisture_0_to_1cm >= soil_wet ~ 0.0,
    soil_moisture_0_to_1cm <= soil_dry ~ 1.0,
    TRUE ~ (soil_wet - soil_moisture_0_to_1cm) / (soil_wet - soil_dry)
  )

  # ---- B: Atmospheric instability bonus [0.0, instability_max] ------------ #

  B <- dplyr::if_else(
    is_day == 1,
    instability_max * pmin(1.0, cape / cape_ref),
    0.0
  )

  # ---- M: Directional attenuation factor {offwind_factor, 1.0} ------------ #

  M <- dplyr::if_else(D == 1L, 1.0, offwind_factor)

  # ---- Composite LRI ------------------------------------------------------- #

  lri_raw <- f_G * S_W * R_gate * V * (1 + B) * M

  lri <- pmin(lri_raw, 100)

  lri
}


# ---- Internal helpers ------------------------------------------------------ #

# Compass label → degrees lookup (meteorological convention: direction wind
# blows FROM). Eight principal compass points only.
LRI_COMPASS_DEGREES <- c(
  "N"  =   0,
  "NE" =  45,
  "E"  =  90,
  "SE" = 135,
  "S"  = 180,
  "SW" = 225,
  "W"  = 270,
  "NW" = 315
)

# Test whether a wind direction theta (degrees) falls within a receptor arc
# [alpha, beta] (degrees, clockwise), expanded by ±delta_theta on each edge.
.lri_direction_in_arc <- function(theta, alpha, beta, delta_theta) {
  alpha_exp <- (alpha - delta_theta) %% 360
  beta_exp  <- (beta  + delta_theta) %% 360

  if (alpha_exp <= beta_exp) {
    theta >= alpha_exp & theta <= beta_exp
  } else {
    theta >= alpha_exp | theta <= beta_exp
  }
}


#' Litter Risk Index for landfill operations (tibble API)
#'
#' Computes an hourly Litter Risk Index in the range \[0, 100\] for each timestep
#' in a meteorological forecast tibble. Wraps \code{litter_risk_index()} and
#' accepts a tibble with named columns rather than individual vectors.
#'
#' @param met_data A tibble (or data frame) with one row per hourly timestep,
#'   containing at minimum:
#'   \describe{
#'     \item{\code{wind_gusts_10m}}{Peak gust at 10 m (km/h).}
#'     \item{\code{wind_speed_10m}}{Mean wind speed at 10 m (km/h).}
#'     \item{\code{wind_direction_10m}}{Wind direction at 10 m (degrees, FROM).}
#'     \item{\code{precipitation}}{Hourly precipitation (mm).}
#'     \item{\code{soil_moisture_0_to_1cm}}{Volumetric soil moisture (m³/m³).}
#'     \item{\code{cape}}{Convective Available Potential Energy (J/kg).}
#'     \item{\code{is_day}}{Daytime flag (0 or 1).}
#'   }
#'
#' @param receptor_arcs A list of \code{character(2)} vectors specifying the
#'   compass arcs toward sensitive receptors. See \code{\link{litter_risk_index}}
#'   for full details.
#'
#' @param ... Additional arguments passed to \code{\link{litter_risk_index}}
#'   (e.g. \code{gust_threshold}, \code{soil_wet}, \code{offwind_factor}).
#'
#' @return Numeric vector of length \code{nrow(met_data)} giving the Litter
#'   Risk Index in \eqn{[0, 100]} for each forecast hour.
#'
#' @export
generate_litter_risk_index <- function(met_data, receptor_arcs, ...) {
  checkmate::assert_data_frame(met_data, min.rows = 1)

  required_cols <- c(
    "wind_gusts_10m", "wind_speed_10m", "wind_direction_10m",
    "precipitation", "soil_moisture_0_to_1cm", "cape", "is_day"
  )
  missing_cols <- setdiff(required_cols, names(met_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "{.arg met_data} is missing required columns: {.val {missing_cols}}.",
      "i" = paste0(
        "Required: wind_gusts_10m (km/h), wind_speed_10m (km/h), ",
        "wind_direction_10m (degrees), precipitation (mm), ",
        "soil_moisture_0_to_1cm (m\u00b3/m\u00b3), cape (J/kg), is_day (0/1)."
      )
    ))
  }

  litter_risk_index(
    wind_gusts_10m         = met_data$wind_gusts_10m,
    wind_speed_10m         = met_data$wind_speed_10m,
    wind_direction_10m     = met_data$wind_direction_10m,
    precipitation          = met_data$precipitation,
    soil_moisture_0_to_1cm = met_data$soil_moisture_0_to_1cm,
    cape                   = met_data$cape,
    is_day                 = met_data$is_day,
    receptor_arcs          = receptor_arcs,
    ...
  )
}
