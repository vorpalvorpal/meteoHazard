#' Litter hazard index for windblown litter at landfill sites
#'
#' Computes an hourly **meteorological litter hazard**: given the weather in a
#' forecast hour, how strong is the propensity for loose litter to be entrained
#' from the working face and moved. This is a point-source **relative** index,
#' **not** a dispersion model, and it is direction-agnostic — wind *direction*,
#' site geometry, and barriers are handled separately by [litter_exposure()].
#'
#' @section Status — relative index, no fixed operational scale (issue #11):
#' The hard `min(., 100)` cap and the package's fixed `categorise_litter()`
#' tiers were removed: this returns the relative entrainment x transport index
#' directly. With the default `entrainment_max = 50` and `transport_max = 2` the
#' index spans ~0-100 by construction, but that scaling is a modelling choice,
#' not a calibrated operational scale. Mapping it onto a site-specific index /
#' tiers is a calibration step delivered by forthcoming calibration tooling
#' (issues #11/#8).
#'
#' The function does not query any API. The caller fetches the required hourly
#' variables from the Open-Meteo `/v1/forecast` endpoint (winds in m/s,
#' `&wind_speed_unit=ms`) and passes them as numeric vectors, one element per
#' forecast hour.
#'
#' @section Units:
#' The dimensional inputs (`wind_gusts_10m`, `wind_speed_10m`, `precipitation`)
#' may be supplied either as bare numerics in the documented unit or as
#' \pkg{units} objects, which are converted automatically (a dimensionally
#' incompatible unit is an error). `soil_moisture_0_to_1cm` and `wetness` are
#' dimensionless ratios and are taken as-is. The returned index is
#' dimensionless and is a plain numeric.
#'
#' @section Model:
#' The index is a multiplicative combination of entrainment, transport
#' potential, and a rainfall gate:
#' \deqn{LRI = E \cdot T \cdot R}
#'
#' \describe{
#'   \item{Entrainment `E`}{A **phenomenological mobilization-threshold CDF**,
#'     not a soil-physics-derived flux: as the gust rises past a threshold, an
#'     increasing *fraction* of the litter population on the working face is
#'     recruited into motion (dry film/bags first, at modest gusts; heavier
#'     items such as bottles need much stronger gusts, if they move at all).
#'     The AP-42 (EPA 2006) / Bagnold-type excess-power form
#'     \eqn{E = a\,\min(1, ((G-G_t)/\Delta G)^n)} is borrowed for its *shape*
#'     only — a threshold-gated, saturating ramp is the right qualitative
#'     behaviour for a mobilization CDF — and is evaluated directly on the
#'     10 m **gust** `G` (m/s) rather than on a soil friction velocity `u*`.
#'     Macro-litter is a different transport regime from mineral dust/sand
#'     (item shape, area-to-mass ratio, and lofting mechanics all differ), so
#'     `entrainment_max`, `gust_threshold`, and `gust_reference` are
#'     litter-calibration placeholders, not soil-erosion constants transferred
#'     from AP-42/Shao & Lu: Mellink et al. (2024) report near-100% mobility
#'     for bags vs. near-0% for bottles at the *same* ~2.3 m/s wind speed —
#'     material identity, not a single physical threshold, dominates: see
#'     Valyrakis et al. (2024) on energetic-flow-structure-driven entrainment
#'     of sediment and plastic debris, and [litter_wetness()] / `material` for
#'     how this package's threshold/saturation terms attempt to track it.}
#'   \item{Moisture-raised threshold}{A damp surface needs a stronger gust to
#'     release litter (Fecan et al. 1999 supplies the *shape*): the threshold
#'     rises with normalised surface wetness, \eqn{G_t = G_{t0}(1 + \gamma
#'     s^{\beta})}, where `s` is `wetness` directly when supplied, else the
#'     Fecan-style clamp of `soil_moisture_0_to_1cm` (see
#'     `soil_dry`/`soil_wet`).}
#'   \item{Material-aware graded saturation (R-B3a)}{At or above the
#'     saturation mark (`s == 1`) the surface is treated as fully wet. Loose
#'     film/thin plastic (`material = "film"`, the default — Mellink et al.
#'     2024 identify film/bags as the dominant landfill-escape material) can
#'     still be partially entrained off a saturated surface (residual
#'     entrainment `E0 * (1 - saturation_penalty)`), whereas absorbent paper
#'     (`material = "paper"`) is hard-vetoed to zero once saturated — it has
#'     soaked up water rather than sitting on top of a film. The graded
#'     penalty only touches the saturated region; the unsaturated
#'     threshold-rise ramp below `s == 1` is identical for both materials.}
#'   \item{Transport potential `T`}{Driven by the **mean wind** directly (not
#'     the gust; transport is flight-height advection). A linear ramp from 1 to
#'     `transport_max` between `wind_transport_onset` and `wind_transport_ref`
#'     m/s — the "how far is it moved" penalty. Mellink et al. (2024) observed
#'     macroplastic transport velocity rising with wind speed (~141 to 219 m/h
#'     over 2.3--3.2 m/s), supporting a monotone transport-vs-mean-wind ramp;
#'     the onset/reference wind speeds are uncalibrated engineering placeholders
#'     (carried from v2), not fitted constants, and the reach/geometry mapping
#'     of this potential is deferred to [litter_exposure()].}
#'   \item{Rainfall gate `R`}{A smooth ramp (R-B4), not a hard cutoff: `R = 1`
#'     at or below `rain_onset`, `R = 0` at or above `rain_threshold`, linear
#'     between. Windblown/aeolian transport only — this package does not model
#'     rain-driven runoff transport of litter (Mellink et al. 2024 show runoff
#'     is a materially different, non-aeolian mobilisation pathway).}
#' }
#' Any single suppressor (rain, saturated + paper, sub-threshold gust) drives
#' the index to zero. The maximum attainable value is
#' `entrainment_max * transport_max` (= 100 with the defaults 50 and 2); this
#' is not a fixed cap (R-A3) — an inflated `entrainment_max` legitimately
#' exceeds it.
#'
#' @param wind_gusts_10m Numeric vector. Peak wind gust at 10 m (m/s),
#'   Open-Meteo `wind_gusts_10m` (fetch with `&wind_speed_unit=ms`). Drives
#'   entrainment.
#' @param wind_speed_10m Numeric vector. Mean wind speed at 10 m (m/s),
#'   Open-Meteo `wind_speed_10m`. Drives transport potential.
#' @param precipitation Numeric vector. Hourly precipitation (mm), Open-Meteo
#'   `precipitation`. Feeds the smooth rainfall gate.
#' @param soil_moisture_0_to_1cm Numeric vector or `NULL` (default). Volumetric
#'   water content of the 0--1 cm layer (m^3/m^3), Open-Meteo
#'   `soil_moisture_0_to_1cm`. Required **unless** `wetness` is supplied; used
#'   as a relative surface-wetness surrogate that raises the entrainment
#'   threshold and supplies the saturation mark.
#' @param wetness Numeric vector or `NULL` (default). Precomputed
#'   litter-surface wetness in `[0, 1]` from [litter_wetness_vec()] (a faster,
#'   litter-specific alternative to the gridded soil-moisture proxy — see
#'   `use_wetness_state` in [litter_hazard()]). When supplied it **supersedes**
#'   `soil_moisture_0_to_1cm` (ADJ-6); supplying both warns
#'   (`meteoHazard_litter_wetness_supersedes`).
#' @param material Character. `"film"` (default) or `"paper"`. Governs how a
#'   saturated surface is treated — see @section Model, R-B3a.
#' @param gust_threshold Dry threshold gust (m/s) below which entrainment is
#'   zero. Default `3.9737355`. Replaces the v3.0 `kappa`/`z0`/`ustar_t0`
#'   friction-velocity parameterisation: `u*` is a linear rescaling of the
#'   gust under the neutral log wind profile, so those four knobs collapse to
#'   two identifiable gust speeds without losing any behaviour (needed for the
#'   calibration workflow in issue #26). The default reproduces the old
#'   `ustar_t0 = 0.30` m/s threshold to rounding
#'   (`0.30 / (0.40 / log(200))` = 3.97373802; pinned as 3.9737355).
#' @param gust_reference Gust (m/s) at which entrainment saturates. Default
#'   `13.9080743`. Must exceed `gust_threshold`. Reproduces the old
#'   `ustar_ref = 1.05` m/s to rounding
#'   (`1.05 / (0.40 / log(200))` = 13.90808309).
#' @param entrainment_max Maximum entrainment score. Default 50.
#' @param excess_exponent Power-law exponent on the gust excess. Default 2
#'   (AP-42-shaped); calibratable in `[2, 3]`.
#' @param moisture_gain Maximum fractional increase of the threshold as the
#'   surface approaches wet. Default 2.0. High calibration uncertainty.
#' @param moisture_curve Curvature of the moisture-threshold rise. Default 0.5
#'   (concave, Fecan-type). High calibration uncertainty.
#' @param soil_dry Soil moisture at or below which the surface is fully dry
#'   (m^3/m^3). Default 0.05. Only used on the `soil_moisture_0_to_1cm` path.
#' @param soil_wet Soil moisture at or above which the normalised wetness `s`
#'   reaches 1 (m^3/m^3); also the "saturated" mark on the soil-moisture path.
#'   Default 0.20. Must exceed `soil_dry`.
#' @param saturation_penalty Fractional entrainment loss on a saturated `film`
#'   surface, so residual entrainment there is `E0 * (1 - saturation_penalty)`.
#'   Default 0.7 (i.e. 70% reduction, 30% residual).
#' @param paper_veto_wetness Normalised wetness (`wetness`, when supplied) at
#'   or above which a `paper` surface is hard-vetoed (`E = 0`). Default 0.8.
#' @param wind_transport_onset Mean wind below which transport adds nothing
#'   (m/s). Default 5.5 (~20 km/h).
#' @param wind_transport_ref Mean wind at which transport saturates (m/s).
#'   Default 15 (~54 km/h). Must exceed `wind_transport_onset`.
#' @param transport_max Maximum transport multiplier. Default 2.0.
#' @param rain_onset Hourly precipitation (mm) at or below which the rainfall
#'   gate is fully open (`R = 1`). Default 0.2.
#' @param rain_threshold Hourly precipitation (mm) at or above which the
#'   rainfall gate is fully closed (`R = 0`). Default 0.5. Must exceed
#'   `rain_onset`.
#'
#' @return Numeric vector of length equal to the inputs, the relative litter
#'   hazard index for each forecast hour (dimensionless; ~0-100 by construction
#'   with the default `entrainment_max`/`transport_max`, but not a calibrated
#'   operational scale — see the Status section).
#'
#' @references
#' EPA (2006). AP-42, Section 13.2.5: Industrial Wind Erosion. Supplies the
#' threshold-gated excess-power *shape* of the entrainment CDF (not its
#' numbers — see @section Model).
#'
#' Fecan, F., Marticorena, B. and Bergametti, G. (1999). Parametrization of the
#' increase of the aeolian erosion threshold wind friction velocity due to soil
#' moisture for arid and semi-arid areas. \emph{Annales Geophysicae}, 17,
#' 149--157. \doi{10.1007/s00585-999-0149-7}. Basis for the moisture-raised
#' threshold *shape*.
#'
#' Mellink, Y. A. M., van Emmerik, T. H. M. & Mani, T. (2024). Wind- and
#' rain-driven macroplastic mobilization and transport on land.
#' \emph{Scientific Reports} 14, 3898. \doi{10.1038/s41598-024-53971-8}.
#'
#' Valyrakis, M. et al. (2024). The role of energetic flow structures on the
#' aeolian transport of sediment and plastic debris. \emph{Acta Mechanica
#' Sinica}. \doi{10.1007/s10409-024-24467-x}.
#'
#' @seealso [litter_exposure()] for the direction- and geometry-aware exposure
#'   layer that sits on top of this hazard index; [litter_wetness_vec()] for
#'   the litter-specific wetness state consumed via `wetness`.
#' @export
litter_hazard_vec <- function(
  wind_gusts_10m,
  wind_speed_10m,
  precipitation,
  soil_moisture_0_to_1cm = NULL,
  wetness              = NULL,
  material             = c("film", "paper"),
  gust_threshold       = 3.9737355,
  gust_reference       = 13.9080743,
  entrainment_max      = 50,
  excess_exponent      = 2,
  moisture_gain        = 2.0,
  moisture_curve       = 0.5,
  soil_dry             = 0.05,
  soil_wet             = 0.20,
  saturation_penalty   = 0.7,
  paper_veto_wetness   = 0.8,
  wind_transport_onset = 5.5,
  wind_transport_ref   = 15.0,
  transport_max        = 2.0,
  rain_onset           = 0.2,
  rain_threshold       = 0.5
) {

  # ---- Normalise dimensional inputs (bare = documented unit; units = converted) #
  # A units-tagged input is converted to the canonical unit (erroring on a
  # dimensional mismatch); a bare numeric is assumed already in that unit. The
  # gate/excess arithmetic below then runs on plain canonical-unit doubles.
  # soil_moisture_0_to_1cm and wetness are dimensionless ratios and stay plain.
  wind_gusts_10m <- .drop_to(wind_gusts_10m, "m/s", arg = "wind_gusts_10m")
  wind_speed_10m <- .drop_to(wind_speed_10m, "m/s", arg = "wind_speed_10m")
  precipitation  <- .drop_to(precipitation,  "mm",  arg = "precipitation")

  # ---- Validate meteorological inputs (complete, non-negative, aligned) ---- #
  n <- length(wind_gusts_10m)
  checkmate::assert_numeric(wind_gusts_10m, lower = 0, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(wind_speed_10m, lower = 0, any.missing = FALSE, len = n)
  checkmate::assert_numeric(precipitation, lower = 0, any.missing = FALSE, len = n)

  material <- match.arg(material)

  # ---- ADJ-4/ADJ-6: soil_moisture_0_to_1cm / wetness contract seam --------- #
  # Exactly one of the two surface-wetness inputs is required. wetness is the
  # faster litter-specific state (R-B3b); soil_moisture_0_to_1cm is the
  # gridded (multi-day-memory) proxy used historically. If both are supplied,
  # wetness wins -- it is the more physically appropriate driver for a litter
  # surface -- and the caller is warned so a stray leftover argument doesn't
  # silently get ignored.
  have_sm  <- !is.null(soil_moisture_0_to_1cm)
  have_wet <- !is.null(wetness)
  if (!have_sm && !have_wet) {
    cli::cli_abort(
      c("Exactly one of {.arg soil_moisture_0_to_1cm} or {.arg wetness} must be supplied.",
        "i" = "Use {.arg soil_moisture_0_to_1cm} (basic mode) or {.arg wetness} (litter-wetness-state mode, see litter_wetness_vec())."),
      class = "meteoHazard_input_error"
    )
  }
  if (have_sm && have_wet) {
    cli::cli_warn(
      c("Both {.arg soil_moisture_0_to_1cm} and {.arg wetness} were supplied.",
        "i" = "{.arg wetness} supersedes {.arg soil_moisture_0_to_1cm} for this call (ADJ-6)."),
      class = "meteoHazard_litter_wetness_supersedes"
    )
    have_sm <- FALSE
  }

  if (have_wet) {
    checkmate::assert_numeric(wetness, lower = 0, upper = 1, any.missing = FALSE, len = n)
  } else {
    checkmate::assert_numeric(soil_moisture_0_to_1cm, lower = 0, upper = 1,
                              any.missing = FALSE, len = n)
  }

  # ---- Validate parameters and the ordering constraints --------------------- #
  checkmate::assert_number(gust_threshold, lower = 0)
  checkmate::assert_number(gust_reference)
  if (gust_reference <= gust_threshold) {
    cli::cli_abort(c(
      "{.arg gust_reference} ({gust_reference}) must be greater than {.arg gust_threshold} ({gust_threshold}).",
      "i" = "The entrainment excess range (gust_reference - gust_threshold) must be positive."
    ), class = "meteoHazard_input_error")
  }
  checkmate::assert_number(entrainment_max, lower = 0)
  checkmate::assert_number(excess_exponent, lower = 0)
  checkmate::assert_number(moisture_gain, lower = 0)
  checkmate::assert_number(moisture_curve, lower = 0)
  checkmate::assert_number(soil_dry, lower = 0)
  checkmate::assert_number(soil_wet, lower = 0)
  if (soil_dry >= soil_wet) {
    cli::cli_abort(c(
      "{.arg soil_dry} ({soil_dry}) must be less than {.arg soil_wet} ({soil_wet}).",
      "i" = "The moisture ramp denominator (soil_wet - soil_dry) must be positive."
    ), class = "meteoHazard_input_error")
  }
  checkmate::assert_number(saturation_penalty, lower = 0, upper = 1)
  checkmate::assert_number(paper_veto_wetness, lower = 0, upper = 1)
  checkmate::assert_number(wind_transport_onset, lower = 0)
  checkmate::assert_number(wind_transport_ref)
  if (wind_transport_ref <= wind_transport_onset) {
    cli::cli_abort(c(
      "{.arg wind_transport_ref} ({wind_transport_ref}) must be greater than {.arg wind_transport_onset} ({wind_transport_onset}).",
      "i" = "The transport ramp denominator must be positive."
    ), class = "meteoHazard_input_error")
  }
  checkmate::assert_number(transport_max, lower = 1)
  checkmate::assert_number(rain_onset, lower = 0)
  checkmate::assert_number(rain_threshold)
  if (rain_threshold <= rain_onset) {
    cli::cli_abort(c(
      "{.arg rain_threshold} ({rain_threshold}) must be greater than {.arg rain_onset} ({rain_onset}).",
      "i" = "The rain ramp denominator (rain_threshold - rain_onset) must be positive."
    ), class = "meteoHazard_input_error")
  }

  # ---- Normalised surface wetness `s` and the saturation mark -------------- #
  # R-B2/entrainment computed directly in gust units (see @param gust_threshold):
  # u* is a linear rescaling of the gust under the neutral log wind profile
  # (u* = kappa*G/ln(z/z0)), so re-deriving u*t/u*ref from kappa/z0/ustar_t0/
  # ustar_ref and then dividing back out is equivalent to working in gust
  # units throughout -- two identifiable knobs instead of four (issue #26).
  #
  # w_norm = wetness directly when supplied (already normalised [0,1]);
  # otherwise the Fecan-style clamp of soil_moisture_0_to_1cm onto [0,1]
  # between soil_dry and soil_wet. "saturated" marks the region where the
  # material-aware graded-saturation treatment (R-B3a) applies: wetness >=
  # paper_veto_wetness on the wetness path (a dedicated, independently
  # calibratable mark), or SM >= soil_wet on the soil-moisture path (ADJ-2:
  # reusing soil_wet keeps the pre-existing soil-moisture worked oracles
  # unchanged).
  if (have_wet) {
    w_norm    <- wetness
    saturated <- wetness >= paper_veto_wetness
  } else {
    w_norm    <- pmin(1, pmax(0, (soil_moisture_0_to_1cm - soil_dry) / (soil_wet - soil_dry)))
    saturated <- soil_moisture_0_to_1cm >= soil_wet
  }

  # ---- Moisture-raised entrainment threshold (Fecan et al. 1999 shape) ----- #
  # A damp surface raises the threshold gust; a dry surface (w_norm == 0)
  # sits at gust_threshold itself.
  gust_t <- gust_threshold * (1 + moisture_gain * w_norm^moisture_curve)

  # ---- Entrainment: phenomenological mobilization-threshold CDF ------------ #
  # E0 = a * min(1, (max(0, G - Gt) / (Gref - Gt0))^n). The excess is floored
  # at 0 before exponentiation so a fractional exponent never meets a negative
  # base. See @section Model for why this is litter-calibrated, not a
  # soil-transport flux (Valyrakis 2024; Mellink 2024).
  denom <- gust_reference - gust_threshold
  excess <- pmax(0, wind_gusts_10m - gust_t)
  E0 <- entrainment_max * pmin(1, (excess / denom)^excess_exponent)

  # ---- R-B3a/ADJ-2/ADJ-3: material-aware graded saturation ----------------- #
  # Only touches the saturated region (s == 1); the unsaturated threshold-rise
  # ramp above is identical for both materials. paper: hard veto (it has
  # absorbed the water, not just sat under a film). film (default -- ADJ-3,
  # the dominant landfill-escape material per Mellink et al. 2024): a graded
  # residual, entrainment_max scaled down by saturation_penalty rather than
  # zeroed, because a film can still be picked up off a wet surface.
  E <- E0
  if (material == "paper") {
    E[saturated] <- 0
  } else {
    E[saturated] <- E0[saturated] * (1 - saturation_penalty)
  }

  # ---- Transport potential from the mean wind ------------------------------ #
  # Linear ramp [1, transport_max] over the mean wind (m/s). Transport is
  # flight-height advection, so it is driven by the mean wind, not the gust.
  # Mellink et al. (2024, Sci. Rep. 14:3898, doi:10.1038/s41598-024-53971-8)
  # measured macroplastic transport velocity rising with wind speed (~141->219
  # m/h over 2.3-3.2 m/s), supporting a monotone transport-vs-mean-wind ramp; the
  # onset/reference wind speeds themselves are uncalibrated v2 placeholders.
  transport <- 1 + (transport_max - 1) *
    pmin(1, pmax(0, wind_speed_10m - wind_transport_onset) /
           (wind_transport_ref - wind_transport_onset))

  # ---- R-B4/ADJ-5: smooth rainfall gate ------------------------------------- #
  # A linear ramp from fully open (R=1, at/below rain_onset) to fully closed
  # (R=0, at/above rain_threshold), replacing the v3.0 hard cutoff at a single
  # rain_threshold -- a light shower does not instantaneously suppress litter
  # movement the way a genuinely wetting rain event does. Scope: windblown/
  # aeolian transport only; rain-driven runoff transport of litter (a
  # materially different pathway -- Mellink et al. 2024) is out of scope.
  rain_gate <- 1 - pmin(1, pmax(0, (precipitation - rain_onset) / (rain_threshold - rain_onset)))

  # ---- Composite relative index -------------------------------------------- #
  # E in [0, entrainment_max], transport in [1, transport_max], rain_gate in
  # [0, 1]: a bounded-by-construction relative screening index (with the
  # default entrainment_max=50 / transport_max=2 it spans ~0-100), but NOT a
  # fixed 0-100 cap (R-A3/issue #11) -- an inflated entrainment_max legitimately
  # exceeds 100. Site-specific tiers come from calibration tooling (issues
  # #11/#8).
  E * transport * rain_gate
}


#' Litter hazard index for a forecast tibble
#'
#' Computes the hourly litter hazard index ([litter_hazard_vec()]) for each row
#' of a meteorological forecast tibble. Wraps the vector API, accepting a tibble
#' with named columns rather than individual vectors. Two column contracts are
#' supported, selected by `use_wetness_state`.
#'
#' @param met_data A tibble (or data frame) with one row per hourly timestep.
#'   Always requires `wind_gusts_10m` (m/s), `wind_speed_10m` (m/s), and
#'   `precipitation` (mm). When `use_wetness_state = FALSE` (default), also
#'   requires `soil_moisture_0_to_1cm` (m^3/m^3). When `use_wetness_state =
#'   TRUE`, also requires the five columns consumed by [litter_wetness_vec()]:
#'   `temperature_2m` (degC), `relative_humidity_2m` (%), `shortwave_radiation`
#'   (W/m^2), `wind_speed_10m`, and `precipitation` (the latter two already
#'   required above).
#' @param use_wetness_state Logical. Default `FALSE` (uses
#'   `soil_moisture_0_to_1cm` directly, as in earlier package versions). When
#'   `TRUE`, the litter-specific wetness state ([litter_wetness_vec()], run
#'   with its own defaults -- see ADJ-8) is computed from `met_data` first and
#'   passed through as `wetness`, in place of `soil_moisture_0_to_1cm`.
#' @param ... Additional calibration parameters forwarded to
#'   [litter_hazard_vec()] ONLY (e.g. `rain_threshold`, `soil_wet`,
#'   `gust_threshold`, `material`). ADJ-8: `...` is not forwarded to the
#'   internal [litter_wetness_vec()] call in the `use_wetness_state = TRUE`
#'   path -- that call always uses its own defaults, so a hazard-side
#'   calibration override (e.g. a custom `rain_threshold` for the hazard's
#'   rain gate) cannot accidentally perturb the wetness-reset dynamics.
#'
#' @return Numeric vector of length `nrow(met_data)`, the relative litter hazard
#'   index for each forecast hour (see [litter_hazard_vec()]). Issue #11 removed
#'   the fixed 0-100 scale across hazards in favour of physical/relative outputs.
#'
#' @seealso [litter_hazard_vec()], [litter_wetness_vec()], [litter_exposure()].
#' @export
litter_hazard <- function(met_data, use_wetness_state = FALSE, ...) {
  checkmate::assert_data_frame(met_data, min.rows = 1)
  checkmate::assert_flag(use_wetness_state)

  always_required <- c("wind_gusts_10m", "wind_speed_10m", "precipitation")

  if (use_wetness_state) {
    wetness_required <- c(
      "temperature_2m", "relative_humidity_2m", "shortwave_radiation",
      "wind_speed_10m", "precipitation"
    )
    required_cols <- union(always_required, wetness_required)
    .assert_required_cols(
      met_data, required_cols, arg = "met_data",
      info = paste0(
        "use_wetness_state = TRUE requires wind_gusts_10m (m/s), ",
        "wind_speed_10m (m/s), precipitation (mm), temperature_2m (degC), ",
        "relative_humidity_2m (%), and shortwave_radiation (W/m²)."
      )
    )

    # ADJ-8: litter_wetness_vec() runs with its own defaults here; `...` below
    # is forwarded to litter_hazard_vec() only.
    wetness <- litter_wetness_vec(
      precipitation         = met_data$precipitation,
      temperature_2m        = met_data$temperature_2m,
      relative_humidity_2m  = met_data$relative_humidity_2m,
      wind_speed_10m        = met_data$wind_speed_10m,
      shortwave_radiation   = met_data$shortwave_radiation
    )

    litter_hazard_vec(
      wind_gusts_10m = met_data$wind_gusts_10m,
      wind_speed_10m = met_data$wind_speed_10m,
      precipitation  = met_data$precipitation,
      wetness        = wetness,
      ...
    )
  } else {
    required_cols <- c(always_required, "soil_moisture_0_to_1cm")
    .assert_required_cols(
      met_data, required_cols, arg = "met_data",
      info = paste0(
        "Required: wind_gusts_10m (m/s), wind_speed_10m (m/s), ",
        "precipitation (mm), soil_moisture_0_to_1cm (m³/m³)."
      )
    )

    litter_hazard_vec(
      wind_gusts_10m         = met_data$wind_gusts_10m,
      wind_speed_10m         = met_data$wind_speed_10m,
      precipitation          = met_data$precipitation,
      soil_moisture_0_to_1cm = met_data$soil_moisture_0_to_1cm,
      ...
    )
  }
}
