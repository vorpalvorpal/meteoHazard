#' Calculate Thermal Work Limit (TWL)
#'
#' Calculates the maximum sustainable metabolic rate (W/m^2) that
#' well-hydrated, acclimatised individuals can maintain in a specific
#' thermal environment. Based on the Brake & Bates (2002) methodology.
#'
#' Weather parameters (`temp`, `wind_speed`, `RH`, `direct_solar`,
#' `diffuse_solar`, `pressure`) that are not supplied are automatically
#' retrieved from the [Open-Meteo API](https://open-meteo.com/) using the
#' provided `datetime`, `latitude` and `longitude`. An internet connection
#' is required for auto-fetching.
#'
#' Open-Meteo returns wind speed at 10 m height (`wind_speed_10m`). This is
#' corrected to approximately 1 m (body level) using a logarithmic wind
#' profile with roughness length z0 = 0.01 m (open terrain), giving a
#' correction factor of `ln(1/z0) / ln(10/z0)` ≈ 0.667. Supply `wind_speed`
#' directly (already at body level) to skip this correction.
#'
#' @param datetime POSIXct datetime vector (required).
#' @param latitude Latitude in decimal degrees (required).
#' @param longitude Longitude in decimal degrees (required).
#' @param temp Dry bulb air temperature in degrees Celsius, or `NULL` to
#'   fetch from Open-Meteo.
#' @param wind_speed Wind speed in m/s, or `NULL` to fetch from Open-Meteo.
#' @param RH Relative humidity in percent, or `NULL` to fetch from
#'   Open-Meteo.
#' @param direct_solar Direct beam solar radiation in W/m^2, or `NULL` to
#'   fetch from Open-Meteo.
#' @param diffuse_solar Diffuse sky solar radiation in W/m^2, or `NULL` to
#'   fetch from Open-Meteo.
#' @param pressure Barometric pressure in hPa (or kPa if `convert_pressure =
#'   FALSE`), or `NULL` to fetch from Open-Meteo.
#' @param wet_bulb Natural (unventilated) wet bulb temperature in degrees
#'   Celsius, or `NULL` (the default) to estimate it from `temp`, `RH`,
#'   `pressure`, `wind_speed`, and `globe_temp`.
#' @param albedo Ground albedo (0--1). Default 0.12 for asphalt.
#' @param Icl Intrinsic clothing thermal resistance in clo. Default 0.6.
#' @param icl Clothing vapour permeation efficiency (0--1). Default 0.45.
#' @param max_core_temp Maximum acceptable core temperature in degrees Celsius.
#'   Default 38.2.
#' @param max_sweat_rate Maximum acceptable sweat rate in kg/(m^2.hr).
#'   Default 0.67.
#' @param convert_pressure If `TRUE` (the default), converts **user-supplied**
#'   pressure from hPa to kPa. API-fetched pressure (Open-Meteo
#'   `surface_pressure`) is always converted hPa→kPa regardless of this flag.
#' @param verbose If `TRUE` (the default), show progress and diagnostic
#'   messages.
#'
#' @return Numeric vector of TWL values in W/m^2.
#'
#' @details
#' TWL categories (from Brake & Bates 2002):
#' \itemize{
#'   \item > 220 W/m^2: Unrestricted work
#'   \item 140--220 W/m^2: Acclimatisation zone
#'   \item 115--140 W/m^2: Buffer zone
#'   \item < 115 W/m^2: Withdrawal required
#' }
#'
#' @references
#' Brake, D.J. and Bates, G.P. (2002) Limiting Metabolic Rate (Thermal Work
#' Limit) as an Index of Thermal Stress. Applied Occupational and Environmental
#' Hygiene, 17:3, 176-186.
#'
#' @export
generate_twl <- function(datetime,
                         latitude,
                         longitude,
                         temp = NULL,
                         wind_speed = NULL,
                         RH = NULL,
                         direct_solar = NULL,
                         diffuse_solar = NULL,
                         pressure = NULL,
                         wet_bulb = NULL,
                         albedo = 0.12,
                         Icl = 0.6,
                         icl = 0.45,
                         max_core_temp = 38.2,
                         max_sweat_rate = 0.67,
                         convert_pressure = TRUE,
                         verbose = TRUE) {

  # --- Input validation (runs BEFORE any API call) ---
  # All range checks use na.rm = TRUE so NA values propagate rather than error.
  n <- length(datetime)

  if (!inherits(datetime, "POSIXct") || n < 1L) {
    cli::cli_abort(
      "`datetime` must be a POSIXct vector with at least one element.",
      class = "meteoHazard_input_error"
    )
  }

  # latitude / longitude: numeric, length 1 or n, valid range
  for (arg_name in c("latitude", "longitude")) {
    val  <- get(arg_name)
    vlen <- length(val)
    if (!is.numeric(val) || !(vlen == 1L || vlen == n)) {
      cli::cli_abort(
        "`{arg_name}` must be numeric with length 1 or length(datetime) ({n}).",
        class = "meteoHazard_input_error"
      )
    }
  }
  if (any(latitude  < -90  | latitude  > 90,  na.rm = TRUE)) {
    cli::cli_abort(
      "`latitude` values must be in [-90, 90].",
      class = "meteoHazard_input_error"
    )
  }
  if (any(longitude < -180 | longitude > 180, na.rm = TRUE)) {
    cli::cli_abort(
      "`longitude` values must be in [-180, 180].",
      class = "meteoHazard_input_error"
    )
  }

  # Per-observation weather args: numeric, length 1 or n, range constraints
  weather_args <- list(
    temp          = temp,
    wind_speed    = wind_speed,
    RH            = RH,
    direct_solar  = direct_solar,
    diffuse_solar = diffuse_solar,
    pressure      = pressure,
    wet_bulb      = wet_bulb
  )
  for (arg_name in names(weather_args)) {
    val <- weather_args[[arg_name]]
    if (is.null(val)) next  # NULL means fetch from API; validated later
    vlen <- length(val)
    if (!is.numeric(val) || !(vlen == 1L || vlen == n)) {
      cli::cli_abort(
        "`{arg_name}` must be numeric with length 1 or length(datetime) ({n}).",
        class = "meteoHazard_input_error"
      )
    }
  }
  if (!is.null(RH)           && any(RH           < 0   | RH           > 100, na.rm = TRUE))
    cli::cli_abort("`RH` values must be in [0, 100].",           class = "meteoHazard_input_error")
  if (!is.null(wind_speed)   && any(wind_speed   < 0,                         na.rm = TRUE))
    cli::cli_abort("`wind_speed` values must be >= 0.",           class = "meteoHazard_input_error")
  if (!is.null(direct_solar) && any(direct_solar < 0,                         na.rm = TRUE))
    cli::cli_abort("`direct_solar` values must be >= 0.",         class = "meteoHazard_input_error")
  if (!is.null(diffuse_solar) && any(diffuse_solar < 0,                        na.rm = TRUE))
    cli::cli_abort("`diffuse_solar` values must be >= 0.",        class = "meteoHazard_input_error")
  if (!is.null(pressure)     && any(pressure     <= 0,                         na.rm = TRUE))
    cli::cli_abort("`pressure` values must be > 0.",              class = "meteoHazard_input_error")

  # Scalar parameters
  if (!is.numeric(albedo) || length(albedo) != 1L || albedo < 0 || albedo > 1)
    cli::cli_abort("`albedo` must be a single numeric value in [0, 1].",       class = "meteoHazard_input_error")
  if (!is.numeric(icl)    || length(icl)    != 1L || icl    < 0 || icl    > 1)
    cli::cli_abort("`icl` must be a single numeric value in [0, 1].",          class = "meteoHazard_input_error")
  if (!is.numeric(Icl)    || length(Icl)    != 1L || Icl    < 0)
    cli::cli_abort("`Icl` must be a single numeric value >= 0.",               class = "meteoHazard_input_error")
  if (!is.numeric(max_sweat_rate) || length(max_sweat_rate) != 1L || max_sweat_rate <= 0)
    cli::cli_abort("`max_sweat_rate` must be a single numeric value > 0.",     class = "meteoHazard_input_error")
  if (!is.numeric(max_core_temp)  || length(max_core_temp)  != 1L || !is.finite(max_core_temp))
    cli::cli_abort("`max_core_temp` must be a single finite numeric value.",   class = "meteoHazard_input_error")

  # --- Record whether pressure will be API-fetched (affects unit conversion) ---
  pressure_from_api <- is.null(pressure)

  # --- Determine which weather fields need fetching ---
  # Note: wind_speed_10m from Open-Meteo is at 10 m height; it is corrected
  # to body level (~1 m) below using a log-profile wind correction.
  field_map <- c(
    temp          = "temperature_2m",
    wind_speed    = "wind_speed_10m",
    RH            = "relative_humidity_2m",
    direct_solar  = "direct_radiation",
    diffuse_solar = "diffuse_radiation",
    pressure      = "surface_pressure"
  )

  needs_fetch <- vapply(
    list(temp, wind_speed, RH, direct_solar, diffuse_solar, pressure),
    is.null, logical(1)
  )
  names(needs_fetch) <- names(field_map)
  fields_needed <- field_map[needs_fetch]

  # Track whether wind_speed came from the API (needs height correction)
  wind_from_api <- is.null(wind_speed)

  if (length(fields_needed) > 0L) {
    if (verbose) {
      cli_alert_info(
        "Fetching missing weather data: {paste(names(fields_needed), collapse = ', ')}"
      )
    }
    api_data <- fetch_openmeteo(
      datetime, latitude, longitude, unname(fields_needed), verbose = verbose
    )
    if (is.null(temp))          temp          <- api_data[["temperature_2m"]]
    if (is.null(wind_speed))    wind_speed    <- api_data[["wind_speed_10m"]]
    if (is.null(RH))            RH            <- api_data[["relative_humidity_2m"]]
    if (is.null(direct_solar))  direct_solar  <- api_data[["direct_radiation"]]
    if (is.null(diffuse_solar)) diffuse_solar <- api_data[["diffuse_radiation"]]
    if (is.null(pressure))      pressure      <- api_data[["surface_pressure"]]
  }

  # Recycle any length-1 user scalars to n_obs so all vectors are coherent.
  # API-fetched vectors are already length n; rep_len() is a no-op for them.
  n_obs <- length(datetime)
  temp          <- rep_len(temp,          n_obs)
  wind_speed    <- rep_len(wind_speed,    n_obs)
  RH            <- rep_len(RH,            n_obs)
  direct_solar  <- rep_len(direct_solar,  n_obs)
  diffuse_solar <- rep_len(diffuse_solar, n_obs)
  pressure      <- rep_len(pressure,      n_obs)
  if (!is.null(wet_bulb)) wet_bulb <- rep_len(wet_bulb, n_obs)

  if (verbose) {
    cli_h1("TWL Calculation")
    cli_alert_info("Processing {n_obs} observation{?s}")
  }

  # Constants
  LR     <- 16.5
  lambda <- TWL_CONSTANTS$LATENT_HEAT_TWL_KJ   # 2430 kJ/kg at skin temp ~30 °C [ASHRAE Eq. 14]

  # Unit conversions
  # API-fetched pressure is always in hPa (Open-Meteo surface_pressure) and
  # must be divided by 10 to get kPa.  User-supplied pressure is divided by 10
  # only when convert_pressure = TRUE (i.e. the user passed hPa, not kPa).
  if (pressure_from_api) {
    pressure <- pressure / 10
  } else if (convert_pressure) {
    pressure <- pressure / 10
  }

  # Apply log-profile wind height correction for API-sourced wind speed.
  # Open-Meteo supplies wind_speed_10m (at 10 m); Brake & Bates use wind at
  # body level (~1 m). Log-profile with roughness length z0 = 0.01 m gives:
  #   v_1m = v_10m * ln(1/z0) / ln(10/z0)  =  v_10m * 0.667
  if (wind_from_api) {
    WIND_HEIGHT_FACTOR <- log(1 / 0.01) / log(10 / 0.01)  # ≈ 0.667
    wind_speed <- wind_speed * WIND_HEIGHT_FACTOR
    if (verbose) {
      cli_alert_info(
        "Applied 10 m -> 1 m wind height correction (factor {round(WIND_HEIGHT_FACTOR, 3)})"
      )
    }
  }

  # Constrain wind speed to valid range (Brake & Bates recommend 0.2-4.0 m/s)
  wind_speed_orig <- wind_speed
  wind_speed <- pmax(0.2, pmin(4.0, wind_speed))
  n_capped <- sum(wind_speed_orig != wind_speed, na.rm = TRUE)
  if (verbose && n_capped > 0) {
    cli_alert_warning(
      "{n_capped} wind speed value{?s} constrained to [0.2, 4.0] m/s range"
    )
  }

  # Calculate solar position
  if (verbose) cli_alert("Calculating solar position...")
  solar_pos <- calculate_solar_position(datetime, latitude, longitude)

  # Calculate globe temperature
  if (verbose) cli_alert("Calculating globe temperature...")
  globe_temp <- calculate_globe_temp(
    temp, wind_speed, direct_solar, diffuse_solar,
    solar_pos$zenith, albedo, verbose = verbose
  )

  # Mean radiant temperature from globe temperature — ISO 7726 formula.
  # For a 150 mm black globe (emissivity 0.95):
  #   trad = ((Tg+273.15)^4 + 1.1e8 * V^0.6 / (emiss * D^0.4) * (Tg - ta))^0.25 - 273.15
  # where D = 0.15 m, emiss = 0.95, V = wind speed (m/s), ta = air temp (°C).
  trad <- ((globe_temp + 273.15)^4 +
    (1.1e8 / (0.95 * 0.15^0.4)) * wind_speed^0.6 * (globe_temp - temp))^0.25 - 273.15

  # Psychrometric calculations
  if (verbose) cli_alert("Calculating psychrometric variables...")
  es <- calc_sat_vp(temp)
  pa <- (RH / 100) * es
  temp_dewpoint <- calc_dew_point(temp, RH)

  # Natural wet bulb temperature
  if (is.null(wet_bulb)) {
    if (verbose) cli_alert("Calculating natural wet bulb temperature...")
    wet_bulb <- calculate_natural_wet_bulb(
      temp, RH, pressure, wind_speed, globe_temp,
      verbose = verbose, show_progress = FALSE
    )
  } else {
    if (verbose) cli_alert("Using supplied wet bulb temperature.")
  }

  # Solve for TWL
  if (verbose) cli_alert("Computing TWL...")

  pb_id <- NULL
  if (verbose && n_obs > 100) {
    pb_id <- cli_progress_bar("Computing TWL", total = n_obs)
  }

  TWL <- pmap_dbl(
    list(
      temp, wind_speed, RH, direct_solar, diffuse_solar,
      pressure, pa, temp_dewpoint,
      wet_bulb, globe_temp, trad, seq_len(n_obs)
    ),
    function(temp, wind_speed, RH, direct_solar, diffuse_solar,
             pressure, pa, temp_dewpoint,
             wet_bulb, globe_temp, trad, idx) {

      if (verbose && n_obs > 100) {
        cli_progress_update(id = pb_id, set = idx)
      }

      # Return NA if inputs are NA
      if (is.na(temp) || is.na(wind_speed) || is.na(RH) ||
          is.na(direct_solar) || is.na(diffuse_solar) || is.na(pressure)) {
        return(NA_real_)
      }

      # Solve for this point
      result <- solve_twl_single(
        temp, wind_speed, RH, direct_solar, diffuse_solar,
        pressure, pa, temp_dewpoint,
        wet_bulb, globe_temp, trad,
        max_core_temp, max_sweat_rate,
        Icl, icl, LR, lambda,
        index = if (verbose) idx else NULL
      )

      # Apply withdrawal limits
      if (!is.na(result)) {
        # DB > 44 degrees C withdrawal limit
        if (temp > 44) {
          result <- min(result, 115)
        }

        # WB > 32 degrees C withdrawal limit
        if (!is.na(wet_bulb) && wet_bulb > 32) {
          result <- min(result, 115)
        }

        # Upper-bound clamp only: the solver guarantees result >= TWL_FLOOR
        # (60 W/m^2), so the lower branch is unreachable and has been removed.
        if (verbose && result > 400) {
          cli_alert_warning(
            "Observation {idx}: TWL ({round(result, 1)} W/m\u00b2) above 400, constrained to ceiling"
          )
        }
        result <- min(result, TWL_CONSTANTS$TWL_CEILING)
      }

      result
    }
  )

  if (verbose && n_obs > 100) {
    cli_progress_done(id = pb_id)
  }

  # Summary statistics
  if (verbose) {
    n_valid <- sum(!is.na(TWL))
    n_invalid <- sum(is.na(TWL))

    if (n_valid > 0) {
      twl_range <- range(TWL, na.rm = TRUE)
      twl_mean <- mean(TWL, na.rm = TRUE)

      cli_alert_success(
        "TWL calculation complete: {n_valid} valid result{?s}, {n_invalid} invalid"
      )
      cli_alert_info(
        "TWL range: [{round(twl_range[1], 1)}, {round(twl_range[2], 1)}] W/m\u00b2, mean: {round(twl_mean, 1)} W/m\u00b2"
      )

      # Categorise results
      n_withdrawal <- sum(TWL < 115, na.rm = TRUE)
      n_buffer <- sum(TWL >= 115 & TWL < 140, na.rm = TRUE)
      n_acclimatisation <- sum(TWL >= 140 & TWL < 220, na.rm = TRUE)
      n_unrestricted <- sum(TWL >= 220, na.rm = TRUE)

      cli_h2("TWL Categories")
      cli_ul(c(
        "Withdrawal (<115 W/m\u00b2): {n_withdrawal} ({round(100*n_withdrawal/n_valid, 1)}%)",
        "Buffer (115-140 W/m\u00b2): {n_buffer} ({round(100*n_buffer/n_valid, 1)}%)",
        "Acclimatisation (140-220 W/m\u00b2): {n_acclimatisation} ({round(100*n_acclimatisation/n_valid, 1)}%)",
        "Unrestricted (>=220 W/m\u00b2): {n_unrestricted} ({round(100*n_unrestricted/n_valid, 1)}%)"
      ))

      if (n_withdrawal > 0) {
        cli_alert_warning(
          "{n_withdrawal} observation{?s} in withdrawal zone"
        )
      }
    } else {
      cli_alert_danger("No valid TWL results calculated")
    }
  }

  return(TWL)
}
