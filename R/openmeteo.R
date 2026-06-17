#' Build an Open-Meteo httr2 request (without performing it)
#'
#' Constructs an httr2 request object for the given URL with the standard
#' meteoHazard user-agent, a 30-second timeout, and up to three retries with
#' a fixed 2-second back-off. The request is NOT performed here; call
#' [om_perform()] to execute it.
#'
#' @param url Character scalar. The fully-constructed Open-Meteo URL.
#' @return An `httr2_request` object.
#' @keywords internal
om_request <- function(url) {
  httr2::request(url) |>
    httr2::req_user_agent(
      "meteoHazard R package (https://github.com/vorpalvorpal/meteoHazard)"
    ) |>
    httr2::req_timeout(30) |>
    httr2::req_retry(max_tries = 3, backoff = ~2)
}

#' Perform an Open-Meteo httr2 request and return the parsed body
#'
#' Executes the request and deserialises the JSON response body. Raises a
#' non-call error if the HTTP request fails or if the API returns an error
#' payload (`{"error": true, "reason": "..."}`).
#'
#' @param req An `httr2_request` object produced by [om_request()].
#' @return A list (the parsed JSON body), with at least `$hourly$time` and
#'   `$hourly[[field]]` elements for a successful variable request.
#' @keywords internal
om_perform <- function(req) {
  resp <- tryCatch(
    {
      raw <- httr2::req_perform(req)
      httr2::resp_body_json(raw, simplifyVector = TRUE)
    },
    error = function(e) {
      stop(
        "Failed to fetch data from Open-Meteo: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (isTRUE(resp$error)) {
    stop("Open-Meteo API error: ", resp$reason, call. = FALSE)
  }

  resp
}

#' Fetch weather data from the Open-Meteo API
#'
#' Queries the Open-Meteo forecast or historical archive API to retrieve
#' hourly weather observations for the given location and times. Only the
#' variables listed in `fields` are requested.
#'
#' Input `datetime` values must be in UTC (or a UTC-equivalent timezone).
#' Sub-hourly timestamps are truncated to the hour before matching against
#' API data.
#'
#' @param datetime POSIXct datetime vector (UTC).
#' @param latitude Latitude in decimal degrees.
#' @param longitude Longitude in decimal degrees.
#' @param fields Character vector of Open-Meteo hourly variable names to
#'   request (e.g. `"temperature_2m"`, `"wind_speed_10m"`).
#' @param verbose If `TRUE`, report progress.
#'
#' @return A named list of numeric vectors, one per requested field, aligned
#'   to the input `datetime` values. Values that cannot be matched are `NA`.
#' @keywords internal
fetch_openmeteo <- function(
  datetime,
  latitude,
  longitude,
  fields,
  verbose = FALSE
) {
  if (length(fields) == 0L) {
    return(list())
  }

  # Always work in UTC to avoid session-timezone mismatches
  datetime_utc <- as.POSIXct(format(datetime, tz = "UTC"), tz = "UTC")
  dates <- as.Date(datetime_utc, tz = "UTC")
  date_min <- min(dates, na.rm = TRUE)
  date_max <- max(dates, na.rm = TRUE)
  today <- Sys.Date()

  fields_csv <- paste(fields, collapse = ",")

  # Boundary between the forecast and archive endpoints.
  # Forecast covers [today-92, today+16]; archive covers everything older.
  B <- today - 92L

  if (verbose) {
    cli_alert("Fetching weather data from Open-Meteo...")
  }

  # ---------------------------------------------------------------------------
  # Endpoint selection via capability flags
  # ---------------------------------------------------------------------------
  # need_archive  : date range extends before the boundary B (archive-only data)
  # need_forecast : date range reaches B or beyond (forecast endpoint needed)
  #
  # Four combinations:
  #   (F=F, A=T) Entirely before B        -> single archive request
  #   (F=T, A=F) Entirely within forecast -> single forecast request
  #   (F=T, A=T) Straddles boundary       -> archive [date_min, B) +
  #                                          forecast [B, date_max]; merge with
  #                                          forecast first so it wins on overlap
  #   (F=F, A=F) Cannot happen (non-empty date range); guard returns empty list
  # ---------------------------------------------------------------------------

  need_archive <- date_min < B
  need_forecast <- date_max >= B

  arch <- NULL
  fc <- NULL

  if (need_archive) {
    # Archive sub-request.
    # end_date = min(date_max, B - 1) to avoid requesting into the forecast
    # window; because need_archive is TRUE we know date_min < B, so
    # start_date < end_date always holds.
    arch_end <- min(date_max, B - 1L)
    arch_url <- sprintf(
      "https://archive-api.open-meteo.com/v1/archive?latitude=%.6f&longitude=%.6f&hourly=%s&wind_speed_unit=ms&start_date=%s&end_date=%s&timeformat=iso8601&timezone=UTC",
      latitude, longitude, fields_csv,
      format(date_min, "%Y-%m-%d"),
      format(arch_end, "%Y-%m-%d")
    )
    arch <- om_perform(om_request(arch_url))
  }

  if (need_forecast) {
    # Forecast sub-request.
    # past_days reaches back to max(date_min, B) so we cover the full range
    # requested without asking for archive data the forecast endpoint cannot
    # provide.
    past <- max(0L, min(92L, as.integer(today - max(date_min, B))))
    fore <- max(1L, min(16L, as.integer(date_max - today) + 1L))
    fc_url <- sprintf(
      "https://api.open-meteo.com/v1/forecast?latitude=%.6f&longitude=%.6f&hourly=%s&wind_speed_unit=ms&past_days=%d&forecast_days=%d&timeformat=iso8601&timezone=UTC",
      latitude, longitude, fields_csv, past, fore
    )
    fc <- om_perform(om_request(fc_url))
  }

  if (need_forecast && need_archive) {
    # Both requests ran: merge with FORECAST first so its values win on any
    # overlapping boundary hour.
    api_times_raw <- c(fc$hourly$time, arch$hourly$time)

    n_fc <- length(fc$hourly$time)
    n_arch <- length(arch$hourly$time)

    hourly_data <- list(time = api_times_raw)
    for (f in fields) {
      fc_vals <- fc$hourly[[f]]
      arch_vals <- arch$hourly[[f]]
      # Guard against a field missing from one response (e.g. a gap in archive).
      if (is.null(fc_vals)) fc_vals <- rep(NA_real_, n_fc)
      if (is.null(arch_vals)) arch_vals <- rep(NA_real_, n_arch)
      hourly_data[[f]] <- c(fc_vals, arch_vals)
    }
  } else if (need_forecast) {
    api_times_raw <- fc$hourly$time
    hourly_data <- fc$hourly
  } else if (need_archive) {
    api_times_raw <- arch$hourly$time
    hourly_data <- arch$hourly
  } else {
    # Unreachable with a non-empty date range, but guard defensively.
    return(stats::setNames(
      lapply(fields, function(f) rep(NA_real_, length(datetime))),
      fields
    ))
  }

  # ---------------------------------------------------------------------------
  # Align API timestamps to input datetimes (truncated to the hour)
  # ---------------------------------------------------------------------------
  api_times <- as.POSIXct(
    api_times_raw,
    format = "%Y-%m-%dT%H:%M",
    tz = "UTC"
  )
  input_hours <- as.POSIXct(
    format(datetime_utc, "%Y-%m-%d %H:00:00", tz = "UTC"),
    tz = "UTC"
  )
  idx <- match(input_hours, api_times)

  n_matched <- sum(!is.na(idx))
  n_total <- length(datetime)
  n_unmatched <- sum(is.na(idx))

  if (verbose) {
    cli_alert_success(
      "Matched {n_matched}/{n_total} observation{?s} to API data"
    )
  }

  if (n_unmatched > 0L) {
    cli_alert_warning(
      "{n_unmatched} observation{?s} could not be matched to API data and will be NA"
    )
  }

  result <- lapply(fields, function(f) {
    vals <- hourly_data[[f]]
    if (is.null(vals)) {
      rep(NA_real_, n_total)
    } else {
      raw_vals <- as.numeric(vals[idx])
      n_na <- sum(is.na(raw_vals))
      if (n_na > 0L && verbose) {
        cli_alert_warning(
          "Field '{f}': {n_na} matched value{?s} are NA (may be at archive boundary)"
        )
      }
      raw_vals
    }
  })
  names(result) <- fields

  result
}
