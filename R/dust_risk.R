#' Dust emission potential for exposed landfill surfaces
#'
#' Computes an hourly dust emission potential (g/m²) for an exposed landfill
#' surface using the US EPA AP-42 Section 13.2.5 industrial wind erosion
#' methodology, with a dynamic soil moisture correction replacing AP-42's
#' static climatological moisture factor.
#'
#' The function is designed to operate on a single forecast hour, accepting
#' the meteorological variables provided by the Open-Meteo forecast API.
#' It is intended as the core engine of an operational dust warning system
#' for landfill sites regulated under the NSW EPA framework.
#'
#' @section Calculation overview:
#' The computation proceeds in six steps:
#' \enumerate{
#'   \item Convert Tyler Sieve number to aggregate diameter (m).
#'   \item Compute dry threshold friction velocity via Shao & Lu (2000).
#'   \item Adjust threshold for soil moisture via Fécan et al. (1999).
#'   \item Derive effective friction velocity from hourly wind data.
#'   \item Compute AP-42 erosion potential P (g/m²).
#'   \item Scale by particle-size multiplier k to give E (g/m²).
#' }
#'
#' @param tyler_sieve_no Integer. Tyler Standard Sieve number representing the
#'   mode of the surface aggregate size distribution of the most erodible
#'   exposed surface. Collect a surface sample from the active cell or
#'   intermediate cover, dry-sieve using Tyler sieves, and use the sieve number
#'   corresponding to the modal size class. Valid values:
#'   3, 4, 5, 6, 8, 9, 10, 14, 20, 24, 28, 32, 35, 42, 48, 60, 65, 80, 100,
#'   115, 150, 170, 200, 250, 270, 325, 400.
#'
#' @param clay_percent Numeric. Clay fraction of the surface material (< 2 µm
#'   particles), expressed as percent by mass (e.g., \code{15} for 15% clay).
#'
#' @param wind_speed_10m Numeric. Hourly mean wind speed at 10 m (km/h).
#'
#' @param wind_gusts_10m Numeric. Maximum 3-second gust in the preceding hour
#'   at 10 m height (km/h).
#'
#' @param soil_moisture Numeric. Volumetric water content of the top 0–1 cm
#'   soil layer (m³/m³).
#'
#' @param z0 Numeric. Aerodynamic roughness length (m). Default: 0.005 m.
#'
#' @param bulk_density Numeric. Dry bulk density (Mg/m³). Default: 1.6 Mg/m³.
#'
#' @param k Numeric. AP-42 particle size multiplier. Default: 0.5 (PM₁₀).
#'
#' @return Numeric scalar. Estimated dust emission potential E (g/m²). Returns
#'   0 when conditions do not exceed the erosion threshold.
#'
#' @references
#' US EPA (2006) AP-42 Section 13.2.5: Industrial Wind Erosion.
#' Shao, Y. and Lu, H. (2000). doi:10.1029/2000JD900304.
#' Fécan, F. et al. (1999). doi:10.1007/s00585-999-0149-7.
#'
#' @export
dust_emission_potential <- function(
  tyler_sieve_no,
  clay_percent,
  wind_speed_10m,
  wind_gusts_10m,
  soil_moisture,
  z0 = 0.005,
  bulk_density = 1.6,
  k = 0.5
) {
  # ---- Input validation --------------------------------------------------- #

  valid_sieves <- as.integer(names(TYLER_SIEVE_DIAMETERS_M))

  checkmate::assert_int(tyler_sieve_no)
  if (!as.character(tyler_sieve_no) %in% names(TYLER_SIEVE_DIAMETERS_M)) {
    cli::cli_abort(
      c(
        "Invalid {.arg tyler_sieve_no}: {tyler_sieve_no}.",
        "i" = "Valid Tyler Sieve numbers are: {valid_sieves}."
      )
    )
  }

  checkmate::assert_number(clay_percent, lower = 0, upper = 100)
  checkmate::assert_number(wind_speed_10m, lower = 0)
  checkmate::assert_number(wind_gusts_10m, lower = 0)
  checkmate::assert_number(soil_moisture, lower = 0, upper = 1)
  checkmate::assert_number(z0, lower = 1e-6)
  checkmate::assert_number(bulk_density, lower = 0)
  checkmate::assert_number(k, lower = 0)

  # ---- Physical constants ------------------------------------------------- #

  A_N   <- 0.0123
  rho_p <- 2650.0
  rho_a <- 1.225
  g     <- 9.81
  gamma <- 3.0e-4
  kappa <- 0.40
  z     <- 10.0

  # ---- Step 1: Tyler Sieve number → aggregate diameter (m) --------------- #

  d <- TYLER_SIEVE_DIAMETERS_M[[as.character(tyler_sieve_no)]]

  # ---- Step 2: Dry threshold friction velocity (m/s) --------------------- #

  u_star_t_dry <- sqrt(A_N * (rho_p / rho_a * g * d + gamma / (rho_a * d)))

  # ---- Step 3: Moisture-adjusted threshold friction velocity (m/s) -------- #

  w_prime <- 0.0014 * clay_percent^2 + 0.17 * clay_percent
  w <- (soil_moisture / bulk_density) * 100

  if (w > w_prime) {
    u_star_t <- u_star_t_dry * sqrt(1 + 1.21 * (w - w_prime)^0.68)
  } else {
    u_star_t <- u_star_t_dry
  }

  # ---- Step 4: Effective friction velocity from forecast wind data (m/s) -- #

  U_10  <- wind_speed_10m / 3.6
  U_gust <- wind_gusts_10m / 3.6
  U_fm  <- max(U_10, U_gust * 0.70)
  u_star <- (kappa * U_fm) / log(z / z0)

  # ---- Step 5: Erosion potential P (g/m²) -------------------------------- #

  excess <- u_star - u_star_t

  if (excess > 0) {
    P <- 58 * excess^2 + 25 * excess
  } else {
    P <- 0.0
  }

  # ---- Step 6: Apply particle-size multiplier ----------------------------- #

  E <- k * P

  return(E)
}


#' Dust Risk Index for landfill operations
#'
#' Computes an hourly Dust Risk Index in the range \[0, 100\] for each timestep
#' in a meteorological forecast tibble. Wraps \code{dust_emission_potential()}
#' and normalises the raw emission potential (g/m²) to a bounded 0–100 index
#' suitable for operational strip plots and trend dashboards.
#'
#' @param met_data A tibble (or data frame) with one row per hourly timestep,
#'   containing at minimum:
#'   \describe{
#'     \item{\code{wind_speed_10m}}{Mean wind speed at 10 m (km/h).}
#'     \item{\code{wind_gusts_10m}}{Peak gust at 10 m (km/h).}
#'     \item{\code{soil_moisture_0_to_1cm}}{Volumetric soil moisture (m³/m³).}
#'   }
#'
#' @param tyler_sieve_no Integer. Tyler sieve number for the modal aggregate
#'   size of the erodible surface. Default: \code{20L} (0.833 mm, coarse sand).
#'
#' @param clay_percent Numeric. Clay fraction (% by mass). Default: \code{10}.
#'
#' @param z0 Numeric. Aerodynamic roughness length (m). Default: \code{0.005}.
#'
#' @param bulk_density Numeric. Dry bulk density (Mg/m³). Default: \code{1.6}.
#'
#' @param k Numeric. AP-42 particle-size multiplier. Default: \code{0.5} (PM₁₀).
#'
#' @param scale_max Numeric. Emission potential (g/m²) mapped to index = 100.
#'   Default: \code{2.0} g/m². Values at or above this are capped at 100.
#'   Tier interpretation (default scale):
#'   \tabular{rll}{
#'     \strong{Index} \tab \strong{Tier} \tab \strong{E (g/m²)} \cr
#'     0–24   \tab LOW      \tab 0.0–0.5 \cr
#'     25–49  \tab MODERATE \tab 0.5–1.0 \cr
#'     50–74  \tab HIGH     \tab 1.0–1.5 \cr
#'     75–100 \tab EXTREME  \tab 1.5–2.0+
#'   }
#'
#' @return Numeric vector of length \code{nrow(met_data)} giving the Dust
#'   Risk Index in \eqn{[0, 100]} for each forecast hour.
#'
#' @export
generate_dust_risk_index <- function(
  met_data,
  tyler_sieve_no = 20L,
  clay_percent   = 10,
  z0             = 0.005,
  bulk_density   = 1.6,
  k              = 0.5,
  scale_max      = 2.0
) {
  checkmate::assert_data_frame(met_data, min.rows = 1)

  required_cols <- c("wind_speed_10m", "wind_gusts_10m", "soil_moisture_0_to_1cm")
  missing_cols  <- setdiff(required_cols, names(met_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "{.arg met_data} is missing required columns: {.val {missing_cols}}.",
      "i" = "Required: wind_speed_10m (km/h), wind_gusts_10m (km/h), soil_moisture_0_to_1cm (m\u00b3/m\u00b3)."
    ))
  }

  emissions <- mapply(
    FUN = dust_emission_potential,
    wind_speed_10m = met_data$wind_speed_10m,
    wind_gusts_10m = met_data$wind_gusts_10m,
    soil_moisture  = met_data$soil_moisture_0_to_1cm,
    MoreArgs = list(
      tyler_sieve_no = tyler_sieve_no,
      clay_percent   = clay_percent,
      z0             = z0,
      bulk_density   = bulk_density,
      k              = k
    ),
    SIMPLIFY = TRUE
  )

  pmin(emissions / scale_max * 100, 100)
}


# Tyler Standard Sieve series: sieve number → nominal opening in metres.
#
# Source: Tyler Industrial Products sieve series, as referenced in
# US EPA AP-42 Section 13.2.5, Table 13.2.5-1.
TYLER_SIEVE_DIAMETERS_M <- c(
  "3" = 0.006680,
  "4" = 0.004699,
  "5" = 0.003962,
  "6" = 0.003327,
  "8" = 0.002362,
  "9" = 0.001981,
  "10" = 0.001651,
  "14" = 0.001168,
  "20" = 0.000833,
  "24" = 0.000701,
  "28" = 0.000589,
  "32" = 0.000495,
  "35" = 0.000417,
  "42" = 0.000351,
  "48" = 0.000295,
  "60" = 0.000246,
  "65" = 0.000208,
  "80" = 0.000175,
  "100" = 0.000147,
  "115" = 0.000124,
  "150" = 0.000104,
  "170" = 0.000088,
  "200" = 0.000074,
  "250" = 0.000063,
  "270" = 0.000053,
  "325" = 0.000044,
  "400" = 0.000037
)
