#' Categorise TWL values into work restriction zones
#'
#' Maps TWL values to the four standard zones defined by Brake & Bates (2002).
#'
#' @param twl TWL values in W/m^2, either a bare numeric or a \pkg{units} object
#'   (e.g. the [generate_twl()] return value); a tagged value is converted.
#' @return Character vector of categories: `"Withdrawal"`, `"Buffer"`,
#'   `"Acclimatisation"`, or `"Unrestricted"`.
#'
#' @examples
#' categorise_twl(c(100, 130, 180, 250))
#'
#' @export
categorise_twl <- function(twl) {
  # Accept either a bare numeric (assumed W/m^2) or a units object from
  # generate_twl() (converted); the band logic runs on plain W/m^2 doubles.
  twl <- .drop_to(twl, "W/m^2", arg = "twl")
  case_when(
    is.na(twl) ~ NA_character_,
    twl < 115 ~ "Withdrawal",
    twl < 140 ~ "Buffer",
    twl < 220 ~ "Acclimatisation",
    TRUE ~ "Unrestricted"
  )
}

#' Get TWL category colour
#'
#' Returns hex colour codes corresponding to each TWL zone for use in
#' visualisations.
#'
#' @param twl TWL values in W/m^2, either a bare numeric or a \pkg{units} object
#'   (e.g. the [generate_twl()] return value); a tagged value is converted.
#' @return Character vector of hex colour codes:
#' \describe{
#'   \item{`"#D32F2F"`}{Red -- Withdrawal}
#'   \item{`"#FF9800"`}{Orange -- Buffer}
#'   \item{`"#FFC107"`}{Amber -- Acclimatisation}
#'   \item{`"#4CAF50"`}{Green -- Unrestricted}
#'   \item{`"#CCCCCC"`}{Grey -- NA}
#' }
#'
#' @examples
#' twl_colour(c(100, 130, 180, 250, NA))
#'
#' @export
twl_colour <- function(twl) {
  # Accept a bare numeric (W/m^2) or a units object from generate_twl().
  twl <- .drop_to(twl, "W/m^2", arg = "twl")
  case_when(
    is.na(twl) ~ "#CCCCCC",
    twl < 115 ~ "#D32F2F",
    twl < 140 ~ "#FF9800",
    twl < 220 ~ "#FFC107",
    TRUE ~ "#4CAF50"
  )
}

# The 0-100 hazard tiers below share a common four-band palette (green ->
# amber -> orange -> red, grey for NA), the same scheme as the TWL zones above,
# so a dashboard can colour every hazard consistently. The band cut-points are
# provisional and tracked for calibration (see the odour calibration issue and
# the litter/dust tier notes in their specs).

#' Categorise a litter hazard index into operational tiers
#'
#' Maps the 0-100 litter hazard ([litter_hazard()]) to the four operational
#' tiers documented in `specs/Litter_v3.md`. The cut-points (20/45/70) are
#' pre-calibration estimates.
#'
#' @param hazard Numeric vector of litter hazard index values in `[0, 100]`.
#' @return Character vector of tiers: `"LOW"`, `"MODERATE"`, `"HIGH"`,
#'   `"EXTREME"`, or `NA`.
#'
#' @examples
#' categorise_litter(c(10, 30, 55, 85))
#'
#' @seealso [litter_hazard()], [litter_colour()].
#' @export
categorise_litter <- function(hazard) {
  case_when(
    is.na(hazard) ~ NA_character_,
    hazard < 20 ~ "LOW",
    hazard < 45 ~ "MODERATE",
    hazard < 70 ~ "HIGH",
    TRUE ~ "EXTREME"
  )
}

#' Get litter hazard tier colour
#'
#' Returns hex colour codes for the litter hazard tiers (see
#' [categorise_litter()]), on the package's shared green/amber/orange/red
#' hazard palette.
#'
#' @param hazard Numeric vector of litter hazard index values in `[0, 100]`.
#' @return Character vector of hex colour codes (`"#4CAF50"` LOW, `"#FFC107"`
#'   MODERATE, `"#FF9800"` HIGH, `"#D32F2F"` EXTREME, `"#CCCCCC"` NA).
#'
#' @examples
#' litter_colour(c(10, 30, 55, 85, NA))
#'
#' @seealso [categorise_litter()].
#' @export
litter_colour <- function(hazard) {
  case_when(
    is.na(hazard) ~ "#CCCCCC",
    hazard < 20 ~ "#4CAF50",
    hazard < 45 ~ "#FFC107",
    hazard < 70 ~ "#FF9800",
    TRUE ~ "#D32F2F"
  )
}

#' Categorise a dust hazard index into operational tiers
#'
#' Maps the 0-100 dust hazard ([dust_hazard()]) to the four operational tiers
#' documented in `specs/Dust.md`. The cut-points (25/50/75) are pre-calibration
#' estimates.
#'
#' @param hazard Numeric vector of dust hazard index values in `[0, 100]`.
#' @return Character vector of tiers: `"LOW"`, `"MODERATE"`, `"HIGH"`,
#'   `"EXTREME"`, or `NA`.
#'
#' @examples
#' categorise_dust(c(10, 35, 60, 90))
#'
#' @seealso [dust_hazard()], [dust_colour()].
#' @export
categorise_dust <- function(hazard) {
  case_when(
    is.na(hazard) ~ NA_character_,
    hazard < 25 ~ "LOW",
    hazard < 50 ~ "MODERATE",
    hazard < 75 ~ "HIGH",
    TRUE ~ "EXTREME"
  )
}

#' Get dust hazard tier colour
#'
#' Returns hex colour codes for the dust hazard tiers (see [categorise_dust()]),
#' on the package's shared green/amber/orange/red hazard palette.
#'
#' @param hazard Numeric vector of dust hazard index values in `[0, 100]`.
#' @return Character vector of hex colour codes (`"#4CAF50"` LOW, `"#FFC107"`
#'   MODERATE, `"#FF9800"` HIGH, `"#D32F2F"` EXTREME, `"#CCCCCC"` NA).
#'
#' @examples
#' dust_colour(c(10, 35, 60, 90, NA))
#'
#' @seealso [categorise_dust()].
#' @export
dust_colour <- function(hazard) {
  case_when(
    is.na(hazard) ~ "#CCCCCC",
    hazard < 25 ~ "#4CAF50",
    hazard < 50 ~ "#FFC107",
    hazard < 75 ~ "#FF9800",
    TRUE ~ "#D32F2F"
  )
}

#' Categorise an odour exposure into operational tiers
#'
#' Maps the 0-100 odour exposure ([odour_exposure()] / [odour_risk()]) to the
#' four operational tiers documented in the README. The cut-points (15/40/70)
#' are provisional and tunable through `map_c50`; see the odour calibration
#' issue.
#'
#' @param exposure Numeric vector of odour exposure values in `[0, 100]`.
#' @return Character vector of tiers: `"LOW"`, `"MODERATE"`, `"HIGH"`,
#'   `"VERY HIGH"`, or `NA`.
#'
#' @examples
#' categorise_odour(c(10, 30, 55, 85))
#'
#' @seealso [odour_exposure()], [odour_risk()], [odour_colour()].
#' @export
categorise_odour <- function(exposure) {
  case_when(
    is.na(exposure) ~ NA_character_,
    exposure < 15 ~ "LOW",
    exposure < 40 ~ "MODERATE",
    exposure < 70 ~ "HIGH",
    TRUE ~ "VERY HIGH"
  )
}

#' Get odour exposure tier colour
#'
#' Returns hex colour codes for the odour exposure tiers (see
#' [categorise_odour()]), on the package's shared green/amber/orange/red hazard
#' palette.
#'
#' @param exposure Numeric vector of odour exposure values in `[0, 100]`.
#' @return Character vector of hex colour codes (`"#4CAF50"` LOW, `"#FFC107"`
#'   MODERATE, `"#FF9800"` HIGH, `"#D32F2F"` VERY HIGH, `"#CCCCCC"` NA).
#'
#' @examples
#' odour_colour(c(10, 30, 55, 85, NA))
#'
#' @seealso [categorise_odour()].
#' @export
odour_colour <- function(exposure) {
  case_when(
    is.na(exposure) ~ "#CCCCCC",
    exposure < 15 ~ "#4CAF50",
    exposure < 40 ~ "#FFC107",
    exposure < 70 ~ "#FF9800",
    TRUE ~ "#D32F2F"
  )
}
