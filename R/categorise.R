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

# NOTE (#11): the litter/dust/odour 0-100 operational tier helpers
# (categorise_litter/dust/odour and their *_colour partners) were removed.
# Those hazards now expose physical quantities only (odour: relative
# concentration from odour_exposure(); dust: dust_flux() / crust-adjusted
# relative flux from dust_hazard(); litter: the relative index from
# litter_hazard()). Mapping a physical value onto a site-specific operational
# index and tiers is a calibration step delivered by forthcoming calibration
# tooling, not by fixed package cut-points. TWL keeps its tiers above because
# its W/m^2 zones are physically grounded (Brake & Bates 2002), not a guessed
# 0-100 scale.
