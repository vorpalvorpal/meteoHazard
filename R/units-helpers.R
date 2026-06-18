# ===========================================================================
# Shared units helpers (built on the `units` package).
#
# The package's public functions accept dimensional inputs either as bare
# numerics (assumed to already be in the documented unit) or as `units` objects
# (which are converted, erroring on a dimensional mismatch). These helpers
# centralise that contract so a wrong-unit input can never be silently misread
# -- the class of bug behind the historical km/h-vs-m/s litter/dust error.
#
# Internally the dimensional inputs are normalised to a canonical unit; the
# gate / case_when / coalesce logic then runs on plain doubles (units objects
# cannot be compared against bare literals, and ifelse() strips units), while
# genuinely physical OUTPUTS (e.g. the TWL in W/m^2, temperatures in degC) are
# re-tagged as `units` objects on the way out. Dimensionless scores (the 0-100
# indices, the relative odour hazard) and percentage / ratio fields (cloud
# cover, relative humidity, soil moisture) stay plain numeric throughout.
# ===========================================================================

# Coerce `x` to a `units` object in `unit`. A bare numeric is assumed to already
# be in `unit` (and is tagged with it); an existing `units` object is converted
# to `unit`, erroring if the dimensions are incompatible. NULL passes through.
.as_units <- function(x, unit, arg = deparse(substitute(x))) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be numeric (bare or a {.cls units} object), not {.cls {class(x)}}.",
      class = "meteoHazard_input_error"
    )
  }
  tryCatch(
    units::set_units(x, unit, mode = "standard"),
    error = function(e) {
      cli::cli_abort(
        c("{.arg {arg}} carries units incompatible with the expected {.val {unit}}.",
          "x" = conditionMessage(e)),
        class = "meteoHazard_input_error"
      )
    }
  )
}

# As .as_units(), but return a plain double in `unit` (units dropped) for the
# internal arithmetic and gate logic. Bare numerics pass through their values
# unchanged (assumed already in `unit`); units objects are converted first.
.drop_to <- function(x, unit, arg = deparse(substitute(x))) {
  u <- .as_units(x, unit, arg = arg)
  if (is.null(u)) return(NULL)
  units::drop_units(u)
}
