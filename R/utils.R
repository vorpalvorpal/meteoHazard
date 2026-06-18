# ===========================================================================
# Shared internal utilities used across the hazard / exposure functions.
# Single-sourcing the input-validation idioms and the NA-coalescing helper so
# the four hazard families validate identically and cannot silently diverge.
# ===========================================================================

# Abort (classed) if `data` is missing any of `required` columns. `arg` is the
# argument name used in the message. Returns invisibly on success.
.assert_required_cols <- function(data, required, arg = "met_data", info = NULL) {
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    msg <- c("{.arg {arg}} is missing required columns: {.val {missing_cols}}.")
    if (!is.null(info)) msg <- c(msg, "i" = info)
    cli::cli_abort(msg, class = "meteoHazard_input_error")
  }
  invisible(data)
}

# Abort (classed) if any of `cols` in `data` is non-numeric. Assumes the columns
# are present (run .assert_required_cols first). Returns invisibly on success.
.assert_numeric_cols <- function(data, cols, arg = "met_data") {
  for (col in cols) {
    if (!is.numeric(data[[col]])) {
      cli::cli_abort(
        "{.arg {arg}} column {.val {col}} must be numeric, not {.cls {class(data[[col]])}}.",
        class = "meteoHazard_input_error"
      )
    }
  }
  invisible(data)
}

# Replace NA in `x` with `default` (scalar). Vectorised; preserves length.
.na_fill <- function(x, default) {
  ifelse(is.na(x), default, x)
}
