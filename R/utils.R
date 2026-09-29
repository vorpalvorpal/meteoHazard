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


# ---- Missing-input (NA) policy --------------------------------------------- #
# Shared by dust_hazard(), litter_risk(), odour_risk(), generate_twl() and the
# layers under them: a row whose required inputs include an NA returns NA for
# that row only (the other rows are computed as normal), with ONE classed
# summary warning per call (`meteoHazard_missing_input`). Genuinely invalid,
# non-missing values (negative wind, RH > 100, ...) still error.

# Logical vector of length `n`: TRUE where any of `...` (each length 1 or `n`,
# NULL ignored) is NA at that row.
.missing_rows <- function(n, ...) {
  miss <- logical(n)
  for (x in list(...)) {
    if (is.null(x)) next
    miss <- miss | rep_len(is.na(x), n)
  }
  miss
}

# As .missing_rows(), over the `cols` of data frame `data` (absent columns are
# ignored; presence is validated by the caller).
.missing_cols_rows <- function(data, cols) {
  do.call(.missing_rows, c(list(nrow(data)), unname(lapply(cols, function(col) data[[col]]))))
}

# Emit the single summary warning for rows returned as NA because of missing
# inputs. Silent when nothing is missing. `fn` names the calling function.
.warn_missing_rows <- function(miss, fn) {
  rows <- which(miss)
  if (length(rows) == 0L) return(invisible(FALSE))
  shown <- rows[seq_len(min(10L, length(rows)))]
  more  <- if (length(rows) > length(shown)) {
    sprintf(" (and %d more)", length(rows) - length(shown))
  } else {
    ""
  }
  one <- length(rows) == 1L
  head_msg <- sprintf(
    "%d of %d row%s %s a missing input; returned NA for %s.",
    length(rows), length(miss), if (length(miss) == 1L) "" else "s",
    if (one) "has" else "have", if (one) "that row" else "those rows"
  )
  cli::cli_warn(
    c(
      paste0("{.fn {fn}}: ", head_msg),
      "i" = paste0("Row", if (one) "" else "s", ": ",
                   paste(shown, collapse = ", "), more, "."),
      "i" = "Warning class {.cls meteoHazard_missing_input} (muffle it by class to silence)."
    ),
    class = "meteoHazard_missing_input"
  )
  invisible(TRUE)
}

# Evaluate `expr` with inner missing-input warnings muffled, so a wrapper can
# emit one summary warning of its own instead of one per layer.
.muffle_missing_input <- function(expr) {
  withCallingHandlers(
    expr,
    meteoHazard_missing_input = function(w) invokeRestart("muffleWarning")
  )
}
